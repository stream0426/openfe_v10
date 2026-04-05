#!/usr/bin/env python3
"""Build an OpenFE RBFE alchemical network from an external YAML config."""

from __future__ import annotations

import argparse
import functools
import json
import os
import pathlib
import shutil
from collections.abc import Iterable, Mapping
from typing import Any

import yaml


def _ensure_writable_cache_dirs(root: pathlib.Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    cache_dirs = {
        "MPLCONFIGDIR": root / "matplotlib",
        "NUMBA_CACHE_DIR": root / "numba",
        "XDG_CACHE_HOME": root / "xdg",
    }
    for env_var, path in cache_dirs.items():
        path.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault(env_var, str(path))


def _deep_merge(base: dict[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = _deep_merge(dict(merged[key]), value)
        else:
            merged[key] = value
    return merged


def _load_yaml(path: pathlib.Path) -> dict[str, Any]:
    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise TypeError(f"Expected top-level mapping in {path}, got {type(data).__name__}")
    return data


def _normalize_protocol_config(protocol_config: Mapping[str, Any]) -> dict[str, Any]:
    normalized = _deep_merge({}, protocol_config)

    lower_string_paths = [
        ("solvation_settings", "solvent_model"),
        ("solvation_settings", "box_shape"),
        ("partial_charge_settings", "partial_charge_method"),
        ("partial_charge_settings", "off_toolkit_backend"),
        ("simulation_settings", "sampler_method"),
        ("simulation_settings", "sams_flatness_criteria"),
        ("engine_settings", "compute_platform"),
        ("alchemical_settings", "softcore_LJ"),
    ]
    upper_string_paths = [
        ("forcefield_settings", "nonbonded_method"),
    ]

    for section, key in lower_string_paths:
        section_map = normalized.get(section)
        if isinstance(section_map, Mapping) and isinstance(section_map.get(key), str):
            normalized[section][key] = section_map[key].lower()

    for section, key in upper_string_paths:
        section_map = normalized.get(section)
        if isinstance(section_map, Mapping) and isinstance(section_map.get(key), str):
            normalized[section][key] = section_map[key].upper()

    ligand_forcefield = normalized.pop("ligand_forcefield", None)
    if ligand_forcefield is not None:
        ff_settings = dict(normalized.get("forcefield_settings", {}))
        ff_settings["small_molecule_forcefield"] = ligand_forcefield
        normalized["forcefield_settings"] = ff_settings

    return normalized


def _water_forcefield_files(solvent_model: str) -> list[str]:
    water_ff = {
        "tip3p": [
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml",
        ],
        "spce": [
            "amber/spce_standard.xml",
            "amber/spce_HFE_multivalent.xml",
        ],
        "tip4pew": [
            "amber/tip4pew_standard.xml",
            "amber/tip4pew_HFE_multivalent.xml",
        ],
    }
    try:
        return water_ff[solvent_model]
    except KeyError as exc:
        raise ValueError(
            "Automatic forcefield synchronization is only supported for "
            f"tip3p, spce, and tip4pew. Got solvent_model={solvent_model!r}. "
            "For other models, set protocol.forcefield_settings.forcefields explicitly."
        ) from exc


def _sync_forcefield_settings_with_solvent(protocol_config: dict[str, Any]) -> dict[str, Any]:
    solvation_settings = protocol_config.get("solvation_settings")
    if not isinstance(solvation_settings, Mapping):
        return protocol_config

    solvent_model = solvation_settings.get("solvent_model")
    if not isinstance(solvent_model, str):
        return protocol_config

    forcefield_settings = dict(protocol_config.get("forcefield_settings", {}))
    if "forcefields" in forcefield_settings:
        return protocol_config

    current_forcefields = forcefield_settings.get(
        "forcefields",
        [
            "amber/ff14SB.xml",
            "amber/tip3p_standard.xml",
            "amber/tip3p_HFE_multivalent.xml",
            "amber/phosaa10.xml",
            "amber/lipid17_merged.xml",
        ],
    )
    preserved = [
        ff for ff in current_forcefields
        if "tip3p_standard.xml" not in ff
        and "tip3p_HFE_multivalent.xml" not in ff
        and "spce_standard.xml" not in ff
        and "spce_HFE_multivalent.xml" not in ff
        and "tip4pew_standard.xml" not in ff
        and "tip4pew_HFE_multivalent.xml" not in ff
    ]
    forcefield_settings["forcefields"] = preserved[:1] + _water_forcefield_files(solvent_model) + preserved[1:]
    protocol_config["forcefield_settings"] = forcefield_settings
    return protocol_config


def _resolve_path(path_str: str, base_dir: pathlib.Path) -> pathlib.Path:
    path = pathlib.Path(path_str)
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()


def _resolve_optional_paths(
    value: None | str | Iterable[str],
    base_dir: pathlib.Path,
) -> list[pathlib.Path]:
    if value is None:
        return []
    if isinstance(value, str):
        return [_resolve_path(value, base_dir)]
    return [_resolve_path(item, base_dir) for item in value]


def _load_molecules(paths: Iterable[pathlib.Path]) -> list[Any]:
    from gufe import SmallMoleculeComponent
    from rdkit import Chem

    molecules = []
    for input_path in paths:
        if input_path.is_dir():
            entries = sorted(input_path.iterdir())
        else:
            entries = [input_path]

        for entry in entries:
            suffix = entry.suffix.lower()
            if suffix == ".sdf":
                supplier = Chem.SDMolSupplier(str(entry), removeHs=False)
                molecules.extend(SmallMoleculeComponent(mol) for mol in supplier if mol is not None)
            elif suffix == ".mol2":
                mol = Chem.MolFromMol2File(str(entry), removeHs=False)
                if mol is not None:
                    molecules.append(SmallMoleculeComponent.from_rdkit(mol))

    if not molecules:
        joined = ", ".join(str(path) for path in paths)
        raise ValueError(f"No ligand molecules could be loaded from: {joined}")
    return molecules


def _load_protein(path: pathlib.Path) -> Any:
    from gufe import ProteinComponent

    suffix = path.suffix.lower()
    if suffix == ".pdb":
        return ProteinComponent.from_pdb_file(str(path))
    if suffix in {".cif", ".pdbx"}:
        return ProteinComponent.from_pdbx_file(str(path))
    raise ValueError(f"Unsupported protein file type for {path}")


def _build_mapper(mapper_config: Mapping[str, Any]) -> Any:
    from kartograf import KartografAtomMapper
    from lomap import LomapAtomMapper

    method = str(mapper_config.get("method", "kartograf")).lower()
    settings = dict(mapper_config.get("settings", {}))

    mapper_choices = {
        "kartograf": KartografAtomMapper,
        "kartografatommappper": KartografAtomMapper,
        "kartografatommapper": KartografAtomMapper,
        "lomap": LomapAtomMapper,
        "lomapatommapper": LomapAtomMapper,
    }
    try:
        mapper_cls = mapper_choices[method]
    except KeyError as exc:
        supported = ", ".join(sorted(mapper_choices))
        raise ValueError(f"Unsupported mapper '{method}'. Choose from: {supported}") from exc

    return mapper_cls(**settings)


def _build_network_planner(network_config: Mapping[str, Any], n_cores: int) -> Any:
    from openfe.setup.ligand_network_planning import (
        generate_lomap_network,
        generate_maximal_network,
        generate_minimal_redundant_network,
        generate_minimal_spanning_network,
        generate_radial_network,
    )

    method = str(network_config.get("method", "minimal_spanning")).lower()
    settings = dict(network_config.get("settings", {}))
    settings.setdefault("n_processes", n_cores)

    planner_choices = {
        "mst": generate_minimal_spanning_network,
        "minimal_spanning": generate_minimal_spanning_network,
        "generate_minimal_spanning_network": generate_minimal_spanning_network,
        "minimal_redundant": generate_minimal_redundant_network,
        "generate_minimal_redundant_network": generate_minimal_redundant_network,
        "maximal": generate_maximal_network,
        "generate_maximal_network": generate_maximal_network,
        "radial": generate_radial_network,
        "generate_radial_network": generate_radial_network,
        "lomap": generate_lomap_network,
        "generate_lomap_network": generate_lomap_network,
    }
    try:
        planner = planner_choices[method]
    except KeyError as exc:
        supported = ", ".join(sorted(planner_choices))
        raise ValueError(f"Unsupported network planner '{method}'. Choose from: {supported}") from exc

    return functools.partial(planner, **settings)


def _read_network_edges(path: pathlib.Path) -> list[tuple[str, str]]:
    suffix = path.suffix.lower()

    if suffix in {".yaml", ".yml"}:
        data = _load_yaml(path)
        edges = data.get("edges", data.get("pairs", data))
    elif suffix == ".json":
        with path.open() as handle:
            data = json.load(handle)
        edges = data.get("edges", data.get("pairs", data))
    else:
        raw_lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        edges = []
        for line in raw_lines:
            if line.startswith("#"):
                continue
            fields = line.replace(",", " ").split()
            if len(fields) == 5 and fields[1] == "#" and fields[3] == "->":
                edges.append((fields[2], fields[4]))
                continue
            if len(fields) == 3 and fields[1] in {"->", ">>"}:
                edges.append((fields[0], fields[2]))
                continue
            if len(fields) == 2:
                edges.append((fields[0], fields[1]))
                continue
            raise ValueError(f"Unsupported edge format in {path}: '{line}'")

    if not isinstance(edges, list):
        raise TypeError(f"Expected a list of edge pairs in {path}")

    parsed_edges = []
    for edge in edges:
        if not isinstance(edge, (list, tuple)) or len(edge) != 2:
            raise ValueError(f"Each edge in {path} must be a pair, got {edge!r}")
        parsed_edges.append((str(edge[0]), str(edge[1])))

    if not parsed_edges:
        raise ValueError(f"No edges found in {path}")
    return parsed_edges


def _build_explicit_network(
    ligands: list[Any],
    mapper: Any,
    network_config: Mapping[str, Any],
    base_dir: pathlib.Path,
) -> Any | None:
    from openfe.setup.ligand_network_planning import generate_network_from_names

    explicit_map = network_config.get("fep_map")
    if explicit_map is None:
        explicit_map = network_config.get("map_file")
    if explicit_map is None:
        explicit_map = network_config.get("edges_file")
    if explicit_map is None:
        return None

    edge_path = _resolve_path(str(explicit_map), base_dir)
    edge_pairs = _read_network_edges(edge_path)
    return generate_network_from_names(
        ligands=ligands,
        mapper=mapper,
        names=edge_pairs,
    )


def _build_solvent(solvent_config: Mapping[str, Any]) -> Any:
    from gufe import SolventComponent
    from openff.units import unit

    solvent_kwargs = dict(solvent_config.get("component", {}))
    ion_concentration = solvent_kwargs.get("ion_concentration")
    if isinstance(ion_concentration, str):
        solvent_kwargs["ion_concentration"] = unit.Quantity(ion_concentration)
    return SolventComponent(**solvent_kwargs)


def _build_protocol(protocol_config: Mapping[str, Any]) -> Any:
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

    defaults = RelativeHybridTopologyProtocol.default_settings()
    normalized = _normalize_protocol_config(protocol_config)
    normalized = _sync_forcefield_settings_with_solvent(normalized)
    merged = _deep_merge(defaults.model_dump(), normalized)
    settings = type(defaults).model_validate(merged)
    return RelativeHybridTopologyProtocol(settings=settings)


def _assign_partial_charges(
    molecules: list[Any],
    protocol: Any,
    n_cores: int,
    overwrite_charges: bool,
) -> list[Any]:
    from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges

    charge_settings = protocol.settings.partial_charge_settings
    return bulk_assign_partial_charges(
        molecules=molecules,
        overwrite=overwrite_charges,
        method=charge_settings.partial_charge_method,
        toolkit_backend=charge_settings.off_toolkit_backend,
        generate_n_conformers=charge_settings.number_of_conformers,
        nagl_model=charge_settings.nagl_model,
        processors=n_cores,
    )


class ConfigurableRBFEPlanner:
    def __init__(
        self,
        name: str,
        mapper: Any,
        network_planner: Any,
        explicit_network: Any | None,
        protocol: Any,
        adaptive_protocol: bool,
    ) -> None:
        from openfe.setup.atom_mapping.lomap_scorers import default_lomap_score

        self.name = name
        self.mapper = mapper
        self.network_planner = network_planner
        self.explicit_network = explicit_network
        self.protocol = protocol
        self.mapping_scorer = default_lomap_score
        self.adaptive_protocol = adaptive_protocol

    def __call__(
        self,
        ligands: list[Any],
        solvent: Any,
        protein: Any,
        cofactors: list[Any],
    ) -> tuple[Any, Any]:
        from gufe import AlchemicalNetwork, Transformation
        from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
        from openfe.setup.chemicalsystem_generator.easy_chemicalsystem_generator import (
            EasyChemicalSystemGenerator,
        )

        ligand_network = self.explicit_network
        if ligand_network is None:
            ligand_network = self.network_planner(
                ligands=ligands,
                mappers=[self.mapper],
                scorer=self.mapping_scorer,
            )
        chemical_system_generator = EasyChemicalSystemGenerator(
            solvent=solvent,
            protein=protein,
            cofactors=cofactors,
        )

        transformations = []
        nodes = []
        for ligand_mapping in ligand_network.edges:
            for stateA, stateB in zip(
                chemical_system_generator(ligand_mapping.componentA),
                chemical_system_generator(ligand_mapping.componentB),
            ):
                protocol_settings = self.protocol.settings.unfrozen_copy()
                if self.adaptive_protocol and isinstance(self.protocol, RelativeHybridTopologyProtocol):
                    protocol_settings = self.protocol._adaptive_settings(
                        stateA=stateA,
                        stateB=stateB,
                        mapping=ligand_mapping,
                        initial_settings=protocol_settings,
                    )

                transformation = Transformation(
                    stateA=stateA,
                    stateB=stateB,
                    mapping=ligand_mapping,
                    name=f"{self.name}_{stateA.name}_{stateB.name}",
                    protocol=self.protocol.__class__(settings=protocol_settings),
                )
                transformations.append(transformation)
                nodes.extend([stateA, stateB])

        alchemical_network = AlchemicalNetwork(
            nodes=nodes,
            edges=transformations,
            name=self.name,
        )
        return alchemical_network, ligand_network


def _write_outputs(
    alchemical_network: Any,
    ligand_network: Any,
    output_dir: pathlib.Path,
    config_path: pathlib.Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    network_json = output_dir / f"{output_dir.name}.json"
    alchemical_network.to_json(network_json)

    graphml_path = output_dir / "ligand_network.graphml"
    graphml_path.write_text(ligand_network.to_graphml())

    transformations_dir = output_dir / "transformations"
    transformations_dir.mkdir(parents=True, exist_ok=True)
    for transformation in alchemical_network.edges:
        transformation_name = transformation.name or transformation.key
        transformation.to_json(transformations_dir / f"{transformation_name}.json")

    shutil.copy2(config_path, output_dir / "config.used.yaml")


def _build_from_config(config_path: pathlib.Path) -> pathlib.Path:
    config = _load_yaml(config_path)
    base_dir = config_path.parent

    cache_dir = _resolve_path(
        config.get("cache_dir", ".openfe-cache"),
        base_dir,
    )
    _ensure_writable_cache_dirs(cache_dir)

    inputs_cfg = config.get("inputs", {})
    planning_cfg = config.get("planning", {})
    solvent_cfg = config.get("solvent", {})
    protocol_cfg = config.get("protocol", {})

    ligand_paths = _resolve_optional_paths(inputs_cfg.get("ligands"), base_dir)
    if not ligand_paths:
        raise ValueError("Config must define inputs.ligands")

    protein_path = inputs_cfg.get("protein")
    if protein_path is None:
        raise ValueError("Config must define inputs.protein")

    protein = _load_protein(_resolve_path(protein_path, base_dir))
    ligands = _load_molecules(ligand_paths)
    cofactors = _load_molecules(_resolve_optional_paths(inputs_cfg.get("cofactors"), base_dir)) if inputs_cfg.get("cofactors") else []

    n_cores = int(planning_cfg.get("n_cores", 1))
    overwrite_charges = bool(planning_cfg.get("overwrite_charges", False))
    adaptive_protocol = bool(planning_cfg.get("adaptive_protocol", True))
    campaign_name = str(planning_cfg.get("name", "rbfe_campaign"))

    mapper = _build_mapper(planning_cfg.get("mapper", {}))
    network_planner = _build_network_planner(planning_cfg.get("network", {}), n_cores=n_cores)
    solvent = _build_solvent(solvent_cfg)
    protocol = _build_protocol(protocol_cfg)

    ligands = _assign_partial_charges(
        molecules=ligands,
        protocol=protocol,
        n_cores=n_cores,
        overwrite_charges=overwrite_charges,
    )
    if cofactors:
        cofactors = _assign_partial_charges(
            molecules=cofactors,
            protocol=protocol,
            n_cores=n_cores,
            overwrite_charges=overwrite_charges,
        )

    explicit_network = _build_explicit_network(
        ligands=ligands,
        mapper=mapper,
        network_config=planning_cfg.get("network", {}),
        base_dir=base_dir,
    )

    planner = ConfigurableRBFEPlanner(
        name=campaign_name,
        mapper=mapper,
        network_planner=network_planner,
        explicit_network=explicit_network,
        protocol=protocol,
        adaptive_protocol=adaptive_protocol,
    )
    alchemical_network, ligand_network = planner(
        ligands=ligands,
        solvent=solvent,
        protein=protein,
        cofactors=cofactors,
    )

    output_dir = _resolve_path(config.get("output_dir", "rbfe_output"), base_dir)
    _write_outputs(
        alchemical_network=alchemical_network,
        ligand_network=ligand_network,
        output_dir=output_dir,
        config_path=config_path,
    )
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate OpenFE RBFE input JSON files from a YAML config.",
    )
    parser.add_argument(
        "config",
        type=pathlib.Path,
        help="Path to the RBFE YAML config.",
    )
    parser.add_argument(
        "--print-output-dir",
        action="store_true",
        help="Print only the generated output directory path.",
    )
    args = parser.parse_args()

    output_dir = _build_from_config(args.config.resolve())
    if args.print_output_dir:
        print(output_dir)
        return

    result = {
        "output_dir": str(output_dir),
        "network_json": str(output_dir / f"{output_dir.name}.json"),
        "transformations_dir": str(output_dir / "transformations"),
        "config_snapshot": str(output_dir / "config.used.yaml"),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
