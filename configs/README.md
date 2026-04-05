# RBFE Config Workflow

This repository now includes a config-driven RBFE input generator:

```bash
conda run -n openfe_v10_env python scripts/build_rbfe_network.py configs/rbfe_config.example.yaml
```

The script writes:

- one alchemical network JSON
- one `ligand_network.graphml`
- one JSON per transformation under `transformations/`
- one `config.used.yaml` snapshot for provenance

## Generate `*.json` Inputs

The OpenFE execution inputs are generated from a YAML config by
`scripts/build_rbfe_network.py`.

### 1. Prepare a config YAML

Start from the example:

```bash
cp configs/rbfe_config.example.yaml configs/my_rbfe.yaml
```

Edit at least these fields:

- `inputs.ligands`
- `inputs.protein`
- `planning.network`
- `protocol`

If you already have an external FEP edge map, set:

```yaml
planning:
  network:
    fep_map: path/to/fep_map_edges.txt
```

### 2. Generate the network JSON and transformation JSON files

Local conda environment:

```bash
conda run -n openfe_v10_env python scripts/build_rbfe_network.py configs/my_rbfe.yaml
```

Apptainer:

```bash
apptainer exec --bind "$PWD":"$PWD" --bind /tmp:/tmp openfe.sif \
  /opt/conda/envs/openfe/bin/python \
  "$PWD/scripts/build_rbfe_network.py" \
  "$PWD/configs/my_rbfe.yaml"
```

### 3. Generated files

If `output_dir: ../outputs/tyk2_rbfe` is set in the YAML, the generator writes:

- `outputs/tyk2_rbfe/tyk2_rbfe.json`
  The full alchemical network JSON.
- `outputs/tyk2_rbfe/ligand_network.graphml`
  The ligand network graph.
- `outputs/tyk2_rbfe/transformations/*.json`
  One JSON file per solvent/complex transformation. These are the files passed
  to `openfe quickrun`.
- `outputs/tyk2_rbfe/config.used.yaml`
  A snapshot of the exact config used.

### 4. Run one generated transformation

```bash
openfe quickrun \
  outputs/tyk2_rbfe/transformations/<transformation>.json \
  -d run-<transformation>
```

## Apptainer / HPC

For HPC runs, build the container locally and copy the resulting `.sif` file
to the cluster. This avoids solving the environment on the HPC side.

### 1. Build the container locally

The repository includes an Apptainer definition file at
`apptainer/openfe.def`. It copies this repository into the image, creates the
OpenFE conda environment, and installs the local checkout.

Build command:

```bash
apptainer build openfe.sif apptainer/openfe.def
```

Notes:

- This definition uses `Bootstrap: docker`, so the build machine needs network
  access to pull the base image.
- `apptainer build` often requires root, `sudo`, or a site-supported fakeroot
  setup depending on the local machine.
- The `%files` section currently points at this checkout:
  `/home/hiroyuki/code/openfe_v10`.
  If you build from a different path, update `apptainer/openfe.def` first.

After the build, transfer `openfe.sif` and your input files to the cluster.

### 2. Generate RBFE input JSONs inside the container

```bash
apptainer exec --bind "$PWD":"$PWD" --bind /tmp:/tmp openfe.sif \
  /opt/conda/envs/openfe/bin/python \
  "$PWD/scripts/build_rbfe_network.py" \
  "$PWD/configs/rbfe_config.example.yaml"
```

This creates the network JSON and the per-transformation JSON files from the
external YAML config.

### 3. Run one OpenFE transformation inside the container

```bash
apptainer exec --nv --bind "$PWD":"$PWD" --bind /tmp:/tmp openfe.sif \
  /opt/conda/envs/openfe/bin/openfe quickrun \
  "$PWD/outputs/tyk2_rbfe/transformations/<transformation>.json" \
  -d "$PWD/run-<transformation>"
```

Use `--nv` on GPU nodes. On CPU-only systems, omit it.

### 4. Submit as a batch job

A Slurm example is included at `apptainer/run_rbfe_job.sh`.

Typical usage:

```bash
mkdir -p logs
sbatch \
  --export=ALL,SIF_PATH=/path/to/openfe.sif,PROJECT_DIR=$PWD,CONFIG_PATH=$PWD/configs/rbfe_config.example.yaml,TRANSFORM_JSON=$PWD/outputs/tyk2_rbfe/transformations/<transformation>.json \
  apptainer/run_rbfe_job.sh
```

That script does two steps:

- builds the alchemical network JSONs from the YAML config
- runs `openfe quickrun` for one transformation JSON

The config is designed around OpenFE's RBFE tutorials and protocol docs:

- RBFE Python tutorial: <https://docs.openfree.energy/en/stable/tutorials/rbfe_python_tutorial.html>
- CLI tutorial: <https://docs.openfree.energy/en/stable/tutorials/rbfe_cli_tutorial.html>
- Relative Hybrid Topology protocol guide: <https://docs.openfree.energy/en/stable/guide/protocols/relativehybridtopology.html>

Sections:

- `inputs`
  Paths to ligands, protein, and optional cofactors.
- `planning`
  Campaign name, mapper choice, network planner choice, CPU count, and whether
  to keep OpenFE's adaptive protocol expansion for charge-changing transforms.
  If `planning.network.fep_map` is set, that explicit edge list is used in
  preference to the auto-generated network.
- `solvent`
  `SolventComponent` construction settings such as ions and concentration.
- `protocol`
  Nested overrides applied on top of
  `RelativeHybridTopologyProtocol.default_settings()`.
  `protocol.ligand_forcefield` is accepted as a short alias for
  `protocol.forcefield_settings.small_molecule_forcefield`.

Practical notes:

- Input and output paths are resolved relative to the config file location.
- `mapper.method` currently supports `kartograf` and `lomap`.
- `network.method` supports `minimal_spanning`, `minimal_redundant`,
  `maximal`, `radial`, and `lomap`.
- `planning.network.fep_map` accepts txt/json/yaml edge lists. This is the
  right place to inject an externally prepared FEP map.
- Quantity-like protocol values can be given as strings such as
  `1.0 nanosecond` or `0.15 molar`.
- Enum-like protocol strings are normalized before validation, so values like
  `TIP4PEW`, `CUDA`, and `PME` are accepted.
- If `protocol.forcefield_settings.forcefields` is not explicitly set, the
  script now synchronizes the Amber water/ion ffxml files to match
  `protocol.solvation_settings.solvent_model` for `tip3p`, `spce`, and
  `tip4pew`.
- Ligand force fields can be changed via `protocol.ligand_forcefield` or
  `protocol.forcefield_settings.small_molecule_forcefield`, for example
  `openff-2.2.1`, `gaff-2.11`, or `espaloma-0.3.2`.
- If `protocol.solvation_settings.solvent_model` is `TIP4PEW`, the script also
  synchronizes the Amber water ffxml files to `tip4pew_standard.xml` and
  `tip4pew_HFE_multivalent.xml` unless you explicitly override
  `protocol.forcefield_settings.forcefields`.
- If you want to compare FEP conditions, copy the example YAML and vary only
  the sections of interest, for example:
  `planning.mapper`, `protocol.solvation_settings.solvent_model`,
  `protocol.simulation_settings`, or `protocol.lambda_settings`.
