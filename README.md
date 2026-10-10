# Underwater coverage: code and reproducibility

This repository contains computational source, parameters, genuine public ocean inputs,
archived numerical results, and reproduction instructions.

## Run offline

Python 3.12 is recommended. From the repository root:

```text
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python code/reproduce.py --mode verify
```

The entry point safely extracts `frozen_data.zip` when `data/summary.json` is absent.
No ocean-service account or LaTeX installation is needed. Checks cover physical
distributions, spherical integration, union probabilities, input hashes, all raw
Monte Carlo outcomes, chronological train/test separation, and numerical refinement.

```text
# Rebuild plots, numerical tables and formula/source mappings
python code/reproduce.py --mode figures
# Recompute all experiments, then verify and regenerate outputs
python code/reproduce.py --mode full
```

Generated plots are written to `results/figures`; numerical table exports to
`results/tables`. These are generated outputs excluded from version control.
The full run comprises 320 optimizer runs and 21 million target--array Monte
Carlo trials. Runtime depends on hardware. Seeds and exact settings are retained.

## Real ocean inputs

The input is Mercator Ocean International's explicitly public OceanBench
CMEMS/GLORYS12V1 reanalysis derivative: daily, 1-degree data from three ocean
cells in 2018--2019, totaling 2,190 daily records. Two original depth layers
are retained and linearly interpolated to 150 m. This input is reanalysis;
displacement and coverage estimates are computations under declared assumptions.

```text
python code/download_public_cmems.py
```

`DATA_LICENSE_AND_PROVENANCE.md` records original URLs, provider attribution,
resolution, terms, depth layers and stations. Exact requests, metadata and
SHA-256 hashes are in `data/cmems`. Data-source attribution is required for
redistribution and reproducibility.

`download_cmems.py` is an optional native-resolution alternative. Install
`copernicusmarine==2.5.0` separately and use your own account if running it.
Its output goes to a separate directory and is not the frozen historical input.

## Source and frozen data

- `code/`: current experiments, verification, downloading, plotting and table export.
- `legacy_support/`: original supplied Python code, archived separately from the corrected calculations.
- `FORMULAS_AND_PSEUDOCODE.md` and `formula_registry.json`: 27 definitions, assumptions, pseudocode and executable function locations.
- `data/parameters.json`, `data/ocean/parameters.json`: complete experimental settings.
- `data/optimizer_runs.csv`, `data/ocean/optimizer_runs.csv`: repeated runs and recorded timings.
- `data/mc_raw`, `data/ocean/mc_raw`: original computed binary trial outcomes.
- `data/cmems`: real downloaded inputs and their extraction provenance.
- `PUBLIC_FILE_MANIFEST.json`: source/archive integrity checks.
- `.github/workflows/reproducibility.yml`: automatic offline verification.

The frozen archive contains computational inputs and numerical records only.
Mathematical benchmarks and generated Monte Carlo outcomes are identified as
computations; they are not presented as environmental observations.

New computational code is MIT licensed. External data and legacy source retain
their original terms, as stated in the applicable attribution files.
