# Underwater coverage reproducibility

Complete source and frozen inputs for **Probabilistic Coverage Optimization
for Underwater Sensor Networks under Three-Dimensional Localization
Uncertainty**.

Repository: https://github.com/houfence-design/underwater-coverage-reproducibility.

The project distinguishes a declared synthetic benchmark, real Copernicus
Marine GLORYS12V1 currents, and model-derived current-conditioned replay. It
does not report sensor sea trials, trained CFD--LSTM coupling, real-time control
or measured energy savings. Their research interfaces and evidence requirements
are retained in the manuscript.

`FIGURE_DATA_PROVENANCE.md` lists the exact source and interpretation of every
figure and table. No environmental record is manufactured. Mathematical-check
samples and Monte Carlo outcomes are expressly labeled computed quantities.

## Reproduce offline

Python 3.12 is recommended. From this directory:

```text
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python code/reproduce.py --mode verify
```

The entry point extracts `frozen_data.zip` automatically when the data directory
is absent. Every raw Monte Carlo outcome, input subset, seed, numerical result
and source hash is included. The archive also contains generated vector/raster
figures and readable manuscript, supplement and formula-companion PDFs.
No Copernicus login is needed to reproduce frozen
results. Verification checks the physical distributions, spherical integration,
union order, independent raw outcomes, chronological split, nominal constraints,
training-only selection, numerical refinement and Monte Carlo agreement.

Regenerate figures and LaTeX tables:

```text
python code/reproduce.py --mode figures
```

Rerun every numerical experiment, including all optimizer starts and 21 million
full target--array Monte Carlo trials:

```text
python code/reproduce.py --mode full
```

Runtime depends on hardware; original timings are recorded separately for the
synthetic and ocean-replay environments. Timing is not an energy or real-time
control metric. The ocean-replay script also supports `--resume` for an
interrupted run only when all input hashes and parameters match.

Compile the manuscript, Supplementary Information and formula companion with a
TeX Live installation exposing `pdflatex` and `bibtex`:

```text
python code/reproduce.py --mode figures --compile
```

Confidential editor/reviewer reports and their verbatim response documents are
kept in the author's local delivery package, outside the public repository.
Historical source supplied with the original is retained in `legacy_support`.
It is explicitly excluded from the corrected computational evidence; its
original notebook images and source files remain in the local revision archive.

## Redownload the real input

```text
python code/download_public_cmems.py
```

This uses Mercator Ocean International's explicitly public OceanBench archive.
It downloads only the two required original depth layers, extracts three cells,
then interpolates to 150 m. The downloaded derivative is daily at **1 degree**,
although its underlying GLORYS12V1 product is native 1/12 degree. Use
`DATA_LICENSE_AND_PROVENANCE.md` and the extraction manifest for exact sources,
original metadata and attribution. The upstream release version is not stated
by the public archive; object hashes provide exact identity.

`download_cmems.py` is an optional authenticated native-resolution alternative.
It writes a separate `cmems_native_optional` directory and is not the input to
the published historical replay. Install `copernicusmarine==2.5.0` separately
if using that optional path. Never commit account credentials.

## Formula-to-code inventory

`FORMULAS_AND_PSEUDOCODE.md`, `formula_registry.json`, and
`paper/formula_code_companion.pdf` link **every displayed equation and
computational inline expression** to pseudocode, units, executable functions
and source positions. The energy objective is a prospective interface that
rejects missing calibrated energy values.

## Experimental design

- Target: 100 × 20 × 25 m at 137.5--162.5 m absolute depth; sensing radius 20 m.
- Physical box: [-500,500]² × [120,270] m; nine-node structured grid, minimum nominal spacing 40 m.
- Parent residual scales: horizontal 120 m, vertical 40 m; actual coordinates physically truncated.
- Synthetic benchmark: four methods × 20 starts, maximum 600 calls each; 30 × 50,000 MC trials for five layouts.
- Real-current replay: three sites, 2018 planning and 2019 held out; four methods × 20 starts per site, maximum 300 calls each.
- Current displacement is an assumed 600 s frozen-flow proxy shared across sensors.
- Union probabilities are computed conditional on common forcing before spatial and temporal averaging.
- Replay MC: 30 × 50,000 trials for three layouts at each site; temporal intervals use block resampling separately.

## Persistent release

GitHub provides public source access. The code-release DOI required by the
journal remains pending deposition in a recognized DOI-assigning repository.
`ZENODO_METADATA.json` and the prepared release archive support that action.
The Copernicus product DOI identifies the ocean input and is not a DOI for this
project. No unissued code DOI is represented as assigned.

Original Python source is MIT licensed; external data and template components
retain their provider terms.
