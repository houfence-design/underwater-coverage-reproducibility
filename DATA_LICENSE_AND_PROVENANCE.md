# Data attribution and provenance

The environmental input is the historical GLORYS12V1 derivative openly
distributed by Mercator Ocean International through its official OceanBench
public EDITO archive. It is Copernicus Marine ocean reanalysis, not direct
sensor-position observations.

- Original product: Global Ocean Physics Reanalysis, GLORYS12V1.
- Historical product alias: GLOBAL_REANALYSIS_PHY_001_030.
- Current product identifier: GLOBAL_MULTIYEAR_PHY_001_030.
- Product DOI: https://doi.org/10.48670/moi-00021.
- Public archive documentation: https://oceanbench.readthedocs.io/en/latest/one-degree-track.html.
- Exact access implementation: https://github.com/mercator-ocean/oceanbench/blob/main/oceanbench/datasets/reference.py.
- Downloaded archive: 1-degree horizontal derivative, daily 2018--2019 records.
- Native product resolution: 1/12 degree; this resolution is **not** assigned to the downloaded derivative.
- Retained original depths: 130.66600036621094 and 155.85069274902344 m.
- Analysis depth: 150 m, linear interpolation in depth only.
- Stations: 115.5 E / 20.5 N; 130.5 E / 25.5 N; 44.5 W / 40.5 N.
- Planning year: 2018; held-out year: 2019; 730 finite daily records per site.
- Upstream release version is not stated in the archive. Exact URL and SHA-256 records identify every downloaded object instead.

The provider's terms remain applicable to these data. The source-code MIT
license does not relicense external ocean data. See the provider's
[licence](https://marine.copernicus.eu/user-corner/service-commitments-and-licence)
and [citation guidance](https://help.marine.copernicus.eu/en/articles/4444611-citing-copernicus-marine-products-and-services).

Redistribution acknowledgement: E.U. Copernicus Marine Service Information,
https://doi.org/10.48670/moi-00021.

Source: European Union, Copernicus Marine Service Information or Data
(2018--2019), ©Mercator Ocean. Public historical derivative distributed by
Mercator Ocean International through OceanBench.

Publication acknowledgement: This study has been conducted using E.U.
Copernicus Marine Service Information: https://doi.org/10.48670/moi-00021.

`data/cmems/provenance.json` and `public_download_requests.json` freeze extraction
details and checksums. Three small NetCDF files retain both original depth
layers; `currents_150m.csv` retains all interpolated daily velocities. Reproduction
from the frozen subset works offline and requires no ocean-service account.

Velocities are externally produced real reanalysis data. The frozen-current
displacement rule, residual Gaussian position scales and replay coverage are
explicit model assumptions/derived calculations. Daily one-degree data cannot
resolve flow within the 100 m target or calibrate sensor positioning error.
