# Data directory

The real observational catalogue is not redistributed here. It is downloaded
by `code/data_source.py` from Chae's corrected-RUWE release (Zenodo record
10986733) and stored under the ignored `data/external/` directory.

Authoritative observational file:

- filename: `gaia_dr3_MSMS_d200pc_ruwe.csv`
- rows: 81,880 pairs
- MD5: `1b6c5063163a4e6c07043d13aeb70f55`

Files whose names begin with `legacy_newton_mock_v5_` are reductions of Chae's
official `Newton_dr3_MSMS_d200pc_5.csv` (Zenodo record 10652994, MD5
`997069b5635200853896d005a238dde9`). They are synthetic Newtonian data. They
are retained solely so that the v1/v2 provenance error can be reproduced; they
must never be used for an observational inference.

`legacy_phase_G_results.csv` is the old synthetic grid. The corrected grid is
`../results/corrected_phase_G_results.csv`.
