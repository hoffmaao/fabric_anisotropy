# fabric_anisotropy

Ice crystal orientation fabric from the polarimetric data of the CReSIS
ground-based accumulation radar.

The repo runs one method, the **node method**. From the four channels the
radar records (HH, VV, HV, VH) it returns, for every 2 km of track, the
fabric's slow axis and the horizontal eigenvalue difference dlam(z) along
it, with its uncertainty and its resolution. It counts the whole fringes
of the birefringent delay from the half-wave nodes of the cross-polarised
power and takes the slope between them from the carrier phase, so nothing
is unwrapped.

**Start with [`docs/cresis_tutorial.md`](docs/cresis_tutorial.md).** It
covers cloning on the CReSIS servers, the test suite, one Ridge A frame
with the numbers you should get, a whole survey, the figures, and the
mistakes to avoid. A fresh clone runs as-is under any account.

The relation between delay and fabric is that of Rathmann (2026),
"Inferring glacier ice-crystal orientation fabrics from oblique
polarimetric radar" (Proc. R. Soc. A, RSPA-2026-0001; reference Python
code: https://github.com/nicholasmr/ptt), at nadir.

## Layout

| path | what it holds |
|---|---|
| `+ptt/` | the method: `nodeFabric` (the estimator) and the six functions it is built from |
| `opr_fabric/server/` | what runs on the CReSIS machines: `coreg_one.sh` and `run_coreg_frame.m` (step 1), `run_node_frame.m` (step 2), `coreg_batch.sh` and `node_survey.sh` (a whole survey) |
| `opr_fabric/test/` | the MATLAB test suite |
| `scripts/figures/` | the figure scripts, Python |
| `docs/` | the tutorial and two reference pages (below) |
| `figs/` | figure output, gitignored |

- [`docs/method.md`](docs/method.md) - what each function does and why
  the method is built the way it is. Each function also has a full header
  of its own (`help ptt.<name>`).
- [`docs/scripts.md`](docs/scripts.md) - every script, test and figure,
  and how the work root on the CReSIS machines is laid out.

## Data and figures

Neither products nor figures are committed.

- **Products** copied from the server go under `~/data/opr/scar/node`;
  the `SCAR_DATA` environment variable moves `~/data/opr/scar`.
- **Figures** go to `figs/`, which `.gitignore` excludes. Each figure
  script takes an `<out_dir>` argument and defaults to `figs/`.

## Requirements

MATLAB with the Signal Processing and Statistics toolboxes (R2024b on the
CReSIS machines, where it is `/opt/sw/matlab/2024b/bin/matlab`; the tests
also pass on R2025b). The OPR toolbox is not needed. The figure scripts
need Python 3 with numpy, scipy, h5py and matplotlib.
