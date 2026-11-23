# Scripts

Every script, test and figure around the method code in
[`method.md`](method.md). [`cresis_tutorial.md`](cresis_tutorial.md) walks
through them on Ridge A.

## The node method, end to end

On the server (`opr_fabric/server/`):

- `run_coreg_frame.m` - step 1 for one frame. Reads the range window and
  the surface from the season's polarimetric product, loads the four
  channels `CSARP_standardphase_{HH,VV,HV,VH}` over it, checks that HH is
  the image the product was built from, fits one range delay per channel
  over the firn (20-150 m) and writes the coregistered channels, down to
  1500 m, to `stages/quadpol_cf/coreg_cache/creg_<tag>.mat` (0.3-3.2 GB a
  frame at Ridge A, 1-3 min). A cache made with the same settings is left
  alone.
- `coreg_one.sh <site_root> <day_seg> <frm>` - that script under a lock,
  with a log in `coreg_logs/`. `coreg_batch.sh <site> [K]` runs it over
  every frame of a survey, K (3) at a time, and skips frames that have a
  cache or no polarimetric product.
- `run_node_frame.m` - step 2 for one frame. Reads the cache, cuts it into
  blocks of `blk_m` (2000 m) with `quadpol_frame_blocks.m`, leaves out
  traces heading more than 10 deg from their block's mean (a turn), and
  runs `ptt.nodeFabric` on each block. Writes
  `stages/node/node_<tag>_B<blk_m><tagx>.mat`: `R`, one struct per block
  (the estimator's output plus `theta_geo_deg`), and `meta`, the block
  table (along-track position, lat, lon, heading, trace count). Its header
  lists the optional variables; `z_max` can cut the profile shallower than
  the cache's 1500 m, and a deeper one is an error.
- `quadpol_frame_blocks.m` - the block cut: channels from the cache,
  positions from the season's polarimetric product, the per-trace heading
  as a great-circle azimuth smoothed on the doubled-angle phasor.
- `node_survey.sh <site> [max_jobs] [threads] [pattern] [blk_m] [tagx]
  [extra] [cache_dir]` - `run_node_frame` over every cached frame of a
  site, `max_jobs` (4) at a time. Skips frames whose product exists, locks
  the frame it is fitting, and reports how many products it ended with.
  `tagx` names a variant pass so it never replaces the standard products;
  `extra` is MATLAB run before the frame, e.g. `"psi_step_deg=2.5;"`.
- `fabric_paths.m`, `fabric_paths.sh` - where the code and the work root
  are (next section).

Sites known to the two launchers: `ridge_a`, `taylor_dome`, `thwaites`,
`eastwind`. The method has been run at Ridge A only.

## Running on the CReSIS machines

Processing happens in a scratch work root, laid out as

    <work>/code          this repo (it may also live anywhere else)
    <work>/stages        products: stages/quadpol_cf/coreg_cache/, stages/node/
    <work>/coreg_logs, <work>/logs, <work>/locks

Nothing names that root. `opr_fabric/server/fabric_paths.m` derives both it
and the repo root by walking up from its own location until it finds the
`+ptt` toolbox - looking for the toolbox rather than counting directories,
so it stays right if the tree is nested differently and fails loudly rather
than putting a wrong directory on the path.

The shell launchers split the same two roles the MATLAB side does, so one
name cannot mean two things across the two languages. CODE location comes
from the file: a launcher reads its worker from its own directory
(`"$FAB_DIR/<worker>.sh"`), because a worker is a sibling in the repo and
stays one wherever the tree is copied, and it passes that same directory to
MATLAB as an explicit `addpath` in every `-batch` string, so no run depends
on the operator's cwd or on a personally saved MATLAB path.

The PRODUCT root is derived by `opr_fabric/server/fabric_paths.sh` on the
SAME `+ptt` anchor, so the two sides cannot disagree about which tree a run
belongs to: the repo root is the directory holding `+ptt`, and the work root
is its parent. That is ONE candidate, never a list, and `stages/` is
confirmation on it. A missing `stages/` is fatal and the walk stops there
rather than continuing to an ancestor - binding to an ancestor's work root
is worse than not resolving at all, since a second checkout would then take
its locks in the LIVE tree while MATLAB, on the same `+ptt` anchor, wrote
products into the test tree.

`FABRIC_ROOT` overrides the product root on **both** sides and moves neither
the code nor the workers. It **picks the candidate; it does not exempt it**:
a work root is a work root by the same test however it was obtained, so the
override goes through the identical checks the derived value does - exists,
can be entered, resolved to an absolute path, and (in the shell) holds
`stages/`. An override is more likely to be wrong than a derived path, not
less: `FABRIC_ROOT=$HOME`, or a path off by one component, is a typo nothing
else would catch. It is rejected rather than created, because a typo that is
created on first write quietly collects the products. A relative value
resolves against the working directory on both sides.

Every failure is loud on both sides - the shell returns non-zero so the
launchers' `|| exit 1` fires, and MATLAB errors. Neither can hand back an
empty path: `set -u` does not catch a set-but-empty variable, and an empty
`$FAB` would resolve every path against `/`, miss every skip-if-cached
check, and recompute the survey.

The two helpers differ in exactly two ways. Both are claims about the PAIR
rather than about either side, so they are owned here and the two headers
point at this section instead of restating them.

**`stages/`.** The launchers require it to pre-exist; `fabric_paths.m` does
not. The reason is silent-versus-loud failure. The two MATLAB scripts create
the directories they write under `stages/`, so against a missing `stages/`
they simply make it. The launchers have no such backstop: their
skip-if-cached checks would match nothing and their `mkdir` locks would land
in a fresh tree, so a run pointed at the wrong root would quietly recompute
every frame. Pre-checking is what turns that silence into an error.

**Symlinks.** `bash`'s `pwd` is logical, so a work root reached through a
symlink comes back from `fabric_paths.sh` as the path written;
`fabric_paths.m` goes through MATLAB's `pwd`, which reports the OS `getcwd`
and returns the resolved target. Different strings for the same directory.
This is benign today - nothing compares them, both sides only join them onto
further path components, so the difference reaches log and error text only.
It would stop being benign if anything compared a shell-derived root against
a MATLAB-derived one, keyed a cache or lock name on the string, or recorded
it in a product for a later run to match.

So the launchers run in place out of the checkout - `bash
<clone>/opr_fabric/server/coreg_batch.sh ridge_a` - **there is no deploy
step**, and **nothing carries a username**, so the same tree runs from any
user's scratch and a second checkout can sit beside the live one for
testing.

DATA paths stay absolute - the `/cresis/...` season roots are real mount
points, not part of the work tree, and are named in the two launchers.

### Standing up a work root

The repo is public, so the CReSIS machines clone it over HTTPS with no
GitHub credentials:

    mkdir -p <work>/stages
    git clone https://github.com/hoffmaao/fabric_anisotropy.git <work>/code

A clone that lives somewhere else, as in
[`cresis_tutorial.md`](cresis_tutorial.md), serves just as well: export
`FABRIC_ROOT=<work>` and the scripts take the work root from there instead
of from the clone's parent. For a branch that has not been pushed, ship the
history as a bundle:

    git bundle create /tmp/fabric.bundle <branch>            # locally
    cat /tmp/fabric.bundle | ssh mem1 'cat > ~/fabric.bundle' # ship
    git clone --branch <branch> ~/fabric.bundle <work>/code   # on mem1

The `mkdir` is what makes `<work>` a work root: it is the `stages/` the
launchers confirm before they will run, so create it first.

## Tests

`opr_fabric/test/`. Every test is a script that errors if any of its checks
fails. `run_all_tests.m` runs them as one suite and names every failure at
the end; `run_all_tests.sh` is the same from a shell: it finds MATLAB and
exits non-zero if any test failed. The suite is the gate's test command in
`.no-mistakes.yaml`.

```sh
bash opr_fabric/test/run_all_tests.sh quick   # 2 min on a laptop, about 10 on mem1
bash opr_fabric/test/run_all_tests.sh test_node_fabric
```

On a shared CReSIS node, prefix `MATLAB_THREADS=8 nice`: uncapped, MATLAB
takes most of the node's cores.

- `test_node_fabric.m` - `ptt.nodeFabric` on four-channel synthetic columns
  (`synth_quadpol_column.m`: speckle, an off-centre band, uncalibrated
  gains, a layer-locked reflection phase, and depolarised backscatter that
  goes through the column as real echoes do). Ten verdicts: the delay
  scale, both members of the 90-deg alias pair with the slow one named,
  node depths and dlam(z) on a clean column, the same with uncalibrated
  channels at the coherence of Ridge A, a weak column with no nodes, a
  missing trace, nodes counted in the band's frequency and not the LO's,
  a reflection phase that must not come back as fabric, and independence
  of the correlation window and of the grid.
- `test_xcorr_traveltime.m` - the delay measurement: no fringe slip, with
  and without depolarisation, calibrated uncertainty, the slope of the
  delay giving the strength, and what a pair rotated off the fabric axes
  reads.
- `test_coregister_constant.m` - the coregistration: the three chain
  delays recovered to 0.02 samples while a deeper birefringent lag is left
  in place.

## Figures

`scripts/figures/`, Python with numpy, scipy, h5py and matplotlib. Each
reads a directory of node products (first argument; default
`$SCAR_DATA/node`, `$SCAR_DATA` defaulting to `~/data/opr/scar`) and writes
into `figs/` (second argument), which git ignores. Nothing is drawn on the
panels but their letters; each script's header holds the captions.

- `node_export.py` - reads every product; writes `node_blocks.csv` (one
  row per block), `node_profiles.npz` and `node_table.mat` (the profiles
  as one flat table, a row per block and 25 m of depth); prints the
  survey's axis and band strengths by heading family; draws
  `node_heading_test.png`, `node_ladders.png` and `node_crossings.png`.
- `node_sections.py` - every frame as a section with its nodes:
  `node_sections_ratio.png` on the cross-pol ratio they were picked from,
  `node_sections_dlam.png` on the stored dlam(z).
- `node_solve.py` - `ptt.nodeFabric`'s linear solve in Python, from the
  measurements a product stores, so another depth model (a finer one, or
  none) can be tried over a whole survey without the radar data.
  Reproduces the stored profile to 1e-10.
- `node_representation.py` - re-solves every block under the depth model
  and reports what a column resolves with depth, whether the 10 m
  structure accumulates, whether the differences between blocks repeat on
  independent passes, and the horizon near 480 m. Draws
  `node_depth_model.png` and `node_lateral_field.png`; writes
  `node_resolved.npz`.
- `morrow_style.py` - the figure style those scripts share.
