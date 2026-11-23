# Tutorial: ice fabric from quad-pol radar on the CReSIS servers

This page takes you from nothing to a fabric product on a CReSIS machine
using your own account, and then to the figures. The repo runs one method,
the node method, and this page uses the Ridge A survey
(2024_Antarctica_Ground2, January 2025) as the worked example, with the
numbers you should get. Nothing in the repo names a user, so every step
works the same for any account.

## What the method measures

Ice with a horizontal crystal fabric is birefringent. A radio wave
polarised along the fabric's slow axis travels slower than one polarised
across it, so the two echoes of one reflector arrive apart by a delay that
grows with depth. The node method reads that delay from the four channels
the accumulation radar records (HH, VV, HV, VH) and returns, for every
along-track block of 2 km:

- the **slow axis**: one azimuth per block, geographic, clockwise from
  north, mod 180;
- **dlam(z)**: the difference between the two horizontal fabric
  eigenvalues along that axis, as a depth profile with its uncertainty and
  its resolution.

It has three steps.

| step | what happens | run by |
|---|---|---|
| 1. Coregister | One range delay per channel, fitted in the firn, puts the four channels on one range axis. | `coreg_one.sh` (`run_coreg_frame`) |
| 2. Fit each block | Channel gains from the block itself (reciprocity and the firn). The slow axis from the delay between the two synthesised polarisations, read against synthesis azimuth. Then dlam(z) from two kinds of data in one linear solve: the half-wave **nodes** of the cross-polarised power, which count the whole fringes, and the carrier **phase** between neighbouring depth windows, which gives the slope between nodes. | `run_node_frame` (`ptt.nodeFabric`) |
| 3. Look | A table of blocks, sections, the depth model and the lateral field. | `scripts/figures/node_*.py` |

Nothing is unwrapped anywhere: the nodes give the fringe count, and the
phase is only ever differenced between neighbouring windows.

What it does not give, so you do not look for it:

- **An axis that turns with depth.** One axis per block. That holds at
  Ridge A; it would be wrong in a shear margin.
- **Reflection anisotropy.** Power is used only to find the node depths.
- **Structure finer than about 130 m in depth.** Finer stripes are a phase
  the reflectors add, not fabric. The product states its own resolution.
- **The top 100 m.** The firn (30-90 m) is the reference.

The estimator's header (`help ptt.nodeFabric`) has the full account;
[`method.md`](method.md) maps the functions behind it.

## 1. Clone and test

The repo is public, so an HTTPS clone works on the CReSIS machines without
GitHub credentials:

```sh
ssh <you>@mem1.cresis.ku.edu
cd /kucresis/scratch/$USER/scripts        # or wherever you keep code
git clone https://github.com/hoffmaao/fabric_anisotropy.git
```

That is the whole install. Do NOT copy scripts out of the clone: they find
their own checkout, and a stale copy elsewhere on your path is the one
thing that could run old code. To update later, `git -C <clone> pull`.

Check the clone before touching real data. The test suite builds synthetic
radar data with a known fabric, so it needs no data. It takes about 10 min
on mem1 (2 min on a recent laptop) and should end with
`3 passed, 0 failed`. `MATLAB_THREADS=8` keeps it to 8 cores; uncapped,
MATLAB takes most of a shared node:

```sh
MATLAB_THREADS=8 nice bash fabric_anisotropy/opr_fabric/test/run_all_tests.sh quick
```

## 2. Make a work root

Products and caches go in a **work root**, a scratch directory of your own
that holds `stages/`. Make one and point `FABRIC_ROOT` at it (put the
`export` in your `~/.bashrc` so every shell has it):

```sh
export FABRIC_ROOT=/kucresis/scratch/$USER/fabric
mkdir -p $FABRIC_ROOT/stages
```

The two lines below are used in every command on this page:

```sh
C=/kucresis/scratch/$USER/scripts/fabric_anisotropy      # your clone
ROOT=/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2 # the season's data
```

## 3. One frame

The example is Ridge A frame `20250111_01_003`: 10 km driven north-south,
five blocks of 2 km.

**Coregister.** About 3 minutes:

```sh
bash $C/opr_fabric/server/coreg_one.sh $ROOT 20250111_01 3
```

It prints `start 20250111_01_003` and then `done  20250111_01_003 rc=0 ok`.
The channels are now in
`$FABRIC_ROOT/stages/quadpol_cf/coreg_cache/creg_20250111_01_003.mat`
(1.6 GB for this frame), and near the end of the log
`$FABRIC_ROOT/coreg_logs/creg_20250111_01_003.log` are the delays it found,
in samples:

    source check: max rel diff 0 (IDENTICAL)
    row offsets: VV +0.204  HV +0.076  VH +0.127

Running the same line again prints `skip  20250111_01_003 (cache exists)`.
It is one delay per channel on purpose. A coregistration that varies with
depth also follows the delay between the slow and fast echoes and removes
it, and that delay is what the next step measures.

**Fit.** About half an hour for this frame (6 to 9 minutes a block), as a
detached job:

```sh
cd $FABRIC_ROOT
nohup nice /opt/sw/matlab/2024b/bin/matlab -batch \
  "addpath('$C/opr_fabric/server'); day_seg='20250111_01'; frm=3; run_node_frame" \
  > node_20250111_01_003.log 2>&1 &
```

The product is `$FABRIC_ROOT/stages/node/node_20250111_01_003_B2000.mat`.

## 4. Check what you got

The log has two lines per block. For block 1 of this frame:

    block 1: heading 3.9, gains |gr| 0.95 +21 deg |gt| 0.61 -92 deg, slow axis 110.0 geo
      (plateau ratio 2.61, contrast 45/0 1.37/0.50 dB), 6 nodes (phase) at 734.7 MHz,
      chi2/dof 1.10, abstain below Inf m, dlam 900-1400 0.069
    window phase 0.0095 ns, 12.4 numbers resolved, one value a mean over 124 m,
      residual lag-one -0.26

and for the whole frame you should get these numbers (the last digit can
move with the MATLAB release):

| block | heading | slow axis | plateau ratio | nodes | chi2/dof | dlam 900-1400 m |
|---|---|---|---|---|---|---|
| 1 | 3.9 | 110.0 | 2.61 | 6 | 1.10 | 0.069 |
| 2 | 179.4 | 112.6 | 2.60 | 7 | 1.10 | 0.070 |
| 3 | 178.4 | 110.6 | 2.63 | 7 | 1.06 | 0.071 |
| 4 | 176.5 | 112.2 | 3.06 | 7 | 1.18 | 0.073 |
| 5 | skipped: `heading spread 22.9 deg, 889 of 2055 traces within 10 deg` | | | | | |

Block 5 holds the turn at the end of the line. A block has to have one
antenna frame, so it is left out and recorded as skipped.

What to look at, per block:

| in the log | what it is | at Ridge A |
|---|---|---|
| `slow axis ... geo` | the block's slow axis, geographic | 100 to 120, most within 102 to 116 |
| `plateau ratio` | how cleanly the axis is determined: the separation of the two halves of the azimuth scan over their scatter | 2 to 15: about 3 on north-south lines, 6 on east-west ones |
| `N nodes` | half-wave nodes found down to 1500 m | 5 to 9, most often 6 or 7 |
| `chi2/dof` | misfit of the linear solve against its stated errors | 0.96 to 1.30 |
| `abstain below` | depth below which it returns NaN because the node count and the phase disagree; `Inf` means it never did | `Inf` in all 133 blocks |
| `dlam 900-1400` | mean strength over 900-1400 m | 0.057 to 0.094, median 0.069 |
| `one value a mean over` | the vertical resolution | about 130 m |

In MATLAB, the product holds a struct array `R` with one element per block
and the block table `meta`:

```matlab
P = load('node_20250111_01_003_B2000.mat');
b = 2;
plot(P.R(b).dlam, P.R(b).zm); axis ij                 % dlam(z), z in metres
P.R(b).theta_geo_deg                                  % slow axis, deg from north
[P.meta.lat(b), P.meta.lon(b)]                        % where the block is
```

| field of `R(b)` | meaning |
|---|---|
| `theta_geo_deg` | slow axis, geographic, clockwise from north, mod 180 |
| `zm`, `dlam`, `sigma_dlam` | depth [m], strength and its one-sigma, on a 10 m grid; NaN below `abstain_below_m` |
| `resolution_m` | the depth range one value of `dlam` averages over |
| `intervals` | the answer with no depth model at all: mean dlam over 125 m spans (`z1`, `z2`, `dlam`, `sigma`), from the summed phase. Where it and `dlam` differ, believe this one. |
| `nodes` | the node ladder: depth `z`, count `n`, residual `resid_ns` |
| `theta_scan` | the azimuth scan behind the axis, with `plateau_ratio` |
| `gains` | the channel gains applied |
| `skipped`, `error` | a block left out (a turn inside it) or one that failed, and why |

For this frame the mean of `dlam` over 300-1300 m is 0.057, 0.058, 0.059
and 0.061 in blocks 1 to 4, and the nodes of block 1 are at 423, 675, 840,
1006, 1173 and 1329 m.

## 5. A whole survey

Two launchers, both idempotent: a frame that is already done is skipped,
so after any failure you launch the same line again.

```sh
cd $FABRIC_ROOT
nohup bash $C/opr_fabric/server/coreg_batch.sh ridge_a > coreg_ridge_a.log 2>&1 &
# when that has finished:
nohup bash $C/opr_fabric/server/node_survey.sh ridge_a > node_ridge_a.log 2>&1 &
```

The first builds the cache of every frame, three at a time. The second
fits every cached frame, four at a time with four threads each.
`node_survey.sh ridge_a 4 4 20250111_01` does one segment only; the
script's header lists its other arguments.

For the 36 Ridge A frames (133 blocks) expect:

- caches: 42 GB in all, 0.3 to 3.2 GB a frame. That is a lot on a
  nearly full filesystem. If someone has built them already, do not build
  them again: symlink their `creg_*.mat` files into your
  `$FABRIC_ROOT/stages/quadpol_cf/coreg_cache/` and go straight to
  `node_survey.sh`.
- the fit: about 8 minutes a block, 18 hours of fitting in all. It took
  3 h 07 min at eight frames at a time (`node_survey.sh ridge_a 8`);
  allow about twice that at the default four.
- products: 27 MB in all.

## 6. Figures

The figure scripts are Python (numpy, scipy, h5py, matplotlib). mem1's
system Python has no matplotlib, and the products are small, so copy them
to your own machine and work there:

```sh
# on your machine, from a clone of this repo
mkdir -p ~/data/opr/scar/node
rsync -av <you>@mem1.cresis.ku.edu:/kucresis/scratch/<you>/fabric/stages/node/ ~/data/opr/scar/node/
python scripts/figures/node_export.py
python scripts/figures/node_sections.py
python scripts/figures/node_representation.py
```

Each reads `~/data/opr/scar/node` (or `$SCAR_DATA/node`, or a directory
given as its first argument) and writes into the repo's `figs/`, which git
ignores. Nothing is drawn on the panels but their letters; each script's
header holds the captions.

| script | writes | shows |
|---|---|---|
| `node_export.py` | `node_blocks.csv`, `node_profiles.npz`, `node_table.mat`, `node_heading_test.png`, `node_ladders.png`, `node_crossings.png` | one row per block; the axis against the heading it was driven at; the same ice driven in two directions |
| `node_sections.py` | `node_sections_ratio.png`, `node_sections_dlam.png` | every frame as a section, with its nodes on the data they were picked from |
| `node_representation.py` | `node_depth_model.png`, `node_lateral_field.png`, `node_resolved.npz` | what a column resolves with depth; whether the differences between blocks are real |

The numbers they print for the 36 Ridge A frames:

| | |
|---|---|
| blocks | 134 read, 133 fitted, one skipped for a turn, none failed |
| slow axis, blocks driven north-south | 108.0 (circular sd 4.3, 95 blocks) |
| slow axis, blocks driven east-west | 110.1 (circular sd 5.0, 38 blocks) |
| median dlam, north-south blocks, by band 200-400 / 400-600 / 600-900 / 900-1200 / 1200-1450 m | 0.022 / 0.051 / 0.060 / 0.066 / 0.075 |
| the same, east-west blocks | 0.024 / 0.048 / 0.061 / 0.065 / 0.076 |
| the same ice driven in two directions (27 pairs within 1.5 km) | axis differs by +0.6 deg in the median (rms 2.6); node depth by -3 m (rms 25); dlam by under 0.002 in every band |
| the same ice on another pass (24 pairs within 0.6 km) | dlam agrees to 0.0010 over 125 m, the axis to 0.7 deg |
| one value of dlam | a mean over 130 m; 11.8 independent values in a column; stated error 0.0021 |
| the layers' phase per window | 0.011 ns, 3.0 deg of carrier |
| the change near 480 m | depth 480 m (460 to 500 in most blocks), step +0.037 (0.028 to 0.048) |

Two of these are the tests to repeat at any new site before believing a
result there:

- **The axis does not turn with the antennas.** Blocks driven north-south
  and east-west give the same geographic axis (`node_heading_test.png`).
- **The same ice gives the same answer on another pass.** Blocks within
  1.5 km of each other, driven in different directions, agree in node
  depth, strength and axis (`node_crossings.png`).

## How fine the profile can go

One value of `dlam` is a mean over about 130 m of depth. That comes from
the depth model's correlation length (500 m), not from the grid the
profile is written on, and it can be changed. On the server, as a variant
pass beside the standard one:

```sh
bash $C/opr_fabric/server/node_survey.sh ridge_a 4 4 20250111_01_003 2000 _L250 "depth_model=struct('len_m',250);"
```

or on your own machine with no rerun at all, because a product stores the
measurements the solve was given: `node_solve.resolve(row, dict(node_solve.SMOOTH, len_m=250.0))`
on a block from `node_export.read_product`.

What each length gives, over the 133 Ridge A blocks:

| correlation length | independent values in 1.4 km | one value averages over | stated error | difference on another pass | log evidence against the default |
|---|---|---|---|---|---|
| 1000 m | 8 | 196 m | 0.0015 | 0.0007 | -113 |
| 500 m (default) | 12 | 130 m | 0.0021 | 0.0010 | 0 |
| 250 m | 18 | 85 m | 0.0031 | 0.0012 | -410 |
| 125 m | 27 | 53 m | 0.0044 | 0.0016 | -1293 |
| 60 m | 44 | 36 m | 0.0057 | 0.0025 | -2582 |

"Another pass" is 24 pairs of blocks under 0.6 km apart on different
frames. Three things to take from the table:

- **A finer profile is not a sharper picture.** The detail it adds is
  noise that the smoothness constraint lets through and shapes. Its width
  follows the setting, not the ice (40, 30 and 20 m at half height for the
  three finer rows), and a synthetic column with a perfectly smooth fabric,
  the measured coherence and the measured layer phase produces as much of
  it as a real block does (0.0010 and 0.0023 against 0.0014 and 0.0027 at
  250 and 125 m, away from the step near 470 m). The survey's evidence
  prefers the default.
- **That noise repeats, so repeating is not proof.** Most of it is a phase
  the reflecting layers put between the two polarisations. It belongs to
  the ice, which is why it comes back on another pass far inside the
  stated error, but it does not add delay to what lies below, as a change
  of fabric does.
- **The sampling does not matter.** On a 20 m, 5 m or 2.5 m grid
  (`model_dz`) the profile is the same to 0.0001, and with correlation
  windows of 100, 25 or 12.5 m in place of 50 m (`win_m`, `step_m`) to
  0.001, a quarter of the stated error.

For detail under 130 m, use the means the product carries with no
constraint at all (`intervals`: the phase summed over a span, 125 m by
default, `interval_m=62.5` for half that, at an error of about 0.004), and
believe a feature there only if the level below it changes too.

## Other sites

The launchers know the other accumulation-radar surveys (`taylor_dome`,
`thwaites`, `eastwind`), but the node method has been run only at Ridge A.
At a new site run one frame first and read its log as in step 4, then
repeat the two tests above. Where the fabric axis turns with depth, as it
does at Thwaites, one axis per block is the wrong model, and the method
has not been tested there.

## Do not do this

Each of these is a mistake that has been made with this radar's data.

- **Do not take the channels from the polarimetric product.** Take them
  from the cache (`creg_*.mat`). In the product, `sec` (VV) is not
  registered to HH at all: the HH-VV coherence on `20250111_01_003` is
  0.15 at 1100 m where the registered pair gives 0.32. And `sec_reg` is
  registered with a field of offsets, which follows the birefringent delay
  and removes it.
- **Do not unwrap the phase yourself and differentiate it.** A 1D unwrap
  slips wherever the coherence dips. On the same frame that gave dlam
  0.014 where the answer is 0.06.
- **Do not gate on coherence at 0.35-0.4.** Below the firn the coherence
  is 0.3-0.4 everywhere, so that threshold is a coin flip. The method
  keeps every window above 0.2 and weights it by its own error.
- **Do not read depth structure finer than `resolution_m`.** Bands of
  tens of metres are the reflectors' phase. Use `intervals` if in doubt.
- **Do not call layer flattening coregistration.** Shifting both channels
  along track by the same amount cannot change the phase between them.
- **Do not set colour limits from percentiles.** Limits of +-0.5 hide a
  0.06 signal. The figures here use 0 to 0.10.
- **Do not expect a rerun to change a finished frame.** The launchers skip
  products that exist. For a variant give it a suffix (`tagx`, the sixth
  argument of `node_survey.sh`), and check an option's name in the
  script's header: a misspelt one is silently ignored.

## Etiquette on the shared nodes

- Run long jobs under `nohup` and `nice`. The launchers already do, and
  cap each MATLAB at 4 or 8 threads.
- `/kucresis/scratch` is shared and nearly full. Keep outputs under your
  own directory, share caches by symlink, and delete test runs.
- Never delete someone else's cache or product.

## When something goes wrong

| symptom | cause and fix |
|---|---|
| `no stages/ in ...` from a launcher | set `FABRIC_ROOT` and `mkdir -p $FABRIC_ROOT/stages` (step 2). |
| `skip  <frame> (no polarimetric product)` | the season has no `CSARP_polarimetric` product for that frame; the range window, the surface and the trace positions come from it. Ridge A has 47 frames of channels and 36 with that product. |
| `skip  <frame> (locked)` | another job holds the frame, or one was killed. If nothing is running, remove `.lock_<frame>` in the cache directory (coregistration) or `$FABRIC_ROOT/locks/node_<frame>` (fit). |
| `no coregistered channels for <frame>` | step 3's first command has not run for that frame; the message prints the line to run. |
| `block N: ... skipped` in the node log | the track turns inside the block; that block has no single antenna frame and is left out. |
| `block N: FAILED: ...` | the estimator's own message follows; the block is recorded with its `error` and the rest of the frame is unaffected. |
| `Undefined function 'ptt.<name>'` | the script was copied out of the clone. Run it from the clone. |
| figures: `No module named 'h5py'` | `pip install numpy scipy h5py matplotlib` in the Python you run them with. |

## What used to be here

Until November 2026 the repo also held the estimators and pipelines that
came before the node method: a quad-pol least-squares fit of the HH-VV
coherence phase (`run_season_*`, `run_quadpol_pipeline`), the OPR dual-pol
module (`run_fabric`, `fabric_task`), the Ershadi et al. (2022) chain, a
Fujita-model fit, delta-k, a port of Rathmann's CMP inversion, and the
figure scripts of each. They were removed so that there is one way to
process the data. The last commit that has them is `fe54276`
(`git checkout fe54276` in a second clone brings them back to look at).
