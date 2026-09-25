# Tutorial: running the fabric inversion on the CReSIS servers

This page takes you from nothing to a fabric product on a CReSIS machine
using your own account. It assumes you already process radar data with the
OPR toolbox there, meaning your `startup.m` sets up `gRadar`. Nothing in this
repo names a user, so every step below works the same for any account.

## Which methods this repo runs

There are two production chains. They need different inputs and measure
different things.

| | **Dual-pol traveltime** (`run_fabric`) | **Quad-pol coherence** (`run_season_*`) |
|---|---|---|
| Input | an OPR `CSARP_polarimetric` product (HH and a rotated VV, coregistered, SNAPHU-unwrapped phase) | the four raw channels `CSARP_standardphase_{HH,VV,HV,VH}` |
| Measures | the traveltime difference between the two polarizations, dtau(depth) | the complex HH-VV coherence at every synthesized antenna azimuth |
| Method | dtau from the unwrapped interferogram phase, with the coregistration offsets fixing the sign and the whole-fringe ambiguity (`ptt.blendTraveltime`), averaged over blocks of 1000 traces (`ptt.blockAverage`). It is then inverted for dlam in one smoothness-regularized joint solve through the Maxwell-Garnett firn column of Rathmann (2026) (`ptt.invertBlocks`). | a least-squares fit of the whole coherence field to the exact single-column birefringence model (`ptt.quadpolFabricLS`). It uses CRB coherence weights, fits a nuisance for the antenna-frame pedestal, and never unwraps; where the data cannot decide, it returns NaN. The fabric axis comes from a frame pass over ~2 km segments (`ptt.quadpolFrameTheta`), and each ~125 m block is fitted on that axis. |
| Gives | dlam = lam_sec - lam_ref, the horizontal fabric contrast between the product's two FIXED synthesized axes, per depth interval and block. It does NOT estimate orientation. | the fabric orientation theta0 (geographic, mod 180) AND the contrast dlam, per depth window and block |
| Output | `CSARP_fabric_joint/<day_seg>/Data_<frame>.mat` plus `_dlam.jpg` and `_dtau.jpg`, in the OPR season tree | `stages/quadpol/quadpol_section_<frame>.mat` under your work root |
| Cost | ~15 s per frame (debug); a few minutes under slurm, including the ~3 min compile on every launch | ~1 h per Ridge A frame from its coregistration cache, plus ~50 min to build the cache the first time |

The quad-pol pipeline also runs the published Ershadi et al. (2022) chain
(`ptt.ershadiFabric`) on the same data and saves it alongside, **for
comparison only**. On this radar a flat ~-3.6 dB cross-pol pedestal locks its
cross-pol minimum to the antenna, so do not report its orientation.

`opr_fabric/server/` also holds research runners that are not production
chains, including `run_ershadi_r` (reflection ratio), `run_nymand_step1`
(an independent power-based check of theta0) and `run_deltak_*`
(split-spectrum dtau). [`scripts.md`](scripts.md) lists what each one
concluded, and [`method.md`](method.md) covers every `ptt.*` function.

## 1. Clone

The repo is public, so an HTTPS clone works on the CReSIS machines without
GitHub credentials. Put it next to your other OPR scripts:

```sh
ssh <you>@mem1.cresis.ku.edu
cd /kucresis/scratch/$USER/scripts        # or wherever you keep opr/ and run_opr/
git clone https://github.com/hoffmaao/fabric_anisotropy.git
```

That is the whole install. Do NOT copy anything into `opr/matlab` or
`run_opr`. The scripts find their own checkout, and a stale copy elsewhere
on your path is the one thing that could run old code.

To check the clone before touching real data, run the fast half of the test
suite. It builds synthetic radar data with a known fabric and stubs out OPR,
so it needs neither data nor your OPR setup. It takes about 4 min on mem1 and
should end with `14 passed, 0 failed`. `MATLAB_THREADS=8` keeps it to 8
cores; uncapped, MATLAB takes most of a shared node:

```sh
MATLAB_THREADS=8 nice bash fabric_anisotropy/opr_fabric/test/run_all_tests.sh quick
```

## 2. Run the dual-pol example

Start MATLAB the way you normally do for OPR, so your `startup.m` runs
(`/opt/sw/matlab/2024b/bin/matlab`). Then run the script from any directory:

```matlab
run('/kucresis/scratch/<you>/scripts/fabric_anisotropy/opr_fabric/run_fabric.m')
```

The shipped example inverts Ridge A frame `20250108_02_009` from the
2024_Antarctica_Ground2 season. It puts the checkout on the path, checks
that OPR is set up, creates the task, runs it, and prints
`Chain 1 succeeded`. Output lands under your `gRadar.out_path`:

    <out_path>/accum/2024_Antarctica_Ground2/CSARP_fabric_joint/20250108_02/
        Data_20250108_02_009.mat        dlam, dlam_top_depth, dlam_bot_depth,
                                        dlam_quality, dlam_ref_degenerate, dtau_blk, ...
        Data_20250108_02_009_dlam.jpg   dlam by depth interval and block
        Data_20250108_02_009_dtau.jpg   block-averaged dtau against two-way time

For this frame, dlam rises from ~0.03 near the surface to ~0.06-0.07 below
1000 m. The first interval carries `dlam_ref_degenerate = true`, so quote
fabric from ~200 m down.

## 3. Point it at your own data

Edit the **User Setup** section of `run_fabric.m`. You normally change four
lines:

```matlab
params = read_param_xls(opr_filename_param('accum_param_<season>.xlsx'));
params = opr_set_params(params,'cmd.generic',1,'day_seg','<YYYYMMDD_SS>');
params = opr_set_params(params,'cmd.frms',[<frames>]);          % [] = all
params = opr_set_params(params,'fabric.in_path','<product>');   % CSARP_<product>
```

`fabric.in_path` is the polarimetric product to read, WITHOUT the `CSARP_`
prefix. Products are named after whoever processed them
(`polarimetric_unwrap_dlilien`, `polarimetric_ndh`, ...). If you name one
that does not exist, the run stops before submitting anything and lists the
season's `CSARP_polarimetric*` products for you to choose from. Frames
missing from the product you chose are skipped, each with a message.

Keep `fabric.out_path` distinct from other people's products in the shared
season tree. `fabric_joint` is the convention. An absolute path such as
`/kucresis/scratch/$USER/fabric_out/CSARP_fabric_joint` writes there instead.

The remaining knobs (coherence threshold, block size, number of depth
intervals, the firn column) are documented where they are set, in the
header of `opr_fabric/fabric.m`. The defaults are the ones every published
section used.

## 4. Scale up

`run_fabric.m` runs in `debug` mode by default, one frame after another in
your MATLAB session, which is fine for a season at ~15 s a frame. To use the
cluster instead, switch the `if 0` above `cluster.type = 'slurm'` to
`if 1`. OPR then compiles `fabric_task` into your `cluster_job` as it does
for any other task. It recompiles on every slurm launch (about 3 min),
because `cluster_job` is shared with every other OPR task. The compile runs
`mcc` through the shell, and a login shell on mem1 does not have MATLAB's
`bin/` on `PATH`, so start MATLAB with it there, or the run stops at
`mcc: command not found`:

```sh
export PATH=/opt/sw/matlab/2024b/bin:$PATH
matlab
```

To only create the tasks and run them later (OPR's
usual behavior), set `run_chain = false` and run
`cluster_run(cluster_load_chain(<id>))` afterwards.

## 5. Run the quad-pol pipeline

The quad-pol scripts write intermediate products and ~0.7 GB coregistration
caches, so they need a **work root**, a scratch directory holding `stages/`.
Make one and point `FABRIC_ROOT` at it:

```sh
export FABRIC_ROOT=/kucresis/scratch/$USER/fabric
mkdir -p $FABRIC_ROOT/stages/quadpol
```

(Without `FABRIC_ROOT` the work root is the directory above the clone,
which here is your `scripts/` directory. The shell launchers refuse to run
there unless it has a `stages/`, rather than guessing.)

One frame, as a detached job on the shared node:

```sh
C=/kucresis/scratch/$USER/scripts/fabric_anisotropy
cd $FABRIC_ROOT
nohup nice /opt/sw/matlab/2024b/bin/matlab -batch \
  "maxNumCompThreads(8); addpath('$C/opr_fabric/server'); day_seg='20250108_02'; frm=9; run_season_ridge_a" \
  > quadpol_20250108_02_009.log 2>&1 &
```

The `run_season_*.m` files are one per site: `ridge_a`, `taylor_dome`,
`thwaites_wais_mcmurdo`, `eastwind` and `eastgrip`. Each declares its
season's data root and quirks in its header. The result is
`$FABRIC_ROOT/stages/quadpol/quadpol_section_20250108_02_009.mat`, with one
struct `res`. The fabric axis theta0 is per segment and depth window
(`ls_theta_seg`, with `ls_dlam_seg` and their jackknife errors
`ls_se_*_seg`). The LS block section is `sec_dlam_ls` with `sec_se_dlam_ls`
and `sec_resid_ls` per window and block, each block fitted on its segment's
axis. The frame pedestal is `ls_pedestal`, and the Ershadi comparison is
`sec_theta` / `sec_dlam`. Where the season carries `CSARP_layer` bed picks
the record is blanked below each trace's bed before anything is fitted, and
`sec_bed` holds each block's median pick (NaN where there is none, as at
Ridge A, whose record ends in ice).

A whole survey runs in two stages. First build the coregistration caches,
the expensive pass. Then invert, which follows the cache builder frame by
frame:

```sh
nohup bash $C/opr_fabric/server/coreg_batch.sh ridge_a 3 > coreg_ridge_a.log 2>&1 &
nohup bash $C/opr_fabric/server/invert_batch.sh ridge_a 2 > invert_ridge_a.log 2>&1 &
```

Both are idempotent: frames already done are skipped, and failures are
retried a capped number of times. If a run dies, just launch it again.

## Etiquette on the shared nodes

- Run long jobs under `nohup` and `nice`, and cap MATLAB at
  `maxNumCompThreads(8)`. The launchers already do this.
- Coreg caches cost ~50 min each to rebuild. Don't delete them casually, and
  do share them: symlink another user's `creg_<frame>.mat` into your own
  `stages/quadpol/coreg_cache/` to rerun an inversion in minutes.
- `/kucresis/scratch` is a shared, nearly full filesystem. Put outputs under
  your own directory, and clean up test runs.

## When something goes wrong

| symptom | cause and fix |
|---|---|
| `OPR is not set up in this MATLAB session` | your `startup.m` did not run. Start MATLAB normally, or run `startup` first. |
| `No requested frame of ... has an input product` | wrong `fabric.in_path`. Pick one of the products the error lists. |
| `Chain 1 succeeded` but no `.mat` | a task-side failure. Look for a `warning` in the task log, above the chain summary. |
| `Undefined function 'ptt.<name>'` | the script was copied out of the clone. Run it from the clone. |
| shell launcher: `no stages/ in ...` | set `FABRIC_ROOT` and `mkdir -p $FABRIC_ROOT/stages/quadpol` (step 5). |
| `mcc: command not found` / `mcc failed to compile` (slurm) | put `/opt/sw/matlab/2024b/bin` on `PATH` before starting MATLAB (step 4). |
| `matlab -batch` sits at 100% CPU after an error | an old copy of `run_fabric.m` with an unconditional `dbstop if error`. Update the clone (`git pull`). |

To update later, run `git -C <clone> pull`. Nothing else needs redeploying.
