# opr_fabric: fabric inversion module for the OPR toolbox

Drop-in processing module for the OPR toolbox
(https://gitlab.com/openpolarradar/opr) that turns the existing
`CSARP_polarimetric` product into profiles of the horizontal ice-fabric
contrast dlam = lam_x - lam_y versus depth and along-track position,
using the theory of Rathmann (2026) implemented in the `+ptt` package
(one directory up).

## Pipeline

1. Per-polarization echograms (existing OPR flow): qlook/sar/array
   processing of each polarization channel into separate products, e.g.
   `CSARP_standardphase_HH`, `..._VV`, `..._HV`, `..._VH` (see
   `run_polarimetric.m` in the toolbox; the 2024_Antarctica_Ground2 season
   is set up this way).
2. `polarimetric.m` (existing): synthesizes a rotated basis
   (`synth_rot_deg`), coregisters ref/sec, forms the multilooked
   interferogram + coherence, optionally SNAPHU-unwraps the phase. Run
   with `coregistration.en = true` and ideally `snaphu_en = true`.
3. `fabric.m` / `fabric_task.m` (this module): `fabric_task.m` is a thin
   OPR adapter (file discovery, product loading, output/figure
   conventions) around the pure numerical chain in `+ptt` -
   `ptt.blendTraveltime` (or `ptt.deltakTraveltime`) -> `ptt.blockAverage`
   -> `ptt.invertBlocks` - which reads its options directly from the
   `param.fabric` struct and is equally callable from standalone scripts
   and tests. Per frame,
   - traveltime differences dtau(twtt, x) = t_sec - t_ref from one of
     three estimators (`param.fabric.dtau_source`):
     - `'phase'` (default): blend of coregistration `row_offset * dt`
       (fixes the sign and the integer 1/fc fringe ambiguity) with the
       interferogram phase (sub-ns precision; SNAPHU-unwrapped when
       present, otherwise wrapped phase with per-pixel fringe resolution),
     - `'coreg'`: row offsets alone (detected-power products),
     - `'deltak'`: split-spectrum ladder over the ref/sec SLC spectra
       (`ptt.deltakTraveltime`): sub-band interferograms cross-multiplied
       per pixel then multilooked, coarse-to-fine in synthetic wavelength,
       giving ABSOLUTE dtau with no phase unwrapping and no fringe
       blending. Needs `ref` and the unregistered `sec` in the product
       (~2x task memory). CAUTION: `sec_reg` is envelope-shifted by
       coregistration and carries no group delay - delta-k must use the
       raw `sec`. Motivated by SNAPHU region errors and coreg-anchored
       blend corrections at the low-coherence Thwaites margin (the two
       failure modes it eliminates); at high-coherence sites it mainly
       serves as an unwrap-free cross-check.
   - dtau referenced to zero just below the surface return (removes
     channel timing/phase biases and the unwrapping constant),
   - coherence-weighted averaging into along-track blocks,
   - inversion through the Maxwell-Garnett firn model for
     piecewise-constant dlam over `num_intervals` depth intervals:
     smoothness-regularized joint solve
     (ptt.invertHorizontalFabricJoint, `inversion = 'joint'`, the module
     default, robust to noisy data) or exact layer stripping
     (ptt.invertHorizontalFabric, `inversion = 'stripping'`),
   - output `CSARP_<out_path>/<day_seg>/Data_*.mat` + overview images
     (`out_path` defaults to `fabric`, but the drivers override it to
     `fabric_joint` so joint-solver results sit beside the earlier
     `CSARP_fabric` stripping outputs for comparison).

## Deployment on the CReSIS servers

- `fabric.m`, `fabric_task.m` -> `opr/matlab/processing/` (or keep on your
  personal path); `run_fabric.m` -> your `run_opr` repo
  (`gRadar.path_override`), edited per season.
- `+ptt` (from the project root) must be on the MATLAB path, e.g. copy to
  `opr/matlab/+ptt` or your `run_opr` repo.
- For compiled cluster modes, add `{'fabric_task.m' 2}` to
  `gRadar.cluster.hidden_depend_funs` in startup.m and re-run
  cluster_compile. `cluster.type = 'debug'` needs none of that.
- Param spreadsheet: add a `fabric` worksheet (row 1 field names, row 2
  type codes, one row per segment matching the `cmd` sheet order), e.g.
  `out_path`(t), `in_path`(t), `fc`(r), `block_size`(r),
  `num_intervals`(r), `ptt.H`(r), `ptt.bco_depth`(r). Enable per segment
  via the `cmd` sheet `generic` column `{{'fabric','fabric'}}`, or just
  drive everything from `run_fabric.m` overrides (current default).
  Note master.m only propagates ctrl_chains for `analysis` generic steps,
  so run via `run_fabric.m` until that one-line change is upstreamed.
- For whole-season batches that bypass master.m and the param spreadsheets
  entirely, `server/run_fabric_scratch.m` runs `fabric_task` (regularized
  joint inversion) over every frame of each input product in its
  per-product season table (EAGER traverse seasons 2022-2026, including
  both independent 2024-25 polarimetric processings) through the
  `test/stubs` opr_* shims, writing `CSARP_<out_path>` outputs
  (`fabric_joint`, or `fabric_joint_jp` for the second 2024-25 processing)
  to a user scratch tree instead of the shared season tree; see its header
  for the table, paths, and the `matlab -batch` launch line.
  `server/run_deltak_scratch.m` is the same kind of batch with
  `dtau_source = 'deltak'`, restricted to the products that carry the
  complex ref/sec SLCs, writing `CSARP_fabric_deltak(_jp)` beside the
  joint results for comparison; see its header.
  `server/run_deltak_stages.m` is a single-frame diagnostic, not a
  batch: it runs all three estimators over one polarimetric frame and
  saves a depth profile for every rung of the delta-k ladder beside the
  phase estimators, to locate where the Ridge A delta-k amplitude is
  lost; `scripts/figures/deltak_stages.py` plots the result (see both
  headers).
- The rest of `server/` is single-purpose runners for the SCAR work. Like
  the batches above they carry absolute mem1 paths by design and are
  launched by hand with `matlab -batch`; each one's header states its
  inputs, its launch line and what it concluded.
  - `run_sections.m` - the shipped high-resolution 2D sections for all
    four profiles (Ridge A, Thwaites, NEGIS, and the EastGRIP borehole
    line) at the retuned settings.
    The per-frame coherence gate is swept and chosen on ADJACENT-BLOCK
    AGREEMENT, not misfit rms - rms is a post-regularization residual that
    always improves with more free parameters and would endorse an
    overfit - and among gates that agree about as well as the best one,
    the deepest-reaching wins, because the figures shade by coherence so
    weak cells arrive visibly faded rather than silently absent. Writes
    the `fabric_sections.mat` both section figures read.
  - `ridge_a_resolution.m` - the sweep of depth resolution, block size and
    regularization that established those settings.
  - `ridge_a_edge_test.m` - is the deep Ridge A downturn physical or a
    domain edge? Sweeps the interval count, the column-model bed depth
    `par.H`, and the record truncation depth (the decisive test: an edge
    effect propagates inward, so values well above the cut must not move
    when the cut moves).
  - `fringe_check.m` - end-to-end validation of the chain against the raw
    fringe density, which is proportional to dlam with no inversion
    involved. Calibrates against `ptt.twttDifference` directly rather than
    a hand-derived constant, on purpose: hand-deriving it drops the
    two-way factor that lives inside `twttDifference` and produces a
    spurious factor of 2 that looks exactly like a real bug.
  - `negis_reposition_tracks.m` - rebuilds the 2024_Greenland_Ground2
    trajectories from the GPS files, since the records radar-time-to-GPS
    sync failed for every segment but 20240618_01.
  - `negis_interferogram.m` - forms and extracts the HH x conj(VV)
    interferograms for the candidate NEGIS frames, plus the lines that
    pass within 2 km of the EastGRIP borehole, which
    `scripts/figures/egrip_azimuthal.py` fits together for the horizontal
    fabric ellipse; its header states why that second group must stay in
    step with the figure's frame list. Culls the traces the
    traverse stopped for BEFORE multilooking (they multilook to high
    coherence with random range phase, so no downstream coherence test
    catches them) and ships the un-normalised power sums so the look count
    can be raised locally.
  - `run_negis_fabric.m` - repackages each NEGIS qlook HH/VV pair in its
    `targets` list (the shear-margin line and the EastGRIP borehole line)
    into `CSARP_polarimetric` layout and runs the ordinary `fabric_task`
    delta-k chain over it, so the inversion is identical to the one that
    produced the other sites rather than a reimplementation. A target that
    loses too many traces to the cull, or whose inversion throws, warns
    and is skipped so the rest of the batch still runs.

## Margin display extracts

For local interferogram QC of the Thwaites line segments
(`scripts/figures/thwaites_interferogram.py`; its docstring says which
segments cross the eastern shear margin), a scratch-side
`extract_margin.py` on mem1 condenses each raw ~1.9 GB
`CSARP_polarimetric` frame into a compact `margin_<day_seg>.mat` display
extract: decimated power in dB for both channels, coherence, wrapped and
SNAPHU-unwrapped phase, coregistration row offsets, and geometry
(Time/Latitude/Longitude). The raw ref/sec complex images never leave
the server; rsync the extracts to `~/data/opr/margin/` locally.

## Ground accum radar (accum3 / EAGER) channel mapping

From the mission defaults and lever_arm.m of recent ground seasons
(2023-2025, config `psc_eager_configHV*`): 2 Tx x 2 Rx colocated crossed
bowties on the sled, rx path 1 = H = ALONG-TRACK polarization, rx path 2 =
V = CROSS-TRACK, zero baseline. The 4 waveforms are 2 Tx pols x 2 pulse
lengths (wf1/wf3 short 0.1-1 us, wf2/wf4 long 1-8 us), so [wf adc] pairs
map to (deep waveforms in bold for fabric work):

| [wf adc] | pol | pulse |    | [wf adc] | pol | pulse |
|---|---|---|---|---|---|---|
| [1 1] | HH | short |    | **[2 1]** | **HH** | long |
| [1 2] | HV | short |    | [2 2] | HV | long |
| [3 1] | VH | short |    | [4 1] | VH | long |
| [3 2] | VV | short |    | **[4 2]** | **VV** | long |

With `synth_rot_deg = 0`, ref = HH and sec = VV, so this module's
dlam = lam_cross-track - lam_along-track.

Season/data caveats to check before interpreting results:
- Phase-preserved products require `array.method = 'standardphase'` (no
  multilook); the default spreadsheets and the public portal's
  CSARP_standard_* use power-detected 'standard' (hence dtau_source =
  'coreg' for those).
- `radar.chan_equal_dB/deg` are all zero in every recent ground season: no
  channel equalization has been applied. Fine for single-pair HH/VV dtau
  (surface referencing absorbs biases), but REQUIRED before trusting
  rotated-basis synthesis (synth_rot_deg ~= 0) or HV/VH use; note
  polarimetric_task.m itself substitutes VH for HV due to a known HV
  amplitude scaling issue.
- Some 2024_Antarctica_Ground2 segments are marked "All polarizations bad
  except VV deep waveform. Failed connector on H channel? Do not process."
  in the cmd sheet notes; 2025_Antarctica_Ground2 notes "polarimetric
  phase unwrapping hangs". Check the notes column per segment.
- There is no `polarimetric` worksheet reader in the toolbox; its params
  (and this module's, unless you add a `fabric` sheet) are set from run
  scripts via opr_set_params.

## Physics caveats

- Common-offset data constrains ONLY the horizontal contrast along the
  synthesized axes: dlam > 0 means more c-axis concentration along the
  secondary (rotated V) axis than the reference (rotated H) axis. The
  vertical eigenvalue lam_z and bubble close-off depth are assumed
  (`param.fabric.ptt`); the synthetic tests show dlam is insensitive to
  those assumptions at small offsets.
- The synthesized basis should be aligned with the horizontal fabric
  principal axes (choose `synth_rot_deg` in polarimetric.m, e.g. from the
  rotation movie or minimum cross-pol energy); misalignment mixes in
  polarization rotation that this scalar-traveltime model does not
  capture. Azimuth scanning via multiple `synth_rot_deg` runs is a
  natural extension.
- With `dtau_source = 'phase'` the inversion takes the depth GRADIENT of
  dtau from the interferogram phase but its absolute LEVEL from
  coregistration: `ptt.blockAverage` shifts each block by a WHOLE number
  of fringes, `round(median(dtau_coreg - dtau)*fc)/fc`. Testing showed
  dlam is invariant to a whole-fringe coregistration offset but sensitive
  to a fractional one - a quarter fringe is enough to move the result.
  This is the likely mechanism behind the weak Thwaites result, where
  coregistration is poorest; correcting it means changing how
  `ptt.invertBlocks` anchors the level, which has not been done.
- Exact layer stripping (`inversion = 'stripping'`) amplifies noise
  between depth intervals; the default joint solve suppresses this with
  its smoothness penalty at the cost of some depth resolution. If
  profiles still oscillate, increase `block_size` / `mlook_window`,
  reduce `num_intervals`, or (joint mode) raise `reg`.

## Test

`test/test_fabric_task.m` builds a synthetic CSARP_polarimetric frame from
a known fabric (with noise, wrong-sign convention, unwrapping constant,
channel timing bias, decaying coherence, and a waveform image-combination
seam declared in `param.array.img_comb` that only the seam mask can
remove), runs the real `fabric_task` with
stubbed OPR support functions (`test/stubs/`), and asserts the inferred
dlam matches the truth in all three `dtau_source` modes - `'phase'`,
`'coreg'`, and `'deltak'` (the last from synthetic band-limited ref/sec
SLCs with the true delay applied spectrally). Runs in MATLAB or Octave:

```sh
docker run --rm --platform linux/amd64 -v "$PWD/../..":/work \
  -w /work/opr_fabric/test gnuoctave/octave:latest \
  octave --no-gui test_fabric_task.m
```
