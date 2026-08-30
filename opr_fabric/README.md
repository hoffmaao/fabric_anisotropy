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
   - dtau referenced to zero over a coherent BAND just below the surface
     return (removes channel timing/phase biases and the unwrapping
     constant; a single-bin reference injects its own error as a constant
     on every node, and a surface-anchored solve can only absorb a
     constant into its shallowest interval's dlam - the mechanism behind
     the spurious near-surface fabric; see `param.fabric.ref_band_twtt`
     in `fabric.m` and `test/test_copol_surface.m`),
   - coherence-weighted averaging into along-track blocks,
   - inversion through the Maxwell-Garnett firn model for
     piecewise-constant dlam over `num_intervals` depth intervals:
     smoothness-regularized joint solve
     (ptt.invertHorizontalFabricJoint, `inversion = 'joint'`, the module
     default, robust to noisy data) or exact layer stripping
     (ptt.invertHorizontalFabric, `inversion = 'stripping'`); the joint
     solve differences its forward model against the same reference and
     carries the residual reference error as an explicit offset nuisance,
     so the first interval is reference-degenerate and quotable fabric
     starts at interval 2,
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
- The rest of `server/` is single-purpose runners for the SCAR work, the
  comparison against published methods, and the quad-pol scattering-matrix
  work. Like the batches above they keep the `/cresis` DATA roots absolute
  by design, but DERIVE the scratch work root and the repo root from their
  own location (`server/fabric_paths.m`, and `server/fabric_paths.sh` for
  the shell launchers), so no script names a user's scratch. That layout,
  the `FABRIC_ROOT` override and how to stand a work root up are owned by
  `docs/scripts.md` ("Running on the CReSIS machines"). They
  are launched by hand with `matlab -batch`; each one's header states its
  inputs, its launch line and what it concluded - including where the
  conclusion is NEGATIVE, which several of the quad-pol ones are. Runtimes
  are long by nature - tens of minutes a frame for coregistration, longer
  again for a tiled SNAPHU unwrap, over 47-frame surveys - so a batch warns
  and continues per frame rather than letting one bad frame cost the run.
  Each script's header carries its own measured runtime.
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
    can be raised locally. Like `run_negis_fabric.m` below, a frame that
    fails extraction warns and is skipped rather than ending the run, and
    the tail of the log names every frame that went missing - the figures
    say only 'missing extract', which cannot tell a failed frame from one
    that was never requested.
  - `run_negis_fabric.m` - repackages each NEGIS qlook HH/VV pair in its
    `targets` list (by default the shear-margin line and the EastGRIP
    borehole line) into `CSARP_polarimetric` layout and runs the ordinary
    `fabric_task` chain over it, so the inversion is identical to the one
    that produced the other sites rather than a reimplementation. dtau comes
    from the SNAPHU phase, writing `CSARP_fabric_snaphu_negis`. This season
    never went through `polarimetric.m`, so the Goldstein filtering
    (`ptt.goldsteinFilter`) and the SNAPHU unwrap both happen HERE, in the
    repackaging, against the same SNAPHU binary and the byte layouts and
    command line taken verbatim from `polarimetric_task.m`, so this season is
    unwrapped by the same means as the products it is compared against. It
    ran on delta-k
    until the EastGRIP core comparison measured what that cost the joint
    inversion (RMS 0.411 against SNAPHU's 0.282); the header states the
    numbers and where they came from. Two things differ from the Antarctic
    path and both are argued in the header: SNAPHU is run TILED, because a
    single tile on this season's 13501 x ~2000 grid ran over an hour without
    finishing, and an existing product is reused only if it carries a
    NON-EMPTY `snaphu_out_phase`, so a pre-SNAPHU product is rebuilt rather
    than silently falling back to wrapped phase and reporting a 'phase' run
    that never saw one. `targets` is overridable by the caller
    (`matlab -batch "addpath('<code>/opr_fabric/server'); targets = {...};
    run_negis_fabric"`) so another frame set runs through the SAME
    repackaging, cull and surface pick instead of a fork that could drift
    in any of them. A target that fails at
    repackaging, the trace cull or the inversion warns and is skipped so the
    rest of the batch still runs; a SNAPHU failure alone is narrower - the
    product is still saved and still usable via delta-k, so it warns and
    stores an empty phase, which the reuse test then treats as unbuilt.
  - `negis_full_inventory.m` / `negis_full_batch.m` - the inventory measures,
    per frame, whether the qlook product is complex (written with
    `inc_dec = 0`, so the HH/VV phase difference survives) and whether its
    day GPS file covers it. Neither fact is recorded anywhere, and an
    11-hour batch is the wrong place to discover them. The batch is then
    just that frame list driven through `run_negis_fabric.m`: 55 distinct
    frames, the whole season except the 12 real-valued ones of
    `20240618_01`. Frames the earlier 20 km batch already built are reused
    rather than re-filtered and re-unwrapped.
  - `egrip_zeising.m` - the Zeising et al. (2023, TC 17, 1097) phase
    co-registration estimator, implemented so THE METHOD WE ARE COMPARED
    AGAINST runs on our own EastGRIP lines rather than on its published
    numbers from another site. Runs server-side because its lagged product
    s_hh(j) conj(s_vv(j+l)) needs the separate complex channels at full
    range resolution, which the staged zero-lag extracts cannot give.
    Its frame list is the one both sides of the comparison start from.
  - `egrip_collect_inversion.m` - concatenates our own inversion's dlam over
    those same nine lines into one small file to mirror back. `in_name` /
    `out_name` are overridable so the SNAPHU and the superseded delta-k runs
    can both be collected without editing the file, which is what makes the
    switch between them measurable side by side.
  - `run_survey_fabric.m` - the same chain as `run_sections.m` but over
    every frame of a survey, so the result can be MAPPED rather than
    sectioned. Settings are `run_sections.m`'s verbatim, because a map drawn
    at a different block size or gate would disagree with the published
    section along the leg they share. Its header tabulates the azimuth
    diversity of each product set, since orientation needs two legs ~20 deg
    apart over the same ground: Thwaites is listed as ribbons-only rather
    than left out.
  - `run_quadpol_frame.m` - the quad-pol chain (`ptt.quadpolMoments` ->
    `quadpolAzimuth` -> `quadpolFabric`) on one raw four-channel frame.
    Calibration is checked BEFORE the science and reported either way:
    HV/VH reciprocity, and whether the cross-pol ratio oscillates with
    depth or sits at a leakage floor.
  - `quadpol_coreg_frame.m` - coregisters VV, HV and VH onto HH with the OPR
    toolbox, on the window and settings the shipped product itself recorded,
    and first checks the product's own `ref` against the standardphase HH so
    a renamed source cannot silently be different data. Offset fields are
    saved DECIMATED 32x: at full resolution they were 499 MB a frame, larger
    than the images they describe, and carry nothing finer than the tiling
    could resolve.
  - `run_quadpol_pipeline.m` - coregistration and both inversions
    (`ptt.ershadiFabric` and `ptt.quadpolFabricLS`, the second as a
    laterally-segmented frame theta0 pass and a per-block section with
    theta0 pinned to the block's segment profile)
    for one profile in a single pass. Deliberately one script: the
    coregistered images are ~370 MB a frame, so a split would cost more in
    I/O than the inversion, and would invite the halves disagreeing about
    the window - which is the mistake behind the retracted result in
    03d292e, where the inversion ran on channels coregistration never
    touched. ~71 min on a 6601 x 4346 frame the first time; reruns load
    the coreg_cache (keyed to the product's window and tiling settings
    and to the cached array's OWN trace count, so a changed cull rule
    cannot reuse a cache whose trace axis no longer matches) and take
    minutes. `coreg_only = true` stops after building the cache,
    and the 2022/2023 seasons - which shipped the polarimetric product
    only in its SNAPHU-unwrapped form - fall back to
    `CSARP_polarimetric_unwrap` for the window and coregistration
    settings. The FRAME PASS - the antenna-frame pedestal (an
    instrument constant, one estimate per frame), the frame-pooled
    theta0(z) profile, and the geographic-frame machinery for CURVING
    frames (on a curve the fabric rotates through the antennas, so
    antenna-frame moment averaging smears it; the p95 heading-spread
    dispatch trips on GPS jitter by design, so most frames take the
    geographic path, and `test/test_quadpol_curved.m` is its
    regression) - lives in `ptt.quadpolFrameTheta`, whose header owns
    the detail. It re-fits theta0 per ~2 km along-track segment
    (`seg_len_m` overridable, default 2000 m; frames shorter than two
    segments keep the single frame fit) so lateral fabric variation is
    resolved instead of pooled away - the pooled handoff was the
    two-pass design's single point of failure on laterally-varying
    frames (`test/test_egrip_blocks.m` pins what the pooled handoff
    costs, `test/test_quadpol_segmented.m` the rescue). Each section block
    inherits its segment's geographic profile through
    `ptt.thetaProfileAt` and converts it through its OWN heading; the
    per-segment fits are saved beside the frame fields as
    `ls_theta_seg`/`ls_q_seg`/`ls_dlam_seg`/`ls_resid_seg`/`ls_seg_x`/
    `ls_nseg`.
    `z_max` (default 1500 m) is overridable like `day_seg`/`frm` for
    special runs that need the full record; any non-default depth gets a
    `_z<depth>` suffix on both the coreg cache and the section output
    name, so such runs can never clobber the standard batch products.
    Where NO polarimetric product exists at all - the 2024 Greenland
    (EastGRIP) season ships only `CSARP_qlook_{HH,VV,HV,VH}`, with
    all-NaN Surface and a third of the traces stationary - the pipeline
    runs in QLOOK MODE, supplying each missing piece explicitly (the
    header argues every choice): coregistration defaults copied from
    what the Antarctic products recorded, the range window from a
    leading-edge surface pick (the `negis_interferogram.m` rule) plus
    `z_max`, a stationary-trace cull BEFORE any block geometry, and a
    loud error on a real-valued channel (segment 20240618_01 is a
    real-only setup day). Antarctic behavior is unchanged.
    Shipped coordinates that place the frame OFF THE ICE SHEET (below
    60 deg of latitude on most of its positioned traces - the 2024
    Greenland qlook products put the traverse in the Gulf of Guinea) are
    rebuilt by interpolating `CSARP_reference_trajectory/ref_<seg>.mat`
    onto the frame's own GPS times, because along-track distance and
    track azimuth - hence block geometry and every `theta0_geo` - are
    computed from them. The rebuild is not optional: a missing reference,
    or one a frame overhangs by more than the one-second tolerance the
    range check accepts, errors the frame rather than positioning it
    wrongly - inside that tolerance the query is clamped onto the end, so
    a second of slop cannot abort a frame on a gap that does not exist.
    The trigger is the defect, not the product family, so a qlook season
    whose coordinates are sound is left alone.
    The RANGE WINDOW and the coregistration TILING are both fitted to
    the frame that exists: a recorded `max_rbin` past the end of the
    data is clamped (the thin-ice seasons record windows from deeper
    configurations, and one 19 km line was being discarded over it), and
    a frame too short to carry the four along-track tiles
    `coregistration` needs gets a proportionally finer tiling instead of
    failing. The tiling is fitted AFTER the cull, on the trace count
    coregistration is actually handed; frames that already satisfy the
    condition are untouched and bit-identical. The LS dlam
    search ceiling is site-aware and overridable per run via `dlam_max`:
    the estimator default 0.25 for the Antarctic sites, 0.45 in qlook
    mode because the EGRIP contrast (~0.3) sits above the default cap
    (`test/test_egrip_cap.m` pins both sides). `nblk_tr` is likewise
    overridable, and EVERY season auto-defaults it to Ridge A's ~125 m
    block LENGTH from the measured trace spacing, not to a trace count:
    125 TRACES is 125 m only where traces sit ~1 m apart, so at
    EastGRIP's ~2.8 m and Eastwind's 2.43 m spacing it would average
    over 2-3x the ice every other site does, which block-level
    comparisons between surveys cannot survive. The re-cut carries a
    25% length tolerance so the ~1 m seasons keep the exact block
    boundaries their finished frames were cut on.
  - `run_season_*.m` over `run_season_frame.m` - declarative per-season
    fact sheets for that pipeline: each season script declares what is
    different about its acquisition (site root, plus the overrides
    auto-detection cannot know, like EastGRIP's `dlam_max` 0.45) and
    carries the hard-won season knowledge in its header - pedestal sign
    family, real-only setup days, product quirks. The deliberately
    behavior-free `run_season_frame.m` applies the declaration to the
    pipeline's variable interface; variables set at the call site
    always win, so special runs stay one-liners
    (`z_max=4000; run_season_ridge_a`).
  - `coreg_batch.sh` / `coreg_one.sh` and `invert_batch.sh` /
    `invert_one.sh` - survey-scale drivers for that pipeline on the
    shared node: the coreg batch builds the caches (the one genuinely
    expensive pass), and the invert batch FOLLOWS it, picking frames up
    as their caches appear and ending only when no work remains and the
    cache builder has exited. Idempotent skip-if-done workers with mkdir
    locks, nice and 8-thread caps, and a per-frame retry cap; each
    header states its usage and nohup launch line.
  - `site_chain.sh` - queues additional survey sites through that
    coreg-batch/invert-follower pair end to end, one site at a time,
    starting only after every currently running batch has gone fully
    quiet; its header carries the launch line.
  - `egrip_chain.sh` - the EastGRIP season driver for the pipeline's
    qlook mode: gates on the already-revalidated section of the
    validation frame, and only past the gate batches every complex frame
    of the season with the same idempotent lock/retry discipline as the
    other drivers. The gate stopping the batch is the feature: a failing
    estimator must not burn days of compute. What it gates on was
    re-stated once the residual floor was diagnosed as the 2psi
    chan_equal-family calibration signature rather than an estimator
    failure - finite fraction > 0.60 AND block self-consistency (median
    adjacent-block dlam step < 0.05) over 300-1100 m, with the residual
    RECORDED as the documented site floor instead of gated on. The
    script's header carries the measured values behind that decision.
  - `extract_sweep_egrip.m` - pulls the measured C(psi,z) and P(psi,z)
    azimuth sweeps for that validation frame, for exactly that open
    diagnosis: the suspect is anisotropic reflectivity from the EGRIP
    girdle, an even cos-2psi structure in the co-pol power that the LS
    nuisance basis cannot absorb, and that hypothesis is testable in
    Phh(psi) directly.
  - `mech_experiment_009.m` - split-sample tests (interleaved azimuth
    halves, interleaved trace halves, reflectivity correlation) of the
    LS profile's residual 50-150 m depth wiggles, run from the frame-009
    coreg cache; established they are medium structure shared between
    estimators, not estimator noise.
  - `run_quadpol_survey.m` / `run_ershadi_survey.m` - the same two
    inversions over every frame of a survey, the second in 200-trace blocks
    of near-constant heading because the published method assumes a
    stationary sounding and 45 of 47 Ridge A frames wander by more than
    5 deg. Both exist for the heading test: if orientation really is
    recoverable from a single line, theta expressed geographically must not
    depend on which way the vehicle drove. A raster survey supplies that
    control for free. They also record the calibration diagnostics per frame,
    which is what separates a fixed antenna-isolation pedestal from
    depth- and site-varying volume scattering - one frame cannot.
  - `quadpol_coherence_test.m` / `quadpol_coreg_check.m` - the two
    single-frame measurements that settled why |C_HHVV| looked too low to
    pass the Ershadi et al. gate. The first separates genuine along-track
    decorrelation from destruction by the 200-trace average; the second
    reads the coherence of the already-coregistered `ref`/`sec_reg` pair
    straight out of the shipped product, which is the number the raw
    standardphase estimate has to be compared against.

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

The quad-pol side of `+ptt` has its own round-trip tests beside it, each run
the same way:

- `test_quadpol.m` - synthetic column with a known principal azimuth and a
  depth-growing contrast through `ptt.quadpolMoments` -> `quadpolAzimuth` ->
  `quadpolFabric`, asserting the two properties a co-polarized pair does not
  have: orientation recovered from a SINGLE antenna azimuth, and a contrast
  that is the true lam_max - lam_min rather than the projection, so it does
  not vary with the azimuth it was measured from.
- `test_ershadi.m` - `ptt.ershadiFabric` on that same synthetic and the same
  truth, so any difference between the two implementations is theirs rather
  than the test's. Passes `deramped = false` because the synthetic is a model
  and not radar data, which is the distinction the paper itself draws.
- `test_rotations.m` - `ptt.rotateMoments` against `ptt.rotatePolarization`.
  The 4x4 shortcut is only sound if it is the same linear map, and the two
  are written out independently, so a transcription slip would rotate the
  survey's geographic average by the wrong angle and surface only as an
  inflated circular spread - the very quantity the heading test keys on.
- `test_quadpol_ls.m` - `ptt.quadpolFabricLS` on the same truth, PLUS the
  failure mode that motivated it: an antenna-fixed reciprocal leakage term
  added to the cross-polarized channels at the level the real system shows.
  Under that leakage ershadiFabric locks (theta error ~40 deg, dlam
  collapsed to the cos-projection; reported by the test, not asserted)
  while the LS fit must stay within a few degrees and a few percent. Also
  asserts the isotropic case reports dlam ~0 with theta0 ABSTAINING - the
  folded-noise floor and the leakage-shaped axis are both regressions this
  estimator exists to avoid - and that the theta0-pinned two-pass path the
  pipeline section uses reproduces the free fit's contrast.
- `test_quadpol_curved.m` - the curving-line adaptation
  `run_quadpol_pipeline.m` dispatches to. A synthetic 90-deg arc with the
  fabric fixed geographically: the standard antenna-frame path must
  REPRODUCE the smearing collapse of dlam (asserted - the observed symptom
  on turning profiles), while the geographic-frame path - sub-block
  moments rotated north-referenced, the survey-calibrated antenna-fixed
  pedestal pinned as a precomputed field over the measured heading
  distribution - must recover axis and contrast, with and without that
  pedestal.
- `test_egrip_cap.m` - the dlam search ceiling (its header carries the
  `matlab -batch` line, as does the next one). The estimator's default
  `dlam_max` = 0.25 sits BELOW the EGRIP core's ~0.2-0.35 contrast, and
  a window whose true rate exceeds the cap rails and abstains or locks
  an aliased branch: at the default a constant 0.32 synthetic fabric
  abstains on 100% of windows, and at the pipeline's qlook-mode 0.45 it
  is recovered exactly at both 125- and 14-trace blocks. This is the
  regression for the site-aware `DLAM_MAX` in `run_quadpol_pipeline.m`.
- `test_egrip_blocks.m` - the two-pass coupling, on the season's REAL
  2.83 m trace spacing (from the `CSARP_reference_trajectory` rebuild;
  the shipped qlook coordinates inflated it 3.3x). Three verdicts. The
  pooled frame pass ABSTAINS rather than lies as lateral dlam variation
  grows - over 0.010 -> 0.020 dlam/km coverage in the 300-1100 m band
  falls 35% -> 5% and live windows over the whole profile fall 16 -> 8 of
  34, while the band axis error stays within 6.8 deg - and the test
  asserts the pair, so no rate may combine a profile the pipeline would
  hand on with an inaccurate axis. The two halves deliberately take
  different window populations: liveness over the full grid, because that
  is what `ptt.quadpolFrameTheta` gates its `min_seg_windows` on, and the
  axis error over the 300-1100 m band, because `sigma = gpd*ddlam*z` grows
  with depth and shallow windows would dilute a deep-only lie below the
  bound. A rate whose band holds too few live windows to average is exempt
  from the axis bound - that is abstention, not a lie - so the mildest
  rate is required to be both handed on and scored, which is what keeps
  the bound from going hollow. The pooled handoff is nonetheless load-bearing: at a
  modest 0.010 dlam/km the frame axis is only ~3 deg off over that same
  300-1100 m band, yet blocks
  handed it recover dlam ~12x worse than blocks handed the true axis -
  which is what the segmented frame pass (`ptt.quadpolFrameTheta`)
  exists to fix. And block length has margin: with the axis held true the
  collapse is governed by `sigma = gpd*ddlam*z` alone, so `L_max`
  scales as 1/L and the 125 m re-cut buys exactly the length ratio
  354/125 = 2.83x - 14.5x headroom over EastGRIP's expected ~0.01
  dlam/km instead of 5.1x. Its header records why it was rebuilt on the
  real geometry, which invalidated the block-size premise it was written
  for rather than just rescaling it.
- `test_quadpol_segmented.m` - the segmented frame pass
  (`ptt.quadpolFrameTheta` + `ptt.thetaProfileAt`) against the two
  properties the design must have: a SUPERSET (a laterally-uniform frame
  reproduces the pooled architecture's block dlam) and the RESCUE (its
  own killer ramp - 0.10 dlam along the frame, enough to kill the pooled
  pass outright rather than merely degrade it - recovers: block dlam error
  0.225 -> 0.002, correlation 1.00 - and a mid-frame 45 deg axis step,
  the shear-margin case, resolves to under 3 deg in the pure segments
  with block dlam unbiased).

`test/test_copol_surface.m` covers the co-polarized chain's surface
reference at the solver level: the OLD surface-anchored call must
REPRODUCE the spurious near-surface dlam ~ 0.06 that 0.3 ns of single-bin
reference error manufactures, and the call with `obs.zref` must recover
the offset as the nuisance, honor the isotropic cap from the first
quotable interval down, and stay unbiased at depth.

`test/test_goldstein.m` covers the unwrapping side, `ptt.goldsteinFilter`
(its header carries the `matlab -batch` line, as the quad-pol tests do). It
measures RESIDUES rather than asserting on the filtered waveform, because
residues are what SNAPHU pays for in branch cuts and are the whole reason
`run_negis_fabric.m` filters before it unwraps - a filter that reduced
nothing would leave the unwrap as slow as it was and nothing downstream
would say so. The residue drop is paired with a phase-error check, since
flattening the signal would cut residues too. It also pins the properties
the caller depends on: alpha = 0 is exactly the identity (so the
overlap-add and taper normalization reconstruct the input rather than
approximate it), single in gives single out for the float32 SNAPHU file, a
grid smaller than one window warns and returns the input instead of zeros
that would read as a coherence collapse stages later, and an alpha outside
[0,1] is rejected rather than silently clamped.
