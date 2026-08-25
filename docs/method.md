# Method code (`+ptt` package)

Reference detail for every function in the `+ptt` package. Each function
carries a full header of its own (`help ptt.<name>` in MATLAB); this page
is the map and the reasoning that connects them.

Coordinate convention: z is height above the bed, zhat = z/H in [0, 1].
Units: meters and nanoseconds in the theory functions; the interferometric
chain uses seconds to match OPR products (each function's header states its
units).

- `ptt.columnProfiles` - closed-form firn-ice column: Herron-Langway
  density, power-law bubble eccentricity, fabric eigenvalue profiles, and
  the resulting eigenpermittivities/velocities/slownesses via anisotropic
  Maxwell-Garnett mixing (paper eqs. 2.5, 3.1-3.5, 4.2).
- `ptt.twttDifference` - forward model of the TWTT difference t_xz - t_y
  for a CMP geometry (half-offset L, reflector height z), including
  refractive ray bending and the refractive shadow zone (eqs. 2.11-2.13, A1).
- `ptt.invertFabric` - full 5-parameter CMP inversion
  (zhat_bco, lam_x_sfc, lam_x_bed, lam_z_sfc, lam_z_bed) via SQP (eq. 4.1+).
  Needs varying-offset data.
- `ptt.invertHorizontalFabric` - common-offset workflow: layer-stripping
  inversion for the horizontal fabric contrast lam_x - lam_y per reflector
  interval. This is the method for standard accumulation-radar profiling,
  where the Tx-Rx offset is fixed; lam_z and BCO depth must be assumed but
  the contrast is insensitive to them at small offsets.
- `ptt.invertHorizontalFabricJoint` - smoothness-regularized joint variant
  of the common-offset inversion: solves all depth intervals
  simultaneously (Gauss-Newton with a first-difference Tikhonov penalty,
  coherence-weighted misfit), robust where per-interval dtau increments
  are below the noise level, at the cost of some depth resolution. Given
  `obs.zref` (the height the dtau field was zeroed at) it differences the
  forward model against that reference and estimates the residual
  reference error as a near-free constant-offset nuisance, flagging the
  shallowest interval reference-degenerate - the fix for the single-bin
  reference error that manufactured a spurious near-surface dlam
  (`opr_fabric/test/test_copol_surface.m` reproduces both sides). The
  fabric products save the flag as `dlam_ref_degenerate` and the
  recovered offset as `ref_offset`, so consumers read the flag instead
  of hard-coding "skip interval 1"; the figure scripts consume it
  through the `dropped_intervals` predicate in
  `scripts/figures/fabric_qc.py`.
- Interferometric processing chain (pure numerics, no OPR dependencies;
  options structs use the same field names as the OPR fabric worksheet):
  - `ptt.goldsteinFilter` - Goldstein-Werner adaptive spectral filter of the
    COMPLEX interferogram, run before unwrapping. Not cosmetic: SNAPHU's
    cost is set by residues, which noise-driven excursions between adjacent
    pixels create, so filtering first is what keeps the solution from being
    dominated by branch cuts. Never applied to the wrapped phase, which
    would average across the +-pi branch cut and invent fringes. Defaults
    match `scripts/prototypes/branch_cuts.py`, where the alpha/window/stride
    choice was measured against residue counts.
  - `ptt.blendTraveltime` - dtau map from interferogram phase blended with
    coregistration offsets (surface referencing, sign detection, fringe
    ambiguity resolution)
  - `ptt.deltakTraveltime` - absolute dtau map by split-spectrum delta-k
    over the ref/sec SLC spectra (coarse-to-fine sub-band ladder, no
    phase unwrapping; alternative to the phase/coreg blend)
  - `ptt.surfaceReference` - the surface-referencing convention
    (coherence mask, reference bins) shared by the dtau estimators
  - `ptt.imgCombSeam` - fast-time mask of the waveform image-combination
    seams, where the trace is a crossfade of a short and a long pulse
    rather than an ice property. The boundaries are read per frame from
    the product's own `param.array.img_comb` / `img_comb_mult` (0.9 us
    after the surface return with a 0.1 us image-1 guard on the accum3
    grids, 8 us with a 1 us guard on the deeper settings, none on
    single-image frames) instead of being hardcoded, placed with the same
    surface-relative formula `img_combine.m` uses, and the masked band
    runs forward from the boundary because that is the way OPR blends.
    Applied by `ptt.surfaceReference`, so every dtau estimator inherits
    it; disable with `param.fabric.seam_mask_en = false` or retune the
    guard with `param.fabric.seam_mask_win`. Delta-k
    gets a wider band: it resolves its ladder integers from a tau_A
    smoothed over analysis cells, so every cell whose smoothing window
    touched the seam is dropped too (`ptt.deltakDefaults` is the one
    definition of that cell geometry, and `ptt.deltakTraveltime` keeps
    those cells out of its shallow referencing median as well). The two
    bands cost very different amounts of record - measured on the
    6601-bin accum3 frame 20250108_02_009 at the default guard:

    | `dtau_source`    | 0.9 us / 0.1 us blend        | 8 us / 1 us blend              |
    | ---------------- | ---------------------------- | ------------------------------ |
    | `phase`, `coreg` | 0.8-1.1 us, 90 bins (1.4%)   | 7-10 us, 900 bins (13.6%)      |
    | `deltak`         | 0.5-1.4 us, 270 bins (4.1%)  | 6.7-10.3 us, 1080 bins (16.4%) |

    Budget a delta-k run against the `deltak` row: the reach triples the
    cost on the 0.9 us settings.
  - `ptt.blockAverage` - coherence-weighted along-track block averaging
    with per-block fringe correction
  - `ptt.invertBlocks` - twtt-to-depth mapping and per-block inversion
    (exact layer stripping or the regularized joint solve, selected by
    `opts.inversion`). Where the seam mask (or an incoherent run) leaves a
    gap, the node's dtau is interpolated across it so the joint solve
    keeps a continuous chain, and the node is flagged in
    `inv.interpolated`, saved as `dlam_interpolated` next to
    `dlam_quality`/`dlam_clipped`. Layer stripping differences consecutive
    nodes, so a fabricated node corrupts its own interval and the one
    below it; `scripts/figures/fabric_qc.py` is the single definition of
    that predicate, and the figure scripts drop both the way they drop
    pegged ones - the unmasked coherence in `dlam_quality` cannot tell
    them apart from measured nodes (`deltak_vs_joint.py` compares two
    chains, so it drops an interval flagged in either one).
  - `ptt.twttDepthMap` - vertical twtt vs depth from the column model
- Quad-polarimetric chain (the full 2x2 scattering matrix, where a
  co-polarized HH/VV pair measures one azimuth and cannot be rotated at
  all - so orientation is in principle recoverable along a single line
  instead of only at survey crossings):
  - `ptt.rotatePolarization` - the response at a synthetic antenna azimuth,
    T = R' S R. Reciprocity (S_hv == S_vh) is NOT assumed, so finite antenna
    isolation stays visible downstream as a calibration diagnostic instead of
    being averaged away here.
  - `ptt.quadpolMoments` - the [Nt x 4 x 4] Hermitian moment matrix over
    (hh, vv, hv, vh). Every observable used here is quadratic in S and
    rotation is linear, so a whole azimuth sweep is a fixed trig combination
    of those 16 numbers per range bin rather than a re-synthesis of the grid.
  - `ptt.rotateMoments` - the same rotation applied to a moment matrix,
    Mr = W M W', for a 4x4 multiply per block instead of another pass over
    the images. Needed because trace averaging is only valid in a frame the
    signal is stationary in: a vehicle-fixed antenna frame turns with the
    heading, so averaging there biases the answer toward an antenna-fixed
    one by itself. Rotating each block into a common GEOGRAPHIC frame first
    inverts that preference, and running both is what separates real
    leakage from an artifact of the averaging.
  - `ptt.quadpolAzimuth` - the synthetic sweep of co-pol power, cross-pol
    power and the HH-VV coherence from that moment matrix. For a
    birefringent column the cross-pol power factors exactly as
    sin^2(2(psi - theta)) sin^2(delta(z)/2), so azimuth fixes orientation
    and depth fixes birefringence independently.
  - `ptt.quadpolFabric` - inverts the sweep for the horizontal principal
    azimuth (projected onto the cos/sin 4psi harmonic, not an argmin, so a
    near-isotropic layer reports a flat sweep instead of a confident angle)
    plus two deliberately independent contrast estimates: the coherence
    phase gradient, and the cross-polarized node spacing.
  - `ptt.ershadiFabric` - Ershadi et al. (2022, TC 16, 1719) implemented
    faithfully and kept SEPARATE from `ptt.quadpolFabric`, so the two can be
    run on the same data and scored against each other. Every departure from
    that file is a place the paper differs and is marked (E1)-(E5) in the
    header, including why the birefringence coefficient carries sqrt(eps')
    rather than the printed eps'. Note this is the paper's DIRECT stage
    only; their published profiles additionally pass through a constrained
    nonlinear fit of the Fujita model (their Sect. 3.5) drawn on a
    continuous depth parameterization, which is why they show no gaps.
  - `ptt.quadpolFabricLS` - the least-squares replacement for the direct
    chain's weak point. ershadiFabric takes the fabric axis from the
    cross-polarized minimum; on this system that minimum is antenna-locked
    (89.6 +- 1.9 deg across 1790 Ridge A blocks spanning all headings,
    because the cross-pol channels sit on a flat ~-3.6 dB instrument
    pedestal), so Psi gets evaluated ~25 deg off-axis: dlam scales down by
    cos 2*offset and the odd-pi coherence nulls turn into banded dropouts.
    This estimator instead fits the measured complex coherence C(psi, z)
    over ALL azimuths to the exact single-column birefringence model in a
    sliding depth window (theta0, delta0, ddelta/dz, plus a closed-form
    coherence scale), so the axis comes from the coherence field the
    co-pol channels dominate, the nulls are modelled rather than gated,
    dlam is signed with no folded-noise floor, and theta0 abstains where
    the birefringence is unresolvable. Two-pass use: a frame theta0
    pass, laterally segmented (`ptt.quadpolFrameTheta` below), then
    per-block dlam with theta0 pinned.
    KNOWN SYSTEMATIC (measured 11 Aug 2026, unresolved, and CONFIRMED not
    to be the heading-wrap bug - see below): dlam carries a
    heading-family bias of order 0.01 (~17%) - same-ice crossing
    pairs between the N-S lines and the rows at Ridge A differ by +0.009
    median with 93% sign consistency (26 of n = 28 line crossings, one
    pair each, `scripts/figures/crossing_pairs.py`), consistent with
    residual pedestal-fabric coupling that grows with the axis-to-antenna
    angle (N-S lines hold the axis ~6.5 deg from an antenna, rows ~18.5,
    the NW-SE connectors ~36.5 - nearest the 45-deg degenerate geometry).
    The combinations involving the NW-SE connectors are not quotable in
    either direction: their crossings are one cluster of end-of-row stubs
    and survive deduplication as n = 2 (NW-SE minus row) and n = 3 (N-S
    minus NW-SE). Orientation is unaffected. Until fixed, quote deep
    dlam with a +-0.005-0.01 family systematic; near-aligned lines are
    least coupled and read high (0.067-0.074 deep at Ridge A). THE
    PLANNED FIX is raw-channel calibration (chan_equal): fit the complex
    HV/VH leakage once per system at the channel level and correct
    BEFORE synthesis, removing the pedestal at its source; the per-frame
    ls_pedestal record and the crossing-pair test (paired same-ice
    difference -> 0) are its inputs and acceptance criterion.
    RULED OUT as the cause: the 0/180 heading-wrap interpolation bug
    fixed in 0e793d5 lived exactly on the N-S family, so all 17 N-S and
    curved Ridge A frames were reswept with the fixed code, and all 17
    still have a surviving pre-fix copy to be compared block for block
    against. Over that whole set the deep dlam moved by a median of
    +0.00000 and a mean of +0.00000 across 1574 blocks (p10/p90
    -0.00003/+0.00004, 99.9% under 0.001, two blocks over it, none over
    0.005, largest single block 0.0048), and frame theta0 was unchanged
    to 0.12 deg. No block's usability flipped either way across the whole
    2025 series (gain 0, lost 0), so nothing sits outside that n: the
    null is not hiding a block that gained or lost a usable contrast.
    The bug was real but the per-block circular mean over 125-200 traces
    absorbed the handful of mis-interpolated headings it produced. The
    family systematic is therefore instrument physics, not that defect,
    and chan_equal is the right thing to build against it.
    `scripts/compare_sections.py` is how that measurement is made and
    will re-make it, but only where both product generations exist on
    disk: the pre-fix inputs survive solely in a local mirror
    (`all_sections_final.tar` under `SCAR_DATA`), which is not in this
    repo and cannot be regenerated from it, so the numbers above are a
    recorded one-time measurement, not something CI can re-derive.
    On the timeline, which reads as contradictory until the zones are
    lined up: 0e793d5 is stamped 09:05 +0200, i.e. 02:05 US/Central, and
    the resweep ran 03:09-10:40 CDT, so the fixed code was in place for
    all of it and every 2025-series member of the mirror (9 Aug 23:40 to
    10 Aug 02:09) is pre-fix.
  - `ptt.quadpolFrameTheta` - the frame pass of that two-pass use,
    laterally segmented: one antenna-frame pedestal per frame (an
    instrument constant) plus a geographic theta0(z) profile re-fitted
    per ~2 km along-track segment, because one pooled theta0 handoff is
    the two-pass design's single point of failure on laterally-varying
    frames. Curving frames take the geographic sub-block machinery
    validated on a 90-deg arc; the header owns the detail, and
    `opr_fabric/test/test_quadpol_segmented.m` /
    `test_quadpol_curved.m` are its regressions.
  - `ptt.thetaProfileAt` - the per-block handoff from that frame pass:
    the segment profiles interpolated at a block's along-track position
    on the doubled-angle phasor, with robust q-weighted end rows because
    the estimator clamps out-of-range depth windows to the terminal row.
  - `ptt.coregisterChannels` - aligns every channel onto the reference by
    CALLING the OPR toolbox `coregistration`, deliberately with no
    implementation of its own (it errors if the toolbox is absent), so the
    quad-pol path is aligned by the same code and defaults that built the
    shipped `CSARP_polarimetric` products. A hand-rolled global-delay fit
    was tried first and was worse on both pairs; the header has the numbers.

## Cross-survey consistency

Audited over 130 frames and seven surveys (24 Aug 2026), because
block-level and depth-level quantities are only comparable between sites
if every step is computed the same way. What follows is what that audit
established - the checks that PASSED are recorded too, so they are not
re-litigated.

**Consistent, verified:**

- **Depth axis.** One formula everywhere, `z = (t - t_surf)*c/(2*sqrt(3.171))`,
  with the surface from the product's own `Surface` (or, in qlook mode
  where that field is NaN, a leading-edge pick at 5% of the trace peak).
- **Permittivity.** The 3.171 in that formula is the same `eps_bar` that
  `ptt.constants` gives `ptt.quadpolFabricLS` for the phase-rate to dlam
  conversion, so depth and contrast use one ice model. (Worth restating
  because they are declared in different files and could drift apart.)
- **Estimator settings.** Azimuth grid, short/fit windows, CRB weighting
  and the two-pass structure take their defaults at every site.
- **Segment length.** 2 km everywhere.
- **Block length.** ~125 m everywhere, enforced in every season by
  re-cutting the trace count where the measured spacing demands it
  (`opr_fabric/server/run_quadpol_pipeline.m` owns the rule and the
  tolerance that keeps the ~1 m sites' existing block boundaries).

**Known differences, deliberate:**

- **dlam search ceiling.** 0.45 in qlook mode against 0.25 elsewhere,
  because the EGRIP core sits near 0.3. Checked for bias: no site rails
  against its ceiling - p99.9 is below the cap at all seven, and under
  0.1% of cells come within 2% of it - so the difference does not
  distort any comparison.
- **Coregistration tiling.** Each season is aligned with the settings its
  own product recorded (Eastwind `Tt 51 Tx 101 ov 25/25`, the others
  `101/301/50/150`), which is what reproduces the shipped `sec_reg`.
  Short frames additionally get a finer tiling so they can be tiled at
  all.

**Caveats to carry:**

- **No firn correction in the quad-pol depth axis.** It uses solid-ice
  velocity from the surface down, so depths in the upper ~100 m are
  compressed by roughly 10 m. This is consistent BETWEEN surveys, so it
  biases no site comparison, but it does mean quad-pol depths and the
  co-polarized chain - which models firn densification through
  `ptt.columnProfiles` - disagree slightly near the surface. `z` is saved
  per frame, so a remap is possible after the fact if the two ever need
  to be overlaid precisely.
- **The per-frame pedestal is unreliable on low-coherence frames**, so
  the `ls_pedestal` record needs a quality cut before it is used as
  chan_equal calibration input. One season spans a3 -0.264 to +0.438
  across its sites while one site inside it holds sd 0.011. A frame whose
  frame-pass pedestal did not converge saves `ls_pedestal` as exactly
  `[0 0 0]`, which is the marker to drop, not a measurement. The marker is
  `ptt.pedestalMarker()` and the test for it is `ptt.pedestalFailed`,
  which every consumer uses instead of a finiteness check - the marker is
  finite, so `isfinite` reads it as a measurement. The pipeline itself
  tests it that way: a marked frame fits each section block's own
  pedestal rather than anchoring every block at zero.
