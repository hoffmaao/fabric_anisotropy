# Scripts

Every script that validates or applies the method code in
[`method.md`](method.md), end to end.

Figure scripts write into the repository's `figs/` directory, which is
gitignored - see the "Figure inputs and outputs" section at the end.

- `scripts/synthetic_experiments.m` - reproduces the paper's synthetic CMP
  experiments 1-3 (noisy forward data, uninformed initial guess).
- `scripts/synthetic_common_offset.m` - synthetic validation of the
  common-offset layer-stripping method.
- `scripts/accum_inversion_template.m` - template for real accumulation
  radar picks (fill in section 1).
- `scripts/extract_swath.py` - run on mem1 (not locally): reduces a full
  rds `CSARP_music3D` tomographic frame to a compact npz that can be
  rsynced off the server - a multilooked uint8 dB cube, deliberately
  quantized because it only ever becomes movie frames, plus the
  look-angle energy reductions that are used quantitatively, kept at
  full precision. Look angles are written in degrees under `theta_deg`.
- `scripts/reduce_ifg_movie.py` - run on mem1 like `extract_swath.py`,
  making the same trade for the wrapped-phase movie: each ~690 MB Ridge A
  `CSARP_polarimetric` frame becomes a ~1.5 MB npz of display-resolution
  uint8 phase and coherence (NaN handled BEFORE the cast, because numpy
  casts NaN to an arbitrary uint8), while the geometry - read
  quantitatively by the map marker and distance axis - stays float.
  `scripts/figures/ridge_a_ifg_movie.py` animates the result locally.
- `scripts/compare_sections.py` - differences two directories of
  `quadpol_section_*.mat` for the same frames: per-block deep dlam
  change over 1150-1500 m (same >= 5 finite cells rule as
  `crossing_pairs.py`) and the frame theta0 change as a doubled-angle
  phasor difference, per frame and pooled. Only blocks usable in BOTH
  generations have a delta, so the blocks whose usability flipped are
  counted beside n in the `gain` and `lost` columns rather than dropped
  silently - a null result must not be able to hide a block that gained
  or lost a usable deep contrast. This is what measured the
  heading-wrap resweep above; it needs both product generations staged
  locally, so it re-makes that measurement only where they are.
- `scripts/figures/` - Python (matplotlib + scipy; cartopy required by
  the map figures) figure scripts that reproduce the analysis figures
  from the `opr_fabric/server/run_fabric_scratch.m` and
  `run_deltak_scratch.m` batch outputs mirrored locally:
  season transects, per-season depth profiles, the 2024-25
  cross-validation of the two independent polarimetric processings, the
  Ridge A azimuthal inversion (per-depth fit over the grid's drive
  headings for the horizontal fabric orientation and principal contrast
  lam1 - lam2), the Thwaites flow-frame analysis (ITS_LIVE velocities
  sampled at each block rotate the measured contrast into the flow frame
  across the eastern shear margin, compared against the Ridge A P(z)
  baseline), the interferogram-chain QC of the Thwaites line segments
  (one figure per segment stacking power, coherence, wrapped/unwrapped
  phase, and coregistration offsets from compact server-side extracts;
  the script's docstring says which segments cross the eastern shear
  margin and which sit in the fast band),
  single-frame scene figures of the eastern-shear-margin crossing and of
  a Ridge A grid leg (map plus stacked along-track panels),
  an all-survey summary (coverage map, median profiles, drive orientations;
  the map uses cartopy if installed, else a plain lat/lon scatter),
  along-profile depth sections of the fabric solution per survey region,
  the delta-k vs joint-chain validation (per-interval dlam scatter and
  difference histogram at Ridge A, printed stats for Thwaites as well),
  the per-stage delta-k ladder diagnostic (`deltak_stages.py`, from the
  single-frame extract written by
  `opr_fabric/server/run_deltak_stages.m`: dtau profiles and retained
  depth increments per ladder rung against the phase estimators, to
  locate where the Ridge A amplitude is lost),
  the matching delta-k eigenvalue-difference depth sections along both
  traverses (the Ridge A panel is flagged unvalidated there because
  delta-k is amplitude-suppressed against the joint chain at that site),
  and zoomed regional maps of each survey area with reference-site markers
  and an Antarctica locator. The maps overlay LIMA/MODIS MOA satellite imagery
  when rasterio and the locally downloaded mosaics are present (download
  commands in `scripts/figures/antarctic_basemap.py`), and fall back to
  coastline-only otherwise. See each script's docstring for usage; the
  inversion chain itself stays MATLAB/Octave.
- The SCAR talk figures are a subset of `scripts/figures/` with their own
  conventions (see "Figure inputs and outputs" below):
  - `scar_two_panel.py` - survey overview + wrapped interferogram for
    Thwaites and Ridge A: the survey over log-scale ITS_LIVE speed with
    all lines white and the focused profile black, red start / white end
    markers repeated on the map inset and at the interferogram's left and
    right edges, no legends, plain black scale bars. Thwaites is framed on
    the WHOLE staged ITS_LIVE window rather than on a box padded round its
    tracks, so the glacier trunk the survey sits beside stays in frame.
    Each site also gets a FABRIC COMPANION slide on the wrapped-phase
    panel's exact geometry - the quad-pol LS dlam section placed on the
    same distance axis by nearest trace, on the same absolute-TWTT axis
    and limits, with the same end markers and colorbar placement, cells
    fading toward grey by their own coherence - so the two slides overlay;
    it prefers a staged deep `_z*` section variant (the overridable-`z_max`
    reruns) so the fabric panel can span the wrapped-phase record.
  - `negis_two_panel.py` - the same slide for the NEGIS onset at EGRIP,
    on a Greenland polar stereographic projection, from the
    `negis_interferogram.m` extract; it can raise the look count locally
    because that extract ships the un-normalised power sums.
  - `egrip_two_panel.py` - the same slide for the EastGRIP borehole,
    windowed on the drill site rather than on the survey (so it needs no
    zoom inset) and with the interferogram plotted against DEPTH, so the
    panel can be set directly beside core measurements. Its docstring says
    why the line drawn is the third-nearest one.
  - `egrip_azimuthal.py` - the horizontal fabric ellipse at EastGRIP from
    the nine phase-preserving lines that pass within 2 km of the borehole
    at different headings (`negis_interferogram.m` extracts them). Fits
    dlam = -P cos 2(alpha - theta) per depth band, so orientation as well
    as strength is recovered without the core's missing azimuth. Works
    from the FRINGE RATE rather than the layer stripping, and accepts a
    band only where its two independent rate estimators agree - they fail
    in opposite directions, which is what limits the solve to above 650 m.
  - `egrip_core_compare.py` - that radar solve against the EastGRIP core's
    own eigenvalues with depth. The core gives magnitudes but no
    eigenvectors, so all three candidate horizontal differences are drawn;
    the line staging, the acceptance rule and the whole per-band solve are
    imported from `egrip_azimuthal.py` rather than restated, so the two
    EastGRIP figures cannot disagree about which lines entered a solve or
    about the P it produced. The horizontal bars are the same solve run on
    each rate estimator on its own, which is not the same thing as
    averaging their two P values - the fit is linear in the estimates, P
    is their norm.
  - `scar_style.py` - the one definition of the styling constants,
    end markers, panel, zoom-inset and continental-locator geometry and
    the track-azimuth/geometry helpers those figures share, so the site
    slides cannot drift apart. Scale bars: the plain bottom-right bar for
    the EPSG:3031 panels is `scale_bar_br` in `antarctic_basemap.py`,
    shared by the two-panel stills and the movie frames (two private
    copies of it had drifted once already); `scar_style` keeps the
    projection-parameterized variant the Greenland-projection figures
    pass their own rounding rule to. Every two-panel figure draws
    that locator - Antarctica for Thwaites and Ridge A, Greenland for
    NEGIS and EastGRIP - from cached Natural Earth land, so it needs no
    network; the study area is the map panel's OWN extent drawn on the
    continent as a black box, rather than a symbol placed near it, and
    only the locator's corner is chosen per figure, by whatever free
    space that panel has. A survey is a small enough fraction of its ice
    sheet to vanish when drawn true to scale, so that box carries a
    legibility floor, and every one of these four sites hits it: read
    the square as "here", not as a scale bar. The `continental_inset`
    docstring gives the numbers. Also used by the EastGRIP solve above
    and the two section figures below.
  - `fabric_sections.py` - per-site 2D dlam sections from
    `opr_fabric/server/run_sections.m`, each cell blended toward grey by
    its node coherence (a diverging ramp cannot reuse the interferogram's
    fade-to-black), with a two-dimensional key for that mapping.
  - `fabric_three_sites.py` - the cross-site depth-profile comparison.
    Magnitudes are NOT directly comparable: dlam is the difference of the
    two horizontal eigenvalues on each profile's own axis and its
    perpendicular, and those azimuths differ per site, which is why both
    figures name them.
- The published-method comparison and quad-pol figures share those same
  conventions (input staged under `SCAR_DATA`, output to the `<out_dir>`
  argument rather than `figs/`) and import `scar_style.py` for them:
  - `egrip_method_compare.py` - our joint inversion against Zeising et al.
    (2023) and against our own bare fringe rate, all three on the SAME nine
    borehole-proximal EastGRIP lines, aperture, surface pick, depth bands and
    azimuthal solve, so the remaining difference is the estimator and not the
    staging or the constants. RMS against the core: 0.282 joint, 0.280 fringe
    rate, 0.356 Zeising. The docstring records that the same joint inversion
    scored 0.411 on delta-k, why the per-frame dtau misfit ROSE when SNAPHU
    replaced it (so misfit alone would have preferred the worse answer), and
    which core eigenvalue pair is horizontal at which depth.
  - `fabric_map.py` - depth-averaged horizontal fabric in map view for every
    survey (Ridge A, Taylor Dome, NEGIS, Eastwind, Thwaites), after Nymand et
    al. (2025) fig. 3: ribbons coloured by each leg's own PROJECTION
    dlam_obs = -P cos 2(alpha - theta), and orientation crosses drawn only in
    grid cells holding two legs far enough apart in azimuth to separate P
    from theta. Thwaites is a single-azimuth transect, so it gets one panel
    rather than an empty second one that would read as a weak fabric.
  - `quadpol_diagnostic.py` - the azimuth sweep itself for one frame, because
    the method's validity is visible there and not in the numbers it returns:
    vertical nulls fix orientation, horizontal nulls fix birefringence, and a
    flat panel means the cross-pol is antenna leakage rather than ice.
  - `quadpol_section.py` - the quad-pol dlam/theta section beside the
    co-polarized one, each panel labelled with what it actually measures
    (projected onto this line's axes vs at the principal axes) rather than
    both with "dlam". Its theta panel is drawn as a DIAGNOSTIC and says so:
    on Ridge A the recovered orientation tracks the antenna frame, at 2.0 deg
    circular spread there against 41.7 geographic. When the pipeline output
    carries the `ptt.quadpolFabricLS` fields it adds a third row - the LS
    coherence-field fit above the published direct chain - on a shared
    color scale; older .mat files without them still draw the two-row
    layout.
  - `quadpol_sites.py` - the one shared definition of survey identity and
    drawing parameters for the multi-site quad-pol maps and depth movies.
    Sites are selected by POSITION, not tag prefix - the 2024 tags alone
    span Thwaites, WAIS Divide and McMurdo, so a prefix map would draw
    three sites on one axis - with Eastwind and McMurdo, which overlap
    spatially, split by SEASON (a 2022 single-spot rotation experiment vs
    a 2024 15 km transect). Each site carries its own depth bands, colour
    ceiling, movie depth range (stopped where the site's median dlam
    crosses the estimator's 0.01 abstention threshold; the measured
    per-window tables are in its comments) and window/step (thin-ice
    sites get 40 m windows at 5 m steps - 150 m is half the ice at a
    300 m shelf, and McMurdo's earlier 450 m cut was detecting the ice
    base, not a fabric limit).
  - `quadpol_fabric_maps.py` - the multi-site quad-pol maps drawn from
    those definitions: one script, three modes (principal contrast; the
    eigenvalue difference projected onto GRID north; both horizontal
    axes as crosses whose ARM-LENGTH DIFFERENCE carries dlam, the mean
    length deliberately meaningless because common-offset polarimetry
    cannot constrain absolute eigenvalue magnitudes). Grid north because
    per-block meridian convergence, measured numerically, reaches ~69 deg
    at Ridge A - there the true-north and grid-north projections disagree
    IN SIGN. Extents are square with span-scaled markers, so fixed-km
    bars sized for Ridge A's ~25 km grid cannot sprawl across Taylor
    Dome's few-km-wide strip.
  - `quadpol_depth_movie.py` - the survey map redrawn per sliding depth
    window at a fixed colour scale and extent, so the movie shows the
    fabric strengthening down the column rather than every depth
    autoscaling to look equally anisotropic; its docstring carries the
    measured Ridge A profile and the 100%-coverage fact that makes the
    pale shallow frames measurements of weak fabric, not gaps.
  - `ridge_a_ifg_movie.py` - the wrapped-phase movie: a fixed basemap
    (REMA over grey) with a walking segment marker beside that segment's
    wrapped interferogram, no panel titles and no zoom inset (either
    would rewrite or re-frame itself every segment), colorbar below the
    map. `--still` renders one segment through the identical code path -
    extent and TWTT window still computed over ALL segments - so a slide
    cuts to the movie with zero on-screen movement. Reads the
    `scripts/reduce_ifg_movie.py` npz extracts.
  - `quadpol_heading_test.py` and `ershadi_heading_test.py` - the decisive
    ice-or-antennas test, run on an existing raster survey at no extra
    acquisition cost: theta expressed geographically must be independent of
    the driving heading. The failure mode drawn in panel (a) is the diagonal
    OFFSET by the measured theta_ant, not the 1:1 line, because both
    pipelines form theta_geo as theta_ant + heading (mod 180); it is drawn
    wrapped into two branches and labelled with the theta_ant behind it.
    The first works per frame, the second per 200-trace heading block, which
    puts the test WITHIN frames as well as between them - same ice, same
    calibration, only the heading differing. Both screen frames on HV/VH
    reciprocity first, since the test means nothing where the
    cross-polarized channels are measuring the system.
  - `crossing_pairs.py` - the same-ice crossing-pair test of the
    heading-family dlam systematic, and the acceptance criterion for the
    planned chan_equal raw-channel calibration (the paired difference
    should go to ~0 once the cross-pol pedestal is removed at source).
    Where two lines of different heading families cross, their nearest
    blocks see the same ice, so their difference isolates what the
    acquisition geometry adds and not the survey's real NW-SE gradient.
    Blocks are classified by their OWN sec_az, since curved connectors
    change family mid-frame, and pairs are deduped to one per crossing so
    a cluster meeting at a single crossing cannot vote repeatedly. That
    dedup clusters midpoints across the whole family combo rather than
    within a frame pair, because a crossing is a PLACE: scoped per frame
    pair it let one crossing vote once per frame covering it, which on
    this grid inflated the N-S/row count from 28 to 45. All 17 of those
    removals are one physical line re-flown under a second frame tag, not
    a lost independent sample - `merge_audit` establishes that on the
    CANDIDATE set rather than on the survivors (whose separations exceed
    the cluster radius by construction, so they could never have shown
    over-merging): no cluster is more than 101 m wide across track where
    two distinct parallel lines are ~1.5 km apart, and the closest
    candidate midpoint in another cluster is 1487 m against a 458 m
    radius. It re-derives both on every run, on the same clustering the
    reported pairs came from, and it runs BEFORE each combo's headline:
    a combo whose audit names a fused cluster is not reported at all (no
    median, no n, and it is left out of the npz) and the run exits
    non-zero, so a survey with tighter line spacing cannot quietly
    report too few crossings. Clusters it could not measure, because
    MIN_CELLS holes left their blocks with no index-adjacent sibling to
    take a direction from, are counted as unverified rather than passed,
    while a side that was measured over the tolerance fails its cluster
    whatever the other side did - evidence of a fusion outranks the
    absence of evidence. The npz records that verdict per combo
    (verified/failed/unverified cluster counts, the across-track width,
    the candidate gap, the radius) with a run-level `audit_pass`, so a
    later comparison can tell a fully audited run from one it should not
    lean on.
    Its reach is bounded: projecting across track separates two lines
    only insofar as they are near-parallel, so two converging
    same-family lines or one line curving back over itself are outside
    what it can certify.
    `test_crossing_pairs.py` beside it pins the invariant on synthetic
    geometry (a line split across consecutive frames, a curved connector
    contributing to both families, a repeat pass, two crossings that must
    stay apart, a deliberately fused pair the audit must catch, a cluster
    nothing could measure that must read unverified rather than clean,
    and a cluster with one unmeasurable side and one measured fusion that
    must read failed); it reads no data files and runs under pytest or
    directly.
- `scripts/prototypes/deltak_remedy.py` - a PARKED investigation into the
  delta-k stage-A suppression at Ridge A (deramped coherent multilook
  before the cross products). It runs and reports, but it does not resolve
  the suppression and its deep profile still anti-correlates with SNAPHU;
  the docstring states the numbers. Not wired into `+ptt`.
- `scripts/prototypes/extract_polpair.py` (run on mem1) and
  `scripts/prototypes/crossover_check.py` - the PoRaPy comparison. The
  first synthesizes along- and across-axis co-pol traces from one frame's
  full scattering matrix, so no crossing geometry or two-line
  coregistration enters; the second runs PoRaPy's sliding-bin correlation
  on them (`PORAPY_SRC` points at Lilien's package, which is not
  vendored). Their docstrings deliberately carry the measured results as
  the record: the raw differential phase confirms the LS chain band for
  band, and the envelope correlation recovers ~a third of the delay -
  explicitly not a criticism of PoRaPy, whose real method correlates on a
  bed reflection this radar never reaches.
- `scripts/figures/ghost2_swath_movie.py` sits outside that batch-output
  family altogether (as `ridge_a_ifg_movie.py` above does, from its own
  mem1 extracts): it renders a look-angle sweep movie (one
  along-track radargram per steering angle, TWTT axis, rotating beam icon)
  and an energy-vs-look-angle plot from the rds GHOST2 swath data that
  rode the same 2024-25 traverse train. It reads either a raw
  `CSARP_music3D` .mat (on the server, via h5py) or the `extract_swath.py`
  npz (locally), writes .mp4 when ffmpeg is available and .gif otherwise,
  and `--selftest` exercises it on a synthetic cube where the real,
  server-only data is unavailable.

### Figure inputs and outputs

Figures follow this project's usual convention of keeping mirrored
products out of git, on both ends.

INPUT data - ITS_LIVE velocity windows streamed from S3, the NEGIS track
and interferogram extracts pulled off mem1, the `fabric_sections.mat`
stage - is staged under `~/data/opr/scar`; override with the `SCAR_DATA`
environment variable. The EastGRIP figures additionally read
`egrip_fabric_949248.tab`, the core's own fabric eigenvalues (Weikusat et
al. 2022, PANGAEA.949248, CC-BY-4.0); it is third-party data, so it is
staged there like the rest rather than committed.

OUTPUT goes to the repository's `figs/` directory, which `.gitignore`
excludes, so figures live beside the code that made them without ever
being pushed. Scripts take an `<out_dir>` argument and default to `figs/`,
so `python scripts/figures/<name>.py` with no arguments does the right
thing from anywhere in the tree; pass a path only when you deliberately
want the output somewhere else.
