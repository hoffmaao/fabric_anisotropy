# Scripts

Every script that validates or applies the method code in
[`method.md`](method.md), end to end.

Figure scripts write into the repository's `figs/` directory, which is
gitignored - see the "Figure inputs and outputs" section at the end.

## Running on the CReSIS machines

Processing happens in a scratch work root, laid out as

    <work>/code          this repo, including the batch launchers
    <work>/stages        products (stages/quadpol/, coreg_cache/, ...)
    <work>/invert_logs   per-frame logs

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
its locks and write its fail markers in the LIVE tree while MATLAB, on the
same `+ptt` anchor, wrote products into the test tree. A launcher copied out
of the repo has no `+ptt` above it and is its own work root.

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
point at this section instead of restating them - each header describing
only what its own side does.

**`stages/`.** The launchers require it to pre-exist; `fabric_paths.m` does
not. The reason is silent-versus-loud failure, not who creates the
directory: of the 17 MATLAB consumers that write under `stages/`, only three
create it (`run_quadpol_pipeline`, `quadpol_coreg_frame`,
`run_deltak_stages`) and the other fourteen neither create nor check it,
calling `save()` straight in. So MATLAB against a missing `stages/` still
fails loudly - just **late**, those fourteen dying at the closing `save()`
with "Cannot create file" after a full survey or section pass, discarding
the compute. The launchers have no such backstop: their skip-if-cached
checks would match nothing and their `mkdir` locks would land in a fresh
tree, so the run would quietly recompute every frame at ~50 min each.
Pre-checking is what turns that silence into an error; `fabric_paths.m`
needs no such pre-check because it cannot be silent.

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
<work>/code/opr_fabric/server/coreg_batch.sh ridge_a` - and equally from a
copy in `<work>`. **There is no deploy step**, and **nothing carries a
username**, so the same tree runs from any user's scratch and a second
checkout can sit beside the live one for testing.

Deliberately NOT relative to the working directory: the batch launchers cd
to `<work>`, but the one-liner form needs the script's directory on the
path to find the script at all, so cwd is not a reliable anchor for both.

DATA paths stay absolute - `/cresis/...` season roots and `gps_dir` are real
mount points, not part of the work tree, and each season declares its own
`site_root` in its `run_season_*` driver. So does the shared
`/kucresis/scratch/software/snaphu`.

### Standing up a work root

The repo is public, so the CReSIS machines clone it over HTTPS with no
GitHub credentials (verified from mem1, 22 Sep 2026):

    mkdir -p <work>/stages/quadpol
    git clone https://github.com/hoffmaao/fabric_anisotropy.git <work>/code

A clone that lives somewhere else, such as beside your OPR scripts as in
[`cresis_tutorial.md`](cresis_tutorial.md), serves just as well: export
`FABRIC_ROOT=<work>` and the scripts take the work root from there instead
of from the clone's parent. For a branch that has not been pushed, ship the
history as a bundle:

    git bundle create /tmp/fabric.bundle <branch>            # locally
    cat /tmp/fabric.bundle | ssh mem1 'cat > ~/fabric.bundle' # ship
    git clone --branch <branch> ~/fabric.bundle <work>/code   # on mem1

The `mkdir` is what makes `<work>` a work root: it is the `stages/` the
launchers confirm before they will run, so create it first - without it they
refuse rather than reaching for the enclosing directory's. Nothing else is
deployed: the launchers stay in `<work>/code/opr_fabric/server` and find
their workers, `<work>` and MATLAB's path from there.

That is a real clone at a known commit, which a copied tree is not: the
live `code/` was hand-copied file by file and `git -C code rev-parse HEAD`
fails on it, so there is no way to say which commit produced a product.

### Reproducing a known result

A frame reruns from its coregistration cache in ~30 min instead of ~100,
so stage just the cache for the frame you want:

    cp <live>/stages/quadpol/coreg_cache/creg_<tag>.mat \
       <work>/stages/quadpol/coreg_cache/
    cd <work> && matlab -batch "addpath('<work>/code/opr_fabric/server'); \
      day_seg='20250108_02'; frm=9; run_season_ridge_a"

Verified 29 Aug 2026 across three seasons - `20250108_02_009` (Ridge A),
`20260106_02_001` (Taylor Dome) and `20221206_02_001` (Eastwind, which
exercises the short-frame window clamp, 2.43 m spacing and length-based
block sizing). A fresh clone into an empty work root reproduced the staged
products **bit-identically, 52 of 52 leaves** including the nested
`row_med` struct, with no configuration beyond the cache.

The one field that did not match was `day_seg` on Eastwind, and in the
direction that matters: the LIVE product holds the unreadable string
marker and the reproduction holds the char text. 21 of 141 staged products
carry that marker (all 17 Eastwind, 2 from 2024, 2 from 2026) because they
were produced before `char(day_seg)` was added. Nothing is lost - `tag`
goes through `sprintf` and is correct in every product, so `day_seg` is
`tag` minus the frame suffix - and every new run writes it readable.

Compare with h5py rather than by file size - two identical products differ
in bytes - and dereference cell fields such as `pairs` rather than
comparing HDF5 object references, which never match.

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
- `scripts/compare_const_theta.py` - run on mem1 (h5py): differences
  every `_ct` (constant-orientation) product against its default twin
  over the quotable 200-1500 m band, and states per site whether the
  assumption holds. It compares the SEGMENT profiles the blocks inherit
  (`ls_theta_seg`) and the BLOCK fits (`sec_dlam_ls`, `sec_resid_ls`),
  not the frame-pooled pass, which constant mode never touches on a
  multi-segment frame - an earlier version read `ls_theta0_geo` and
  reported two products identical while the block residual had moved by
  0.05. The verdict is the contrast-weighted per-window axis spread
  before holding (`th_spread_seg`), with the residual cost as the second
  opinion; the pooled contrast is printed but never read as evidence.
  Its thresholds are provisional until the three Thwaites control frames
  in `opr_fabric/server/const_theta_batch.sh` (the known rotating column,
  run precisely because the assumption is wrong there) calibrate them
  against the divides. That batch writes `_ct` products for every cached
  frame with a season root except the rest of Thwaites (`DRY=1` only
  builds the work list), runs from the repo directory and makes MATLAB
  prove which pipeline it loaded, because a stale work-root copy once
  silently dropped the `_ct` tag and overwrote a default product.
- `scripts/prototypes/coreg_lag_diag.py`, `coreg_alongtrack_diag.py`,
  `coreg_ramp_diag.py` - the 2 Sep 2026 measurement of WHY the Thwaites
  margin loses HH-VV coherence, run locally on a shipped
  `CSARP_polarimetric` product (`~/projects/polarimetry_thwaites/data/`):
  raw vs registered coherence with the applied offset field and the depth
  phase gradient; the HH-VV cross-correlation as a function of range LAG
  in 60 m windows (the correlation the tile pick actually faces);
  coherence against along-track window length with the lateral phase
  gradient; and block coherence with and without a local along-track
  phase-ramp compensation. Verdict on 20240108_01_001: a single lag peak
  everywhere, dead at all lags where coherence is dead, so range
  coregistration is not the loss; the western fringe fan is depolarised
  at the single-look level (true |C| below ~0.2 at 40-60 dB SNR), the east
  is the 45-deg blind geometry the quad-pol synthesis handles, and below
  ~1400 m it is SNR. `ptt.coregistration` (a package copy of the toolbox
  function, reachable through `ptt.coregisterChannels` with `impl='ptt'`)
  keeps the one real defect found - a quadratic peak refinement that
  leaves its patch and applies -21.8 bin shifts - behind `peak_clamp_en`.
- `scripts/figures/power_pattern_figs.py` - the power-extinction figure
  set: Ridge A 009's pattern P(psi, z) with axis and nodes, its node series
  and dlam against the coherence estimator; three Thwaites margin blocks'
  patterns; and the margin context panel - raw |C|, power-only node
  visibility kappa, power |dlam| where kappa allows over the coherence
  chain's blocks, and the birefringent delay predicted from the fitted
  fabric drawn over the toolbox's measured tile offsets. Inputs are the
  scratch products of `power_figdata.m` (power fits on the moment dumps)
  and `coreg_lag_diag.py`.
- `opr_fabric/server/dump_block_moments.m` - per-frame along-track block
  moment matrices from the coreg cache, antenna AND geographic frame,
  with block position and heading (positions from the polarimetric
  product, or from the section product's blocks for the qlook seasons).
  A frame becomes tens of MB that every moment-level estimator runs on
  in minutes on a laptop, on identical blocks;
  `scripts/prototypes/block_estimator_compare.m` runs the coherence LS
  (free and held axis) and the power extinction fit on them, equalising
  the channels first (`ptt.equaliseChannels`) and rotating the power
  fit's pedestal basis by the block heading.
- `opr_fabric/server/egrip_records_resync.m` - re-syncs the
  2024_Greenland_Ground2 records to the final GPS with the toolbox's
  `records_update`, in a PRIVATE support tree (~/scratch/opr_support_egrip,
  gps symlinked to the group directory), driven by each records file's
  own stored parameters. Every day's gps file is cresis-final_2025 now,
  but the records of 20240620_01 onward were still arena-field with
  positions in the Gulf of Guinea (median shift on re-sync 8,725 km).
  CAUTION recorded in the header: records_update also regenerates the
  reference trajectory under the DATA output path, and its first run
  (before that path was overridden) rewrote the production
  CSARP_reference_trajectory/ref_20240620_02.mat; the override now sends
  everything to <private>/opr_data. The channel products for those
  segments are still the stale ones until they are reprocessed from the
  re-synced records.
- `opr_fabric/server/gnss_despike.m` - repairs isolated bad fixes in a
  ground-traverse GNSS solution, writing to a private support tree with
  `gps_source` suffixed `+despike`. A sample is a spike when it departs
  from a running-median local trend by more than `dev_max` metres, NOT
  when the speed implied from the previous fix is high: that speed test is
  sample-rate dependent and fails silently, flagging 14229 non-spikes in
  EastGRIP's 20 Hz solution while missing its maximum entirely, and
  working at Kamb's 1 Hz. Runs longer than `max_run` are reported and left,
  since a sustained excursion is more likely real than a glitch. Used for
  Kamb (2022_Antarctica_Ground), which has no post-processed solution to
  re-sync to; EastGRIP needed `egrip_records_resync.m` instead and is
  clean by this test.
- `scripts/make_bed_by_block.py` - builds the `bed_by_block.json` that
  `quadpol_depth_movie.py` masks its blocks with, from the `CSARP_layer`
  picks, given one or more season roots holding `CSARP_layer/<day_seg>/`.
  Bed depth is derived as `(twtt_bottom - twtt_surface) * c /
  (2*sqrt(3.171))` - the SAME permittivity the pipeline's depth axis uses,
  so the mask and the fabric it masks sit on one ice model, and with no
  firn correction for the same reason. The surface and bottom layers are
  chosen BY NAME, never by position: these files carry internal reflectors
  too, and differencing whatever happens to be layer 2 would give a
  systematically shallow bed. The bed preference is (1) the per-trace
  median of whichever of `bottom_HH`/`bottom_VV`/`bottom_HV`/`bottom_VH`
  the file carries, (2) `bottom`, (3) `bottom_mc`, (4) a uniquely
  bottom-named layer. There is no positional tier at all: `lyr_id` is the
  row's ordinal, so `id 2` is `bottom_mc` at McMurdo and `surface_dem` at
  Eastwind, and binding it differenced a 224 m Eastwind shelf against a
  0.6 m DEM and wrote a 1 m bed. THE NAMES COME FROM A DIFFERENT FILE THAN
  THE PICKS: the per-frame `Data_<seg>_<frm>.mat` carries `twtt` but no
  `lyr_name`, and the catalogue sits once per segment in `layer_<seg>.mat`,
  so row i of `twtt` is joined to entry i of `lyr_name` after checking the
  two counts agree. With that catalogue in hand the ROWS of `twtt` are its
  layers, never whichever axis is longer: a frame with fewer traces than
  layers - the 2022 sites are 13 soundings at one spot - would otherwise
  read as one bogus layer per trace. A frame whose layers cannot be named,
  or whose row count does not match the catalogue, is SKIPPED with a line
  in the log, never bound by position and never transposed until the counts
  agree. `bottom_mc` ranks LAST on measurement, not taste: over the 9
  frames of 2022/2023_Antarctica_Ground carrying both, it sits a median
  44.6 m ABOVE the polarimetric bed while `bottom_mc_bot` sits 44.4 m
  below, so the pair brackets a basal zone rather than picking it - and
  preferring it would grey out ~45 m of real ice at Eastwind and McMurdo,
  the two seasons where it is the only non-polarimetric option. The SURFACE
  order is the reverse - exact `surface` first, per-channel second -
  because all 107 layer files across the three seasons carry a plain
  `surface` and 69 carry a `surface_dem` beside it that is a DEM, not a
  radar pick, and must never win; `_dem` layers are barred from the fuzzy
  tier for the same reason. Anything ambiguous or unnamed is skipped with a
  line in the log. The layers each frame's bed was differenced from are
  recorded under the output's reserved `_layers` key, so a wrong binding is
  auditable after the run. Each section block takes the nearest pick WITHIN
  ITS OWN FRAME and within half a block (62.5 m), because the trace axes
  differ once the qlook frames are culled and an index-for-index match
  would be wrong; measured over 56 Taylor Dome blocks the nearest pick sits
  a median 0.7 m away, so the radius turns nothing away. The picks
  themselves pass through UNTOUCHED - no thickness plausibility filter, no
  substituted frame medians, because thin ice is real at these sites and a
  small pick is not evidence of a bad one. The one validity check is on the
  SIGN: a bottom at or above the surface is a malformed record rather than
  thin ice, and becomes null with a warning naming the frame and the count,
  because a negative depth is finite and would grey that block for the
  whole movie where null leaves it drawn (measured, this fires on nothing -
  0 of 66010 finite picks). Blocks with no pick within the radius are
  written as null and drawn unmasked. The JSON is written to a temp name
  and renamed into place, as the pipeline does for its coreg cache, so a
  killed run cannot leave a truncated file that the movie then treats as
  fatal. The output is a JSON object keyed by frame tag (`20250108_02_009`,
  matching `quadpol_section_<tag>.mat`), each value a list of metres below
  the surface in block order whose length is that frame's block count; the
  movie warns and skips masking for any frame where the length disagrees,
  which is the signal to rerun this after a block-size change. Input
  sections and the output file both live under `SCAR_DATA`. The write is a
  full replacement and REFUSES to drop a season: sections are globbed from
  all of `SCAR_DATA` while layer files are only sought under the roots
  given, so a one-site rerun would leave every other season's frames out -
  and an absent tag is drawn unmasked, so those movies would stop masking
  on nothing more than the per-frame warning the movie prints. If the
  existing file carries frames this run produced none for, the script
  exits non-zero listing them with the reason it logged for
  each - no layer file under the roots given, an unreadable one, a
  catalogue that does not align - or noting that the frame was never
  reached. It does not assert a cause: a narrowed root list is only one
  explanation, and the stricter name/row alignment can legitimately skip a
  frame an older, looser run bound. Rerun with every site root, or pass
  `--replace` to write this run's frames alone. Stale entries are never
  merged forward, which would hide which run produced what.
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
  conventions (input staged under `SCAR_DATA`, output to `figs/` unless an
  `<out_dir>` argument says otherwise) and import `scar_style.py` for them:
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
    ceiling, movie depth range and window/step, and the comments record
    the measurement behind each. The COLOUR CEILING is measured, not
    chosen: the p95 of the quantity the movie actually draws (per-block
    median dlam over one depth window) over well-fit windows only, which
    reproduces Ridge A's long-standing 0.08 and un-clips the sites where
    a site-median ceiling was saturating half the map. The MOVIE RANGE
    stops where the site stops being measurable - the estimator's 0.01
    abstention threshold at most sites, but the BED where a `CSARP_layer`
    pick puts it shallower than a fabric reading would suggest (Taylor
    Dome), since a cut read off decaying dlam can end the movie inside
    real ice. WINDOW/STEP are per site: thin-ice sites get 40 m windows
    at 5 m steps - 150 m is half the ice at a 300 m shelf, and McMurdo's
    earlier 450 m cut was detecting the ice base, not a fabric limit.
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
    pale shallow frames measurements of weak fabric, not gaps. Segments
    are blended toward grey by their LS residual and drop out entirely
    once the bed enters the window, so the key carries a grey swatch:
    grey is OFF the ramp, not the bottom of it. The bed comes from
    `bed_by_block.json` (below); with no such file the movie says so and
    draws every block unmasked.
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
    channel calibration (`ptt.calibrateChannels`, issue #23: the paired
    difference should go to ~0 once the channel gains are removed at
    source, and its co-pol phase is the one term the pairs settle).
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
want the output somewhere else, or set `FABRIC_FIGS` to move that default
for a whole session. The directory is resolved from `scar_style.py`'s own
location and created on import, so nothing has to exist beforehand.
