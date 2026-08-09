# fabric_anisotropy

Inferring ice crystal orientation fabric from polarimetric radar traveltimes,
implementing the theory of Rathmann (2026), "Inferring glacier ice-crystal
orientation fabrics from oblique polarimetric radar" (Proc. R. Soc. A,
RSPA-2026-0001; reference Python code: https://github.com/nicholasmr/ptt).
Target application: polarimetric data from the CReSIS ground-based
accumulation radar.

## Method code (`+ptt` package)

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
  are below the noise level, at the cost of some depth resolution.
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
    rather than the printed eps'.
  - `ptt.coregisterChannels` - aligns every channel onto the reference by
    CALLING the OPR toolbox `coregistration`, deliberately with no
    implementation of its own (it errors if the toolbox is absent), so the
    quad-pol path is aligned by the same code and defaults that built the
    shipped `CSARP_polarimetric` products. A hand-rolled global-delay fit
    was tried first and was worse on both pairs; the header has the numbers.

Scripts (each validates or applies the above end to end):

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
  conventions (see "SCAR figure inputs and outputs" below):
  - `scar_two_panel.py` - survey overview + wrapped interferogram for
    Thwaites and Ridge A: the survey over log-scale ITS_LIVE speed with
    all lines white and the focused profile black, red start / white end
    markers repeated on the map inset and at the interferogram's left and
    right edges, no legends, plain black scale bars. Thwaites is framed on
    the WHOLE staged ITS_LIVE window rather than on a box padded round its
    tracks, so the glacier trunk the survey sits beside stays in frame.
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
    end markers, scale bar, panel, zoom-inset and continental-locator
    geometry and the track-azimuth/geometry helpers those figures share,
    so the site slides cannot drift apart. Every two-panel figure draws
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
    circular spread there against 41.7 geographic.
  - `quadpol_heading_test.py` and `ershadi_heading_test.py` - the decisive
    ice-or-antennas test, run on an existing raster survey at no extra
    acquisition cost: theta expressed geographically must be independent of
    the driving heading, and the diagonal of panel (a) is the failure mode.
    The first works per frame, the second per 200-trace heading block, which
    puts the test WITHIN frames as well as between them - same ice, same
    calibration, only the heading differing. Both screen frames on HV/VH
    reciprocity first, since the test means nothing where the
    cross-polarized channels are measuring the system.
- `scripts/prototypes/deltak_remedy.py` - a PARKED investigation into the
  delta-k stage-A suppression at Ridge A (deramped coherent multilook
  before the cross products). It runs and reports, but it does not resolve
  the suppression and its deep profile still anti-correlates with SNAPHU;
  the docstring states the numbers. Not wired into `+ptt`.
- `scripts/figures/ghost2_swath_movie.py` is the one figure script outside
  that batch-output family: it renders a look-angle sweep movie (one
  along-track radargram per steering angle, TWTT axis, rotating beam icon)
  and an energy-vs-look-angle plot from the rds GHOST2 swath data that
  rode the same 2024-25 traverse train. It reads either a raw
  `CSARP_music3D` .mat (on the server, via h5py) or the `extract_swath.py`
  npz (locally), writes .mp4 when ffmpeg is available and .gif otherwise,
  and `--selftest` exercises it on a synthetic cube where the real,
  server-only data is unavailable.

### SCAR figure inputs and outputs

The SCAR talk figures follow this project's usual convention of keeping
mirrored products out of git. Their INPUT data - ITS_LIVE velocity
windows streamed from S3, the NEGIS track and interferogram extracts
pulled off mem1, the `fabric_sections.mat` stage - is staged under
`~/data/opr/scar`, override with the `SCAR_DATA` environment variable.
The EastGRIP figures additionally read `egrip_fabric_949248.tab`, the
core's own fabric eigenvalues (Weikusat et al. 2022, PANGAEA.949248,
CC-BY-4.0); it is third-party data, so it is staged there like the rest
rather than committed.
Their OUTPUT goes wherever the `<out_dir>` argument points, which for the
talk is `~/presentations/SCAR_figures`, outside the repo; unlike the rest
of `scripts/figures/`, they do not write to `figs/`.

## OPR toolbox integration (`opr_fabric/`)

`opr_fabric/` contains a drop-in OPR processing module (`fabric.m`,
`fabric_task.m`, `run_fabric.m`) that consumes the toolbox's existing
`CSARP_polarimetric` product (interferogram + coherence + SNAPHU phase +
coregistration offsets; in delta-k mode the complex ref/sec SLCs
instead) and produces `CSARP_fabric_joint` profiles of the
horizontal fabric contrast. See `opr_fabric/README.md` for the pipeline,
output naming, and deployment notes, and `opr_fabric/test/` for its
end-to-end test.

Figures land in `figs/`. To run without a MATLAB license (e.g. CI or quick
checks), the scripts also work in Octave:

```sh
docker run --rm --platform linux/amd64 \
  -v "$PWD":/work -w /work/scripts gnuoctave/octave:latest \
  octave --no-gui synthetic_experiments.m
```

## MATLAB environment

Local MATLAB for testing before running on the CReSIS servers
(mem1/mem2/gpu1). Uses the same browser-based MATLAB (`matlab-proxy`) as
the OPR hackathon containers.

## Usage

```sh
docker compose up -d
open http://localhost:8888
```

Sign in with your MathWorks account (online licensing) when prompted. If you
have access to a network license server instead, set `MLM_LICENSE_FILE` in
`docker-compose.yml`.

This directory is mounted inside the container at
`/home/matlab/fabric_anisotropy`, so edits on the host are visible in MATLAB
immediately and anything MATLAB writes there survives container restarts.

Stop with `docker compose down`.

## Notes

- The image is `fabric-matlab-custom:r2024b`: `mathworks/matlab:r2024b`
  (linux/amd64, emulated via Rosetta on Apple Silicon) with the Optimization
  Toolbox added via `mpm` and committed locally with `docker commit` (see the
  comment in `docker-compose.yml`). On a machine without that committed
  image, or to match the MATLAB release on the CReSIS servers, start from the
  matching `mathworks/matlab` tag, install the toolboxes as below, and
  `docker commit` it under the name in `docker-compose.yml`.
- To add more toolboxes, run e.g.:

  ```sh
  docker exec -it --user root fabric-matlab \
    mpm install --release=r2024b --destination=/opt/matlab/R2024b \
    Signal_Processing_Toolbox Image_Processing_Toolbox
  docker commit fabric-matlab fabric-matlab-custom:r2024b
  ```

  (The `docker commit` persists the toolboxes into the image; without it they
  are lost when the container is removed. If a fixed toolbox set emerges,
  bake it into a Dockerfile.)
