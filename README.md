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
    the product's own `param.array.img_comb` (0.9 us with a 0.1 us blend
    on the accum3 grids, 8 us with a 1 us blend on the deeper settings,
    none on single-image frames) instead of being hardcoded, and the
    masked band runs forward from the boundary because that is the way
    OPR crossfades. Applied by `ptt.surfaceReference`, so every dtau
    estimator inherits it; disable with `param.fabric.seam_mask_en =
    false` or retune the guard with `param.fabric.seam_mask_win`. Delta-k
    gets a wider band: it resolves its ladder integers from a tau_A
    smoothed over analysis cells, so every cell whose smoothing window
    touched the seam is dropped too (`ptt.deltakDefaults` is the one
    definition of that cell geometry).
  - `ptt.blockAverage` - coherence-weighted along-track block averaging
    with per-block fringe correction
  - `ptt.invertBlocks` - twtt-to-depth mapping and per-block inversion
    (exact layer stripping or the regularized joint solve, selected by
    `opts.inversion`). Where the seam mask (or an incoherent run) leaves a
    gap, the node's dtau is interpolated across it so the joint solve
    keeps a continuous chain, and the node is flagged in
    `inv.interpolated`, saved as `dlam_interpolated` next to
    `dlam_quality`/`dlam_clipped`. The figure scripts drop those intervals
    the way they drop pegged ones - the unmasked coherence in
    `dlam_quality` cannot tell them apart from measured nodes.
  - `ptt.twttDepthMap` - vertical twtt vs depth from the column model

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
- `scripts/figures/` - Python (matplotlib + scipy; cartopy required only
  by the regional maps) figure scripts that reproduce the analysis figures
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
  the matching delta-k eigenvalue-difference depth sections along both
  traverses (the Ridge A panel is flagged unvalidated there because
  delta-k is amplitude-suppressed against the joint chain at that site),
  and zoomed regional maps of each survey area with reference-site markers
  and an Antarctica locator. The maps overlay LIMA/MODIS MOA satellite imagery
  when rasterio and the locally downloaded mosaics are present (download
  commands in `scripts/figures/antarctic_basemap.py`), and fall back to
  coastline-only otherwise. See each script's docstring for usage; the
  inversion chain itself stays MATLAB/Octave.
- `scripts/figures/ghost2_swath_movie.py` is the one figure script outside
  that batch-output family: it renders a look-angle sweep movie (one
  along-track radargram per steering angle, TWTT axis, rotating beam icon)
  and an energy-vs-look-angle plot from the rds GHOST2 swath data that
  rode the same 2024-25 traverse train. It reads either a raw
  `CSARP_music3D` .mat (on the server, via h5py) or the `extract_swath.py`
  npz (locally), writes .mp4 when ffmpeg is available and .gif otherwise,
  and `--selftest` exercises it on a synthetic cube where the real,
  server-only data is unavailable.

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
