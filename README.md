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
  - `ptt.blockAverage` - coherence-weighted along-track block averaging
    with per-block fringe correction
  - `ptt.invertBlocks` - twtt-to-depth mapping and per-block inversion
    (exact layer stripping or the regularized joint solve, selected by
    `opts.inversion`)
  - `ptt.twttDepthMap` - vertical twtt vs depth from the column model

Scripts (each validates or applies the above end to end):

- `scripts/synthetic_experiments.m` - reproduces the paper's synthetic CMP
  experiments 1-3 (noisy forward data, uninformed initial guess).
- `scripts/synthetic_common_offset.m` - synthetic validation of the
  common-offset layer-stripping method.
- `scripts/accum_inversion_template.m` - template for real accumulation
  radar picks (fill in section 1).
- `scripts/figures/` - Python (matplotlib + scipy; cartopy required only
  by the regional maps) figure scripts that reproduce the analysis figures
  from `opr_fabric/server/run_fabric_scratch.m` batch outputs mirrored
  locally: season transects, per-season depth profiles, the 2024-25
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
  an all-survey summary (coverage map, median profiles, drive orientations;
  the map uses cartopy if installed, else a plain lat/lon scatter),
  along-profile depth sections of the fabric solution per survey region,
  and zoomed regional maps of each survey area with reference-site markers
  and an Antarctica locator. The maps overlay LIMA/MODIS MOA satellite imagery
  when rasterio and the locally downloaded mosaics are present (download
  commands in `scripts/figures/antarctic_basemap.py`), and fall back to
  coastline-only otherwise. See each script's docstring for usage; the
  inversion chain itself stays MATLAB/Octave.

## OPR toolbox integration (`opr_fabric/`)

`opr_fabric/` contains a drop-in OPR processing module (`fabric.m`,
`fabric_task.m`, `run_fabric.m`) that consumes the toolbox's existing
`CSARP_polarimetric` product (interferogram + coherence + SNAPHU phase +
coregistration offsets) and produces `CSARP_fabric_joint` profiles of the
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
