# fabric_anisotropy

Inferring ice crystal orientation fabric from polarimetric radar traveltimes,
implementing the theory of Rathmann (2026), "Inferring glacier ice-crystal
orientation fabrics from oblique polarimetric radar" (Proc. R. Soc. A,
RSPA-2026-0001; reference Python code: https://github.com/nicholasmr/ptt).
Target application: polarimetric data from the CReSIS ground-based
accumulation radar.

## Method code (`+ptt` package)

Coordinate convention: z is height above the bed, zhat = z/H in [0, 1].
Units: meters and nanoseconds.

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
- Interferometric processing chain (pure numerics, no OPR dependencies;
  options structs use the same field names as the OPR fabric worksheet):
  - `ptt.blendTraveltime` - dtau map from interferogram phase blended with
    coregistration offsets (surface referencing, sign detection, fringe
    ambiguity resolution)
  - `ptt.blockAverage` - coherence-weighted along-track block averaging
    with per-block fringe correction
  - `ptt.invertBlocks` - twtt-to-depth mapping and per-block
    layer-stripping inversion
  - `ptt.twttDepthMap` - vertical twtt vs depth from the column model

Scripts (each validates or applies the above end to end):

- `scripts/synthetic_experiments.m` - reproduces the paper's synthetic CMP
  experiments 1-3 (noisy forward data, uninformed initial guess).
- `scripts/synthetic_common_offset.m` - synthetic validation of the
  common-offset layer-stripping method.
- `scripts/accum_inversion_template.m` - template for real accumulation
  radar picks (fill in section 1).

## OPR toolbox integration (`opr_fabric/`)

`opr_fabric/` contains a drop-in OPR processing module (`fabric.m`,
`fabric_task.m`, `run_fabric.m`) that consumes the toolbox's existing
`CSARP_polarimetric` product (interferogram + coherence + SNAPHU phase +
coregistration offsets) and produces `CSARP_fabric` profiles of the
horizontal fabric contrast. See `opr_fabric/README.md` for the pipeline
and deployment notes, and `opr_fabric/test/` for its end-to-end test.

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

- Image is `mathworks/matlab:r2024b` (linux/amd64, emulated via Rosetta on
  Apple Silicon). To match the MATLAB release on the CReSIS servers, change
  the tag in `docker-compose.yml` and `docker compose up -d` again.
- The base image ships MATLAB only. To add toolboxes, run e.g.:

  ```sh
  docker exec -it --user root fabric-matlab \
    mpm install --release=r2024b --destination=/opt/matlab/R2024b \
    Signal_Processing_Toolbox Image_Processing_Toolbox
  ```

  (Toolbox installs go into the container layer; they persist across
  `stop`/`start` but are lost if the container is removed. If a fixed toolbox
  set emerges, bake it into a Dockerfile.)
