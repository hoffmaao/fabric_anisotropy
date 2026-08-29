# fabric_anisotropy

To Infer ice crystal orientation fabric from polarimetric radar traveltime and refelctivity, we implement the theory of Rathmann (2026), "Inferring glacier ice-crystal orientation fabrics from oblique polarimetric radar" (Proc. R. Soc. A, RSPA-2026-0001; reference Python code: https://github.com/nicholasmr/ptt) and Ershadi (2022) and use these methods to interpret polarimetric data from the CReSIS ground-based accumulation radar.

## Layout

| path | what it holds |
|---|---|
| `+ptt/` | the method: column physics, CMP and common-offset inversions, and interferometric and quad-pol estimator chain |
| `opr_fabric/` | a drop-in OPR processing module plus the server drivers and the MATLAB test suite |
| `scripts/` | synthetic validations, batch drivers, and the figure scripts |
| `docs/` | reference detail (below) |
| `figs/` | figure output, gitignored |

Two reference pages carry the detail that used to live in this file:

- [`docs/method.md`](docs/method.md) - every `+ptt` function, and choice that informed the fabric algorthim. Each function also has a full header of its own (`help ptt.<name>`).
- [`docs/scripts.md`](docs/scripts.md) - every script and figure, with information on what it validates or produces, and where its inputs come from, plus how to stand up a scratch work root and run the analysis on the CReSIS machines.

## Data and figures

Neither mirrored products nor figures are committed.

- **Input** is staged under `~/data/opr/scar`; override with the
  `SCAR_DATA` environment variable.
- **Output** goes to `figs/`, which `.gitignore` excludes. Figure scripts
  take an `<out_dir>` argument and default to `figs/`, so running one with
  no arguments puts the figure there from anywhere in the tree.

## OPR toolbox integration (`opr_fabric/`)

`opr_fabric/` contains a drop-in OPR processing module (`fabric.m`,
`fabric_task.m`, `run_fabric.m`) that consumes the toolbox's existing
`CSARP_polarimetric` product (interferogram + coherence + SNAPHU phase +
coregistration offsets; in delta-k mode the complex ref/sec SLCs instead)
and produces `CSARP_fabric_joint` profiles of the horizontal fabric
contrast. See `opr_fabric/README.md` for the pipeline, output naming and
deployment notes, and `opr_fabric/test/` for its end-to-end test.

To run without a MATLAB license (CI or quick checks), the scripts also work
in Octave:

```sh
docker run --rm --platform linux/amd64 \
  -v "$PWD":/work -w /work/scripts gnuoctave/octave:latest \
  octave --no-gui synthetic_experiments.m
```

## MATLAB environment

Local MATLAB for testing before running on the CReSIS servers
(mem1/mem2/gpu1), using the same browser-based MATLAB (`matlab-proxy`) as
the OPR hackathon containers:

```sh
docker compose up -d
open http://localhost:8888
```

Sign in with your MathWorks account (online licensing) when prompted; for a
network license server set `MLM_LICENSE_FILE` in `docker-compose.yml`
instead. The repository is mounted at `/home/matlab/fabric_anisotropy`, so
host edits appear in MATLAB immediately and anything MATLAB writes survives
container restarts. Stop with `docker compose down`.

The image is `fabric-matlab-custom:r2024b`: `mathworks/matlab:r2024b`
(linux/amd64, emulated via Rosetta on Apple Silicon) with the Optimization
Toolbox added via `mpm` and committed locally with `docker commit` (see the
comment in `docker-compose.yml`). On a machine without that committed image,
or to match the MATLAB release on the CReSIS servers, start from the
matching `mathworks/matlab` tag, install the toolboxes as below, and
`docker commit` it under the name in `docker-compose.yml`:

```sh
docker exec -it --user root fabric-matlab \
  mpm install --release=r2024b --destination=/opt/matlab/R2024b \
  Signal_Processing_Toolbox Image_Processing_Toolbox
docker commit fabric-matlab fabric-matlab-custom:r2024b
```

The `docker commit` is what persists the toolboxes; without it they are lost
when the container is removed. If a fixed toolbox set emerges, bake it into
a Dockerfile.
