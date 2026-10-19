# LIS Patch Testbed

This directory contains a reusable testbed for fast, real-LIS data assimilation
screening before launching full CONUS experiments.

The design goal is conservative:

- never modify existing production run directories;
- never modify the baseline LISF source tree for experimental code changes;
- cap every generated MPI run at 12 ranks;
- preserve the same LIS executable, DA code path, CDF files, observation readers,
  perturbation files, and restart logic as the parent experiment;
- make short real runs easy first, and only enable spatial patch subsetting after
  restart/CDF tile-order checks pass.

## Recommended Workflow

1. If the method requires source changes, clone the source tree first:

```bash
python tools/patch_testbed/clone_lis_source.py \
  --source /land1/user/gychoi/LIS/LISF_MMDA \
  --name LISF_MMDA_test_DOEE \
  --out-root /tmp/lis_patch_testbed_sources
```

Only edit and compile the cloned source tree. Do not patch the baseline source
tree when testing an experimental method.

2. Build a short isolated experiment from an existing run:

```bash
python tools/patch_testbed/build_patch_experiment.py \
  --template-run /land1/user/gychoi/LIS/MMDA/lis/run_pooled_MMDA \
  --template-config lis.config.super_pooled_anomcdf_3hrcdf_2019_2021_ycdf \
  --exe LIS_MMDA_0816 \
  --experiment smoke_mmda_pooled_ycdf \
  --start 2021-01-01T00:00:00 \
  --end 2021-01-03T23:59:00 \
  --out-root /tmp/lis_patch_testbed_runs
```

3. Preflight the generated run:

```bash
python tools/patch_testbed/run_preflight.py \
  /tmp/lis_patch_testbed_runs/smoke_mmda_pooled_ycdf
```

4. Run only after preflight passes:

```bash
cd /tmp/lis_patch_testbed_runs/smoke_mmda_pooled_ycdf
./run_lis.sh
```

For a one-rank smoke test, the builder also writes `run_lis_direct.sh`. This
uses the same executable and config but bypasses `mpirun`, which is useful on
interactive nodes where MPI launch sockets are restricted.

Some LIS executables still initialize MPI/ESMF even in one-rank direct mode.
If a run fails immediately with `UCX ERROR socket create failed`, the failure is
from the execution environment, not from the generated config. Run the generated
`run_lis.sh` on a normal compute node, or run `run_lis_direct.sh` directly from
a shell where MPI initialization is allowed.

5. For spatial patch runs that use SMAP L2, create a patch-specific SMAP
   observation directory before running a multi-day or longer window:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/make_smap_patch_obsdir.py \
  --config /tmp/lis_patch_testbed_runs/<run>/lis.config.<run> \
  --out-dir /tmp/lis_patch_testbed_runs/<run>/patch_inputs/smap_l2_patch_qc \
  --write-config /tmp/lis_patch_testbed_runs/<run>/lis.config.<run>.smap_patch_obsdir \
  --start 2021-01-01T00:00 \
  --end 2021-02-01T23:59 \
  --require-da-qc
```

Then run the `.smap_patch_obsdir` config. This does not change the DA physics:
it only points LIS at a directory containing symlinks to SMAP L2 files that
overlap the patch after the SMAP DA QC screen. Files outside the patch are not
opened by LIS.

This step is important for actual period tests. The unmodified SMAP L2 reader
searches by UTC hour and then reads full swath arrays before spatial filtering.
On a tiny patch, most full-swath reads can produce zero usable local
observations. The patch observation directory removes those expensive no-op
reads without modifying the baseline source tree.

6. Diagnose the output:

```bash
python tools/patch_testbed/diagnose_patch_run.py \
  /tmp/lis_patch_testbed_runs/smoke_mmda_pooled_ycdf
```

## Screening Suites

Methods are registered in `methods.example.json`. Time windows are registered
in `windows.example.json`. Build a method x window suite with:

```bash
python tools/patch_testbed/build_window_suite.py \
  --registry tools/patch_testbed/methods.example.json \
  --windows tools/patch_testbed/windows.example.json \
  --tag trial01 \
  --out-root /tmp/lis_patch_testbed_runs \
  --max-ranks 12
```

Use `--method` and `--window` to build a smaller subset. The generated run
directories are isolated from the source run directories and write outputs under
their own `output/` directories.

## ISMN-Centered Windows

First select station candidates:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/select_ismn_patches.py \
  --n-each 2 \
  --min-distance-deg 2.0
```

Then find actual assimilation windows where those station-grid cells have finite
DA diagnostics in an existing full run:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/find_station_assim_windows.py \
  /land4/LSM_OUTPUT/MMDA/<full_run> \
  --template-config /tmp/lis_patch_testbed_runs/<generated_run>/lis.config.<generated_run> \
  --patches tools/patch_testbed/patches.generated.json \
  --out tools/patch_testbed/windows.station_assim.json \
  --not-before 2021-01-01T00:00:00
```

Build runs from those station windows:

```bash
python tools/patch_testbed/build_window_suite.py \
  --registry tools/patch_testbed/methods.example.json \
  --windows tools/patch_testbed/windows.station_assim.json \
  --tag trial01 \
  --patches tools/patch_testbed/patches.generated.json \
  --max-ranks 12
```

After a run finishes, sample the DA diagnostics back at the selected station
grid cells:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/station_da_report.py \
  /tmp/lis_patch_testbed_runs/<run>
```

`run_preflight.py` checks that restart timestamps are not after the generated
run start time. This prevents invalid tests that use a future restart to replay
an earlier assimilation window.

`diagnose_patch_run.py` scans both LIS logs and root-level run logs such as
`patch_testbed_run.log`. If you intentionally re-run a case with a different
stdout file name, old failed `patch_testbed_run.log` files may still be reported
until they are overwritten or removed.

## Verification Tools

Basic config/run checks:

```bash
python tools/patch_testbed/run_preflight.py /tmp/lis_patch_testbed_runs/<run>
```

CDF inventory:

```bash
python tools/patch_testbed/cdf_inventory.py /tmp/lis_patch_testbed_runs/<run>/lis.config.<run>
```

Patch-subsetting safety check, using the same Python environment used for
notebook diagnostics:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/verify_patch_inputs.py \
  /tmp/lis_patch_testbed_runs/<run>/lis.config.<run>
```

`verify_patch_inputs.py` returning `spatial_patch_subset_safe = NO` means the
run is still usable as a short-duration real-LIS screening run, but its
parameter/CDF/restart files should not be spatially cropped yet.

Restart-to-domain grid map:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/make_restart_grid_map.py \
  /tmp/lis_patch_testbed_runs/<run>/lis.config.<run> \
  --out-nc /tmp/lis_patch_testbed_runs/<run>_restart_grid_map.nc \
  --out-csv /tmp/lis_patch_testbed_runs/<run>_restart_grid_map.csv
```

This writes a 2D map of the restart compressed grid index and tile start
position. Use it before writing any spatial subsetter; it tells whether restart
grid cells fall on `LANDMASK` or only on the broader `DOMAINMASK`.

Spatial patch inputs:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/make_spatial_patch_inputs.py \
  /tmp/lis_patch_testbed_runs/<run>/lis.config.<run> \
  --patches tools/patch_testbed/patches.generated.json \
  --patch-name improved_filtered_01_SCAN_Uvalde \
  --out-dir /tmp/lis_patch_testbed_runs/<run>/patch_inputs \
  --write-config /tmp/lis_patch_testbed_runs/<run>/lis.config.<run> \
  --align-mask-to-restart
```

GMAO perturbation restart subset:

```bash
/home/gychoi/anaconda3/envs/kite/bin/python \
  tools/patch_testbed/subset_dapert_restart.py \
  /tmp/lis_patch_testbed_runs/<run>/lis.config.<run> \
  --full-config /path/to/full/lis.config \
  --secondary \
  --update-config
```

The primary `LIS_DAPERT_*.bin` is required when `Perturbations start mode` is
`restart`; LIS treats a missing primary DAPERT as fatal. Secondary MM-EnKF
DAPERT files are optional for correctness because the secondary reader can
cold-start, but subsetting them avoids unnecessary random-field generation when
the file dimensions match the patch grid.

## CDF Policy

For method screening, use `subset_existing_cdf` semantics: keep the same CDF
files used by the full-domain experiment. This keeps the CDF transform identical
to the full CONUS run.

For experiments that change the CDF algorithm itself, first rebuild the CDF on
the full CONUS or a large stable domain, then subset only after validating the
`ngrid` and restart tile ordering. Do not rebuild tiny-patch CDFs for the main
comparison because their sample population can be unstable.

## Patch Subsetting Status

Current verification result for
`/tmp/lis_patch_testbed_runs/spatial_mmda_uvalde_20210105_v6`:

- patch grid: 15 x 15 around SCAN Uvalde;
- common active cells after reconciling LIS domain, DA obs domain, model
  restarts, and all CDF domains: 183;
- model restarts: `ntiles = 1830 = 183 x 10 members`;
- CDF files: `ngrid = 183`;
- `verify_patch_inputs.py` reports `spatial_patch_subset_safe = YES`;
- 1-rank real LIS run completed with no obvious log errors;
- output includes DAOBS, MMENKF innov/mean/spread/stddev diagnostics,
  SURFACEMODEL history/restart, and a new primary DAPERT restart.
- SMAP patch-observation directory test:
  - original SMAP candidate files for the 2021-01-05 12:00 to 14:59 window: 4;
  - patch+DA-QC overlapping files: 1;
  - rerun with the patch SMAP directory completed in 34 s;
  - core innov diagnostics matched the original run exactly
    (`maxabsdiff = 0` for `obs_std_01`, `innov_combined_01`, `ninnov_01`,
    `K_super_01`, `Pf_super_inter_01`, and `Pf_super_intra_01`).

Important constraints:

- CDF domain files can store the active mask in `SoilMoist_domain`, not only
  `LANDMASK`; spatial input generation handles this explicitly.
- Full-run restart active cells and CDF active cells are not necessarily
  identical. The patch active mask must be their coordinate intersection.
- Fortran sequential unformatted DAPERT files are big-endian in this build.
- Do not create a spatial crop by slicing compressed `ngrid` or `ntiles`
  indices directly. Use coordinate joins.
