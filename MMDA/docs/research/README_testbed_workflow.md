# LIS Patch Testbed Workflow

This testbed is for fast, real LIS data assimilation screening on small CONUS
patches before launching full-domain experiments.

## Directory Policy

Use these locations consistently:

```bash
TESTBED_ROOT=/land1/user/gychoi/LIS/testbed
OUTPUT_ROOT=/land4/LSM_OUTPUT/Noah-MP5.0/ERA5/025deg/DA_testbed
TOOLS=$TESTBED_ROOT/tools/patch_testbed
PY=/home/gychoi/anaconda3/envs/kite/bin/python
```

Patch inputs, configs, run scripts, and executables are kept under:

```bash
$TESTBED_ROOT/runs/<run_name>
```

LIS model and DA outputs are written under:

```bash
$OUTPUT_ROOT/<run_name>
```

Do not use `/tmp/lis_patch_testbed_runs` for persistent runs.

## 1. CDF Preparation

For method screening, do not rebuild tiny-patch CDFs. Use full-domain CDFs from
the parent experiment and subset them spatially with the patch tool. This keeps
the CDF transform consistent with the full CONUS experiment.

The EnKF 30-member base config currently uses:

```bash
$TESTBED_ROOT/templates/lis.config.enkf_30m_base_2019_2020
```

Check CDF paths with:

```bash
rg -n "CDF file|CDF read option|number of bins" \
  $TESTBED_ROOT/templates/lis.config.enkf_30m_base_2019_2020
```

If the CDF algorithm itself changes, rebuild the CDF on full CONUS or a large
stable domain first, then use this patch workflow to subset it. Do not derive
scientific conclusions from CDFs generated only on a tiny patch.

## 2. Generate Patch Inputs

To build all three standard EnKF 30-member patches:

```bash
cd $TESTBED_ROOT
$TOOLS/make_enkf_regional_2019_2020_30m_inputs.sh all
```

To build one region only:

```bash
$TOOLS/make_enkf_regional_2019_2020_30m_inputs.sh west
$TOOLS/make_enkf_regional_2019_2020_30m_inputs.sh central
$TOOLS/make_enkf_regional_2019_2020_30m_inputs.sh east
```

The standard regions are:

```text
west     : lon=-110.800, lat=39.533, half_width=7
central  : lon= -98.803, lat=30.421, half_width=7
east     : lon= -87.450, lat=34.183, half_width=7
```

Each run directory contains:

```text
$TESTBED_ROOT/runs/enkf_<region>_2019_2020_30m/
  LIS
  run_lis_direct_patch.sh
  lis.config.enkf_<region>_2019_2020_30m.patch
  lis.config.enkf_<region>_2019_2020_30m.patch.smap_patch_obsdir
  patch_inputs_v3/
```

## 3. What The Patch Input Script Does

`make_enkf_regional_2019_2020_30m_inputs.sh` runs two steps.

First, it spatially subsets LIS inputs:

```bash
$PY $TOOLS/make_spatial_patch_inputs.py \
  $TESTBED_ROOT/templates/lis.config.enkf_30m_base_2019_2020 \
  --center-lon <lon> \
  --center-lat <lat> \
  --half-width-grid 7 \
  --out-dir $TESTBED_ROOT/runs/<run_name>/patch_inputs_v3 \
  --write-config $TESTBED_ROOT/runs/<run_name>/lis.config.<run_name>.patch \
  --align-mask-to-restart \
  --drop-perturb-restart
```

This subsets:

- LIS domain and DA observation domain
- CDF files and CDF domain files
- Noah-MP restart file

The restart subset is coordinate-based. It does not blindly slice `ntiles`.
For a 30-member restart, the expected size is:

```text
ntiles = patch_active_cells x 30
```

Check it with:

```bash
ncdump -h $TESTBED_ROOT/runs/<run_name>/patch_inputs_v3/restart/LIS_EnRST_NOAHMP50_201812312330.d01_v2_30_conus.nc | sed -n '1,20p'
```

Second, it creates a patch-specific SMAP L2 symlink cache:

```bash
$PY $TOOLS/make_smap_patch_obsdir.py \
  --config $TESTBED_ROOT/runs/<run_name>/lis.config.<run_name>.patch \
  --src-dir /land2/SMAP/level2 \
  --out-dir $TESTBED_ROOT/runs/<run_name>/patch_inputs_v3/smap_l2_patch_bbox_2019_2020 \
  --write-config $TESTBED_ROOT/runs/<run_name>/lis.config.<run_name>.patch.smap_patch_obsdir \
  --start 2018-12-31T23:30 \
  --end 2020-12-31T23:59
```

The final config to run is always the `.smap_patch_obsdir` file.

## 4. Run LIS

Run from the individual run directory:

```bash
cd $TESTBED_ROOT/runs/enkf_west_utah_plateau_2019_2020_30m
./run_lis_direct_patch.sh > run_$(date +%y%m%d_%H%M).log 2>&1 &
```

For central and east:

```bash
cd $TESTBED_ROOT/runs/enkf_central_texas_2019_2020_30m
./run_lis_direct_patch.sh > run_$(date +%y%m%d_%H%M).log 2>&1 &

cd $TESTBED_ROOT/runs/enkf_east_alabama_tennessee_2019_2020_30m
./run_lis_direct_patch.sh > run_$(date +%y%m%d_%H%M).log 2>&1 &
```

The run script uses 12 MPI ranks:

```bash
mpirun -np 12 ./LIS -f ./lis.config.<run_name>.patch.smap_patch_obsdir
```

## 5. Check Run Status

Check logs:

```bash
tail -f $TESTBED_ROOT/runs/<run_name>/log/lislog.0000
```

Check outputs:

```bash
find $OUTPUT_ROOT/<run_name> -maxdepth 3 -type f | sed -n '1,40p'
```

Expected output subdirectories:

```text
DAOBS
DAPERT
EnKF
SURFACEMODEL
```

## 6. Validate Against ISMN

Use the patch validation script after outputs exist:

```bash
cd /land1/user/gychoi/LIS/LISF_MMDA

$PY tools/patch_testbed/validate_patch_vs_full_ismn.py \
  --patch-run $OUTPUT_ROOT/<run_name> \
  --ol-run /land4/LSM_OUTPUT/MMDA/MP50_OL_era \
  --enkf-run /land4/LSM_OUTPUT/MMDA/MP50_EnKF30_era_std_NoStatePert_2018_2022_smaptime \
  --out-dir $OUTPUT_ROOT/<run_name>/ismn_validation_vs_ol_enkf
```

This produces station metrics, time series figures, precipitation bars, and SMAP
assimilation-event markers from DAOBS.

## 7. Common Checks

Confirm the final config does not point to `/tmp`:

```bash
rg -n "/tmp/lis_patch_testbed_runs" $TESTBED_ROOT/runs/<run_name> -g '*.config*'
```

Confirm the output directory:

```bash
rg -n "Output directory" $TESTBED_ROOT/runs/<run_name>/lis.config.<run_name>.patch.smap_patch_obsdir
```

Confirm the patch SMAP cache:

```bash
find $TESTBED_ROOT/runs/<run_name>/patch_inputs_v3/smap_l2_patch_bbox_2019_2020 -maxdepth 1 -type l | wc -l
```

Confirm restart member count:

```bash
ncdump -h $TESTBED_ROOT/runs/<run_name>/patch_inputs_v3/restart/LIS_EnRST_NOAHMP50_201812312330.d01_v2_30_conus.nc | sed -n '1,12p'
```
