# mask-intensity-analysis-repo

A repository containing all the masked intensity analysis approaches used in the lab. Quantifies fluorescence intensity in masked subcellular regions with local background subtraction.

**Last updated:** 2026-04-01

---

## Installation

Requires **Python 3.9 or higher**.

### Mac / Linux

```bash
git clone https://github.com/marcusjoshm/mask-intensity-analysis-repo.git
cd mask-intensity-analysis-repo
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Windows

```cmd
git clone https://github.com/marcusjoshm/mask-intensity-analysis-repo.git
cd mask-intensity-analysis-repo
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

---

## Scripts Overview

| Script | Analysis Type | Description |
|--------|--------------|-------------|
| `per_particle_analysis.py` | Per-particle donut | Measures fluorescence intensity in individual P-bodies and/or Stress Granules using a donut-shaped annulus around each particle for local background estimation. |
| `whole_field_analysis.py` | Whole-field aggregate | Measures aggregate mNG and Halo channel intensities across all P-body pixels vs. all dilute-mask pixels in each field of view. |

---

## `per_particle_analysis.py` CLI Reference

### Required Arguments

| Flag | Description |
|------|-------------|
| `--data-dir` | Directory containing .tif/.tiff images |
| `--output-pbody` | Output CSV for P-body results *(required when not using --preset)* |
| `--output-sg` | Output CSV for Stress Granule results *(required when not using --preset)* |
| `--output` | Base output CSV path *(required when using --preset; generates `<name>_p-body.csv` and `<name>_sg.csv` automatically)* |

### Preset

| Flag | Description |
|------|-------------|
| `--preset` | Named parameter preset. Locks all analysis parameters. Only `--data-dir`, `--output` / `--output-pbody` / `--output-sg`, and `--single-cell` can be specified alongside a preset. Available: `m7g-cap-v1` |

### Analysis Parameters

These are locked when a `--preset` is active. When no preset is used, they fall back to their defaults.

| Flag | Type | Default (no preset) | Description |
|------|------|---------------------|-------------|
| `--buffer` | int | 4 | Buffer zone gap in pixels between particle edge and donut inner edge |
| `--donut` | int | 5 | Donut ring width in pixels |
| `--bg-mode` | choice | `donut` | Background subtraction mode: `donut` (per-particle median), `donut-mean` (per-particle mean), `flat` (single fixed value) |
| `--bg-value` | int | 1 | Flat background value when `--bg-mode=flat` |
| `--exclude-cap-zero` / `--no-exclude-cap-zero` | bool | True | Exclude zero-valued pixels from Cap channel when estimating donut background |
| `--min-size` | int | 10 | Minimum particle size in pixels; smaller particles are excluded |
| `--bgsub-k` | float | 2.5 | Global background subtraction threshold: mu + k*sigma of the Gaussian fit to the Cap image background peak |
| `--no-bgsub` | flag | False | Disable automatic global background subtraction of the Cap image |
| `--export-donuts` | flag | False | Export binary donut mask TIFFs for overlay visualization |

### Other Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--single-cell` | flag | False | Aggregate per-particle results by cell using a `cp_mask` segmentation file. Output has one row per cell per image set. |

### Usage Examples

**Without preset (custom parameters):**

```bash
python per_particle_analysis.py \
    --data-dir /path/to/data \
    --buffer 5 --donut 5 --bg-mode donut-mean --min-size 10 \
    --output-pbody /path/to/pbody_results.csv \
    --output-sg /path/to/sg_results.csv
```

**With preset:**

```bash
python per_particle_analysis.py \
    --preset m7g-cap-v1 \
    --data-dir /path/to/data \
    --output /path/to/results.csv
# Creates: /path/to/results_p-body.csv and /path/to/results_sg.csv
```

**With preset and single-cell mode:**

```bash
python per_particle_analysis.py \
    --preset m7g-cap-v1 \
    --data-dir /path/to/data \
    --output /path/to/results.csv \
    --single-cell
```

---

## `whole_field_analysis.py` CLI Reference

### Required Arguments

| Flag | Description |
|------|-------------|
| `--data-dir` | Directory containing .tif/.tiff images |
| `--output` | Output CSV path |

### Preset

| Flag | Description |
|------|-------------|
| `--preset` | Named parameter preset. Locks all analysis parameters. Only `--data-dir`, `--output`, and `--single-cell` can be specified alongside a preset. Available: `decapping-sensor-v1` |

### Background Modes

These are locked when a `--preset` is active.

| Flag | Type | Default (no preset) | Description |
|------|------|---------------------|-------------|
| `--mng-bg-mode` | str or int | `mean` | mNG background subtraction mode. Options: `mean`, `median`, `mode`, `top_quintile`, `top_quartile`, `top_decile`, or an integer for a manual flat value |
| `--halo-bg-mode` | str or int | `median` | Halo background subtraction mode. All options from mng-bg-mode plus: `mng-nan` (mean of Halo where mNG is NaN in dilute mask), `mng-nan-median`, `mng-nan-mode`, `mng-nan-top_quintile`, `mng-nan-top_quartile`, `mng-nan-top_decile`, `mng-nan-max`, `SiR_mean` (mean of Halo inside SiR_mask from As condition, applied to both As and UT pairs) |
| `--exclude-halo-zero` / `--no-exclude-halo-zero` | bool | True | Exclude zero values from Halo channel when estimating background |
| `--exclude-halo-one` / `--no-exclude-halo-one` | bool | False | Exclude 0 and 1 values from Halo channel when estimating background (supersedes --exclude-halo-zero) |

### Filter Options

These are locked when a `--preset` is active.

| Flag | Type | Default (no preset) | Description |
|------|------|---------------------|-------------|
| `--min-size` | int | 10 | Minimum P-body mask component size in pixels |
| `--SiR-subtract` | choice | None | Filter Halo channel using SiR mask before all other filters. Options: `zero`, `NaN`. Cannot be combined with `--halo-bg-mode SiR_mean`. |
| `--FLIM-filter` | choice | None | Filter Halo channel using FLIM interaction mask. Options: `zero`, `NaN` |
| `--mNG-filter` | choice | None | Filter mNG channel using Dcp2_mask. Options: `zero`, `NaN` |
| `--mNG-in-FLIM` | choice | `no` | Restrict both mNG and Halo to where both interaction_mask and Dcp2_mask are positive. Options: `yes`, `no` |
| `--percent` | flag | False | Add columns reporting percent overlap of Halo in mNG regions |
| `--save-processed` | flag | False | Save processed mNG and Halo TIFFs after background subtraction and masking |

### Other Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--single-cell` | flag | False | Per-cell analysis using a `cp_mask` segmentation file. Output has one row per cell per image set. |

### Usage Examples

**Without preset (custom parameters):**

```bash
python whole_field_analysis.py \
    --data-dir /path/to/data \
    --mng-bg-mode mean --halo-bg-mode median --min-size 10 \
    --output /path/to/results.csv
```

**With preset:**

```bash
python whole_field_analysis.py \
    --preset decapping-sensor-v1 \
    --data-dir /path/to/data \
    --output /path/to/results.csv
```

**With preset and single-cell mode:**

```bash
python whole_field_analysis.py \
    --preset decapping-sensor-v1 \
    --data-dir /path/to/data \
    --output /path/to/results.csv \
    --single-cell
```

---

## File Naming Conventions

Both scripts group TIFF files by detecting **channel keywords** in the filename. Everything before the channel keyword (minus trailing underscores) becomes the **group key** that links files from the same image set.

**Pattern:** `{group_key}_{channel_keyword}.tif`

### Channel Keywords

**`per_particle_analysis.py`:** `P-body_mask`, `SG_mask`, `Cap`, `pnorm`, `sgnorm`, `cp_mask`

**`whole_field_analysis.py`:** `P-body_mask`, `dilute_mask`, `interaction_mask`, `Dcp2_mask`, `SiR_mask`, `Halo`, `mNG`, `cp_mask`

### Example Data Directory

```
my_data/
├── WT_Rep1_NaAsO2_Cap.tif
├── WT_Rep1_NaAsO2_P-body_mask.tif
├── WT_Rep1_NaAsO2_SG_mask.tif
├── WT_Rep1_NaAsO2_pnorm.tif
├── WT_Rep1_NaAsO2_sgnorm.tif
├── WT_Rep2_NaAsO2_Cap.tif
├── WT_Rep2_NaAsO2_P-body_mask.tif
├── WT_Rep2_NaAsO2_SG_mask.tif
├── WT_Rep2_NaAsO2_pnorm.tif
└── WT_Rep2_NaAsO2_sgnorm.tif
```

In this example, the group keys are `WT_Rep1_NaAsO2` and `WT_Rep2_NaAsO2`. Each group has 5 files (Cap + P-body_mask + SG_mask + pnorm + sgnorm), enabling both P-body and SG analysis.

### SiR_mean Mode File Naming

When using `--halo-bg-mode SiR_mean`, group keys must start with `As_` (arsenite-treated) or `UT_` (untreated). The script pairs groups by condition and computes the Halo background from the SiR_mask region of the As condition.

```
my_data/
├── As_sample1_P-body_mask.tif
├── As_sample1_dilute_mask.tif
├── As_sample1_Halo.tif
├── As_sample1_mNG.tif
├── As_sample1_SiR_mask.tif
├── UT_sample1_P-body_mask.tif
├── UT_sample1_dilute_mask.tif
├── UT_sample1_Halo.tif
└── UT_sample1_mNG.tif
```

---

## Assay: Decapping Sensor Analysis

**Script:** `whole_field_analysis.py`
**Preset:** `decapping-sensor-v1`

This assay measures the enrichment of a HaloTag-labeled decapping sensor in P-bodies relative to the dilute cytoplasm, using mNG as a normalization channel. The `SiR_mean` background mode computes Halo background from the SiR (DNA dye) mask region in the arsenite-treated condition and applies it to both treated and untreated samples.

### Required Input Files

Each image set needs: `P-body_mask`, `dilute_mask`, `Halo`, `mNG`, `Dcp2_mask`

For the `SiR_mean` background mode, the arsenite (As) condition also needs: `SiR_mask`

Group keys must be prefixed with `As_` or `UT_` for condition pairing.

### Example Command

```bash
python whole_field_analysis.py \
    --preset decapping-sensor-v1 \
    --data-dir /path/to/decapping_data/exports/tiff \
    --output /path/to/decapping_results.csv
```

### Key Output Columns

| Column | Description |
|--------|-------------|
| `mNG_pbody_mean` | Mean mNG intensity in P-bodies (background-subtracted) |
| `mNG_dilute_mean` | Mean mNG intensity in dilute phase (background-subtracted) |
| `halo_pbody_mean` | Mean Halo intensity in P-bodies (background-subtracted, restricted to mNG-valid pixels) |
| `halo_dilute_mean` | Mean Halo intensity in dilute phase (background-subtracted, restricted to mNG-valid pixels) |
| `mNG_pbody_over_dilute` | mNG enrichment ratio (P-body / dilute) |
| `halo_pbody_over_dilute` | Halo enrichment ratio (P-body / dilute) |
| `halo_over_mNG_pbody` | Halo-to-mNG ratio in P-bodies |
| `halo_over_mNG_dilute` | Halo-to-mNG ratio in dilute phase |
| `pct_halo_in_mNG_pbody` | Percent of Dcp2-positive P-body pixels with Halo signal |
| `pct_halo_in_mNG_dilute` | Percent of Dcp2-positive dilute pixels with Halo signal |

---

## Assay: m7G Cap Immunofluorescence Analysis

**Script:** `per_particle_analysis.py`
**Preset:** `m7g-cap-v1`

This assay measures the enrichment of m7G Cap immunofluorescence signal in P-bodies and/or Stress Granules relative to the local cytoplasmic background. Each particle is individually analyzed using a donut-shaped annulus for background estimation.

When both P-body and SG masks are present, Cap pixels inside SG_mask are set to NaN before P-body analysis because P-bodies embedded in stress granules cannot be reliably quantified.

### Required Input Files

**P-body analysis:** `P-body_mask`, `Cap`, `pnorm`
**SG analysis:** `SG_mask`, `Cap`, `sgnorm`
**Both (recommended):** All five files above
**Single-cell mode:** Also requires `cp_mask`

### Example Command

```bash
python per_particle_analysis.py \
    --preset m7g-cap-v1 \
    --data-dir /path/to/m7G_cap_data \
    --output /path/to/m7g_results.csv
# Creates: /path/to/m7g_results_p-body.csv and /path/to/m7g_results_sg.csv
```

### Key Output Columns (P-body)

| Column | Description |
|--------|-------------|
| `pbody_id` | Unique particle ID within the image |
| `pbody_area_px` | Particle area in pixels |
| `donut_area_px` | Donut background region area in pixels |
| `bg_value` | Estimated background intensity for this particle |
| `cap_pbody_mean` | Mean Cap intensity in the P-body (background-subtracted) |
| `cap_dilute_mean` | Mean Cap intensity in the donut (background-subtracted) |
| `cap_pbody_over_dilute` | Cap enrichment ratio (P-body / donut) |
| `pnorm_pbody_over_dilute` | Normalization channel enrichment ratio |
| `cap_over_pnorm_pbody` | Cap-to-normalization ratio in P-body |
| `cap_over_pnorm_dilute` | Cap-to-normalization ratio in donut |

SG output columns follow the same pattern with `sg` replacing `pbody` and `sgnorm` replacing `pnorm`.

---

## Preset Reference

### `m7g-cap-v1` (per_particle_analysis.py)

| Parameter | Value |
|-----------|-------|
| `--buffer` | 5 |
| `--donut` | 5 |
| `--bg-mode` | donut-mean |
| `--min-size` | 10 |
| `--bgsub-k` | 2.5 |
| `--no-bgsub` | False (bgsub enabled) |
| `--bg-value` | 1 |
| `--exclude-cap-zero` | True |
| `--export-donuts` | False |

### `decapping-sensor-v1` (whole_field_analysis.py)

| Parameter | Value |
|-----------|-------|
| `--min-size` | 5 |
| `--mng-bg-mode` | 0 (manual flat value) |
| `--halo-bg-mode` | SiR_mean |
| `--mNG-filter` | NaN |
| `--percent` | True |
| `--exclude-halo-zero` | True |
| `--exclude-halo-one` | False |
| `--SiR-subtract` | None (disabled) |
| `--FLIM-filter` | None (disabled) |
| `--mNG-in-FLIM` | no |
| `--save-processed` | False |

---

## Adding New Presets

To add a new preset for a new assay or a revised protocol:

1. Open the relevant script (`per_particle_analysis.py` or `whole_field_analysis.py`)
2. Find the `PRESETS` dictionary near the top of the file
3. Add a new entry with a versioned name (e.g., `m7g-cap-v2`, `my-new-assay-v1`)
4. **Never modify an existing preset** — always create a new version to preserve reproducibility

Example:

```python
PRESETS = {
    'm7g-cap-v1': { ... },  # DO NOT MODIFY
    'm7g-cap-v2': {          # New version with updated parameters
        'buffer': 6,
        'donut': 6,
        'bg_mode': 'donut-mean',
        'min_size': 15,
        'bgsub_k': 3.0,
        'no_bgsub': False,
        'bg_value': 1,
        'exclude_cap_zero': True,
        'export_donuts': False,
    },
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
