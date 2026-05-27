# mask-intensity-analysis-repo

A repository containing all the masked intensity analysis approaches used in the lab. Quantifies fluorescence intensity in masked subcellular regions with local background subtraction.

**Last updated:** 2026-05-27

---

## Table of Contents

- [Installation](#installation)
- [Scripts Overview](#scripts-overview)
- [File Naming Conventions](#file-naming-conventions)
  - [Channel Keywords](#channel-keywords)
  - [Example Data Directory](#example-data-directory)
  - [The cp_mask File](#the-cp_mask-file)
  - [SiR_mean Mode File Naming](#sir_mean-mode-file-naming)
- [`per_particle_analysis.py` CLI Reference](#per_particle_analysispy-cli-reference)
  - [Required Arguments](#required-arguments)
  - [Preset](#preset)
  - [Analysis Parameters](#analysis-parameters)
  - [Other Options](#other-options)
  - [Usage Examples](#usage-examples)
  - [Full `--help` output](#per_particle_analysispy-full---help-output)
- [`whole_field_analysis.py` CLI Reference](#whole_field_analysispy-cli-reference)
  - [Required Arguments](#required-arguments-1)
  - [Preset](#preset-1)
  - [Background Modes](#background-modes)
  - [Filter Options](#filter-options)
  - [Three-Region Options (v4 / v5)](#three-region-options-v4--v5)
  - [Other Options](#other-options-1)
  - [Usage Examples](#usage-examples-1)
  - [Full `--help` output](#whole_field_analysispy-full---help-output)
- [Assay: m7G Cap Immunofluorescence Analysis](#assay-m7g-cap-immunofluorescence-analysis)
- [Assay: Decapping Sensor Analysis](#assay-decapping-sensor-analysis)
  - [Two-region presets: v1, v2, v3](#two-region-presets-v1-v2-v3)
  - [Three-region presets: v4, v5](#three-region-presets-v4-v5)
  - [Required input files per preset](#required-input-files-per-preset)
  - [Output columns](#output-columns)
  - [Percent column semantics](#percent-column-semantics)
- [Generating cp_mask Files with ImageJ](#generating-cp_mask-files-with-imagej)
- [Single-Cell Output](#single-cell-output)
- [Preset Reference](#preset-reference)
- [Adding New Presets](#adding-new-presets)
- [License](#license)

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
| `whole_field_analysis.py` | Whole-field aggregate | Measures aggregate mNG and Halo channel intensities across all P-body pixels vs. all dilute-mask pixels (two-region) or across P-body, intermediate, and dilute regions (three-region: v4 / v5) in each field of view. |
| `ROI_to_CellposeLabels.ijm` | ImageJ macro | Generates a `cp_mask` cell segmentation label image from ROI Manager selections, used with the `--single-cell` flag. |

---

## File Naming Conventions

Both scripts group TIFF files by detecting **channel keywords** in the filename. Everything before the channel keyword (minus trailing underscores) becomes the **group key** that links files from the same image set.

**Pattern:** `{group_key}_{channel_keyword}.tif`

### Channel Keywords

**`per_particle_analysis.py`:** `P-body_mask`, `SG_mask`, `Cap`, `pnorm`, `sgnorm`, `cp_mask`

**`whole_field_analysis.py`:** `P-body_mask`, `dilute_mask`, `interaction_mask`, `interaction_mask_2`, `Dcp2_mask`, `Dcp2_mask_2`, `SiR_mask`, `Halo`, `mNG`, `cp_mask`

The two newer `*_2` channels (`Dcp2_mask_2`, `interaction_mask_2`) are required only when running the three-region intermediate-assemblies presets (v4 / v5) or passing `--intermediate-assemblies` directly.

### Example Data Directory

```
my_data/
├── WT_Rep1_NaAsO2_Cap.tif
├── WT_Rep1_NaAsO2_P-body_mask.tif
├── WT_Rep1_NaAsO2_SG_mask.tif
├── WT_Rep1_NaAsO2_pnorm.tif
├── WT_Rep1_NaAsO2_sgnorm.tif
├── WT_Rep1_NaAsO2_cp_mask.tif      ← optional, for --single-cell mode
├── WT_Rep2_NaAsO2_Cap.tif
├── WT_Rep2_NaAsO2_P-body_mask.tif
├── WT_Rep2_NaAsO2_SG_mask.tif
├── WT_Rep2_NaAsO2_pnorm.tif
├── WT_Rep2_NaAsO2_sgnorm.tif
└── WT_Rep2_NaAsO2_cp_mask.tif      ← optional, for --single-cell mode
```

In this example, the group keys are `WT_Rep1_NaAsO2` and `WT_Rep2_NaAsO2`. Each group has 5 files (Cap + P-body_mask + SG_mask + pnorm + sgnorm), enabling both P-body and SG analysis. The `cp_mask` files are optional and only needed when using `--single-cell`.

### The `cp_mask` File

The `cp_mask` (cellpose mask) file is a 16-bit TIFF label image used for single-cell analysis. Each pixel value represents a unique cell ID (1, 2, 3, ...), and background pixels are 0. When `--single-cell` is passed, the script uses this mask to assign each particle (or each pixel region) to the cell that contains it, then aggregates measurements per cell rather than per particle or per field.

- **Format:** 16-bit TIFF with integer labels (supports up to 65,535 cells)
- **Naming:** Must follow the standard convention: `{group_key}_cp_mask.tif`
- **Generation:** Use the included `ROI_to_CellposeLabels.ijm` ImageJ macro (see [Generating cp_mask Files with ImageJ](#generating-cp_mask-files-with-imagej) below), or any cell segmentation tool (e.g., Cellpose) that produces integer-labeled masks
- **Behavior:** If `--single-cell` is passed but a group is missing its `cp_mask` file, that group is skipped with a warning

### SiR_mean Mode File Naming

When using `--halo-bg-mode SiR_mean` (the default for `decapping-sensor-v1`), group keys must start with `As_` (arsenite-treated) or `UT_` (untreated). The script pairs groups by condition and computes the Halo background from the SiR_mask region of the As condition. The computed background is applied to both members of the pair.

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

Other presets (v2 onwards) do not require the `As_` / `UT_` naming and process each group independently. v2 requires a `SiR_mask` for every image set; v3, v4, and v5 do not use `SiR_mask` at all.

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

### `per_particle_analysis.py` Full `--help` output

```
usage: per_particle_analysis.py [-h] --data-dir DATA_DIR
                                [--preset {m7g-cap-v1}] [--buffer BUFFER]
                                [--donut DONUT]
                                [--bg-mode {donut,donut-mean,flat}]
                                [--bg-value BG_VALUE]
                                [--exclude-cap-zero | --no-exclude-cap-zero]
                                [--min-size MIN_SIZE] [--bgsub-k BGSUB_K]
                                [--no-bgsub] [--export-donuts]
                                [--output OUTPUT | --output-pbody OUTPUT_PBODY]
                                [--output-sg OUTPUT_SG] [--single-cell]

Per-particle donut background subtraction analysis of fluorescence intensities
in masked subcellular regions (P-bodies and/or Stress Granules).

options:
  -h, --help            show this help message and exit
  --data-dir DATA_DIR   Directory containing .tif/.tiff images
  --preset {m7g-cap-v1}
                        Named parameter preset for a specific assay. Locks all
                        analysis parameters; only --data-dir,
                        --output/--output-pbody/--output-sg, and --single-cell
                        can be specified alongside a preset. Available:
                        m7g-cap-v1
  --buffer BUFFER       Buffer zone dilation in pixels (default without
                        preset: 4)
  --donut DONUT         Donut ring width in pixels (default without preset: 5)
  --bg-mode {donut,donut-mean,flat}
                        Background subtraction mode: "donut" estimates per-
                        particle from donut median, "donut-mean" from donut
                        mean, "flat" uses a single fixed value (default
                        without preset: donut)
  --bg-value BG_VALUE   Flat background value when --bg-mode=flat (default
                        without preset: 1)
  --exclude-cap-zero, --no-exclude-cap-zero
                        Exclude zero values from Cap channel before estimating
                        background (default without preset: True, use --no-
                        exclude-cap-zero to disable)
  --min-size MIN_SIZE   Only analyze particles larger than this many pixels
                        (default without preset: 10)
  --bgsub-k BGSUB_K     Background subtraction threshold expressed as mu +
                        k*sigma of the background Gaussian fit. Higher k is
                        more permissive (default without preset: 2.5)
  --no-bgsub            Disable automatic background subtraction of the Cap
                        image
  --export-donuts       Export binary donut mask TIFFs (one per image set) for
                        overlay visualization
  --output OUTPUT       Base output CSV path (use with --preset). Generates
                        <name>_p-body.csv and <name>_sg.csv automatically.
  --output-pbody OUTPUT_PBODY
                        Output CSV for P-body results (use without --preset,
                        requires --output-sg)
  --output-sg OUTPUT_SG
                        Output CSV for Stress Granule results (use without
                        --preset, requires --output-pbody)
  --single-cell         Aggregate all particle measurements by single cell
                        using a _cp_mask.tiff segmentation file. Each integer
                        value in the mask represents one cell. Particles are
                        assigned to the cell containing the majority of their
                        pixels. Output will have one row per cell per image
                        set.

The script auto-detects which analyses to run based on available files per
image set. When both P-body and SG files are present, Cap pixels inside
SG_mask are set to NaN before P-body analysis (P-bodies embedded in stress
granules cannot be reliably quantified).
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
| `--preset` | Named parameter preset. Locks all analysis parameters. Only `--data-dir`, `--output`, and `--single-cell` can be specified alongside a preset. Available: `decapping-sensor-v1`, `decapping-sensor-v2`, `decapping-sensor-v3`, `decapping-sensor-v4`, `decapping-sensor-v5` |

### Background Modes

These are locked when a `--preset` is active.

| Flag | Type | Default (no preset) | Description |
|------|------|---------------------|-------------|
| `--mng-bg-mode` | str or int | `mean` | mNG background subtraction mode. Options: `mean`, `median`, `mode`, `top_quintile`, `top_quartile`, `top_decile`, or an integer for a manual flat value |
| `--halo-bg-mode` | str or int | `median` | Halo background subtraction mode. All options from `--mng-bg-mode` plus: `mng-nan` (mean of Halo where mNG is NaN in dilute mask), `mng-nan-median`, `mng-nan-mode`, `mng-nan-top_quintile`, `mng-nan-top_quartile`, `mng-nan-top_decile`, `mng-nan-max`, `SiR_mean` (mean of Halo inside SiR_mask from As condition, applied to both As and UT pairs) |
| `--exclude-halo-zero` / `--no-exclude-halo-zero` | bool | True | Exclude zero values from Halo channel when estimating background |
| `--exclude-halo-one` / `--no-exclude-halo-one` | bool | False | Exclude 0 and 1 values from Halo channel when estimating background (supersedes `--exclude-halo-zero`) |

### Filter Options

These are locked when a `--preset` is active.

| Flag | Type | Default (no preset) | Description |
|------|------|---------------------|-------------|
| `--min-size` | int | 10 | Minimum P-body mask component size in pixels |
| `--SiR-subtract` | choice | None | Filter Halo channel using SiR mask before all other filters. Options: `zero`, `NaN`. Cannot be combined with `--halo-bg-mode SiR_mean`. |
| `--SiR-filter` | flag | False | Restrict Halo measurements to within `SiR_mask`: condensed Halo = `P-body_mask & SiR_mask`, dilute Halo = `dilute_mask & SiR_mask`. mNG measurements are unaffected. Requires a `SiR_mask` file in every image set. Cannot be combined with `--SiR-subtract` or `--halo-bg-mode SiR_mean`. |
| `--FLIM-filter` | choice | None | Filter Halo channel using the FLIM `interaction_mask`. Options: `zero` (sets Halo outside `interaction_mask` to 0), `NaN` (sets to NaN). |
| `--mNG-filter` | choice | None | Filter mNG channel using `Dcp2_mask`. Options: `zero`, `NaN`. |
| `--mNG-in-FLIM` | choice | `no` | Restrict both mNG and Halo to where both `interaction_mask` and `Dcp2_mask` are positive. Options: `yes`, `no`. |
| `--percent` | flag | False | Add per-region columns reporting the percent of mNG-mask area in each compartment that also has Halo signal. See [Percent column semantics](#percent-column-semantics). |
| `--save-processed` | flag | False | Save processed mNG and Halo TIFFs after background subtraction and masking, for troubleshooting. |

### Three-Region Options (v4 / v5)

These enable the three-region "intermediate assemblies" analysis used by `decapping-sensor-v4` and `decapping-sensor-v5`. See [Three-region presets: v4, v5](#three-region-presets-v4-v5) for the full description.

| Flag | Type | Default (no preset) | Description |
|------|------|---------------------|-------------|
| `--intermediate-assemblies` | flag | False | Switch from the two-region (P-body / dilute) layout to a three-region (P-body / intermediate / dilute) layout. The intermediate region is defined by `Dcp2_mask_2` (mNG) and `interaction_mask_2` (Halo); the dilute region is redefined to exclude the intermediate. Requires `Dcp2_mask_2` and `interaction_mask_2` files in every image set, plus a global `mNG_filter` and `FLIM_filter`. Compatible with `--single-cell`. Not compatible with any SiR option. |
| `--intermediate-zero-fill` | flag | False | Modifier on `--intermediate-assemblies`: Halo means in the intermediate and dilute regions use v3-style zero-fill (pixels outside the inner interaction mask are zeroed and INCLUDED in the mean) instead of being excluded. mNG dilute region is also broadened to `dilute_mask & ~Dcp2_mask_2`. Percent columns are unchanged — only Halo means differ. Requires `--intermediate-assemblies` to be on. |

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

**With a two-region preset (v1, v2, or v3):**

```bash
python whole_field_analysis.py \
    --preset decapping-sensor-v3 \
    --data-dir /path/to/data \
    --output /path/to/results.csv
```

**With a three-region preset (v4 or v5):**

```bash
python whole_field_analysis.py \
    --preset decapping-sensor-v4 \
    --data-dir /path/to/data \
    --output /path/to/results.csv
```

**With preset and single-cell mode:**

```bash
python whole_field_analysis.py \
    --preset decapping-sensor-v4 \
    --data-dir /path/to/data \
    --output /path/to/results.csv \
    --single-cell
```

### `whole_field_analysis.py` Full `--help` output

```
usage: whole_field_analysis.py [-h] --data-dir DATA_DIR --output OUTPUT
                               [--preset {decapping-sensor-v1,decapping-sensor-v2,decapping-sensor-v3,decapping-sensor-v4,decapping-sensor-v5}]
                               [--mng-bg-mode MNG_BG_MODE]
                               [--halo-bg-mode HALO_BG_MODE]
                               [--exclude-halo-zero | --no-exclude-halo-zero]
                               [--exclude-halo-one | --no-exclude-halo-one]
                               [--min-size MIN_SIZE]
                               [--SiR-subtract {zero,NaN}] [--SiR-filter]
                               [--FLIM-filter {zero,NaN}]
                               [--mNG-filter {zero,NaN}]
                               [--mNG-in-FLIM {yes,no}] [--percent]
                               [--save-processed] [--intermediate-assemblies]
                               [--intermediate-zero-fill] [--single-cell]

Whole-field analysis of fluorescence channel intensities using P-body and
dilute masks. Background subtraction is computed per field from the dilute
mask region.

options:
  -h, --help            show this help message and exit
  --data-dir DATA_DIR   Directory containing .tif/.tiff images
  --output OUTPUT       Output CSV path
  --preset {decapping-sensor-v1,decapping-sensor-v2,decapping-sensor-v3,decapping-sensor-v4,decapping-sensor-v5}
                        Named parameter preset for a specific assay. Locks all
                        analysis parameters; only --data-dir, --output, and
                        --single-cell can be specified alongside a preset.
                        Available: decapping-sensor-v1, decapping-sensor-v2,
                        decapping-sensor-v3, decapping-sensor-v4, decapping-
                        sensor-v5
  --mng-bg-mode MNG_BG_MODE
                        mNG background subtraction: "mean" (default without
                        preset), "median", "mode", "top_quintile",
                        "top_quartile", "top_decile", or an integer for a
                        manual flat value
  --halo-bg-mode HALO_BG_MODE
                        Halo background subtraction: "median" (default without
                        preset), "mean", "mode", "top_quintile",
                        "top_quartile", "top_decile", "mng-nan" (mean of Halo
                        where mNG is NaN in dilute mask), "mng-nan-median",
                        "mng-nan-mode", "mng-nan-top_quintile", "mng-nan-
                        top_quartile", "mng-nan-top_decile", "mng-nan-max",
                        "SiR_mean" (mean of Halo inside SiR_mask from As
                        condition, applied to both As and UT pairs; cannot be
                        combined with --SiR-subtract), or an integer for a
                        manual flat value
  --exclude-halo-zero, --no-exclude-halo-zero
                        Exclude zero values from Halo channel before
                        estimating background (default without preset: True,
                        use --no-exclude-halo-zero to disable)
  --exclude-halo-one, --no-exclude-halo-one
                        Exclude 0 and 1 values from Halo channel before
                        estimating background (default without preset: False;
                        when enabled, supersedes --exclude-halo-zero)
  --min-size MIN_SIZE   Only include P-body mask components larger than this
                        many pixels (default without preset: 10)
  --SiR-subtract {zero,NaN}
                        Apply SiR mask filtering to Halo channel before all
                        other filters. "zero" sets Halo pixels where SiR_mask
                        > 0 to 0; "NaN" sets them to NaN. Requires SiR_mask
                        file in the image set. Cannot be combined with --halo-
                        bg-mode SiR_mean.
  --SiR-filter          Restrict Halo measurements to within SiR_mask:
                        condensed Halo = P-body_mask & SiR_mask, dilute Halo =
                        dilute_mask & SiR_mask. mNG measurements are
                        unaffected. Requires a SiR_mask file in every image
                        set. Cannot be combined with --SiR-subtract or --halo-
                        bg-mode SiR_mean.
  --FLIM-filter {zero,NaN}
                        Apply FLIM mask filtering to Halo channel before
                        analysis. "zero" sets Halo pixels inside
                        interaction_mask to 0; "NaN" sets them to NaN.
                        Requires interaction_mask file in the image set.
  --mNG-filter {zero,NaN}
                        Apply Dcp2_mask filtering to mNG channel before
                        analysis. "zero" sets mNG pixels inside Dcp2_mask to
                        0; "NaN" sets them to NaN. Requires Dcp2_mask file in
                        the image set.
  --mNG-in-FLIM {yes,no}
                        How to handle mNG channel when --FLIM-filter is used.
                        "yes": only analyze mNG and Halo where both
                        interaction_mask and Dcp2_mask > 0; "no" (default
                        without preset): do not change how mNG channel is
                        handled.
  --percent             Add columns reporting the percent of mNG-valid pixels
                        occupied by Halo, computed separately for P-body and
                        dilute masks
  --save-processed      Save processed mNG and Halo .tif files (after bg
                        subtraction and masking) for troubleshooting
  --intermediate-assemblies
                        Three-region measurement (P-body, intermediate,
                        dilute) using Dcp2_mask_2 and interaction_mask_2 to
                        define the intermediate region. The dilute region is
                        redefined as dilute_mask & Dcp2_mask & ~Dcp2_mask_2
                        (mNG) and dilute_mask & interaction_mask &
                        ~interaction_mask_2 (Halo). Requires Dcp2_mask_2 and
                        interaction_mask_2 files in every image set.
                        Compatible with --single-cell (each cell gets per-
                        region rows). Not compatible with any SiR option.
  --intermediate-zero-fill
                        Modifies --intermediate-assemblies to use v3-style
                        halo handling: Halo means in the intermediate and
                        dilute regions are dragged down by zeros from pixels
                        outside the inner interaction mask, instead of those
                        pixels being excluded from the mean. mNG dilute region
                        is also broadened to dilute_mask & ~Dcp2_mask_2
                        (without intersecting Dcp2_mask). Percent columns are
                        unchanged — only Halo means differ from
                        --intermediate-assemblies. Requires --intermediate-
                        assemblies to be on.
  --single-cell         Group all measurements by single cells using a
                        _cp_mask.tiff segmentation file. Each integer value in
                        the mask represents one cell. Output will have one row
                        per cell per image set.
```

---

## Assay: m7G Cap Immunofluorescence Analysis

**Script:** `per_particle_analysis.py`
**Preset:** `m7g-cap-v1`

This assay measures the enrichment of m7G Cap immunofluorescence signal in P-bodies and/or Stress Granules relative to the local cytoplasmic background. Each particle is individually analyzed using a donut-shaped annulus for background estimation.

When both P-body and SG masks are present, Cap pixels inside `SG_mask` are set to NaN before P-body analysis because P-bodies embedded in stress granules cannot be reliably quantified.

### Required Input Files

| Analysis | Required channels |
|----------|-------------------|
| P-body analysis | `P-body_mask`, `Cap`, `pnorm` |
| SG analysis | `SG_mask`, `Cap`, `sgnorm` |
| Both (recommended) | all five above |
| Single-cell mode | additionally `cp_mask` |

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

## Assay: Decapping Sensor Analysis

**Script:** `whole_field_analysis.py`
**Presets:** `decapping-sensor-v1` … `decapping-sensor-v5`

This assay measures the enrichment of a HaloTag-labeled decapping sensor in P-bodies (and, for v4 / v5, in a separate "intermediate assemblies" compartment) relative to the dilute cytoplasm. mNG is used as a normalization channel and `Dcp2_mask` defines where the mNG signal is real.

The preset versions differ in how the Halo background is handled and in whether the analysis is two-region (P-body / dilute) or three-region (P-body / intermediate / dilute). Once a preset is published, it is **never modified** — new behavior is added as a new version.

### Two-region presets: v1, v2, v3

These produce one row per image set (or per cell with `--single-cell`) with measurements over two compartments:

- **P-body** — `P-body_mask`
- **Dilute** — `dilute_mask`

| Preset | Halo background source | Notes |
|--------|------------------------|-------|
| `decapping-sensor-v1` | `SiR_mean` (mean Halo inside `SiR_mask` of the paired As field, applied to both As and UT) | Requires `As_*` / `UT_*` group key prefixes and a `SiR_mask` in the As condition. |
| `decapping-sensor-v2` | Disabled (`halo_bg_mode=0`) | Uses `SiR_mask` as an intersection filter on Halo measurements instead of for background. Requires a `SiR_mask` file in **every** image set. Percent columns are dropped. |
| `decapping-sensor-v3` | Disabled (`halo_bg_mode=0`) | Drops SiR entirely. Uses `interaction_mask` via `FLIM_filter='zero'` so Halo outside `interaction_mask` is set to 0 before measurement. Requires an `interaction_mask` file in every image set. Percent columns kept. |

### Three-region presets: v4, v5

These produce one row per image set (or per cell with `--single-cell`) with measurements over three compartments:

- **P-body** — `P-body_mask` (same as v3)
- **Intermediate** — `Dcp2_mask_2` (mNG) intersected with `interaction_mask_2` (Halo); represents "intermediate assemblies" sized between P-bodies and the dilute phase
- **Dilute** — `dilute_mask & Dcp2_mask & ~Dcp2_mask_2` (mNG) and `dilute_mask & interaction_mask & ~interaction_mask_2` (Halo); the true dilute phase with the intermediate region subtracted out

Both presets share the global v3-style processing (`mNG_filter='NaN'`, `FLIM_filter='zero'`, no SiR, no Halo background subtraction). The only thing that differs between v4 and v5 is how the **Halo means** in the intermediate and dilute regions are computed:

| Preset | Halo handling in intermediate / dilute regions |
|--------|------------------------------------------------|
| `decapping-sensor-v4` | **Intersect**: Halo means use only pixels inside the inner interaction mask. Clean signal mean — represents typical intensity *where Halo is actually detected*. |
| `decapping-sensor-v5` | **v3-style zero-fill**: pixels outside the inner interaction mask are zeroed and included in the mean, dragging it down. Represents the average intensity *across the whole compartment*, with non-detection areas counted as zero. Also broadens the mNG dilute region (and thus `dilute_area_px`) to `dilute_mask & ~Dcp2_mask_2` to match v3 convention. |

mNG means, integrals, areas, and the percent columns are identical between v4 and v5. Only `halo_intermediate_mean`, `halo_dilute_mean`, and (in v5) `dilute_area_px` differ.

### Required input files per preset

| File | v1 | v2 | v3 | v4 | v5 |
|------|----|----|----|----|----|
| `P-body_mask` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `dilute_mask` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `Halo` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `mNG` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `Dcp2_mask` | ✅ | ✅ | ✅ | ✅ | ✅ |
| `SiR_mask` (As only, As/UT pairing required) | ✅ | — | — | — | — |
| `SiR_mask` (every image set) | — | ✅ | — | — | — |
| `interaction_mask` | — | — | ✅ | ✅ | ✅ |
| `Dcp2_mask_2` | — | — | — | ✅ | ✅ |
| `interaction_mask_2` | — | — | — | ✅ | ✅ |
| `cp_mask` (for `--single-cell`) | optional | optional | optional | optional | optional |

### Output columns

Two-region output (v1, v2, v3) — common columns:

| Column | Description |
|--------|-------------|
| `pbody_area_px` | Filtered P-body mask area (pixels) |
| `dilute_area_px` | Dilute mask area (pixels) |
| `mng_bg_value` | mNG background value subtracted |
| `halo_bg_value` | Halo background value subtracted |
| `mNG_pbody_mean` / `mNG_dilute_mean` | Mean mNG in each compartment (background-subtracted) |
| `mNG_pbody_integ` / `mNG_dilute_integ` | Integrated mNG in each compartment |
| `halo_pbody_mean` / `halo_dilute_mean` | Mean Halo in each compartment (background-subtracted, restricted to mNG-valid pixels) |
| `halo_pbody_integ` / `halo_dilute_integ` | Integrated Halo in each compartment |
| `mNG_pbody_over_dilute` | mNG enrichment ratio (P-body / dilute) |
| `halo_pbody_over_dilute` | Halo enrichment ratio (P-body / dilute) |
| `halo_over_mNG_pbody` | Halo / mNG ratio in P-bodies |
| `halo_over_mNG_dilute` | Halo / mNG ratio in dilute phase |
| `pct_halo_in_mNG_pbody` | (with `--percent`) See [Percent column semantics](#percent-column-semantics) |
| `pct_halo_in_mNG_dilute` | (with `--percent`) See [Percent column semantics](#percent-column-semantics) |

Three-region output (v4, v5) — adds an `intermediate` set:

Every column above with `pbody` / `dilute` is mirrored by an `intermediate` version: `intermediate_area_px`, `mNG_intermediate_mean`, `mNG_intermediate_integ`, `halo_intermediate_mean`, `halo_intermediate_integ`, `mNG_intermediate_over_dilute`, `halo_intermediate_over_dilute`, `halo_over_mNG_intermediate`, and (with `--percent`) `pct_halo_in_mNG_intermediate`.

Single-cell mode (`--single-cell`) prepends `cell_id`, `cell_area_px`, `particle_count`, and appends `mNG_cell_mean`.

### Percent column semantics

The `pct_halo_in_mNG_*` columns (active with `--percent`, which is on by default for `decapping-sensor-v*`) answer one question per compartment:

> **Within compartment X, what fraction of the mNG-mask area also has Halo interaction signal?**

Each percent is computed as:

```
pct_halo_in_mNG_X = (compartment X's mNG-mask area ∩ relevant interaction mask) / (compartment X's mNG-mask area) × 100
```

For v4 / v5 specifically:

```
pct_halo_in_mNG_pbody        = (P-body ∩ Dcp2_mask ∩ interaction_mask)                                         / (P-body ∩ Dcp2_mask)
pct_halo_in_mNG_intermediate = (Dcp2_mask_2 ∩ interaction_mask_2)                                              / Dcp2_mask_2
pct_halo_in_mNG_dilute       = (dilute ∩ Dcp2_mask ∩ ~Dcp2_mask_2 ∩ interaction_mask ∩ ~interaction_mask_2)    / (dilute ∩ Dcp2_mask ∩ ~Dcp2_mask_2)
```

In single-cell mode, both numerator and denominator are additionally intersected with the cell region.

Each percent is bounded in [0%, 100%] and is an **independent** indicator of how much of each compartment is occupied by Halo interaction signal. The three percents are NOT three slices of a shared whole — they will not generally sum to 100% and have no reason to.

---

## Generating cp_mask Files with ImageJ

The included `ROI_to_CellposeLabels.ijm` macro creates a cellpose-style cell segmentation mask from manually drawn ROIs in ImageJ/FIJI. This mask is used with the `--single-cell` flag in either Python script.

### Prerequisites

1. Open a TIFF image of the correct dimensions in ImageJ (e.g., one of the channel images from the image set you want to analyze)
2. Draw cell boundary ROIs using any selection tool (freehand, polygon, etc.) and add each one to the ROI Manager (`T` to add, or Analyze > Tools > ROI Manager)

### Running the Macro

1. Open the macro: Plugins > Macros > Run... > select `ROI_to_CellposeLabels.ijm`
2. The macro creates a new 16-bit label image where each ROI is filled with a unique integer (1, 2, 3, ...)
3. Save the result as TIFF: File > Save As > Tiff...
4. **Rename the file** to match the naming convention: `{group_key}_cp_mask.tif` (e.g., `WT_Rep1_NaAsO2_cp_mask.tif`)
5. Place the file in the same data directory as the other channel images for that group

---

## Single-Cell Output

When `--single-cell` is used, the output changes from one row per particle (or per field) to **one row per cell per image set**:

- **`per_particle_analysis.py`:** Particles are assigned to the cell containing the majority of their pixels. Per-cell means are area-weighted across all particles in that cell, and integrated intensities are summed. Ratios are recomputed from the aggregated means. Additional columns: `cell_id`, `cell_area_px`, `n_pbodys` (or `n_sgs`), `total_pbody_area_px` (or `total_sg_area_px`).
- **`whole_field_analysis.py`:** Each compartment mask is intersected with the cell region. Intensities are measured within each cell's portion of each compartment independently. Additional columns: `cell_id`, `cell_area_px`, `particle_count`, `mNG_cell_mean`.

Cells that contain no particles (in per-particle mode) or no mask pixels (in whole-field mode) still appear in the output with NaN for all intensity metrics, so every segmented cell is represented.

`--single-cell` is compatible with every `decapping-sensor-v*` preset, including v4 / v5.

---

## Preset Reference

### `m7g-cap-v1` (per_particle_analysis.py)

| Parameter | Value |
|-----------|-------|
| `--buffer` | 5 |
| `--donut` | 5 |
| `--bg-mode` | donut-mean |
| `--min-size` | 4 |
| `--bgsub-k` | 2.5 |
| `--no-bgsub` | False (bgsub enabled) |
| `--bg-value` | 1 |
| `--exclude-cap-zero` | True |
| `--export-donuts` | False |

### `decapping-sensor-v*` (whole_field_analysis.py)

| Parameter | v1 | v2 | v3 | v4 | v5 |
|-----------|----|----|----|----|----|
| `--min-size` | 2 | 2 | 2 | 2 | 2 |
| `--mng-bg-mode` | 0 | 0 | 0 | 0 | 0 |
| `--halo-bg-mode` | `SiR_mean` | 0 | 0 | 0 | 0 |
| `--mNG-filter` | NaN | NaN | NaN | NaN | NaN |
| `--FLIM-filter` | None | None | zero | zero | zero |
| `--SiR-subtract` | None | None | None | None | None |
| `--SiR-filter` | False | True | False | False | False |
| `--mNG-in-FLIM` | no | no | no | no | no |
| `--percent` | True | False | True | True | True |
| `--exclude-halo-zero` | True | True | True | True | True |
| `--exclude-halo-one` | False | False | False | False | False |
| `--save-processed` | False | False | False | False | False |
| `--intermediate-assemblies` | False | False | False | True | True |
| `--intermediate-zero-fill` | False | False | False | False | True |

---

## Adding New Presets

To add a new preset for a new assay or a revised protocol:

1. Open the relevant script (`per_particle_analysis.py` or `whole_field_analysis.py`)
2. Find the `PRESETS` dictionary near the top of the file
3. Add a new entry with a versioned name (e.g., `m7g-cap-v2`, `my-new-assay-v1`)
4. **Never modify an existing preset** — always create a new version to preserve reproducibility. Published results that cite a preset name must remain reproducible. Adding a *new* preset-controlled key requires adding that key (with an inert default) to every existing preset so `validate_presets` passes, but the preset's *behavior* must not change.

Example:

```python
PRESETS = {
    'decapping-sensor-v5': { ... },  # DO NOT MODIFY
    'decapping-sensor-v6': {          # New version with updated parameters
        'min_size': 2,
        'mng_bg_mode': 0,
        'halo_bg_mode': 0,
        # ... all preset-controlled keys ...
    },
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
