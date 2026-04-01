# Brainstorm: Prepare Repository for GitHub Hosting

**Date:** 2026-04-01
**Status:** Draft

## What We're Building

Preparing the `decapping_sensor_bs_scripts` repository for GitHub hosting so it can be installed and used consistently across lab computers (Mac, Linux, Windows). This involves renaming, cross-platform compatibility, documentation, and standardized assay configuration presets.

## Key Decisions

### 1. Repository Name: `mask-intensity-analysis-repo`

The repo will be renamed from `decapping_sensor_bs_scripts` to **`mask-intensity-analysis-repo`**. The name describes what the code does (intensity analysis using image masks) rather than a specific assay. The README will describe it as "a repository containing all the masked intensity analysis approaches used in the lab."

### 2. Script Renaming

| Current Name | New Name |
|---|---|
| `analyze_m7G_cap.py` | `per_particle_analysis.py` |
| `analyze_whole_field.py` | `whole_field_analysis.py` |

Names describe the analysis approach (per-particle donut vs. whole-field aggregate) rather than the specific assay.

### 3. Installation Method: Clone + venv + requirements.txt

Users will:
1. Clone the repo from GitHub
2. Create a Python virtual environment
3. Install dependencies via `pip install -r requirements.txt`
4. Run scripts directly with `python per_particle_analysis.py ...`

Platform-specific instructions (Mac/Linux/Windows) will be provided in the README for venv creation and activation (different shell commands).

**Dependencies for requirements.txt:**
- numpy
- pandas
- scipy
- scikit-image
- tifffile

### 4. Assay Configuration Presets via `--preset` Flag

A `--preset <name>` flag will be added to each script to load standardized analysis parameters for a specific assay. This ensures methodological consistency across the lab.

**Behavior rules:**
- If `--preset` is passed alongside any individual analysis flag that the preset controls, the script errors out (no silent overrides allowed).
- Presets only control analysis parameters. Users always provide `--data-dir` and `--output` paths themselves.
- `--single-cell` is never part of a preset; users add it if needed.

#### Preset: `m7g-cap-v1` (for `per_particle_analysis.py`)

| Parameter | Value |
|---|---|
| `--buffer` | 5 |
| `--donut` | 5 |
| `--bg-mode` | donut-mean |
| `--min-size` | 10 |
| `--bgsub-k` | 2.5 (default, bgsub enabled) |

**Usage:**
```
python per_particle_analysis.py --preset m7g-cap-v1 --data-dir /path/to/data --output-pbody /path/to/pbody.csv --output-sg /path/to/sg.csv
```

#### Preset: `decapping-sensor-v1` (for `whole_field_analysis.py`)

| Parameter | Value |
|---|---|
| `--min-size` | 5 |
| `--mng-bg-mode` | 0 |
| `--halo-bg-mode` | SiR_mean |
| `--mNG-filter` | NaN |
| `--percent` | True |

**Usage:**
```
python whole_field_analysis.py --preset decapping-sensor-v1 --data-dir /path/to/data --output /path/to/results.csv
```

### 5. README Structure

The README will include:
1. **Repo description** — what `mask-intensity-analysis-repo` is and what it does
2. **Installation instructions** — platform-specific (Mac, Linux, Windows) for clone + venv + requirements.txt
3. **CLI usage for `per_particle_analysis.py`** — complete flag/option reference table
4. **CLI usage for `whole_field_analysis.py`** — complete flag/option reference table
5. **File naming conventions** — how input TIFFs must be named (channel keywords, group key prefix)
6. **Assay-specific sections:**
   - **Decapping Sensor Analysis** — which script, which preset, expected inputs, how to interpret outputs
   - **m7G Cap Immunofluorescence Analysis** — which script, which preset, expected inputs, how to interpret outputs
7. **Preset reference** — all available presets and their locked parameters

### 6. Windows Compatibility

Key changes needed:
- Replace any hardcoded `/` path separators with `os.path.join()` or `pathlib.Path`
- Ensure TIFF file reading works with Windows-style paths
- Test that `argparse` boolean flags work correctly on Windows Python
- Provide Windows-specific venv instructions (`python -m venv venv` + `venv\Scripts\activate`)

## Why This Approach

- **Clone + venv** is the simplest installation method that all lab members can follow regardless of Python packaging experience.
- **`--preset` with strict conflict errors** ensures methodological consistency — no one accidentally runs a slightly different configuration.
- **Generic script names** make the tools reusable for future assays without renaming.
- **Detailed README** reduces the need for in-person onboarding.

## Open Questions

None remaining.

## Resolved Questions

1. **Archive directory** — Exclude via .gitignore. Old scripts are superseded and would confuse new users.
2. **Image data & results** — Exclude all data (TIFFs), results (CSVs, XLSX), and Prism files via .gitignore. Only code, docs, and requirements.txt go on GitHub.
3. **Python version** — Require Python 3.9+ for broadest compatibility across lab machines.

## .gitignore Plan

Exclude:
- `archive/`
- `venv/`
- `*.tif`, `*.tiff`
- `*.csv`, `*.xlsx`
- `*.prism`
- `*_images/`, `*_results/`, `*_cap_images/`, `*_cap_results/`
- `__pycache__/`, `*.pyc`
- `.DS_Store`
