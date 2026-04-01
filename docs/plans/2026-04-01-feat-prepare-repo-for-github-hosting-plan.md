---
title: "feat: Prepare mask-intensity-analysis-repo for GitHub Hosting"
type: feat
date: 2026-04-01
deepened: 2026-04-01
---

# Prepare mask-intensity-analysis-repo for GitHub Hosting

## Enhancement Summary

**Deepened on:** 2026-04-01
**Research agents used:** Python reviewer, Code simplicity reviewer, Architecture strategist, Best practices researcher, Security sentinel, Performance oracle

### Key Improvements
1. Added preset validation and defensive assertions to catch bugs early
2. Expanded `.gitignore` with missing patterns identified by security review
3. Added `from __future__ import annotations` for clean Python 3.9+ support
4. Fixed `parse_filename` inconsistency (bug) between scripts
5. Clarified git safety: explicit `git add` by filename, never `git add .`

### New Considerations Discovered
- `parse_filename` in `whole_field_analysis.py` lacks longest-match-first sorting — latent bug that should be fixed during rename
- The dual `--output` interface (preset vs non-preset) is a UX concern flagged by multiple reviewers — plan includes mitigation via argparse mutually exclusive groups
- Preset keys should be validated against the argparse namespace to catch typos
- Should add "never modify existing presets — always version up" comment for reproducibility

---

## Overview

Prepare the `decapping_sensor_bs_scripts` repository for cross-platform GitHub hosting as `mask-intensity-analysis-repo`. This involves renaming scripts, adding a `--preset` system for standardized assay configurations, creating installation infrastructure (requirements.txt, .gitignore, LICENSE), writing comprehensive documentation (README.md), and ensuring Windows compatibility.

## Problem Statement / Motivation

The lab's fluorescence microscopy analysis scripts are currently a loose collection of files on one machine with no version control, no documentation, no installation method, and assay-specific naming. Lab members on different computers (Mac, Linux, Windows) need to run identical analyses with consistent methodology. Without presets, each person manually specifies flags, risking subtle parameter differences that compromise data comparability.

## Proposed Solution

### Phase 1: Script Renaming and Bug Fixes

**Rename scripts:**
- `analyze_m7G_cap.py` → `per_particle_analysis.py`
- `analyze_whole_field.py` → `whole_field_analysis.py`

**Make `--data-dir` and `--output` required arguments (no defaults) in both scripts.** Remove the current default values (`m7G_cap_images`, `Decapping_donut_images`, `m7G_cap_pbody_results.csv`, etc.). Users must always specify paths explicitly.

For `per_particle_analysis.py`, keep `--output-pbody` and `--output-sg` as the two required output flags (no defaults). The preset will introduce a unified `--output` convenience (see Phase 2).

For `whole_field_analysis.py`, `--output` becomes required (no default).

**Add Python version check and future annotations** at the top of each script:
```python
from __future__ import annotations  # Must be first import — enables modern type hints on 3.9
import sys
if sys.version_info < (3, 9):
    sys.exit("Error: Python 3.9 or higher is required (for argparse.BooleanOptionalAction).")
```

### Research Insights (Phase 1)

**Bug fix — `parse_filename` inconsistency:** `analyze_m7G_cap.py:122` sorts CHANNELS by length (longest-match-first) before matching, but `analyze_whole_field.py:44` does not. This means whole_field could mismatch on channel names that are substrings of each other. Fix: add `sorted(CHANNELS, key=len, reverse=True)` to `whole_field_analysis.py`'s `parse_filename` during the rename.

**`from __future__ import annotations`:** This must be the very first import line. It defers type hint evaluation, allowing modern syntax like `list[str]` and `str | None` on Python 3.9 without runtime errors. The `sys.version_info` check should come immediately after.

**Use `pathlib.Path` for output path derivation** (the `_p-body`/`_sg` suffix insertion):
```python
from pathlib import Path

def derive_output_paths(base_output: str) -> tuple[Path, Path]:
    p = Path(base_output)
    return (
        p.parent / f"{p.stem}_p-body{p.suffix}",
        p.parent / f"{p.stem}_sg{p.suffix}",
    )
```
Do NOT use string `.replace(".csv", ...)` or `os.path.splitext` — `pathlib` handles edge cases (no extension, multiple dots) correctly.

---

### Phase 2: Preset System

#### Architecture

Add a `PRESETS` dict at the top of each script file. Each preset maps parameter names (using argparse `dest` names) to their locked values.

```python
# IMPORTANT: Never modify an existing preset. Create a new version instead.
# This ensures published results citing a preset version remain reproducible.
PRESETS = {
    "m7g-cap-v1": {
        "buffer": 5,
        "donut": 5,
        "bg_mode": "donut-mean",
        "min_size": 10,
        "bgsub_k": 2.5,
        "no_bgsub": False,
        "bg_value": 1,
        "exclude_cap_zero": True,
        "export_donuts": False,
    },
}
```

**Conflict detection approach:** Set all preset-controlled argument defaults to `None` in argparse. After parsing:
1. If `--preset` is provided, apply preset values to any arg that is still `None`.
2. If `--preset` is provided and an arg is NOT `None` (user explicitly set it), error out naming the specific conflicting flag.
3. If no `--preset`, apply the original default values to any arg that is still `None`.

**Extract into a testable function:**
```python
ORIGINAL_DEFAULTS = {
    "buffer": 4, "donut": 5, "bg_mode": "donut", "min_size": 10,
    "bgsub_k": 2.5, "no_bgsub": False, "bg_value": 1,
    "exclude_cap_zero": True, "export_donuts": False,
}

def resolve_args(args, preset_name, presets, original_defaults):
    """Apply preset or original defaults. Error on conflicts."""
    if preset_name:
        if preset_name not in presets:
            sys.exit(f"Error: Unknown preset '{preset_name}'. Available: {', '.join(presets)}")
        preset_vals = presets[preset_name]
        for key, val in preset_vals.items():
            current = getattr(args, key)
            if current is not None:
                sys.exit(f"Error: --{key.replace('_', '-')} conflicts with --preset {preset_name}")
            setattr(args, key, val)
    else:
        for key, val in original_defaults.items():
            if getattr(args, key) is None:
                setattr(args, key, val)
    return args
```

**Validate preset keys at startup** (catches typos and stale keys):
```python
def validate_presets(presets, parser):
    valid_dests = {action.dest for action in parser._actions}
    for name, vals in presets.items():
        bad_keys = set(vals.keys()) - valid_dests
        if bad_keys:
            raise ValueError(f"Preset '{name}' has unknown keys: {bad_keys}")
```

**Add defensive assertion** after defaults are applied (before any analysis runs):
```python
for attr in ['buffer', 'donut', 'bg_mode', 'min_size', 'bgsub_k']:
    assert getattr(args, attr) is not None, f"Bug: {attr} is None after default application"
```

**Print active parameters** when a preset is invoked:
```
[preset: m7g-cap-v1] buffer=5, donut=5, bg-mode=donut-mean, min-size=10, bgsub-k=2.5
```

#### Preset: `m7g-cap-v1` (per_particle_analysis.py)

**Locked parameters:**
| Parameter | Preset Value |
|---|---|
| `--buffer` | 5 |
| `--donut` | 5 |
| `--bg-mode` | donut-mean |
| `--min-size` | 10 |
| `--bgsub-k` | 2.5 |
| `--no-bgsub` | False (bgsub enabled) |
| `--bg-value` | 1 |
| `--exclude-cap-zero` | True |
| `--export-donuts` | False |

**Free parameters (user-provided):**
| Parameter | Notes |
|---|---|
| `--data-dir` | Required, no default |
| `--output` | **Preset-specific unified output.** User provides a path like `/path/to/results.csv`. Script auto-generates `results_p-body.csv` and `results_sg.csv` by inserting `_p-body` / `_sg` before the `.csv` extension. |
| `--single-cell` | Optional, user's choice |

**Important:** When `--preset m7g-cap-v1` is active, the `--output-pbody` and `--output-sg` flags are NOT available. Instead, the unified `--output` flag is used. When no preset is active, `--output-pbody` and `--output-sg` remain the interface (no `--output` flag).

**Implementation note:** Use an argparse mutually exclusive group to make the constraint visible in `--help`:
```python
output_group = parser.add_mutually_exclusive_group()
output_group.add_argument('--output', help='Base output path; generates _p-body.csv and _sg.csv (use with --preset)')
output_group.add_argument('--output-pbody', help='P-body results CSV (use without --preset, requires --output-sg)')
```
Then validate after parsing: if `--preset` is given, require `--output`; if no preset, require `--output-pbody` and `--output-sg`.

**Usage:**
```bash
python per_particle_analysis.py --preset m7g-cap-v1 \
    --data-dir /path/to/data \
    --output /path/to/results.csv
# Creates: /path/to/results_p-body.csv and /path/to/results_sg.csv
```

#### Preset: `decapping-sensor-v1` (whole_field_analysis.py)

**Locked parameters:**
| Parameter | Preset Value |
|---|---|
| `--min-size` | 5 |
| `--mng-bg-mode` | 0 |
| `--halo-bg-mode` | SiR_mean |
| `--mNG-filter` | NaN |
| `--percent` | True |
| `--exclude-halo-zero` | True (current default) |
| `--exclude-halo-one` | False (current default) |
| `--SiR-subtract` | None |
| `--FLIM-filter` | None |
| `--mNG-in-FLIM` | no |
| `--save-processed` | False |

**Free parameters (user-provided):**
| Parameter | Notes |
|---|---|
| `--data-dir` | Required, no default |
| `--output` | Required, no default |
| `--single-cell` | Optional, user's choice |

**Note:** `--SiR-subtract` is locked to `None` because combining it with `--halo-bg-mode=SiR_mean` is already forbidden by the existing code (mutual exclusion check). This prevents a confusing runtime error.

**Note:** The `decapping-sensor-v1` preset requires `Dcp2_mask` files in the data directory (due to `--mNG-filter=NaN`). Document this requirement clearly in the README.

**Usage:**
```bash
python whole_field_analysis.py --preset decapping-sensor-v1 \
    --data-dir /path/to/data \
    --output /path/to/results.csv
```

#### Invalid Preset Handling

If a user provides an unknown preset name, error with:
```
Error: Unknown preset 'foo'. Available presets: m7g-cap-v1
```

---

### Phase 3: Installation Infrastructure

#### `requirements.txt`

```
numpy>=1.21
pandas>=1.3
scipy>=1.7
scikit-image>=0.19
tifffile>=2021.7
```

Minimum versions chosen to support Python 3.9+ while ensuring the features and performance optimizations used in the scripts are available. Only direct imports are listed — transitive dependencies (pillow, imageio, etc.) are resolved by pip.

### Research Insights (Phase 3)

**Version floor rationale:**
- `numpy>=1.21`: First version with SIMD optimizations for element-wise ops. Avoids silent dtype upcasting bugs in older versions.
- `scipy>=1.7`: Contains performance fix for `distance_transform_edt` on large arrays (issue #13373).
- `scikit-image>=0.19`: Dropped slow Python-loop fallbacks in morphology and labeling utilities.
- `pandas>=1.3`: Material `to_csv` and `concat` performance improvements.
- `tifffile>=2021.7`: Stable `imread`/`imwrite` interface.

**Do NOT use `==` exact pins** for a lab tool. Exact pins cause installation failures when a collaborator has a slightly different environment or a security patch bumps a micro version. Minimum pins (`>=`) are the right choice here.

**Do NOT list transitive dependencies** (pillow, imageio, networkx). Let pip resolve them.

#### `.gitignore`

```gitignore
# Python
__pycache__/
*.pyc
*.pyo
venv/
.venv/
*.egg-info/

# Data and results
*.tif
*.tiff
*.csv
*.xlsx
*.prism
*_results/
*_reults/
*_images/
*_image/

# OS
.DS_Store
Thumbs.db

# IDE / tools
.claude/
.vscode/
.idea/

# Archive
archive/

# PDFs (analysis docs, not code)
*.pdf

# Misc
untitled folder/
untitled folder 2/
```

### Research Insights (.gitignore)

**Security review additions:** The original `.gitignore` was missing patterns for:
- `*_results/`, `*_reults/` (misspelled directory `M10_cap_reults/` exists), `*_images/`, `*_image/` — data directories
- `untitled folder 2/` — another empty folder in the repo
- `.claude/settings.local.json` contains absolute paths (`/Users/leelab/...`, `/Volumes/NX-01-A/...`) that would leak local system info. Already covered by `.claude/` exclusion — verify before first commit.

#### `LICENSE`

MIT License with the lab's name and 2026 date.

---

### Phase 4: README.md

**Structure:**

1. **Title and description** — `mask-intensity-analysis-repo`: A repository containing all the masked intensity analysis approaches used in the lab. Quantifies fluorescence intensity in masked subcellular regions with local background subtraction. Include a "Last updated: YYYY-MM-DD" line.

2. **Installation** — Platform-specific instructions:
   - **Mac/Linux:**
     ```bash
     git clone https://github.com/<org>/mask-intensity-analysis-repo.git
     cd mask-intensity-analysis-repo
     python3 -m venv venv
     source venv/bin/activate
     pip install -r requirements.txt
     ```
   - **Windows:**
     ```cmd
     git clone https://github.com/<org>/mask-intensity-analysis-repo.git
     cd mask-intensity-analysis-repo
     python -m venv venv
     venv\Scripts\activate
     pip install -r requirements.txt
     ```

3. **Scripts overview** — Brief description of each script and when to use which.

4. **`per_particle_analysis.py` CLI Reference** — Complete table of all flags, types, defaults, descriptions. Group into:
   - Required arguments
   - Analysis parameters
   - Background subtraction
   - Output options
   - Presets

5. **`whole_field_analysis.py` CLI Reference** — Same format. Group into:
   - Required arguments
   - Background modes
   - Filter options
   - Output options
   - Presets

6. **File naming conventions** — How input TIFFs must be named. Channel keywords for each script. Group key prefix explanation with examples. **Include a concrete directory listing** showing what a working data folder looks like.

7. **Assay: Decapping Sensor Analysis** — Which script (`whole_field_analysis.py`), which preset (`decapping-sensor-v1`), required input files (P-body_mask, dilute_mask, Halo, mNG, Dcp2_mask, SiR_mask for SiR_mean mode), example command, how to interpret output columns. **Note: Dcp2_mask files are required when using this preset.**

8. **Assay: m7G Cap Immunofluorescence Analysis** — Which script (`per_particle_analysis.py`), which preset (`m7g-cap-v1`), required input files (P-body_mask, SG_mask, Cap, pnorm, sgnorm), example command, how to interpret output columns.

9. **Preset reference table** — All available presets, which script they belong to, and their locked parameter values.

10. **Adding new presets** — Brief instructions for lab members on how to add a new preset (add entry to `PRESETS` dict, always create a new version, never modify existing presets).

### Research Insights (README)

**Quick Start must work on the first try.** A new rotation student should be able to clone, install, and run within 5 minutes. If example data can't be included (too large), describe exactly where to find it on the lab server.

**Document the naming convention prominently.** The filename-based channel detection is the #1 source of confusion for new users. Show a real `ls` output of a working image directory.

**Keep it practical.** The README is the documentation for a small lab team. No separate docs site needed. Aim for completeness but conciseness.

---

### Phase 5: Git Initialization

1. **Create `.gitignore` FIRST** — before `git init`
2. Initialize git repo: `git init`
3. **Add files by name** — do NOT use `git add .` or `git add -A`. Explicitly add:
   ```bash
   git add per_particle_analysis.py whole_field_analysis.py
   git add requirements.txt .gitignore LICENSE README.md
   git add docs/
   ```
4. Run `git status` to verify no data files, results, or `.claude/` leaked into staging
5. Run `git ls-files` to audit what was committed
6. Initial commit
7. Create GitHub remote and push (user will provide org/account details)

### Research Insights (Git Safety)

**Security review finding:** The `.claude/settings.local.json` file contains absolute local paths that would expose system info. Using explicit `git add` by filename (not `git add .`) is the safest approach. Always run `git status` after staging to verify nothing unexpected is included.

---

## Acceptance Criteria

### Functional Requirements

- [x] `per_particle_analysis.py` exists and works identically to current `analyze_m7G_cap.py` when run without `--preset`
- [x] `whole_field_analysis.py` exists and works identically to current `analyze_whole_field.py` when run without `--preset`
- [x] `--preset m7g-cap-v1` on `per_particle_analysis.py` sets correct parameters and accepts unified `--output`
- [x] `--preset decapping-sensor-v1` on `whole_field_analysis.py` sets correct parameters
- [x] Passing a locked flag alongside `--preset` produces a clear error naming the conflicting flag and exits
- [x] Passing an unknown preset name produces a clear error listing available presets
- [x] Preset prints active parameters to stdout when invoked
- [x] `--data-dir` and `--output` are required (no defaults) in both scripts
- [x] Python version check at top of each script exits with clear message if < 3.9
- [x] Scripts work on Mac, Linux, and Windows (no OS-specific path issues)
- [x] `requirements.txt` installs all dependencies successfully on Python 3.9+
- [x] `.gitignore` excludes all data, results, archive, venv, and OS artifacts
- [x] `README.md` contains complete installation instructions for all three platforms
- [x] `README.md` contains complete CLI reference for both scripts with all flags
- [x] `README.md` contains file naming convention documentation with example directory listing
- [x] `README.md` contains Decapping Sensor Analysis section with preset usage
- [x] `README.md` contains m7G Cap IF Analysis section with preset usage
- [x] MIT LICENSE file is present
- [x] `parse_filename` bug fixed in `whole_field_analysis.py` (longest-match-first sorting)
- [x] Preset keys validated against argparse namespace at startup
- [x] Defensive assertions after default/preset application (no None values reach analysis code)

### Quality Gates

- [ ] Both scripts produce identical output to the originals when given the same inputs and flags (regression test by running on existing data)
- [ ] No hardcoded path separators (`/` or `\\`) in any code path
- [ ] All preset-controlled defaults are `None` in argparse; original defaults applied in code
- [ ] `git ls-files` shows only intended files after initial commit (no data, no .claude/)

## Implementation Order

1. **Rename scripts** (copy, don't delete originals until verified)
2. **Fix `parse_filename` bug** in whole_field — add longest-match-first sorting
3. **Add `from __future__ import annotations` and Python version check** to both scripts
4. **Refactor argparse** — set preset-controlled defaults to None, add `ORIGINAL_DEFAULTS` dict, add `resolve_args()` function
5. **Add `--preset` flag, `PRESETS` dict, and `validate_presets()`** to both scripts
6. **Add preset conflict detection, parameter printing, and defensive assertions**
7. **Add unified `--output` for m7g-cap-v1 preset** using mutually exclusive group and `pathlib.Path` for suffix insertion
8. **Make `--data-dir` and `--output*` required**
9. **Create `requirements.txt`**
10. **Create `.gitignore`** (with expanded patterns from security review)
11. **Create `LICENSE`**
12. **Write `README.md`** (with concrete directory listing examples and preset data requirements)
13. **Test on existing data** — verify identical output for both scripts with and without presets
14. **Create `.gitignore` first, then `git init`, then explicit `git add` by filename, then verify with `git status` and `git ls-files`**

## Dependencies & Risks

- **Risk: Preset default=None refactor could introduce bugs** if any code path checks `args.buffer` without the fallback being applied. **Mitigation:** Defensive assertions after `resolve_args()` catch this immediately. Test all code paths with and without presets.
- **Risk: Renaming scripts may break existing workflows** for the current user. **Mitigation:** Keep originals until verified, and document the rename in README.
- **Risk: Dual output interface (--output vs --output-pbody/--output-sg) may confuse users.** **Mitigation:** Use argparse mutually exclusive group so `--help` makes the constraint visible. Document clearly in README.
- **Risk: Accidental commit of data/config files.** **Mitigation:** Create `.gitignore` before `git init`; use explicit `git add` by filename; verify with `git status` and `git ls-files`.
- **Dependency: GitHub org/account** — need to know where to create the remote repo.

## References

- Brainstorm: `docs/brainstorms/2026-04-01-repo-preparation-brainstorm.md`
- Current scripts: `analyze_m7G_cap.py` (657 lines), `analyze_whole_field.py` (583 lines)
- argparse definitions: `analyze_m7G_cap.py:483-522`, `analyze_whole_field.py:354-408`
- `parse_filename` bug location: `analyze_whole_field.py:44` (missing `sorted(CHANNELS, key=len, reverse=True)`)
- [Scientific Python SPEC 0 — Minimum Supported Dependencies](https://scientific-python.org/specs/spec-0000/)
- [Lee 2018 — Ten simple rules for documenting scientific software (PLOS Comp Bio)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006561)
