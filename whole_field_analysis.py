#!/usr/bin/env python3
"""
Whole-field analysis of fluorescence channel intensities using P-body and
dilute masks. Background subtraction is computed per field from the dilute
mask region.
"""
from __future__ import annotations

import sys
if sys.version_info < (3, 9):
    sys.exit("Error: Python 3.9 or higher is required (for argparse.BooleanOptionalAction).")

import argparse
import glob
import math
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import tifffile
from scipy import stats
from skimage import measure


CHANNELS = ['P-body_mask', 'dilute_mask', 'interaction_mask', 'Dcp2_mask', 'SiR_mask', 'Halo', 'mNG', 'cp_mask']


def parse_bg_mode(value):
    """Parse a bg mode argument: 'mean', 'median', 'mode', 'top_quartile', 'mng-nan', or an integer."""
    if value in ('mean', 'median', 'mode', 'top_quintile', 'top_quartile', 'top_decile',
                 'mng-nan', 'mng-nan-median', 'mng-nan-mode', 'mng-nan-top_quintile',
                 'mng-nan-top_quartile', 'mng-nan-top_decile', 'mng-nan-max',
                 'SiR_mean'):
        return value
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"'{value}' is not 'mean', 'median', 'mode', 'top_quintile', 'top_quartile', 'top_decile', "
            f"'mng-nan', 'mng-nan-median', 'mng-nan-mode', 'mng-nan-top_quintile', "
            f"'mng-nan-top_quartile', 'mng-nan-top_decile', 'mng-nan-max', 'SiR_mean', or an integer"
        )


def parse_filename(filepath):
    """Parse a .tif/.tiff filename by detecting a known channel keyword."""
    basename = os.path.basename(filepath)
    name_no_ext = os.path.splitext(basename)[0]
    for ch in sorted(CHANNELS, key=len, reverse=True):
        idx = name_no_ext.find(ch)
        if idx != -1:
            group_key = name_no_ext[:idx].strip('_')
            return {'group_key': group_key, 'channel': ch}
    return None


def group_image_sets(data_dir):
    """Group .tif/.tiff files by shared name prefix."""
    files = (glob.glob(os.path.join(data_dir, '*.tif'))
             + glob.glob(os.path.join(data_dir, '*.tiff')))
    groups = defaultdict(dict)
    for f in files:
        parsed = parse_filename(f)
        if parsed is None:
            print(f"Warning: could not parse filename (no channel keyword found): "
                  f"{os.path.basename(f)}")
            continue
        groups[parsed['group_key']][parsed['channel']] = f
    return groups


def filter_mask_by_size(mask_img, min_size):
    """Label connected components and remove those <= min_size pixels."""
    binary = mask_img > 0
    labeled = measure.label(binary)
    if min_size <= 0:
        return binary
    filtered = np.zeros_like(binary)
    for prop in measure.regionprops(labeled):
        if prop.area > min_size:
            filtered[labeled == prop.label] = True
    return filtered


def compute_bg_value(image, dilute_mask, bg_mode, exclude_zeros=False, exclude_ones=False):
    """Compute or return the background value for a channel.

    Background values are always returned as integers (rounded up),
    since pixel intensities represent photon counts.
    """
    if isinstance(bg_mode, int):
        return bg_mode
    pixels = image[dilute_mask]
    if exclude_ones:
        pixels = pixels[pixels > 1]
    elif exclude_zeros:
        pixels = pixels[pixels != 0]
    # Remove NaN values so they don't affect calculations
    pixels = pixels[~np.isnan(pixels)]
    if len(pixels) == 0:
        return 0
    if bg_mode == 'mean':
        raw = np.mean(pixels)
    elif bg_mode == 'median':
        raw = np.median(pixels)
    elif bg_mode == 'mode':
        # For integer-valued images, scipy.stats.mode finds the most frequent value.
        result = stats.mode(pixels.astype(int), keepdims=False)
        raw = float(result.mode)
    elif bg_mode == 'top_quintile':
        raw = np.percentile(pixels, 80)
    elif bg_mode == 'top_quartile':
        raw = np.percentile(pixels, 75)
    elif bg_mode == 'top_decile':
        raw = np.percentile(pixels, 90)
    else:
        raw = np.mean(pixels)
    return math.ceil(raw)


def _measure_region(mng_sub, halo_sub, mng_valid, pbody_mask, dilute_mask,
                    compute_percent=False, flim_mask_img=None, mng_mask_img=None):
    """Compute mNG/Halo metrics for given P-body and dilute mask regions."""
    mng_pbody_pixels = mng_sub[pbody_mask]
    mng_dilute_pixels = mng_sub[dilute_mask]

    mng_pbody_mean = np.nanmean(mng_pbody_pixels) if mng_pbody_pixels.size > 0 and np.any(~np.isnan(mng_pbody_pixels)) else np.nan
    mng_pbody_integ = np.nansum(mng_pbody_pixels) if mng_pbody_pixels.size > 0 and np.any(~np.isnan(mng_pbody_pixels)) else np.nan
    mng_dilute_mean = np.nanmean(mng_dilute_pixels) if mng_dilute_pixels.size > 0 and np.any(~np.isnan(mng_dilute_pixels)) else np.nan
    mng_dilute_integ = np.nansum(mng_dilute_pixels) if mng_dilute_pixels.size > 0 and np.any(~np.isnan(mng_dilute_pixels)) else np.nan

    halo_pbody_valid = pbody_mask & mng_valid
    halo_dilute_valid = dilute_mask & mng_valid

    halo_pbody_pixels = halo_sub[halo_pbody_valid]
    halo_dilute_pixels = halo_sub[halo_dilute_valid]

    halo_pbody_mean = np.nanmean(halo_pbody_pixels) if halo_pbody_pixels.size > 0 and np.any(~np.isnan(halo_pbody_pixels)) else np.nan
    halo_pbody_integ = np.nansum(halo_pbody_pixels) if halo_pbody_pixels.size > 0 and np.any(~np.isnan(halo_pbody_pixels)) else np.nan
    halo_dilute_mean = np.nanmean(halo_dilute_pixels) if halo_dilute_pixels.size > 0 and np.any(~np.isnan(halo_dilute_pixels)) else np.nan
    halo_dilute_integ = np.nansum(halo_dilute_pixels) if halo_dilute_pixels.size > 0 and np.any(~np.isnan(halo_dilute_pixels)) else np.nan

    def _ratio(a, b):
        if b and not np.isnan(a) and not np.isnan(b) and b != 0:
            return a / b
        return np.nan

    result = {
        'pbody_area_px': int(pbody_mask.sum()),
        'dilute_area_px': int(dilute_mask.sum()),
        'mNG_pbody_mean': mng_pbody_mean,
        'mNG_dilute_mean': mng_dilute_mean,
        'mNG_pbody_integ': mng_pbody_integ,
        'mNG_dilute_integ': mng_dilute_integ,
        'halo_pbody_mean': halo_pbody_mean,
        'halo_dilute_mean': halo_dilute_mean,
        'halo_pbody_integ': halo_pbody_integ,
        'halo_dilute_integ': halo_dilute_integ,
        'mNG_pbody_over_dilute': _ratio(mng_pbody_mean, mng_dilute_mean),
        'halo_pbody_over_dilute': _ratio(halo_pbody_mean, halo_dilute_mean),
        'halo_over_mNG_pbody': _ratio(halo_pbody_mean, mng_pbody_mean),
        'halo_over_mNG_dilute': _ratio(halo_dilute_mean, mng_dilute_mean),
    }

    if compute_percent:
        overlap_mask = np.ones(halo_sub.shape, dtype=bool)
        if flim_mask_img is not None:
            overlap_mask &= flim_mask_img > 0
        if mng_mask_img is not None:
            overlap_mask &= mng_mask_img > 0

        dcp2_mask = mng_mask_img > 0 if mng_mask_img is not None else np.ones(halo_sub.shape, dtype=bool)

        dcp2_in_pbody = int((pbody_mask & dcp2_mask).sum())
        if dcp2_in_pbody > 0:
            pbody_overlap = int((pbody_mask & overlap_mask).sum())
            result['pct_halo_in_mNG_pbody'] = (pbody_overlap / dcp2_in_pbody) * 100
        else:
            result['pct_halo_in_mNG_pbody'] = np.nan

        dcp2_in_dilute = int((dilute_mask & dcp2_mask).sum())
        if dcp2_in_dilute > 0:
            dilute_overlap = int((dilute_mask & overlap_mask).sum())
            result['pct_halo_in_mNG_dilute'] = (dilute_overlap / dcp2_in_dilute) * 100
        else:
            result['pct_halo_in_mNG_dilute'] = np.nan

    return result


def analyze_image_set(pbody_mask_path, dilute_mask_path, halo_path, mng_path,
                      mng_bg_mode, halo_bg_mode, min_size, exclude_halo_zero=True,
                      exclude_halo_one=False, sir_mask_path=None,
                      sir_subtract_mode=None, flim_mask_path=None,
                      flim_filter_mode=None, mng_in_flim=False,
                      mng_mask_path=None, mng_filter_mode=None,
                      compute_percent=False,
                      save_processed=False, save_dir=None, group_key=None,
                      halo_bg_override=None,
                      single_cell=False, cp_mask_path=None):
    """Analyze one image set and return whole-field measurements."""
    pbody_mask_raw = tifffile.imread(pbody_mask_path)
    dilute_mask_img = tifffile.imread(dilute_mask_path)
    halo_img = tifffile.imread(halo_path).astype(np.float64)
    mng_img = tifffile.imread(mng_path).astype(np.float64)

    # --- SiR subtract: apply to raw Halo data before all other filters ---
    if sir_subtract_mode is not None and sir_mask_path is not None:
        sir_mask_img = tifffile.imread(sir_mask_path)
        sir_inside = sir_mask_img > 0
        if sir_subtract_mode == 'zero':
            halo_img[sir_inside] = 0.0
        elif sir_subtract_mode == 'NaN':
            halo_img[sir_inside] = np.nan
        print(f"  SiR subtract applied: Halo pixels where SiR_mask > 0 set to {sir_subtract_mode}")

    # --- FLIM filter: apply to raw data before background subtraction ---
    if flim_filter_mode is not None and flim_mask_path is not None:
        flim_mask_img = tifffile.imread(flim_mask_path)
        flim_outside = flim_mask_img == 0
        if flim_filter_mode == 'zero':
            halo_img[flim_outside] = 0.0
        elif flim_filter_mode == 'NaN':
            halo_img[flim_outside] = np.nan
        print(f"  FLIM filter applied: Halo pixels outside interaction_mask set to {flim_filter_mode}")

    # --- mNG filter: apply to mNG channel before background subtraction ---
    if mng_filter_mode is not None and mng_mask_path is not None:
        mng_mask_img = tifffile.imread(mng_mask_path)
        mng_outside = mng_mask_img == 0
        if mng_filter_mode == 'zero':
            mng_img[mng_outside] = 0.0
        elif mng_filter_mode == 'NaN':
            mng_img[mng_outside] = np.nan
        print(f"  mNG filter applied: mNG pixels outside Dcp2_mask set to {mng_filter_mode}")

    # --- mNG-in-FLIM: restrict both channels to where interaction_mask AND Dcp2_mask are both > 0 ---
    if mng_in_flim and flim_mask_path is not None and mng_mask_path is not None:
        flim_for_joint = tifffile.imread(flim_mask_path) if flim_filter_mode is None else flim_mask_img
        mng_for_joint = tifffile.imread(mng_mask_path) if mng_filter_mode is None else mng_mask_img
        outside_both = (flim_for_joint == 0) | (mng_for_joint == 0)
        mng_img[outside_both] = np.nan
        halo_img[outside_both] = np.nan
        print(f"  mNG-in-FLIM: mNG and Halo set to NaN outside (interaction_mask & Dcp2_mask)")

    # Filter P-body mask by min size
    pbody_mask = filter_mask_by_size(pbody_mask_raw, min_size)
    dilute_mask = dilute_mask_img > 0

    # --- mNG background subtraction first (Halo 'mng-nan' mode depends on it) ---
    mng_bg = compute_bg_value(mng_img, dilute_mask, mng_bg_mode)
    mng_sub = mng_img - mng_bg
    mng_sub[mng_sub <= 0] = np.nan

    # --- Halo background subtraction ---
    if halo_bg_override is not None:
        halo_bg = halo_bg_override
    elif isinstance(halo_bg_mode, str) and halo_bg_mode.startswith('mng-nan'):
        # Use Halo intensity where mNG is NaN within the dilute mask.
        # These pixels have no detectable mNG signal, so Halo intensity there
        # represents background fluorescence rather than true signal.
        mng_nan_in_dilute = dilute_mask & np.isnan(mng_sub)
        halo_bg_pixels = halo_img[mng_nan_in_dilute]
        if exclude_halo_one:
            halo_bg_pixels = halo_bg_pixels[halo_bg_pixels > 1]
        elif exclude_halo_zero:
            halo_bg_pixels = halo_bg_pixels[halo_bg_pixels != 0]
        if len(halo_bg_pixels) > 0:
            if halo_bg_mode == 'mng-nan':
                halo_bg = math.ceil(np.mean(halo_bg_pixels))
            elif halo_bg_mode == 'mng-nan-median':
                halo_bg = math.ceil(np.median(halo_bg_pixels))
            elif halo_bg_mode == 'mng-nan-mode':
                result = stats.mode(halo_bg_pixels.astype(int), keepdims=False)
                halo_bg = int(result.mode)
            elif halo_bg_mode == 'mng-nan-top_quintile':
                halo_bg = math.ceil(np.percentile(halo_bg_pixels, 80))
            elif halo_bg_mode == 'mng-nan-top_quartile':
                halo_bg = math.ceil(np.percentile(halo_bg_pixels, 75))
            elif halo_bg_mode == 'mng-nan-top_decile':
                halo_bg = math.ceil(np.percentile(halo_bg_pixels, 90))
            elif halo_bg_mode == 'mng-nan-max':
                halo_bg = math.ceil(np.max(halo_bg_pixels))
        else:
            # No NaN mNG pixels in dilute mask — no region qualifies as
            # background, so no correction can be applied.
            halo_bg = 0
            print(f"  Warning: no NaN mNG pixels in dilute mask; "
                  f"Halo background set to 0 (no subtraction)")
    else:
        halo_bg = compute_bg_value(halo_img, dilute_mask, halo_bg_mode,
                                   exclude_zeros=exclude_halo_zero, exclude_ones=exclude_halo_one)

    halo_sub = halo_img - halo_bg
    halo_sub = np.maximum(halo_sub, 0)

    # --- Halo measurements depend on mNG validity ---
    mng_valid = ~np.isnan(mng_sub)

    # Save processed images for troubleshooting
    if save_processed and save_dir is not None:
        prefix = group_key or 'unknown'
        tifffile.imwrite(os.path.join(save_dir, f'{prefix}_mNG_processed.tiff'),
                         mng_sub.astype(np.float32))
        halo_processed = halo_sub.copy()
        halo_processed[~mng_valid] = np.nan
        tifffile.imwrite(os.path.join(save_dir, f'{prefix}_Halo_processed.tiff'),
                         halo_processed.astype(np.float32))
        print(f"  Saved processed images: {prefix}_mNG_processed.tiff, {prefix}_Halo_processed.tiff")

    # Pre-load masks needed for percent computation
    flim_mask_loaded = tifffile.imread(flim_mask_path) if compute_percent and flim_mask_path else None
    mng_mask_loaded = tifffile.imread(mng_mask_path) if compute_percent and mng_mask_path else None

    mng_bg_label = f"manual={mng_bg}" if isinstance(mng_bg_mode, int) else f"{mng_bg_mode}={mng_bg}"
    if halo_bg_override is not None:
        halo_bg_label = f"SiR_mean={halo_bg}"
    elif isinstance(halo_bg_mode, int):
        halo_bg_label = f"manual={halo_bg}"
    else:
        halo_bg_label = f"{halo_bg_mode}={halo_bg}"
    print(f"  P-body mask area: {pbody_mask.sum()} px | Dilute mask area: {dilute_mask.sum()} px")
    print(f"  mNG bg: {mng_bg_label} | Halo bg: {halo_bg_label}")

    # --- Single-cell mode: compute metrics per cell ---
    if single_cell and cp_mask_path is not None:
        cp_mask_img = tifffile.imread(cp_mask_path)
        cell_ids = np.unique(cp_mask_img)
        cell_ids = cell_ids[cell_ids != 0]

        results = []
        for cell_id in cell_ids:
            cell_region = cp_mask_img == cell_id
            cell_pbody = pbody_mask & cell_region
            cell_dilute = dilute_mask & cell_region

            metrics = _measure_region(mng_sub, halo_sub, mng_valid,
                                      cell_pbody, cell_dilute,
                                      compute_percent, flim_mask_loaded, mng_mask_loaded)
            metrics['cell_id'] = int(cell_id)
            metrics['cell_area_px'] = int(cell_region.sum())
            metrics['mng_bg_value'] = mng_bg
            metrics['halo_bg_value'] = halo_bg
            results.append(metrics)

        print(f"  Single-cell: {len(cell_ids)} cells analyzed")
        return results

    # --- Whole-field measurements ---
    result = _measure_region(mng_sub, halo_sub, mng_valid,
                             pbody_mask, dilute_mask,
                             compute_percent, flim_mask_loaded, mng_mask_loaded)
    result['mng_bg_value'] = mng_bg
    result['halo_bg_value'] = halo_bg
    return [result]


# IMPORTANT: Never modify an existing preset. Create a new version instead.
# This ensures published results citing a preset version remain reproducible.
PRESETS = {
    'decapping-sensor-v1': {
        'min_size': 5,
        'mng_bg_mode': 0,
        'halo_bg_mode': 'SiR_mean',
        'mNG_filter': 'NaN',
        'percent': True,
        'exclude_halo_zero': True,
        'exclude_halo_one': False,
        'SiR_subtract': None,
        'FLIM_filter': None,
        'mNG_in_FLIM': 'no',
        'save_processed': False,
    },
}

ORIGINAL_DEFAULTS = {
    'min_size': 10,
    'mng_bg_mode': 'mean',
    'halo_bg_mode': 'median',
    'mNG_filter': None,
    'percent': False,
    'exclude_halo_zero': True,
    'exclude_halo_one': False,
    'SiR_subtract': None,
    'FLIM_filter': None,
    'mNG_in_FLIM': 'no',
    'save_processed': False,
}

PRESET_CONTROLLED_ARGS = set(ORIGINAL_DEFAULTS.keys())


def validate_presets(presets, parser):
    """Verify that all preset keys correspond to valid argparse destinations."""
    valid_dests = {action.dest for action in parser._actions}
    for name, vals in presets.items():
        bad_keys = set(vals.keys()) - valid_dests
        if bad_keys:
            raise ValueError(f"Preset '{name}' references unknown parameters: {bad_keys}")


def resolve_args(args, presets, original_defaults, preset_controlled):
    """Apply preset or original defaults. Error on conflicts with preset."""
    if args.preset:
        preset_name = args.preset
        if preset_name not in presets:
            sys.exit(f"Error: Unknown preset '{preset_name}'. "
                     f"Available presets: {', '.join(presets)}")
        preset_vals = presets[preset_name]
        for key, val in preset_vals.items():
            current = getattr(args, key)
            if current is not None:
                sys.exit(f"Error: --{key.replace('_', '-')} conflicts with "
                         f"--preset {preset_name}. Presets lock all analysis "
                         f"parameters; only --data-dir, --output, and "
                         f"--single-cell can be specified alongside a preset.")
            setattr(args, key, val)
        # Print active preset parameters
        param_str = ', '.join(f'{k}={v}' for k, v in preset_vals.items())
        print(f"[preset: {preset_name}] {param_str}")
    else:
        for key, val in original_defaults.items():
            if getattr(args, key) is None:
                setattr(args, key, val)

    # Defensive assertion: args that should have a non-None value must be set
    # (Some args like SiR_subtract, FLIM_filter legitimately default to None)
    for attr in preset_controlled:
        if original_defaults.get(attr) is not None:
            assert getattr(args, attr) is not None, \
                f"Bug: {attr} is None after default/preset application"

    return args


def main():
    parser = argparse.ArgumentParser(
        description='Whole-field analysis of fluorescence channel intensities using '
                    'P-body and dilute masks. Background subtraction is computed per '
                    'field from the dilute mask region.'
    )
    parser.add_argument('--data-dir', required=True,
                        help='Directory containing .tif/.tiff images')
    parser.add_argument('--output', required=True,
                        help='Output CSV path')

    # --- Preset ---
    parser.add_argument('--preset', choices=list(PRESETS.keys()), default=None,
                        help='Named parameter preset for a specific assay. Locks all analysis '
                             'parameters; only --data-dir, --output, and --single-cell can be '
                             'specified alongside a preset. '
                             f'Available: {", ".join(PRESETS.keys())}')

    # --- Analysis parameters (defaults=None; resolved by preset or ORIGINAL_DEFAULTS) ---
    parser.add_argument('--mng-bg-mode', type=parse_bg_mode, default=None,
                        help='mNG background subtraction: "mean" (default without preset), "median", '
                             '"mode", "top_quintile", "top_quartile", "top_decile", '
                             'or an integer for a manual flat value')
    parser.add_argument('--halo-bg-mode', type=parse_bg_mode, default=None,
                        help='Halo background subtraction: "median" (default without preset), "mean", '
                             '"mode", "top_quintile", "top_quartile", "top_decile", '
                             '"mng-nan" (mean of Halo where mNG is NaN in dilute mask), '
                             '"mng-nan-median", "mng-nan-mode", "mng-nan-top_quintile", '
                             '"mng-nan-top_quartile", "mng-nan-top_decile", "mng-nan-max", '
                             '"SiR_mean" (mean of Halo inside SiR_mask from As condition, '
                             'applied to both As and UT pairs; cannot be combined with --SiR-subtract), '
                             'or an integer for a manual flat value')
    parser.add_argument('--exclude-halo-zero', action=argparse.BooleanOptionalAction, default=None,
                        help='Exclude zero values from Halo channel before estimating background '
                             '(default without preset: True, use --no-exclude-halo-zero to disable)')
    parser.add_argument('--exclude-halo-one', action=argparse.BooleanOptionalAction, default=None,
                        help='Exclude 0 and 1 values from Halo channel before estimating background '
                             '(default without preset: False; when enabled, supersedes --exclude-halo-zero)')
    parser.add_argument('--min-size', type=int, default=None,
                        help='Only include P-body mask components larger than this many pixels '
                             '(default without preset: 10)')
    parser.add_argument('--SiR-subtract', choices=['zero', 'NaN'], default=None,
                        help='Apply SiR mask filtering to Halo channel before all other filters. '
                             '"zero" sets Halo pixels where SiR_mask > 0 to 0; '
                             '"NaN" sets them to NaN. Requires SiR_mask file in the image set. '
                             'Cannot be combined with --halo-bg-mode SiR_mean.')
    parser.add_argument('--FLIM-filter', choices=['zero', 'NaN'], default=None,
                        help='Apply FLIM mask filtering to Halo channel before analysis. '
                             '"zero" sets Halo pixels inside interaction_mask to 0; '
                             '"NaN" sets them to NaN. Requires interaction_mask file in the image set.')
    parser.add_argument('--mNG-filter', choices=['zero', 'NaN'], default=None,
                        help='Apply Dcp2_mask filtering to mNG channel before analysis. '
                             '"zero" sets mNG pixels inside Dcp2_mask to 0; '
                             '"NaN" sets them to NaN. Requires Dcp2_mask file in the image set.')
    parser.add_argument('--mNG-in-FLIM', choices=['yes', 'no'], default=None,
                        help='How to handle mNG channel when --FLIM-filter is used. '
                             '"yes": only analyze mNG and Halo where both interaction_mask and Dcp2_mask > 0; '
                             '"no" (default without preset): do not change how mNG channel is handled.')
    parser.add_argument('--percent', action='store_true', default=None,
                        help='Add columns reporting the percent of mNG-valid pixels '
                             'occupied by Halo, computed separately for P-body and dilute masks')
    parser.add_argument('--save-processed', action='store_true', default=None,
                        help='Save processed mNG and Halo .tif files (after bg subtraction and masking) '
                             'for troubleshooting')

    # --- Other ---
    parser.add_argument('--single-cell', action='store_true', default=False,
                        help='Group all measurements by single cells using a _cp_mask.tiff '
                             'segmentation file. Each integer value in the mask represents one cell. '
                             'Output will have one row per cell per image set.')
    args = parser.parse_args()

    # --- Validate and resolve presets ---
    validate_presets(PRESETS, parser)
    args = resolve_args(args, PRESETS, ORIGINAL_DEFAULTS, PRESET_CONTROLLED_ARGS)

    groups = group_image_sets(args.data_dir)
    if not groups:
        print(f"No image sets found in {args.data_dir}")
        return

    sir_subtract_mode = getattr(args, 'SiR_subtract', None)
    flim_filter_mode = getattr(args, 'FLIM_filter', None)
    mng_filter_mode = getattr(args, 'mNG_filter', None)
    mng_in_flim = getattr(args, 'mNG_in_FLIM', 'no') == 'yes'

    # --- Check for incompatible --SiR-subtract and --halo-bg-mode SiR_mean ---
    if args.halo_bg_mode == 'SiR_mean' and sir_subtract_mode is not None:
        print("WARNING: --halo-bg-mode SiR_mean cannot be combined with --SiR-subtract.")
        print("Please remove --SiR-subtract when using --halo-bg-mode SiR_mean.")
        return

    def _check_channels(group_key, channels):
        """Validate that a group has all required channels. Returns True if valid."""
        required = {'P-body_mask', 'dilute_mask', 'Halo', 'mNG'}
        if not required.issubset(channels.keys()):
            print(f"Skipping {group_key}: missing channels {required - channels.keys()}")
            return False
        if sir_subtract_mode is not None and 'SiR_mask' not in channels:
            print(f"Skipping {group_key}: --SiR-subtract requires SiR_mask file but none found")
            return False
        if flim_filter_mode is not None and 'interaction_mask' not in channels:
            print(f"Skipping {group_key}: --FLIM-filter requires interaction_mask file but none found")
            return False
        if mng_filter_mode is not None and 'Dcp2_mask' not in channels:
            print(f"Skipping {group_key}: --mNG-filter requires Dcp2_mask file but none found")
            return False
        if mng_in_flim and ('interaction_mask' not in channels or 'Dcp2_mask' not in channels):
            print(f"Skipping {group_key}: --mNG-in-FLIM yes requires both interaction_mask and Dcp2_mask")
            return False
        if args.single_cell and 'cp_mask' not in channels:
            print(f"Skipping {group_key}: --single-cell requires cp_mask file but none found")
            return False
        return True

    def _analyze_group(group_key, channels, halo_bg_override=None):
        """Run analyze_image_set for a single group and return a list of result dicts."""
        print(f"Analyzing: {group_key}")
        results = analyze_image_set(
            channels['P-body_mask'], channels['dilute_mask'],
            channels['Halo'], channels['mNG'],
            args.mng_bg_mode, args.halo_bg_mode, args.min_size, args.exclude_halo_zero,
            args.exclude_halo_one,
            sir_mask_path=channels.get('SiR_mask'),
            sir_subtract_mode=sir_subtract_mode,
            flim_mask_path=channels.get('interaction_mask'),
            flim_filter_mode=flim_filter_mode,
            mng_in_flim=mng_in_flim,
            mng_mask_path=channels.get('Dcp2_mask'),
            mng_filter_mode=mng_filter_mode,
            compute_percent=args.percent,
            save_processed=args.save_processed,
            save_dir=args.data_dir,
            group_key=group_key,
            halo_bg_override=halo_bg_override,
            single_cell=args.single_cell,
            cp_mask_path=channels.get('cp_mask'),
        )
        for r in results:
            r['group'] = group_key
        return results

    all_results = []

    if args.halo_bg_mode == 'SiR_mean':
        # --- SiR_mean mode: pair groups by UT/As condition ---
        pairs = {}  # pair_key -> {'As': (group_key, channels), 'UT': (group_key, channels)}
        for group_key, channels in groups.items():
            first_underscore = group_key.index('_') if '_' in group_key else -1
            if first_underscore == -1:
                print(f"Skipping {group_key}: cannot determine UT/As condition (no underscore in name)")
                continue
            condition = group_key[:first_underscore]
            pair_key = group_key[first_underscore + 1:]
            if condition not in ('As', 'UT'):
                print(f"Skipping {group_key}: first token '{condition}' is not 'As' or 'UT'")
                continue
            pairs.setdefault(pair_key, {})[condition] = (group_key, channels)

        for pair_key in sorted(pairs.keys()):
            pair = pairs[pair_key]
            # Check both members exist
            if 'As' not in pair:
                print(f"Skipping UT_{pair_key}: no matching As condition found")
                continue
            if 'UT' not in pair:
                print(f"Skipping As_{pair_key}: no matching UT condition found")
                continue

            as_group_key, as_channels = pair['As']
            ut_group_key, ut_channels = pair['UT']

            # Validate channels for both
            if not _check_channels(as_group_key, as_channels):
                continue
            if not _check_channels(ut_group_key, ut_channels):
                continue

            # Require SiR_mask in the As group
            if 'SiR_mask' not in as_channels:
                print(f"Skipping pair {pair_key}: --halo-bg-mode SiR_mean requires "
                      f"SiR_mask file in As condition but none found for {as_group_key}")
                continue

            # Compute SiR_mean background from As condition's raw Halo
            as_halo_img = tifffile.imread(as_channels['Halo']).astype(np.float64)
            as_sir_mask = tifffile.imread(as_channels['SiR_mask'])
            sir_pixels = as_halo_img[as_sir_mask > 0]
            if len(sir_pixels) == 0:
                print(f"Skipping pair {pair_key}: no pixels where SiR_mask > 0 in {as_group_key}")
                continue
            halo_bg_sir = math.ceil(np.mean(sir_pixels))
            print(f"  SiR_mean background for pair {pair_key}: {halo_bg_sir} "
                  f"(from {as_group_key}, {len(sir_pixels)} SiR pixels)")

            # Analyze both As and UT with the computed background
            all_results.extend(_analyze_group(as_group_key, as_channels, halo_bg_override=halo_bg_sir))
            all_results.extend(_analyze_group(ut_group_key, ut_channels, halo_bg_override=halo_bg_sir))
    else:
        # --- Standard mode: process each group independently ---
        for group_key, channels in sorted(groups.items()):
            if not _check_channels(group_key, channels):
                continue
            all_results.extend(_analyze_group(group_key, channels))

    if not all_results:
        print("No complete image sets found.")
        return

    df = pd.DataFrame(all_results)
    if args.single_cell:
        col_order = ['group', 'cell_id', 'cell_area_px',
                     'pbody_area_px', 'dilute_area_px',
                     'mng_bg_value', 'halo_bg_value',
                     'mNG_pbody_mean', 'mNG_dilute_mean',
                     'mNG_pbody_integ', 'mNG_dilute_integ',
                     'halo_pbody_mean', 'halo_dilute_mean',
                     'halo_pbody_integ', 'halo_dilute_integ',
                     'mNG_pbody_over_dilute', 'halo_pbody_over_dilute',
                     'halo_over_mNG_pbody', 'halo_over_mNG_dilute']
    else:
        col_order = ['group', 'pbody_area_px', 'dilute_area_px',
                     'mng_bg_value', 'halo_bg_value',
                     'mNG_pbody_mean', 'mNG_dilute_mean',
                     'mNG_pbody_integ', 'mNG_dilute_integ',
                     'halo_pbody_mean', 'halo_dilute_mean',
                     'halo_pbody_integ', 'halo_dilute_integ',
                     'mNG_pbody_over_dilute', 'halo_pbody_over_dilute',
                     'halo_over_mNG_pbody', 'halo_over_mNG_dilute']
    if args.percent:
        col_order += ['pct_halo_in_mNG_pbody', 'pct_halo_in_mNG_dilute']
    df = df[col_order]
    df.to_csv(args.output, index=False)
    print(f"\nResults saved to {args.output}")
    if args.single_cell:
        print(f"Total cells analyzed: {len(df)}")
        print(f"\nSummary by group:")
        print(df.groupby('group')[['mNG_pbody_mean', 'mNG_dilute_mean',
                                    'halo_pbody_mean', 'halo_dilute_mean']].mean().to_string())
    else:
        print(f"Total fields analyzed: {len(df)}")
        print(f"\nSummary:")
        print(df[['group', 'mNG_pbody_mean', 'mNG_dilute_mean',
                  'halo_pbody_mean', 'halo_dilute_mean']].to_string(index=False))


if __name__ == '__main__':
    main()
