#!/usr/bin/env python3
"""
Universal group-wise cfDNA size analysis from detailed 5 bp bin files.

This script restores the original cfDNA_group_ks_mwu_analysis-style outputs while
keeping universal group assignment (metadata/group-column/inference):
1) KS-based line plot of group-average size distributions.
2) Boxplots for percent cfDNA by fragment class (grid + combined).
3) Boxplots for mean fragment size by fragment class (grid + combined).
4) CSV summaries for KS and pairwise MWU tests.
"""

import argparse
import glob
import os
import re
import warnings
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import ks_2samp, mannwhitneyu

# ============================================================================
# USER SETTINGS
# ============================================================================
DATA_DIRECTORY = "/mnt/d/coronary_cfDNA/cfDNA_detailed_size_analysis"
OUTPUT_DIRECTORY = "/mnt/d/coronary_cfDNA/group_statistics_results"

FRAGMENT_TYPES = ["mononucleosomal", "dinucleosomal", "trinucleosomal", "hmw"]
FILE_SEARCH_PATTERNS = ["**/{fragment_type}_{suffix}", "**/*{fragment_type}*{suffix}"]

DEFAULT_SAMPLE_METADATA_CSV = None
DEFAULT_METADATA_SAMPLE_COL = "sample"
DEFAULT_METADATA_GROUP_COL = "group"

GROUP_COLUMN_CANDIDATES = ["group", "condition", "cohort", "class", "label"]
SAMPLE_NAME_COLUMN_CANDIDATES = ["bam_file", "sample", "sample_name", "filename", "file", "id"]
GROUP_RULES = [
    # {"label": "brca2", "pattern": r"^brca2[-_]"},
    # {"label": "control", "pattern": r"^(ctrl|control)[-_]"},
    {"label": "stemi", "pattern": r"^cor_.*barcode[_-]?(0[1-5])\b"},
    {"label": "recovery", "pattern": r"^rec_.*barcode[_-]?(0[6-9]|10)\b"},
    {"label": "control", "pattern": r"^con_.*barcode[_-]?(1[1-5])\b"},
]
GROUP_ORDER = ["brca2", "ctrl", "stemi", "recovery", "control", "unknown"]
GROUP_COLOR_MAP = {
    # "brca2": "#d62728",
    # "ctrl": "#1f77b4",
    "stemi": "#d62728",
    "recovery": "#1f77b4",
    "control": "#0C9B30",
    "unknown": "#7f7f7f",
}
DROP_UNKNOWN_GROUP = False
SIZE_MIN = 25
SIZE_MAX = 1000
# ============================================================================


def normalize_runtime_path(path_str: str) -> str:
    path_str = str(path_str)
    if os.name == "nt" and path_str.startswith("/mnt/"):
        parts = path_str.split("/")
        if len(parts) >= 4 and len(parts[2]) == 1:
            drive = parts[2].upper() + ":"
            tail = parts[3:]
            return os.path.join(drive + os.sep, *tail)
    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", path_str):
        drive = path_str[0].lower()
        rest = path_str[2:].replace("\\", "/").lstrip("/")
        return f"/mnt/{drive}/{rest}"
    return path_str


def parse_args():
    parser = argparse.ArgumentParser(
        description="Universal KS + MWU cfDNA group comparison from size-bin tables"
    )
    parser.add_argument("--data-dir", default=DATA_DIRECTORY)
    parser.add_argument("--output-dir", default=OUTPUT_DIRECTORY)
    parser.add_argument("--fragments", default=",".join(FRAGMENT_TYPES))
    parser.add_argument("--drop-unknown", action="store_true", default=DROP_UNKNOWN_GROUP)
    parser.add_argument("--size-min", type=int, default=SIZE_MIN)
    parser.add_argument("--size-max", type=int, default=SIZE_MAX)
    parser.add_argument("--sample-metadata-csv", default=DEFAULT_SAMPLE_METADATA_CSV)
    parser.add_argument("--metadata-sample-col", default=DEFAULT_METADATA_SAMPLE_COL)
    parser.add_argument("--metadata-group-col", default=DEFAULT_METADATA_GROUP_COL)
    return parser.parse_args()


def ensure_output_directory(path):
    os.makedirs(path, exist_ok=True)


def find_matching_file(data_dir, fragment_type, suffix):
    matches = []
    for pattern_template in FILE_SEARCH_PATTERNS:
        pattern = pattern_template.format(fragment_type=fragment_type, suffix=suffix)
        matches.extend(glob.glob(os.path.join(data_dir, pattern), recursive=True))
    matches = sorted(set(matches))
    return matches[0] if matches else None


def get_sample_names(df):
    for col in SAMPLE_NAME_COLUMN_CANDIDATES:
        if col in df.columns:
            return df[col].astype(str).values, col
    return df.index.astype(str).values, "index"


def get_ordered_groups(groups):
    observed = {g for g in groups if pd.notna(g)}
    ordered = [group for group in GROUP_ORDER if group in observed]
    remaining = sorted([group for group in observed if group not in GROUP_ORDER])
    return ordered + remaining


def build_group_color_map(groups):
    color_map = dict(GROUP_COLOR_MAP)
    extras = [g for g in get_ordered_groups(groups) if g not in color_map]
    if extras:
        palette = plt.cm.tab20(np.linspace(0, 1, max(3, len(extras))))
        for idx, group in enumerate(extras):
            rgb = palette[idx % len(palette)][:3]
            color_map[group] = "#%02x%02x%02x" % tuple(int(255 * c) for c in rgb)
    return color_map


def infer_group_from_sample_name(sample_name):
    normalized = os.path.splitext(os.path.basename(str(sample_name)))[0].lower()

    # Handle BRCA2/ctrl cohorts by prefix + barcode/sample-number ranges.
    # brca2: 1-6 and 13-18
    # ctrl: 7-12 and 19-24
    if normalized.startswith(("brca2-", "brca2_")):
        matches = re.findall(r"(?:barcode[_-]?)?(\d{1,2})", normalized)
        if matches:
            num = int(matches[-1])
            if (1 <= num <= 6) or (13 <= num <= 18):
                return "brca2"
        return "brca2"

    if normalized.startswith(("ctrl-", "ctrl_", "control-", "control_")):
        matches = re.findall(r"(?:barcode[_-]?)?(\d{1,2})", normalized)
        if matches:
            num = int(matches[-1])
            if (7 <= num <= 12) or (19 <= num <= 24):
                return "ctrl"
        return "ctrl"

    for rule in GROUP_RULES:
        if re.search(rule["pattern"], normalized):
            return rule["label"]
    return "unknown"


def format_group_label(group: str) -> str:
    if group.lower() == "brca2":
        return "BRCA2-999del5"
    if group.lower() in {"ctrl", "control"}:
        return "control"
    if group.lower() == "stemi":
        return "STEMI"
    return group.capitalize()


def load_group_map(metadata_csv, sample_col, group_col):
    if metadata_csv is None:
        return None
    if not os.path.exists(metadata_csv):
        raise FileNotFoundError(f"Metadata CSV not found: {metadata_csv}")
    meta_df = pd.read_csv(metadata_csv)
    if sample_col not in meta_df.columns or group_col not in meta_df.columns:
        raise ValueError(f"Metadata CSV must contain '{sample_col}' and '{group_col}'")
    out = {}
    for _, row in meta_df[[sample_col, group_col]].dropna().iterrows():
        out[str(row[sample_col]).strip()] = str(row[group_col]).strip().lower()
    return out


def assign_groups(df, sample_names, group_map=None):
    if group_map is not None:
        groups = np.array([group_map.get(str(s).strip(), "unknown") for s in sample_names])
        if len(pd.unique(groups)) > 1:
            return groups, "metadata"

    for col in GROUP_COLUMN_CANDIDATES:
        if col in df.columns:
            groups = df[col].astype(str).str.strip().str.lower().values
            if len(pd.unique(groups)) > 1:
                return groups, f"column:{col}"

    groups = np.array([infer_group_from_sample_name(name) for name in sample_names])
    return groups, "inferred"


def benjamini_hochberg(p_values):
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adjusted[i] = prev
    q = np.empty(n, dtype=float)
    q[order] = np.clip(adjusted, 0.0, 1.0)
    return q


def parse_size_columns(df: pd.DataFrame, fragment_type: str, size_min=None, size_max=None):
    pattern = re.compile(rf"^{re.escape(fragment_type)}_(\d+)bp$")
    col_to_size = {}

    for col in df.columns:
        m = pattern.match(col)
        if not m:
            continue
        size_bp = int(m.group(1))
        if size_min is not None and size_bp < size_min:
            continue
        if size_max is not None and size_bp > size_max:
            continue
        col_to_size[col] = size_bp

    if not col_to_size:
        raise ValueError(f"No size columns found for {fragment_type}")

    ordered_cols = sorted(col_to_size, key=lambda c: col_to_size[c])
    ordered_sizes = np.array([col_to_size[c] for c in ordered_cols], dtype=float)
    return ordered_cols, ordered_sizes


def assign_groups_from_df(df, group_map=None, drop_unknown=False):
    sample_names, sample_source = get_sample_names(df)
    groups, group_source = assign_groups(df, sample_names, group_map=group_map)

    work_df = df.copy()
    work_df["_sample"] = sample_names
    work_df["_group"] = groups

    if drop_unknown:
        work_df = work_df[work_df["_group"] != "unknown"].copy()
        if work_df.empty:
            return None, None, None, None

    unique_groups = [g for g in get_ordered_groups(work_df["_group"].unique()) if pd.notna(g)]
    if len(unique_groups) < 2:
        return None, None, None, None

    return work_df, unique_groups, group_source, sample_source


def load_fragment_tables(data_dir, fragments, suffix, group_map=None, drop_unknown=False):
    tables = {}
    group_source_any = None
    sample_source_any = None

    for fragment_type in fragments:
        file_path = find_matching_file(data_dir, fragment_type, f"size_{suffix}.csv")
        if not file_path:
            raise FileNotFoundError(f"Missing file for {fragment_type}: *size_{suffix}.csv")

        df = pd.read_csv(file_path)
        if df.empty:
            raise ValueError(f"Empty table for {fragment_type}: {file_path}")

        work_df, groups, group_source, sample_source = assign_groups_from_df(
            df,
            group_map=group_map,
            drop_unknown=drop_unknown,
        )
        if work_df is None:
            raise ValueError(f"No usable groups for {fragment_type}: {file_path}")

        work_df = work_df.copy()
        work_df["bam_file"] = work_df["_sample"].astype(str)
        work_df["group"] = work_df["_group"].astype(str)
        tables[fragment_type] = work_df.drop(columns=["_sample", "_group"])
        if group_source_any is None:
            group_source_any = group_source
        if sample_source_any is None:
            sample_source_any = sample_source

    return tables, group_source_any, sample_source_any


def build_distribution_long_for_ks(dist_tables, size_min, size_max):
    parts = []
    for fragment_type, df in dist_tables.items():
        size_cols, size_vals = parse_size_columns(df, fragment_type, size_min=size_min, size_max=size_max)
        sub = df[["bam_file", "group"] + size_cols].copy()
        map_col_size = {col: int(bp) for col, bp in zip(size_cols, size_vals)}
        long_df = sub.melt(
            id_vars=["bam_file", "group"],
            value_vars=size_cols,
            var_name="size_col",
            value_name="percent_cfDNA",
        )
        long_df["size_bp"] = long_df["size_col"].map(map_col_size).astype(int)
        parts.append(long_df.drop(columns=["size_col"]))

    combined = pd.concat(parts, ignore_index=True)
    combined = (
        combined.groupby(["bam_file", "group", "size_bp"], as_index=False)["percent_cfDNA"]
        .sum()
        .sort_values(["group", "bam_file", "size_bp"])
    )
    return combined


def run_ks_test(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ks_stat, p_value = ks_2samp(x, y, alternative="two-sided", mode="exact")

    if any("Exact calculation unsuccessful" in str(w.message) for w in caught):
        ks_stat, p_value = ks_2samp(x, y, alternative="two-sided", mode="asymp")

    if p_value == 0.0:
        p_value = np.finfo(float).tiny
    return ks_stat, p_value


def compute_peak_ks_tests(dist_long, group_order):
    peak_rows = []
    for (bam_file, group), sub in dist_long.groupby(["bam_file", "group"]):
        idx = sub["percent_cfDNA"].idxmax()
        peak_rows.append(
            {
                "bam_file": bam_file,
                "group": group,
                "peak_size_bp": int(sub.loc[idx, "size_bp"]),
                "peak_percent_cfDNA": float(sub.loc[idx, "percent_cfDNA"]),
            }
        )
    peak_df = pd.DataFrame(peak_rows)

    ks_records = []
    for g1, g2 in combinations(group_order, 2):
        x = peak_df.loc[peak_df["group"] == g1, "peak_size_bp"].dropna()
        y = peak_df.loc[peak_df["group"] == g2, "peak_size_bp"].dropna()
        stat, p = run_ks_test(x, y)
        ks_records.append(
            {
                "test": "KS_on_peak_size_distribution",
                "group1": g1,
                "group2": g2,
                "n1": int(len(x)),
                "n2": int(len(y)),
                "statistic": stat,
                "p_value": p,
            }
        )
    return peak_df, pd.DataFrame(ks_records)


def _pick_local_maxima(size_vals, pct_vals, min_peak_bp, top_n, min_sep_bp, min_height_ratio=0.02):
    global_max = float(np.nanmax(pct_vals)) if len(pct_vals) > 0 else 0.0
    height_floor = global_max * min_height_ratio

    candidates = []
    for i in range(1, len(pct_vals) - 1):
        if (
            size_vals[i] >= min_peak_bp
            and pct_vals[i] > pct_vals[i - 1]
            and pct_vals[i] >= pct_vals[i + 1]
            and pct_vals[i] >= height_floor
        ):
            candidates.append((int(size_vals[i]), float(pct_vals[i])))

    if not candidates and len(pct_vals) > 0:
        valid_idx = [i for i in range(len(pct_vals)) if size_vals[i] >= min_peak_bp]
        if not valid_idx:
            valid_idx = list(range(len(pct_vals)))
        top_idx = [valid_idx[i] for i in np.argsort(pct_vals[valid_idx])[-max(1, top_n):]]
        candidates = [(int(size_vals[i]), float(pct_vals[i])) for i in top_idx]

    candidates = sorted(candidates, key=lambda t: t[1], reverse=True)
    selected = []
    for size_bp, pct in candidates:
        if all(abs(size_bp - s) >= min_sep_bp for s, _ in selected):
            selected.append((size_bp, pct))
        if len(selected) >= top_n:
            break
    return selected


def compute_major_peaks(
    dist_long,
    group_order,
    peaks_per_group=5,
    max_total_peaks=5,
    min_sep_bp=70,
    min_peak_bp=25,
    merge_tol_bp=25,
    source_prominence_ratio=0.25,
):
    mean_curve = (
        dist_long.groupby(["group", "size_bp"], as_index=False)["percent_cfDNA"]
        .mean()
        .sort_values(["group", "size_bp"])
    )

    candidates = []
    for group in group_order:
        sub = mean_curve[mean_curve["group"] == group]
        if sub.empty:
            continue

        sizes = sub["size_bp"].to_numpy(dtype=float)
        pcts = sub["percent_cfDNA"].to_numpy(dtype=float)
        global_max = float(np.nanmax(pcts)) if pcts.size else 0.0
        peaks = _pick_local_maxima(
            sizes,
            pcts,
            min_peak_bp=min_peak_bp,
            top_n=peaks_per_group,
            min_sep_bp=min_sep_bp,
        )
        for size_bp, pct in peaks:
            candidates.append(
                {
                    "group": group,
                    "peak_size_bp": int(size_bp),
                    "peak_percent_cfDNA": float(pct),
                }
            )

        # Force-check canonical broad regions so lower-amplitude di/tri humps
        # are preserved when mononucleosomal peaks dominate absolute height.
        window_specs = [
            ("dinuc", 250, 420, max(0.10, 0.015 * global_max)),
            ("trinuc", 430, 650, max(0.06, 0.0075 * global_max)),
        ]
        for _, low_bp, high_bp, min_height in window_specs:
            mask = (sizes >= low_bp) & (sizes <= high_bp)
            if not np.any(mask):
                continue

            window_idx = np.where(mask)[0]
            local_idx = window_idx[int(np.argmax(pcts[window_idx]))]
            local_height = float(pcts[local_idx])
            if local_height < min_height:
                continue

            if 0 < local_idx < (len(pcts) - 1):
                if not (pcts[local_idx] >= pcts[local_idx - 1] and pcts[local_idx] >= pcts[local_idx + 1]):
                    continue

            candidates.append(
                {
                    "group": group,
                    "peak_size_bp": int(sizes[local_idx]),
                    "peak_percent_cfDNA": local_height,
                }
            )

    if not candidates:
        return pd.DataFrame(columns=["peak_id", "peak_size_bp", "peak_percent_cfDNA", "source_groups"])

    clusters = []
    for cand in sorted(candidates, key=lambda r: r["peak_size_bp"]):
        placed = False
        for cluster in clusters:
            if abs(cand["peak_size_bp"] - cluster["center_bp"]) <= merge_tol_bp:
                cluster["members"].append(cand)
                sizes = np.array([m["peak_size_bp"] for m in cluster["members"]], dtype=float)
                weights = np.array([m["peak_percent_cfDNA"] for m in cluster["members"]], dtype=float)
                cluster["center_bp"] = float(np.average(sizes, weights=weights))
                cluster["max_pct"] = max(cluster["max_pct"], cand["peak_percent_cfDNA"])
                cluster["groups"].add(cand["group"])
                cluster["group_max_pct"][cand["group"]] = max(
                    cluster["group_max_pct"].get(cand["group"], 0.0),
                    cand["peak_percent_cfDNA"],
                )
                placed = True
                break
        if not placed:
            clusters.append(
                {
                    "center_bp": float(cand["peak_size_bp"]),
                    "max_pct": float(cand["peak_percent_cfDNA"]),
                    "groups": {cand["group"]},
                    "group_max_pct": {cand["group"]: float(cand["peak_percent_cfDNA"])},
                    "members": [cand],
                }
            )

    def _region(center_bp: float) -> str:
        if center_bp < 240:
            return "mono"
        if center_bp < 430:
            return "dinuc"
        if center_bp < 650:
            return "trinuc"
        return "other"

    clusters = sorted(clusters, key=lambda c: c["max_pct"], reverse=True)
    selected = []

    def add_top_from_region(region_name: str, max_from_region: int):
        region_items = [c for c in clusters if _region(c["center_bp"]) == region_name]
        for cluster in region_items[:max_from_region]:
            if cluster not in selected:
                selected.append(cluster)

    # Preserve dominant mono structure while guaranteeing representation of
    # secondary nucleosomal periodicity when signal exists.
    add_top_from_region("mono", 2)
    add_top_from_region("dinuc", 1)
    add_top_from_region("trinuc", 1)

    for cluster in clusters:
        if len(selected) >= max_total_peaks:
            break
        if cluster not in selected:
            selected.append(cluster)

    clusters = sorted(selected[:max_total_peaks], key=lambda c: c["center_bp"])

    rows = []
    for i, cluster in enumerate(clusters, start=1):
        dominant_groups = sorted(
            [
                g
                for g, g_peak in cluster["group_max_pct"].items()
                if g_peak >= source_prominence_ratio * cluster["max_pct"]
            ]
        )
        if not dominant_groups:
            dominant_groups = sorted(cluster["groups"])

        rows.append(
            {
                "peak_id": f"peak_{i}",
                "peak_size_bp": int(round(cluster["center_bp"])),
                "peak_percent_cfDNA": float(cluster["max_pct"]),
                "source_groups": "/".join(dominant_groups),
            }
        )

    return pd.DataFrame(rows)


def compute_peak_window_ks_tests(dist_long, peak_df, group_order, window_bp=25):
    local_peak_rows = []
    ks_rows = []

    for _, peak_row in peak_df.iterrows():
        peak_id = peak_row["peak_id"]
        center = int(peak_row["peak_size_bp"])

        window = dist_long[
            (dist_long["size_bp"] >= center - window_bp)
            & (dist_long["size_bp"] <= center + window_bp)
        ]

        for (bam_file, group), sub in window.groupby(["bam_file", "group"]):
            idx = sub["percent_cfDNA"].idxmax()
            local_area = float(sub["percent_cfDNA"].sum())
            local_peak_rows.append(
                {
                    "peak_id": peak_id,
                    "window_center_bp": center,
                    "window_bp": int(window_bp),
                    "bam_file": bam_file,
                    "group": group,
                    "local_peak_size_bp": int(sub.loc[idx, "size_bp"]),
                    "local_peak_percent_cfDNA": float(sub.loc[idx, "percent_cfDNA"]),
                    "local_peak_area_percent_cfDNA": local_area,
                }
            )

        local_peak_df = pd.DataFrame(
            [
                row
                for row in local_peak_rows
                if row["peak_id"] == peak_id and row["window_center_bp"] == center
            ]
        )

        metrics = {
            "local_peak_height_percent": "local_peak_percent_cfDNA",
            "local_peak_position_bp": "local_peak_size_bp",
            "local_peak_area_percent": "local_peak_area_percent_cfDNA",
        }

        for metric_name, metric_col in metrics.items():
            for g1, g2 in combinations(group_order, 2):
                x = local_peak_df.loc[local_peak_df["group"] == g1, metric_col].dropna()
                y = local_peak_df.loc[local_peak_df["group"] == g2, metric_col].dropna()
                ks_stat, p_value = run_ks_test(x, y)
                ks_rows.append(
                    {
                        "test": "KS_on_local_peak_feature_distribution",
                        "metric": metric_name,
                        "peak_id": peak_id,
                        "window_center_bp": center,
                        "window_bp": int(window_bp),
                        "group1": g1,
                        "group2": g2,
                        "n1": int(len(x)),
                        "n2": int(len(y)),
                        "statistic": ks_stat,
                        "p_value": p_value,
                    }
                )

    return pd.DataFrame(local_peak_rows), pd.DataFrame(ks_rows)


def build_percent_by_class(dist_tables):
    rows = []
    for fragment_type, df in dist_tables.items():
        size_cols, _ = parse_size_columns(df, fragment_type)
        percent_in_class = df[size_cols].sum(axis=1).to_numpy(dtype=float)
        rows.append(
            pd.DataFrame(
                {
                    "bam_file": df["bam_file"].values,
                    "group": df["group"].values,
                    "fragment_type": fragment_type,
                    "percent_cfDNA": percent_in_class,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def build_mean_size_by_class(count_tables):
    rows = []
    for fragment_type, df in count_tables.items():
        size_cols, size_vals = parse_size_columns(df, fragment_type)
        counts = df[size_cols].to_numpy(dtype=float)
        totals = counts.sum(axis=1)
        weighted = counts @ size_vals
        mean_size = np.divide(
            weighted,
            totals,
            out=np.full(weighted.shape, np.nan, dtype=float),
            where=totals > 0,
        )
        rows.append(
            pd.DataFrame(
                {
                    "bam_file": df["bam_file"].values,
                    "group": df["group"].values,
                    "fragment_type": fragment_type,
                    "mean_size_bp": mean_size,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def compute_mwu_by_fragment(metric_df, value_col, metric_name, group_order):
    records = []
    for fragment_type in FRAGMENT_TYPES:
        sub = metric_df[metric_df["fragment_type"] == fragment_type]
        for g1, g2 in combinations(group_order, 2):
            x = sub.loc[sub["group"] == g1, value_col].dropna()
            y = sub.loc[sub["group"] == g2, value_col].dropna()
            if len(x) == 0 or len(y) == 0:
                stat, p_val = np.nan, np.nan
            else:
                stat, p_val = mannwhitneyu(x, y, alternative="two-sided")
            records.append(
                {
                    "test": f"MWU_{metric_name}",
                    "metric": metric_name,
                    "fragment_type": fragment_type,
                    "group1": g1,
                    "group2": g2,
                    "n1": int(len(x)),
                    "n2": int(len(y)),
                    "statistic": stat,
                    "p_value": p_val,
                }
            )
    out = pd.DataFrame(records)
    if not out.empty:
        out["q_value_bh"] = benjamini_hochberg(out["p_value"].fillna(1.0).values)
    return out


def _format_pvalue_text(p, sci_threshold: float = 1e-5):
    if pd.isna(p):
        return "nan"

    p = float(p)
    if p < sci_threshold:
        return f"{p:.2e}"

    return f"{p:.5f}".rstrip("0").rstrip(".")


def _pvalue_to_stars(p_value: float) -> str:
    if pd.isna(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def add_significance_bars(ax, comparisons, x_positions, data_max, data_min):
    entries = []
    for g1, g2, p_value in comparisons:
        if pd.isna(p_value):
            label = "ns (p=nan)"
        elif p_value < 0.05:
            stars = _pvalue_to_stars(float(p_value))
            p_txt = _format_pvalue_text(p_value)
            label = f"{stars} (p={p_txt})" if stars else f"p={p_txt}"
        else:
            label = f"ns (p={_format_pvalue_text(p_value)})"
        entries.append((g1, g2, label))

    if not entries:
        return

    value_range = max(float(data_max - data_min), 1.0)
    step = 0.11 * value_range
    base_y = float(data_max) + 0.06 * value_range
    text_offset = 0.016 * value_range
    top_y = base_y

    for level, (g1, g2, label) in enumerate(entries):
        x1 = x_positions[g1]
        x2 = x_positions[g2]
        y = base_y + level * step
        h = 0.28 * step
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c="#333333")
        ax.text(
            (x1 + x2) / 2,
            y + h + text_offset,
            label,
            ha="center",
            va="bottom",
            fontsize=12,
            color="#222222",
            fontweight="bold",
        )
        top_y = max(top_y, y + h + 2 * text_offset)

    ax.set_ylim(top=max(ax.get_ylim()[1], top_y + 0.04 * value_range))


def add_pairwise_text_box(ax, stats_df):
    lines = []
    for _, row in stats_df.iterrows():
        lines.append(
            f"{format_group_label(row['group1'])} vs {format_group_label(row['group2'])}: p={_format_pvalue_text(row['p_value'])}"
        )

    if not lines:
        return

    text = "MWU p-values\n" + "\n".join(lines)
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8, "edgecolor": "gray"},
    )


def plot_ks_distribution(
    dist_long,
    major_peak_df,
    peakwise_ks_df,
    out_path,
    size_min,
    size_max,
    group_order,
    palette,
):
    sns.set_theme(style="whitegrid")

    mean_curve = (
        dist_long.groupby(["group", "size_bp"], as_index=False)["percent_cfDNA"]
        .mean()
        .sort_values(["group", "size_bp"])
    )
    sample_counts = (
        dist_long[["bam_file", "group"]]
        .drop_duplicates()
        .groupby("group")
        .size()
        .to_dict()
    )

    fig, ax = plt.subplots(figsize=(13, 7))
    for group in group_order:
        sub = mean_curve[mean_curve["group"] == group]
        if sub.empty:
            continue
        ax.plot(
            sub["size_bp"],
            sub["percent_cfDNA"],
            color=palette.get(group, "#555555"),
            linewidth=2.8,
            label=f"{format_group_label(group)} (n={sample_counts.get(group, 0)})",
        )

    y_max = float(mean_curve["percent_cfDNA"].max()) if not mean_curve.empty else 1.0
    y_top = max(y_max * 1.35, y_max + 1.0)
    ax.set_ylim(bottom=0, top=y_top)

    placed_label_positions = []

    for i, peak_row in major_peak_df.reset_index(drop=True).iterrows():
        peak_id = peak_row["peak_id"]
        center = int(peak_row["peak_size_bp"])
        peak_val = float(peak_row["peak_percent_cfDNA"])
        color = "#444444"

        ax.axvline(center, color=color, linestyle="--", linewidth=1.4, alpha=0.45)

        sub_stats = peakwise_ks_df[
            (peakwise_ks_df["peak_id"] == peak_id)
            & (peakwise_ks_df["window_center_bp"] == center)
            & (peakwise_ks_df["metric"] == "local_peak_height_percent")
        ]
        stat_lines = [
            f"{format_group_label(row['group1'])} vs {format_group_label(row['group2'])}: p={_format_pvalue_text(row['p_value'])}"
            for _, row in sub_stats.iterrows()
        ]
        stat_text = "\n".join(stat_lines) if stat_lines else "No test data"
        source_groups = peak_row.get("source_groups", "")
        source_text = f" [from {source_groups}]" if source_groups else ""

        if center < 170:
            x_text = float(center)
            y_text = float(min(peak_val + 0.32 + 0.12 * i, y_top - 0.10 * y_max))
            ha = "center"
        else:
            x_text = float(center + (18 if center < (size_min + size_max) / 2 else -225))
            y_text = float(min(peak_val + 0.28 + 0.15 * i, y_top - 0.12 * y_max))
            ha = "left"

        x_text = float(np.clip(x_text, size_min + 12, size_max - 220))
        min_dx = 185.0
        min_dy = max(0.60, 0.07 * y_top)

        attempts = 0
        while attempts < 12:
            has_conflict = any(
                abs(x_text - px) < min_dx and abs(y_text - py) < min_dy
                for px, py in placed_label_positions
            )
            if not has_conflict:
                break

            y_text += 0.45
            if y_text > (y_top - 0.08 * y_max):
                y_text = max(peak_val + 0.18, 0.20)
                x_text += 110.0 if (attempts % 2 == 0) else -110.0
                x_text = float(np.clip(x_text, size_min + 12, size_max - 220))
            attempts += 1

        ax.text(
            x_text,
            y_text,
            f"{peak_id} ({center} bp){source_text}\n{stat_text}",
            fontsize=10,
            ha=ha,
            va="bottom",
            linespacing=1.25,
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "alpha": 0.92,
                "edgecolor": color,
            },
        )
        placed_label_positions.append((x_text, y_text))

    ax.set_xlim(size_min, size_max)
    ax.set_xlabel("Fragment size (bp)", fontsize=14, fontweight="bold")
    ax.set_ylabel("% of cfDNA at each size bin", fontsize=14, fontweight="bold")
    ax.set_title(
        "Kolmogorov-Smirnov Group Peak Comparison\n"
        "Average cfDNA size distributions with peak-window KS p-values",
        fontsize=17,
        fontweight="bold",
        pad=12,
    )
    ax.tick_params(axis="x", labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")

    legend = ax.legend(loc="upper right", fontsize=12, title="Group")
    legend.get_title().set_fontweight("normal")
    legend.get_title().set_fontsize(13)
    for text in legend.get_texts():
        txt = text.get_text()
        text.set_fontweight("normal")
        if txt.startswith("BRCA2-999del5"):
            text.set_fontstyle("italic")

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_metric_grid(metric_df, value_col, ylabel, title, out_path, stats_df, group_order, palette):
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))

    for ax, fragment_type in zip(axes.flatten(), FRAGMENT_TYPES):
        sub = metric_df[metric_df["fragment_type"] == fragment_type]
        frag_stats = stats_df[stats_df["fragment_type"] == fragment_type]

        sns.boxplot(
            data=sub,
            x="group",
            y=value_col,
            hue="group",
            order=group_order,
            hue_order=group_order,
            palette=palette,
            showfliers=False,
            width=0.62,
            ax=ax,
            legend=False,
        )
        sns.stripplot(
            data=sub,
            x="group",
            y=value_col,
            order=group_order,
            color="black",
            alpha=0.55,
            jitter=0.15,
            size=4,
            ax=ax,
        )

        ax.set_title(fragment_type, fontsize=14, fontweight="bold")
        ax.set_xlabel("Group", fontsize=13, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")
        ax.set_xticks(np.arange(len(group_order)))
        ax.set_xticklabels([format_group_label(g) for g in group_order])
        ax.tick_params(axis="x", labelrotation=0, labelsize=12)
        ax.tick_params(axis="y", labelsize=12)
        for tick in ax.get_xticklabels():
            tick.set_fontweight("bold")
            if tick.get_text() == "BRCA2-999del5":
                tick.set_fontstyle("italic")
        for tick in ax.get_yticklabels():
            tick.set_fontweight("bold")

        comparisons = [(row["group1"], row["group2"], row["p_value"]) for _, row in frag_stats.iterrows()]
        sub_values = sub[value_col].dropna()
        if not sub_values.empty:
            add_significance_bars(
                ax,
                comparisons,
                {group: idx for idx, group in enumerate(group_order)},
                float(sub_values.max()),
                float(sub_values.min()),
            )

    fig.suptitle(title, fontsize=18, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_metric_combined(metric_df, value_col, ylabel, title, out_path, stats_df, group_order, palette):
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(14, 7))

    sns.boxplot(
        data=metric_df,
        x="fragment_type",
        y=value_col,
        hue="group",
        order=FRAGMENT_TYPES,
        hue_order=group_order,
        palette=palette,
        showfliers=False,
        ax=ax,
    )
    sns.stripplot(
        data=metric_df,
        x="fragment_type",
        y=value_col,
        hue="group",
        order=FRAGMENT_TYPES,
        hue_order=group_order,
        dodge=True,
        alpha=0.3,
        size=2.8,
        palette={group: "black" for group in group_order},
        ax=ax,
        legend=False,
    )

    handles, labels = ax.get_legend_handles_labels()
    if len(handles) >= len(group_order):
        legend = ax.legend(
            handles[: len(group_order)],
            [format_group_label(x) for x in labels[: len(group_order)]],
            title="Group",
            fontsize=12,
            title_fontsize=13,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
        )
        legend.get_title().set_fontweight("bold")
        for text in legend.get_texts():
            text.set_fontweight("bold")
            if text.get_text() == "BRCA2-999del5":
                text.set_fontstyle("italic")

    ax.set_xlabel("Fragment class", fontsize=14, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=14, fontweight="bold")
    ax.set_title(title, fontsize=17, fontweight="bold")
    ax.tick_params(axis="x", labelrotation=0, labelsize=12)
    ax.tick_params(axis="y", labelsize=12)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight("bold")

    # Intentionally no bracket overlays in combined panels to avoid visual clutter.

    fig.tight_layout(rect=(0.0, 0.0, 0.84, 1.0))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    global DATA_DIRECTORY, OUTPUT_DIRECTORY, FRAGMENT_TYPES, DROP_UNKNOWN_GROUP

    args = parse_args()
    DATA_DIRECTORY = normalize_runtime_path(args.data_dir)
    OUTPUT_DIRECTORY = normalize_runtime_path(args.output_dir)
    FRAGMENT_TYPES = [frag.strip() for frag in args.fragments.split(",") if frag.strip()]
    DROP_UNKNOWN_GROUP = bool(args.drop_unknown)
    metadata_csv = normalize_runtime_path(args.sample_metadata_csv) if args.sample_metadata_csv else None

    print("=== Universal cfDNA Group KS + MWU Analysis ===")
    print(f"Data directory: {DATA_DIRECTORY}")
    print(f"Output directory: {OUTPUT_DIRECTORY}")
    print(f"Fragments: {FRAGMENT_TYPES}")
    print(f"Drop unknown groups: {DROP_UNKNOWN_GROUP}")
    print(f"Size range for KS plot: {args.size_min}-{args.size_max} bp")
    print(f"Metadata CSV: {metadata_csv}")

    if not os.path.exists(DATA_DIRECTORY):
        raise FileNotFoundError(f"Data directory not found: {DATA_DIRECTORY}")

    ensure_output_directory(OUTPUT_DIRECTORY)

    group_map = load_group_map(
        metadata_csv,
        sample_col=args.metadata_sample_col,
        group_col=args.metadata_group_col,
    )

    dist_tables, group_source, sample_source = load_fragment_tables(
        DATA_DIRECTORY,
        FRAGMENT_TYPES,
        "distributions",
        group_map=group_map,
        drop_unknown=DROP_UNKNOWN_GROUP,
    )
    count_tables, _, _ = load_fragment_tables(
        DATA_DIRECTORY,
        FRAGMENT_TYPES,
        "counts",
        group_map=group_map,
        drop_unknown=DROP_UNKNOWN_GROUP,
    )

    observed_groups = set()
    for df in dist_tables.values():
        observed_groups.update(df["group"].dropna().unique().tolist())
    group_order = get_ordered_groups(observed_groups)
    if len(group_order) < 2:
        raise ValueError("Need at least two groups with data")
    palette = build_group_color_map(group_order)

    dist_long = build_distribution_long_for_ks(dist_tables, args.size_min, args.size_max)
    peak_df, ks_df = compute_peak_ks_tests(dist_long, group_order)
    major_peak_df = compute_major_peaks(
        dist_long,
        group_order,
        peaks_per_group=5,
        max_total_peaks=5,
        min_sep_bp=70,
        min_peak_bp=25,
        merge_tol_bp=25,
    )
    local_peak_df, peakwise_ks_df = compute_peak_window_ks_tests(
        dist_long,
        major_peak_df,
        group_order,
        window_bp=25,
    )

    percent_df = build_percent_by_class(dist_tables)
    mean_size_df = build_mean_size_by_class(count_tables)

    mwu_percent_df = compute_mwu_by_fragment(percent_df, "percent_cfDNA", "percent_cfDNA", group_order)
    mwu_size_df = compute_mwu_by_fragment(mean_size_df, "mean_size_bp", "mean_size_bp", group_order)

    ks_plot_path = os.path.join(OUTPUT_DIRECTORY, "ks_group_size_distribution_25_1000.png")
    pct_grid_path = os.path.join(OUTPUT_DIRECTORY, "boxplot_percent_cfDNA_by_fragment_class.png")
    pct_combined_path = os.path.join(OUTPUT_DIRECTORY, "boxplot_percent_cfDNA_all_fragment_classes.png")
    size_grid_path = os.path.join(OUTPUT_DIRECTORY, "boxplot_mean_size_by_fragment_class.png")
    size_combined_path = os.path.join(OUTPUT_DIRECTORY, "boxplot_mean_size_all_fragment_classes.png")

    plot_ks_distribution(
        dist_long,
        major_peak_df,
        peakwise_ks_df,
        ks_plot_path,
        args.size_min,
        args.size_max,
        group_order,
        palette,
    )
    plot_metric_grid(
        percent_df,
        value_col="percent_cfDNA",
        ylabel="% of total cfDNA",
        title="Group comparison: % of fragments by class",
        out_path=pct_grid_path,
        stats_df=mwu_percent_df,
        group_order=group_order,
        palette=palette,
    )
    plot_metric_combined(
        percent_df,
        value_col="percent_cfDNA",
        ylabel="% of total cfDNA",
        title="All fragment classes together: % of total cfDNA",
        out_path=pct_combined_path,
        stats_df=mwu_percent_df,
        group_order=group_order,
        palette=palette,
    )
    plot_metric_grid(
        mean_size_df,
        value_col="mean_size_bp",
        ylabel="Mean fragment size (bp)",
        title="Group comparison: mean fragment size by class",
        out_path=size_grid_path,
        stats_df=mwu_size_df,
        group_order=group_order,
        palette=palette,
    )
    plot_metric_combined(
        mean_size_df,
        value_col="mean_size_bp",
        ylabel="Mean fragment size (bp)",
        title="All fragment classes together: mean fragment size",
        out_path=size_combined_path,
        stats_df=mwu_size_df,
        group_order=group_order,
        palette=palette,
    )

    ks_df.to_csv(os.path.join(OUTPUT_DIRECTORY, "ks_peak_size_tests.csv"), index=False)
    peakwise_ks_df.to_csv(os.path.join(OUTPUT_DIRECTORY, "ks_peak_window_tests.csv"), index=False)
    major_peak_df.to_csv(os.path.join(OUTPUT_DIRECTORY, "detected_peak_centers.csv"), index=False)
    peak_df.to_csv(os.path.join(OUTPUT_DIRECTORY, "peak_sizes_by_sample.csv"), index=False)
    local_peak_df.to_csv(os.path.join(OUTPUT_DIRECTORY, "local_peak_sizes_by_sample.csv"), index=False)
    mwu_percent_df.to_csv(os.path.join(OUTPUT_DIRECTORY, "mwu_percent_fragment_class_tests.csv"), index=False)
    mwu_size_df.to_csv(os.path.join(OUTPUT_DIRECTORY, "mwu_mean_size_fragment_class_tests.csv"), index=False)

    summary_df = pd.DataFrame(
        [
            {
                "n_groups": len(group_order),
                "groups": "|".join(group_order),
                "group_source": group_source,
                "sample_name_source": sample_source,
                "n_samples": int(dist_long[["bam_file", "group"]].drop_duplicates().shape[0]),
            }
        ]
    )
    summary_csv = os.path.join(OUTPUT_DIRECTORY, "universal_group_statistics_summary.csv")
    summary_df.to_csv(summary_csv, index=False)

    print("\nGenerated files:")
    print(f"  - {ks_plot_path}")
    print(f"  - {pct_grid_path}")
    print(f"  - {pct_combined_path}")
    print(f"  - {size_grid_path}")
    print(f"  - {size_combined_path}")
    print(f"  - {os.path.join(OUTPUT_DIRECTORY, 'ks_peak_size_tests.csv')}")
    print(f"  - {os.path.join(OUTPUT_DIRECTORY, 'ks_peak_window_tests.csv')}")
    print(f"  - {os.path.join(OUTPUT_DIRECTORY, 'detected_peak_centers.csv')}")
    print(f"  - {os.path.join(OUTPUT_DIRECTORY, 'peak_sizes_by_sample.csv')}")
    print(f"  - {os.path.join(OUTPUT_DIRECTORY, 'local_peak_sizes_by_sample.csv')}")
    print(f"  - {os.path.join(OUTPUT_DIRECTORY, 'mwu_percent_fragment_class_tests.csv')}")
    print(f"  - {os.path.join(OUTPUT_DIRECTORY, 'mwu_mean_size_fragment_class_tests.csv')}")
    print(f"  - {summary_csv}")
    print("\n=== Universal Group Statistics Complete ===")


if __name__ == "__main__":
    main()
