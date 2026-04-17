#!/usr/bin/env python3
"""
Universal group comparison (MWU / Kruskal-Wallis) for cfDNA size-bin tables.

Reads *_size_distributions.csv and compares each size bin across groups.
- If exactly 2 groups: Mann-Whitney U
- If >2 groups: Kruskal-Wallis
"""

import argparse
import glob
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu
from itertools import combinations

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
    {"label": "stemi", "pattern": r"^cor_.*barcode[_-]?(0[1-5])\b"},
    {"label": "recovery", "pattern": r"^rec_.*barcode[_-]?(0[6-9]|10)\b"},
    {"label": "control", "pattern": r"^con_.*barcode[_-]?(1[1-5])\b"},
]
GROUP_ORDER = ["stemi", "recovery", "control", "unknown"]
GROUP_COLOR_MAP = {
    "stemi": "#d62728",
    "recovery": "#1f77b4",
    "control": "#2ca02c",
    "unknown": "#7f7f7f",
}
DROP_UNKNOWN_GROUP = False
TOP_BINS_FOR_PLOTS = 20
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
    parser = argparse.ArgumentParser(description="Universal MWU/Kruskal for cfDNA size-bin groups")
    parser.add_argument("--data-dir", default=DATA_DIRECTORY)
    parser.add_argument("--output-dir", default=OUTPUT_DIRECTORY)
    parser.add_argument("--fragments", default=",".join(FRAGMENT_TYPES))
    parser.add_argument("--drop-unknown", action="store_true", default=DROP_UNKNOWN_GROUP)
    parser.add_argument("--top-bins", type=int, default=TOP_BINS_FOR_PLOTS)
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
    for rule in GROUP_RULES:
        if re.search(rule["pattern"], normalized):
            return rule["label"]
    return "unknown"


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


def run_group_tests(df, fragment_type, group_map=None, drop_unknown=False):
    size_columns = [col for col in df.columns if col.endswith("bp")]
    if len(size_columns) == 0:
        return None, None, None, None, None

    sample_names, sample_source = get_sample_names(df)
    groups, group_source = assign_groups(df, sample_names, group_map=group_map)

    work_df = df.copy()
    work_df["_group"] = groups

    if drop_unknown:
        work_df = work_df[work_df["_group"] != "unknown"].copy()
        if work_df.empty:
            return None, None, None, None, None

    unique_groups = [g for g in get_ordered_groups(work_df["_group"].unique()) if pd.notna(g)]
    if len(unique_groups) < 2:
        return None, None, None, None, None

    mwu_rows = []
    kruskal_rows = []
    for size_bin in size_columns:
        vectors = [work_df.loc[work_df["_group"] == g, size_bin].dropna().values for g in unique_groups]
        if not all(len(v) > 0 for v in vectors):
            continue

        # Always run pairwise MWU for every group pair.
        for i, j in combinations(range(len(unique_groups)), 2):
            stat, p_val = mannwhitneyu(vectors[i], vectors[j], alternative="two-sided")
            med_a = float(np.median(vectors[i]))
            med_b = float(np.median(vectors[j]))
            effect = med_a - med_b
            mwu_rows.append(
                {
                    "fragment_type": fragment_type,
                    "size_bin": size_bin,
                    "test": "mwu",
                    "group_1": unique_groups[i],
                    "group_2": unique_groups[j],
                    "statistic": float(stat),
                    "p_value": float(p_val),
                    "effect_median_diff": effect,
                    "n_group_1": int(len(vectors[i])),
                    "n_group_2": int(len(vectors[j])),
                }
            )

        # Also output a Kruskal-Wallis row when there are more than two groups.
        if len(unique_groups) > 2:
            stat, p_val = kruskal(*vectors)
            group_medians = {g: float(np.median(vectors[i])) for i, g in enumerate(unique_groups)}
            spread = max(group_medians.values()) - min(group_medians.values())
            kruskal_rows.append(
                {
                    "fragment_type": fragment_type,
                    "size_bin": size_bin,
                    "test": "kruskal",
                    "group_1": "|".join(unique_groups),
                    "group_2": "",
                    "statistic": float(stat),
                    "p_value": float(p_val),
                    "effect_median_range": float(spread),
                    "n_total": int(sum(len(v) for v in vectors)),
                }
            )

    if not mwu_rows and not kruskal_rows:
        return None, None, None, None, None

    mwu_df = pd.DataFrame(mwu_rows)
    if not mwu_df.empty:
        mwu_df["q_value_bh"] = benjamini_hochberg(mwu_df["p_value"].values)
        mwu_df = mwu_df.sort_values(["p_value", "size_bin", "group_1", "group_2"]).reset_index(drop=True)

    kruskal_df = pd.DataFrame(kruskal_rows)
    if not kruskal_df.empty:
        kruskal_df["q_value_bh"] = benjamini_hochberg(kruskal_df["p_value"].values)
        kruskal_df = kruskal_df.sort_values(["p_value", "size_bin"]).reset_index(drop=True)

    return mwu_df, kruskal_df, unique_groups, group_source, sample_source


def make_summary_plot(df, stats_df, fragment_type, out_dir, top_bins=20, group_map=None):
    ensure_output_directory(out_dir)

    if stats_df is None or stats_df.empty:
        return

    # Collapse pairwise rows to one p-value per bin using the strongest pairwise signal.
    top_stats = (
        stats_df.groupby("size_bin", as_index=False)["p_value"]
        .min()
        .sort_values("p_value", ascending=True)
        .head(top_bins)
    )
    top_bins_list = top_stats["size_bin"].tolist()
    if not top_bins_list:
        return

    sample_names, _ = get_sample_names(df)
    groups, _ = assign_groups(df, sample_names, group_map=group_map)
    plot_df = df.copy()
    plot_df["_group"] = groups

    top_bins_list = [b for b in top_bins_list if b in plot_df.columns]
    if len(top_bins_list) == 0:
        return

    melted = plot_df[["_group"] + top_bins_list].melt(id_vars="_group", var_name="size_bin", value_name="percentage")
    melted = melted.rename(columns={"_group": "group"})

    ordered_groups = get_ordered_groups(melted["group"].unique())
    palette = build_group_color_map(ordered_groups)

    plt.figure(figsize=(max(10, 0.45 * len(top_bins_list)), 6.5))
    sns.boxplot(data=melted, x="size_bin", y="percentage", hue="group", order=top_bins_list, hue_order=ordered_groups, palette=palette, showfliers=False)
    plt.xticks(rotation=90)
    plt.title(f"{fragment_type}: Top {len(top_bins_list)} significant bins")
    plt.xlabel("Size Bin")
    plt.ylabel("Percentage")
    plt.tight_layout()

    out_png = os.path.join(out_dir, f"{fragment_type}_top_significant_bins.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved plot: {out_png}")


def main():
    global DATA_DIRECTORY, OUTPUT_DIRECTORY, FRAGMENT_TYPES, DROP_UNKNOWN_GROUP, TOP_BINS_FOR_PLOTS

    args = parse_args()
    DATA_DIRECTORY = normalize_runtime_path(args.data_dir)
    OUTPUT_DIRECTORY = normalize_runtime_path(args.output_dir)
    FRAGMENT_TYPES = [frag.strip() for frag in args.fragments.split(",") if frag.strip()]
    DROP_UNKNOWN_GROUP = bool(args.drop_unknown)
    TOP_BINS_FOR_PLOTS = max(1, int(args.top_bins))
    metadata_csv = normalize_runtime_path(args.sample_metadata_csv) if args.sample_metadata_csv else None

    print("=== Universal Group Statistics Analysis ===")
    print(f"Data directory: {DATA_DIRECTORY}")
    print(f"Output directory: {OUTPUT_DIRECTORY}")
    print(f"Fragments: {FRAGMENT_TYPES}")
    print(f"Drop unknown groups: {DROP_UNKNOWN_GROUP}")
    print(f"Metadata CSV: {metadata_csv}")

    if not os.path.exists(DATA_DIRECTORY):
        raise FileNotFoundError(f"Data directory not found: {DATA_DIRECTORY}")

    ensure_output_directory(OUTPUT_DIRECTORY)

    group_map = load_group_map(
        metadata_csv,
        sample_col=args.metadata_sample_col,
        group_col=args.metadata_group_col,
    )

    all_mwu_results = []
    all_kruskal_results = []
    summary_rows = []

    for fragment_type in FRAGMENT_TYPES:
        print("\n" + "=" * 60)
        print(f"Analyzing fragment: {fragment_type}")
        print("=" * 60)

        file_path = find_matching_file(DATA_DIRECTORY, fragment_type, "size_distributions.csv")
        if not file_path:
            print(f"Missing file for {fragment_type}: *_size_distributions.csv")
            continue

        df = pd.read_csv(file_path)
        if df.empty:
            print(f"Skipping {fragment_type}: empty data")
            continue

        mwu_df, kruskal_df, groups, group_source, sample_source = run_group_tests(
            df,
            fragment_type,
            group_map=group_map,
            drop_unknown=DROP_UNKNOWN_GROUP,
        )
        if mwu_df is None or mwu_df.empty:
            print(f"No valid tests for {fragment_type}")
            continue

        mwu_csv = os.path.join(OUTPUT_DIRECTORY, f"{fragment_type}_pairwise_mwu_statistics.csv")
        mwu_df.to_csv(mwu_csv, index=False)
        print(f"Saved pairwise MWU stats: {mwu_csv}")

        if kruskal_df is not None and not kruskal_df.empty:
            kruskal_csv = os.path.join(OUTPUT_DIRECTORY, f"{fragment_type}_kruskal_statistics.csv")
            kruskal_df.to_csv(kruskal_csv, index=False)
            print(f"Saved Kruskal stats: {kruskal_csv}")

        make_summary_plot(
            df,
            mwu_df,
            fragment_type,
            OUTPUT_DIRECTORY,
            top_bins=TOP_BINS_FOR_PLOTS,
            group_map=group_map,
        )

        n_sig_pairwise = int((mwu_df["q_value_bh"] < 0.05).sum())
        n_sig_kruskal = int((kruskal_df["q_value_bh"] < 0.05).sum()) if kruskal_df is not None and not kruskal_df.empty else 0
        summary_rows.append(
            {
                "fragment_type": fragment_type,
                "n_pairwise_mwu_tests": int(len(mwu_df)),
                "n_pairwise_mwu_significant_q_lt_0_05": n_sig_pairwise,
                "n_kruskal_tests": int(len(kruskal_df)) if kruskal_df is not None else 0,
                "n_kruskal_significant_q_lt_0_05": n_sig_kruskal,
                "group_source": group_source,
                "sample_name_source": sample_source,
                "groups": "|".join(groups),
            }
        )
        all_mwu_results.append(mwu_df)
        if kruskal_df is not None and not kruskal_df.empty:
            all_kruskal_results.append(kruskal_df)

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = os.path.join(OUTPUT_DIRECTORY, "universal_group_statistics_summary.csv")
        summary_df.to_csv(summary_csv, index=False)
        print(f"Saved summary: {summary_csv}")

    if all_mwu_results:
        merged_mwu = pd.concat(all_mwu_results, ignore_index=True)
        merged_mwu_csv = os.path.join(OUTPUT_DIRECTORY, "universal_pairwise_mwu_statistics_all_fragments.csv")
        merged_mwu.to_csv(merged_mwu_csv, index=False)
        print(f"Saved merged pairwise MWU stats: {merged_mwu_csv}")

    if all_kruskal_results:
        merged_kruskal = pd.concat(all_kruskal_results, ignore_index=True)
        merged_kruskal_csv = os.path.join(OUTPUT_DIRECTORY, "universal_kruskal_statistics_all_fragments.csv")
        merged_kruskal.to_csv(merged_kruskal_csv, index=False)
        print(f"Saved merged Kruskal stats: {merged_kruskal_csv}")

    print("\n=== Universal Group Statistics Complete ===")


if __name__ == "__main__":
    main()
