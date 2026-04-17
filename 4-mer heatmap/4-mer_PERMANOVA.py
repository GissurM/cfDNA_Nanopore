from __future__ import annotations

import argparse
import itertools
import json
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.spatial.distance import pdist, squareform
from scipy.stats import mannwhitneyu


DEFAULT_GROUP_PREFIXES = "cor:stemi,con:control,rec:recovery"
DEFAULT_GROUP_ORDER = "stemi,control,recovery"


def parse_mapping(raw: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid mapping '{item}'. Use prefix:group format.")
        key, value = item.split(":", 1)
        out[key.strip().lower()] = value.strip().lower()
    return out


def parse_group_order(raw: str) -> list[str]:
    return [x.strip().lower() for x in raw.split(",") if x.strip()]


def canonical_motifs() -> list[str]:
    return ["".join(motif) for motif in itertools.product("ACGT", repeat=4)]


def bh_fdr(p_values: list[float] | NDArray[np.float64]) -> NDArray[np.float64]:
    p = np.asarray(p_values, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    cumulative = 1.0
    for rank in range(n - 1, -1, -1):
        value = ranked[rank] * n / (rank + 1)
        cumulative = min(cumulative, value)
        adjusted[rank] = cumulative
    output = np.empty(n, dtype=float)
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Pairwise PERMANOVA for coronary 4-mer motif data")
    parser.add_argument("--motif-dir", default="/mnt/d/coronary_cfDNA/fragmentomics_coronary_fixed")
    parser.add_argument("--reference", default="/mnt/c/Users/gissu/Documents/hg38/hg38.fa")
    parser.add_argument("--outdir", default="/mnt/d/coronary_cfDNA/fragmentomics_coronary_fixed/permanova_output")
    parser.add_argument("--permutations", type=int, default=100000)
    parser.add_argument("--random-seed", type=int, default=42)
    parser.add_argument(
        "--group-prefixes",
        default=DEFAULT_GROUP_PREFIXES,
        help="Comma-separated prefix:group mapping, e.g. cor:stemi,con:control,rec:recovery",
    )
    parser.add_argument(
        "--group-order",
        default=DEFAULT_GROUP_ORDER,
        help="Comma-separated preferred group order for reporting",
    )
    parser.add_argument(
        "--focus-group",
        default="stemi",
        help="Group used to select focused shared-motif and driver heatmap summaries",
    )
    return parser.parse_args()


def build_expected_frequencies(reference_path: Path, cache_path: Path) -> pd.DataFrame:
    motifs = canonical_motifs()
    if cache_path.exists():
        cached = pd.read_csv(cache_path)
        if set(cached["motif"]) == set(motifs):
            return cached.sort_values("motif").reset_index(drop=True)

    encoding = {"A": 0, "C": 1, "G": 2, "T": 3}
    mask = (1 << 8) - 1
    counts = np.zeros(256, dtype=np.int64)
    rolling = 0
    valid = 0

    with reference_path.open("r", encoding="ascii", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip().upper()
            if not line:
                continue
            if line.startswith(">"):
                rolling = 0
                valid = 0
                continue
            for base in line:
                value = encoding.get(base)
                if value is None:
                    rolling = 0
                    valid = 0
                    continue
                rolling = ((rolling << 2) | value) & mask
                valid += 1
                if valid >= 4:
                    counts[rolling] += 1

    total = counts.sum()
    if total == 0:
        raise RuntimeError("No valid 4-mers were counted from reference genome")

    expected = pd.DataFrame(
        {
            "motif": motifs,
            "expected_count": counts,
            "expected_frequency": counts / total,
        }
    )
    expected.to_csv(cache_path, index=False)
    return expected


def infer_group(sample_name: str, prefix_map: dict[str, str]) -> str:
    lowered = sample_name.lower()
    for prefix, group in prefix_map.items():
        if lowered.startswith(prefix):
            return group
    return "unknown"


def read_motif_table(path: Path, motif_order: list[str], prefix_map: dict[str, str]) -> pd.DataFrame:
    table = pd.read_csv(path)
    table.columns = [c.strip().lower() for c in table.columns]
    required = {"motif", "count", "frequency"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")

    sample = path.name.replace(".motif.csv", "")
    group = infer_group(sample, prefix_map=prefix_map)

    table["motif"] = table["motif"].str.upper()
    table = table.groupby("motif", as_index=False)[["count", "frequency"]].sum()
    frame = pd.DataFrame({"motif": motif_order})
    merged = frame.merge(table, on="motif", how="left").fillna({"count": 0.0, "frequency": 0.0})
    merged["sample"] = sample
    merged["group"] = group
    return merged


def discover_count_matrix(
    motif_dir: Path,
    motif_order: list[str],
    prefix_map: dict[str, str],
    group_order: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    files = sorted(motif_dir.glob("*.motif.csv"))
    if not files:
        raise FileNotFoundError(f"No .motif.csv files found in {motif_dir}")

    tables = [read_motif_table(path, motif_order, prefix_map=prefix_map) for path in files]
    long_df = pd.concat(tables, ignore_index=True)

    metadata = long_df[["sample", "group"]].drop_duplicates().copy()
    allowed = set(group_order)
    metadata = metadata[metadata["group"].isin(allowed)].copy()
    if metadata.empty:
        raise ValueError("No samples mapped to the configured groups by filename prefix")

    counts = (
        long_df.pivot(index="sample", columns="motif", values="count")
        .loc[metadata["sample"], motif_order]
        .copy()
    )
    return metadata.reset_index(drop=True), counts


def normalize_by_expected(count_matrix: pd.DataFrame, expected: pd.DataFrame) -> pd.DataFrame:
    expected_series = expected.set_index("motif").loc[count_matrix.columns, "expected_frequency"]
    observed = count_matrix.div(count_matrix.sum(axis=1), axis=0)
    oe = observed.div(expected_series, axis=1)
    normalized = oe.div(oe.sum(axis=1), axis=0)
    return normalized * 100.0


def clr_transform(composition: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    proportions = composition / 100.0
    adjusted = proportions + pseudocount
    logged = np.log(adjusted)
    return logged.sub(logged.mean(axis=1), axis=0)


def one_way_permanova(
    distance_matrix: NDArray[np.float64],
    groups: NDArray[np.object_],
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    n = distance_matrix.shape[0]
    unique_groups = np.unique(groups)
    if unique_groups.size < 2:
        raise ValueError("PERMANOVA requires at least two groups")

    condensed = squareform(distance_matrix, checks=False)
    total_ss = (condensed ** 2).sum() / n

    def within_group_ss(labels: NDArray[np.object_]) -> float:
        value = 0.0
        for group in np.unique(labels):
            indices = np.where(labels == group)[0]
            if indices.size < 2:
                continue
            sub = distance_matrix[np.ix_(indices, indices)]
            sub_condensed = squareform(sub, checks=False)
            value += (sub_condensed ** 2).sum() / indices.size
        return value

    observed_within = within_group_ss(groups)
    observed_between = total_ss - observed_within
    df_between = unique_groups.size - 1
    df_within = n - unique_groups.size
    observed_f = (observed_between / df_between) / (observed_within / df_within)

    permutation_stats = np.empty(permutations, dtype=float)
    for i in range(permutations):
        shuffled = rng.permutation(groups)
        perm_within = within_group_ss(shuffled)
        perm_between = total_ss - perm_within
        permutation_stats[i] = (perm_between / df_between) / (perm_within / df_within)

    p_value = (1 + np.sum(permutation_stats >= observed_f)) / (permutations + 1)
    r_squared = observed_between / total_ss
    return {
        "pseudo_f": float(observed_f),
        "p_value": float(p_value),
        "r_squared": float(r_squared),
        "df_between": float(df_between),
        "df_within": float(df_within),
    }


def _one_way_anova_f(values: NDArray[np.float64], labels: NDArray[np.object_]) -> float:
    unique_groups = np.unique(labels)
    n = values.shape[0]
    k = unique_groups.size
    if k < 2 or n <= k:
        return float("nan")

    grand_mean = float(values.mean())
    ss_between = 0.0
    ss_within = 0.0
    for group in unique_groups:
        group_vals = values[labels == group]
        group_mean = float(group_vals.mean())
        ss_between += len(group_vals) * (group_mean - grand_mean) ** 2
        ss_within += float(np.sum((group_vals - group_mean) ** 2))

    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n - k)
    return float(np.inf if ms_within == 0 else ms_between / ms_within)


def permdisp_test(
    clr_values: NDArray[np.float64],
    groups: NDArray[np.object_],
    permutations: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    unique_groups = np.unique(groups)
    distances = np.zeros(clr_values.shape[0], dtype=float)

    for group in unique_groups:
        idx = np.where(groups == group)[0]
        centroid = clr_values[idx].mean(axis=0)
        distances[idx] = np.linalg.norm(clr_values[idx] - centroid, axis=1)

    observed_f = _one_way_anova_f(distances, groups)
    permutation_stats = np.empty(permutations, dtype=float)
    for i in range(permutations):
        shuffled = rng.permutation(groups)
        permutation_stats[i] = _one_way_anova_f(distances, shuffled)

    p_value = float((1 + np.sum(permutation_stats >= observed_f)) / (permutations + 1))
    result = {
        "permdisp_f": float(observed_f),
        "permdisp_p_value": p_value,
    }
    for group in unique_groups:
        idx = np.where(groups == group)[0]
        result[f"permdisp_mean_dist_{group}"] = float(np.mean(distances[idx]))
    return result


def pairwise_permanova_and_motif_counts(
    normalized: pd.DataFrame,
    clr_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    group_order: list[str],
    permutations: int,
    seed: int,
) -> pd.DataFrame:
    pairs = list(itertools.combinations(group_order, 2))
    rows = []

    for i, (group_a, group_b) in enumerate(pairs):
        subset_meta = metadata[metadata["group"].isin([group_a, group_b])].copy()
        subset_samples = subset_meta["sample"].tolist()

        subset_clr = clr_matrix.loc[subset_samples]
        subset_dist = squareform(pdist(subset_clr.to_numpy(dtype=float), metric="euclidean"))
        subset_groups = subset_meta["group"].to_numpy(dtype=object)
        rng = np.random.default_rng(seed + i + 100)
        stat = one_way_permanova(subset_dist, subset_groups, permutations, rng)
        disp = permdisp_test(subset_clr.to_numpy(dtype=float), subset_groups, permutations, np.random.default_rng(seed + i + 500))

        # Motif-level MWU with BH FDR for this pair
        pvals = []
        for motif in normalized.columns:
            a_vals = normalized.loc[subset_meta[subset_meta["group"] == group_a]["sample"], motif]
            b_vals = normalized.loc[subset_meta[subset_meta["group"] == group_b]["sample"], motif]
            _, p = mannwhitneyu(a_vals, b_vals, alternative="two-sided")
            pvals.append(p)
        qvals = bh_fdr(np.asarray(pvals, dtype=float))

        rows.append(
            {
                "comparison": f"{group_a} vs {group_b}",
                "group_a": group_a,
                "group_b": group_b,
                "n_group_a": int((subset_meta["group"] == group_a).sum()),
                "n_group_b": int((subset_meta["group"] == group_b).sum()),
                "pseudo_f": stat["pseudo_f"],
                "p_value": stat["p_value"],
                "r2": stat["r_squared"],
                "motifs_fdr_lt_0_05": int(np.sum(qvals < 0.05)),
                "motifs_fdr_lt_0_10": int(np.sum(qvals < 0.10)),
                "motifs_fdr_lt_0_20": int(np.sum(qvals < 0.20)),
                "motifs_fdr_lt_0_50": int(np.sum(qvals < 0.50)),
                "total_motifs": int(normalized.shape[1]),
                "df_between": stat["df_between"],
                "df_within": stat["df_within"],
                "permdisp_f": disp["permdisp_f"],
                "permdisp_p_value": disp["permdisp_p_value"],
                f"permdisp_mean_dist_{group_a}": disp.get(f"permdisp_mean_dist_{group_a}", np.nan),
                f"permdisp_mean_dist_{group_b}": disp.get(f"permdisp_mean_dist_{group_b}", np.nan),
            }
        )

    out = pd.DataFrame(rows)
    out["p_value_adjusted"] = bh_fdr(out["p_value"].to_numpy(dtype=float))
    return out


def pairwise_motif_results(normalized: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    discovered_groups = sorted(metadata["group"].unique().tolist())
    pairs = list(itertools.combinations(discovered_groups, 2))
    rows = []

    for group_a, group_b in pairs:
        subset_meta = metadata[metadata["group"].isin([group_a, group_b])].copy()
        sample_a = subset_meta[subset_meta["group"] == group_a]["sample"].tolist()
        sample_b = subset_meta[subset_meta["group"] == group_b]["sample"].tolist()
        pvals = []
        motif_stats = []

        for motif in normalized.columns:
            a_vals = normalized.loc[sample_a, motif]
            b_vals = normalized.loc[sample_b, motif]
            stat, p_value = mannwhitneyu(a_vals, b_vals, alternative="two-sided")
            mean_a = float(a_vals.mean())
            mean_b = float(b_vals.mean())
            delta = mean_a - mean_b
            pvals.append(p_value)
            motif_stats.append(
                {
                    "comparison": f"{group_a} vs {group_b}",
                    "group_a": group_a,
                    "group_b": group_b,
                    "motif": motif,
                    "mean_group_a_pct": mean_a,
                    "mean_group_b_pct": mean_b,
                    "delta_pct": delta,
                    "abs_delta_pct": abs(delta),
                    "u_statistic": float(stat),
                    "p_value": float(p_value),
                }
            )

        qvals = bh_fdr(np.asarray(pvals, dtype=float))
        for row, q_value in zip(motif_stats, qvals):
            row["q_value"] = float(q_value)
            row["higher_in"] = row["group_a"] if row["delta_pct"] > 0 else row["group_b"] if row["delta_pct"] < 0 else "tie"
            rows.append(row)

    results = pd.DataFrame(rows)
    return results.sort_values(["comparison", "q_value", "p_value", "abs_delta_pct"], ascending=[True, True, True, False]).reset_index(drop=True)


def build_overlap_summary(motif_results: pd.DataFrame) -> pd.DataFrame:
    comparisons = motif_results["comparison"].drop_duplicates().tolist()
    thresholds = [0.05, 0.10, 0.20, 0.50]
    rows = []

    for i, comparison_a in enumerate(comparisons):
        motifs_a = motif_results[motif_results["comparison"] == comparison_a]
        for comparison_b in comparisons[i + 1 :]:
            motifs_b = motif_results[motif_results["comparison"] == comparison_b]
            set_a_all = set(motifs_a["motif"])
            set_b_all = set(motifs_b["motif"])
            for threshold in thresholds:
                set_a = set(motifs_a.loc[motifs_a["q_value"] < threshold, "motif"])
                set_b = set(motifs_b.loc[motifs_b["q_value"] < threshold, "motif"])
                shared = set_a & set_b
                rows.append(
                    {
                        "comparison_a": comparison_a,
                        "comparison_b": comparison_b,
                        "q_threshold": threshold,
                        "motifs_in_a": int(len(set_a)),
                        "motifs_in_b": int(len(set_b)),
                        "shared_motifs": int(len(shared)),
                        "jaccard_index": float(len(shared) / len(set_a | set_b)) if (set_a or set_b) else 0.0,
                        "shared_motif_list": ";".join(sorted(shared)),
                        "total_motifs": int(len(set_a_all & set_b_all)),
                    }
                )

    return pd.DataFrame(rows)


def build_shared_motif_table(
    motif_results: pd.DataFrame,
    focus_group: str,
    threshold: float = 0.20,
) -> pd.DataFrame:
    focus_rows = motif_results[(motif_results["group_a"] == focus_group) | (motif_results["group_b"] == focus_group)]
    focus_comparisons = focus_rows["comparison"].drop_duplicates().tolist()
    if len(focus_comparisons) < 2:
        return pd.DataFrame(
            columns=[
                "motif",
                "comparison_left",
                "q_value_left",
                "delta_pct_left",
                "higher_in_left",
                "comparison_right",
                "q_value_right",
                "delta_pct_right",
                "higher_in_right",
                "same_direction",
            ]
        )

    left_name, right_name = focus_comparisons[:2]
    left = motif_results[(motif_results["comparison"] == left_name) & (motif_results["q_value"] < threshold)].copy()
    right = motif_results[(motif_results["comparison"] == right_name) & (motif_results["q_value"] < threshold)].copy()

    if left.empty or right.empty:
        return pd.DataFrame(
            columns=[
                "motif",
                "comparison_left",
                "q_value_left",
                "delta_pct_left",
                "higher_in_left",
                "comparison_right",
                "q_value_right",
                "delta_pct_right",
                "higher_in_right",
                "same_direction",
            ]
        )

    merged = left[["motif", "q_value", "delta_pct", "higher_in"]].merge(
        right[["motif", "q_value", "delta_pct", "higher_in"]],
        on="motif",
        how="inner",
        suffixes=("_left", "_right"),
    )
    merged["comparison_left"] = left_name
    merged["comparison_right"] = right_name
    merged["same_direction"] = merged["higher_in_left"] == merged["higher_in_right"]
    merged["mean_abs_delta"] = (
        merged["delta_pct_left"].abs() + merged["delta_pct_right"].abs()
    ) / 2.0
    return merged.sort_values(["same_direction", "mean_abs_delta"], ascending=[False, False]).reset_index(drop=True)


def build_driver_heatmap_table(motif_results: pd.DataFrame, focus_group: str, top_n: int = 40) -> pd.DataFrame:
    focus_rows = motif_results[(motif_results["group_a"] == focus_group) | (motif_results["group_b"] == focus_group)].copy()
    if focus_rows.empty:
        return pd.DataFrame(columns=["motif"])

    ranking = (
        focus_rows.groupby("motif", as_index=False)
        .agg(min_q_value=("q_value", "min"), max_abs_delta=("abs_delta_pct", "max"))
        .sort_values(["min_q_value", "max_abs_delta"], ascending=[True, False])
        .head(top_n)
    )
    selected = ranking["motif"].tolist()

    heatmap = motif_results[motif_results["motif"].isin(selected)].pivot(
        index="motif", columns="comparison", values="delta_pct"
    )
    heatmap = heatmap.reindex(selected)
    heatmap = heatmap.loc[:, sorted(heatmap.columns)]
    return heatmap.reset_index()


def plot_driver_heatmap(heatmap_table: pd.DataFrame, outpath: Path) -> None:
    if heatmap_table.empty:
        return

    plot_df = heatmap_table.set_index("motif")
    values = plot_df.to_numpy(dtype=float)
    vmax = float(np.nanmax(np.abs(values))) if values.size else 0.1
    if vmax == 0:
        vmax = 0.1

    fig_h = max(6.0, 0.34 * len(plot_df) + 1.8)
    fig, ax = plt.subplots(figsize=(7.5, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(plot_df.shape[1]))
    ax.set_xticklabels(plot_df.columns, rotation=20, ha="right")
    ax.set_yticks(np.arange(plot_df.shape[0]))
    ax.set_yticklabels(plot_df.index)
    ax.set_xlabel("Pairwise comparison")
    ax.set_ylabel("Motif")
    ax.set_title("Main 4-mer Motif Drivers Across Pairwise Comparisons")

    for i in range(plot_df.shape[0]):
        for j in range(plot_df.shape[1]):
            value = values[i, j]
            text_color = "white" if abs(value) > 0.55 * vmax else "black"
            ax.text(j, i, f"{value:.02f}", ha="center", va="center", fontsize=7.5, color=text_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Signed delta in normalized motif frequency (%)")

    fig.tight_layout()
    fig.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_driver_heatmap_horizontal(
    heatmap_table: pd.DataFrame,
    outpath: Path,
    max_motifs: int = 24,
) -> None:
    if heatmap_table.empty:
        return

    plot_df = heatmap_table.set_index("motif").copy()
    values = plot_df.to_numpy(dtype=float)
    ranking = np.nanmax(np.abs(values), axis=1)
    top_idx = np.argsort(ranking)[::-1][:max_motifs]
    compact = plot_df.iloc[top_idx].copy()

    # Horizontal layout: comparisons as rows, motifs as columns.
    compact_t = compact.T
    vals = compact_t.to_numpy(dtype=float)
    vmax = float(np.nanmax(np.abs(vals))) if vals.size else 0.1
    if vmax == 0:
        vmax = 0.1

    fig_w = max(12.5, 0.42 * compact_t.shape[1] + 3.8)
    fig_h = max(4.8, 0.9 * compact_t.shape[0] + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(vals, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)

    ax.set_xticks(np.arange(compact_t.shape[1]))
    ax.set_xticklabels(compact_t.columns, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(np.arange(compact_t.shape[0]))
    ax.set_yticklabels(compact_t.index, fontsize=9)
    ax.set_xlabel("Motif")
    ax.set_ylabel("Pairwise comparison")
    ax.set_title(f"Main 4-mer Motif Drivers (Horizontal Compact Top {compact_t.shape[1]})")

    for i in range(compact_t.shape[0]):
        for j in range(compact_t.shape[1]):
            value = vals[i, j]
            text_color = "white" if abs(value) > 0.55 * vmax else "black"
            ax.text(j, i, f"{value:.02f}", ha="center", va="center", fontsize=6.9, color=text_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.03)
    cbar.set_label("Signed delta in normalized motif frequency (%)")

    fig.tight_layout()
    fig.savefig(outpath, dpi=240, bbox_inches="tight")
    plt.close(fig)


def build_pairwise_permanova_permdisp_summary(df: pd.DataFrame) -> pd.DataFrame:
    for col in ["permdisp_f", "permdisp_p_value"]:
        if col not in df.columns:
            df[col] = np.nan

    out = df[
        [
            "comparison",
            "n_group_a",
            "n_group_b",
            "pseudo_f",
            "r2",
            "p_value",
            "p_value_adjusted",
            "permdisp_f",
            "permdisp_p_value",
        ]
    ].copy()

    interpretation = []
    for _, row in out.iterrows():
        p_perm = float(row["p_value"])
        p_disp = float(row["permdisp_p_value"])
        if p_perm < 0.05 and p_disp < 0.05:
            interpretation.append("PERMANOVA sig.; dispersion also differs")
        elif p_perm < 0.05:
            interpretation.append("PERMANOVA sig.; no strong dispersion evidence")
        else:
            interpretation.append("No PERMANOVA significance")
    out["interpretation"] = interpretation
    return out


def plot_pairwise_permanova_permdisp_summary(df: pd.DataFrame, outpath: Path, title_suffix: str) -> None:
    cols = [
        "comparison",
        "n_group_a",
        "n_group_b",
        "pseudo_f",
        "r2",
        "p_value",
        "p_value_adjusted",
        "permdisp_p_value",
        "interpretation",
    ]
    labels = [
        "Comparison",
        "n A",
        "n B",
        "Pseudo-F",
        "R²",
        "PERMANOVA p",
        "Adj. p",
        "PERMDISP p",
        "Interpretation",
    ]

    data = []
    for _, row in df[cols].iterrows():
        data.append(
            [
                str(row["comparison"]),
                str(int(row["n_group_a"])),
                str(int(row["n_group_b"])),
                f"{float(row['pseudo_f']):.3f}",
                f"{float(row['r2']):.3f}",
                f"{float(row['p_value']):.4f}",
                f"{float(row['p_value_adjusted']):.4f}",
                f"{float(row['permdisp_p_value']):.4f}",
                textwrap.fill(str(row["interpretation"]), width=44),
            ]
        )

    fig, ax = plt.subplots(figsize=(21, 5.4))
    ax.axis("off")
    col_widths = [0.12, 0.06, 0.06, 0.09, 0.07, 0.10, 0.08, 0.08, 0.19]
    table = ax.table(cellText=data, colLabels=labels, colWidths=col_widths, cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 2.0)

    for j in range(len(labels)):
        cell = table[0, j]
        cell.set_facecolor("#2c4770")
        cell.set_text_props(color="white", fontweight="bold")

    p_cols = [labels.index("PERMANOVA p"), labels.index("Adj. p"), labels.index("PERMDISP p")]
    interpretation_col = labels.index("Interpretation")
    for i in range(len(data)):
        bg = "#f7f7f7" if i % 2 == 0 else "#ffffff"
        for j in range(len(labels)):
            c = table[i + 1, j]
            c.set_facecolor(bg)
            if j == interpretation_col:
                c._text.set_ha("left")
            if j in p_cols:
                try:
                    pv = float(data[i][j])
                    if pv < 0.05:
                        c.set_facecolor("#c6efce")
                    elif pv < 0.10:
                        c.set_facecolor("#ffeb9c")
                except ValueError:
                    pass

    ax.set_title(f"Pairwise PERMANOVA + PERMDISP Summary — {title_suffix}", fontsize=13, fontweight="bold", pad=14)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_pairwise_permanova_permdisp_summary_html(df: pd.DataFrame, outpath: Path, title_suffix: str) -> None:
    cols = [
        "comparison",
        "n_group_a",
        "n_group_b",
        "pseudo_f",
        "r2",
        "p_value",
        "p_value_adjusted",
        "permdisp_p_value",
        "interpretation",
    ]
    headers = {
        "comparison": "Comparison",
        "n_group_a": "n A",
        "n_group_b": "n B",
        "pseudo_f": "Pseudo-F",
        "r2": "R&sup2;",
        "p_value": "PERMANOVA <i>p</i>",
        "p_value_adjusted": "Adj. <i>p</i>",
        "permdisp_p_value": "PERMDISP <i>p</i>",
        "interpretation": "Interpretation",
    }

    head_html = "".join(f"<th>{headers[c]}</th>" for c in cols)
    body = ""
    for i, (_, row) in enumerate(df.iterrows()):
        bg = "#f7f7f7" if i % 2 == 0 else "#ffffff"
        body += f'<tr style="background:{bg}">'
        for c in cols:
            v = row[c]
            style = ""
            if c in {"p_value", "p_value_adjusted", "permdisp_p_value"}:
                pv = float(v)
                if pv < 0.05:
                    style = ' style="background:#c6efce;font-weight:bold"'
                elif pv < 0.10:
                    style = ' style="background:#ffeb9c"'
                txt = f"{pv:.4f}"
            elif c in {"pseudo_f", "r2"}:
                txt = f"{float(v):.3f}"
            elif c in {"n_group_a", "n_group_b"}:
                txt = str(int(v))
            else:
                txt = str(v)
            body += f"<td{style}>{txt}</td>"
        body += "</tr>\n"

    html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
<meta charset=\"UTF-8\">
<title>Pairwise PERMANOVA + PERMDISP</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #222; }}
  h2 {{ color: #2c4770; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
  th {{ background: #2c4770; color: white; padding: 10px 12px; text-align: center; }}
  td {{ padding: 8px 12px; text-align: center; border-bottom: 1px solid #ddd; }}
  td:first-child, td:last-child {{ text-align: left; }}
  tr:hover {{ background: #eef2ff !important; }}
</style>
</head>
<body>
<h2>Pairwise PERMANOVA + PERMDISP Summary — {title_suffix}</h2>
<table>
  <thead><tr>{head_html}</tr></thead>
  <tbody>
{body}  </tbody>
</table>
</body>
</html>
"""
    outpath.write_text(html, encoding="utf-8")


def plot_table(df: pd.DataFrame, outpath: Path, title_suffix: str) -> None:
    cols = [
        "comparison",
        "pseudo_f",
        "p_value",
        "p_value_adjusted",
        "permdisp_p_value",
        "r2",
        "motifs_fdr_lt_0_05",
        "motifs_fdr_lt_0_10",
        "motifs_fdr_lt_0_20",
        "motifs_fdr_lt_0_50",
        "total_motifs",
    ]
    labels = [
        "Comparison",
        "Pseudo-F",
        "p-value",
        "Adj. p-value",
        "PERMDISP p",
        "R²",
        "Motifs\nq<0.05",
        "Motifs\nq<0.10",
        "Motifs\nq<0.20",
        "Motifs\nq<0.50",
        "Total\nmotifs",
    ]

    data = []
    for _, row in df[cols].iterrows():
        data.append(
            [
                row["comparison"],
                f"{row['pseudo_f']:.3f}",
                f"{row['p_value']:.4f}",
                f"{row['p_value_adjusted']:.4f}",
                f"{row['permdisp_p_value']:.4f}",
                f"{row['r2']:.3f}",
                str(int(row["motifs_fdr_lt_0_05"])),
                str(int(row["motifs_fdr_lt_0_10"])),
                str(int(row["motifs_fdr_lt_0_20"])),
                str(int(row["motifs_fdr_lt_0_50"])),
                str(int(row["total_motifs"])),
            ]
        )

    fig, ax = plt.subplots(figsize=(17, 4.6))
    ax.axis("off")
    table = ax.table(cellText=data, colLabels=labels, cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 2.1)

    for j in range(len(labels)):
        cell = table[0, j]
        cell.set_facecolor("#2c4770")
        cell.set_text_props(color="white", fontweight="bold")

    p_cols = [cols.index("p_value"), cols.index("p_value_adjusted"), cols.index("permdisp_p_value")]
    for i in range(len(data)):
        bg = "#f7f7f7" if i % 2 == 0 else "#ffffff"
        for j in range(len(labels)):
            c = table[i + 1, j]
            c.set_facecolor(bg)
            if j in p_cols:
                pv = float(data[i][j])
                if pv < 0.05:
                    c.set_facecolor("#c6efce")
                elif pv < 0.10:
                    c.set_facecolor("#ffeb9c")

    ax.set_title(f"Pairwise PERMANOVA Results — {title_suffix}", fontsize=13, fontweight="bold", pad=14)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_table_html(df: pd.DataFrame, outpath: Path, title_suffix: str) -> None:
    cols = [
        "comparison",
        "pseudo_f",
        "p_value",
        "p_value_adjusted",
        "permdisp_p_value",
        "r2",
        "motifs_fdr_lt_0_05",
        "motifs_fdr_lt_0_10",
        "motifs_fdr_lt_0_20",
        "motifs_fdr_lt_0_50",
        "total_motifs",
    ]
    headers = {
        "comparison": "Comparison",
        "pseudo_f": "Pseudo-F",
        "p_value": "<i>p</i>-value",
        "p_value_adjusted": "Adj. <i>p</i>-value",
        "permdisp_p_value": "PERMDISP <i>p</i>",
        "r2": "R&sup2;",
        "motifs_fdr_lt_0_05": "Motifs q&lt;0.05",
        "motifs_fdr_lt_0_10": "Motifs q&lt;0.10",
        "motifs_fdr_lt_0_20": "Motifs q&lt;0.20",
        "motifs_fdr_lt_0_50": "Motifs q&lt;0.50",
        "total_motifs": "Total motifs",
    }

    head_html = "".join(f"<th>{headers[c]}</th>" for c in cols)
    body = ""
    for i, (_, row) in enumerate(df.iterrows()):
        bg = "#f7f7f7" if i % 2 == 0 else "#ffffff"
        body += f'<tr style="background:{bg}">'
        for c in cols:
            v = row[c]
            style = ""
            if c in {"p_value", "p_value_adjusted", "permdisp_p_value"}:
                pv = float(v)
                if pv < 0.05:
                    style = ' style="background:#c6efce;font-weight:bold"'
                elif pv < 0.10:
                    style = ' style="background:#ffeb9c"'
                txt = f"{pv:.4f}"
            elif c in {"pseudo_f", "r2"}:
                txt = f"{float(v):.3f}"
            elif c == "comparison":
                txt = str(v)
            else:
                txt = str(int(v))
            body += f"<td{style}>{txt}</td>"
        body += "</tr>\n"

    html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
<meta charset=\"UTF-8\">
<title>Pairwise PERMANOVA</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 40px; color: #222; }}
  h2 {{ color: #2c4770; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
  th {{ background: #2c4770; color: white; padding: 10px 12px; text-align: center; }}
  td {{ padding: 8px 12px; text-align: center; border-bottom: 1px solid #ddd; }}
  td:first-child {{ text-align: left; font-weight: 500; }}
  tr:hover {{ background: #eef2ff !important; }}
</style>
</head>
<body>
<h2>Pairwise PERMANOVA Results — {title_suffix}</h2>
<table>
  <thead><tr>{head_html}</tr></thead>
  <tbody>
{body}  </tbody>
</table>
</body>
</html>
"""
    outpath.write_text(html, encoding="utf-8")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prefix_map = parse_mapping(args.group_prefixes)
    preferred_group_order = parse_group_order(args.group_order)
    if not preferred_group_order:
        raise ValueError("--group-order must contain at least two groups")

    motif_order = canonical_motifs()
    expected = build_expected_frequencies(Path(args.reference), outdir / "hg38_expected_4mer_frequencies.csv")
    metadata, counts = discover_count_matrix(
        Path(args.motif_dir),
        motif_order,
        prefix_map=prefix_map,
        group_order=preferred_group_order,
    )

    discovered = sorted(metadata["group"].unique().tolist())
    group_order = [g for g in preferred_group_order if g in discovered]
    for g in discovered:
        if g not in group_order:
            group_order.append(g)
    if len(group_order) < 2:
        raise RuntimeError("Need at least 2 groups for pairwise comparisons")

    title_suffix = ", ".join(group_order)

    normalized = normalize_by_expected(counts, expected)
    clr = clr_transform(normalized)

    pairwise = pairwise_permanova_and_motif_counts(
        normalized=normalized,
        clr_matrix=clr,
        metadata=metadata,
        group_order=group_order,
        permutations=args.permutations,
        seed=args.random_seed,
    )

    base_cols = [
        "comparison",
        "group_a",
        "group_b",
        "n_group_a",
        "n_group_b",
        "pseudo_f",
        "p_value",
        "p_value_adjusted",
        "permdisp_f",
        "permdisp_p_value",
        "r2",
        "motifs_fdr_lt_0_05",
        "motifs_fdr_lt_0_10",
        "motifs_fdr_lt_0_20",
        "motifs_fdr_lt_0_50",
        "total_motifs",
        "df_between",
        "df_within",
    ]
    dist_cols = sorted([c for c in pairwise.columns if c.startswith("permdisp_mean_dist_")])
    pairwise = pairwise[base_cols[:10] + dist_cols + base_cols[10:]].copy()

    pairwise_summary = build_pairwise_permanova_permdisp_summary(pairwise)

    motif_results = pairwise_motif_results(normalized, metadata)
    overlap_summary = build_overlap_summary(motif_results)
    shared_focus_motifs = build_shared_motif_table(motif_results, focus_group=args.focus_group.lower(), threshold=0.20)
    driver_heatmap_table = build_driver_heatmap_table(motif_results, focus_group=args.focus_group.lower(), top_n=40)

    top_motif_rows = []
    for comparison, sub in motif_results.groupby("comparison"):
        top = sub.sort_values(["q_value", "abs_delta_pct"], ascending=[True, False]).head(25).copy()
        top_motif_rows.append(top)
    top_motifs = pd.concat(top_motif_rows, ignore_index=True)

    metadata.to_csv(outdir / "sample_group_mapping.csv", index=False)
    normalized.to_csv(outdir / "normalized_motif_matrix_pct.csv")
    pairwise.to_csv(outdir / "pairwise_permanova_results.csv", index=False)
    pairwise_summary.to_csv(outdir / "pairwise_permanova_permdisp_summary.csv", index=False)
    motif_results.to_csv(outdir / "pairwise_motif_results.csv", index=False)
    top_motifs.to_csv(outdir / "pairwise_motif_top25_per_comparison.csv", index=False)
    overlap_summary.to_csv(outdir / "pairwise_motif_overlap_summary.csv", index=False)
    shared_focus_motifs.to_csv(outdir / f"shared_motifs_{args.focus_group.lower()}_q0_20.csv", index=False)
    driver_heatmap_table.to_csv(outdir / "pairwise_motif_driver_heatmap_table.csv", index=False)
    plot_table(pairwise, outdir / "pairwise_permanova_results.png", title_suffix=title_suffix)
    write_table_html(pairwise, outdir / "pairwise_permanova_results.html", title_suffix=title_suffix)
    plot_pairwise_permanova_permdisp_summary(
        pairwise_summary,
        outdir / "pairwise_permanova_permdisp_summary.png",
        title_suffix=title_suffix,
    )
    write_pairwise_permanova_permdisp_summary_html(
        pairwise_summary,
        outdir / "pairwise_permanova_permdisp_summary.html",
        title_suffix=title_suffix,
    )
    plot_driver_heatmap(driver_heatmap_table, outdir / "pairwise_motif_driver_heatmap.png")
    plot_driver_heatmap_horizontal(
        driver_heatmap_table,
        outdir / "pairwise_motif_driver_heatmap_horizontal_top24.png",
        max_motifs=24,
    )

    payload = {
        "samples": int(len(metadata)),
        "group_counts": metadata["group"].value_counts().to_dict(),
        "group_order": group_order,
        "pairwise_permanova_results": pairwise.to_dict(orient="records"),
        "pairwise_permanova_permdisp_summary": pairwise_summary.to_dict(orient="records"),
        "pairwise_motif_top25_per_comparison": top_motifs.to_dict(orient="records"),
        "pairwise_motif_overlap_summary": overlap_summary.to_dict(orient="records"),
        "pairwise_motif_driver_heatmap_table": driver_heatmap_table.to_dict(orient="records"),
    }
    (outdir / "pairwise_permanova_results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
