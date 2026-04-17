#!/usr/bin/env python3
"""
Genome-wide promoter methylation PERMANOVA/PERMDISP for coronary cfDNA.

Groups are inferred from BED file name prefixes:
- con-* -> control
- rec-* -> recovery
- cor-* -> STEMI

Pipeline
1) Build gene-level promoter regions from GENCODE transcript TSS windows.
2) Aggregate methylated/total CpG counts per promoter per sample via bedtools intersect.
3) Build sample x promoter methylation matrix (percent methylation).
4) Filter promoters by per-group minimum observed samples, median-impute remaining sparse NAs.
5) Run global and pairwise PERMANOVA (Euclidean on arcsin-sqrt transformed methylation).
6) Run global and pairwise PERMDISP (distance-to-centroid permutation ANOVA).
7) Export tables + simple QC plots.

Notes
- The SVD-based coordinate compression is exact for Euclidean distances between samples,
  and makes permutation testing much faster when promoters >> samples.
- PERMANOVA significance should be interpreted alongside PERMDISP.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.stats import kruskal, mannwhitneyu


DEFAULT_GENCODE_GTF = "/mnt/c/Users/gissu/Downloads/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
DEFAULT_BED_DIR = "/mnt/d/coronary_cfDNA/mod_BED_coronary"
DEFAULT_OUTPUT_DIR = "/mnt/d/coronary_cfDNA/mod_BED_coronary/permanova_genomewide_output"

DEFAULT_MIN_CPGS = 3
DEFAULT_MIN_SAMPLES = 4
DEFAULT_PERMUTATIONS = 10000
DEFAULT_RANDOM_SEED = 42

DEFAULT_GROUP_ALIASES = {
    "con": "control",
    "rec": "recovery",
    "cor": "stemi",
}
DEFAULT_GROUP_LABELS = {
    "control": "Control",
    "recovery": "Recovery",
    "stemi": "STEMI",
}
DEFAULT_GROUP_ORDER = ["control", "recovery", "stemi"]


def format_comparison_label(comp: str, group_labels: Dict[str, str]) -> str:
    parts = comp.split("_vs_")
    if len(parts) < 2:
        return comp
    return " vs ".join(group_labels.get(p, p) for p in parts)


def parse_key_value_arg(raw: str, lower_values: bool = True) -> Dict[str, str]:
    """Parse comma-separated key:value entries into a dict."""
    out: Dict[str, str] = {}
    if not raw:
        return out
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid mapping item '{item}'. Expected key:value format.")
        key, value = item.split(":", 1)
        key = key.strip().lower()
        value = value.strip().lower() if lower_values else value.strip()
        if not key or not value:
            raise ValueError(f"Invalid mapping item '{item}'. Empty key or value.")
        out[key] = value
    return out


def parse_group_order_arg(raw: str) -> List[str]:
    if not raw:
        return []
    return [x.strip().lower() for x in raw.split(",") if x.strip()]


@dataclass
class SampleInfo:
    sample_id: str
    barcode: str
    group: str
    bed_file: str


def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_field.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip().strip('"')
    return attrs


class GencodePromoterCatalog:
    """Build merged gene-level promoter coordinates from transcript TSS windows."""

    def __init__(self, gtf_path: str):
        self.gtf_path = Path(gtf_path)
        self._catalog: Optional[Dict[str, Dict[str, object]]] = None

    def build_catalog(self) -> Dict[str, Dict[str, object]]:
        if self._catalog is not None:
            return self._catalog

        if not self.gtf_path.exists():
            raise FileNotFoundError(f"GENCODE GTF not found: {self.gtf_path}")

        opener = gzip.open if self.gtf_path.suffix == ".gz" else open
        gene_data: Dict[str, Dict[str, object]] = {}

        with opener(self.gtf_path, "rt") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue

                chrom, _, feature, start, end, _, strand, _, attrs = parts
                if feature != "transcript":
                    continue

                attr_map = parse_gtf_attributes(attrs)
                gene_name = attr_map.get("gene_name")
                if not gene_name:
                    continue

                try:
                    start_i = int(start)
                    end_i = int(end)
                except ValueError:
                    continue

                tss = start_i if strand == "+" else end_i
                if strand == "+":
                    promoter_start = max(1, tss - 2000)
                    promoter_end = tss + 500
                else:
                    promoter_start = max(1, tss - 500)
                    promoter_end = tss + 2000

                chrom_name = chrom if chrom.startswith("chr") else f"chr{chrom}"
                if gene_name not in gene_data:
                    gene_data[gene_name] = {
                        "gene_name": gene_name,
                        "gene_id": attr_map.get("gene_id", "N/A"),
                        "chromosome": chrom_name,
                        "strand": strand,
                        "promoter_start": promoter_start,
                        "promoter_end": promoter_end,
                        "transcript_count": 1,
                    }
                else:
                    gene_data[gene_name]["promoter_start"] = min(
                        int(gene_data[gene_name]["promoter_start"]), promoter_start
                    )
                    gene_data[gene_name]["promoter_end"] = max(
                        int(gene_data[gene_name]["promoter_end"]), promoter_end
                    )
                    gene_data[gene_name]["transcript_count"] = int(
                        gene_data[gene_name]["transcript_count"]
                    ) + 1

        self._catalog = gene_data
        return self._catalog


def aggregate_sample_with_bedtools(task: Tuple[str, str, int]) -> Tuple[np.ndarray, np.ndarray]:
    """Aggregate per-promoter methylated and total CpGs for one sample using bedtools."""
    bed_file, promoter_bed, gene_count = task
    meth_counts = np.zeros(gene_count, dtype=np.int64)
    total_counts = np.zeros(gene_count, dtype=np.int64)

    cmd = [
        "bedtools",
        "intersect",
        "-a",
        bed_file,
        "-b",
        promoter_bed,
        "-wa",
        "-wb",
        "-f",
        "1.0",
    ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 17:
            continue
        try:
            total_cpg = int(float(parts[9]))
            meth_cpg = int(float(parts[11]))
            gene_idx = int(parts[-1])
        except ValueError:
            continue

        if 0 <= gene_idx < gene_count:
            meth_counts[gene_idx] += meth_cpg
            total_counts[gene_idx] += total_cpg

    stderr_text = proc.stderr.read() if proc.stderr is not None else ""
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"bedtools intersect failed for {bed_file}: {stderr_text.strip()}")

    return meth_counts, total_counts


def discover_samples(bed_dir: Path, group_aliases: Dict[str, str]) -> List[SampleInfo]:
    samples: List[SampleInfo] = []
    bed_files = sorted(bed_dir.glob("*barcode*.bed"))
    if not bed_files:
        raise FileNotFoundError(f"No .bed files found in {bed_dir}")

    for bed_file in bed_files:
        m = re.match(r"([A-Za-z]+)[_\-]barcode(\d+)(?:[_\-]\d+)?$", bed_file.stem, re.IGNORECASE)
        if not m:
            continue
        prefix = m.group(1).lower()
        barcode = m.group(2)
        group = group_aliases.get(prefix, prefix)

        sample_id = f"{group}-barcode{barcode}"
        samples.append(SampleInfo(sample_id=sample_id, barcode=barcode, group=group, bed_file=str(bed_file)))

    if not samples:
        raise RuntimeError("No files matched expected prefixes con/rec/cor with barcode naming.")

    return samples


def write_promoter_bed(outdir: Path, promoter_items: List[Tuple[str, Dict[str, object]]]) -> Path:
    out_path = outdir / "_tmp_promoters_permanova.bed"
    with open(out_path, "w", encoding="utf-8") as handle:
        for idx, (_, info) in enumerate(promoter_items):
            chrom = str(info["chromosome"])
            start = int(info["promoter_start"])
            end = int(info["promoter_end"])
            handle.write(f"{chrom}\t{start}\t{end}\tgene_{idx}\t{idx}\n")
    return out_path


def aggregate_all_samples(
    samples: Sequence[SampleInfo],
    promoter_bed: Path,
    gene_count: int,
    max_workers: int,
) -> Tuple[np.ndarray, np.ndarray]:
    tasks = [(s.bed_file, str(promoter_bed), gene_count) for s in samples]
    meth_matrix = np.zeros((len(samples), gene_count), dtype=np.int64)
    total_matrix = np.zeros((len(samples), gene_count), dtype=np.int64)

    workers = max(1, min(max_workers, len(samples)))
    print(f"Aggregating {len(samples)} samples with {workers} workers")
    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_idx = {
            executor.submit(aggregate_sample_with_bedtools, task): i for i, task in enumerate(tasks)
        }
        done = 0
        for future in as_completed(future_to_idx):
            i = future_to_idx[future]
            meth_counts, total_counts = future.result()
            meth_matrix[i, :] = meth_counts
            total_matrix[i, :] = total_counts
            done += 1
            print(f"  aggregated {done}/{len(samples)}: {Path(samples[i].bed_file).name}")

    return meth_matrix, total_matrix


def filter_promoters_by_group_coverage(
    total_matrix: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    min_cpgs: int,
    min_samples: int,
) -> np.ndarray:
    groups_arr = np.asarray(groups)
    keep = np.ones(total_matrix.shape[1], dtype=bool)
    for group in group_order:
        idx = np.where(groups_arr == group)[0]
        if len(idx) == 0:
            continue
        observed = (total_matrix[idx, :] >= min_cpgs).sum(axis=0)
        keep &= observed >= min_samples
    return keep


def median_impute_by_promoter(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Median-impute NaNs per promoter; return imputed matrix and per-promoter NA fraction."""
    x = values.copy()
    na_mask = np.isnan(x)
    na_fraction = na_mask.mean(axis=0)

    medians = np.nanmedian(x, axis=0)
    medians = np.where(np.isnan(medians), 0.0, medians)
    rows, cols = np.where(na_mask)
    x[rows, cols] = medians[cols]
    return x, na_fraction


def to_euclidean_preserving_coords(x: np.ndarray) -> np.ndarray:
    """Project samples to low dimension that preserves pairwise Euclidean distances exactly."""
    xc = x - x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(xc, full_matrices=False)
    tol = 1e-12
    r = int(np.sum(s > tol))
    if r == 0:
        raise RuntimeError("All samples are numerically identical after preprocessing.")
    return u[:, :r] * s[:r]


def _ss_from_labels(coords: np.ndarray, labels: np.ndarray) -> Tuple[float, float, int, int]:
    unique = np.unique(labels)
    n = coords.shape[0]
    k = len(unique)
    if k < 2:
        raise ValueError("Need at least two groups.")
    if n <= k:
        raise ValueError("Need more samples than groups for pseudo-F denominator.")

    grand_centroid = coords.mean(axis=0)
    ss_between = 0.0
    ss_within = 0.0
    for g in unique:
        idx = labels == g
        group_coords = coords[idx]
        centroid = group_coords.mean(axis=0)
        ss_between += group_coords.shape[0] * float(np.sum((centroid - grand_centroid) ** 2))
        ss_within += float(np.sum((group_coords - centroid) ** 2))

    return ss_between, ss_within, n, k


def permanova(coords: np.ndarray, labels: Sequence[str], permutations: int, rng: np.random.Generator) -> Dict[str, float]:
    labels_arr = np.asarray(labels)
    ss_between, ss_within, n, k = _ss_from_labels(coords, labels_arr)
    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n - k)
    pseudo_f = np.inf if ms_within == 0 else ms_between / ms_within
    r2 = ss_between / (ss_between + ss_within) if (ss_between + ss_within) > 0 else 0.0

    if permutations <= 0:
        p_value = 1.0
    else:
        count = 0
        for _ in range(permutations):
            perm_labels = rng.permutation(labels_arr)
            p_ss_between, p_ss_within, p_n, p_k = _ss_from_labels(coords, perm_labels)
            p_ms_between = p_ss_between / (p_k - 1)
            p_ms_within = p_ss_within / (p_n - p_k)
            p_f = np.inf if p_ms_within == 0 else p_ms_between / p_ms_within
            if p_f >= pseudo_f:
                count += 1
        p_value = (count + 1) / (permutations + 1)

    return {
        "pseudo_f": float(pseudo_f),
        "r_squared": float(r2),
        "p_value": float(p_value),
        "n_samples": int(n),
        "n_groups": int(k),
        "permutations": int(permutations),
    }


def _one_way_anova_f(values: np.ndarray, labels: np.ndarray) -> Tuple[float, int, int]:
    unique = np.unique(labels)
    n = values.shape[0]
    k = len(unique)
    if k < 2:
        raise ValueError("Need at least two groups.")
    if n <= k:
        raise ValueError("Need more samples than groups for ANOVA denominator.")

    grand_mean = float(np.mean(values))
    ss_between = 0.0
    ss_within = 0.0
    for g in unique:
        group_vals = values[labels == g]
        group_mean = float(np.mean(group_vals))
        ss_between += len(group_vals) * (group_mean - grand_mean) ** 2
        ss_within += float(np.sum((group_vals - group_mean) ** 2))

    ms_between = ss_between / (k - 1)
    ms_within = ss_within / (n - k)
    f_stat = np.inf if ms_within == 0 else ms_between / ms_within
    return float(f_stat), n, k


def permdisp(coords: np.ndarray, labels: Sequence[str], permutations: int, rng: np.random.Generator) -> Dict[str, float]:
    labels_arr = np.asarray(labels)
    unique = np.unique(labels_arr)

    distances = np.zeros(coords.shape[0], dtype=float)
    mean_distance_by_group: Dict[str, float] = {}
    for g in unique:
        idx = labels_arr == g
        centroid = coords[idx].mean(axis=0)
        d = np.linalg.norm(coords[idx] - centroid, axis=1)
        distances[idx] = d
        mean_distance_by_group[str(g)] = float(np.mean(d))

    observed_f, n, k = _one_way_anova_f(distances, labels_arr)

    if permutations <= 0:
        p_value = 1.0
    else:
        count = 0
        for _ in range(permutations):
            perm_labels = rng.permutation(labels_arr)
            perm_f, _, _ = _one_way_anova_f(distances, perm_labels)
            if perm_f >= observed_f:
                count += 1
        p_value = (count + 1) / (permutations + 1)

    result: Dict[str, float] = {
        "f_stat": float(observed_f),
        "p_value": float(p_value),
        "n_samples": int(n),
        "n_groups": int(k),
        "permutations": int(permutations),
    }
    for g in unique:
        result[f"mean_distance_{str(g)}"] = mean_distance_by_group[str(g)]
    return result


def bh_fdr(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([], dtype=float)
    order = np.argsort(p)
    sorted_p = p[order]
    adj = sorted_p * n / np.arange(1, n + 1)
    for i in range(n - 2, -1, -1):
        adj[i] = min(adj[i], adj[i + 1])
    adj = np.minimum(adj, 1.0)
    out = np.zeros(n, dtype=float)
    out[order] = adj
    return out


def _group_palette(group_order: Sequence[str], group_labels: Dict[str, str]) -> Dict[str, str]:
    fixed = {
        "control": "#1f77b4",
        "recovery": "#2ca02c",
        "stemi": "#d62728",
    }
    cmap = plt.get_cmap("tab20")
    palette: Dict[str, str] = {}
    for i, g in enumerate(group_order):
        palette[g] = fixed.get(g, cmap(i % 20))
    return palette


def plot_sample_scatter(
    coords: np.ndarray,
    groups: Sequence[str],
    sample_ids: Sequence[str],
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    out_png: Path,
) -> None:
    groups_arr = np.asarray(groups)
    if coords.shape[1] == 1:
        xs = coords[:, 0]
        ys = np.zeros(coords.shape[0], dtype=float)
        xlabel = "Axis 1"
        ylabel = "Axis 2 (constant)"
    else:
        xs = coords[:, 0]
        ys = coords[:, 1]
        xlabel = "Axis 1"
        ylabel = "Axis 2"

    palette = _group_palette(group_order, group_labels)
    plt.figure(figsize=(8, 6))
    for g in group_order:
        idx = np.where(groups_arr == g)[0]
        if len(idx) == 0:
            continue
        plt.scatter(
            xs[idx],
            ys[idx],
            s=80,
            alpha=0.9,
            label=group_labels.get(g, g),
            color=palette[g],
            edgecolors="black",
            linewidths=0.4,
        )
        for i in idx:
            plt.text(xs[i], ys[i], sample_ids[i], fontsize=7, alpha=0.8)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title("Genome-wide promoter methylation sample embedding")
    plt.grid(alpha=0.2)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def plot_global_boxplot(
    sample_df: pd.DataFrame,
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    out_png: Path,
) -> None:
    palette = _group_palette(group_order, group_labels)
    order = [g for g in group_order if g in set(sample_df["group"])]

    plt.figure(figsize=(7, 5))
    positions = np.arange(len(order))
    data = [sample_df.loc[sample_df["group"] == g, "global_methylation_pct"].to_numpy() for g in order]
    labels = [group_labels.get(g, g) for g in order]

    bp = plt.boxplot(data, positions=positions, widths=0.6, patch_artist=True)
    for patch, g in zip(bp["boxes"], order):
        patch.set_facecolor(palette[g])
        patch.set_alpha(0.35)

    for i, g in enumerate(order):
        vals = sample_df.loc[sample_df["group"] == g, "global_methylation_pct"].to_numpy()
        jitter = np.random.default_rng(42 + i).normal(0, 0.04, size=len(vals))
        plt.scatter(np.full(len(vals), i) + jitter, vals, s=45, color=palette[g], edgecolors="black", linewidths=0.3)

    plt.xticks(positions, labels)
    plt.ylabel("Coverage-weighted global promoter methylation (%)")
    plt.title("Global methylation by cohort")
    plt.grid(axis="y", alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def build_readable_summary_table(
    permanova_df: pd.DataFrame,
    permdisp_df: pd.DataFrame,
    group_labels: Dict[str, str],
) -> pd.DataFrame:
    merged = permanova_df.merge(
        permdisp_df[["comparison", "test_scope", "f_stat", "p_value"]],
        on=["comparison", "test_scope"],
        how="left",
        suffixes=("_permanova", "_permdisp"),
    )

    rows: List[Dict[str, object]] = []
    for _, row in merged.iterrows():
        p_perm = float(row["p_value_permanova"])
        p_disp = float(row["p_value_permdisp"]) if pd.notna(row["p_value_permdisp"]) else np.nan
        q_perm = float(row["p_adj_bh"]) if pd.notna(row.get("p_adj_bh", np.nan)) else np.nan

        if p_perm < 0.05 and pd.notna(p_disp) and p_disp < 0.05:
            interpretation = "PERMANOVA sig.; dispersion also differs (interpret centroid shift cautiously)"
        elif p_perm < 0.05 and (pd.isna(p_disp) or p_disp >= 0.05):
            interpretation = "PERMANOVA sig.; no strong dispersion evidence"
        else:
            interpretation = "No PERMANOVA significance"

        rows.append(
            {
                "comparison": format_comparison_label(str(row["comparison"]), group_labels),
                "scope": "Global" if row["test_scope"] == "global" else "Pairwise",
                "n_samples": int(row["n_samples"]),
                "pseudo_f": float(row["pseudo_f"]),
                "r_squared": float(row["r_squared"]),
                "p_permanova": p_perm,
                "q_permanova_bh_pairwise_only": q_perm,
                "f_permdisp": float(row["f_stat"]) if pd.notna(row["f_stat"]) else np.nan,
                "p_permdisp": p_disp,
                "interpretation": interpretation,
            }
        )

    out = pd.DataFrame(rows)
    scope_order = {"Global": 0, "Pairwise": 1}
    out["_scope_order"] = out["scope"].map(scope_order).fillna(9)
    out = out.sort_values(["_scope_order", "comparison"]).drop(columns=["_scope_order"]).reset_index(drop=True)
    return out


def plot_readable_summary_table(table_df: pd.DataFrame, out_png: Path) -> None:
    display_df = table_df.copy()
    for col in [
        "pseudo_f",
        "r_squared",
        "p_permanova",
        "q_permanova_bh_pairwise_only",
        "f_permdisp",
        "p_permdisp",
    ]:
        if col in display_df.columns:
            display_df[col] = display_df[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.4g}")

    fig_h = max(2.5, 0.55 * (len(display_df) + 1))
    fig, ax = plt.subplots(figsize=(18, fig_h))
    ax.axis("off")

    col_widths = [0.16, 0.07, 0.06, 0.07, 0.07, 0.08, 0.11, 0.07, 0.07, 0.24]
    tbl = ax.table(
        cellText=display_df.values,
        colLabels=display_df.columns,
        colWidths=col_widths,
        cellLoc="center",
        loc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    tbl.scale(1, 1.35)

    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#e9ecef")
            cell.set_text_props(weight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#f8f9fa")

    ax.set_title("Genome-wide Methylation PERMANOVA/PERMDISP Summary", fontsize=13, fontweight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def write_readable_summary_html(table_df: pd.DataFrame, out_html: Path) -> None:
    html = [
        "<!doctype html>",
        "<html><head><meta charset='utf-8'><title>PERMANOVA/PERMDISP Summary</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;margin:24px;background:#f6f8fa;color:#222;}",
        "h2{margin-top:0;}",
        "table{border-collapse:collapse;width:100%;background:white;}",
        "th,td{border:1px solid #d0d7de;padding:8px 10px;font-size:13px;}",
        "th{background:#e9ecef;text-align:center;}",
        "tr:nth-child(even){background:#f8f9fa;}",
        "td:last-child{text-align:left;}",
        "</style></head><body>",
        "<h2>Genome-wide Methylation PERMANOVA/PERMDISP Summary</h2>",
        table_df.to_html(index=False, border=0, escape=False),
        "</body></html>",
    ]
    out_html.write_text("\n".join(html), encoding="utf-8")


def compute_group_weighted_promoter_means(
    meth_matrix: np.ndarray,
    total_matrix: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    min_cpgs: int,
) -> Dict[str, np.ndarray]:
    """Compute per-promoter coverage-weighted mean methylation per group."""
    group_means: Dict[str, np.ndarray] = {}
    groups_arr = np.asarray(groups)
    for group in group_order:
        idx = np.where(groups_arr == group)[0]
        n_promoters = meth_matrix.shape[1]
        means = np.full(n_promoters, np.nan, dtype=float)
        if len(idx) == 0:
            group_means[group] = means
            continue

        g_meth = meth_matrix[idx, :].astype(float)
        g_total = total_matrix[idx, :].astype(float)
        valid = g_total >= float(min_cpgs)

        weighted_meth = np.where(valid, g_meth, 0.0).sum(axis=0)
        weighted_total = np.where(valid, g_total, 0.0).sum(axis=0)
        np.divide(weighted_meth, weighted_total, out=means, where=weighted_total > 0)
        means *= 100.0
        group_means[group] = means

    return group_means


def plot_spearman_pair(
    mean_a: np.ndarray,
    mean_b: np.ndarray,
    label_a: str,
    label_b: str,
    out_png: Path,
    out_summary: Path,
) -> None:
    valid = np.isfinite(mean_a) & np.isfinite(mean_b)
    x = mean_a[valid]
    y = mean_b[valid]

    if x.size > 1:
        rho, pvalue = scipy_stats.spearmanr(x, y)
        rho = float(rho)
        pvalue = float(pvalue)
    else:
        rho = float("nan")
        pvalue = float("nan")

    plt.figure(figsize=(8, 8))
    plt.scatter(x, y, alpha=0.25, s=10, color="#1f77b4", edgecolors="none")

    if x.size > 0 and y.size > 0:
        min_val = float(np.nanmin([np.nanmin(x), np.nanmin(y)]))
        max_val = float(np.nanmax([np.nanmax(x), np.nanmax(y)]))
        plt.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="black", linewidth=1)

    plt.xlabel(f"{label_a} weighted mean promoter methylation (%)")
    plt.ylabel(f"{label_b} weighted mean promoter methylation (%)")
    plt.title(
        f"Promoter methylation concordance: {label_a} vs {label_b}\n"
        f"Spearman rho = {rho:.4f}, N = {int(valid.sum())}"
    )
    plt.text(
        0.05,
        0.95,
        f"Spearman rho = {rho:.4f}\nN = {int(valid.sum())}",
        transform=plt.gca().transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    with open(out_summary, "w", encoding="utf-8") as f:
        f.write(f"Spearman analysis: {label_a} vs {label_b}\n")
        f.write(f"N promoters: {int(valid.sum())}\n")
        f.write(f"Spearman rho: {rho:.8f}\n")
        f.write(f"Spearman p-value: {pvalue:.8e}\n")


def write_spearman_outputs(
    meth_matrix: np.ndarray,
    total_matrix: np.ndarray,
    groups: Sequence[str],
    group_order: Sequence[str],
    group_labels: Dict[str, str],
    min_cpgs: int,
    outdir: Path,
) -> None:
    group_means = compute_group_weighted_promoter_means(
        meth_matrix=meth_matrix,
        total_matrix=total_matrix,
        groups=groups,
        group_order=group_order,
        min_cpgs=min_cpgs,
    )

    for g1, g2 in itertools.combinations(group_order, 2):
        label_a = group_labels.get(g1, g1)
        label_b = group_labels.get(g2, g2)
        slug = f"{g1}_vs_{g2}"
        out_png = outdir / f"spearman_{slug}.png"
        out_txt = outdir / f"spearman_{slug}_summary.txt"
        plot_spearman_pair(group_means[g1], group_means[g2], label_a, label_b, out_png, out_txt)


def run_analysis(args: argparse.Namespace) -> None:
    np.random.seed(args.random_seed)
    rng = np.random.default_rng(args.random_seed)

    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found in PATH. Install bedtools first.")

    bed_dir = Path(args.bed_dir)
    if not bed_dir.exists():
        raise FileNotFoundError(f"Input BED directory not found: {bed_dir}")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    group_aliases = parse_key_value_arg(args.group_aliases, lower_values=True)
    if not group_aliases:
        group_aliases = DEFAULT_GROUP_ALIASES.copy()

    group_labels = parse_key_value_arg(args.group_labels, lower_values=False)
    if not group_labels:
        group_labels = DEFAULT_GROUP_LABELS.copy()

    samples = discover_samples(bed_dir, group_aliases=group_aliases)
    sample_df = pd.DataFrame([s.__dict__ for s in samples])
    sample_df = sample_df.sort_values(["group", "sample_id"]).reset_index(drop=True)

    # Keep order deterministic and match sample arrays.
    samples = [SampleInfo(**row) for row in sample_df.to_dict(orient="records")]

    group_counts = sample_df["group"].value_counts().to_dict()
    discovered_groups = sorted(group_counts.keys())
    requested_order = parse_group_order_arg(args.group_order) if args.group_order else []

    group_order: List[str] = []
    for g in requested_order if requested_order else DEFAULT_GROUP_ORDER:
        if g in discovered_groups and g not in group_order:
            group_order.append(g)
    for g in discovered_groups:
        if g not in group_order:
            group_order.append(g)

    for g in discovered_groups:
        if g not in group_labels:
            group_labels[g] = g.upper() if len(g) <= 5 else g.title()

    print("Discovered samples:")
    for g in group_order:
        print(f"  {group_labels.get(g, g):<12}: {group_counts.get(g, 0)}")

    if len(group_order) < 2:
        raise RuntimeError("Need at least 2 groups for between-group tests.")

    min_group_size = min(group_counts.get(g, 0) for g in group_order)
    if min_group_size < 2:
        raise RuntimeError("Need at least 2 samples per group for stable group-level tests.")

    min_samples = args.min_samples
    if min_samples > min_group_size:
        print(
            f"WARNING: min-samples={min_samples} exceeds smallest group size={min_group_size}. "
            f"Using {min_group_size}."
        )
        min_samples = min_group_size

    catalog = GencodePromoterCatalog(args.gencode_gtf).build_catalog()
    promoter_items = list(catalog.items())
    gene_count = len(promoter_items)
    print(f"Loaded promoters: {gene_count:,}")

    promoter_bed = write_promoter_bed(outdir, promoter_items)
    try:
        meth_matrix, total_matrix = aggregate_all_samples(
            samples=samples,
            promoter_bed=promoter_bed,
            gene_count=gene_count,
            max_workers=args.workers if args.workers is not None else max(1, os.cpu_count() or 1),
        )
    finally:
        try:
            promoter_bed.unlink(missing_ok=True)
        except Exception:
            pass

    groups = [s.group for s in samples]
    sample_ids = [s.sample_id for s in samples]

    keep = filter_promoters_by_group_coverage(total_matrix, groups, group_order, args.min_cpgs, min_samples)
    kept_genes = np.where(keep)[0]
    if len(kept_genes) < 20:
        raise RuntimeError(
            "Too few promoters left after filtering. Try lowering --min-cpgs or --min-samples."
        )

    meth_pct = np.full(meth_matrix.shape, np.nan, dtype=float)
    valid = total_matrix >= args.min_cpgs
    meth_pct[valid] = (meth_matrix[valid] / total_matrix[valid]) * 100.0
    meth_pct = meth_pct[:, keep]
    total_kept = total_matrix[:, keep]

    meth_imputed, na_fraction = median_impute_by_promoter(meth_pct)

    # Arcsin-sqrt transform for beta-like methylation proportions.
    eps = 1e-6
    beta = np.clip(meth_imputed / 100.0, eps, 1.0 - eps)
    transformed = np.arcsin(np.sqrt(beta))

    coords = to_euclidean_preserving_coords(transformed)

    # Global weighted methylation per sample (interpretability companion metric).
    total_all = total_kept.sum(axis=1).astype(float)
    meth_all = meth_matrix[:, keep].sum(axis=1).astype(float)
    global_meth = np.divide(meth_all, total_all, out=np.zeros_like(meth_all), where=total_all > 0) * 100.0

    sample_out = sample_df.copy()
    sample_out["global_methylation_pct"] = global_meth
    sample_out.to_csv(outdir / "sample_global_methylation.csv", index=False)

    embedding_df = pd.DataFrame({"sample_id": sample_ids, "group": groups, "axis1": coords[:, 0]})
    if coords.shape[1] > 1:
        embedding_df["axis2"] = coords[:, 1]
    else:
        embedding_df["axis2"] = 0.0
    embedding_df.to_csv(outdir / "sample_embedding.csv", index=False)

    # Global (all-group) tests
    global_perm = permanova(coords, groups, args.permutations, rng)
    global_disp = permdisp(coords, groups, args.permutations, rng)
    global_comparison_name = "_vs_".join(group_order)

    permanova_rows: List[Dict[str, object]] = [
        {
            "comparison": global_comparison_name,
            "test_scope": "global",
            **global_perm,
        }
    ]
    permdisp_rows: List[Dict[str, object]] = [
        {
            "comparison": global_comparison_name,
            "test_scope": "global",
            **global_disp,
        }
    ]

    # Pairwise tests
    pairwise = list(itertools.combinations(group_order, 2))
    for g1, g2 in pairwise:
        idx = [i for i, g in enumerate(groups) if g in (g1, g2)]
        sub_coords = coords[idx, :]
        sub_groups = [groups[i] for i in idx]

        p_res = permanova(sub_coords, sub_groups, args.permutations, rng)
        d_res = permdisp(sub_coords, sub_groups, args.permutations, rng)

        permanova_rows.append(
            {
                "comparison": f"{g1}_vs_{g2}",
                "test_scope": "pairwise",
                **p_res,
            }
        )
        permdisp_rows.append(
            {
                "comparison": f"{g1}_vs_{g2}",
                "test_scope": "pairwise",
                **d_res,
            }
        )

    permanova_df = pd.DataFrame(permanova_rows)
    permdisp_df = pd.DataFrame(permdisp_rows)

    # FDR only on the three pairwise PERMANOVA p-values
    pair_mask = permanova_df["test_scope"] == "pairwise"
    permanova_df.loc[pair_mask, "p_adj_bh"] = bh_fdr(permanova_df.loc[pair_mask, "p_value"].tolist())
    permanova_df.loc[~pair_mask, "p_adj_bh"] = np.nan

    permanova_df.to_csv(outdir / "genomewide_permanova_results.csv", index=False)
    permdisp_df.to_csv(outdir / "genomewide_permdisp_results.csv", index=False)

    readable_table = build_readable_summary_table(permanova_df, permdisp_df, group_labels)
    readable_table.to_csv(outdir / "genomewide_permanova_readable_table.csv", index=False)
    plot_readable_summary_table(readable_table, outdir / "genomewide_permanova_readable_table.png")
    write_readable_summary_html(readable_table, outdir / "genomewide_permanova_readable_table.html")

    # Global methylation univariate tests.
    group_arrays = {
        g: sample_out.loc[sample_out["group"] == g, "global_methylation_pct"].to_numpy(dtype=float)
        for g in group_order
    }
    kw_stat, kw_p = kruskal(*(group_arrays[g] for g in group_order))

    uw_rows: List[Dict[str, object]] = [
        {
            "test": f"kruskal_global_methylation_{len(group_order)}group",
            "statistic": float(kw_stat),
            "p_value": float(kw_p),
        }
    ]
    for g1, g2 in pairwise:
        stat, p = mannwhitneyu(group_arrays[g1], group_arrays[g2], alternative="two-sided")
        uw_rows.append({"test": f"mannwhitney_{g1}_vs_{g2}", "statistic": float(stat), "p_value": float(p)})

    uw_df = pd.DataFrame(uw_rows)
    mw_mask = uw_df["test"].str.startswith("mannwhitney_")
    uw_df.loc[mw_mask, "p_adj_bh"] = bh_fdr(uw_df.loc[mw_mask, "p_value"].tolist())
    uw_df.to_csv(outdir / "global_methylation_univariate_tests.csv", index=False)

    # Promoter coverage QC summary.
    qc_df = pd.DataFrame(
        {
            "n_promoters_total": [gene_count],
            "n_promoters_kept": [len(kept_genes)],
            "fraction_promoters_kept": [len(kept_genes) / gene_count],
            "mean_missing_fraction_before_impute": [float(np.mean(na_fraction))],
            "median_missing_fraction_before_impute": [float(np.median(na_fraction))],
        }
    )
    qc_df.to_csv(outdir / "promoter_matrix_qc.csv", index=False)

    plot_sample_scatter(
        coords,
        groups,
        sample_ids,
        group_order,
        group_labels,
        outdir / "sample_embedding.png",
    )
    plot_global_boxplot(sample_out, group_order, group_labels, outdir / "global_methylation_boxplot.png")
    write_spearman_outputs(
        meth_matrix=meth_matrix[:, keep],
        total_matrix=total_matrix[:, keep],
        groups=groups,
        group_order=group_order,
        group_labels=group_labels,
        min_cpgs=args.min_cpgs,
        outdir=outdir,
    )

    # Short report.
    lines = [
        "# Genome-wide promoter methylation PERMANOVA/PERMDISP",
        "",
        "## Inputs",
        f"- BED directory: {bed_dir}",
        f"- GENCODE GTF: {args.gencode_gtf}",
        f"- Group aliases: {group_aliases}",
        f"- Group order used: {group_order}",
        "",
        "## Matrix",
        f"- Samples: {len(samples)}",
        f"- Promoters total: {gene_count:,}",
        f"- Promoters kept: {len(kept_genes):,}",
        f"- Mean missing fraction before impute: {float(np.mean(na_fraction)):.4f}",
        "",
        "## Key outputs",
        "- genomewide_permanova_results.csv",
        "- genomewide_permdisp_results.csv",
        "- genomewide_permanova_readable_table.csv",
        "- genomewide_permanova_readable_table.png",
        "- genomewide_permanova_readable_table.html",
        "- global_methylation_univariate_tests.csv",
        "- sample_global_methylation.csv",
        "- sample_embedding.csv",
        "- sample_embedding.png",
        "- global_methylation_boxplot.png",
        "- spearman_<group1>_vs_<group2>.png",
        "- spearman_<group1>_vs_<group2>_summary.txt",
    ]
    (outdir / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")

    print("Analysis complete")
    print(f"Results directory: {outdir}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Genome-wide promoter methylation PERMANOVA/PERMDISP for coronary cfDNA"
    )
    parser.add_argument("--bed-dir", default=DEFAULT_BED_DIR, help="Directory with bedMethyl BED files")
    parser.add_argument("--gencode-gtf", default=DEFAULT_GENCODE_GTF, help="GENCODE GTF path (.gz supported)")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR, help="Output directory")
    parser.add_argument("--min-cpgs", type=int, default=DEFAULT_MIN_CPGS, help="Minimum CpGs required per promoter/sample")
    parser.add_argument(
        "--min-samples",
        type=int,
        default=DEFAULT_MIN_SAMPLES,
        help="Minimum observed samples per group for a promoter to be retained",
    )
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS, help="Permutation count")
    parser.add_argument("--workers", type=int, default=None, help="Max worker processes for bedtools aggregation")
    parser.add_argument("--random-seed", type=int, default=DEFAULT_RANDOM_SEED, help="Random seed")
    parser.add_argument(
        "--group-aliases",
        type=str,
        default="",
        help=(
            "Optional prefix:group mapping (comma-separated), e.g. "
            "'con:control,rec:recovery,cor:stemi'. Unmapped prefixes use themselves."
        ),
    )
    parser.add_argument(
        "--group-labels",
        type=str,
        default="",
        help=(
            "Optional group:display_label mapping (comma-separated), e.g. "
            "'control:Control,recovery:Recovery,stemi:STEMI'."
        ),
    )
    parser.add_argument(
        "--group-order",
        type=str,
        default=",".join(DEFAULT_GROUP_ORDER),
        help=(
            "Optional comma-separated group order for reporting/plotting, e.g. "
            "'control,recovery,stemi'. Discovered groups not listed are appended."
        ),
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
