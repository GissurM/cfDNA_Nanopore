#!/usr/bin/env python3
"""Generic PCA analysis for motif tables (*.motif.csv).

Expected columns in each input file: motif, count, frequency.
Groups are inferred from filename prefixes or explicit mappings.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import LabelEncoder, StandardScaler


def parse_mapping(raw: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not raw:
        return out
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid mapping '{item}'. Use key:value format.")
        key, value = item.split(":", 1)
        out[key.strip().lower()] = value.strip()
    return out


def infer_group(sample_name: str, prefix_map: Dict[str, str]) -> str:
    lowered = sample_name.lower()
    for prefix, group in prefix_map.items():
        if lowered.startswith(prefix):
            return group
    return "unknown"


def draw_group_circle(ax, x_values: np.ndarray, y_values: np.ndarray, color: str) -> None:
    from matplotlib.patches import Circle

    if x_values.size == 0:
        return
    center_x = float(np.mean(x_values))
    center_y = float(np.mean(y_values))
    distances = np.sqrt((x_values - center_x) ** 2 + (y_values - center_y) ** 2)
    base_radius = float(np.max(distances)) if distances.size > 0 else 0.0
    padding = max(0.12, base_radius * 0.15)
    circle = Circle(
        (center_x, center_y),
        base_radius + padding,
        fill=False,
        edgecolor=color,
        linewidth=1.8,
        linestyle="--",
        alpha=0.85,
    )
    ax.add_patch(circle)
    ax.autoscale_view()


def make_palette(groups: List[str]) -> Dict[str, str]:
    cmap = plt.get_cmap("tab20")
    return {g: cmap(i % 20) for i, g in enumerate(groups)}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="PCA analysis for motif frequency tables")
    parser.add_argument("--input-dir", default="/home/gissu/UXM_deconv/results2", help="Directory with *.motif.csv files")
    parser.add_argument("--output-dir", default="/mnt/c/Users/gissu/Documents/pca_analysis", help="Output directory")
    parser.add_argument(
        "--group-prefixes",
        default="brca2:BRCA2,ctrl:ctrl",
        help="Comma-separated prefix:group mapping used to infer group from sample filename",
    )
    parser.add_argument(
        "--include-unknown",
        action="store_true",
        help="Include samples that do not match any --group-prefixes mapping",
    )
    parser.add_argument("--top-motifs", type=int, default=25, help="Number of most abundant motifs to include")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    print("=== Motif PCA Analysis ===")
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    prefix_map = parse_mapping(args.group_prefixes)
    if not prefix_map:
        raise ValueError("At least one --group-prefixes mapping is required.")

    csv_files = sorted(input_dir.glob("*.motif.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No .motif.csv files found in {input_dir}")

    all_data: List[pd.DataFrame] = []
    for file_path in csv_files:
        sample_name = file_path.name.replace(".motif.csv", "")
        group = infer_group(sample_name, prefix_map)
        if group == "unknown" and not args.include_unknown:
            continue

        df = pd.read_csv(file_path)
        cols = {c.lower().strip(): c for c in df.columns}
        required = {"motif", "count", "frequency"}
        missing = required.difference(cols.keys())
        if missing:
            raise ValueError(f"{file_path.name} missing required columns: {sorted(missing)}")

        frame = df[[cols["motif"], cols["count"], cols["frequency"]]].copy()
        frame.columns = ["motif", "count", "frequency"]
        frame["motif"] = frame["motif"].astype(str).str.upper()
        frame["sample"] = sample_name
        frame["group"] = group
        all_data.append(frame)

    if not all_data:
        raise RuntimeError("No samples available after group filtering.")

    combined_data = pd.concat(all_data, ignore_index=True)
    group_counts = combined_data[["sample", "group"]].drop_duplicates()["group"].value_counts().to_dict()
    if len(group_counts) < 2:
        raise RuntimeError("At least two groups are required for comparative PCA plotting.")

    print("Discovered groups:")
    for group, count in sorted(group_counts.items()):
        print(f"  {group}: {count}")

    print(f"\n=== Using Top {args.top_motifs} Most Common Motifs ===")
    motif_counts = combined_data.groupby("motif")["count"].sum().sort_values(ascending=False)
    top_motifs = motif_counts.head(args.top_motifs).index.tolist()

    pca_data = combined_data[combined_data["motif"].isin(top_motifs)].copy()
    pca_matrix = pca_data.pivot(index="sample", columns="motif", values="frequency").fillna(0.0)
    sample_info = pca_data[["sample", "group"]].drop_duplicates().set_index("sample").reindex(pca_matrix.index)

    print(f"PCA matrix shape: {pca_matrix.shape} (samples x motifs)")
    pca_scaled = StandardScaler().fit_transform(pca_matrix.to_numpy(dtype=float))

    pca = PCA()
    pca_result = pca.fit_transform(pca_scaled)
    explained_var = pca.explained_variance_ratio_
    cumulative_var = np.cumsum(explained_var)

    n_components = pca_result.shape[1]
    if n_components < 2:
        raise RuntimeError("Need at least 2 principal components. Add more samples or motifs.")

    print(f"PC1 explains {explained_var[0]:.1%} of variance")
    print(f"PC2 explains {explained_var[1]:.1%} of variance")
    print(f"PC1+PC2 explain {cumulative_var[1]:.1%} of variance")

    groups = sorted(sample_info["group"].unique().tolist())
    colors = make_palette(groups)

    fig = plt.figure(figsize=(20, 12))

    plt.subplot(2, 4, 1)
    plot_k = min(10, len(explained_var))
    plt.plot(range(1, plot_k + 1), explained_var[:plot_k], "bo-", linewidth=2, markersize=8)
    plt.xlabel("Principal Component", fontsize=12, fontweight="bold")
    plt.ylabel("Explained Variance Ratio", fontsize=12, fontweight="bold")
    plt.title("Scree Plot", fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3)

    plt.subplot(2, 4, 2)
    plt.plot(range(1, plot_k + 1), cumulative_var[:plot_k], "ro-", linewidth=2, markersize=8)
    plt.axhline(y=0.8, color="gray", linestyle="--", alpha=0.7, label="80% variance")
    plt.xlabel("Principal Component", fontsize=12, fontweight="bold")
    plt.ylabel("Cumulative Variance Ratio", fontsize=12, fontweight="bold")
    plt.title("Cumulative Variance Explained", fontsize=14, fontweight="bold")
    plt.legend()
    plt.grid(True, alpha=0.3)

    def scatter_panel(subplot_idx: int, pc_x: int, pc_y: int, title: str) -> None:
        plt.subplot(2, 4, subplot_idx)
        for group in groups:
            mask = sample_info["group"].to_numpy() == group
            plt.scatter(
                pca_result[mask, pc_x],
                pca_result[mask, pc_y],
                c=[colors[group]],
                label=f"{group} (n={int(mask.sum())})",
                s=100,
                alpha=0.7,
            )
        ax = plt.gca()
        for group in groups:
            mask = sample_info["group"].to_numpy() == group
            if int(mask.sum()) > 0:
                draw_group_circle(ax, pca_result[mask, pc_x], pca_result[mask, pc_y], colors[group])

        plt.xlabel(f"PC{pc_x + 1} ({explained_var[pc_x]:.1%} variance)", fontsize=12, fontweight="bold")
        plt.ylabel(f"PC{pc_y + 1} ({explained_var[pc_y]:.1%} variance)", fontsize=12, fontweight="bold")
        plt.title(title, fontsize=14, fontweight="bold")
        plt.legend(fontsize=8)
        plt.grid(True, alpha=0.3)

    scatter_panel(3, 0, 1, "PCA: PC1 vs PC2")
    if n_components >= 3:
        scatter_panel(4, 0, 2, "PCA: PC1 vs PC3")
        scatter_panel(5, 1, 2, "PCA: PC2 vs PC3")
    else:
        plt.subplot(2, 4, 4)
        plt.text(0.5, 0.5, "PC3 unavailable", ha="center", va="center")
        plt.axis("off")
        plt.subplot(2, 4, 5)
        plt.text(0.5, 0.5, "PC3 unavailable", ha="center", va="center")
        plt.axis("off")

    motif_names = pca_matrix.columns.to_numpy()
    plt.subplot(2, 4, 6)
    loadings_pc1 = pca.components_[0]
    idx1 = np.argsort(np.abs(loadings_pc1))[-10:]
    plt.barh(range(len(idx1)), loadings_pc1[idx1], color="steelblue")
    plt.yticks(range(len(idx1)), [motif_names[i] for i in idx1])
    plt.xlabel("PC1 Loading", fontsize=12, fontweight="bold")
    plt.title("Top 10 Motifs for PC1", fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3)

    plt.subplot(2, 4, 7)
    loadings_pc2 = pca.components_[1]
    idx2 = np.argsort(np.abs(loadings_pc2))[-10:]
    plt.barh(range(len(idx2)), loadings_pc2[idx2], color="darkgreen")
    plt.yticks(range(len(idx2)), [motif_names[i] for i in idx2])
    plt.xlabel("PC2 Loading", fontsize=12, fontweight="bold")
    plt.title("Top 10 Motifs for PC2", fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3)

    plt.subplot(2, 4, 8)
    for group in groups:
        mask = sample_info["group"].to_numpy() == group
        plt.scatter(pca_result[mask, 0], pca_result[mask, 1], c=[colors[group]], label=group, s=80, alpha=0.6)
    ax = plt.gca()
    for group in groups:
        mask = sample_info["group"].to_numpy() == group
        if int(mask.sum()) > 0:
            draw_group_circle(ax, pca_result[mask, 0], pca_result[mask, 1], colors[group])

    top_idx = np.argsort(np.abs(loadings_pc1) + np.abs(loadings_pc2))[-8:]
    scale_factor = 3.0
    for i in top_idx:
        plt.arrow(
            0,
            0,
            loadings_pc1[i] * scale_factor,
            loadings_pc2[i] * scale_factor,
            head_width=0.1,
            head_length=0.1,
            fc="black",
            ec="black",
            alpha=0.7,
        )
        plt.text(
            loadings_pc1[i] * scale_factor * 1.1,
            loadings_pc2[i] * scale_factor * 1.1,
            motif_names[i],
            fontsize=9,
            ha="center",
            va="center",
        )

    plt.xlabel(f"PC1 ({explained_var[0]:.1%} variance)", fontsize=12, fontweight="bold")
    plt.ylabel(f"PC2 ({explained_var[1]:.1%} variance)", fontsize=12, fontweight="bold")
    plt.title("PCA Biplot", fontsize=14, fontweight="bold")
    plt.legend(fontsize=8)
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "4mer_PCA_comprehensive_analysis.png", dpi=300, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=(12, 8))
    for group in groups:
        mask = sample_info["group"].to_numpy() == group
        plt.scatter(
            pca_result[mask, 0],
            pca_result[mask, 1],
            c=[colors[group]],
            label=f"{group} (n={int(mask.sum())})",
            s=150,
            alpha=0.8,
            edgecolors="black",
            linewidth=1,
        )
    ax = plt.gca()
    for group in groups:
        mask = sample_info["group"].to_numpy() == group
        if int(mask.sum()) > 0:
            draw_group_circle(ax, pca_result[mask, 0], pca_result[mask, 1], colors[group])

    for i, sample in enumerate(sample_info.index.tolist()):
        plt.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.7)

    plt.xlabel(f"PC1 ({explained_var[0]:.1%} variance explained)", fontsize=14, fontweight="bold")
    plt.ylabel(f"PC2 ({explained_var[1]:.1%} variance explained)", fontsize=14, fontweight="bold")
    plt.title("PCA of 4-mer Motif Frequencies", fontsize=16, fontweight="bold")
    plt.legend(fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "4mer_PCA_detailed_scatter.png", dpi=300, bbox_inches="tight")
    plt.close()

    print("\n=== PCA Analysis Results ===")
    if n_components >= 3:
        print(f"Total variance explained by first 3 PCs: {cumulative_var[2]:.1%}")
    else:
        print(f"Total variance explained by available PCs: {cumulative_var[-1]:.1%}")

    sil_score_msg = "NA"
    if len(groups) >= 2 and pca_result.shape[0] > len(groups):
        le = LabelEncoder()
        group_labels = le.fit_transform(sample_info["group"].to_numpy())
        sil_score = silhouette_score(pca_result[:, :2], group_labels)
        sil_score_msg = f"{sil_score:.3f}"
    print(f"Silhouette score (PC1+PC2): {sil_score_msg}")

    n_save = min(5, pca_result.shape[1])
    pca_results_df = pd.DataFrame(
        pca_result[:, :n_save],
        columns=[f"PC{i + 1}" for i in range(n_save)],
        index=pca_matrix.index,
    )
    pca_results_df["group"] = sample_info["group"]
    pca_results_df.to_csv(output_dir / "4mer_PCA_scores.csv")

    loadings_df = pd.DataFrame(
        pca.components_[:n_save].T,
        columns=[f"PC{i + 1}" for i in range(n_save)],
        index=pca_matrix.columns,
    )
    loadings_df.to_csv(output_dir / "4mer_PCA_loadings.csv")

    print("\n=== Top Contributing Motifs ===")
    top_pc_count = min(3, pca.components_.shape[0])
    for pc_num in range(top_pc_count):
        loadings = pca.components_[pc_num]
        top_idx = np.argsort(np.abs(loadings))[-5:][::-1]
        print(f"\nPC{pc_num + 1} top 5 contributors:")
        for i, idx in enumerate(top_idx):
            print(f"  {i + 1}. {motif_names[idx]}: {loadings[idx]:.3f}")

    sample_info.reset_index().to_csv(output_dir / "4mer_sample_group_mapping.csv", index=False)

    print("\n=== Files Saved ===")
    print("  - 4mer_PCA_comprehensive_analysis.png")
    print("  - 4mer_PCA_detailed_scatter.png")
    print("  - 4mer_PCA_scores.csv")
    print("  - 4mer_PCA_loadings.csv")
    print("  - 4mer_sample_group_mapping.csv")
    print(f"  - All files saved in: {output_dir}")
    print("\n=== Analysis Complete ===")


if __name__ == "__main__":
    main()
