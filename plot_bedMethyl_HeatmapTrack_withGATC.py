#!/usr/bin/env python3
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

REFERENCE_FASTA = "XXXXX/AdaptiveSequencing/GCF_030685395.1_ASM3068539v1_genomic.fna"

DATASET_FILES = [
    "NoAdaptive.filtered.bedMethyl",
    "Pmand_enriched.filtered.bedMethyl",
    "Stub_depleted.filtered.bedMethyl",
]

DATASET_LABELS = {
    "NoAdaptive.filtered.bedMethyl": "NoAdaptive",
    "Pmand_enriched.filtered.bedMethyl": "Pmand_enriched",
    "Stub_depleted.filtered.bedMethyl": "Stub_depleted",
}

GLOBAL_TICK_STEP = 500_000
BIN_SIZE = 10_000
SMOOTH_SIGMA_BINS = 2.0
GROUP_SIZE = 3


def read_fasta_sequences(fasta_path):
    sequences = {}
    current_name = None
    seq_chunks = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    sequences[current_name] = "".join(seq_chunks)

                current_name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.upper())
        if current_name:
            sequences[current_name] = "".join(seq_chunks)

    return sequences


def find_gatc_positions(sequence):
    positions = []
    i = sequence.find("GATC")
    while i != -1:
        positions.append(i)
        i = sequence.find("GATC", i + 1)

    return np.array(positions, dtype=float)


def find_bedmethyl_file(folder, base_name):
    candidates = [
        folder / base_name,
        folder / f"{base_name}.bedMethyl",
        folder / f"{base_name}.tsv",
        folder / f"{base_name}.txt",
    ]
    for c in candidates:
        if c.exists():
            return c
            
    return None


def read_bedmethyl_high_only(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    required = {"chrom", "start", "end", "classification"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["classification"] = (df["classification"].astype(str).str.lower().str.strip())
    df["pos"] = (df["start"] + df["end"]) / 2.0
    df = df[df["classification"] == "high"]
    df = df[df["pos"].notna()]

    return df


def gaussian_kernel1d(sigma_bins, truncate=4.0):
    radius = int(truncate * sigma_bins + 0.5)
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-(x**2) / (2 * sigma_bins**2))
    kernel = kernel / kernel.sum()
    return kernel


def smooth_counts(counts, sigma_bins):
    kernel = gaussian_kernel1d(sigma_bins)
    return np.convolve(counts, kernel, mode="same")


def compute_density(pos, chrom_len, bin_size, sigma_bins):
    bins = np.arange(0, chrom_len + bin_size, bin_size)
    if len(bins) < 2:
        bins = np.array([0, chrom_len], dtype=float)
    counts, _ = np.histogram(pos, bins=bins)
    smoothed = smooth_counts(counts.astype(float), sigma_bins)

    return smoothed, bins


def plot_tracks(
    input_folders,
    output_file,
    chrom_filter=None,
    dpi=600,
    bin_size=BIN_SIZE,
    sigma_bins=SMOOTH_SIGMA_BINS,
):
    sequences = read_fasta_sequences(REFERENCE_FASTA)
    chrom_names = list(sequences.keys())
    if chrom_filter:
        chrom_names = [c for c in chrom_names if c in chrom_filter]
    if not chrom_names:
        raise ValueError("No chromosomes selected after applying --chroms filter.")

    chrom_lengths = {c: len(sequences[c]) for c in chrom_names}
    width_ratios = [chrom_lengths[c]for c in chrom_names]
    panel_order = []
    panel_data = {}

    #read in methylation data
    for folder_str in input_folders:
        folder = Path(folder_str)
        if not folder.exists():
            raise FileNotFoundError(f"Input folder does not exist: {folder}")
        folder_label = folder.name

        for dataset_file in DATASET_FILES:
            path = find_bedmethyl_file(folder, dataset_file)
            if path is None:
                raise FileNotFoundError(f"Could not find file '{dataset_file}' in folder: {folder}")

            dataset_label = DATASET_LABELS[dataset_file]
            row_label = f"{folder_label} | {dataset_label}"
            panel_order.append((folder_label, dataset_label, row_label))
            df = read_bedmethyl_high_only(path)
            for chrom in chrom_names:
                chrom_len = chrom_lengths[chrom]
                sub = df[df["chrom"] == chrom]
                if sub.empty:
                    pos = np.array([], dtype=float)
                else:
                    pos = sub["pos"].values

                smoothed, bins = compute_density(pos,chrom_len,bin_size,sigma_bins)
                panel_data[(folder_label, dataset_label, chrom)] = (smoothed, bins)

    #get GATC density
    gatc_data = {}
    for chrom in chrom_names:
        seq = sequences[chrom]
        pos = find_gatc_positions(seq)
        smoothed, bins = compute_density( pos,len(seq),bin_size,sigma_bins)
        gatc_data[chrom] = (smoothed, bins)
        
    # normalize 
    gatc_norm_data = {}
    methyl_gatc_scaled_data = {}
    
    for chrom in chrom_names:
        gatc_smoothed, gatc_bins = gatc_data[chrom]
    
        gatc_max = float(gatc_smoothed.max()) if gatc_smoothed.size > 0 else 0.0
        if gatc_max == 0:
            gatc_max = 1.0
        gatc_norm = np.clip(gatc_smoothed / gatc_max, 0, 1)
        gatc_norm_data[chrom] = (gatc_norm, gatc_bins)
    
        for folder_label, dataset_label, row_label in panel_order:
            methyl_smoothed, methyl_bins = panel_data[(folder_label, dataset_label, chrom)]
            methyl_scaled = methyl_smoothed / gatc_max
            methyl_plot_norm = np.clip(methyl_scaled, 0, 1)
            methyl_gatc_scaled_data[(folder_label, dataset_label, chrom)] = (methyl_plot_norm, methyl_bins)
            
    # get source data frame 
    source_rows = []
    for chrom in chrom_names:
        gatc_norm, gatc_bins = gatc_norm_data[chrom]
        bin_starts = gatc_bins[:-1]
        bin_ends = gatc_bins[1:]
        bin_mids = (bin_starts + bin_ends) / 2.0
    
        for i in range(len(gatc_norm)):
            row = {
                "chrom": chrom,
                "start": bin_starts[i],
                "end": bin_ends[i],
                "mid": bin_mids[i],
                "GATC_density": gatc_norm[i],
            }
    
            for folder_label, dataset_label, row_label in panel_order:
                methyl_plot_norm, methyl_bins = methyl_gatc_scaled_data[(folder_label, dataset_label, chrom)]
                col_name = f"{folder_label} | {dataset_label}"
                row[col_name] = (methyl_plot_norm[i] if i < len(methyl_plot_norm) else np.nan)
            source_rows.append(row)

    source_df = pd.DataFrame(source_rows)
    source_df.to_csv("XXXXX/AdaptiveSequencing/PaperFigures_IK82/Figure3C_SourceData.csv",index=False)

    #plot
    n_rows = len(panel_order) + 1
    n_cols = len(chrom_names)
    fig_width = max(14,sum(width_ratios) / max(width_ratios) * 3)
    fig_height = max(2.5, n_rows * 0.4)

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(n_rows,n_cols,figure=fig,width_ratios=width_ratios,hspace=0.0,wspace=0.08)
    axes = [[None] * n_cols for _ in range(n_rows)    ]

    #Methylation rows
    for row_idx, (folder_label, dataset_label, row_label) in enumerate(panel_order):
        for col_idx, chrom in enumerate(chrom_names):
            sharex_ax = axes[0][col_idx] if row_idx > 0 else None
            ax = fig.add_subplot(
                gs[row_idx, col_idx],
                sharex=sharex_ax
            )
            axes[row_idx][col_idx] = ax
            norm, bins = methyl_gatc_scaled_data[(folder_label, dataset_label, chrom)]
            if norm.size > 0:
                x0 = bins[0]
                x1 = min(bins[-1], chrom_lengths[chrom]) - 1e-6

                ax.imshow(
                    norm[np.newaxis, :],
                    cmap="Greys",
                    aspect="auto",
                    interpolation="none",
                    extent=[x0, x1, 0, 1],
                    vmin=0,
                    vmax=1,
                    origin="lower",
                    zorder=2,
                )

            ax.set_xlim(0, chrom_lengths[chrom])
            ax.set_xbound(0, chrom_lengths[chrom])
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.margins(x=0, y=0)

            if row_idx == 0:
                ax.set_title(
                    chrom,
                    fontsize=8,
                    pad=2
                )

            if col_idx == 0:
                ax.set_ylabel(
                    row_label,
                    rotation=0,
                    ha="right",
                    va="center",
                    fontsize=8,
                    labelpad=16,
                )

            if row_idx == 0 or (
                row_idx % GROUP_SIZE == 0
                and row_idx != 0
            ):
                ax.axhline(
                    1,
                    color="black",
                    lw=2.0,
                    clip_on=False,
                    zorder=5,
                )

            ax.axhline(0,color="black",lw=0.35,alpha=0.25,clip_on=False,zorder=5)
            ax.tick_params(axis="x",bottom=False,labelbottom=False)
            ax.spines["bottom"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)

    #GATC density
    gatc_row = n_rows - 1
    for col_idx, chrom in enumerate(chrom_names):
        ax = fig.add_subplot(
            gs[gatc_row, col_idx],
            sharex=axes[0][col_idx]
        )

        axes[gatc_row][col_idx] = ax
        norm, bins = gatc_norm_data[chrom]
        if norm.size > 0:
            x0 = bins[0]
            x1 = min(bins[-1], chrom_lengths[chrom]) - 1e-6

            ax.imshow(
                norm[np.newaxis, :],
                cmap="Greys",
                aspect="auto",
                interpolation="none",
                extent=[x0, x1, 0, 1],
                vmin=0,
                vmax=1,
                origin="lower",
                zorder=2,
            )

        ax.set_xlim(0, chrom_lengths[chrom])
        ax.set_xbound(0, chrom_lengths[chrom])
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.margins(x=0, y=0)

        if col_idx == 0:
            ax.set_ylabel(
                "GATC density",
                rotation=0,
                ha="right",
                va="center",
                fontsize=8,
                labelpad=16,)

        ax.axhline(1,color="black",lw=2.0,clip_on=False,zorder=5)
        ax.axhline(0,color="black",lw=2.0,alpha=1.0,clip_on=False,zorder=6)
        ticks = np.arange(0,chrom_lengths[chrom],GLOBAL_TICK_STEP)
        ax.set_xticks(ticks)
        ax.ticklabel_format(style="sci",axis="x",scilimits=(0, 0))
        ax.xaxis.get_offset_text().set_size(7)
        ax.tick_params(axis="x",labelsize=9,pad=1,bottom=True, colors="black",)

        ax.spines["bottom"].set_color("black")
        ax.spines["bottom"].set_linewidth(2.0)
        ax.spines["bottom"].set_alpha(1.0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    #colorbar
    plt.subplots_adjust(
        left=0.12,
        right=0.90,
        top=0.95,
        bottom=0.10
    )
    
    cax = fig.add_axes([0.915, 0.16, 0.018, 0.72])
    norm_shared = Normalize(vmin=0, vmax=1)
    sm_shared = ScalarMappable(norm=norm_shared,cmap="Greys")
    sm_shared.set_array([])
    cbar = fig.colorbar(sm_shared,cax=cax)
    cbar.set_label("Density relative to max GATC density",fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_ticks([0, 1])

    fig.supxlabel("Position (bp)",fontsize=9,y=0.03)
    fig.savefig(output_file,dpi=dpi,bbox_inches="tight")
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_folders",nargs="+")
    parser.add_argument("-o","--output",default="heatmap_tracks.pdf")
    parser.add_argument("--chroms",nargs="+",default=None)
    parser.add_argument("--dpi",type=int,default=600)
    parser.add_argument( "--bin-size",type=int, default=BIN_SIZE)
    parser.add_argument("--sigma-bins", type=float, default=SMOOTH_SIGMA_BINS)
    
    args = parser.parse_args()

    plot_tracks(
        args.input_folders,
        args.output,
        chrom_filter=args.chroms,
        dpi=args.dpi,
        bin_size=args.bin_size,
        sigma_bins=args.sigma_bins,
    )

if __name__ == "__main__":
    main()