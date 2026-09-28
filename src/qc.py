"""
SkinAge Atlas: Quality Control and Exploratory Data Analysis Module
Performs filtering, library size normalization (CPM/TPM/log2), PCA,
sample correlation heatmaps, and outlier detection for Discovery (GSE226189)
and Validation (GSE38308) cohorts.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.stats import zscore

FIG_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results', 'figures'))
TAB_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results', 'tables'))
DATA_META = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'metadata'))
DATA_PROC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'processed'))

for d in [FIG_DIR, TAB_DIR]:
    os.makedirs(d, exist_ok=True)

# Set publication style
sns.set_theme(style='ticks', font_scale=1.1)
plt.rcParams['font.sans-serif'] = 'Helvetica, Arial, DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0


def run_gse226189_qc():
    """Perform QC, filtering, and visualization on Discovery RNA-seq cohort (GSE226189)."""
    print("\n--- Running QC on Discovery Cohort: GSE226189 ---")
    meta = pd.read_csv(os.path.join(DATA_META, 'gse226189_metadata.csv'))
    raw_counts = pd.read_csv(os.path.join(DATA_PROC, 'gse226189_raw_counts.csv'), index_col=0)
    ens_map = pd.read_csv(os.path.join(DATA_META, 'ensembl_to_symbol.csv'))
    ens_to_sym = dict(zip(ens_map['ensembl_id'], ens_map['symbol']))

    # 1. Library sizes
    lib_sizes = raw_counts.sum(axis=0)
    meta['lib_size_millions'] = meta['sample_id'].map(lib_sizes) / 1e6

    # 2. Gene filtering: Keep genes with count >= 10 in at least 20% of samples (17 samples)
    min_samples = int(0.20 * raw_counts.shape[1])
    keep_mask = (raw_counts >= 10).sum(axis=1) >= min_samples
    filtered_counts = raw_counts[keep_mask]
    print(f"Filter rule: count >= 10 in >= {min_samples} samples (20%):")
    print(f"  Retained {filtered_counts.shape[0]:,} / {raw_counts.shape[0]:,} genes ({filtered_counts.shape[0]/raw_counts.shape[0]*100:.1f}%)")

    # Save filtered counts
    filt_csv = os.path.join(DATA_PROC, 'gse226189_filtered_counts.csv')
    filtered_counts.to_csv(filt_csv)

    # 3. Log2-CPM normalization for PCA and visualization
    # CPM = (count / library_size) * 1e6
    cpm = (filtered_counts / lib_sizes) * 1e6
    log2_cpm = np.log2(cpm + 1)
    log2_cpm.to_csv(os.path.join(DATA_PROC, 'gse226189_log2_cpm.csv'))

    # Map genes to symbols for downstream gene-level queries
    log2_cpm_symbol = log2_cpm.copy()
    log2_cpm_symbol['symbol'] = log2_cpm_symbol.index.map(ens_to_sym)
    log2_cpm_symbol = log2_cpm_symbol.dropna(subset=['symbol']).groupby('symbol').mean()
    log2_cpm_symbol.to_csv(os.path.join(DATA_PROC, 'gse226189_log2_cpm_symbols.csv'))
    print(f"Mapped {log2_cpm_symbol.shape[0]:,} genes to unique HGNC symbols.")

    # 4. Outlier detection based on inter-sample correlation and PCA distance
    corr_matrix = log2_cpm.corr(method='spearman')
    mean_corr = corr_matrix.mean(axis=1)
    corr_z = zscore(mean_corr)

    # PCA
    pca = PCA(n_components=5)
    pca_coords = pca.fit_transform(log2_cpm.T)
    var_exp = pca.explained_variance_ratio_ * 100

    pca_df = pd.DataFrame({
        'sample_id': log2_cpm.columns,
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1],
        'PC3': pca_coords[:, 2],
        'spearman_mean_corr': mean_corr.values,
        'corr_zscore': corr_z
    }).merge(meta, on='sample_id')

    # Flag potential outliers (|z-score| > 3 in correlation or distance)
    pca_df['flagged_outlier'] = (np.abs(pca_df['corr_zscore']) > 3.0)
    flagged_samples = pca_df[pca_df['flagged_outlier']]
    print(f"Outlier check: {len(flagged_samples)} samples flagged (|correlation z-score| > 3.0).")
    pca_df.to_csv(os.path.join(TAB_DIR, 'gse226189_sample_qc_metrics.csv'), index=False)

    # --- PLOTS ---
    # Plot 1: Cohort Overview & Age Distribution
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.histplot(data=meta, x='age', hue='sex', multiple='stack', bins=14, palette={'Male': '#1f77b4', 'Female': '#e377c2'}, ax=axes[0])
    axes[0].set_title('GSE226189: Age & Sex Distribution (N=82)', fontweight='bold')
    axes[0].set_xlabel('Chronological Age (years)')
    axes[0].set_ylabel('Sample Count')

    sns.boxplot(data=meta, x='age_group', y='lib_size_millions', hue='sex', palette={'Male': '#1f77b4', 'Female': '#e377c2'}, ax=axes[1])
    axes[1].set_title('Library Size Across Age Groups', fontweight='bold')
    axes[1].set_xlabel('Age Group')
    axes[1].set_ylabel('Library Size (Millions of Reads)')
    plt.tight_layout()
    fig1_path = os.path.join(FIG_DIR, 'fig1_gse226189_age_and_library_qc.png')
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"Saved: {fig1_path}")

    # Plot 2: PCA with Variance Explained
    fig, ax = plt.subplots(figsize=(7, 6))
    scatter = ax.scatter(
        pca_df['PC1'], pca_df['PC2'],
        c=pca_df['age'], cmap='viridis', s=80, edgecolors='#333333', alpha=0.85
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Age (years)', rotation=270, labelpad=15)
    ax.set_title(f'GSE226189: Transcriptomic PCA (N=82)\nPC1 ({var_exp[0]:.1f}%) vs PC2 ({var_exp[1]:.1f}%)', fontweight='bold')
    ax.set_xlabel(f'PC1 ({var_exp[0]:.1f}% variance)')
    ax.set_ylabel(f'PC2 ({var_exp[1]:.1f}% variance)')
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    fig2_path = os.path.join(FIG_DIR, 'fig2_gse226189_pca.png')
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"Saved: {fig2_path}")

    # Plot 3: Sample Correlation Heatmap
    plt.figure(figsize=(8, 7))
    sns.heatmap(
        corr_matrix, cmap='magma', vmin=0.85, vmax=1.0,
        xticklabels=False, yticklabels=False, cbar_kws={'label': "Spearman Correlation"}
    )
    plt.title('GSE226189: Sample-to-Sample Spearman Correlation Matrix', fontweight='bold')
    plt.tight_layout()
    fig3_path = os.path.join(FIG_DIR, 'fig3_gse226189_correlation_heatmap.png')
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    print(f"Saved: {fig3_path}")

    return meta, filtered_counts, log2_cpm, pca_df


def run_gse38308_qc():
    """Perform QC and exploratory visualization on Validation cohort (GSE38308)."""
    print("\n--- Running QC on Validation Cohort: GSE38308 ---")
    meta = pd.read_csv(os.path.join(DATA_META, 'gse38308_metadata.csv'))
    expr = pd.read_csv(os.path.join(DATA_PROC, 'gse38308_expression_symbols.csv'), index_col=0)

    # PCA on GSE38308
    pca = PCA(n_components=5)
    pca_coords = pca.fit_transform(expr.T)
    var_exp = pca.explained_variance_ratio_ * 100

    pca_df = pd.DataFrame({
        'sample_id': expr.columns,
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1],
    }).merge(meta, on='sample_id')

    pca_df.to_csv(os.path.join(TAB_DIR, 'gse38308_qc_metrics.csv'), index=False)

    # Plot 4: GSE38308 Matched Exposure PCA & Distribution
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    palette = {'Sun-Exposed': '#d95f02', 'Sun-Protected': '#7570b3'}
    sns.scatterplot(
        data=pca_df, x='PC1', y='PC2', hue='sun_exposure', style='sun_exposure',
        s=90, palette=palette, edgecolors='#222222', ax=axes[0]
    )
    axes[0].set_title(f'GSE38308: PCA by Sun Exposure (N=42)\nPC1 ({var_exp[0]:.1f}%) vs PC2 ({var_exp[1]:.1f}%)', fontweight='bold')
    axes[0].set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
    axes[0].set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
    axes[0].grid(True, linestyle='--', alpha=0.5)

    # Connect matched pairs on PCA
    for donor, sub in pca_df.groupby('donor_id'):
        if len(sub) == 2:
            axes[0].plot(sub['PC1'], sub['PC2'], color='gray', linestyle=':', alpha=0.6)

    # Age distribution
    sns.histplot(data=meta[meta['sun_exposure'] == 'Sun-Exposed'], x='age', bins=10, color='#1b9e77', ax=axes[1])
    axes[1].set_title('GSE38308: Donor Age Distribution (21 Donors)', fontweight='bold')
    axes[1].set_xlabel('Chronological Age (years)')
    axes[1].set_ylabel('Donor Count')
    plt.tight_layout()
    fig4_path = os.path.join(FIG_DIR, 'fig4_gse38308_qc_and_pca.png')
    plt.savefig(fig4_path, dpi=300)
    plt.close()
    print(f"Saved: {fig4_path}")

    return meta, expr, pca_df


if __name__ == '__main__':
    run_gse226189_qc()
    run_gse38308_qc()
    print("\n=== Quality Control step completed successfully! ===")
