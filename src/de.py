"""
SkinAge Atlas: Differential Expression and Association Analysis Module
Models gene expression against chronological age (continuous and dichotomous young vs old)
in the Discovery cohort (GSE226189), adjusting for sex.
Performs intra-individual paired differential expression for Sun-Exposed vs Sun-Protected skin
in the Validation cohort (GSE38308).
Generates publication-quality Volcano and MA plots.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

FIG_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results', 'figures'))
TAB_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results', 'tables'))
DATA_META = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'metadata'))
DATA_PROC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'processed'))

for d in [FIG_DIR, TAB_DIR]:
    os.makedirs(d, exist_ok=True)

sns.set_theme(style='ticks', font_scale=1.1)
plt.rcParams['font.sans-serif'] = 'Helvetica, Arial, DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0


def run_gse226189_de():
    """
    Model expression against continuous age and young vs old groups, adjusted for sex.
    """
    print("\n--- Running Differential Expression on GSE226189 (Discovery Cohort) ---")
    meta = pd.read_csv(os.path.join(DATA_META, 'gse226189_metadata.csv'))
    log2_cpm_sym = pd.read_csv(os.path.join(DATA_PROC, 'gse226189_log2_cpm_symbols.csv'), index_col=0)
    
    # Align samples
    common_samples = [s for s in meta['sample_id'] if s in log2_cpm_sym.columns]
    meta = meta.set_index('sample_id').loc[common_samples].reset_index()
    expr = log2_cpm_sym[common_samples]

    # 1. Continuous Age association model: Expression ~ Age + Sex (Male=1, Female=0)
    print("Fitting linear model for continuous age association (adjusted for sex)...")
    sex_binary = (meta['sex'] == 'Male').astype(int).values
    age_vals = meta['age'].values

    # Design matrix: Intercept, Age, Sex
    X = np.column_stack([np.ones(len(common_samples)), age_vals, sex_binary])

    results = []
    for gene, y in zip(expr.index, expr.values):
        beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        dof = len(y) - X.shape[1]
        res = y - X @ beta
        sigma2 = np.sum(res**2) / dof
        cov = sigma2 * np.linalg.pinv(X.T @ X)
        se_age = np.sqrt(np.maximum(cov[1, 1], 1e-12))
        t_stat = beta[1] / se_age
        p_val = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=dof))

        results.append({
            'symbol': gene,
            'mean_expr': np.mean(y),
            'slope_per_year': beta[1],
            'log2FC_per_decade': beta[1] * 10,
            'sex_effect': beta[2],
            't_stat': t_stat,
            'p_value': p_val
        })

    de_df = pd.DataFrame(results)
    reject, q_vals, _, _ = multipletests(de_df['p_value'], method='fdr_bh')
    de_df['fdr_qval'] = q_vals
    de_df['neg_log10_padj'] = -np.log10(np.maximum(de_df['fdr_qval'], 1e-30))
    de_df['neg_log10_pval'] = -np.log10(np.maximum(de_df['p_value'], 1e-30))
    de_df['significant_raw'] = (de_df['p_value'] < 0.01) & (np.abs(de_df['log2FC_per_decade']) > 0.1)

    de_df = de_df.sort_values('p_value').reset_index(drop=True)
    out_csv = os.path.join(TAB_DIR, 'gse226189_de_continuous_age.csv')
    de_df.to_csv(out_csv, index=False)

    # 2. Dichotomous Young (<35y, n=19) vs Old (>60y, n=33) Model
    print("Fitting dichotomous Young (<35) vs Old (>60) model...")
    young_samples = meta[meta['age'] < 35]['sample_id'].values
    old_samples = meta[meta['age'] > 60]['sample_id'].values
    yo_samples = list(young_samples) + list(old_samples)
    yo_meta = meta.set_index('sample_id').loc[yo_samples].reset_index()
    yo_expr = expr[yo_samples]

    yo_group = (yo_meta['age'] > 60).astype(int).values  # 1 = Old, 0 = Young
    yo_sex = (yo_meta['sex'] == 'Male').astype(int).values
    X_yo = np.column_stack([np.ones(len(yo_samples)), yo_group, yo_sex])

    yo_results = []
    for gene, y in zip(yo_expr.index, yo_expr.values):
        beta, _, _, _ = np.linalg.lstsq(X_yo, y, rcond=None)
        dof = len(y) - X_yo.shape[1]
        res = y - X_yo @ beta
        sigma2 = np.sum(res**2) / dof
        cov = sigma2 * np.linalg.pinv(X_yo.T @ X_yo)
        se_grp = np.sqrt(np.maximum(cov[1, 1], 1e-12))
        t_stat = beta[1] / se_grp
        p_val = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=dof))

        yo_results.append({
            'symbol': gene,
            'mean_expr': np.mean(y),
            'log2FC_old_vs_young': beta[1],
            't_stat': t_stat,
            'p_value': p_val
        })

    yo_df = pd.DataFrame(yo_results)
    _, yo_q_vals, _, _ = multipletests(yo_df['p_value'], method='fdr_bh')
    yo_df['fdr_qval'] = yo_q_vals
    yo_df['neg_log10_padj'] = -np.log10(np.maximum(yo_df['fdr_qval'], 1e-30))
    yo_df['neg_log10_pval'] = -np.log10(np.maximum(yo_df['p_value'], 1e-30))
    yo_df['significant_raw'] = (yo_df['p_value'] < 0.01) & (np.abs(yo_df['log2FC_old_vs_young']) > 0.4)
    yo_df = yo_df.sort_values('p_value').reset_index(drop=True)
    yo_csv = os.path.join(TAB_DIR, 'gse226189_de_young_vs_old.csv')
    yo_df.to_csv(yo_csv, index=False)

    print(f"GSE226189 DE Results:")
    print(f"  Continuous age model: {de_df['significant_raw'].sum()} genes with p < 0.01 and |log2FC/decade| > 0.1")
    print(f"  Young vs Old model: {yo_df['significant_raw'].sum()} genes with p < 0.01 and |log2FC| > 0.4")

    # Volcano Plot: Continuous Age Model
    plt.figure(figsize=(8.5, 6.5))
    non_sig = de_df[~de_df['significant_raw']]
    sig_up = de_df[(de_df['significant_raw']) & (de_df['log2FC_per_decade'] > 0)]
    sig_down = de_df[(de_df['significant_raw']) & (de_df['log2FC_per_decade'] < 0)]

    plt.scatter(non_sig['log2FC_per_decade'], non_sig['neg_log10_pval'], c='#cccccc', s=15, alpha=0.5, label='Not Significant')
    plt.scatter(sig_up['log2FC_per_decade'], sig_up['neg_log10_pval'], c='#d95f02', s=25, alpha=0.85, label=f'Age-Positive (p<0.01, n={len(sig_up)})')
    plt.scatter(sig_down['log2FC_per_decade'], sig_down['neg_log10_pval'], c='#1f78b4', s=25, alpha=0.85, label=f'Age-Negative (p<0.01, n={len(sig_down)})')

    plt.axhline(-np.log10(0.01), color='#666666', linestyle='--', linewidth=0.9, alpha=0.7, label='p = 0.01')
    plt.axvline(-0.1, color='#666666', linestyle=':', linewidth=0.9, alpha=0.7)
    plt.axvline(0.1, color='#666666', linestyle=':', linewidth=0.9, alpha=0.7)

    # Label top 8 positive and top 8 negative
    top_label = pd.concat([sig_up.head(7), sig_down.head(7)])
    for _, row in top_label.iterrows():
        plt.annotate(
            row['symbol'],
            (row['log2FC_per_decade'], row['neg_log10_pval']),
            xytext=(5, 4), textcoords='offset points',
            fontsize=8.5, fontweight='bold', alpha=0.95
        )

    plt.title('GSE226189: Volcano Plot of Chronological Age Associations\n(Linear Regression Adjusted for Sex, N=82)', fontweight='bold')
    plt.xlabel('log2 Fold Change per Decade of Age')
    plt.ylabel('-log10(p-value)')
    plt.legend(frameon=True, loc='upper left')
    plt.tight_layout()
    fig_volc = os.path.join(FIG_DIR, 'fig5_gse226189_age_volcano.png')
    plt.savefig(fig_volc, dpi=300)
    plt.close()
    print(f"Saved: {fig_volc}")

    # MA Plot (log2FC vs Mean Expression)
    plt.figure(figsize=(8, 6))
    plt.scatter(non_sig['mean_expr'], non_sig['log2FC_per_decade'], c='#cccccc', s=12, alpha=0.4, label='Not Significant')
    plt.scatter(sig_up['mean_expr'], sig_up['log2FC_per_decade'], c='#d95f02', s=22, alpha=0.8, label='Age-Positive')
    plt.scatter(sig_down['mean_expr'], sig_down['log2FC_per_decade'], c='#1f78b4', s=22, alpha=0.8, label='Age-Negative')
    plt.axhline(0, color='black', linestyle='--', linewidth=0.9)
    plt.title('GSE226189: MA Plot (Effect Size vs Mean Abundance)', fontweight='bold')
    plt.xlabel('Mean Expression (log2 CPM)')
    plt.ylabel('log2 Fold Change per Decade')
    plt.legend(frameon=True)
    plt.tight_layout()
    fig_ma = os.path.join(FIG_DIR, 'fig6_gse226189_age_ma_plot.png')
    plt.savefig(fig_ma, dpi=300)
    plt.close()
    print(f"Saved: {fig_ma}")

    return de_df, yo_df


def run_gse38308_de():
    """
    Intra-individual paired differential expression for Sun-Exposed vs Sun-Protected skin.
    21 biological donors, each contributing 1 exposed (A) and 1 protected (B) sample.
    """
    print("\n--- Running Paired DE on GSE38308 (Sun-Exposed vs Sun-Protected) ---")
    meta = pd.read_csv(os.path.join(DATA_META, 'gse38308_metadata.csv'))
    expr = pd.read_csv(os.path.join(DATA_PROC, 'gse38308_expression_symbols.csv'), index_col=0)

    # Collect pairs (1A vs 1B, 2A vs 2B, ..., 21A vs 21B)
    donors = sorted(list(set(meta['donor_id'])))
    pairs_exposed = []
    pairs_protected = []
    donor_list = []

    for d in donors:
        sub = meta[meta['donor_id'] == d]
        exp_sample = sub[sub['sun_exposure'] == 'Sun-Exposed']['sample_id'].values
        prot_sample = sub[sub['sun_exposure'] == 'Sun-Protected']['sample_id'].values
        if len(exp_sample) == 1 and len(prot_sample) == 1:
            if exp_sample[0] in expr.columns and prot_sample[0] in expr.columns:
                pairs_exposed.append(exp_sample[0])
                pairs_protected.append(prot_sample[0])
                donor_list.append(d)

    print(f"Found {len(donor_list)} complete intra-individual matched pairs.")

    mat_exp = expr[pairs_exposed].values
    mat_prot = expr[pairs_protected].values
    diff_mat = mat_exp - mat_prot  # delta = exposed - protected

    # Paired t-test per gene
    mean_diff = np.mean(diff_mat, axis=1)
    t_stats, p_vals = stats.ttest_rel(mat_exp, mat_prot, axis=1)

    reject, q_vals, _, _ = multipletests(p_vals, method='fdr_bh')

    de_val = pd.DataFrame({
        'symbol': expr.index,
        'mean_expr': np.mean((mat_exp + mat_prot) / 2, axis=1),
        'mean_diff_log2FC': mean_diff,
        't_stat': t_stats,
        'p_value': p_vals,
        'fdr_qval': q_vals,
        'neg_log10_padj': -np.log10(np.maximum(q_vals, 1e-30)),
        'neg_log10_pval': -np.log10(np.maximum(p_vals, 1e-30)),
    })
    de_val['significant'] = (de_val['fdr_qval'] < 0.05) & (np.abs(de_val['mean_diff_log2FC']) > 0.15)
    de_val = de_val.sort_values('p_value').reset_index(drop=True)

    out_csv = os.path.join(TAB_DIR, 'gse38308_de_sun_exposure.csv')
    de_val.to_csv(out_csv, index=False)

    num_sig = de_val['significant'].sum()
    num_up = ((de_val['significant']) & (de_val['mean_diff_log2FC'] > 0)).sum()
    num_down = ((de_val['significant']) & (de_val['mean_diff_log2FC'] < 0)).sum()
    print(f"Sun-Exposure DE summary (Paired t-test, FDR < 0.05, |log2FC| > 0.15):")
    print(f"  Total tested genes: {len(de_val):,}")
    print(f"  Significant photoaging genes: {num_sig:,} ({num_up} higher in sun-exposed, {num_down} higher in sun-protected)")

    # Volcano Plot for Sun Exposure
    plt.figure(figsize=(8.5, 6.5))
    non_sig = de_val[~de_val['significant']]
    sig_up = de_val[(de_val['significant']) & (de_val['mean_diff_log2FC'] > 0)]
    sig_down = de_val[(de_val['significant']) & (de_val['mean_diff_log2FC'] < 0)]

    plt.scatter(non_sig['mean_diff_log2FC'], non_sig['neg_log10_padj'], c='#cccccc', s=15, alpha=0.5, label='Not Significant')
    plt.scatter(sig_up['mean_diff_log2FC'], sig_up['neg_log10_padj'], c='#e7298a', s=25, alpha=0.85, label=f'Sun-Exposed Up (FDR<0.05, n={len(sig_up)})')
    plt.scatter(sig_down['mean_diff_log2FC'], sig_down['neg_log10_padj'], c='#66a61e', s=25, alpha=0.85, label=f'Sun-Protected Up (FDR<0.05, n={len(sig_down)})')

    plt.axhline(-np.log10(0.05), color='#555555', linestyle='--', linewidth=0.9, alpha=0.7, label='FDR = 0.05')
    plt.axvline(-0.15, color='#555555', linestyle=':', linewidth=0.9, alpha=0.7)
    plt.axvline(0.15, color='#555555', linestyle=':', linewidth=0.9, alpha=0.7)

    # Annotate top genes
    top_label = pd.concat([sig_up.head(6), sig_down.head(6)])
    for _, row in top_label.iterrows():
        plt.annotate(
            row['symbol'],
            (row['mean_diff_log2FC'], row['neg_log10_padj']),
            xytext=(5, 4), textcoords='offset points',
            fontsize=8.5, fontweight='bold', alpha=0.95
        )

    plt.title('GSE38308: Volcano Plot (Sun-Exposed vs Sun-Protected Skin)\n(Intra-Individual Matched Pairs, N=21 Donors)', fontweight='bold')
    plt.xlabel('Mean Difference in Expression [log2(Exposed) - log2(Protected)]')
    plt.ylabel('-log10(FDR Adjusted P-value)')
    plt.legend(frameon=True, loc='upper right')
    plt.tight_layout()
    fig_volc = os.path.join(FIG_DIR, 'fig7_gse38308_sun_exposure_volcano.png')
    plt.savefig(fig_volc, dpi=300)
    plt.close()
    print(f"Saved: {fig_volc}")

    return de_val


if __name__ == '__main__':
    de_disc, de_yo = run_gse226189_de()
    de_val = run_gse38308_de()
    print("\n=== Differential Expression step completed successfully! ===")
