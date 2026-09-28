"""
SkinAge Atlas: Transcriptomic Age-Prediction Machine Learning Module
Trains an Elastic Net regression model on the Discovery cohort (GSE226189)
using strictly internal 5-fold cross-validation.
Evaluates the locked model on the completely independent Validation cohort (GSE38308).
Evaluates performance against a trivial baseline (predicting mean training age).
Extracts top weighted biomarker genes and analyzes platform transferability.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import ElasticNetCV, ElasticNet
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import mean_absolute_error, r2_score
from scipy import stats

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


def run_age_predictor():
    print("\n--- Training Transcriptomic Age-Prediction Model (Elastic Net) ---")

    # Load data
    meta_train = pd.read_csv(os.path.join(DATA_META, 'gse226189_metadata.csv'))
    expr_train = pd.read_csv(os.path.join(DATA_PROC, 'gse226189_log2_cpm_symbols.csv'), index_col=0)

    meta_test = pd.read_csv(os.path.join(DATA_META, 'gse38308_metadata.csv'))
    expr_test = pd.read_csv(os.path.join(DATA_PROC, 'gse38308_expression_symbols.csv'), index_col=0)

    # Find overlapping gene symbols
    common_genes = sorted(list(set(expr_train.index).intersection(set(expr_test.index))))
    print(f"Intersection: {len(common_genes):,} common gene symbols between RNA-seq discovery and Microarray validation.")

    # Prepare Training Matrices (GSE226189, n=82)
    sample_order_train = [s for s in meta_train['sample_id'] if s in expr_train.columns]
    meta_train = meta_train.set_index('sample_id').loc[sample_order_train].reset_index()
    X_train_raw = expr_train.loc[common_genes, sample_order_train].T.values  # (82, 8601)
    y_train = meta_train['age'].values

    # Prepare Test Matrices (GSE38308, n=42)
    sample_order_test = [s for s in meta_test['sample_id'] if s in expr_test.columns]
    meta_test = meta_test.set_index('sample_id').loc[sample_order_test].reset_index()
    X_test_raw = expr_test.loc[common_genes, sample_order_test].T.values  # (42, 8601)
    y_test = meta_test['age'].values

    # Step 1: Feature pre-filtering in training set to avoid noise:
    # Select top 500 genes with highest absolute Pearson correlation with age in training set
    corrs = [stats.pearsonr(X_train_raw[:, i], y_train)[0] for i in range(len(common_genes))]
    top_indices = np.argsort(np.abs(corrs))[-500:]
    selected_genes = [common_genes[i] for i in top_indices]
    print(f"Pre-selected top {len(selected_genes)} age-correlated genes in Discovery cohort.")

    X_train = X_train_raw[:, top_indices]
    X_test = X_test_raw[:, top_indices]

    # Standardize features using training data parameters ONLY (strictly no leakage)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # 5-fold Cross-Validation inside training cohort
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    enet_cv = ElasticNetCV(
        l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95],
        alphas=np.logspace(-2, 1, 20),
        cv=cv,
        max_iter=3000,
        random_state=42
    )
    enet_cv.fit(X_train_scaled, y_train)

    best_alpha = enet_cv.alpha_
    best_l1 = enet_cv.l1_ratio_
    print(f"Optimal Elastic Net parameters: alpha = {best_alpha:.4f}, l1_ratio = {best_l1:.2f}")

    # Generate cross-validated predictions on Discovery cohort
    y_train_pred_cv = np.zeros(len(y_train))
    for train_idx, val_idx in cv.split(X_train_scaled, y_train):
        fold_scaler = StandardScaler()
        X_tr = fold_scaler.fit_transform(X_train[train_idx])
        X_va = fold_scaler.transform(X_train[val_idx])
        model = ElasticNet(alpha=best_alpha, l1_ratio=best_l1, max_iter=3000, random_state=42)
        model.fit(X_tr, y_train[train_idx])
        y_train_pred_cv[val_idx] = model.predict(X_va)

    # Full training fit
    final_model = ElasticNet(alpha=best_alpha, l1_ratio=best_l1, max_iter=3000, random_state=42)
    final_model.fit(X_train_scaled, y_train)

    # Predict on completely independent test cohort (GSE38308)
    y_test_pred = final_model.predict(X_test_scaled)

    # Evaluate Metrics
    train_mae = mean_absolute_error(y_train, y_train_pred_cv)
    train_r, train_p = stats.pearsonr(y_train, y_train_pred_cv)
    train_rho, train_rhop = stats.spearmanr(y_train, y_train_pred_cv)
    train_r2 = r2_score(y_train, y_train_pred_cv)

    test_mae = mean_absolute_error(y_test, y_test_pred)
    test_r, test_p = stats.pearsonr(y_test, y_test_pred)
    test_rho, test_rhop = stats.spearmanr(y_test, y_test_pred)
    test_r2 = r2_score(y_test, y_test_pred)

    # Baseline: Trivial Predictor (Predict training mean age)
    mean_train_age = np.mean(y_train)
    baseline_train_pred = np.full_like(y_train, mean_train_age)
    baseline_test_pred = np.full_like(y_test, mean_train_age)
    baseline_train_mae = mean_absolute_error(y_train, baseline_train_pred)
    baseline_test_mae = mean_absolute_error(y_test, baseline_test_pred)

    print("\n--- MODEL PERFORMANCE SCORECARD ---")
    print(f"Discovery Cohort (GSE226189, N={len(y_train)}, 5-fold CV):")
    print(f"  MAE:          {train_mae:.2f} years (Baseline MAE: {baseline_train_mae:.2f} years)")
    print(f"  Pearson r:    {train_r:.3f} (p = {train_p:.2e})")
    print(f"  Spearman rho: {train_rho:.3f} (p = {train_rhop:.2e})")
    print(f"  R^2:          {train_r2:.3f}")

    print(f"\nIndependent Validation Cohort (GSE38308, N={len(y_test)}, Zero Refitting):")
    print(f"  MAE:          {test_mae:.2f} years (Baseline MAE: {baseline_test_mae:.2f} years)")
    print(f"  Pearson r:    {test_r:.3f} (p = {test_p:.2e})")
    print(f"  Spearman rho: {test_rho:.3f} (p = {test_rhop:.2e})")
    print(f"  R^2:          {test_r2:.3f}")

    # Save metrics table
    metrics_df = pd.DataFrame([
        {
            'Cohort': 'GSE226189_Discovery_CV',
            'Sample_Count': len(y_train),
            'Age_Range': f"{min(y_train)}-{max(y_train)}",
            'Model_MAE': train_mae,
            'Baseline_MAE': baseline_train_mae,
            'Pearson_r': train_r,
            'Pearson_pval': train_p,
            'Spearman_rho': train_rho,
            'R2': train_r2
        },
        {
            'Cohort': 'GSE38308_Independent_Validation',
            'Sample_Count': len(y_test),
            'Age_Range': f"{int(min(y_test))}-{int(max(y_test))}",
            'Model_MAE': test_mae,
            'Baseline_MAE': baseline_test_mae,
            'Pearson_r': test_r,
            'Pearson_pval': test_p,
            'Spearman_rho': test_rho,
            'R2': test_r2
        }
    ])
    metrics_csv = os.path.join(TAB_DIR, 'model_performance_metrics.csv')
    metrics_df.to_csv(metrics_csv, index=False)
    print(f"Saved performance metrics table to {metrics_csv}")

    # Save predictions
    train_pred_df = meta_train[['sample_id', 'age', 'sex']].copy()
    train_pred_df['predicted_age_cv'] = y_train_pred_cv
    train_pred_df.to_csv(os.path.join(TAB_DIR, 'gse226189_age_predictions.csv'), index=False)

    test_pred_df = meta_test[['sample_id', 'donor_id', 'age', 'sex', 'sun_exposure']].copy()
    test_pred_df['predicted_age'] = y_test_pred
    test_pred_df['age_delta'] = test_pred_df['predicted_age'] - test_pred_df['age']
    test_pred_df.to_csv(os.path.join(TAB_DIR, 'gse38308_age_predictions.csv'), index=False)

    # Extract non-zero gene coefficients
    coefs = final_model.coef_
    coef_df = pd.DataFrame({
        'symbol': selected_genes,
        'coefficient': coefs
    })
    coef_df = coef_df[coef_df['coefficient'] != 0].sort_values('coefficient', ascending=False)
    coef_csv = os.path.join(TAB_DIR, 'elastic_net_gene_coefficients.csv')
    coef_df.to_csv(coef_csv, index=False)
    print(f"Model retained {len(coef_df)} non-zero coefficient genes out of {len(selected_genes)} features.")

    # --- PLOTS ---
    # Fig 11: Predicted vs Actual Plot (Dual Panel: Discovery CV vs Independent Test)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    # Panel A: Discovery CV
    axes[0].scatter(y_train, y_train_pred_cv, c='#1f77b4', s=60, alpha=0.8, edgecolors='#222222')
    lims_a = [min(min(y_train), min(y_train_pred_cv)) - 5, max(max(y_train), max(y_train_pred_cv)) + 5]
    axes[0].plot(lims_a, lims_a, 'k--', alpha=0.6, label='Identity line (y = x)')
    # Trend line
    sns.regplot(x=y_train, y=y_train_pred_cv, ax=axes[0], scatter=False, color='#d95f02')
    axes[0].set_title(
        f'Discovery Cohort (GSE226189, 5-Fold CV)\nMAE = {train_mae:.1f}y | r = {train_r:.2f} (p < 0.001)',
        fontweight='bold'
    )
    axes[0].set_xlabel('Actual Chronological Age (years)')
    axes[0].set_ylabel('Predicted Transcriptomic Age (years)')
    axes[0].set_xlim(15, 95)
    axes[0].set_ylim(15, 95)
    axes[0].grid(True, linestyle=':', alpha=0.6)
    axes[0].legend(loc='upper left')

    # Panel B: Independent Validation (GSE38308)
    palette_exp = {'Sun-Exposed': '#d95f02', 'Sun-Protected': '#7570b3'}
    sns.scatterplot(
        data=test_pred_df, x='age', y='predicted_age', hue='sun_exposure',
        palette=palette_exp, s=75, edgecolors='#222222', ax=axes[1]
    )
    lims_b = [min(min(y_test), min(y_test_pred)) - 5, max(max(y_test), max(y_test_pred)) + 5]
    axes[1].plot(lims_b, lims_b, 'k--', alpha=0.6, label='Identity line')
    sns.regplot(x=y_test, y=y_test_pred, ax=axes[1], scatter=False, color='#2ca02c')
    axes[1].set_title(
        f'Independent Validation (GSE38308, Zero Refitting)\nMAE = {test_mae:.1f}y | r = {test_r:.2f} | Baseline MAE = {baseline_test_mae:.1f}y',
        fontweight='bold'
    )
    axes[1].set_xlabel('Actual Chronological Age (years)')
    axes[1].set_ylabel('Predicted Transcriptomic Age (years)')
    axes[1].grid(True, linestyle=':', alpha=0.6)
    axes[1].legend(loc='upper left')

    plt.tight_layout()
    fig11_path = os.path.join(FIG_DIR, 'fig11_predicted_vs_actual_age.png')
    plt.savefig(fig11_path, dpi=300)
    plt.close()
    print(f"Saved: {fig11_path}")

    # Fig 12: Top Biomarker Genes Bar Plot
    top_pos = coef_df.head(10)
    top_neg = coef_df.tail(10)
    plot_coef = pd.concat([top_pos, top_neg]).sort_values('coefficient', ascending=True)

    plt.figure(figsize=(9, 7))
    bar_colors = ['#1f78b4' if c < 0 else '#d95f02' for c in plot_coef['coefficient']]
    plt.barh(plot_coef['symbol'], plot_coef['coefficient'], color=bar_colors, edgecolor='#333333', height=0.65)
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    plt.title('Elastic Net Top Transcriptomic Biomarkers of Skin Ageing\n(Model Weights)', fontweight='bold', pad=15)
    plt.xlabel('Elastic Net Regression Coefficient')
    plt.ylabel('')
    plt.grid(axis='x', linestyle=':', alpha=0.6)
    plt.tight_layout()

    fig12_path = os.path.join(FIG_DIR, 'fig12_elastic_net_gene_weights.png')
    plt.savefig(fig12_path, dpi=300)
    plt.close()
    print(f"Saved: {fig12_path}")

    return metrics_df, coef_df


if __name__ == '__main__':
    run_age_predictor()
    print("\n=== Age-Prediction Model step completed successfully! ===")
