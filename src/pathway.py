"""
SkinAge Atlas: Pathway and Functional Enrichment Module
Performs Gene Set Enrichment Analysis (GSEA) on Hallmark, Reactome, and GO Biological Process gene sets
using GSEA PreRank on chronological aging (GSE226189) and photoaging (GSE38308).
Generates publication-quality enrichment bubble / bar plots.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp

FIG_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results', 'figures'))
TAB_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'results', 'tables'))

for d in [FIG_DIR, TAB_DIR]:
    os.makedirs(d, exist_ok=True)

sns.set_theme(style='ticks', font_scale=1.1)
plt.rcParams['font.sans-serif'] = 'Helvetica, Arial, DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 1.0


def run_prerank_gsea(rnk_series: pd.Series, gene_set: str, out_prefix: str) -> pd.DataFrame:
    """Run GSEA prerank on a ranked series of genes and save results."""
    print(f"Running GSEA PreRank for {out_prefix} against {gene_set} ({len(rnk_series):,} genes)...")
    res = gp.prerank(
        rnk=rnk_series,
        gene_sets=gene_set,
        min_size=15,
        max_size=500,
        permutation_num=1000,
        seed=42,
        verbose=False
    )
    df = res.res2d.copy()
    df = df.sort_values('NES', ascending=False).reset_index(drop=True)
    out_csv = os.path.join(TAB_DIR, f"{out_prefix}_{gene_set}.csv")
    df.to_csv(out_csv, index=False)
    print(f"  Saved {len(df)} pathway results to {out_csv}")
    return df


def plot_gsea_enrichment(df: pd.DataFrame, title: str, filename: str, top_n: int = 8):
    """Generate dual bar/lollipop plot of top positive and negative enriched pathways."""
    # Top positive (up in age/exposure) and top negative (down in age/exposure)
    pos = df[df['NES'] > 0].sort_values('NES', ascending=False).head(top_n)
    neg = df[df['NES'] < 0].sort_values('NES', ascending=True).head(top_n)
    plot_df = pd.concat([pos, neg]).copy()
    plot_df = plot_df.sort_values('NES', ascending=True)

    # Clean term names (remove library prefixes or suffixes if present)
    plot_df['Clean_Term'] = plot_df['Term'].apply(lambda x: x.split(' (GO:')[0] if ' (GO:' in x else x)

    plt.figure(figsize=(10, 7.5))
    colors = ['#1f78b4' if nes < 0 else '#d95f02' for nes in plot_df['NES']]

    bars = plt.barh(plot_df['Clean_Term'], plot_df['NES'], color=colors, edgecolor='#333333', height=0.65)
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')

    # Annotate FDR
    for bar, qval in zip(bars, plot_df['FDR q-val']):
        x = bar.get_width()
        offset = 0.05 if x > 0 else -0.05
        align = 'left' if x > 0 else 'right'
        fdr_str = f"FDR={qval:.3f}" if qval >= 0.001 else "FDR<0.001"
        plt.text(x + offset, bar.get_y() + bar.get_height() / 2, fdr_str,
                 va='center', ha=align, fontsize=8.5, color='#444444')

    plt.title(title, fontweight='bold', pad=15)
    plt.xlabel('Normalized Enrichment Score (NES)')
    plt.ylabel('')
    xlim = max(abs(plot_df['NES'].min()), abs(plot_df['NES'].max())) * 1.35
    plt.xlim(-xlim, xlim)
    plt.grid(axis='x', linestyle=':', alpha=0.6)
    plt.tight_layout()

    out_path = os.path.join(FIG_DIR, filename)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def run_pathway_analysis():
    """Run pathway analysis for both Discovery (age) and Validation (sun exposure)."""
    # 1. Discovery Cohort: GSE226189 Age-ranked genes
    de_disc = pd.read_csv(os.path.join(TAB_DIR, 'gse226189_de_continuous_age.csv'))
    rnk_age = de_disc[['symbol', 't_stat']].dropna().drop_duplicates(subset=['symbol']).set_index('symbol')['t_stat']
    rnk_age = rnk_age.sort_values(ascending=False)

    hallmark_age = run_prerank_gsea(rnk_age, 'MSigDB_Hallmark_2020', 'gse226189_age')
    plot_gsea_enrichment(
        hallmark_age,
        'GSE226189: Hallmark Pathways Associated with Chronological Ageing\n(Ranked by Linear Model t-statistic)',
        'fig8_gse226189_hallmark_gsea.png'
    )

    reactome_age = run_prerank_gsea(rnk_age, 'Reactome_2022', 'gse226189_age')
    plot_gsea_enrichment(
        reactome_age,
        'GSE226189: Top Reactome Pathways Associated with Chronological Ageing',
        'fig9_gse226189_reactome_gsea.png'
    )

    # 2. Validation Cohort: GSE38308 Sun-exposure ranked genes
    de_val = pd.read_csv(os.path.join(TAB_DIR, 'gse38308_de_sun_exposure.csv'))
    rnk_sun = de_val[['symbol', 't_stat']].dropna().drop_duplicates(subset=['symbol']).set_index('symbol')['t_stat']
    rnk_sun = rnk_sun.sort_values(ascending=False)

    hallmark_sun = run_prerank_gsea(rnk_sun, 'MSigDB_Hallmark_2020', 'gse38308_sun_exposure')
    plot_gsea_enrichment(
        hallmark_sun,
        'GSE38308: Hallmark Pathways Enriched in Sun-Exposed vs Sun-Protected Skin\n(Ranked by Paired t-statistic)',
        'fig10_gse38308_hallmark_gsea.png'
    )

    # Compile plain-language biological explanations
    explanations = [
        {
            'Pathway': 'Protein Secretion',
            'Cohort': 'GSE226189 (Chronological Age)',
            'NES': hallmark_age.loc[hallmark_age['Term'] == 'Protein Secretion', 'NES'].values[0],
            'FDR_qval': hallmark_age.loc[hallmark_age['Term'] == 'Protein Secretion', 'FDR q-val'].values[0],
            'Biological_Direction': 'Enriched in Older Skin',
            'Plain_Language_Summary': (
                "Aged fibroblasts exhibit hyperactive protein secretion, reflecting endoplasmic reticulum stress "
                "and the Senescence-Associated Secretory Phenotype (SASP) where senescent cells secrete inflammatory cytokines, "
                "chemokines, and matrix metalloproteinases."
            ),
            'Literature_Lookup_Recommendation': (
                "Search PubMed for: 'fibroblast senescence SASP secretome skin aging' (e.g., Coppé et al., 2008; Campisi et al.)."
            )
        },
        {
            'Pathway': 'Epithelial Mesenchymal Transition (EMT)',
            'Cohort': 'GSE226189 (Chronological Age)',
            'NES': hallmark_age.loc[hallmark_age['Term'] == 'Epithelial Mesenchymal Transition', 'NES'].values[0],
            'FDR_qval': hallmark_age.loc[hallmark_age['Term'] == 'Epithelial Mesenchymal Transition', 'FDR q-val'].values[0],
            'Biological_Direction': 'Enriched in Older Skin',
            'Plain_Language_Summary': (
                "Includes genes involved in extracellular matrix remodeling, collagen degradation, and cellular plasticity. "
                "As dermal fibroblasts age, they shift into a more myofibroblastic, fibrotic state with altered matrix synthesis."
            ),
            'Literature_Lookup_Recommendation': (
                "Search PubMed for: 'dermal fibroblast aging extracellular matrix collagen EMT plasticity'."
            )
        },
        {
            'Pathway': 'Cholesterol Homeostasis',
            'Cohort': 'GSE226189 (Chronological Age)',
            'NES': hallmark_age.loc[hallmark_age['Term'] == 'Cholesterol Homeostasis', 'NES'].values[0],
            'FDR_qval': hallmark_age.loc[hallmark_age['Term'] == 'Cholesterol Homeostasis', 'FDR q-val'].values[0],
            'Biological_Direction': 'Downregulated in Older Skin',
            'Plain_Language_Summary': (
                "Lipid and sterol synthesis pathways decline significantly in aged dermal and epidermal cells. "
                "This correlates with the clinically documented breakdown of skin barrier function, reduced sebum production, and increased dryness."
            ),
            'Literature_Lookup_Recommendation': (
                "Search PubMed for: 'skin lipid metabolism aging cholesterol barrier function' (e.g., Zouboulis et al.)."
            )
        },
        {
            'Pathway': 'Wnt-beta Catenin Signaling',
            'Cohort': 'GSE226189 (Chronological Age)',
            'NES': hallmark_age.loc[hallmark_age['Term'] == 'Wnt-beta Catenin Signaling', 'NES'].values[0],
            'FDR_qval': hallmark_age.loc[hallmark_age['Term'] == 'Wnt-beta Catenin Signaling', 'FDR q-val'].values[0],
            'Biological_Direction': 'Downregulated in Older Skin',
            'Plain_Language_Summary': (
                "Wnt signaling maintains skin cell identity and regenerative priming. Its loss with age leads to reduced fibroblast density "
                "and compromised dermal replenishment."
            ),
            'Literature_Lookup_Recommendation': (
                "Search PubMed for: 'dermal fibroblast loss of priming Wnt signaling skin aging' (e.g., Solé-Boldo et al., Cell Reports 2020)."
            )
        }
    ]
    expl_df = pd.DataFrame(explanations)
    expl_df.to_csv(os.path.join(TAB_DIR, 'pathway_plain_language_explanations.csv'), index=False)
    print("Saved plain-language pathway explanations to results/tables/pathway_plain_language_explanations.csv")
    return hallmark_age, reactome_age, hallmark_sun


if __name__ == '__main__':
    run_pathway_analysis()
    print("\n=== Pathway analysis step completed successfully! ===")
