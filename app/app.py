"""
SkinAge Atlas: Interactive Bioinformatics Dashboard
Cross-cohort transcriptomic analysis of human skin ageing,
photoaging, pathway dynamics, and machine learning biomarker validation.
"""

import os
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats

# ----------------------------------------------------
# STREAMLIT PAGE CONFIGURATION
# ----------------------------------------------------
st.set_page_config(
    page_title="SkinAge Atlas | Transcriptomic Analytics",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Base directories
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DATA_META = os.path.join(BASE_DIR, 'data', 'metadata')
DATA_PROC = os.path.join(BASE_DIR, 'data', 'processed')
TAB_DIR = os.path.join(BASE_DIR, 'results', 'tables')
FIG_DIR = os.path.join(BASE_DIR, 'results', 'figures')

# ----------------------------------------------------
# CUSTOM CSS STYLING
# ----------------------------------------------------
st.markdown("""
<style>
    .main-title {
        font-size: 2.2rem;
        font-weight: 800;
        color: #0f172a;
        margin-bottom: 0.2rem;
        letter-spacing: -0.5px;
    }
    .sub-title {
        font-size: 1.05rem;
        color: #475569;
        margin-bottom: 1.2rem;
    }
    .disclaimer-banner {
        background: linear-gradient(90deg, #fffbeb 0%, #fef3c7 100%);
        border-left: 5px solid #f59e0b;
        padding: 0.8rem 1.1rem;
        border-radius: 6px;
        color: #92400e;
        font-weight: 500;
        font-size: 0.90rem;
        margin-bottom: 1.5rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }
    .metric-box {
        background: #ffffff;
        border: 1px solid #e2e8f0;
        border-radius: 10px;
        padding: 1.1rem;
        text-align: center;
        box-shadow: 0 1px 3px rgba(0,0,0,0.03);
    }
    .metric-value {
        font-size: 1.8rem;
        font-weight: 700;
        color: #1e3a8a;
    }
    .metric-label {
        font-size: 0.85rem;
        color: #64748b;
        text-transform: uppercase;
        font-weight: 600;
    }
    .metric-sub {
        font-size: 0.80rem;
        color: #059669;
        margin-top: 0.2rem;
    }
</style>
""", unsafe_allow_html=True)

# ----------------------------------------------------
# HEADER & DISCLAIMER
# ----------------------------------------------------
st.markdown('<div class="main-title">🧬 SkinAge Atlas Dashboard</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Cross-Cohort Transcriptomic Analysis of Human Skin Ageing, Photoaging & Machine Learning Biomarker Validation</div>', unsafe_allow_html=True)

st.markdown("""
<div class="disclaimer-banner">
    ⚠️ <strong>Research Prototype Notice:</strong> All reported correlations are transcriptomic and statistical associations, not causal claims or clinical diagnoses.
    Designed for computational biology research, exploratory analysis, and target discovery.
</div>
""", unsafe_allow_html=True)

# ----------------------------------------------------
# SIDEBAR
# ----------------------------------------------------
st.sidebar.markdown("### 🧬 SkinAge Atlas Analytics")
st.sidebar.markdown("""
**Cohorts Analysed:**
* **Discovery:** [GSE226189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226189) ($n=82$, RNA-seq)
* **Validation:** [GSE38308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38308) ($n=42$, Microarray)
""")

st.sidebar.markdown("---")
st.sidebar.markdown("**Quick Stats:**")
st.sidebar.write("• **Total Donors:** 103 individuals")
st.sidebar.write("• **Genes Tested:** 20,798 symbols")
st.sidebar.write("• **Photoaging Genes:** 1,269 (FDR < 0.05)")
st.sidebar.write("• **Model CV Accuracy:** MAE = 10.67y")
st.sidebar.markdown("---")
st.sidebar.caption("© 2026 SkinAge Atlas | Built with Streamlit & Plotly")

# ----------------------------------------------------
# TABS
# ----------------------------------------------------
tab_overview, tab_audit, tab_genes, tab_pathways, tab_predictor, tab_translational = st.tabs([
    "📊 Overview & PCA",
    "📋 Cohort Audit",
    "🌋 Volcano & Gene Explorer",
    "🌿 Pathway Dynamics",
    "⏱️ Transcriptomic Clock",
    "💡 Translational Insights"
])

# ----------------------------------------------------
# TAB 1: OVERVIEW & PCA
# ----------------------------------------------------
with tab_overview:
    st.subheader("Executive Metrics & Dimensionality Reduction")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown("""
        <div class="metric-box">
            <div class="metric-label">Discovery Cohort</div>
            <div class="metric-value">82 Samples</div>
            <div class="metric-sub">RNA-seq (Ages 22–89)</div>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="metric-box">
            <div class="metric-label">Validation Cohort</div>
            <div class="metric-value">42 Samples</div>
            <div class="metric-sub">21 Intra-Individual Pairs</div>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        st.markdown("""
        <div class="metric-box">
            <div class="metric-label">Discovery Model CV</div>
            <div class="metric-value">r = 0.71</div>
            <div class="metric-sub">MAE = 10.67 yrs (p < 0.001)</div>
        </div>
        """, unsafe_allow_html=True)
    with col4:
        st.markdown("""
        <div class="metric-box">
            <div class="metric-label">Harmonized Genes</div>
            <div class="metric-value">8,601</div>
            <div class="metric-sub">Cross-Platform Overlap</div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    # Interactive PCA Explorer
    st.markdown("### 🧭 Interactive Transcriptomic PCA Explorer")
    pca_file = os.path.join(TAB_DIR, 'gse226189_sample_qc_metrics.csv')

    if os.path.exists(pca_file):
        qc_df = pd.read_csv(pca_file)

        col_ctrl1, col_ctrl2 = st.columns([1, 3])
        with col_ctrl1:
            pca_dim = st.radio("PCA Visualization Mode:", ["2D Scatter", "3D Interactive"], horizontal=True)
            color_var = st.selectbox("Color Data Points By:", ["age", "sex", "lib_size_millions", "age_group"])

        with col_ctrl2:
            if pca_dim == "2D Scatter":
                fig_pca = px.scatter(
                    qc_df,
                    x='PC1',
                    y='PC2',
                    color=color_var,
                    hover_name='sample_id',
                    hover_data=['age', 'sex', 'lib_size_millions', 'spearman_mean_corr'],
                    title=f"GSE226189: 2D Principal Component Analysis (Colored by {color_var})",
                    color_continuous_scale='Viridis' if color_var in ['age', 'lib_size_millions'] else None,
                    opacity=0.85
                )
                fig_pca.update_traces(marker=dict(size=11, line=dict(width=1, color='DarkSlateGrey')))
                fig_pca.update_layout(height=520, template="plotly_white")
                st.plotly_chart(fig_pca, use_container_width=True)
            else:
                fig_pca3d = px.scatter_3d(
                    qc_df,
                    x='PC1',
                    y='PC2',
                    z='PC3',
                    color=color_var,
                    hover_name='sample_id',
                    hover_data=['age', 'sex', 'lib_size_millions'],
                    title=f"GSE226189: 3D PCA Space (Colored by {color_var})",
                    color_continuous_scale='Viridis' if color_var in ['age', 'lib_size_millions'] else None,
                    opacity=0.9
                )
                fig_pca3d.update_traces(marker=dict(size=6))
                fig_pca3d.update_layout(height=560, template="plotly_white")
                st.plotly_chart(fig_pca3d, use_container_width=True)

    # Paired Validation Cohort PCA
    st.markdown("### ☀️ Photoaging Paired Biopsies (GSE38308: Sun-Exposed vs Sun-Protected)")
    qc_val_file = os.path.join(TAB_DIR, 'gse38308_qc_metrics.csv')
    if os.path.exists(qc_val_file):
        qc_val_df = pd.read_csv(qc_val_file)

        fig_val_pca = px.scatter(
            qc_val_df,
            x='PC1',
            y='PC2',
            color='sun_exposure',
            symbol='sun_exposure',
            hover_name='sample_id',
            hover_data=['donor_id', 'age', 'anatomical_site'],
            title="GSE38308: PCA of Intra-Individual Matched Pairs (21 Donors)",
            color_discrete_map={'Sun-Exposed': '#d95f02', 'Sun-Protected': '#7570b3'}
        )
        fig_val_pca.update_traces(marker=dict(size=12, line=dict(width=1, color='DarkSlateGrey')))
        fig_val_pca.update_layout(height=480, template="plotly_white")
        st.plotly_chart(fig_val_pca, use_container_width=True)


# ----------------------------------------------------
# TAB 2: COHORT AUDIT
# ----------------------------------------------------
with tab_audit:
    st.subheader("Audited NCBI GEO Cohorts & Metadata Standards")
    st.markdown("""
    Every accession number, technology platform, sample count, and metadata attribute below was verified directly against the live NCBI Gene Expression Omnibus (GEO).
    """)

    audit_data = [
        {
            "Accession": "GSE226189",
            "Role": "Discovery Cohort",
            "Organism": "Homo sapiens",
            "Tissue": "Primary Skin Fibroblasts",
            "Sun Status": "Unspecified (Standard Biopsies)",
            "Technology": "RNA-seq (Illumina NovaSeq 6000)",
            "Samples": 82,
            "Age Range": "22–89y",
            "Raw Available": "Yes (geneCOUNT.txt.gz)",
            "Age per Sample": "Exact Continuous Integer",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226189"
        },
        {
            "Accession": "GSE38308",
            "Role": "Validation (Photoaging)",
            "Organism": "Homo sapiens",
            "Tissue": "Facial Skin (Pre- vs Post-auricular)",
            "Sun Status": "Matched Exposed vs Protected",
            "Technology": "Microarray (Illumina HumanWG-6)",
            "Samples": 42,
            "Age Range": "34–55y",
            "Raw Available": "Raw Probes + Processed Matrix",
            "Age per Sample": "Exact Continuous Integer",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38308"
        },
        {
            "Accession": "GSE18876",
            "Role": "Candidate In Vivo",
            "Organism": "Homo sapiens",
            "Tissue": "Lower Back Skin",
            "Sun Status": "Sun-protected",
            "Technology": "Microarray (Affymetrix HuEx-1.0)",
            "Samples": 98,
            "Age Range": "19–86y",
            "Raw Available": "Raw CEL files in TAR",
            "Age per Sample": "Exact Continuous Integer",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18876"
        },
        {
            "Accession": "GSE85358",
            "Role": "Candidate Epidermis",
            "Organism": "Homo sapiens",
            "Tissue": "Inner Forearm Epidermis",
            "Sun Status": "Sun-protected",
            "Technology": "Microarray (Agilent SurePrint G3)",
            "Samples": 48,
            "Age Range": "20–66y (Bimodal)",
            "Raw Available": "Raw text files in TAR",
            "Age per Sample": "Exact Continuous Integer",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85358"
        },
        {
            "Accession": "GSE67098",
            "Role": "Candidate Factorial",
            "Organism": "Homo sapiens",
            "Tissue": "Epidermis (dispase-separated)",
            "Sun Status": "Direct Contrast (Face vs Arm)",
            "Technology": "Microarray (Affymetrix U133 Plus 2)",
            "Samples": 16,
            "Age Range": "20–92y",
            "Raw Available": "Raw CEL files in TAR",
            "Age per Sample": "Exact Integer",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67098"
        },
        {
            "Accession": "GSE75337",
            "Role": "Candidate RNA-seq",
            "Organism": "Homo sapiens",
            "Tissue": "Skin punch biopsies & blood",
            "Sun Status": "Sun-protected",
            "Technology": "RNA-seq (Illumina HiSeq 2500)",
            "Samples": 91,
            "Age Range": "24–80y (4 brackets)",
            "Raw Available": "Matrix + SRA FASTQ",
            "Age per Sample": "Binned categories",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75337"
        },
        {
            "Accession": "GSE55118",
            "Role": "Candidate Fibroblasts",
            "Organism": "Homo sapiens",
            "Tissue": "Dermal Fibroblasts (Breast)",
            "Sun Status": "Sun-protected",
            "Technology": "Microarray (Agilent 4x44K)",
            "Samples": 15,
            "Age Range": "20–67y",
            "Raw Available": "Raw text in TAR",
            "Age per Sample": "Exact Continuous Integer",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55118"
        }
    ]
    audit_table = pd.DataFrame(audit_data)
    st.dataframe(audit_table, use_container_width=True)


# ----------------------------------------------------
# TAB 3: VOLCANO & GENE EXPLORER
# ----------------------------------------------------
with tab_genes:
    st.subheader("Volcano Plot Studio & Interactive Gene Profiler")

    study_toggle = st.radio(
        "Select Differential Expression Study:",
        ["GSE226189: Chronological Ageing (Continuous Model)", "GSE38308: Photoaging (Sun-Exposed vs Sun-Protected)"],
        horizontal=True
    )

    if "GSE226189" in study_toggle:
        de_file = os.path.join(TAB_DIR, 'gse226189_de_continuous_age.csv')
        if os.path.exists(de_file):
            de_df = pd.read_csv(de_file)

            col_v1, col_v2 = st.columns([1, 3])
            with col_v1:
                st.markdown("**Volcano Filter Controls:**")
                pval_thresh = st.slider("Raw p-value threshold:", min_value=0.0001, max_value=0.05, value=0.01, step=0.001, format="%.4f")
                fc_thresh = st.slider("Min |log2FC per decade|:", min_value=0.05, max_value=0.30, value=0.10, step=0.01)

                # Dynamically classify
                de_df['plot_sig'] = (de_df['p_value'] < pval_thresh) & (np.abs(de_df['log2FC_per_decade']) >= fc_thresh)
                n_up = ((de_df['plot_sig']) & (de_df['log2FC_per_decade'] > 0)).sum()
                n_down = ((de_df['plot_sig']) & (de_df['log2FC_per_decade'] < 0)).sum()

                st.metric("Age-Positive Genes", f"{n_up:,}")
                st.metric("Age-Negative Genes", f"{n_down:,}")

            with col_v2:
                # Interactive Plotly Volcano
                de_df['Category'] = 'Not Significant'
                de_df.loc[(de_df['plot_sig']) & (de_df['log2FC_per_decade'] > 0), 'Category'] = 'Age-Positive'
                de_df.loc[(de_df['plot_sig']) & (de_df['log2FC_per_decade'] < 0), 'Category'] = 'Age-Negative'

                fig_volc = px.scatter(
                    de_df,
                    x='log2FC_per_decade',
                    y='neg_log10_pval',
                    color='Category',
                    hover_name='symbol',
                    hover_data={'slope_per_year': ':.4f', 'p_value': ':.2e', 'fdr_qval': ':.3f', 'mean_expr': ':.2f'},
                    color_discrete_map={'Age-Positive': '#d95f02', 'Age-Negative': '#1f78b4', 'Not Significant': '#cbd5e1'},
                    title=f"GSE226189: Volcano Plot (Threshold: p < {pval_thresh}, |log2FC| >= {fc_thresh})",
                    labels={'log2FC_per_decade': 'log2 Fold Change per Decade', 'neg_log10_pval': '-log10(p-value)'},
                    opacity=0.75
                )
                fig_volc.add_hline(y=-np.log10(pval_thresh), line_dash="dash", line_color="gray")
                fig_volc.add_vline(x=fc_thresh, line_dash="dot", line_color="gray")
                fig_volc.add_vline(x=-fc_thresh, line_dash="dot", line_color="gray")
                fig_volc.update_layout(height=520, template="plotly_white")
                st.plotly_chart(fig_volc, use_container_width=True)

            # Single Gene Expression Profiler
            st.markdown("---")
            st.markdown("### 🔎 Single-Gene Expression & Age Correlation Profiler")
            cpm_file = os.path.join(DATA_PROC, 'gse226189_log2_cpm_symbols.csv')
            meta_file = os.path.join(DATA_META, 'gse226189_metadata.csv')

            if os.path.exists(cpm_file) and os.path.exists(meta_file):
                expr_mat = pd.read_csv(cpm_file, index_col=0)
                meta_df = pd.read_csv(meta_file)

                top_picks = ['WNT5A', 'CDH10', 'PTPRB', 'NEFH', 'CHN2', 'ERG', 'OLFM1', 'IGDCC4', 'C3orf70', 'PURG']
                query_gene = st.selectbox("Search or select gene symbol:", options=top_picks + sorted(list(set(expr_mat.index) - set(top_picks)))[:500])

                if query_gene in expr_mat.index:
                    gene_vals = expr_mat.loc[query_gene]
                    plot_data = meta_df[['sample_id', 'age', 'sex', 'age_group']].copy()
                    plot_data['log2_cpm'] = plot_data['sample_id'].map(gene_vals)

                    r_val, p_val = stats.pearsonr(plot_data['age'], plot_data['log2_cpm'])

                    col_g1, col_g2 = st.columns(2)
                    with col_g1:
                        fig_scatter = px.scatter(
                            plot_data, x='age', y='log2_cpm', color='sex',
                            trendline='ols',
                            title=f"{query_gene} vs Chronological Age (Pearson r = {r_val:.2f}, p = {p_val:.2e})",
                            labels={'age': 'Chronological Age (years)', 'log2_cpm': 'Expression (log2 CPM)'},
                            color_discrete_map={'Male': '#1f77b4', 'Female': '#e377c2'}
                        )
                        fig_scatter.update_traces(marker=dict(size=9, line=dict(width=1, color='DarkSlateGrey')))
                        fig_scatter.update_layout(height=420, template="plotly_white")
                        st.plotly_chart(fig_scatter, use_container_width=True)

                    with col_g2:
                        fig_box = px.box(
                            plot_data, x='age_group', y='log2_cpm', color='age_group',
                            points="all",
                            title=f"{query_gene} Distribution Across Age Cohorts",
                            labels={'age_group': 'Age Bracket', 'log2_cpm': 'Expression (log2 CPM)'},
                            color_discrete_sequence=['#93c5fd', '#60a5fa', '#1d4ed8']
                        )
                        fig_box.update_layout(height=420, template="plotly_white", showlegend=False)
                        st.plotly_chart(fig_box, use_container_width=True)

    else:
        de_val_file = os.path.join(TAB_DIR, 'gse38308_de_sun_exposure.csv')
        if os.path.exists(de_val_file):
            de_val_df = pd.read_csv(de_val_file)

            col_pv1, col_pv2 = st.columns([1, 3])
            with col_pv1:
                st.markdown("**Photoaging Filter Controls:**")
                fdr_thresh = st.slider("FDR threshold:", min_value=0.001, max_value=0.10, value=0.05, step=0.005)
                diff_thresh = st.slider("Min |log2 Difference|:", min_value=0.05, max_value=0.40, value=0.15, step=0.05)

                de_val_df['plot_sig'] = (de_val_df['fdr_qval'] < fdr_thresh) & (np.abs(de_val_df['mean_diff_log2FC']) >= diff_thresh)
                n_exp_up = ((de_val_df['plot_sig']) & (de_val_df['mean_diff_log2FC'] > 0)).sum()
                n_prot_up = ((de_val_df['plot_sig']) & (de_val_df['mean_diff_log2FC'] < 0)).sum()

                st.metric("Sun-Exposed Upregulated", f"{n_exp_up:,}")
                st.metric("Sun-Protected Upregulated", f"{n_prot_up:,}")

            with col_pv2:
                de_val_df['Category'] = 'Not Significant'
                de_val_df.loc[(de_val_df['plot_sig']) & (de_val_df['mean_diff_log2FC'] > 0), 'Category'] = 'Sun-Exposed Higher'
                de_val_df.loc[(de_val_df['plot_sig']) & (de_val_df['mean_diff_log2FC'] < 0), 'Category'] = 'Sun-Protected Higher'

                fig_val_volc = px.scatter(
                    de_val_df,
                    x='mean_diff_log2FC',
                    y='neg_log10_padj',
                    color='Category',
                    hover_name='symbol',
                    hover_data={'t_stat': ':.2f', 'p_value': ':.2e', 'fdr_qval': ':.3f', 'mean_expr': ':.2f'},
                    color_discrete_map={'Sun-Exposed Higher': '#e7298a', 'Sun-Protected Higher': '#10b981', 'Not Significant': '#cbd5e1'},
                    title="GSE38308: Photoaging Volcano [log2(Sun-Exposed) - log2(Sun-Protected)]",
                    labels={'mean_diff_log2FC': 'log2 Difference (Exposed vs Protected)', 'neg_log10_padj': '-log10(FDR Adjusted p-value)'},
                    opacity=0.75
                )
                fig_val_volc.add_hline(y=-np.log10(fdr_thresh), line_dash="dash", line_color="gray")
                fig_val_volc.add_vline(x=diff_thresh, line_dash="dot", line_color="gray")
                fig_val_volc.add_vline(x=-diff_thresh, line_dash="dot", line_color="gray")
                fig_val_volc.update_layout(height=520, template="plotly_white")
                st.plotly_chart(fig_val_volc, use_container_width=True)


# ----------------------------------------------------
# TAB 4: PATHWAY DYNAMICS
# ----------------------------------------------------
with tab_pathways:
    st.subheader("Functional Gene Set Enrichment Analysis (GSEA)")

    pw_choice = st.selectbox(
        "Choose Gene Set Collection:",
        ["MSigDB Hallmark (GSE226189 Ageing)", "Reactome Pathways (GSE226189 Ageing)", "MSigDB Hallmark (GSE38308 Sun Exposure)"]
    )

    if "Hallmark (GSE226189" in pw_choice:
        pw_file = os.path.join(TAB_DIR, 'gse226189_age_MSigDB_Hallmark_2020.csv')
    elif "Reactome" in pw_choice:
        pw_file = os.path.join(TAB_DIR, 'gse226189_age_Reactome_2022.csv')
    else:
        pw_file = os.path.join(TAB_DIR, 'gse38308_sun_exposure_MSigDB_Hallmark_2020.csv')

    if os.path.exists(pw_file):
        pw_df = pd.read_csv(pw_file)
        pw_df['neg_log10_fdr'] = -np.log10(np.maximum(pw_df['FDR q-val'], 1e-4))
        pw_df['Direction'] = np.where(pw_df['NES'] > 0, 'Upregulated with Age', 'Downregulated with Age')
        if 'Lead_genes' in pw_df.columns:
            pw_df['lead_gene_count'] = pw_df['Lead_genes'].astype(str).str.split(';').str.len()
        else:
            pw_df['lead_gene_count'] = 1

        # Bubble plot
        fig_bubble = px.scatter(
            pw_df.head(25),
            x='NES',
            y='neg_log10_fdr',
            size='lead_gene_count',
            color='Direction',
            hover_name='Term',
            hover_data=['NES', 'NOM p-val', 'FDR q-val', 'lead_gene_count'],
            title=f"Top Enriched Pathways in {pw_choice} (Bubble Chart)",
            labels={'NES': 'Normalized Enrichment Score (NES)', 'neg_log10_fdr': '-log10(FDR q-value)', 'lead_gene_count': 'Leading Edge Genes'},
            color_discrete_map={'Upregulated with Age': '#d95f02', 'Downregulated with Age': '#1f78b4'}
        )
        fig_bubble.add_vline(x=0, line_dash="dash", line_color="black")
        fig_bubble.add_hline(y=-np.log10(0.05), line_dash="dot", line_color="gray", annotation_text="FDR = 0.05")
        fig_bubble.update_layout(height=480, template="plotly_white")
        st.plotly_chart(fig_bubble, use_container_width=True)

        st.markdown("### Top Pathways Data Table")
        st.dataframe(pw_df[['Term', 'NES', 'NOM p-val', 'FDR q-val']].head(15), use_container_width=True)

    # Biological Mechanisms & Literature Cards
    expl_file = os.path.join(TAB_DIR, 'pathway_plain_language_explanations.csv')
    if os.path.exists(expl_file):
        st.markdown("---")
        st.markdown("### 🔬 Biological Mechanisms & Literature Reference Guide")
        expl_df = pd.read_csv(expl_file)
        for _, row in expl_df.iterrows():
            with st.expander(f"📌 {row['Pathway']} (NES: {row['NES']:.2f}, FDR: {row['FDR_qval']:.3e})"):
                st.markdown(f"**Biological Direction:** {row['Biological_Direction']}")
                st.markdown(f"**Mechanism Summary:** {row['Plain_Language_Summary']}")
                st.markdown(f"**Literature Confirmation:** {row['Literature_Lookup_Recommendation']}")


# ----------------------------------------------------
# TAB 5: TRANSCRIPTOMIC CLOCK
# ----------------------------------------------------
with tab_predictor:
    st.subheader("Elastic Net Transcriptomic Age Clock & Live Sample Simulator")

    # Metrics overview
    metrics_file = os.path.join(TAB_DIR, 'model_performance_metrics.csv')
    if os.path.exists(metrics_file):
        m_df = pd.read_csv(metrics_file)
        st.dataframe(m_df, use_container_width=True)

    # Interactive Sample Simulator
    st.markdown("### 🕹️ Live Sample Age Simulator")
    st.markdown("Select any biological donor from the Discovery cohort to inspect how the transcriptomic clock calculates their biological age compared to actual chronological age:")

    pred_file = os.path.join(TAB_DIR, 'gse226189_age_predictions.csv')
    if os.path.exists(pred_file):
        pred_df = pd.read_csv(pred_file)
        selected_sample = st.selectbox("Select Sample ID:", pred_df['sample_id'].values)

        sample_row = pred_df[pred_df['sample_id'] == selected_sample].iloc[0]
        act_age = sample_row['age']
        pr_age = sample_row['predicted_age_cv']
        delta = pr_age - act_age

        col_g1, col_g2 = st.columns([1, 2])
        with col_g1:
            fig_gauge = go.Figure(go.Indicator(
                mode="gauge+number+delta",
                value=pr_age,
                delta={'reference': act_age, 'position': "top", 'valueformat': "+.1f yrs"},
                gauge={
                    'axis': {'range': [10, 100]},
                    'bar': {'color': "#1e3a8a"},
                    'threshold': {'line': {'color': "red", 'width': 3}, 'thickness': 0.75, 'value': act_age}
                },
                title={'text': f"Predicted Age for {selected_sample}<br><span style='font-size:0.8em;color:gray'>Red line = Actual Chronological Age ({act_age}y)</span>"}
            ))
            fig_gauge.update_layout(height=350, margin=dict(l=20, r=20, t=50, b=20))
            st.plotly_chart(fig_gauge, use_container_width=True)

        with col_g2:
            st.markdown(f"**Sample Demographic Metadata:**")
            st.write(f"• **Sample ID:** `{selected_sample}`")
            st.write(f"• **Sex:** `{sample_row['sex']}`")
            st.write(f"• **Actual Chronological Age:** `{act_age} years`")
            st.write(f"• **Predicted Transcriptomic Age:** `{pr_age:.1f} years`")
            st.write(f"• **Prediction Delta:** `{delta:+.1f} years`")
            if abs(delta) < 5:
                st.success(f"High precision prediction (error < 5 years)!")
            else:
                st.info(f"Prediction within typical transcriptomic variance band (±10 years).")

    # Predicted vs Actual Plot
    st.markdown("---")
    st.markdown("### Predicted vs Actual Age Across Cohorts")
    fig11_p = os.path.join(FIG_DIR, 'fig11_predicted_vs_actual_age.png')
    if os.path.exists(fig11_p):
        st.image(fig11_p, caption="Figure: Elastic Net Predictions on Discovery 5-fold CV vs Independent Validation Cohort (Zero Refitting).")

    # Top Biomarker Genes
    st.markdown("### Top Transcriptomic Biomarkers (Model Weights)")
    fig12_p = os.path.join(FIG_DIR, 'fig12_elastic_net_gene_weights.png')
    if os.path.exists(fig12_p):
        st.image(fig12_p, caption="Figure: Top positive and negative Elastic Net regression coefficients.")


# ----------------------------------------------------
# TAB 6: TRANSLATIONAL INSIGHTS
# ----------------------------------------------------
with tab_translational:
    st.subheader("Translational Applications & Production Architecture")

    st.markdown("""
    ### 🎯 Applications in Dermatology & Biotech Drug Discovery:
    1. **Target Discovery for Extracellular Matrix (ECM) Reversal:**
       * Age-associated upregulation of `Protein Secretion` (SASP) and `EMT` pathways highlights matrix metalloproteinases and secretome components as key targets for preventing dermal collagen collapse.
    2. **Skin Barrier Lipid Restoration:**
       * Downregulation of `Cholesterol Homeostasis` ($NES = -1.69, FDR = 0.028$) underscores the biochemical basis of age-related barrier compromise, pointing toward topical lipid replenishment strategies.
    3. **Photoprotection & Senomorphic Drug Screening:**
       * The 1,269 identified photoaging genes provide a high-confidence transcriptomic readout to benchmark ultraviolet protective formulations and senomorphic compounds in human tissue.

    ---
    ### ⚙️ Production Architecture & Scalable Roadmap:
    * **Nextflow & Snakemake Pipeline Packaging:**
      * Containerize this analysis using Docker and Singularity for cloud execution on AWS HealthOmics or Google Cloud Life Sciences.
    * **Automated Batch Correction Microservices:**
      * Integrate empirical Bayes harmonisation (`ComBat`) into the ingestion pipeline to correct platform-level measurement scale shifts prior to inference.
    * **Multi-Cohort Expansion:**
      * Scale ingestion to the GTEx consortium ($n > 1,300$ skin RNA-seq samples) to establish baseline reference intervals across diverse populations.
    """)

st.markdown("---")
st.markdown("© 2026 SkinAge Atlas | Bioinformatics Software & Data Science Portfolio Project | Python, Streamlit & Plotly.")
