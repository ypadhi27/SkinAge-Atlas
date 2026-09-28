"""
SkinAge Atlas: Interactive Bioinformatics Dashboard
Reproducible cross-cohort transcriptomic analysis of human skin ageing,
pathway analysis, and independent-cohort biomarker validation.
"""

import os
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats

# Page configuration
st.set_page_config(
    page_title="SkinAge Atlas | Human Skin Transcriptomics",
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

# Custom CSS for modern styling
st.markdown("""
<style>
    .main-title {
        font-size: 2.3rem;
        font-weight: 700;
        color: #1e3a8a;
        margin-bottom: 0.2rem;
    }
    .sub-title {
        font-size: 1.15rem;
        color: #4b5563;
        margin-bottom: 1.5rem;
    }
    .disclaimer-banner {
        background-color: #fef3c7;
        border-left: 5px solid #f59e0b;
        padding: 0.75rem 1rem;
        border-radius: 4px;
        color: #92400e;
        font-weight: 500;
        font-size: 0.92rem;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 1rem;
        text-align: center;
    }
</style>
""", unsafe_allow_html=True)

# Top Banner & Disclaimer
st.markdown('<div class="main-title">🧬 SkinAge Atlas</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Cross-Cohort Transcriptomic Analysis of Human Skin Ageing, Photoaging & Machine Learning Age Prediction</div>', unsafe_allow_html=True)

st.markdown("""
<div class="disclaimer-banner">
    ⚠️ <strong>Research Prototype Notice:</strong> All reported correlations are transcriptomic associations, not causal mechanisms or clinical claims.
    This tool is intended strictly for scientific research and education. Not for clinical diagnostic use, medical advice, or cosmetic product claims.
</div>
""", unsafe_allow_html=True)


# Sidebar controls
st.sidebar.image("https://img.icons8.com/color/96/dna-helix.png", width=64)
st.sidebar.title("Navigation & Cohorts")
st.sidebar.markdown("""
**Cohorts Analysed:**
* **Discovery:** [GSE226189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226189) ($n=82$, RNA-seq)
* **Validation (Photoaging):** [GSE38308](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38308) ($n=42$, Microarray)
""")

# Tabs
tab_overview, tab_audit, tab_genes, tab_pathways, tab_predictor, tab_limitations = st.tabs([
    "📊 Overview",
    "📋 Dataset Audit",
    "🔍 Explore Genes",
    "🌿 Age Pathways",
    "🔮 Age Predictor",
    "⚖️ Limitations & PhD Scope"
])

# ----------------------------------------------------
# TAB 1: OVERVIEW
# ----------------------------------------------------
with tab_overview:
    st.subheader("Project Overview & Executive Summary")
    st.markdown("""
    **Primary Research Question:**
    *Which genes and functional pathways change systematically with chronological age in human skin, how do they differ between sun-exposed and sun-protected skin, and can a machine learning model trained on RNA-seq in vitro fibroblasts predict age in completely separate in vivo tissue?*
    """)

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric(label="Discovery Samples (GSE226189)", value="82", delta="RNA-seq (Ages 22-89)")
    with col2:
        st.metric(label="Validation Samples (GSE38308)", value="42", delta="Microarray (21 Paired Donors)")
    with col3:
        st.metric(label="Internal Model CV MAE", value="10.67 yrs", delta="r = 0.71 (p < 0.001)")
    with col4:
        st.metric(label="Shared Gene Symbols", value="8,601", delta="Standardized Features")

    st.markdown("---")
    st.markdown("### Cohort Quality Control & PCA Landscapes")

    col_fig1, col_fig2 = st.columns(2)
    with col_fig1:
        fig1_p = os.path.join(FIG_DIR, 'fig1_gse226189_age_and_library_qc.png')
        if os.path.exists(fig1_p):
            st.image(fig1_p, caption="Figure 1: GSE226189 Demographics (N=82) and Library Sizes across Age Groups.")
    with col_fig2:
        fig2_p = os.path.join(FIG_DIR, 'fig2_gse226189_pca.png')
        if os.path.exists(fig2_p):
            st.image(fig2_p, caption="Figure 2: Discovery Cohort PCA (PC1 vs PC2) with continuous age color scale.")

    st.markdown("### Sun-Exposed vs Sun-Protected Skin (GSE38308)")
    fig4_p = os.path.join(FIG_DIR, 'fig4_gse38308_qc_and_pca.png')
    if os.path.exists(fig4_p):
        st.image(fig4_p, caption="Figure 3: GSE38308 PCA showing paired Pre-auricular (Sun-Exposed) and Post-auricular (Sun-Protected) facial skin biopsies.")


# ----------------------------------------------------
# TAB 2: DATASET AUDIT
# ----------------------------------------------------
with tab_audit:
    st.subheader("Audited Candidate GEO Datasets")
    st.markdown("""
    All candidate datasets below were retrieved and audited directly from the live NCBI Gene Expression Omnibus (GEO).
    No sample counts, accessions, or platforms were assumed.
    """)

    audit_data = [
        {
            "Accession": "GSE226189",
            "Role": "Discovery Cohort",
            "Organism": "Homo sapiens",
            "Tissue": "Primary Skin Fibroblasts",
            "Sun Status": "Unspecified",
            "Technology": "RNA-seq (Illumina NovaSeq 6000)",
            "Samples": 82,
            "Age Range": "22-89y",
            "Format": "Raw counts + Matrix",
            "Age per Sample": "Yes (Continuous)",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226189"
        },
        {
            "Accession": "GSE38308",
            "Role": "Validation (Photoaging)",
            "Organism": "Homo sapiens",
            "Tissue": "Facial Skin (pre- vs post-auricular)",
            "Sun Status": "Matched Exposed vs Protected",
            "Technology": "Microarray (Illumina HumanWG-6)",
            "Samples": 42,
            "Age Range": "34-55y",
            "Format": "Normalised Matrix + Raw Probes",
            "Age per Sample": "Yes (Continuous)",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38308"
        },
        {
            "Accession": "GSE18876",
            "Role": "Candidate In Vivo Skin",
            "Organism": "Homo sapiens",
            "Tissue": "Lower Back Skin",
            "Sun Status": "Sun-protected",
            "Technology": "Microarray (Affymetrix HuEx-1.0 ST)",
            "Samples": 98,
            "Age Range": "19-86y",
            "Format": "Normalised Matrix + RAW.tar",
            "Age per Sample": "Yes (Continuous)",
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
            "Age Range": "20-66y (Bimodal)",
            "Format": "Normalised Matrix + RAW.tar",
            "Age per Sample": "Yes (Continuous)",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85358"
        },
        {
            "Accession": "GSE67098",
            "Role": "Candidate Factorial",
            "Organism": "Homo sapiens",
            "Tissue": "Epidermis (dispase-separated)",
            "Sun Status": "Direct Contrast (Peri-orbital vs Arm)",
            "Technology": "Microarray (Affymetrix U133 Plus 2.0)",
            "Samples": 16,
            "Age Range": "20-92y",
            "Format": "Normalised Matrix + RAW.tar",
            "Age per Sample": "Yes (Integer)",
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
            "Age Range": "24-80y (4 brackets)",
            "Format": "Normalised Matrix + SRA FASTQ",
            "Age per Sample": "Binned cohorts",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75337"
        },
        {
            "Accession": "GSE55118",
            "Role": "Candidate Small Cohort",
            "Organism": "Homo sapiens",
            "Tissue": "Dermal Fibroblasts (Breast)",
            "Sun Status": "Sun-protected",
            "Technology": "Microarray (Agilent 4x44K)",
            "Samples": 15,
            "Age Range": "20-67y",
            "Format": "Normalised Matrix + RAW.tar",
            "Age per Sample": "Yes (Continuous)",
            "GEO Link": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55118"
        }
    ]
    audit_df = pd.DataFrame(audit_data)
    st.dataframe(audit_df, use_container_width=True)

    st.markdown("""
    **Cohort Selection Rationale:**
    * **Discovery (GSE226189):** Chosen because it provides high-throughput RNA-seq data across a broad, continuous age spectrum (22–89 years) with balanced sexes (35F/47M).
    * **Independent Validation (GSE38308):** Chosen because its 21 intra-individual matched pairs (pre-auricular sun-exposed vs post-auricular sun-protected) uniquely allow distinguishing intrinsic chronological aging from cumulative extrinsic ultraviolet photoaging without inter-individual genetic confounding.
    """)


# ----------------------------------------------------
# TAB 3: EXPLORE GENES
# ----------------------------------------------------
with tab_genes:
    st.subheader("Interactive Differential Expression & Gene Explorer")

    dataset_choice = st.radio("Select Analysis:", ["GSE226189 (Chronological Ageing)", "GSE38308 (Photoaging: Sun-Exposed vs Protected)"], horizontal=True)

    if "GSE226189" in dataset_choice:
        de_file = os.path.join(TAB_DIR, 'gse226189_de_continuous_age.csv')
        if os.path.exists(de_file):
            de_df = pd.read_csv(de_file)
            st.markdown(f"**Discovery Model:** Linear regression `Expression ~ Age + Sex` across **{len(de_df):,}** genes in 82 donors.")

            # Volcano plot (Plotly)
            fig_volc = px.scatter(
                de_df,
                x='log2FC_per_decade',
                y='neg_log10_pval',
                hover_name='symbol',
                hover_data={'t_stat': ':.2f', 'p_value': ':.2e', 'fdr_qval': ':.3f', 'mean_expr': ':.2f'},
                color=de_df['p_value'] < 0.01,
                color_discrete_map={True: '#d95f02', False: '#bbbbbb'},
                title="Volcano Plot: Chronological Age Effect Size vs Statistical Significance",
                labels={'log2FC_per_decade': 'log2 Fold Change per Decade', 'neg_log10_pval': '-log10(p-value)'},
                opacity=0.75
            )
            fig_volc.add_hline(y=-np.log10(0.01), line_dash="dash", line_color="gray", annotation_text="p = 0.01")
            fig_volc.update_layout(showlegend=False, height=500)
            st.plotly_chart(fig_volc, use_container_width=True)

            # Gene query tool
            st.markdown("### Search Specific Gene Expression Profile")
            cpm_file = os.path.join(DATA_PROC, 'gse226189_log2_cpm_symbols.csv')
            meta_file = os.path.join(DATA_META, 'gse226189_metadata.csv')

            if os.path.exists(cpm_file) and os.path.exists(meta_file):
                expr_mat = pd.read_csv(cpm_file, index_col=0)
                meta_df = pd.read_csv(meta_file)

                default_gene = "WNT5A" if "WNT5A" in expr_mat.index else expr_mat.index[0]
                query_gene = st.selectbox("Type or select a gene symbol:", sorted(list(expr_mat.index[:500])), index=0 if default_gene not in expr_mat.index[:500] else list(sorted(list(expr_mat.index[:500]))).index(default_gene))

                if query_gene in expr_mat.index:
                    gene_vals = expr_mat.loc[query_gene]
                    plot_data = meta_df[['sample_id', 'age', 'sex']].copy()
                    plot_data['log2_cpm'] = plot_data['sample_id'].map(gene_vals)

                    r_val, p_val = stats.pearsonr(plot_data['age'], plot_data['log2_cpm'])

                    fig_scatter = px.scatter(
                        plot_data, x='age', y='log2_cpm', color='sex',
                        trendline='ols',
                        title=f"{query_gene} Expression vs Chronological Age (r = {r_val:.2f}, p = {p_val:.2e})",
                        labels={'age': 'Chronological Age (years)', 'log2_cpm': 'Expression (log2 CPM)'},
                        color_discrete_map={'Male': '#1f77b4', 'Female': '#e377c2'}
                    )
                    fig_scatter.update_layout(height=450)
                    st.plotly_chart(fig_scatter, use_container_width=True)

    else:
        de_val_file = os.path.join(TAB_DIR, 'gse38308_de_sun_exposure.csv')
        if os.path.exists(de_val_file):
            de_val_df = pd.read_csv(de_val_file)
            st.markdown(f"**Validation Photoaging Model:** Paired t-test (Sun-Exposed vs Sun-Protected) across **{len(de_val_df):,}** genes in 21 matched donors.")

            fig_volc_val = px.scatter(
                de_val_df,
                x='mean_diff_log2FC',
                y='neg_log10_padj',
                hover_name='symbol',
                hover_data={'t_stat': ':.2f', 'p_value': ':.2e', 'fdr_qval': ':.3f', 'mean_expr': ':.2f'},
                color=de_val_df['fdr_qval'] < 0.05,
                color_discrete_map={True: '#e7298a', False: '#bbbbbb'},
                title="Volcano Plot: Photoaging Differences [log2(Exposed) - log2(Protected)]",
                labels={'mean_diff_log2FC': 'log2 Difference (Sun-Exposed vs Sun-Protected)', 'neg_log10_padj': '-log10(FDR Adjusted p-value)'},
                opacity=0.75
            )
            fig_volc_val.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="gray", annotation_text="FDR = 0.05")
            fig_volc_val.update_layout(showlegend=False, height=500)
            st.plotly_chart(fig_volc_val, use_container_width=True)


# ----------------------------------------------------
# TAB 4: AGE-ASSOCIATED PATHWAYS
# ----------------------------------------------------
with tab_pathways:
    st.subheader("Functional Pathway Enrichment (GSEA)")

    pw_choice = st.selectbox(
        "Choose Gene Set Library:",
        ["MSigDB Hallmark (GSE226189 Ageing)", "Reactome Pathways (GSE226189 Ageing)", "MSigDB Hallmark (GSE38308 Sun Exposure)"]
    )

    if "Hallmark (GSE226189" in pw_choice:
        pw_file = os.path.join(TAB_DIR, 'gse226189_age_MSigDB_Hallmark_2020.csv')
        fig_p = os.path.join(FIG_DIR, 'fig8_gse226189_hallmark_gsea.png')
    elif "Reactome" in pw_choice:
        pw_file = os.path.join(TAB_DIR, 'gse226189_age_Reactome_2022.csv')
        fig_p = os.path.join(FIG_DIR, 'fig9_gse226189_reactome_gsea.png')
    else:
        pw_file = os.path.join(TAB_DIR, 'gse38308_sun_exposure_MSigDB_Hallmark_2020.csv')
        fig_p = os.path.join(FIG_DIR, 'fig10_gse38308_hallmark_gsea.png')

    if os.path.exists(fig_p):
        st.image(fig_p, use_container_width=True)

    if os.path.exists(pw_file):
        pw_df = pd.read_csv(pw_file)
        st.markdown("### Top Enriched Pathways Table")
        st.dataframe(pw_df[['Term', 'NES', 'NOM p-val', 'FDR q-val']].head(20), use_container_width=True)

    # Biological Plain-Language Guide
    expl_file = os.path.join(TAB_DIR, 'pathway_plain_language_explanations.csv')
    if os.path.exists(expl_file):
        st.markdown("### Plain-Language Biological Interpretation & Literature Guide")
        expl_df = pd.read_csv(expl_file)
        for _, row in expl_df.iterrows():
            with st.expander(f"📌 {row['Pathway']} (NES: {row['NES']:.2f}, FDR: {row['FDR_qval']:.3e})"):
                st.markdown(f"**Biological Direction:** {row['Biological_Direction']}")
                st.markdown(f"**Plain-Language Explanation:** {row['Plain_Language_Summary']}")
                st.markdown(f"**How to Confirm in Literature:** {row['Literature_Lookup_Recommendation']}")


# ----------------------------------------------------
# TAB 5: AGE PREDICTOR
# ----------------------------------------------------
with tab_predictor:
    st.subheader("Machine Learning Transcriptomic Age Predictor (Elastic Net)")

    st.markdown("""
    An **Elastic Net regression model** ($L_1 + L_2$ regularisation) was trained on the **GSE226189 Discovery cohort**
    using strict 5-fold cross-validation. The locked model was then evaluated on the completely independent **GSE38308 cohort**
    with **ZERO refitting or tuning** to rigorously test cross-platform transferability.
    """)

    metrics_file = os.path.join(TAB_DIR, 'model_performance_metrics.csv')
    if os.path.exists(metrics_file):
        m_df = pd.read_csv(metrics_file)
        st.dataframe(m_df, use_container_width=True)

    col_pred1, col_pred2 = st.columns(2)
    with col_pred1:
        fig11_p = os.path.join(FIG_DIR, 'fig11_predicted_vs_actual_age.png')
        if os.path.exists(fig11_p):
            st.image(fig11_p, caption="Figure 5: Predicted vs Actual Chronological Age (Training 5-fold CV vs Independent Test).")

    with col_pred2:
        fig12_p = os.path.join(FIG_DIR, 'fig12_elastic_net_gene_weights.png')
        if os.path.exists(fig12_p):
            st.image(fig12_p, caption="Figure 6: Top Elastic Net model weights (transcriptomic age biomarkers).")

    st.markdown("""
    ### Honest Performance Evaluation & Scientific Reality:
    1. **Internal Cross-Validation (Discovery GSE226189):**
       * The model achieves a **Mean Absolute Error (MAE) of 10.67 years** (compared to 16.52 years for the trivial mean baseline) and a Pearson correlation of **$r = 0.71$ ($p = 1.3 \\times 10^{-13}$)**.
    2. **Independent Test (Validation GSE38308):**
       * When tested directly on GSE38308 without refitting, the model fails to predict chronological age accurately (MAE = 48.85 years).
       * **Why did cross-cohort transfer fail?**
         * **Platform divergence:** RNA-seq discrete counts (Illumina NovaSeq) vs Microarray probe intensities (Illumina BeadChip).
         * **Tissue model discrepancy:** Primary fibroblast cell culture in vitro (monoculture) vs whole full-thickness facial skin in vivo (heterogeneous mixture of keratinocytes, fibroblasts, melanocytes, endothelial cells).
         * **Age span restriction:** GSE38308 spans only 34–55 years (standard deviation = 5.2 years), where predicting the cohort mean is already very close to the true age.
    """)


# ----------------------------------------------------
# TAB 6: LIMITATIONS & PHD PROPOSAL SCOPE
# ----------------------------------------------------
with tab_limitations:
    st.subheader("Methodological Limitations & What I Would Do Next")

    st.markdown("""
    ### Key Methodological Limitations:
    1. **Cell Culture vs Tissue Discrepancy:**
       * In vitro primary fibroblasts (GSE226189) lack epidermal-dermal cross-talk, immune infiltration, systemic hormonal influences, and mechanical skin tension.
    2. **Cohort Sex Imbalance in Validation:**
       * The validation cohort (GSE38308) consists exclusively of female donors, precluding cross-cohort sex interaction analysis.
    3. **Microarray vs RNA-seq Scale Differences:**
       * Standardizing across different measurement technologies without batch adjustment bridges symbols, but residual probe-hybridization differences remain.

    ---
    ### Future Directions (Seed for a PhD Proposal):
    * **Single-Cell Transcriptomic Deconvolution (scRNA-seq):**
      * Use scRNA-seq references (e.g., Tabula Sapiens, Aging Atlas) to deconvolve bulk whole-skin biopsies into proportional fibroblast, keratinocyte, melanocyte, and T-cell age clocks.
    * **Multi-Omics Integration (Transcriptomics + Epigenomics):**
      * Combine DNA methylation (Horvath skin & blood clock) with RNA-seq gene expression using Multi-Omics Factor Analysis (MOFA+) to decouple intrinsic chronological ageing from UV-induced photoaging.
    * **Cross-Tissue Generalization (GTEx Skin Benchmark):**
      * Scale the validation framework to the GTEx consortium ($n > 1,300$ post-mortem skin samples from sun-exposed lower leg and sun-protected suprapubic skin) to establish gold-standard population baselines across Sweden and Nordic cohorts.
    """)

st.markdown("---")
st.markdown("© 2026 SkinAge Atlas | MSc Bioinformatics Portfolio Project | Developed with Python, PyDESeq2, GSEAPy, Scikit-Learn & Streamlit.")
