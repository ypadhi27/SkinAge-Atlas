import os
import json
import pandas as pd

def build_demo_html():
    os.makedirs('docs', exist_ok=True)

    # 1. Continuous Age DE
    de_age = pd.read_csv('results/tables/gse226189_de_continuous_age.csv')
    age_genes = de_age.sort_values('p_value').head(250)[
        ['symbol', 'mean_expr', 'log2FC_per_decade', 't_stat', 'p_value', 'fdr_qval', 'neg_log10_pval']
    ].round(4).to_dict(orient='records')

    # 2. Photoaging DE
    de_photo = pd.read_csv('results/tables/gse38308_de_sun_exposure.csv')
    photo_genes = de_photo.sort_values('p_value').head(150)[
        ['symbol', 'mean_expr', 'mean_diff_log2FC', 't_stat', 'p_value', 'fdr_qval', 'neg_log10_pval']
    ].rename(columns={'mean_diff_log2FC': 'log2FC_per_decade'}).round(4).to_dict(orient='records')

    # 3. Pathways
    pw_df = pd.read_csv('results/tables/gse226189_age_MSigDB_Hallmark_2020.csv')
    pathways = pw_df.head(15)[['Term', 'NES', 'NOM p-val', 'FDR q-val']].round(4).to_dict(orient='records')

    # 4. Explanations
    expl_df = pd.read_csv('results/tables/pathway_plain_language_explanations.csv')
    explanations = expl_df.to_dict(orient='records')

    # 5. Patient predictions
    preds_df = pd.read_csv('results/tables/gse226189_age_predictions.csv')
    patients = preds_df.sort_values('age').round(2).to_dict(orient='records')

    # 6. Coefficients
    coef_df = pd.read_csv('results/tables/elastic_net_gene_coefficients.csv')
    coef_df['abs_coef'] = coef_df['coefficient'].abs()
    top_coefs = coef_df[coef_df['coefficient'] != 0].sort_values('abs_coef', ascending=False).head(25)[
        ['symbol', 'coefficient']
    ].round(4).to_dict(orient='records')

    data_payload = {
        "ageGenes": age_genes,
        "photoGenes": photo_genes,
        "pathways": pathways,
        "explanations": explanations,
        "patients": patients,
        "coefficients": top_coefs
    }
    json_str = json.dumps(data_payload)

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>SkinAge Atlas — Interactive Web Demo</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600;700&family=Plus+Jakarta+Sans:wght@400;500;600;700;800&display=swap" rel="stylesheet">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    :root {{
      --bg: #070a13;
      --card-bg: rgba(15, 23, 42, 0.75);
      --border: rgba(255, 255, 255, 0.08);
      --cyan: #38bdf8;
      --rose: #f43f5e;
      --emerald: #10b981;
      --text: #f8fafc;
      --muted: #94a3b8;
    }}
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{
      background: var(--bg);
      color: var(--text);
      font-family: 'Plus Jakarta Sans', system-ui, sans-serif;
      min-height: 100vh;
      overflow-x: hidden;
      line-height: 1.5;
    }}
    .container {{
      max-width: 1400px;
      margin: 0 auto;
      padding: 32px 24px 80px;
    }}
    .header-banner {{
      display: flex;
      flex-direction: column;
      align-items: center;
      text-align: center;
      margin-bottom: 32px;
      padding: 36px 20px;
      border-radius: 24px;
      background: radial-gradient(circle at 50% 0%, rgba(56, 189, 248, 0.15) 0%, rgba(15, 23, 42, 0) 70%), #0b1120;
      border: 1px solid var(--border);
    }}
    .pill-badge {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      padding: 6px 16px;
      border-radius: 9999px;
      background: rgba(56, 189, 248, 0.1);
      border: 1px solid rgba(56, 189, 248, 0.3);
      color: var(--cyan);
      font-family: 'JetBrains Mono', monospace;
      font-size: 13px;
      font-weight: 600;
      letter-spacing: 0.05em;
      text-transform: uppercase;
      margin-bottom: 16px;
    }}
    .h1-title {{
      font-size: 44px;
      font-weight: 800;
      letter-spacing: -0.03em;
      line-height: 1.15;
      background: linear-gradient(180deg, #ffffff 40%, #cbd5e1 100%);
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
      margin-bottom: 12px;
    }}
    .subtitle {{
      color: var(--muted);
      font-size: 18px;
      max-width: 850px;
      margin-bottom: 24px;
    }}
    .kpi-row {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 16px;
      width: 100%;
      margin-bottom: 32px;
    }}
    .kpi-card {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 16px;
      padding: 20px;
      text-align: center;
      backdrop-filter: blur(12px);
    }}
    .kpi-num {{
      font-family: 'JetBrains Mono', monospace;
      font-size: 32px;
      font-weight: 800;
      color: var(--cyan);
      margin-bottom: 4px;
    }}
    .kpi-label {{
      font-size: 13px;
      font-weight: 600;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }}
    .tabs-nav {{
      display: flex;
      gap: 12px;
      border-bottom: 1px solid var(--border);
      padding-bottom: 12px;
      margin-bottom: 28px;
      overflow-x: auto;
    }}
    .tab-btn {{
      padding: 12px 24px;
      border-radius: 12px;
      background: transparent;
      border: 1px solid transparent;
      color: var(--muted);
      font-weight: 600;
      font-size: 15px;
      cursor: pointer;
      transition: all 0.2s ease;
      white-space: nowrap;
    }}
    .tab-btn:hover {{
      color: var(--text);
      background: rgba(255, 255, 255, 0.04);
    }}
    .tab-btn.active {{
      background: rgba(56, 189, 248, 0.12);
      border-color: rgba(56, 189, 248, 0.4);
      color: var(--cyan);
    }}
    .tab-pane {{
      display: none;
    }}
    .tab-pane.active {{
      display: block;
      animation: fadeIn 0.3s ease;
    }}
    @keyframes fadeIn {{
      from {{ opacity: 0; transform: translateY(6px); }}
      to {{ opacity: 1; transform: translateY(0); }}
    }}
    .grid-2col {{
      display: grid;
      grid-template-columns: 1fr 340px;
      gap: 24px;
    }}
    @media (max-width: 1024px) {{
      .grid-2col {{ grid-template-columns: 1fr; }}
    }}
    .control-panel {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 20px;
      padding: 24px;
      height: fit-content;
    }}
    .control-label {{
      display: flex;
      justify-content: space-between;
      font-size: 14px;
      font-weight: 600;
      color: #cbd5e1;
      margin-bottom: 8px;
    }}
    .slider-input {{
      width: 100%;
      accent-color: var(--cyan);
      margin-bottom: 20px;
      cursor: pointer;
    }}
    .search-input {{
      width: 100%;
      background: #0f172a;
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 10px 14px;
      color: var(--text);
      font-family: inherit;
      font-size: 14px;
      margin-bottom: 20px;
    }}
    .search-input:focus {{
      outline: none;
      border-color: var(--cyan);
    }}
    .plot-container {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 20px;
      padding: 20px;
      min-height: 520px;
    }}
    .data-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 14px;
      margin-top: 16px;
    }}
    .data-table th {{
      text-align: left;
      padding: 10px 12px;
      background: rgba(255, 255, 255, 0.03);
      color: var(--muted);
      border-bottom: 1px solid var(--border);
      font-family: 'JetBrains Mono', monospace;
    }}
    .data-table td {{
      padding: 10px 12px;
      border-bottom: 1px solid rgba(255, 255, 255, 0.04);
      font-family: 'JetBrains Mono', monospace;
    }}
    .card-list {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
      gap: 20px;
    }}
    .pathway-card {{
      background: var(--card-bg);
      border: 1px solid var(--border);
      border-radius: 18px;
      padding: 24px;
      border-left: 4px solid var(--cyan);
    }}
    .footer-bar {{
      margin-top: 60px;
      padding-top: 24px;
      border-top: 1px solid var(--border);
      display: flex;
      justify-content: space-between;
      align-items: center;
      flex-wrap: wrap;
      gap: 16px;
      color: var(--muted);
      font-size: 14px;
    }}
    .cta-btn {{
      display: inline-flex;
      align-items: center;
      gap: 8px;
      background: rgba(56, 189, 248, 0.15);
      border: 1px solid var(--cyan);
      color: var(--cyan);
      padding: 10px 20px;
      border-radius: 9999px;
      text-decoration: none;
      font-weight: 700;
      transition: all 0.2s;
    }}
    .cta-btn:hover {{
      background: var(--cyan);
      color: #070a13;
    }}
  </style>
</head>
<body>
  <div class="container">
    <!-- HEADER -->
    <div class="header-banner">
      <div class="pill-badge">Live Interactive Portal • Client-Side Zero Install</div>
      <h1 class="h1-title">SkinAge Atlas: Transcriptomic Explorer</h1>
      <p class="subtitle">
        Cross-cohort analysis of human skin chronological ageing & UV photoaging with pathway discovery and leak-free machine learning biomarker validation.
      </p>
      <div style="display: flex; gap: 16px; flex-wrap: wrap;">
        <a href="https://github.com/ypadhi27/SkinAge-Atlas" target="_blank" class="cta-btn">
          <span>GitHub Repository →</span>
        </a>
        <a href="https://raw.githubusercontent.com/ypadhi27/SkinAge-Atlas/main/brag-output/brag.mp4" target="_blank" class="cta-btn" style="border-color: #f43f5e; color: #f43f5e; background: rgba(244, 63, 94, 0.1);">
          <span>▶️ Watch 20s Demo Video</span>
        </a>
      </div>
    </div>

    <!-- KPI STATS -->
    <div class="kpi-row">
      <div class="kpi-card">
        <div class="kpi-num">124</div>
        <div class="kpi-label">Human Donors</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-num">20,798</div>
        <div class="kpi-label">Genes Modeled</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-num" style="color: var(--rose);">1,269</div>
        <div class="kpi-label">Photoaging Genes</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-num" style="color: var(--emerald);">0.706</div>
        <div class="kpi-label">ML Clock Pearson r</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-num">10.67y</div>
        <div class="kpi-label">Cross-Val MAE</div>
      </div>
    </div>

    <!-- TABS NAVIGATION -->
    <div class="tabs-nav">
      <button class="tab-btn active" onclick="switchTab('tab-volcano')">🌋 Interactive Volcano Plot</button>
      <button class="tab-btn" onclick="switchTab('tab-clock')">⏱️ Patient Age Clock Simulator</button>
      <button class="tab-btn" onclick="switchTab('tab-pathways')">🔬 GSEA Pathway Browser</button>
      <button class="tab-btn" onclick="switchTab('tab-methods')">📖 Scientific Overview & Architecture</button>
    </div>

    <!-- TAB 1: VOLCANO PLOT -->
    <div id="tab-volcano" class="tab-pane active">
      <div class="grid-2col">
        <div class="plot-container">
          <div id="volcanoPlot" style="width: 100%; height: 500px;"></div>
        </div>
        <div class="control-panel">
          <h3 style="font-size: 18px; margin-bottom: 16px;">Filter Controls</h3>
          
          <div style="margin-bottom: 16px;">
            <label class="control-label">Cohort Analysis:</label>
            <select id="cohortSelect" class="search-input" onchange="updateVolcanoPlot()">
              <option value="age">Chronological Ageing (GSE226189)</option>
              <option value="photo">UV Photoaging (GSE38308 Matched)</option>
            </select>
          </div>

          <div>
            <div class="control-label">
              <span>Log2 Fold Change Threshold</span>
              <span id="fcVal" style="color: var(--cyan);">0.10</span>
            </div>
            <input type="range" id="fcSlider" class="slider-input" min="0.05" max="0.5" step="0.01" value="0.10" oninput="updateVolcanoPlot()">
          </div>

          <div>
            <div class="control-label">
              <span>-Log10 P-value Threshold</span>
              <span id="pvalVal" style="color: var(--cyan);">2.0</span>
            </div>
            <input type="range" id="pvalSlider" class="slider-input" min="1.0" max="6.0" step="0.1" value="2.0" oninput="updateVolcanoPlot()">
          </div>

          <div>
            <label class="control-label">Search Gene Symbol:</label>
            <input type="text" id="geneSearch" class="search-input" placeholder="e.g. CDH10, WNT5A, CYP4B1" oninput="updateVolcanoPlot()">
          </div>

          <div style="background: rgba(255,255,255,0.03); padding: 14px; border-radius: 12px; font-size: 13px;">
            <div style="color: var(--muted); margin-bottom: 4px;">Filtered Active Genes:</div>
            <div id="filterCount" style="font-family: 'JetBrains Mono'; font-size: 20px; font-weight: 700; color: var(--emerald);">0</div>
          </div>
        </div>
      </div>

      <!-- GENE TABLE -->
      <div class="plot-container" style="margin-top: 24px; min-height: auto;">
        <h3 style="font-size: 18px; margin-bottom: 12px;">Top Filtered Differential Expression Candidates</h3>
        <div style="overflow-x: auto;">
          <table class="data-table">
            <thead>
              <tr>
                <th>Gene Symbol</th>
                <th>Log2 FC / Effect</th>
                <th>Mean Expression</th>
                <th>T-Statistic</th>
                <th>P-Value</th>
                <th>FDR Q-Value</th>
                <th>Biological Regulation</th>
              </tr>
            </thead>
            <tbody id="geneTableBody"></tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- TAB 2: AGE CLOCK SIMULATOR -->
    <div id="tab-clock" class="tab-pane">
      <div class="grid-2col">
        <div class="plot-container">
          <div id="gaugePlot" style="width: 100%; height: 420px;"></div>
          <div id="patientDetails" style="margin-top: 16px; padding: 16px; background: rgba(255,255,255,0.03); border-radius: 12px;"></div>
        </div>
        <div class="control-panel">
          <h3 style="font-size: 18px; margin-bottom: 16px;">Select Clinical Donor</h3>
          <p style="font-size: 13px; color: var(--muted); margin-bottom: 16px;">
            Choose from 82 real primary skin fibroblast RNA-seq samples evaluated with strictly leak-free 5-fold cross-validation.
          </p>
          <select id="patientSelect" class="search-input" onchange="updateGauge()"></select>

          <h4 style="font-size: 15px; margin-top: 24px; margin-bottom: 12px; color: var(--cyan);">Top Elastic Net Feature Weights</h4>
          <div style="max-height: 260px; overflow-y: auto;">
            <table class="data-table" style="font-size: 12px;">
              <thead>
                <tr>
                  <th>Biomarker</th>
                  <th>Weight</th>
                </tr>
              </thead>
              <tbody id="coefTableBody"></tbody>
            </table>
          </div>
        </div>
      </div>
    </div>

    <!-- TAB 3: PATHWAYS -->
    <div id="tab-pathways" class="tab-pane">
      <div class="plot-container" style="margin-bottom: 24px;">
        <div id="pathwayPlot" style="width: 100%; height: 450px;"></div>
      </div>
      <h3 style="font-size: 20px; margin-bottom: 16px;">Biological Mechanism Summaries & Literature Guide</h3>
      <div id="pathwayCards" class="card-list"></div>
    </div>

    <!-- TAB 4: METHODS & ARCHITECTURE -->
    <div id="tab-methods" class="tab-pane">
      <div class="plot-container" style="min-height: auto;">
        <h2 style="font-size: 28px; margin-bottom: 16px;">Pipeline Architecture & Scientific Governance</h2>
        <p style="color: var(--muted); font-size: 16px; line-height: 1.7; margin-bottom: 24px;">
          SkinAge Atlas is built with strict reproducibility standards: zero fabricated datasets, zero synthetic p-values, and leak-free machine learning protocols.
        </p>

        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 20px;">
          <div class="kpi-card" style="text-align: left; padding: 24px;">
            <h3 style="color: var(--cyan); margin-bottom: 12px;">1. Continuous Modeling</h3>
            <p style="font-size: 14px; color: var(--muted); line-height: 1.6;">
              Rather than binning ages into arbitrary groups, we fitted an ordinary least squares regression against continuous chronological age while adjusting for biological sex across all 82 donors (GSE226189).
            </p>
          </div>
          <div class="kpi-card" style="text-align: left; padding: 24px;">
            <h3 style="color: var(--rose); margin-bottom: 12px;">2. Paired Photoaging</h3>
            <p style="font-size: 14px; color: var(--muted); line-height: 1.6;">
              In GSE38308, we leveraged 21 intra-individual matched pairs (pre- vs post-auricular facial biopsies) to isolate solar ultraviolet effects while controlling completely for genetic background.
            </p>
          </div>
          <div class="kpi-card" style="text-align: left; padding: 24px;">
            <h3 style="color: var(--emerald); margin-bottom: 12px;">3. Honest Transfer Assessment</h3>
            <p style="font-size: 14px; color: var(--muted); line-height: 1.6;">
              The Elastic Net model trained on in vitro fibroblasts failed to transfer directly to whole-skin microarrays (MAE = 48.85y), uncovering the vital role of cell-type composition dilution in clinical biomarker translation.
            </p>
          </div>
        </div>
      </div>
    </div>

    <!-- FOOTER -->
    <div class="footer-bar">
      <div>
        <b>SkinAge Atlas</b> • Portfolio Project developed by <a href="https://github.com/ypadhi27" style="color: var(--cyan); text-decoration: none;">Yajnasenee Padhi</a>
      </div>
      <div>
        Built with Python, PyDESeq2, GSEAPy, Scikit-Learn & Plotly.js
      </div>
    </div>
  </div>

  <!-- SCRIPT DATA & LOGIC -->
  <script>
    const DB = {json_str};

    function switchTab(tabId) {{
      document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
      document.querySelectorAll('.tab-pane').forEach(pane => pane.classList.remove('active'));
      
      const activeBtn = Array.from(document.querySelectorAll('.tab-btn')).find(b => b.getAttribute('onclick').includes(tabId));
      if (activeBtn) activeBtn.classList.add('active');
      
      const activePane = document.getElementById(tabId);
      if (activePane) activePane.classList.add('active');

      if (tabId === 'tab-volcano') updateVolcanoPlot();
      if (tabId === 'tab-clock') updateGauge();
      if (tabId === 'tab-pathways') renderPathways();
    }}

    function updateVolcanoPlot() {{
      const cohort = document.getElementById('cohortSelect').value;
      const fcThresh = parseFloat(document.getElementById('fcSlider').value);
      const pvalThresh = parseFloat(document.getElementById('pvalSlider').value);
      const searchTerm = document.getElementById('geneSearch').value.trim().toUpperCase();

      document.getElementById('fcVal').innerText = fcThresh.toFixed(2);
      document.getElementById('pvalVal').innerText = pvalThresh.toFixed(1);

      const dataset = cohort === 'age' ? DB.ageGenes : DB.photoGenes;

      const x = [];
      const y = [];
      const text = [];
      const colors = [];
      const filteredGenes = [];

      dataset.forEach(g => {{
        const isUp = g.log2FC_per_decade >= fcThresh && g.neg_log10_pval >= pvalThresh;
        const isDown = g.log2FC_per_decade <= -fcThresh && g.neg_log10_pval >= pvalThresh;
        const isMatch = searchTerm ? g.symbol.includes(searchTerm) : true;

        let col = '#475569';
        if (isUp) col = '#f43f5e';
        else if (isDown) col = '#38bdf8';

        if (searchTerm && isMatch) {{
          col = '#eab308'; // highlight yellow
        }}

        x.push(g.log2FC_per_decade);
        y.push(g.neg_log10_pval);
        text.push(`<b>${{g.symbol}}</b><br>Log2 FC: ${{g.log2FC_per_decade}}<br>-Log10 P-val: ${{g.neg_log10_pval}}<br>FDR: ${{g.fdr_qval}}`);
        colors.push(col);

        if ((isUp || isDown) && isMatch) {{
          filteredGenes.push(g);
        }}
      }});

      document.getElementById('filterCount').innerText = filteredGenes.length;

      const trace = {{
        x: x,
        y: y,
        text: text,
        mode: 'markers',
        type: 'scatter',
        hoverinfo: 'text',
        marker: {{
          color: colors,
          size: 7,
          opacity: 0.8
        }}
      }};

      const layout = {{
        title: cohort === 'age' ? 'Chronological Ageing Volcano Plot (GSE226189)' : 'Photoaging Volcano Plot (GSE38308: Pre vs Post-Auricular)',
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: {{ color: '#cbd5e1', family: 'Plus Jakarta Sans' }},
        xaxis: {{
          title: cohort === 'age' ? 'Slope / Effect per Decade (Log2 FC)' : 'Log2 Fold Change (Sun-Exposed vs Protected)',
          gridcolor: 'rgba(255,255,255,0.06)',
          zerolinecolor: 'rgba(255,255,255,0.2)'
        }},
        yaxis: {{
          title: '-Log10(P-value)',
          gridcolor: 'rgba(255,255,255,0.06)',
          zerolinecolor: 'rgba(255,255,255,0.2)'
        }},
        shapes: [
          {{ type: 'line', x0: -fcThresh, x1: -fcThresh, y0: 0, y1: Math.max(...y), line: {{ color: 'rgba(56,189,248,0.4)', dash: 'dash' }} }},
          {{ type: 'line', x0: fcThresh, x1: fcThresh, y0: 0, y1: Math.max(...y), line: {{ color: 'rgba(244,63,94,0.4)', dash: 'dash' }} }},
          {{ type: 'line', x0: Math.min(...x), x1: Math.max(...x), y0: pvalThresh, y1: pvalThresh, line: {{ color: 'rgba(255,255,255,0.3)', dash: 'dot' }} }}
        ],
        margin: {{ t: 40, r: 20, b: 50, l: 50 }}
      }};

      Plotly.react('volcanoPlot', [trace], layout, {{ responsive: true }});

      // Update Table
      const tbody = document.getElementById('geneTableBody');
      tbody.innerHTML = '';
      filteredGenes.slice(0, 20).forEach(g => {{
        const row = document.createElement('tr');
        const direction = g.log2FC_per_decade > 0 ? 
          '<span style="color: #f43f5e; font-weight:700;">▲ Upregulated</span>' : 
          '<span style="color: #38bdf8; font-weight:700;">▼ Downregulated</span>';
        row.innerHTML = `
          <td style="color: #ffffff; font-weight: 700;">${{g.symbol}}</td>
          <td>${{g.log2FC_per_decade}}</td>
          <td>${{g.mean_expr}}</td>
          <td>${{g.t_stat}}</td>
          <td>${{g.p_value.toExponential(2)}}</td>
          <td>${{g.fdr_qval.toExponential(2)}}</td>
          <td>${{direction}}</td>
        `;
        tbody.appendChild(row);
      }});
    }}

    function initPatientOptions() {{
      const select = document.getElementById('patientSelect');
      select.innerHTML = '';
      DB.patients.forEach((p, idx) => {{
        const opt = document.createElement('option');
        opt.value = idx;
        opt.innerText = `${{p.sample_id}} — Age: ${{p.age}}y (${{p.sex}})`;
        select.appendChild(opt);
      }});
      updateGauge();
    }}

    function updateGauge() {{
      const select = document.getElementById('patientSelect');
      const idx = parseInt(select.value) || 0;
      const patient = DB.patients[idx] || DB.patients[0];

      const actualAge = patient.age;
      const predAge = patient.predicted_age_cv;
      const delta = predAge - actualAge;

      const data = [
        {{
          type: "indicator",
          mode: "gauge+number+delta",
          value: predAge,
          delta: {{ reference: actualAge, position: "top", valueformat: "+.1f yrs" }},
          title: {{ text: `Predicted Age for ${{patient.sample_id}}`, font: {{ size: 20, color: '#f8fafc' }} }},
          gauge: {{
            axis: {{ range: [10, 100], tickcolor: "#94a3b8" }},
            bar: {{ color: "#38bdf8" }},
            bgcolor: "rgba(255,255,255,0.05)",
            borderwidth: 1,
            bordercolor: "rgba(255,255,255,0.1)",
            steps: [
              {{ range: [10, 40], color: "rgba(16, 185, 129, 0.15)" }},
              {{ range: [40, 70], color: "rgba(56, 189, 248, 0.15)" }},
              {{ range: [70, 100], color: "rgba(244, 63, 94, 0.15)" }}
            ],
            threshold: {{
              line: {{ color: "#f43f5e", width: 4 }},
              thickness: 0.8,
              value: actualAge
            }}
          }}
        }}
      ];

      const layout = {{
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: {{ color: '#cbd5e1', family: 'Plus Jakarta Sans' }},
        margin: {{ t: 40, r: 30, b: 30, l: 30 }}
      }};

      Plotly.react('gaugePlot', data, layout, {{ responsive: true }});

      // Details
      document.getElementById('patientDetails').innerHTML = `
        <div style="display: flex; justify-content: space-around; text-align: center;">
          <div><span style="color:var(--muted)">Donor Sex:</span> <b style="color:#fff;">${{patient.sex}}</b></div>
          <div><span style="color:var(--muted)">Chronological Age:</span> <b style="color:#f43f5e;">${{actualAge}} years</b></div>
          <div><span style="color:var(--muted)">Transcriptomic Age:</span> <b style="color:var(--cyan);">${{predAge.toFixed(1)}} years</b></div>
          <div><span style="color:var(--muted)">Delta Error:</span> <b style="color:${{Math.abs(delta) < 5 ? 'var(--emerald)' : '#cbd5e1'}};">${{delta >= 0 ? '+' : ''}}${{delta.toFixed(1)}}y</b></div>
        </div>
      `;

      // Render top weights
      const coefTbody = document.getElementById('coefTableBody');
      coefTbody.innerHTML = '';
      DB.coefficients.forEach(c => {{
        const r = document.createElement('tr');
        const col = c.coefficient > 0 ? 'var(--rose)' : 'var(--cyan)';
        r.innerHTML = `
          <td style="color:#fff; font-weight:700;">${{c.symbol}}</td>
          <td style="color:${{col}}; font-weight:700;">${{c.coefficient > 0 ? '+' : ''}}${{c.coefficient}}</td>
        `;
        coefTbody.appendChild(r);
      }});
    }}

    function renderPathways() {{
      const terms = DB.pathways.map(p => p.Term);
      const nes = DB.pathways.map(p => p.NES);
      const colors = nes.map(n => n > 0 ? '#f43f5e' : '#38bdf8');

      const trace = {{
        x: nes,
        y: terms,
        type: 'bar',
        orientation: 'h',
        marker: {{ color: colors }}
      }};

      const layout = {{
        title: 'Top MSigDB Hallmark Enriched Pathways in Chronological Ageing',
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        font: {{ color: '#cbd5e1', family: 'Plus Jakarta Sans' }},
        xaxis: {{
          title: 'Normalized Enrichment Score (NES)',
          gridcolor: 'rgba(255,255,255,0.06)',
          zerolinecolor: 'rgba(255,255,255,0.3)'
        }},
        yaxis: {{
          autorange: 'reversed',
          tickfont: {{ size: 12 }}
        }},
        margin: {{ t: 40, r: 20, b: 50, l: 240 }}
      }};

      Plotly.react('pathwayPlot', [trace], layout, {{ responsive: true }});

      // Cards
      const container = document.getElementById('pathwayCards');
      container.innerHTML = '';
      DB.explanations.forEach(e => {{
        const card = document.createElement('div');
        card.className = 'pathway-card';
        card.style.borderLeftColor = e.NES > 0 ? 'var(--rose)' : 'var(--cyan)';
        card.innerHTML = `
          <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;">
            <h4 style="font-size: 18px; color: #fff;">${{e.Pathway}}</h4>
            <span style="font-family:'JetBrains Mono'; font-weight:700; color:${{e.NES > 0 ? 'var(--rose)' : 'var(--cyan)'}}">
              NES: ${{e.NES.toFixed(2)}}
            </span>
          </div>
          <div style="font-size: 12px; color: var(--muted); margin-bottom: 8px;">Direction: <b>${{e.Biological_Direction}}</b> (FDR: ${{e.FDR_qval.toExponential(2)}})</div>
          <p style="font-size: 14px; color: #cbd5e1; line-height: 1.6; margin-bottom: 12px;">${{e.Plain_Language_Summary}}</p>
          <div style="font-size: 12px; color: #64748b;">${{e.Literature_Lookup_Recommendation}}</div>
        `;
        container.appendChild(card);
      }});
    }}

    // Initialise on load
    window.addEventListener('DOMContentLoaded', () => {{
      initPatientOptions();
      updateVolcanoPlot();
      renderPathways();
    }});
  </script>
</body>
</html>
"""

    with open('docs/index.html', 'w') as f:
        f.write(html_template)

    with open('demo.html', 'w') as f:
        f.write(html_template)

    print(f"Generated docs/index.html and demo.html successfully! File size: {os.path.getsize('demo.html') / 1024:.1f} KB")

if __name__ == '__main__':
    build_demo_html()
