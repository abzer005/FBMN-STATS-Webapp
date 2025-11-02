import streamlit as st
import pandas as pd
import numpy as np
import pingouin as pg
import plotly.express as px
import plotly.graph_objects as go
from scipy.stats import kruskal
import scikit_posthocs as sp

def gen_kruskal_data(group_data, _progress_callback=None):
    total = len(group_data[0].columns)
    for idx, col in enumerate(group_data[0].columns):
        try:
            statistic, p = kruskal(*[df[col] for df in group_data])
            if _progress_callback is not None:
                _progress_callback(idx + 1, total, max(0, total - (idx + 1)))
            yield col, p, statistic
        except ValueError:
            continue

def add_p_correction_to_kruskal(df, correction):
    # add Bonferroni corrected p-values for multiple testing correction
    if "p-corrected" not in df.columns:
        df.insert(2, "p-corrected",
                  pg.multicomp(df["p"].astype(float), method=correction)[1])
    # add significance
    if "significant" not in df.columns:
        df.insert(3, "significant", df["p-corrected"] < 0.05)
    # sort by p-value
    df.sort_values("p", inplace=True)
    return df

def kruskal_wallis(df, attribute, correction, elements, _progress_callback=None):
    combined = pd.concat([df, st.session_state.md], axis=1)
    if elements is not None:
        combined = combined[combined[attribute].isin(elements)]
    groups = combined[attribute].unique()
    metabolite_cols = list(st.session_state.data.columns)
    group_data = [combined[combined[attribute] == group].loc[:, metabolite_cols] for group in groups]

    df = pd.DataFrame(
        np.fromiter(
            gen_kruskal_data(group_data, _progress_callback=_progress_callback),
            dtype=[("metabolite", "U100"), ("p", "f"), ("statistic", "f")],
        )
    )
    df = df.dropna()
    df = add_p_correction_to_kruskal(df, correction)
    df = df[df["metabolite"] != attribute]
    return df

def _get_feature_name_map():
    """Return a mapping metabolite_id -> feature name. Look for common name columns
    in st.session_state.ft_gnps. If nothing found, return None."""
    ft = st.session_state.get("ft_gnps", pd.DataFrame())
    if ft is None or ft.empty:
        return None
    candidates = ["metabolite_name", "name", "feature_name", "compound_name", "compound"]
    for c in candidates:
        if c in ft.columns:
            return ft[c].to_dict()
    return None

@st.cache_resource
def get_kruskal_plot(kruskal):

    kruskal_clean = kruskal[kruskal["metabolite"].notna()].copy()
    if "kruskal_attribute" in st.session_state:
        kr_attr = st.session_state.kruskal_attribute
        kruskal_clean = kruskal_clean[kruskal_clean["metabolite"] != kr_attr]

    if "significant" in kruskal_clean.columns:
        kruskal_clean["significant"] = kruskal_clean["significant"].fillna(False).astype(bool)
    else:
        kruskal_clean["significant"] = False

    total_points = len(kruskal_clean)
    n_significant = int(kruskal_clean["significant"].sum())
    n_insignificant = total_points - n_significant

    st.write(f"Significant: {n_significant}")
    st.write(f"Insignificant: {n_insignificant}")
    st.write(f"Total data points: {total_points}")

    insig = kruskal_clean[~kruskal_clean["significant"]]
    sig = kruskal_clean[kruskal_clean["significant"]]

    fig = go.Figure()
    eps = 1e-12
    def safe_log10_series(s):
        s = pd.to_numeric(s, errors="coerce")
        return s.where(s > 0, np.nan).apply(np.log10)

    def safe_neglog10p(pseries):
        p = pd.to_numeric(pseries, errors="coerce").fillna(1.0)
        p = p.clip(lower=eps)
        return -np.log10(p)

    fig.add_trace(go.Scatter(
        x=safe_log10_series(insig["statistic"]),
        y=safe_neglog10p(insig["p"]),
        mode="markers",
        marker=dict(color="#696880"),
        name="insignificant",
        text=insig["metabolite"],
        hovertemplate="%{text}",
        showlegend=True
    ))
    fig.add_trace(go.Scatter(
        x=safe_log10_series(sig["statistic"]),
        y=safe_neglog10p(sig["p"]),
        mode="markers",
        marker=dict(color="#ef553b"),
        name="significant",
        text=sig["metabolite"],
        hovertemplate="%{text}",
        showlegend=True
    ))

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": f"Kruskal Wallis - {st.session_state.kruskal_attribute.upper()}",
            "font_color": "#3E3D53"
        },
        xaxis_title="log10(H)",
        yaxis_title="-log10(p)",
        legend=dict(title="Legend"),
        width=600,
        height=600
    )
    return fig

@st.cache_resource
def get_metabolite_boxplot(kruskal, metabolite):
    attribute = st.session_state.kruskal_attribute
    p_value = kruskal.set_index("metabolite")._get_value(metabolite, "p-corrected")
    df = pd.concat([st.session_state.data, st.session_state.md], axis=1)[[attribute, metabolite]].copy()

    # Add filename column if available
    df = df.reset_index().rename(columns={"index": "filename"})
    if df.columns[0] == "filename" and st.session_state.data.index.name:
        df.rename(columns={"filename": st.session_state.data.index.name}, inplace=True)

    # Get metabolite name from feature map if available
    def _get_feature_name_map():
        ft = st.session_state.get("ft_gnps", pd.DataFrame())
        if ft is None or ft.empty:
            return None
        candidates = ["metabolite_name", "name", "feature_name", "compound_name", "compound"]
        for c in candidates:
            if c in ft.columns:
                return ft[c].to_dict()
        return None
    feature_map = _get_feature_name_map()
    metabolite_name = feature_map.get(metabolite, metabolite) if feature_map else metabolite
    df["metabolite_name"] = metabolite_name

    df["intensity"] = df[metabolite]
    df["hovertext"] = df.apply(lambda row: f"filename: {row['filename']}<br>attribute&group: {attribute}, {row[attribute]}<br>metabolite: {row['metabolite_name']}<br>intensity: {row['intensity']}", axis=1)

    try:
        p_value_float = float(p_value)
        p_value_str = f"{p_value_float:.2e}"
    except Exception:
        p_value_float = None
        p_value_str = str(p_value)

    # Determine significance
    is_significant = p_value_float is not None and p_value_float < 0.05
    significance_text = "Significant" if is_significant else "Insignificant"
    title = f"{significance_text} Metabolite: {metabolite_name}<br>Corrected p-value: {p_value_str}"
    fig = px.box(
        df,
        x=attribute,
        y=metabolite,
        template="plotly_white",
        width=800,
        height=600,
        points="all",
        color=attribute,
        hover_data=None,
        custom_data=[df["hovertext"]],
    )

    fig.update_traces(hovertemplate="%{customdata[0]}")

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={"text": title, "font_color": "#3E3D53"},
        xaxis_title=attribute,
        yaxis_title="intensity",
    )
    return fig


def gen_pairwise_dunn(group_data, _progress_callback=None):
    """Yield results for pairwise dunn test for all metabolites between two options within the attribute."""
    total = len(group_data[0].columns)
    for idx, col in enumerate(group_data[0].columns):
        p = sp.posthoc_dunn([df[col] for df in group_data]).iloc[0, 1]
        if _progress_callback is not None:
            _progress_callback(idx + 1, total, max(0, total - (idx + 1)))
        yield (col, p)


def add_p_value_correction_to_dunns(dunn, correction):
    if "p-corrected" not in dunn.columns:
        # add Bonferroni corrected p-values
        dunn.insert(
            2, "p-corrected", pg.multicomp(
                dunn["p"].astype(float), method=correction)[1]
        )
        # add significance
        dunn.insert(3, "stats_significant", dunn["p-corrected"] < 0.05)
        # sort by p-value
        dunn.sort_values("p", inplace=True)
    return dunn

def dunn(df, attribute, elements, correction, _progress_callback=None):
    # Use all metabolites from KW results, not just significant ones
    all_metabolites = df["metabolite"]
    valid_metabolites = [m for m in all_metabolites if m in st.session_state.data.columns]
    data = pd.concat(
        [
            st.session_state.data.loc[:, valid_metabolites],
            st.session_state.md[attribute],
        ],
        axis=1,
    )
    data = data[data[attribute].isin(elements)]
    # Calculate means for each group
    groupA, groupB = elements[0], elements[1]
    meanA = data[data[attribute] == groupA][valid_metabolites].mean()
    meanB = data[data[attribute] == groupB][valid_metabolites].mean()

    dunn = pd.DataFrame(
        np.fromiter(
            gen_pairwise_dunn([data[data[attribute] == element].drop(columns=[attribute]) for element in elements], _progress_callback=_progress_callback),
            dtype=[
                ("stats_metabolite", "U100"),
                ("p", "f")
            ],
        )
    )
    dunn = dunn.dropna()
    dunn = add_p_value_correction_to_dunns(dunn, correction)
    # Add mean(A) and mean(B) columns for volcano plot
    dunn["mean(A)"] = dunn["stats_metabolite"].map(meanA)
    dunn["mean(B)"] = dunn["stats_metabolite"].map(meanB)
    dunn["diff"] = dunn["mean(B)"] - dunn["mean(A)"]
    
    # Store the original numeric DataFrame before formatting
    dunn_original = dunn.copy()
    
    # Create display version with formatted p-values
    dunn_display = dunn.copy()
    if "p" in dunn_display.columns:
        dunn_display["p"] = dunn_display["p"].apply(lambda x: f"{x:.2e}" if pd.notnull(x) else x)
    if "p-corrected" in dunn_display.columns:
        dunn_display["p-corrected"] = dunn_display["p-corrected"].apply(lambda x: f"{x:.2e}" if pd.notnull(x) else x)
    
    # Attach the original numeric DataFrame as an attribute
    dunn_display._original = dunn_original
    
    return dunn_display

def _get_dunn_feature_map(df_dunn):
    """Look up feature name mapping for stats_metabolite similar to kruskal."""
    return _get_feature_name_map()

@st.cache_resource
def get_dunn_teststat_plot(df):
    feature_map = _get_dunn_feature_map(df)
    fig = go.Figure()

    def make_hovertext(metabolites):
        htext = []
        for m in metabolites:
            met_name = feature_map[m] if feature_map and m in feature_map else str(m)
            htext.append(f"metabolite: {met_name}")
        return htext

    # Use original DataFrame if available (contains numeric p-values)
    df_numeric = getattr(df, '_original', df)
    
    sig_col = "stats_significant" if "stats_significant" in df_numeric.columns else "significant"
    met_col = "stats_metabolite" if "stats_metabolite" in df_numeric.columns else "metabolite"
    diff_col = "diff" if "diff" in df_numeric.columns else None

    if diff_col is None:
        st.error("Diff column not found in Dunn's results. Please rerun the analysis.")
        return fig

    # Insignificant points
    ins = df_numeric[df_numeric[sig_col] == False]
    if not ins.empty:
        fig.add_trace(
            go.Scatter(
                x=ins[diff_col],
                y=-np.log(ins["p"].astype(float)),
                mode="markers",
                marker=dict(color="#696880"),
                name="insignificant",
                hovertext=make_hovertext(ins[met_col]),
                hoverinfo="text",
            )
        )

    # Significant points
    sig = df_numeric[df_numeric[sig_col] == True]
    if not sig.empty:
        fig.add_trace(
            go.Scatter(
                x=sig[diff_col],
                y=-np.log(sig["p"].astype(float)),
                mode="markers+text",
                marker=dict(color="#ef553b"),
                text=["" for _ in sig[met_col]],
                textposition="top right",
                textfont=dict(color="#ef553b", size=12),
                name="significant",
                hovertext=make_hovertext(sig[met_col]),
                hoverinfo="text",
            )
        )

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": f"DUNN - {st.session_state.kruskal_attribute.upper()}: {st.session_state.dunn_elements[0]} - {st.session_state.dunn_elements[1]} (test-statistic)",
            "font_color": "#3E3D53",
        },
        xaxis_title="diff (mean B - mean A)",
        yaxis_title="-log(p)",
        template="plotly_white",
        legend=dict(
            title="Legend",
            itemsizing='trace',
            font=dict(size=12),
            orientation="v",
            x=1.02,
            y=1,
            xanchor="left",
            yanchor="top"
        ),
        width=600,
        height=600,
    )
    fig.update_yaxes(title_standoff=10)
    return fig

@st.cache_resource
def get_dunn_volcano_plot(df):
    """Volcano plot for Dunn's: x = log2 fold change (mean(B)/mean(A)), y = -log10(p-value).
    Adds metabolite/feature hover labels.
    """
    feature_map = _get_dunn_feature_map(df)
    
    # Use original DataFrame if available (contains numeric values)
    df_numeric = getattr(df, '_original', df).copy()
    
    # compute log2 fold change (B relative to A). avoiding zeros by using a small epsilon
    eps = 1e-9
    meanA = df_numeric["mean(A)"].astype(float) + eps
    meanB = df_numeric["mean(B)"].astype(float) + eps
    
    df_numeric["log2FC"] = np.log2(meanB / meanA)
    df_numeric["neglog10p"] = -np.log10(df_numeric["p"].astype(float) + eps)

    fig = go.Figure()

    def make_hovertext(metabolites):
        htext = []
        for m in metabolites:
            met_name = feature_map[m] if feature_map and m in feature_map else str(m)
            htext.append(f"metabolite: {met_name}")
        return htext

    sig_col = "stats_significant" if "stats_significant" in df_numeric.columns else "significant"
    met_col = "stats_metabolite" if "stats_metabolite" in df_numeric.columns else "metabolite"
    
    ins = df_numeric[df_numeric[sig_col] == False]
    if not ins.empty:
        fig.add_trace(
            go.Scatter(
                x=ins["log2FC"],
                y=ins["neglog10p"],
                mode="markers",
                marker=dict(color="#696880"),
                name="insignificant",
                hovertext=make_hovertext(ins[met_col]),
                hoverinfo="text",
                showlegend=True,
            )
        )

    sig = df_numeric[df_numeric[sig_col] == True]
    if not sig.empty:
        # Group A blue (higher in A, so negative log2FC)
        sig_A = sig[sig["log2FC"] < 0]
        if not sig_A.empty:
            fig.add_trace(
                go.Scatter(
                    x=sig_A["log2FC"],
                    y=sig_A["neglog10p"],
                    mode="markers+text",
                    marker=dict(color="#1f77b4"), 
                    text=["" for _ in sig_A[met_col]],
                    textposition="top right",
                    textfont=dict(color="#1f77b4", size=12),
                    name=f"Significant: {st.session_state.dunn_elements[0]} > {st.session_state.dunn_elements[1]}",
                    hovertext=make_hovertext(sig_A[met_col]),
                    hoverinfo="text",
                    showlegend=True,
                )
            )
        # Group B red (higher in B, so positive log2FC)
        sig_B = sig[sig["log2FC"] > 0]
        if not sig_B.empty:
            fig.add_trace(
                go.Scatter(
                    x=sig_B["log2FC"],
                    y=sig_B["neglog10p"],
                    mode="markers+text",
                    marker=dict(color="#ef553b"),  #red
                    text=["" for _ in sig_B[met_col]],
                    textposition="top right",
                    textfont=dict(color="#ef553b", size=12),
                    name=f"Significant: {st.session_state.dunn_elements[1]} > {st.session_state.dunn_elements[0]}",
                    hovertext=make_hovertext(sig_B[met_col]),
                    hoverinfo="text",
                    showlegend=True,
                )
            )

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": f"DUNN - {st.session_state.kruskal_attribute.upper()}: {st.session_state.dunn_elements[0]} - {st.session_state.dunn_elements[1]} (volcano)",
            "font_color": "#3E3D53",
        },
        xaxis_title="log2(mean B) - log2(mean A)",
        yaxis_title="-log10(p)",
        template="plotly_white",
        legend=dict(
            title="Legend",
            itemsizing='trace',
            font=dict(size=12),
            orientation="v",
            x=1.02,
            y=1,
            xanchor="left",
            yanchor="top"
        ),
        width=700,
        height=600,
    )
    return fig