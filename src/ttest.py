import streamlit as st
import pandas as pd
import pingouin as pg
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import time

def gen_ttest_data(ttest_attribute, target_groups, paired, alternative, correction, p_correction, _progress_callback=None):
    df = pd.concat([st.session_state.data, st.session_state.md], axis=1)
    ttest = []
    columns = list(st.session_state.data.columns)
    total = len(columns)
    start_time = None
    for idx, col in enumerate(columns):
        if idx == 0:
            start_time = time.time()
        
        group1 = df[col][df[ttest_attribute] == target_groups[0]]
        group2 = df[col][df[ttest_attribute] == target_groups[1]]
        
        # Calculate means for volcano plot
        mean1 = group1.mean()
        mean2 = group2.mean()

        # Determine which t-test to use
        # The 'correction' variable is "auto", "True", or "False" (strings)
        correction_param = correction
        if correction == "True":
            correction_param = True
        elif correction == "False":
            correction_param = False
        
        # Calculate t-test
        result = pg.ttest(group1, group2, paired, alternative, correction=correction_param)

        # Determine which test was *actually* used for DF calculation
        was_welch = False
        if correction_param == True:
            was_welch = True
        elif correction_param == 'auto':
            was_welch = (len(group1) != len(group2))
        # else: was_welch = False (already set)

        # Calculate degrees of freedom
        if was_welch:
            # Welch's df calculation
            s1 = np.var(group1, ddof=1)
            s2 = np.var(group2, ddof=1)
            n1 = len(group1)
            n2 = len(group2)
            numerator = (s1/n1 + s2/n2)**2
            denominator = ((s1/n1)**2)/(n1-1) + ((s2/n2)**2)/(n2-1)
            df_welch = numerator / denominator if denominator != 0 else np.nan
            result["df"] = df_welch
            result["ttest_type"] = "Welch"
        else:
            # Student's t-test df
            n1 = len(group1)
            n2 = len(group2)
            # For paired test, df is n-1
            if paired:
                result["df"] = n1 - 1
                result["ttest_type"] = "Paired Student"
            else:
                result["df"] = n1 + n2 - 2
                result["ttest_type"] = "Student"
        
        result["metabolite"] = col
        result["mean(A)"] = mean1
        result["mean(B)"] = mean2
        ttest.append(result)
        
        # Progress callback with estimated time left
        if _progress_callback is not None:
            elapsed = time.time() - start_time if start_time else 0
            done = idx + 1
            est_left = (elapsed / done) * (total - done) if done > 0 else 0
            _progress_callback(done, total, est_left)

    ttest = pd.concat(ttest).set_index("metabolite")
    ttest = ttest.dropna(subset=['p-val']) # Only drop if p-val is NaN

    ttest.insert(8, "p-corrected", pg.multicomp(ttest["p-val"].astype(float), method=p_correction)[1])
    # add significance
    ttest.insert(9, "significance", ttest["p-corrected"] < 0.05)
    ttest.insert(10, "st.session_state.ttest_attribute", ttest_attribute)
    ttest.insert(11, "A", target_groups[0])
    ttest.insert(12, "B", target_groups[1])

    return ttest.sort_values("p-corrected")


@st.cache_resource
def plot_ttest(df):

    # Use the correct t-statistic column name from pingouin output
    t_col = None
    for candidate in ["T", "T-val", "t", "tval"]:
        if candidate in df.columns:
            t_col = candidate
            break
    if t_col is None:
        st.error("No t-statistic column found in t-test results. Columns: " + str(df.columns))
        return go.Figure()


    # Add a column for -log(p-corrected) and significance label
    df = df.copy()
    df["-log_p_corrected"] = df["p-corrected"].apply(lambda x: -np.log(x + 1e-300)) # Add epsilon
    df["sig_label"] = df["significance"].apply(lambda x: "significant" if x else "insignificant")
    df["metabolite_name"] = df.index


    if df.empty:
        st.warning("No t-test results to display. Please check your data or selection.")
        return go.Figure()

    fig = px.scatter(
        df,
        x=t_col,
        y="-log_p_corrected",
        color="sig_label",
        color_discrete_map={"significant": "#ef553b", "insignificant": "#696880"},
        custom_data=["metabolite_name"],
        template="plotly_white",
        width=600,
        height=600,
    )

    # Custom hovertemplate: metabolite&name: <metabolite name>
    fig.update_traces(hovertemplate="metabolite&name: %{customdata[0]}<extra></extra>")

    xlim = [df[t_col].min(), df[t_col].max()]
    x_padding = abs(xlim[1] - xlim[0]) / 5 if (xlim[1] != xlim[0] and pd.notnull(xlim[0]) and pd.notnull(xlim[1])) else 1
    fig.update_layout(xaxis=dict(range=[xlim[0] - x_padding, xlim[1] + x_padding]))

    # Defensive: Only set title if df has at least one row
    if len(df) > 0:
        title_text = f"t-test - FEATURE SIGNIFICANCE - {str(df.iloc[0, 10]).upper()}: {df.iloc[0, 11]} - {df.iloc[0, 12]}"
    else:
        title_text = "t-test - FEATURE SIGNIFICANCE"
    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": title_text,
            "font_color": "#3E3D53",
        },
        xaxis_title="T-statistic",
        yaxis_title="-Log(p-corrected)",
        showlegend=True,  # Enable legend
        legend_title_text="Significance",
    )
    return fig

def _get_ttest_feature_map():
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

@st.cache_resource(show_spinner="Creating volcano plot...")
def get_ttest_volcano_plot(df):
    """Volcano plot for t-test: x = log2 fold change (mean(B)/mean(A)), y = -log10(p-corrected).
    Adds metabolite/feature hover labels.
    """
    feature_map = _get_ttest_feature_map()
    df = df.copy()
    eps = 1e-9
    
    # Ensure means are numeric, handle potential NaNs
    meanA = pd.to_numeric(df["mean(A)"], errors='coerce').fillna(0) + eps
    meanB = pd.to_numeric(df["mean(B)"], errors='coerce').fillna(0) + eps
    p_values = pd.to_numeric(df["p-corrected"], errors='coerce').fillna(1.0) + eps

    df["log2FC"] = np.log2(meanB / meanA)
    df["neglog10p"] = -np.log10(p_values)

    fig = go.Figure()

    def make_hovertext(metabolites):
        htext = []
        for m in metabolites:
            met_name = feature_map.get(m, str(m)) if feature_map else str(m)
            htext.append(f"metabolite&name: {met_name}")
        return htext

    ins = df[df["significance"] == False]
    if not ins.empty:
        fig.add_trace(
            go.Scatter(
                x=ins["log2FC"],
                y=ins["neglog10p"],
                mode="markers",
                marker=dict(color="#696880"),
                name="insignificant",
                hovertext=make_hovertext(ins.index),
                hoverinfo="text",
            )
        )

    sig = df[df["significance"] == True]
    if not sig.empty:
        # Group A blue (log2FC < 0)
        sig_A = sig[sig["log2FC"] < 0]
        if not sig_A.empty:
            fig.add_trace(
                go.Scatter(
                    x=sig_A["log2FC"],
                    y=sig_A["neglog10p"],
                    mode="markers",
                    marker=dict(color="#1f77b4"), 
                    name=f"Significant: {st.session_state.ttest_options[0]} > {st.session_state.ttest_options[1]}",
                    hovertext=make_hovertext(sig_A.index),
                    hoverinfo="text",
                )
            )
        # Group B red (log2FC > 0)
        sig_B = sig[sig["log2FC"] > 0]
        if not sig_B.empty:
            fig.add_trace(
                go.Scatter(
                    x=sig_B["log2FC"],
                    y=sig_B["neglog10p"],
                    mode="markers",
                    marker=dict(color="#ef553b"),
                    name=f"Significant: {st.session_state.ttest_options[1]} > {st.session_state.ttest_options[0]}",
                    hovertext=make_hovertext(sig_B.index),
                    hoverinfo="text",
                )
            )

    fig.update_layout(
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": f"t-test - VOLCANO PLOT - {st.session_state.ttest_attribute.upper()}: {st.session_state.ttest_options[0]} - {st.session_state.ttest_options[1]}",
            "font_color": "#3E3D53",
        },
        xaxis_title="log2(mean B / mean A)",
        yaxis_title="-log10(p-corrected)",
        template="plotly_white",
        legend=dict(
            title="Legend",
            itemsizing='trace'
        ),
        width=700,
        height=600,
    )
    return fig


@st.cache_resource
def ttest_boxplot(df_ttest, metabolite):
    df = pd.concat([st.session_state.md, st.session_state.data], axis=1)
    
    # Filter for the two selected groups
    df_filtered = df[df[st.session_state.ttest_attribute].isin(st.session_state.ttest_options)]
    
    # Create the plot from the filtered dataframe
    fig = px.box(
        df_filtered,
        x=st.session_state.ttest_attribute,
        y=metabolite,
        color=st.session_state.ttest_attribute,
        width=350,
        height=400,
        points="all",
        custom_data=[st.session_state.ttest_attribute, metabolite],
    )

    # Set custom hovertemplate
    attribute = st.session_state.ttest_attribute if hasattr(st, "session_state") and hasattr(st.session_state, "ttest_attribute") else "Attribute"
    hovertemplate = (
        f"attribute: {attribute}<br>"
        "option: %{customdata[0]}<br>"
        f"metabolite&name: {metabolite}<br>"
        "intensity: %{customdata[1]:.3g}<extra></extra>"
    )
    fig.update_traces(hovertemplate=hovertemplate)
    
    # Determine significance label
    is_significant = False
    pvalue = 1.0
    if metabolite in df_ttest.index:
        is_significant = df_ttest.loc[metabolite, "significance"] if "significance" in df_ttest.columns else False
        pvalue = df_ttest.loc[metabolite, "p-corrected"]
        
    sig_label = "Significant Metabolite:" if is_significant else "Insignificant Metabolite:"
    
    fig.update_layout(
        showlegend=False,
        xaxis_title=st.session_state.ttest_attribute.replace("st.session_state.ttest_attribute_", ""),
        yaxis_title="intensity",
        template="plotly_white",
        font={"color": "grey", "size": 12, "family": "Sans"},
        title={
            "text": f"{sig_label} {metabolite}",
            "font_color": "#3E3D53",
        },
    )
    fig.update_yaxes(title_standoff=10)
    
    if pvalue >= 0.05:
        symbol = "ns"
    elif pvalue >= 0.01:
        symbol = "*"
    elif pvalue >= 0.001:
        symbol = "**"
    else:
        symbol = "***"

    # Get y-max from the filtered data
    y_max = df_filtered[metabolite].max()
    top_y = y_max * 1.2 if pd.notnull(y_max) else 1

    # Define x-coordinates for the line/annotation
    x0 = st.session_state.ttest_options[0]
    x1 = st.session_state.ttest_options[1]

    # horizontal line
    fig.add_shape(
        type="line",
        x0=x0,
        y0=top_y,
        x1=x1,
        y1=top_y,
        line=dict(width=1, color="#000000"),
    )
    
    if symbol == "ns":
        y_margin = y_max * 0.05 if pd.notnull(y_max) else 0.05
    else:
        y_margin = y_max * 0.1 if pd.notnull(y_max) else 0.1
        
    # Center the annotation
    # This works for categorical x-axis
    fig.add_annotation(
        x=0.5,  # Center point for 2 categories
        xref="paper", # Use paper reference for x
        y=top_y + y_margin,
        text=f"<b>{symbol}</b>",
        showarrow=False,
        font_color="#555555",
    )
    return fig
