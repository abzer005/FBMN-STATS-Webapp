import streamlit as st
import pandas as pd
import pingouin as pg
import plotly.express as px
import numpy as np

@st.cache_data(show_spinner="Calculating t-test results...")
def gen_ttest_data(ttest_attribute, target_groups, paired, alternative, correction, p_correction, _progress_callback=None):
    df = pd.concat([st.session_state.data, st.session_state.md], axis=1)
    ttest = []
    columns = list(st.session_state.data.columns)
    total = len(columns)
    start_time = None
    for idx, col in enumerate(columns):
        if idx == 0:
            import time
            start_time = time.time()
        group1 = df[col][df[ttest_attribute] == target_groups[0]]
        group2 = df[col][df[ttest_attribute] == target_groups[1]]
        # Determine which t-test to use
        use_welch = False
        if correction == "auto":
            use_welch = (len(group1) != len(group2))
        elif correction == "True":
            use_welch = True
        elif correction == "False":
            use_welch = False
        # Calculate t-test
        result = pg.ttest(group1, group2, paired, alternative, use_welch)
        # Calculate degrees of freedom
        if use_welch:
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
            result["df"] = n1 + n2 - 2
            result["ttest_type"] = "Student"
        result["metabolite"] = col
        ttest.append(result)
        # Progress callback with estimated time left
        if _progress_callback is not None:
            elapsed = time.time() - start_time if start_time else 0
            done = idx + 1
            est_left = (elapsed / done) * (total - done) if done > 0 else 0
            _progress_callback(done, total, est_left)

    ttest = pd.concat(ttest).set_index("metabolite")
    ttest = ttest.dropna()

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
        raise KeyError("No t-statistic column found in t-test results. Columns: " + str(df.columns))


    # Add a column for -log(p-corrected) and significance label
    df = df.copy()
    df["-log_p_corrected"] = df["p-corrected"].apply(lambda x: -np.log(x))
    df["sig_label"] = df["significance"].apply(lambda x: "significant" if x else "insignificant")
    df["metabolite_name"] = df.index


    if df.empty:
        st.warning("No t-test results to display. Please check your data or selection.")
        import plotly.graph_objects as go
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
    x_padding = abs(xlim[1] - xlim[0]) / 5 if xlim[1] != xlim[0] else 1
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
        xaxis_title="T",
        yaxis_title="-Log(p)",
        showlegend=True,  # Enable legend
        legend_title_text="Significance",
    )
    return fig


@st.cache_resource
def ttest_boxplot(df_ttest, metabolite):
    df = pd.concat([st.session_state.md, st.session_state.data], axis=1)
    df1 = pd.DataFrame(
        {
            metabolite: df[df[st.session_state.ttest_attribute] == st.session_state.ttest_options[0]].loc[:, metabolite],
            "option": st.session_state.ttest_options[0],
        }
    )
    df2 = pd.DataFrame(
        {
            metabolite: df[df[st.session_state.ttest_attribute] == st.session_state.ttest_options[1]].loc[:, metabolite],
            "option": st.session_state.ttest_options[1],
        }
    )
    df = pd.concat([df1, df2])
    fig = px.box(
        df,
        x="option",
        y=metabolite,
        color="option",
        width=350,
        height=400,
        points="all",
        custom_data=["option", metabolite],
    )

    # Set custom hovertemplate
    attribute = st.session_state.ttest_attribute if hasattr(st.session_state, "ttest_attribute") else "Attribute"
    hovertemplate = (
        f"attribute: {attribute}<br>"
        "option: %{customdata[0]}<br>"
        f"metabolite&name: {metabolite}<br>"
        "intensity: %{customdata[1]:.3g}<extra></extra>"
    )
    fig.update_traces(hovertemplate=hovertemplate)
    # Determine significance label
    is_significant = df_ttest.loc[metabolite, "significance"] if "significance" in df_ttest.columns else False
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
    pvalue = df_ttest.loc[metabolite, "p-corrected"]
    if pvalue >= 0.05:
        symbol = "ns"
    elif pvalue >= 0.01:
        symbol = "*"
    elif pvalue >= 0.001:
        symbol = "**"
    else:
        symbol = "***"

    top_y = max(df[metabolite]) * 1.2

    if isinstance(df["option"][0], str) and isinstance(df["option"][1], str):
        x0, x1 = 0, 1
    else:
        x0, x1 = st.session_state.ttest_options[0], st.session_state.ttest_options[1]
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
        y_margin = max(df[metabolite]) * 0.05
    else:
        y_margin = max(df[metabolite]) * 0.1
    fig.add_annotation(
        x=(x1-x0)/2,
        y=top_y + y_margin,
        text=f"<b>{symbol}</b>",
        showarrow=False,
        font_color="#555555",
    )
    return fig
