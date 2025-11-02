import streamlit as st
import pandas as pd

from src.common import *
from src.ttest import *

page_setup()

st.markdown("# Student's t-test")

with st.expander("📖 About"):
    st.markdown(
        "The Student's t-test is a statistical test used to determine whether the means of two groups of data are significantly different from each other. The t-test is a parametric test that assumes the data are normally distributed and the variances of the two groups are equal. It is commonly used in hypothesis testing to determine whether there is a significant difference between the means of two populations or to compare the means of two samples. The t-test is a powerful tool in statistical analysis and is often the go-to test for researchers and analysts when analyzing data.")

# Ensure st.session_state.df_ttest is initialized
if "df_ttest" not in st.session_state:
    st.session_state.df_ttest = pd.DataFrame()

if st.session_state.data is not None and not st.session_state.data.empty:
    c1, c2 = st.columns(2)
    def clear_ttest_data():
        st.session_state.df_ttest = pd.DataFrame()
        # Remove tab-related session state if present
        for key in list(st.session_state.keys()):
            if key.startswith("ttest_metabolite") or key.startswith("_page_tab_"):
                del st.session_state[key]

    c1.selectbox(
        "select attribute of interest",
        options=[c for c in st.session_state.md.columns if len(set(st.session_state.md[c])) > 1],
        key="ttest_attribute",
        on_change=clear_ttest_data
    )

    if st.session_state.ttest_attribute is not None:
        attribute_options = list(
            set(st.session_state.md[st.session_state.ttest_attribute].dropna())
        )
        attribute_options.sort()
    else:
        attribute_options = []
    c2.multiselect(
        "select **two** options from the attribute for comparison",
        options=attribute_options,
        default=attribute_options[:2],
        key="ttest_options",
        max_selections=2,
        help="Select two options.",
        on_change=clear_ttest_data
    )
    c1, c2, c3 = st.columns(3)
    v_space(2, c1)
    c1.checkbox(
        "paired",
        False,
        key="ttest_paired",
        help="Specify whether the two observations are related (i.e. repeated measures) or independent.",
        on_change=clear_ttest_data
    )

    correction_options = {
        "auto": "auto",
        "Welch's": "True",
        "Student's": "False"
    }
    c2.selectbox(
        "T-test type",
        options=list(correction_options.keys()),
        key="ttest_correction_label",
        help="Choose Welch's (recommended) to use a version of the test that does not assume equal variances between groups. Choose Student's if you are confident that the group variances are equal. Choose Auto to use Welch's if the group sizes are unequal, or Student's if they are equal.",
        on_change=clear_ttest_data
    )
    c3.selectbox(
        "alternative",
        options=["two-sided", "greater", "less"],
        key="ttest_alternative",
        help="Defines the alternative hypothesis, or tail of the test.",
        on_change=clear_ttest_data
    )

    if c1.button("Run t-test", type="primary", disabled=(len(st.session_state.ttest_options) != 2)):
        # Map label to value for correction
        correction_value = correction_options[st.session_state.ttest_correction_label]
        
        # Add progress bar
        progress_placeholder = st.empty()
        time_placeholder = st.empty()
        
        def progress_callback(done, total, est_left):
            progress = done / total
            progress_placeholder.progress(progress, text=f"Running t-test: {done}/{total}")
            time_placeholder.info(f"Estimated time left: {int(est_left)} seconds")

        st.session_state.df_ttest = gen_ttest_data(
            st.session_state.ttest_attribute,
            st.session_state.ttest_options,
            st.session_state.ttest_paired,
            st.session_state.ttest_alternative,
            correction_value,
            corrections_map[st.session_state.p_value_correction],
            _progress_callback=progress_callback
        )
        
        progress_placeholder.empty()
        time_placeholder.empty()
        st.rerun()

    # Only show tabs if t-test results exist (button pressed and results generated)
    ttest_stat_cols = {"T", "T-val", "t", "tval"}
    df = st.session_state.df_ttest
    if df is not None and not df.empty and any(col in df.columns for col in ttest_stat_cols):
        # Added Volcano plot tab

        tabs = st.tabs(["📈 Feature significance", "📊 Single metabolite plots", "📈 Volcano plot", "📁 Data"])

        with tabs[0]:
            fig = plot_ttest(df)
            show_fig(fig, "t-test")

        with tabs[1]:
            metabolite_options = list(df.index)
            def metabolite_label(m):
                return m.split("&")[0] if "&" in m else m
            st.selectbox("metabolite", metabolite_options, key="ttest_metabolite", format_func=metabolite_label)
            # Only plot if the selected metabolite is in the data columns
            if st.session_state.ttest_metabolite in st.session_state.data.columns:
                fig = ttest_boxplot(df, st.session_state.ttest_metabolite)
                if fig is not None:
                    show_fig(fig, f"ttest-boxplot-{st.session_state.ttest_metabolite}", container_width=True)
            else:
                st.warning(f"Selected metabolite not found in data columns. Please select a valid metabolite.")

        with tabs[2]:
            # Volcano Plot tab (now third)
            if "mean(A)" in df.columns and "mean(B)" in df.columns:
                fig_volcano = get_ttest_volcano_plot(df)
                show_fig(fig_volcano, "ttest-volcano")
            else:
                st.warning("Could not generate volcano plot. Mean values are missing from t-test results. Please re-run the t-test.")

        with tabs[3]:
            # Format p-val and p-corrected columns in scientific notation, but allow full value on cell click
            df_display = df.copy()
            def sci_notation_or_plain(x):
                try:
                    if pd.isnull(x):
                        return x
                    if float(x) == 0:
                        return 0
                    return f"{x:.2e}"
                except Exception:
                    return x
            style_dict = {}
            for col in ["p-val", "p-corrected"]:
                if col in df_display.columns:
                    style_dict[col] = sci_notation_or_plain
            if style_dict:
                styled = df_display.style.format(style_dict)
                st.dataframe(styled, use_container_width=True)
            else:
                st.dataframe(df_display, use_container_width=True)

else:
    st.warning("⚠️ Please complete data preparation step first!")