import streamlit as st
from src.common import *

try:
    from src.pcoa import *

    page_setup()

    st.markdown("# Multivariate Statistics")
    st.markdown("### PERMANOVA & Principal Coordinate Analysis (PCoA)")

    with st.expander("📖 About"):
        st.markdown(
            """
            **Principal Coordinate Analysis (PCoA)** is an **unsupervised** ordination technique that visualizes relationships among samples based on a **distance (dissimilarity) matrix**. 
            Unlike PCA, which relies on Euclidean distance, PCoA can use **different distance metrics** such as Bray–Curtis, Jaccard, or Euclidean, making it suitable for non-normal or compositional data. 
            The method projects samples into new coordinate axes (PCo1, PCo2, etc.) that explain the greatest variation in distances. 
            Typically, the **first 10 coordinates** capture most of the variance, and users can select any two among these to visualize group separation.

            **PERMANOVA (Permutational Multivariate Analysis of Variance)** tests whether the **centroids of groups** differ significantly in multivariate space. 
            It uses **permutation-based resampling** (usually 999 permutations) to compute:
            - **Pseudo-F (test statistic):** measures the ratio of between-group to within-group variation  
            - **R²:** indicates the proportion of total variance explained by the grouping variable (metadata attribute)  
            - **p-value:** shows whether observed group differences are statistically significant under permutation  

            In this app, PERMANOVA results help quantify whether the separation observed in PCoA plots is statistically meaningful, rather than just visual. 
            """
)

        st.image("assets/figures/pcoa.png")


    if st.session_state.data is not None and not st.session_state.data.empty:
        c1, c2 = st.columns(2)
        with c1: 
            st.selectbox(
                "attribute for multivariate analysis",
                [c for c in st.session_state.md.columns if len(set(st.session_state.md[c])) > 1 and len(set(st.session_state.md[c])) != st.session_state.md.shape[0]],
                key="pcoa_attribute",
            )
        with c2: 
            st.selectbox(
                "distance matrix",
                ["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine", "euclidean", "hamming", "jaccard", "matching", "minkowski", "seuclidean", "sqeuclidean"],
                key="pcoa_distance_matrix",
                index = 6
            )

        max_allowed = min(st.session_state.data.shape[0], st.session_state.data.shape[1])
        n_component = min(10, max_allowed)
        st.session_state["n_components"] = n_component


        # Check if pcoa_attribute is valid before accessing DataFrame
        if (st.session_state.pcoa_attribute is None or st.session_state.pcoa_attribute not in st.session_state.md.columns):
            st.warning("Please select a valid attribute for PCoA.\n\nA valid attribute is required to group your samples for multivariate analysis. Make sure to choose a column from your metadata that contains categorical groupings (e.g., treatment, condition, or sample type) with at least two unique values. If no options appear, check that your metadata table is loaded and contains appropriate columns.")
            st.stop()

        attribute_series = st.session_state.md[st.session_state.pcoa_attribute].dropna()
        unique_categories = attribute_series.nunique()

        max_pcoa = 10
        pcoa_labels = [f"PC{i+1}" for i in range(max_pcoa)]

        col1, col2 = st.columns(2)
        with col1:
            pcoa_x_axis = st.selectbox("Interested X-axis for plot", pcoa_labels)
        with col2:
            pcoa_y_axis = st.selectbox("Interested Y-axis for plot", pcoa_labels, index = 1)

        att_col = st.session_state.pcoa_attribute
        
        options = sorted(st.session_state.md[att_col].dropna().unique())
        allowed_categories = st.multiselect(
            f"Filter categories in '{att_col}'", 
            options = options,
            default = options, 
            key = "pcoa_category_filter"
        )
        selected_categories = st.session_state.pcoa_category_filter
        filtered_md = st.session_state.md[st.session_state.md[att_col].isin(selected_categories)]
        filtered_data = st.session_state.data.loc[filtered_md.index]


        if len(filtered_md[att_col].unique()) < 2:
            st.warning("⚠️ PERMANOVA cannot be calculated for this group because there is only one category in the selected attribute. You need at least two categories to perform statistical testing.")
        else:
            permanova, pcoa_result = permanova_pcoa(
                filtered_data,
                st.session_state.pcoa_distance_matrix,
                st.session_state.md[st.session_state.pcoa_attribute],
            )

            if pcoa_x_axis == pcoa_y_axis:
                st.warning("⚠️ X-axis and Y-axis cannot be the same. Please choose different axes to view results.")
            else: 
                if not permanova.empty:
                    t1, t2, t3, t4 = st.tabs(["📁 PERMANOVA statistics", "📈 Principal Coordinate Analysis", "📊 Explained variance", "📁 Data"])
                    with t1:
                        show_table(permanova, "PERMANOVA-statistics")
                    with t2:
                        fig = get_pcoa_scatter_plot(
                            pcoa_result,
                            filtered_md,
                            st.session_state.pcoa_attribute,
                            pcoa_x_axis,
                            pcoa_y_axis, 
                            st.session_state.pcoa_category_filter
                        )
                        show_fig(fig, "principal-coordinate-analysis")
                    with t3:
                        fig = get_pcoa_variance_plot(pcoa_result)
                        show_fig(fig, "pcoa-variance")
                    with t4:
                        show_table(pcoa_result.samples.iloc[:, :10], "principal-coordinates")

    else:
        st.warning("⚠️ Please complete data preparation step first!")

except ModuleNotFoundError:
    st.error("This page requires the `skbio` package, which is not available in the Windows app.")
