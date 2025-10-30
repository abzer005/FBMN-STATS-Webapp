import streamlit as st

from src.common import *
from src.randomforest import *
def clear_rf_outputs():
    for key in [
        'df_oob', 'df_important_features', 'log', 'class_report', 'label_mapping',
        'test_confusion_df', 'train_confusion_df', 'test_accuracy', 'train_accuracy']:
        if key in st.session_state:
            del st.session_state[key]

page_setup()

st.title("Random Forest")

with st.expander("📖 About"):
    st.markdown(
        """Get the most important features explaining the selected attribute with supervised learning via random forest model."""
    )
    st.image("assets/figures/random-forest.png")



if st.session_state.data is not None and not st.session_state.data.empty:
    # Preserve original data and metadata
    if 'data_full' not in st.session_state:
        st.session_state['data_full'] = st.session_state.data.copy()
    if 'md_full' not in st.session_state:
        st.session_state['md_full'] = st.session_state.md.copy()


    use_random_seed = st.checkbox(
        'Use a fixed random seed for reproducibility',
        True,
        key='use_random_seed',
        on_change=clear_rf_outputs
    )
    c1, c2 = st.columns(2)
    c1.selectbox(
        "attribute for supervised learning feature classification",
        options=[c for c in st.session_state.md_full.columns if len(set(st.session_state.md_full[c])) > 1],
        key="rf_attribute",
        on_change=clear_rf_outputs
    )

    # Always use the full metadata for category options
    if st.session_state.rf_attribute is not None and st.session_state.rf_attribute in st.session_state.md_full.columns:
        rf_categories_options = sorted(set(st.session_state.md_full[st.session_state.rf_attribute].dropna()))
    else:
        rf_categories_options = []


    c2.multiselect(
        "select at least 2 categories to include (optional)",
        options=rf_categories_options,
        default=rf_categories_options,
        key="rf_categories",
        help="If you want to include only specific categories for classification, select them here. Otherwise, all categories will be used.",
        on_change=clear_rf_outputs
    )

    # Disable the button if less than two categories are selected
    selected_categories = st.session_state.get("rf_categories", [])
    button_disabled = len(selected_categories) < 2


    c1.number_input(
        "number of trees", 1, 500, 100, 50,
        key = "rf_n_trees",
        help="number of trees for random forest, check the OOB error plot and select a number of trees where the error rate is low and flat",
        on_change=clear_rf_outputs
    )
    
    random_seed = 123 if use_random_seed else None

    if c2.button("Run supervised learning", type="primary", disabled=button_disabled):
        try:
            # Filter data and metadata to only include selected categories, but do NOT overwrite originals
            selected_categories = st.session_state.rf_categories
            if selected_categories:
                mask = st.session_state.md_full[st.session_state.rf_attribute].isin(selected_categories)
                data_filtered = st.session_state.data_full[mask]
                md_filtered = st.session_state.md_full[mask]
            else:
                data_filtered = st.session_state.data_full.copy()
                md_filtered = st.session_state.md_full.copy()

            # Temporarily set filtered data for model
            st.session_state.data = data_filtered
            st.session_state.md = md_filtered

            df_oob, df_important_features, log, class_report, label_mapping, test_confusion_df, train_confusion_df, test_accuracy, train_accuracy = run_random_forest(st.session_state.rf_attribute, st.session_state.rf_n_trees, random_seed)
            st.session_state['df_oob'] = df_oob
            st.session_state['df_important_features'] = df_important_features
            st.session_state['log'] = log
            st.session_state['class_report'] = class_report
            st.session_state['label_mapping'] = label_mapping
            st.session_state['test_confusion_df'] = test_confusion_df
            st.session_state['train_confusion_df'] = train_confusion_df
            st.session_state['test_accuracy'] = test_accuracy
            st.session_state['train_accuracy'] = train_accuracy

            # Restore full data/metadata after model run
            st.session_state.data = st.session_state.data_full.copy()
            st.session_state.md = st.session_state.md_full.copy()
        except Exception as e:
            st.error(f"Failed to run model due to: {str(e)}")
else:
    st.warning("⚠️ Please complete data preparation step first!")

if ('df_important_features' in st.session_state and st.session_state.df_important_features is not None and not st.session_state.df_important_features.empty):
    tabs = st.tabs(["📈 Analyze optimum number of trees", 
                    "📁 Feature ranked by importance", 
                    "📋 Classification Report",
                    "🔍 Confusion Matrix"])
    with tabs[0]:
        fig = get_oob_fig(st.session_state.df_oob)
        show_fig(fig, "oob-error")
    with tabs[1]:
        show_table(st.session_state.df_important_features)
    with tabs[2]:  # Classification Report
        if 'log' in st.session_state:
            st.subheader("Log Messages")
            st.text(st.session_state.log)

        if 'class_report' in st.session_state and 'label_mapping' in st.session_state:
            st.subheader("Classification Report")
        
            # Convert the classification report string to DataFrame
            class_report_df = classification_report_to_df(st.session_state.class_report)
        
            # Convert the label mapping string to DataFrame
            label_mapping_df = label_mapping_to_df(st.session_state.label_mapping)
           
            # Ensure class_report_df's index is set correctly for merging
            class_report_df['class'] = class_report_df['class'].astype(str)
        
            # Merge the DataFrames on 'Class Index'
            merged_df = pd.merge(class_report_df, label_mapping_df, on='class')
            merged_df.set_index('Label', inplace=True)
            st.dataframe(merged_df)
    with tabs[3]:
        st.subheader("Test Set Confusion Matrix")
        st.dataframe(st.session_state.test_confusion_df)
        st.write(f"Test Set Accuracy: {st.session_state.test_accuracy:.2%}")

        st.subheader("Train Set Confusion Matrix")
        st.dataframe(st.session_state.train_confusion_df)
        st.write(f"Train Set Accuracy: {st.session_state.train_accuracy:.2%}")
