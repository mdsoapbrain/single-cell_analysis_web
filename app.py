import streamlit as st
import analysis_functions as af

# Define title
st.title("Single Cell RNA-Seq Data Analysis")

# Sidebar with action buttons and input fields
st.sidebar.header("Control Panel")
if st.sidebar.button("Load libraries"):
    af.load_libraries()

file = st.sidebar.file_uploader("Upload data file")
if file is not None:
    adata = af.load_nsclc_data(file)

if st.sidebar.button("Initialize Seurat object"):
    adata = af.initialize_seurat_object(adata)
    st.write("Seurat object initialized.")

if st.sidebar.button("Perform QC"):
    adata = af.perform_qc(adata)
    st.write("Quality control performed.")

if st.sidebar.button("Filter data"):
    adata = af.filter_data(adata)
    st.write("Data filtered.")

if st.sidebar.button("Normalize data"):
    af.normalize_data(adata)
    st.write("Data normalized.")

if st.sidebar.button("Identify highly variable features"):
    af.identify_highly_variable_features(adata)
    st.write("Highly variable features identified.")

if st.sidebar.button("Scale data"):
    af.scale_data(adata)
    st.write("Data scaled.")

if st.sidebar.button("Perform linear dimensionality reduction"):
    af.perform_linear_dimensionality_reduction(adata)
    st.write("Linear dimensionality reduction performed.")

if st.sidebar.button("Cluster data"):
    resolution = st.sidebar.number_input("Resolution", min_value=0.1, max_value=1.0, value=0.5)
    af.cluster_data(adata, resolution)
    st.write("Data clustered.")

if st.sidebar.button("Perform non-linear dimensionality reduction"):
    af.perform_non_linear_dimensionality_reduction(adata)
    st.write("Non-linear dimensionality reduction performed.")

# Main panel can be populated by writing to the streamlit display with st.write(), st.dataframe(), etc.
