import streamlit as st
import analysis_functions as af
import matplotlib.pyplot as plt
import scanpy as sc
# Define title
st.title("Single Cell RNA-Seq Data Analysis")

# Sidebar with action buttons and input fields
st.sidebar.header("Control Panel")
if st.sidebar.button("Load libraries"):
    af.load_libraries()
    st.write("scanpy, anndata")

file = st.sidebar.file_uploader("Upload data file")
if file is not None:
    adata = af.load_nsclc_data(file)

if st.sidebar.button("Initialize Seurat object"):
    adata = af.initialize_seurat_object(adata)
    st.write("Seurat object initialized.")

if st.sidebar.button("Perform QC"):
    adata, fig = af.perform_qc(adata, plt, sc)
    st.write("Quality control performed.")
    st.dataframe(adata.obs)
    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.pyplot(fig)

min_genes = st.sidebar.number_input("Minimum genes per cell", min_value=1, max_value=5000, value=200)
min_cells = st.sidebar.number_input("Minimum cells per gene", min_value=1, max_value=5000, value=3)
max_genes = st.sidebar.number_input("Maximum genes per cell", min_value=1, max_value=5000, value=2500)
max_mt = st.sidebar.number_input("Maximum MT percent", min_value=1, max_value=100, value=5)

if st.sidebar.button("Filter data"):
    adata = af.filter_data(adata, min_genes, min_cells, max_genes, max_mt)
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
