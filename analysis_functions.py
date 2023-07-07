import scanpy as sc
import anndata
import matplotlib.pyplot as plt

def load_libraries():
    # Scanpy and other libraries are imported at the start, no need for a specific function
    pass

def load_nsclc_data(uploaded_file):
    # Save uploaded file to a temporary file
    with open("temp.h5", "wb") as f:
        f.write(uploaded_file.getbuffer())
    # Now you can pass the path of the temporary file to read_10x_h5()
    adata = sc.read_10x_h5("temp.h5")
    return adata

def initialize_seurat_object(adata):
    # In Scanpy, this step can be mostly achieved while loading the data
    # Additional preprocessing can be performed as needed
    return adata


def perform_qc(adata, plt, sc):
    adata.var_names_make_unique()  # add this line to make var_names unique
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # Calculate n_genes_by_counts
    adata.obs['n_genes_by_counts'] = adata.X.sum(axis=1).A1    
    # Create the violin plots
    fig, axs = plt.subplots(1, 3, figsize=(12, 4))
    fig = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], multi_panel=True)
    
    # plt.savefig("qc.png")
    # plt.close()
    return adata, fig


def filter_data(adata, min_genes, min_cells, max_genes, max_mt):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs.n_genes_by_counts < max_genes, :]
    adata = adata[adata.obs.pct_counts_mt < max_mt, :]
    return adata

def normalize_data(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    return adata

def identify_highly_variable_features(adata):
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    return adata

def scale_data(adata):
    sc.pp.scale(adata, max_value=10)
    return adata

def perform_linear_dimensionality_reduction(adata):
    sc.tl.pca(adata)
    return adata

def cluster_data(adata, resolution=0.5):
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata, resolution=resolution)
    return adata

def perform_non_linear_dimensionality_reduction(adata):
    sc.tl.umap(adata)
    return adata
