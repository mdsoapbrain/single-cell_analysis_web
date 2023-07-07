import scanpy as sc
import anndata

def load_libraries():
    # Scanpy and other libraries are imported at the start, no need for a specific function
    pass

def load_nsclc_data(filename):
    adata = sc.read_10x_h5(filename)
    return adata

def initialize_seurat_object(adata):
    # In Scanpy, this step can be mostly achieved while loading the data
    # Additional preprocessing can be performed as needed
    return adata

def perform_qc(adata):
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    return adata

def filter_data(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
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
