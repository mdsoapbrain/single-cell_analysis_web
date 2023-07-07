# Single Cell RNA-Seq Data Analysis

This project utilizes Python and several key libraries for single-cell RNA-Seq data analysis, including `scanpy`, `numpy`, `scipy`, `pandas`, `matplotlib`, `seaborn`, `scikit-learn`, `natsort`, `networkx`, `h5py`, and two graph clustering libraries, `louvain` and `leidenalg`.

## Environment Setup and Installation

1. It is recommended to use [Anaconda](https://www.anaconda.com/products/distribution) to manage your Python and package environment, but you can also use the Python environment already installed on your machine.

2. Create a new conda environment (optional):

    ```
    conda create -n single_cell_analysis python=3.8
    conda activate single_cell_analysis
    ```

3. Install packages:

    ```
    pip install -r requirements.txt
    ```

## Single Cell RNA-Seq Data Analysis Dashboard

We have provided an interactive dashboard using Streamlit, which makes the single-cell RNA-Seq data analysis process more visual and intuitive.

### Starting the Dashboard

Run the dashboard using the following command in your terminal:

    ```
    streamlit run app.py
    ```


Once the dashboard is running, you will see the title "Single Cell RNA-Seq Data Analysis" at the top of the page.

### Using the Control Panel

On the left-hand side of the dashboard, there's a sidebar titled "Control Panel". This contains buttons and input fields for the analysis.

- **Load libraries**: This button loads necessary Python libraries for the analysis.
- **Upload data file**: Here you can upload your data file. It should be in `.h5` format.
- **Initialize Seurat object**: Initializes the Seurat object using the uploaded data.
- **Perform QC**: Performs quality control checks on the data. Results will be displayed on the main page.
- **Filter data**: Use the four input fields under this section to specify the parameters for filtering the data. These parameters include "Minimum genes per cell", "Minimum cells per gene", "Maximum genes per cell", and "Maximum MT percent".
- **Normalize data**: Normalizes the data.
- **Identify highly variable features**: Identifies highly variable features in the data.
- **Scale data**: Scales the data.
- **Perform linear dimensionality reduction**: Performs PCA on the data.
- **Cluster data**: Uses the Leiden algorithm to cluster the data. The "Resolution" input field specifies the resolution for clustering.
- **Perform non-linear dimensionality reduction**: Performs UMAP on the data for visualization.

### Output and Results

The results of each step of the analysis will be displayed in the main area of the dashboard. This includes the results of the QC step, which are displayed as a dataframe, and a violin plot showing the distribution of gene counts, total counts, and MT percent counts in the cells. After each step, the state of the `Anndata` object is updated, and this updated state can be used as the input for the next step.

**Note**: Ensure you have set the option `'deprecation.showPyplotGlobalUse'` to `False` in the Streamlit settings to avoid deprecation warnings.
