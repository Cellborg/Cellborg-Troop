from fastapi import FastAPI
from pydantic import BaseModel, Field
import os
import boto3
import scanpy as sc
import anndata as ad
import pandas as pd

app = FastAPI()

class QCRequest(BaseModel):
    user: str
    project: str
    dataset: str
    min: int
    max: int
    mt: int

class QCResponse(BaseModel):
    sucess: bool
    message: str
    cell_count: int
    gene_count: int

# variables in Global context
user_environment = {}
adata = ad.AnnData()
workspace_path = r'/tmp'

# set bucket values depending on the environment
def set_user_env():
    global user_environment
    global workspace_path
    # Create a dictionary-like object (similar to R's new.env())
    # Get the environment variable
    environment = os.getenv("ENVIRONMENT", default="dev")
    # Print a message
    print(f"Cellborg QCProcessing R container running in environment: {environment}")
    # Set bucket names based on the environment
    if environment == "dev":
        DATASET_BUCKET = "cellborgdatasetuploadbucket"
        QC_DATASET_BUCKET = "cellborgqcdatasetbucket"
    else:
        DATASET_BUCKET = f"cellborg-{environment}-datasetupload-bucket"
        QC_DATASET_BUCKET = f"cellborg-{environment}-qcdataset-bucket"

    # Assign bucket names to the user_environment dictionary
    user_environment["dataset_bucket"] = DATASET_BUCKET
    user_environment["qc_dataset_bucket"] = QC_DATASET_BUCKET

    #create temp folder for keeping files in the workspace
    #handle windows
    # if os.name == 'nt':
    #     workspace_path = r'C:\temp' 
    # if not os.path.exists(workspace_path):
    #     os.makedirs(workspace_path)

# load dataset files for QC from S3 bucket to the local file system under /tmp
# TODO: change /tmp to something else as this might cause problems.

def load_dataset(qcreq):
    # Load dataset and create Seurat object
    print("Loading dataset and creating Seurat object")
    prefix = f"{qcreq.user}/{qcreq.project}/{qcreq.dataset}/"
    print(prefix)
    
    s3 = boto3.client("s3")

    # List files in the dataset bucket
    response = s3.list_objects_v2(Bucket=user_environment["dataset_bucket"], Prefix=prefix)
    files = response.get("Contents", [])

    for file in files:
        file_name = os.path.basename(file["Key"])
        print(file_name)
        # Save the object locally (you can customize this part)
        local_path = os.path.join(workspace_path, file_name)
        s3.download_file(user_environment["dataset_bucket"], file["Key"], local_path)
    # return value , exception handling - TODO

# load dataset files create AnnData object
def read_10x_mtx():
    adata = sc.read_10x_mtx(
    workspace_path,  # the directory with the `.mtx` file
    var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
    cache=True,  # write a cache file for faster subsequent reading
    )
    return adata

def calculate_qc_metrics():
    global adata
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    # The scanpy function {func}`~scanpy.pp.calculate_qc_metrics` calculates common quality control (QC) metrics, which are largely based on `calculateQCMetrics` from scater {cite}`McCarthy2017`. One can pass specific gene population to {func}`~scanpy.pp.calculate_qc_metrics` in order to calculate proportions of counts for these populations. Mitochondrial, ribosomal and hemoglobin genes are defined by distinct prefixes as listed below. 

    # mitochondrial genes, "MT-" for human, "mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

def voilin_plot():
    # One can now inspect violin plots of some of the computed QC metrics:
    # 
    # * the number of genes expressed in the count matrix
    # * the total counts per cell
    # * the percentage of counts in mitochondrial genes
    # saves violin image file as violin.png, copy this to S3 under plots folder.
    print("------ voilin_plot begins -------")
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save=".png",
    )

    png_file = "./figures/violin.png"

    # Create S3 object key for quality control data
    s3_key = f"{s3_plots_dir}/QcViolinPlot.png"
    upload_plot_to_s3(s3_key,png_file)

def scatter_plot():
    print("------scatter_plot begins -------")
    adata1 = adata[adata.obs.pct_counts_mt < 8, :].copy()
    # Additionally, it is useful to consider QC metrics jointly by inspecting a scatter plot colored by `pct_counts_mt`. 
    sc.pl.scatter(
        adata1, 
        "total_counts", 
        "n_genes_by_counts", 
        color="pct_counts_mt",
        save=".png",
    )
    png_file = "./figures/scatter.png"
        # Create S3 object key for quality control data
    s3_key = f"{s3_plots_dir}/QcScatter.png"
    upload_plot_to_s3(s3_key,png_file)
 
def doublet_detection():
    print("------- doublet detection begins ------")
    sc.pp.scrublet(adata)

def normalize():
    print("-----normalize begins----")
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    print ("----normalize completed-----")

def feature_selection():
    # ## Feature selection
    # 
    # As a next step, we want to reduce the dimensionality of the dataset and only include the most informative genes. This step is commonly known as feature selection. The scanpy function `pp.highly_variable_genes` annotates highly variable genes by reproducing the implementations of Seurat {cite}`Satija2015`, Cell Ranger {cite}`Zheng2017`, and Seurat v3 {cite}`stuart2019comprehensive` depending on the chosen `flavor`. 
    print("-------- feature selection begins --------")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pl.highly_variable_genes(adata)

def dimentionality_reduction():
    # ## Dimensionality Reduction
    # Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
    print("------- dimentionality reduction begins ------")
    sc.tl.pca(adata)

    # Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function {func}`~scanpy.tl.leiden` or {func}`~scanpy.tl.tsne`. In our experience, there does not seem to be signifigant downside to overestimating the numer of principal components.
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save=".png")
    png_file1 = "./figures/pca_variance_ratio.png"
        # Create S3 object key for quality control data
    s3_key = f"{s3_plots_dir}/QCpca_variance_ratio.png"
    upload_plot_to_s3(s3_key,png_file1)
    
    # You can also plot the principal components to see if there are any potentially undesired features (e.g. batch, QC metrics) driving signifigant variation in this dataset. In this case, there isn't anything too alarming, but it's a good idea to explore this.
    sc.pl.pca(
        adata,
        color=["pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (1, 2)],
        ncols=2,
        size=2,
        save=".png",
    )
    # upload plots to s3
    png_file2 = "./figures/pca.png"
        # Create S3 object key for quality control data
    s3_key = f"{s3_plots_dir}/QCpca.png"
    upload_plot_to_s3(s3_key,png_file2)


def nearest_neighbor_graph():
# ## Nearest neighbor graph constuction and visualization
# 
# Let us compute the neighborhood graph of cells using the PCA representation of the data matrix.
    print("----- nearest_neighbor_graph begins -----")
    sc.pp.neighbors(adata)
    # This graph can then be embedded in two dimensions for visualiztion with UMAP (McInnes et al., 2018):
    sc.tl.umap(adata)
    # We can now visualize the UMAP according to the `sample`. 
    sc.pl.umap(
        adata,
        # Setting a smaller point size to get prevent overlap
        size=2,
        save="1.png"
    )
    targetfile=f"{s3_plots_dir}/QCnearest_neighbor.png"
    upload_plot_to_s3(targetfile,"./figures/umap1.png")

def clustering():
# ## Clustering
# 
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) {cite}`traag2019louvain`. Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    print("------- clustering begins --------")
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
    sc.pl.umap(adata, color=["leiden"], save="2.png")
    targetfile=f"{s3_plots_dir}/QCclustering.png"
    upload_plot_to_s3(targetfile,"./figures/umap2.png")

def reassess_qc_and_filtering():
# ## Re-assess quality control and cell filtering 
# 
# As indicated before, we will now re-assess our filtering strategy by visualizing different QC metrics using UMAP. 
    sc.pl.umap(
        adata,
        color=["leiden", "predicted_doublet", "doublet_score"],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
        save="3.png"
    )
    targetfile=f"{s3_plots_dir}/QCre-assess.png"
    upload_plot_to_s3(targetfile,"./figures/umap3.png")

    sc.pl.umap(
        adata,
        color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
        save="4.png"
    )
    targetfile=f"{s3_plots_dir}/QCre-assess-cell-filtering.png"
    upload_plot_to_s3(targetfile,"./figures/umap4.png")

def cell_type_annotation():
    # ## Manual cell-type annotation

# :::{note}
# This section of the tutorial is expanded upon using prior knowledge resources like automated assignment and gene enrichment in the scverse tutorial [here](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html#cell-type-annotation)
# :::

# Cell type annotation is laborous and repetitive task, one which typically requires multiple rounds of subclustering and re-annotation. It's difficult to show the entirety of the process in this tutorial, but we aim to show how the tools scanpy provides assist in this process.

# We have now reached a point where we have obtained a set of cells with decent quality, and we can proceed to their annotation to known cell types. Typically, this is done using genes that are exclusively expressed by a given cell type, or in other words these genes are the marker genes of the cell types, and are thus used to distinguish the heterogeneous groups of cells in our data. Previous efforts have collected and curated various marker genes into available resources, such as [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/), [TF-Marker](http://bio.liclab.net/TF-Marker/), and [PanglaoDB](https://panglaodb.se/). The [cellxgene gene expression tool](https://cellxgene.cziscience.com/gene-expression) can also be quite useful to see which cell types a gene has been expressed in across many existing datasets.

# Commonly and classically, cell type annotation uses those marker genes subsequent to the grouping of the cells into clusters. So, let's generate a set of clustering solutions which we can then use to annotate our cell types. Here, we will use the Leiden clustering algorithm which will extract cell communities from our nearest neighbours graph.
    for res in [0.02, 0.5, 2.0]:
        sc.tl.leiden(
            adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
        )
    # Notably, the number of clusters that we define is largely arbitrary, and so is the `resolution` parameter that we use to control for it. As such, the number of clusters is ultimately bound to the stable and biologically-meaningful groups that we can ultimately distringuish, typically done by experts in the corresponding field or by using expert-curated prior knowledge in the form of markers.
    sc.pl.umap(
        adata,
        color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
        legend_loc="on data",
        save="5.png",
    )
    targetfile=f"{s3_plots_dir}/QCcelltype-annotation-1.png"
    upload_plot_to_s3(targetfile,"./figures/umap5.png")

    sc.pl.umap(
        adata,
        color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
        legend_loc="on data",
        save="6.png",
    )
    targetfile=f"{s3_plots_dir}/QCcelltype-annotation-2.png"
    upload_plot_to_s3(targetfile,"./figures/umap6.png")



def upload_plot_to_s3(s3_key, localfile):
    # Upload the JSON data to S3
    s3.upload_file(localfile, user_environment['qc_dataset_bucket'], s3_key, Callback=print)
    print(f"Uploaded plot png to S3: {s3_key}")

#----- main -------
app = FastAPI()

@app.post("/qc_endpoint", status_code=201)
async def do_qc(qcreq: QCRequest):
    s3_plots_dir = f"{qcreq.user}/{qcreq.project}/{qcreq.dataset}/plots/"
    set_user_env()
    load_dataset(qcreq)
    adata = read_10x_mtx()
    calculate_qc_metrics()
    voilin_plot()
    scatter_plot()
    doublet_detection()
    normalize()
    feature_selection()
    dimentionality_reduction()
    nearest_neighbor_graph()
    clustering()
    reassess_qc_and_filtering()
    cell_type_annotation()
    return fastapi.Response(status_code=200, content='QC completed...')

@app.post("/shutdown")
async def shutdown():
    os.kill(os.getpid(), signal.SIGTERM)
    return fastapi.Response(status_code=200, content='Server shutting down...')