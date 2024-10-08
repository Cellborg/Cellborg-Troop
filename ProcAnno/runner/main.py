import signal
from fastapi import FastAPI, Response
from pydantic import BaseModel, Field
import os
import boto3
import scanpy as sc
import anndata as ad
import pandas as pd
import hdf5plugin


app = FastAPI()

class initializeProjectRequest(BaseModel):
    user:str
    project: str
    datasets: list[str]

class clusteringRequest(BaseModel):
    user: str
    project: str
    resolution: int

# variables in Global context
user_environment = {}
adata = None
s3_client = boto3.client("s3") #may need to add access_key_id and secret_access_key

# set bucket values depending on the environment
def set_user_env():
    global user_environment

    # Get the environment variable  
    environment = os.getenv("ENVIRONMENT", default="dev")
    print(f"Cellborg Processing and Annotations Python container running in environment: {environment}")
    # Set bucket names based on the environment
    if environment == "dev":
        DATASET_BUCKET = "cellborgdatasetuploadbucket"
        QC_DATASET_BUCKET = "cellborgqcdatasetbucket"
        PCA_BUCKET = "cellborgpcabucket"
    else:
        DATASET_BUCKET = f"cellborg-{environment}-datasetupload-bucket"
        QC_DATASET_BUCKET = f"cellborg-{environment}-qcdataset-bucket"
        PCA_BUCKET = f"cellborg-{environment}-pca-bucket"

    # Assign bucket names to the user_environment dictionary
    user_environment["dataset_bucket"] = DATASET_BUCKET
    user_environment["qc_dataset_bucket"] = QC_DATASET_BUCKET
    user_environment["pca_bucket"] = PCA_BUCKET

def upload_to_s3(type, s3_key, localfile):
    # Upload the JSON data to S3
    if type=="pca":
        s3_client.upload_file(localfile, user_environment['pca_bucket'], s3_key, Callback=print)
    print(f"Uploaded plot png to S3: {s3_key}")

def initializeAdata(s3_singlets_path: str, datasets: list[str]):
    global adata
    global user_environment

    samples={datasetId:s3_singlets_path+f"/{datasetId}/singlets.h5ad" for datasetId in datasets} 
    adatas={}
    for sample_id, filepath in samples.items():
        response = s3_client.get_object(Bucket= user_environment["qc_dataset_bucket"], Key=filepath)
        sample_adata = sc.read_text(response['Body'].read().decode('utf-8'))
        #sample_adata.var_names_make_unique()
        adata[sample_id] = sample_adata
    
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()

def normalize():
    global adata
    print("-----normalize begins----")
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    print ("----normalize completed-----")

def feature_selection():
    global adata
    # ## Feature selection
    # 
    # As a next step, we want to reduce the dimensionality of the dataset and only include the most informative genes. This step is commonly known as feature selection. The scanpy function `pp.highly_variable_genes` annotates highly variable genes by reproducing the implementations of Seurat {cite}`Satija2015`, Cell Ranger {cite}`Zheng2017`, and Seurat v3 {cite}`stuart2019comprehensive` depending on the chosen `flavor`. 
    print("-------- feature selection begins --------")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pl.highly_variable_genes(adata)

def dimentionality_reduction(s3_path):
    global adata
    # ## Dimensionality Reduction
    # Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
    print("------- dimentionality reduction begins ------")
    sc.tl.pca(adata)

    #find ideal number of principle components
    df = adata.uns["pca"]['variance_ratio'].cumsum(axis=0)
    n_pcs = 50
    if df[-1] < 0.95:
        print("use 50")
    else:
        for num in range(len(df)):
            if df[num] >=0.95:
                n_pcs = num+1
                break
                
    # Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function {func}`~scanpy.tl.leiden` or {func}`~scanpy.tl.tsne`. In our experience, there does not seem to be signifigant downside to overestimating the numer of principal components.
    sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, save=".png")
    png_file1 = "./figures/pca_variance_ratio.png"
        # Create S3 object key for quality control data
    s3_key = f"{s3_path}/QCpca_variance_ratio.png"
    upload_to_s3("pca",s3_key,png_file1)
    
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
    s3_key = f"{s3_path}/QCpca.png"
    upload_to_s3("pca",s3_key,png_file2)

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
    targetfile=f"{s3_path}/QCnearest_neighbor.png"
    upload_to_s3(targetfile,"./figures/umap1.png")
def clustering(s3_path, resolution):
    global adata

# ## Clustering
# 
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) {cite}`traag2019louvain`. Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    print("------- clustering begins --------")
    sc.tl.leiden(adata, resolution = resolution, flavor="igraph", n_iterations=-1)
    sc.pl.umap(adata, color=["leiden"], save="2.png")
    targetfile=f"{s3_path}/QCclustering.png"
    upload_to_s3(targetfile,"./figures/umap2.png")
    return (
        adata.var_names, adata.obs['leiden'].cat.categories
    )


def cell_type_annotation(s3_path):
    global adata
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
    targetfile=f"{s3_path}/QCcelltype-annotation-1.png"
    upload_to_s3(targetfile,"./figures/umap5.png")

    sc.pl.umap(
        adata,
        color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
        legend_loc="on data",
        save="6.png",
    )
    targetfile=f"{s3_path}/QCcelltype-annotation-2.png"
    upload_to_s3(targetfile,"./figures/umap6.png")

@app.post("/init_endpoint", status_code=200)
async def initialize_project(initReq: initializeProjectRequest):
    try:
        global adata
        global user_environment

        s3_path = f"{initReq.user}/{initReq.project}"
        set_user_env()
        initializeAdata(s3_path, initReq.datasets)
        print('Successfully concatonated all datasets')
        normalize()
        print('Successfully normalized concatonated dataset')
        feature_selection()
        print('Successfully selected features')
        dimentionality_reduction(s3_path)
        print('Successfully reduced dimensions')


        return{
            "success": True,
            "message": "ProcAnno successfully initialized"
        }

    except Exception as err:
        print(err)
        return { 
            "success": False,
            "message": err
        }
@app.post("/clustering", status_code=200)
async def do_clustering(clustReq: clusteringRequest):
    try:
        s3_path = f"{clustReq.user}/{clustReq.project}"
        gene_names, clustering = clustering(s3_path, clustReq.resolution)
        return {
            "success": True,
            "message": "Clustering successfully finished",
            "gene_names": gene_names,
            "clusters": clustering
        }
    except Exception as err:
        print(err)
        return {
            "success": False,
            "message": err
        }
        

