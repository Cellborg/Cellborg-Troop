import signal
from fastapi import FastAPI, Response
from pydantic import BaseModel, Field
import os
import boto3
import scanpy as sc
import anndata as ad
import pandas as pd
import hdf5plugin
import json


app = FastAPI()

class initializeProjectRequest(BaseModel):
    user:str
    project: str
    datasets: list[str]

class clusteringRequest(BaseModel):
    user: str
    project: str
    resolution: float

class annoRequest(BaseModel):
    user: str
    project: str
    resolution: float
    annotations: object

# variables in Global context
user_environment = {}
adata = None
s3_client = boto3.client("s3") #may need to add access_key_id and secret_access_key
resolution_global= 0.0

principle_components=50

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
        s3_client.upload_file(localfile, user_environment['qc_dataset_bucket'], s3_key, Callback=print)
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
    global principle_components
    # ## Dimensionality Reduction
    # Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
    print("------- dimentionality reduction begins ------")
    sc.tl.pca(adata)

    #find ideal number of principle components
    df = adata.uns["pca"]['variance_ratio'].cumsum(axis=0)
    if df[-1] < 0.95:
        print("use 50")
    else:
        for num in range(len(df)):
            if df[num] >=0.95:
                principle_components = num+1
                break
                
    # Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function {func}`~scanpy.tl.leiden` or {func}`~scanpy.tl.tsne`. In our experience, there does not seem to be signifigant downside to overestimating the numer of principal components.
    sc.pl.pca_variance_ratio(adata, n_pcs=principle_components, log=True, save=".png")
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
    global resolution_global

# ## Clustering
# 
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) {cite}`traag2019louvain`. Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    print("------- clustering begins --------")
    sc.tl.leiden(adata, resolution = resolution, flavor="igraph", n_iterations=-1)
    sc.pl.umap(adata, color=["leiden"], save="2.png")
    targetfile=f"{s3_path}/QCclustering.png"
    upload_to_s3(targetfile,"./figures/umap2.png")
    resolution_global = resolution
    return (
        adata.var_names, adata.obs[f"leiden_res_${resolution}"].cat.categories
    )


def cell_type_annotation(annotations):
    global adata
    global resolution_global

    adata.obs["cell_type_lvl1"] = adata.obs[f"leiden_res_${resolution_global}"].map(annotations)

def upload_plot_to_s3(s3_key, localfile):
    # Upload the JSON data to S3
    s3_client.upload_file(localfile, user_environment['qc_dataset_bucket'], s3_key, Callback=print)
    print(f"Uploaded plot png to S3: {s3_key}")

@app.post("/init_endpoint", status_code=200)
async def initialize_project(initReq: initializeProjectRequest):
    try:
        global adata
        global user_environment
        global principle_components

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
        
        #create project_values.json in proj dir
        numpcs = {
            "num_pcs":principle_components
        }
        with open("project_values.json", "w") as outputfile:
            json.dump(numpcs, outputfile)
        print("created project_values.json file")

        input_var_path = f"{initReq.user}/{initReq.project}/project_values.json"
        upload_plot_to_s3(input_var_path, 'project_values.json')
        print("uploaded project_values.json file to s3")

        os.remove("project_values.json")
        print("project_values.json file successfully deleted")

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
        s3_path = f"{clustReq.user}/{clustReq.project}/project_values.json"
        gene_names, clustering = clustering(s3_path, clustReq.resolution)

        #download old project_values file from s3
        s3_client.Bucket(user_environment["qc_dataset_bucket"]).download_file('project_values.json', s3_path)

        #add clustering resolution
        with open('project_values.json') as f:
            data = json.load(f)
            data['clust_resolution'] = clustReq.resolution
            json.dump(data, f)
        print("Added clustering resolution to project values")
        
        #upload updated file to s3
        upload_plot_to_s3(f's3_path+/project_values.json','project_values.json')
        print("uploaded new project values to s3")
        
        #remove temp file locally
        os.remove("project_values.json")


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

@app.post("/annotations", status_code = 200)
async def annotations(annoRequest: annoRequest):
    global adata
    cell_type_annotation(annoRequest.annotations)
    #used to verify that annotations did work
    sc.pl.umap(
    adata,
    color=["cell_type_lvl1"],
    legend_loc="on data",
    )

