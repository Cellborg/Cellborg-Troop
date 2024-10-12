import signal
from fastapi import FastAPI, Response
from pydantic import BaseModel, Field
import os
import boto3
import scanpy as sc
import anndata as ad
import pandas as pd
import hdf5plugin
from datetime import datetime, timedelta
import json

app = FastAPI()

class QCPrePlotRequest(BaseModel):
    user: str
    project: str
    dataset: str
    mt: str


class QCRequest(BaseModel):
    user: str
    project: str
    dataset: str
    min: int
    max: int
    mt: str

class QCDoublets(BaseModel):
   user:str
   project: str
   dataset: str
   countMax: int
   countMin: int
   geneMax:  int
   geneMin:int
   mitoMax: int
   mitoMin: int


class QCFinish(BaseModel):
    user: str
    project: str
    dataset: str
    doubletScore: float

class QCResponse(BaseModel):
    success: bool
    message: str
    cell_count: int
    gene_count: int

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
adata = ad.AnnData()
workspace_path = r'/tmp'
s3_plots_dir = ""
s3 = boto3.client("s3")
resolution_global= 0.0

principle_components=50

# set bucket values depending on the environment
def set_user_env():
    global user_environment
    global workspace_path
    # Create a dictionary-like object (similar to R's new.env())
    # Get the environment variable
    environment = os.getenv("ENVIRONMENT", default="dev")
    # Print a message
    print(f"Cellborg QCProcessing Python container running in environment: {environment}")
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

def calculate_qc_metrics(mt: str):
    global adata
    min_genes=200
    min_cells=3
    sc.pp.filter_cells(adata, min_genes)
    sc.pp.filter_genes(adata, min_cells)
    # The scanpy function {func}`~scanpy.pp.calculate_qc_metrics` calculates common quality control (QC) metrics, which are largely based on `calculateQCMetrics` from scater {cite}`McCarthy2017`. One can pass specific gene population to {func}`~scanpy.pp.calculate_qc_metrics` in order to calculate proportions of counts for these populations. Mitochondrial, ribosomal and hemoglobin genes are defined by distinct prefixes as listed below. 

    # mitochondrial genes, "MT-" for human, "mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith(mt)
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

def voilin_plot():
    global adata
    global s3_plots_dir
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
    global adata

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
    global adata
    print("------- doublet detection begins ------")
    sc.pp.scrublet(adata)
    print("------doublet detection is finished------")

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
    s3_key = f"{s3_path}/QCpca.png"
    upload_plot_to_s3(s3_key,png_file2)

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
    upload_plot_to_s3(targetfile,"./figures/umap1.png")


def nearest_neighbor_graph():
    global adata
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
    upload_plot_to_s3(targetfile,"./figures/umap2.png")
    resolution_global = resolution
    return (
        adata.var_names, adata.obs["leiden"].cat.categories
    )


def gating_adata(countMx, countMn, geneMx, geneMn, mitoMx, mitoMn):
    global adata
    adata = adata[adata.obs.n_genes_by_counts < geneMx, :]
    adata = adata[adata.obs.n_genes_by_counts > geneMn, :]
    adata = adata[adata.obs.pct_counts_mt < mitoMx, :]
    adata = adata[adata.obs.pct_counts_mt > mitoMn, :]
    adata = adata[adata.obs.total_counts < countMx, :]
    adata = adata[adata.obs.total_counts > countMn, :]
    return adata

def reassess_qc_and_filtering():
    global adata
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

def cell_type_annotation(annotations):
    global adata
    global resolution_global

    adata.obs["cell_type_lvl1"] = adata.obs[f"leiden_res_{resolution_global}"].map(annotations)



def upload_plot_to_s3(s3_key, localfile):
    # Upload the JSON data to S3
    s3.upload_file(localfile, user_environment['qc_dataset_bucket'], s3_key, Callback=print)
    print(f"Uploaded plot png to S3: {s3_key}")

def print_time(msg):
    time_now = datetime.now()
    date_time = time_now.strftime("%m/%d/%Y, %H:%M:%S")
    print("===%s: %s" % (msg, time_now))

def initializeAdata(s3_singlets_path: str, datasets: list[str]):
    global adata
    global user_environment

    samples={datasetId:s3_singlets_path+f"/{datasetId}/singlets.h5ad" for datasetId in datasets} 
    adatas={}
    for sample_id, filepath in samples.items():
        s3.download_file(
            Bucket = user_environment["qc_dataset_bucket"], 
            Key = filepath, 
            Filename= 'singlets.h5ad') 
        sample_adata = sc.read_h5ad("singlets.h5ad")
        adatas[sample_id] = sample_adata
        os.remove('singlets.h5ad')
        print("removed file from local")
        #response = s3.get_object(Bucket= user_environment["qc_dataset_bucket"], Key=filepath)
        #print(f"Pulled {sample_id} from s3: ")
        #file_pulled = response['Body'].read()
        #print(file_pulled)
        #print(type(file_pulled))
        #sample_adata = sc.read_h5ad(file_pulled)
        #sample_adata.var_names_make_unique()
        #adatas[sample_id] = sample_adata
    
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()

#----- main -------
app = FastAPI()

@app.post("/qc_pre_plot_endpoint", status_code = 200)
async def do_pre_plot_qc(qcreq: QCPrePlotRequest):
    try:
        t1=datetime.now()
        global adata
        global s3_plots_dir

        s3_plots_dir = f"{qcreq.user}/{qcreq.project}/{qcreq.dataset}/plots"
        print_time("[load_dataset_from_s3]")
        set_user_env()
        load_dataset(qcreq)
        adata = read_10x_mtx()
        print_time("[calculate_qc_metrics]")
        calculate_qc_metrics(qcreq.mt)
        print_time("[generate_plots]")
        voilin_plot()
        
        t2=datetime.now()
        time_elapsed=t2-t1
        #date_time = time_elapsed.strftime("%m/%d/%Y, %H:%M:%S")
        print("[time_elapsed]: %s" % time_elapsed)
        return {"success": True,
                "message": "QC Pre-Plot Completed Successfully"
                }
    except Exception as err:
        print('ERROR: ',err)
        return {"success": False,
                "message": err}

@app.post("/qc_doublet_endpoint", status_code=200)
async def do_doublet_plot_qc(qcreq: QCDoublets):
    try:
        t1=datetime.now()
        #gate adata
        global adata
        countMx = qcreq.countMax
        countMn = qcreq.countMin
        geneMx = qcreq.geneMax
        geneMn = qcreq.geneMin
        mitoMx = qcreq.mitoMax
        mitoMn = qcreq.mitoMin

        print_time("[doublet_detection]")

        gating_adata(countMx, countMn, geneMx, geneMn, mitoMx, mitoMn)
        doublet_detection()

        #doublet detection
        singlets = adata[adata.obs.predicted_doublet == False]
        #create temperary file
        singlets.write_h5ad(
            "tmp.h5ad",
        )
        #upload file to s3
        print_time("[upload_plot_to_s3]")
        s3_key = f"{qcreq.user}/{qcreq.project}/{qcreq.dataset}/singlets.h5ad"
        upload_plot_to_s3(s3_key, 'tmp.h5ad')
        print("Successfully uploaded singlets to s3!!!")
        #del temp file
        os.remove("tmp.h5ad")
        print("temp file successfully deleted")

        counts = {
            "countMax": qcreq.countMax,
            "countMin": qcreq.countMin,
            "geneMax":qcreq.geneMax,
            "geneMin": qcreq.geneMin,
            "mitoMax": qcreq.mitoMax,
            "mitoMin":qcreq.mitoMin
        }
        with open("input_vars.json", "w") as outputfile:
            json.dump(counts, outputfile)
        print("created input_vars.json file")

        input_var_path = f"{qcreq.user}/{qcreq.project}/{qcreq.dataset}/input_vars.json"
        upload_plot_to_s3(input_var_path, 'input_vars.json')
        print("uploaded input_vars.json file to s3")

        os.remove("input_vars.json")
        print("input_vars.json file successfully deleted")

        t2=datetime.now()
        time_elapsed=t2-t1
        print("[time_elapsed]: %s" % time_elapsed)
        return{"success": True,
            "message": "QC Completed Successfully",
            }
    except Exception as err:
        print('ERROR: ',err)
        return{"success": False,
               "message": err}
    
#------------------- Processing & Annotations-----------------

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
        print('ERROR: ',err)
        return { 
            "success": False,
            "message": err
        }
    
@app.post("/clustering", status_code=200)
async def do_clustering(clustReq: clusteringRequest):
    try:
        s3_path = f"{clustReq.user}/{clustReq.project}"
        gene_names, clusters = clustering(s3_path, clustReq.resolution)

        #download old project_values file from s3
        s3.download_file(
            Bucket = user_environment["qc_dataset_bucket"], 
            Key = s3_path+'/project_values.json', 
            Filename= 'project_values.json') 

        #add clustering resolution
        with open('project_values.json', "r+") as f:
            data = json.load(f)
            data['clust_resolution'] = clustReq.resolution
            json.dump(data, f)
        print("Added clustering resolution to project values")
        
        #upload updated file to s3
        upload_plot_to_s3(f"{s3_path}/project_values.json",'project_values.json')
        print("uploaded new project values to s3")
        
        #remove temp file locally
        os.remove("project_values.json")

        print(f"gene_names: {type(gene_names)}", gene_names)
        print(f"cluster: {type(clusters)}", clusters)
        return {
            "success": True,
            "message": "Clustering successfully finished",
            "gene_names": gene_names,
            "clusters": clusters
        }
    except Exception as err:
        print('ERROR: ',err)
        return {
            "success": False,
            "message": err
        }
    
@app.post("/annotations", status_code = 200)
async def annotations(annoRequest: annoRequest):
    global adata
    try:
        cell_type_annotation(annoRequest.annotations)
        #used to verify that annotations did work
        sc.pl.umap(
        adata,
        color=["cell_type_lvl1"],
        legend_loc="on data",
        )
        return{
            "success":True,
            "message":"Annotatons successfully completed"
        }
    except Exception as err:
        print('ERROR: ',err)
        return{
            "success":False,
            "message": err
        }

#@app.post("/qc_finish_doublet_endpoint", status_code = 200)
#async def finish_doublet(qcreq: QCFinish):
#    doubletScore = qcreq.doubletScore
#    singlets = adata[adata.obs.doublet_score < doubletScore]
#    #save adata to s3
#    try:
#        singlets.write_h5ad(
#            "tmp.h5ad",
#        )
#        s3_key = f"{qcreq.user}/{qcreq.project}/{qcreq.dataset}/singlets.h5ad"
#        upload_plot_to_s3(s3_key, 'tmp.h5ad')
#        print("Successfully uploaded singlets to s3!!!")
#        os.remove("temp.h5ad")
#        print("temp file successfully deleted")
#        return{"success": True,
#                "message": "QC Completed Successfully"}
#    except Exception as err:
#        print(err)
#        return{"success": False,
#               "message": err}


#@app.post("/qc_endpoint", status_code=200)
#async def do_qc(qcreq: QCRequest):
#    global adata
#    global s3_plots_dir
#
#    s3_plots_dir = f"{qcreq.user}/{qcreq.project}/{qcreq.dataset}/plots"
#    set_user_env()
#    load_dataset(qcreq)
#    adata = read_10x_mtx()
#    calculate_qc_metrics(qcreq.mt)
#    voilin_plot()
#    scatter_plot()
#
#   #doublets
#    doublet_detection()
#    #everything below this needs to be in analysis
#    normalize()
#    feature_selection()
#    dimentionality_reduction()
#    nearest_neighbor_graph()
#    clustering()
#    reassess_qc_and_filtering()
#    cell_type_annotation()
#    return {"success":       "TRUE", 
#            "message":      "QC Completed Successfully",
#            "cell_count":   adata.n_obs,
#            "gene_count":   adata.n_vars
#    }

@app.post("/shutdown")
async def shutdown():
    os.kill(os.getpid(), signal.SIGTERM)
    return Response(status_code=200, content='Server shutting down...')