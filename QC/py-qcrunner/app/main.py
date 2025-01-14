import signal
from fastapi import FastAPI, Response
from pydantic import BaseModel
import os
import boto3
import scanpy as sc
import anndata as ad
import pandas as pd
from datetime import datetime
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
    global workspace_path
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
    global workspace_path
    try:
        adata = sc.read_10x_mtx(
        workspace_path,  # the directory with the `.mtx` file
        var_names="gene_symbols",  # use gene symbols for the variable names (variables-axis index)
        cache=False,  # write a cache file for faster subsequent reading
        )
        return adata
    except:
        adata = sc.read_mtx(f'{workspace_path}/matrix.mtx')
        adata_bc=pd.read_csv(f'{workspace_path}/barcodes.tsv',header=None)
        adata_features=pd.read_csv(f'{workspace_path}/features.tsv',header=None, delimiter = '\t')
        adata= adata.T
        adata.obs['cell_id']= adata_bc
        adata.var['gene_name'] = adata_features[1].values
        adata.var.index= adata.var['gene_name']

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

    print("Creating violin df")
    print(adata.obs.columns)
    data_df = adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]]

    print("Creating violin json")
    data_for_highcharts = {
        index:{
                "n_genes":row['n_genes_by_counts'],
                "total_counts": row["total_counts"],
                "pct_counts_mt":row["pct_counts_mt"]
        }
        for index, row in data_df.iterrows()
    }
    print("uploading json file")
    with open("highcharts_data.json", "w") as f:
        json.dump(data_for_highcharts, f, indent=4)
    upload_plot_to_s3(f"{s3_plots_dir}/QCViolinPlot.json", "highcharts_data.json")

    print("removing temp json file")
    os.remove("highcharts_data.json")

    #sc.pl.violin(
    #    adata,
    #    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    #    jitter=0.4,
    #    multi_panel=True,
    #    save=".png",
    #)

    #png_file = "./figures/violin.png"

    # Create S3 object key for quality control data
    #s3_key = f"{s3_plots_dir}/QcViolinPlot.png"
    #upload_plot_to_s3(s3_key,png_file)

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




def upload_plot_to_s3(s3_key, localfile):
    # Upload the JSON data to S3
    s3.upload_file(localfile, user_environment['qc_dataset_bucket'], s3_key, Callback=print)
    print(f"Uploaded plot png to S3: {s3_key}")

def print_time(msg):
    time_now = datetime.now()
    date_time = time_now.strftime("%m/%d/%Y, %H:%M:%S")
    print("===%s: %s" % (msg, time_now))


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
        print('ERROR: ',str(err))
        return {
            "success": False,
            "message": str(err)
        }

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

        print('---adata shape after gating----')
        print(adata.shape)
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
        print('ERROR: ',str(err))
        return {
            "success": False,
            "message": str(err)
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
#               "message": str(err)}

@app.post("/shutdown")
async def shutdown():
    os.kill(os.getpid(), signal.SIGTERM)
    return Response(status_code=200, content='Server shutting down...')

@app.get("/health", status_code = 200)
async def health():
    return {"status": "ok"}