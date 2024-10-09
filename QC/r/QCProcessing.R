library(Seurat)
library(Matrix)
library(dplyr)
library(jsonlite)
library(DoubletFinder)
library(aws.s3)

user_environment <- new.env()
environment <- Sys.getenv("ENVIRONMENT")
if (environment == "") {
  environment <- "dev"
}
message("Cellborg QCProcessing R container running in environment: ", environment)

if (environment == "dev") {
  DATASET_BUCKET <- "cellborgdatasetuploadbucket"
  QC_DATASET_BUCKET <- "cellborgqcdatasetbucket"
} else {
  DATASET_BUCKET <- paste0("cellborg-", environment, "-datasetupload-bucket")
  QC_DATASET_BUCKET <- paste0("cellborg-", environment, "-qcdataset-bucket")
}
user_environment$dataset_bucket <- DATASET_BUCKET
user_environment$qc_dataset_bucket <- QC_DATASET_BUCKET

perform_qc <- function(user, project, dataset, min, max, mt) {

  #Load dataset and create Seurat object
  message("Loading dataset and creating seurat object")
  prefix <- paste0(user, "/", project, "/", dataset, "/")
  message(prefix)
  temp_dir <- tempdir()
  files <- get_bucket(
    bucket = user_environment$dataset_bucket,
    prefix = prefix
  )
 

  for (file in files) {
    file_name <- basename(file$Key)
    message(file_name)
    save_object(object = file$Key, bucket = user_environment$dataset_bucket,
                file = file.path(temp_dir, file_name))
  }
  
  data <- Read10X(data.dir = temp_dir)
  seurat_obj <- CreateSeuratObject(
    counts = data, 
    project = dataset, 
    min.cells = 3, 
    min.features = 200)

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  message("Subsetting based on input: ", min, max, mt)
  seurat_obj2 <- subset(seurat_obj, subset = nFeature_RNA > min & nFeature_RNA < max & percent.mt < mt)  # nolint
  with(seurat_obj2, {
    nrna <- seurat_obj2$nFeature_RNA
    nrn <- seurat_obj2$nCount_RNA
    pmt <- seurat_obj2$percent.mt
  })

   rm(seurat_obj)

  ### Write data in JSON format
  json_data <- list(
    nCount_RNA = nrn,
    nFeature_RNA = nrna,
    percent.mt = pmt
  )
  json_file <- tempfile(fileext = ".json")
  plot_data <- toJSON(json_data)
  write(plot_data, json_file)
 

 
  cat("Uploading the json data for qc violin plot to S3...")
  ### Create S3 Object key for quality control data
  s3_key <- paste0(user, "/", project, "/", dataset, "/plots/QcViolinPlotData.json")
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$qc_dataset_bucket,
    show_progress = TRUE,
    multipart = TRUE
  )
  message("Uploaded qc violin plot data to S3: ", s3_key)
  file.remove(json_file)
  
  rm(json_data,plot_data,nrn,nrna,pmt)
  do_gc()

  message("Performing pre-processing of dataset")
  data1 <- NormalizeData(seurat_obj2, 
            normalization.method = "LogNormalize",
            scale.factor = 10000)
  
  rm(seurat_obj2)
  
  data2 <- FindVariableFeatures(data1, selection.method = "vst", nfeatures = 2000)

  rm(data1)

  message("Found variable features")
  json_data <- variable_feature_data(data2)

  plot_data <- toJSON(json_data)

  rm(json_data)
  do_gc()

  json_file <- tempfile(fileext = ".json")
  write(plot_data, json_file)
  ### Create S3 Object key for variable feature data
  s3_key <- paste0(user, "/", project, "/", dataset, "/plots/variableFeatureData.json") # nolint: line_length_linter.
  message(s3_key)

  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$qc_dataset_bucket,
    show_progress = TRUE,
    multipart = TRUE
  )
  message("Uploaded var feature plot data to S3")
  file.remove(json_file)

  rm(plot_data)
  
  all.genes <- rownames(data2)
  data3 <- ScaleData(data2, features = all.genes)
  
  rm(data2)
  do_gc(5)

  varfeat <- VariableFeatures(object = data3)
  data4 <- RunPCA(data3, features = varfeat)
 
  rm(data3)
  do_gc(5)
  
  data5 <- FindNeighbors(data4, dims = 1:10)
  data6 <- FindClusters(data5, resolution = 0.5)

  data7 <- RunUMAP(data6, dims = 1:10)
 
  rm(data4,data3,data5,data6)
  do_gc(10)

  gene_count <- dim(data7)[1]
  cell_count <- dim(data7)[2]
  
  message("Preprocessing complete")
  message("Performing pK identification")
  sweep.res.young <- paramSweep(data7, PCs = 1:10, sct = FALSE)
  sweep.stats_young <- summarizeSweep(sweep.res.young, GT = FALSE)
  bcmvn_young <- find.pK(sweep.stats_young)
  attributes(bcmvn_young)
  m <- bcmvn_young[bcmvn_young$BCmetric == max(bcmvn_young$BCmetric), ]

  result <- m$pK
  result <- as.numeric(result)
 

  doubletrate <- (length(data7$orig.ident) / 1000) * 0.008
  annotations <- data7@meta.data$ClusteringResults


  message("Performing Homotypic Doublet Proportion Estimation")
  homotypic.prop <- modelHomotypic(annotations)

  rm(annotations)
  do_gc()

  nExp_poi <- round(doubletrate*nrow(data7@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  all_singlets <- doubletFinder(data7, PCs = 1:10, pN = 0.25, pK = result, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE,10000) # nolint: line_length_linter.
  
  rm(data7)
  do_gc()
  
  column <- grep(paste0("^","DF.classifications"), colnames(all_singlets@meta.data), value = TRUE) # nolint
  Idents(object = all_singlets) <- column

  singlets <- subset(all_singlets, idents = c("Singlet"), invert = FALSE)
  doublets <- subset(all_singlets, idents = c("Doublet"), invert = FALSE)
  
  rm(all_singlets)
  do_gc()

  singlet_data <- Embeddings(object = singlets, reduction = "umap")
  df_singlets <- as.data.frame(singlet_data)
  df_singlets$class <- "Singlet"
  rm(singlet_data) 
  do_gc()
  

  doublet_data <- Embeddings(object = doublets, reduction = "umap")
  df_doublets <- as.data.frame(doublet_data)
  df_doublets$class <- "Doublet"
  
  rm(doublets,doublet_data)
  do_gc()
  
  combined_df <- rbind(df_singlets, df_doublets)
  
  rm(df_singlets,df_doublets)
  do_gc()

  json_file <- tempfile(fileext = ".json")
  write_json(combined_df, json_file)
  s3_key <- paste0(user, "/", project, "/", dataset, "/plots/doubletplot.json")
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$qc_dataset_bucket,
  )
  message("Uploaded doublet dim plot to s3: ", s3_key)
  unlink(json_file)
  do_gc()
  tmp_file <- tempfile(fileext = ".rds")
  saveRDS(singlets, file = tmp_file)
  s3_key <- paste0(user, "/", project, "/", dataset, "/", dataset, ".rds")
  put_object(
    file = tmp_file,
    object = s3_key,
    bucket = user_environment$qc_dataset_bucket,
    multipart = TRUE,
    verbose = TRUE
  )

  message("Uploaded QC processed dataset to S3")
  unlink(tmp_file)
  message("COMPLETE")
  return_list <- list("cell_count" = cell_count, "gene_count" = gene_count)
  return(return_list)
}

variable_feature_data <- function(
  object,
  selection.method = NULL,
  assay = NULL
) {
  hvf_info <- HVFInfo(
    object = object,
    assay = assay,
    selection.method = selection.method,
    status = TRUE
  )
  var_status <- c('no', 'yes')[unlist(x = hvf_info[, ncol(x = hvf_info)]) + 1] # nolint
  if (colnames(x = hvf_info)[3] == "dispersion.scaled") {
    hvf_info <- hvf_info[, c(1, 2)]
  } else {
    hvf_info <- hvf_info[, c(1, 3)]
  }
  axis_labels <- switch( # nolint
    EXPR = colnames(x = hvf_info)[2],
    "variance.standardized" = c("Average Expression", "Standardized Variance"),
    "dispersion" = c("Average Expression", "Dispersion"),
    "residual_variance" = c("Geometric Mean of Expression", "Residual Variance")
  )

  # Extract the data you need for the plot
  mean_expression <- hvf_info[, 1]
  standardized_variance <- hvf_info[, 2]

  # Create a data frame with the extracted data
  variable_feature_data <- data.frame(
    Mean_Expression = mean_expression,
    Standardized_Variance = standardized_variance
  )
  rm(hvf_info)
  do_gc()

  return(variable_feature_data)
}

do_gc <- function(sleep_seconds=10){
  message("performing gc...")
  print("========objects in memory=======")
  sort( sapply(ls(),function(x){object.size(get(x))}))
  print("========gc begin")

  gc()
  Sys.sleep(sleep_seconds)

  print("========gc done")
  print("========objects in memory=======")
  sort( sapply(ls(),function(x){object.size(get(x))}))
}

#* @post /qc_endpoint
#* @parser json
function(req, res) {
  req_data <- req$body
  dataset <- req_data$dataset
  user <- req_data$user
  project <- req_data$project
  min <- req_data$min
  max <- req_data$max 
  mt <- req_data$mt
  message("Performing QC...")
  return_data <- perform_qc(user, project, dataset, min, max, mt)
  cell_count <- return_data$cell_count
  gene_count <- return_data$gene_count
  message("Successfully performed QC on the dataset")
  res$status <- 200
  return(list(
    success = TRUE,
    message = "QC Completed Successfully",
    cell_count = cell_count,
    gene_count = gene_count
  ))
}

#* @post /shutdown
#* @parser json
function(req, res) {
  message("Shutting down the QC API server...")
  quit(save = "no", status = 0, runLast = TRUE)
}

