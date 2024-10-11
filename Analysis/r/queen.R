library(Seurat)
library(Matrix)
library(dplyr)
library(aws.s3)
library(jsonlite)
library(harmony)
library(RProtoBuf)
library(monocle3)
library(SeuratWrappers)

user_environment <- new.env()
Sys.setenv("AWS_ACCESS_KEY_ID" = "AKIA4TIFZQZ254WVY4YO",
           "AWS_SECRET_ACCESS_KEY" = "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",
           "AWS_DEFAULT_REGION" = "us-west-2")
environment <- Sys.getenv("ENVIRONMENT")
if (environment == "") {
  environment <- "dev"
}
message("Cellborg Queen R container running in environment: ", environment)

if (environment == "dev") {
  DATASET_BUCKET <- "cellborgdatasetuploadbucket"
  QC_DATASET_BUCKET <- "cellborgqcdatasetbucket" #Seurat object
  VAR_FEAT_BUCKET <- "cellborgvariablefeaturebucket"
  CLUSTER_BUCKET <- "cellborgdatasetclusterplotbucket"
  FEATURE_PLOT_BUCKET <- "cellborggenefeatureplotbucket"
  GENE_NAMES_BUCKET <- "cellborgprojectgenenamesbucket"
  PCA_BUCKET <- "cellborgpcabucket"
  HEATMAP_BUCKET <- "cellborgheatmapbucket"
  PSUEDOTIME_BUCKET <- "cellborgpsuedotimebucket"
  DOTPLOT_BUCKET <- "cellborgdotplotbucket"
  VLNPLOTS_BUCKET <- "cellborgviolinplotsbucket"
} else {
  DATASET_BUCKET <- paste0("cellborg-", environment, "-datasetupload-bucket")
  QC_DATASET_BUCKET <- paste0("cellborg-", environment, "-qcdataset-bucket")
  VAR_FEAT_BUCKET <- paste0("cellborg-", environment, "-variablefeature-bucket")
  CLUSTER_BUCKET <- paste0("cellborg-", environment, "-datasetclusterplot-bucket")
  FEATURE_PLOT_BUCKET <- paste0("cellborg-", environment, "-genefeatureplot-bucket")
  GENE_NAMES_BUCKET <- paste0("cellborg-", environment, "-projectgenenames-bucket")
  PCA_BUCKET <- paste0("cellborg-", environment, "-pca-bucket")
  HEATMAP_BUCKET <- paste0("cellborg-", environment, "-heatmap-bucket")
  PSUEDOTIME_BUCKET <- paste0("cellborg-", environment, "-psuedotime-bucket")
  DOTPLOT_BUCKET <- paste0("cellborg-", environment, "-dotplot-bucket")
  VLNPLOTS_BUCKET <- paste0("cellborg-", environment, "-violinplots-bucket")
}
user_environment$dataset_bucket <- DATASET_BUCKET
user_environment$qc_dataset_bucket <- QC_DATASET_BUCKET
user_environment$var_feature_bucket <- VAR_FEAT_BUCKET
user_environment$cluster_bucket <- CLUSTER_BUCKET
user_environment$feature_plot_bucket <- FEATURE_PLOT_BUCKET
user_environment$gene_names_bucket <- GENE_NAMES_BUCKET
user_environment$pca_bucket <- PCA_BUCKET
user_environment$heatmap_bucket <- HEATMAP_BUCKET
user_environment$psuedotime_bucket <- PSUEDOTIME_BUCKET
user_environment$dotplot_bucket <- DOTPLOT_BUCKET
user_environment$vlnplots_bucket <- VLNPLOTS_BUCKET

init_dataset_seurat_object <- function(user, project, datasets, analysis) {

  message("Initializing the seurat object...")

  message("Retrieving datasets from S3 for: ", project)
  # Create an empty Seurat object to store merged data
  merged_seurat <- list()

  # Download datasets from S3 and merge
  for (dataset in datasets) {
    key <- paste0(user, "/", project, "/", dataset, "/", dataset, ".rds")
    message("Retrieving: ", key)
    temp_file <- tempfile()
    save_object(
      object = key,
      bucket = user_environment$qc_dataset_bucket,
      file = temp_file,
      verbose = TRUE
    )
    # Load and Merge the .rds Files
    seurat_obj <- readRDS(temp_file)
    merged_seurat[[length(merged_seurat) + 1]] <- seurat_obj

    # Remove the temporary file
    unlink(temp_file)
  }

  # Merge the Seurat objects
  if (length(merged_seurat) > 1) {
    merged_object <- merge(merged_seurat[[1]], y = merged_seurat[-1])
  } else {
    merged_object <- merged_seurat[[1]]
  }
  message("Datasets received and merged!")
  return(merged_object)
}

variable_feature_search <- function(data, user, project, analysis, nFeatures) {

  variable_data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = nFeatures)

  message("Found variable features")
  json_data <- variable_feature_data(variable_data)
  plot_data <- jsonlite::toJSON(json_data)
  json_file <- tempfile(fileext = ".json")
  write(plot_data, json_file)

  ### Create S3 Object key for variable feature data
  s3_key <- paste0(user, "/", project, "/", analysis, "/variableFeatures.json")
  message(s3_key)

  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$var_feature_bucket
  )
  message("Uploaded to S3")
  file.remove(json_file)
  scaled_data <- perform_pca(variable_data, user, project, analysis)
  return(scaled_data)
}

perform_pca <- function(data, user, project, analysis) {
  data <- ScaleData(data)
  data <- RunPCA(data, features = VariableFeatures(object = data))

  s3_key <- paste0(user, "/", project, "/", analysis, "/embeddings.json")
  json_file <- tempfile(fileext = ".json")
  pca_scores <- Embeddings(data, reduction = "pca")
  jsonlite::write_json(pca_scores, json_file) #UPLOAD TO S3
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$pca_bucket
  )
  message("Uploaded pca embeddings to s3: ", s3_key)

  s3_key <- paste0(user, "/", project, "/", analysis, "/loadings.json")
  tmp_file <- tempfile(fileext = ".json")
  loadings <- Loadings(data[["pca"]])[, 1:2]
  df <- as.data.frame(loadings)
  df$gene_name <- rownames(df)
  jsonlite::write_json(df, tmp_file) #UPLOAD TO S3
  put_object(
    file = tmp_file,
    object = s3_key,
    bucket = user_environment$pca_bucket
  )
  message("Uploaded pca loadings to s3: ", s3_key)
  file.remove(json_file)
  file.remove(tmp_file)
  return(data)
}

perform_clustering <- function(
    data,
    user,
    project,
    analysis,
    cluster_neighbors,
    cluster_clusters,
    cluster_dimensions,
    cluster_reduction) {

  cluster_neighbors <- cluster_neighbors
  cluster_clusters <- as.double(cluster_clusters)
  cluster_dimensions <- cluster_dimensions
  cluster_reduction <- as.character(cluster_reduction)

  ### PROCESS DATA
  s3_key <- paste0(
                   user,
                   "/",
                   project,
                   "/",
                   analysis,
                   "/",
                   "n=", cluster_neighbors,
                   "&c=", cluster_clusters,
                   "&d=", cluster_dimensions,
                   "&r=", cluster_reduction,
                   ".json",
                   sep = "")

  message(s3_key)

  if (length(unique(data@meta.data$orig.ident)) > 1) { #need to use harmony
    data <- RunHarmony(data, group.by.vars = "orig.ident")
    data <- FindNeighbors(data, dims = 1:cluster_neighbors, reduction = "harmony") # nolint: line_length_linter.
    data <- FindClusters(data, resolution = cluster_clusters, reduction = "harmony") # nolint: line_length_linter.
    if (cluster_reduction == "umap") {
      data <- RunUMAP(data, dims = 1:cluster_dimensions, reduction = "harmony")
    } else if (cluster_reduction == "pca") {
      data <- RunPCA(data, dims = 1:cluster_dimensions, reduction = "harmony")
    } else if (cluster_reduction == "tsne") {
      data <- RunTSNE(data, dims = 1:cluster_dimensions, reduction = "harmony")
    } else {
      message("Error: Invalid Reduction type provided. ", cluster_reduction)
    }
  } else { #Harmony is not used
    data <- FindNeighbors(data, dims = 1:cluster_neighbors)
    data <- FindClusters(data, resolution = cluster_clusters) # nolint: line_length_linter.
    if (cluster_reduction == "umap") {
      data <- RunUMAP(data, dims = 1:cluster_dimensions)
    } else if (cluster_reduction == "pca") {
      data <- RunPCA(data, dims = 1:cluster_dimensions)
    } else if (cluster_reduction == "tsne") {
      data <- RunTSNE(data, dims = 1:cluster_dimensions)
    } else {
      message("Error: Invalid Reduction type provided. ", cluster_reduction)
    }
  }

  ### Extracting plotting data from processed Seurat object
  results <- data[[cluster_reduction]]
  cluster_cells <- as.numeric(table(data$seurat_clusters))

  # Extract cell coordinates (first two dimensions for a 2D plot)
  cell_coordinates <- results@cell.embeddings[, 1:2]
  cell_metadata <- data@meta.data
  cluster_assignments <- cell_metadata$seurat_clusters
  num_clusters <- length(unique(cluster_assignments))

  data_for_json <- list(
    total_clusters = num_clusters,
    clusterCounts = cluster_cells,
    cluster = as.numeric(cluster_assignments),
    X = cell_coordinates[, 1],
    Y = cell_coordinates[, 2]
  )

  ### UPLOAD JSON TO S3
  temp_file <- tempfile(fileext = ".json")
  jsonlite::write_json(data_for_json, temp_file)

  put_object(
    file = temp_file,
    object = s3_key,
    bucket = user_environment$cluster_bucket
  )

  file.remove(temp_file)
  getGeneNamesAndUploadToS3(data, user, project, analysis)

  data.cds <- as.cell_data_set(data)
  data.cds <- cluster_cells(cds = data.cds, reduction_method = "UMAP")
  data.cds <- learn_graph(data.cds, use_partition = TRUE, close_loop = FALSE)
  order_cells_json <- order_up_cells(data.cds, reduction_method = "UMAP")

  json_temp <- tempfile(fileext = ".json")
  write(order_cells_json, json_temp)
  ocs3_key <- paste0(user, "/", project, "/", analysis, "/order_cells.json", sep = "")
  message("Uploading to s3: ", ocs3_key)
  put_object(
    file = json_temp,
    object = ocs3_key,
    bucket = user_environment$psuedotime_bucket
  )
  message("Uploaded order cells plot to s3")
  file.remove(json_temp)
  return_list <- list("seurat" = data, "monocle" = data.cds)
  return(return_list)
}

get_gene_feature_plot_data <- function(data, gene_name, user, project, analysis) { # nolint: line_length_linter.
  message("Processing feature plot")
  plot_data <- feature_plot_data(data, features = gene_name)
  json_file <- tempfile(fileext = ".json")
  jsonlite::write_json(plot_data, json_file)

  s3_key <- paste0(user, "/", project, "/", analysis, "/", gene_name, "/", "featurePlot.json") # nolint: line_length_linter.
  message(s3_key)
  message("Uploading feature plot data to s3")
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$feature_plot_bucket
  )
  message("Feature Plot data uploaded to S3")
  file.remove(json_file)
}

feature_plot_data <- function(
  object,
  features,
  dims = c(1, 2),
  cells = NULL,
  cols = if (blend) {
    c('lightgrey', '#ff0000', '#00ff00')
  } else {
    c('lightgrey', 'blue')
  },
  pt.size = NULL,
  order = FALSE,
  min.cutoff = NA,
  max.cutoff = NA,
  reduction = NULL,
  split.by = NULL,
  keep.scale = "feature",
  shape.by = NULL,
  slot = 'data',
  blend = FALSE,
  blend.threshold = 0.5,
  label = FALSE,
  label.size = 4,
  label.color = "black",
  repel = FALSE,
  ncol = NULL,
  coord.fixed = FALSE,
  by.col = TRUE,
  sort.cell = NULL,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
) {
  # TODO: deprecate fully on 3.2.0
  if (!is.null(x = sort.cell)) {
    warning(
      "The sort.cell parameter is being deprecated. Please use the order ",
      "parameter instead for equivalent functionality.",
      call. = FALSE,
      immediate. = TRUE
    )
    if (isTRUE(x = sort.cell)) {
      order <- sort.cell
    }
  }
  # Check keep.scale param for valid entries
  if (!(is.null(x = keep.scale)) && !(keep.scale %in% c("feature", "all"))) {
    stop("`keep.scale` must be set to either `feature`, `all`, or NULL")
  }
  # Get the DimReduc to use
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  # Name the reductions
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident', features),
    cells = cells,
    slot = slot
  )
  return(data)
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

  return(variable_feature_data)
}

plot_cells_and_upload_to_s3 <- function(cds,
                       x=1,
                       y=2,
                       reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                       color_cells_by="cluster",
                       group_cells_by="cluster",
                       genes=NULL,
                       show_trajectory_graph=TRUE,
                       user,
                       project,
                       analysis,
                       trajectory_graph_color="grey28",
                       trajectory_graph_segment_size=0.75,
                       norm_method = c("log", "size_only"),
                       label_cell_groups = TRUE,
                       label_groups_by_cluster=TRUE,
                       group_label_size=2,
                       labels_per_group=1,
                       label_branch_points=TRUE,
                       label_roots=TRUE,
                       label_leaves=TRUE,
                       graph_label_size=2,
                       cell_size=0.35,
                       cell_stroke= I(cell_size / 2),
                       alpha = 1,
                       min_expr=0.1,
                       rasterize=FALSE,
                       scale_to_range=TRUE,
                       label_principal_points = FALSE) {
  feature_label <- min_val_for_feature <- max_val_for_feature <- NULL # no visible binding
  cell_group <- cell_color <- cells_in_cluster <- per <- NULL # no visible binding
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  low_dim_coords <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >=max(x,y),
                          msg = paste("x and/or y is too large. x and y must",
                                      "be dimensions in reduced dimension",
                                      "space."))
  if(!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                  "pseudotime") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must one of",
                                        "'cluster', 'partition', 'pseudotime,",
                                        "or a column in the colData table."))

    if(color_cells_by == "pseudotime") {
      tryCatch({pseudotime(cds, reduction_method = reduction_method)},
               error = function(x) {
                 stop("No pseudotime for ", reduction_method,
                            " calculated. Please run order_cells with ",
                            "reduction_method = ", reduction_method,
                            " before attempting to color by pseudotime.")})

    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))

  norm_method = match.arg(norm_method)

  assertthat::assert_that((group_cells_by %in% c("cluster", "partition")) ||
                          (group_cells_by %in% names(SummarizedExperiment::colData(cds))),
                          msg = paste("group_cells_by must be one of",
                          "'cluster', 'partition', or a column in the colData table."))

  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))

  if (show_trajectory_graph &&
      is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }
  if (label_principal_points &&
      is.null(principal_graph(cds)[[reduction_method]])) {
    message("Cannot label principal points when no trajectory to plot. Has learn_graph() been called yet?")
    label_principal_points = FALSE
  }

  if (label_principal_points) {
    label_branch_points <- FALSE
    label_leaves <- FALSE
    label_roots <- FALSE
  }


  gene_short_name <- NA
  sample_name <- NA
  #sample_state <- colData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA

  S_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[,c(x,y)])

  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)

  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster") {
    data_df$cell_group <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    data_df$cell_group <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else{
    data_df$cell_group <- SummarizedExperiment::colData(cds)[data_df$sample_name, group_cells_by]
  }

  if (color_cells_by == "cluster"){
    data_df$cell_color <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "pseudotime") {
    data_df$cell_color <-
      tryCatch({pseudotime(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]}, error = function(e) {NULL})
  } else{
    data_df$cell_color <- colData(cds)[data_df$sample_name, color_cells_by]
  }

  ## Graph info
  if (show_trajectory_graph) {

    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))

    dp_mst <- cds@principal_graph[[reduction_method]]

    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(
                           source="sample_name",
                           source_prin_graph_dim_1="prin_graph_dim_1",
                           source_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(
                           target="sample_name",
                           target_prin_graph_dim_1="prin_graph_dim_1",
                           target_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "target")
  }

  ## Marker genes
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes)[[2]] >= 2){
      markers = unlist(genes[,1], use.names=FALSE)
    } else {
      markers = genes
    }
	markers_rowData <- rowData(cds)[(rowData(cds)$gene_short_name %in% markers) |
							        (row.names(rowData(cds)) %in% markers),,drop=FALSE]
	markers_rowData <- as.data.frame(markers_rowData)
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
#      assertthat::assert_that(!is.null(size_factors(cds_exprs)))
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))

      if (!is.null(dim(genes)) && dim(genes)[[2]] >= 2){
        #genes = as.data.frame(genes)
        #row.names(genes) = genes[,1]
        #genes = genes[row.names(cds_exprs),]

        agg_mat = as.matrix(aggregate_gene_expression(cds, genes, norm_method=norm_method, gene_agg_fun="mean", scale_agg_values=FALSE))
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        if (is.factor(genes[,2]))
          markers_exprs$feature_id = factor(markers_exprs$feature_id,
                                            levels=levels(genes[,2]))

        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method = "size_only"
        expression_legend_label = "Expression score"
      } else {
        cds_exprs@x = round(10000*cds_exprs@x)/10000
        markers_exprs = matrix(cds_exprs, nrow=nrow(markers_rowData))
        colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_rowData,
                               by.x = "feature_id", by.y="row.names")
        if (is.null(markers_exprs$gene_short_name)) {
          markers_exprs$feature_label <-
            as.character(markers_exprs$feature_id)
        } else {
          markers_exprs$feature_label <-
            as.character(markers_exprs$gene_short_name)
        }

        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | !as.character(markers_exprs$feature_label) %in% markers,
                                              as.character(markers_exprs$feature_id),
                                              as.character(markers_exprs$feature_label))

        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                              levels = markers)
        if (norm_method == "size_only")
          expression_legend_label = "Expression"
        else
          expression_legend_label = "log10(Expression)"
      }

      if (scale_to_range){
				markers_exprs = dplyr::group_by(markers_exprs, feature_label) %>%
				  dplyr::mutate(max_val_for_feature = max(value),
								min_val_for_feature = min(value)) %>%
				  dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
				expression_legend_label = "% Max"
			  }
			}
		  }

		  if (label_cell_groups && is.null(color_cells_by) == FALSE){
			if (is.null(data_df$cell_color)){
			  if (is.null(genes)){
				message(color_cells_by, " not found in colData(cds), cells will ",
							  "not be colored")
			  }
			  text_df = NULL
			  label_cell_groups = FALSE
			}else{
			  if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {

				if (label_groups_by_cluster && is.null(data_df$cell_group) == FALSE){
				  text_df = data_df %>%
					dplyr::group_by(cell_group) %>%
					dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
					dplyr::group_by(cell_color, .add=TRUE) %>%
					dplyr::mutate(per=dplyr::n()/cells_in_cluster)
				  median_coord_df = text_df %>%
					dplyr::summarize(fraction_of_group = dplyr::n(),
									 text_x = stats::median(x = data_dim_1),
									 text_y = stats::median(x = data_dim_2))
				  text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
											   dplyr::distinct())
				  text_df = suppressMessages(dplyr::inner_join(text_df,
															   median_coord_df))
				  text_df = text_df %>% dplyr::group_by(cell_group) %>%
					dplyr::top_n(labels_per_group, per)
				} else {
				  text_df = data_df %>% dplyr::group_by(cell_color) %>%
					dplyr::mutate(per=1)
				  median_coord_df = text_df %>%
					dplyr::summarize(fraction_of_group = dplyr::n(),
									 text_x = stats::median(x = data_dim_1),
									 text_y = stats::median(x = data_dim_2))
				  text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
											   dplyr::distinct())
				  text_df = suppressMessages(dplyr::inner_join(text_df,
															   median_coord_df))
				  text_df = text_df %>% dplyr::group_by(cell_color) %>%
					dplyr::top_n(labels_per_group, per)
				}
				text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
			  } else {
				message("Cells aren't colored in a way that allows them to ",
                        "be grouped.")
				text_df = NULL
				label_cell_groups = FALSE
			  }
			}
		  }
		  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0){
			data_df <- merge(data_df, markers_exprs, by.x="sample_name",
							 by.y="cell_id")
			data_df$value <- with(data_df, ifelse(value >= min_expr, value, NA))
			ya_sub <- data_df[!is.na(data_df$value),]
			na_sub <- data_df[is.na(data_df$value),]
			if(norm_method == "size_only"){
        print("------------------1----------------")
    } else {
      print("------------------2----------------")
      
      if (scale_to_range){
        print("------------------3----------------")
        
      }else{
        print("------------------4----------------")

      }
      print("------------------5----------------")
      
    }
  } else {
    print("------------------6----------------")
    
    if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        print("------------------7----------------")
        
        message("cluster_cells() has not been called yet, can't ",
                      "color cells by cluster")
      } else{
        print("------------------8----------------")
      }
      print("------------------9----------------")

    } else if (methods::is(data_df$cell_color, "numeric")) {
      print("------------------10----------------")

    } else {
      print("------------------11----------------")

    }

  }
  if (show_trajectory_graph){
    print("------------------12----------------")

    if (label_principal_points) {
      print("------------------13----------------")
    }
    if (label_branch_points){
      print("------------------14----------------")
    }
    if (label_leaves){
      print("------------------15----------------")
    }
    if (label_roots){
      print("------------------16----------------")
    }
  }

  if(label_cell_groups) {
    print("------------------17----------------")

    # If we're coloring by gene expression, don't hide the legend
    if (is.null(markers_exprs))
    print("------------------18----------------")

  }
  print("------------------19----------------")

  json_file <- tempfile(fileext = ".json")
  data_list <- list(
    data = data_df,
    edges = edge_df
  )
  jsonlite::write_json(data_list, json_file)
  s3_key <- paste0(user, "/", project, "/", analysis, "/psuedotime.json")
  message(s3_key)
  message("Uploading psuedotime plot data to s3")
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$psuedotime_bucket
  )
  message("upload to S3 completed")
  file.remove(json_file)
}

select_trajectory_roots <- function(cds, x=1, y=2, # nocov start
                                    reduction_method) {
  prin_graph_dim_1 <- prin_graph_dim_2 <- V1 <- V2 <- NULL # no visible binding
  source_prin_graph_dim_1 <- target_prin_graph_dim_1 <- NULL # no visible binding
  source_prin_graph_dim_2 <- target_prin_graph_dim_2 <- NULL # no visible binding
  reduced_dim_coords <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst)

  ica_space_df <- as.data.frame(reduced_dim_coords)
  reduced_dims <- as.data.frame(SingleCellExperiment::reducedDims(cds)[[reduction_method]])

  colnames(reduced_dims)<-vapply(seq_along(colnames(reduced_dims)),function(i){paste0('V',i)},c('a'))
  if (is.null(ica_space_df) || nrow(ica_space_df) == 0) {
    stop("ica_space_df is NULL or has no rows.")
  }
  num_reduced_dim <- ncol(ica_space_df)
  print(num_reduced_dim)
  if( num_reduced_dim >= 3 )
  {
    if( num_reduced_dim > 3 )
    {
      message( reduction_method, ' space has ', num_reduced_dim, ' dimensions but is shown in three.' )
    }
    use_3d=TRUE
  }
  else
  {
    use_3d=FALSE
  }
  print(use_3d)
  colnames(ica_space_df) <- vapply(seq_along(ica_space_df),function(i){paste0('prin_graph_dim_',i)},c('a'))

  ica_space_df$sample_name <- row.names(ica_space_df)
  ica_space_df$sample_state <- row.names(ica_space_df)

  dp_mst <- principal_graph(cds)[[reduction_method]]

  if (is.null(dp_mst)){
    stop("You must first call orderCells() before using this function")
  }

  if (use_3d){
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(source="sample_name",
                                    source_prin_graph_dim_1="prin_graph_dim_1",
                                    source_prin_graph_dim_2="prin_graph_dim_2",
                                    source_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(target="sample_name",
                                    target_prin_graph_dim_1="prin_graph_dim_1",
                                    target_prin_graph_dim_2="prin_graph_dim_2",
                                    target_prin_graph_dim_3="prin_graph_dim_3"),
                       by = "target")
  }else{
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(source="sample_name",
                                    source_prin_graph_dim_1="prin_graph_dim_1",
                                    source_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select(target="sample_name",
                                    target_prin_graph_dim_1="prin_graph_dim_1",
                                    target_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "target")
  }

  num_roots <- nrow(ica_space_df)
  sel <- rep(FALSE, nrow(ica_space_df))

  if (use_3d){
    ica_space_df$keep <- FALSE
    keeprows <- rep(TRUE, nrow(ica_space_df))
    ica_space_df[keeprows,]$keep <- TRUE
    
    if(sum(!ica_space_df$keep) == 0) {
      cols <- c("black")
    } else {
      cols <- c("red", "black")
    }
    data_list <- list(
        ica_space_data = list(
            x = ica_space_df$prin_graph_dim_1,
            y = ica_space_df$prin_graph_dim_2,
            z = ica_space_df$prin_graph_dim_3,
            key = ica_space_df$sample_name,
            color = ica_space_df$keep,
            colors = cols
        ),
        reduced_dims_data = list(
            x = reduced_dims$V1,
            y = reduced_dims$V2,
            z = reduced_dims$V3
        )
    )
    json_data <- jsonlite::toJSON(data_list, dataframe = "rows")
    return(json_data)

  } else {

    keeprows <- rep(TRUE, nrow(ica_space_df))
    keep    <- ica_space_df[keeprows, , drop = FALSE]
    exclude <- ica_space_df[!keeprows, , drop = FALSE]
    data_list <- list(
        keep = keep,
        exclude = exclude,
        reduced_dims = reduced_dims,
        edge_df = edge_df
    )
    json_data <- jsonlite::toJSON(data_list, dataframe = "rows")
    return(json_data)
  }
}

order_up_cells <- function(cds,
                        reduction_method = "UMAP",
                        root_pr_nodes=NULL,
                        root_cells=NULL,
                        verbose = FALSE){

  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::are_equal("UMAP", reduction_method),
                          msg = paste("Currently only 'UMAP' is accepted as a",
                                      "reduction_method."))
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste0("No dimensionality reduction for ",
                                      reduction_method, " calculated. ",
                                      "Please run reduce_dimension with ",
                                      "reduction_method = ", reduction_method,
                                      ", cluster_cells, and learn_graph ",
                                      "before running order_cells."))
  assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]),
                          msg = paste("No cell clusters for",
                                      reduction_method, "calculated.",
                                      "Please run cluster_cells with",
                                      "reduction_method =", reduction_method,
                                      "and run learn_graph before running",
                                      "order_cells."))
  assertthat::assert_that(!is.null(principal_graph(cds)[[reduction_method]]),
                          msg = paste("No principal graph for",
                                      reduction_method, "calculated.",
                                      "Please run learn_graph with",
                                      "reduction_method =", reduction_method,
                                      "before running order_cells."))
  assertthat::assert_that(
    igraph::vcount(principal_graph(cds)[[reduction_method]]) < 10000,
    msg = paste("principal graph is too large. order_cells doesn't support",
                "more than 10 thousand centroids."))
  if(!is.null(root_pr_nodes)) {
    assertthat::assert_that(
      all(root_pr_nodes %in%
            igraph::V(principal_graph(cds)[[reduction_method]])$name),
      msg = paste("All provided root_pr_nodes must be present in the",
                  "principal graph."))
  }

  if(!is.null(root_cells)) {
    assertthat::assert_that(all(root_cells %in% row.names(colData(cds))),
                            msg = paste("All provided root_cells must be",
                                        "present in the cell data set."))
  }
  assertthat::assert_that(!all(c(!is.null(root_cells),
                                 !is.null(root_pr_nodes))),
                          msg = paste("Please specify either root_pr_nodes",
                                      "or root_cells, not both."))

  if (is.null(root_pr_nodes) & is.null(root_cells)){
    root_pr_nodes <- select_trajectory_roots(cds, reduction_method = reduction_method)
    return (root_pr_nodes)
  } else if (!is.null(root_cells)) {
    closest_vertex <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex
    root_pr_nodes <- unique(paste("Y_", closest_vertex[root_cells,], sep=""))
    return (root_pr_nodes)
  }
}

getGeneNamesAndUploadToS3 <- function(data, user, project, analysis) { # nolint
  gene_names <- rownames(data)
  gene_names_data <- list(geneNames = gene_names)
  json_file <- tempfile(fileext = ".json")
  jsonlite::write_json(gene_names_data, json_file)

  s3_key <- paste0(user, "/", project, "/", analysis, "/geneNames.json")
  message(s3_key)
  message("Uploading to s3")
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$gene_names_bucket
  )
  message("upload to S3 completed")
  file.remove(json_file)
}

DoTheHeatmap <- function(
  object,
  features = NULL,
  cells = NULL,
  group.by = 'ident',
  slot = 'scale.data',
  assay = NULL
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  
  # make sure features are present
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot,
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(
    object = object,
    slot = slot)[features, cells, drop = FALSE])))
  object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
  group.by <- group.by %||% 'ident'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  data$Cluster <- groups.use

  heatmap_data_list <- list(
    heatmapValues = data
  )
  return(heatmap_data_list)
}

DotPlotData <- function(
  object,
  assay = NULL,
  features,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = 'radius',
  scale.min = NA,
  scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))

  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  return(data.plot)
}

ExIPlot <- function(
  object,
  features,
  type = 'violin',
  idents = NULL,
  ncol = NULL,
  sort = FALSE,
  assay = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  adjust = 1,
  cols = NULL,
  pt.size = 0,
  group.by = NULL,
  split.by = NULL,
  log = FALSE,
  slot = 'data',
  stack = FALSE,
  combine = TRUE,
  fill.by = NULL,
  flip = FALSE,
  add.noise = TRUE,
  raster = NULL
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  if (isTRUE(x = stack)) {
    if (!is.null(x = ncol)) {
      warning(
        "'ncol' is ignored with 'stack' is TRUE",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    if (!is.null(x = y.max)) {
      warning(
        "'y.max' is ignored when 'stack' is TRUE",
        call. = FALSE,
        immediate. = TRUE
      )
    }
  } else {
    ncol <- ncol %||% ifelse(
      test = length(x = features) > 9,
      yes = 4,
      no = min(length(x = features), 3)
    )
  }
  data <- FetchData(object = object, vars = features, slot = slot)
  pt.size <- pt.size %||% AutoPointSize(data = object)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    } else if (length(x = cols) == 1 && cols == 'interaction') {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- levels(x = split)
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    message("WHHHHHHHHHHHHHHHHHHHHAST 1111111111111111")
    return(MultiExIPlot(
      type = type,
      data = data,
      idents = idents,
      split = split,
      sort = sort,
      same.y.lims = same.y.lims,
      adjust = adjust,
      cols = cols,
      pt.size = pt.size,
      log = log,
      fill.by = fill.by,
      add.noise = add.noise,
      flip = flip
    ))
  }

  # Return the updated data
  return(list(data = data, idents = idents))
}

uploadDotPlotDataToS3 <- function(data, user, project, analysis, genes) {
  gene_features <- jsonlite::fromJSON(genes)
  message("getting dot plot data for: ", gene_features)
  dotplot_data <- DotPlotData(object = data, features = gene_features)
  json_file <- tempfile(fileext = ".json")
  jsonlite::write_json(dotplot_data, json_file)
  s3_key <- paste0(user, "/", project, "/", analysis, "/dotplot.json")
  message(s3_key)
  message("Uploading dot plot data to s3")
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$dotplot_bucket
  )
  message("upload to S3 completed")
  file.remove(json_file)
}

uploadViolinPlotDataToS3 <- function(data, user, project, analysis, genes) {
  gene_features <- jsonlite::fromJSON(genes)
  message("getting violin plot data for: ", gene_features)
  vln_plots <- ExIPlot(object = data, features = gene_features)
  idents_vector <- unlist(vln_plots$idents)
  vln_plots$data$cluster_id <- idents_vector[rownames(vln_plots$data)]
  json_file <- tempfile(fileext = ".json")
  jsonlite::write_json(vln_plots$data, json_file)
  s3_key <- paste0(user, "/", project, "/", analysis, "/vlnplots.json")
  message(s3_key)
  message("Uploading violin plot data to s3")
  put_object(
    file = json_file,
    object = s3_key,
    bucket = user_environment$vlnplots_bucket
  )
  message("upload to S3 completed")
  file.remove(json_file)
}

gatherHeatmapDataAndUploadToS3 <- function(data, user, project, analysis, gene_names) {
  gene_features <- jsonlite::fromJSON(gene_names)
  heatmap_data <- DoTheHeatmap(data, features = gene_features)
  readProtoFiles("/app/heatmap.proto")
  heatmap_proto <- new(HeatmapData)
  heatmap_proto$featureNames <- colnames(heatmap_data[[1]])[-ncol(heatmap_data[[1]])]
  for (i in 1:nrow(heatmap_data[[1]])) {
    row <- new(DataRow)
    row$cellBarcode <- rownames(heatmap_data[[1]])[i]
    row$values <- as.numeric(heatmap_data[[1]][i, -ncol(heatmap_data[[1]])]) # Exclude the "ident" column
    row$cluster <- as.character(heatmap_data[[1]][i, ncol(heatmap_data[[1]])])
    heatmap_proto$dataRows[i] <- row
  }
  temp <- tempfile()
  s3key <- paste0(user, "/", project, "/", analysis, "/", "heatmap_data.bin")
  writeBin(RProtoBuf::serialize(heatmap_proto, NULL), temp)
  message("Uploading to S3: ", s3key)

  put_object(
    file = temp,
    object = s3key,
    bucket = user_environment$heatmap_bucket
  )
  message("Uploaded to S3")
  file.remove(temp)
}

find_nearest_vertex <- function(data_matrix, target_points, block_size=50000,
                                process_targets_in_blocks=FALSE){
  closest_vertex = c()
  if (process_targets_in_blocks == FALSE){
    num_blocks = ceiling(ncol(data_matrix) / block_size)
    if(num_blocks < 1) warning('bad loop: num_blocks < 1')
    for (i in 1:num_blocks){
      if (i < num_blocks){
        block = data_matrix[,((((i-1) * block_size)+1):(i*block_size))]
      }else{
        block = data_matrix[,((((i-1) * block_size)+1):(ncol(data_matrix)))]
      }
      distances_Z_to_Y <- proxy::dist(t(block), t(target_points))
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1,
                                        function(z) { which.min(z) } )
      closest_vertex = append(closest_vertex, closest_vertex_for_block)
    }
  }else{
    num_blocks = ceiling(ncol(target_points) / block_size)
    dist_to_closest_vertex = rep(Inf, length(ncol(data_matrix)))
    closest_vertex = rep(NA, length(ncol(data_matrix)))
    if(num_blocks < 1) warning('bad loop: num_blocks < 1')
    for (i in 1:num_blocks){
      if (i < num_blocks){
        block = target_points[,((((i-1) * block_size)+1):(i*block_size))]
      }else{
        block = target_points[,((((i-1) * block_size)+1):(ncol(target_points)))]
      }
      distances_Z_to_Y <- proxy::dist(t(data_matrix), t(block))
      closest_vertex_for_block <- apply(distances_Z_to_Y, 1,
                                        function(z) { which.min(z) } )
      if(nrow(distances_Z_to_Y) < 1) warning('bad loop: nrow(distances_Z_to_Y) < 1')
      new_block_distances <- distances_Z_to_Y[cbind(1:nrow(distances_Z_to_Y),
                                                    closest_vertex_for_block)]
      updated_nearest_idx <- which(new_block_distances < dist_to_closest_vertex)
      closest_vertex[updated_nearest_idx] <-
        closest_vertex_for_block[updated_nearest_idx] + (i-1) * block_size
      dist_to_closest_vertex[updated_nearest_idx] <-
        new_block_distances[updated_nearest_idx]
    }
  }
  stopifnot(length(closest_vertex) == ncol(data_matrix))
  return (closest_vertex)
}

extract_general_graph_ordering <- function(cds,
                                           root_pr_nodes,
                                           verbose=TRUE,
                                           reduction_method) {
  Z <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])
  Y <- cds@principal_graph_aux[[reduction_method]]$dp_mst
  pr_graph <- principal_graph(cds)[[reduction_method]]

  parents <- rep(NA, length(igraph::V(pr_graph)))
  states <- rep(NA, length(igraph::V(pr_graph)))

  if(any(is.na(igraph::E(pr_graph)$weight))) {
    igraph::E(pr_graph)$weight <- 1
  }

  # do pseudotime calculation on the cell-wise graph
  # 1. identify nearest cells to the selected principal node
  # 2. build a cell-wise graph for each Louvain group
  # 3. run the distance function to assign pseudotime for each cell
  closest_vertex <- find_nearest_vertex(Y[, root_pr_nodes, drop = FALSE], Z)
  closest_vertex_id <- colnames(cds)[closest_vertex]

  cell_wise_graph <-
    cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree
  cell_wise_distances <- igraph::distances(cell_wise_graph,
                                           v = closest_vertex_id)

  if (length(closest_vertex_id) > 1){
    node_names <- colnames(cell_wise_distances)
    pseudotimes <- apply(cell_wise_distances, 2, min)
  }else{
    node_names <- names(cell_wise_distances)
    pseudotimes <- cell_wise_distances
  }

  names(pseudotimes) <- node_names

  ordering_df <- data.frame(sample_name = igraph::V(cell_wise_graph)$name,
                            pseudo_time = as.vector(pseudotimes)
  )
  row.names(ordering_df) <- ordering_df$sample_name
  return(ordering_df)
}

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

branch_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points = branch_points[branch_points %in% root_nodes(cds, reduction_method) == FALSE]
  return(branch_points)
}

leaf_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  leaves <- which(igraph::degree(g) == 1)
  leaves = leaves[leaves %in% root_nodes(cds, reduction_method) == FALSE]
  return(leaves)
}

root_nodes <- function(cds, reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  root_pr_nodes <- which(names(igraph::V(g)) %in%
                    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes
  return(root_pr_nodes)
}

order_down_cells <- function(cds, root_pr_nodes,
                        reduction_method = "UMAP",
                        root_cells=NULL,
                        verbose = FALSE) {

  print(root_pr_nodes)
  cds@principal_graph_aux[[reduction_method]]$root_pr_nodes <- root_pr_nodes

  cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes, verbose,
                                                reduction_method)
  cds@principal_graph_aux[[reduction_method]]$pseudotime <-
    cc_ordering[row.names(colData(cds)), ]$pseudo_time
  names(cds@principal_graph_aux[[reduction_method]]$pseudotime) <-
    row.names(colData(cds))

  return(cds)
}

uploadPsuedotimePlotDataToS3 <- function(data, user, project, analysis, points) {
  data_points <- unlist(strsplit(points, " "))
  message(data_points)
  data.cds <- order_down_cells(data, data_points)
  message("Gathering psuedotime plot data and uploading it to s3")
  plot_cells_and_upload_to_s3(cds = data.cds, 
                                color_cells_by = "pseudotime", 
                                show_trajectory_graph = TRUE, 
                                user = user, 
                                project = project, 
                                analysis = analysis)
}

annotateSueratObjectClusters <- function(data, annotations) {
  rename_list <- jsonlite::fromJSON(annotations, simplifyVector = FALSE)
  args <- c(list(object = data), rename_list)

  message("Annotating the clusters: ", rename_list)
  renamed_data <- do.call(Seurat::RenameIdents, args)
  user_environment$cluster_object <- renamed_data  
  user_environment$scaled_object <- renamed_data           
  message("Annotated and updated seurat objects")
}

getAllMarkersDataAndUploadToS3 <- function(data, user, project, analysis) {
  message("Running find all markers..")
  markers <- FindAllMarkers(data)
  message("All Markers found")
  temp <- tempfile(fileext = ".csv")
  s3key <- paste0(user, "/", project, "/", analysis, "/", "allMarkersData.csv")
  write.csv(markers, temp)
  message("Uploading to S3: ", s3key)

  put_object(
    file = temp,
    object = s3key,
    bucket = user_environment$heatmap_bucket
  )
  message("Uploaded all markers to S3")
  file.remove(temp)
}

getFindMarkersDataAndUploadToS3 <- function(data, user, project, analysis, cluster1, cluster2) {
  message("Running find markers..", cluster1, cluster2)
  markers <- FindMarkers(data, ident.1 = cluster1, ident.2 = cluster2)
  message("Markers were found")
  temp <- tempfile(fileext = ".csv")
  s3key <- paste0(user, "/", project, "/", analysis, "/", "findMarkersData.csv")
  write.csv(markers, temp)
  message("Uploading to S3: ", s3key)

  put_object(
    file = temp,
    object = s3key,
    bucket = user_environment$heatmap_bucket
  )
  message("Uploaded found markers to S3")
  file.remove(temp)
}

#* @post /init_endpoint
#* @parser json
function(req, res) {
  req_data <- req$body
  datasets <- req_data$datasets
  user <- req_data$user
  project <- req_data$project
  analysis <- req_data$analysis
  message("Initializing...")
  message(datasets)
  data <- init_dataset_seurat_object(
    user = user,
    project = project,
    datasets = datasets,
    analysis = analysis)
  message("Successfully initialized dataset seurat object")
  user_environment$init_object <- data
  res$status <- 200
  return(list(success = TRUE, message = "Initialized project successfully"))
}

#* @post /variableFeatures
#* @parser json
function(req, res) {
  user <- req$body$user
  project <- req$body$project
  analysis <- req$body$analysis
  nFeatures <- req$body$nFeatures
  qc_data <- user_environment$init_object
  data <- variable_feature_search(qc_data, user, project, analysis, as.numeric(nFeatures))
  message("Successfully found variable features and performed PCA.")
  user_environment$scaled_object <- data
  res$status <- 200
  return(list(success = TRUE, message = "Found variable features and performed PCA successfully"))
}

#* @post /cluster
#* @parser json
function(req, res) {
  cluster_data <- req$body
  return_data <- perform_clustering(
    user_environment$scaled_object,
    cluster_data$user,
    cluster_data$project,
    cluster_data$analysis,
    cluster_data$cluster_neighbors,
    cluster_data$cluster_clusters,
    cluster_data$cluster_dimensions,
    cluster_data$cluster_reduction
  )
  message("Successfully clustered the data.")
  user_environment$cluster_object <- return_data$seurat
  user_environment$monocle_object <- return_data$monocle
  res$status <- 200
  return(list(success = TRUE, message = "Performed clustering successfully"))
}

#* @post /psuedotime
#* @parser json
function(req, res) {
  psuedo_data <- req$body
  uploadPsuedotimePlotDataToS3(
    user_environment$monocle_object,
    psuedo_data$user,
    psuedo_data$project,
    psuedo_data$analysis,
    psuedo_data$points
  )
  message("Successfully uploaded the psuedotime plot data.")
  res$status <- 200
  return(list(success = TRUE, message = "Performed psuedotime successfully"))
}

#* @post /genefeature
#* @parser json
function(req, res) {
  gene_name <- req$body$gene_name
  user <- req$body$user
  project <- req$body$project
  analysis <- req$body$analysis
  message(gene_name)
  get_gene_feature_plot_data(
    user_environment$cluster_object,
    gene_name,
    user,
    project,
    analysis
  )
  message("Successfully processed feature plot.")
  res$status <- 200
  return(list(success = TRUE, message = "Processed the gene feature plot successfully"))
}

#* @post /heatmap_analysis
#* @parser json
function(req, res) {
  gene_names <- req$body$geneNames
  user <- req$body$user
  project <- req$body$project
  analysis <- req$body$analysis
  message(gene_names)
  gatherHeatmapDataAndUploadToS3(user_environment$cluster_object, user, project, analysis, gene_names)
  message("Successfully processed the heatmap.")
  res$status <- 200
  return(list(success = TRUE, message = "Processed the gene feature plot successfully"))
}

#* @post /dotplot
#* @parser json
function(req, res) {
  dotplot_data <- req$body
  uploadDotPlotDataToS3(
    user_environment$cluster_object,
    dotplot_data$user,
    dotplot_data$project,
    dotplot_data$analysis,
    dotplot_data$genes
  )
  message("Successfully uploaded the dot plot data.")
  res$status <- 200
  return(list(success = TRUE, message = "Gathered dot plots successfully"))
}

#* @post /violinplots
#* @parser json
function(req, res) {
  violin_data <- req$body
  uploadViolinPlotDataToS3(
    user_environment$cluster_object,
    violin_data$user,
    violin_data$project,
    violin_data$analysis,
    violin_data$genes
  )
  message("Successfully uploaded the violin plot data.")
  res$status <- 200
  return(list(success = TRUE, message = "Gathered violin plots successfully"))
}

#* @post /annotations
#* @parser json
function(req, res) {
  annotations_data <- req$body
  message("Annotating the Seurat object clusters: ", annotations_data$annotations)
  annotateSueratObjectClusters(
    user_environment$cluster_object,
    annotations_data$annotations
  )
  message("Successfully annotated the Suerat clusters.")
  res$status <- 200
  return(list(success = TRUE, message = "Performed annotations successfully"))
}
#* @post /allMarkers
#* @parser json
function(req, res) {
  allMarkers_data <- req$body
  getAllMarkersDataAndUploadToS3(
    user_environment$cluster_object,
    allMarkers_data$user,
    allMarkers_data$project,
    allMarkers_data$analysis
  )
  message("Successfully uploaded the all markers data.")
  res$status <- 200
  return(list(success = TRUE, message = "Gathered all markers data successfully"))
}

#* @post /findMarkers
#* @parser json
function(req, res) {
  findMarkers_data <- req$body
  getFindMarkersDataAndUploadToS3(
    user_environment$cluster_object,
    findMarkers_data$user,
    findMarkers_data$project,
    findMarkers_data$analysis,
    findMarkers_data$cluster1,
    findMarkers_data$cluster2
  )
  message("Successfully uploaded the find markers data.")
  res$status <- 200
  return(list(success = TRUE, message = "Gathered find markers data successfully"))
}

#* @post /shutdown
#* @parser json
function(req, res) {
  message("Shutting down the Cellborg Analysis Troop server...")
  quit(save = "no", status = 0, runLast = TRUE)
}