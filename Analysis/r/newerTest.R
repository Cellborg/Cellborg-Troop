library(Seurat)
library(Matrix)
library(dplyr)
library(jsonlite)
library(CellChat)

netVisual_circle_data <-function(net, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                            weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                            edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                            edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                            arrow.width=1,arrow.size = 0.2){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0

  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0


  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }

  return(net)
}
netVisual_aggregate_data <- function(object, signaling, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                layout = c("circle","hierarchy","chord","spatial"),
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                alpha.image = 0.15, point.size = 1.5,
                                group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20,
                                ...) {
  layout <- match.arg(layout)
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }


  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  # prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
  }  else if (layout == "spatial") {
    prob.sum <- apply(prob, c(1,2), sum)
    if (vertex.weight == "incoming"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$indeg
    }
    if (vertex.weight == "outgoing"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$outdeg
    }
    coordinates <- object@images$coordinates
    labels <- object@idents
    gg <- netVisual_spatial(prob.sum, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)

  } else if (layout == "chord") {
    prob.sum <- apply(prob, c(1,2), sum)
    
    plot_data <- list(
      prob_sum = prob.sum,
      color_use = color.use,
      sources_use = sources.use,
      targets_use = targets.use,
      remove_isolate = remove.isolate,
      group = group,
      cell_order = cell.order,
      lab_cex = vertex.label.cex,
      small_gap = small.gap,
      big_gap = big.gap,
      scale = scale,
      reduce = reduce,
      signaling_name = signaling.name,
      show_legend = show.legend,
      legend_pos_x = legend.pos.x,
      legend_pos_y = legend.pos.y
    )
  }
  return(plot_data)
}


#Pre-Processing (Seurat)
cat("Loading dataset and creating seurat object")
pbmc <- Read10X(data.dir = "C:\\Users\\Pranik2\\Downloads\\filtered_gene_bc_matrices\\hg19") # nolint: line_length_linter.
pbmc <- CreateSeuratObject(counts = pbmc, min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

cat("Performing pre-processing of dataset")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RenameIdents(object = pbmc, c("0" = "Zero", "1" = "One", "2" = "Two", "3" = "Three", "4" = "Four", "5" = "Five", "6" = "Six", "7" = "Seven", "8" = "Eight" ))
#Compute CellChat
cellchat <- createCellChat(object = pbmc, group.by = "ident")
message("Cellchat obj created")
CellChatDB <- CellChatDB.human #SELECT SPECIES USER INPUT NEEDED, SEE ALL SPECIES
CellChatDB$interaction["H2-BI_KIR3DL1","interaction_name"] = "H2-BL_KIR3DL1"
CellChatDB$interaction["H2-BI_KIR3DL1","ligand"] = "H2-Bl"
CellChatDB$interaction["H2-BI_KIR3DL1","interaction_name_2"] = "H2-bl - Kir3dl1"
rownames(CellChatDB$interaction)[rownames(CellChatDB$interaction) == "H2-BI_KIR3DL1"] = "H2-BL_KIR3DL1"

CellChatDB$interaction["H2-EA-PS_CD4","interaction_name"] = "H2-EA_CD4"
CellChatDB$interaction["H2-EA-PS_CD4","ligand"] = "H2-Ea"
CellChatDB$interaction["H2-EA-PS_CD4","interaction_name_2"] = "H2-ea - Cd4"
rownames(CellChatDB$interaction)[rownames(CellChatDB$interaction) == "H2-EA-PS_CD4"] = "H2-EA_CD4"
dplyr::glimpse(CellChatDB$interaction)
message("WE HAVE GLUMPEDS")
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) # necessary even if using the whole database
message("Cellchat subsetted")
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
message("IDENTITY HAS BEEN EXPRESSEDNIGAG")
cellchat <- projectData(cellchat, PPI.human) #NEED TO CHANGE BASED ON USER INPUT
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
message("Aggregating")
groupSize <- as.numeric(table(cellchat@idents))
#Select Pathways
pathways.show <- c("ICAM") #USER INPUT REQUIRED

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
ppp <- netVisual_aggregate_data(cellchat, signaling = pathways.show, layout = "chord")
