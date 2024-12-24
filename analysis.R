
library(data.table)
library(limma)
library(Seurat)
library(dplyr)
library(PCAtools)
library(ggplot2)
library(ggridges)
library(plyr)
library(GSVA)
library(cowplot)
library(harmony)
library(ggforce)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(Scissor)
library(pheatmap)
library(DoubletFinder)
library(survival)
library(survminer)

###########################################functions###########################################
filter_10x<<-function(scExp,outdir,gene_num_min=500,gene_num_max=8000,MT_frac_min=10){
  system(paste("mkdir -p ", outdir, sep =""))
  exp <- scExp
  exp[["percent.mt"]] <- PercentageFeatureSet(exp, pattern = "^MT-")
  
  pdf(paste(outdir, "/Features.pdf", sep = ""), height = 7, width = 18)
  plot1<-VlnPlot(exp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = F)
  plot2<-VlnPlot(exp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = T)
  print(plot1)
  print(plot2)
  dev.off()
  
  
  
  ## Filter
  exp <- subset(exp, subset = nFeature_RNA > as.numeric(gene_num_min) & nFeature_RNA < as.numeric(gene_num_max) & percent.mt < as.numeric(MT_frac_min))
  # Normalization
  if (ncol(exp)>100){
    exp <- NormalizeData(exp, verbose = FALSE)
    exp <- FindVariableFeatures(exp, selection.method = "vst", nfeatures = 2000)
    pdf(paste(outdir, "/Variable_genes.pdf", sep = ""), height = 7, width = 18)
    top10 <- head(VariableFeatures(exp), 10)
    plot1_exp <- VariableFeaturePlot(exp)
    plot2_exp <- LabelPoints(plot = plot1_exp, points = top10, repel = TRUE)
    #CombinePlots(plots = list(plot1_exp, plot2_exp))
    print(plot2_exp)
    dev.off()
    ## Pre-process
    exp <- ScaleData(exp)
    exp <- RunPCA(exp)
    exp <- RunUMAP(exp, dims = 1:10)
    exp <- RunTSNE(exp, reduction = "pca", dims = 1:20, check_duplicates = FALSE)
    ## pK Identification
    sweep.res.list_exp <- paramSweep_v3(exp, PCs = 1:10, sct = FALSE)
    sweep.stats_exp <- summarizeSweep(sweep.res.list_exp, GT = FALSE)
    bcmvn_exp <- find.pK(sweep.stats_exp)
    mpK<-as.numeric(as.vector(bcmvn_exp$pK[which.max(bcmvn_exp$BCmetric)]))
    
    ## Clustering
    exp <- FindNeighbors(exp, reduction = "pca", dims = 1:20)
    exp <- FindClusters(exp, resolution = 0.2)
    pdf(paste(outdir, "/Cluster_merge_umap.pdf", sep = ""), height = 7, width = 18)
    p1 <- DimPlot(exp, reduction = "umap")
    p2 <- DimPlot(exp, reduction = "umap", label = TRUE)
    print(plot_grid(p1, p2))
    dev.off()
    pdf(paste(outdir, "/Cluster_merge_tsne.pdf", sep = ""), height = 7, width = 18)
    p1 <- DimPlot(exp, reduction = "tsne")
    p2 <- DimPlot(exp, reduction = "tsne", label = TRUE)
    print(plot_grid(p1, p2))
    dev.off()
    
    ## Homotypic Doublet Proportion Estimate
    annotations <- exp@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
    nExp_poi <- round(0.075*nrow(exp@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    exp <- doubletFinder_v3(exp, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    exp <- doubletFinder_v3(exp, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = paste("pANN","0.25",mpK,nExp_poi,sep="_"), sct = FALSE)
    
    name <- paste("DF.classifications", 0.25, mpK, nExp_poi.adj, sep = "_")
    pdf(paste(outdir, "/Doublet_umap.pdf", sep = ""), height = 7, width = 9)
    print(DimPlot(exp, reduction = "umap", group.by = name))
    dev.off()
    colnames(exp@meta.data)[ncol(exp@meta.data)] = "DF.classifications"
    exp <- subset(exp, subset = DF.classifications == "Singlet")
    ## Save
    saveRDS(exp, file = paste(outdir, "/Exp.pre_precessed.rds", sep = ""))
  }else{
    print("Not enough cells for QC")
  }
}
survplot <- function(dat = mydata, type = "OS", fit = fit, pval = pval){
  p <- ggsurvplot(fit,
                  linetype = 1,
                  #censor.shape=45,
                  data = dat,
                  size = 1, # change line size
                  palette = c("#6495EDFF", "#FF4500FF"),# custom color palettes
                  #conf.int = TRUE, # Add confidence interval
                  pval = paste('p = ', round(pval, 3)), # Add p-value
                  risk.table = T, # Add risk table
                  #tables.theme = theme_survminer(font.main = 10),
                  #risk.table.col = "strata",# Risk table color by groups
                  legend = "right",
                  #legend.labs = c("G1 (n = 7)", "G2 (n = 26)", "G3 (n = 24)"), # Change legend labels
                  risk.table.height = 0.25, # Useful to change when you have multiple groups
                  ggtheme = theme_classic2(), # Change ggplot2 theme
                  xlab = "Time (years)",
                  ylab = paste0("Probability of ", type))
  return(p)
}
surv_bestcut <- function(mtr, gene, cli_info, num = 20){
  tmp_mtr <- mtr[which(rownames(mtr) == gene), ]
  if (length(tmp_mtr[is.na(tmp_mtr)]) == 0) {
    tmp_mtr <- tmp_mtr
  }else{
    tmp_mtr <- tmp_mtr[-which(is.na(tmp_mtr))]
  }
  common_samples <- intersect(names(tmp_mtr), rownames(cli_info))
  cluster_surv <- cli_info[common_samples, ]
  tmp_mtr <- as.data.frame(t(tmp_mtr))[common_samples, ]
  sevalue <- as.numeric(tmp_mtr)
  values <- c()
  hr_rfs <- c()
  hr_rfs_high <- c()
  hr_rfs_low <- c()
  pva_rfs <- c()
  n_high <- c()
  n_low <- c()
  for (i in c(round(length(sevalue)/10):round(length(sevalue)/10*9))) {
    cluster_surv$Type = ifelse(sevalue > sort(sevalue)[i], "1.High", "0.Low")
    values <- c(values, sort(sevalue)[i])
    ##PFS
    tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
    hr_rfs <- c(hr_rfs, tmp$conf.int[[1]])
    pva_rfs <- c(pva_rfs, tmp$logtest[[3]])
    n_high <- c(n_high, (length(sevalue)-i))
    n_low <- c(n_low, i)
    hr_rfs_high <- c(hr_rfs_high, tmp$conf.int[[4]])
    hr_rfs_low <- c(hr_rfs_low, tmp$conf.int[[3]])
    if (i == num) {
      tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
      fit <- survfit(Surv(OS.time, OS) ~ Type, data = cluster_surv)
      pfs <- survplot(cluster_surv, type = "OS-LIHC", fit = fit, pval = tmp$logtest[3])
      pdf(paste("../3.Analysis/3.X.singlegene.", gene, ".surv.pdf", sep = ""), width = 6.5, height = 5)
      print(pfs)
      dev.off()
    }
  }
  res <- data.frame(ID = gene,
                    cutoff = values,
                    HR_OS = hr_rfs,
                    HR_OS_H = hr_rfs_high,
                    HR_OS_L = hr_rfs_low,
                    Pvalue_OS = pva_rfs,
                    n_high = n_high,
                    n_low = n_low)
  
}


as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
CountCellFreq <- function(loc = "Lauren", metadata, metacluster = "CompClass"){
  res <- data.frame()
  for (clusters in unique(metadata[, loc])){
    tmpdf <- metadata[metadata[, loc] == clusters, ]
    tmpres <- as.data.frame(prop.table(table(tmpdf[, c("Patient", metacluster)])))
    tmpres$Freq <- 100*tmpres$Freq
    tmpres$Class <- clusters
    res <- rbind(res, tmpres)
  }
  return(res)
}


do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  return(count.dist.melt.ext.tb)
}


###########################################1.MergeDatasets###########################################
####################1.1.generate RDS####################
#######GSE125449#######
allfiles <- list.files("../../LiverTME/0.data/GSE125449/", "barcodes.tsv.gz")
allfiles <- paste(strsplit2(allfiles, split = "_")[, 1], strsplit2(allfiles, split = "_")[, 2], sep = "_")

for (i in 1:length(allfiles)) {
  files = allfiles[i]
  dir.create(paste("../../LiverTME/0.data/GSE125449/", files, sep = ""))
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE125449/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE125449/", pattern = "barcodes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE125449/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE125449/", files, "/barcodes.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE125449/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE125449/", pattern = "genes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE125449/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE125449/", files, "/features.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE125449/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE125449/", pattern = "matrix.mtx.gz")
  file.copy(paste("../../LiverTME/0.data/GSE125449/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE125449/", files, "/matrix.mtx.gz", sep = ""), overwrite = T)
  GSE125449.mtr.tmp <- Read10X(paste("../../LiverTME/0.data/GSE125449/", files, sep = ""))
  GSE125449.saminfo.tmp <- data.frame(ID = colnames(GSE125449.mtr.tmp),
                                      GSEID = files)
  rownames(GSE125449.saminfo.tmp) <- GSE125449.saminfo.tmp$ID
  GSE125449.obj.tmp <- CreateSeuratObject(GSE125449.mtr.tmp, meta.data = GSE125449.saminfo.tmp)
  if (i == 1) {
    GSE125449.obj <- GSE125449.obj.tmp
  }else{
    GSE125449.obj <- merge(GSE125449.obj, y = GSE125449.obj.tmp)
  }
}

GSE125449.obj$ID <- rownames(GSE125449.obj@meta.data)
saminfo1 <- fread("../../LiverTME/0.data/GSE125449/GSE125449_Set1_samples.txt.gz", header = T, stringsAsFactors = F, data.table = F)
saminfo2 <- fread("../../LiverTME/0.data/GSE125449/GSE125449_Set2_samples.txt.gz", header = T, stringsAsFactors = F, data.table = F)
saminfo <- rbind(saminfo1, saminfo2)
saminfo$Patient <- saminfo$Sample
saminfo$Tissue <- "Tumor"
saminfo$Metas <- NA
saminfo$Lauren <- NA
rownames(saminfo) <- saminfo$`Cell Barcode`
GSE125449.obj <- AddMetaData(GSE125449.obj, saminfo[, c(4, 5)])
GSE125449.obj$GSEID <- "GSE125449"
GSE125449.obj$Dataset <- "GSE125449"
saveRDS(GSE125449.obj, file = "../0.data/GSE125449.obj.rds")

#######GSE140228#######
dir.create("../../LiverTME/0.data/GSE140228/GSE140228")
file.copy("../../LiverTME/0.data/GSE140228/GSE140228_UMI_counts_Droplet_barcodes.tsv.gz", "../../LiverTME/0.data/GSE140228/GSE140228/barcodes.tsv.gz", overwrite = T)
file.copy("../../LiverTME/0.data/GSE140228/GSE140228_UMI_counts_Droplet_genes.tsv.extra.gz", "../../LiverTME/0.data/GSE140228/GSE140228/features.tsv.gz", overwrite = T)
file.copy("../../LiverTME/0.data/GSE140228/GSE140228_UMI_counts_Droplet.mtx.gz", "../../LiverTME/0.data/GSE140228/GSE140228/matrix.mtx.gz", overwrite = T)

GSE140228.mtr.tmp <- Read10X("../../LiverTME/0.data/GSE140228/GSE140228/")
GSE140228.saminfo.tmp <- data.frame(ID = colnames(GSE140228.mtr.tmp),
                                    GSEID = "GSE140228")
rownames(GSE140228.saminfo.tmp) <- GSE140228.saminfo.tmp$ID
GSE140228.obj.tmp <- CreateSeuratObject(GSE140228.mtr.tmp, meta.data = GSE140228.saminfo.tmp)
GSE140228.obj <- GSE140228.obj.tmp
saminfo <- fread("../../LiverTME/0.data/GSE140228/GSE140228_UMI_counts_Droplet_cellinfo.tsv.gz", header = T, stringsAsFactors = F, data.table = F)
saminfo$Patient <- saminfo$Donor
saminfo$Tissue <- saminfo$Tissue
saminfo$Metas <- NA
saminfo$Lauren <- NA
rownames(saminfo) <- saminfo$Barcode
GSE140228.obj <- AddMetaData(GSE140228.obj, saminfo[, c(10, 3)])
GSE140228.obj$GSEID <- "GSE140228"
GSE140228.obj$Dataset <- "GSE140228"
saveRDS(GSE140228.obj, file = "../0.data/GSE140228.obj.rds")

#######GSE112271#######
allfiles <- list.files("../../LiverTME/0.data/GSE112271/", "barcodes.tsv.gz")
allfiles <- strsplit2(allfiles, split = "_")[, 1]

for (i in 1:length(allfiles)) {
  files = allfiles[i]
  dir.create(paste("../../LiverTME/0.data/GSE112271/", files, sep = ""))
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE112271/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE112271/", pattern = "barcodes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE112271/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE112271/", files, "/barcodes.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE112271/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE112271/", pattern = "genes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE112271/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE112271/", files, "/features.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE112271/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE112271/", pattern = "matrix.mtx.gz")
  file.copy(paste("../../LiverTME/0.data/GSE112271/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE112271/", files, "/matrix.mtx.gz", sep = ""), overwrite = T)
  GSE112271.mtr.tmp <- Read10X(paste("../../LiverTME/0.data/GSE112271/", files, sep = ""))
  GSE112271.saminfo.tmp <- data.frame(ID = colnames(GSE112271.mtr.tmp),
                                      GSEID = files)
  rownames(GSE112271.saminfo.tmp) <- GSE112271.saminfo.tmp$ID
  GSE112271.obj.tmp <- CreateSeuratObject(GSE112271.mtr.tmp, meta.data = GSE112271.saminfo.tmp)
  if (i == 1) {
    GSE112271.obj <- GSE112271.obj.tmp
  }else{
    GSE112271.obj <- merge(GSE112271.obj, y = GSE112271.obj.tmp)
  }
}

GSE112271.obj$ID <- rownames(GSE112271.obj@meta.data)
saminfo <- fread("../../LiverTME/0.data/GSE112271/GSE112271.txt", header = F, stringsAsFactors = F, data.table = F)
saminfo$Patient <- paste("GSE112271", "_", gsub(" ", "_", strsplit2(saminfo$V2, "\\.")[, 1]), sep = "")
saminfo$Tissue <- "Tumor"
saminfo$Metas <- NA
saminfo$Lauren <- NA
rownames(saminfo) <- saminfo$V1
saminfo <- merge(saminfo, GSE112271.obj@meta.data, by.x = "V1", by.y = "GSEID", all.y = T)
rownames(saminfo) <- saminfo$ID
GSE112271.obj <- AddMetaData(GSE112271.obj, saminfo[, c(3, 4)])
GSE112271.obj$GSEID <- "GSE112271"
GSE112271.obj$Dataset <- "GSE112271"
saveRDS(GSE112271.obj, file = "../0.data/GSE112271.obj.rds")

#######GSE178318#######
allfiles <- list.files("../../LiverTME/0.data/GSE178318/", "barcodes.tsv.gz")
allfiles <- strsplit2(allfiles, split = "_")[, 1]

for (i in 1:length(allfiles)) {
  files = allfiles[i]
  dir.create(paste("../../LiverTME/0.data/GSE178318/", files, sep = ""))
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE178318/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE178318/", pattern = "barcodes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE178318/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE178318/", files, "/barcodes.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE178318/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE178318/", pattern = "genes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE178318/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE178318/", files, "/features.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE178318/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE178318/", pattern = "matrix.mtx.gz")
  file.copy(paste("../../LiverTME/0.data/GSE178318/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE178318/", files, "/matrix.mtx.gz", sep = ""), overwrite = T)
  GSE178318.mtr.tmp <- Read10X(paste("../../LiverTME/0.data/GSE178318/", files, sep = ""))
  GSE178318.saminfo.tmp <- data.frame(ID = colnames(GSE178318.mtr.tmp),
                                      GSEID = files)
  rownames(GSE178318.saminfo.tmp) <- GSE178318.saminfo.tmp$ID
  GSE178318.obj.tmp <- CreateSeuratObject(GSE178318.mtr.tmp, meta.data = GSE178318.saminfo.tmp)
  if (i == 1) {
    GSE178318.obj <- GSE178318.obj.tmp
  }else{
    GSE178318.obj <- merge(GSE178318.obj, y = GSE178318.obj.tmp)
  }
}

GSE178318.obj$ID <- rownames(GSE178318.obj@meta.data)

saminfo <- data.frame(ID = colnames(GSE178318.obj),
                      Patient = strsplit2(colnames(GSE178318.obj), "_")[, 2],
                      Tissue = strsplit2(colnames(GSE178318.obj), "_")[, 3])

saminfo$Metas <- NA
saminfo$Lauren <- NA
rownames(saminfo) <- saminfo$ID
GSE178318.obj <- AddMetaData(GSE178318.obj, saminfo[, c(2, 3)])
GSE178318.obj$GSEID <- "GSE178318"
GSE178318.obj$Dataset <- "GSE178318"
saveRDS(GSE178318.obj, file = "../0.data/GSE178318.obj.rds")

#######GSE151530#######
allfiles <- list.files("../../LiverTME/0.data/GSE151530/", "barcodes.tsv.gz")
allfiles <- strsplit2(allfiles, split = "_")[, 1]

for (i in 1:length(allfiles)) {
  files = allfiles[i]
  dir.create(paste("../../LiverTME/0.data/GSE151530/", files, sep = ""))
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE151530/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE151530/", pattern = "barcodes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE151530/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE151530/", files, "/barcodes.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE151530/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE151530/", pattern = "genes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE151530/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE151530/", files, "/features.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE151530/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE151530/", pattern = "matrix.mtx.gz")
  file.copy(paste("../../LiverTME/0.data/GSE151530/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE151530/", files, "/matrix.mtx.gz", sep = ""), overwrite = T)
  GSE151530.mtr.tmp <- Read10X(paste("../../LiverTME/0.data/GSE151530/", files, sep = ""))
  GSE151530.saminfo.tmp <- data.frame(ID = colnames(GSE151530.mtr.tmp),
                                      GSEID = files)
  rownames(GSE151530.saminfo.tmp) <- GSE151530.saminfo.tmp$ID
  GSE151530.obj.tmp <- CreateSeuratObject(GSE151530.mtr.tmp, meta.data = GSE151530.saminfo.tmp)
  if (i == 1) {
    GSE151530.obj <- GSE151530.obj.tmp
  }else{
    GSE151530.obj <- merge(GSE151530.obj, y = GSE151530.obj.tmp)
  }
}

GSE151530.obj$ID <- rownames(GSE151530.obj@meta.data)
saminfo <- fread("../../LiverTME/0.data/GSE151530/GSE151530_Info.txt.gz", header = T, stringsAsFactors = F, data.table = F)
saminfo$Patient <- saminfo$Sample
saminfo$Tissue <- "Tumor"
saminfo$Metas <- NA
saminfo$Lauren <- NA
rownames(saminfo) <- saminfo$Cell
GSE151530.obj <- AddMetaData(GSE151530.obj, saminfo[, c(5, 6)])
GSE151530.obj$GSEID <- "GSE151530"
GSE151530.obj$Dataset <- "GSE151530"
saveRDS(GSE151530.obj, file = "../0.data/GSE151530.obj.rds")

#######GSE166635#######
allfiles <- list.files("../../LiverTME/0.data/GSE166635/", "barcodes.tsv.gz")
allfiles <- strsplit2(allfiles, split = "_")[, 1]

for (i in 1:length(allfiles)) {
  files = allfiles[i]
  dir.create(paste("../../LiverTME/0.data/GSE166635/", files, sep = ""))
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE166635/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE166635/", pattern = "barcodes.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE166635/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE166635/", files, "/barcodes.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE166635/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE166635/", pattern = "features.tsv.gz")
  file.copy(paste("../../LiverTME/0.data/GSE166635/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE166635/", files, "/features.tsv.gz", sep = ""), overwrite = T)
  sepfiles1 <- list.files("../../LiverTME/0.data/GSE166635/", pattern = files)
  sepfiles2 <- list.files("../../LiverTME/0.data/GSE166635/", pattern = "matrix.mtx.gz")
  file.copy(paste("../../LiverTME/0.data/GSE166635/", intersect(sepfiles1, sepfiles2), sep = ""), 
            paste("../../LiverTME/0.data/GSE166635/", files, "/matrix.mtx.gz", sep = ""), overwrite = T)
  GSE166635.mtr.tmp <- Read10X(paste("../../LiverTME/0.data/GSE166635/", files, sep = ""))
  GSE166635.saminfo.tmp <- data.frame(ID = colnames(GSE166635.mtr.tmp),
                                      GSEID = files)
  rownames(GSE166635.saminfo.tmp) <- GSE166635.saminfo.tmp$ID
  GSE166635.obj.tmp <- CreateSeuratObject(GSE166635.mtr.tmp, meta.data = GSE166635.saminfo.tmp)
  if (i == 1) {
    GSE166635.obj <- GSE166635.obj.tmp
  }else{
    GSE166635.obj <- merge(GSE166635.obj, y = GSE166635.obj.tmp)
  }
}

GSE166635.obj$ID <- rownames(GSE166635.obj@meta.data)
GSE166635.obj$Patient <- GSE166635.obj$GSEID
GSE166635.obj$Tissue <- "Tumor"
GSE166635.obj$GSEID <- "GSE166635"
GSE166635.obj$Dataset <- "GSE166635"
saveRDS(GSE166635.obj, file = "../0.data/GSE166635.obj.rds")



#######GSE164522#######
saminfo <- fread("../../LiverTME/0.data/GSE164522/GSE164522_CRLM_metadata.csv.gz", header = T, data.table = F)
saminfo$V1 <- gsub("-", ".", saminfo$V1)
GSE164522_CRLM_PT.mtr <- fread("../../LiverTME/0.data/GSE164522/GSE164522_CRLM_PT_expression.csv.gz", header = T, data.table = F)
GSE164522_CRLM_MT.mtr <- fread("../../LiverTME/0.data/GSE164522/GSE164522_CRLM_MT_expression.csv.gz", header = T, data.table = F)
GSE164522_CRLM_MN.mtr <- fread("../../LiverTME/0.data/GSE164522/GSE164522_CRLM_MN_expression.csv.gz", header = T, data.table = F)
GSE164522_CRLM_PN.mtr <- fread("../../LiverTME/0.data/GSE164522/GSE164522_CRLM_PN_expression.csv.gz", header = T, data.table = F)
GSE164522_CRLM_LN.mtr <- fread("../../LiverTME/0.data/GSE164522/GSE164522_CRLM_LN_expression.csv.gz", header = T, data.table = F)
GSE164522_CRLM_PBMC.mtr <- fread("../../LiverTME/0.data/GSE164522/GSE164522_CRLM_PBMC_expression.csv.gz", header = T, data.table = F)
rownames(GSE164522_CRLM_PT.mtr) <- GSE164522_CRLM_PT.mtr$V1
rownames(GSE164522_CRLM_MT.mtr) <- GSE164522_CRLM_MT.mtr$V1
rownames(GSE164522_CRLM_MN.mtr) <- GSE164522_CRLM_MN.mtr$V1
rownames(GSE164522_CRLM_PN.mtr) <- GSE164522_CRLM_PN.mtr$V1
rownames(GSE164522_CRLM_LN.mtr) <- GSE164522_CRLM_LN.mtr$V1
rownames(GSE164522_CRLM_PBMC.mtr) <- GSE164522_CRLM_PBMC.mtr$V1
GSE164522_CRLM_PT.mtr <- GSE164522_CRLM_PT.mtr[, -1]
GSE164522_CRLM_MT.mtr <- GSE164522_CRLM_MT.mtr[, -1]
GSE164522_CRLM_MN.mtr <- GSE164522_CRLM_MN.mtr[, -1]
GSE164522_CRLM_PN.mtr <- GSE164522_CRLM_PN.mtr[, -1]
GSE164522_CRLM_LN.mtr <- GSE164522_CRLM_LN.mtr[, -1]
GSE164522_CRLM_PBMC.mtr <- GSE164522_CRLM_PBMC.mtr[, -1]

GSE164522_CRLM.mtr <- cbind(GSE164522_CRLM_PT.mtr, GSE164522_CRLM_MT.mtr, GSE164522_CRLM_MN.mtr,
                            GSE164522_CRLM_PN.mtr, GSE164522_CRLM_LN.mtr, GSE164522_CRLM_PBMC.mtr)
GSE164522_CRLM.saminfo <- data.frame(ID = c(names(GSE164522_CRLM_PT.mtr), names(GSE164522_CRLM_MT.mtr), names(GSE164522_CRLM_MN.mtr),
                                            names(GSE164522_CRLM_PN.mtr), names(GSE164522_CRLM_LN.mtr), names(GSE164522_CRLM_PBMC.mtr)),
                                     Tissue = c(rep("PT", ncol(GSE164522_CRLM_PT.mtr)), rep("MT", ncol(GSE164522_CRLM_MT.mtr)),
                                                rep("MN", ncol(GSE164522_CRLM_MN.mtr)), rep("PN", ncol(GSE164522_CRLM_PN.mtr)),
                                                rep("LN", ncol(GSE164522_CRLM_LN.mtr)), rep("PBMC", ncol(GSE164522_CRLM_PBMC.mtr))))
GSE164522_CRLM.saminfo <- merge(GSE164522_CRLM.saminfo, saminfo[, c(1, 7)], by.x = "ID", by.y = "V1")
names(GSE164522_CRLM.saminfo) <- c("ID", "Tissue", "Patient")
rownames(GSE164522_CRLM.saminfo) <- GSE164522_CRLM.saminfo$ID
GSE164522_CRLM.obj <- CreateSeuratObject(GSE164522_CRLM.mtr, meta.data = GSE164522_CRLM.saminfo)
GSE164522_CRLM.obj@meta.data$GSEID <- "GSE164522"
GSE164522_CRLM.obj$Dataset <- "GSE164522"
saveRDS(GSE164522_CRLM.obj, file = "../0.data/GSE164522.obj.rds")


#######scCRLM#######
load("../../LiverTME/0.data/scCRLM/metadata.rda")
load("../../LiverTME/0.data/scCRLM/exprmatrix.rda")

scCRLM.mtr <- exprmatrix
scCRLM.saminfo <- data.frame(ID = rownames(metadata),
                             Patient = metadata$patient,
                             Tissue = metadata$tissue)
rownames(scCRLM.saminfo) <- scCRLM.saminfo$ID
scCRLM.obj <- CreateSeuratObject(scCRLM.mtr, meta.data = scCRLM.saminfo)
scCRLM.obj@meta.data$GSEID <- "scCRLM"
scCRLM.obj$Dataset <- "scCRLM"
saveRDS(scCRLM.obj, file = "../0.data/scCRLM.obj.rds")

#######GSE132465#######
GSE132465.obj <- readRDS("../../../../../other_works/wangyun/20230403CRCLiver/0.data/GSE132465.rds")
saveRDS(GSE132465.obj, file = "../0.data/GSE132465.obj.rds")


####################1.2.patientRDSfilter####################
GSE125449 <- readRDS("../0.data/GSE125449.obj.rds")
GSE140228 <- readRDS("../0.data/GSE140228.obj.rds")
GSE112271 <- readRDS("../0.data/GSE112271.obj.rds")
GSE178318 <- readRDS("../0.data/GSE178318.obj.rds")
GSE151530 <- readRDS("../0.data/GSE151530.obj.rds")
GSE166635 <- readRDS("../0.data/GSE166635.obj.rds")
GSE164522 <- readRDS("../0.data/GSE164522.obj.rds")
scCRLM <- readRDS("../0.data/scCRLM.obj.rds")
GSE132465 <- readRDS("../0.data/GSE132465.obj.rds")

dim(GSE125449)
dim(GSE140228)
dim(GSE112271)
dim(GSE178318)
dim(GSE151530)
dim(GSE166635)
dim(GSE164522)
dim(scCRLM)
dim(GSE132465)


GSE125449_list<-SplitObject(GSE125449,split.by = "Patient")
GSE140228_list<-SplitObject(GSE140228,split.by = "Patient")
GSE112271_list<-SplitObject(GSE112271,split.by = "Patient")
GSE178318_list<-SplitObject(GSE178318,split.by = "Patient")
GSE151530_list<-SplitObject(GSE151530,split.by = "Patient")
GSE166635_list<-SplitObject(GSE166635,split.by = "Patient")
GSE164522_list<-SplitObject(GSE164522,split.by = "Patient")
scCRLM_list<-SplitObject(scCRLM,split.by = "Patient")
GSE132465_list<-SplitObject(GSE132465,split.by = "Patient")

total_seu_list<-list(
  "GSE125449"=GSE125449_list,
  "GSE140228"=GSE140228_list,
  "GSE112271"=GSE112271_list,
  "GSE178318"=GSE178318_list,
  "GSE151530"=GSE151530_list,
  "GSE166635"=GSE166635_list,
  "GSE164522"=GSE164522_list,
  "scCRLM"=scCRLM_list,
  "GSE132465"=GSE132465_list
)
saveRDS(total_seu_list,"../0.data/combined_seu_list.RDS")

library(Seurat)
total_seu_list<-readRDS("../0.data/combined_seu_list.RDS")

dir.create("../0.data/1.QC")
QC_fun<-function(i){
  .libPaths(c("/home/yukai/R/4.0", 
              "/home/rstudio_cu02/anaconda3/envs/R4.0/lib/R/library", 
              "/home/yukai/miniconda3/envs/R41/lib/R/library"))
  library(Seurat)
  seu_obj<-data.lst[[i]]
  Patient<-unique(seu_obj$Patient)
  filter_10x(seu_obj,paste("../0.data/1.QC/",Patient,"/",sep=""))
}

library(doSNOW)
cl <- makeCluster(20)
registerDoSNOW(cl)
for (k in 1:length(total_seu_list)){
  data.lst<<-total_seu_list[[k]]
  iterations <- length(data.lst)
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  result <- foreach(i = 1:iterations, .combine = rbind, 
                    .options.snow = opts,
                    .packages=c("PCAtools","data.table","ggplot2","ggridges",
                                "dplyr","cowplot","harmony","RColorBrewer","DoubletFinder")
  ) %dopar% QC_fun(i)
  close(pb)
}
stopCluster(cl) 

####################1.3.integrate####################
library(PCAtools)
library(Seurat)
library(data.table)
library(ggplot2)
library(ggridges)
library(dplyr)
library(cowplot)
library(harmony)
library(RColorBrewer)
library(DoubletFinder)

if (T){
  system("ls ../0.data/1.QC/ > name.list")
  ##intergrate data
  file_info<-data.frame(data.table::fread("name.list",header = F))
  colnames(file_info)<-c("Sample_ID")
  sample_name_list<-file_info$Sample_ID
  k=0
  data.lst<-list()
  for (i in 1:length(sample_name_list)){
    #i=1
    Sample_ID<-sample_name_list[i]
    path<-paste("../0.data/1.QC/",Sample_ID,"/Exp.pre_precessed.rds",sep="")
    if (file.exists(path)){
      Sc_Exp<-readRDS(path)
      Sc_Exp@meta.data<-Sc_Exp@meta.data
      #Sc_Exp$Anatomic.region[which(Sc_Exp$Anatomic.region %in% c("Ascending colon","Caecum","Descending colon","Distal Descending colon","Sigmoid colon","Transverse colon"))]<-"Colon"
      #Patient_ID=unique(Sc_Exp$Patient)
      #Study<-unique(Sc_Exp$Dataset)
      #MSI<-unique(Sc_Exp$MSI)
      #Sample.type=unique(Sc_Exp$Type)
      print(ncol(Sc_Exp))
      #Sc_Exp$Type<-Sc_Exp$Sample.type
      #print(unique(Sc_Exp$Sample.type))
      if ((ncol(Sc_Exp)>=500)){
        k=k+1
        data.lst[[k]]<-Sc_Exp
      }
    }
    
  }
  
  dir.create("../1.integrate/")
  #saveRDS(data.lst, file = "2.intergrate/scExp_filter_list.rds")
  
  k=0
  result_table<-data.frame()
  for (i in data.lst){
    k=k+1
    print(paste(k,length(i$Patient),unique(i$Patient)))
    tmp_table<-data.frame(order=k,
                          Dataset= unique(i$Dataset),
                          #Type=unique(i$Type),
                          Patient=unique(i$Patient),
                          number=length(i$Patient))
    result_table<-rbind(result_table,tmp_table)
  }
  library(dplyr)
  
  Study_list<-unique(result_table$Dataset)
  data.lst1<-list( "GSE125449" =list()  ,    
                   "GSE140228"    =list()  ,   
                   "GSE112271" =list() , 
                   "GSE178318"       =list()     ,    
                   "GSE151530"  =list() ,   
                   "GSE166635"   =list()   ,   
                   "GSE164522"    =list() ,
                   "scCRLM"    =list()  ,
                   "GSE132465"    =list()
  )
  k=0
  for (i in data.lst){
    Dataset=unique(i$Dataset)
    data.lst1[[Dataset]][length(data.lst1[[Dataset]])+1]<-i
  }
  data.lst2<-list()
  for (Study in Study_list){
    print(Study)
    print(length(data.lst1[[Study]]))
    tmp_data.lst<-data.lst1[[Study]]
    if (length(tmp_data.lst) > 1) {
      data.lst2[[Study]]<-merge(tmp_data.lst[[1]],tmp_data.lst[2:length(tmp_data.lst)])
    }else{
      data.lst2[[Study]]<- tmp_data.lst[[1]]
    }
  }
  data.lst<-data.lst2
  
  ## Normalization
  warning("############## Normalization ##############")
  #data.lst <- list(exp_P1,exp_P2,exp_P3,exp_P4,exp_P5,exp_P6,exp_P7)
  data.lst <- lapply(X = data.lst, FUN = function(x){
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  features <- SelectIntegrationFeatures(object.list = data.lst,nfeatures =3000 )
  data.lst <- lapply(X = data.lst, FUN = function(x){
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  #data.lst<-list( data.lst[1:29],data.lst[31:length(data.lst)])
  ## Integration
  warning("############## Integration ##############")
  #data.lst<-data.lst[-285]
  anchors <- FindIntegrationAnchors(object.list = data.lst, anchor.features = features, reduction = "rpca",
                                    reference = 7,
                                    dims = 1:50)
  combined <- IntegrateData(anchorset = anchors,
                            dims = 1:50)
  DefaultAssay(combined) <- "integrated"
  saveRDS(combined, file = "../1.integrate/integrated.rds")
  
}


####################1.4.overview####################
##
combined<-readRDS("../1.integrate/integrated.rds")
combined <- subset(combined, subset = Tissue %in% c("MT", "MN", "Tumor", "CRC", "LM", "Normal", 
                                                    "Colon_T", "Liver_P", "Liver_T", "PT", "NC"))
table(combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "MT", "LM", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "MN", "NL", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "Tumor", "PL", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "CRC", "PC", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "LM", "LM", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "Normal", "NL", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "Colon_T", "PC", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "Liver_P", "NL", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "Liver_T", "LM", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "PT", "PC", combined$Tissue)
combined$Tissue <- ifelse(combined$Tissue == "NC", "CRC", combined$Tissue)

#combined$Type<-combined$Anatomic.region
combined <- ScaleData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
ElbowPlot(combined)
chosen.elbow <- findElbowPoint(Stdev(object = combined, reduction = "pca"))
combined <-RunHarmony(combined,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:15)
combined <- RunTSNE(combined, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:15)
combined <- FindClusters(combined, resolution = 0.6)
#combined<-readRDS("../1.integrate/Combined.rds")
DimPlot(combined, reduction = "umap", label = TRUE, raster=FALSE)

## Save
warning("############## Save data ##############")
saveRDS(combined, file = "../1.integrate/Combined.rds")


## Visualization
warning("############## Visualization ##############")
combined <- readRDS("../1.integrate/Combined.rds")
pbmc<-combined
#p0 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", raster=FALSE)
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "Patient", raster=FALSE)
p3 <- DimPlot(pbmc, reduction = "umap", group.by = "Dataset", raster=FALSE)
pdf("../1.integrate/Merged.features.umap.pdf", width = 22, height = 6)
p2|p3
dev.off()

p2 <- DimPlot(pbmc, reduction = "tsne", group.by = "Patient", raster=FALSE)
p3 <- DimPlot(pbmc, reduction = "tsne", group.by = "Dataset", raster=FALSE)
pdf("../1.integrate/Merged.features.tsne.pdf", width = 22, height = 6)
p2|p3
dev.off()

umap.out = cbind("Barcode" = rownames(Embeddings(object = combined, reduction = "umap")), Embeddings(object = combined, reduction = "umap"))
#dir.create("../1.integrate/")
write.table(umap.out, file="../1.integrate/Umap_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
tsne.out = cbind("Barcode" = rownames(Embeddings(object = combined, reduction = "tsne")), Embeddings(object = combined, reduction = "tsne"))
write.table(tsne.out, file="../1.integrate/Tsne_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
cluster.out = cbind("Barcode" = rownames(combined@meta.data), "Cluster" = combined@meta.data$seurat_clusters)
write.table(cluster.out, file="../1.integrate/Cluster_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)

combined_down <- subset(x = combined, downsample = 1000)
markers <- FindAllMarkers(combined_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
write.table(top10, "../1.integrate/Combined.markers.top10.txt",row.names = F, col.names = T, sep = "\t", quote = F)
write.table(markers, "../1.integrate/Combined.markers.all.txt",row.names = F, col.names = T, sep = "\t", quote = F)

#markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
#write.table(top10, "../1.integrate/Combined.markers.top10.txt",row.names = F, col.names = T, sep = "\t", quote = F)
#write.table(markers, "../1.integrate/Combined.markers.all.txt",row.names = F, col.names = T, sep = "\t", quote = F)

#pdf("../1.integrate/Combined_Cluster_marker_heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
#print(DoHeatmap(combined, features = top10$gene,label = F))
#dev.off()

#pdf ("../1.integrate/Cluter_DotPlot.pdf", height = 8, width = 40)
#print(DotPlot(combined, features = unique(top10$gene),
#              cols = c("grey","blue"))+RotatedAxis()+
#        scale_x_discrete("")+scale_y_discrete(""))
#dev.off()


###########################################2.Annotation###########################################
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(sscVis)
library(ggplot2)
library(harmony)
library(cowplot)

####################2.1.MajorCluster####################
dir.create("../2.Annote/")
combined <- readRDS("../1.integrate/Combined.rds")
DimPlot(combined, reduction = "umap", label = TRUE, raster=TRUE)
pdf("../2.Annote/MajorCluster.markers.pdf",width = 20,height = 25)
FeaturePlot(combined,features = c("MS4A1", "CD19", "CD79A", ##B cells
                                  "SDC1", "MZB1", ##Plasma cells
                                  "CD3D", ## T  cells
                                  "CD4", ## CD4+ T cells
                                  "FOXP3", "IL2RA", "CTLA4", ## CD4+ Treg cells
                                  "CD8A", ## CD8+ T cells
                                  "GNLY", "KLRD1", ## NK
                                  "IDO1", "LAMP3", "CLEC9A", "XCR1", "CLNK", ##cDC1 DC
                                  "CD1C", "CLEC10A", ##cDC2
                                  "CCR7", "CCL22", ##DC mature
                                  "LILRA4", "GPR183", "IL3RA", "CLEC4C", ##pDC DC cells
                                  "CD68", "APOE", "CD5L", "MARCO", "C1QB", ##Macrophage
                                  "MS4A2", "TPSB2", ##Mast cells,
                                  "CD14", "VCAN", "FCN1", ##Monocyte
                                  "CSF3R", "FCGR3B", ##Neutrophils
                                  "PECAM1", "VWF", "CDH5", "SELE", ##Endothelial  cells
                                  "CCL21", ##Endothelial  cells Lymphatic
                                  "PDGFRA", "FAP", "COL1A1", ##Fibroblast
                                  "ACTG2", ##MyoFibroblast
                                  "MSLN", "CALB2", ##Mesothelial cells,
                                  "KRT8", "EPCAM", "CDH1", "CLDN6", "CLDN7", ##Epithial cells,
                                  "MKI67", "CDK1",  ##prolification cells,
                                  "FYXD2", "TM4SF4", "ANXA4", ##cholangiocytes
                                  "APOC3", "FABP1", "APOA1", ##hepatocytes,
                                  "CD133", "SOX9"  ##Hepatic progenitor
),ncol = 6,min.cutoff = 0,raster = T)
dev.off()

FeaturePlot(combined,features = c("CD3D"),min.cutoff = 0,raster = T)

markers <- read.delim("../1.integrate/Combined.markers.all.txt", header = T)
top20genes <- top_n(group_by(markers,cluster), n = 20, wt = avg_log2FC)
panglaomarkers <- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv", 
                             header = T)
i = 0
cluster2celltype <- data.frame()
for (i in unique(top20genes$cluster)) {
  newres <- merge(panglaomarkers, top20genes[top20genes$cluster == i, ], by.x = "official.gene.symbol", by.y = "gene")
  newtmp <- as.data.frame(table(newres$cell.type))
  newtmp <- newtmp[order(newtmp$Freq, decreasing = T), ]
  tmpres <- data.frame(Cluster = i,
                       Cells = newtmp[1:5, ]$Var1,
                       COunts = newtmp[1:5, ]$Freq)
  cluster2celltype <- rbind(cluster2celltype, tmpres)
}
View(cluster2celltype)

#combined <- subset(combined, subset = seurat_clusters != 10)

majorcluster <- c("NK", "T", "T", "T", "NK", 
                  "Macro", "T", "Epi", "Plasma", "B", 
                  "Mono", "T", "Epi", "Prolif", "DC",
                  "Epi", "Fibro", "Endo", "DC", "Mast", 
                  "Unknown", "Macro", "Prolif", "T", "Epi", 
                  "Unknown", "Unknown", "Epi", "Prolif")

names(majorcluster)<-levels(combined)
combined<-RenameIdents(combined,majorcluster)
combined$MajorCluster<-Idents(combined)
combined <- subset(combined, subset = MajorCluster != "Unknown")
DimPlot(combined, reduction = "umap", label = TRUE, raster=TRUE)

saveRDS(combined,"../2.Annote/Merged.pca.major.rds")     

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
pdf("../2.Annote/Vlnplot.pdf",width = 24,height = 30)
VlnPlot(obj = combined , features = c("MS4A1", "CD19", "CD79A", ##B cells
                                      "SDC1", "MZB1", ##Plasma cells
                                      "CD3D", ## T  cells
                                      "CD4", ## CD4+ T cells
                                      "FOXP3", "IL2RA", "CTLA4", ## CD4+ Treg cells
                                      "CD8A", ## CD8+ T cells
                                      "GNLY", "KLRD1", ## NK
                                      "IDO1", "LAMP3", "CLEC9A", "XCR1", "CLNK", ##cDC1 DC
                                      "CD1C", "CLEC10A", ##cDC2
                                      "CCR7", "CCL22", ##DC mature
                                      "LILRA4", "GPR183", "IL3RA", "CLEC4C", ##pDC DC cells
                                      "CD68", "APOE", "CD5L", "MARCO", "C1QB", ##Macrophage
                                      "MS4A2", "TPSB2", ##Mast cells,
                                      "CD14", "VCAN", "FCN1", ##Monocyte
                                      "CSF3R", "FCGR3B", ##Neutrophils
                                      "PECAM1", "VWF", "CDH5", "SELE", ##Endothelial  cells
                                      "CCL21", ##Endothelial  cells Lymphatic
                                      "PDGFRA", "FAP", "COL1A1", ##Fibroblast
                                      "ACTG2", ##MyoFibroblast
                                      "MSLN", "CALB2", ##Mesothelial cells,
                                      "KRT8", "EPCAM", "CDH1", "CLDN6", "CLDN7", ##Epithial cells,
                                      "MKI67", "CDK1",  ##prolification cells,
                                      "FYXD2", "TM4SF4", "ANXA4", ##cholangiocytes
                                      "APOC3", "FABP1", "APOA1", ##hepatocytes,
                                      "CD133", "SOX9"  ##Hepatic progenitor
),ncol = 6,pt.size = 0)
dev.off()


p <- DimPlot(combined, reduction = "umap", label = TRUE, raster=TRUE)
p
ggsave("../2.Annote/MajorCluster_merge_umap.pdf",p,width = 5,height = 4 )

p <- DimPlot(combined, reduction = "umap", group.by = "Dataset", label = TRUE, raster=TRUE)
p
ggsave("../2.Annote/MajorCluster_merge_umap.study.pdf",p,width = 5,height = 4 )

#markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.table(markers,"../2.Annote/MajorCluster_merge.markers.txt",sep="\t",row.names=T,col.names=T,quote=F)
#top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
#pdf("../2.Annote/MajorCluster_merge.markers.heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
#DoHeatmap(combined, features = top10$gene,label = F)
#dev.off()

#pdf ("../2.Annote/MajorCluster_merge.dotplot.pdf", height = 8, width = 30)
#DotPlot(combined, features = unique(top10$gene),
#        cols = c("grey","blue"))+RotatedAxis()+
#  scale_x_discrete("")+scale_y_discrete("")
#dev.off()

####################2.2.Minorcluster####################
###########2.2.1.T and NK cells###########
gene_list<-unique(
  c("CD3D","CD3E","AC016831.7","AC017002.1","AC069363.1","AC092580.4","AC133644.2","ACTG1",
    "ADAM19","ALOX5AP","AMICA1","ANXA1","ANXA2","APOBEC3C","AQP3","AREG",
    "ARID5B","ATXN1","BATF","BATF3","BAX","BIRC3","BNIP3",
    "BOLA2B","C12orf75","C16orf54","C1orf56","CALR","CCL17",
    "CCL20","CCL3","CCL3L3","CCL4","CCL4L2","CCL5","CCR7","CD2",
    "CD4","CD40LG","CD52","CD63","CD69","CD7","CD74","CD84","CD8A",
    "CD8B","CD9","CDC42SE1","CDK2AP2","CEBPD","CHST12","CLC","CLDND1",
    "CLEC2B","CLIC1","CLIC3","CMC1","CORO1A","COTL1","CREM","CRIP1","CRTAM",
    "CSF2","CST7","CTC-425F1.4","CTLA4","CTNNB1","CTSC","CTSH","CTSL","CTSW",
    "CXCL10","CXCL13","CXCR4","CXCR5","CXCR6","CYBA","CYCS","CYSLTR1","DACT1","DDIT4","DDX21","DUSP2",
    "DUSP4","EFHD2","EGR1","EGR2","EGR3","EIF5A","EMP3","ENO1","EOMES","EVI2A","EVI2B",
    "F2R","FABP5","FCGR3A","FGFBP2","FLJ21408","GADD45B","GAS5","GCNT4","GIMAP4",
    "GIMAP7","GLRX","GNG2","GNG8","GNLY","GPR155","GPR171","GPR183","GZMA","GZMB","GZMH",
    "GZMK","HAVCR2","HBA1","HBA2","HBB","HCST","HILPDA","HNRNPH1","HOPX","HSP90AB1","HSPA8","HSPE1","ICOS","ID2",
    "IER3","IFI27","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","IFITM2","IFNG","IL10",
    "IL13","IL17A","IL17F","IL2","IL22","IL2RA","IL3","IL31","IL32","IL4","IL4R","IL5",
    "IL7R","IL9","INSIG1","IQCG","IRF1","IRF4","ISG15","ISG20","ITGB1","JUN",
    "JUNB","KLF2","KLF6","KLRB1","KLRC1","KLRC2","KLRD1","KLRG1","LAG3","LDHA","LDHB",
    "LGALS1","LGALS3","LINC00152","LINC00861","LPAR6","LSP1","LTA","LTB","LY6E","LYAR"
    ,"LYST","MAF","MAL","MDM4","MIAT","MIF","MIR155HG","MIR4435-2HG","MT2A","MTA2",
    "MTCO1P12","MTCO2P12","MT-ND5","MT-RNR1","MT-TE","MT-TM","MX1","MX2","MYC","MYO1F",
    "MYO1G","NCL","NEAT1","NFKBIA","NFKBID","NKG7","NME1","NOSIP","NR3C1","NR4A1","NR4A2",
    "OAS1","OAS3","ODC1","PABPC1","PASK","PATL2","PDCD4","PEBP1P3","PGK1","PHKG1","PHLDA1",
    "PIK3IP1","PIM2","PIM3","PKM","PLAC8","PMCH","PPP1R14B","PRDM1","PRF1","PTGER2","PTMS",
    "PYHIN1","RARRES3","RCAN2","RGCC","RGS1","RGS16","RNF213","RP11-138A9.1","RP11-138A9.2",
    "RP11-1399P15.1","RP11-291B21.2","RP11-347P5.1","RP11-693N9.2","RPL22L1","RPS27L","RSAD2",
    "RTKN2","S100A10","S100A11","S100A4","S100A6","SAMD9","SAMD9L","SAT1","SDF2L1","SEC11C",
    "SELL","SESN3","SET","SGK1","SH3BGRL3","SLFN5","SOCS1","SOCS2","SRGN","SRM","STAT1","SYNE2",
    "TAGAP","TC2N","TIGIT","TMEM173","TNF","TNFRSF18","TNFRSF4","TNFRSF9","TNFSF10","TNIP3","TNRC6B",
    "TPI1","TPM4","TPR","TRAC","TRBC1","TRBC2","TRDC","TRGC1","TRGC2","TSC22D3","TSHZ2","TUBA1B","TXN",
    "TYROBP","USP10","VIM","VTRNA1-3","XAF1","XBP1","XCL1","XCL2","YBX3","YPEL3","YWHAH","ZBED2","ZEB2",
    "ZFP36","ZFP36L2","ZNF683","PDCD1",
    "CTLA4","FOXP3", #Treg
    "TBX21","CXCR3", #Th1, T-bet
    "RORC","RORA","IL23R","CCR6", #Th17,RORγt
    "CCR10", #Th22
    "CXCR5","ICOS","PDCD1", #Tfh
    "ITGA1","CXCR6","CCL5", #Trm
    "AXNA1","CCR7","TCF7","IL7R", #Tcm
    "LTB","LEF1","MAL", #Tn
    "GZMK","GZMA","GZMH","HLA-DPA1","HLA-DPB1","HLA-DRB1","HLA-DRB6", ##Tem
    "NKG7","PRF1","GNLY",  ##Temra
    "LAG3","CXCL13","PDCD1","TIGIT","HAVCR2", #Tex
    "SLC4A10","CEBPD", ##MAIT
    "MKI67","TUBA1B","STMN1","HMGB2","HMGN2" ##Prolif
  )
)

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
Idents(combined)<-combined$MajorCluster
sub_sce<-subset(combined,idents=c("T", "NK"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"
dir.create("../2.Annote/2.1.Tcells")
saveRDS(sub_sce, file = "../2.Annote/2.1.Tcells/Merged.pca.major.Tcells.rds")

sub_sce <- ScaleData(sub_sce, verbose = FALSE,features = gene_list)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 1000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE,features = gene_list)

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 1.5)
DimPlot(sub_sce, reduction = "umap", label = TRUE)
DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
umap.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "umap")), 
                 Embeddings(object = sub_sce, reduction = "umap"))

write.table(umap.out, file="../2.Annote/2.1.Tcells/Umap_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
tsne.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "tsne")), 
                 Embeddings(object = sub_sce, reduction = "tsne"))
write.table(tsne.out, file="../2.Annote/2.1.Tcells/Tsne_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
cluster.out = cbind("Barcode" = rownames(sub_sce@meta.data), 
                    "Cluster" = sub_sce@meta.data$seurat_clusters)
write.table(cluster.out, file="../2.Annote/2.1.Tcells/Cluster_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
## Save
warning("############## Save data ##############")
saveRDS(sub_sce, file = "../2.Annote/2.1.Tcells/Merged.pca.major.Tcells.pca.rds")


## Visualization
warning("############## Visualization ##############")
sub_sce<-readRDS("../2.Annote/2.1.Tcells/Merged.pca.major.Tcells.pca.rds")
pdf("../2.Annote/2.1.Tcells/Cluster_merge_umap.pdf", height = 7, width = 18)
p1<-DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
p2 <- DimPlot(sub_sce, reduction = "umap", label = TRUE)
print(plot_grid(p1,p2))
dev.off()
pdf("../2.Annote/2.1.Tcells/Cluster_merge_umap2.pdf", height = 7, width = 18)
p1 <- DimPlot(sub_sce, reduction = "umap", split.by = "MajorCluster",label = TRUE)
print(plot_grid(p1))
dev.off()

sub_sce_down <- subset(x = sub_sce, downsample = 1000)
markers <- FindAllMarkers(sub_sce_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../2.Annote/2.1.Tcells/sub_sce_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)

## Visualization
warning("############## Visualization ##############")
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
pdf("../2.Annote/2.1.Tcells/sub_sce_Cluster_marker_heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
print(DoHeatmap(sub_sce, features = top10$gene,label = F))
dev.off()

pdf ("../2.Annote/2.1.Tcells/Cluter_DotPlot.pdf", height = 8, width = 30)
print(DotPlot(sub_sce, features = unique(c(top10$gene,"ICOS")),
              cols = c("grey","blue"))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete(""))
dev.off()

## minor level clustering
warning("############## Annotation genes ##############")
sub_sce<-readRDS("../2.Annote/2.1.Tcells/Merged.pca.major.Tcells.pca.rds")
Idents(sub_sce)<-sub_sce$seurat_clusters

DimPlot(sub_sce, reduction = "umap", label = TRUE)
pdf("../2.Annote/2.1.Tcells/replot.markers.pdf",width = 10,height = 8)
FeaturePlot(sub_sce,features = c("CD3D", ## T  cells
                                 "CD4", ## CD4+ T cells
                                 "FOXP3", "IL2RA", "CTLA4", ## CD4+ Treg cells
                                 "CD8A", ## CD8+ T cells
                                 "GNLY", "KLRD1" ## NK
),ncol = 3,min.cutoff = 0,raster = T)
dev.off()

markers <- read.delim("../2.Annote/2.1.Tcells/sub_sce_Cluster_markers.txt", header = T)
top20genes <- top_n(group_by(markers,cluster), n = 40, wt = avg_log2FC)
panglaomarkers <- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv", header = T)
i = 0
cluster2celltype <- data.frame()
for (i in unique(top20genes$cluster)) {
  newres <- merge(panglaomarkers, top20genes[top20genes$cluster == i, ], by.x = "official.gene.symbol", by.y = "gene")
  newtmp <- as.data.frame(table(newres$cell.type))
  newtmp <- newtmp[order(newtmp$Freq, decreasing = T), ]
  tmpres <- data.frame(Cluster = i,
                       Cells = newtmp[1:5, ]$Var1,
                       COunts = newtmp[1:5, ]$Freq)
  cluster2celltype <- rbind(cluster2celltype, tmpres)
}
View(cluster2celltype)

top10genes <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
makgs <- c("CTLA4","FOXP3", #Treg
           "TBX21","CXCR3", #Th1, T-bet
           "RORC","RORA","IL23R","CCR6", #Th17,RORγt
           "CCR10", #Th22
           "CXCR5","ICOS","PDCD1", #Tfh
           "ITGA1","CXCR6","CCL5", #Trm
           "AXNA1","CCR7","TCF7","IL7R", #Tcm
           "LTB","LEF1","MAL", #Tn
           "GZMK","GZMA","GZMH","HLA-DPA1","HLA-DPB1","HLA-DRB1","HLA-DRB6", ##Tem
           "NKG7","PRF1","GNLY",  ##Temra
           "LAG3","CXCL13","PDCD1","TIGIT","HAVCR2", #Tex
           "SLC4A10","CEBPD", ##MAIT
           "MKI67","TUBA1B","STMN1","HMGB2","HMGN2" ##Prolif
)
intersect(top20genes[top20genes$cluster == 23, ]$gene, makgs)
new.cluster.id<-c("CD8T_MAIT_CEBPD", "CD8T_Tem_GZMK", "CD8T_Trm_ITGA1", "CD8T_Tfh_PDCD1",
                  "NK_NK_GNLY", "NK_NK_GNLY", "CD4T_Treg_FOXP3", "CD8T_Tem_GZMK", 
                  "CD4T_Tcm_IL7R", "CD4T_Tcm_IL7R", "CD8T_Tem_GZMK", "CD4T_Tn_MAL", 
                  "CD8T_Trm_ITGA1", "CD8T_Tcm_CCR7", "CD4T_Treg_FOXP3", "CD4T_Tcm_IL7R", 
                  "CD4T_Trm_CXCR6", "CD8T_Tem_GZMK", "CD4T_Treg_CTLA4", "CD8T_NKT_GNLY", 
                  "CD4T_Tcm_CCR7", "UN_UN_UN", "NK_NK_NKG7", "CD4T_Tn_MAL", 
                  "UN_UN_UN", "UN_UN_UN") 

names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)

sub_sce$MinorCluster <- Idents(sub_sce)
sub_sce$Celltype <- strsplit2(sub_sce$MinorCluster, "_")[, 1]
sub_sce$CellClass <- strsplit2(sub_sce$MinorCluster, "_")[, 2]
sub_sce$Genes <- strsplit2(sub_sce$MinorCluster, "_")[, 3]

sub_sce <- subset(sub_sce, subset = CellClass != "UN")

tmpdf <- data.frame(MinorCluster = sub_sce$MinorCluster,
                    Celltype = sub_sce$Celltype,
                    CellClass = sub_sce$CellClass)
tmpdf <- tmpdf[order(tmpdf$Celltype, tmpdf$CellClass), ]

Idents(sub_sce)<-factor(Idents(sub_sce),levels= as.character(unique(tmpdf$MinorCluster)))   
saveRDS(sub_sce, file = "../2.Annote/2.1.Tcells/Merged.pca.major.Tcells.pca.minor.rds")

## see annotations
warning("############## see annotations ##############")
p <- DimPlot(sub_sce, reduction = "umap", group.by = "MinorCluster", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.1.Tcells/minor_Cluster_merge_umap.pdf",p,width = 8,height = 5 )

p <- DimPlot(sub_sce, reduction = "umap", group.by = "Dataset", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.1.Tcells/minor_Cluster_study_umap.pdf",p,width = 12,height = 9 )

p1 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.1.Tcells/minor_Cluster_cellclass_umap.pdf",p1,width = 8,height = 7 )

###########2.2.2.B cells###########
gene_list<-c(
  ##Class switch
  "CD27","TNFRSF13B","EPS15","TAF6","UBE2I","SUB1","SCRN1","AFTPH",
  "TRRAP","ZBTB32","GPR25","SNED1","HLA-DQB2","SPIB","CD80","BLK",
  "NDUFA9","FCRL2","UBE2N","NARFL","SP140","PTPN6","CPSF4","TRAF3",
  "SLC12A3","UBE2G1","PRDM10","KHDRBS2","MS4A1","TNFRSF17","COL19A1",
  "PAX5","CXCR5","RAD17","DEPDC5","NGLY1","RAPGEF1","ABI1","BAIAP3",
  "PIKFYVE","CR1","SEC24A","ADAMDEC1",
  ##Memory
  "GABRA4","GPX5","SLC17A1","CYLC2","DSCR4","SHISA6","NPHS2","CXCL13",
  "TSHB","PRKCB","ANKRD34C","CPB1","KIF5A","KCNA5","CD3EAP","DDX4","PLIN1",
  "SIGLEC6","SLC12A3","GPRC5D","MS4A1","PNOC","ZIC3","CCR6","PAX5","S100G",
  "SLC17A7","SERPINA4","CXCR5","KIAA0125","GAD2","PIKFYVE","NPY5R","ADAM20",
  "ADAM30","SYN2","ZNF548","RIC3","MYOZ3","CD72","SELP","GNRHR","SCRN1","ZBTB32",
  "MAP3K9","CHRNA2","CD180","UNC5C","ULK4","VPREB3","SLCO1C1","FCRL2","DSP","POU4F2",
  "PRPH2","TSPAN13","CD22","MS4A5","CD37","CHRM2","GRIN2B","AP1M2","CSN1S1","STAP1",
  "MBD4","RNGTT","TNFRSF17","AQP8","CPA2","ART1","HLA-DPB1","CYP2C19","GRM6","HSD3B2",
  "CER1","MC4R","GK2","CD19","CNR2","DPP6","CHST5","TAS2R14","OTC","ZNF747","HNRNPL",
  "GGA2","SPIB","TRPM3","CLDN17","BLK","DCC","SSX3","FMO6P","GLYAT","MIOS","NTRK3",
  "MEFV","CD1C","CAPN3","NCAN","KRT75","SP140","VN1R1","MBL2","CYP2A7","HECW1","OBSCN",
  "FSHR","CHP2","KALRN","TCTN2","COLEC10","TRMT61A","HTN3","MOGAT2","HRH4","CETP","WNT16",
  "ADAMTS12","LECT2","SLN","INHBC","CCR9","TNFRSF13B","TMPRSS11D","FMO1","MGAT5","CD79B",
  "RRH","CASQ2","ADCY2","NT5C","QRSL1","HTR3A","ADAM21","LY86","C4BPA","SLC30A10","KRT2",
  "PNLIPRP1","PROZ","HCRTR2","PTH1R","KHDRBS2","KCNJ10","SLC24A2","GNAT2","WNT2","AICDA",
  "MEP1B","CRB1","SLC5A7","FSCN2","PGLYRP4","ODC1","CD79A","BAIAP3","ACRV1","SYPL1","FSCN3","COX6A2","AIPL1",
  ##Naive
  "CD24","MARCH8","DAZL","PRKCB","FRS2","GMFB","TRA2B","BCL2L11","HLA-DOA","GCM1",
  "CIITA","STAG3","MS4A1","PNOC","CCR6","COL19A1","PGAM2","ZNF154","PAX5","CXCR5",
  "PWP1","SLC30A4","CSNK1G3","CD1A","AKAP6","PIKFYVE","FCER2","ADAM20","CD72","TCL1B",
  "USP7","MAP3K9","POU2F1","CUBN","CD180","TCL1A","VPREB3","FCRL2","DSP","RB1","TCL6",
  "TSPAN13","CD22","CD37","SNX2","STAP1","CDK13","TBC1D5","HSPA4","GNG3","MBD4","C10orf76",
  "AP3B1","MCM9","MFN1","DEF8","CD19","SIPA1L3","RBM15","GGA2","P2RY10","TREML2","WDR74",
  "SPIB","MMP17","BLK","RERE","CACNA1F","MYBPC2","NOC3L","PTCH2","CAPN3","SP140","SDK2",
  "PRDM4","RRAS2","LY9","GH1","MATN1","EGOT","PHKG1","N4BP3","SNTG2","PYGM","MYO3A","MGAT5",
  "CD79B","USP6","TRAPPC9","PRDM2","BMP3","BCL2L10","GPR18","KHDRBS2","SMC6","CD79A","SYN3",
  "UTP6","UBE2O","CXCL1","CXCL2","CXCL6","CXCL8","CXCL9","CXCL10","CXCL13","CXCL14",
  ##Plasma
  "MTDH","PRDX4","CCNC","GPR37L1","B4GALT3","DNAJC4","RAB3A","IMP4","MYL2","CAMP",
  "LTB4R2","VPREB1","ISCU","TIMM17B","C19orf73","CDH15","FNDC3A","DDOST","CCDC40",
  "GOLGA4","SRPRB","MGAT2","RGS13","STMN4","NOS2","SAP30BP","SEC61A1","PPIB","IDE",
  "ALG5","ZBP1","LMAN2","CACNA1S","WNT1","EBAG9","R3HCC1","CA7","GRM4","BMP8B","WDR45",
  "CHRNA4","GOLGB1","PDE6A","UBA5","UBE2G1","HDLBP","GDF2","AVIL","PNOC","TP73","RPN2",
  "RAX","SCFD1","GPRC5D","UGGT1","PICK1","PPIL2","TERT","PRX","RFX2","C6orf25","SLC13A2",
  "HIST1H2BB","NGLY1","ZDHHC4","C16orf58","SPATS2","KIAA0125","CSHL1","BFSP2","MANF","MIA3",
  "KDELR2","NPPC","TNNT3","ACBD3","RASIP1","OGFOD2","UBXN4","HOOK2","PTGER1","GUCA1A","UTF1",
  "SLC35B1","TRABD","CD180","TSHR","CYP11A1","USP48","RNF113A","CELA2B","NTRK1","VPREB3","ATF6",
  "FCRL2","CNKSR1","KLF15","GRIN1","ADM2","FTCD","LBX1","UFSP2","CCDC88A","HSF4","RALY","NEUROG1",
  "FBP2","ACBD4","CRYGC","RNGTT","APOA1","SIX5","TNFRSF17","DDN","SEC61G","EPO","CHRNG","PCDHA5",
  "CUX2","HSPA6","CRYBB3","PDIA6","TBL2","CUL7","GABRR2","FGF6","PDIA2","CASP10","KRT10","HAND2",
  "IRGC","ENTPD1","SLCO5A1","CD19","IGF1","DEF8","GOLGA3","MCTS1","CAV1","SNAPC4","APOC3","SIPA1L3",
  "RGS1","DPAGT1","SLC35C2","ARL1","SLC6A13","PABPC4","CLCNKB","GH2","MAPK8IP3","ARSA","GNB3","TM9SF1",
  "ALPI","MATN4","GLT8D1","SLC38A10","SHANK1","T","CRYBA4","PRM1","DRD5","SEC62","ABCB9","DNASE1L2",
  "CRYBB1","IRF4","SSR4","RNF103","SERP1","SEC63","TMEM39A","CHST8","TMEM59","SEC61B","POMC","DKKL1",
  "DRD4","GMPPA","SSR1","FN3K","SYT5","RAD17","GRWD1","POU3F3","MYH13","VSX1","TMED10","TBX4","ITGA8",
  "MYL7","MAST1","TCF3","SEMA6C","SEC24A","NPPA","NDOR1","AVP","PHOX2A","RPN1","CNR1","RNF208",
  "ZNF133","ZPBP","YIPF1","ZNF37A","CCL25","CNTD2","KCNQ4","DMTF1","MARS","CLINT1","P2RY4","TNFRSF13B",
  "PREB","CNPY2","PTPRS","TMED9","CSPP1","IFT52","CD79B","HEYL","SMPD2","SLC5A2","TREH","GORASP2","LAX1",
  "SRP54","MBTPS1","AUP1","SHBG","YIPF2","CCDC33","KNTC1","C21orf2","AMPD1","KCNN3","CCDC121","ALG9","ELL",
  "GP9","LMAN1L","CYBA","HSP90B1","GPLD1","SURF1","CEACAM21","PRDM14","MRPS31","CD27","TSSK2","MIS12","TG",
  "SS18","ARHGEF16","CD79A","FKBP2","DAD1","ZNF142","MAGEF1","LEFTY2","KLC2","DOK3","SPCS1","THAP4",
  "ERGIC3","NPAS1","COX6A2","AIPL1",
  "IGHG1","IGHG2","IGHG4","IGHA1","IGHA2","IGHM","IGHD","IGHE",
  "MKI67","TUBA1B","STMN1","HMGB2","HMGN2", ##Prolif
  "IFI27","IFI44L","IFI6","IFIT1","IFIT2","IFIT3","IFITM2","IFNG","ISG15","ISG20",
  "GPR183","XIST")  

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
Idents(combined)<-combined$MajorCluster
sub_sce<-subset(combined,idents=c("B","Plasma"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"
dir.create("../2.Annote/2.2.Bcells")
saveRDS(sub_sce, file = "../2.Annote/2.2.Bcells/Merged.pca.major.Bcells.rds")

sub_sce <- ScaleData(sub_sce, verbose = FALSE,features = gene_list)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 3000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE,features = gene_list )

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 0.7)
DimPlot(sub_sce, reduction = "umap", label = TRUE)
DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
umap.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "umap")), 
                 Embeddings(object = sub_sce, reduction = "umap"))

write.table(umap.out, file="../2.Annote/2.2.Bcells/Umap_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
tsne.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "tsne")), 
                 Embeddings(object = sub_sce, reduction = "tsne"))
write.table(tsne.out, file="../2.Annote/2.2.Bcells/Tsne_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
cluster.out = cbind("Barcode" = rownames(sub_sce@meta.data), 
                    "Cluster" = sub_sce@meta.data$seurat_clusters)
write.table(cluster.out, file="../2.Annote/2.2.Bcells/Cluster_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
## Save
warning("############## Save data ##############")
saveRDS(sub_sce, file = "../2.Annote/2.2.Bcells/Merged.pca.major.Bcells.pca.rds")

## Visualization
warning("############## Visualization ##############")
sub_sce<-readRDS("../2.Annote/2.2.Bcells/Merged.pca.major.Bcells.pca.rds")
pdf("../2.Annote/2.2.Bcells/Cluster_merge_umap.pdf", height = 7, width = 18)
p1<-DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
p2 <- DimPlot(sub_sce, reduction = "umap", label = TRUE)
print(plot_grid(p1,p2))
dev.off()
pdf("../2.Annote/2.2.Bcells/Cluster_merge_umap2.pdf", height = 7, width = 18)
p1 <- DimPlot(sub_sce, reduction = "umap", split.by = "MajorCluster",label = TRUE)
print(plot_grid(p1))
dev.off()

sub_sce_down <- subset(x = sub_sce, downsample = 1000)
markers <- FindAllMarkers(sub_sce_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../2.Annote/2.2.Bcells/sub_sce_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)

## Visualization
warning("############## Visualization ##############")
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
pdf("../2.Annote/2.2.Bcells/sub_sce_Cluster_marker_heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
print(DoHeatmap(sub_sce, features = top10$gene,label = F))
dev.off()

pdf ("../2.Annote/2.2.Bcells/Cluter_DotPlot.pdf", height = 8, width = 30)
print(DotPlot(sub_sce, features = unique(c(top10$gene,"ICOS")),
              cols = c("grey","blue"))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete(""))
dev.off()

## minor level clustering
warning("############## Annotation genes ##############")
sub_sce<-readRDS("../2.Annote/2.2.Bcells/Merged.pca.major.Bcells.pca.rds")
Idents(sub_sce)<-sub_sce$seurat_clusters

DimPlot(sub_sce, reduction = "umap", label = TRUE)
pdf("../2.Annote/2.2.Bcells/replot.markers.pdf",width = 10,height = 6)
FeaturePlot(sub_sce,features = c("MS4A1", "CD19", "CD79A", ##B cells
                                 "SDC1", "MZB1", "VPREB1",
                                 "CD83" ##Plasma cells
),ncol = 3,min.cutoff = 0,raster = T)
dev.off()

markers <- read.delim("../2.Annote/2.2.Bcells/sub_sce_Cluster_markers.txt", header = T)
top20genes <- top_n(group_by(markers,cluster), n = 20, wt = avg_log2FC)
panglaomarkers <- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv", header = T)
i = 0
cluster2celltype <- data.frame()
for (i in unique(top20genes$cluster)) {
  newres <- merge(panglaomarkers, top20genes[top20genes$cluster == i, ], by.x = "official.gene.symbol", by.y = "gene")
  newtmp <- as.data.frame(table(newres$cell.type))
  newtmp <- newtmp[order(newtmp$Freq, decreasing = T), ]
  tmpres <- data.frame(Cluster = i,
                       Cells = newtmp[1:5, ]$Var1,
                       COunts = newtmp[1:5, ]$Freq)
  cluster2celltype <- rbind(cluster2celltype, tmpres)
}
View(cluster2celltype)

top10genes <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

new.cluster.id<-c("B_Bn_LTB", "Plasma_Plasma_MZB1", "B_Bn_CD72", "Plasma_Plasma_RRBP1", 
                  "Plasma_Plasma_PRDX4", "B_Bn_BANK1", "Plasma_Plasma_RRBP1", "B_Bm_CCR7", 
                  "B_Bm_CXCL2", "Plasma_Plasma_PTPRS", "Plasma_Plasma_SDF2L1", "B_Bn_IFI27", 
                  "B_Bm_RGS13", "Plasma_Plasma_APOA1", "B_Bm_CCR7", "UN_UN_UN", 
                  "UN_UN_UN") 

names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)

sub_sce$MinorCluster <- Idents(sub_sce)
sub_sce$Celltype <- strsplit2(sub_sce$MinorCluster, "_")[, 1]
sub_sce$CellClass <- strsplit2(sub_sce$MinorCluster, "_")[, 2]
sub_sce$Genes <- strsplit2(sub_sce$MinorCluster, "_")[, 3]

sub_sce <- subset(sub_sce, subset = CellClass != "UN")

tmpdf <- data.frame(MinorCluster = sub_sce$MinorCluster,
                    Celltype = sub_sce$Celltype,
                    CellClass = sub_sce$CellClass)
tmpdf <- tmpdf[order(tmpdf$Celltype, tmpdf$CellClass), ]

Idents(sub_sce)<-factor(Idents(sub_sce),levels= as.character(unique(tmpdf$MinorCluster)))                
saveRDS(sub_sce, file = "../2.Annote/2.2.Bcells/Merged.pca.major.Bcells.pca.minor.rds")
## see annotations
warning("############## see annotations ##############")
p <- DimPlot(sub_sce, reduction = "umap", group.by = "MinorCluster", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.2.Bcells/minor_Cluster_merge_umap.pdf",p,width = 8,height = 5 )

p <- DimPlot(sub_sce, reduction = "umap", group.by = "Dataset", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.2.Bcells/minor_Cluster_study_umap.pdf",p,width = 12,height = 9 )

p1 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.2.Bcells/minor_Cluster_cellclass_umap.pdf",p1,width = 8,height = 7 )

###########2.2.3.Myeloid cells###########
##Macrophage (CD68)
##M1 macro (CD38, GPR18, FPR2, HLA-DRA, ITGAX, CD86, NOS2)
##M2 macro (EGR2, MYC, CD163, VEGFA, MAF)
##cDC1 (IDO1, LAMP3)
##cDC2 (CD1C)
##pDC (LILRA4, GPR183)
##monocytes (FCGR3A, CD14, S100A9)
##neutrophils (CDK11B, FCGR3A, and MME)

allgenes <- c("CD68", "CD38", "GPR18", "FPR2", "HLA-DRA", "ITGAX", "CD86", "NOS2",
              "EGR2", "MYC", "CD163", "VEGFA", "MAF", "IDO1", "LAMP3", "CD1C",
              "LILRA4", "GPR183", "FCGR3A", "CD14", "CDK11B", "MME")


gene_list<- read.table("/home/yukai/projects/sclearn/NPC2023fromMei/0.data/Myeloid_markers.txt",header=F)
gene_list <- unique(c(gene_list$V1, allgenes))

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
Idents(combined)<-combined$MajorCluster
sub_sce<-subset(combined,idents=c("Macro", "Mono", "DC"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"
dir.create("../2.Annote/2.3.Myeloid")
saveRDS(sub_sce, file = "../2.Annote/2.3.Myeloid/Merged.pca.major.Myeloidcells.rds")

sub_sce <- ScaleData(sub_sce, verbose = FALSE,features = gene_list)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 3000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE,features = gene_list)

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 0.8)
DimPlot(sub_sce, reduction = "umap", label = TRUE)
DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
umap.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "umap")), 
                 Embeddings(object = sub_sce, reduction = "umap"))

write.table(umap.out, file="../2.Annote/2.3.Myeloid/Umap_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
tsne.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "tsne")), 
                 Embeddings(object = sub_sce, reduction = "tsne"))
write.table(tsne.out, file="../2.Annote/2.3.Myeloid/Tsne_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
cluster.out = cbind("Barcode" = rownames(sub_sce@meta.data), 
                    "Cluster" = sub_sce@meta.data$seurat_clusters)
write.table(cluster.out, file="../2.Annote/2.3.Myeloid/Cluster_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
## Save
warning("############## Save data ##############")
saveRDS(sub_sce, file = "../2.Annote/2.3.Myeloid/Merged.pca.major.Myeloidcells.pca.rds")

## Visualization
warning("############## Visualization ##############")
sub_sce<-readRDS("../2.Annote/2.3.Myeloid/Merged.pca.major.Myeloidcells.pca.rds")
pdf("../2.Annote/2.3.Myeloid/Cluster_merge_umap.pdf", height = 7, width = 18)
p1<-DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
p2 <- DimPlot(sub_sce, reduction = "umap", label = TRUE)
print(plot_grid(p1,p2))
dev.off()
pdf("../2.Annote/2.3.Myeloid/Cluster_merge_umap2.pdf", height = 7, width = 18)
p1 <- DimPlot(sub_sce, reduction = "umap", split.by = "MajorCluster",label = TRUE)
print(plot_grid(p1))
dev.off()

sub_sce_down <- subset(x = sub_sce, downsample = 1000)
markers <- FindAllMarkers(sub_sce_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../2.Annote/2.3.Myeloid/sub_sce_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)

## Visualization
warning("############## Visualization ##############")
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
pdf("../2.Annote/2.3.Myeloid/sub_sce_Cluster_marker_heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
print(DoHeatmap(sub_sce, features = top10$gene,label = F))
dev.off()

pdf ("../2.Annote/2.3.Myeloid/Cluter_DotPlot.pdf", height = 8, width = 30)
print(DotPlot(sub_sce, features = unique(c(top10$gene,"ICOS")),
              cols = c("grey","blue"))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete(""))
dev.off()

## minor level clustering
warning("############## Annotation genes ##############")
sub_sce<-readRDS("../2.Annote/2.3.Myeloid/Merged.pca.major.Myeloidcells.pca.rds")
Idents(sub_sce)<-sub_sce$seurat_clusters

##Macrophage (CD68)
##M1 macro (CD38, GPR18, FPR2, HLA-DRA, ITGAX, CD86, NOS2)
##M2 macro (EGR2, MYC, CD163, VEGFA, MAF)
##cDC1 (IDO1, LAMP3)
##cDC2 (CD1C)
##pDC (LILRA4, GPR183)
##monocytes (FCGR3A, CD14, S100A9)
##neutrophils (CDK11B, FCGR3A, and MME)

DimPlot(sub_sce, reduction = "umap", label = TRUE)
pdf("../2.Annote/2.3.Myeloid/replot.markers.pdf",width = 10,height = 18)
FeaturePlot(sub_sce,features = c("IDO1", "LAMP3", "CLEC9A", "XCR1", "CLNK", ##cDC1 DC
                                 "CD1C", "CLEC10A", ##cDC2
                                 "CCR7", "CCL22", ##DC mature
                                 "LILRA4", "GPR183", "IL3RA", "CLEC4C", ##pDC DC cells
                                 "CD68", "APOE", "CD5L", "MARCO", "C1QB", ##Macrophage
                                 "CD14", "VCAN", "FCN1", ##Monocyte
                                 "CSF3R", "FCGR3B" ##Neutrophils
)
,ncol = 3,min.cutoff = 0,raster = T)
dev.off()
FeaturePlot(sub_sce,features = c("MYL9")
            ,min.cutoff = 0,raster = T)
markers <- read.delim("../2.Annote/2.3.Myeloid/sub_sce_Cluster_markers.txt", header = T)
top20genes <- top_n(group_by(markers,cluster), n = 20, wt = avg_log2FC)
panglaomarkers <- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv", header = T)
i = 0
cluster2celltype <- data.frame()
for (i in unique(top20genes$cluster)) {
  newres <- merge(panglaomarkers, top20genes[top20genes$cluster == i, ], by.x = "official.gene.symbol", by.y = "gene")
  newtmp <- as.data.frame(table(newres$cell.type))
  newtmp <- newtmp[order(newtmp$Freq, decreasing = T), ]
  tmpres <- data.frame(Cluster = i,
                       Cells = newtmp[1:5, ]$Var1,
                       COunts = newtmp[1:5, ]$Freq)
  cluster2celltype <- rbind(cluster2celltype, tmpres)
}
View(cluster2celltype)

top10genes <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

new.cluster.id<-c("Macro_Macro_TREM2", "cDC_cDC2_CD1C", "Mono_Mono_VCAN", "Macro_Macro_GPR183",
                  "Macro_Macro_CXCL12", "Macro_Macro_SPP1", "Macro_Macro_CCL3", "UN_UN_UN",
                  "Mono_Mono_CCL20", "Mono_Mono_SPP1", "Mono_Mono_FCN1", "cDC_cDC1_IDO1",
                  "Mono_Mono_FCN1", "pDC_pDC_LILRA4", "Macro_Macro_C1QC", "Mono_Mono_S100A8", 
                  "pDC_pDC_LILRA4")

names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)

sub_sce$MinorCluster <- Idents(sub_sce)
sub_sce$Celltype <- strsplit2(sub_sce$MinorCluster, "_")[, 1]
sub_sce$CellClass <- strsplit2(sub_sce$MinorCluster, "_")[, 2]
sub_sce$Genes <- strsplit2(sub_sce$MinorCluster, "_")[, 3]

sub_sce <- subset(sub_sce, subset = CellClass != "UN")

tmpdf <- data.frame(MinorCluster = sub_sce$MinorCluster,
                    Celltype = sub_sce$Celltype,
                    CellClass = sub_sce$CellClass)
tmpdf <- tmpdf[order(tmpdf$Celltype, tmpdf$CellClass), ]

Idents(sub_sce)<-factor(Idents(sub_sce),levels= as.character(unique(tmpdf$MinorCluster)))                
saveRDS(sub_sce, file = "../2.Annote/2.3.Myeloid/Merged.pca.major.Myeloidcells.pca.minor.rds")
## see annotations
warning("############## see annotations ##############")
p <- DimPlot(sub_sce, reduction = "umap", group.by = "MinorCluster", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.3.Myeloid/minor_Cluster_merge_umap.pdf",p,width = 8,height = 5 )

p <- DimPlot(sub_sce, reduction = "umap", group.by = "Dataset", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.3.Myeloid/minor_Cluster_study_umap.pdf",p,width = 12,height = 9 )

p1 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.3.Myeloid/minor_Cluster_cellclass_umap.pdf",p1,width = 8,height = 7 )


###########2.2.4.Stromal cells###########
allmarkers<- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv",header=T)
gene_list <- allmarkers[allmarkers$cell.type %in% c("Fibroblasts", "Endothelial cells"), ]$official.gene.symbol

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
Idents(combined)<-combined$MajorCluster
sub_sce<-subset(combined,idents=c("Fibro","Endo"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"
dir.create("../2.Annote/2.4.Stromal")
saveRDS(sub_sce, file = "../2.Annote/2.4.Stromal/Merged.pca.major.Stromalcells.rds")

sub_sce <- ScaleData(sub_sce, verbose = FALSE,features = gene_list)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 3000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE,features = gene_list )

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 0.8)
DimPlot(sub_sce, reduction = "umap", label = TRUE)
DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
umap.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "umap")), 
                 Embeddings(object = sub_sce, reduction = "umap"))

write.table(umap.out, file="../2.Annote/2.4.Stromal/Umap_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
tsne.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "tsne")), 
                 Embeddings(object = sub_sce, reduction = "tsne"))
write.table(tsne.out, file="../2.Annote/2.4.Stromal/Tsne_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
cluster.out = cbind("Barcode" = rownames(sub_sce@meta.data), 
                    "Cluster" = sub_sce@meta.data$seurat_clusters)
write.table(cluster.out, file="../2.Annote/2.4.Stromal/Cluster_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
## Save
warning("############## Save data ##############")
saveRDS(sub_sce, file = "../2.Annote/2.4.Stromal/Merged.pca.major.Stromalcells.pca.rds")

## Visualization
warning("############## Visualization ##############")
sub_sce<-readRDS("../2.Annote/2.4.Stromal/Merged.pca.major.Stromalcells.pca.rds")
pdf("../2.Annote/2.4.Stromal/Cluster_merge_umap.pdf", height = 7, width = 18)
p1<-DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
p2 <- DimPlot(sub_sce, reduction = "umap", label = TRUE)
print(plot_grid(p1,p2))
dev.off()
pdf("../2.Annote/2.4.Stromal/Cluster_merge_umap2.pdf", height = 7, width = 18)
p1 <- DimPlot(sub_sce, reduction = "umap", split.by = "MajorCluster",label = TRUE)
print(plot_grid(p1))
dev.off()

sub_sce_down <- subset(x = sub_sce, downsample = 1000)
markers <- FindAllMarkers(sub_sce_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../2.Annote/2.4.Stromal/sub_sce_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)

## Visualization
warning("############## Visualization ##############")
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
pdf("../2.Annote/2.4.Stromal/sub_sce_Cluster_marker_heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
print(DoHeatmap(sub_sce, features = top10$gene,label = F))
dev.off()

pdf ("../2.Annote/2.4.Stromal/Cluter_DotPlot.pdf", height = 8, width = 30)
print(DotPlot(sub_sce, features = unique(c(top10$gene,"ICOS")),
              cols = c("grey","blue"))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete(""))
dev.off()

## minor level clustering
warning("############## Annotation genes ##############")
sub_sce<-readRDS("../2.Annote/2.4.Stromal/Merged.pca.major.Stromalcells.pca.rds")
Idents(sub_sce)<-sub_sce$seurat_clusters

DimPlot(sub_sce, reduction = "umap", label = TRUE)
pdf("../2.Annote/2.4.Stromal/replot.markers.pdf",width = 10,height = 18)
FeaturePlot(sub_sce,features = c("PECAM1", "VWF", "CDH5", "SELE", ##Endothelial  cells
                                 "CCL21", ##Endothelial  cells Lymphatic
                                 "PDGFRA", "FAP", "COL1A1", ##Fibroblast
                                 "ACTG2", ##MyoFibroblast
                                 "MSLN", "CALB2" ##Mesothelial cells,
)
,ncol = 3,min.cutoff = 0,raster = T)
dev.off()

markers <- read.delim("../2.Annote/2.4.Stromal/sub_sce_Cluster_markers.txt", header = T)
top20genes <- top_n(group_by(markers,cluster), n = 20, wt = avg_log2FC)
panglaomarkers <- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv", header = T)
i = 0
cluster2celltype <- data.frame()
for (i in unique(top20genes$cluster)) {
  newres <- merge(panglaomarkers, top20genes[top20genes$cluster == i, ], by.x = "official.gene.symbol", by.y = "gene")
  newtmp <- as.data.frame(table(newres$cell.type))
  newtmp <- newtmp[order(newtmp$Freq, decreasing = T), ]
  tmpres <- data.frame(Cluster = i,
                       Cells = newtmp[1:5, ]$Var1,
                       COunts = newtmp[1:5, ]$Freq)
  cluster2celltype <- rbind(cluster2celltype, tmpres)
}
View(cluster2celltype)

top10genes <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

new.cluster.id<-c("Fibro_Fibro_APOD", "Endo_Endo_PLVAP", "Fibro_Fibro_MYH11", "Fibro_Fibro_RGS5",
                  "Endo_Endo_SELE", "Fibro_Fibro_CST1", "Endo_Endo_CLDN5", "UN_UN_UN",
                  "Endo_Endo_RBP7", "Endo_Endo_FCN3", "Endo_Endo_CCL21") 

names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)

sub_sce$MinorCluster <- Idents(sub_sce)
sub_sce$Celltype <- strsplit2(sub_sce$MinorCluster, "_")[, 1]
sub_sce$CellClass <- strsplit2(sub_sce$MinorCluster, "_")[, 2]
sub_sce$Genes <- strsplit2(sub_sce$MinorCluster, "_")[, 3]

sub_sce <- subset(sub_sce, subset = CellClass != "UN")

tmpdf <- data.frame(MinorCluster = sub_sce$MinorCluster,
                    Celltype = sub_sce$Celltype,
                    CellClass = sub_sce$CellClass)
tmpdf <- tmpdf[order(tmpdf$Celltype, tmpdf$CellClass), ]

Idents(sub_sce)<-factor(Idents(sub_sce),levels= as.character(unique(tmpdf$MinorCluster)))                
saveRDS(sub_sce, file = "../2.Annote/2.4.Stromal/Merged.pca.major.Stromalcells.pca.minor.rds")
## see annotations
warning("############## see annotations ##############")
p <- DimPlot(sub_sce, reduction = "umap", group.by = "MinorCluster", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.4.Stromal/minor_Cluster_merge_umap.pdf",p,width = 8,height = 5 )

p <- DimPlot(sub_sce, reduction = "umap", group.by = "Dataset", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.4.Stromal/minor_Cluster_study_umap.pdf",p,width = 12,height = 9 )

p1 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.4.Stromal/minor_Cluster_cellclass_umap.pdf",p1,width = 8,height = 7 )

###########2.2.5.inferCNV###########
library(infercnv)

dir.create("../2.Annote/2.5.inferCNV")

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
Idents(combined)<-combined$MajorCluster
sub_sce<-subset(combined,idents=c("Epi", "Endo"))

##select cancer and CAF to perform infercnv
genecode <- fread("/home/yukai/projects/sclearn/NPC2023fromMei/0.data/gencode.v25.annotation4infercnv.bed", header = F, data.table = F)
genecode <- genecode %>% distinct(V1, .keep_all = TRUE)
genecode <- genecode[genecode$V1 %in% sub_sce@assays[["RNA"]]@counts@Dimnames[[1]], ]
write.table(genecode, 
            '../2.Annote/2.5.inferCNV/gencode.v25.annotation4infercnv.unique.bed', 
            sep = '\t', row.names = F, col.names = F, quote = F)
tmpmtr <- as_matrix(sub_sce@assays$RNA@counts)
fwrite(as.data.frame(tmpmtr), 
       '../2.Annote/2.5.inferCNV/epi.count.txt', 
       sep = '\t', row.names = T, col.names = T, quote = F)

res <- as.data.frame(sub_sce$MajorCluster)
rownames(res) <- rownames(sub_sce@meta.data)
write.table(res, 
            '../2.Annote/2.5.inferCNV/epi.sampleinfo', 
            sep = '\t', row.names = T, col.names = F, quote = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../2.Annote/2.5.inferCNV/epi.count.txt",
                                    annotations_file="../2.Annote/2.5.inferCNV/epi.sampleinfo",
                                    delim="\t",
                                    gene_order_file="../2.Annote/2.5.inferCNV/gencode.v25.annotation4infercnv.unique.bed",
                                    ref_group_names=c("Endo")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="../2.Annote/2.5.inferCNV/", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             resume_mode=FALSE)

##deal results
infercnv_obj = readRDS("../2.Annote/2.5.inferCNV/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$Endo
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- c(test_loc$`Hepatic progenitor`, test_loc$Epi, test_loc$Hepatocytes,
              test_loc$Cholangiocytes, test_loc$`Hepa prolif`)

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)

gn <- rownames(expr)
geneFile <- read.table("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/2.Annote/2.5.inferCNV/gencode.v25.annotation4infercnv.unique.bed",
                       header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]

#cluster, extract results
set.seed(20210418)
kmeans.result <- kmeans(t(expr), 10)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #merge
kmeans_df_s=arrange(kmeans_df,kmeans_class) #sort
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #transfer to factor
head(kmeans_df_s)

#annotation heatmap and color
library(ComplexHeatmap)
library(circlize)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:length(unique(kmeans_df_s$kmeans_class))]
names(color_v)=as.character(1:length(unique(kmeans_df_s$kmeans_class)))
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))

#heatmap
ranselect <- sort(sample(1:ncol(expr), 4000))
pdf("../2.Annote/2.5.inferCNV/infercnv.res.cluster.heatmap.pdf",width = 15,height = 10)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),][ranselect, ], #绘图数据的CB顺序和注释CB顺序保持一致
             col = colorRamp2(c(0.8,1,1.2), c("#377EB8","#F0F0F0","#E41A1C")), #10x(0.8, 1, 1.2)   smartseq(0.4, 1, 1.6)
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), #sort chr
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno[ranselect, ], #add annotation
             row_title = NULL,column_title = NULL,
             use_raster = T, raster_device = "CairoPNG")
draw(ht, heatmap_legend_side = "right")
dev.off()

#save CB for each class
fwrite(kmeans_df_s, file = "../2.Annote/2.5.inferCNV/infercnv.res.cluster.kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)

##visualize cluster CNV score
expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
kmeans_df_s$CB=rownames(kmeans_df_s)
CNV_score=CNV_score%>%inner_join(kmeans_df_s,by="CB")

p <- CNV_score%>%ggplot(aes(kmeans_class,CNV_score))+geom_violin(aes(fill=kmeans_class),color="NA")+
  scale_fill_manual(values = color_v)+
  theme_bw()
p
ggsave("../2.Annote/2.5.inferCNV/infercnv.res.cluster.CNV level.pdf",p,width = 10,height = 6,units = "cm")
fwrite(CNV_score, file = "../2.Annote/2.5.inferCNV/infercnv.res.cluster.CNV.score.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)


###########2.2.6.Epi cells###########
allmarkers<- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv",header=T)
gene_list <- allmarkers[allmarkers$cell.type %in% c("Epithelial cells"), ]$official.gene.symbol

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
Idents(combined)<-combined$MajorCluster
sub_sce<-subset(combined,idents=c("Epi"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"
dir.create("../2.Annote/2.6.Epi")
saveRDS(sub_sce, file = "../2.Annote/2.6.Epi/Merged.pca.major.Epicells.rds")

sub_sce <- ScaleData(sub_sce, verbose = FALSE,features = gene_list)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 3000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE,features = gene_list )

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 0.6)
DimPlot(sub_sce, reduction = "umap", label = TRUE)
DimPlot(sub_sce, reduction = "umap", group.by = "Tissue", label = TRUE)
DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
umap.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "umap")), 
                 Embeddings(object = sub_sce, reduction = "umap"))

write.table(umap.out, file="../2.Annote/2.6.Epi/Umap_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
tsne.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "tsne")), 
                 Embeddings(object = sub_sce, reduction = "tsne"))
write.table(tsne.out, file="../2.Annote/2.6.Epi/Tsne_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
cluster.out = cbind("Barcode" = rownames(sub_sce@meta.data), 
                    "Cluster" = sub_sce@meta.data$seurat_clusters)
write.table(cluster.out, file="../2.Annote/2.6.Epi/Cluster_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
## Save
warning("############## Save data ##############")
saveRDS(sub_sce, file = "../2.Annote/2.6.Epi/Merged.pca.major.Epicells.pca.rds")

## Visualization
warning("############## Visualization ##############")
sub_sce<-readRDS("../2.Annote/2.6.Epi/Merged.pca.major.Epicells.pca.rds")
pdf("../2.Annote/2.6.Epi/Cluster_merge_umap.pdf", height = 7, width = 18)
p1<-DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
p2 <- DimPlot(sub_sce, reduction = "umap", label = TRUE)
print(plot_grid(p1,p2))
dev.off()
pdf("../2.Annote/2.6.Epi/Cluster_merge_umap2.pdf", height = 7, width = 18)
p1 <- DimPlot(sub_sce, reduction = "umap", split.by = "MajorCluster",label = TRUE)
print(plot_grid(p1))
dev.off()

sub_sce_down <- subset(x = sub_sce, downsample = 1000)
markers <- FindAllMarkers(sub_sce_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../2.Annote/2.6.Epi/sub_sce_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)

## Visualization
warning("############## Visualization ##############")
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
pdf("../2.Annote/2.6.Epi/sub_sce_Cluster_marker_heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
print(DoHeatmap(sub_sce, features = top10$gene,label = F))
dev.off()

pdf ("../2.Annote/2.6.Epi/Cluter_DotPlot.pdf", height = 8, width = 30)
print(DotPlot(sub_sce, features = unique(c(top10$gene,"ICOS")),
              cols = c("grey","blue"))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete(""))
dev.off()

## minor level clustering
warning("############## Annotation genes ##############")
sub_sce<-readRDS("../2.Annote/2.6.Epi/Merged.pca.major.Epicells.pca.rds")
Idents(sub_sce)<-sub_sce$seurat_clusters

##CNV identify
CNVscore <- fread("../2.Annote/2.5.inferCNV/infercnv.res.cluster.CNV.score.txt", header = T, data.table = F)
cellinfo <- data.frame(ID = rownames(sub_sce@meta.data),
                       seurat_clusters = sub_sce$seurat_clusters)
cellinfo <- merge(cellinfo, CNVscore, by.x = "ID", by.y = "CB", all.y = T)
cellinfo$seurat_clusters <- ifelse(cellinfo$class == "normal", "42", cellinfo$seurat_clusters)
cellinfo$CNV_score <- ifelse(cellinfo$CNV_score > 0.01, 0.01, cellinfo$CNV_score)
p <- cellinfo%>%ggplot(aes(seurat_clusters,CNV_score))+geom_violin(aes(fill=seurat_clusters),color="NA")+
  scale_fill_manual(values = color_v)+
  theme_bw()
p
ggsave("../2.Annote/2.6.Epi/infercnv.cluster.CNVscore.pdf",p,width = 20,height = 6,units = "cm")

##annote
DimPlot(sub_sce, reduction = "umap", label = TRUE)
pdf("../2.Annote/2.6.Epi/replot.markers.pdf",width = 10,height = 8)
FeaturePlot(sub_sce,features = c("KRT8", "EPCAM", "CDH1", "CLDN6", "CLDN7", ##Epithial cells,
                                 "KDM4B", "KLK15", "PLIN2", "GLCCI1", "PGA3"
),ncol = 3,min.cutoff = 0,raster = T)
dev.off()

markers <- read.delim("../2.Annote/2.6.Epi/sub_sce_Cluster_markers.txt", header = T)
top20genes <- top_n(group_by(markers,cluster), n = 20, wt = avg_log2FC)
panglaomarkers <- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv", header = T)
i = 0
cluster2celltype <- data.frame()
for (i in unique(top20genes$cluster)) {
  newres <- merge(panglaomarkers, top20genes[top20genes$cluster == i, ], by.x = "official.gene.symbol", by.y = "gene")
  newtmp <- as.data.frame(table(newres$cell.type))
  newtmp <- newtmp[order(newtmp$Freq, decreasing = T), ]
  tmpres <- data.frame(Cluster = i,
                       Cells = newtmp[1:5, ]$Var1,
                       COunts = newtmp[1:5, ]$Freq)
  cluster2celltype <- rbind(cluster2celltype, tmpres)
}
View(cluster2celltype)

top10genes <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

new.cluster.id<-c("Epi_Chol_CEACAM6", "Epi_Hepa_APOA1", "Epi_Malig_DPEP1", "Epi_Malig_ICAM1", 
                  "Epi_Hepa_APOA1", "UN_UN_UN", "Epi_Malig_MUC1", "Epi_Hepa_PRG4", 
                  "UN_UN_UN", "Epi_Malig_EPCAM", "Epi_Hepa_ITGA1", "Epi_Malig_IL1R1", 
                  "Epi_Hepa_ITGA5", "Epi_Malig_NRP2", "Epi_Malig_KLK10", "UN_UN_UN", 
                  "Epi_Malig_CCND1", "UN_UN_UN", "Epi_Hepa_IGFBP2", "Epi_Malig_CXCL3", 
                  "Epi_Hepa_MT1H") 
DimPlot(sub_sce, reduction = "umap", group.by = "Tissue", label = TRUE)

names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)

sub_sce$MinorCluster <- Idents(sub_sce)
sub_sce$Celltype <- strsplit2(sub_sce$MinorCluster, "_")[, 1]
sub_sce$CellClass <- strsplit2(sub_sce$MinorCluster, "_")[, 2]
sub_sce$Genes <- strsplit2(sub_sce$MinorCluster, "_")[, 3]

sub_sce <- subset(sub_sce, subset = CellClass != "UN")

tmpdf <- data.frame(MinorCluster = sub_sce$MinorCluster,
                    Celltype = sub_sce$Celltype,
                    CellClass = sub_sce$CellClass)
tmpdf <- tmpdf[order(tmpdf$Celltype, tmpdf$CellClass), ]

Idents(sub_sce)<-factor(Idents(sub_sce),levels= as.character(unique(tmpdf$MinorCluster)))                
saveRDS(sub_sce, file = "../2.Annote/2.6.Epi/Merged.pca.major.Epicells.pca.minor.rds")
## see annotations
warning("############## see annotations ##############")
p <- DimPlot(sub_sce, reduction = "umap", group.by = "MinorCluster", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.6.Epi/minor_Cluster_merge_umap.pdf",p,width = 8,height = 5 )

p <- DimPlot(sub_sce, reduction = "umap", group.by = "Dataset", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.6.Epi/minor_Cluster_study_umap.pdf",p,width = 12,height = 9 )

p1 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.6.Epi/minor_Cluster_cellclass_umap.pdf",p1,width = 8,height = 7 )


###########2.2.7.Mast cells###########
gene_list<- read.table("/home/yukai/projects/sclearn/NPC2023fromMei/0.data/Myeloid_markers.txt",header=F)
gene_list <- unique(c(gene_list$V1))

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
Idents(combined)<-combined$MajorCluster
sub_sce<-subset(combined,idents=c("Mast"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"
dir.create("../2.Annote/2.7.Mast")
saveRDS(sub_sce, file = "../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.rds")

sub_sce <- ScaleData(sub_sce, verbose = FALSE,features = gene_list)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 3000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE,features = gene_list)

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 0.8)
DimPlot(sub_sce, reduction = "umap", label = TRUE)
DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
umap.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "umap")), 
                 Embeddings(object = sub_sce, reduction = "umap"))

write.table(umap.out, file="../2.Annote/2.7.Mast/Umap_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
tsne.out = cbind("Barcode" = rownames(Embeddings(object = sub_sce, reduction = "tsne")), 
                 Embeddings(object = sub_sce, reduction = "tsne"))
write.table(tsne.out, file="../2.Annote/2.7.Mast/Tsne_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
cluster.out = cbind("Barcode" = rownames(sub_sce@meta.data), 
                    "Cluster" = sub_sce@meta.data$seurat_clusters)
write.table(cluster.out, file="../2.Annote/2.7.Mast/Cluster_out.csv", sep = "\t", quote = F, row.names = F, col.names = T)
## Save
warning("############## Save data ##############")
saveRDS(sub_sce, file = "../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.pca.rds")

## Visualization
warning("############## Visualization ##############")
sub_sce<-readRDS("../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.pca.rds")
pdf("../2.Annote/2.7.Mast/Cluster_merge_umap.pdf", height = 7, width = 18)
p1<-DimPlot(sub_sce, reduction = "umap", group.by = "MajorCluster", label = TRUE)
p2 <- DimPlot(sub_sce, reduction = "umap", label = TRUE)
print(plot_grid(p1,p2))
dev.off()
pdf("../2.Annote/2.7.Mast/Cluster_merge_umap2.pdf", height = 7, width = 18)
p1 <- DimPlot(sub_sce, reduction = "umap", split.by = "MajorCluster",label = TRUE)
print(plot_grid(p1))
dev.off()

sub_sce_down <- subset(x = sub_sce, downsample = 1000)
markers <- FindAllMarkers(sub_sce_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../2.Annote/2.7.Mast/sub_sce_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)

## Visualization
warning("############## Visualization ##############")
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
pdf("../2.Annote/2.7.Mast/sub_sce_Cluster_marker_heatmap.pdf", height = 20, width = 20)
#png("2.Intergration/Cluter_marker_heatmap.png", height = 4000, width = 5000)
print(DoHeatmap(sub_sce, features = top10$gene,label = F))
dev.off()

pdf ("../2.Annote/2.7.Mast/Cluter_DotPlot.pdf", height = 8, width = 30)
print(DotPlot(sub_sce, features = unique(c(top10$gene,"ICOS")),
              cols = c("grey","blue"))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete(""))
dev.off()

## minor level clustering
warning("############## Annotation genes ##############")
sub_sce<-readRDS("../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.pca.rds")
Idents(sub_sce)<-sub_sce$seurat_clusters

##Macrophage (CD68)

DimPlot(sub_sce, reduction = "umap", label = TRUE)
pdf("../2.Annote/2.7.Mast/replot.markers.pdf",width = 10,height = 18)
FeaturePlot(sub_sce,features = c("IDO1", "LAMP3", "CLEC9A", "XCR1", "CLNK", ##cDC1 DC
                                 "CD1C", "CLEC10A", ##cDC2
                                 "CCR7", "CCL22", ##DC mature
                                 "LILRA4", "GPR183", "IL3RA", "CLEC4C", ##pDC DC cells
                                 "CD68", "APOE", "CD5L", "MARCO", "C1QB", ##Macrophage
                                 "CD14", "VCAN", "FCN1", ##Monocyte
                                 "CSF3R", "FCGR3B" ##Neutrophils
)
,ncol = 3,min.cutoff = 0,raster = T)
dev.off()
FeaturePlot(sub_sce,features = c("MYL9")
            ,min.cutoff = 0,raster = T)
markers <- read.delim("../2.Annote/2.7.Mast/sub_sce_Cluster_markers.txt", header = T)
top20genes <- top_n(group_by(markers,cluster), n = 20, wt = avg_log2FC)
panglaomarkers <- read.delim("/home/yukai/work/projects/XXLYDDAQ/wangruiqi/LiverTME/0.data/PanglaoDB_markers_27_Mar_2020.tsv", header = T)
i = 0
cluster2celltype <- data.frame()
for (i in unique(top20genes$cluster)) {
  newres <- merge(panglaomarkers, top20genes[top20genes$cluster == i, ], by.x = "official.gene.symbol", by.y = "gene")
  newtmp <- as.data.frame(table(newres$cell.type))
  newtmp <- newtmp[order(newtmp$Freq, decreasing = T), ]
  tmpres <- data.frame(Cluster = i,
                       Cells = newtmp[1:5, ]$Var1,
                       COunts = newtmp[1:5, ]$Freq)
  cluster2celltype <- rbind(cluster2celltype, tmpres)
}
View(cluster2celltype)

top10genes <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

new.cluster.id<-c("Mast_Mast_KRT86", "Mast_Mast_TPSAB1", "Mast_Mast_CPA3", "Mast_Mast_S100A8",
                  "Mast_Mast_HSPA1A", "Mast_Mast_ITGA2B", "Mast_Mast_IL32", "UN_UN_UN")

names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)

sub_sce$MinorCluster <- Idents(sub_sce)
sub_sce$Celltype <- strsplit2(sub_sce$MinorCluster, "_")[, 1]
sub_sce$CellClass <- strsplit2(sub_sce$MinorCluster, "_")[, 2]
sub_sce$Genes <- strsplit2(sub_sce$MinorCluster, "_")[, 3]

sub_sce <- subset(sub_sce, subset = CellClass != "UN")

tmpdf <- data.frame(MinorCluster = sub_sce$MinorCluster,
                    Celltype = sub_sce$Celltype,
                    CellClass = sub_sce$CellClass)
tmpdf <- tmpdf[order(tmpdf$Celltype, tmpdf$CellClass), ]

Idents(sub_sce)<-factor(Idents(sub_sce),levels= as.character(unique(tmpdf$MinorCluster)))                
saveRDS(sub_sce, file = "../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.pca.minor.rds")
## see annotations
warning("############## see annotations ##############")
p <- DimPlot(sub_sce, reduction = "umap", group.by = "MinorCluster", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.7.Mast/minor_Cluster_merge_umap.pdf",p,width = 8,height = 5 )

p <- DimPlot(sub_sce, reduction = "umap", group.by = "Dataset", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.7.Mast/minor_Cluster_study_umap.pdf",p,width = 12,height = 9 )

p1 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass", label = TRUE, raster=TRUE)
ggsave("../2.Annote/2.7.Mast/minor_Cluster_cellclass_umap.pdf",p1,width = 8,height = 7 )




####################2.3.MergeAnnotations####################
library(Seurat)
library(data.table)
library(sscVis)
library(ggplot2)
library(ggridges)
library(dplyr)
library(cowplot)
library(PCAtools)

combined<-readRDS("../2.Annote/Merged.pca.major.rds")
combined$MinorCluster<-"Unknown"

name_list<-c("../2.Annote/2.1.Tcells/Merged.pca.major.Tcells.pca.minor.rds",
             "../2.Annote/2.2.Bcells/Merged.pca.major.Bcells.pca.minor.rds",
             "../2.Annote/2.3.Myeloid/Merged.pca.major.Myeloidcells.pca.minor.rds",
             "../2.Annote/2.4.Stromal/Merged.pca.major.Stromalcells.pca.minor.rds",
             "../2.Annote/2.6.Epi/Merged.pca.major.Epicells.pca.minor.rds",
             "../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.pca.minor.rds")

for (file_name in name_list){
  sub_sce<-readRDS(file_name)
  ann_table<-data.frame(Cell_ID=rownames(sub_sce@meta.data),
                        MinorCluster=sub_sce$MinorCluster)
  Cluster_list<-unique(ann_table$MinorCluster)
  for (Cluster in Cluster_list){
    Cell_ID_list<-ann_table$Cell_ID[which(ann_table$MinorCluster==Cluster)]
    combined$MinorCluster[which(rownames(combined@meta.data) %in%Cell_ID_list)]<-Cluster
  }
}

combined$Celltype <- strsplit2(combined$MinorCluster, "_")[, 1]
combined$CellClass <- strsplit2(combined$MinorCluster, "_")[, 2]
combined$Genes <- strsplit2(combined$MinorCluster, "_")[, 3]
DimPlot(combined, reduction = "umap", group.by = "Celltype", raster=FALSE)

##Subcluster
combined <- subset(combined, subset = Celltype != "Unknown")
DimPlot(combined, reduction = "umap", group.by = "Celltype", raster=FALSE)

combined <- ScaleData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 3000)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
ElbowPlot(combined)
chosen.elbow <- findElbowPoint(Stdev(object = combined, reduction = "pca"))
combined <-RunHarmony(combined,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:15)
combined <- RunTSNE(combined, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:15)
combined <- FindClusters(combined, resolution = 0.6)
DimPlot(combined, reduction = "umap", group.by = "Celltype", raster=FALSE)

## Visualization
dir.create("../2.Annote/MergeAnno")
saveRDS(combined,"../2.Annote/MergeAnno/Merged.pca.major.minor.rds") 

warning("############## Visualization ##############")
pbmc<-combined
#p0 <- DimPlot(pbmc, reduction = "umap", group.by = "orig.ident", raster=FALSE)
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "Celltype", raster=FALSE)
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "MinorCluster", raster=FALSE)
p3 <- DimPlot(pbmc, reduction = "umap", group.by = "Dataset", raster=FALSE)
pdf("../2.Annote/MergeAnno/Merged.features.umap.pdf", width = 24, height = 6)
p1|p2|p3
dev.off()


#newrm <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.rds")
newrm <- combined
newrm$TopCluster <- as.character(newrm$Celltype)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("CD4T"), "CD4+ T cell", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("CD8T"), "CD8+ T cell", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("NK"), "NK cell", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("B"), "B cell", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("Endo"), "Endothelial cell", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("Fibro"), "Fibroblast", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("Epi"), "Epithelial cell", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("Macro", "Mono"), "Macrophage/Monocyte", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("Plasma"), "Plasma cell", newrm$TopCluster)
newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("Mast"), "Mast cell", newrm$TopCluster)
#newrm$TopCluster <- ifelse(newrm$TopCluster %in% c("Neu"), "Neutrophils", newrm$TopCluster)

#newrm$TopCluster <- ifelse(newrm$CellClass %in% c("Tprolif", "Bprolif", "NKprolif"), "Prolif cell", newrm$TopCluster)


newrm$TopCluster<-factor(newrm$TopCluster,levels=unique(sort(newrm$TopCluster)))
DimPlot(newrm, reduction = "umap", group.by = "TopCluster", raster=FALSE)
DimPlot(newrm, reduction = "umap", group.by = "seurat_clusters", raster=FALSE, label = T)

#newrm <- ScaleData(newrm, verbose = FALSE)
#newrm <- FindVariableFeatures(newrm, selection.method = "vst", nfeatures = 3000)
#newrm <- RunPCA(newrm, npcs = 30, verbose = FALSE)
#ElbowPlot(newrm)
#chosen.elbow <- findElbowPoint(Stdev(object = newrm, reduction = "harmony"))
#newrm <-RunHarmony(newrm,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)
#newrm <- RunUMAP(newrm, reduction = "harmony", dims = 1:15)
#newrm <- RunTSNE(newrm, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
#newrm <- FindNeighbors(newrm, reduction = "harmony", dims = 1:15)
#DimPlot(newrm, reduction = "umap", group.by = "TopCluster", raster=FALSE)
saveRDS(newrm,"../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.rds") 

Idents(newrm) <- newrm$TopCluster
newrm_down <- subset(x = newrm, downsample = 1000)
markers <- FindAllMarkers(newrm_down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../2.Annote/MergeAnno/Combined_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)
top10 <-  top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

pdf ("../2.Annote/MergeAnno/Cluter_DotPlot.pdf", height = 20, width = 40)
DotPlot(combined, features = unique(top10$gene),
        cols = c("grey","blue"))+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
dev.off()

###########################################3.Analysis###########################################
library(Seurat)
library(dplyr)
library(SingleCellExperiment)
library(sscVis)
library(ggplot2)
library(cowplot)

dir.create("../3.Analysis")
pbmc <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.rds")
pbmc$CompClass <- as.character(pbmc$TopCluster)
pbmc$CompClass <- ifelse(pbmc$Celltype == "CD4T", "T cell CD4", pbmc$CompClass)
pbmc$CompClass <- ifelse(pbmc$Celltype == "CD8T", "T cell CD8", pbmc$CompClass)

####################3.1.Overview####################
###########3.1.1.all cells###########
p3 <- DimPlot(pbmc, reduction = "umap", group.by = "Dataset", raster=TRUE)
p3
pdf("../3.Analysis/3.1.overview.allcell.dataset.pdf", width = 4.5, height = 3)
p3
dev.off()


my_cols <- c('CD4+ T cell'='#F68282', 'CD8+ T cell'='#F78247','NK cell'='#ff9a36','Plasma cell'='#31C53F','B cell'='#1FA195',
             'Macrophage/Monocyte'='#B95FBB','cDC'='#D4D915','pDC'='#2FF18B',
             'Mast cell'='#25aff5','Epithelial cell'='#4B4BF7','Endothelial cell'='#AC8F14',
             'Fibroblast'='#E6C122')
my_cols2 <- my_cols[order(names(my_cols))]
p2 <- DimPlot(pbmc, reduction = "umap", group.by = "TopCluster", cols = my_cols2, raster=TRUE, label = TRUE) + NoLegend()
p2
pdf("../3.Analysis/3.1.overview.allcell.cellclass.pdf", width = 5, height = 5)
p2
dev.off()

###########3.1.2.immune cells###########
Idents(pbmc)<-pbmc$TopCluster
sub_sce<-subset(pbmc,idents=c("B cell", "cDC", "Macrophage/Monocyte", "Mast cell", 
                              "NK cell", "pDC", "Plasma cell", "CD4+ T cell", "CD8+ T cell"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"

sub_sce <- ScaleData(sub_sce, verbose = FALSE)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 1000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE)

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 1.5)
DimPlot(sub_sce, reduction = "umap", label = TRUE)

my_cols <- c('CD4+ T cell'='#F68282', 'CD8+ T cell'='#F78247','NK cell'='#ff9a36','Plasma cell'='#31C53F','B cell'='#1FA195',
             'Macrophage/Monocyte'='#B95FBB','cDC'='#D4D915','pDC'='#2FF18B',
             'Mast cell'='#25aff5','Epithelial cell'='#4B4BF7','Endothelial cell'='#AC8F14',
             'Fibroblast'='#E6C122')

p2 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass", raster=TRUE, label = TRUE) + NoLegend()
p2
pdf("../3.Analysis/3.1.overview.immunecells.cellclass.pdf", width = 5, height = 5)
p2
dev.off()


saveRDS(sub_sce, file = "../3.Analysis/Merged.pca.major.immunecells.rds")


###########3.1.3.epithelial cells###########
sub_sce<-readRDS("../2.Annote/2.6.Epi/Merged.pca.major.Epicells.pca.minor.rds")
#combined<-readRDS("2.intergrate/Combined.rds")
DimPlot(sub_sce, reduction = "umap", label = TRUE)

my_cols <- c('Hepa'='#4B4BF7', 'Chol'='#2FF18B', 'Malig'='#AC8F14')

p2 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass",cols =  my_cols, raster=TRUE, label = TRUE) + NoLegend()
p2
pdf("../3.Analysis/3.1.overview.epicells.cellclass.pdf", width = 5, height = 5)
p2
dev.off()

saveRDS(sub_sce, file = "../3.Analysis/Merged.pca.major.epicells.rds")


###########3.1.4.stromal cells###########
Idents(pbmc)<-pbmc$TopCluster
sub_sce<-subset(pbmc,idents=c("Endothelial cell", "Fibroblast"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"

sub_sce <- ScaleData(sub_sce, verbose = FALSE)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 1000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE)

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 1.5)
DimPlot(sub_sce, reduction = "umap", label = TRUE)

my_cols <- c('Endo'='#AC8F14',
             'Fibro'='#E6C122')

p2 <- DimPlot(sub_sce, reduction = "umap", group.by = "CellClass",cols =  my_cols, raster=TRUE, label = TRUE) + NoLegend()
p2
pdf("../3.Analysis/3.1.overview.stromalcells.cellclass.pdf", width = 5, height = 5)
p2
dev.off()

saveRDS(sub_sce, file = "../3.Analysis/Merged.pca.major.stromalcells.rds")


###########3.1.5.difference in cell comp###########
library(data.table)
library(Seurat)
library(dplyr)
library(plyr)
meta.tb <- pbmc@meta.data
out.prefix <- "../3.Analysis/3.1.cellcount"

cellInfo.tb <- meta.tb[meta.tb$Tissue %in% c("PL", "LM", "PC", "CRC"), ]
OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$CompClass,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.TopCluster.pl.lm.pc.crc",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.1.cellcount.STARTRAC.dist.TopCluster.pl.lm.pc.crc.OR.dist.txt", 
       row.names = T, quote = F)

cellInfo.tb <- meta.tb[meta.tb$Tissue %in% c("PL", "NL"), ]
OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$CompClass,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.TopCluster.pl.nl",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.1.cellcount.STARTRAC.dist.TopCluster.pl.nl.OR.dist.txt", 
       row.names = T, quote = F)

###########3.1.6.featureplot of marker genes###########
pdf("../3.Analysis/3.1.MajorCluster.markers.pdf",width = 16,height = 8)
FeaturePlot(pbmc,features = c("MS4A1", ##B cells
                              "MZB1", ##Plasma cells
                              "CD3D", ## T  cells
                              "CD4", ## CD4+ T cells
                              "CD8A", ## CD8+ T cells
                              "GNLY", ## NK
                              "LAMP3", ##cDC1 DC
                              "CD1C", ##cDC2
                              "LILRA4", ##pDC DC cells
                              "CD68", ##Macrophage
                              "MS4A2", ##Mast cells,
                              "FCN1", ##Monocyte
                              "VWF", ##Endothelial  cells
                              "COL1A1", ##Fibroblast
                              "KRT8" ##Epithial cells,
),ncol = 5,min.cutoff = 0,raster = T)
dev.off()

###########3.1.7.featureplot of gastric genes###########
#many gastric markers (ANXA10, VSIG1, CLDN18, CTSE, TFF2, MUC5AC and MUC6) 
#PGC, GIF, GAST, and ATP4A
subsce <- readRDS("../3.Analysis/Merged.pca.major.epicells.rds")
Idents(subsce) <- subsce$CellClass
subsce <- subset(subsce, subset = CellClass %in% c("Hepa", "Malig"))

tmpdf <- data.frame(ID = subsce$ID,
                    CellClass = subsce$CellClass,
                    Exp = FetchData(subsce, vars = "ALB")$ALB)
p <- tmpdf%>%ggplot(aes(CellClass,Exp))+geom_violin(aes(fill=CellClass),color="NA")+
  scale_fill_manual(values = c("#377EB8","#E41A1C"))+
  theme_bw()
p
ggsave("../3.Analysis/3.1.MajorCluster.gastric.markers.pdf",p,width = 4,height = 3)

###########3.1.8.cell distribution###########
res <- as.data.frame(prop.table(table(pbmc$CompClass)))
res <- res[order(res$Freq, decreasing = F), ]
res$Var1 <- factor(res$Var1, res$Var1)

p1 <- ggplot(res, aes(x = Freq, y = Var1, fill = Var1))+
  geom_bar(stat="identity", position=position_dodge())+
  xlab("fraction")+
  #stat_compare_means(label = "p.signif")+
  scale_x_continuous(expand = c(0, 0))+
  ylab("")+
  theme_classic2()+
  theme(legend.position = "none")

p1
ggsave("../3.Analysis/3.1.cell.distribution.pdf", p1, width = 5, height = 3)

table(as.factor(as.character(combined$orig.ident)))
prop.table(table(Idents(pbmc)))
cellratio <- prop.table(table(pbmc$CompClass, as.factor(as.character(pbmc$Dataset))), margin = 2)

cellratio <- as.data.frame(cellratio)

cellratio$Var1 <- factor(cellratio$Var1, res$Var1)
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  coord_flip()+
  #scale_fill_manual(values = cell_type_cols)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(legend.position = "none")
p
pdf("../3.Analysis/3.1.dataset.cellcount.pdf", width = 4, height = 3)
p
dev.off()

###########3.1.9.cnvscore###########
CNV_score <- fread("../2.Annote/2.5.inferCNV/infercnv.res.cluster.CNV.score.txt",
                   header = T, stringsAsFactors = F, data.table = F)
saminfo <- data.frame(ID = pbmc$ID,
                      Group = pbmc$CellClass)
saminfo <- saminfo[saminfo$Group %in% c("Hepa", "Malig", "Chol"), ]
CNV_score <- merge(CNV_score, saminfo, by.x = "CB", by.y = "ID")

CNV_score$CNV_score <- ifelse(CNV_score$CNV_score > 0.015, 0.015, CNV_score$CNV_score)

p <- CNV_score%>%ggplot(aes(Group,CNV_score))+geom_violin(aes(fill=Group),color="NA")+
  scale_fill_manual(values = c("#377EB8","#E41A1C"))+
  theme_bw()
p
ggsave("../3.Analysis/3.1.infercnv.epi.CNV.pdf",p,width = 4,height = 3)


###########3.1.10.dataset.summary###########
cellratio <- as.data.frame(table(pbmc@meta.data[, c("Dataset", "Patient", "CompClass")]))
cellratio <- cellratio[cellratio$Freq > 0, ]

fwrite(cellratio, file = "../3.Analysis/3.1.dataset.summary.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)


####################3.2.TME subtype####################
###########3.2.1.subtypes###########
pbmc.tumor <- subset(pbmc, subset = Tissue %in% c("LM", "PL"))
cellratio <- prop.table(table(pbmc.tumor$CompClass, as.factor(as.character(pbmc.tumor$Patient))), margin = 2)
cellratio <- as.data.frame(cellratio)
cellratio <- cellratio[!(cellratio$Var1 %in% c("Endothelial cell", "Fibroblast")), ]
res <- acast(cellratio,  Var1~Var2, value.var="Freq")

res <- t(scale(t(res)))
pheatmap(as.matrix(res),
         scale = "row",
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = T,
         show_colnames = T,
         color =colorRampPalette(c("blue", "white","red"))(50),
         fontsize = 10
)

resdf <- as.data.frame(t(res))
get_max_column <- function(row) {
  max_value <- max(row[-1])  # 忽略第一列（ID列）
  max_column <- names(row)[-1][row[-1] == max_value]
  return(max_column)
}
resdf$Max_Column <- apply(resdf, 1, get_max_column)
resdf$Class <- NA
resdf$Class <- ifelse(resdf$Max_Column %in% c("Macrophage/Monocyte", "pDC"), "M", resdf$Class)
resdf$Class <- ifelse(resdf$Max_Column %in% c("Epithelial cell"), "D", resdf$Class)
resdf$Class <- ifelse(resdf$Max_Column %in% c("T cell CD4", "T cell CD8", "NK cell"), "T", resdf$Class)
resdf$Class <- ifelse(resdf$Max_Column %in% c("Mast cell", "cDC", "Plasma cell", "B cell"), "B", resdf$Class)

classdf <- data.frame(ID = rownames(resdf),
                      Class = resdf$Class)

classdf$Class <- factor(classdf$Class, c("B", "T", "M", "D"))
classdf <- classdf[order(classdf$Class, decreasing = F), ]

cells <- as.data.frame(prop.table(table(pbmc$CompClass)))
cells <- cells[!(cells$Var1 %in% c("Endothelial cell", "Fibroblast")), ]
cells <- cells[order(cells$Freq, decreasing = F), ]
res <- as.data.frame(res)[as.character(cells$Var1), classdf$ID]

#prepare clinical info
cliinfo <- pbmc.tumor@meta.data[, c("Dataset", "Tissue", "Patient")]
cliinfo <- unique(cliinfo)
rownames(cliinfo) <- cliinfo$Patient
cliinfo <- cliinfo[classdf$ID, ]

#plot
library(ComplexHeatmap)
ha = HeatmapAnnotation(
  Dataset = cliinfo$Dataset,
  Tissue = cliinfo$Tissue,
  Subtype = classdf$Class,
  col = list(Subtype = c(B = "#FCE38AFF", T = "#4DBBD5FF", M = "#E64B35FF", D = "#D4B499FF")
  ),
  na_col = "gray80"
)

newres <- scale(t(scale(t(res))))
Heatmap(newres, cluster_columns = F, top_annotation = ha, 
        cluster_rows = F, show_column_names = F)

pdf("../3.Analysis/3.2.immune.subtype.heatmap.pdf", width = 8, height = 3)
Heatmap(newres, cluster_columns = F, top_annotation = ha, 
        cluster_rows = F, show_column_names = F)
dev.off()


#save pbmc.tumor with immune subtype info
immunesub <- data.frame(ID = rownames(pbmc.tumor@meta.data),
                        Patient = pbmc.tumor$Patient)
rownames(immunesub) <- immunesub$ID
immunesub <- merge(immunesub, classdf, by.x = "Patient", by.y = "ID", all.x = T)
rownames(immunesub) <- immunesub$ID

pbmc.tumor <- AddMetaData(pbmc.tumor, immunesub)
saveRDS(pbmc.tumor, file = "../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.rds")


###########3.2.2.pathways###########
pbmc.tumor <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.rds")
Idents(pbmc.tumor)<-pbmc.tumor$Class
pbmc.tumor.down <- subset(x = pbmc.tumor, downsample = 2000)

pbmc.tumor.mtr <- pbmc.tumor.down@assays$RNA@counts

library(clusterProfiler)
library(msigdbr)
library(dplyr)
hall <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  dplyr::select(gs_description, gene_symbol)
tfs <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>%
  dplyr::select(gs_name, gene_symbol)

hallgmt <- split(hall$gene_symbol, hall$gs_name)
kegggmt <- split(kegg$gene_symbol, kegg$gs_description)
tfsgmt <- split(tfs$gene_symbol, tfs$gs_name)

res_hall <- gsva(as.matrix(pbmc.tumor.mtr), hallgmt, kcdf="Gaussian",method = "ssgsea",parallel.sz=10)
res_kegg <- gsva(as.matrix(pbmc.tumor.mtr), kegggmt, kcdf="Gaussian",method = "ssgsea",parallel.sz=10)
res_tfs <- gsva(as.matrix(pbmc.tumor.mtr), tfsgmt, kcdf="Gaussian",method = "ssgsea",parallel.sz=10)

dim(res_hall)
class(res_hall)

pbmc.tumor.down.hall<-CreateSeuratObject(res_hall)
Idents(pbmc.tumor.down.hall)<-pbmc.tumor.down$Class

pbmc.tumor.down.kegg<-CreateSeuratObject(res_kegg)
Idents(pbmc.tumor.down.kegg)<-pbmc.tumor.down$Class

pbmc.tumor.down.tfs<-CreateSeuratObject(res_tfs)
Idents(pbmc.tumor.down.tfs)<-pbmc.tumor.down$Class

saveRDS(pbmc.tumor.down.hall, file = "../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.pathactivity.hall.rds")
saveRDS(pbmc.tumor.down.kegg, file = "../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.pathactivity.kegg.rds")
saveRDS(pbmc.tumor.down.tfs, file = "../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.pathactivity.tfs.rds")

##differential analysis
markers <- FindAllMarkers(pbmc.tumor.down.hall, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.)
write.table(markers,"../3.Analysis/3.2.diffpath.hall.txt",sep="\t",row.names=F,col.names=T,quote=F)

markers <- FindAllMarkers(pbmc.tumor.down.kegg, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.)
write.table(markers,"../3.Analysis/3.2.diffpath.kegg.txt",sep="\t",row.names=F,col.names=T,quote=F)

markers <- FindAllMarkers(pbmc.tumor.down.tfs, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.)
write.table(markers,"../3.Analysis/3.2.diffpath.tfs.txt",sep="\t",row.names=F,col.names=T,quote=F)

##plot
markers <- fread("../3.Analysis/3.2.diffpath.hall.txt", header = T, stringsAsFactors = F, data.table = F)
markers$logFDR <- -log10(markers$p_val_adj)
selmarkers <- markers[markers$gene %in% c("HALLMARK-INFLAMMATORY-RESPONSE", 
                                          "HALLMARK-IL2-STAT5-SIGNALING", "HALLMARK-FATTY-ACID-METABOLISM"), ]
p1 <- ggplot(markers, aes(x = avg_log2FC, y = logFDR, color = cluster))+
  geom_point()+
  xlab("avg log2FC")+
  ylab("-log(FDR)")+
  scale_color_manual(values = c("#FCE38AFF", "#D4B499FF", "#E64B35FF", "#4DBBD5FF"))+
  geom_label_repel(data = selmarkers, aes(x = avg_log2FC, y = logFDR, label = gene))+
  theme_classic2()
p1
ggsave("../3.Analysis/3.2.diffpath.hall.pdf", p1, width = 7, height = 5.5)

markers <- fread("../3.Analysis/3.2.diffpath.kegg.txt", header = T, stringsAsFactors = F, data.table = F)
markers$logFDR <- -log10(markers$p_val_adj)
selmarkers <- markers[markers$gene %in% c("Fatty acid metabolism", 
                                          "Jak-STAT signaling pathway", "Toll-like receptor signaling pathway"), ]
p1 <- ggplot(markers, aes(x = avg_log2FC, y = logFDR, color = cluster))+
  geom_point()+
  xlab("avg log2FC")+
  ylab("-log(FDR)")+
  scale_color_manual(values = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF", "#1a9641"))+
  geom_label_repel(data = selmarkers, aes(x = avg_log2FC, y = logFDR, label = gene))+
  theme_classic2()
p1
ggsave("../3.Analysis/3.2.diffpath.kegg.pdf", p1, width = 7, height = 5.5)

markers <- fread("../3.Analysis/3.2.diffpath.tfs.txt", header = T, stringsAsFactors = F, data.table = F)
markers$logFDR <- -log10(markers$p_val_adj)
markers$gene <- gsub("-TARGET-GENES", "", markers$gene)
selmarkers <- markers[markers$gene %in% c("NR1H4", 
                                          "ZNF597", "MAP2K1"), ]
p1 <- ggplot(markers, aes(x = avg_log2FC, y = logFDR, color = cluster))+
  geom_point()+
  xlab("avg log2FC")+
  ylab("-log(FDR)")+
  scale_color_manual(values = c("#FCE38AFF", "#4DBBD5FF", "#E64B35FF", "#1a9641"))+
  geom_label_repel(data = selmarkers, aes(x = avg_log2FC, y = logFDR, label = gene))+
  theme_classic2()
p1
ggsave("../3.Analysis/3.2.diffpath.tfs.pdf", p1, width = 7, height = 5.5)

###########3.2.3.subtype cell count###########
cellratio <- as.data.frame(table(pbmc.tumor@meta.data[, c("Class", "CompClass")]))
cellratio <- cellratio[cellratio$Freq > 0, ]

p <- ggplot()+
  geom_tile(data = cellratio, aes(x = CompClass, y = Class, fill = Freq))+
  geom_text(data = cellratio, aes(x = CompClass, y = Class, label = Freq))+
  theme_classic2()+
  scale_fill_gradient2(low="blue", high="red",mid = "white", midpoint = 0, na.value = "gray80")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p

ggsave("../3.Analysis/3.2.subtype.cell.count.pdf", p, width = 7, height = 4)

####################3.3.CrossTalk####################
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringsAsFactors = F)
pbmc <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.rds")
pbmc$CompClass <- as.character(pbmc$TopCluster)
pbmc$CompClass <- ifelse(pbmc$Celltype == "CD4T", "T cell CD4", pbmc$CompClass)
pbmc$CompClass <- ifelse(pbmc$Celltype == "CD8T", "T cell CD8", pbmc$CompClass)

pbmc3k.final <- subset(pbmc, subset = CellClass %in% c("Hepa", "Chol"), invert = TRUE)
pbmc3k.final <- subset(pbmc3k.final, subset = Tissue %in% c("LM", "PL"))

###########3.3.1.cellchat for topcluster###########
cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "CompClass")
summary(cellchat)
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
groupSize

CellChatDB <- CellChatDB.human
str(CellChatDB)
# 包含interaction、complex、cofactor和geneInfo这四个dataframe
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
#CellChatDB.use <- subsetData(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB

# 在矩阵的所有基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个)
cellchat <- subsetData(cellchat)
# 相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat <- identifyOverExpressedInteractions(cellchat)
# 上一步运行的结果存储在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)
# 找到配体受体关系后，projectData将配体受体对的表达质投射到PPI上，来对@data.signaling中的表达值进行校正。
# 结果保存在@data.project

# 根据表达值推测细胞互作的概率(cellphonedb是用平均表达值代表互作强度)
# 如果不想用上一步PPI矫正的结果，raw.use = TRUE即可
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
fwrite(df.net, "../3.Analysis/3.3.cellchat.malig2all.net_lr.txt", sep = "\t")
saveRDS(cellchat, "../3.Analysis/3.3.cellchat.malig2all.net_lr.rds")

# 统计细胞和细胞之间通信的数量(有多少个配体-受体对)和强度(概率)
cellchat <- aggregateNet(cellchat)
# 计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
pdf("../3.Analysis/3.3.cellchat.malig2all.net_lr.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
pdf("../3.Analysis/3.3.cellchat.malig2all.net_lr.sep.pdf")
for(i in 1:nrow(mat)){
  if (rownames(mat)[i] == "Epithelial cell") {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                     title.name = rownames(mat)[i])
  }
}
dev.off()

cellcommu <- subsetCommunication(cellchat)
cellcommu <- cellcommu[cellcommu$source == "Epithelial cell" & cellcommu$target != "Epithelial cell", ]
cellcommu <- top_n(group_by(cellcommu,target), n = 6, wt = prob)
cellcommu$ID <- paste(cellcommu$source, cellcommu$target, sep = " -> ")

p <- ggplot(cellcommu, aes(x = ID, y = interaction_name_2, color = prob))+
  geom_point(size = 3)+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", na.value = "gray80")+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p
ggsave("../3.Analysis/3.3.cellchat.malig2all.net_lr.pairs.pdf", p, width = 8, height = 4)


###########3.3.2.cellchat for subcluster###########
pbmc.tumor <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.rds")
Idents(pbmc.tumor)<-pbmc.tumor$Class
pbmc.tumor.down <- subset(x = pbmc.tumor, downsample = 1000)
markers <- FindAllMarkers(pbmc.tumor.down, only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
fwrite(markers, "../3.Analysis/3.3.cellchat.subtype2all.markers.txt")

type_d <- pbmc.tumor@meta.data[pbmc.tumor$Class == "D", ]$Patient
type_b <- pbmc.tumor@meta.data[pbmc.tumor$Class == "B", ]$Patient
type_t <- pbmc.tumor@meta.data[pbmc.tumor$Class == "T", ]$Patient
type_m <- pbmc.tumor@meta.data[pbmc.tumor$Class == "M", ]$Patient

pbmc3k.final$CompClass <- ifelse(pbmc3k.final$Patient %in% type_d & pbmc3k.final$CompClass == "Epithelial cell", "D", pbmc3k.final$CompClass)
pbmc3k.final$CompClass <- ifelse(pbmc3k.final$Patient %in% type_b & pbmc3k.final$CompClass == "Epithelial cell", "B", pbmc3k.final$CompClass)
pbmc3k.final$CompClass <- ifelse(pbmc3k.final$Patient %in% type_t & pbmc3k.final$CompClass == "Epithelial cell", "T", pbmc3k.final$CompClass)
pbmc3k.final$CompClass <- ifelse(pbmc3k.final$Patient %in% type_m & pbmc3k.final$CompClass == "Epithelial cell", "M", pbmc3k.final$CompClass)

cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "CompClass")
summary(cellchat)
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
groupSize

CellChatDB <- CellChatDB.human
str(CellChatDB)
# 包含interaction、complex、cofactor和geneInfo这四个dataframe
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
#CellChatDB.use <- subsetData(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB

# 在矩阵的所有基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个)
cellchat <- subsetData(cellchat)
# 相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat <- identifyOverExpressedInteractions(cellchat)
# 上一步运行的结果存储在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)
# 找到配体受体关系后，projectData将配体受体对的表达质投射到PPI上，来对@data.signaling中的表达值进行校正。
# 结果保存在@data.project

# 根据表达值推测细胞互作的概率(cellphonedb是用平均表达值代表互作强度)
# 如果不想用上一步PPI矫正的结果，raw.use = TRUE即可
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
fwrite(df.net, "../3.Analysis/3.3.cellchat.subtype2all.net_lr.txt", sep = "\t")
saveRDS(cellchat, "../3.Analysis/3.3.cellchat.subtype2all.net_lr.rds")

# 统计细胞和细胞之间通信的数量(有多少个配体-受体对)和强度(概率)
cellchat <- aggregateNet(cellchat)
# 计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))

mat <- cellchat@net$weight
pdf("../3.Analysis/3.3.cellchat.subtype2all.net_lr.sep.pdf", width = 20, height = 4)
par(mfrow = c(1,4), xpd = T)
for(i in 1:nrow(mat)){
  if (rownames(mat)[i] %in% c("B", "D", "M", "T")) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                     title.name = rownames(mat)[i])
  }
}
dev.off()

##heatmap show markers
tmpcellchat <- readRDS("../3.Analysis/3.3.cellchat.malig2all.net_lr.rds")
tmpcellcommu <- subsetCommunication(tmpcellchat)
tmpcellcommu <- tmpcellcommu[tmpcellcommu$source == "Epithelial cell" & tmpcellcommu$target != "Epithelial cell", ]
tmpcellcommu <- top_n(group_by(tmpcellcommu,target), n = 20, wt = prob)

markers <- fread("../3.Analysis/3.3.cellchat.subtype2all.markers.txt", data.table = F, stringsAsFactors = F)

res4plot <- markers[markers$gene %in% tmpcellcommu$ligand, ]
res4meta <- res4plot

res4meta$Class <- "Not"
res4meta$Class <- ifelse(res4meta$p_val_adj < 0.01, "Sig", res4meta$Class)
res4meta$gene <- factor(res4meta$gene, unique(res4meta$gene))

p <- ggplot()+
  geom_tile(data = res4meta, aes(x = gene, y = cluster, fill = avg_log2FC))+
  geom_text(data = res4meta[res4meta$Class == "Sig", ], aes(x = gene, y = cluster, label = "*"))+
  theme_classic2()+
  scale_fill_gradient2(low="blue", high="red",mid = "white", midpoint = 0, na.value = "gray80")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p

ggsave("../3.Analysis/3.3.cellchat.subtype2all.net_lr.topmarkers.pdf", p, width = 12, height = 2)

##dot plot
library(scales)
cellcommu <- subsetCommunication(cellchat)
cellcommu <- cellcommu[cellcommu$source %in% c("B", "D", "M", "T") & !(cellcommu$target %in% c("B", "D", "M", "T")), ]
cellcommu <- top_n(group_by(cellcommu,target), n = 15, wt = prob)
cellcommu$ID <- paste(cellcommu$source, cellcommu$target, sep = " -> ")

p <- ggplot(cellcommu, aes(x = ID, y = interaction_name_2, colour = prob))+
  geom_point(size = 3)+
  theme_bw()+
  scale_colour_gradient(low="blue", high="red")+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p
ggsave("../3.Analysis/3.3.cellchat.subtype2all.net_lr.pairs.pdf", p, width = 10, height = 5)



####################3.4.TCGA cross analysis####################
pbmc.tumor <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.rds")
Idents(pbmc.tumor) <- "CompClass"
pbmc.tumor <- subset(pbmc.tumor, subset = CompClass != "Epithelial cell")
pbmc.tumor <- subset(x = pbmc.tumor, downsample = 3000)

library(Scissor)
#Prepare the scRNA-seq data
sc_dataset <- pbmc.tumor@assays$RNA@counts

dim(sc_dataset)
sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)
names(sc_dataset)

cds.embed <- sc_dataset@reductions$umap
int.embed <- pbmc.tumor@reductions$umap

sc_dataset@reductions$umap <- int.embed

DimPlot(sc_dataset, reduction = 'umap', label = T, label.size = 10)
#Prepare the bulk data and phenotype
tcgaexp <- fread("/home/yukai/data/Datasets/TCGA20201022/htseq/fpkm_gene/TCGA-LIHC.htseq_fpkm.tsv",
                 header = T, stringsAsFactors = F, data.table = F)
rownames(tcgaexp) <- tcgaexp$Ensembl_ID
tcgaexp <- tcgaexp[, -1]

survinfo <- fread("/home/yukai/data/Datasets/TCGA20201022/survival_TCGA/LIHC_survival.txt",
                  header = T, stringsAsFactors = F, data.table = F)
survinfo$sample <- paste(survinfo$sample, "A", sep = "")
rownames(survinfo) <- survinfo$sample

comsams <- intersect(rownames(survinfo), names(tcgaexp))
tcgaexp <- tcgaexp[, comsams]
survinfo <- survinfo[comsams, ]

###########3.4.0.test.os###########
for (selnum in c(1:20)) {
  bulk_dataset <- as.matrix(tcgaexp)
  set.seed(selnum)  # 设置随机数种子，确保每次运行结果相同  
  random_integers <- sample(1:421, 200, replace = FALSE)  
  bulk_dataset <- bulk_dataset[, random_integers]
  phenotype <- survinfo[random_integers, c("OS.time", "OS")]
  names(phenotype) <- c("time", "status")
  infos2 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.0001, 
                    family = "cox", Save_file = paste("../3.Analysis/3.4.Scissor_LIHC_survival.os.", selnum, ".RData"))
  
  Scissor_select <- rep(0, ncol(sc_dataset))
  names(Scissor_select) <- colnames(sc_dataset)
  Scissor_select[infos2$Scissor_pos] <- 1
  Scissor_select[infos2$Scissor_neg] <- 2
  sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
  DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1))
  
  pdf(paste("../3.Analysis/3.4.Scissor_LIHC_survival.os.", selnum, ".pdf"), width = 6, height = 5)
  print(DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1)))
  dev.off()
  
  fwrite(cbind(sc_dataset@meta.data, pbmc.tumor@meta.data[, c("Patient", "CompClass")]), 
         file = paste("../3.Analysis/3.4.Scissor_LIHC_survival.os.", selnum, ".txt"), 
         quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
}
###########3.4.0.test.pfs###########
for (selnum in c(1:20)) {
  bulk_dataset <- as.matrix(tcgaexp)
  set.seed(selnum)  # 设置随机数种子，确保每次运行结果相同  
  random_integers <- sample(1:421, 200, replace = FALSE)  
  bulk_dataset <- bulk_dataset[, random_integers]
  phenotype <- survinfo[random_integers, c("PFI.time", "PFI")]
  names(phenotype) <- c("time", "status")
  infos2 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.0001, 
                    family = "cox", Save_file = paste("../3.Analysis/3.4.Scissor_LIHC_survival.pfs.", selnum, ".RData"))
  
  Scissor_select <- rep(0, ncol(sc_dataset))
  names(Scissor_select) <- colnames(sc_dataset)
  Scissor_select[infos2$Scissor_pos] <- 1
  Scissor_select[infos2$Scissor_neg] <- 2
  sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
  DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1))
  
  pdf(paste("../3.Analysis/3.4.Scissor_LIHC_survival.pfs.", selnum, ".pdf"), width = 6, height = 5)
  print(DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1)))
  dev.off()
  
  fwrite(cbind(sc_dataset@meta.data, pbmc.tumor@meta.data[, c("Patient", "CompClass")]), 
         file = paste("../3.Analysis/3.4.Scissor_LIHC_survival.pfs.", selnum, ".txt"), 
         quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
}
###########3.4.1.survival.os###########
bulk_dataset <- as.matrix(tcgaexp)
phenotype <- survinfo[, c("OS.time", "OS")]
names(phenotype) <- c("time", "status")
#Execute Scissor to select the informative cells
infos2 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.0001, 
                  family = "cox", Save_file = '../3.Analysis/3.4.Scissor_LIHC_survival.os.RData')

Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos2$Scissor_pos] <- 1
Scissor_select[infos2$Scissor_neg] <- 2
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1))

pdf('../3.Analysis/3.4.Scissor_LIHC_survival.os.pdf', width = 6, height = 5)
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1))
dev.off()

fwrite(cbind(sc_dataset@meta.data, pbmc.tumor@meta.data[, c("Patient", "CompClass")]), 
       file = "../3.Analysis/3.4.Scissor_LIHC_survival.os.txt", 
       quote = FALSE, sep = '\t', row.names = T, col.names = T)

###########3.4.2.survival.pfs###########
bulk_dataset <- as.matrix(tcgaexp)
phenotype <- survinfo[, c("PFI.time", "PFI")]
names(phenotype) <- c("time", "status")
#Execute Scissor to select the informative cells
infos2 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.0001, 
                  family = "cox", Save_file = '../3.Analysis/3.4.Scissor_LIHC_survival.pfs.RData')

Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos2$Scissor_pos] <- "1"
Scissor_select[infos2$Scissor_neg] <- "2"
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1))

pdf('../3.Analysis/3.4.Scissor_LIHC_survival.pfs.pdf', width = 6, height = 5)
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1))
dev.off()

fwrite(cbind(sc_dataset@meta.data, pbmc.tumor@meta.data[, c("Patient", "CompClass")]), 
       file = "../3.Analysis/3.4.Scissor_LIHC_survival.pfs.txt", 
       quote = FALSE, sep = '\t', row.names = T, col.names = T)


###########3.4.3.mutation landscape###########
library(maftools)
mafile <- read.delim("/home/yukai/data/TCGA/mutation_maf/LIHC.Mutation_filter.txt", header = T, stringsAsFactors = F, sep = '\t')

oncotsg <- read.delim("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/EMT.Meta.Metabolic.oncoTSG.list", 
                      header =F, sep = '\t')
mafile = mafile[mafile$Hugo_Symbol %in% oncotsg$V3,]

flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
          "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
          "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")

laml1 = read.maf(maf = '/home/yukai/data/TCGA/mutation_maf/LIHC.Mutation_filter.txt')

col = ggsci::pal_npg("nrc")(10)
names(col) = c('Frame_Shift_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Multi_Hit', 
               'Nonsense_Mutation', 'In_Frame_Del', 'Splice_Site', 'Frame_Shift_Ins',
               'Nonstop_Mutation', 'Translation_Start_Site')
#oncoplot(maf = laml1, top = 15, removeNonMutated = F, colors = col)
pdf("../3.Analysis/3.4.LIHC.mutation.landscape.pdf", width = 5, height = 4)
oncoplot(maf = laml1, top = 15, removeNonMutated = F, colors = col)
dev.off()

###########3.4.4.mutation and scdata###########
allgenes <- c("TP53", "CTNNB1", "APOB", "TTN", "MUC16")

for (selgene in allgenes) {
  #selgene <- "TP53"
  #Apply Scissor with logistic regression
  mafile <- read.delim("/home/yukai/data/TCGA/mutation_maf/LIHC.Mutation_filter.txt", header = T, stringsAsFactors = F, sep = '\t')
  allmut <- unique(mafile[mafile$Hugo_Symbol == selgene, ]$Tumor_Sample_Barcode)
  allnon <- setdiff(unique(mafile$Tumor_Sample_Barcode), allmut)
  
  allsams <- data.frame(ID = c(allmut, allnon),
                        Class = c(rep(1, length(allmut)), rep(0, length(allnon))))
  rownames(allsams) <- allsams$ID
  
  tmptcgaexp <- tcgaexp
  comsams <- intersect(allsams$ID, names(tmptcgaexp))
  tmptcgaexp <- tmptcgaexp[, comsams]
  allsams <- allsams[comsams, ]
  
  bulk_dataset <- as.matrix(tmptcgaexp)
  TP53_mutation <-  allsams$Class
  TP53_mutation <- setNames(TP53_mutation, allsams$ID)
  table(TP53_mutation)
  phenotype <- TP53_mutation
  tag <- c('wild-type', 'TP53 mutant')
  infos4 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = 0.0001, 
                    family = "binomial", Save_file = paste("../3.Analysis/3.4.Scissor_LIHC_", selgene, ".mut.RData"))
  
  Scissor_select <- rep(0, ncol(sc_dataset))
  names(Scissor_select) <- colnames(sc_dataset)
  Scissor_select[infos4$Scissor_pos] <- 1
  Scissor_select[infos4$Scissor_neg] <- 2
  sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
  DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
  
  pdf(paste("../3.Analysis/3.4.Scissor_LIHC_", selgene, ".mut.pdf"), width = 6, height = 5)
  print(DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), order = c(2,1)))
  dev.off()
  
  fwrite(cbind(sc_dataset@meta.data, pbmc.tumor@meta.data[, c("Patient", "CompClass")]), 
         file = paste("../3.Analysis/3.4.Scissor_LIHC_", selgene, ".mut.txt"), 
         quote = FALSE, sep = '\t', row.names = T, col.names = T)
  
  
}

###########3.4.5.cell group difference for each data###########
for (selnum in c(1:20)) {
  #os
  meta.tb <- fread(paste("../3.Analysis/3.4.Scissor_LIHC_survival.os.", selnum, ".txt"), header = T, stringsAsFactors = F, data.table = F)
  out.prefix <- "../3.Analysis/3.4.cellcount"
  
  cellInfo.tb <- meta.tb[meta.tb$scissor != 0, ]
  OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                               meta.cluster = cellInfo.tb$CompClass,
                               colname.patient = "Patient",
                               loc = cellInfo.tb$scissor,
                               out.prefix=sprintf(paste("%s.STARTRAC.dist.TopCluster.os.", selnum) ,out.prefix),
                               pdf.width=4,pdf.height=6,verbose=1)
  
  fwrite(as.data.frame(OR.all.list$OR.dist.mtx), paste("../3.Analysis/3.4.cellcount.STARTRAC.dist.TopCluster.os.OR.dist.", selnum, ".txt"), 
         row.names = T, quote = F)
  
  #pfs
  meta.tb <- fread(paste("../3.Analysis/3.4.Scissor_LIHC_survival.pfs.", selnum, ".txt"), header = T, stringsAsFactors = F, data.table = F)
  out.prefix <- "../3.Analysis/3.4.cellcount"
  
  cellInfo.tb <- meta.tb[meta.tb$scissor != 0, ]
  OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                               meta.cluster = cellInfo.tb$CompClass,
                               colname.patient = "Patient",
                               loc = cellInfo.tb$scissor,
                               out.prefix=sprintf(paste("%s.STARTRAC.dist.TopCluster.pfs.", selnum) ,out.prefix),
                               pdf.width=4,pdf.height=6,verbose=1)
  
  fwrite(as.data.frame(OR.all.list$OR.dist.mtx), paste("../3.Analysis/3.4.cellcount.STARTRAC.dist.TopCluster.pfs.OR.dist.", selnum, ".txt"), 
         row.names = T, quote = F)
}

#mutations
allgenes <- c("TP53", "CTNNB1", "APOB", "TTN", "MUC16")

for (selgene in allgenes) {
  meta.tb <- fread(paste("../3.Analysis/3.4.Scissor_LIHC_", selgene, ".mut.txt"), header = T, stringsAsFactors = F, data.table = F)
  out.prefix <- "../3.Analysis/3.4.cellcount"
  
  cellInfo.tb <- meta.tb[meta.tb$scissor != 0, ]
  OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                               meta.cluster = cellInfo.tb$CompClass,
                               colname.patient = "Patient",
                               loc = cellInfo.tb$scissor,
                               out.prefix=sprintf(paste("%s.STARTRAC.dist.TopCluster.", selgene, sep = ""),out.prefix),
                               pdf.width=4,pdf.height=6,verbose=1)
  
  fwrite(as.data.frame(OR.all.list$OR.dist.mtx), paste("../3.Analysis/3.4.cellcount.STARTRAC.dist.TopCluster.", selgene, ".OR.dist.txt", sep = ""), 
         row.names = T, quote = F)
  
}


###########3.4.6.tcga xcell###########
tcgaexp <- fread("/home/yukai/data/Datasets/TCGA20201022/htseq/fpkm_gene/TCGA-LIHC.htseq_fpkm.tsv",
                 header = T, stringsAsFactors = F, data.table = F)
rownames(tcgaexp) <- tcgaexp$Ensembl_ID
tcgaexp <- tcgaexp[, -1]
library(xCell)
res <- xCellAnalysis(tcgaexp, signatures = NULL, genes = NULL, spill = NULL,
                     rnaseq = TRUE, file.name = NULL, scale = TRUE, alpha = 0.5,
                     save.raw = FALSE, parallel.sz = 4, parallel.type = "SOCK",
                     cell.types.use = NULL)
res <- res[-c(nrow(res)-2, nrow(res)-1, nrow(res)), ]
write.table(res, "../3.Analysis/3.4.tcga.scores.immunecell.xcell.txt", row.names = T, col.names = NA, sep = "\t", quote = F)

###########3.4.7.xcell and mutation###########
selgene <- "TP53"
#NK, mast, B, plasma
#Apply Scissor with logistic regression
mafile <- read.delim("/home/yukai/data/TCGA/mutation_maf/LIHC.Mutation_filter.txt", header = T, stringsAsFactors = F, sep = '\t')
allmut <- unique(mafile[mafile$Hugo_Symbol == selgene, ]$Tumor_Sample_Barcode)
allnon <- setdiff(unique(mafile$Tumor_Sample_Barcode), allmut)

allsams <- data.frame(ID = c(allmut, allnon),
                      Class = c(rep("Mut", length(allmut)), rep("WT", length(allnon))))
rownames(allsams) <- allsams$ID

xcellinfo <- fread("../3.Analysis/3.4.tcga.scores.immunecell.xcell.txt",
                   header = T, stringsAsFactors = F, data.table = F)
rownames(xcellinfo) <- xcellinfo$V1
xcellinfo <- xcellinfo[, -1]

xcellinfo <- as.data.frame(t(xcellinfo))
xcellinfo$ID <- rownames(xcellinfo)
xcellinfo <- merge(xcellinfo, allsams, by.x = "ID", by.y = "ID")

selcells <- c("CD8+ T-cells", "Mast cells", "B-cells", "Endothelial cells")
xcellinfo <- xcellinfo[, c(selcells, "Class")]
res <- melt(xcellinfo, id.vars = "Class")

p2 <- ggplot(res, aes(x = Class, y = value, color = Class))+
  #geom_violin(aes(fill = Type), scale = "width")+
  geom_boxplot(outlier.shape = NA, width = 0.7, )+
  #geom_sina()+
  #stat_boxplot(aes(ymin=..lower..,ymax=..upper..), outlier.shape = NA)+
  stat_compare_means(label = "p.signif", step.increase = 0)+
  ylab("fraction")+
  scale_color_manual(values = c("#6495EDFF", "#FF4500FF"))+
  theme_bw()
p2
q = facet(p2, facet.by = "variable", scales = "free_y", nrow = 1)
q
pdf("../3.Analysis/3.4.tcga.scores.immunecell.xcell.mut.pdf", width = 10, height = 3)
q
dev.off()



###########3.4.8.xcell and survival###########
xcellinfo <- fread("../3.Analysis/3.4.tcga.scores.immunecell.xcell.txt",
                   header = T, stringsAsFactors = F, data.table = F)
rownames(xcellinfo) <- xcellinfo$V1
xcellinfo <- xcellinfo[, -1]
cli_info <- fread("/home/yukai/data/Datasets/TCGA20201022/survival_TCGA/LIHC_survival.txt",
                  header = T, sep = '\t', data.table  = F)
rownames(cli_info) <- paste(cli_info$sample, "A", sep = "")
cli_info <- cli_info[as.numeric(substr(rownames(cli_info), 14, 15)) < 10, ]


for (gene in rownames(xcellinfo)) {
  survres <- surv_bestcut(xcellinfo, gene, cli_info, num = 20)
  print(survres[order(survres$Pvalue_OS), ][1, ])
  survres <- surv_bestcut(xcellinfo, gene, cli_info, num = survres[order(survres$Pvalue_OS), ][1, 8])
}

##Fibroblast (Monocytes)
##pDC (Th1)
##Mast cell (Mast cell)
##NK cell (pDC)
##CD8+ Tcell (CD8+ T-cells)


####################3.5.Key subgroups####################
#CAF：FAP, ACTA2
#fibroblasts：COL1A1, FN1
###########3.5.1.class TA and NA / featureplot###########
pbmc.tumor <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.rds")
Idents(pbmc.tumor)<-pbmc.tumor$TopCluster
sub_sce<-subset(pbmc.tumor,idents=c("Fibroblast"))
#combined<-readRDS("2.intergrate/Combined.rds")

DefaultAssay(sub_sce) <- "integrated"

sub_sce <- ScaleData(sub_sce, verbose = FALSE)
sub_sce <- FindVariableFeatures(sub_sce, selection.method = "vst", nfeatures = 1000)
sub_sce <- RunPCA(sub_sce, npcs = 30, verbose = FALSE)

ElbowPlot(sub_sce)
#chosen.elbow <- findElbowPoint(Stdev(object = sub_sce, reduction = "harmony"))
sub_sce <-RunHarmony(sub_sce,"Dataset",reduction = "pca",assay.use = "integrated",plot_convergence = TRUE)

sub_sce <- RunUMAP(sub_sce, reduction = "harmony", dims = 1:15)
sub_sce <- RunTSNE(sub_sce, reduction = "harmony", dims = 1:15, check_duplicates = FALSE)
sub_sce <- FindNeighbors(sub_sce, reduction = "harmony", dims = 1:15)

sub_sce <- FindClusters(sub_sce, resolution = 1.2)
DimPlot(sub_sce, reduction = "umap", label = TRUE)

markers <- FindAllMarkers(sub_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../3.Analysis/3.5.sub_sce_Cluster_markers.txt",sep="\t",row.names=T,col.names=T,quote=F)


top10 <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
pdf ("../3.Analysis/3.5.sub_sce_Cluter_DotPlot.pdf", height = 8, width = 30)
print(DotPlot(sub_sce, features = unique(c(top10$gene,"ICOS")),
              cols = c("grey","blue"))+RotatedAxis()+
        scale_x_discrete("")+scale_y_discrete(""))
dev.off()

## minor level clustering
DimPlot(sub_sce, reduction = "umap", label = TRUE)
pdf("../3.Analysis/3.5.sub_sce_Cluter_replot.markers.pdf",width = 10,height = 5)
FeaturePlot(sub_sce,features = c("ACTA2", "FAP", ##CAFs
                                 "ACTG2", ##MyoFibroblast
                                 "COL1A1", "FN1" ##Fibroblast,
)
,ncol = 3,min.cutoff = 0,raster = T)
dev.off()
pdf("../3.Analysis/3.5.sub_sce_Cluter_replot.markers.pdf",width = 10,height = 5)
FeaturePlot(sub_sce,features = c("ACTA2", "FAP", ##CAFs
                                 "ACTG2", ##MyoFibroblast
                                 "COL1A1", "FN1" ##Fibroblast,
)
,ncol = 3,min.cutoff = 0,raster = T)

pdf("../3.Analysis/3.5.sub_sce_Cluter_replot.markers.vlnplot.pdf",width = 6,height = 4)
VlnPlot(sub_sce, c("ACTA2", "FAP", ##CAFs
                   "ACTG2", ##MyoFibroblast
                   "COL1A1", "FN1" ##Fibroblast,
                   ), stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()


markers <- read.delim("../3.Analysis/3.5.sub_sce_Cluster_markers.txt", header = T)
top10genes <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)

new.cluster.id<-c("NAF_PTP4A3", "NAF_ACTG2", "NAF_SOD3", "CAF_MMP2",
                  "NAF_FABP5", "CAF_CXCL12", "NAF_JUNB", "NAF_RGS5",
                  "NAF_JUN", "CAF_CXCL14", "CAF_ASPN", "CAF_CD74",
                  "CAF_CXCL12", "NAF_PTPRC") 

names(new.cluster.id)<-levels(sub_sce)
sub_sce<-RenameIdents(sub_sce,new.cluster.id)

sub_sce$SelCellCluster <- Idents(sub_sce)
sub_sce$SelCellCelltype <- strsplit2(sub_sce$SelCellCluster, "_")[, 1]
sub_sce$SelCellGenes <- strsplit2(sub_sce$SelCellCluster, "_")[, 2]

tmpdf <- data.frame(SelCellCluster = sub_sce$SelCellCluster,
                    SelCellCelltype = sub_sce$SelCellCelltype,
                    SelCellGenes = sub_sce$SelCellGenes)
tmpdf <- tmpdf[order(tmpdf$SelCellCelltype, tmpdf$SelCellCluster), ]

Idents(sub_sce)<-factor(Idents(sub_sce),levels= as.character(unique(tmpdf$SelCellCluster)))                
saveRDS(sub_sce, file = "../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")

my_cols <- c('NAF'='#31C53F',
             'CAF'='#ff9a36')

p2 <- DimPlot(sub_sce, reduction = "umap", group.by = "SelCellCelltype",cols =  my_cols, raster=TRUE, label = TRUE) + NoLegend()
p2
pdf("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.pdf", width = 5, height = 5)
p2
dev.off()

###########3.5.2.class in lauren###########
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")

meta.tb <- sub_sce@meta.data
out.prefix <- "../3.Analysis/3.5.cellcount"

cellInfo.tb <- meta.tb
OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SelCellCelltype,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.TopCluster.lauren",out.prefix),
                             pdf.width=4,pdf.height=3,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.5.cellcount.STARTRAC.dist.TopCluster.lauren.OR.dist.txt", 
       row.names = T, quote = F)

###########3.5.3.marker expression###########
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")
Idents(sub_sce) <- "SelCellCelltype"

pdf("../3.Analysis/3.5.sub_sce.marker.exp.pdf",width = 8,height = 3)
FeaturePlot(sub_sce,features = c("FAP", "COL1A1"
),ncol = 3,min.cutoff = 0,raster = T)
dev.off()

selgene <- c("FAP", "COL1A1")
pdf("../3.Analysis/3.5.sub_sce.marker.vlnplot.exp.pdf",width = 6,height = 4)
VlnPlot(sub_sce, selgene, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()


###########3.5.4.diff transcriptors###########
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")
Idents(sub_sce) <- "SelCellCelltype"
sub_sce.down <- subset(x = sub_sce, downsample = 2000)

pbmc.tumor.mtr <- sub_sce.down@assays$RNA@counts

library(clusterProfiler)
library(msigdbr)
library(dplyr)
tfs <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>%
  dplyr::select(gs_name, gene_symbol)

tfsgmt <- split(tfs$gene_symbol, tfs$gs_name)

res_tfs <- gsva(as.matrix(pbmc.tumor.mtr), tfsgmt, kcdf="Gaussian",method = "ssgsea",parallel.sz=10)

pbmc.tumor.down.tfs<-CreateSeuratObject(res_tfs)
Idents(pbmc.tumor.down.tfs)<-sub_sce.down$SelCellCelltype
saveRDS(pbmc.tumor.down.tfs, file = "../3.Analysis/3.5.sub_sce.pathactivity.tfs.rds")

##differential analysis
markers <- FindAllMarkers(pbmc.tumor.down.tfs, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.)
write.table(markers,"../3.Analysis/3.5.diffpath.tfs.txt",sep="\t",row.names=F,col.names=T,quote=F)

##plot
markers <- fread("../3.Analysis/3.5.diffpath.tfs.txt", header = T, stringsAsFactors = F, data.table = F)
markers$gene <- gsub("-TARGET-GENES", "", markers$gene)
markers$gene <- strsplit2(markers$gene, "-")[, 1]
selmarkers <- top_n(group_by(markers,cluster), n = 10, wt = avg_log2FC)
selmarkers$avg_log2FC <- ifelse(selmarkers$cluster == "NAF", -selmarkers$avg_log2FC, selmarkers$avg_log2FC)
selmarkers <- selmarkers[order(selmarkers$avg_log2FC, decreasing = T), ]
selmarkers$gene <- factor(selmarkers$gene, selmarkers$gene)

p2 <- ggplot(selmarkers, aes(x = gene, y = avg_log2FC, fill = cluster))+
  geom_bar(stat = "identity", color = "black", width = 0.7)+
  theme_classic()+
  #scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_manual(values = c('#ff9a36', '#31C53F'))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
p2
ggsave("../3.Analysis/3.5.diffpath.tfs.pdf", p2, width = 6, height = 3)

####################3.6.detail subclassification####################
###########3.6.1.detail clusters###########
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")
Idents(sub_sce) <- "SelCellCluster"
Idents(sub_sce)<-factor(as.character(Idents(sub_sce)),levels= sort(as.character(unique(sub_sce$SelCellCluster))))  

allgenes <- c(unique(unique(sub_sce$SelCellGenes)), "ACTA2")
pdf("../3.Analysis/3.6.sub_sce.marker.vlnplot.exp.pdf",width = 7,height = 5)
VlnPlot(sub_sce, allgenes, stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()


Idents(sub_sce) <- "SelCellCluster"
markers <- FindAllMarkers(sub_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../3.Analysis/3.5.sub_sce_Cluster_markers.SelCellCluster.txt",sep="\t",row.names=T,col.names=T,quote=F)


Idents(sub_sce) <- "SelCellCelltype"
markers <- FindAllMarkers(sub_sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../3.Analysis/3.5.sub_sce_Cluster_markers.SelCellCelltype.txt",sep="\t",row.names=T,col.names=T,quote=F)

###########3.6.2.monocyte###########
#dyn.load("/home/rstudio_cu02/anaconda3/envs/R4.0/lib/./libspatialite.so.7")
library(monocle3)
library(Seurat)
library(tidyverse)
library(patchwork)

dir.create("../3.Analysis/Monocle3")
combined <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")
Idents(combined) <- 'SelCellCluster'
##创建CDS对象并预处理数据
pbmc <- combined
data <- GetAssayData(pbmc, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#UMAP
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "SelCellCluster")

#import umap from seurat
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]

cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "SelCellCluster")
p <- p1|p2
ggsave("../3.Analysis/Monocle3/reduction.umap.pdf", p, width = 6, height = 3)

#visualize some genes
ciliated_genes <- c("ACTA2")
plot_cells(cds, 
           genes = ciliated_genes,
           label_cell_groups =FALSE,
           show_trajectory_graph = F)

#cluster your cells for trajectory
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "SelCellCluster",
               label_groups_by_cluster =FALSE,
               show_trajectory_graph = TRUE)
ggsave("../3.Analysis/Monocle3/reduction.umap.trajectory.pdf", p, width = 5, height = 4)

#pseudotime
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime",
               label_groups_by_cluster =FALSE,
               show_trajectory_graph = TRUE)
p
ggsave("../3.Analysis/Monocle3/reduction.umap.trajectory.pseudotime.pdf", p, width = 5, height = 4)

###########3.6.3.cell number among different datasets###########
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")

cellratio <- as.data.frame(table(sub_sce@meta.data[, c("Tissue", "SelCellCluster")]))
cellratio <- cellratio[cellratio$Freq > 0, ]

p <- ggplot()+
  geom_tile(data = cellratio, aes(x = SelCellCluster, y = Tissue, fill = Freq))+
  geom_text(data = cellratio, aes(x = SelCellCluster, y = Tissue, label = Freq))+
  theme_classic()+
  scale_fill_gradient2(low="blue", high="red",mid = "white", midpoint = 0, na.value = "gray80")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p

ggsave("../3.Analysis/3.6.subtype.cell.count.pdf", p, width = 7, height = 4)


cellratio <- prop.table(table(sub_sce$SelCellCluster, as.factor(as.character(sub_sce$Dataset))), margin = 2)
cellratio <- as.data.frame(cellratio)
cellratio$Var1 <- factor(as.character(cellratio$Var1), sort(unique(as.character(cellratio$Var1))))
p <- ggplot(cellratio) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1), stat = "sum", width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y ='Ratio')+
  coord_flip()+
  #scale_fill_manual(values = cell_type_cols)+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))
p
pdf("../3.Analysis/3.6.dataset.cellcount.pdf", width = 8, height = 5)
p
dev.off()

###########3.6.4.cell number count among subtypes###########
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")

meta.tb <- sub_sce@meta.data
out.prefix <- "../3.Analysis/3.6.cellcount"

cellInfo.tb <- meta.tb
OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$SelCellCluster,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.TopCluster.lauren",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.6.cellcount.STARTRAC.dist.TopCluster.lauren.OR.dist.txt", 
       row.names = T, quote = F)


###########3.6.5.crosstalk###########
pbmc.tumor <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.rds")
Idents(pbmc.tumor) <- "CompClass"
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")
pbmc3k.immune <- subset(pbmc.tumor, subset = CompClass %in% c("T cell CD8", "Epithelial cell"))
sub_sce$CompClass <- sub_sce$SelCellCluster

pbmc3k.final <- merge(pbmc3k.immune, y = sub_sce, add.cell.ids = c("4K", "8K"), project = "PBMC12K")

cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "CompClass")
summary(cellchat)
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
groupSize

CellChatDB <- CellChatDB.human
str(CellChatDB)
# 包含interaction、complex、cofactor和geneInfo这四个dataframe
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
#CellChatDB.use <- subsetData(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB

# 在矩阵的所有基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个)
cellchat <- subsetData(cellchat)
# 相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
cellchat <- identifyOverExpressedInteractions(cellchat)
# 上一步运行的结果存储在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human)
# 找到配体受体关系后，projectData将配体受体对的表达质投射到PPI上，来对@data.signaling中的表达值进行校正。
# 结果保存在@data.project

# 根据表达值推测细胞互作的概率(cellphonedb是用平均表达值代表互作强度)
# 如果不想用上一步PPI矫正的结果，raw.use = TRUE即可
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = T)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
fwrite(df.net, "../3.Analysis/3.6.cellchat.malig2all.net_lr.txt", sep = "\t")
saveRDS(cellchat, "../3.Analysis/3.6.cellchat.malig2all.net_lr.rds")

# 统计细胞和细胞之间通信的数量(有多少个配体-受体对)和强度(概率)
cellchat <- aggregateNet(cellchat)
# 计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
pdf("../3.Analysis/3.6.cellchat.malig2all.net_lr.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge = F, title.name = "Interaction weights/strength")
dev.off()

mat <- cellchat@net$weight
pdf("../3.Analysis/3.6.cellchat.malig2all.net_lr.sep.pdf", width = 12, height = 12)
par(mfrow = c(4,4), xpd=TRUE)
for(i in 1:nrow(mat)){
  if (substr(rownames(mat)[i], 1, 3) %in% c("CAF", "NAF")) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i,] <- mat[i,]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                     title.name = rownames(mat)[i])
  }
}
dev.off()

cellcommu <- subsetCommunication(cellchat)
cellcommu <- cellcommu[substr(cellcommu$source, 1, 3) == c("CAF", "NAF") & cellcommu$target %in% c("Epithelial cell", "T cell CD8"), ]
cellcommu <- top_n(group_by(cellcommu,source), n = 10, wt = prob)
cellcommu$target <- as.character(cellcommu$target)
cellcommu$target <- ifelse(cellcommu$target == "Epithelial cell", "Cancer cell", cellcommu$target)
cellcommu$ID <- paste(cellcommu$source, cellcommu$target, sep = " -> ")

p <- ggplot(cellcommu, aes(x = ID, y = interaction_name_2, color = prob))+
  geom_point(size = 3)+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", na.value = "gray80")+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p
ggsave("../3.Analysis/3.6.cellchat.malig2all.net_lr.pairs.pdf", p, width = 10, height = 6)




####################3.7.signature application####################
###########3.7.1.marker gene selection###########
pbmc.tumor <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.tumor.rds")
sub_sce <- readRDS("../3.Analysis/3.5.sub_sce_Cluster.pca.minor.rds")

selsamples <- rownames(sub_sce@meta.data[sub_sce$SelCellCelltype == "CAF", ])
pbmc.tumor$SelCellCelltype <- ifelse(rownames(pbmc.tumor@meta.data) %in% selsamples, "CAF", "Others")

pbmc.tumor.down <- subset(x = pbmc.tumor, downsample = 5000)
Idents(pbmc.tumor.down) <- "SelCellCelltype"
markers <- FindAllMarkers(pbmc.tumor.down, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"../3.Analysis/3.7.sub_sce_Cluster_markers.SelCellCluster.txt",sep="\t",row.names=T,col.names=T,quote=F)
top40 <-  top_n(group_by(markers,cluster), n = 40, wt = avg_log2FC)
write.table(top40,"../3.Analysis/3.7.sub_sce_Cluster_markers.SelCellCluster.top40.txt",sep="\t",row.names=T,col.names=T,quote=F)


top10 <-  top_n(group_by(markers,cluster), n = 40, wt = avg_log2FC)
top10 <- top10[top10$cluster == "CAF", ]
pdf("../3.Analysis/3.7.sub_sce_Cluster_markers.heatmap.pdf",width = 12,height = 3)
VlnPlot(pbmc.tumor.down, top10$gene, stack = TRUE, sort = TRUE, flip = FALSE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()


###########3.7.2.build score###########
surv_bestcut <- function(mtr, gene, cli_info, num = 20){
  tmp_mtr <- mtr[which(rownames(mtr) == gene), ]
  if (length(tmp_mtr[is.na(tmp_mtr)]) == 0) {
    tmp_mtr <- tmp_mtr
  }else{
    tmp_mtr <- tmp_mtr[-which(is.na(tmp_mtr))]
  }
  common_samples <- intersect(names(tmp_mtr), rownames(cli_info))
  cluster_surv <- cli_info[common_samples, ]
  tmp_mtr <- as.data.frame(tmp_mtr)[, common_samples]
  sevalue <- as.numeric(tmp_mtr)
  values <- c()
  hr_os <- c()
  pva_os <- c()
  n_high <- c()
  n_low <- c()
  for (i in c(round(length(sevalue)/4):round(length(sevalue)/4*3))) {
    cluster_surv$Type = ifelse(sevalue > sort(sevalue)[i], "1.High", "0.Low")
    values <- c(values, sort(sevalue)[i])
    tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
    hr_os <- c(hr_os, tmp$conf.int[[1]])
    pva_os <- c(pva_os, tmp$logtest[[3]])
    n_high <- c(n_high, (length(sevalue)-i))
    n_low <- c(n_low, i)
    if (i == num) {
      ##OS
      tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
      fit <- survfit(Surv(OS.time, OS) ~ Type, data = cluster_surv)
      os <- survplot(cluster_surv, type = paste("OS", gene, sep = "-"), fit = fit, pval = tmp$logtest[3])
      ##PFS
      print(os)
    }
  }
  res <- data.frame(ID = gene,
                    cutoff = values,
                    HR_OS = hr_os,
                    Pvalue_OS = pva_os,
                    n_high = n_high,
                    n_low = n_low)
  
}

##### lasso cox model #####
desurv <- read.delim("../3.Analysis/3.7.sub_sce_Cluster_markers.SelCellCluster.top40.txt", header = T, sep = '\t', stringsAsFactors = F)
genelist <- desurv[desurv$cluster == "CAF", ]$gene
tcgaMtr <- fread("/home/yukai/data/Datasets/TCGA20201022/htseq/fpkm_gene/TCGA-LIHC.htseq_fpkm.tsv", header = T, sep = '\t', check.names = F, stringsAsFactors = F, data.table = F)
rownames(tcgaMtr) <- tcgaMtr$Ensembl_ID
tcgaMtr <- tcgaMtr[, -1]
colData <- data.frame(Type = ifelse(substr(colnames(tcgaMtr), 14, 15) < 10, 'Tumor', 'Normal'),
                      ID = colnames(tcgaMtr))
colData$ID <- as.character(colData$ID)
colData <- colData[colData$Type == 'Tumor', ]
tcgaMtr <- tcgaMtr[, colData$ID]
tcgaMtr4meta <- tcgaMtr[genelist, ]
tcgaMtr4meta <- na.omit(tcgaMtr4meta)

exprSet=tcgaMtr4meta
used_genes = rownames(exprSet)
survinfo <- fread("/home/yukai/data/Datasets/TCGA20201022/survival_TCGA/LIHC_survival.txt",
                  header = T, stringsAsFactors = F, data.table = F)
survinfo$sample <- paste(survinfo$sample, "A", sep = "")
rownames(survinfo) <- survinfo$sample
meta <- na.omit(survinfo[, c(1:4)])
common_sample = intersect(meta$sample, colnames(exprSet))
rownames(meta) <- meta$sample
meta <- meta[common_sample, ]
exprSet <- exprSet[, common_sample]
identical(colnames(exprSet),rownames(meta))
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
x=t(exprSet)
#y = as.numeric(meta$OS)
y=as.matrix(data.frame(time = meta$OS.time,
                       status = meta$OS))
library(glmnet)
choose_gene_min <- c()
start_num <- 0
fam <- "cox"
while(length(choose_gene_min) < 4){
  if (start_num > 20) {
    break;
  }
  alpha = 1
  model_lasso <- glmnet(x, y, family = fam, nlambda = 1000, alpha=alpha)
  cv_fit <- cv.glmnet(x=x, y=y, family = fam, nlambda = 1000,alpha = alpha)
  pdf("../3.Analysis/3.7.lasso.features.pdf", width = 5, height = 5)
  plot(cv_fit)
  plot(cv_fit$glmnet.fit,xvar ="lambda")
  abline(v=log(c(cv_fit$lambda.min,cv_fit$lambda.1se)),lty=2)
  dev.off()
  model_lasso_min <- glmnet(x=x, y=y, family = fam, alpha = alpha, lambda=cv_fit$lambda.min)
  model_lasso_1se <- glmnet(x=x, y=y, family = fam, alpha = alpha, lambda=cv_fit$lambda.1se)
  choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
  choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
  length(choose_gene_min)  #70
  length(choose_gene_1se)  #40
  start_num = start_num+1
}
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
#如上得到根据模型预测每个样本的生存概率的单列矩阵

select_genes <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
coefs <- as.numeric(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
cors <- c()
pvas <- c()
for (gen in select_genes) {
  corre <- cor.test(as.numeric(lasso.prob[, 1]), as.numeric(x[, which(colnames(x) == gen)]))
  cors <- c(cors, corre$estimate[[1]])
  pvas <- c(pvas, corre$p.value)
}

res_mtr <- data.frame(ID = select_genes,
                      Lasso_Coef = coefs,
                      Corr = cors,
                      Pvalue = pvas,
                      stringsAsFactors = F)

write.table(res_mtr, file = "../3.Analysis/3.7.lasso.lncRNA.coef.txt",
            sep = "\t", quote = FALSE, row.names = F, col.names = T)


re=cbind(y ,lasso.prob)
#合并样本预测值与真实值
head(re)
re=as.data.frame(re)
colnames(re)=c('time', 'event','prob_min','prob_1se')
re$event=as.factor(re$event)
library(ggpubr)
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
pdf("../3.Analysis/3.7.lasso.riskscore.boxplot.pdf", width = 5, height = 5)
p1
dev.off()
library(ROCR)
pred_min <- prediction(re[,3], re[,2])
auc_min = performance(pred_min,"auc")@y.values[[1]]
#求得AUC值
perf_min <- performance(pred_min,"tpr","fpr")
pdf("../3.Analysis/3.7.lasso.riskscore.roc.pdf", width = 5, height = 5)
plot(perf_min,colorize=FALSE, col="blue")
#绘图
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# y=x
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
# 加AUC值
dev.off()
##plot survival
meta$Value <- re$prob_min
meta$Type <- ifelse(meta$Value > median(meta$Value), "0.High", "1.Low")
write.table(meta, file = '../3.Analysis/3.7.risk.score.logistic.txt',
            sep = "\t", quote = FALSE, row.names = F, col.names = T)

restcga <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
restcga[order(restcga$Pvalue_OS, decreasing = F), ][1, ]
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = restcga[order(restcga$Pvalue_OS, decreasing = F), ][1, 6])
pdf("../3.Analysis/3.7.lasso.riskscore.survival.pdf", width = 6, height = 5)
surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = restcga[order(restcga$Pvalue_OS, decreasing = F), ][1, 6])
dev.off()

##### lasso cox model evaluated #####
lncscoef <- fread("../3.Analysis/3.7.lasso.lncRNA.coef.txt", header = T, stringsAsFactors = F, data.table = F)

#CHCC
gsesam <- read.delim("../0.data/CHCC.os.txt", header = T, row.names = 1)
gsemtr1 <- read.delim("../0.data/CHCC.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr1))
gsesam <- gsesam[common_sample, ]
gsemtr1 <- gsemtr1[, common_sample]

#gsemtr <- gsemtr[lncscoef$ID, ]
gsemtr <- gsemtr1[sample(1:40000, 9), ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$OS.time,
                   OS = gsesam$OS,
                   CHCC = colSums(exp_sl))
meta <- na.omit(meta)
write.table(meta, "../3.Analysis/3.7.lasso.sample.score.CHCC.txt", row.names = F, col.names = T, sep = "\t", quote = F)

res1 <- surv_bestcut(as.data.frame(t(meta)), "CHCC", meta, num = 20)
res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6]
res <- surv_bestcut(as.data.frame(t(meta)), "CHCC", meta, num = res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6])
pdf("../3.Analysis/3.7.lasso.riskscore.survival.CHCC.pdf", width = 6, height = 5)
surv_bestcut(as.data.frame(t(meta)), "CHCC", meta, num = res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6])
dev.off()

#GSE14520
gsesam <- read.delim("../0.data/GSE14520.os.txt", header = T, row.names = 1)
gsemtr1 <- read.delim("../0.data/GSE14520.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr1))
gsesam <- gsesam[common_sample, ]
gsemtr1 <- gsemtr1[, common_sample]

#gsemtr <- gsemtr1[lncscoef$ID, ]
gsemtr <- gsemtr1[sample(1:40000, 9), ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$OS.time,
                   OS = gsesam$OS,
                   GSE14520 = colSums(exp_sl))
meta <- na.omit(meta)
write.table(meta, "../3.Analysis/3.7.lasso.sample.score.GSE14520.txt", row.names = F, col.names = T, sep = "\t", quote = F)

res1 <- surv_bestcut(as.data.frame(t(meta)), "GSE14520", meta, num = 20)
res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6]
res <- surv_bestcut(as.data.frame(t(meta)), "GSE14520", meta, num = res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6])
pdf("../3.Analysis/3.7.lasso.riskscore.survival.GSE14520.pdf", width = 6, height = 5)
surv_bestcut(as.data.frame(t(meta)), "GSE14520", meta, num = res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6])
dev.off()

#ICGC
gsesam <- read.delim("../0.data/ICGC.os.txt", header = T, row.names = 1)
gsemtr1 <- read.delim("../0.data/ICGC.matrix", header = T, row.names = 1)
gsemtr1 <- as.data.frame(scale(t(scale(t(gsemtr1)))))
common_sample = intersect(rownames(gsesam), colnames(gsemtr1))
gsesam <- gsesam[common_sample, ]
gsemtr1 <- gsemtr1[, common_sample]

#gsemtr <- gsemtr1[lncscoef$ID, ]
gsemtr <- gsemtr1[sample(1:40000, 9), ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$OS.time,
                   OS = gsesam$OS,
                   ICGC = colSums(exp_sl))
meta <- na.omit(meta)
write.table(meta, "../3.Analysis/3.7.lasso.sample.score.ICGC.txt", row.names = F, col.names = T, sep = "\t", quote = F)

res1 <- surv_bestcut(as.data.frame(t(meta)), "ICGC", meta, num = 20)
res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6]
res <- surv_bestcut(as.data.frame(t(meta)), "ICGC", meta, num = res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6])
pdf("../3.Analysis/3.7.lasso.riskscore.survival.ICGC.pdf", width = 6, height = 5)
surv_bestcut(as.data.frame(t(meta)), "ICGC", meta, num = res1[order(res1$Pvalue_OS, decreasing = F), ][1, 6])
dev.off()


###########3.7.3.score and timwROC###########
library(timeROC)
library(survival)
library(ggplot2)
tr <- read.delim("../3.Analysis/3.7.risk.score.logistic.txt", row.names = 1)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res1 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res1$Class <- "TCGA"
# ROC.b的设置与其他两个相比，marker设置不同，可以看到tr$b前后多一个负号(“-”)
tr <- fread("../3.Analysis/3.7.lasso.sample.score.CHCC.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$CHCC,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res2 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res2$Class <- "CHCC"

tr <- fread("../3.Analysis/3.7.lasso.sample.score.GSE14520.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$GSE14520,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res3 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res3$Class <- "GSE14520"

tr <- fread("../3.Analysis/3.7.lasso.sample.score.ICGC.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$ICGC,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res4 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res4$Class <- "ICGC"

res <- rbind(res1, res2, res3, res4)

p <- ggplot(res, aes(x=time, y=auc, color = Class)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 0.5, linetype = 2)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_y_continuous(breaks = seq(0.4, 1.0, 0.1), expand = c(0,0.1))+
  scale_x_continuous(breaks = seq(0, 365*5, 365*1))+
  theme(axis.title.y = element_text(size = 15, colour="black", angle=90),
        axis.title.x = element_text(size = 15, colour="black", angle=0))+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 0),
        axis.text.y = element_text(size = 10, colour = "black", angle = 0))+
  xlab("time (days)") + ylab("AUC")

p


ggsave("../3.Analysis/3.7.timeROC.pdf", p, width = 5, height = 3)


###########3.7.4.score and calibration and nomogram###########
library(rms)  ##绘制列线图
library(survival)  ##生存分析包
scores=read.table("../3.Analysis/3.7.risk.score.logistic.txt",header = T)
cliniinfo <- fread("/home/yukai/data/Datasets/TCGA20201022/phenotype/TCGA-LIHC.GDC_phenotype.tsv", header = T, stringsAsFactors = F)
luadph <- data.frame(ID = cliniinfo$submitter_id.samples,
                     Age = cliniinfo$age_at_initial_pathologic_diagnosis,
                     Gender = cliniinfo$gender.demographic,
                     TNMs = cliniinfo$tumor_stage.diagnoses,
                     stringsAsFactors = F)

luadph$Gender <- ifelse(luadph$Gender == "", NA, luadph$Gender)
luadph$TNMs <- ifelse(luadph$TNMs == "", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "not reported", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ia", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ib", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iia", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iib", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiia", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiib", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiic", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iva", "stage iv", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ivb", "stage iv", luadph$TNMs)
res <- merge(scores, luadph, by.x = "sample", by.y = "ID")

d <- na.omit(res)
###添加变量标签，在列线图上展示分类标签，作图用到
d$Age <- as.numeric(d$Age)
d$Survival <- d$OS
d$time <- d$OS.time
dd<-datadist(d) #设置工作环境变量，将数据整合
options(datadist='dd') #设置工作环境变量，将数据整合
##
coxm <- cph(Surv(time,Survival)~Age+Gender+TNMs+Value,x=T,y=T,data=d,surv=T)
###绘制cox回归生存概率的nomogram图
## 构建Nomo图的对象只能是rms保重d额cph()函数
survival = Survival(coxm)
survival1 = function(x)survival(1*365,x)
survival3 = function(x)survival(3*365,x)
survival5 = function(x)survival(5*365,x)
nom <- nomogram(coxm,fun=list(survival1,survival3,survival5), ##算出不同时间节点生存率值，显示在列线图上
                funlabel = c('1-year probability',
                             '3-year probability',
                             '5-year probability'),
                lp=F,
                fun.at=c('0.9','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1'))
par(mar=c(2,5,3,2),cex=0.8)##mar 图形空白边界  cex 文本和符号大小
pdf("../3.Analysis/3.7.risk.score.normagram.pdf", width = 5, height = 5)
plot(nom,xfrac=0.6)
dev.off()
###这里演示，采用age sex 和ph.ecog 构建coxph模型
##构建校准曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
pdf("../3.Analysis/3.7.risk.score.calibration.pdf", width = 5, height = 5)
coxm_1 <- cph(Surv(time,Survival)~Age+Gender+TNMs+Value,data=d,surv=T,x=T,y=T,time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=50,B=1000)
##绘制1年生存期校准曲线
par(mar=c(7,4,4,3),cex=1.0)
plot(cal_1,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-year OS',#便签
     ylab='Actual 1-year OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)))#设置一个颜色
##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_2 <- cph(Surv(time,Survival)~Age+Gender+TNMs+Value,data=d,surv=T,x=T,y=T,time.inc = 3*365)
cal_2<-calibrate(coxm_2,u=3*365,cmethod='KM',m=50,B=1000)
plot(cal_2,lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 3-year OS',#便签
     ylab='Actual 3-year OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)))#设置一个颜色
##绘制5年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(Surv(time,Survival)~Age+Gender+TNMs+Value,data=d,surv=T,x=T,y=T,time.inc = 5*365)
cal_3<-calibrate(coxm_3,u=5*365,cmethod='KM',m=50,B=1000)
plot(cal_3,lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5-year OS',#便签
     ylab='Actual 5-year OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)))#设置一个颜色
dev.off()

####################3.8.CRC Liver metastasis####################
library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")

pbmc <- readRDS("../2.Annote/MergeAnno/Merged.pca.major.minor.rmoutlier.rds")
pbmc$CompClass <- as.character(pbmc$TopCluster)
pbmc$CompClass <- ifelse(pbmc$Celltype == "CD4T", "T cell CD4", pbmc$CompClass)
pbmc$CompClass <- ifelse(pbmc$Celltype == "CD8T", "T cell CD8", pbmc$CompClass)

###########3.8.1.cell count###########
out.prefix <- "../3.Analysis/3.8.cellcount"
meta.tb <- pbmc@meta.data
meta.tb <- meta.tb[meta.tb$CompClass != "Epithelial cell", ]

##CRC
cellInfo.tb <- meta.tb[meta.tb$Tissue %in% c("PC", "CRC"), ]
OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$TopCluster,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.TopCluster.colon",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.cellcount.STARTRAC.dist.TopCluster.colon.OR.dist.txt", 
       row.names = T, quote = F)

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$CellClass,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.CellClass.colon",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.cellcount.STARTRAC.dist.CellClass.colon.OR.dist.txt", 
       row.names = T, quote = F)

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$MinorCluster,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.MinorCluster.colon",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.cellcount.STARTRAC.dist.MinorCluster.colon.OR.dist.txt", 
       row.names = T, quote = F)

##liver
cellInfo.tb <- meta.tb[meta.tb$Tissue %in% c("LM", "PL"), ]

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$TopCluster,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.TopCluster.liver",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.cellcount.STARTRAC.dist.TopCluster.liver.OR.dist.txt", 
       row.names = T, quote = F)

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$CellClass,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.CellClass.liver",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.cellcount.STARTRAC.dist.CellClass.liver.OR.dist.txt", 
       row.names = T, quote = F)

OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$MinorCluster,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.MinorCluster.liver",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.cellcount.STARTRAC.dist.MinorCluster.liver.OR.dist.txt", 
       row.names = T, quote = F)



###########3.8.2.cell cluster###########
library(SingleCellExperiment)
##visualiza Mast
sub_sce <- readRDS("../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.pca.minor.rds")
sub_sce <- subset(sub_sce, subset = Tissue != "NL")
sce<-as.SingleCellExperiment(sub_sce)
reducedDim(sce, "umap")<-Embeddings(object = sub_sce, reduction = "umap")
reducedDim(sce, "tsne")<-Embeddings(object = sub_sce, reduction = "tsne")

Cluster<-levels(Idents(sub_sce))
g.colSet <- readRDS("/home/yukai/projects/sclearn/NPC2023fromMei/0.data/colSet.list.rds")
g.colSet1<-g.colSet$Set8[1:length(Cluster)]
names(g.colSet1)<-Cluster
g.colSet1<-list(g.colSet1)

p1 <- ssc.plot.tsne(sce, columns = "MinorCluster", 
                    reduced.name = "umap",
                    colSet=g.colSet1,
                    size=0.5,
                    label=2,
                    par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                    show.legend = T, 
                    par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                    base_aspect_ratio = 1.35)

ggsave("../3.Analysis/3.8.cellcount.sel_Cluster_merge_umap.MinorCluster.pdf",p1,width = 7,height = 5 )

p2 <- ssc.plot.tsne(sce, columns = "MinorCluster", 
                    splitBy = "Tissue",
                    reduced.name = "umap",
                    colSet=g.colSet1,
                    size=0.5,
                    label=2,
                    par.repel = list(force = 1,bg.color="white",bg.r=0.15),
                    show.legend = T, 
                    par.geneOnTSNE=list(scales="free",pt.order="random",pt.alpha=0.8),
                    base_aspect_ratio = 1.35)

ggsave("../3.Analysis/3.8.cellcount.sel_Cluster_merge_umap.Tissue.pdf",p2,width = 12,height = 9 )


out.prefix <- "../3.Analysis/3.8.mastcount"
meta.tb <- pbmc@meta.data
meta.tb <- meta.tb[meta.tb$CompClass == "Mast cell", ]

##CRC
cellInfo.tb <- meta.tb[meta.tb$Tissue %in% c("PC", "CRC"), ]
OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$MinorCluster,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.MinorCluster.colon",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.mastcount.STARTRAC.dist.MinorCluster.colon.OR.dist.txt", 
       row.names = T, quote = F)

##liver
cellInfo.tb <- meta.tb[meta.tb$Tissue %in% c("LM", "PL"), ]
OR.all.list <- do.tissueDist(cellInfo.tb=cellInfo.tb,
                             meta.cluster = cellInfo.tb$MinorCluster,
                             colname.patient = "Patient",
                             loc = cellInfo.tb$Tissue,
                             out.prefix=sprintf("%s.STARTRAC.dist.MinorCluster.liver",out.prefix),
                             pdf.width=4,pdf.height=6,verbose=1)

fwrite(as.data.frame(OR.all.list$OR.dist.mtx), "../3.Analysis/3.8.mastcount.STARTRAC.dist.MinorCluster.liver.OR.dist.txt", 
       row.names = T, quote = F)

sub_sce <- subset(pbmc, subset = TopCluster %in% c("Mast cell"))
Idents(sub_sce) <- "MinorCluster"
selgene <- c("IL32", "S100A8", "KRT86", "TPSAB1", "ITGA2B", "CPA3", "HSPA1A")
pdf("../3.Analysis/3.8.sub_sce.marker.vlnplot.exp.pdf",width = 6,height = 4)
VlnPlot(sub_sce, selgene, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()


###########3.8.3.monocle###########
library(monocle3)
library(Seurat)
library(tidyverse)
library(patchwork)

dir.create("../3.Analysis/Monocle4Mast")
combined <- readRDS("../2.Annote/2.7.Mast/Merged.pca.major.Mastcells.pca.minor.rds")
Idents(combined) <- 'MinorCluster'
##创建CDS对象并预处理数据
pbmc11 <- combined
data <- GetAssayData(pbmc11, assay = 'RNA', slot = 'counts')
cell_metadata <- pbmc11@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#UMAP
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "MinorCluster")

#import umap from seurat
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(pbmc, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]

cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "MinorCluster")
p <- p1|p2
p
ggsave("../3.Analysis/Monocle4Mast/reduction.umap.pdf", p, width = 6, height = 3)

#visualize some genes
ciliated_genes <- c("ACTA2")
plot_cells(cds, 
           genes = ciliated_genes,
           label_cell_groups =FALSE,
           show_trajectory_graph = F)

#cluster your cells for trajectory
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

cds <- learn_graph(cds)
p = plot_cells(cds, color_cells_by = "MinorCluster",
               label_groups_by_cluster =FALSE,
               show_trajectory_graph = TRUE)
ggsave("../3.Analysis/Monocle4Mast/reduction.umap.trajectory.pdf", p, width = 5, height = 4)

#pseudotime
cds <- order_cells(cds)
p = plot_cells(cds, color_cells_by = "pseudotime",
               label_groups_by_cluster =FALSE,
               show_trajectory_graph = TRUE)
p
ggsave("../3.Analysis/Monocle4Mast/reduction.umap.trajectory.pseudotime.pdf", p, width = 5, height = 4)



###########3.8.4.pathway###########
dir.create("../3.Analysis/3.8.scmetabolism")
dir.create("../3.Analysis/3.8.scmetabolism/sepfiles")
pbmc <- subset(pbmc, subset = Tissue != "NL")
pbmc.down <- subset(x = pbmc, downsample = 1000)
library(scMetabolism)
library(ggplot2)
library(rsvd)

allclass <- as.character(unique(pbmc.down$CellClass))
for (subcells in allclass) {
  print(paste("Starting computing ", subcells, " !!!!!!!!!!!", sep = ""))
  tmppbmc <- subset(pbmc.down, subset = CellClass == subcells)
  print(paste(ncol(tmppbmc), " Cells", sep = ""))
  countexp.Seurat<-sc.metabolism.Seurat(obj = tmppbmc, method = "AUCell", imputation = F, ncores = 20, metabolism.type = "KEGG")
  saveRDS(countexp.Seurat, file = paste("../3.Analysis/3.8.scmetabolism/sepfiles/", subcells, ".scmetabolism.rds", sep = ""))
}

metacomb <- data.frame()
for (subcells in allclass) {
  scmeta <- readRDS(paste("../3.Analysis/3.8.scmetabolism/sepfiles/", subcells, ".scmetabolism.rds", sep = ""))
  metamtr <- scmeta@assays$METABOLISM$score
  metamtr <- as.data.frame(t(metamtr))
  metacomb <- rbind(metacomb, metamtr)
}

sams <- pbmc.down@meta.data
rownames(sams) <- gsub("-", ".", rownames(sams))
allcomb <- metacomb[rownames(sams), ]
scmetabolism.obj <- CreateSeuratObject(t(allcomb), meta.data = sams)
saveRDS(scmetabolism.obj, file = "../3.Analysis/3.8.scmetabolism/combined.scmetabolism.rds")

##de pathways
scmetabolism.obj$SelCluster <- scmetabolism.obj$CompClass
scmetabolism.obj$SelCluster <- ifelse(scmetabolism.obj$MinorCluster %in% c("Mast_Mast_IL32", 
                                                                         "Mast_Mast_S100A8", 
                                                                         "Mast_Mast_KRT86"), "Mast_1", scmetabolism.obj$SelCluster)
scmetabolism.obj$SelCluster <- ifelse(scmetabolism.obj$MinorCluster %in% c("Mast_Mast_TPSAB1", 
                                                                         "Mast_Mast_ITGA2B", 
                                                                         "Mast_Mast_CPA3",
                                                                         "Mast_Mast_HSPA1A"), "Mast_2", scmetabolism.obj$SelCluster)
Idents(scmetabolism.obj) <- "SelCluster"
markers <- FindAllMarkers(scmetabolism.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.)
write.table(markers,"../3.Analysis/3.8.scmetabolism/combined.scmetabolism.de.metas.txt",sep="\t",row.names=T,col.names=T,quote=F)


##show pathways
##
selpathways <- c("Inositol phosphate metabolism", 
                 "Phenylalanine, tyrosine and tryptophan biosynthesis", 
                 "Glycosaminoglycan biosynthesis - keratan sulfate")

tmppbmc <- scmetabolism.obj
pdf("../3.Analysis/3.8.scmetabolism/selpath.Vlnplot.immunecells.pdf",width = 10,height = 5)
VlnPlot(tmppbmc, selpathways, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()

tmppbmc <- subset(scmetabolism.obj, subset = TopCluster %in% c("Mast cell"))
Idents(tmppbmc) <- tmppbmc$MinorCluster
pdf("../3.Analysis/3.8.scmetabolism/selpath.Vlnplot.mastcells.minor.pdf",width = 10,height = 5)
VlnPlot(tmppbmc, selpathways, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()

Idents(tmppbmc) <- tmppbmc$Tissue
pdf("../3.Analysis/3.8.scmetabolism/selpath.Vlnplot.mastcells.tissue.pdf",width = 10,height = 5)
VlnPlot(tmppbmc, selpathways, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
dev.off()

###########3.8.5.tcga###########
library(GSVA)
library(msigdbr)
library(GSEABase)
tcgaMtr <- fread("/home/yukai/data/TCGA_flowchart/TCGA_new_VERSION/htseq_fpkm/gene/TCGA-LIHC.htseq_fpkm.tsv.cv.txt", header = T, sep = '\t', check.names = F, stringsAsFactors = F, data.table = F)
rownames(tcgaMtr) <- tcgaMtr$Ensembl_ID
tcgaMtr <- tcgaMtr[, -1]
colData <- data.frame(Type = ifelse(substr(colnames(tcgaMtr), 14, 15) < 10, 'Tumor', 'Normal'),
                      ID = colnames(tcgaMtr))
colData$ID <- as.character(colData$ID)
colData <- colData[colData$Type == 'Tumor', ]
tcgaMtr <- tcgaMtr[, colData$ID]
tcgaMtr <- scale(t(scale(t(tcgaMtr))))

##signature
markers <- fread("../2.Annote/MergeAnno/Combined_Cluster_markers.txt", header = F, stringsAsFactors = F, data.table = F)
markers <- markers[markers$V7 == "Mast cell", ]
markers <- markers[order(markers$V6), ]

top10genes <- markers$V8[1:50]
write.table(t(c("Markers", "www", top10genes)), "../3.Analysis/3.8.signature.gmt", row.names = F, col.names = F, sep = "\t", quote = F)

immuneCirc <- getGmt("../3.Analysis/3.8.signature.gmt")
tcgaMtr <- as.matrix(tcgaMtr)
mast_signature <- gsva(tcgaMtr, immuneCirc, method = "ssgsea")
mast_signature <- as.data.frame(t(mast_signature))
mast_signature$ColID <- rownames(mast_signature)

##xcells
xcell <- fread("../3.Analysis/3.4.tcga.scores.immunecell.xcell.txt", header = T, stringsAsFactors = F, data.table = F)
rownames(xcell) <- xcell$V1
xcell <- xcell[, -1]
xcell <- as.data.frame(t(xcell))
##KEGG pathway
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(plyr)
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_description, gene_symbol)
kegggmt <- split(kegg$gene_symbol, kegg$gs_description)
keggpath <- gsva(as.matrix(tcgaMtr), kegggmt, kcdf="Gaussian",method = "gsva",parallel.sz=10)
keggpath <- as.data.frame(t(keggpath))
keggpath$ColID <- rownames(keggpath)

##metasis
metas <- fread("/home/yukai/work/gc_data/gc_combine_analysis/20210118_datas/0.database/EMT.Meta.Metabolic.oncoTSG.list",
               header = F, stringsAsFactors = F, data.table = F)
metas <- metas[metas$V2 == "Meta", ]
metagmt <- split(metas$V3, metas$V2)
metapath <- gsva(as.matrix(tcgaMtr), metagmt, kcdf="Gaussian",method = "gsva",parallel.sz=10)
metapath <- as.data.frame(t(metapath))
metapath$ColID <- rownames(metapath)
##merge
commonsams <- intersect(rownames(mast_signature), rownames(xcell))
mergedf <- data.frame(ID = commonsams,
                      Metas = metapath[commonsams, ]$Meta,
                      Mast = xcell[commonsams, ]$`Mast cells`,
                      Signature = mast_signature[commonsams, ]$Markers,
                      IPM = keggpath[commonsams, ]$`Inositol phosphate metabolism`,
                      PM = keggpath[commonsams, ]$`Phenylalanine metabolism`,
                      GB = keggpath[commonsams, ]$`Glycosaminoglycan biosynthesis - keratan sulfate`)

corre <- cor.test(mergedf$Metas,mergedf$Mast,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p1 <- ggplot(mergedf, aes(x = Metas, y = Mast))+
  geom_point()+
  ggtitle(plottitle)+
  geom_smooth(method="lm",color="#1a9641") + 
  xlab("Metastasis Index")+
  ylab("Estimated Mast cell fraction")+
  theme_classic2()
p1

corre <- cor.test(mergedf$Metas,mergedf$Signature,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p2 <- ggplot(mergedf, aes(x = Metas, y = Signature))+
  geom_point()+
  ggtitle(plottitle)+
  geom_smooth(method="lm",color="#1a9641") + 
  xlab("Metastasis Index")+
  ylab("Inositol phosphate metabolism")+
  theme_classic2()
p2

corre <- cor.test(mergedf$Metas,mergedf$IPM,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p3 <- ggplot(mergedf, aes(x = Metas, y = IPM))+
  geom_point()+
  ggtitle(plottitle)+
  geom_smooth(method="lm",color="#1a9641") + 
  xlab("Metastasis Index")+
  ylab("Estimated Mast cell fraction")+
  theme_classic2()
p3

corre <- cor.test(mergedf$Metas,mergedf$PM,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p4 <- ggplot(mergedf, aes(x = Metas, y = PM))+
  geom_point()+
  ggtitle(plottitle)+
  geom_smooth(method="lm",color="#1a9641") + 
  xlab("Metastasis Index")+
  ylab("Phenylalanine metabolism")+
  theme_classic2()
p4

corre <- cor.test(mergedf$Metas,mergedf$GB,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p5 <- ggplot(mergedf, aes(x = Metas, y = GB))+
  geom_point()+
  ggtitle(plottitle)+
  geom_smooth(method="lm",color="#1a9641") + 
  xlab("Metastasis Index")+
  ylab("Glycosaminoglycan biosynthesis")+
  theme_classic2()
p5
pdf("../3.Analysis/3.8.mastcount.sig.tcga.cor.pdf", width = 10, height = 6)
grid.arrange(p1,p2,p3,p4,p5, nrow = 2)
dev.off()


