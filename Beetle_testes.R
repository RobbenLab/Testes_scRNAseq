#########################################################
###########     Setup         ###########################
#########################################################
#load all required libraries
require(Seurat)
require(venn)
require(monocle3)
library(SeuratWrappers)
require(org.Dm.eg.db)
require(AnnotationDbi)
require(rtracklayer)
require(dplyr)
require(tibble)
require(pheatmap)
library(scales)
require(Matrix)
require(escape)
library(dittoSeq)
require(plyr)
require(data.table)
require(ggplot2)
require(sparsepca)
require(umap)
require(qlcMatrix)
# devtools::install_github("kieranrcampbell/ouija")
require(ouija)

#First set wd
wd <- "/home/robbenm/LuberLab/Beetle/"
setwd(wd)
#Next load in annotation file
anno <- read.csv(file = "./Annotation/tcas3/annotation.tsv",header = T,sep = "\t")

#read in data
#Put the output of the cellranger run into ./CellRanger_out/feature
tri <- CreateSeuratObject(Read10X(data.dir = "./CellRanger_out/feature"),min.cells = 3, min.features = 100) 
head(tri@meta.data)
head(GetAssayData(tri))

#########################################################
########### Data Exploration and Filter  ################
#########################################################

#First find the mitochondrial genes using annotation
tcasmt <- subset(anno, D.Chrom == "mitochondrion MT")$tcas_name
tri[["percent.mt"]] <- PercentageFeatureSet(tri, features = tcasmt[1]) #can't determine mitochondrial because not human, we would need to look up the mitochondrial genes in tri which is too much work
VlnPlot(tri, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(tri, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tri, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#Now we can subset for cells that express features between 100-3000 count and less than 5% mitochondrial
tri <- subset(tri, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 0.05)
nrow(tri@meta.data) 
#We removed about 400 doublet and high mt expressing cells

#Now we will normalize and scale the data
tri <- NormalizeData(tri)
tri <- FindVariableFeatures(tri, selection.method = "vst")
tri <- ScaleData(tri, features = rownames(tri))

# We want to filter out cells that are most likely not germline or cyst related
epi = anno[anno$D.Symbol %in% "shg",]$tcas_name
mhc <- anno[anno$D.Symbol %in% "Mhc",]$tcas_name
mhc <- anno[anno$D.Symbol %in% "Prm",]$tcas_name
ppn <- anno[anno$D.Symbol %in% c("lz","NimC1","Hml","atilla","odd","He"),]$tcas_name 
neur <- c(anno[anno$D.Symbol %in% c("para","mmd"),]$tcas_name,"TcasGA2-TC030021")
pig <- c(anno[anno$D.Symbol %in% c("Sox100B","ems"),]$tcas_name,"TcasGA2-TC006235","TcasGA2-TC014650","TcasGA2-TC033720","TcasGA2-TC009280")
fb <- anno[anno$D.Symbol %in% c("Lsp2","Fbp1","Fbp2"),]$tcas_name
cy <- anno[anno$D.Symbol %in% c("tj","stg","esg"),]$tcas_name
germ <- c("TcasGA2-TC009035","TcasGA2-TC006703","TcasGA2-TC011729")

#Plot germline cells
venn(list(Somatic = which(as.numeric(rowSums(FetchData(tri, c(mhc,neur,fb,pig,ppn,epi)))) > 0.1),
          Germ_cyst = rownames(GetAssayData(tri))))

#Remove non-germ cells
somcells <- unique(c(which(as.numeric(rowSums(FetchData(tri, mhc))) > 0.1),which(as.numeric(rowSums(FetchData(tri, neur))) > 0.1),which(as.numeric(rowSums(FetchData(tri, fb))) > 0.1),which(as.numeric(rowSums(FetchData(tri, cy))) > 0.1),which(as.numeric(rowSums(FetchData(tri, pig))) > 0.1),which(as.numeric(rowSums(FetchData(tri, ppn))) > 0.1),which(as.numeric(rowSums(FetchData(tri, epi))) > 0.1)))
remove <- rownames(tri@meta.data[somcells,])
trigerm <- tri[,!colnames(tri) %in% remove]
trigerm
trigerm <- FindVariableFeatures(object = trigerm)

#Calculate number of cells in which a certain gene is expressed
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}
calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

num_exp <- PrctCellExpringGene(tri400germ,genes = rownames(tri400germ),group.by = "all")
over10 <- num_exp[num_exp$Cell_proportion > .1,]

#We can also predict stage of mitotic activity

g2m <- blast_ann_name[blast_ann_name$D.Symbol %in% c("cana","aurA","CycB","fzy","twe","pbl","RanGAP","stg","sub","Su(var)205","Cdk1","Klp61F","pav","glu","msps","pigs","Nek2","HP1b","Bub1","mars","mapmodulin","LBR","CTCF","Cks85A","HP1c","Cap-D2",
                                                     "cmet","scra","BubR1","Det","vih","Dsp1",""),]$tcas_name 
s <- blast_ann_name[blast_ann_name$D.Symbol %in% c("Blm","spn-A","PCNA","CycE","RnrS","DNApol-α50","l(2)dtl","Mcm2","spel1","tos","dpa","Mcm5","Ts","Mcm6","Fen1","Cdc45","Usp1","CG15141","CG10336","PCNA2","RPA2","Caf1-105",
                                                   "CG11788","Cdc6","Slbp","Claspin","DNApol-α180","RfC4","Psf2"),]$tcas_name 

trigerm <- CellCycleScoring(trigerm,g2m.features = g2m,s.features = s)

#########################################################
###########       Marker set up          ################
#########################################################

#General markers of spermatogenesis in fly
markers <- read.table(file = "./Markers/spermmarkers.csv",header = T,sep = ",")
head(markers)
b.markers <- anno[anno$D.Symbol %in% markers$Gene.symbol,]
markers[!markers$Gene.symbol %in% unique(b.markers$D.Symbol),]
#Only 18 markers from the original 43 exist in beetle
DoHeatmap(tri, features = b.markers)

#B2 Tubulin, Rad50 and Enolase can be used to stage beetle spermatogenesis
venn(list(B2tubulin = which(as.numeric(rowSums(FetchData(tri,"TcasGA2-TC009035" ))) > 0.1),
          Rad50 = which(as.numeric(rowSums(FetchData(tri, "TcasGA2-TC006703"))) > 0.1),
          enolase = which(as.numeric(rowSums(FetchData(tri, "TcasGA2-TC011729"))) > 0.1),
          all = rownames(GetAssayData(tri))))

heatmap(t(as.matrix(FetchData(tri,germ))),labCol = F,labRow = c("B2tubulin","Rad50","Enolase"))

#Markers indicated for germline vs cyst expression from fly data
flyup <- read.table("./Drosophila_results/Ovary/fxljbc_run2_wilcox_seurat_de_table_up.txt",header = F,sep = "\t"
flydown <- read.table("./Drosophila_results/Ovary/fxljbc_run2_wilcox_seurat_de_table_down.txt",header = F,sep = "\t")
#flyup <- read.table("./Drosophila_results/fxljbc_run1_wilcox_seurat_de_table_up.txt",header = F,sep = "\t") #alternatively use testes
#flydown <- read.table("./Drosophila_results/fxljbc_run1_wilcox_seurat_de_table_down.txt",header = F,sep = "\t")
flygermup <- anno[anno$D.Symbol %in% c(flyup[,3]),]$tcas_name 
flygermdown <- anno[anno$D.Symbol %in% c(flydown[,3]),]$tcas_name
#Check the markers expressed in over 25% of cells
highflyup <- over10[over10$Markers %in% flygermup,]
FeaturePlot(trigermplus,features = highflyup$Markers)
DotPlot(trigermplus,features = highflyup$Markers,cols = c("yellow","blue"))
#Now do down 
highflydown <- over10[over10$Markers %in% flygermdown,]
FeaturePlot(trigermplus,features = highflydown$Markers)
DotPlot(trigermplus,features = highflydown$Markers,cols = c("yellow","blue"))


#########################################################
###########     Umap and clustering       ################
#########################################################


#Now we can generate a umap and cell clusters
#First we need to generate principle components weighted to biologically significant markers
fly <- unique(c(highflydown$Markers,highflyup$Markers))
plus <- blast_ann_name[blast_ann_name$D.Symbol %in% c("spz","bru1"),]$tcas_nam
nvar <- 400
match <- rownames(GetAssayData(trigerm))[as.numeric(na.omit(match(c(germ,epi,mhc,ppn,neur,pig,cy,fb,plus,fly),rownames(GetAssayData(trigerm)))))]
cdpr <- prcomp(x = GetAssayData(trigerm)[c(rep(match,each = 10),VariableFeatures(object = trigerm)[1:nvar]),])
#Biologically significant markers are weighted 10 times against the 400 most variable genes expressed
trigerm[["pca"]] <- CreateDimReducObject(embeddings = cdpr$rotation,stdev = cdpr$sdev,key = "PCA_", assay = DefaultAssay(trigerm))
trigerm <- ProjectDim(object = trigerm)
nPC = 8;
#We use 8 dimensions of Principal components which show the greatest variance
triflygerm <- RunUMAP(object = trigerm, reduction = "pca", dims = 1:nPC);
#Plot out the umap
DimPlot(triflygerm,reduction = "umap")
#We demonstrate the expression of markers known to stage beetle spermatogenesis
FeaturePlot(triflygerm,features = germ,cols = c("yellow","blue"))

#Now we can cluster
triflygerm <- FindNeighbors(triflygerm,reduction = "pca",dims = c(1:8))
triflygerm <- FindClusters(triflygerm,resolution = c(0.2,seq(0.25,4,0.25)))
sapply(grep("res",colnames(triflygerm@meta.data),value = TRUE),
       function(x) length(unique(triflygerm@meta.data[,x])))

#We choose the resolution 0.5 because it shows the most granularity between cell types
DimPlot(triflygerm,group.by = "RNA_snn_res.0.5",label = T)
Idents(triflygerm) <- "RNA_snn_res.0.5"

#Now we can rename the cells to match what we have found 
cell_types <- c("Spermatogonia 2",#0
                "Secondary Spermatocyte 2",#1
                "Secondary Spermatocyte 1",#2
                "Spermatogonia 1",#3
                "GSC",#4
                "Primary Spermatocyte 1",#5
                "Early Spermatid 1",#6
                "Early Spermatid 2",#7
                "Primary Spermatoocyte 4",#8
                "Spermatozoa 1",#9
                "Late Spermatid 1",#10
                "Primary Spermatocyte 3",#11
                "Spermatozoa 2",#12
                "Late Spermatid 2",#13
                "Spermatogonia 3",#14
                "Primary Spermatocyte 2"#15
)
#Cell types were named based on expression of staging markers
names(cell_types) <- levels(triflygerm)
triflygerm <- RenameIdents(triflygerm, cell_types)
cols <- c(4,3,0,14,5,15,11,8,2,1,6,7,10,13,9,12) #Stage based ordering for columns
levels(triflygerm) <- cell_types[cols+1]

DimPlot(triflygerm,reduction = "umap", label = T)
# save(triflygerm,file = "./Clustering/Triflygerm.Rdata")


#########################################################
###########     Find Markers             ################
#########################################################

#We can use the following code to predict markers based on differential expression between group and all
tfgmarker <- FindAllMarkers(triflygerm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tfgmarkeranno <- cbind(tfgmarker,anno[tfgmarker$gene,])
lapply(split(tfgmarkeranno,tfgmarkeranno$cluster),function(x) head(x[,c(2,3,4,5,10,11,19)]))

write.table(tfgmarkeranno,file = "./Markers/allmarkers.tsv",sep = "\t",quote = F)

#########################################################
##############     Ouija Psuedotime      ################
#########################################################


#First get a sample of 20 cells per cluster
ls <- lapply(as.list(cell_types), function(x) which(as.character(Idents(triflygerm)) == x))
sample <- unlist(lapply(ls,function(x) x[sample(seq_len(length(x)), 20)]))
germ_mat <- t(as.matrix(GetAssayData(triflygerm)[germ,sample]))
colnames(germ_mat) <- c("Tubulin","Rad50","Enolase")
oui <- ouija(germ_mat,iter = 500)
#With that oui should have everything
print(oui)
plot_diagnostics(oui)
plot_expression(oui,ncol = 1)
plot_switch_times(oui)
cmo <- consistency_matrix(oui)
plot_consistency(oui)
pseudotimes <- map_pseudotime(oui)

test.df <- data.frame(germ_mat,pseudotimes,Cluster = rep(cell_types,each = 20))
long.df <- test.df %>% gather(Gene,Expression,-pseudotimes,-Cluster)
long.df$Cluster <- factor(long.df$Cluster,levels = cell_types[cols+1])

ggplot(long.df,aes(x = pseudotimes,y = Expression,color = Cluster)) +
  geom_point() +
  facet_wrap(~Gene)


ggplot(long.df,aes(x = Cluster,y = pseudotimes)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tfgsub <- triflygerm[,sample]
tfgsub[["Psuedotime"]] <- pseudotimes
FeaturePlot(tfgsub,features = "Psuedotime",reduction = "umap")

#########################################################
###########     Monocle trajectory       ################
#########################################################


#We will also calculate monocle trajectory based on gene expression among clusters
cds <- as.cell_data_set(triflygerm)
rowData(cds)$gene_short_name <- rownames(GetAssayData(triflygerm))
cds
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(triflygerm[["RNA"]])
cds <- preprocess_cds(cds, num_dim = 100)
cds <- cluster_cells(cds, reduction_method = "UMAP",resolution=1e-2)
# plot_cells(cds,color_cells_by = "partition")
# cds <- learn_graph(cds)
partitions(cds)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
#This shows the monocle trajectory for clusters

#We will also generate co-expression modules in monocle
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=Idents(triflygerm))
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# agg_mat <- agg_mat[,cols+1]
#Plot out clusters as a heatmep
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=F,
                   scale="column", clustering_method="ward.D2",
                   fontsize=12,angle_col = 45)


#########################################################
###########     Marker expression        ################
#########################################################

#We can also visualize marker expression by cluster
#First lets look at expression of all markers identified in findallmarkers
# DoHeatmap( #Takes a long time to run
#   triflygerm,
#   features = unique(rownames(tgpmarkeranno)),
#   cells = NULL,
#   group.by = "ident",
#   group.bar = TRUE,
#   group.colors = NULL,
#   disp.min = -2.5,
#   disp.max = NULL,
#   slot = "scale.data",
#   assay = NULL,
#   label = TRUE,
#   size = 5.5,
#   hjust = 0,
#   angle = 45,
#   raster = TRUE,
#   draw.lines = TRUE,
#   lines.width = NULL,
#   group.bar.height = 0.02,
#   combine = TRUE
# )

exp <- AverageExpression(
  triflygerm,
  assays = NULL,
  features = unique(rownames(tfgmarkeranno)),
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)
exp <- exp$RNA

pheatmap(as.matrix(exp),cluster_rows=TRUE, cluster_cols=F,
         scale="column", show_rownames = F,clustering_method="ward.D2",
         fontsize=12,angle_col = 45)

#Next let's plot out markers that we identified using monocle modules

exp <- AverageExpression( #First get expression by cluster
  triflygerm,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)
exp <- exp$RNA

GSC_markers <- gene_module_df[gene_module_df$module == 10,]
GSC_markers <- within(GSC_markers, Symbol <- anno[GSC_markers$id,]$D.Symbol)
Spermatogonia_markers <- gene_module_df[gene_module_df$module %in% c(9,7),]#9,7
Spermatogonia_markers <- within(Spermatogonia_markers, Symbol <- anno[Spermatogonia_markers$id,]$D.Symbol)
Primary_markers <- gene_module_df[gene_module_df$module %in% c(5),]#5,13,14
Primary_markers <- within(Primary_markers, Symbol <- anno[Primary_markers$id,]$D.Symbol)
Secondary_markers <- gene_module_df[gene_module_df$module %in% c(18,3,2),]#18,3,2
Secondary_markers <- within(Secondary_markers, Symbol <- anno[Secondary_markers$id,]$D.Symbol)
Spermatid_markers <- gene_module_df[gene_module_df$module %in% c(1,8,12),]#8,12,1
Spermatid_markers <- within(Spermatid_markers, Symbol <- anno[Spermatid_markers$id,]$D.Symbol)
Spz_markers <- gene_module_df[gene_module_df$module == 6,] #maybe 3 and 2?
Spz_markers <- within(Spz_markers, Symbol <- anno[Spz_markers$id,]$D.Symbol)

gm <- GSC_markers[GSC_markers$Symbol %in% c("CG17377","His1:CG33801","kto","ctp","Mical","CLIP-190"),]
sgm <- Spermatogonia_markers[Spermatogonia_markers$Symbol %in% c("CG31792","CG30356","tn","S-Lap7","BBP1","Pkn","Jupiter","S-Lap1","NA","Topors","CG5903","Sbxn1-3"),]
pm <- rbind(Primary_markers[Primary_markers$Symbol %in% c("unc-13","CG1786","CG2127","Udh","CG14893","CG14331","Moe","knon","LManIII","Gbs-70E"),], Primary_markers[Primary_markers$id == "TcasGA2-TC032252",])
sm <- Secondary_markers[Secondary_markers$Symbol %in% c("CG15203","CG31676","EbpIII","yellow-c","mgl","fabp","GILT1","PHGPx","RPs20","Fer2LCH","RpLP2","Fer1HCH"),]
spm <- Spermatid_markers[Spermatid_markers$Symbol %in% c("RpLP0","CG9822","CG6770","kdn","CG9826","CG33178","Hsp83","p23","COX5B","CG31638","COX5B","CG15449","PrptA","14-3-3epsilon","CG14342"),]

spzm <- Spz_markers[Spz_markers$Symbol %in% c("CG14762","spz"),]

exp_markers <- exp[c(gm$id,sgm$id,pm$id,sm$id,spm$id),]
names <- c(gm$Symbol,sgm$Symbol,pm$Symbol,sm$Symbol,spm$Symbol)
rowanno <- c(rep("Germline",nrow(gm)),rep("Spermatogonia",nrow(sgm)),rep("Primary Spermatocyte",nrow(pm)),rep("Secondary Spermatocyte",nrow(sm)),rep("Spermatid",nrow(spm)))
colanno <- c("Germline", rep("Spermatogonia",3),rep("Spermatocyte",6),rep("Spermatid",4),rep("Spermatozoa",2))

pheatmap(t(apply(exp_markers, 1, rescale)),labels_row = names,angle_col = 45,fontsize_row = 10,cluster_cols = F,cluster_rows = F)#,annotation_col = colanno,annotation_row = rowanno)

#Also make a heatmap of germline markers to benchmark
pheatmap(t(apply(exp[germ,], 1, rescale)),labels_row = c("B2tubulin","Rad50","Enolase"),angle_col = 45,fontsize_row = 10,cluster_cols = F,cluster_rows = F,show_colnames = F)


#We also want to make plots for expression of interesting markers
markernames <- c("spz","ProtA iso b1","ProtA iso b2","PknA","S-Lap7","CG17816","Unc-13","CG14762")
ids <- paste("TcasGA2-TC",c("031906","007670","007828","032948","001175","016377","034988","035032"),sep="")
i <- 8
FeaturePlot(triflygerm,features = ids[i],cols = c("yellow","blue")) + labs(title = markernames[i])

#########################################################
###########     Ploidy                   ################
#########################################################

#Look at mitotic stage by cluster
DimPlot(triflygerm,reduction = "umap",group.by= "Phase")

#in order to phase haploid cells we need to first get snps so that we can find cells that are 100% homozygous
cov <- readMM(file = "./SNP/scAlleleCount/covmat.mtx")
alt <- readMM(file = "./SNP/scAlleleCount/altmat.mtx")
ref <- readMM(file = "./SNP/scAlleleCount/refmat.mtx")
#So ref has the number of reads of a reference allele (column) per barcode cell (row), alt has the alternative allele and cov has the total number of reads for the site
#So if I want to know if it is homozygous for either than I can divide ref and alt by cov, if either is 1 then it is homo for either and if both are < 1 then it is hetero
#Lets first cut the cells to the ones that I am interested in
cells <- match(rownames(trigerm@meta.data),cellBarcodes)
covsub <- cov[cells,]
altsub <- alt[cells,]
refsub <- ref[cells,]
#To divide first get summary in triplet format
scov <- summary(covsub)
salt <- summary(altsub)
sref <- summary(refsub)
#Combine 2 triplets into a single table
saltcov <- merge(salt, scov, by=c("i", "j"))
altcov <- sparseMatrix(i=saltcov[,1], j=saltcov[,2], x=saltcov[,3]/saltcov[,4]) #Turn back into sparse matrix and divide matching values

srefcov <- merge(sref, scov, by=c("i", "j"))
refcov <- sparseMatrix(i=srefcov[,1], j=srefcov[,2], x=srefcov[,3]/srefcov[,4])

#So I want to get two things from this, finding out the % het of each cell to get which cells are diploid v haploid and to look at the x and y to find which of the haploid cells are y only.
sac <- summary(altcov)
src <- summary(refcov)

marc <- merge(sac, src, by=c("i", "j"))
het <- sparseMatrix(i=marc[,1], j=marc[,2], x=rep(1,nrow(marc)))
hetnum <- rowSums(het) #Number of heterozygous
covnum <- rowSums(sparseMatrix(i=scov[,1], j=scov[,2], x=rep(1,nrow(scov))))
hetperc <- hetnum/covnum
#Now I can make a vector to add to the metadata for trigermplus
dip <- rep("diploid",nrow(trigermplus@meta.data))
dip[hetperc < 0.025 ] <- "haploid"
#Now add it to the seurat object
triflygerm[['Ploidy']] <- dip
DimPlot(triflygerm,reduction = "umap",group.by= "Ploidy")

#########################################################
###########     Unpooling                ################
#########################################################

#We will cluster the allele table in order to determine which cells originated in the same cluster
#First use the data to make a single sparse matrix with 0,1,2 for homozygous reference, heterozygous, or homozygous alt
hetra <- data.frame(marc[,c(1,2)], x = rep(1,nrow(marc)))
homr <- src[src$x == 1,]
homr$x <- rep(0,nrow(homr))
homa <- sac[sac$x == 1,]
homa$x <- rep(2,nrow(homa))
geno <- rbind(hetra,homr,homa)
genotype <- sparseMatrix(i = geno[,1], j = geno[,2],x = geno[,3]) 
#Calculate a distance matrix over the genotype matrix
dist <- cosSparse(t(genotype))
clust <- hclust(as.dist(dist))
#There were 4 individuals so set k = 4 
mycl <- cutree(clust,k = 4)
#Add population back into the dataset
triflygerm[["population"]] <- as.character(mycl)
DimPlot(triflygerm,reduction = "umap",group.by = "population")
#Run statistics
datbl <- table(Idents(triflygerm),triflygerm@meta.data$population)
stats <- data.frame(Sum = apply(datbl,1,sum), Mean = apply(datbl,1,mean), SD = apply(datbl,1,sd))
poptbl <- data.frame(as.data.frame.matrix(as.matrix(datbl)),stats)
write.table(poptbl,file = "./poptable.tsv",quote = F,sep = "\t")
ggplot(triflygerm@meta.data,aes(group = "RNA_snn_res.0.5",x = "population")) +
  geom_bar(stat = "count")


#########################################################
###########       GSEA                   ################
#########################################################

#In order to find patterns of expressions among clusters we need to perform gene set enrichment
#We can read in biological process GO data for drosophila and then connect it back to beetle through annotation
gs <- readLines(con = "./Annotation/gene_set_Biological_Process_2018.txt")
gs <- strsplit(gs,"\t\t")
gsnames <- unlist(lapply(gs,function(x) x[[1]]))
gsval <- strsplit(unlist(lapply(gs,function(x) x[[2]])),"\t")
names(gsval) <- gsnames
gsval_tri <- lapply(gsval,function(x) rownames(subset(anno, D.Symbol %in% x)))
#Now we can perform a per cell level enrichment of the gene sets
ES.triflygerm <- enrichIt(obj = triflygerm,
                          gene.sets = gsval_tri,
                          groups = 1000, cores = 12,
                          min.size = 5)
write.table(ES.triflygerm,file = "./GSEA/Triflygerm_BiologicalProcess_fly_results.tsv",quote = F,sep = "\t")
triflygerm_gsea <- AddMetaData(triflygerm, ES.triflygerm)
colors <- colorRampPalette(c("#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF"))


#########################################################
###########       X Linked expression    ################
#########################################################

#We should look at how genes on the x axis are expressed in each cluster

exp <- AverageExpression(
  triflygerm,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)
exp <- exp$RNA
#Next we attach position information to each gene and transform it by cluster
exp_df <- cbind(anno[(rownames(exp)),c(2:4)],exp)
exp_df_long <- gather(exp_df,Cluster,Expression,-Chrom,-start,-ID,,factor_key = T)
exp_df_long <- exp_df_long[!is.na(exp_df_long$Chrom),]
exp_df_long$Chrom <- factor(exp_df_long$Chrom,levels = c("X","2","3","4","5","6","7","8","9","10"))
exp_df_long$Cluster <- factor(exp_df_long$Cluster, levels = cell_types[cols+1])
# rownames(exp_df_long) <- exp_df_long$ID

ggplot(exp_df_long,aes(x = start,y = Expression))+
  geom_line() +
  facet_grid(Cluster~Chrom,scales = "free_x")

ggplot(subset(exp_df_long, Chrom == "X"),aes(x = start,y = Expression))+
  geom_line() +
  facet_grid(Cluster~Chrom) +
  coord_cartesian(ylim=c(0,.1))+
  theme(strip.text.y.right = element_text(angle = 0))

xexp <- scale(colMeans(GetAssayData(triflygerm)[rownames(subset(exp_df,Chrom == "X")),]))
xexp[is.infinite(xexp)] <- 0
xexp <- pmin(xexp,1.5)
hist(xexp)
triflygerm[["logX_expression"]] <- xexp
FeaturePlot(triflygerm,features = "logX_expression",max.cutoff = 1,col = c("yellow","blue"))


ggplot(triflygerm@meta.data, aes(x = as.factor(RNA_snn_res.0.5), y = logX_expression)) +
  geom_violin() +
  facet_wrap(~Ploidy)
