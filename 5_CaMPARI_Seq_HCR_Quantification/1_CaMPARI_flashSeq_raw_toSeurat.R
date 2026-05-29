library(dplyr)
library(Seurat)
library(patchwork)
library (tidyverse)


OM_20240412_CaMPARI_FS.data <- read.table("U:/Scientific Data/RG-AS04-Data01/Oded_Mayseless/sequencing_data/20240412_CampariFS_raw/processed_data/featureCounts/count_gene.txt", sep = "", header = T, row.names = 1)
OM_20240412_CaMPARI_FS_metadata <- OM_20240412_CaMPARI_FS.data[,1:5]
OM_20240412_CaMPARI_FS.data = OM_20240412_CaMPARI_FS.data[,-c(1:5)]
cellnames=colnames(OM_20240412_CaMPARI_FS.data)
cellnames=unlist(lapply(cellnames, function(x) strsplit(x,"STAR.")[[1]][2]))
cellnames=unlist(lapply(cellnames, function(x) strsplit(x,"_Aligned")[[1]][1]))
cellnames=sub("\\.", "_", cellnames)
colnames(OM_20240412_CaMPARI_FS.data)=cellnames
OM_20240412_CaMPARI_FS.data$ensembl_gene_id<-rownames(OM_20240412_CaMPARI_FS.data)
OM_20240412_CaMPARI_FS.data<-OM_20240412_CaMPARI_FS.data[,c(3073,1:3072)]  

#Extracting gene_names from ENSEMBL IDs -version 103- and adding them as another column
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl <- useEnsembl(biomart = 'genes', dataset = 'drerio_gene_ensembl', version = 103)

mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('drerio_gene_ensembl', mart)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'external_gene_name',
    'ensembl_gene_id',
    'gene_biotype'),
  uniqueRows = TRUE)

head(annotLookup)



#number of rows in each dataset
nrow(annotLookup)
nrow(OM_20240412_CaMPARI_FS.data)

##### change gene name if it's corresponding several ensembl ids #####
## to rename dup gene names
##1. list genes that are duplicated
##2. calculate #total reads for each ENS
##3. Rank ENS for each gene base on #total counts
##4. Rename only the ones with lower counts (keep the name for the one with highest count)
##5. Rename in sequence of counts, with special characters not overlapping with normal gene names (gene_name;2 gene_name#2 /)

annotLookup<-annotLookup[,c(2,1,3)]

colnames(annotLookup)<-c("ensembl_gene_id","external_gene_name","gene.type")
annotLookup<-subset(annotLookup,annotLookup$external_gene_name!="")

## only use the genes from count matrix
annotLookup<-annotLookup[annotLookup$ensembl_gene_id%in%OM_20240412_CaMPARI_FS.data$ensembl_gene_id,]

## check how many genes with one gene name to many ensembl id
annotLookup<-data.frame(annotLookup)
#gene.freq<-count(annotLookup, external_gene_name)
library(dplyr)

gene.freq <- annotLookup %>%
  filter(!is.na(external_gene_name), external_gene_name != "") %>%
  dplyr::count(external_gene_name, name = "n")

colnames(gene.freq)[2]<-"freq"
gene.to.many.ens<-subset(gene.freq,gene.freq$freq>1)
gene.to.many.ens<-annotLookup[annotLookup$external_gene_name%in%gene.to.many.ens$external_gene_name,]
gene.to.one.ens<-annotLookup[!annotLookup$external_gene_name%in%gene.to.many.ens$external_gene_name,]

## for each ensembl gene id, calculate the mean expression 
dup.gene.name<-unique(gene.to.many.ens$external_gene_name)

modified.gene.name<-data.frame(ensembl_gene_id=NA,mean.exp=NA,external_gene_name=NA)
for (i in c(1:length(dup.gene.name))) {
  temp1<-gene.to.many.ens[gene.to.many.ens$external_gene_name%in%dup.gene.name[i],]
  temp2<-OM_20240412_CaMPARI_FS.data[OM_20240412_CaMPARI_FS.data$ensembl_gene_id%in%temp1$ensembl_gene_id,]
  temp2<-data.frame(temp2)
  temp2$mean.exp<-apply(temp2[,-1],1, mean)
  temp2<-temp2[,c(1,ncol(temp2))]
  temp2$external_gene_name<-rep(dup.gene.name[i],nrow(temp2))
  # rank by expression
  temp2<-temp2[order(temp2$mean.exp,decreasing=T),]
  for (j in c(2:nrow(temp2))) {
    k<-j-1
    temp2$external_gene_name[j]<-paste0(temp2$external_gene_name[j],"_Alt_",k)
  }
  modified.gene.name<-rbind(modified.gene.name,temp2)
}
## merge gene.to.one.ens and modified.gene.name
gene.to.many.ens<-merge(gene.to.many.ens,modified.gene.name,by="ensembl_gene_id")
gene.to.many.ens<-gene.to.many.ens[,c(1,5,3)]
colnames(gene.to.many.ens)[2]<-"external_gene_name"
gene.name.id.new<-rbind(gene.to.one.ens,gene.to.many.ens)
## merge with expression count
OM_20240412_CaMPARI_FS.data<-merge(OM_20240412_CaMPARI_FS.data,gene.name.id.new,by="ensembl_gene_id")
rownames(OM_20240412_CaMPARI_FS.data)<-OM_20240412_CaMPARI_FS.data$external_gene_name

#debugging
rownames(OM_20240412_CaMPARI_FS.data) <- OM_20240412_CaMPARI_FS.data$external_gene_name
OM_20240412_CaMPARI_FS.data[] <- lapply(OM_20240412_CaMPARI_FS.data, as.numeric)
library(Matrix)
sparse_data <- Matrix(as.matrix(OM_20240412_CaMPARI_FS.data), sparse = TRUE)
OM_20240412_CaMPARI_FS <- CreateSeuratObject(counts = sparse_data, project = "OM_20240412_CaMPARI_FS", min.cells = 3, min.features = 200)
saveRDS(OM_20240412_CaMPARI_FS, file = "C:/Oded_data/single_cell_seq/data/20241104_OM_20240412_CaMPARI_FS.rds")
#OM_20240412_CaMPARI_FS <- readRDS(file = "C:/Oded_data/single_cell_seq/data/20241104_OM_20240412_CaMPARI_FS.rds")

### start seurat pipeline

OM_20240412_CaMPARI_FS <- CreateSeuratObject(counts = OM_20240412_CaMPARI_FS.data, project = "OM_20240412_CaMPARI_FS", min.cells = 3, min.features = 200, row.names = rownames(OM_20240412_CaMPARI_FS.data) )
OM_20240412_CaMPARI_FS

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
OM_20240412_CaMPARI_FS[["percent.mt"]] <- PercentageFeatureSet(OM_20240412_CaMPARI_FS, pattern = "^mt-")

# Visualize QC metrics as a violin plot
Campari_QC_violinPlots<-VlnPlot(OM_20240412_CaMPARI_FS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

save_dir<-"U:/Scientific Data/RG-AS04-Data01/Oded_Mayseless/innate_behavior/figures/20241115_single_cell/20250509"
# Visualize QC metrics as a violin plot
ggsave(filename = file.path(save_dir,"Campari_QC_violinPlots.svg"),
       plot = Campari_QC_violinPlots, device = "svg",
       width = 8,
       height = 8)
print(paste0("plot saved to ", file.path(save_dir, "Campari_QC_violinPlots.svg")))
# Extracting the metrics from metadata
cell_metrics <- OM_20240412_CaMPARI_FS@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]
# Calculating summary statistics for each metric
summary_stats <- summary(cell_metrics)
print(summary_stats)


#Subsetting after QC with >200 - < 12000 Feature counts, between 1750000 and 200k Count_RNA and <4% mtRNA (basically all cells)
OM_20240412_CaMPARI_FS <- subset(OM_20240412_CaMPARI_FS, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 10)
# Add a new column to the metadata that identifies the sample based on the cell names
OM_20240412_CaMPARI_FS$sample <- ifelse(grepl("E3", rownames(OM_20240412_CaMPARI_FS@meta.data)), "E3",
                                  ifelse(grepl("TCA", rownames(OM_20240412_CaMPARI_FS@meta.data)), "TCA",
                                         ifelse(grepl("ATP", rownames(OM_20240412_CaMPARI_FS@meta.data)), "ATP",
                                                ifelse(grepl("Cytidine", rownames(OM_20240412_CaMPARI_FS@meta.data)), "Cytidine", "Other"))))

# Extract the replicate number based on the pattern after the sample name
extract_replicate <- function(rownames) {
  replicate <- sub(".*_(E3|TCA|ATP|Cytidine|Cad)\\.(\\d).*", "\\2", rownames)
  return(replicate)
}

OM_20240412_CaMPARI_FS$replicate <- extract_replicate(rownames(OM_20240412_CaMPARI_FS@meta.data))

# Add combined sample and replicate information
OM_20240412_CaMPARI_FS$sample_replicate <- paste(OM_20240412_CaMPARI_FS$sample, OM_20240412_CaMPARI_FS$replicate, sep = "_")

# Create a data frame to count the number of each replicate
replicate_counts <- OM_20240412_CaMPARI_FS@meta.data %>%
  group_by(sample_replicate) %>%
  summarize(count = n()) %>%
  arrange(desc(count))

# Print the replicate counts
print(replicate_counts)


OM_20240412_CaMPARI_FS
OM_20240412_CaMPARI_FS <- NormalizeData(OM_20240412_CaMPARI_FS)
OM_20240412_CaMPARI_FS <-  FindVariableFeatures(OM_20240412_CaMPARI_FS,selection.method = 'vst')
OM_20240412_CaMPARI_FS <-  ScaleData(OM_20240412_CaMPARI_FS,features = rownames(OM_20240412_CaMPARI_FS))
OM_20240412_CaMPARI_FS <-  RunPCA(OM_20240412_CaMPARI_FS)
OM_20240412_CaMPARI_FS <-  FindNeighbors(OM_20240412_CaMPARI_FS,dims = 1:30,verbose = FALSE)


