ref_mapping<-function(ref_data,query_data,features_to_anchor,ref_dataset_title,query_dataset_title,cell_lbl ="seurat_clusters")
  {
#source('C:/Oded_data/single_cell_seq/dev/scicore_scripts/cluster_comparison_2.R')
lable_transfer_ref<-ref_data
lable_transfer_query<-query_data

if (is.null(features_to_anchor)) {

  anchors <- FindTransferAnchors(reference = lable_transfer_ref, query = lable_transfer_query,
                               dims = 1:30 ,reference.reduction = "pca",scale = FALSE)
} else {
  anchors <- FindTransferAnchors(reference = lable_transfer_ref, query = lable_transfer_query, features = features_to_anchor,
                                 dims = 1:30 ,reference.reduction = "pca",scale = FALSE)
}
predictions <- TransferData(anchorset = anchors, refdata = lable_transfer_ref$seurat_clusters,
                            dims = 1:30)
lable_transfer_query <- AddMetaData(lable_transfer_query, metadata = predictions)

lable_transfer_ref <- RunUMAP(lable_transfer_ref, dims = 1:30, reduction = "pca", return.model = TRUE)
lable_transfer_query <- MapQuery(anchorset = anchors, reference = lable_transfer_ref, query = lable_transfer_query,
                                 refdata = list(celltype = cell_lbl), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(lable_transfer_ref, reduction = "umap", label = TRUE, label.size = 3,
              repel = TRUE) + ggtitle(ref_dataset_title) + 
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank()) +
              NoLegend()

p2 <- DimPlot(lable_transfer_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + ggtitle(query_dataset_title)+ 
              theme_minimal() +
              theme(plot.title = element_text(hjust = 0.5),
                    axis.text = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks = element_blank()) +
              NoLegend()
p1 + p2

temp<-as.data.frame(lable_transfer_query$predicted.id)
temp2<-temp$`lable_transfer_query$predicted.id`

a<-cluster_comparison_3(lable_transfer_ref$seurat_clusters,lable_transfer_query$predicted.celltype,ref_dataset_title,query_dataset_title,do_bootstrap = TRUE, B_boot = 1000)

result <- list(
  plots = list(p1, p2, a$plot),
  mapped_query = lable_transfer_query,
  comparison = a   # <--- Save all cluster_comparison_2 output (which includes proportions)
)

return(result)
}
