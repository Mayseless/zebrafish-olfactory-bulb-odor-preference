ChooseClusterResolutionDownsample <- function(
    input.srobj, n.pcs, sample.name =  format(paste(Sys.Date(), "SilhouetteScore", sep="_")),
    res.low = .01, res.high=10, res.n = 40, bias = "over", subsample = FALSE,print_fig = TRUE, figdir=F) {
  
  ######## step 1: save the input seurat object as a new temporary object, 
  ########  we dont want to overwrite or change the original one with all of the parameter scans
  
  srobj.tmp = input.srobj 
  # in case there have been other things calculated in the metadata, just cut down to simplify/avoid errors
  srobj.tmp@meta.data = srobj.tmp@meta.data[,c(2:3)] # should just be the nUMI and nGene	
  
  
  ######## step 2: calculate the FindClusters over a large range of resolutions
  print("Performing parameter scan over multiple resolutions...")
  
  set.res = round(exp(seq(log(res.low), log(res.high), length.out=res.n)), digits=3)
  srobj.tmp = FindClusters(srobj.tmp, dims.use = n.pcs, k.param=20,
                           resolution=set.res[1], save.SNN=T, 
                           verbose=FALSE)
  
  for(i in 2:length(set.res)){
    srobj.tmp = FindClusters(
      srobj.tmp, resolution=set.res[i], verbose=FALSE)
    print(paste("          ", round(100*i/length(set.res)), "% done with parameter scan", sep=""))
  }
  
  
  ######## step 3: output plot of how the resolution changes the number of clusters you get
  n.clusters = vector(mode="numeric", length=length(set.res))
  names(n.clusters) = set.res
  for(i in 1:length(n.clusters)){
    n.clusters[i] = length(table(as.vector(srobj.tmp@meta.data[,paste("RNA_snn_res.", names(n.clusters)[i], sep="")])))
  }
  
  ######## step 4: calculate the silhouette width for each resolution
  print("Computing a silhouette width for each cell, for each resolution...")
  require(cluster)
  
  dist.temp = cor(t(srobj.tmp@reductions$pca@cell.embeddings[,1:n.pcs]), method="pearson")
  # Subsample based on the input parameter
  if (subsample) {
    random.cells.choose = sample(1:nrow(dist.temp), round(nrow(dist.temp) / 10, digits = 0))  # Subsample 10%
  } else {
    random.cells.choose = 1:nrow(dist.temp)  # Use all cells
  }
  
  dist.temp.downsample = dist.temp[random.cells.choose, random.cells.choose]
  sil.all.matrix = matrix(data=NA, nrow=nrow(dist.temp.downsample), ncol=0)
  
  for(i in 1:length(set.res)){
    clusters.temp = as.numeric(as.vector(
      srobj.tmp@meta.data[random.cells.choose,paste("RNA_snn_res.", set.res[i], sep="")]))
    if(length(table(clusters.temp))>1){
      sil.out = silhouette(clusters.temp, as.dist(1-as.matrix(dist.temp.downsample)))
      sil.all.matrix = cbind(sil.all.matrix, sil.out[,3])
    }
    if(length(table(clusters.temp))==1){
      sil.all.matrix = cbind(sil.all.matrix, rep(0, length(clusters.temp)))
    }
    print(paste("          ", round(100*i/length(set.res)), "% done with silhouette calculation", sep=""))
    
  }
  
  ######## step 5: calculate summary metric to compare the silhouette distributions,
  ########  average has worked well so far... could be different
  
  print("Identifying a best resolution to maximize silhouette width")
  sil.average = colMeans(sil.all.matrix)
  names(sil.average) = set.res
  
  
  ######## step 6: automate choosing resolution that maximizes the silhouette 
  hist.out = hist(sil.average, length(sil.average)/1.2,  plot=FALSE)
  
  #  take the ones that fall into the top bin, 
  #  and the max OR MIN of those  ******* can change this to under vs over cluster
  if(bias=="over"){
    resolution.choice = as.numeric(max(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  if(bias=="under"){
    resolution.choice = as.numeric(min(
      names(sil.average[which(sil.average>hist.out$breaks[length(hist.out$breaks)-1])])))
  }
  
  # get the silhouette of the best resolution: 
  silhouette.best = as.numeric(sil.average[paste(resolution.choice)])
  
  print(paste("Best Resolution Choice: ", resolution.choice, ", with average silhouette score of: ",
              round(silhouette.best, digits=3), ", giving ", as.numeric(n.clusters[paste(resolution.choice)]),
              " clusters", sep=""))
  
  
  ######### step 7: output plot and data 
  
  if (print_fig) {
    setwd(figdir)
    
    print(paste0("Outputting summary statistics and returning seurat object... ",
                 "This will create a pdf in your output directory,",
                 " and will return your input seurat object amended with the best choice",
                 " for clusters (found as Best.Clusters in the meta.data matrix, and set to your new ident)..."))
    
    pdf(paste(sample.name, ".pdf", sep=""),
        width=10, height=4, useDingbats=FALSE)
    par(mfrow=c(1,3))
    # Resolution vs # of Clusters
    plot(set.res, n.clusters, col="black", pch=19,
         type="p", xlab="Resolution", ylab="# Clusters",
         main="Resolution vs. # Clusters")
    # Resolution vs Average Silhouette
    plot(set.res, sil.average, col="black", pch=19,
         type="p", xlab="Resolution", ylab="Average Silhouette",
         main="Resolution vs. Average Silhouette")
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=resolution.choice, col="dodgerblue2", lty=2)
    
    # N Clusters vs Average Silhouette
    plot(n.clusters, sil.average, col="black", pch=19,
         type="p", xlab="# Clusters", ylab="Average Silhouette",
         main="# Clusters vs. Average Silhouette")
    abline(h=hist.out$breaks[length(hist.out$breaks)-1], col="firebrick3", lty=2)
    abline(v=as.numeric(n.clusters[paste(resolution.choice)]), col="dodgerblue2", lty=2)
    dev.off()
  }
  
  ######## step 8: return the original seurat object, with the metadata containing a 
  ########         concatenated vector with the clusters defined by the best choice here,
  ########         as well as the ident set to this new vector
  
  Best.Clusters = srobj.tmp@meta.data[,paste("RNA_snn_res.", resolution.choice, sep="")]
  
  input.srobj$Best.Clusters = Best.Clusters
  Idents(input.srobj) = input.srobj$Best.Clusters
  input.srobj@misc$resolution.choice <- resolution.choice
  return(input.srobj)
}


# dir = "Path_to_output_dir"
# 
# scData.combined <- ChooseClusterResolutionDownsample(PN_subset, 30,figdir = dir)
