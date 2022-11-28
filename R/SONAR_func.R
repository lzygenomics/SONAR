#source("helper_preprocess.R")
#' SONAR.preprocess
#'
#' @param sc_count scRNA-seq expression matrix of reference(Details and formats are in vignettes)
#' @param sc_cell_type cell type annotation for reference(Details and formats are in vignettes)
#' @param sc_nUMI cell library size(Details and formats are in vignettes)
#' @param sp_coords spot expression coordinate(Details and formats are in vignettes)
#' @param sp_count spot expression matrix(Details and formats are in vignettes)
#' @param sp_nUMI spot library size(Details and formats are in vignettes)
#' @param cores the number of CPU threads that you want to use
#' @param type_min_cell Specifies a threshold to filter out cell types with a smaller number of cells
#' @param spot_min_UMI Specifies a threshold to filter out spots with fewer UMIs
#' @param preclus_resolution Specifies the resolution of the preclustering
#' @import methods Seurat utils Matrix
#' @return preprocessed data
#' @export

SONAR.preprocess<-function(sc_count,sc_cell_type,sc_nUMI,sp_coords,sp_count,sp_nUMI,cores=8,type_min_cell=0,spot_min_UMI=100,preclus_resolution=0.7){
  reference <- Reference(sc_count,sc_cell_type,sc_nUMI)
  puck <- SpatialRNA(sp_coords,sp_count,sp_nUMI)
  SONAR_obj <- create.SONAR(puck, reference, max_cores = cores, CELL_MIN_INSTANCE = type_min_cell,UMI_min = spot_min_UMI)
  SONAR_obj <-fitBulk(SONAR_obj)
  #u[j,t] is preprocessed type t,gene j 's expression
  u <- SONAR_obj@cell_type_info$renorm[[1]]
  #y[j,n]is spot n,gene j 's expression
  y <- as.matrix(SONAR_obj@spatialRNA@counts)
  #N[n] spot n library size
  N <- SONAR_obj@spatialRNA@nUMI
  coord<-SONAR_obj@spatialRNA@coords
  #preclustering
  print("####Part 4: ")
  print("Begin: Pre-clustering")
  yy<-CreateSeuratObject(y)
  yy<-SCTransform(yy)
  yy <- RunPCA(yy, assay = "SCT", verbose = FALSE)
  yy <- FindNeighbors(yy,reduction = "pca", dims=1:30)
  yy <- FindClusters(object = yy, verbose = T, resolution = preclus_resolution)
  yy <- RunUMAP(yy, reduction = "pca", dims = 1:10, label = T)
  #DimPlot(yy, reduction = "umap")
  precluster_label<-as.numeric(Idents(yy))
  print("End: Pre-clustering")
  p<-list(u=u,y=y,N=N,coord=coord,precluster_label=precluster_label)
  print("All preprocessing is complete!")
  return(p)
}
#' SONAR.deliver
#'
#' @param processed_data the data from preprocessed steps
#' @param path the path of delivering the preprocessed data to SONAR
#' @import utils
#' @return this part is to save the files, return the complete state
#' @export

SONAR.deliver<-function(processed_data,path){
  write.table(processed_data$u,file = paste(path,'u.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'u.txt',sep=""))) {
    stop("reference file error")
  }
  write.table(processed_data$y,file = paste(path,'y.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'y.txt',sep=""))) {
    stop("spatial file error")
  }
  write.table(processed_data$N,file = paste(path,'N.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'N.txt',sep=""))) {
    stop("spots library file error")
  }
  write.table(processed_data$coord,file=paste(path,'coord.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'coord.txt',sep=""))) {
    stop("coord file ready error")
  }
  write.table(processed_data$precluster_label,file=paste(path,'label.txt',sep=""),sep = ",")
  if (!file.exists(paste(path,'label.txt',sep=""))) {
    stop("precluster file error")
  }
  return(print("deliver complete!"))
}
#' SONAR.deconvolute
#'
#' @param fname the path to SONAR's MATLAB core codes
#' @param path the path to preprocessed data
#' @param h bandwidth
#' @param verbose show the command in MATLAB
#' @param desktop Determines whether the matlab run process is displayed
#' @param splash Determines whether the matlab run process is displayed
#' @param display Determines whether the matlab run process is displayed
#' @param wait wait the SONAR's deconvolution complete
#' @param single_thread Determine whether to use only single threads
#' @import matlabr
#' @return the state of deconvolution
#' @export

SONAR.deconvolute<-function (fname,path,h,verbose = TRUE, desktop = FALSE, splash = FALSE,
                             display = FALSE, wait = FALSE, single_thread = FALSE)
{
  stopifnot(file.exists(fname))
  matcmd = get_matlab(desktop = desktop, splash = splash,
                      display = display, wait = wait, single_thread = single_thread)
  cmd=paste0("userpath('",path,"');","thepath='",path,"';h=",h,";","SONAR_main;","exit")
  cmd = paste0(matcmd, '\"',cmd, '\"')
  if (verbose) {
    message("Command run is:")
    message(cmd)
  }
  x <- system(cmd, wait = wait)
  return(x)
}
