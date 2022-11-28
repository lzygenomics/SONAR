#In this visualization section, we follow part of the CARD's visualization code(Ying Ma, Xiang Zhou. Nature Biotechnology)
#pie plot
#' SONAR.visualize.pie
#'
#' @param proportion deconvolution results matrix from SONAR(Details and formats are in vignettes)
#' @param spatial_location spatial coordinates matrix(Details and formats are in vignettes)
#' @param colors color vector that you want to use
#' @import ggplot2 gtools scatterpie
#' @return ggplot object
#' @export

SONAR.visualize.pie <- function(proportion,spatial_location,colors = colors){
  res_SONAR = as.data.frame(proportion)
  res_SONAR = res_SONAR[,mixedsort(colnames(res_SONAR))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_SONAR)==rownames(location))!= nrow(res_SONAR)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  data = cbind(res_SONAR,location)
  ct.select = colnames(res_SONAR)
  p = suppressMessages(ggplot() + geom_scatterpie(aes(x=x, y=y,r = 200),data=data,
                                                  cols=ct.select,color=NA) + coord_fixed(ratio = 1) +
                         scale_fill_manual(values =  colors)+
                         theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 15,face="bold"),
                               legend.text=element_text(size = 15),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'in'),
                               strip.text = element_text(size = 15,face="bold",vjust = 2),
                               legend.position="bottom"
                         )+
                         scale_y_reverse()+
                         guides(fill=guide_legend(title="Cell Type")))
  return(p)
}
#proportion plot
#' SONAR.prop
#'
#' @param scale if TRUE, the proportion will be scaled by Min-Max Scaling; else use the original predicted proportion
#' @param proportion deconvolution results matrix from SONAR(Details and formats are in vignettes)
#' @param spatial_location spatial coordinates matrix(Details and formats are in vignettes)
#' @param select_type the types you want to visualize
#' @param colors color vector that you want to use
#' @param NumCols decide the col-numbers you want to arrange
#' @import ggplot2 gtools reshape2
#' @return ggplot object
#' @export

SONAR.prop <- function(scale=F,proportion,spatial_location,select_type,colors = c("grey","red"),NumCols){
  res_SONAR = as.data.frame(proportion)
  res_SONAR = res_SONAR[,order(colnames(res_SONAR))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_SONAR)==rownames(location))!= nrow(res_SONAR)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  select_type = select_type
  res_SONAR = res_SONAR[,select_type]
  if(scale==F){
    res_SONAR_scale = as.data.frame(apply(res_SONAR,2,function(x){
      x}))
  }else{
    res_SONAR_scale = as.data.frame(apply(res_SONAR,2,function(x){
      (x - min(x)) / (max(x) - min(x))}))
  }
  res_SONAR_scale$x = as.numeric(location$x)
  res_SONAR_scale$y = as.numeric(location$y)
  mData = melt(res_SONAR_scale,id.vars = c("x","y"))
  colnames(mData)[3] <- "Cell_Type"
  b = c(0,1)
  p = suppressMessages(ggplot(mData, aes(x, y)) +
                         geom_point(aes(colour = value),size = 3.0) +
                         scale_color_gradientn(colours = colors) +
                         #scale_color_viridis_c(option = 2)+
                         #scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+
                         facet_wrap(~Cell_Type,ncol = NumCols)+
                         coord_fixed(ratio = 1)+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "in"),
                               #legend.position=c(0.14,0.76),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 15,face="bold",vjust = 1),
                               legend.text=element_text(size = 14),
                               strip.text = element_text(size = 13,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'in'))+
                         labs(color="proportion")+
                         scale_y_reverse())
  return(p)
}
