#### Thank you for using our code, wish you all the best! ####
#### Key parameters to generate different scenario of Homo-Area can be found in line 36——42 ####  

#### Clear the work space
rm(list=ls())
library(Matrix)
library(CMP)
options(warn =-1)
set.seed(110)

#### Prepare the material for Simulation
# May need to modify the path of se_quartz.RDS
se_quartz <- readRDS("./se_quartz.RDS")
counts <- as(se_quartz@assays$RNA@counts,"dgCMatrix")
cell_types <- se_quartz@meta.data$nnet2
names(cell_types) <- colnames(counts)# create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type


#### Setting for Homo-Area simulation
# Cell number for dominant types, we default to two dominant types per spot
D1 <- 6
D2 <- 2
# Maximum type number of spots in this area 
T <- 4
# Specify the cell type names which have been selected
types<-sample(1:length(unique(cell_types)),T,replace = FALSE)
# Show the cell type information
cat("Dominant type has",N_s,"cells,","including type index",types,",the first denotes the dominant type")

#### Begin sampling, here we default area with size of 20*20 spots
type_name<-c("B cells","CD14+ Monocytes","CD4 T cells","CD8 T cells","Dendritic cells","FCGR3A+ Monocytes","HEK cells","NK cells")
# Keep true cell type proportion and expression matrix
true_prob<-matrix(0,nrow = 400,ncol = 8)
spot_counts<-matrix(0,nrow(counts),1)
for (m in 1:400) {
  # Sample the specific cell numbers for each cell types in each spot
  # rCMP can generate number follow general Poisson distribution, which can control the dispersion
  # mu define the expectation,nu define the dispersion,bigger nu will make less dispersion,nu = 1 is traditional Poisson distribution
  # x_max define the max number you want to generate
  # If you want to make the area more continuous or homogeneous, please set a big nu
  # If you want to make the area more random or heterogeneous, you can set a small nu
  # You can adjust different parameters for 4 types to get Homo-Areas with different property
  Dominant1<-rCMP(n=1,mu = D1,nu=20,x_max =7)
  Dominant2<-rCMP(n=1,mu = D2,nu=20,x_max =3)
  sub1<-rCMP(n=1,mu = 0.1,nu=2,x_max =1)
  sub2<-rCMP(n=1,mu = 0.1,nu=2,x_max =1)
  cell_pertype<-c(Dominant1,Dominant2,sub1,sub2)
  spot_num<-sum(cell_pertype)
  spot_types<-types
  spot_cell_pertype<-cell_pertype
  spot_counts_i <-matrix(0,nrow(counts),1)
  for (a in 1:length(spot_types)) {
    # prepare the submatrix of selected cell types
    type1_counts <- counts[,which(cell_types == type_name[spot_types[a]])]
    type1_counts <- as.matrix(type1_counts)
    spot_type1 <- type1_counts[,sample(ncol(type1_counts),spot_cell_pertype[a],replace = FALSE)]
    if (spot_cell_pertype[a]==1) {
      spot_counts_i <- cbind(spot_counts_i,spot_type1)
    }else{
      spot_counts_i <- cbind(spot_counts_i,rowSums(spot_type1))
    }
  }
  spot_counts <- cbind(spot_counts,rowSums(spot_counts_i))
  # save the true cell type proportions
  for (r in 1:length(spot_types)) {
    true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
  }
}
spot_counts <- spot_counts[,-1]
colnames(spot_counts) <- c(1:ncol(spot_counts))

#### define the coordinate
coords<-matrix(0,nrow=1,ncol=2)
for(cd in 1:20){
  xcoord <- c(1:20)
  xcoord <- 10*xcoord
  ycoord <- rep(10*cd,times=20)
  temp_ords <- cbind(xcoord,ycoord)
  coords <- rbind(coords,temp_ords)
}
coords <-coords[-1,]
coords <- as.data.frame(coords)
rownames(coords)<-c(1:nrow(coords))

#### Save the data
#save(spot_counts,coords,true_prob,file="./results.RData")