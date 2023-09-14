#### Thank you for using our code, wish you all the best! ####
# Note: In this script, our code comments are not very detailed,
# please see the simulation code of Homo-Area for detailed comments. 
# The main difference between Compo-Area code and Homo-Area code is to assign-
# different properties to each sub-region, which is reflected in the loop of different regions to generate spot
# To make each sub-region more homogeneous, in this code, we set a big nu parameters for dominant cell type in rCMP to control the dispersion
# You can adjust the nu parameter according to your requirement

#### Clear the work space
rm(list=ls())
library(Matrix)
library(CMP)
options(warn =-1)

#### Prepare the material for Simulation
# May need to modify the path of se_quartz.RDS
se_quartz <- readRDS("./se_quartz.RDS")
counts <- as(se_quartz@assays$RNA@counts,"dgCMatrix")
cell_types <- se_quartz@meta.data$nnet2
names(cell_types) <- colnames(counts)# create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type


#### Setting for Background of Compo-Area simulation
## Sub-region below the bottom small window
set.seed(43)
# Cell number for dominant types, we default to 1 dominant types per spot
N_s<-6
# Maximum type number of spots in this area 
T<-4
# Specify the cell type names which have been selected
types<-sample(1:length(unique(cell_types)),T,replace = FALSE)
# Show the cell type information
cat("Dominant type has",N_s,"cells,","including type index",types,",the first denotes the dominant type")

type_name<-c("B cells","CD14+ Monocytes","CD4 T cells","CD8 T cells","Dendritic cells","FCGR3A+ Monocytes","HEK cells","NK cells")

true_prob<-matrix(0,nrow = 800,ncol = 8)
spot_counts<-matrix(0,nrow(counts),1)
for (m in 1:140) {
  Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
  sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  cell_pertype<-c(Dominant,sub1,sub2,sub3)
  spot_num<-sum(cell_pertype)
  spot_types<-types
  spot_cell_pertype<-cell_pertype
  spot_counts_i <-matrix(0,nrow(counts),1)

  for (a in 1:length(spot_types)) {
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
  for (r in 1:length(spot_types)) {
    true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
  }
}

## Sub-region near the bottom small window 

for(kk in 0:5)
{
  # left part of small window
  for (m in (141+20*kk):(147+20*kk)) {
    Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
    sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    cell_pertype<-c(Dominant,sub1,sub2,sub3)
    spot_num<-sum(cell_pertype)
    spot_types<-types
    spot_cell_pertype<-cell_pertype
    spot_counts_i <-matrix(0,nrow(counts),1)
    
    for (a in 1:length(spot_types)) {
      
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
    
    for (r in 1:length(spot_types)) {
      true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
    }
  }
  # small window part
  set.seed(20)
  N_s<-6
  T<-4
  types<-sample(1:length(unique(cell_types)),T,replace = FALSE)
  cat("Dominant type has",N_s,"cells,","specifically type",types,",first is dominant")
  for (m in (148+20*kk):(153+20*kk)) {
    Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
    sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    cell_pertype<-c(Dominant,sub1,sub2,sub3)

    spot_num<-sum(cell_pertype)
    spot_types<-types
    spot_cell_pertype<-cell_pertype

    spot_counts_i <-matrix(0,nrow(counts),1)

    for (a in 1:length(spot_types)) {

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

    for (r in 1:length(spot_types)) {
      true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
    }
  }
  # right part of small window
  set.seed(43)
  N_s<-6
  T<-4
  types<-sample(1:length(unique(cell_types)),T,replace = FALSE)
  cat("Dominant type has",N_s,"cells,","specifically type",types,",first is dominant")
  for (m in (154+20*kk):(160+20*kk)) {
    Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
    sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    cell_pertype<-c(Dominant,sub1,sub2,sub3)
    spot_num<-sum(cell_pertype)
    spot_types<-types
    spot_cell_pertype<-cell_pertype
    spot_counts_i <-matrix(0,nrow(counts),1)
    for (a in 1:length(spot_types)) {
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
    for (r in 1:length(spot_types)) {
      true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
    }
  }
}
# small window to median stripe
for (m in 261:380) {

  Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
  sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  cell_pertype<-c(Dominant,sub1,sub2,sub3)

  spot_num<-sum(cell_pertype)
  spot_types<-types
  spot_cell_pertype<-cell_pertype

  spot_counts_i <-matrix(0,nrow(counts),1)

  for (a in 1:length(spot_types)) {

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

  for (r in 1:length(spot_types)) {
    true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
  }
}
## median stripe
set.seed(20)

N_s<-6

T<-4

types<-sample(1:length(unique(cell_types)),T,replace = FALSE)

cat("Dominant type has",N_s,"cells,","specifically type",types,",first is dominant")

for (m in 381:420) {
  
  Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
  sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  cell_pertype<-c(Dominant,sub1,sub2,sub3)
  spot_num<-sum(cell_pertype)
  spot_types<-types
  spot_cell_pertype<-cell_pertype
  spot_counts_i <-matrix(0,nrow(counts),1)
  for (a in 1:length(spot_types)) {
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
  for (r in 1:length(spot_types)) {
    true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
  }
}
## region between the median stripe and the next(above) small window  
set.seed(43)

N_s<-6

T<-4

types<-sample(1:length(unique(cell_types)),T,replace = FALSE)

cat("Dominant type has",N_s,"cells,","specifically type",types,",first is dominant")

for (m in 421:540) {
 
  Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
  sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  cell_pertype<-c(Dominant,sub1,sub2,sub3)
 
  spot_num<-sum(cell_pertype)
  spot_types<-types
  spot_cell_pertype<-cell_pertype

  spot_counts_i <-matrix(0,nrow(counts),1)

  for (a in 1:length(spot_types)) {
    
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
 
  for (r in 1:length(spot_types)) {
    true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
  }
}

## Sub-region near the top small window 
for(kk in 0:5)
{
  # left part of window
  for (m in (541+20*kk):(547+20*kk)) {
    
    Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
    sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    cell_pertype<-c(Dominant,sub1,sub2,sub3)
    
    spot_num<-sum(cell_pertype)
    spot_types<-types
    spot_cell_pertype<-cell_pertype
    
    spot_counts_i <-matrix(0,nrow(counts),1)
    
    for (a in 1:length(spot_types)) {
      
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
    
    for (r in 1:length(spot_types)) {
      true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
    }
  }
  # window
  set.seed(20)
  N_s<-6
  T<-4
  types<-sample(1:length(unique(cell_types)),T,replace = FALSE)
  cat("Dominant type has",N_s,"cells,","specifically type",types,",first is dominant")
  for (m in (548+20*kk):(553+20*kk)) {
    Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
    sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    cell_pertype<-c(Dominant,sub1,sub2,sub3)
    spot_num<-sum(cell_pertype)
    spot_types<-types
    spot_cell_pertype<-cell_pertype
    spot_counts_i <-matrix(0,nrow(counts),1)
    for (a in 1:length(spot_types)) {
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
    for (r in 1:length(spot_types)) {
      true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
    }
  }
  # right part
  set.seed(43)
  N_s<-6
  T<-4
  types<-sample(1:length(unique(cell_types)),T,replace = FALSE)
  cat("Dominant type has",N_s,"cells,","specifically type",types,",first is dominant")
  for (m in (554+20*kk):(560+20*kk)) {
    Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
    sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
    cell_pertype<-c(Dominant,sub1,sub2,sub3)
    spot_num<-sum(cell_pertype)
    spot_types<-types
    spot_cell_pertype<-cell_pertype
    spot_counts_i <-matrix(0,nrow(counts),1)
    for (a in 1:length(spot_types)) {
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
    for (r in 1:length(spot_types)) {
      true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
    }
  }
}
# rest part above the top window
for (m in 661:800) {
  
  Dominant<-rCMP(n=1,mu = N_s,nu=20,x_max =8)
  sub1<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub2<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  sub3<-rCMP(n=1,mu = 0.10,nu=2.3,x_max =2)
  cell_pertype<-c(Dominant,sub1,sub2,sub3)
  spot_num<-sum(cell_pertype)
  spot_types<-types
  spot_cell_pertype<-cell_pertype
  spot_counts_i <-matrix(0,nrow(counts),1)
  for (a in 1:length(spot_types)) {
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
  for (r in 1:length(spot_types)) {
    true_prob[m,spot_types[r]]<-spot_cell_pertype[r]/spot_num
  }
}
spot_counts <- spot_counts[,-1]
colnames(spot_counts) <- c(1:ncol(spot_counts))

#### define the coordinate
coords<-matrix(0,nrow=1,ncol=2)
for(cd in 1:40){
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