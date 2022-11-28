#In this preprocess part, we cite RCTD's preprocess code(Dylan M. Cable, et al. Nature Biotechnology.)
#Citation in our Methods parts.

setClass("Reference",
         slots = c(
           cell_types = "factor",
           counts = "dgCMatrix",
           nUMI = "numeric"
         ),
         prototype = list(
           cell_types = NULL,
           counts = NULL,
           nUMI = NA_integer_
         )
)
Reference <- function(counts, cell_types, nUMI = NULL, require_int = TRUE) {
  counts <- check_counts(counts, 'Reference', require_2d = T, require_int = require_int)
  if(is.null(nUMI)) {
    nUMI = Matrix::colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'Reference', require_2d = T, require_int = require_int)
  }
  check_cell_types(cell_types)
  barcodes <- intersect(intersect(names(nUMI), names(cell_types)), colnames(counts))
  if(length(barcodes) == 0)
    stop('Reference: cell_types, counts, and nUMI do not share any barcode names. Please ensure that names(cell_types)
         matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),length(cell_types),dim(counts)[2]))
    warning('Reference: some barcodes in nUMI, cell_types, or counts were not mutually shared. Such barcodes were removed.')
  if(sum(nUMI[barcodes] != Matrix::colSums(counts[,barcodes])) > 0)
    warning('Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
            is intended, there is no problem.')
  missing_cell_types <- names(which(table(cell_types[barcodes]) == 0))
  if(length(missing_cell_types) > 0)
    warning(paste('Reference: missing cell types with no occurences: ',paste(missing_cell_types,collapse=', ')))
  new("Reference", cell_types = cell_types[barcodes], counts = counts[,barcodes], nUMI = nUMI[barcodes])
}

check_cell_types <- function(cell_types) {
  if(class(cell_types) != 'factor')
    stop('Reference: cell_types is not a factor. Please format cell_types as a factor.')
  if(length(cell_types) < 2)
    stop('Reference: length(cell_types) < 2. cell_types needs to be a factor with length equal to the number of cells.')
  if(length(levels(cell_types)) < 2)
    stop('Reference: length(levels(cell_types)) < 2. cell_types needs to be a factor with multiple levels for each cell type.')
  if(is.null(names(cell_types)))
    stop('Reference: names(cell_types) is null. Please enter cell barcodes as names')
}

convert_old <- function(old_reference) {
  cell_types <- old_reference@meta.data$liger_ident_coarse
  nUMI <- old_reference@meta.data$nUMI
  names(cell_types) <- rownames(old_reference@meta.data);
  names(nUMI) <- rownames(old_reference@meta.data);
  Reference(old_reference@assays$RNA@counts, cell_types, nUMI)
}

restrict_reference <- function(reference, barcodes) {
  reference@counts <- reference@counts[,barcodes]
  reference@nUMI <- reference@nUMI[barcodes]
  reference@cell_types <- reference@cell_types[barcodes]
  return(reference)
}

restrict_reference_cell_types <- function(reference, cell_type_list) {
  new_ref <- (restrict_reference(reference, names(reference@cell_types)[reference@cell_types %in% cell_type_list]))
  new_ref@cell_types <- droplevels(new_ref@cell_types)
  return(new_ref)
}
check_counts <- function(counts, f_name, require_2d = F, require_int = T) {
  if(class(counts) != 'dgCMatrix') {
    if(class(counts) != 'matrix')
      tryCatch({
        counts <- as(counts,'matrix')
      }, error = function(e) {
        stop(paste0(f_name,': could not convert counts to matrix using as(counts,\'matrix\'). Please check that
             counts is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.'))
      })
    counts <- as(counts,"dgCMatrix")
  }
  if(dim(counts)[1] == 1) #check more than one gene
    stop(paste0(f_name,': the first dimension of counts is 1, indicating only one gene present. Please format counts so that
           the first dimension is greater than 1.'))
  if(dim(counts)[2] == 1)
    if(require_2d)
      stop(paste0(f_name,': the second dimension of counts is 1, indicating only one cell present. Please format counts so that
           the second dimension is greater than 1.'))
  else
    warning(paste0(f_name,': the second dimension of counts is 1, indicating only one cell/pixel present. If this is unintended,
        please format counts so that the second dimension is greater than 1.'))
  if(!is.numeric(counts[1,1]))
    stop(paste0(f_name,': elements of counts are not numeric'))
  if(require_int) {
    if(max(abs(counts %% 1)) > 1e-6)
      stop(paste0(f_name,': counts does not contain integers'))
  }
  if(is.null(rownames(counts)))
    stop(paste0(f_name,': rownames(counts) is null. Please enter gene names as rownames'))
  if(is.null(colnames(counts)))
    stop(paste0(f_name,': colnames(counts) is null. Please enter barcodes as colnames'))
  return(counts)
}
check_UMI <- function(nUMI, f_name, require_2d = F, require_int = T) {
  if(!is.atomic(nUMI))
    stop(paste0(f_name,': nUMI is not an atomic vector. Please format nUMI as an atomic vector.'))
  if(!is.numeric(nUMI))
    stop(paste0(f_name,': nUMI is not numeric'))
  if(require_int) {
    if(max(abs(nUMI %% 1)) > 1e-6)
      stop(paste0(f_name,': nUMI does not contain integers'))
  }
  if(is.null(names(nUMI)))
    stop(paste0(f_name,': names(nUMI) is null. Please enter barcodes as names'))
  if(length(nUMI) == 1)
    if(require_2d)
      stop(paste0(f_name,': the length of nUMI is 1, indicating only one cell present. Please format nUMI so that
           the length is greater than 1.'))
  else
    warning(paste0(f_name,': the length of nUMI is 1, indicating only one cell present. If this is unintended,
        please format nUMI so that the length is greater than 1.'))
}
setClass("SpatialRNA",
         slots = c(
           coords = "data.frame",
           counts = "dgCMatrix",
           nUMI = "numeric"
         ),
         prototype = list(
           coords = data.frame(NULL),
           counts = NULL,
           nUMI = NA_integer_
         )
)

fake_coords <- function(counts) {
  coords <- data.frame(Matrix(0,nrow=dim(counts)[2],ncol=2))
  colnames(coords) <- c('x','y')
  rownames(coords) <- colnames(counts)
  return(coords)
}

read.VisiumSpatialRNA <- function (datadir)
{
  coords <- readr::read_csv(file = paste(datadir, "spatial/tissue_positions_list.csv",
                                         sep = "/"),
                            col_names = c("barcodes", "in_tissue", "x", "y", "pxl_col_in_fullres", "pxl_row_in_fullres"))
  coords = tibble::column_to_rownames(coords, var = "barcodes")
  counts <- Seurat::Read10X_h5(paste0(datadir, "/filtered_feature_bc_matrix.h5"))
  puck = SpatialRNA(coords[,c('x','y')], counts)
  restrict_puck(puck, colnames(puck@counts))
}


SpatialRNA <- function(coords, counts, nUMI = NULL, use_fake_coords = FALSE, require_int = FALSE) {
  counts <- check_counts(counts, 'SpatialRNA', require_int = require_int)
  if(use_fake_coords)
    coords <- fake_coords(counts)
  else
    coords <- check_coords(coords)
  if(is.null(nUMI)) {
    nUMI = Matrix::colSums(counts)
    names(nUMI) <- colnames(counts)
  } else {
    check_UMI(nUMI, 'SpatialRNA', require_int = require_int)
  }
  barcodes <- intersect(intersect(names(nUMI), rownames(coords)), colnames(counts))
  if(length(barcodes) == 0)
    stop('SpatialRNA: coords, counts, and nUMI do not share any barcode names. Please ensure that rownames(coords)
         matches colnames(counts) and names(nUMI)')
  if(length(barcodes) < max(length(nUMI),dim(coords)[1],dim(counts)[2]))
    warning('SpatialRNA: some barcodes in nUMI, coords, or counts were not mutually shared. Such barcodes were removed.')
  if(sum(nUMI[barcodes] != Matrix::colSums(counts[,barcodes])) > 0)
    warning('SpatialRNA: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this
            is intended, there is no problem.')
  new("SpatialRNA", coords = coords[barcodes,], counts = counts[,barcodes], nUMI = nUMI[barcodes])
}

check_UMI <- function(nUMI, f_name, require_2d = F, require_int = T) {
  if(!is.atomic(nUMI))
    stop(paste0(f_name,': nUMI is not an atomic vector. Please format nUMI as an atomic vector.'))
  if(!is.numeric(nUMI))
    stop(paste0(f_name,': nUMI is not numeric'))
  if(require_int) {
    if(max(abs(nUMI %% 1)) > 1e-6)
      stop(paste0(f_name,': nUMI does not contain integers'))
  }
  if(is.null(names(nUMI)))
    stop(paste0(f_name,': names(nUMI) is null. Please enter barcodes as names'))
  if(length(nUMI) == 1)
    if(require_2d)
      stop(paste0(f_name,': the length of nUMI is 1, indicating only one cell present. Please format nUMI so that
           the length is greater than 1.'))
  else
    warning(paste0(f_name,': the length of nUMI is 1, indicating only one cell present. If this is unintended,
        please format nUMI so that the length is greater than 1.'))
}

check_counts <- function(counts, f_name, require_2d = F, require_int = T) {
  if(class(counts) != 'dgCMatrix') {
    if(class(counts) != 'matrix')
      tryCatch({
        counts <- as(counts,'matrix')
      }, error = function(e) {
        stop(paste0(f_name,': could not convert counts to matrix using as(counts,\'matrix\'). Please check that
             counts is coercible to matrix, such as a matrix, dgCmatrix, or data.frame.'))
      })
    counts <- as(counts,"dgCMatrix")
  }
  if(dim(counts)[1] == 1) #check more than one gene
    stop(paste0(f_name,': the first dimension of counts is 1, indicating only one gene present. Please format counts so that
           the first dimension is greater than 1.'))
  if(dim(counts)[2] == 1)
    if(require_2d)
      stop(paste0(f_name,': the second dimension of counts is 1, indicating only one cell present. Please format counts so that
           the second dimension is greater than 1.'))
  else
    warning(paste0(f_name,': the second dimension of counts is 1, indicating only one cell/pixel present. If this is unintended,
        please format counts so that the second dimension is greater than 1.'))
  if(!is.numeric(counts[1,1]))
    stop(paste0(f_name,': elements of counts are not numeric'))
  if(require_int) {
    if(max(abs(counts %% 1)) > 1e-6)
      stop(paste0(f_name,': counts does not contain integers'))
  }
  if(is.null(rownames(counts)))
    stop(paste0(f_name,': rownames(counts) is null. Please enter gene names as rownames'))
  if(is.null(colnames(counts)))
    stop(paste0(f_name,': colnames(counts) is null. Please enter barcodes as colnames'))
  return(counts)
}

check_coords <- function(coords) {
  if(class(coords) != 'data.frame') {
    tryCatch({
      coords <- as(coords,'data.frame')
    }, error = function(e) {
      stop('SpatialRNA: could not convert coords to data.frame using as(coords,\'data.frame\'). Please check that
           coords is coercible to data.frame, such as a matrix or data.frame .')
    })
  }
  if(dim(coords)[2] != 2) #check more than one gene
    stop('SpatialRNA: the second dimension of coords is not 2. Please enforce that dim(coords)[2] == 2 (x and y coordinates).')
  colnames(coords) <- c('x','y')
  if(!(is.numeric(coords$x) & is.numeric(coords$y)))
    stop('SpatialRNA: coords is not numeric')
  if(is.null(rownames(coords)))
    stop('SpatialRNA: rownames(coords) is null. Please enter barcodes as rownames')
  return(coords)
}

restrict_counts <- function(puck, gene_list, UMI_thresh = 1, UMI_max = 20000) {
  keep_loc = (puck@nUMI >= UMI_thresh) & (puck@nUMI <= UMI_max)
  puck@counts = puck@counts[gene_list,keep_loc]
  puck@nUMI = puck@nUMI[keep_loc]
  return(puck)
}

restrict_puck <- function(puck, barcodes) {
  barcodes = intersect(colnames(puck@counts), barcodes)
  puck@counts = puck@counts[,barcodes]
  puck@nUMI = puck@nUMI[barcodes]
  puck@coords = puck@coords[barcodes,]
  return(puck)
}


coerce_old <- function(puck) {
  new("SpatialRNA", coords = puck@coords, counts = puck@counts, nUMI = puck@nUMI)
}

create.SONAR <- function(spatialRNA, reference, max_cores = 8, test_mode = FALSE, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, UMI_min = 100, UMI_max = 20000000, UMI_min_sigma = 300,
                         class_df = NULL, CELL_MIN_INSTANCE = 5, cell_type_names = NULL) {

  config <- list(gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff, gene_cutoff_reg = gene_cutoff_reg, fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min, UMI_min_sigma = UMI_min_sigma, max_cores = max_cores,
                 N_epoch = 8, N_X = 50000, K_val = 100, N_fit = 1000, N_epoch_bulk = 30, MIN_CHANGE_BULK = 0.0001, MIN_CHANGE_REG = 0.001, UMI_max = UMI_max, MIN_OBS = 3)
  if(test_mode)
    config <- list(gene_cutoff = .00125, fc_cutoff = 0.5, gene_curoff_reg = 0.002, fc_cutoff_reg = 0.75, UMI_min = 1000,
                   N_epoch = 1, N_X = 50000, K_val = 100, N_fit = 50, N_epoch_bulk = 4, MIN_CHANGE_BULK = 1, MIN_CHANGE_REG = 0.001, UMI_max = 200000, MIN_OBS = 3, max_cores = 1, UMI_min_sigma = 300)
  if(is.null(cell_type_names))
    cell_type_names <- levels(reference@cell_types)
  cell_type_info <- list(info = process_cell_type_info(reference, cell_type_names = cell_type_names, CELL_MIN = CELL_MIN_INSTANCE), renorm = NULL)
  puck.original = restrict_counts(spatialRNA, rownames(spatialRNA@counts), UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  print("####Part 2: ")
  print('Find DE genes for regression: ')
  #puckMeans <- rowMeans(sweep(puck@counts, 2 , puck@nUMI, '/'))
  gene_list_reg = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff_reg, expr_thresh = config$gene_cutoff_reg, MIN_OBS = config$MIN_OBS)
  if(length(gene_list_reg) == 0)
    stop("Error: 0 regression differentially expressed genes found")
  print('Find DE genes for platform effect normalization: ')
  gene_list_bulk = get_de_genes(cell_type_info$info, puck.original, fc_thresh = config$fc_cutoff, expr_thresh = config$gene_cutoff, MIN_OBS = config$MIN_OBS)
  if(length(gene_list_bulk) == 0)
    stop("Error: 0 bulk differentially expressed genes found")
  puck = restrict_counts(puck.original, gene_list_bulk, UMI_thresh = config$UMI_min, UMI_max = config$UMI_max)
  puck = restrict_puck(puck, colnames(puck@counts))
  if(is.null(class_df))
    class_df <- data.frame(cell_type_info$info[[2]], row.names = cell_type_info$info[[2]]); colnames(class_df)[1] = "class"
  internal_vars <- list(gene_list_reg = gene_list_reg, gene_list_bulk = gene_list_bulk, proportions = NULL, class_df = class_df)
  new("SONAR", spatialRNA = puck, originalSpatialRNA = puck.original, reference = reference, config = config, cell_type_info = cell_type_info, internal_vars = internal_vars)
}
process_cell_type_info <- function(reference, cell_type_names, CELL_MIN = 5) {
  print("####Part 1: ")
  print("Start: Normalize reference matrix")
  print(paste("number of cells in reference:", dim(reference@counts)[2]))
  print(paste("number of genes in reference:", dim(reference@counts)[1]))
  cell_counts = table(reference@cell_types)
  print(cell_counts)

  if(min(cell_counts) < CELL_MIN)
    stop(paste0("Normalize reference matrix error: need a minimum of ",CELL_MIN, " cells for each cell type in the reference"))
  cell_type_info <- get_cell_type_info(reference@counts, reference@cell_types, reference@nUMI
                                       , cell_type_names = cell_type_names)
  print("End: Normalize reference matrix")
  return(cell_type_info)
}
get_cell_type_info <- function(raw.data, cell_types, nUMI, cell_type_names = NULL) {
  if(is.null(cell_type_names))
    cell_type_names = levels(cell_types)

  n_cell_types = length(cell_type_names)

  get_cell_mean <- function(cell_type) {
    cell_type_data = raw.data[,cell_types == cell_type]
    cell_type_umi = nUMI[cell_types == cell_type]
    normData = sweep(cell_type_data,2,cell_type_umi,`/`)
    return(rowSums(normData) / dim(normData)[2])
  }

  cell_type = cell_type_names[1]
  cell_type_means <- data.frame(get_cell_mean(cell_type))
  colnames(cell_type_means)[1] = cell_type
  for (cell_type in cell_type_names[2:length(cell_type_names)]) {
    cell_type_means[cell_type] = get_cell_mean(cell_type)
  }
  return(list(cell_type_means, cell_type_names, n_cell_types))
}
get_de_genes <- function(cell_type_info, puck, fc_thresh = 1.25, expr_thresh = .00015, MIN_OBS = 3) {
  total_gene_list = c()
  epsilon = 1e-9
  bulk_vec = rowSums(puck@counts)
  gene_list = rownames(cell_type_info[[1]])
  if(length(grep("mt-",gene_list)) > 0)
    gene_list = gene_list[-grep("mt-",gene_list)]
  gene_list = intersect(gene_list,names(bulk_vec))
  if(length(gene_list) == 0)
    stop("Find DE genes: Error: 0 common genes between SpatialRNA and Reference objects. Please check for gene list nonempty intersection.")
  gene_list = gene_list[bulk_vec[gene_list] >= MIN_OBS]
  for(cell_type in cell_type_info[[2]]) {
    if(cell_type_info[[3]] > 2)
      other_mean = rowMeans(cell_type_info[[1]][gene_list,cell_type_info[[2]] != cell_type])
    else {
      other_mean <- cell_type_info[[1]][gene_list,cell_type_info[[2]] != cell_type]
      names(other_mean) <- gene_list
    }
    logFC = log(cell_type_info[[1]][gene_list,cell_type] + epsilon) - log(other_mean + epsilon)
    type_gene_list = which((logFC > fc_thresh) & (cell_type_info[[1]][gene_list,cell_type] > expr_thresh)) #| puck_means[gene_list] > expr_thresh)
    print(paste0("Find DE genes: ", cell_type, " found ",length(type_gene_list)," DE genes"))
    total_gene_list = union(total_gene_list, type_gene_list)
  }
  total_gene_list = gene_list[total_gene_list]
  print(paste0("Find DE genes: total DE genes: ",length(total_gene_list)))
  return(total_gene_list)
}
setClass("SONAR",
         slots = c(
           spatialRNA = "SpatialRNA",
           originalSpatialRNA = "SpatialRNA",
           reference = "Reference",
           config = "list",
           cell_type_info = 'list',
           internal_vars = 'list',
           results = 'list',
           de_results = 'list',
           internal_vars_de = 'list'
         ),
         prototype = list(
           spatialRNA = NULL,
           originalSpatialRNA = NULL,
           reference = NULL,
           config = list(),
           cell_type_info = list(info = NULL, renorm = NULL),
           internal_vars = list(),
           results = list(),
           de_results = list(),
           internal_vars_de = list()
         )
)
#下面开始矫正平台效应
fitBulk <- function(SONAR) {
  #prepareBulkData这个函数在utils.R中，通过这一步，准备好Utj和所有空间转录组加起来的Ni
  bulkData <- prepareBulkData(SONAR@cell_type_info$info[[1]], SONAR@spatialRNA, SONAR@internal_vars$gene_list_bulk)
  print("####Part 3: ")
  print('Begin: platform effects normalization')
  #decompose_full在SONAR_helper.R中，注意以这里的传递的参数为准，具体的函数体用到了IRWLS那个R脚本里的函数
  #他现在应该是得到了空间转录组的pixel混合成的bulk的基因表达值；以及从单细胞里面得到的每一类type的gene的表达平均比例；我感觉他是要做回归
  #问题来了，他要用什么IRWSL的什么做回归呢？有没有指定类型
  #我估计他应该是普通的最小二乘回归，多元线性的，只不过用了迭代加权最小二乘回归更鲁棒了，可以对离群点有一些处理，constrain指的是权重和是否为一
  decompose_results <- decompose_full(bulkData$X, sum(SONAR@spatialRNA@nUMI),
                                      bulkData$b, verbose = F, constrain = F, MIN_CHANGE = SONAR@config$MIN_CHANGE_BULK,
                                      n.iter = 100, bulk_mode = T)
  SONAR@internal_vars$proportions <- decompose_results$weights
  SONAR@cell_type_info$renorm = SONAR@cell_type_info$info
  SONAR@cell_type_info$renorm[[1]] = get_norm_ref(SONAR@spatialRNA, SONAR@cell_type_info$info[[1]], SONAR@internal_vars$gene_list_bulk, decompose_results$weights)
  #######这边我自己也改造了一下，return之前的代码是我给加上的
  #SONAR@cell_type_info$lzynorm = SONAR@cell_type_info$info
  #SONAR@cell_type_info$lzynorm[[1]] = lzy_get_norm_ref(SONAR@spatialRNA, SONAR@cell_type_info$info[[1]], SONAR@internal_vars$gene_list_bulk, decompose_results$weights)
  print('End: platform effects normalization')
  return(SONAR)
}
prepareBulkData <- function(cell_type_means, puck, gene_list, MIN_OBS = 10) {
  bulk_vec = rowSums(puck@counts)
  gene_list <- intersect(names(which(bulk_vec >= MIN_OBS)),gene_list)
  nUMI = sum(puck@nUMI)
  X = as.matrix(cell_type_means[gene_list,] * nUMI)
  b = bulk_vec[gene_list]
  return(list(X=X, b=b))
}
decompose_full <- function(cell_type_profiles, nUMI, bead, constrain = TRUE, OLS = FALSE, verbose = F, n.iter = 50, MIN_CHANGE = 0.001, bulk_mode = F) {
  results = solveIRWLS.weights(cell_type_profiles,bead,nUMI,OLS = OLS, constrain = constrain,
                               verbose = verbose, n.iter = n.iter, MIN_CHANGE = MIN_CHANGE, bulk_mode = bulk_mode)
  return(results)
}
solveIRWLS.weights <-function(S,B,nUMI, OLS=FALSE, constrain = TRUE, verbose = FALSE,
                              n.iter = 50, MIN_CHANGE = .001, bulk_mode = F, solution = NULL){
  if(!bulk_mode)
    B[B > K_val] <- K_val
  if(OLS) {
    solution<-solveOLS(S,B, constrain = constrain) #first solve OLS, use this solution to find a starting point for the weights
    return(list(weights = solution, converged = T))
  }
  if(is.null(solution)) {
    solution <- numeric(dim(S)[2])
    solution[] <- 1/length(solution) #actually, just ignore the OLS solution
  }
  #solution <- runif(length(solution))*2 / length(solution) # random initialization
  names(solution) <- colnames(S)

  S_mat <<- matrix(0,nrow = dim(S)[1],ncol = dim(S)[2]*(dim(S)[2] + 1)/2)
  counter = 1
  for(i in 1:dim(S)[2])
    for(j in i:dim(S)[2]) {
      S_mat[,counter] <<- S[,i] * S[,j] # depends on n^2
      counter <- counter + 1
    }

  iterations<-0 #now use dampened WLS, iterate weights until convergence
  changes<-c()
  change<-1;
  while(change > MIN_CHANGE && iterations<n.iter){
    new_solution<-solveWLS(S,B,solution, nUMI,constrain=constrain, bulk_mode = bulk_mode)
    change<-norm(as.matrix(new_solution-solution))
    if(verbose) {
      print(paste("Change:",change))
      print(solution)
    }
    solution <- new_solution
    iterations<-iterations+1
  }
  return(list(weights = solution, converged = (change <= MIN_CHANGE)))
}
solveWLS<-function(S,B,initialSol, nUMI, bulk_mode = F, constrain = F){
  solution<-pmax(initialSol,0)
  prediction = abs(S%*%solution)
  threshold = max(1e-4, nUMI * 1e-7)
  prediction[prediction < threshold] <- threshold
  gene_list = rownames(S)
  derivatives <- get_der_fast(S, B, gene_list, prediction, bulk_mode = bulk_mode)
  d_vec <- -derivatives$grad
  D_mat <- psd(derivatives$hess)
  norm_factor <- norm(D_mat,"2")
  D_mat <- D_mat / norm_factor
  d_vec <- d_vec / norm_factor
  epsilon <- 1e-7; D_mat <- D_mat + epsilon * diag(length(d_vec))
  A<-cbind(diag(dim(S)[2]))
  bzero<- (-solution)
  alpha = 0.3
  if(constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1 - sum(solution),bzero)
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A_const,b_const,meq=1)$solution
  } else {
    solution <- solution + alpha*quadprog::solve.QP(D_mat,d_vec,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
}
solveOLS<-function(S,B, constrain = T){
  D<-t(S)%*%S
  d<-t(S)%*%B
  norm_factor <- norm(D,"2")
  D <- D / norm_factor
  d <- d / norm_factor
  epsilon <- 1e-7; D <- D + epsilon * diag(length(d))
  A<-cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  if(constrain) {
    A_const = t(rbind(1,A))
    b_const <-c(1 - sum(solution),bzero)
    solution <- quadprog::solve.QP(D,d,A_const,b_const,meq=1)$solution
  } else {
    solution <- quadprog::solve.QP(D,d,A,bzero,meq=0)$solution
  }
  names(solution)<-colnames(S)
  return(solution)
}
get_der_fast <- function(S, B, gene_list, prediction, bulk_mode = F) {
  if(bulk_mode) {
    #d1_vec <- -t(log(prediction) - log(B))
    #d2_vec <- -t(1/prediction)
    d1_vec <- -2*t((log(prediction) - log(B))/prediction)
    d2_vec <- -2*t((1 - log(prediction) + log(B))/prediction^2)
  } else {
    d1_d2 <- get_d1_d2(B, prediction)
    d1_vec <- d1_d2$d1_vec
    d2_vec <- d1_d2$d2_vec
  }
  grad = -d1_vec %*% S;
  hess_c <- -d2_vec %*% S_mat
  hess <- matrix(0,nrow = dim(S)[2], ncol = dim(S)[2])
  counter = 1
  for(i in 1:dim(S)[2]) {
    l <- dim(S)[2] - i
    hess[i,i:dim(S)[2]] <- hess_c[counter:(counter+l)]
    hess[i,i] <- hess[i,i] / 2
    counter <- counter + l + 1
  }
  hess <- hess + t(hess)
  return(list(grad=grad, hess=hess))
}
psd <- function(H) {
  eig <- eigen(H); epsilon = 1e-3
  if(length(H) == 1)
    P <- eig$vectors %*% pmax(eig$values,epsilon) %*% t(eig$vectors)
  else
    P <- eig$vectors %*% diag(pmax(eig$values,epsilon)) %*% t(eig$vectors)
  return(P)
}
get_norm_ref <- function(puck, cell_type_means, gene_list, proportions) {
  #bulk_vec这一步是将空间转录组的数据的各个点按照gene进行加和，把spot们变成一个bulk，是对列进行加和处理
  bulk_vec = rowSums(puck@counts)
  #weight_avg这一步：cell_type_means是R中create.SONAR之后得到的对象cell_type_info的第一列，也就是每一类的每个gene的平均表达比例
  ########然后cell_type_means数据结构为行名是gene，列名是类别名
  #还是weight_avg这一步：原本cell_type_means是各个类别每个基因的表达比例；现在先通过线性回归求出bulk中各细胞类型的比例，然后对这个基准的类别基因表达比例进行加权之后所有类别求和，得到的是这个bulk的对细胞类型加权后的，针对本bulk的每个基因的表达比例，表示的是通过ref表示出来的gene的比例情况
  #注意：2是对列的操作
  weight_avg = rowSums(sweep(cell_type_means[gene_list,],2,proportions / sum(proportions),'*'))
  #这个target_means是空间转录组数据简单相加成bulk后，gene_list中的基因的表达比例，表示的是空间转录组观测到的gene的比例情况
  target_means = bulk_vec[gene_list]/sum(puck@nUMI)
  #1是对行的操作
  #最后一行这个操作的意思是因为weight_avg表示的是bulk的通过根据ref的回归，所得到的由ref表示出来的bulk的基因们的表达水平，而target_means表示的是直接通过bulk简单的直观的看到的该基因的表达水平；因此这两个的比值其实就是ref和bulk的平台效应。举个例子，如果bulk简单直观的表达水平高，ref表示出来的gene表达水平低，也就是相除小于1，这时候让这个cell_type_means除以小于1的数，所以就变大了，所以相当于对ref的平均表达水平往高处调了，这就是对平台效应的矫正
  cell_type_means_renorm = sweep(cell_type_means[gene_list,],1,weight_avg / target_means,'/')
}



