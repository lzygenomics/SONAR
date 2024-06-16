# SONAR
SONAR is the algorithm of cell-type deconvolution for spatial transcriptomics

## Installation
#library(devtools)

devtools:: install_github("lzygenomics/SONAR")

## Dependence

- R version >= 4.0.5.
- R packages: 
  - this.path>=0.5.1; Matrix>=1.3.4; data.table>=1.14.0; Seurat>=4.0.3;  matlabr=1.5.2; R.matlab=3.7.0
- MATLAB version >= R2019a 
  - Just install the MATLAB, subsequent operations do not require you to use the MATLAB language

## Run SONAR

1. Install the dependence

2. Install SONAR

3. Download the SONAR files: **Follow the Example/SONAR-entrance.Rmd** IS ALL YOU NEED

## Files Annotation

- **/inst/extdata/**	

  *These files store the required input information, See SONAR-entrance.Rmd for specific format requirements.*


- **/result/**

  *These files store the output results*.

  Proportions of all cell types in each spot (**SONAR.results.txt**);

  Spatial pie plots (**pie.pdf**);

  Spatial distribution of specific cell types proportion (**abs_prop.pdf** / **scaled_prop.pdf**);

  Colocalization(correlation along the spatial) for pairs of cell types (**colocalization.pdf**).


- **/core-code/**	

  *no need to operate.*

  *These files store the core code, and store the preprocessed data that delivered to SONAR*.
  
## A brif Example

Please Follow the Example/SONAR-entrance.Rmd

## Publication Link

  https://www.nature.com/articles/s41467-023-40458-9

## How to cite SONAR

Liu, Z., Wu, D., Zhai, W. et al. SONAR enables cell type deconvolution with spatially weighted Poisson-Gamma model for spatial transcriptomics. Nat Commun 14, 4727 (2023). https://doi.org/10.1038/s41467-023-40458-9

## Thank you, and happy researching!
