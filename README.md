# SONAR: A Probabilistic Framework for Cell-Type Deconvolution in Spatial Transcriptomics

> **Important update:** This repository is kept as the original MATLAB-backed SONAR archive. A MATLAB-free R/Rcpp version is available at **[SONAR_R](https://github.com/lzygenomics/SONAR_R)**.

SONAR is an algorithm developed for cell-type deconvolution in spatial transcriptomics. It integrates spatial information in a balanced way to enhance performance and robustness.

## Key Features

- **Enhanced Signal-to-Noise Ratio**  
  Utilizes the similarity of spatial locations to boost the signal-to-noise ratio, helping to extract meaningful patterns from the data.

- **Local Spatial Heterogeneity Consideration**  
  Incorporates local spatial heterogeneity to avoid over-reliance on spatial information, ensuring that the inherent biological diversity is accurately captured.

- **Robust Probabilistic Framework**  
  Built on a probabilistic framework, SONAR delivers more robust and reliable results compared to deep learning–based approaches.
  
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

1. Install the dependence (pay attention that you need to install the **MATLAB**)

2. Install SONAR

3. Download the SONAR files (This **files structure** will help you run on the custom data)

4. Open **Example/SONAR-entrance.Rmd** , you could run and get the example results.

For running the custom dataset, you could modify the data preparation stage in **Example/SONAR-entrance.Rmd** with the same format, 
and substitude the single cell data and spatial data in this files structure

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
