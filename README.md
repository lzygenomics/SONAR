# SONAR
SONAR is the algorithm of cell-type deconvolution for spatial transcriptomics

## Dependence

- R version >= 4.0.5.
- R packages: 
  - this.path>=0.5.1; Matrix>=1.3.4; data.table>=1.14.0; Seurat>=4.0.3;  matlabr=1.5.2; R.matlab=3.7.0
- MATLAB version >= R2019a 
  - Just install the MATLAB, subsequent operations do not require you to use the MATLAB language

## Run SONAR

1. Install the dependence
2. Download the SONAR-code files
3. **Follow the  ./SONAR-code/core-code/SONAR.R** IS ALL YOU NEED

## Files Annotation

- **./SONAR-code/core-code/**	

  *These files store the core code*.
  
  Follow the ./SONAR-code/core-code/SONAR.R. 
  
  Other files are called automatically and generally require no additional action.

- **./SONAR-code/input/**	

  *These files store the required input information, See SONAR.R for specific format requirements.*


- **./SONAR-code/result/**

  *These files store the output results*.

  Proportions of all cell types in each spot (**SONAR.results.txt**);

  Spatial pie plots (**pie.pdf**);

  Spatial distribution of specific cell types proportion (**abs_prop.pdf** / **scaled_prop.pdf**);

  Colocalization(correlation along the spatial) for pairs of cell types (**colocalization.pdf**).


- **./SONAR-code/data/**

  *no need to operate.* 

  *These files store the preprocessed data that delivered to SONAR*.

