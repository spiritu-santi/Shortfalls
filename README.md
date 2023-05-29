# Digital Knowledge Shortfalls

Files needed to run the results from:

- Ramírez-Barahona S, Cuervo-Robayo AP, Magallón S. 2022. Assessing digital accessible botanical knowledge and priorities for (re)exploration and discovery of plant diversity across Mesoamerica.

## R script
The repository follows the file system of the R project (NewPhytol_KEW.Rproj). All directories are needed to properly run the functions.  
All analyses were performed in R.  
The script is fully functional, but some of the functions are not yet fully automated, some options are 'hard-coded', and pipelines have not been benchmarked for efficiency.  

- The 'data' directory contains some of the datasets (< 100MB) necessary to conduct the analyses (other data are generated within the script).  
- The 'code' directory contains the R script to perform the analyses and plot the results.  
- The 'figures' direrctory contains output plots generated within the script (NOTE: the final figures were obtained by editing figures in Adobe Illustrator).
- The 'output' directory contains some of the output files (< 100MB) of the code (other outputs are generated within the script).
- The 'interim' directory is empty and intermediate outputs will be saved here.

The R code was last updated March 10, 2023 (KewCleaning.R) and contains all the functions necessary to process the geographic data as obtained
from GBIF and harmonize it against Kew's World Checklist of Vascular Plants (https://powo.science.kew.org/).

Code and data are available at Zenodo:
[![DOI](https://zenodo.org/badge/194142746.svg)](https://doi.org/10.5281/zenodo.7982903)

When using the code please cite:
- Ramírez-Barahona S, Cuervo-Robayo AP, Magallón S. 2022. Assessing digital accessible botanical knowledge and priorities for (re)exploration and discovery of plant diversity across Mesoamerica.
