# COVID-19 Project

This repository contains source code used in analysis of patients infected by SARS-CoV-2.

Manuscript: "in progress". 

By: The Independent Clinical Epigenetics Laboratory, Poland

### Start
To start work with code make sure that R>=4.0.3 (https://cran.r-project.org/), 
Python >=3.8, <3.10 (https://www.python.org/), pip (https://packaging.python.org/tutorials/installing-packages/) and poetry (https://python-poetry.org/) are installed.

To download repository:
	
		git clone https://github.com/EpiGenMed/covid-project.git

### R scripts
        
Repository contains 1 file:
        
* **ModifiedRefBase.R** -> containing  **modified_refBase** function what is a copy of refBase (https://rdrr.io/bioc/ChAMP/man/champ.refbase.html) from ChAMP package. Note, that part of code used to estimate white-blood cell fractions was re-implemented using RPC method from EpiDISH (https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html) package.

* **install_packages.R** -> file containing requirements.
        
Before run install required packages type:
        
        Rscript install_packages.R
        

#### Examples

##### **ModifiedRefBase.R**

Input *.csv* file must be a 2-dimensional frame containing beta-values, with **CpGs** as rows and **samples** as columns. 

Function modified_refBase return list of three obejcts: 

1) Estimated white blood cell fractions for each sample
2) Corrected mynorm
3) Estimated white blood cell fractions for each sample after correction for WBC fractions

        setwd(<path_to_dir>)
        source(<path_to_file.R>)
        mynorm <- data.table::fread(<path_to_mynorm_file>, data.table = F)
        mynorm <- data.frame(mynorm, row.names = 1)

        cfc <- modified_refBase(mynorm)

### Python notebooks

This part of repository contains notebooks and a source code used in the main part of the analysis.

To make code transferable we used *Poetry* (https://python-poetry.org/) dependency menager. To create virtual environment and install all required dependencies type:

        pip install poetry [if needed] 
        cd <project directory> && poetry install 
        
Then to run notebook use:

        poetry run jupyter-lab
        

*src/* directory contains files:

   
    * figures.py -> with function to generate plots (col_pallete.py contains colors used per sample group)
    
    * utils.py -> additional functions.py
    

*statistics/src* directory contains file:

    * stats.py -> class implementing methods to identify DMPs
