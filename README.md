# COVID-19 Project

This repository contains source code used in analysis of patients infected by SARS-CoV-2.

Manuscript: "Coherent response of blood methylome to SARS-CoV-2 infection between three independent populations of patients". 

Team: The Independent Clinical Epigenetics Laboratory, Poland

### Start
To start work with code make sure that R>=4.0.3 (https://cran.r-project.org/), 
Python >=3.8, <3.10 (https://www.python.org/), pip (https://packaging.python.org/tutorials/installing-packages/) and poetry (https://python-poetry.org/) are installed.

To download repository:
	
		git clone https://github.com/EpiGenMed/covid-project.git

### R scripts
        
Repository contains 3 files:
        
* **FAnalysis.R** -> containing **enrich** function implementing gene-set-enrichment-analysis using methylGSA (https://bioconductor.org/packages/release/bioc/html/methylGSA.html) package.

* **RawDataProcessing.R** -> containing code used to process *IDAT files using ChAMP (https://www.bioconductor.org/packages/release/bioc/html/ChAMP.html) package.

* **ModifiedRefBase.R** -> containing  **modified_refBase** function what is a copy of refBase (https://rdrr.io/bioc/ChAMP/man/champ.refbase.html) from ChAMP package. Note, that part of code used to estimate white-blood cell fractions was re-implemented using RPC method from EpiDISH (https://www.bioconductor.org/packages/release/bioc/html/EpiDISH.html) package.

* **install_packages.R** -> file containing requirements.
        
Before run install required packages type:
        
        Rscript install_packages.R
        

#### Examples

##### **FAnalysis.R**

Input *.csv* file must contain two columns **CpG** with CpG ID and **p-value** with corresponding value. Input file is an output from StatsAnalysis class described below.

	setwd(<path_to_dir>)
        source(<path_to_file.R>)
        enrich(result_dir=<path_to_output_directory>,
               input_data=<path_to_input_csv_file>,
               ) # Other arguments described in methylGSA package documentation


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
        

Additionally *src/* directory contains two files:

    * stats.py -> class implementing methods to identify DMPs.
   
    * enrichment_analysis.py -> class to perform enrichment analysis of identfied DMPs.
    
    
#### Examples


##### **stats.py**

Implementation of statistical analysis process to identfiy DMPs.
    
    import pandas as pd
    from stats import StatsAnalysis
    
    control_mynorm = pd.read_csv(<path_to_csv_file>, index_col=0) # Control samples
    target_mynorm = pd.read_csv(<path_to_csv_file>, index_col=0) # Target samples
    
    epic = pd.read_csv("<path_to_EPIC_manifest>", index_col=0, low_memory=False)
    # Load EPIC manifest 
    
    stats = StatsAnalysis(target_mynorm, control_mynorm, epic) # Initialize object
    # stats.extract_probes(regions_ra=<region>) # Optional extract only CpGs in specific regions for example: regions_ra="TSS1500|TSS200" will extract probes in 
    # TSS1500 or TSS200
    results = stats.run()
    
    results.to_csv(<path_to_save_file>)
    

##### **enrichment_analysis.py**

Enrichment analysis for set of CpGs identifed using StatsAnalysis implementation.

    import pandas as pd
    from src.enrichemnt_analysis import EnrichmentAnalysis
    
    bg = pd.read_csv(<path_to_bg_file>, index_col=0) # Background is set of CpGs from mynorm file used to identfied DMPs 
    report = pd.read_csv(<path_to_report_file>, index_col=0) # Report file containg DMPs is filtered output from StatsAnalysis
    
    ea = EnrichmentAnalysis(
                            mynorm=bg,
                            report=report,
                            manifest_path=<path_to_EPIC_manifest>,
                            )
    ea.prepare_bg()
    ea.calculate_frequency()
    ea.estimate_fc()
    ea.vis() # or ea.vis(<path_to_output_dir>) to export results
