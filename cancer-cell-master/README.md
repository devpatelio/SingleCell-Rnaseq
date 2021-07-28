# LYMPHOMA MICROENVIRONMENT ANALYSIS - README

Welcome to my sc-RNAseq code, designed to be a straightforward and user-friendly pipeline that helps interpret and analyze key cell environments in a tissue sample. 
It's 


### Running into Errors / Troubleshooting
1. Python Package Error 
    If running in Jupyter Lab or notebook, you can take a series of steps to
    make sure that your libraries are properly installed. For one, make sure 
    that your conda environment has all of the required packages. 

    To check your libraries, enter the following in your terminal or folder:
    ```

    >> conda env list
    ##If your environment is there, you can just delete it and rewrite it with an updated yml file

    >> conda env remove --name _name_of_conda_environment_
    >> conda env create -f sc-rnaseq.yml
    ##Once installed, you can activate the environment in terminal
        ##Note that the name of the environment is in the yml file, so 
        ##you can go in and manually edit that name to activate it below


    >> conda activate sc-rnaseq
    ##Now, you can switch the jupyter kernel to the conda environment
    ##If all else fails, install all error-libraries in the conda env 
    
    >> (sc-rnaseq) conda install _package_name_
    >> (sc-rnaseq) pip install _package_name_ 
    ```


2. Dataset formatting errors
    If working with 10x Genomics Data, you must make sure that it follows the specified format in 
    [10xDatasets](https://support.10xgenomics.com/single-cell-gene-expression/datasets). If you would like to 
    test a dataset, you can go into any of these single-cell expression datasets and do the following: 
    
    - [ ] Click on *Gene Expression - Feature/cell matrix* (raw or filtered)
    - [ ] Click on *Gene Expression' - Feature/cell matrix HDF5* (raw or filtered)
    - [ ] OPTIONAL Download any metadata such as cells per tag, dataset summary, etc.

    - [ ] Unzip tsv.gz files, if ERROR: 79, install [The Unarchiver](https://theunarchiver.com/) 
          and use the application to unzip the files
    - [ ] You can now input the absolute path of the folder that contains the matrix.mx and other
          files as the downloadable root path


3. Divide By 0s: 
    Although this should never be the case thanks to error handling, check if your mitochondrial gene
    plot is at 0 and if the error is being throw in cell 7. If so, you must make sure that your Error 
    Threshold Input is greater than 0.

    Otherwise, you must make an informed visual choice to determine the most effective way to filter 
    the data. You can always visit it back to make changes. 


4. Stuck in an awkward loop
    In some cases, this can probably be a memory error or if you've modified the code and certain flags are
    misplaced. To solve this, add a `break` to any of the loops or just remove them if you do not want a dynamic
    integration. 

    For other instances, just interrupt the cell and rerun it by adding manual inputs to the parameters -> you 
    can trace this back based on where the input variable is used. 


5. TSNE, UMAP, or PCA Stuck
    In most cases, this will take more than 6-10 seconds depending on your computer internals. To maximize
    your calculations, you can add the flag `n_pcs=14` to `sc.tl__ ()` change the number of pcs being used.

    In other cases, it would be a good idea to reduce the total data you've passed through or it's just too 
    big. You can track how much is completed by changing the learning_rate parameter, but should always be 
    between 100 to 1000. Refer to [Documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.tsne.html)
