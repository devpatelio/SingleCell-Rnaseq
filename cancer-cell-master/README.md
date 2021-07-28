# :microscope: LYMPHOMA MICROENVIRONMENT ANALYSIS - README

Welcome to my sc-RNAseq code, designed to be a straightforward and user-friendly pipeline that helps interpret and analyze key cell environments in a tissue sample. 
The example provided here is meant to give you an intuition of how and why I designed this code and the different biological underlying mechanisms at play
that can be interpreted through microenvironment analysis. 

> This informatics pipeline primarily utilized scanpy, a powerful bioinformatics library built by the [Theis Lab](https://github.com/theislab/scanpy). 

## Project Outline
This project is designed to provide a complete overview of the entire sc-RNAseq pipeline. The tool has become widely popular for a variety of reasons including: 
- mapping genotypes to phenotypes is challenging in biology and medicine → powerful method is currently transcriptome analysis, but transcriptome information in a cell reflects activity of a small subset of genes
body's cells each express unique transcriptome, and increasing evidence shows that gene expression is heterogeneous even in similar cell types
- majority of transcriptome analysis is based on assumption that cells from given tissue are homogeneous and these studies likely miss important cell-to-cell variability
- most biological processes are stochastic, thus we need more precise understanding of transcriptome in individual cells for finding their role in cell functions and understanding how gene expression can promote certain genes
using cDNA of individual cells, single-cell RNA sequencing was born → this helped provide high-resolution views of single-cell heterogeneity on global scale and allowed for finding differences in gene expression between individual cells has potential to identify rare populations that cannot be detected from analysis of pooled cells
    - eg. finding and characterizing outlier cells within population has potential implications for furthering understanding of drug resistance and relapse in cancer treatment with bioinformatics pipelines, can now learn more about highly diverse immune cell populations in healthy and diseased states
- scRNA-seq is being utilized to delineate cell lineage relationships in early development, myoblast differentiation, and lymphocyte fate determination

These key features make it very promising, primarily because it has been very useful for understanding single-cell identity, involving the different gene expressions to their response and interaction with external stimuli such as neighbouring cells or diseases like tumors, which is the primary focus of the demo. 

This pipeline is an interactive, user-friendly visualization and analysis toolset which can be used to better analyze and interpret single-cell RNA sequencing data of a set of tissues, and then analyze different features of these cells such as their gene expression activity (transcriptome), protein interactions, and cell localization. 

The program is designed to visualize and analyze sc-RNAseq data which ultimately looks at gene and cell activity for the individual types of cells in a microenvironment. This dataset is based off of Hodgkin's Lymphoma samples sequenced by 10x Genomics, and analyzed in the pipline which looks at a tumor microenvironment and the different kinds of cells there are (T-cells, CD4 vs CD8+ T Cells vs NKCs in a tissue, FOXP3+/CD3+/SIGLEC 1 T-cell -> the primary goal will be to find patterns and essentially ID different kinds of cells and try to identify whether or not the immune system of the body is actually killing off the tumor. 

Further Readings: 
> [Analysis of scRNA-seq Data](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html)
> [Best Practices in scRNAseq Analysis](https://www.embopress.org/doi/full/10.15252/msb.20188746)
> [StatQuest Playlist](https://www.youtube.com/user/joshstarmer)
> [scanpy Documentation](https://scanpy.readthedocs.io/en/stable/usage-principles.html#workflow)

## Biological Significance

**Brief history of transcriptional profiling**: 

- 2000s, microarrays enabled systematic measurements of transcriptional changes → mRNAs labeled, loaded onto chip where each spot on chip was coated to catch mRNA of specific gene → color of spot indicates ratio of transcript
- Bulk RNA sequencing made it much more feasible which transcriptional differences could be measured
    - quantification of gene expression
    - comparitive transcriptomics
- limitations involve the inability to resolve heterogeneity
- single cell also gives you dynamics of gene expression and cell identity, primarily because it can help you look at cell topology
- goals of scRNA-seq are 1) measure distribution of expression levels for each gene across population of cells 2) measure transcirptional differnces across and within groups of cells 3) resolve single-cell heterogenity
- applications involve characterization of cell and transcriptional composition of tissues (cell type changes and their correlatoin)
- evaluate developmental processes (see expression changes while undergoing differentiation)
- evaluate how cancer evolution leads to mechanisms of therapeutic resistance

**Droplet-based approach**: 

1. Tissue disassociated and samples created to separate individual cells, need to make sure tissue maintains integrity to prevent digestion of extracellular layers
2. samples collected, separated into individual droplets (use cytometry to evaluate cell viability through FACS to measure cell composition by pre-binding cell with fluorescent or antibodies to find expectation of proportion of cell types)
3. oil and water reaction has a single bead with barcode that bind to the specific mRNA of a single cell + done through cell where water and oils are displaced → you have a poisson distribution of beads and cells (proportion of concentration is based on concentration of cells and beads)
4. cells are lysed, reverse transcriptase synthesizes cDNA and applied molecular identity, use pCR and then add sequencing adapters for sequencing and assembly (3'-5')
    - lysis part is flawed right now because mRNA molecules don't always bind to the beads, RT enzymes are known to be inefficient, and RT + PCR depend on specific transcript, making this more variable
- onyl about 5% of mRNA is captured → UMI: total count of unique molculare identifies which are labelled on each mRNA with a unique barcode at the RT step
    - UMI reduces amplification error → evaluation of performance found that 10x is one of the better ones → SUPeR-seq wth Accuracy of 0.95, 4 million copies of mRNA sequence in every cell

**Limitations:** 

- low capture probability and data sparsity is a big problem → you can use BDRhapsody to limit the scope of experiment to a few hundred mRNAs
- inability to analyze full-length transcriptomes → can use smartSEQ2
- inability to resovle spatial information → where computation comes into play
- integrating with measurements in Microscopy, FACs, total_Seq, hashtag tech, etc.

- preparing scRNAseq data for clustering → preprocessing: align and count UMIS
- finding feature selection and do dimensionality reduction → you can visualize hereogeneity with t-SNE and UMAP
- **HVGs for preprocessing**:
    - distinguish technical noise by looking at distribution of population of cells → HVGs are identified because they are interesting biological markers
    - PCA is then applied as a result
        - visualized using UMAP and t-SNE
- Quality control → normalization → feauture selection → dimensionality reduction → cell-cell distances → unsupervised clustering

**Computational challenges in scRNA-seq:** 

- computational pipelines for handling raw data files is limited → not many tools and is still in infancy
- first step is pre-processing data → once reads are obtained from well-designed scRNA-seq experiments, quality control performed
- read alignment is next step and tools available for this procedure are used for bulk RNA and can be used here → when adding transcripts of known quantity and sequence for calibration and QC, low-mapping ratio of endogenous RNA to spike-ins is indication of low -quality library caused by RNA degradation

    ![https://s3-us-west-2.amazonaws.com/secure.notion-static.com/168c4f07-b933-42da-b0e9-c5bddc0eaa78/Screen_Shot_2021-07-16_at_12.07.13.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/168c4f07-b933-42da-b0e9-c5bddc0eaa78/Screen_Shot_2021-07-16_at_12.07.13.png)

- following alignment, reads allocated to exonic, intronic, or intergenic features using transcript annotation in format **General Transcript**
    - only reads that map to exonic loci with high mapping quality considered for generation of gene expression matrix → N x m where N=cells, m=genes
- there is a presence of zero-inflated counts due to dropout or transient gene expression, which requires normalization to remove cell-specific bias
    - read count of gene in each cell expected to be proportional to gene-specific expression level and cell-specific scaling factors

**My Project:** Single-Cell RNA Sequence Analysis (Downstream) w. Interactive UI for Hodgkin's Lymphoma, Dissociated Tumor: Targeted, Pan-Cancer Panel
[Single-Cell Transcriptome Analysis Reveals Disease-Defining T-cell Subsets in the Tumor Microenvironment of Classic Hodgkin Lymphoma](https://pubmed.ncbi.nlm.nih.gov/31857391/)
- cHL characterized by extensive microneivornment composed of different types of noncancerous normal immune cells → severla types of T cells, B cells, eosinophils, and marcophages, and rare populations of clonal malignant Hodgkin and Reed-Sternberg cells
    - some findings support concept that HRS cells recruit immune cells to form tumor-supporting, regulatory tumor microenvironment (TME) with limited antitumor activity in cHL
    - complex interactions between HRS cells and TME remains partially understood → looking at symbiotic cellular cross-talk may lead to development of novel biomarkers and therapeutic approaches
- immune-checkpoint inhibitors like PD-1 have shown dramatic efficacy in relapsed or refractory cHL, with overall response rate of 65-87% and durable remissions of 1.5 years
    - remains unclear which cells are most important targets of immune-checkpoint inhibitors and which components are most relevant for immune-escpae phenotype in cHL
- goal of this project is to characterize immune phenotype of the TME in cHL and identify important associations between immune cell types and respective clinilca outcome
- key differentiator between lymphomas and tumors is that they are dervied from lymphocytes that progessionaly interact with other immune cells in ecosystem of microenvironment


### System Requirements
Make sure you have Python3.5 or greater installed, along with anaconda setup with pip as your default package manager. 
Also requires that you have experience with python and package installation + some basic biology background to understand
some of the biological comparisons drawn here with the data. 

Package specifications and versions in the [.yml file](sc-rnaseq.yml) for reference.

### The Dataset and Exploration
The dataset is titled Hodgkin's Lymphoma, Dissociated Tumor: Targeted, Pan-Cancer Panel which is a [Single Cell Gene Expression Dataset by Cell Ranger 4.0.0](https://support.10xgenomics.com/single-cell-gene-expression/datasets/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer)

The sample is collected from a disassociated lymph node tumor of a 19-year-old male. The sample was obtained from 10x Genomics and the whole transcriptome in the dataset is generated with Chromium GEM Single Cell 3' Reagent Kits. There are approximately 11,332 reads per cell (not necessarily important unless you're interested in genome assembly metrics). Out of the tumor microenvironment, 3049 cells were detected with an 97.2% confidence score for the reads being mapped to the targeted transcriptome. The total targeted genes detected are 1165.

Download links: 
[Feature / cell matrix HDF5 (filtered)](https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer_filtered_feature_bc_matrix.h5)
[Feature / cell matrix (filtered)](https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer_filtered_feature_bc_matrix.tar.gz)


### The UI and Getting Started
**File Downloads and Setting up Anaconda**
You can download the files using the links above, and save them to the directory we will create below. If you do not have anaconda already installed, you can refer to the [documentation](https://docs.anaconda.com/anaconda/install/index.html). Once you've set up anaconda, enter the following lines in your terminal: 

```
>> mkdir cancer-target-main
>> cd cancer-target-main

##You can manually install this by clicking the link, but for anyone with terminal experience
>> wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer_filtered_feature_bc_matrix.h5
>> wget https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer/Targeted_NGSC3_DI_HodgkinsLymphoma_Pan_Cancer_filtered_feature_bc_matrix.h5

##You can also unzip manually, but this is also convenient
>> unzip cancer-target-main/filtered_feature_bc_matrix.h5
```

From here, you can download the `targeted-cancer-lymphoma.ipynb` file in the same directory. Then, run the notebook in VSCode or jupyter lab. To install jupyter lab, visit this [link](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html)

** Advanced Startup from Github **
You can simply clone the repository and run the primary *targeted-cancer-lymphoma.ipynb* file. 


> Now you're ready to follow along and understand your chosen tissue microenvironment. 

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
