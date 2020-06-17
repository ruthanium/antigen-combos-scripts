### _In silico_ Screening of Antigen Combinations
This is a repository for the scripts used to do the analysis in the Dannenfelser et al. paper, _Discriminatory power of combinatorial antigen recognition in cancer T cell therapies_. Top predicted combinations can be viewed on our [antigen explorer](http://antigen.princeton.edu) webserver.

### Citation
> Discriminatory power of combinatorial antigen recognition in cancer T cell therapies.
Dannenfelser R, Allen G, VanderSluis B, Koegel AK, Levinson S, Stark SR, Yao V, Tadych A, Troyanskaya OG, Lim WA. FULL CITATION TO COME

<!-- (DOI badge for later?[![DOI](https://zenodo.org/badge/126377943.svg)](https://zenodo.org/badge/latestdoi/126377943)) -->

### Usage
Due to size contraints, this repository contains expression data for a subset of the genes analyzed in the paper. However it is easy to extend these scripts by processing an input gene expression matrix using the same format of the sample matrix `data/test-normalized-matrix.txt`, with genes as columns and samples as rows. In addition to genes, be sure to have a `dataset` column with unique identifiers, a `type` column with values of either `cancer` or `normal`, and a `tissue.cancer` column with a more specific sample label such as the type of cancer (e.g., `Glioblastoma Multiforme`) or the tissue type (e.g., `Adipose Tissue`).

To build the version used in the paper, download the 2015 release of RNAseq data from [TCGA](https://portal.gdc.cancer.gov/) and version 6 of [GTEx](https://www.gtexportal.org/home/). Correct for batch effects using your favorite normalization scheme (we used COMBAT) and make sure the expression values are in log_2 space.

#### Instructions
1. Build the geometric sketches needed for the training and test sets:
   - setup a virtual environment and install relevant python packages (here we used conda)
   ```bash
   conda create -n antigen python=3.7
   conda activate antigen
   pip install -r requirements.txt
   ```
   - run `CreateSketches.py`
   ```bash
   python src/CreateSketches.py
   ```
   - once complete, sketch files should be constructed in `results/sketches` and the virtual environment is no longer necessary.
2. The remaining steps will be in R. First, load all R functions in the  `Functions.R` file.
3. Calculate clustering-based scores for singles, doubles, and triples with `CalculateClusteringScores.R`. Depending on your computational resources running on the sample file may take a while to complete. Considering setting the `tr` variable to `True` to further reduce the number of double and triple combinations computed for testing purposes.
4. Run the decision tree evaluations in `Evaluate.R` to calculate performance metrics (F1, precision, recall). Like cluster score calculations, this script may take a long time to run based on your computing resources and the size of the starting dataset.
5. Use `Analyze.R` to process the top results and conduct the analysis carried out in the paper.

### User Agreement
This code is free to use and modify for non-commerical use.
