### _In silico_ Screening of Antigen Combinations
This is a repository for the scripts used to do the analysis in the Dannenfelser et al. paper, _Discriminatory power of combinatorial antigen recognition in cancer T cell therapies_. Top predicted combination can be viewed on our [antigen explorer](http://antigen.princeton.edu) webserver.

### Citation
> Discriminatory power of combinatorial antigen recognition in cancer T cell therapies.
Dannenfelser R, Allen G, VanderSluis B, Koegel AK, Levinson S, Stark SR, Yao V, Tadych A, Troyanskaya OG, Lim WA. FULL CITATION TO COME

<!-- (DOI badge for later?[![DOI](https://zenodo.org/badge/126377943.svg)](https://zenodo.org/badge/latestdoi/126377943)) -->

### Usage
1. Need to gather TCGA and GTEx RSEM counts
2. Combine using COMBAT or favorite batch correction scheme (matrix of sampleIDs by genes)
   - first column (dataset (dataset specific identifier), last columns: batch (GTEx or TCGA), type (tissue or cancer), tissue.cancer( specific sample type e.g., colon cancer or adipose tissue))
   - should be log transformed values
3. Build sketches with `CreateSketches.py` (need to make sure env correct - python 3, requirements.txt etc)
4. load all R functions - `Functions.R`
5. Calculate clustering-based scores for singles, doubles, and triplets `CalculateClusteringScores.R`
6.
