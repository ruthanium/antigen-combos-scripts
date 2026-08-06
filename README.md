# Antigen Combos

Patient-aware single-cell discovery and validation of combinatorial
tumor-associated antigen gates. The primary workflow supports:

- biomarker-anchored partner discovery, including measured protein-positive
  cell calls or an explicitly labeled RNA surrogate;
- patient-separated discovery and held-out testing;
- AND, OR, NOT, and three-input Boolean antigen gates;
- malignant/normal epithelial coverage and patient prevalence;
- cell-cluster decoding and marker expression-level reports;
- independent validation in the ICBI colorectal cancer atlas; and
- corrected comparison with the original PCASA gate calculations.

The main CRC results and validation reports are versioned under `results/`.
Downloaded single-cell matrices and API caches remain excluded because of
their size.

This repository builds on the scripts used for the Dannenfelser et al. paper,
_Discriminatory power of combinatorial antigen recognition in cancer T cell
therapies_. The original TCGA/GTEx bulk-expression workflow is retained as
legacy/reference functionality, while the newer development focuses on
patient-aware single-cell analysis.

### Citation
> Discriminatory power of combinatorial antigen recognition in cancer T cell therapies.
Dannenfelser R, Allen G, VanderSluis B, Koegel AK, Levinson S, Stark SR, Yao V, Tadych A, Troyanskaya OG, Lim WA. Cell Systems. 2020. [https://doi.org/10.1016/j.cels.2020.08.002](https://doi.org/10.1016/j.cels.2020.08.002)
<!-- (DOI badge for later?[![DOI](https://zenodo.org/badge/126377943.svg)](https://zenodo.org/badge/latestdoi/126377943)) -->

### Usage
Due to size contraints, this repository contains expression data for a subset of the genes analyzed in the paper. The instructions below will work off this data. However it is easy to extend these scripts by processing an input gene expression matrix using the same format of the sample matrix `data/test-normalized-matrix.txt`, with genes as columns and samples as rows. In addition to genes, be sure to have a `dataset` column with unique identifiers, a `type` column with values of either `cancer` or `normal`, and a `tissue.cancer` column with a more specific sample label such as the type of cancer (e.g., `Glioblastoma Multiforme`) or the tissue type (e.g., `Adipose Tissue`).

To build the version used in the paper, download the 2015 release of RNAseq data from [TCGA](https://portal.gdc.cancer.gov/) and version 7 of [GTEx](https://www.gtexportal.org/home/). Correct for batch effects using your favorite normalization scheme (we used COMBAT) and make sure the expression values are in log_2 space.

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

### Python pipeline

The computational R stages have also been ported to Python. The Python port
uses CSV intermediates so that results can be inspected without R:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements-python.txt

python src/CreateSketches.py
python src/CalculateClusteringScores.py
python src/Evaluate.py
python src/Analyze.py
```

The complete doubles and triples search can be expensive. To validate the
pipeline quickly, limit doubles and triples to 50 combinations per cancer:

```bash
python src/CalculateClusteringScores.py --test-run
python src/Evaluate.py
python src/Analyze.py
```

For a still smaller smoke test, process one cancer:

```bash
python src/CalculateClusteringScores.py --test-run --cancer "Acute Myeloid Leukemia"
```

Python outputs are written to `results/`:

- `singles.csv`, `doubles.csv`, and `triples.csv`: clustering metrics
- `*.total.csv`: clustering metrics plus held-out precision, recall, and F1
- `top_*.csv` and `top_combinations.csv`: ranked, diversity-filtered results
- `performance_summary.csv` and `f1_by_gate_size.png`: summary analysis

The custom Davies–Bouldin calculation and the AND-gate evaluation logic are
ported directly. Decision stumps use scikit-learn with settings matching the
important `rpart` defaults, but exact thresholds can differ slightly between
the two tree implementations.

### Biomarker-anchored single-cell partner discovery

`FindBiomarkerPartners.py` accepts an anchor biomarker and finds pairs of RNA
partners for a three-input AND gate. It reports patient-level tumor and normal
coverage and decodes each gate across the published epithelial cell clusters.
Partner selection uses discovery patients only; complete patients are held out
for testing. The top held-out gates are then tested in the independent ICBI CRC
atlas after excluding the overlapping Pelka/GSE178341 dataset. Patient
bootstrap confidence intervals are included for held-out coverage.

```bash
python src/FindBiomarkerPartners.py \
  --biomarker CEACAM5 \
  --extra-genes GPR35 \
  --top-partners 20 \
  --max-normal 0.01 \
  --output-dir results/biomarker_partners/CEACAM5
```

By default, detected biomarker RNA is explicitly treated as a surrogate.
Measured protein-positive cells can instead be supplied as a delimited file
with a `cell_id`, `barcode`, or `cellID` column and an optional numeric column:

```bash
python src/FindBiomarkerPartners.py \
  --biomarker CEACAM5 \
  --protein-calls CEACAM5_protein_calls.tsv \
  --protein-threshold 1.0
```

Outputs include `partner_ranking.csv`, `two_partner_gates.csv`,
`two_partner_cluster_decode.csv`, `marker_expression_by_split.csv`,
`heldout_marker_expression_by_cluster.csv`,
`external_marker_expression.csv`, and `run_summary.csv`.

### User Agreement
This code is released under a [Clear BSD License](https://github.com/ruthanium/antigen-combos-scripts/blob/master/LICENSE) and is intended for research/academic use only. For commercial use, please contact: Laurie Tzodikov (Assistant Director, Office of Technology Licensing), Princeton University, 609-258-7256.
