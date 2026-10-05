# Rumen nitrogen and amino-acid metabolism in growing male buffaloes

Processed data and code supporting the manuscript
"Metagenomic characteristics of rumen nitrogen and amino-acid metabolism
in growing male buffaloes".

This package contains processed data for 30 end-of-trial rumen-fluid samples from
individually housed male buffaloes, five diets with six animals each, and a
150-day feeding trial. The individual buffalo is the experimental unit.
The code starts from the processed KO abundance matrix, not raw FASTQ files.

## Contents

- `data/`: KO abundance, all-gene denominators, fixed function definitions,
  and EC/KO memberships.
- `data_processed/`: sample identifiers, diet group and initial body weight.
- `analysis/`: offline score preparation, R statistics, R plotting and
  independent Python/R validation.
- `results/nitrogen_function/`: per-animal scores and complete statistical outputs.
- `source_data/`: data underlying each figure panel.
- `figures/nitrogen_tahp_revised/`: two main and four supplementary figures,
  with vector PDF/SVG and raster PNG/TIFF exports.
- `tables/`: manuscript table data.
- `environment/`: software versions used for verification.
- `MANIFEST.csv`: file sizes and SHA-256 checksums of the packaged contents.

## Reproduction

Python 3 with the packages in `requirements.txt`, and R 4.5.3 with the packages
listed in `analysis/install_dependencies.R`, are required. Arial was used for
the supplied plots; install that font for matching labels. R plotting needs
Cairo support. Install dependencies once:

```sh
python -m pip install -r requirements.txt
Rscript analysis/install_dependencies.R
```

From the repository root, run in this order:

```sh
python analysis/prepare_nitrogen_function.py
Rscript analysis/run_nitrogen_function.R
python analysis/validate_nitrogen_function.py
Rscript analysis/render_release_figures.R
Rscript analysis/validate_release_figures.R
python analysis/clean_release_metadata.py
```

The last command strips document metadata and normalizes archive timestamps;
it does not change numeric data or plotted pixels. It refreshes checksums.
An existing R library may be selected using `TAHP_R_LIBRARY`.
No network calls are made by the analysis itself.

## Abundance and models

The KO matrix sums gene-level relative abundance for genes assigned to a
single KO. Function memberships are fixed in `data/function_definitions.json`.
Duplicate KOs within a function are counted once. The primary denominator is
the sum of relative abundance across all catalog genes for that sample;
figure values are percentages of that denominator. The alternative denominator
is the total abundance of single-KO-assigned genes.

Scores use natural logarithms. When a zero occurs, it is replaced with half
the smallest positive fraction for that function. The three ratios use
differences of log scores. KO contributions use summed abundance across all
30 animals divided by the summed abundance of the corresponding function.
Functions may overlap and must not be summed as disjoint categories.

Diet models adjust for initial body weight, use HC3 covariance and a
four-degree-of-freedom Wald F test. BH families contain 23 overall diet tests
and 230 pairwise contrasts. The 80 carbon-function models adjust for diet,
initial body weight and log KO annotation load. Each pair excludes shared
KOs. Primary inference uses HC3 with BH correction across 80 pairs.
Residual permutation, leave-one-animal-out and diet-stratified bootstrap
analyses are sensitivity checks, not replacements for the primary tests.
Seeds and resampling counts are supplied in the function definitions.

## Data access

The manuscript lists NCBI BioProject `PRJNA1321136` for sequencing data.

This package contains study-derived data and function identifier mappings.
It does not contain the full KEGG database, downloaded database entry records,
third-party article PDFs, raw sequencing reads or private project notes.

## Version and reuse

Version `v1.0.0`, prepared 2026-10-05.
Project-authored code uses MIT; study-derived data, figures and documentation
use CC BY 4.0. See `LICENSE` and `LICENSES/README.md` for license scope.
