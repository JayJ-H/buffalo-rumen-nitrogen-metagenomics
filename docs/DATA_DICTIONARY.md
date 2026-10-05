# Data dictionary and provenance

## Samples

`sample_id`: laboratory identifier of the rumen-fluid sequencing sample.
`animal_id`: identifier of the individually housed buffalo; read as a string.
`diet_group`: A, B, C, D or E, six animals per group.
`initial_BW_kg`: baseline body weight in kilograms.

`metadata_correction_log.csv` records two confirmed transcription corrections:
animal 2319 belongs to B and animal 2353 belongs to C. Neither correction
represents an experimental reassignment. The analysis inputs already use the
confirmed groups. Original records were not overwritten.

## KO data

`data/KO_reconstructed_gene_relative.tsv` has KOs in rows and the 30 rumen
samples in columns. Values sum relative abundance of genes with a single KO
assignment. They are not raw read counts. The scale is inherited from the
gene-level abundance table; all-gene column totals are approximately 1,000,000.
Exact sample totals are supplied in `gene_definition_evidence.json`.

The supplied KO matrix is a processed analysis input, not a full gene catalog.
The analysis package does not rebuild sequence assembly, gene prediction or
read-to-gene abundance estimation. Those processing steps are described in
the manuscript methods; raw sequencing data use the NCBI accession route.

`function_definitions.json` supplies the fixed 20 base functions, three ratios,
four carbon-processing sets, EC categories, resampling counts, zero-replacement
rules and multiple-testing family sizes. Memberships preserve the definitions
used in the manuscript. KO identifiers within each function are unique;
different functions can share KOs.

## Scores and statistics

`scores_per_animal.csv` contains natural-log scores for 23 responses,
unlogged fractions (`__fraction` suffix), all-gene and KO-gene denominators,
and the natural log of the KO/all-gene abundance ratio.
`within_KO_scores_per_animal.csv` uses the KO-gene denominator.
`carbon_pair_scores.csv` contains the log nitrogen response and log carbon
predictor for each pair, in sample order. `carbon_pair_membership.csv`
documents overlap removal and zero replacement.

`P_HC3` and `q_HC3` are the primary raw and BH-adjusted probabilities.
`P_permutation` and `q_permutation` refer to sensitivity tests.
`Beta` is the regression coefficient in natural-log score units.
`Partial_r` is the correlation of nitrogen and carbon residuals after adjustment.
Confidence intervals use the model residual degrees of freedom.
`LOO` files contain leave-one-animal-out estimates.
Permutation and bootstrap index files use R's one-based row numbering.
`CV` in `functional_description.csv` is a fraction, whereas `CV_percent` in
Table 2 is that value multiplied by 100.

## Tables and figures

`tables/chinese_draft/Table_1.csv` contains diet ingredients and calculated
nutrient levels on a dry-matter basis (%), using Chinese row labels.
`tables/Table2_enzyme_abundance.csv` contains the four enzyme gene collections,
detected KO counts, mean abundance as percent of all genes, CV (%) and
detected-animal counts. It is calculated from `functional_description.csv`.

Figure panel source files identify every plotted quantity. Supplementary
violin panels use reversible power-of-ten display scaling, recorded in
`FigureS4_display_scaling.csv`. This scaling is for plotting, not model fitting.
