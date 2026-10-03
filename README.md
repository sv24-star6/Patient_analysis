# Patient analysis: neem limonoids in oral cancer

Scripts and processed inputs for patient expression and survival analyses.

## Run
Install Python dependencies:
```
pip install pandas numpy scipy statsmodels lifelines matplotlib
```
From the extracted folder:
```
python analyse.py
python plot_figure9_expression.py
python plot_figure10_survival.py
```
The analysis script produces integrated_target_evidence.csv (the full-precision basis for Table 5). The plotting scripts write PNG (300 dpi) and vector PDF figures into figure9_expression/ and figure10_survival/.

## Data and interpretation
GSE30784: 167 tumours and 45 normal samples in the primary comparison; dysplasia excluded. GSE9844: 26 tumours and 12 normal controls. Processed gene expression was summarised by the median across eligible probes. Expression tests use Welch tests with BH correction across 14 targets per cohort. Replication requires matching direction and corrected P < 0.05 in both cohorts.

TCGA-HNSC PanCancer Atlas expression and clinical data were restricted to oral sites using anatomical annotations. Adjusted Cox models include age, sex and stage (286 cases, 129 deaths); expression is standardised using the eligible 307-case cohort. No primary adjusted survival association survived BH correction across 14 targets. These results assess disease association, not limonoid efficacy or target engagement.

## Files
analyse.py: statistical analyses from included processed inputs.
prepare.py: original preparation script, included for provenance; requires original source files in data/ and is not needed to reproduce figures from included inputs. Required files: gpl570.gz, gse30784.gz, gse9844.gz, clinical_patient.tsv, clinical_sample.tsv, tcga_expression.tsv, xena_clinical.tsv. Raw downloads are not included.
results/: processed analysis inputs and saved results used for the manuscript figures. These are derived data, not raw source matrices.

Data sources: GEO GSE30784 and GSE9844; cBioPortal hnsc_tcga_pan_can_atlas_2018; UCSC Xena TCGA clinical anatomical annotations. Cite the source studies and resources in the manuscript. No software licence has been selected for this package.
