Python scripts and processed data for expression and survival analysis of 14 predicted targets in oral cancer.
**Quick start**
Run these commands from the repository folder:
pip install -r requirements.txt
python analyse.py
python plot_figure9_expression.py
python plot_figure10_survival.py

**Datasets**
- GSE30784: 167 tumours and 45 normal samples; dysplasia excluded.
- GSE9844: 26 tumours and 12 normal controls.
- TCGA-HNSC: Oral-site cases; adjusted survival analysis includes 286 patients and 129 deaths.
Expression comparisons use Welch tests. Survival models adjust for age, sex and stage. Both analyses apply Benjamini–Hochberg correction across 14 targets.
**Files and outputs**
- analyse.py: Runs statistical analyses and generates results/integrated_target_evidence.csv for Table 5.
- plot_figure9_expression.py: Creates expression plots in figure9_expression/.
- plot_figure10_survival.py: Creates the survival forest plot in figure10_survival/.
- results/: Contains processed inputs and saved results.
- prepare.py: Prepares data from original source files, which are not included. This step is unnecessary when using the supplied processed inputs.
Figures - exported as 300-dpi PNG and vector PDF files.
**Findings**
EGFR, MMP9 and MMP7 showed replicated tumour upregulation. No target retained a significant adjusted survival association after multiple-testing correction. These findings support target prioritisation for further investigation; they do not establish limonoid anticancer activity.
