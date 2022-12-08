# SGI_cancer_recurrence_NIMO
## Requirements
- python (2.7)
- pandas (v 0.24.2)
- matplotlib (v 2.0.0)
- numpy (v 1.16.6)
- scipy (v 1.2.2)

## Installation
All python packages can be installed via pip (https://pypi.org/project/pip/).

## Description
Codes to identify NIMO signature genes by using transcriptome- and methylome-based cell type markers.

## Extract NIMO signatures
- Codes to identify NIMO signature genes is provided under the './code/scripts/transcriptome_methylome_signature_comparison'
- run 'return_network_expansion_of_deconvolution_signature_genes.py'

## Patients recurrence prediction
- Code for reproducing cross-validation of patients recurrence prediction is provided under the './code/scripts/recur_ML_prediction/'
- For SGI cohorts prediction, run 'test_ML_predict_recur_SGI_cohort.py'
- For GSE107442 cohorts prediction, run 'test_ML_predict_recur_GSE107422.py'
