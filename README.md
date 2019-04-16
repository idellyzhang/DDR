# DDR

## Requirements:
R (packages: edgeR)
Python (packages: pandas)

## Running the tool for finding differentially expressed genes 
```
Rscript DDR_Ref.R
```
If you have another python version that has pandas:
```
Rscript DDR_Ref.R PYTHON_PATH
```

## Output files
* final_out.csv: Normalized count data
* ref_cpm_BRCA.csv: Expression of reference genes
* overlap_test_fdr_1_BRCA_DEC.csv: Differentially expressed genes with fdr < 0.1
* overlap_test_fdr_05_BRCA.csv: Differentially expressed genes with fdr < 0.05
