# DDR

## Requirements:
* R (packages: edgeR)
* Python (packages: pandas)

## Running step wise

### Step0 Preprocess
```
Rscript step0_preprocess.R
```
### Step1 Calculate Stats
```
Rscript step1_calculateStats.R 
```
### Step2 Find Reference set
```
Rscript step2_findRef.R
```
### Step3 Find overlap using Fisher's Exact test
Since the example dataset has 115 samples in group1.
```
Rscript step3_overlapFisher.R 115
```


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
