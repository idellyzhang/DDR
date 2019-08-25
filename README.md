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
#### Enter RNASeq or microarray as input argument.
```
Rscript step1_calculateStats.R RNASeq
Rscript step1_calculateStats.R microarray 

```
### Step2 Find Reference set
#### Enter RNASeq or microarray as input argument.

```
Rscript step2_findRef.R RNASeq
Rscript step2_findRef.R microarray

```
### Step3 Find overlap using Fisher's Exact test
#### Enter number of samples in first group as first input argument.
#### Enter RNASeq or microarray as second input argument.

Since the example dataset has 115 samples in group1.
```
Rscript step3_overlapFisher.R 115 RNASeq
Rscript step3_overlapFisher.R 115 microarray

```


## Running all the steps together 
```
Rscript DDR_Ref.R
```
If you have another python version that has pandas:
```
Rscript DDR_Ref.R PYTHON_PATH
```

## Output files
* final_out.csv: Normalized count data
* ref_cpm.csv: Expression of reference genes
* overlap_test_fdr_1_RNASeq.csv".csv: Differentially expressed genes with fdr < 0.1
* overlap_test_fdr_05_RNASeq.csv: Differentially expressed genes with fdr < 0.05
