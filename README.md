# DDR

## Requirements:
* R (packages: edgeR)
* Python (packages: pandas)

## RNASeq Pipeline: 
Following series of steps show how to run the DDR method on RNASeq data. 

### Step0 Preprocess
Step0 should be used to format the count table such that the input table has the gene expression table with first n columns as group1 samples (Triple negative:TN) and remaining columns as samples from group2 (other:OT). Rows represent the ENSG ids.
```
Rscript step0_preprocess.R
```
### Step1 Calculate Stats
In this step, the count table is normalized and the covariance, standard deviation, mean and MFC are being calculated.
```
Rscript step1_calculateStats.R RNASeq

```
### Step2 Find Reference set
In this step, reference set of genes are being determined. The output file 'ref_cpm.csv' stores expression level of these reference genes.
```
Rscript step2_findRef.R RNASeq

```
### Step3 Find overlap using Fisher's Exact test
#### Enter number of samples in first group as first input argument.
Since the example dataset has 115 samples in group1.
```
Rscript step3_overlapFisher.R 115 RNASeq
```

## Output files
* final_out.csv: Normalized count data
* ref_cpm.csv: Expression of reference genes
* overlap_test_fdr_1_[RNASeq|microarray].csv or : Differentially expressed genes with fdr < 0.1
* overlap_test_fdr_05_[RNASeq|microarray].csv: Differentially expressed genes with fdr < 0.05

## Pipeline for microarray data. Note that the steps are similar to for RNASeq data. 
Rscript step0_preprocess_microarray.R
Rscript step1_calculateStats.R microarray 
Rscript step2_findRef.R microarray
Rscript step3_overlapFisher.R 115 microarray

## Running all the steps together 
```
Rscript DDR_Ref.R
```
If you have another python version that has pandas:
```
Rscript DDR_Ref.R PYTHON_PATH
```
