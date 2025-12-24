# JOBS
JOint model viewing Bk-eQTLs as a weighted sum of Sc-eQTLs [JOBS](https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00076-X)


## Table of contents
* [Introduction](#Introduction)
* [Installation](#Installation)
* [Quick tutorial](#Quick_tutorial)
* [Input files](#Input_files)
* [Contact](#Contact)

## Introduction
Here, we propose a JOint model viewing Bk-eQTLs as a weighted sum of Sc-eQTLs (JOBS) from constituent cell types. JOBS borrows strength from large sample sizes of bk-eQTLs to improve sc-eQTLs analysis, with improvements bigger for more common cell types. With more accurate sc-eQTLs effect estimates, all integrative analyses can be improved. 

JObS summary statistics (analyzing eQTLGen and OneK1K) are available for download from the [Zenodo repo](https://zenodo.org/records/15097677). It is developed and maintained by [Lida Wang](https://lidawangpsu.github.io/Lida/) at [Dajiang Liu's Group](https://dajiangliu.blog).

## Citation
If you use this dataset or JOBS package, please cite the following publication:

Wang et al., An atlas of single-cell eQTLs dissects autoimmune disease genes and identifies novel drug classes for treatment, Cell Genomics (2025), https://doi.org/10.1016/j.xgen.2025.100820

## Installation
JOBS is hosted on GitHub, making it easy to install and update. Before installing. The following R packages are required: MASS, dplyr, data.table, tidyr, R.utils, pracma, and quadprog.

Steps to Install:
```r
install.packages("devtools")
library(devtools)
```

Install JOBS from the GitHub repository:

```
devtools::install_github("LidaWangPSU/JOBS/JOBS")
library(JOBS)
```


## Quick tutorial
### 1. Prepare bulk and single cell eQTLs summary statistics (effect size and s.e.). 

#### Effect Size Matrix
This matrix contains eQTL effect sizes for both bulk and single-cell data. We recommend using standardized effect sizes and standard errors, for example:

$\beta=z score/\sqrt{N}, s.e.=1/\sqrt{N}$

- **Columns**:
  - **Column 1**: Gene-SNP pair identifiers.
  - **Column 2**: Bulk effect size.
  - **Columns 3+**: Cell type-specific eQTL effect sizes.
 
**Example**:

|    Gene-snp pair    |      Bulk     | Cell type 1  |  Cell type 2 | ...... |  Cell type k |
| ------------------- |      ----     | -----------  |  ----------- | ------ |  ----------- |
| ENSG00000XXXXX-snp1 |      0.65     |      0.55    |      0.30    | ...... |      0.75    |  
| ENSG00000XXXXX-snp2 |     -0.50     |     -0.62    |     -0.35    | ...... |      0.02    |


#### Standard Error Matrix (S.E.)
This matrix has the same dimensions as the effect size matrix and represents the standard errors associated with the eQTL effect sizes.

- **Columns**:
  - **Column 1**: Gene-SNP pair identifiers (should match the effect size matrix).
  - **Column 2**: Bulk effect size standard error.
  - **Columns 3+**: Standard errors for the cell type-specific eQTL effect sizes.

**Example**:

|    Gene-snp pair    |      Bulk     | Cell type 1  |  Cell type 2 | ...... |  Cell type k |
| ------------------- |      ----     | -----------  |  ----------- | ------ |  ----------- |
| ENSG00000XXXXX-snp1 |    0.0071     |    0.0316    |    0.0316    | ...... |    0.0316    |  
| ENSG00000XXXXX-snp2 |    0.0071     |    0.0316    |    0.0316    | ...... |    0.0316    |

- It is acceptable to have missing data in the input. We will analyze all gene-SNP pairs unless sc-eQTL data is entirely missing for all cell types. In other cases, we will simply ignore the missing values and keep "NA" in the output.
   
We have a small example input data [here](https://github.com/LidaWangPSU/JOBS/blob/main/example_data/). Data were subsetted from brain bulk and single cell eQTLs as an example to run the script.

```r
library(data.table)
beta <- as.data.frame(fread("~/example_beta_chr22.txt.gz"))
se <- as.data.frame(fread("~/example_beta_se_chr22.txt.gz"))
```
  
### 2. Run JOBS
#### JOBS Step 1: cell type weights
Here, we used Non-negative least squares to estimate cell type weights.
You can also specify the weights by yourself, e.g. estimated from scRNAseq data or other methods.
```r
weight <- jobs.nnls.weights(beta,se)
```


#### JOBS Step 2: refine eQTLs estimation 
- **weight**: K numeric numbers for K cell types, add up to 1.
- **COR**: Whether consider correlation across cell types, default is FALSE.
```r
jobs_eqtls <- jobs.eqtls(beta,se,weight,COR=F)
```

### 3. Output results
The output includes:
* Refined effect size: in the same format as effect size input
* Refined S.E: in the same format as s.e. input

We uploaded pre-computed gene expression prediction weights [here](https://github.com/LidaWangPSU/JOBS/blob/main/example_data/)

## Contact
Lida Wang [lida.wang.96@gmail.com](lida.wang.96@gmail.com)
