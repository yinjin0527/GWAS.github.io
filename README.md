# GWAS step by step

## Install PLINK using conda 
```
conda install -c bioconda plink
conda install -c bioconda/label/cf201901 plink
```

## Using a example dataset from plink example
```
mkdir -p ~/GWAS
cd GWAS
wget http://zzz.bwh.harvard.edu/plink/dist/example.zip
unzip example.zip
cd example/

```
## Check Ped file maping information and change to binary (fam bed bai) files
```
$ plink --file wgas1
$ plink --file wgas1 --make-bed --out wgas1
# --mind 0.05" if you want to sent quality more than 95% snp
```
![filemap](https://www.researchgate.net/publication/281588338/figure/fig6/AS:281417826226190@1444106646899/Genome-wide-association-data-files-GWA-data-files-are-typically-organized-into_W640.jpg) 







## PLink quanlity control,command and function.
| Step                                   | Command                                                      | Function                                                     |
| -------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| 1: Missingness of SNPs and individuals | ‐‐geno                                                       | Excludes SNPs that are missing in a large proportion of the subjects. In this step, SNPs with low genotype calls are removed. |
| ‐‐mind                                 | Excludes individuals who have high rates of genotype missingness. In this step, individual with low genotype calls are removed. |                                                              |
| 2: Sex discrepancy                     | ‐‐check‐sex                                                  | Checks for discrepancies between sex of the individuals recorded in the dataset and their sex based on X chromosome heterozygosity/homozygosity rates. |
| 3: Minor allele frequency (MAF)        | ‐‐maf                                                        | Includes only SNPs above the set MAF threshold.              |
| 4: Hardy–Weinberg equilibrium (HWE)    | ‐‐hwe                                                        | Excludes markers which deviate from Hardy–Weinberg equilibrium. |
| 5: Heterozygosity                      | For an example script see https://github.com/MareesAT/GWA_tutorial/ | Excludes individuals with high or low heterozygosity rates   |
| 6: Relatedness                         | ‐‐genome                                                     | Calculates identity by descent (IBD) of all sample pairs.    |
| ‐‐min                                  | Sets threshold and creates a list of individuals with relatedness above the chosen threshold. Meaning that subjects who are related at, for example, pi‐hat >0.2 (i.e., second degree relatives) can be detected. |                                                              |
| 7: Population stratification           | ‐‐genome                                                     | Calculates identity by descent (IBD) of all sample pairs.    |
| ‐‐cluster ‐‐mds‐plot *k*               | Produces a *k*‐dimensional representation of any substructure in the data, based on IBS. |                                                              |
