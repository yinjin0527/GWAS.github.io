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


## quality control
Remove Id with following properties:
call rate <0.98
absval inbreeding coefficient >0.2
if trio dataset: mendel errors>10000
sex check errors
if reported sex =male Fhet<0.5
if reported dex=female Fhet>0.5

Remove SNPs with following properties:
call rate <0.98
call rate difference between case/control >0.2
if trio dataset: mendel error>4
HWE P-value:
if trio P<1e-6
if case only, P<1e-6
if control only, P<1e-6
if both cases and control, P calculated from cases <1e-10 or from control <1e-6
invariant

| # per ID   | mind=0.02              | include only  IDs with missing-rate < NUM                    |
| ---------- | ---------------------- | ------------------------------------------------------------ |
| #  per SNP | geno=0.02              | include only SNPs with missing-rate < NUM                    |
| #  per SNP | maf=0.00               | include only SNPs with MAF > NUM                             |
| #  per SNP | midi=0.02              | include only SNPs with case/control  missing-rate-difference < NUM |
| #  per SNP | pre_geno=0.05          | include only SNPs with missing-rate < NUM  (before ID filter) |
| #  per SNP | lmend=4                | number of mendelian errors per ID                            |
| #  per ID  | imend =10000           | max number of mendelian errors per SNP                       |
| #  per ID  | Fhet_th=.2             | NUM -NUM < FHET < NUM                                        |
|            | hwe_th_co=0.000001     | NUM  HWE_controls < NUM                                      |
|            | hwe_th_ca=0.0000000001 | NUM HWE_cases < NUM                                          |
|            | withpna=0              | inlcude SNPs with p=NA (monomorph)                           |
|            | sexmin = 10            | INT minimum number of chrX SNPs to perform  sexcheck, default: $sexcheck_min |

```
#first step
#create two types of missing files id missing .imiss file, snp .lmiss file
$ plink --bfile wgas1 --missing --out miss_stat
#check for id missing

#missing individuals (N MISS) and the proportion of individuals missing (F MISS)
head miss_stat.imiss

#Filtered out SNPs with genotyping efficiency below 95%
plink --noweb --bfile wgas1 --geno 0.05 --make-bed --out geno_wgas1



#MAF test
plink --bfile wgas1 --freq --out freq_stat
head freq_stat.frq

plink --bfile test --geno 005 --hwe 0.000001 --mad 0.5 --mind 0.1 --make-bed --out test_qc1

