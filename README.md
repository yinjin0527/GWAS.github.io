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


## Quality control
professor Paschou's pipeline:
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

|            | functioln              | explanation                                                  |
| ---------- | ---------------------- | ------------------------------------------------------------ |
| #  per SNP | geno                   | include only SNPs with missing-rate < NUM                    |
| #  per SNP | maf                    | include only SNPs with MAF > NUM                             |
| #  per SNP | test-missing           | include only SNPs with case/control  missing-rate-difference < NUM |
| #  per SNP | pre_geno               | include only SNPs with missing-rate < NUM  (before ID filter) |
| #  per SNP | lmend                  | number of mendelian errors per ID                            |
| #  per ID  | mind                   | include only  IDs with missing-rate < NUM                    |
| #  per ID  | imend                  | max number of mendelian errors per SNP                       |
| #  per ID  | Fhet_th                | NUM -NUM < FHET < NUM                                        |
|            | hwe_th_co              | NUM  HWE_controls < NUM                                      |
|            | hwe_th_ca              | NUM HWE_cases < NUM                                          |
|            | withpna                | inlcude SNPs with p=NA (monomorph)                           |
|            | check-sex              | INT minimum number of chrX SNPs to perform  sexcheck, default: $sexcheck_min |


### Before quality control
```
#1.check missing value
#create two types of missing files id missing .imiss file, snp .lmiss file
$ plink --bfile wgas1 --missing --out miss_stat
#check for id missing

#2.check sample gender
plink --bfile wgas1  --check-sex --out gender

#3.missing individuals (N MISS) and the proportion of individuals missing (F MISS)
head miss_stat.imiss
head miss_stat.lmiss
```
 | FID     | IID MISS_PHENO | N_MISS   N | _GENO | F_MISS |          |
| ------- | -------------- | ---------- | ----- | ------ | :------: |
| CH18526 | NA18526        | N          | 1057  | 228694 | 0.004622 |
| CH18524 | NA18524        | N          | 967   | 228694 | 0.004228 |
| CH18529 | NA18529        | N          | 638   | 228694 | 0.00279  |
| CH18558 | NA18558        | N          | 2342  | 228694 | 0.01024  |
| CH18532 | NA18532        | N          | 1408  | 228694 | 0.006157 |
| CH18561 | NA18561        | N          | 479   | 228694 | 0.002095 |
| CH18562 | NA18562        | N          | 1263  | 228694 | 0.005523 |
| CH18537 | NA18537        | N          | 573   | 228694 | 0.002506 |
| CH18603 | NA18603        | N          | 453   | 228694 | 0.001981 |
 
 
| CH   | R      SNP | N_MISS | N_GE | NO   F_MISS |
| ---- | ---------- | ------ | ---- | ----------- |
| 1    | rs3094315  | 1      | 90   | 0.01111     |
| 1    | rs6672353  | 1      | 90   | 0.01111     |
| 1    | rs4040617  | 0      | 90   | 0           |
| 1    | rs2905036  | 0      | 90   | 0           |
| 1    | rs4245756  | 0      | 90   | 0           |
| 1    | rs4075116  | 0      | 90   | 0           |
| 1    | rs9442385  | 1      | 90   | 0.01111     |
| 1    | rs6603781  | 1      | 90   | 0.01111     |
| 1    | rs11260562 | 2      | 90   | 0.02222     |
   
   
   ```
#4.check MAF
plink --bfile wgas1 --freq --out freq_stat
head freq_stat.frq
```
| CHR  | SNP        | A1   | A2   | MAF      | NCHROBS |
| ---- | ---------- | ---- | ---- | -------- | ------- |
| 1    | rs3094315  | G    | A    | 0.1236   | 178     |
| 1    | rs6672353  | A    | G    | 0.005618 | 178     |
| 1    | rs4040617  | G    | A    | 0.1167   | 180     |
| 1    | rs2905036  | 0    | T    | 0        | 180     |
| 1    | rs4245756  | 0    | C    | 0        | 180     |
| 1    | rs4075116  | C    | T    | 0.05556  | 180     |
| 1    | rs9442385  | T    | G    | 0.3933   | 178     |
| 1    | rs6603781  | 0    | G    | 0        | 178     |
| 1    | rs11260562 | A    | G    | 0.02841  | 176     |
   
   ```

#5.check HWE
plink --bfile wgas1 --hardy out HWE
head HWE.hwe
```
| CHR  | SNP       | TEST  | A1   | A2   | GENO    | O(HET)  | E(HET)  | P      |
| ---- | --------- | ----- | ---- | ---- | ------- | ------- | ------- | ------ |
| 1    | rs3094315 | ALL   | G    | A    | 0/22/67 | 0.2472  | 0.2166  | 0.3476 |
| 1    | rs3094315 | AFF   | G    | A    | 0/15/33 | 0.3125  | 0.2637  | 0.5771 |
| 1    | rs3094315 | UNAFF | G    | A    | 0/7/34  | 0.1707  | 0.1562  | 1      |
| 1    | rs6672353 | ALL   | A    | G    | 0/1/88  | 0.01124 | 0.01117 | 1      |
| 1    | rs6672353 | AFF   | A    | G    | 0/1/48  | 0.02041 | 0.0202  | 1      |
| 1    | rs6672353 | UNAFF | A    | G    | 0/0/40  | 0       | 0       | 1      |
| 1    | rs4040617 | ALL   | G    | A    | 0/21/69 | 0.2333  | 0.2061  | 0.5994 |
| 1    | rs4040617 | AFF   | G    | A    | 0/14/35 | 0.2857  | 0.2449  | 0.5714 |
| 1    | rs4040617 | UNAFF | G    | A    | 0/7/34  | 0.1707  | 0.1562  | 1      |


```
#6.Heterozygosity
plink --bfile wgas1  --het --out inbreed 
head inbreed.het

```
| FID     | IID     | O(HOM) | E(HOM)   | N(NM)  | F         |
| ------- | ------- | ------ | -------- | ------ | --------- |
| CH18526 | NA18526 | 130790 | 1.31E+05 | 187689 | 0.001195  |
| CH18524 | NA18524 | 131114 | 1.31E+05 | 187759 | 0.005865  |
| CH18529 | NA18529 | 130736 | 1.31E+05 | 188040 | -0.00402  |
| CH18558 | NA18558 | 129964 | 1.30E+05 | 186635 | 0.000776  |
| CH18532 | NA18532 | 129957 | 1.31E+05 | 187376 | -0.00894  |
| CH18561 | NA18561 | 130857 | 1.31E+05 | 188173 | -0.003463 |
| CH18562 | NA18562 | 129830 | 1.31E+05 | 187492 | -0.01329  |
| CH18537 | NA18537 | 130543 | 1.31E+05 | 188095 | -0.008016 |
| CH18603 | NA18603 | 130580 | 1.31E+05 | 188191 | -0.008352 |

```
#7.Relatedness
plink --bfile wgas1  --genome  --out pihat
head pihat.genome
```
| FID1    | IID1    | FID2    | IID2 RT |      | EZ   | Z0     | Z1    | Z2 P   | I_HAT PHE |      | DST      | PPC    | RATIO  |
| ------- | ------- | ------- | ------- | ---- | ---- | ------ | ----- | ------ | --------- | ---- | -------- | ------ | ------ |
| CH18526 | NA18526 | CH18524 | NA18524 | UN   | NA   | 0.8139 | 0     | 0.1861 | 0.1861    | -1   | 0.797465 | 0.4549 | 1.9923 |
| CH18526 | NA18526 | CH18529 | NA18529 | UN   | NA   | 0.8093 | 0     | 0.1907 | 0.1907    | -1   | 0.799983 | 0.1855 | 1.9404 |
| CH18526 | NA18526 | CH18558 | NA18558 | UN   | NA   | 0.8204 | 0     | 0.1796 | 0.1796    | -1   | 0.797216 | 0.1506 | 1.9311 |
| CH18526 | NA18526 | CH18532 | NA18532 | UN   | NA   | 0.8158 | 0.015 | 0.1692 | 0.1767    | -1   | 0.796774 | 0.2359 | 1.9515 |
| CH18526 | NA18526 | CH18561 | NA18561 | UN   | NA   | 0.8164 | 0     | 0.1836 | 0.1836    | -1   | 0.798742 | 0.616  | 2.0202 |
| CH18526 | NA18526 | CH18562 | NA18562 | UN   | NA   | 0.8233 | 0     | 0.1767 | 0.1767    | -1   | 0.796276 | 0.0171 | 1.862  |
| CH18526 | NA18526 | CH18537 | NA18537 | UN   | NA   | 0.8193 | 0     | 0.1807 | 0.1807    | 0    | 0.797589 | 0.1876 | 1.9407 |
| CH18526 | NA18526 | CH18603 | NA18603 | UN   | NA   | 0.8202 | 0     | 0.1798 | 0.1798    | 0    | 0.797009 | 0.2407 | 1.9528 |
| CH18526 | NA18526 | CH18540 | NA18540 | UN   | NA   | 0.8097 | 0     | 0.1903 | 0.1903    | -1   | 0.800521 | 0.105  | 1.9167 |
  
  

### perform quality control
```
#1)Delete SNPs and individuals with high levels of missingness

plink --bfile wgas1 --geno 0.02 --make-bed --out wgas2
plink --bfile wgas2 --mind 0.02 --make-bed --out wgas3
```
23441 variants removed due to missing genotype data (--geno).
205253 variants and 90 people pass filters and QC.

205253 variants and 89 people pass filters and QC.
Among remaining phenotypes, 48 are cases and 41 are controls.

```
# 2) Check for sex discrepancy
# 1. Delete individuals with sex discrepancy.
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
plink --bfile wgas3 --remove sex_discrepancy.txt --make-bed --out wgas4 


# 2. impute-sex.
#plink --bfile wgas4 --impute-sex --make-bed --out wgas5
```
No variants removed due to sex discrepancy for this small dataset

```
3) MAF

# Generate a bfile with autosomal SNPs only and delete SNPs with a low minor allele frequency (MAF).

# Select autosomal SNPs only.
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' wgas5.bim > snp_1_22.txt
plink --bfile wgas5 --extract snp_1_22.txt --make-bed --out wgas6

# Remove SNPs with a low MAF frequency.
plink --bfile wgas6 --maf 0.05 --make-bed --out wgas7
```
60890 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
144363 variants and 89 people pass filters and QC.


```
4) HWE
# Therefore, we use two steps, first we use a stringent HWE threshold for controls, followed by a less stringent threshold for the case data.
plink --bfile wgas7 --hwe 1e-6 --make-bed --out wgas8

# This second HWE step only focusses on cases because in the controls all SNPs with a HWE p-value < hwe 1e-6 were already removed
plink --bfile wgas8 --hwe 1e-10 --hwe-all --make-bed --out wgas9
```
Total genotyping rate is 0.998381.
--hwe: 0 variants removed due to Hardy-Weinberg exact test.
144363 variants and 89 people pass filters and QC.``

```
5) heterozygosity rate 
# remove individuals with a heterozygosity rate deviating more than 3 sd from the mean.
# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high inversion regions  and prune the SNPs using the command --indep-pairwise’.
# The parameters ‘50 5 0.2’ stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.

plink --file wgas9--make-set inversion-ld.txt --write-set --out inversion
plink --bfile wgas9 --exclude inversion-ld.txt --range --indep-pairwise 50 5 0.2 --out indepSNP



# Adapt this file to make it compatible for PLINK, by removing all quotation marks from the file and selecting only the first two columns.
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt

# Remove heterozygosity rate outliers.
plink --bfile wgas9 --remove het_fail_ind.txt --make-bed --out wgas10
```
```
6) heterozygosity rate 
# Assuming a random population sample we are going to exclude all individuals above the pihat threshold of 0.2 in this tutorial.
# Check for relationships between individuals with a pihat > 0.2.

plink --bfile wgas10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2

# The HapMap dataset is known to contain parent-offspring relations. 
# The following commands will visualize specifically these parent-offspring relations, using the z values. 
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
```

```
7) relatedness
# The generated plots show a considerable amount of related individuals (explentation plot; PO = parent-offspring, UN = unrelated individuals) in the Hapmap data, this is expected since the dataset was constructed as such.
# Normally, family based data should be analyzed using specific family based methods. relatedness here  as cryptic relatedness in a random population sample.

plink --bfile wgas10 --filter-founders --make-bed --out wgas11

#  individuals with a pihat >0.2.
plink --bfile wgas11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders

# For each pair of 'related' individuals with a pihat > 0.2, we recommend to remove the individual with the lowest call rate. 
plink --bfile wgas11 --missing

# Delete the individuals with the lowest call rate in 'related' pairs with a pihat > 0.2 
plink --bfile wgas11 --remove 0.2_low_call_rate_pihat.txt --make-bed --out wgas12
```
```
###
8)if trio then Mendel error rate need be considered
plink --bfile wgas12 --me 1 1 --set-me-missing --make-bed --out wgas13
###
```
This dataset is not trio


