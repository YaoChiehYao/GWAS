#!/usr/bin/env bash
"""
In assignment 10, we used the Plink bioinformatic tool from the
Board Institute. This script runs according to the PLINK tutorial
and is automatic for the GWAS data analysis process, from clean
data to trait association analysis.
"""

### Clean data
# check if hapmap1.ped and hapmap1.map are in the current folder
plink --file hapmap1 --noweb
# Step 1: Make a make a binary PED file
# Filter flags --maf 0.01 --geno 0.01 --mind 0.01 --hew 0.001 --hwe-all
plink --file hapmap1 --make-bed --out hapmap1 --noweb
# Ckeck binary hapmap1.ped file 
# Get three bifile hapmap1.fam, hapmap1.bed and hapmap1.bim for analysis
plink --bfile hapmap1 --noweb
# Working with Binary hapmap1.ped file
# USE bfile to generate allele frequency and codes for each SNP
plink --bfile hapmap1 --freq --out freq_stat --noweb
# Perform a stratified analysis use within
plink --bfile hapmap1 --freq --within pop.phe --out freq_stat --noweb
# more freq_stat.frq.strat
# Check specific SNP frequency within population file pop.phe
plink --bfile hapmap1 --snp rs1891905 --freq --within pop.phe --out snp1_frq_stat --noweb
# Check for Hardy-Weinberg disequilibrium (p-value as threshold)
plink --bfile hapmap1 --hardy --out hapmap1_hardy --noweb
# Step 2: Genotyping rate per individual and per biomarker in the file
plink --bfile hapmap1 --missing --out miss_stat --noweb
# Per-individual genotyping/missing rate 
# more miss_stat.imiss
# Per-marker (locus) genotyping/missing rate
# more miss_stat.lmiss

# Check batch effects for differential genotyping rate (case vs control)
# Plink --bfile hapmap1 --test-missing --out hapmap1_test 
# awk '$5<0.001' hapmap1_test.missing

# Check batch effects for differential genotyping rate (case vs control)
# Plink --bfile hapmap_filtered --test-missing --out hapmap1_test_postfilter 
# awk '$5<0.001' hapmap1_test_postfilter.missing

# Check for excess heterozygosity
# plink --bfile hapmap_filtered --het --out hapmap1_het
# awk '{print $6}' hapmap1_het.het | sort -g | tail
# Drop p-value higher 0.03
# awk '$6>0.03' hapmap1_het.het
# remove bad individuals 
# echo "JPT257 1" > bad.samples
# echo "JPT260 1" >> bad.samples

# remove bad SNPs 
# echo "rs10907185 1" > bad.snps
# echo "rs262688" >> bad.snps

# plink --bfile hapmap1_filtered --remove bad.samples --exclude bad.snps --recode --out hapmap1_clean
# Check for related and duplicate sample 
# plink --bfile hapmap1_clean --genome --out hapmap1_genome
# awk '$10>0.05' hapmap1_genome.genome

### Basic association analysis 
plink --bfile hapmap1 --assoc --out as1 --noweb
# Get an association analysis with sorted result
plink --bfile hapmap1 --assoc --adjust --out as2 --noweb
# Get an adjusted significance values that control for multiple testing
plink --bfile hapmap1 --pheno pop.phe --assoc --adjust --out as3 --noweb
# Run --model flag include basic allelic test, the Cochran-Armitage trend 
# test,dominant and recessive models and a genotypic test
plink --bfile hapmap1 --model --snp rs2222162 --out mod1 --noweb
# Run --cell flag to force the genotypic tests when less than 5 observations
plink --bfile hapmap1 --model --cell 0 --snp rs2222162 --out mod2 --noweb
# Stratification analysis(separate samples to clusters)
# --cluster flag produce IBS clustering, --ppc flag is the size of each cluster
plink --bfile hapmap1 --cluster --mc 2 --ppc 0.05 --out str1 --noweb
# Check the clustering result
# more str1.cluster1



### Association analysis, accounting for clusters
# Use the Cochran-Mantel-Haenszel (CMH) association statistic to test the SNP-diease 
# association on clustering group. 
# --within flag to perform association test conditional on the matching
# --adjust flag is to get a sorted list of CMH result
plink --bfile hapmap1 --mh --within str1.cluster2 --adjust --out aac1 --noweb
# check the adjusted result
# more aac1.cmh.adjusted
# Version 2
plink --bfile hapmap1 --cluster --cc --ppc 0.01 --out version2 --noweb
# Repeat association analysis (optional) if the inflation factor is reduced
plink --bfile hapmap1 --mh --within version2.cluster2 --adjust --out aac2 --noweb
# Version 3  
plink --bfile hapmap1 --cluster --K 2 --out version3 --noweb
# Repeat association analysis (optional) if the inflation factor is reduced
plink --bfile hapmap1 --mh --within pop.phe --adjust --out aac3 --noweb
# Cluster visualization 
plink --bfile hapmap1 --cluster --matrix --out ibd_view --noweb
# Then use R to plot the ibd_view.mibs file


### Quantitative trait association analysis
#  Analyse this quantitative trait directly 
plink --bfile hapmap1 --assoc --pheno qt.phe --out quant1 --noweb
#  A clustered permutation approach to test phenotype with quantitative trait  
plink --bfile hapmap1 --assoc --pheno qt.phe --perm --within str1.cluster2 --out quant2 --noweb
#  An enhanced version use --mperm to specify the test replicates, and speed up the process.
plink --bfile hapmap1 --assoc --pheno qt.phe --mperm 1000 --within str1.cluster2 --out quant3 --noweb
# Test the association of phenotype with the two populations
plink --bfile hapmap1 --pheno qt.phe --gxe --covar pop.phe --snp rs2222162 --out quant3 --noweb



### Extracting a SNP of interested 
# Use recoedAD to extract the snp id rs2222162, and output a separated file called rec_snp1
plink --bfile hapmap1 --snp rs2222162 --recodeAD --out rec_snp1 --noweb
# We can read the rec_snp1 file by read.table function in R and pipe to a variable d
d <- read.table("rec_snp1.recode.raw" , header=T)
# Use a linear model glm function in R statistic package to study their coefficient and correlation   
summary(glm(PHENOTYPE-1 ~ rs2222162_A, data=d, family="binomial"))
