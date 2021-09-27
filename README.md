# parsimony_SNPs.py: a program to identify parsimony-informative SNPs in VCF files.

##### Written by Andre E. Moncrieff, 2021.

## Introduction 

**Action**: this program filters a VCF file containing biallelic SNP data such that all sites will have at least two individuals homozygous for the minor allele.

**Purpose**: to select sites that will be parsimony informative even after creating a consensus DNA sequence (e.g., a Phylip file) that represents heterozygous sites (sites with no fixed variants) with IUPAC ambiguity codes. Some programs (e.g., IQ-TREE) are not able to integrate information from these heterozygous sites. Thus, especially for analyses that depend on unlinked SNPs, where a "representative SNP" is chosen for a given sequence interval of DNA, it can be helpful to first remove
heterozygous sites before SNP thinning. This guarantees that any "representative SNP" will provide phylogenetic information. 

## Step-by-step instructions 
#### (*A few quick manual steps*)

- Remove the header lines in your VCF before the line starting with '#CHROM'
- Change your vcf file extension to .txt
- Run: `python parsimony_SNPs.py --vcf_file yourVCF.txt` 
- Output file is csv: 'parsimony_SNPs.csv' for easy visualization
- Finally, if desired, add full VCF header info and change extension back to .vcf

## Citation

[![DOI](https://zenodo.org/badge/411024089.svg)](https://zenodo.org/badge/latestdoi/411024089)

**Moncrieff, A.E.** 2021. Parsimony SNPs v1.0: a program to identify parsimony-informative SNPs in VCF files. DOI:10.5281/zenodo.5532931