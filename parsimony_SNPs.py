# !/usr/bin/env python
# encoding: utf-8

"""
parsimony_SNPs.py: a program to identify parsimony-informative SNPs in VCF files.

Action: this program filters a VCF file containing biallelic SNP data such that 
all sites will have at least two individuals homozygous for the minor allele.

Purpose: to select sites that will be parsimony informative even after creating 
a consensus DNA sequence (e.g., a Phylip file) that represents heterozygous sites 
(sites with no fixed variants) with IUPAC ambiguity codes. Some programs (e.g., IQ-TREE) 
are not able to integrate information from these heterozygous sites. Thus, especially 
for analyses that depend on unlinked SNPs, where a "representative SNP" is chosen for a 
given sequence interval of DNA, it can be helpful to first remove heterozygous sites 
before SNP thinning. This guarantees that any "representative SNP" will provide 
phylogenetic information. 

Copyright 2021 Andre E. Moncrieff. All rights reserved.

"""


import argparse
import pandas
import numpy
import csv
import itertools
import re


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_file", required=True,
                        help="Enter the file name (including txt extension)",
                        type=str)
    args = parser.parse_args()
    return args


def read_in_csv(txt_file_dataframe):
    raw_dataframe = pandas.read_csv(txt_file_dataframe,
                                    sep='\t',
                                    encoding = "ISO-8859-1",
                                    dtype=str)
    return raw_dataframe


def rows_to_list(vcf_df):
    list_of_rows = vcf_df.values.tolist()
    return list_of_rows


def trim_lists(list_of_rows):
    trimmed_list_of_rows = []
    for item in list_of_rows:
        trimmed_row= item[9:]
        trimmed_list_of_rows.append(trimmed_row)
    return trimmed_list_of_rows


def genotype_list(list_of_rows):
    genotypes_by_row = []
    for row in list_of_rows:
        rowlist = []
        for item in row:
            genotype = item.split(':')[0]
            rowlist.append(genotype)
        genotypes_by_row.append(rowlist)
    return genotypes_by_row


def set_of_genotypes(genotypes_by_row):
    list_of_genotypes = []
    for row in genotypes_by_row:
        list_of_genotypes.extend(row)
    set_of_genotypes = set(list_of_genotypes)
    list_of_genotypes2 = sorted(set_of_genotypes)
    return list_of_genotypes2


def allele_list(genotypes_by_row):
    alleles_by_row = []
    for row in genotypes_by_row:
        row_list_alleles = []
        for genotype in row:
            alleles = re.split('\\||/', genotype)
            row_list_alleles.extend(alleles)
        alleles_by_row.append(row_list_alleles)
    return alleles_by_row


def minor_allele_by_row(list_of_alleles_by_row):
    minor_allele_by_row = []
    #row = list_of_alleles_by_row[-1]
    for row in list_of_alleles_by_row:
        allelezero = row.count('0')
        allele1 = row.count('1')
        #allele2 = row.count('2')
        #allele3 = row.count('3')
        #allele4 = row.count('4')
        if allelezero > allele1:
            minor_allele_by_row.append('1')
        else:
            minor_allele_by_row.append('0')
    return minor_allele_by_row


def count_minor_allele_homos(genotypes_by_row, minor_alyl_by_row):
    number_of_minor_allele_homos = []
    for row,mAllele in zip(genotypes_by_row, minor_alyl_by_row):
            homo_min_allele_unphased = row.count(mAllele+'/'+mAllele)
            homo_min_allele_phased = row.count(mAllele+'|'+mAllele)
            total_count = homo_min_allele_unphased + homo_min_allele_phased
            number_of_minor_allele_homos.append(total_count)
    #return count_of_minor_allele_homos
    #return minor_alyl_by_row[5]
    return number_of_minor_allele_homos
    #return
    #print(len(genotypes_by_row))
    #print(len(minor_alyl_by_row))


def count_alternate_individuals(genotypes_by_row):
    homo_count = []
    for row in genotypes_by_row:
        homo1by1 = row.count('1/1')
        het1by2 = row.count('1/2')
        het1by3 = row.count('1/3')
        het1by4 = row.count('1/4')
        homo1p1 = row.count('1|1')
        homo2by2 = row.count('2/2')
        het2by3 = row.count('2/3')
        het2by4 = row.count('2/4')
        homo2p2 = row.count('2|2')
        homo3by3 = row.count('3/3')
        homo3p3 = row.count('3|3')
        homo4by4 = row.count('4/4')

        homosum = homo1by1 + het1by2 + het1by3 + het1by4 + homo1p1 + homo2by2 \
        + het2by3 + het2by4 + homo2p2 + homo3by3 + homo3p3 + homo4by4

        homo_count.append(homosum)
    return homo_count


def index_counts(homo_count_list):
    desirables_index_list = []
    for item in homo_count_list:
        if item >= 2:
            desirables_index_list.append("yes")
        else:
            desirables_index_list.append("no")
    return desirables_index_list


def add_column_pandas(vcf_df, desirables_list):
    vcf_df['Indicator'] = desirables_list
    return vcf_df


def filter_dataframe(new_vcf_df):
    filtered_dataf = new_vcf_df[new_vcf_df['Indicator'].isin(['yes'])]
    return filtered_dataf


def delete_column(filtered_dataf):
    del filtered_dataf['Indicator']
    return filtered_dataf


def main():
    #create args object
    args = parser()
    #read in dataframes
    vcf_df = read_in_csv(args.vcf_file)
    #print(vcf_df)
    list_of_rows = rows_to_list(vcf_df)
    #print(list_of_rows)
    trimmed_list_of_rows = trim_lists(list_of_rows)
    #print(trimmed_list_of_rows)
    list_of_genotypes = genotype_list(trimmed_list_of_rows)
    #print(list_of_genotypes)
    unique_genotypes = set_of_genotypes(list_of_genotypes)
    #print(unique_genotypes)
    list_of_alleles_by_row = allele_list(list_of_genotypes)
    min_allele_by_row = minor_allele_by_row(list_of_alleles_by_row)
    #print(list_of_alleles_by_row)
    min_allele_homos = count_minor_allele_homos(list_of_genotypes, min_allele_by_row)
    list_of_alt_counts = count_alternate_individuals(list_of_genotypes)
    #print(list_of_alt_counts)
    #print(len(list_of_homo_counts))
    #print("\n")
    #desirables_list = index_counts(list_of_alt_counts)
    desirables_list = index_counts(min_allele_homos)
    #print(desirables_list)
    new_vcf_df = add_column_pandas(vcf_df, desirables_list)
    #print(new_vcf_df)
    filtered_datafr = filter_dataframe(new_vcf_df)
    #print(filtered_datafr)
    final_SNP_df = delete_column(filtered_datafr)
    #print(final_SNP_df)
    #print(desirables_list)
    #print(len(desirables_list))
    final_SNP_df.to_csv('parsimony_SNPs.csv', sep='\t', index=False)


if __name__ == '__main__':
    main()
