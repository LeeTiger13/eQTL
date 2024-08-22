# -*- UTF-8 -*-
# !/usr/bin/env python3
# by zhiyuan

import numpy as np
import pandas as pd
import sys, os

vcf_file = sys.argv[1]
prefix = sys.argv[2]

f1 = open(vcf_file, 'r')
for l in f1:
    if l.startswith('#C'):
        samples = l.strip().split('\t')[9:]
        break
f1.close()

print("reading vcf file...")
vcf = pd.read_table(vcf_file, sep='\t', header=None, comment='#')
vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples

print("writing info file...")
# Prepare to make info file
vcf_info = vcf[['ID','CHROM', 'POS']]
vcf_info.columns = ['SNP', 'Chromosome', 'Position']
vcf_info.to_csv(prefix + '.info', sep='\t', index=False, header=True)

print("writing geno file...")
# Prepare to make geno file
vcf_geno = vcf[['ID'] + samples]
vcf_geno = vcf_geno.replace('0/0', '0').replace('0/1', '1').replace('1/1', '2').replace('./.', '1')
vcf_geno = vcf_geno.set_index('ID')
vcf_geno.index.name = 'Taxa'
vcf_geno.to_csv(prefix + '.geno', sep='\t', index=True, header=True)