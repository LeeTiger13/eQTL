# -*- UTF-8 -*-
# !/usr/bin/env python3
# by zhiyuan

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import sys
import networkx as nx
from community import community_louvain

genotype_file = sys.argv[1]
trans_file = sys.argv[2]

genotype = pd.read_csv(genotype_file,header=0,index_col=0,sep='\t')
ck_trans = pd.read_csv(trans_file,header=0,sep='\t')
ck_trans = ck_trans.loc[ck_trans['p-value']<0.05/3333393,:]

for gene in ck_trans['gene'].unique():
    gene_ck_df = ck_trans[ck_trans['gene']==gene]
    gene_ck_df['chr'] = gene_ck_df['SNP'].apply(lambda x: int(x.split('_')[0][2:]))
    gene_ck_df['pos'] = gene_ck_df['SNP'].apply(lambda x: int(x.split('_')[1]))
    gene_ck_df.sort_values(by=['chr','pos'],inplace=True)
    gene_ck_df.reset_index(drop=True,inplace=True)
    gene_ck_df['bin'] = None
    bin_num = 0
    prev_chr = None
    prev_ps = 0
    for index, row in gene_ck_df.iterrows():
        if prev_chr is None or row['chr'] != prev_chr or row['pos'] - prev_ps > 10000:
            bin_num += 1
        bin_label = f'bin{bin_num}'
        prev_chr = row['chr']
        prev_ps = row['pos']
        gene_ck_df.at[index, 'bin'] = bin_label
    gene_ck_df = gene_ck_df.groupby('bin').filter(lambda x: len(x)>=3)
    gene_ck_df.to_csv(f'{gene}_ms.txt',sep='\t',index=False)
    try:
        gene_ck_df = gene_ck_df.groupby('bin').apply(lambda x: x.loc[x['p-value'].idxmin()])
        gene_ck_df.sort_values(by=['chr','pos'],inplace=True)
        genotype_gene_ck = genotype.loc[gene_ck_df['SNP'],:]
        if genotype_gene_ck.shape[0] > 1:
            correlation_matrix = np.corrcoef(genotype_gene_ck.T, rowvar=False)**2
            np.fill_diagonal(correlation_matrix, 0)
            correlation_matrix = pd.DataFrame(correlation_matrix, index=genotype_gene_ck.index, columns=genotype_gene_ck.index)
            upper_triangle_matrix = correlation_matrix.where(correlation_matrix > np.triu(correlation_matrix, k=1).astype(bool))
            adjacency_matrix = upper_triangle_matrix.values
            threshold = 0.1
            adjacency_matrix = (adjacency_matrix > threshold).astype(int)
            G = nx.from_numpy_array(adjacency_matrix)
            partition = community_louvain.best_partition(G)
            cluster_df = pd.DataFrame({'SNP': genotype_gene_ck.index, 'cluster': partition.values()})
            cluster_df['cluster'] = cluster_df['cluster'].apply(lambda x: f'{gene}_{x+1}')
        else:
            cluster_df = pd.DataFrame({'SNP': gene_ck_df['SNP'], 'cluster': f'{gene}_{1}'})
        gene_result = pd.merge(gene_ck_df, cluster_df, on='SNP', how='left')
        gene_result.to_csv(f'{gene}_trans_ld.txt',sep='\t',index=False)
        gene_fin = gene_result.groupby('cluster').apply(lambda x: x.loc[x['p-value'].idxmin()])
        gene_fin.to_csv(f'{gene}_trans_ld_final.txt',sep='\t',index=False)
    except:
        print(f"{gene} has no significant SNPs!")


