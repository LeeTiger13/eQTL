# !/bin/bash
# by zhiyuan

workdir=/work/
datadir=$workdir/data
scriptdir=$workdir/scripts

#####################################################################
# Preparation of genotype files for eQTL analysis
#####################################################################
cd $workdir
mkdir 01.genotype_files
cd 01.genotype_files

python $scriptdir/1_Preparation_of_Genotype.py $datadir/genotype.vcf eQTL_analysis

#####################################################################
# Processing expression data
#####################################################################
cd $workdir
mkdir 02.expression_data
cd 02.expression_data

Rscript $scriptdir/2_Preprocessing_Gene_Expression_Data.R

#####################################################################
# eQTL analysis
#####################################################################
cd $workdir
mkdir 03.eQTL_analysis
cd 03.eQTL_analysis
mkdir 01.result

Rscript $scriptdir/3_Matrixeqtl.R --genelocfile $datadir/gene_loc.txt --expressionfile $workdir/02.expression_data/gene.expressed.normal.txt --genofille $workdir/01.genotype_files/eQTL_analysis.geno --snpposinfo $workdir/01.genotype_files/eQTL_analysis.info --covariate $workdir/02.expression_data/fin_cov_10HF_5PS.txt --outpath $workdir/03.eQTL_analysis/01.result

#####################################################################
# Trans-eQTL filtering
#####################################################################
mkidr 02.trans_eQTL_filtering
cd 02.trans_eQTL_filtering
mkdir 01.tem
cd 01.tem

python $scriptdir/4.Trans-eQTL_LD_Filtering.py $workdir/01.genotype_files/eQTL_analysis.geno $workdir/03.eQTL_analysis/01.result/trans_qqnorm_result.txt

awk 'FNR==1 && NR!=1{next;}{print}' *trans_ld_final.txt > $workdir/03.eQTL_analysis/02.trans_eQTL_filtering/trans_ld_final.txt

#####################################################################
# SMR
#####################################################################
cd $workdir
mkdir 04.SMR
cd 04.SMR

cat $workdir/03.eQTL_analysis/01.result/cis_qqnorm_result.txt $workdir/03.eQTL_analysis/01.result/trans_qqnorm_result.txt > $workdir/04.SMR/mateQTL.txt

smr --eqtl-summary mateQTL.txt --matrix-eqtl-format --make-besd --out mybesd 
smr --bfile /work/data/gene --gwas-summary mygwas.ma --beqtl-summary myeqtl --out mysmr --thread-num 10 --cis-wind 200 --peqtl-smr 1.5e-8





















