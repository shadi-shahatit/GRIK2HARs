#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024

#### Content:
#### 1. GWAS date retrieval  
#### 2. LDSC regression
#### Appendix I - Tools and scripts

#### 1. GWAS date retrieval

## download GWAS summery stats from https://doi.org/10.1093/molbev/msae001 

## ASD
## (n = 18,381 cases; 27,969 controls) (Grove et al. 2019) https://figshare.com/articles/dataset/asd2019/14671989
awk -F'\t' '$1 = 6' ASD_Grove2019.txt > ASD_Grove2019_chr6.txt

## ADHD
## (n = 38,691 cases; 186,843 controls) (Demontis et al. 2023) https://figshare.com/articles/dataset/adhd2022/22564390

## SCZ
## (n = 76,755 cases; 243,649 controls) (Trubetskoy et al. 2022) https://figshare.com/articles/dataset/scz2022/19426775
grep -v '#' SCZ_trubetskoy2022.tsv | awk -F'\t' 'NR > 1 && $1 = 6' > SCZ_trubetskoy2022_chr6.tsv
awk -F'\t' '$11 < 0.05' SCZ_trubetskoy2022_chr6.tsv > SCZ_trubetskoy2022_chr6_p0.05.tsv

## BD
## (n = 41,917 cases; 371,549 controls) (Mullins et al. 2021) https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594
grep -v '#' BD_Mullins2021.tsv | awk -F'\t' '$1 = 6' > BD_Mullins2021_chr6.tsv

## MDD
## (n = 2,074 cases; 2,925 controls) (Meng et al. 2024) https://figshare.com/articles/dataset/mdd2023diverse/24799299
sed 's/,/\t/g' MDD_meng2024.csv | awk -F'\t' '$1 = 6' > MDD_meng2024_chr6.csv
awk 'NR > 1' MDD_meng2024_chr6.csv | awk '$9 < 0.05' >  MDD_meng2024_chr6_p0.05.csv

## AD
## (n = 71,880 cases; 383,378 controls) (Jansen et al. 2019) https://www-nature-com.proxy-ub.rug.nl/articles/s41588-018-0311-9#additional-information

## intelligence
## (n = 269,867 cases;  controls) (Savage et al. 2018) https://www.nature.com/articles/s41588-018-0152-6#Sec23

## intersections with GWAS_snps.R

#### 2. LDSC regression

## install ldsc

git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create --file environment.yml
source activate ldsc

## get the ldsc regression and population data files

## snps list from hapmap project
wget https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/eur_w_ld_chr/w_hm3.snplist
awk '{if ($1!="SNP") {print $1} }' /home/shadi/Desktop/S4_Project/GWAS_snps/w_hm3.snplist > listHM3.txt
# or from https://github.com/perslab/CELLECT/blob/master/data/ldsc/w_hm3.snplist

## population baseline data
mkdir -p eur_w_ld_chr
cd eur_w_ld_chr
wget -r -l1 -H -nd -A .gz -e robots=off https://ibg.colorado.edu/cdrom2023/session/Day-3d%20LD%20Score%20heritability%20+%20Genetic%20Correlation,%20Andrew%20Grotzinger/LDSC_Practical_2023/EUR/eur_w_ld_chr/
wget -r -l1 -H -nd -np -e robots=off -A "*.M_5_50,*.ldscore.gz,README,w_hm3.snplist" "https://ibg.colorado.edu/cdrom2023/session/Day-3d%20LD%20Score%20heritability%20+%20Genetic%20Correlation,%20Andrew%20Grotzinger/LDSC_Practical_2023/EUR/eur_w_ld_chr/"

## download data from https://zenodo.org/records/10515792
## unzip 10515792.zip and then tar -xzvf its files
tar -xzvf 1000G_Phase3_plinkfiles.tgz

## munge the stats summary form GWAS data

## run munge_sumstats.py

## ASD

ldsc/munge_sumstats.py \
	--sumstats ASD_Grove2019.txt \
	--N 18381 \
	--out Munge_ASD_Grove2019 \
	--no-alleles \
	--snp SNP 

## SCZ

## remove # from the file

grep -v '#' SCZ_trubetskoy2022.tsv  > SCZ_trubetskoy2022_edited.tsv 

ldsc/munge_sumstats.py \
	--sumstats SCZ_trubetskoy2022_edited.tsv \
	--N 76755 \
	--out Munge_SCZ_trubetskoy2022 \
	--no-alleles \
	--snp ID 

## ldsc regression of the munged files

## run ldsc.py - h2

ldsc/ldsc.py \
  --h2 Munge_ASD_Grove2019.sumstats.gz \
  --ref-ld-chr 1000G_Phase3_ldscores/LDscore \
  --w-ld-chr 1000G_Phase3_ldscores/LDscore \
  --out ldsc_ASD_Grove2019

ldsc/ldsc.py \
  --h2 Munge_ASD_Grove2019.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out ldsc_ASD_Grove2019

## run ldsc.py - partitioned h2 per functional categories of binary annotation (published online from https://zenodo.org/records/10515792)

# ldsc/ldsc.py \
# 	--h2 Munge_ASD_Grove2019.sumstats.gz \
# 	--ref-ld-chr 1000G_Phase3_baselineLD_v2.3_ldscores/baseline. \
# 	--w-ld-chr 1000G_Phase3_weights_hm3_no_MHC/ \
# 	--overlap-annot \
# 	--frqfile-chr 1000G_Phase3_frq/ \
# 	--out ldsc_ASD_Grove2019_basline

## for ASD
ldsc/ldsc.py \
	--h2 Munge_ASD_Grove2019.sumstats.gz \
	--ref-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/baselineLD. \
	--w-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--frqfile-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out ldsc_ASD_Grove2019_basline
## extract UCSC functional categories 
(head -n 1 ldsc_ASD_Grove2019_basline.results && grep 'UCSC' ldsc_ASD_Grove2019_basline.results) > ASD_h2perUCSC.txt

## for SCZ
ldsc/ldsc.py \
	--h2 Munge_SCZ_trubetskoy2022.sumstats.gz \
	--ref-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/baselineLD. \
	--w-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--frqfile-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out ldsc_SCZ_trubetskoy2022_basline
## extract UCSC functional categories 
(head -n 1 ldsc_SCZ_trubetskoy2022_basline.results && grep 'UCSC' ldsc_SCZ_trubetskoy2022_basline.results) > SCZ_h2perUCSC.txt

## run ldsc.py - partitioned h2 per HAR_all, HAR_grik2, HAR_grik2_flanking of binary annotation

cut -f1,2,3,4 GSE180714_HARs.bed | grep -v -e 'chrX' -e 'chrY' | tail -n +2 > HAR_3100_nosex.bed
grep 'GRIK2' GSE180714_HARs.bed | cut -f1,2,3,4 | grep -v -e 'chrX' -e 'chrY' | tail -n +2 > HAR_grik2.bed 
cut -f 1,2,3,4 HAR_3100_grik2_flanking.bed >  HAR_grik2_flanking.bed

## create the annot file for HAR_all, HAR_grik2, and HAR_grik2_flanking and then calculate LD score
## code ref: https://kevinlkx.github.io/analysis_pipelines/sldsc_pipeline.html#tutorials_for_partitioning_heritability_using_s-ldsc_(stratified_ld_score_regression)

## using unix for HAR_all

for chrom in {1..22}
do
  echo ${chrom}

  ## Step 1: Creating an annot file
  echo "Make ldsc-friendly annotation files for ${ANNOT}.bed"
  python ldsc/make_annot.py \
  --bed-file HAR_ldsc/HAR_3100_nosex.bed \
  --bimfile 10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
  --annot-file HAR_ldsc/HAR_3100.${chrom}.annot.gz

  ## Step 2: Computing LD scores with an annot file
  echo "Computing LD scores with the annot file ${ANNOT}.${chrom}.annot.gz"
  python ldsc/ldsc.py \
  --l2 \
  --bfile 10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
  --print-snps listHM3.txt \
  --ld-wind-cm 1 \
  --annot HAR_ldsc/HAR_3100.${chrom}.annot.gz \
  --thin-annot \
  --out HAR_ldsc/HAR_3100.${chrom}
done

## using R for HAR_all

for chrom in {1..22}
do
  echo ${chrom}

  ## Step 1: Creating an annot file
  echo "Make ldsc-friendly annotation files for bed"
  Rscript ldsc/make_ldsc_binary_annot.R \
  HAR_ldsc/HAR_3100.bed \
  10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
  HAR_ldsc_R/HAR_3100_R.${chrom}.annot.gz "full-annot"
done
for chrom in {1..22}
do
  ## Step 2: Computing LD scores with an annot file
  echo "Computing LD scores with the annot file annot.gz"
  python ldsc/ldsc.py \
  --l2 \
  --bfile 10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
  --print-snps listHM3.txt \
  --ld-wind-cm 1 \
  --annot HAR_ldsc_R/HAR_3100_R.${chrom}.annot.gz \
  --out HAR_ldsc_R/HAR_3100_R.${chrom}
done

## using R for HAR_grik2

for chrom in {1..22}
do
  echo ${chrom}

  ## Step 1: Creating an annot file
  echo "Make ldsc-friendly annotation files for bed"
  Rscript ldsc/make_ldsc_binary_annot.R \
  HAR_ldsc/HAR_grik2.bed \
  10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
  HAR_ldsc_R/HAR_grik2_R.${chrom}.annot.gz "full-annot"
done
for chrom in {1..22}
do
  ## Step 2: Computing LD scores with an annot file
  echo "Computing LD scores with the annot file annot.gz"
  python ldsc/ldsc.py \
  --l2 \
  --bfile 10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
  --print-snps listHM3.txt \
  --ld-wind-cm 1 \
  --annot HAR_ldsc_R/HAR_grik2_R.${chrom}.annot.gz \
  --out HAR_ldsc_R/HAR_grik2_R.${chrom}
done

## using R for HAR_grik2_flanking

for chrom in {1..22}
do
  echo ${chrom}

  ## Step 1: Creating an annot file
  echo "Make ldsc-friendly annotation files for bed"
  Rscript ldsc/make_ldsc_binary_annot.R \
  HAR_ldsc/HAR_grik2_flanking.bed \
  10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom}.bim \
  HAR_ldsc_R/HAR_grik2_flanking_R.${chrom}.annot.gz "full-annot"
done
for chrom in {1..22}
do
  ## Step 2: Computing LD scores with an annot file
  echo "Computing LD scores with the annot file annot.gz"
  python ldsc/ldsc.py \
  --l2 \
  --bfile 10515792/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chrom} \
  --print-snps listHM3.txt \
  --ld-wind-cm 1 \
  --annot HAR_ldsc_R/HAR_grik2_flanking_R.${chrom}.annot.gz \
  --out HAR_ldsc_R/HAR_grik2_flanking_R.${chrom}
done

## ISSUE: Prop._h2 is always 1.0, refer to https://github.com/bulik/ldsc/issues/138

	## ISSUE_FIX_1: add a step 3 to the code which adds a base to thin anno files using
		## cat ${annotfile} | awk '{OFS="\t"; if (NR==1) {print "base",$1} else {print "1",$1} }' > ${annotfile_out}
		## cat HAR_grik2_R.6.annot | awk '{OFS="\t"; if (NR==1) {print "base",$1} else {print "1",$1} }' | gzip > HAR_grik2_R_based.6.annot.gz
	
	## ISSUE_FIX_2: edit code chunk #1 to #2 in ldsc/regressions.py script
		## code chunk #1
			## overlap_matrix_prop = np.zeros([self.n_annot,self.n_annot])
			## for i in range(self.n_annot):`
			##    overlap_matrix_prop[i, :] = overlap_matrix[i, :] / M_annot
		## code chunk #2
			## overlap_matrix_prop = np.zeros([self.n_annot,self.n_annot])
			## for i in range(self.n_annot):
			##    overlap_matrix_prop[i, :] = overlap_matrix[i, :] / M_tot

## run ldsc.py - HAR_all

## for ASD
ldsc/ldsc.py \
	--h2 Munge_ASD_Grove2019.sumstats.gz \
	--ref-ld-chr HAR_ldsc/HAR_3100. \
	--frqfile-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out ASD_HAR_baselineLD

## for SCZ
ldsc/ldsc.py \
	--h2 Munge_SCZ_trubetskoy2022.sumstats.gz \
	--ref-ld-chr HAR_ldsc/HAR_3100. \
	--frqfile-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out SCZ_HAR_baselineLD

## run ldsc.py - chr6

## for ASD
ldsc/ldsc.py \
	--h2 Munge_ASD_Grove2019.sumstats.gz \
	--ref-ld HAR_ldsc/HAR_3100.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out ASD_HAR_baselineLD_chr6_changesscipt

## for SCZ
ldsc/ldsc.py \
	--h2 Munge_SCZ_trubetskoy2022.sumstats.gz \
	--ref-ld HAR_ldsc/HAR_3100.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out SCZ_HAR_baselineLD_chr6_changesscipt 

## run ldsc.py - HAR_grik2

## for ASD
ldsc/ldsc.py \
	--h2 Munge_ASD_Grove2019.sumstats.gz \
	--ref-ld HAR_ldsc_R/HAR_grik2_R.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out ASD_grik2_baselineLD_chr6_changesscipt
	
## for SCZ
ldsc/ldsc.py \
	--h2 Munge_SCZ_trubetskoy2022.sumstats.gz \
	--ref-ld HAR_ldsc_R/HAR_grik2_R.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out SCZ_grik2_baselineLD_chr6_changesscipt

## run ldsc.py - HAR_grik2_flanking

## for ASD
ldsc/ldsc.py \
	--h2 Munge_ASD_Grove2019.sumstats.gz \
	--ref-ld HAR_ldsc_R/HAR_grik2_flanking_R.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out ASD_grik2_flanking_baselineLD_chr6_changesscipt
	
## for SCZ
ldsc/ldsc.py \
	--h2 Munge_SCZ_trubetskoy2022.sumstats.gz \
	--ref-ld HAR_ldsc_R/HAR_grik2_flanking_R.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out SCZ_grik2_flanking_baselineLD_chr6_changesscipt

## Rewrite the code as one chuck for one phenotype GWAS file

## for BD

## data preparation
sed '/^##/d' BD_Mullins2021.tsv | cat - <(grep '^##.*POS' BD_Mullins2021.tsv) > BD_Mullins2021_edited.tsv
## munge
ldsc/munge_sumstats.py \
	--sumstats BD_Mullins2021_edited.tsv \
	--N 41917 \
	--out Munge_BD_Mullins2021 \
	--no-alleles \
	--snp ID 
## h2 partitioned
ldsc/ldsc.py \
	--h2 Munge_BD_Mullins2021.sumstats.gz \
	--ref-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/baselineLD. \
	--w-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--frqfile-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC. \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out ldsc_BD_Mullins2021_basline
(head -n 1 ldsc_BD_Mullins2021_basline.results && grep 'UCSC' ldsc_BD_Mullins2021_basline.results) > BD_h2perUCSC.txt
## for HAR_all
ldsc/ldsc.py \
	--h2 Munge_BD_Mullins2021.sumstats.gz \
	--ref-ld-chr HAR_ldsc/HAR_3100. \
	--frqfile-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out BD_HAR_baselineLD
## for chr6
ldsc/ldsc.py \
	--h2 Munge_BD_Mullins2021.sumstats.gz \
	--ref-ld HAR_ldsc/HAR_3100.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out BD_HAR_baselineLD_chr6_changesscipt
## for HAR_grik2
ldsc/ldsc.py \
	--h2 Munge_BD_Mullins2021.sumstats.gz \
	--ref-ld HAR_ldsc_R/HAR_grik2_R.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out BD_grik2_baselineLD_chr6_changesscipt
## for HAR_grik2_flanking
ldsc/ldsc.py \
	--h2 Munge_BD_Mullins2021.sumstats.gz \
	--ref-ld HAR_ldsc_R/HAR_grik2_flanking_R.6 \
	--frqfile /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_frq/1000G.EUR.QC.6 \
	--w-ld /home/shadi/Desktop/S4_Project/GWAS_snps/10515792/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.6 \
	--overlap-annot --print-cov --print-coefficients --print-delete-vals \
	--out BD_grik2_flanking_baselineLD_chr6_changesscipt

#### Appendix I - Tools and scripts

## conda v24.1.2           
## MView v1.67
## biopython v1.84
## blast v2.14.0
## meme v5.5.5 & v5.0.2
## paml v4.10.7
## plink2 v2.00a5
## python v3.12.3
## ucsc-bigwigtowig v448
## ucsc-fatotwobit v455
## ucsc-fetchchromsizes v357
## vcf2maf v1.6.22
## vcfanno v0.3.3
## vcftools v0.1.17
## bcftools v1.19
## bedtools v2.25.0
## CLUSTER-BUSTER 2020-05-07 10:30:42 -0400  commit: ac1d33c
## RStudio 2024.04.2+764 "Chocolate Cosmos" Release (e4392fc9ddc21961fd1d0efd47484b43f07a4177, 2024-06-05) for Ubuntu Jammy rstudio/2024.04.2+764
## https://github.com/bulik/ldsc



