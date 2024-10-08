#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024

#### Content:
#### 1. positive selection in grik2 promoters in the human branch compared to other primates
#### Appendix I - Tools and scripts

#### Adaptation from adaptiPhy script from Wray lab
#### refs https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6752-4 https://github.com/wodanaz/adaptiPhy

#### 1. positive selection in grik2 promoters in the human branch compared to other primates

## grik2 promoter range = chr6:101391708-101393707

## 30way multiple alignments and hg38 reference - chr6 not genomewide

for chr in chr6; 
do wget --timestamping 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/multiz100way/maf/'$chr.maf.gz; done

for chr in chr6; 
do wget --timestamping 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/'$chr.fa.gz; done

## subset based on species

mkdir -p primates
maf_parse chr6.maf --seqs hg38,ponAbe2,panTro4,rheMac3 > chr6.primates.maf
msa_split chr6.primates.maf --refseq chr6.fa --gap-strip ANY -q --in-format MAF --features chr6.feat.bed --for-features

msa_split chr6.primates.maf --refseq chr6.fa --gap-strip ANY -q --in-format MAF --features chr6.feat_NFR.bed --for-features

## extract query alignments (aka your regulatory features)

mkdir features
for chr in chr6; 
	do grep -w $chr grik2_promoter_human.bed | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat.bed; 
done

for chr in chr6; 
	do grep -w $chr NFRchr6.sorted.bed | awk '{print $1 "\t" $2 "\t" $3 }' | sort -k1,1 -k2,2 -V >  features/$chr.feat_NFR.bed; 
done

for file in *promoter_msa_split.101391709-101393707.fa; do echo $file >> all.list;done
python prunning.py 
for file in *prunned; do echo $file >> all.prunned.list;done

python filtering.py 

cat goodalignments.txt > queries.list
cp queries.list /neutral
cd /neutral

python bfgenerator_global.py query/queries.list

cat /home/shadi/Desktop/S4_Project/HyPhy_grik2reg/query/* > queries.list
cp queries.list /home/shadi/Desktop/S4_Project/HyPhy_grik2reg/query
cd  /home/shadi/Desktop/S4_Project/HyPhy_grik2reg/query

/home/shadi/Desktop/S4_Project/HyPhy_grik2reg/reference

for file in *.null.bf; do echo $file >> null.hg38.list; done
for file in *.alt.bf; do echo $file >> alt.hg38.list; done

for file in *hg38.null.res; do echo $file >> null.hg38.log; done;
for file in *hg38.alt.res; do echo $file >> alt.hg38.log; done;

for filename in `cat alt.hg38.log`; do grep -H "BEST LOG-L:" $filename >> alt.hg38.tab; done;
for filename in `cat null.hg38.log`; do grep -H "BEST LOG-L:" $filename >> null.hg38.tab; done;

for file in *panTro5.null.res; do echo $file >> null.panTro5.log; done;
for file in *panTro5.alt.res; do echo $file >> alt.panTro5.log; done;

for filename in `cat alt.panTro5.log`; do grep -H "BEST LOG-L:" $filename >> alt.panTro5.tab; done;
for filename in `cat null.panTro5.log`; do grep -H "BEST LOG-L:" $filename >> null.panTro5.tab; done;

awk '{print $1 "\t" $3}' null.hg38.tab | sort -k1,1 -V  > nulls.hg38.tab
awk '{print $1 "\t" $3}' alt.hg38.tab | sort -k1,1 -V > alts.hg38.tab

awk '{print $1 "\t" $3}' null.panTro5.tab | sort -k1,1 -V  > nulls.panTro5.tab
awk '{print $1 "\t" $3}' alt.panTro5.tab | sort -k1,1 -V > alts.panTro5.tab

awk -F"." '{print $1 ":" $2 "\t" $3 }'  nulls.hg38.tab > col1.hg38.tab
paste col1.hg38.tab nulls.hg38.tab  | awk '{print $1 "\t" $2 "\t" $4   }' > lnulls.hg38.tab
paste lnulls.hg38.tab alts.hg38.tab | awk '{print $1 "\t" $2  "\t" $3 "\t" $5  }' |  sort -k1,1 -V > likelihoods.hg38.tab
 
awk -F"." '{print $1 ":" $2 "\t" $3 }'  nulls.panTro5.tab > col1.panTro5.tab
paste col1.panTro5.tab nulls.panTro5.tab  | awk '{print $1 "\t" $2 "\t" $4   }'> lnulls.panTro5.tab
paste lnulls.panTro5.tab alts.panTro5.tab | awk '{print $1 "\t" $2  "\t" $3 "\t" $5  }' |  sort -k1,1 -V > likelihoods.panTro5.tab

cat likelihoods.hg38.tab likelihoods.panTro5.tab  |  sort -k1,1 -V | sed 1i"chromosome\tbranch\tlnull\tlalt" > likelihoods.tab

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



