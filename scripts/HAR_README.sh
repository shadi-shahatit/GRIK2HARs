#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024

#### Content:
#### 1. extracting HARs multiple sequence alignment and plasmid order  
#### 	1.1 MView
#### 	1.2 pattern serach for the HAR coordinates
#### 	1.3 twobittofa and orthologous seq extraction
#### 	1.4 bed files of HAR snps
#### 	1.5 order HAR plasmids from idt for human
#### 2. Neanderthal and Denisovan in HARs
#### 3. Transcription factor binding site analysis
#### 	3.1 The meme suite
#### 	3.2 Cluster-Buster
#### 	3.3 motifbreakR
#### Appendix I - Tools and scripts

#### extracting HARs multiple sequence alignment and plasmid order  

#### 1.1 MView

## install MView

tar xvzf mview-1.67.tar.gz
perl install.pl
mview -help

## run MView

mview -in chr6.maf -html head -bold -css on -coloring any chr6.maf > data_chr6.html

mview -in chr6.maf -html head -bold -css on -coloring any -block 2000 chr6.maf  > data_chr6_block_2.html

#### 1.2 pattern serach for the HAR coordinates

grep -E '\b(101296272|10129627[3-9]|101296[3-4][0-9]{2}|1012965[0-3][0-9]|10129654[0-6])\b' chr6.maf 

grep 'aattgtttcaagatgtacagatgcgagtgagaagagtgagtaaaacattgagcttgtgaaagaaggaagatgagtgaatttctcttattgggcttctctccagaaagagtcatctggcgagttggcaaagacaaaatcaaacaaccatgtgccagatgtgcacggctgttgttgggccatctgcagcttgcttttatgccatatgccacttttacatttggtcggctggggtatcagaggcaaaattcatggggcaattgcttttgtattattct' chr6.maf ## Not working!

#### 1.3 twobittofa and orthologous seq extraction

## download twobittofa from http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod 744 twoBitToFa

## extract the seq

./twoBitToFa https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit -bed=HAR_3100_grik2.bed stdout > HAR_human_seq.fa
./twoBitToFa https://hgdownload.soe.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.2bit -bed=chimpliftover_Girskis_3100_grik2.bed stdout > HAR_chimp_liftover_seq.fa
./twoBitToFa https://hgdownload.soe.ucsc.edu/goldenPath/panPan2/bigZips/panPan2.2bit -bed=bonoboliftover_Girskis_3100_grik2.bed stdout > HAR_bonobo_liftover_seq.fa

./twoBitToFa https://hgdownload.soe.ucsc.edu/goldenPath/panPan2/bigZips/panPan2.2bit -bed=macaqueliftover_Girskis_3100_grik2.bed stdout > HAR_macaque_liftover_seq.fa

## sepearte each seq to a fasta file

awk '/^>/{header=$0; gsub("[>\x27]", "", header); if (out) close(out); out=header".txt"; print > out; next} {print > out}' HAR_human_seq.fa 
awk '/^>/{header=$0; gsub("[>\x27]", "", header); if (out) close(out); out=header".txt"; print > out; next} {print > out}' HAR_chimp_liftover_seq.fa
awk '/^>/{header=$0; gsub("[>\x27]", "", header); if (out) close(out); out=header".txt"; print > out; next} {print > out}' HAR_bonobo_liftover_seq.fa

awk '/^>/{header=$0; gsub("[>\x27]", "", header); if (out) close(out); out=header".txt"; print > out; next} {print > out}' HAR_macaque_liftover_seq.fa

## combine the human and ortho seq in jalview
## extract human SNPs

directory="/home/shadi/Desktop/S4_Project/HAR/HAR_maf_vcf"
for file in "$directory"/HARsv2_25*.fa; do
    if [ -f "$file" ] && [ -r "$file" ]; then
    	snp-sites "$file" -v > "${file%.fa}.vcf"
    else
        echo "Error: $file does not exist or is not readable."
    fi
done

directory="/home/shadi/Desktop/S4_Project/HAR/HAR_maf_vcf"
HAR_positions=(101296272 101518375 101931326 102249110 103398754 103468078 103543581)
for i in {14..20}; do
    file="$directory/HARsv2_25${i}.vcf"
    if [ -f "$file" ] && [ -r "$file" ]; then
        num=${HAR_positions[$((i - 14))]}
        awk -v num="$num" 'BEGIN {OFS="\t"} {
            if ($1 ~ /^#/ || $1 == "") {
                print $0;
            } else {
                $1 = "6";
                $2 += num;
                print $0;
            }
        }' "$file" > "${file%.vcf}_modified.vcf"
    else
        echo "Error: $file does not exist or is not readable."
    fi
done

## check whether HARsv2_2520 seq make sense

## download cactus 447-way MAF for the hg38 chr6:103543581-103544084 from table browser and view with mview

mview -in maf447_HARsv2_2520_human.maf -html head -bold -css on -coloring any maf447_HARsv2_2520_human.maf > maf447_HARsv2_2520_human.html
maf_parse maf447_HARsv2_2520_human.maf --order hg38.chr6,Pan_troglodytes.CM009244.2,Pan_paniscus.CM003389.1,Macaca_mulatta.NC_041757.1,Mus_musculus.chr10 > maffiltered_HARsv2_2520_human.maf
mview -in maffiltered_HARsv2_2520_human.maf -html head -bold -css on -coloring mismatch maffiltered_HARsv2_2520_human.maf > maffiltered_HARsv2_2520_human.html

## to get the right coordinates for HARsv2_2520, extract variations from cactus 447-way maf file

## use maf2vcf.pl script from ckandoth
## download hg38 fasta seq
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
perl maf2vcf.pl --input-maf maffiltered_HARsv2_2520_human.maf --output-dir HARsv2_2520_cactus_vcfs --ref-fasta /home/shadi/Desktop/S4_Project/HAR/HAR_maf_vcf/hg38.fa
## use msa2vcf from jvarkit
## install java in your source
sudo apt install default-jre
sudo apt install openjdk-17-jre-headless
grep -v -e '##maf version=1' -e '# maf_parse --order hg38.chr6,Pan_troglodytes.CM009244.2,Pan_paniscus.CM003389.1,Macaca_mulatta.NC_041757.1,Mus_musculus.chr10 maf447_HARsv2_2520_human.maf' -e '#eof' -e 'score' -e 'C 0' maffiltered_HARsv2_2520_human.maf | cut -d' ' -f2- | awk '{gsub(/[0-9+]/, ""); gsub(/(^|[^-])-(?![^-]|$)/, ""); print}' > test.maf
## add "CLUSTAL W (1.81) multiple sequence alignment"
java -jar jvarkit.jar msa2vcf test.maf 

## NOT WORKING

## use PHAST instead to create a fasta file from the maf file and reiterate the snp-sites loop to get a vcf file

maf_parse maffiltered_HARsv2_2520_human.maf --order hg38.chr6,Pan_troglodytes.CM009244.2,Pan_paniscus.CM003389.1 > maf3species_HARsv2_2520_human.maf
maf_parse maf3species_HARsv2_2520_human.maf -p
msa_view maf3species_HARsv2_2520_human.maf -o FASTA > HARsv2_2520_cactus.fa

## be careful of the start site coordinate (103543581), the maf starts with 103543580 so remove manually (or with awk) the first nucleotide which is A (not part of HAR, and also conserved)

## extract human SNPs for HARsv2_2520_cactus.fa

snp-sites HARsv2_2520_cactus.fa -v > HARsv2_2520.vcf

directory="/home/shadi/Desktop/S4_Project/HAR/HAR_maf_vcf"
HAR_position=103543581
file="$directory/HARsv2_2520.vcf"
awk -v num="$HAR_position" 'BEGIN {OFS="\t"} {
    if ($1 ~ /^#/ || $1 == "") {
        print $0;
    } else {
        $1 = "6";
        $2 += num;
        print $0;
    }
}' "$file" > "${file%.vcf}_modified.vcf"

#### 1.4 bed files of HAR snps

## create a bed file from the vcf files of HAR snps
## note that this output is tailored for motifbreakR package and snps.from.file function where start and end are the same. In a UCSC bed file format the start is one base before the seqeunce

output_bed="merged_HAR_snp.bed"
rm -f $output_bed

for i in {2514..2520}
do
    vcf_file="HARsv2_${i}_modified.vcf"
    if [ -f "$vcf_file" ]; then
        awk 'BEGIN {OFS="\t"} 
            !/^#/ {
                start = $2-1;
                end = start+1;
                custom_field = "chr"$1":"start":"$4":"$5;
                score = 0;
                strand = "+";
                print "chr"$1, start, end, custom_field, score, strand
            }' $vcf_file >> $output_bed
    else
        echo "File $vcf_file does not exist."
    fi
done

#### 1.5 order HAR plasmids from idt for human

## create a new HAR_human_seq.fa file with new coordinates that have potential regulatory flanking regions

awk 'BEGIN {
    OFS = "\t";
    print "chr", "start_flanking", "end_flanking", "id", "size", "percentage", "start", "end"
}
NR > 1 {
    size = $3-$2;
    percentage = int(size*0.25)+1;
    start_flanking = $2 - percentage;
    end_flanking = $3 + percentage;
    print $1, int(start_flanking), int(end_flanking), $4, size, percentage, $2, $3;
}' HAR_3100_grik2.bed > HAR_3100_grik2_flanking_head.bed
sed '1d' HAR_3100_grik2_flanking_head.bed > HAR_3100_grik2_flanking.bed
cat HAR_3100_grik2_flanking.bed

## extract the seq
./twoBitToFa https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit -bed=HAR_3100_grik2_flanking.bed stdout > HAR_flanking_human_seq.fa

## check if restriction sites of SpeI (ACTAGT) and SnaBI (TACGTA) are inside the target seq
grep -e 'ACTAGT' -e 'TACGTA' HAR_flanking_human_seq.fa
## ACTAGT is in HAR17 and HAR20
## what about SalI (GTCGAC) and and SnaBI (TACGTA)
grep -e 'GTCGAC' -e 'TACGTA' HAR_flanking_human_seq.fa
## add restriction sites with 4-5 random bps (e.g., CAGTT)
awk '/^>/ {print $0; next} {print "CAGTTGTCGAC"$0"TACGTACAGTT"}' HAR_flanking_human_seq.fa > HAR_flanking_human_seq_restsites.txt

#### 2. Neanderthal and Denisovan in HARs

## Neanderthal

## Neanderthal SNPs that underwent a selective sweep in early humans - 3 data files of snps genotypes, snps +ve selection Z score, lowest 5% score

wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ntSssSnps.txt.gz
gunzip ntSssSnps.txt.gz
awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' ntSssSnps.txt | sort -k1,1 -k2,2n | sed 's/ \+/\t/g' > ntSssSnps_edited.bed
bedtools intersect -a ntSssSnps_edited.bed -b hg19liftover_Girskis_3100_grik2.bed -wa -wb > ntSssSnps_intersected.bed
## could be -wo -wa wb
## or just directly from the ucsc table tool
## output 4 SNPs

wget https://hgdownload.cse.ucsc.edu/gbdb/hg19/neandertal/bbi/ntSssZScorePMVar.bw
chmod a+x bigWigToWig 
./bigWigToWig ntSssZScorePMVar.bw ntSssZScorePMVar.wig
wig2bed < ntSssZScorePMVar.wig > ntSssZScorePMVar.bed
awk '{$2=$2+1; print }' ntSssZScorePMVar.bed | sort -k1,1 -k2,2n | sed 's/ \+/\t/g' > ntSssZScorePMVar_0based.bed
bedtools intersect -a ntSssZScorePMVar_0based.bed -b hg19liftover_Girskis_3100_grik2.bed -wa -wb > ntSssZScorePMVar_intersected.bed
## output 4 SNPs
## first one is a red SNP in which at least four of the six modern human genomes are derived while and all observed Neandertal alleles are ancestral

wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ntSssTop5p.txt.gz
gunzip ntSssTop5p.txt.gz
awk '{print $2,$3,$4,$5,$6}' ntSssTop5p.txt | sort -k1,1 -k2,2n | sed 's/ \+/\t/g' > ntSssTop5p_edited.bed
bedtools intersect -a ntSssTop5p_edited.bed -b hg19liftover_Girskis_3100_grik2.bed -wa -wb > ntSssTop5p_intersected.bed
## output 0 SNPs

## look for neanderthal_snps in HAR.R

## Denisovan

## Modern derived - 3 data files of fixed, fixed and dbSNP, high frequency snps

wget https://hgdownload.cse.ucsc.edu/gbdb/hg19/dhcHumDerDenAnc/dhcHumDerDenAncAllFixed.bb
wget https://hgdownload.cse.ucsc.edu/gbdb/hg19/dhcHumDerDenAnc/dhcHumDerDenAncAllFixedDbSnp.bb
wget https://hgdownload.cse.ucsc.edu/gbdb/hg19/dhcHumDerDenAnc/dhcHumDerDenAncAllHighFreq.bb

chmod a+x bigBedToBed 
./bigBedToBed dhcHumDerDenAncAllFixed.bb dhcHumDerDenAncAllFixed.bed
./bigBedToBed dhcHumDerDenAncAllFixedDbSnp.bb dhcHumDerDenAncAllFixedDbSnp.bed
./bigBedToBed dhcHumDerDenAncAllHighFreq.bb dhcHumDerDenAncAllHighFreq.bed

## focus on fixed and high frequency

bedtools intersect -a dhcHumDerDenAncAllFixed.bed -b hg19liftover_Girskis_3100_grik2.bed -wa -wb > FixedDenisovan_intersected.bed
bedtools intersect -a dhcHumDerDenAncAllHighFreq.bed -b hg19liftover_Girskis_3100_grik2.bed -wa -wb > HighFreqDenisovan_intersected.bed
bedtools intersect -a dhcHumDerDenAncAllFixedDbSnp.bed -b hg19liftover_Girskis_3100_grik2.bed -wa -wb > FixedDbSnpDenisovan_intersected.bed
## output 0 SNPs

## Denisovan variants

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfDenisovaPinky.txt.gz

bedtools intersect -b HARsv2_2514_modified.vcf  HARsv2_2515_modified.vcf  HARsv2_2516_modified.vcf  HARsv2_2517_modified.vcf  HARsv2_2518_modified.vcf  HARsv2_2519_modified.vcf HARsv2_2520_modified.vcf -a HAR_dhcVcfDenisovaPinky.vcf -wa -wb > Denisovan_intersected.bed ## no intersections whatsoever

## keep Denisovan variants file (HAR_dhcVcfDenisovaPinky.vcf) to consider its snps

## Modern human variants 10+1 sample

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfDNK02.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP00456.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP00521.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP00542.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP00665.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP00778.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP00927.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP00998.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP01029.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP01284.txt.gz
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/dhcVcfHGDP01307.txt.gz

## intersect ucsc table tool files with HAR

directory="/home/shadi/Desktop/S4_Project/HAR/Denisovan"

for file in "$directory"/HAR_dhcVcf*.bed; do
    if [ -f "$file" ] && [ -r "$file" ]; then
        bedtools intersect -b "$file" -a hg19liftover_Girskis_3100_grik2.bed -wa -wb > "${file%.bed}_intersected.bed"
    else
        echo "Error: $file does not exist or is not readable."
    fi
done

## intersect ucsc table tool files with HAR at once

bedtools intersect -a hg19liftover_Girskis_3100_grik2.bed -b HAR_dhcVcfHGDP00521.bed HAR_dhcVcfHGDP01284.bed  HAR_dhcVcfHGDP00927.bed  HAR_dhcVcfHGDP00542.bed  HAR_dhcVcfHGDP01307.bed  HAR_dhcVcfHGDP00998.bed  HAR_dhcVcfHGDP00665.bed  HAR_dhcVcfHGDP00456.bed  HAR_dhcVcfHGDP01029.bed  HAR_dhcVcfHGDP00778.bed HAR_dhcVcfDNK02.bed -wa -wb -sorted -names HGDP00521_French HGDP01284_Mandenka  HGDP00927_Yoruba  HGDP00542_Papuan  HGDP01307_Dai HGDP00998_Karitiana  HGDP00665_Sardinian  HGDP00456_Mbuti  HGDP01029_San  HGDP00778_Han DNK02_Dinka | sort -k1,1 -k7,7n > superDenisovan_intersected.bed 

## inclde human diversity into account to put Denisovan variants into perspective for HAR_dhcVcfDenisovaPinky.vcf

awk '{print$6,$7,$8,$9,$5,$4}' superDenisovan_intersected.bed | sort -k1,1 -k2,2n | sed 's/ \+/\t/g' > superDenisovan_intersected_edited.bed  
grep -v '#' HAR_dhcVcfDenisovaPinky.vcf | awk '{print$1,$2-1,$2,$3}' | awk '{ if ($1 == "6") $1 = "chr6"; print }' | sort -k1,1 -k2,2n | sed 's/ \+/\t/g' > HAR_dhcVcfDenisovaPinky_edited.bed 
## revisit the 0-based versus 1-based coordinates

bedtools intersect -b superDenisovan_intersected_edited.bed -a HAR_dhcVcfDenisovaPinky_edited.bed -wa -wb > Denisovan_Pinky_intersected.bed

## used -v to check intervals that do not intersect - output 0 interval

## look for denisovan_snps in HAR.R

#### 3. Transcription factor binding site analysis

## JASPAR2024

## downlaod JASPAR2024 sites from ucsc table tool with hg38 HAR_grik2 coordinates
## filter significant scores (<400) and count
awk '$5 > 400' jaspar24_HARgrik2_38.bed | cut -f7 | sort | uniq -c | sort -nr > jaspar24_HARgrik2_38_counted.bed

#### 3.1 The meme suite

## FIMO

## downlaod JASPAR2024 individual PFMs files with meme format from https://jaspar.elixir.no/downloads/
## downlaod JASPAR2024 single batch file (txt) with meme format from https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt 
## FOR VERTEBRATES

fimo JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt  HAR_human_seq.fa # worked :)
sea --m JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt  --p HAR_human_seq.fa ## not working :(
sea --m JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt  --p HAR_human_seq.fa --no-pgc ## not working :(

## install meme directly from Source ## not working :(

## downlaod meme file tar version from https://meme-suite.org/meme/doc/download.html
tar zxf meme-5.5.5.tar.gz
          cd meme-5.5.5
          ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
          make
          make test
          make install
## downlaod dependencies    
cd scripts/
perl dependencies.pl 
sudo cpan Math::CDF
sudo cpan HTML::Template 
sudo cpan Sys::Info
sudo cpan XML::Simple
sudo cpan Log::Log4perl 
sudo cpan XML::Compile::SOAP11
sudo cpan XML::Compile::WSDL11
sudo cpan XML::Compile::Transport::SOAPHTTP
sudo apt-get install libexpat1-dev
sudo apt-get install libxml2-dev
sudo apt-get install libicui18n.so.58

## install meme older version in separate conda env

conda install -c bioconda/label/cf201901 meme

centrimo HAR_human_seq.fa JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
centrimo -verbosity 1 -oc centrimo_out_2 -score 5.0 -ethresh 10.0 HAR_human_seq.fa JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

dreme-py3 -p HAR_human_seq.fa

meme HAR_human_seq.fa -dna -mod oops -revcomp -nmotifs 1000

## run sea on the meme suite website https://meme-suite.org/meme/tools/sea

## FIMO - chimp

fimo -oc fimo_out_chimp JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt HAR_chimp_liftover_seq.fa

## TOMTOM

meme HAR_human_seq.fa -dna -mod oops -revcomp -nmotifs 1000
meme -oc meme_out_chimp HAR_chimp_liftover_seq.fa -dna -mod oops -revcomp -nmotifs 1000

tomtom JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt ./meme_out_chimp/meme.html ./meme_out/meme.html
tomtom ./meme_out_chimp/meme.html ./meme_out/meme.html

#### 3.2 Cluster-Buster 

## old as 

## install cluster-buster directly from Source

git clone https://github.com/weng-lab/cluster-buster
cd cluster-buster
make cbust
./cbust -h

## run cluster-buster

## downlaod JASPAR2024 single batch file (txt) with jaspar format from https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt

awk '/^>/ {print; next} {gsub(/^[ACGT] /,""); gsub(/\[|\]/,""); print}' JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt > JASPAR2024_pfms_cbust.txt

./cluster-buster/cbust -f 1 -c0 -m0 -l JASPAR2024_pfms_cbust.txt HAR_human_seq.fa > cbust_out_1.txt
./cluster-buster/cbust -f 2 -c0 -m0 -l JASPAR2024_pfms_cbust.txt HAR_human_seq.fa > cbust_out_2.txt
./cluster-buster/cbust -f 3 -c0 -m0 -l JASPAR2024_pfms_cbust.txt HAR_human_seq.fa > cbust_out_3.txt

## other flags: -g3 -c5 - m3 , -c0.01 -g20 , -c0 -m0

grep 'Score' cbust_out_2.txt

## cluster 1
grep -e 'MA1155.1' -e 'MA0884.2' -e 'MA0629.2' JASPAR2024_pfms_cbust.txt 
## cluster 2
grep -e 'MA0649.2' -e 'MA1147.2' -e 'MA1146.2' JASPAR2024_pfms_cbust.txt 
## cluster 3
grep -e 'MA0756.3' JASPAR2024_pfms_cbust.txt 

## add TF id

awk '{if ($1 ~ /^>/) {gsub(/ +/, "", $1); $1=$1"-"$2; $2=""; sub(/^-/, ""); } print}' JASPAR2024_pfms_cbust.txt > JASPAR2024_pfms_cbust_id.txt
./cluster-buster/cbust -f 2 -c0 -m0 -l JASPAR2024_pfms_cbust_id.txt HAR_human_seq.fa > cbust_out_id2.txt

#### 3.3 motifbreakR

## Forge a BSgenome package for the HAR elements
## https://tomguest.netlify.app/tutorial/bsgenome/

## check TF from JASPAR in grik2 co expression module

grep "^MOTIF" TFBS/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt | awk '{print$3}' > TFBS/TF_genesym_JASPAR2024.txt
grep -Fxf grik2_module_gene_names.txt  TFBS/TF_genesym_JASPAR2024.txt

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



