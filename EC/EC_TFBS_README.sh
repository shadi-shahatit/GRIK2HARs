#### Transcription factor binding site analysis of FRMPD2 and FRMPD2B promoters
#### September 2024, for Eirini

#### Content:
#### 1. The meme suite
####	1.1 create cortical synaptogenesis TF list    
#### 	1.2 run FIMO
####    1.3 extract results
####        - top scoring (1% & 5%)
####        - uniqueness (complete & partial)
####    1.4 combine to make life easier
#### 2. Cluster-Buster
####	2.1 create cortical synaptogenesis TF list and cbust file prep
####	2.2 run cbust
####	2.3 extract results (top scoring cluster or cis regulatory module)
#### 3. intersect with genomic regions
#### 	3.1 HAQER

#### 1. The meme suite

####	1.1 create cortical synaptogenesis TF list    

## filter JASPAR2024 based on cortical tissue expression for human and macaque during synaptogenesis. From PsychENCODE data, we filtered genes for specific tissues and developmental time points (aka, synaptogenesis_genes_human.txt, synaptogenesis_genes_macaque.txt). Now, from this genes list, we only consider the TFs labelled from JASPAR to produce JASPAR2024_synaptogenesis_genes_human and JASPAR2024_synaptogenesis_genes_macaque in meme and jaspar format for fimo and cbust

## kept 510 out of 879 TFs in humans
## kept 384 out of 879 TFs in macaques

## for more details, this is the code snipt for how to do it (EXTRA - can be skipped):

## downlaod JASPAR2024 single batch file (txt) with meme format from https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

## JASPAR2024_synaptogenesis_genes in humans

genes=$(awk '{print tolower($0)}' synaptogenesis_genes_human.txt | paste -sd'|' -)
echo "$genes"
awk -v genes="$genes" '
BEGIN {IGNORECASE=1; RS=""; FS="\n"}
/^MOTIF/ {
    split($0, motif_line, " ");  # Split the MOTIF line into an array by spaces
    gene_name = tolower(motif_line[3]);  # Get the third field (gene name)
    if (gene_name ~ genes) {
        print $0;  # If the gene name matches, print the entire motif block
    }
}
' JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt > JASPAR2024_synaptogenesis_meme.txt
## manually add "MEME version 4 ALPHABET= ACGT strands: + - Background letter frequencies A 0.25 C 0.25 G 0.25 T 0.25" JASPAR2024_synaptogenesis_meme.txt > JASPAR2024_synaptogenesis_genes_human_meme.txt
grep "^MOTIF" JASPAR2024_synaptogenesis_meme.txt | awk '{print$3}' > JASPAR2024_synaptogenesis_meme_genesym.txt
grep -A 4 -e "" -e "" -e "" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt > JASPAR2024_synaptogenesis_genes_human.txt
grep ">" JASPAR2024_synaptogenesis_genes_human.txt | wc -l
grep ">" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt | wc -l
# filtered 879 to 510 

## JASPAR2024_synaptogenesis_genes in macaque

genes=$(awk '{print tolower($0)}' synaptogenesis_genes_macaque.txt | paste -sd'|' -)
echo "$genes"
awk -v genes="$genes" '
BEGIN {IGNORECASE=1; RS=""; FS="\n"}
/^MOTIF/ {
    split($0, motif_line, " ");  # Split the MOTIF line into an array by spaces
    gene_name = tolower(motif_line[3]);  # Get the third field (gene name)
    if (gene_name ~ genes) {
        print $0;  # If the gene name matches, print the entire motif block
    }
}
' JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt > JASPAR2024_synaptogenesis_meme_macaque.txt
## manually add "MEME version 4 ALPHABET= ACGT strands: + - Background letter frequencies A 0.25 C 0.25 G 0.25 T 0.25" JASPAR2024_synaptogenesis_meme_macaque.txt > JASPAR2024_synaptogenesis_genes_macaque_meme.txt
grep "^MOTIF" JASPAR2024_synaptogenesis_meme_macaque.txt | awk '{print$3}' > JASPAR2024_synaptogenesis_meme_genesym_macaque.txt
grep -A 4 -e "" -e "" -e "" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt > JASPAR2024_synaptogenesis_genes_macaque.txt
grep ">" JASPAR2024_synaptogenesis_genes_macaque.txt | wc -l
grep ">" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt | wc -l
# filtered 879 to 384

#### 	1.2 run FIMO

## get the promoter seqeunces in fasta format

FRMPD2Bpromoter.txt
6exFRMPD2humanpromoter.txt
6exFRMPD2macaquepromoter.txt 

## For FRMPD2Bpromoter.txt
fimo -oc fimo_out_FRMPD2Bpromoter JASPAR2024_synaptogenesis_genes_human_meme.txt FRMPD2Bpromoter.txt 

## For 6exFRMPD2humanpromoter.txt
fimo -oc fimo_out_6exFRMPD2humanpromoter JASPAR2024_synaptogenesis_genes_human_meme.txt 6exFRMPD2humanpromoter.txt

## For 6exFRMPD2macaquepromoter.txt
fimo -oc fimo_out_6exFRMPD2macaquepromoter JASPAR2024_synaptogenesis_genes_macaque_meme.txt 6exFRMPD2macaquepromoter.txt

####	1.3 extract results (unique and top scoring)
####    1.4 combine to make life easier

## see EC_TFBS.R

#### 2. Cluster-Buster

####	2.1 create cortical synaptogenesis TF list and cbust file prep

## use JASPAR2024_synaptogenesis_genes_human and JASPAR2024_synaptogenesis_genes_macaque in jaspar format for cbust

## cbust prep

sed '/^--$/d' JASPAR2024_synaptogenesis_genes_human.txt > JASPAR2024_synaptogenesis_genes_human_mod.txt
awk '/^>/ {print; next} {gsub(/^[ACGT] /,""); gsub(/\[|\]/,""); print}' JASPAR2024_synaptogenesis_genes_human_mod.txt > JASPAR2024_pfms_synaptogenesis_cbust.txt

sed '/^--$/d' JASPAR2024_synaptogenesis_genes_macaque.txt > JASPAR2024_synaptogenesis_genes_macaque_mod.txt
awk '/^>/ {print; next} {gsub(/^[ACGT] /,""); gsub(/\[|\]/,""); print}' JASPAR2024_synaptogenesis_genes_macaque_mod.txt > JASPAR2024_pfms_synaptogenesis_macaque_cbust.txt

####	2.1 run cbust

./cluster-buster/cbust -f 1 -c0 -m0 JASPAR2024_pfms_synaptogenesis_cbust.txt FRMPD2Bpromoter.txt > cbustout_FRMPD2Bpromoter_1.txt
./cluster-buster/cbust -f 1 -c0 -m0 -l JASPAR2024_pfms_synaptogenesis_cbust.txt 6exFRMPD2humanpromoter.txt  > cbustout_6exFRMPD2humanpromoter_1.txt
./cluster-buster/cbust -f 1 -c0 -m0 -l JASPAR2024_pfms_synaptogenesis_macaque_cbust.txt 6exFRMPD2macaquepromoter.txt  > cbustout_6exFRMPD2macaquepromoter_1.txt

## other flags for cbust: -g3 -c5 - m3 , -c0.01 -g20 , -c0 -m0, -c0 -m0 -l

## this gives cbust output in format 1 (a good format to read the results), but still the files are not clean. So, we use some unix manipulation to clean them and attach proper gene names we can recognize rather than motif IDs (we also make sure that the # are in the unread unneeded lines)

grep "^MOTIF" JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt | awk '{print $2,$3}' > gene_mapping.txt

gene_mapping_file="gene_mapping.txt"
input_file="cbustout_FRMPD2Bpromoter_1.txt"
output_file="cbustout_FRMPD2Bpromoter_1_symbols.txt"
temp_file=$(mktemp)
awk '{gsub(/\/|&/, "\\&"); print "s/" $1 "/" $2 "/g"}' "$gene_mapping_file" > sed_commands.txt
sed -f sed_commands.txt "$input_file" > "$temp_file"
mv "$temp_file" "$output_file"
rm sed_commands.txt

gene_mapping_file="gene_mapping.txt"
input_file="cbustout_6exFRMPD2humanpromoter_1.txt"
output_file="cbustout_6exFRMPD2humanpromoter_1_symbols.txt"
temp_file=$(mktemp)
awk '{gsub(/\/|&/, "\\&"); print "s/" $1 "/" $2 "/g"}' "$gene_mapping_file" > sed_commands.txt
sed -f sed_commands.txt "$input_file" > "$temp_file"
mv "$temp_file" "$output_file"
rm sed_commands.txt

gene_mapping_file="gene_mapping.txt"
input_file="cbustout_6exFRMPD2macaquepromoter_1.txt"
output_file="cbustout_6exFRMPD2macaquepromoter_1_symbols.txt"
temp_file=$(mktemp)
awk '{gsub(/\/|&/, "\\&"); print "s/" $1 "/" $2 "/g"}' "$gene_mapping_file" > sed_commands.txt
sed -f sed_commands.txt "$input_file" > "$temp_file"
mv "$temp_file" "$output_file"
rm sed_commands.txt

####	2.3 extract results (top scoring cluster or cis regulatory module)

## see EC_TFBS.R

#### 3. intersect with genomic regions

#### 	3.1 HAQER

bedtools sort -i HAQER_chr10.bed > HAQER_chr10_sorted.bed

bedtools intersect -a FRMPD2B_hg38ref.bed -b HAQER_chr10_sorted.bed > intersect_FRMPD2B_HAQER.bed
bedtools closest -a FRMPD2B_hg38ref.bed -b HAQER_chr10_sorted.bed > closest_FRMPD2B_HAQER.bed

## no intersections at all but one 867 bp long HAQER is close (197157 bp from TSS)




