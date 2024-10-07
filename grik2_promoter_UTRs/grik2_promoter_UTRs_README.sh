#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024

#### Content:
#### 1. Promoter and 5' 3' UTRs of grik2  
#### 	1.1 Promoter
#### 	1.2 5' 3' UTRs
#### 	1.3 order promoter plasmids from idt
#### 2. Transcription factor binding site analysis for promoter
#### 	2.1 The meme suite
#### 	2.2 Cluster-Buster
#### Appendix I - Tools and scripts

#### 1. Promoter and 5' 3' UTRs of grik2

#### 1.1 Promoter

## USCS table browser

## NOT for grik2 coordinates from ensembl
## BUT for grik2 coordinates from ensembl canonical or ncbi database
## choose seqeunce then promoter/upstream seq of 2000/1000 bp
## filter for ensembl canonical ENS###
## use genecode, transmap ensembl, and ensembl gene  

## grik2 species coordinates

# human - GRCH38
# human_ncbi_NC_000006.12_2898 <- "chr6:101393708-102070083"
# human_ensembl_ENSG00000164418 <- "chr6:100962701-102081622"
# human_ensembl_canonical_ENST00000369134.9 <- "chr6:101393708-102070083"

# chimp - Pantro3.0
# chimp_ncbi_NC_072403.2 <- "chr5:112823618-114006919"
# chimp_ensembl_ENSPTRG00000018449 <- "chr6:104757582-105448610"
# chimp_ensembl_canonical_ENSPTRT00000034097.3 <- "chr6:104757582-105448610"

# macaque - Mmul_10
# macaque_ncbi_NC_041757.1 <- "chr4:71196062-71892145"
# macaque_ensembl_ENSMMUG00000009856 <- "chr4:71197621-71886540"
# macaque_ensembl_canonical_ENSMMUT00000072242.2 <- "chr4:71197621-71886540"

# mouse - GRCm39
# mouse_ncbi_NC_000076.7 <- "chr10:48969776-49666523"
# mouse_ensembl_ENSMUSG00000056073 <- "chr10:48970929-49664862"
# mouse_ensembl_canonical_ENSMUST00000218823.2 <- "chr10:48970929-49664862"

## human

grep 'ENST00000369134.9' grik2_promoter2000_human.fa | sed 's/^>//' > ENST00000369134.9.lst
seqtk subseq grik2_promoter2000_human.fa ENST00000369134.9.lst > grik2_promoter2000_human_canonical.fa
grep 'ENST00000369134.9' grik2_promoter1000_human.fa | sed 's/^>//' > ENST00000369134.9.lst
seqtk subseq grik2_promoter1000_human.fa ENST00000369134.9.lst > grik2_promoter1000_human_canonical.fa

## chimp - Pantro6

grep 'ENSPTRT00000034097.3' grik2_promoter2000_chimp.fa

## macaque

grep 'ENSMMUT00000072242.2' grik2_promoter2000_macaque.fa | sed 's/^>//' > ENSMMUT00000072242.2.lst
seqtk subseq grik2_promoter2000_macaque.fa ENSMMUT00000072242.2.lst > grik2_promoter2000_macaque_canonical.fa
grep 'ENSMMUT00000072242.2' grik2_promoter1000_macaque.fa | sed 's/^>//' > ENSMMUT00000072242.2.lst
seqtk subseq grik2_promoter1000_macaque.fa ENSMMUT00000072242.2.lst > grik2_promoter1000_macaque_canonical.fa

## mouse

grep 'ENSMUST00000218823.2' grik2_promoter2000_mouse.fa | sed 's/^>//' > ENSMUST00000218823.2.lst
seqtk subseq grik2_promoter2000_mouse.fa ENSMUST00000218823.2.lst > grik2_promoter2000_mouse_canonical.fa
grep 'ENSMUST00000218823.2' grik2_promoter1000_mouse.fa | sed 's/^>//' > ENSMUST00000218823.2.lst
seqtk subseq grik2_promoter1000_mouse.fa ENSMUST00000218823.2.lst > grik2_promoter1000_mouse_canonical.fa

## cluster MSA via jalview

## ensembl

## NCBI

## R getpromoter and BiomaRt

## check whether promoter seq make sense

## download cactus 447-way MAF for the human promoter coordinates from table browser
## view with mview

mview -in maf447_promoter1000_human.maf -html head -bold -css on -coloring any maf447_promoter1000_human.maf > maf447_promoter1000_human.html

## mview cannot process by species, use maf_parse from phast  

maf_parse maf447_promoter1000_human.maf --order hg38.chr6,Macaca_mulatta.NC_041757.1,Mus_musculus.chr10 > maffiltered_promoter1000_human.maf

mview -in maffiltered_promoter1000_human.maf -html head -bold -css on -coloring any maffiltered_promoter1000_human.maf > maffiltered_promoter1000_human.html
mview -in maffiltered_promoter1000_human.maf -html head -bold -css on -coloring mismatch maffiltered_promoter1000_human.maf > maffiltered_promoter1000_human.html

#### 1.2 5' 3' UTRs

## USCS table browser

## NOT for grik2 coordinates from ensembl
## BUT for grik2 coordinates from ensembl canonical or ncbi database
## choose seqeunce then 5' 3' UTRs 
## filter for ensembl canonical ENS###

## human

grep 'ENST00000369134.9' grik2_5UTR_human.fa | sed 's/^>//' > ENST00000369134.9.lst
seqtk subseq grik2_5UTR_human.fa ENST00000369134.9.lst > grik2_5UTR_human_canonical.fa
grep 'ENST00000369134.9' grik2_3UTR_human.fa | sed 's/^>//' > ENST00000369134.9.lst
seqtk subseq grik2_3UTR_human.fa ENST00000369134.9.lst > grik2_3UTR_human_canonical.fa

## mouse

grep 'ENSMUST00000218823.2' grik2_5UTR_mouse.fa | sed 's/^>//' > ENSMUST00000218823.2.lst
seqtk subseq grik2_5UTR_mouse.fa ENSMUST00000218823.2.lst > grik2_5UTR_mouse_canonical.fa
grep 'ENSMUST00000218823.2' grik2_3UTR_mouse.fa | sed 's/^>//' > ENSMUST00000218823.2.lst
seqtk subseq grik2_3UTR_mouse.fa ENSMUST00000218823.2.lst > grik2_3UTR_mouse_canonical.fa

## cluster MSA via jalview

#### 1.3 order promoter plasmids from idt

## merge promoter seq for human, macaque, and mouse
cat grik2_promoter1000_human_canonical.fa grik2_promoter1000_macaque_canonical.fa grik2_promoter1000_mouse_canonical.fa > promoter1000_hmm_canonical.txt
## check if restriction sites of MluI (ACGCGT) and AgeI (ACCGGT) are inside the target seq
grep -e 'ACGCGT' -e 'ACCGGT' promoter1000_hmm_canonical.txt 
## add restriction sites with 4-5 random bps (e.g., TAGCG)
awk '/^>/ {print $0; next} {print "TAGCGACGCGT"$0"ACCGGTTAGCG"}' promoter1000_hmm_canonical.txt > promoter1000_hmm_restsites.txt

#### 2. Transcription factor binding site analysis for promoter

#### 2.1 The meme suite

## FIMO

## downlaod JASPAR2024 single batch file (txt) with meme format from https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

## human

fimo -oc fimo_out_1 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt grik2_promoter1000_human_canonical.fa 
fimo -oc fimo_out_2 JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt grik2_promoter2000_human_canonical.fa 

centrimo -verbosity 1 -oc centrimo_out_1 -score 5.0 -ethresh 10.0 grik2_promoter1000_human_canonical.fa JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
centrimo -verbosity 1 -oc centrimo_out_2 -score 5.0 -ethresh 10.0 grik2_promoter2000_human_canonical.fa JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

## macaque and mouse

fimo -oc fimo_mac JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt grik2_promoter1000_macaque_canonical.fa
fimo -oc fimo_mouse JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt grik2_promoter1000_mouse_canonical.fa

## code to count and get uniq - may not work!

# Function to count lines after the 99 percentile of q-value column
count_lines_after_99_percentile() {
  local file=$1
  PVAL_99=$(awk -F '\t' '{print $3}' $file | sort -n | awk 'NR < (0.01 * NR) {print; exit}')
  awk -F '\t' -v pval="$PVAL_99" '$3 < pval' $file | wc -l
}

# Count the number of lines in each file
echo "Line counts:"
wc -l fimo_mac/fimo.tsv
wc -l fimo_mouse/fimo.tsv
wc -l fimo_out_1/fimo.tsv

# Count the number of lines after the 99 percentile of column q-value
echo "Lines after 99 percentile:"
echo "fimo_mac:"
count_lines_after_99_percentile fimo_mac/fimo.tsv
echo "fimo_mouse:"
count_lines_after_99_percentile fimo_mouse/fimo.tsv
echo "fimo_out_1:"
count_lines_after_99_percentile fimo_out_1/fimo.tsv

# Extract the second column and find unique entries compared to the other files
awk -F '\t' '{print $2}' fimo_mac/fimo.tsv | sort | uniq > fimo_mac_col2.txt
awk -F '\t' '{print $2}' fimo_mouse/fimo.tsv | sort | uniq > fimo_mouse_col2.txt
awk -F '\t' '{print $2}' fimo_out_1/fimo.tsv | sort | uniq > fimo_out_1_col2.txt

# Find unique values in each file compared to the others
comm -23 fimo_mac_col2.txt <(sort fimo_mouse_col2.txt fimo_out_1_col2.txt | uniq) > unique_to_fimo_mac.txt
comm -23 fimo_mouse_col2.txt <(sort fimo_mac_col2.txt fimo_out_1_col2.txt | uniq) > unique_to_fimo_mouse.txt
comm -23 fimo_out_1_col2.txt <(sort fimo_mac_col2.txt fimo_mouse_col2.txt | uniq) > unique_to_fimo_out_1.txt

echo "Unique values in fimo_mac compared to the others:"
cat unique_to_fimo_mac.txt
echo "Unique values in fimo_mouse compared to the others:"
cat unique_to_fimo_mouse.txt
echo "Unique values in fimo_out_1 compared to the others:"
cat unique_to_fimo_out_1.txt

#### 2.2 Cluster-Buster 

## filter JASPAR2024 based on brain expression for human and macaque
## human
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
grep "^MOTIF" JASPAR2024_synaptogenesis_meme.txt | awk '{print$3}' > JASPAR2024_synaptogenesis_meme_genesym.txt
grep -A 4 -e "Arnt" -e "PAX6" -e "RORA" -e "RELA" -e "NR1H2::RXRA" -e "Znf423" -e "NFIC::TLX1" -e "ZNF354C" -e "Pou5f1::Sox2" -e "EWSR1-FLI1" -e "Arid3a" -e "INSM1" -e "RARA::RXRA" -e "Tfcp2l1" -e "SREBF1" -e "SREBF2" -e "Atf1" -e "Creb3l2" -e "FOXG1" -e "Foxj2" -e "BHLHE41" -e "HNF1B" -e "HSF1" -e "JDP2" -e "MLX" -e "NEUROG2" -e "Nr2e1" -e "OLIG2" -e "SRF" -e "TFAP4" -e "ZIC1" -e "RARA" -e "BCL6B" -e "KLF16" -e "ZNF143" -e "HSF2" -e "MEF2D" -e "MEIS2" -e "NFKB1" -e "PKNOX2" -e "POU3F2" -e "POU3F3" -e "POU3F4" -e "PROX1" -e "SMAD3" -e "TGIF2" -e "MGA" -e "TFAP2B" -e "TFAP2C" -e "HEY1" -e "OLIG1" -e "OLIG3" -e "CEBPG" -e "FOXC1" -e "FOXO4" -e "FOXO6" -e "RXRG" -e "Rarb" -e "TP53" -e "GMEB2" -e "MTF1" -e "E2F1" -e "SOX21" -e "Sox1" -e "TFAP2A" -e "ZNF24" -e "FOSB::JUN" -e "FOSL1::JUN" -e "FOSB::JUNB" -e "TCF7L1" -e "FERD3L" -e "FOXN3" -e "HES6" -e "NR2F1" -e "RXRG" -e "SMAD5" -e "TFAP4" -e "TGIF2LX" -e "TGIF2LY" -e "ZNF135" -e "ZNF136" -e "ZNF382" -e "ZNF460" -e "ZNF528" -e "ZNF707" -e "ZNF8" -e "POU2F1::SOX2" -e "ZFP14" -e "ZNF418" -e "SP1" -e "RFX1" -e "SP2" -e "BHLHA15" -e "SP4" -e "KLF12" -e "RFX3" -e "BHLHE22" -e "TFE3" -e "KLF10" -e "NR1D2" -e "SOX18" -e "Thap11" -e "ZNF75D" -e "BACH1" -e "ARNT::HIF1A" -e "Ddit3::Cebpa" -e "Foxq1" -e "ZNF524" -e "ZNF75A" -e "ZNF766" -e "ZNF770" -e "ZSCAN16" -e "ZBTB17" -e "FOXP4" -e "FOXS1" -e "ZNF184" -e "ZNF213" -e "BCL11A" -e "EPAS1" -e "IKZF2" -e "ZBTB11" -e "ZBTB24" -e "ZNF157" -e "ZNF175" -e "ZNF35" -e "ZNF558" -e "ZSCAN21" -e "FEZF2" -e "ALX3" -e "ARNT2" -e "ASCL1" -e "ATF2" -e "ATF3" -e "ATF4" -e "ATF6" -e "ATF7" -e "MAFK" -e "BACH2" -e "BARX1" -e "BARX2" -e "BATF" -e "JUN" -e "BATF3" -e "BCL6" -e "BHLHE22" -e "BHLHE40" -e "CEBPA" -e "CEBPB" -e "CEBPD" -e "CEBPE" -e "CLOCK" -e "CREB1" -e "CREB3" -e "CREB3L1" -e "CREB3L4" -e "CREM" -e "CTCF" -e "CUX1" -e "CUX2" -e "DBP" -e "DLX1" -e "E2F3" -e "E2F4" -e "EBF3" -e "EGR1" -e "EGR3" -e "EGR4" -e "ELF1" -e "ELF2" -e "ELF3" -e "ELK1::SREBF2" -e "TEF" -e "ELK3" -e "ELK4" -e "EMX1" -e "EMX2" -e "ERF" -e "ERF::FIGLA" -e "FIGLA" -e "ERF::FOXI1" -e "ERF::FOXO1" -e "ERF::HOXB13" -e "MAX" -e "ERF::NHLH1" -e "NHLH1" -e "ESRRA" -e "ESRRB" -e "ETS1" -e "ETS2" -e "ETV1" -e "ETV2::FIGLA" -e "ETV3" -e "ETV4" -e "ETV5" -e "ETV5::DRGX" -e "ETV5::FIGLA" -e "ETV5::FOXI1" -e "ETV5::FOXO1" -e "ETV5::HOXA2" -e "ETV6" -e "FOS" -e "JUNB" -e "JUND" -e "FOSB::JUNB" -e "FOSL1" -e "FOSL2" -e "FOXF2" -e "FOXJ2::ELF1" -e "FOXK1" -e "FOXK2" -e "FOXO1::ELF1" -e "FOXO1::ELK1" -e "FOXO1::ELK3" -e "FOXO1::FLI1" -e "FOXP1" -e "FOXP2" -e "SOX2" -e "TCF3" -e "HES1" -e "HIC2" -e "HIF1A" -e "HINFP" -e "HLF" -e "HMBOX1" -e "HNF1A" -e "IKZF1" -e "IRF2" -e "IRF3" -e "JDP2" -e "JUN::JUNB" -e "KLF11" -e "KLF15" -e "KLF3" -e "KLF4" -e "KLF6" -e "KLF7" -e "KLF9" -e "LHX2" -e "LHX6" -e "LMX1A" -e "LMX1B" -e "MAF" -e "MAF::NFE2" -e "MAFA" -e "MAFF" -e "MAFG::NFE2L1" -e "MAX::MYC" -e "MYC" -e "MAZ" -e "MEF2A" -e "MEF2C" -e "MEIS2" -e "MEIS3" -e "MGA::EVX1" -e "MITF" -e "MLXIPL" -e "MNT" -e "MSANTD3" -e "MSX1" -e "MXI1" -e "MYCN" -e "MYOD1" -e "MZF1" -e "NEUROD1" -e "NEUROG2" -e "NFATC3" -e "NFATC4" -e "NFIA" -e "NFIB" -e "NFIC" -e "NFIX" -e "NFYA" -e "NFYB" -e "NFYC" -e "NKX6-2" -e "NR1D1" -e "NR2C1" -e "NR2C2" -e "NR2F1" -e "NR2F2" -e "NR3C1" -e "NR3C2" -e "NR4A1" -e "NR4A2" -e "NR4A2::RXRA" -e "NRL" -e "OTX1" -e "PATZ1" -e "PBX1" -e "PBX3" -e "PKNOX1" -e "PLAGL2" -e "POU2F1" -e "POU2F2" -e "POU2F3" -e "POU3F1" -e "POU4F3" -e "POU6F1" -e "PPARD" -e "PPARG" -e "RARA::RXRG" -e "RARB" -e "RBPJ" -e "RELB" -e "RFX4" -e "RFX5" -e "RFX7" -e "RORA" -e "RORB" -e "SATB1" -e "SCRT1" -e "SCRT2" -e "SMAD2" -e "SNAI3" -e "SOX10" -e "SOX12" -e "SOX13" -e "SOX14" -e "SOX15" -e "SOX4" -e "SOX8" -e "SOX9" -e "SP3" -e "SP9" -e "SREBF1" -e "STAT1" -e "STAT1::STAT2" -e "STAT3" -e "TAL1::TCF3" -e "TBR1" -e "TCF12" -e "TCF4" -e "TCF7L2" -e "TCFL5" -e "TEAD1" -e "TEAD2" -e "TFAP2A" -e "TFAP2B" -e "TFAP2C" -e "TFAP2E" -e "TFAP4::ETV1" -e "TFAP4::FLI1" -e "TFCP2" -e "TFDP1" -e "TFEB" -e "TFEC" -e "THAP1" -e "THRA" -e "THRB" -e "TWIST1" -e "USF1" -e "USF2" -e "VEZF1" -e "XBP1" -e "ZBTB12" -e "ZBTB14" -e "ZBTB18" -e "ZBTB26" -e "ZBTB33" -e "ZBTB7A" -e "ZBTB7B" -e "ZEB1" -e "ZKSCAN1" -e "ZKSCAN3" -e "ZKSCAN5" -e "ZNF140" -e "ZNF148" -e "ZNF16" -e "ZNF189" -e "ZNF211" -e "ZNF214" -e "ZNF257" -e "ZNF263" -e "ZNF274" -e "ZNF281" -e "ZNF282" -e "ZNF317" -e "ZNF320" -e "ZNF324" -e "ZNF331" -e "ZNF341" -e "ZNF343" -e "ZNF354A" -e "ZNF384" -e "ZNF410" -e "ZNF416" -e "ZNF417" -e "ZNF449" -e "ZNF454" -e "ZNF549" -e "ZNF574" -e "ZNF610" -e "ZNF652" -e "ZNF667" -e "ZNF675" -e "ZNF680" -e "ZNF682" -e "ZNF692" -e "ZNF701" -e "ZNF708" -e "ZNF740" -e "ZNF76" -e "ZNF768" -e "ZNF784" -e "ZNF816" -e "ZNF85" -e "ZSCAN29" -e "ZSCAN31" -e "Ebf4" -e "Hmga1" -e "Nr1h3" -e "Nanog" -e "Ahr::Arnt" -e "Arid3b" -e "Arntl" -e "Arx" -e "Atf3" -e "Bach1::Mafk" -e "Jun" -e "Bcl11B" -e "Bhlha15" -e "Dlx2" -e "Dlx5" -e "Foxo1" -e "Esrrg" -e "Foxj3" -e "Gmeb1" -e "Hand1::Tcf3" -e "Hmx1" -e "Hnf1A" -e "Ikzf3" -e "Lef1" -e "Mafb" -e "Mafg" -e "Mlxip" -e "Neurod2" -e "Nfat5" -e "Nfe2l2" -e "Npas2" -e "Nr1H2" -e "Nr1h3::Rxra" -e "Nrf1" -e "Olig2" -e "Plagl1" -e "Pparg::Rxra" -e "Prdm15" -e "Prdm4" -e "Ptf1A" -e "Smad4" -e "Sox11" -e "Sox17" -e "Sox5" -e "Spi1" -e "Stat2" -e "Stat4" -e "Stat5a::Stat5b" -e "Stat5b" -e "Stat6" -e "Tcf12" -e "Yy1" -e "Zfp335" -e "Zfx" -e "Zic1::Zic2" -e "Zic2" -e "Zic3" -e "LIN54" -e "TBP" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt > JASPAR2024_synaptogenesis_genes_human.txt

grep ">" JASPAR2024_synaptogenesis_genes_human.txt | wc -l
grep ">" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt | wc -l

#filtered 879 to 510 

sed '/^--$/d' JASPAR2024_synaptogenesis_genes_human.txt > JASPAR2024_synaptogenesis_genes_human_mod.txt
awk '/^>/ {print; next} {gsub(/^[ACGT] /,""); gsub(/\[|\]/,""); print}' JASPAR2024_synaptogenesis_genes_human_mod.txt > JASPAR2024_pfms_synaptogenesis_cbust.txt
awk '/^>/ {print; next} {gsub(/^[ACGT] /,""); gsub(/\[|\]/,""); print}' JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt > JASPAR2024_pfms_cbust.txt

## run cluster-buster

./cluster-buster/cbust -f 1 -c0 -m0 JASPAR2024_pfms_synaptogenesis_cbust.txt grik2_promoter1000_human_canonical.fa > cbustout_humanpro_1.txt
./cluster-buster/cbust -f 1 -c0 -m0 -l JASPAR2024_pfms_synaptogenesis_cbust.txt grik2_promoter1000_human_canonical.fa  > cbustout_humanpro_2.txt
./cluster-buster/cbust -f 1 -c0 -m0 -l JASPAR2024_pfms_synaptogenesis_cbust.txt grik2_promoter1000_human_canonical.fa  > cbustout_humanpro_3.txt

## other flags: -g3 -c5 - m3 , -c0.01 -g20 , -c0 -m0, -c0 -m0 -l

## macaque

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
grep "^MOTIF" JASPAR2024_synaptogenesis_meme_macaque.txt | awk '{print$3}' > JASPAR2024_synaptogenesis_meme_genesym_macaque.txt
grep -A 4 -e "Arnt" -e "PAX6" -e "RORA" -e "RELA" -e "NR1H2::RXRA" -e "Znf423" -e "NFIC::TLX1" -e "ZNF354C" -e "Pou5f1::Sox2" -e "EWSR1-FLI1" -e "Arid3a" -e "INSM1" -e "RARA::RXRA" -e "Tfcp2l1" -e "SREBF1" -e "SREBF2" -e "KLF5" -e "Atf1" -e "Creb3l2" -e "FOXG1" -e "Foxj2" -e "BHLHE41" -e "HNF1B" -e "HSF1" -e "JDP2" -e "MLX" -e "Nr2e1" -e "OLIG2" -e "SRF" -e "TFAP4" -e "ZIC1" -e "RARA" -e "BCL6B" -e "EGR2" -e "KLF16" -e "ZNF143" -e "HSF2" -e "MEF2D" -e "MEIS2" -e "NFKB1" -e "PKNOX2" -e "POU3F2" -e "POU3F3" -e "PROX1" -e "SMAD3" -e "MGA" -e "TFAP2B" -e "TFAP2C" -e "HEY1" -e "OLIG1" -e "OLIG3" -e "CEBPG" -e "FOXC1" -e "FOXO4" -e "FOXO6" -e "FOXP3" -e "RXRG" -e "Rarb" -e "TP53" -e "GMEB2" -e "MTF1" -e "E2F1" -e "SOX21" -e "Sox1" -e "TFAP2A" -e "ZNF24" -e "FOSB::JUN" -e "FOSL1::JUN" -e "FOSB::JUNB" -e "IRF5" -e "FERD3L" -e "FOXN3" -e "HES6" -e "NR2F1" -e "SMAD5" -e "TFAP4" -e "ZNF135" -e "ZNF136" -e "ZNF382" -e "ZNF528" -e "ZNF707" -e "ZNF8" -e "POU2F1::SOX2" -e "ZFP14" -e "ZNF418" -e "SP1" -e "RFX1" -e "SP2" -e "SP4" -e "KLF12" -e "RFX3" -e "BHLHE22" -e "TFE3" -e "KLF10" -e "NR1D2" -e "SOX18" -e "Thap11" -e "ZNF75D" -e "BACH1" -e "ARNT::HIF1A" -e "Ddit3::Cebpa" -e "Foxq1" -e "ZNF524" -e "ZNF75A" -e "ZNF766" -e "ZNF770" -e "ZSCAN16" -e "ZBTB17" -e "FOXP4" -e "FOXS1" -e "ZNF184" -e "ZNF213" -e "BCL11A" -e "EPAS1" -e "IKZF2" -e "ZBTB11" -e "ZBTB24" -e "ZNF157" -e "ZNF175" -e "ZNF35" -e "ZNF547" -e "ZNF558" -e "ZSCAN21" -e "FEZF2" -e "ARNT2" -e "ASCL1" -e "ATF2" -e "ATF3" -e "ATF4" -e "ATF6" -e "ATF7" -e "MAFK" -e "BACH2" -e "BARX1" -e "BARX2" -e "BATF" -e "JUN" -e "BATF3" -e "BCL6" -e "BHLHE22" -e "BHLHE40" -e "CEBPA" -e "CEBPB" -e "CEBPD" -e "CEBPE" -e "CEBPG" -e "CLOCK" -e "CREB1" -e "CREB3" -e "CREB3L1" -e "CREB3L4" -e "CREM" -e "CTCF" -e "CTCFL" -e "CUX1" -e "CUX2" -e "DBP" -e "DLX1" -e "E2F3" -e "E2F4" -e "EBF3" -e "EGR1" -e "EGR3" -e "EGR4" -e "ELF1" -e "ELF2" -e "ELF3" -e "SPDEF" -e "ELK1::SREBF2" -e "TEF" -e "ELK3" -e "ELK4" -e "EMX1" -e "EMX2" -e "ERF" -e "MAX" -e "ESRRA" -e "ESRRB" -e "ETS1" -e "ETS2" -e "ETV1" -e "ETV3" -e "ETV4" -e "ETV5" -e "ETV6" -e "FOS" -e "FOS::JUN" -e "JUNB" -e "FOSL1" -e "FOXF2" -e "FOXJ2::ELF1" -e "FOXK1" -e "FOXK2" -e "FOXP1" -e "FOXP2" -e "SOX2" -e "TCF3" -e "HES1" -e "HEY2" -e "HIC2" -e "HIF1A" -e "HINFP" -e "HLF" -e "HMBOX1" -e "HNF1A" -e "IKZF1" -e "IRF2" -e "IRF3" -e "JUN::JUNB" -e "KLF11" -e "KLF15" -e "KLF3" -e "KLF4" -e "KLF6" -e "KLF7" -e "KLF9" -e "LHX2" -e "LHX6" -e "LMX1A" -e "LMX1B" -e "MAF" -e "MAFA" -e "MAFF" -e "MYC" -e "MAZ" -e "MEF2A" -e "MEF2C" -e "MEIS3" -e "MITF" -e "MLXIPL" -e "MNT" -e "MXI1" -e "MYCN" -e "MYOD1" -e "MZF1" -e "NEUROD1" -e "NFATC3" -e "NFIA" -e "NFIB" -e "NFIC" -e "NFIX" -e "NFKB2" -e "NFYA" -e "NKX6-2" -e "NR1D1" -e "NR2C1" -e "NR2F2" -e "NR3C1" -e "OTX1" -e "OVOL2" -e "PATZ1" -e "PBX1" -e "POU2F1" -e "PPARD" -e "PPARG" -e "RBPJ" -e "RELB" -e "RFX4" -e "RFX7" -e "RORB" -e "SATB1" -e "SOX10" -e "SOX12" -e "SOX4" -e "SOX9" -e "SP3" -e "SP9" -e "STAT1" -e "STAT3" -e "TBR1" -e "TCF4" -e "TEAD1" -e "TFEB" -e "THRA" -e "THRB" -e "TRPS1" -e "USF1" -e "USF2" -e "XBP1" -e "ZBTB7A" -e "ZBTB7B" -e "ZBTB33" -e "ZEB1" -e "ZKSCAN1" -e "ZKSCAN3" -e "ZNF140" -e "ZNF148" -e "ZNF189" -e "ZNF211" -e "ZNF263" -e "ZNF281" -e "ZNF282" -e "ZNF317" -e "ZNF320" -e "ZNF331" -e "ZNF341" -e "ZNF343" -e "ZNF354A" -e "ZNF384" -e "ZNF410" -e "ZNF416" -e "ZNF454" -e "ZNF549" -e "ZNF574" -e "ZNF582" -e "ZNF652" -e "ZNF667" -e "ZNF675" -e "ZNF701" -e "ZNF768" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt > JASPAR2024_synaptogenesis_genes_macaque.txt
grep ">" JASPAR2024_synaptogenesis_genes_macaque.txt | wc -l
grep ">" JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt | wc -l

#filtered 879 to 384

sed '/^--$/d' JASPAR2024_synaptogenesis_genes_macaque.txt > JASPAR2024_synaptogenesis_genes_macaque_mod.txt
awk '/^>/ {print; next} {gsub(/^[ACGT] /,""); gsub(/\[|\]/,""); print}' JASPAR2024_synaptogenesis_genes_macaque_mod.txt > JASPAR2024_pfms_synaptogenesis_macaque_cbust.txt
awk '/^>/ {print; next} {gsub(/^[ACGT] /,""); gsub(/\[|\]/,""); print}' JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt > JASPAR2024_pfms_cbust.txt

## run cluster-buster

./cluster-buster/cbust -f 1 -c0 -m0 -l JASPAR2024_pfms_synaptogenesis_macaque_cbust.txt grik2_promoter1000_macaque_canonical.fa > cbustout_macaquepro_1.txt
./cluster-buster/cbust -f 2 -c0 -m0 -l JASPAR2024_pfms_synaptogenesis_macaque_cbust.txt grik2_promoter1000_macaque_canonical.fa  > cbustout_macaquepro_2.txt
./cluster-buster/cbust -f 3 -c0 -m0 -l JASPAR2024_pfms_synaptogenesis_macaque_cbust.txt grik2_promoter1000_macaque_canonical.fa  > cbustout_macaquepro_3.txt

## get right gene names

grep "^MOTIF" JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt | awk '{print $2,$3}' > gene_mapping.txt

gene_mapping_file="gene_mapping.txt"
input_file="cbustout_humanpro_1.txt"
output_file="cbustout_humanpro_1_gene_symbols.txt"
temp_file=$(mktemp)
awk '{gsub(/\/|&/, "\\&"); print "s/" $1 "/" $2 "/g"}' "$gene_mapping_file" > sed_commands.txt
sed -f sed_commands.txt "$input_file" > "$temp_file"
mv "$temp_file" "$output_file"
rm sed_commands.txt

gene_mapping_file="gene_mapping.txt"
input_file="cbustout_macaquepro_1.txt"
output_file="cbustout_macaquepro_1_gene_symbols.txt"
temp_file=$(mktemp)
awk '{gsub(/\/|&/, "\\&"); print "s/" $1 "/" $2 "/g"}' "$gene_mapping_file" > sed_commands.txt
sed -f sed_commands.txt "$input_file" > "$temp_file"
mv "$temp_file" "$output_file"
rm sed_commands.txt

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



