#### The Roles of GRIK2 and HARs in the Evolution and Development of Synapses
#### Shadi Shahatit, Master's Thesis, IBENS, Paris 2024

#### Content:
#### 1. PsychENCODE data
#### 2. Synaptogenesis cortical TFs in humans
#### 3. WGCNA
#### Appendix I - Tools and scripts

#### 1. PsychENCODE data

## Replot Fig.1d (Liu et al., 2024)
## download data from http://evolution.psychencode.org/#

(head -1 nhp_development_RPKM_rmTechRep.txt && grep 'ENSG00000164418' nhp_development_RPKM_rmTechRep.txt) > grik2_nhp_development_RPKM_rmTechRep.txt

sed "1s/^/gene_name\t/" nhp_development_RPKM_rmTechRep.txt > nhp_development_RPKM_rmTechRep_genename.txt

#### 2. Synaptogenesis cortical TFs in humans

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

#### 3. WGCNA

## filter encode data based on human, tissue, and age

awk 'BEGIN { FS=OFS="\t" } { for (i=1; i<=NF; i++) if ($i !~ /RMB/) output = (output ? output OFS : "") $i; print output; output="" }' nhp_development_RPKM_rmTechRep_genename.txt > human_encode_data.txt

awk -F"\t" 'NR==1 {
    # Initialize arrays for filtering
    for (i=1; i<=NF; i++) {
        if ($i ~ /HSB119|HSB127|HSB130|HSB136|HSB126|HSB145|HSB123|HSB135/) samples[i]=1;
        if ($i ~ /CGE|DTH|LGE|MGE|OC|PC|TC|URL|MSC|VFC_001|STC_GAIIx|A1C_GAIIx|M1C_GAIIx/) tissues[i]=1;
        if (samples[i] || !tissues[i]) keep[i]=1;
    }
    # Print header with filtered columns
    for (i=1; i<=NF; i++) {
        if (i == 1 || keep[i]) printf "%s\t", $i;
    }
    print "";
}
NR>1 {
    # Print data rows with filtered columns, keeping the first column
    printf "%s\t", $1;
    for (i=2; i<=NF; i++) {
        if (keep[i]) printf "%s\t", $i;
    }
    print "";
}' human_encode_data.txt > human_encode_data_filtered.txt

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



