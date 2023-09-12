# Irenaeus Chan
# May 1st, 2022

library(dplyr)
library(stringr)
library(jsonlite)
library(sqldf)

bolton_bick_vars <- "/Users/irenaeuschan/Documents/Irenaeus/data/bick.bolton.vars3.txt"
cosmic_hotspots <- "/Users/irenaeuschan/Documents/Irenaeus/data/COSMIC.heme.myeloid.hotspot.tsv"
gene_list <- "/Users/irenaeuschan/Documents/Irenaeus/UKBB/oncoKB_CGC_pd_table_disparity_KB.csv"
pd_table <- "/Volumes/bolton/Active/Protected/Annotation_Files/pd_table_kbreview_bick_trunc4_oncoKB_SAFE.filtered_genes_oncoKB_CGC.tsv"

#Pilot <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/pilot.archer.combined.FPpass.tsv"
#Pilot <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/pilot.mgi.combined.FPpass.tsv"
Pilot <- "/Volumes/bolton/Active/Projects/ArcherPilot/Pilot/TERRA/pilot.combined.tsv"
#Trios <- "/Users/irenaeuschan/Documents/Irenaeus/MGI_Yizhe/trios.combined.tsv"
#MGI_EXTRA <- "/Volumes/bolton/Active/Projects/MGI_Data/TERRA_EXTRA/MGI_EXTRA.combined.tsv"
#Dilution <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/dilution.combined.tsv"
#Final <- "/Volumes/bolton/Active/Projects/ProstateCancer/TERRA/prostate.final.FPpass.tsv"
#Prostate <- "/Users/irenaeuschan/Documents/Irenaeus/ProstateCancer/prostate.final.combined.FPpass.filtered_KB2.csv"
Final <- "/Volumes/bolton/Active/Projects/GoodCell/TERRA/GoodCell2.combined.FPpass.tsv"
#Final <- "/Users/irenaeuschan/Documents/Irenaeus/Freidman/final.combined.Freidman.FPpass.tsv"
#Final <- "/Users/irenaeuschan/Documents/Irenaeus/ArcherDX/final.combined.FPpass.tsv"
#Orig <- "/Users/irenaeuschan/Documents/Irenaeus/ArcherDX/data/variant_review_IC_31722_KB_complete_updated.csv"
#Alex_Filter <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/alex_filter.csv"
# alex_filter <- read.csv(Alex_Filter, header = TRUE)
#Final <- "/Volumes/bolton/Active/Projects/ElementBio/TERRAFinal/final.combined.tsv"

unfiltered <- read.table(Final, sep='\t', header = TRUE, comment.char = '', quote = '')
final <- unfiltered

#final <- left_join(final, harvard_duplicates, by=c("SN_TAG"="fakeid")) %>% 
#  mutate(id_export = ifelse(is.na(id_export), SN_TAG, id_export)) %>%
#  mutate(SN_TAG = id_export) %>% select(-id_export)

# Remove all variants that fail our general filter "all_fp_pass_XGB"
# These filters are: PoN_pvalue, PoN_Zscore, lowAF, long_indels, bcbio, di_tri_nuc
# This does NOT include fp_filter, that needs to be handled separately
final <- final %>% dplyr::filter(as.logical(all_fp_pass_XGB))

# Remove any weird Lofreq N calls
final <- final %>% dplyr::filter(ifelse(grepl('N', ALT), FALSE, TRUE))
final <- final %>% dplyr::filter(ifelse(grepl('N', REF), FALSE, TRUE))

# Remove Off-Target & Introns          
final <- final %>% dplyr::filter(Consequence_VEP != "")

# Remove Variants that Failed PoN
# Bonferroni Corrected p-value for Fisher's Exact Test (0.05 / Size of BED Panel)
pon_pvalue = 2.114164905e-6
pon_pvalue = 0.000000511148141
pon_pvalue = 8.85582713e-7
final <- final %>% mutate(pon_pvalue_Mutect_Raw = ifelse(is.na(pon_pvalue_Mutect_Raw), 0, pon_pvalue_Mutect_Raw),
                 pon_pvalue_Lofreq_Raw = ifelse(is.na(pon_pvalue_Lofreq_Raw), 0, pon_pvalue_Lofreq_Raw),
                 pon_pvalue_Vardict_Raw = ifelse(is.na(pon_pvalue_Vardict_Raw), 0, pon_pvalue_Vardict_Raw))
final <- final %>% mutate(pon_FP_pass_XGB = ifelse(
  pon_pvalue_Mutect_Raw <= pon_pvalue & pon_pvalue_Lofreq_Raw <= pon_pvalue & pon_pvalue_Vardict_Raw <= pon_pvalue, TRUE, FALSE
))
final <- final %>% filter(pon_FP_pass_XGB == TRUE)

# Remove Duplicated Rows - If large, this will take a while
# final <- final[!(duplicated(final) | duplicated(final, fromLast = TRUE)), ]

# For the Dilution Data Only
#final <- final %>% mutate(Sample = ifelse(grepl("T", SN_TAG, fixed = TRUE), paste0(unlist(lapply(strsplit(SN_TAG, "T", fixed = TRUE), "[[", 1)), "T"), SN_TAG))
#final <- final %>% mutate(sample_key = paste0(Sample, " ", key))
#dilution_final <- left_join(alex_filter %>% dplyr::select(sample_key, Germline), final, by=c("sample_key"="sample_key"))
#final <- final %>% mutate(SN_TAG = unlist(lapply(strsplit(SN_TAG, "_", fixed = TRUE), "[[", 1)),)

# Dealing with fp_filter
# First filter out which samples failed out other FP filters:
# - PoN 
# - Long INDELs (100 bp)
# - DiNuc and TriNuc in Vardict
# - BCBIO
final <- final %>% filter(as.logical(pon_FP_pass_XGB), 
                          as.logical(long100_indel_pass_XGB), 
                          as.logical(long_indel_pass_XGB), 
                          as.logical(di_tri_vard_pass_XGB),
                          as.logical(bcbio_pass_XGB))

# Then filter out samples that failed Varscan's FP filter
final <- final %>% filter(FP_Filter_DETP20_XGB == 0,
                          FP_Filter_MMQS100_XGB == 0,
                          FP_Filter_MMQSD50_XGB == 0,
                          FP_Filter_NRC_XGB == 0,
                          FP_Filter_PB10_XGB == 0,
                          FP_Filter_RLD25_XGB == 0)

# Calculate average VAF and average Alt Depth
final$average_AF <- final %>% dplyr::select(gt_AF_Mutect, gt_AF_Lofreq, gt_AF_Vardict) %>% rowMeans(., na.rm=TRUE)
final$average_AD <- final %>% dplyr::select(gt_AD_alt_Mutect, gt_AD_alt_Lofreq, gt_AD_alt_Vardict) %>% rowMeans(., na.rm=TRUE)

# Filter out variants lower than 0.001 VAF because hard to tell the difference between artifact and real
final <- final %>% filter(average_AF >= 0.001)

# Filter out samples that failed MCV4 if it had an Alt Depth above 5
final <- final %>% filter(ifelse(FP_Filter_MVC4_XGB == 1 & average_AD > 5, FALSE, TRUE))

# Filter out samples that failed SB1 only if it had an Alt Depth above 10
final <- final %>% filter(ifelse(FP_Filter_SB1_XGB == 1 & average_AD > 10, FALSE, TRUE))

# Only one Caller
final <- final %>%
  mutate(CALL_BY_1 = dplyr::case_when(
    FILTER_Mutect != "" & FILTER_Lofreq == "" & FILTER_Vardict == "" ~ TRUE,
    FILTER_Mutect == "" & FILTER_Lofreq != "" & FILTER_Vardict == "" ~ TRUE,
    FILTER_Mutect == "" & FILTER_Lofreq == "" & FILTER_Vardict != "" ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(CALL_BY_CALLER = dplyr::case_when(
    FILTER_Mutect != "" & FILTER_Lofreq == "" & FILTER_Vardict == "" ~ 'mutect',
    FILTER_Mutect == "" & FILTER_Lofreq != "" & FILTER_Vardict == "" ~ 'lofreq',
    FILTER_Mutect == "" & FILTER_Lofreq == "" & FILTER_Vardict != "" ~ 'vardict',
    FILTER_Mutect != "" & FILTER_Lofreq != "" & FILTER_Vardict == "" ~ 'mutect;lofreq',
    FILTER_Mutect != "" & FILTER_Lofreq == "" & FILTER_Vardict != "" ~ 'mutect;vardict',
    FILTER_Mutect == "" & FILTER_Lofreq != "" & FILTER_Vardict != "" ~ 'lofreq;vardict',
    FILTER_Mutect != "" & FILTER_Lofreq != "" & FILTER_Vardict != "" ~ 'mutect;lofreq;vardict',
  )) 

# Filter our non PD variants
final <- final %>% dplyr::filter(ch_pd2 >= 1)
# Has to be PASS by at least ONE caller
final <- final %>% dplyr::filter(as.logical(PASS_BY_1))
# Filter out Low VAF Mutect Calls
final <- final %>% filter(!(CALL_BY_CALLER == "mutect" & average_AF < 0.01))
# Throw out Silent Mutations
final <- final %>% filter(Consequence_VEP != "synonymous_variant")

# If Median VAF >= 35% and NEVER REPORTED in B/B or COSMIC
# Calculate the Median VAF per Variant
final <- final %>% left_join(., final %>% group_by(key) %>% summarise(median_VAF = median(average_AF)))
final <- final %>% filter(!(median_VAF >= 0.35 & sourcetotalsc_XGB == 0 & CosmicCount == 0))

min_alt <- 5

# Strand Bias Filter - 90/10
final <- final %>% rowwise() %>% mutate(pass_strand_bias = ifelse(sum(pass_strand_bias_Mutect, pass_strand_bias_Lofreq, pass_strand_bias_Vardict, na.rm = TRUE) >= 1, TRUE, FALSE))
# Strand Bias Filter - Minimum Counts
final <- final %>% mutate(enough_strand_evidence_Mutect = ifelse(gt_AD_alt_Mutect >= min_alt & (AltFwd_Mutect_Raw == 0 | AltRev_Mutect_Raw == 0), 0, 1),
                          enough_strand_evidence_Lofreq = ifelse(gt_AD_alt_Lofreq >= min_alt & (AltFwd_Lofreq_Raw == 0 | AltRev_Lofreq_Raw == 0), 0, 1),
                          enough_strand_evidence_Vardict = ifelse(gt_AD_alt_Vardict >= min_alt & (AltFwd_Vardict_Raw == 0 | AltRev_Vardict_Raw == 0), 0, 1))
final <- final %>% rowwise() %>% mutate(pass_strand_evidence = ifelse(sum(enough_strand_evidence_Mutect, enough_strand_evidence_Lofreq, enough_strand_evidence_Vardict, na.rm = TRUE) >= 1, TRUE, FALSE))
final <- final %>% filter(pass_strand_evidence == TRUE)
# Min Alt Count
final <- final %>% mutate(min_alt_count_Mutect = ifelse((AltFwd_Mutect_Raw+AltRev_Mutect_Raw) < min_alt, 0, 1))
final <- final %>% mutate(min_alt_count_Lofreq = ifelse((AltFwd_Lofreq_Raw+AltRev_Lofreq_Raw) < min_alt, 0, 1))
final <- final %>% mutate(min_alt_count_Vardict = ifelse((AltFwd_Vardict_Raw+AltRev_Vardict_Raw) < min_alt, 0, 1))
final <- final %>% rowwise() %>% mutate(pass_min_alt_count = ifelse(sum(min_alt_count_Mutect, min_alt_count_Lofreq, min_alt_count_Vardict, na.rm = TRUE) >= 1, TRUE, FALSE))
final <- final %>% filter(pass_min_alt_count == TRUE) %>% dplyr::select(-min_alt_count_Mutect, -min_alt_count_Lofreq, -min_alt_count_Vardict)

# Redefine homopolymer filters from string to boolean
final <- final %>%
  mutate(case_NXXX = as.logical(case_NXXX),
         case_XNXX = as.logical(case_XNXX),
         case_XXNX = as.logical(case_XXNX),
         case_XXXN = as.logical(case_XXXN),
         case_NNXX = as.logical(case_NNXX),
         case_XNNX = as.logical(case_XNNX),
         case_XXNN = as.logical(case_XXNN)) %>%
  mutate(
    # For SNPs or INDELs with length <=2 use the dust_score_10 over max of 5, 3, & 10
    pass_complexity_filter = dplyr::case_when(
      abs(nchar(REF)-nchar(ALT)) <= 2 ~ !(case_NXXX) &
        !(case_XNXX) &
        !(case_XXNX) &
        !(case_XXXN) &
        !(case_NNXX) &
        !(case_XNNX) &
        !(case_XXNN) &
        dust_score_10 < 7,
      abs(nchar(REF)-nchar(ALT)) > 3 ~ dust_score < 7,
      TRUE ~ !(case_NXXX) &
        !(case_XNXX) &
        !(case_XXNX) &
        !(case_XXXN) &
        !(case_NNXX) &
        !(case_XNNX) &
        !(case_XXNN) &
        dust_score < 7
    ))

# Filter Homopolymer that are not "Exact Match in BB
#rescue <- final %>% filter(!as.logical(pass_homopolymer_filter))
#final <- final %>% filter(!(as.logical(pass_homopolymer_filter) == FALSE & sourcetotalsloci_XGB == 0))
#rescue <- rescue %>% filter(sourcetotalsloci_XGB >= 5 | CosmicCount >= 25 | heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5)
#final <- rbind(final, rescue)

# Recurrent Filter Parameters
# 80 Samples have R882H Mutation in ArcherDX Dataset
count_threshold<- max(ceiling(length(unique((final$SN_TAG)))*0.06), 5)
bb_count_threshold<- max(ceiling(length(unique((final$SN_TAG)))*0.03), 2)
n_samples_min_vaf <- 0.001
n_samples_min_count <- 5

# N_samples
final <- final %>% dplyr::left_join(., final %>%
                                      dplyr::distinct(key, SN_TAG) %>%
                                      group_by(key) %>%
                                      mutate(nsamples = dplyr::n()) %>%
                                      dplyr::ungroup() %>%
                                      dplyr::select(key, nsamples) %>%
                                      dplyr::distinct(), by=c("key"="key"))

final <- final %>% dplyr::left_join(., final %>% 
                                      dplyr::filter(average_AF >= n_samples_min_vaf) %>%
                                      dplyr::distinct(key, SN_TAG) %>%
                                      dplyr::group_by(key) %>%
                                      dplyr::mutate(nsamples_min_vaf = dplyr::n()) %>%
                                      dplyr::ungroup() %>%
                                      dplyr::select(key, nsamples_min_vaf) %>%
                                      dplyr::distinct(), by=c("key"="key")) %>%
  mutate(nsamples_min_vaf = ifelse(is.na(nsamples_min_vaf), 0, nsamples_min_vaf))

final <- final %>% dplyr::left_join(., final %>% 
                                      dplyr::filter(average_AD >= n_samples_min_count) %>%
                                      dplyr::distinct(key, SN_TAG) %>%
                                      dplyr::group_by(key) %>%
                                      dplyr::mutate(nsamples_min_count = dplyr::n()) %>%
                                      dplyr::ungroup() %>%
                                      dplyr::select(key, nsamples_min_count) %>%
                                      dplyr::distinct(), by=c("key"="key")) %>%
  mutate(nsamples_min_count = ifelse(is.na(nsamples_min_count), 0, nsamples_min_count))

# We have to save ASXL1 G646W because this is a very recurrent variant that may be removed from our recurrent filter
final <- final %>% dplyr::filter(nsamples_min_vaf < count_threshold | (key == 'chr20 32434638 A>AG' & average_AF >= 0.05) | (key == 'chr20 32434638 A>AGG' & average_AF >= 0.05))        # Recurrent Filter

# Last filter to remove any variants that have a high recurrent count but are not reported inside Kelly's or Bick's dataset
final <- final %>% dplyr::filter(!(nsamples_min_vaf > bb_count_threshold & (sourcetotalsc_XGB <= 25 & CosmicCount <= 50)))

# Gene Stuff
full_gene_list <- read.csv(gene_list, header = TRUE)
TSG_gene_list <- unname(unlist(full_gene_list %>% filter(isTSG == 1) %>% dplyr::select(Gene)))
gene_list <- full_gene_list$Gene
pd_table <- read.table(pd_table, header = TRUE, sep = "\t")
bick_email <- pd_table %>% filter(source == 'Bick_email')

nonsense_mutation <- c("frameshift_variant", "stop_lost", "stop_gained", "transcript_ablation")
missense_mutation <- c("missense_variant", "inframe_deletion", "inframe_insertion")
SF3B1_positions <- c(622, 623, 624, 625, 626, 662, 663, 664, 665, 666, 700, 701, 702, 703, 704, 740, 741, 742)
clinvar_sig_terms <- c("Likely_pathogenic", "Likely_pathogenic&_drug_response", "Pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic/Likely_pathogenic&_drug_response")

# Germline Filters
# Keep variants that are below our gnomAD VAF filter (<= gnomAD Population VAF of 0.0005)
final <- final %>% dplyr::rename(pass_max_sub_gnomAD_AF = max_gnomAD_AF)
final <- final %>% rowwise() %>% mutate(max_sub_gnomAD_AF = max(max_gnomAD_AF_VEP, max_gnomADe_AF_VEP, max_gnomADg_AF_VEP, na.rm = TRUE)) %>% ungroup()
final$max_pop_gnomAD_AF <- final %>% dplyr::select(gnomAD_AF_VEP, gnomADe_AF_VEP, gnomADg_AF_VEP) %>% apply(., 1, function(x){max(x, na.rm = T)})
final <- final %>% mutate(comp_germline = ifelse(max_pop_gnomAD_AF <= 0.0005 | (key == 'chr20 32434638 A>AG' & average_AF >= 0.05) | (key == 'chr20 32434638 A>AGG' & average_AF >= 0.05), 0, 1))
final <- final %>% mutate(comp_germline = ifelse(average_AF >= 0.35 & (sourcetotalsc_XGB == 0 & sourcetotalsp_XGB == 0 & CosmicCount == 0), 1, comp_germline))
final <- final %>% mutate(comp_germline = ifelse(average_AF >= 0.35 & Gene %in% c("DNMT3A", "TET2", "ASXL1", "PPM1D") & VariantClass %in% nonsense_mutation, 0, comp_germline))
#final <- final %>% mutate(comp_germline = ifelse(average_AF >= 0.25 & average_AF < 0.35 & sourcetotalsc_XGB < 5 & CosmicCount < 25 & max_sub_gnomAD_AF > 0.0001, 1, comp_germline))
final <- final %>% mutate(comp_germline = ifelse(median_VAF >= 0.35 & average_AF >= 0.25, 1, comp_germline))

# Adding INDEL lengths
final <- final %>% mutate(alt_len = nchar(ALT), ref_len = nchar(REF))
# Filter out Complex INDELs
final <- final %>% dplyr::filter(ifelse(CALL_BY_1 == TRUE & (nchar(REF) >= 20 | nchar(ALT) >= 20) & PINDEL_MATCH == "", FALSE, TRUE))

# Adding Near Columns
vars <- read.table(bolton_bick_vars, sep = "\t", header = T, comment.char = "")
vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:","",vars$HGVSc_VEP), sep="_")
vars <- vars %>% dplyr::mutate(key = paste0(CHROM, ':', POS, ' ', REF, '>', ALT))

vars$aa.pos <- as.numeric(str_extract(vars$loci.vep, "\\d+"))
vars$CHROM.POS <- with(vars, paste0(CHROM,":",POS))
vars$GENE.AA.POS <- with(vars, paste(SYMBOL_VEP,aa.pos,sep=":"))
vars$gene_aachange <- with(vars, paste(SYMBOL_VEP,AAchange2,sep=":"))
vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:","",vars$HGVSc_VEP), sep="_")

## new lines to handle DNMT3A:NA, ASXL1:NA etc.
vars$GENE.AA.POS <- with(vars, ifelse(is.na(vars$aa.pos),NA,paste(SYMBOL_VEP,aa.pos,sep=":"))) # new line
vars$GENE.AA.POS[vars$GENE.AA.POS=="NA:NA"] = NA # new line
vars$gene_cDNAchange <-gsub("_",":", vars$gene_cDNAchange)
vars$gene_aachange <- with(vars, paste(SYMBOL_VEP,AAchange2,sep=":"))

final$aa.pos <- as.numeric(str_extract(final$AAchange.x, "\\d+"))
final$gene_aachange <- with(final, paste(SYMBOL_VEP, AAchange.x,sep=":"))

final$near.BB.HS <- apply(final[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange", "gene_cDNAchange")], 1, function(x) {
  p = c(-3:0,1:3)
  n = c(-9:0,1:9)
  
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% vars$GENE.AA.POS
  
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% vars$CHROM.POS
  
  if (any(any.in.p)) {
    res <- unique(vars[vars$GENE.AA.POS %in% vector.p, c("gene_aachange","source.totals.p")])
    my_pre_return <- c(x[["gene_aachange"]], paste0(res$gene_aachange, "(", res$source.totals.p,")"))
    if(grepl("Ter", x["gene_aachange"])) {
      return(my_pre_return[grepl("Ter", my_pre_return)])
    } else {
      return(my_pre_return[!grepl("Ter", my_pre_return)])
    }
  } else if (any(any.in.n)) {
    res <- unique(vars[vars$CHROM.POS %in% vector.n, c("CHROM.POS", "source.totals.c", "gene_cDNAchange")])
    my_pre_return <- c(x[["gene_cDNAchange"]], paste0(res$gene_cDNAchange, "(", res$source.totals.c,")"))
    if(grepl("del|ins|dup", x["gene_cDNAchange"])) {
      return(my_pre_return[grepl("del|ins|dup", my_pre_return)])
    } else {
      return(my_pre_return[!grepl("del|ins|dup", my_pre_return)])
    }
  } else {
    return("")
  }
})

# toJSON
final$near.BB.HS <- sapply(final$near.BB.HS, function(x) {
  ifelse (x == "", return(x), return(toJSON(x)))
})
final$near.BB.HS <- ifelse(grepl(",", final$near.BB.HS), final$near.BB.HS, "")

ct <- read.table(cosmic_hotspots,sep = "\t", header = T)
ct$gene <- gsub("_.*", "", ct$Gene_HGVSp_VEP)
ct$aa.pos <- gsub(".*_", "", ct$Gene_HGVSp_VEP)
ct$cDNAchange <- gsub(".*:", "", ct$HGVSC)
ct$aa.pos <- as.numeric(str_extract(ct$aa.pos, "\\d+"))
colnames(ct)[2] <- "key"
ct$CHROM.POS <- unlist(lapply(ct$key, function(x) paste(str_split(x,":")[[1]][1],str_split(x,":")[[1]][2],sep = ":")))
ct$GENE.AA.POS <- with(ct, paste(gene, aa.pos, sep=":"))
ct$Gene_HGVSc_VEP <- with(ct, paste(gene, cDNAchange, sep=":"))

final$near.heme.cosmic.HS <- apply(final[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange","gene_cDNAchange")], 1, function(x) {
  p = c(-3:0,1:3)
  n = c(-9:0,1:9)
  
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% ct$GENE.AA.POS
  
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% ct$CHROM.POS
  
  if (any(any.in.p)) {
    my_pre_return <- c(x[["gene_aachange"]],
                       paste(ct$Gene_HGVSp_VEP[ct$GENE.AA.POS %in% vector.p],
                             paste0(ct$GENOMIC_MUTATION_ID[ct$GENE.AA.POS %in% vector.p],
                                    "(",ct$cosmic_count_chr[ct$GENE.AA.POS %in% vector.p],",",
                                    ct$haematopoietic_and_lymphoid_tissue_count_chr[ct$GENE.AA.POS %in% vector.p],",",
                                    ct$myeloid_count_chr[ct$GENE.AA.POS %in% vector.p],")"), 
                             sep="|"))
    if(grepl("Ter", x["gene_aachange"])) {
      return(my_pre_return[grepl("Ter", my_pre_return) | grepl("del|ins|dup", my_pre_return)])
    } else {
      return(my_pre_return[!grepl("Ter", my_pre_return)])
    }
  } else if (any(any.in.n)) {
    my_pre_return <- c(x[["gene_cDNAchange"]],
                       paste(ct$Gene_HGVSc_VEP[ct$CHROM.POS %in% vector.n],
                             paste0(ct$GENOMIC_MUTATION_ID[ct$CHROM.POS %in% vector.n],
                                    "(",ct$cosmic_count_chr[ct$CHROM.POS %in% vector.n],",",
                                    ct$haematopoietic_and_lymphoid_tissue_count_chr[ct$CHROM.POS %in% vector.n],",",
                                    ct$myeloid_count_chr[ct$CHROM.POS %in% vector.n],")"), 
                             sep="|"))
    if(grepl("del|ins|dup", x["gene_cDNAchange"])) {
      return(my_pre_return[grepl("del|ins|dup", my_pre_return)])
    } else {
      return(my_pre_return[!grepl("del|ins|dup", my_pre_return)])
    }
  } else {
    return("")
  }
}) 

# toJSON
final$near.heme.cosmic.HS <- sapply(final$near.heme.cosmic.HS, function(x) {
  ifelse (x == "", return(x), return(toJSON(x)))
})

final$near.BB.loci.HS <- apply(final[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange", "gene_cDNAchange")], 1, function(x) {
  p = c(-3:0,1:3)
  n = c(-9:0,1:9)
  
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% vars$GENE.AA.POS
  
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% vars$CHROM.POS
  
  if (any(any.in.p)) {
    res <- unique(vars[vars$GENE.AA.POS %in% vector.p, c("gene_loci_vep", "truncating", "source.totals.loci.truncating")])
    if(grepl("Ter", x["gene_aachange"])) {
      res <- res %>% filter(ifelse(truncating == "truncating", TRUE, FALSE))
      if (nrow(res) != 0){
        return(c(x[["gene_aachange"]], paste0(res$gene_loci_vep, "(", res$source.totals.loci.truncating,")")))
      } else {
        return("")
      }
    } else {
      res <- res %>% filter(ifelse(truncating == "not", TRUE, FALSE))
      if (nrow(res) != 0){
        return(c(x[["gene_aachange"]], paste0(res$gene_loci_vep, "(", res$source.totals.loci.truncating,")")))
      } else {
        return("") 
      }
    }
  } else if (any(any.in.n)) {
    res <- unique(vars[vars$CHROM.POS %in% vector.n, c("CHROM.POS", "n.HGVSc", "gene_cDNAchange", "gene_loci_vep")])
    if(grepl("del|ins|dup", x["gene_cDNAchange"])) {
      res <- res %>% filter(ifelse(grepl("del|ins|dup",gene_cDNAchange), TRUE, FALSE))
      res <- res %>% group_by(n.HGVSc, gene_loci_vep) %>% summarise(n = sum(n.HGVSc))
      if (nrow(res) != 0){
        return(c(x[["gene_cDNAchange"]], paste0(res$gene_loci_vep, "(", res$n,")"))) 
      } else {
        return("")
      }
    } else {
      res <- res %>% filter(ifelse(grepl("del|ins|dup",gene_cDNAchange), FALSE, TRUE))
      res <- res %>% group_by(gene_loci_vep) %>% summarise(n = sum(n.HGVSc))
      if (nrow(res) != 0){
        return(c(x[["gene_cDNAchange"]], paste0(res$gene_loci_vep, "(", res$n,")")))
      } else {
        return("")  
      }
    }
  } else {
    return("")
  }
})

# toJSON
final$near.BB.loci.HS <- sapply(final$near.BB.loci.HS, function(x) {
  ifelse (x == "", return(x), return(toJSON(x)))
})
final$near.BB.loci.HS <- ifelse(grepl(",", final$near.BB.loci.HS), final$near.BB.loci.HS, "")

final$near.heme.cosmic.loci.HS <- apply(final[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange","gene_cDNAchange")], 1, function(x) {
  p = c(-3:0,1:3)
  n = c(-9:0,1:9)
  
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% ct$GENE.AA.POS
  
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% ct$CHROM.POS
  
  if (any(any.in.p)) {
    res <- unique(ct[ct$GENE.AA.POS %in% vector.p, c("gene_loci_vep", "truncating", "cosmic_count.loci.truncating", "heme_count.loci.truncating", "myeloid_count.loci.truncating")])
    if(grepl("Ter", x["gene_aachange"])) {
      res <- res %>% filter(ifelse(truncating == TRUE, TRUE, FALSE))
      if (nrow(res) != 0){
        return(c(x[["gene_aachange"]], paste0(res$gene_loci_vep, "(", res$cosmic_count.loci.truncating, " ", res$heme_count.loci.truncating, " ", res$myeloid_count.loci.truncating,")")))
      } else {
        return("")
      }
    } else {
      res <- res %>% filter(ifelse(truncating == FALSE, TRUE, FALSE))
      if (nrow(res) != 0){
        return(c(x[["gene_aachange"]], paste0(res$gene_loci_vep, "(", res$cosmic_count.loci.truncating, " ", res$heme_count.loci.truncating, " ", res$myeloid_count.loci.truncating,")")))
      } else {
        return("")
      }
    }
  } else if (any(any.in.n)) {
    res <- unique(ct[ct$CHROM.POS %in% vector.n, c("CHROM.POS", "cDNAchange", "gene_loci_vep", "cosmic_count.totals.c", "heme_count.totals.c", "myeloid_count.totals.c")])
    if(grepl("del|ins|dup", x["cDNAchange"])) {
      res <- res %>% filter(ifelse(grepl("del|ins|dup",cDNAchange), TRUE, FALSE))
      res <- res %>% group_by(cosmic_count.totals.c, heme_count.totals.c, myeloid_count.totals.c, gene_loci_vep) %>% 
        summarise(cosmic_count.totals.c = sum(cosmic_count.totals.c),
                  heme_count.totals.c = sum(heme_count.totals.c),
                  myeloid_count.totals.c = sum(myeloid_count.totals.c))
      if (nrow(res) != 0){
        return(c(x[["cDNAchange"]], paste0(res$gene_loci_vep, "(", res$cosmic_count.totals.c, " ", res$heme_count.totals.c, " ", res$myeloid_count.totals.c,")"))) 
      } else {
        return("")
      }
    } else {
      res <- res %>% filter(ifelse(grepl("del|ins|dup",cDNAchange), FALSE, TRUE))
      res <- res %>% group_by(cosmic_count.totals.c, heme_count.totals.c, myeloid_count.totals.c, gene_loci_vep) %>% 
        summarise(cosmic_count.totals.c = sum(cosmic_count.totals.c),
                  heme_count.totals.c = sum(heme_count.totals.c),
                  myeloid_count.totals.c = sum(myeloid_count.totals.c))
      if (nrow(res) != 0){
        return(c(x[["cDNAchange"]], paste0(res$gene_loci_vep, "(", res$cosmic_count.totals.c, " ", res$heme_count.totals.c, " ", res$myeloid_count.totals.c,")"))) 
      } else {
        return("")
      }
    }
  } else {
    return("")
  }
})

# toJSON
final$near.heme.cosmic.loci.HS <- sapply(final$near.heme.cosmic.loci.HS, function(x) {
  ifelse (x == "", return(x), return(toJSON(x)))
})


final$near.heme.cosmic.loci.HS.logic <- unlist(lapply(final$near.heme.cosmic.loci.HS, function(x) {
  if (x == "") {
    return(FALSE)
  } else {
    near <- fromJSON(x)
    for (i in 2:length(near)) {
      loci_key <- unlist(lapply(strsplit(near[i], "(", fixed = TRUE), "[[", 1))
      totals <- gsub("\\).*", "", gsub(".*\\(", "", near[i]))
      cosmic<-as.numeric(unlist(lapply(strsplit(totals, " ", fixed = TRUE), "[[", 1)))
      heme<-as.numeric(unlist(lapply(strsplit(totals, " ", fixed = TRUE), "[[", 2)))
      myeloid<-as.numeric(unlist(lapply(strsplit(totals, " ", fixed = TRUE), "[[", 3)))
      
      if (cosmic > 25 | (heme >= 10 | myeloid >= 5)) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
}))

final$near.BB.loci.HS.logic <- unlist(lapply(final$near.BB.loci.HS, function(x) {
  total_source <- 0
  if (x == "") {
    return(FALSE)
  } else {
    near <- fromJSON(x)
    # Start at 2 because the first index is always OUR mutation
    for (i in 2:length(near)){
      #loci_key <- unlist(lapply(strsplit(near[i], "(", fixed = TRUE), "[[", 1)) 
      totals <- unlist(lapply(strsplit(near[i], "(", fixed = TRUE), "[[", 2)) 
      # If there is a comma, it means that it has BOTH bick and bolton information, has to be handled differently
      if (grepl(',', totals)) {
        kelly <- unlist(lapply(strsplit(totals, ",", fixed = TRUE), "[[", 1))
        kelly <- gsub("[^0-9.-]", "", unlist(lapply(strsplit(kelly, ":", fixed = TRUE), "[[", 2)))
        bick <- unlist(lapply(strsplit(totals, ",", fixed = TRUE), "[[", 2))
        bick <- gsub("[^0-9.-]", "", unlist(lapply(strsplit(bick, ":", fixed = TRUE), "[[", 2)))
        source <- as.numeric(kelly) + as.numeric(bick)
        # This means it is still protein counts, not cDNA
      } else if (grepl(':', totals)){
        source <- as.numeric(gsub("[^0-9.-]", "", unlist(lapply(strsplit(totals, ":", fixed = TRUE), "[[", 2))))
        # cDNA
      } else {
        source <- as.numeric(gsub("[^0-9.-]", "", totals))
      }
      total_source <- total_source + source
      if (source >= 5 | total_source >= 5) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
}))

if (!("n.loci.truncating.vep" %in% colnames(final))) {
  final$truncating <- ifelse(grepl("Ter", final$AAchange), "truncating", "not")
  dims <- dim(final)[[1]]
  tmp <- final[final$key %in% vars$key | final$gene_loci_vep %in% vars$gene_loci_vep, ]
  tmp <- sqldf("SELECT l.*, r.`n.loci.truncating.vep`, r.`source.totals.loci.truncating`
             FROM `tmp` as l
             LEFT JOIN `vars` as r
             on (l.key = r.key AND l.truncating = r.truncating) OR (l.gene_loci_vep = r.gene_loci_vep AND l.truncating = r.truncating)")
  tmp <- tmp %>% select(key, gene_loci_vep, truncating, n.loci.truncating.vep)
  tmp <- tmp[!duplicated(tmp), ]
  final <- final %>% left_join(tmp %>% select(key, gene_loci_vep, truncating, n.loci.truncating.vep),
                             by = c("key", "gene_loci_vep", "truncating")
  )
  paste0("dims match after sqldf: ", dim(final)[[1]] == dims)
  final <- final %>% dplyr::select(-truncating)
}

# Adding Recurrent Proportions
min_samples_for_recurrent <- max(ceiling(length(unique((final$SN_TAG)))*0.001),2)  # It's 5 in UKBB
total_samples = length(unique(final$SN_TAG)) # 1934 (for Archer)
total_heme_cosmic = 29234
total_BB = 122691
final <- final %>% mutate(prop_nsamples = nsamples_min_vaf/total_samples, 
                          prop_Cosmic = (0.5+heme_cosmic_count)/total_heme_cosmic, 
                          prop_BB = (0.5+sourcetotalsc_XGB)/total_BB) %>%
  mutate(ratio_to_BB = prop_nsamples/prop_BB, ratio_to_cosmic = prop_nsamples/prop_Cosmic)

# These are considered to be recurrent
final <- final %>% mutate(pass_prop_recurrent = case_when(
  Gene %in% unique(bick_email$Gene) ~ TRUE,
  nsamples_min_vaf > min_samples_for_recurrent & sourcetotalsc_XGB == 0 & heme_cosmic_count != 0 & ratio_to_cosmic > 835 ~ FALSE,
  nsamples_min_vaf > min_samples_for_recurrent & sourcetotalsc_XGB != 0 & heme_cosmic_count == 0 & ratio_to_BB > 1335 ~ FALSE,
  nsamples_min_vaf > min_samples_for_recurrent & sourcetotalsc_XGB != 0 & heme_cosmic_count != 0 & (ratio_to_BB > 1335 & ratio_to_cosmic > 835) ~ FALSE,
  TRUE ~ TRUE
)) %>% filter(pass_prop_recurrent == TRUE | (key == 'chr20 32434638 A>AG' & average_AF >= 0.05) | (key == 'chr20 32434638 A>AGG' & average_AF >= 0.05))

# Autofail Variants with NO support in B/B and Cosmic that have nsamples >20
final <- final %>% filter(ifelse(
  sourcetotalsc_XGB == 0 & heme_cosmic_count == 0 & nsamples_min_vaf > 20 & !(Gene %in% unique(bick_email$Gene)), FALSE, TRUE
)) 

# Putative Driver Rules 
# ---------------------
# 1) TSG + Nonsense Mutation --> PD = 1
# 2) OncoKB is Reviewed by Pathologists --> PD = 1
# 3) If OncoKB Reports as 'Neutral' but A LOT of Support from B/B then B/B takes precedence
# 4) OncoKB No Support --> PD = 0
# 5) Missense Variant + Cosmic Support  --> PD = 1
# 6) Missense Variant + B/B Loci Count + SIFT & PolyPhen Support --> PD = 1
# 7) Missense Variant + Near B/B Hotspot | Near Cosmic Hotspot + SIFT & PolyPhen Support --> PD = 1
# 8a) SRSF2 + Hotspot
# 8b) SRSF2 + OncoKB
# 9) SF3B1 Rules
# 10) IDH1 and IDH2 Hotspots
# 11) JAK2 Rules
# 12) PPM1D Exon 6 Rules
# 13) Missense Variant + B/B AA Support w/ EITHER SIFT | PolyPhen Support --> PD = 1
# 14) Missense Variant + Extension of Termination Codon --> PD = 1
# 15) TSG + Splice Acceptor/Donor Variant --> PD = 1
# 16) TSG + ClinVar Support --> PD == 1
# 17) Synonymous Variant in Splicing Region (Handled with SpliceAI)
final <- final %>% mutate(pd_reason = case_when(
  Gene %in% TSG_gene_list & VariantClass %in% nonsense_mutation ~ "Nonsense Mutation in TSG",
  grepl("Oncogenic", oncoKB) ~ "OncoKB",
  Gene %in% gene_list & VariantClass %in% missense_mutation & (sourcetotalsp_XGB >= 10 | sourcetotalsc_XGB >= 5) ~ "B/B Hotspot >= 10",
  Gene %in% gene_list & grepl("Neutral", oncoKB) ~ "Not PD",
  Gene %in% gene_list & VariantClass %in% missense_mutation & (CosmicCount >= 10 | heme_cosmic_count >= 5 | myeloid_cosmic_count >= 1) ~ "COSMIC",
  Gene %in% gene_list & VariantClass %in% missense_mutation & (n.loci.vep - n.loci.truncating.vep) >= 5 & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) ~ "Loci + SIFT/PolyPhen",
  Gene %in% gene_list & VariantClass %in% missense_mutation & (near.BB.loci.HS.logic == TRUE | near.heme.cosmic.loci.HS.logic == TRUE) & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) ~ "Near Hotspot + SIFT/PolyPhen",
  Gene == "SRSF2" & VariantClass %in% missense_mutation & aa.pos == 95 ~ "SRSF2 Hotspot",
  # Gene == "SRSF2" & grepl("Oncogenic", oncoKB) ~ "SRSF2 OncoKB", # redundant
  Gene == "SF3B1" & VariantClass %in% missense_mutation & aa.pos %in% SF3B1_positions ~ "SF3B1 Hotspot",
  Gene == "IDH1" & VariantClass %in% missense_mutation & aa.pos == 132 ~ "IDH1 Hotspot",
  Gene == "IDH2" & VariantClass %in% missense_mutation & (aa.pos == 140 | aa.pos == 172) ~ "IDH2 Hotspot"
  # Gene == "JAK2" & grepl("Oncogenic", oncoKB) ~ "JAK2 OncoKB", # redundant
  Gene == "PPM1D" & VariantClass %in% nonsense_mutation & EXON_VEP == "6/6" ~ "Nonsense Mutation on Exon 6",
  Gene %in% gene_list & VariantClass %in% missense_mutation & sourcetotalsp_XGB >= 1 & (grepl("deleterious", SIFT_VEP) | grepl("damaging", PolyPhen_VEP)) ~ "B/B Hotspot + SIFT/PolyPhen",
  # Gene %in% gene_list & VariantClass %in% missense_mutation & grepl("extTer", gene_aachange) ~ "Termination Extension", # remove
  Gene %in% TSG_gene_list & VariantClass %in% c("splice_donor_variant", "splice_aceptor_variant", "splice_region_variant") & Consequence_VEP != "splice_region_variant&synonymous_variant" ~ "Splicing Mutation",
  Gene %in% TSG_gene_list & clinvar_CLINSIGN_VEP %in% clinvar_sig_terms ~ "ClinVar",
  Gene %in% gene_list & Consequence_VEP == "splice_region_variant&synonymous_variant" ~ "Not PD",
  TRUE ~ "Not PD"
), putative_driver = ifelse(pd_reason != "Not PD", 1, 0))

# OncoKB API automatically classifies ALL splicing mutations as oncogenic, however if the mutation is synonymous, then we have to change it
final <- final %>% mutate(putative_driver = ifelse(grepl("Oncogenic", oncoKB) & Consequence_VEP == "splice_region_variant&synonymous_variant" & !(clinvar_CLINSIGN_VEP %in% clinvar_sig_terms), 0, putative_driver),
                          pd_reason = ifelse(grepl("Oncogenic", oncoKB) & Consequence_VEP == "splice_region_variant&synonymous_variant" & !(clinvar_CLINSIGN_VEP %in% clinvar_sig_terms), "Not PD", pd_reason))

# Additional Bick's Email Rules
unique(bick_email$Gene)
ZBTB33 <- unname(unlist(bick_email %>% filter(Gene == "ZBTB33") %>% mutate(AAchange = paste0(aa_ref, aa_pos, aa_alt)) %>% dplyr::select(AAchange) %>% dplyr::filter(AAchange != "***")))
final <- final %>% mutate(putative_driver = case_when(
  Gene == "ZBTB33" & VariantClass %in% missense_mutation & AAchange.x %in% ZBTB33 ~ 1, 
  Gene %in% unique(bick_email$Gene) & (VariantClass %in% nonsense_mutation) ~ 1,
  TRUE ~ putative_driver
), pd_reason = case_when(
  Gene == "ZBTB33" & VariantClass %in% missense_mutation & AAchange.x %in% ZBTB33 ~ "Bick's Email", 
  Gene %in% unique(bick_email$Gene) & (VariantClass %in% nonsense_mutation) ~ "Bick's Email",
  TRUE ~ pd_reason
))

# Review Stuff
final <- final %>% mutate(Review = "No Review")

# Review Rule for SRSF2 for Matthew Walters
# Need to split Protein_position_VEP up into beginning and end
final <- final %>% mutate(Protein_position_VEP_start = ifelse(grepl("-", Protein_position_VEP), str_split(Protein_position_VEP,"-")[[1]][1], Protein_position_VEP),
                          Protein_position_VEP_end = ifelse(grepl("-", Protein_position_VEP), str_split(Protein_position_VEP,"-")[[1]][2], Protein_position_VEP))

# If it is an Inframe Insertion or Deletion that OVERLAPS P95, then mark for review
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Gene == "SRSF2" & (VariantClass == "inframe_deletion" | VariantClass == "inframe_insertion") & as.numeric(Protein_position_VEP_start) <= 95 & as.numeric(Protein_position_VEP_end) >= 95 ~ "MW Review",
  Review != "No Review" & Gene == "SRSF2" & (VariantClass == "inframe_deletion" | VariantClass == "inframe_insertion") & as.numeric(Protein_position_VEP_start) <= 95 & as.numeric(Protein_position_VEP_end) >= 95 ~ paste(Review, "MW Review", sep = ";"),
  TRUE ~ Review
))

# For SRSF2, SF3B1, IDH1, IDH2, and JAK2 (If not already identified as PD or review... remove all other variants)
SSIIJ = c("SRSF2", "SF3B1", "IDH1", "IDH2", "JAK2")
final <- final %>% filter(ifelse(
  Gene %in% SSIIJ & putative_driver == 0 & Review == "No Review", FALSE, TRUE
))

final <- final %>% mutate(Review = ifelse(CALL_BY_CALLER == "mutect" & (average_AF >= 0.01 & average_AF < 0.02), "LowVAF Mutect", Review))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & nchar(REF) > 5 | nchar(ALT) > 5 ~ "Long INDEL",
  Review != "No Review" & nchar(REF) > 5 | nchar(ALT) > 5 ~ paste(Review, "Long INDEL", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & alt_len >= 2 & ref_len >= 2 ~ "Complex INDEL",
  Review != "No Review" & alt_len >= 2 & ref_len >= 2 ~ paste(Review, "Complex INDEL", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Consequence_VEP == "splice_region_variant&synonymous_variant" ~ "Splice Region Variant",
  Review != "No Review" & Consequence_VEP == "splice_region_variant&synonymous_variant" ~ paste(Review, "Splice Region Variant", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & average_AF >= 0.2 ~ "High VAF",
  Review != "No Review" & average_AF >= 0.2 ~ paste(Review, "High VAF", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & pass_homopolymer_filter == "false" ~ "Homopolymer Region",
  Review != "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & pass_homopolymer_filter == "false" ~ paste(Review, "Homopolymer Region", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) & putative_driver == 0 & sourcetotalsc_XGB >= 1 ~ "B/B Missense Review",
  Review != "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) & putative_driver == 0 & sourcetotalsc_XGB >= 1 ~ paste(Review, "B/B Missense Review", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) & putative_driver == 0 ~ "S/P Missense Review",
  Review != "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) & putative_driver == 0 ~ paste(Review, "S/P Missense Review", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & (near.BB.loci.HS.logic == TRUE | near.heme.cosmic.loci.HS.logic == TRUE) & (grepl("deleterious", SIFT_VEP) | grepl("damaging", PolyPhen_VEP)) & putative_driver == 0 ~ "NHS Missense Review",
  Review != "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & (near.BB.loci.HS.logic == TRUE | near.heme.cosmic.loci.HS.logic == TRUE) & (grepl("deleterious", SIFT_VEP) | grepl("damaging", PolyPhen_VEP)) & putative_driver == 0 ~ paste(Review, "NHS Missense Review", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & (n.loci.vep - n.loci.truncating.vep) >= 5 & (grepl("deleterious", SIFT_VEP) | grepl("damaging", PolyPhen_VEP)) & putative_driver == 0 ~ "HS Missense Review",
  Review != "No Review" & Gene %in% gene_list & VariantClass %in% missense_mutation & (n.loci.vep - n.loci.truncating.vep) >= 5 & (grepl("deleterious", SIFT_VEP) | grepl("damaging", PolyPhen_VEP)) & putative_driver == 0 ~ paste(Review, "HS Missense Review", sep = ";"),
  TRUE ~ Review
))
final <- final %>% mutate(Review = case_when(
  Review == "No Review" & nsamples_min_vaf > min_samples_for_recurrent ~ "Recurrent",
  Review != "No Review" & nsamples_min_vaf > min_samples_for_recurrent ~ paste(Review, "Recurrent", sep = ";"),
  TRUE ~ Review
))

final <- final %>% mutate(Review = case_when(
  Review == "No Review" & Gene %in% unique(bick_email %>% filter(Gene != 'ZBTB33') %>% dplyr::select(Gene)) ~ "Bick's Email",
  Review != "No Review" & Gene %in% unique(bick_email %>% filter(Gene != 'ZBTB33') %>% dplyr::select(Gene)) ~ paste(Review, "Bick's Email", sep = ";"),
  TRUE ~ Review
))

# If something is ONLY in Review because it was recurrent, we can auto pass it if the recurrence wasn't that significant
final <- final %>% mutate(auto_pass_recurrence = case_when(
  Review == "Recurrent" & sourcetotalsc_XGB == 0 & heme_cosmic_count == 0 & (ratio_to_BB <= 200 | ratio_to_cosmic <= 53) ~ TRUE,
  Review == "Recurrent" & sourcetotalsc_XGB == 0 & heme_cosmic_count != 0 & ratio_to_cosmic <= 53 ~ TRUE,
  Review == "Recurrent" & sourcetotalsc_XGB != 0 & heme_cosmic_count == 0 & ratio_to_BB <= 200 ~ TRUE,
  Review == "Recurrent" & sourcetotalsc_XGB != 0 & heme_cosmic_count != 0 & (ratio_to_BB <= 200 & ratio_to_cosmic <= 53) ~ TRUE,
  TRUE ~ FALSE
)) %>% mutate(Review = ifelse(auto_pass_recurrence, "Was Recurrent", Review))

# If Long INDEL or Complex INDEL have at least SOME PINDEL Support (>= 10) then it should be okay.
final$MAX_PINDEL_MATCH <- unlist(lapply(final$PINDEL_MATCH, function(x) {
  if (x == "") {
    return(0)
  } else {
    hold<-unlist(str_extract_all(x, "\\(\\d+\\, \\d+\\)"))
    if (length(hold) > 1) {
      pindel_count <- 0
      for(i in 1:length(hold)){
        temp<-gsub("[^0-9.-]", "", unlist(lapply(strsplit(unlist(hold[i]), ",", fixed = TRUE), "[[", 2)))
        if (temp > pindel_count) {
          pindel_count <- temp 
        }
      }
    } else {
      pindel_count<- gsub("[^0-9.-]", "", unlist(lapply(strsplit(unlist(hold), ",", fixed = TRUE), "[[", 2)))
    }
    return(pindel_count)
  }
}))
final <- final %>% mutate(Review = case_when(
  grepl("INDEL", Review) & MAX_PINDEL_MATCH >= 10 ~ "Pindel Match",
  TRUE ~ Review
))

# If SpliceAI says a Splice Region Variant is a Splice Donor and Acceptor, we change to Putative Driver
if ("SpliceAI_pred_SYMBOL_VEP" %in% colnames(final)){
  final <- final %>% mutate(putative_driver = ifelse(grepl("Splice Region Variant", Review) & (SpliceAI_pred_DS_AG_VEP >= 0.8 | SpliceAI_pred_DS_AL_VEP >= 0.8 | SpliceAI_pred_DS_DG_VEP >= 0.8 | SpliceAI_pred_DS_DL_VEP >= 0.8), 1, putative_driver))
}

final <- final %>% arrange(CHROM, POS, REF, ALT)

review <- final %>% filter(case_when(
  grepl("LowVAF Mutect", Review) & putative_driver == 1 ~ TRUE,
  grepl("Long INDEL", Review) & putative_driver == 1 ~ TRUE,
  grepl("Complex INDEL", Review) & putative_driver == 1 ~ TRUE,
  grepl("High VAF", Review) & putative_driver == 1 ~ TRUE,
  grepl("Missense Review", Review) & putative_driver == 0 ~ TRUE,
  grepl("^Recurrent$", Review) & putative_driver == 1 ~ TRUE,
  grepl("Homopolymer Region", Review) & putative_driver == 1 ~ TRUE,
  grepl("MW Review", Review) & putative_driver == 0 ~ TRUE,
  grepl("Pindel Match", Review) & putative_driver == 1 ~ TRUE,
  grepl("Splice Region Variant", Review) & putative_driver == 0 ~ TRUE
)) %>% filter(ifelse(comp_germline == 0, TRUE, FALSE))

table(review %>% dplyr::select(Review, putative_driver))

passed <- final %>% filter(case_when(
  putative_driver == 1 & Review == "No Review" ~ TRUE,
  putative_driver == 1 & "Missense Review" %in% Review ~ TRUE,
  putative_driver == 1 & Review == "Was Recurrent" ~ TRUE,
  putative_driver == 1 & grepl("MW Review", Review) ~ TRUE,
  putative_driver == 1 & grepl("Splice Region Variant", Review) ~ TRUE
))

table(passed %>% dplyr::select(Review, putative_driver))

check<-c(review$sample_key, passed$sample_key)
table(final %>% filter(ifelse(sample_key %in% check, FALSE, TRUE)) %>% dplyr::select(putative_driver, Review), useNA = "always")

write.csv(final, "final.combined.FPpass.filtered.csv", row.names = FALSE)
write.csv(passed, "passed.combined.FPpass.filtered.csv", row.names = FALSE)
write.csv(review, "review.combined.FPpass.filtered.csv", row.names = FALSE)
