#!/usr/local/bin/Rscript

# Post-ArCH Pipeline Filters (v2.1.0+)
# Created May 21, 2024
# Author: Irenaeus Chan (chani@wustl.edu)

library(optparse)
library(tidyverse)
library(data.table)
source("/opt/bin/post_pipeline/supportFunctions.R", local = TRUE)
#source("~/Documents/Irenaeus/Scripts/R_docker_IChan/post_pipeline/supportFunctions.R")

option_list = list(
  make_option(c("--tsv"), type="character", default=NULL, 
              help="The final combined TSV after all samples are finished running", metavar="character"),
  make_option("--p_value", type="double", default=2.114164905e-6,
              help="The Bonferroni Corrected P-Value [default = %default]", metavar="double"),
  make_option("--bolton_bick_vars", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/bick.bolton.vars3.txt",
              help="The directory path where the 'bick.bolton.vars' data  is stored. [default_path = %default]", metavar="character"),
  make_option("--cosmic", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/COSMIC.heme.myeloid.hotspot.w_truncating_counts.tsv",
              help="The directory path where the 'COSMIC.heme.myeloid.hotspot.w_truncating_counts.tsv' data is stored. [default_path = %default]", metavar="character"),
  make_option("--gene_list", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/oncoKB_CGC_pd_table_disparity_KB_BW.csv",
              help="The directory path where the 'putative driver annotations for specific genes' data is stored. [default_path = %default]", metavar="character"),
  make_option("--pd_table", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/pd_table_kbreview_bick_trunc4_oncoKB_SAFE.filtered_genes_oncoKB_CGC.tsv",
              help="The directory path where the 'putative driver annotations for specific gene and loci' data is stored. [default_path = %default]", metavar="character"),
  make_option("--prefix", type="character", default="final", 
              help="The output prefix e.g. <prefix>.all.csv [default = %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

message("Reading in final file from ArCH...")
df <- fread(opt$tsv, sep='\t', header=TRUE, na.strings=c("NA", "NaN", "NULL", ""))

message("Reading in supporting files...")
vars <- fread(opt$bolton_bick_vars, sep = "\t", header = T)
ct <- fread(opt$cosmic, sep = "\t", header = TRUE)
gene_list <- fread(opt$gene_list, sep = ",", header = TRUE)
pd_table <- fread(opt$pd_table, sep = "\t", header = TRUE)

# Process Flat Database Files
vars <- prepare_vars_file(vars)
ct <- prepare_cosmic_file(ct)
TSG_gene_list <- gene_list[isTSG == 1]$Gene
gene_list <- gene_list$Gene
bick_email <- pd_table[source == "Bick_email"]

# all_fp_pass summarizes:
# - Panel of Normal P-Value
# - Low Allele Frequency
# - Long Indel without Pindel Support
# - BCBIO Filter
# - Z-Score for Allele Frequency
# - Di/Tri-Nucleotide Indel
# - Long Indel (>100bp)
#
# Does NOT include:
# - FP_Filters: Because of the MVC4 and SB1 and average_AD situation
# - gnomAD Filters
#
df <- df %>% filter(as.logical(all_fp_pass))

# In case we did not filter for PoN P-Value during the workflow for exploration
df <- df %>%
  mutate(
    pon_FP_pass = ifelse( # panel of normal p value. If any of the callers' fisher p-values are above the threshold then fail
      is.na(pon_pvalue_Mutect) & is.na(pon_pvalue_Lofreq) & is.na(pon_pvalue_Vardict),
      FALSE,
      ifelse(
        (is.na(pon_pvalue_Mutect) | pon_pvalue_Mutect <= 2.1e-6) & 
          (is.na(pon_pvalue_Lofreq) | pon_pvalue_Lofreq <= 2.1e-6) & 
          (is.na(pon_pvalue_Vardict) | pon_pvalue_Vardict <= 2.1e-6),
        TRUE,
        FALSE
      ))
  )
df <- df %>% filter(pon_FP_pass)

# Handle Varscan's FP Filters
df <- df %>% filter(FP_Filter_DETP20 == 0,
                    FP_Filter_MMQS100 == 0,
                    FP_Filter_MMQSD50 == 0,
                    FP_Filter_NRC == 0,
                    FP_Filter_PB10 == 0,
                    FP_Filter_RLD25 == 0)
# Filter out samples that failed MCV4 if it had an Alt Depth above 5
df <- df %>% filter(ifelse(FP_Filter_MVC4 == 1 & average_AD > 5, FALSE, TRUE))
# Filter out samples that failed SB1 only if it had an Alt Depth above 10
df <- df %>% filter(ifelse(FP_Filter_SB1 == 1 & average_AD > 10, FALSE, TRUE))

# Adding more caller information
df <- df %>%
  mutate(
    CALL_BY_1 = case_when(
      FILTER_Mutect != "" & FILTER_Lofreq == "" & FILTER_Vardict == "" ~ TRUE,
      FILTER_Mutect == "" & FILTER_Lofreq != "" & FILTER_Vardict == "" ~ TRUE,
      FILTER_Mutect == "" & FILTER_Lofreq == "" & FILTER_Vardict != "" ~ TRUE,
      TRUE ~ FALSE
    ),
    FILTER_Mutect = replace_na(FILTER_Mutect, ""),
    FILTER_Lofreq = replace_na(FILTER_Lofreq, ""),
    FILTER_Vardict = replace_na(FILTER_Vardict, ""),
    CALL_BY_CALLER = case_when(
      FILTER_Mutect != "" & FILTER_Lofreq == "" & FILTER_Vardict == "" ~ 'mutect',
      FILTER_Mutect == "" & FILTER_Lofreq != "" & FILTER_Vardict == "" ~ 'lofreq',
      FILTER_Mutect == "" & FILTER_Lofreq == "" & FILTER_Vardict != "" ~ 'vardict',
      FILTER_Mutect != "" & FILTER_Lofreq != "" & FILTER_Vardict == "" ~ 'mutect;lofreq',
      FILTER_Mutect != "" & FILTER_Lofreq == "" & FILTER_Vardict != "" ~ 'mutect;vardict',
      FILTER_Mutect == "" & FILTER_Lofreq != "" & FILTER_Vardict != "" ~ 'lofreq;vardict',
      FILTER_Mutect != "" & FILTER_Lofreq != "" & FILTER_Vardict != "" ~ 'mutect;lofreq;vardict',
    )
  )


df <- df %>% filter(ch_pd >= 1)                                           # Remove variants not annotated by AnnotatePD
df <- df %>% filter(PASS_BY_1)
df <- df %>% filter(!(CALL_BY_CALLER == "mutect" & average_AF < 0.01))    # Filter out Low VAF Mutect Calls
df <- df %>% filter(Consequence_VEP != "synonymous_variant")              # Throw out Silent Mutations
df <- df %>% filter(pass_strand_bias)                                     # Has to pass the 90/10 Strand Bias
df <- df %>% filter(pass_min_alt_count)                                   # Needs to have at least 5 reads with at least 1 count on both Fwd and Rev
df <- df %>% distinct()

# Adding nsamples and median_vaf
df <- df %>% 
  left_join(
    df %>%
      group_by(key) %>%
      mutate(nsamples = n_distinct(subject), median_AF = median(average_AF)) %>%
      distinct(key, nsamples, median_AF),
    by = "key"
  )

# If median_AF >= 35% and NEVER REPORTED in B/B or COSMIC
df <- df %>% filter(!(median_AF >= 0.35 & nsamples > 1 & n.HGVSc == 0 & CosmicCount == 0))

# Recurrent Filter Parameters - Based on 80 Samples have R882H Mutation in ArcherDX Dataset 
count_threshold<- max(ceiling(length(unique((df$subject)))*0.06), 5)
bb_count_threshold<- max(ceiling(length(unique((df$subject)))*0.03), 2)
# We have to save ASXL1 G646W because this is a very recurrent variant that may be removed from our recurrent filter
df <- df %>% filter(nsamples < count_threshold | (key == 'chr20:32434638:A:AG' & average_AF >= 0.05) | (key == 'chr20:32434638:A:AGG' & average_AF >= 0.05))        # Recurrent Filter
# Last filter to remove any variants that have a high recurrent count but are not reported inside Kelly's or Bick's dataset
df <- df %>% dplyr::filter(!(nsamples > bb_count_threshold & (n.HGVSc <= 25 & CosmicCount <= 50)))

nonsense_mutation <- c("frameshift_variant", "stop_lost", "stop_gained", "transcript_ablation")
# Germline Filters
# - Keep variants that are below our gnomAD VAF filter (<= gnomAD Population VAF of 0.0005)
# - Save Nonsense Mutations in DTAP
# - Average VAF <= 0.25 OR ELSE B/B >= 5 OR ELSE CosmicCount >= 25
df <- df %>% 
  mutate(
    germline = case_when(
      Gene %in% c("DNMT3A", "TET2", "ASXL1", "PPM1D") & VariantClass %in% nonsense_mutation ~ FALSE, 
      (max_gnomADe_AF_VEP < 0.005 & max_gnomADg_AF_VEP < 0.005 & max_pop_gnomAD_AF < 0.0005) | key == "chr20:32434638:A:AG" & average_AF >= 0.05 | key == "chr20:32434638:A:AGG" & average_AF >= 0.05 ~ FALSE,
      average_AF >= 0.35 & n.HGVSc == 0 & n.HGVSp == 0 & CosmicCount == 0 ~ TRUE,
      average_AF <= 0.25 | n.HGVSc >= 5 | CosmicCount >= 25 ~ FALSE,
      median_AF >= 0.35 & average_AF >= 0.25 & nsamples > 1 ~ TRUE
    )
  )
# Adding INDEL lengths
df <- df %>% mutate(alt_len = nchar(ALT), ref_len = nchar(REF))
# Filter out Complex INDELs
df <- df %>% filter(ifelse(CALL_BY_1 == TRUE & (nchar(REF) >= 20 | nchar(ALT) >= 20) & PINDEL_MATCH == "", FALSE, TRUE))
# Autofail Variants with NO support in B/B and Cosmic that have nsamples >20
df <- df %>% filter(ifelse(source.totals.c == 0 & heme_cosmic_count == 0 & nsamples > 20 & !(Gene %in% c("SRCAP", "YLPM1", "ZBTB33", "ZNF318")), FALSE, TRUE))

# Add Near Hotspots
df$aa.pos <- as.numeric(str_extract(df$AAchange, "\\d+"))
df$gene_aachange <- with(df, paste(SYMBOL_VEP, AAchange,sep="_"))
df$n.HGVSc[is.na(df$n.HGVSc)] <- 0
df$n.HGVSp[is.na(df$n.HGVSp)] <- 0
df$n.loci.vep[is.na(df$n.loci.vep)] <- 0
df$n.loci.truncating.vep[is.na(df$n.loci.truncating.vep)] <- 0

df$near.BB.loci.HS <- apply(df, 1, function(x) near_BB_loci_HS(x, vars))
df$near.COSMIC.loci.HS <- apply(df, 1, function(x) near_COSMIC_loci_HS(x, ct))
df$nearBBLogic <- df$near.BB.loci.HS != ''
df$nearCosmicHemeLogic <- df$near.COSMIC.loci.HS != ''

# Adding Recurrent Proportions
min_samples_for_recurrent <- max(ceiling(length(unique((df$subject)))*0.001),2)  # It's 5 in UKBB
total_samples = length(unique(df$subject)) # 1934 (for Archer)
total_heme_cosmic = 29234
total_BB = 122691
df <- df %>% mutate(prop_nsamples = nsamples/total_samples, 
                          prop_Cosmic = (0.5+heme_cosmic_count)/total_heme_cosmic, 
                          prop_BB = (0.5+n.HGVSc)/total_BB) %>%
  mutate(ratio_to_BB = prop_nsamples/prop_BB, ratio_to_cosmic = prop_nsamples/prop_Cosmic)
# These are considered to be recurrent
df <- df %>% mutate(pass_prop_recurrent = case_when(
  Gene %in% unique(bick_email$Gene) ~ TRUE,
  nsamples > min_samples_for_recurrent & n.HGVSc == 0 & heme_cosmic_count != 0 & ratio_to_cosmic > 835 ~ FALSE,
  nsamples > min_samples_for_recurrent & n.HGVSc != 0 & heme_cosmic_count == 0 & ratio_to_BB > 1335 ~ FALSE,
  nsamples > min_samples_for_recurrent & n.HGVSc != 0 & heme_cosmic_count != 0 & (ratio_to_BB > 1335 & ratio_to_cosmic > 835) ~ FALSE,
  TRUE ~ TRUE
)) 
df <- df %>% filter(pass_prop_recurrent == TRUE | (key == 'chr20:32434638:A:AG' & average_AF >= 0.05) | (key == 'chr20:32434638:A:AGG' & average_AF >= 0.05))

# DETERMINE PATHOGENICITY
nonsense_mutation <- c("frameshift_variant", "stop_lost", "stop_gained", "transcript_ablation")
missense_mutation <- c("missense_variant", "inframe_deletion", "inframe_insertion")
SF3B1_positions <- c(622, 623, 624, 625, 626, 662, 663, 664, 665, 666, 700, 701, 702, 703, 704, 740, 741, 742)
clinvar_sig_terms <- c("Likely_pathogenic", "Likely_pathogenic&_drug_response", "Pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic/Likely_pathogenic&_drug_response")
splicingSynonymous = c("splice_donor_5th_base_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", "splice_region_variant,intron_variant", "splice_region_variant,non_coding_transcript_exon_variant", "synonymous_variant", "splice_region_variant,synonymous_variant")
ZBTB33 <- bick_email[Gene == "ZBTB33"] %>% mutate(AAchange = paste0(aa_ref, aa_pos, aa_alt)) %>% pull(AAchange)
ZBTB33 <- ZBTB33[ZBTB33 != "***"]

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
# 18) Bick's Email Rules - ZBTB33
# 19) Bick's Email Rules - Other Genes
df <- df %>% mutate(pd_reason = case_when(
  Gene %in% TSG_gene_list & VariantClass %in% nonsense_mutation ~ "Nonsense Mutation in TSG",
  grepl("Oncogenic", oncoKB) ~ "OncoKB",
  Gene %in% gene_list & VariantClass %in% missense_mutation & (n.HGVSp >= 10 | n.HGVSc >= 5) ~ "B/B Hotspot >= 10",
  Gene %in% gene_list & grepl("Neutral", oncoKB) ~ "Not PD",
  Gene %in% gene_list & VariantClass %in% missense_mutation & (CosmicCount >= 10 | heme_cosmic_count >= 5 | myeloid_cosmic_count >= 1) ~ "COSMIC",
  Gene %in% gene_list & VariantClass %in% missense_mutation & n.loci.truncating.vep >= 5 & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) ~ "Loci + SIFT/PolyPhen",
  Gene %in% gene_list & VariantClass %in% missense_mutation & (nearBBLogic == TRUE | nearCosmicHemeLogic == TRUE) & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) ~ "Near Hotspot + SIFT/PolyPhen",
  Gene == "SRSF2" & VariantClass %in% missense_mutation & aa.pos == 95 ~ "SRSF2 Hotspot",
  # Gene == "SRSF2" & grepl("Oncogenic", oncoKB) ~ "SRSF2 OncoKB", # redundant
  Gene == "SF3B1" & VariantClass %in% missense_mutation & aa.pos %in% SF3B1_positions ~ "SF3B1 Hotspot",
  Gene == "IDH1" & VariantClass %in% missense_mutation & aa.pos == 132 ~ "IDH1 Hotspot",
  Gene == "IDH2" & VariantClass %in% missense_mutation & (aa.pos == 140 | aa.pos == 172) ~ "IDH2 Hotspot",
  # Gene == "JAK2" & grepl("Oncogenic", oncoKB) ~ "JAK2 OncoKB", # redundant
  Gene == "PPM1D" & VariantClass %in% nonsense_mutation & EXON_VEP == "6/6" ~ "Nonsense Mutation on Exon 6",
  Gene %in% gene_list & VariantClass %in% missense_mutation & n.HGVSp >= 1 & (grepl("deleterious", SIFT_VEP) | grepl("damaging", PolyPhen_VEP)) ~ "B/B Hotspot + SIFT/PolyPhen",
  # Gene %in% gene_list & VariantClass %in% missense_mutation & grepl("extTer", gene_aachange) ~ "Termination Extension", # remove
  Gene %in% TSG_gene_list & VariantClass %in% c("splice_donor_variant", "splice_acceptor_variant", "splice_region_variant") & !grepl(paste(splicingSynonymous, collapse = "|"), Consequence_VEP) ~ "Splicing Mutation",
  Gene %in% TSG_gene_list & clinvar_CLNSIG_VEP %in% clinvar_sig_terms ~ "ClinVar",
  Gene %in% gene_list & grepl(paste(splicingSynonymous, collapse = "|"), Consequence_VEP) ~ "Not PD",
  Gene == "ZBTB33" & VariantClass %in% missense_mutation & AAchange %in% ZBTB33 ~ "Bick's Email",
  Gene %in% unique(bick_email$Gene) & VariantClass %in% nonsense_mutation ~ "Bick's Email",
  TRUE ~ "Not PD"
), putative_driver = ifelse(pd_reason != "Not PD", 1, 0))

# OncoKB API automatically classifies ALL splicing mutations as oncogenic, however if the mutation is synonymous, then we have to change it
df <- df %>% mutate(putative_driver = ifelse(grepl("Oncogenic", oncoKB) & grepl(paste(splicingSynonymous, collapse = "|"), Consequence_VEP) & !(clinvar_CLNSIG_VEP %in% clinvar_sig_terms), 0, putative_driver),
                    pd_reason = ifelse(grepl("Oncogenic", oncoKB) & grepl(paste(splicingSynonymous, collapse = "|"), Consequence_VEP) & !(clinvar_CLNSIG_VEP %in% clinvar_sig_terms), "Not PD", pd_reason))

# Start of Review
df$Review <- ""

# Review Rule for SRSF2 for Matthew Walters
# Need to split Protein_position_VEP up into beginning and end
df <- df %>% mutate(Protein_position_start = ifelse(grepl("-", Protein_position_VEP), str_split(Protein_position_VEP,"-")[[1]][1], Protein_position_VEP),
                          Protein_position_end = ifelse(grepl("-", Protein_position_VEP), str_split(Protein_position_VEP,"-")[[1]][2], Protein_position_VEP))
# If it is an Inframe Insertion or Deletion that OVERLAPS P95, then mark for review
df <- df %>% mutate(
  Review = ifelse(Gene == "SRSF2" & VariantClass %in% c("inframe_deletion", "inframe_insertion") & as.numeric(Protein_position_start) <= 95 & as.numeric(Protein_position_end) >= 95, "MW Review", Review)
) 

# For SRSF2, SF3B1, IDH1, IDH2, and JAK2 (If not already identified as PD or review... remove all other variants)
df <- df %>% filter(ifelse(
  Gene %in% c("SRSF2", "SF3B1", "IDH1", "IDH2", "JAK2") & putative_driver == 0 & Review == "", FALSE, TRUE
))

review_conditions <- list(
  list((df$CALL_BY_CALLER == "mutect") & (df$average_AF >= 0.01 & df$average_AF < 0.02), "LowVAF Mutect"),
  list((nchar(df$REF) > 5) | (nchar(df$ALT) > 5), "Long INDEL"),
  list((nchar(df$REF) >= 2) & (nchar(df$ALT) >= 2), "Complex INDEL"),
  list((df$average_af >= 0.2), "High VAF"),
  list((df$Gene %in% gene_list) & (df$VariantClass %in% missense_mutation) & (!is.na(df$homopolymerCase)), "Homopolymer Region"),
  list((df$Gene %in% gene_list) & (df$VariantClass %in% missense_mutation) & grepl("deleterious", df$SIFT_VEP) & grepl("damaging", df$PolyPhen_VEP) & (df$putative_driver == 0) & (df$n.HGVSc >= 1), "B/B Missense Review"),
  list((df$Gene %in% gene_list) & (df$VariantClass %in% missense_mutation) & grepl("deleterious", df$SIFT_VEP) & grepl("damaging", df$PolyPhen_VEP) & (df$putative_driver == 0), "S/P Missense Review"),
  list((df$Gene %in% gene_list) & (df$VariantClass %in% missense_mutation) & (df$nearBBLogic | df$nearCosmicHemeLogic) & (grepl("deleterious", df$SIFT_VEP) | grepl("damaging", df$PolyPhen_VEP)) & (df$putative_driver == 0), "NHS Missense Review"),
  list((df$Gene %in% gene_list) & (df$VariantClass %in% missense_mutation) & ((df$n.loci.vep - df$n.loci.truncating.vep) >= 5) & (grepl("deleterious", df$SIFT_VEP) | grepl("damaging", df$PolyPhen_VEP)) & (df$putative_driver == 0), "HS Missense Review"),
  list((df$nsamples > 5) & (df$n.HGVSc < 25) & (df$CosmicCount < 50), "Recurrent"),
  list((df$Gene %in% c("SRCAP", "YLPM1", "ZNF318")), "Bick's Email")
)
for (condition in review_conditions) {
  df$Review[condition[[1]]] <- paste(df$Review[condition[[1]]], condition[[2]], sep = "|")
}
df$Review[df$Review == ""] <- "No Review"

# If something is ONLY in Review because it was recurrent, we can auto pass it if the recurrence wasn't that significant
df <- df %>% mutate(auto_pass_recurrence = case_when(
  Review == "|Recurrent" & (ratio_to_BB <= 200 | ratio_to_cosmic <= 53) ~ TRUE,
  Review == "|Recurrent" & (n.HGVSc != 0 & heme_cosmic_count == 0 & ratio_to_BB <= 200) ~ TRUE,
  TRUE ~ FALSE
)) %>% mutate(Review = ifelse(auto_pass_recurrence, "Was Recurrent", Review))

# If Long INDEL or Complex INDEL have at least SOME PINDEL Support (>= 10) then it should be okay.
df$pindel_AD <- as.numeric(sapply(strsplit(as.character(df$PINDEL_MATCH), ";"), "[", 2))
df$pindel_AD[is.na(df$pindel_AD)] <- 0
df <- df %>% mutate(Review = ifelse(grepl("INDEL", Review) & pindel_AD >= 10, "Pindel Match", Review))

# If SpliceAI says a Splice Region Variant is a Splice Donor and Acceptor, we change to Putative Driver
if ("SpliceAI_pred_SYMBOL_VEP" %in% colnames(df)){
  df <- df %>% mutate(putative_driver = ifelse(grepl("Splice Region Variant", Review) & (SpliceAI_pred_DS_AG_VEP >= 0.8 | SpliceAI_pred_DS_AL_VEP >= 0.8 | SpliceAI_pred_DS_DG_VEP >= 0.8 | SpliceAI_pred_DS_DL_VEP >= 0.8), 1, putative_driver),
                      pd_reason = ifelse(grepl("Splice Region Variant", Review) & (SpliceAI_pred_DS_AG_VEP >= 0.8 | SpliceAI_pred_DS_AL_VEP >= 0.8 | SpliceAI_pred_DS_DG_VEP >= 0.8 | SpliceAI_pred_DS_DL_VEP >= 0.8), "SpliceAI", pd_reason)
  )
}
df <- df %>% arrange(CHROM, POS, REF, ALT)

review <- df %>% filter(case_when(
  grepl("LowVAF Mutect", Review) & putative_driver == 1 ~ TRUE,
  grepl("Long INDEL", Review) & putative_driver == 1 ~ TRUE,
  grepl("Complex INDEL", Review) & putative_driver == 1 ~ TRUE,
  grepl("High VAF", Review) & putative_driver == 1 ~ TRUE,
  grepl("^\\|Recurrent$", Review) & putative_driver == 1 ~ TRUE,
  grepl("Homopolymer Region", Review) & putative_driver == 1 ~ TRUE,
  grepl("Bick's Email", Review) & putative_driver == 1 ~ TRUE,
  grepl("MW Review", Review) & putative_driver == 0 ~ TRUE,
  grepl("Pindel Match", Review) & putative_driver == 1 ~ TRUE,
  TRUE ~ FALSE
)) %>% filter(ifelse(!germline, TRUE, FALSE))
table(review %>% dplyr::select(Review, putative_driver))

passed <- df %>% filter(case_when(
  putative_driver == 1 & Review == "No Review" ~ TRUE,
  putative_driver == 1 & grepl("Missense Review", Review) ~ TRUE,
  putative_driver == 1 & Review == "Was Recurrent" ~ TRUE,
  putative_driver == 1 & grepl("MW Review", Review) ~ TRUE,
  putative_driver == 1 & grepl("Splice Region Variant", Review) ~ TRUE,
  TRUE ~ FALSE
))
table(passed %>% dplyr::select(Review, putative_driver))

check<-c(review$sample_key, passed$sample_key)
table(df %>% filter(ifelse(sample_key %in% check, FALSE, TRUE)) %>% dplyr::select(putative_driver, Review), useNA = "always")

message("Finished determining pathogenicity...")

write.csv(df, paste0(opt$prefix, ".all.csv"), row.names = FALSE)
write.csv(passed, paste0(opt$prefix, ".pass.csv"), row.names = FALSE)
write.csv(review, paste0(opt$prefix, ".review.csv"), row.names = FALSE)
