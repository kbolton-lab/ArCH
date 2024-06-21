#!/usr/local/bin/Rscript

library(optparse)
library(tidyverse)
library(data.table)
library(vcfR)
library(purrr)
source("/opt/bin/combine/supportFunctions.R", local = TRUE)
#source("/storage1/fs1/bolton/Active/Users/IrenaeusChan/OptimizationsWDL/xgb/supportFunctions.R")

option_list = list(
    make_option(c("--mutect"), type="character", default=NULL, 
              help="Mutect TSV from AnnotatePD e.g. mutect.sample_name.final.annotated.tsv", metavar="character"),
    make_option(c("--lofreq"), type="character", default=NULL, 
              help="Lofreq TSV from AnnotatePD e.g. lofreq.sample_name.final.annotated.tsv", metavar="character"),
    make_option(c("--vardict"), type="character", default=NULL, 
              help="Vardict TSV from AnnotatePD e.g. vardict.sample_name.final.annotated.tsv", metavar="character"),
    make_option(c("--pindel"), type="character", default=NULL, 
              help="Pindel VCF e.g. pindel.sample_name.vcf.gz", metavar="character"),
    make_option(c("--pon"), type="character", default=NULL, 
              help="The PoN Pileup from mskGetBaseCounts e.g. sample_name.pon.total.counts.vcf.gz", metavar="character"),
    make_option("--p_value", type="double", default=2.114164905e-6,
              help="The Bonferroni Corrected P-Value [default = %default]", metavar="double"),
    make_option("--sample_name", type="character", default="TUMOR",
              help="The Sample Name [default = %default]", metavar="character"),
    make_option("--model", type="logical", default=FALSE,
            help="Whether to output a dataframe that can be used in the model or not [default = %default]", metavar="logical")
); 
opt_parser = OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# Reading in Variant Callers
M <- fread(opt$mutect, header=TRUE, sep="\t", na.strings=c("NA", "NaN", "NULL", ""), nrows = -1)
L <- fread(opt$lofreq, header=TRUE, sep="\t", na.strings=c("NA", "NaN", "NULL", ""), nrows = -1)
V <- fread(opt$vardict, header=TRUE, sep="\t", na.strings=c("NA", "NaN", "NULL", ""), nrows = -1)

# Reading in Pindel
P <- vcfR2tidy(read.vcfR(opt$pindel), single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat
P <- P %>% separate(gt_AD, c("gt_AD_ref", "gt_AD_alt"), sep = ",", extra = "merge", fill = "right") %>%
            mutate(
                key = paste0(CHROM, ":", POS, ":", REF, ":", ALT),
                PINDEL_MATCH = paste0(gt_AD_ref, ";", gt_AD_alt),
            ) %>%
            select(key, PINDEL_MATCH) %>%
            distinct(key, .keep_all = TRUE)

# Using the PoN Pileup, we can calculate the mean and standard deviation of the VAFs
pileup <- vcfR2tidy(read.vcfR(opt$pon), single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat
pileup <- pileup %>%
    mutate(key = paste0(CHROM, ":", POS, ":", REF, ":", ALT)) %>%
    group_by(key) %>%
    mutate(
        PONvafmean = mean(gt_VF, na.rm = TRUE),
        PONvafstd = sd(gt_VF, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(key, PONvafmean, PONvafstd) %>% 
    distinct(key, .keep_all = TRUE) 

# Combine into one dataframe
message("Comining all callers into 1...")
final <- combine_all_callers(M, L, V)

final <- final %>%
left_join(P, by = c("key")) %>%
left_join(pileup, by = c("key")) %>%
mutate(
    sample_key = paste0(SAMPLE, " ", key),
    ref_len = nchar(REF),
    alt_len = nchar(ALT),
    gt_AD_ref_Mutect_z = abs(gt_AD_ref_Mutect - mean(gt_AD_ref_Mutect, na.rm = TRUE)) / sd(gt_AD_ref_Mutect, na.rm = TRUE),
    gt_AD_ref_Lofreq_z = abs(gt_AD_ref_Lofreq - mean(gt_AD_ref_Lofreq, na.rm = TRUE)) / sd(gt_AD_ref_Lofreq, na.rm = TRUE),
    gt_AD_ref_Vardict_z = abs(gt_AD_ref_Vardict - mean(gt_AD_ref_Vardict, na.rm = TRUE)) / sd(gt_AD_ref_Vardict, na.rm = TRUE),
    Fwd_StrandBias_Mutect = AltFwd_Mutect/(AltFwd_Mutect+AltRev_Mutect) > 0.9 | AltFwd_Mutect/(AltFwd_Mutect+AltRev_Mutect) < 0.1,
    Rev_StrandBias_Mutect = AltRev_Mutect/(AltFwd_Mutect+AltRev_Mutect) > 0.9 | AltRev_Mutect/(AltFwd_Mutect+AltRev_Mutect) < 0.1,
    Fwd_StrandBias_Lofreq = AltFwd_Lofreq/(AltFwd_Lofreq+AltRev_Lofreq) > 0.9 | AltFwd_Lofreq/(AltFwd_Lofreq+AltRev_Lofreq) < 0.1,
    Rev_StrandBias_Lofreq = AltRev_Lofreq/(AltFwd_Lofreq+AltRev_Lofreq) > 0.9 | AltRev_Lofreq/(AltFwd_Lofreq+AltRev_Lofreq) < 0.1,
    Fwd_StrandBias_Vardict = AltFwd_Vardict/(AltFwd_Vardict+AltRev_Vardict) > 0.9 | AltFwd_Vardict/(AltFwd_Vardict+AltRev_Vardict) < 0.1,
    Rev_StrandBias_Vardict = AltRev_Vardict/(AltFwd_Vardict+AltRev_Vardict) > 0.9 | AltRev_Vardict/(AltFwd_Vardict+AltRev_Vardict) < 0.1,
    fail_strand_bias_Mutect = Fwd_StrandBias_Mutect & Rev_StrandBias_Mutect,
    fail_strand_bias_Lofreq = Fwd_StrandBias_Lofreq & Rev_StrandBias_Lofreq,
    fail_strand_bias_Vardict = Fwd_StrandBias_Vardict & Rev_StrandBias_Vardict,
    pass_strand_bias = !fail_strand_bias_Mutect | !fail_strand_bias_Lofreq | !fail_strand_bias_Vardict,
    pass_min_alt_count_Mutect = gt_AD_alt_Mutect >= 5 & (AltFwd_Mutect != 0 & AltRev_Mutect != 0),
    pass_min_alt_count_Lofreq = gt_AD_alt_Lofreq >= 5 & (AltFwd_Lofreq != 0 & AltRev_Lofreq != 0),
    pass_min_alt_count_Vardict = gt_AD_alt_Vardict >= 5 & (AltFwd_Vardict != 0 & AltRev_Vardict != 0),
    pass_min_alt_count = pass_min_alt_count_Mutect | pass_min_alt_count_Lofreq | pass_min_alt_count_Vardict
)

# If pass_strand_bias and pass_min_alt_count is NA, we set to FALSE
final <- final %>% mutate(
    pass_strand_bias = ifelse(is.na(pass_strand_bias), FALSE, pass_strand_bias),
    pass_min_alt_count = ifelse(is.na(pass_min_alt_count), FALSE, pass_min_alt_count)
)

final$average_AF <- rowMeans(subset(final, select = c("gt_AF_Mutect", "gt_AF_Lofreq", "gt_AF_Vardict")), na.rm = TRUE)
final$average_AD <- rowMeans(subset(final, select = c("gt_AD_alt_Mutect", "gt_AD_alt_Lofreq", "gt_AD_alt_Vardict")), na.rm = TRUE)
final <- final %>% mutate(pon_af_zscore = abs(average_AF - PONvafmean) / PONvafstd)

# TODO: Call Complex - Maybe scrap it completely...

# Split Filters 
final <- split_filters(final)

# Prepare Model Dataframe?
if (opt$model) {
    model_df <- prepare_for_model(final)
    write.table(model_df, file = paste0(opt$sample_name, ".model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Final Set of FP Filters
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
final <- final %>%
    mutate(
        pon_FP_pass = ifelse( # panel of normal p value. If any of the callers' fisher p-values are above the threshold then fail
        is.na(pon_pvalue_Mutect) & is.na(pon_pvalue_Lofreq) & is.na(pon_pvalue_Vardict),
        FALSE,
        ifelse(
            (is.na(pon_pvalue_Mutect) | pon_pvalue_Mutect <= opt$p_value) & 
            (is.na(pon_pvalue_Lofreq) | pon_pvalue_Lofreq <= opt$p_value) & 
            (is.na(pon_pvalue_Vardict) | pon_pvalue_Vardict <= opt$p_value),
            TRUE,
            FALSE
        )),
        low_AF_pass = average_AF >= 0.001,
        long_indel_pass = !(ref_len > 20 | alt_len > 20) | !is.na(PINDEL_MATCH),
        bcbio_pass = FILTER_Vardict_BCBIO != 1,
        zscore_pass = pon_af_zscore >= 1.96,
        di_tri_vard_pass = !(ref_len == alt_len & ref_len > 1) | !is.na(PINDEL_MATCH),
        long100_indel_pass = !(ref_len > 100 | alt_len > 100),
        PASS_BY_1 = FILTER_Mutect_PASS + FILTER_Lofreq_PASS + FILTER_Vardict_PASS >= 1,
        PASS_BY_2 = FILTER_Mutect_PASS + FILTER_Lofreq_PASS + FILTER_Vardict_PASS >= 2,
        PASS_BY_3 = FILTER_Mutect_PASS + FILTER_Lofreq_PASS + FILTER_Vardict_PASS >= 3,
        all_fp_pass = pon_FP_pass & low_AF_pass & long_indel_pass & bcbio_pass & zscore_pass & di_tri_vard_pass & long100_indel_pass
    )

message("Summarizing final dataframe...")
final <- final %>%
select(
    CHROM, POS, ID, REF, ALT, key, sample_key, SAMPLE, subject, Gene, average_AF, average_AD,
    ends_with("_Mutect"), gt_AD_ref_Mutect_z, ends_with("_Lofreq"), gt_AD_ref_Lofreq_z, ends_with("_Vardict"), gt_AD_ref_Vardict_z, PINDEL_MATCH,
    FP_filter, starts_with("FP_Filter_"), pon_FP_pass, low_AF_pass, long_indel_pass, bcbio_pass, zscore_pass, di_tri_vard_pass, long100_indel_pass, all_fp_pass,
    PASS_BY_1, PASS_BY_2, PASS_BY_3,
    ends_with("_VEP"), starts_with("clinvar"), 
    ref_len, alt_len, context_5, context_3, context_5_3, dust_score, homopolymerCase, pass_strand_bias, pass_min_alt_count,
    PON_RefDepth, PON_AltDepth, pon_pvalue_Mutect, pon_pvalue_Lofreq, pon_pvalue_Vardict, PONvafmean, PONvafstd, pon_af_zscore,
    max_gnomADe_AF_VEP, max_gnomADg_AF_VEP, pass_max_sub_gnomAD_AF, max_sub_gnomAD_AF, max_pop_gnomAD_AF,
    VariantClass, AAchange, gene_loci_p, gene_loci_c, gene_aachange, gene_cDNAchange, 
    source.totals.loci, source.totals.loci.truncating, source.totals.p, source.totals.c, n.loci.vep, n.loci.truncating.vep, n.HGVSp, n.HGVSc,
    COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count, oncoKB, isOncogenic, isTSG, isTruncatingHotSpot,
    ch_pd, WHY_CH
)

message("Writing out final annotated file")
write.table(final, file = paste0(opt$sample_name, ".final.annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

