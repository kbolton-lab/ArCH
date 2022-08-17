# Irenaeus Chan
# May 1st, 2022

library(dplyr)
library(stringr)
library(jsonlite)

bolton_bick_vars <- "/Users/irenaeuschan/Documents/Irenaeus/data/bick.bolton.vars3.txt"
cosmic_hotspots <- "/Users/irenaeuschan/Documents/Irenaeus/data/COSMIC.heme.myeloid.hotspot.tsv2"
#Pilot <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/pilot.archer.combined.FPpass.tsv"
#Pilot <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/pilot.mgi.combined.FPpass.tsv"
#Trios <- "/Users/irenaeuschan/Documents/Irenaeus/MGI_Yizhe/trios.combined.tsv"
#Dilution <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/dilution.combined.tsv"
Final <- "/Volumes/bolton/Active/projects/ProstateCancer/TERRA/prostate.final.combined.FPpass.tsv"
Prostate <- "/Users/irenaeuschan/Documents/Irenaeus/ProstateCancer/prostate.final.combined.FPpass.filtered_KB2.csv"
#Final <- "/Users/irenaeuschan/Documents/Irenaeus/ArcherDX/final.combined.FPpass.tsv"
#orig <- "/Users/irenaeuschan/Documents/Irenaeus/ArcherDX/data/variant_review_IC_31722_KB_complete_updated.csv"
#Alex_Filter <- "/Users/irenaeuschan/Documents/Irenaeus/archer_pilot_data/alex_filter.csv"
# alex_filter <- read.csv(Alex_Filter, header = TRUE)
final <- read.table(Final, sep='\t', header = TRUE)

# Remove all variants that fail our general filter "all_fp_pass_XGB"
final <- final %>% dplyr::filter(all_fp_pass_XGB == "true")

# Remove Off-Target & Introns          
final <- final %>% dplyr::filter(Consequence_VEP != "")

# Remove Duplicated Rows
final <- final[!(duplicated(final) | duplicated(final, fromLast = TRUE)), ]

# For the Dilution Data Only
#final <- final %>% mutate(Sample = ifelse(grepl("T", SN_TAG, fixed = TRUE), paste0(unlist(lapply(strsplit(SN_TAG, "T", fixed = TRUE), "[[", 1)), "T"), SN_TAG))
#final <- final %>% mutate(sample_key = paste0(Sample, " ", key))
#dilution_final <- left_join(alex_filter %>% dplyr::select(sample_key, Germline), final, by=c("sample_key"="sample_key"))
#final <- final %>% mutate(SN_TAG = unlist(lapply(strsplit(SN_TAG, "_", fixed = TRUE), "[[", 1)),)

# Dealing with FP_filter
# First filter out which samples failed out other FP filters:
# - PoN 
# - Long INDELs (100 bp)
# - DiNuc and TriNuc in Vardict
# - BCBIO
final <- final %>% filter(pon_FP_pass_XGB == 'true', 
                          long100_indel_pass_XGB == 'true', 
                          long_indel_pass_XGB == 'true', 
                          di_tri_vard_pass_XGB == 'true',
                          bcbio_pass_XGB == 'true')

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

# Calculate the Median VAF per Variant
final <- final %>% left_join(., final %>% group_by(key) %>% summarise(median_VAF = median(average_AF)))

# Filter out samples that failed MCV4 if it had an Alt Depth above 5
final <- final %>% filter(ifelse(FP_Filter_MVC4_XGB == 1 & average_AD > 5, FALSE, TRUE))

# Filter out samples that failed SB1 only if it had an Alt Depth above 10
final <- final %>% filter(ifelse(FP_Filter_SB1_XGB == 1 & average_AD > 10, FALSE, TRUE))

# Recurrent Filter Parameters
# 80 Samples have R882H Mutation in ArcherDX Dataset
count_threshold<-ceiling(length(unique((final$SN_TAG)))*0.06)
bb_count_threshold<-ceiling(length(unique((final$SN_TAG)))*0.03)
n_samples_min_vaf <- 0.001
n_samples_min_count <- 5

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
final <- final %>% dplyr::filter(PASS_BY_1 == 'true')

# Strand Bias Filter - 90/10
final <- final %>% rowwise() %>% 
  mutate(pass_strand_bias = ifelse(sum(pass_strand_bias_Mutect, pass_strand_bias_Lofreq, pass_strand_bias_Vardict, na.rm = TRUE) >= 1, TRUE, FALSE)) %>% ungroup() %>%
  filter(pass_strand_bias)

# Redefine homopolymer filters from string to boolean
final <- final %>%
  mutate(case_NXXX = ifelse(case_NXXX == "true", TRUE, FALSE),
         case_XNXX = ifelse(case_XNXX == "true", TRUE, FALSE),
         case_XXNX = ifelse(case_XXNX == "true", TRUE, FALSE),
         case_XXXN = ifelse(case_XXXN == "true", TRUE, FALSE),
         case_NNXX = ifelse(case_NNXX == "true", TRUE, FALSE),
         case_XNNX = ifelse(case_XNNX == "true", TRUE, FALSE),
         case_XXNN = ifelse(case_XXNN == "true", TRUE, FALSE)) %>%
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

# Filter Homopolymer
rescue <- final %>% filter(!as.logical(pass_homopolymer_filter))
final <- final %>% filter(as.logical(pass_homopolymer_filter))
rescue <- rescue %>% filter(sourcetotalsloci_XGB > 5 | CosmicCount > 25 | heme_cosmic_count > 20 | myeloid_cosmic_count > 10)
final <- rbind(final, rescue)

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
final <- final %>% dplyr::filter(nsamples_min_vaf <= count_threshold | (key == 'chr20 32434638 A>AG' & average_AF >= 0.05) | (key == 'chr20 32434638 A>AGG' & average_AF >= 0.05))        # Recurrent Filter

# Last filter to remove any variants that have a high recurrent count but are not reported inside Kelly's or Bick's dataset
final <- final %>% dplyr::filter(!(nsamples_min_vaf >= bb_count_threshold & source.totals.c == ""))

# Germline Filters
# Keep variants that are below our gnomAD VAF filter (<= gnomAD Population VAF of 0.0005)
final <- final %>% dplyr::rename(pass_max_sub_gnomAD_AF = max_gnomAD_AF)
final <- final %>% rowwise() %>% mutate(max_sub_gnomAD_AF = max(max_gnomAD_AF_VEP, max_gnomADe_AF_VEP, max_gnomADg_AF_VEP, na.rm = TRUE)) %>% ungroup()
final$max_pop_gnomAD_AF <- final %>% dplyr::select(gnomAD_AF_VEP, gnomADe_AF_VEP, gnomADg_AF_VEP) %>% apply(., 1, function(x){max(x, na.rm = T)})
final <- final %>% mutate(germline = ifelse(max_pop_gnomAD_AF <= 0.0005 | (key == 'chr20 32434638 A>AG' & average_AF >= 0.05) | (key == 'chr20 32434638 A>AGG' & average_AF >= 0.05), 0, 1))
final <- final %>% mutate(germline = ifelse(average_AF >= 0.35 & (sourcetotalsc_XGB <= 25 | CosmicCount <= 50), 1, germline))
final <- final %>% mutate(germline = ifelse(average_AF >= 0.25 & average_AF < 0.35 & sourcetotalsc_XGB < 5 & CosmicCount < 25 & max_sub_gnomAD_AF > 0.0001, 1, germline))
final <- final %>% mutate(germline = ifelse(median_VAF >= 0.35 & average_AF >= 0.25, 1, germline))

# Adding INDEL
final <- final %>% mutate(alt_len = nchar(ALT), ref_len = nchar(REF))

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

res <- apply(final[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange", "gene_cDNAchange")], 1, function(x) {
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

if (nrow(final) > 1) {
  final$near.BB.HS <- res 
} else {
  final$near.BB.HS <- list(unlist(as.list(res)))
}

# toJSON
final$near.BB.HS <- sapply(final$near.BB.HS, function(x) {
  if (x == "") {
    return(x)
  } else {
    return(toJSON(x))
  }
})
final$near.BB.HS <- ifelse(grepl(",", final$near.BB.HS), final$near.BB.HS, "")

ct <- read.table(cosmic_hotspots,sep = "\t", quote = "", header = T, comment.char = "")
ct$gene <- gsub("_.*", "", ct$Gene_HGVSp_VEP)
ct$aa.pos <- gsub(".*_", "", ct$Gene_HGVSp_VEP)
ct$cDNAchange <- gsub(".*:", "", ct$HGVSC)
ct$aa.pos <- as.numeric(str_extract(ct$aa.pos, "\\d+"))
colnames(ct)[2] <- "key"
ct$CHROM.POS <- unlist(lapply(ct$key, function(x) paste(str_split(x,":")[[1]][1],str_split(x,":")[[1]][2],sep = ":")))
ct$GENE.AA.POS <- with(ct, paste(gene, aa.pos, sep=":"))
ct$Gene_HGVSp_VEP
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
    return(
      c(x[["gene_aachange"]],
        paste(ct$Gene_HGVSp_VEP[ct$GENE.AA.POS %in% vector.p],
              paste0(ct$GENOMIC_MUTATION_ID[ct$GENE.AA.POS %in% vector.p],
                     "(",ct$cosmic_count_chr[ct$GENE.AA.POS %in% vector.p],",",
                     ct$haematopoietic_and_lymphoid_tissue_count_chr[ct$GENE.AA.POS %in% vector.p],",",
                     ct$myeloid_count_chr[ct$GENE.AA.POS %in% vector.p],")"),
              sep="|"))
    )
  } else if (any(any.in.n)) {
    return(
      c(x[["gene_cDNAchange"]],
        paste(ct$Gene_HGVSc_VEP[ct$CHROM.POS %in% vector.n],
              paste0(ct$GENOMIC_MUTATION_ID[ct$CHROM.POS %in% vector.n],
                     "(",ct$cosmic_count_chr[ct$CHROM.POS %in% vector.n],",",
                     ct$haematopoietic_and_lymphoid_tissue_count_chr[ct$CHROM.POS %in% vector.n],",",
                     ct$myeloid_count_chr[ct$CHROM.POS %in% vector.n],")"),
              sep="|"))
    )
  } else {
    return("")
  }
})

# toJSON
final$near.heme.cosmic.HS <- sapply(final$near.heme.cosmic.HS, function(x) {
  if (x == "") {
    return(x)
  } else {
    return(toJSON(x))
  }
})

# Nonsense Mutations (Stop Gained)
nonsense_gene_list <- c("ASXL1", "CHEK2", "TP53", "DNMT3A", "TET2")
nonsense_mutation <- c("frameshift_variant", "start_lost", "inframe_deletion", "stop_gained", "inframe_insertion")
SF3B1_positions <- c(622, 623, 624, 625, 626, 662, 663, 664, 665, 666, 700, 701, 702, 703, 704, 740, 741, 742)

final <- final %>% mutate(putative_driver = case_when(
  nsamples_min_vaf <= 1 & Gene %in% nonsense_gene_list & VariantClass %in% nonsense_mutation  ~ 1,
  nsamples_min_vaf <= 1 & Gene %in% nonsense_gene_list & VariantClass == "missense_variant" & grepl("Oncogenic", oncoKB) ~ 1,
  nsamples_min_vaf <= 1 & Gene %in% nonsense_gene_list & VariantClass == "missense_variant" & grepl("Neutral", oncoKB) ~ 0,
  nsamples_min_vaf <= 1 & Gene %in% nonsense_gene_list & VariantClass == "missense_variant" & CosmicCount >= 10 & (heme_cosmic_count >= 1 | myeloid_cosmic_count >= 1) ~ 1,
  nsamples_min_vaf <= 1 & Gene %in% nonsense_gene_list & VariantClass == "missense_variant" & n.HGVSp.y >= 5 ~ 1,
  nsamples_min_vaf <= 1 & Gene %in% nonsense_gene_list & VariantClass == "missense_variant" & n.loci.vep >= 1 & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) ~ 1,
  nsamples_min_vaf <= 1 & Gene %in% nonsense_gene_list & VariantClass == "missense_variant" & (near.BB.HS != "" | near.heme.cosmic.HS != "") & grepl("deleterious", SIFT_VEP) & grepl("damaging", PolyPhen_VEP) ~ 1,
  nsamples_min_vaf <= 1 & Gene == "SRSF2" & VariantClass == "missense_variant" & aa.pos == 95 ~ 1,
  nsamples_min_vaf <= 1 & Gene == "SF3B1" & VariantClass == "missense_variant" & aa.pos %in% SF3B1_positions ~ 1,
  nsamples_min_vaf <= 1 & Gene == "JAK2" & grepl("Oncogenic", oncoKB) ~ 1,
  nsamples_min_vaf <= 1 & Gene == "PPM1D" & VariantClass %in% nonsense_mutation & EXON_VEP == "6/6" ~ 1
))

write.csv(final, "final.combined.FPpass.filtered.csv", row.names = FALSE)


