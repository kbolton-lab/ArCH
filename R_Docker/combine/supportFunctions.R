combine_all_callers <- function(M, L, V) {

  #M <- M %>% mutate_if(is.character, as.character) %>%
  #  mutate_if(ends_with("AF_VEP"), as.numeric)
  M <- M %>% mutate(across(where(is.character), as.character),
                  across(ends_with("AF_VEP"), as.numeric))
  #L <- L %>% mutate_if(is.character, as.character) %>%
  #  mutate_if(ends_with("AF_VEP"), as.numeric)
  L <- L %>% mutate(across(where(is.character), as.character),
                  across(ends_with("AF_VEP"), as.numeric))
  #V <- V %>% mutate_if(is.character, as.character) %>%
  #  mutate_if(ends_with("AF_VEP"), as.numeric)
  V <- V %>% mutate(across(where(is.character), as.character),
                  across(ends_with("AF_VEP"), as.numeric))

  # Small checks to account for merging issues with data types
  M$PUBMED_VEP <- as.character(M$PUBMED_VEP)
  L$PUBMED_VEP <- as.character(L$PUBMED_VEP)
  V$PUBMED_VEP <- as.character(V$PUBMED_VEP)

  # We can create a single "annotation dataframe"
  A <- bind_rows(M, L, V) %>% 
  select(
    CHROM, POS, ID, REF, ALT,
    key,
    SAMPLE, Indiv, subject,
    PON_RefDepth, PON_AltDepth,
    FP_filter,
    ends_with("_VEP"),
    starts_with("clinvar"),
    pass_max_sub_gnomAD_AF, max_sub_gnomAD_AF, max_pop_gnomAD_AF,
    context_5, context_3, context_5_3,
    dust_score,
    homopolymerCase,
    AAchange,
    gene_loci_p, gene_loci_c, gene_loci_vep, gene_aachange, gene_cDNAchange,
    n.loci.vep, source.totals.loci, n.loci.truncating.vep, source.totals.loci.truncating, n.HGVSp, source.totals.p, n.HGVSc, source.totals.c,
    COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count,
    VariantClass, Gene, oncoKB, isOncogenic, isTSG, isTruncatingHotSpot,
    ch_pd, WHY_CH
  ) %>% 
  distinct(key, .keep_all = TRUE) 

  A$n.loci.vep <- tidyr::replace_na(A$n.loci.vep, 0)
  A$n.loci.truncating.vep <- tidyr::replace_na(A$n.loci.truncating.vep, 0)

  M <- M %>% select(
      QUAL_mutect, FILTER_mutect, 
      AS_FilterStatus,
      AS_SB_TABLE,
      AS_UNIQ_ALT_READ_COUNT, CONTQ, DP, ECNT, GERMQ, MBQ, MFRL, MMQ, MPOS,NALOD,
      NCount, NLOD, OCM, PON, POPAF, ROQ, RPA, RU, SEQQ, STR, STRANDQ, STRQ, TLOD,
      PON_2AT2_percent, PON_NAT2_percent, PON_MAX_VAF,PON_FISHER,
      ends_with("_mutect"),
      key
  ) %>% rename_with(~paste0(., "_mutect"), -contains("_mutect")) %>%
  rename_with(~sub("_mutect$", "_Mutect", .x), ends_with("_mutect")) %>%
  dplyr::rename(
      key = key_Mutect,
      pon_pvalue_Mutect = PON_FISHER_Mutect,
      RefFwd_Mutect = RDF_Mutect,
      RefRev_Mutect = RDR_Mutect,
      AltFwd_Mutect = ADF_Mutect,
      AltRev_Mutect = ADR_Mutect
  ) %>%
  group_by(key) %>%
  mutate(
    SBF_pvalue_Mutect = fisher.test(matrix(c(RefFwd_Mutect, RefRev_Mutect, AltFwd_Mutect, AltRev_Mutect), nrow = 2), simulate.p.value = TRUE)$p.value
  ) %>%
  ungroup()

  L <- L %>% select(
      QUAL_lofreq, FILTER_lofreq,
      DP, AF, SB, RDF_lofreq, RDR_lofreq, ADF_lofreq, ADR_lofreq, 
      INDEL, CONSVAR, HRUN,
      PON_2AT2_percent, PON_NAT2_percent, PON_MAX_VAF, PON_FISHER,
      ends_with("_lofreq"),
      key
  ) %>% rename_with(~paste0(., "_lofreq"), -contains("_lofreq")) %>%
  rename_with(~sub("_lofreq$", "_Lofreq", .x), ends_with("_lofreq")) %>%
  dplyr::rename(
      key = key_Lofreq,
      pon_pvalue_Lofreq = PON_FISHER_Lofreq,
      RefFwd_Lofreq = RDF_Lofreq,
      RefRev_Lofreq = RDR_Lofreq,
      AltFwd_Lofreq = ADF_Lofreq,
      AltRev_Lofreq = ADR_Lofreq
  ) %>%
  group_by(key) %>%
  mutate(
    SBF_pvalue_Lofreq = fisher.test(matrix(c(RefFwd_Lofreq, RefRev_Lofreq, AltFwd_Lofreq, AltRev_Lofreq), nrow = 2), simulate.p.value = TRUE)$p.value
  ) %>%
  ungroup()

  V <- V %>% select(
      QUAL_vardict, FILTER_vardict,
      TYPE, DP, END, VD, AF, BIAS, REFBIAS, VARBIAS, PMEAN, PSTD, ReadQual,
      QSTD, SBF, ODDRATIO, MQ, SN, HIAF, ADJAF, SHIFT3, MSI, MSILEN, NM,
      LSEQ, RSEQ, GDAMP, TLAMP, NCAMP, AMPFLAG, HICNT, HICOV, SPLITREAD, SPANPAIR,
      SVTYPE, SVLEN, DUPRATE,
      PON_2AT2_percent, PON_NAT2_percent, PON_MAX_VAF,PON_FISHER,
      ends_with("_vardict"),
      key
  ) %>% rename_with(~paste0(., "_vardict"), -contains("_vardict")) %>%
  rename_with(~sub("_vardict$", "_Vardict", .x), ends_with("_vardict")) %>%
  dplyr::rename(
      key = key_Vardict,
      pon_pvalue_Vardict = PON_FISHER_Vardict,
      RefFwd_Vardict = RDF_Vardict,
      RefRev_Vardict = RDR_Vardict,
      AltFwd_Vardict = ADF_Vardict,
      AltRev_Vardict = ADR_Vardict
  ) %>%
  group_by(key) %>%
  mutate(
    SBF_pvalue_Vardict = fisher.test(matrix(c(RefFwd_Vardict, RefRev_Vardict, AltFwd_Vardict, AltRev_Vardict), nrow = 2), simulate.p.value = TRUE)$p.value
  ) %>%
  ungroup()
  
  dfs <- list(M, L, V)
  final <- reduce(dfs, full_join, by = c("key")) %>%
    left_join(A, by = c("key"))
  return(final)
}

prepare_for_model <- function(df) {
  for_model <- df %>%
  separate(MBQ_Mutect, c("MBQ_Mutect_1", "MBQ_Mutect_2"), sep = ",", extra = "merge", fill = "right") %>%
  separate(MFRL_Mutect, c("MFRL_Mutect_1", "MFRL_Mutect_2"), sep = ",", extra = "merge", fill = "right") %>%
  separate(MMQ_Mutect, c("MMQ_Mutect_1", "MMQ_Mutect_2"), sep = ",", extra = "merge", fill = "right") %>%
  dplyr::rename(
    PONRefDepth = PON_RefDepth,
    PONAltDepth = PON_AltDepth,
    FILTER_Vardict_NM5_25 = FILTER_Vardict_NM5.25
  ) %>%
  mutate(
    sourcetotalsloci = n.loci.vep,
    sourcetotalsp = n.HGVSp,
    sourcetotalsc = n.HGVSc,
    hemecosmiccount = heme_cosmic_count,
    myeloidcosmiccount = myeloid_cosmic_count,
    FPpass = ifelse(FP_filter == "PASS", TRUE, FALSE),
    FILTER_Mutect_Not_Detected = ifelse(is.na(FILTER_Mutect), TRUE, FALSE),
    FILTER_Lofreq_Not_Detected = ifelse(is.na(FILTER_Lofreq), TRUE, FALSE),
    FILTER_Vardict_Not_Detected = ifelse(is.na(FILTER_Vardict), TRUE, FALSE),
    passhomopolymerfilter = ifelse(homopolymerCase == "", 1, 0),
    maxgnomADAF = ifelse(pass_max_sub_gnomAD_AF, 1, 0),
    dustscore = dust_score,
    dustscore3 = dust_score,
    dustscore5 = dust_score,
    dustcore10 = dust_score,
    clinvarSCOREVEP = 0
  )
}

split_filters <- function(df) {
  mutect_filters <- c(
      'PASS',
      'FAIL',
      'base_qual',
      'clustered_events',
      'contamination',
      'duplicate',
      'fragment',
      'germline',
      'haplotype',
      'low_allele_frac',
      'map_qual',
      'multiallelic',
      'n_ratio',
      'normal_artifact',
      'orientation',
      'panel_of_normals',
      'position',
      'possible_numt',
      'slippage',
      'strand_bias',
      'strict_strand',
      'weak_evidence'
  )

  vardict_filters <- c(
      'PASS',
      'q22.5',
      'Q10',
      'p8',
      'SN1.5',
      'Bias',
      'pSTD',
      'd3',
      'v2',
      'min_af',
      'MSI12',
      'NM5.25',
      'InGap',
      'InIns',
      'Cluster0bp',
      'LongMSI',
      'AMPBIAS',
      'BCBIO'
  )

  lofreq_filters <- c(
      'PASS',
      'min_snvqual_45',
      'min_indelqual_31',
      'min_dp_10',
      'sb_fdr',
      'min_snvqual_62',
      'min_indelqual_49'
  )

  fp_filters <- c(
    'RLD25',
    'MVF0',
    'SB1',
    'NRC',
    'PB10',
    'IRC',
    'MMQSD50',
    'DETP20',
    'MVC4',
    'MMQS100',
    'MQD30'
  )

  # This works
  filter_dictionary <- list(
    mutect = setNames(rep(0, length(mutect_filters)), paste0('FILTER_Mutect_', mutect_filters)),
    lofreq = setNames(rep(0, length(lofreq_filters)), paste0('FILTER_Lofreq_', lofreq_filters)),
    vardict = setNames(rep(0, length(vardict_filters)), paste0('FILTER_Vardict_', vardict_filters)),
    fp = setNames(rep(0, length(fp_filters)), paste0('FP_Filter_', fp_filters))
  )

  df[, names(filter_dictionary$mutect)] <- rep(0, nrow(df))
  df[, names(filter_dictionary$lofreq)] <- rep(0, nrow(df))
  df[, names(filter_dictionary$vardict)] <- rep(0, nrow(df))
  df[, names(filter_dictionary$fp)] <- rep(0, nrow(df))

  # Set binary values for each filter in the dictionary
  df[, names(filter_dictionary$mutect)] <- lapply(names(filter_dictionary$mutect), function(x) as.integer(grepl(sub("FILTER_Mutect_", "", x), df$FILTER_Mutect)))
  df[, names(filter_dictionary$lofreq)] <- lapply(names(filter_dictionary$lofreq), function(x) as.integer(grepl(sub("FILTER_Lofreq_", "", x), df$FILTER_Lofreq)))
  df[, names(filter_dictionary$vardict)] <- lapply(names(filter_dictionary$vardict), function(x) as.integer(grepl(sub("FILTER_Vardict_", "", x), df$FILTER_Vardict)))
  df[, names(filter_dictionary$fp)] <- lapply(names(filter_dictionary$fp), function(x) as.integer(grepl(sub("FP_Filter_", "", x), df$FILTER_Mutect) | grepl(sub("FP_Filter_", "", x), df$FILTER_Lofreq) | grepl(sub("FP_Filter_", "", x), df$FILTER_Vardict)))

  # Remove the FP_Filters from FILTER_Mutect, FILTER_Lofreq, and FILTER_Vardict
  df <- df %>%
    mutate(FILTER_Mutect = gsub(paste0(fp_filters, collapse = "|"), "", FILTER_Mutect),
           FILTER_Lofreq = gsub(paste0(fp_filters, collapse = "|"), "", FILTER_Lofreq),
           FILTER_Vardict = gsub(paste0(fp_filters, collapse = "|"), "", FILTER_Vardict))

  # Check if ";" is the last character, if it is, then remove it
  df <- df %>%
    mutate(FILTER_Mutect = gsub(";+$", "", FILTER_Mutect),
           FILTER_Lofreq = gsub(";+$", "", FILTER_Lofreq),
           FILTER_Vardict = gsub(";+$", "", FILTER_Vardict))

  # Check if FILTER_Mutect == ";" and replace
  df <- df %>%
    mutate(FILTER_Mutect = if_else(FILTER_Mutect == ";", "", FILTER_Mutect),
           FILTER_Lofreq = if_else(FILTER_Lofreq == ";", "", FILTER_Lofreq),
           FILTER_Vardict = if_else(FILTER_Vardict == ";", "", FILTER_Vardict))

  # If empty string is left, replace with PASS
  df <- df %>%
    mutate(FILTER_Mutect = if_else(FILTER_Mutect == "", "PASS", FILTER_Mutect),
           FILTER_Lofreq = if_else(FILTER_Lofreq == "", "PASS", FILTER_Lofreq),
           FILTER_Vardict = if_else(FILTER_Vardict == "", "PASS", FILTER_Vardict))

  # If FILTER_Mutect is PASS, then set FILTER_Mutect_PASS to 1
  df <- df %>%
    mutate(FILTER_Mutect_PASS = if_else(FILTER_Mutect == "PASS", 1, 0),
           FILTER_Lofreq_PASS = if_else(FILTER_Lofreq == "PASS", 1, 0),
           FILTER_Vardict_PASS = if_else(FILTER_Vardict == "PASS", 1, 0))

  # Fill NA with 0 for Filter_Mutect_PASS
  df <- df %>%
    mutate(FILTER_Mutect_PASS = if_else(is.na(FILTER_Mutect_PASS), 0, FILTER_Mutect_PASS),
           FILTER_Lofreq_PASS = if_else(is.na(FILTER_Lofreq_PASS), 0, FILTER_Lofreq_PASS),
           FILTER_Vardict_PASS = if_else(is.na(FILTER_Vardict_PASS), 0, FILTER_Vardict_PASS))

  return(df)
}

# This is temporarily unused..
complex_variants <- function(full) {
  # Creates the row numbers
  M$Column1 <- 1:nrow(M)
  L$Column1 <- 1:nrow(L)
  V$Column1 <- 1:nrow(V)

  V$IsComplex <- rep(FALSE, nrow(V))
  M$IsComplex <- rep(FALSE, nrow(M))
  L$IsComplex <- rep(FALSE, nrow(L))
  V$IsComplexPart <- rep(FALSE, nrow(V))
  M$IsComplexPart <- rep(FALSE, nrow(M))
  L$IsComplexPart <- rep(FALSE, nrow(L))
  V$ComplexKey <- rep("", nrow(V))
  M$ComplexKey <- rep("", nrow(M))
  L$ComplexKey <- rep("", nrow(L))
  
  sqo_mutect <- NULL
  if (nrow(M) > 0) {
    sqo_mutect <- sqo_main(M,V,paste0("mutect_complex_",SAMPLE_NAME),directoryloc)
  }
  if (!is.null(sqo_mutect)) {
    # check if the complex should " pass ". It should be within 1 order of magnitude from the average in the complex. It should also have less than 150 bp in both REF and ALT. The complex should also have more than 1 variant making it up
    sqo_mutect$pass <- abs(log10(sqo_mutect$VAF / sapply(sqo_mutect$bestVAF,mean))) <= 1 & sapply(strsplit(sapply(strsplit(sqo_mutect$key," "),"[",3),">"),nchar) < 150 & sapply(strsplit(sapply(strsplit(sqo_mutect$key," "),"[",4),">"),nchar) < 150 & sapply(sqo_mutect$calpos_best,length) > 1
    permutation <- lapply(sqo_mutect$bestVAF,sort,decreasing=TRUE)
    chromosomes <- sapply(strsplit(sqo_mutect$key," "),"[",1)
    m_dict <- list()
    for (i in 1:nrow(sqo_mutect)) {
      for (j in 1:length(permutation[[i]])) {
        s <- paste0(chromosomes[i]," ", sqo_mutect[i,"calpos_best"][permutation[[i]][j]]," ",sqo_mutect[i,"calref_best"][permutation[[i]][j]],">",sqo_mutect[i,"calalt_best"][permutation[[i]][j]])
        m_dict[[s]] <- sqo_mutect[i,"key"]
        l <- which(L$key == s)[1]
        M[l,"IsComplexPart"] <- TRUE
        M[l,"ComplexKey"] <- sqo_mutect[i,"key"]
      }
      if (sqo_mutect[i,"pass"]) {
        V[which(V$key == sqo_mutect[i,"key"]),"IsComplex"] <- TRUE
      }
    }
  }
  
  sqo_lofreq <- NULL
  if (nrow(L) > 0) {
    sqo_lofreq <- sqo_main(L,V,paste0("lofreq_complex_",SAMPLE_NAME),directoryloc)
  }
  if (!is.null(sqo_lofreq)) {
    sqo_lofreq$pass <- abs(log10(sqo_lofreq$VAF / sapply(sqo_lofreq$bestVAF,mean))) <= 1 & sapply(strsplit(sapply(strsplit(sqo_lofreq$key," "),"[",3),">"),nchar) < 150 & sapply(strsplit(sapply(strsplit(sqo_lofreq$key," "),"[",4),">"),nchar) < 150 & sapply(sqo_lofreq$calpos_best,length) > 1
    permutation <- lapply(sqo_lofreq$bestVAF,sort,decreasing=TRUE)
    chromosomes <- sapply(strsplit(sqo_lofreq$key," "),"[",1)
    l_dict <- list()
    for (i in 1:nrow(sqo_lofreq)) {
      for (j in 1:length(permutation[[i]])) {
        s <- paste0(chromosomes[i]," ", sqo_lofreq[i,"calpos_best"][permutation[[i]][j]]," ",sqo_lofreq[i,"calref_best"][permutation[[i]][j]],">",sqo_lofreq[i,"calalt_best"][permutation[[i]][j]])
        l_dict[[s]] <- sqo_lofreq[i,"key"]
        l <- which(L$key == s)[1]
        L[l,"IsComplexPart"] <- TRUE
        L[l,"ComplexKey"] <- sqo_lofreq[i,"key"]
      }
      if (sqo_lofreq[i,"pass"]) {
        V[which(V$key == sqo_lofreq[i,"key"]),"IsComplex"] <- TRUE
      }
    }
  }
  
  return(list(L,M,V))
}