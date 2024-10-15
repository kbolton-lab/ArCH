getVAFs <- function(x) {
  mutect_VAF <- grep("Mutect2_gt_AF", colnames(x))
  vardict_VAF <- grep("Vardict_gt_AF", colnames(x))
  VAFS <- c(mutect_VAF, vardict_VAF)
  VAFS
}

readBed <- function(filePath) {
  tmp <- fread(filePath, data.table = FALSE)
  gr <- GRanges(
    seqnames = tmp[, 1],
    ranges = IRanges(
      start = tmp[, 2] + 1, # rtracklayer import start as +1
      end = tmp[, 3]
    )
  )
  if (ncol(tmp) == 4) {
    gr$socre <- tmp[, 4]
  } else if (ncol(tmp) == 5) {
    gr$name <- tmp[, 4]
    gr$socre <- tmp[, 5]
  } else if (ncol(tmp) == 6) {
    gr$name <- tmp[, 4]
    gr$socre <- tmp[, 5]
    strand(gr) <- tmp[, 6]
  }
  gr
}

sameLetters <- function(x) {
  letter <- substr(x, 1, 1)
  for (pos in 1:nchar(x)) {
    if (substr(x, pos, pos) != letter) {
      return(FALSE)
    }
  }
  return(TRUE)
}

saveRDS.gz <- function(object, file) {
  con <- pipe(paste0("gzip -1 > ", file), "wb")
  saveRDS(object, file = con)
  close(con)
}

readVcfFiles <- function(mutect2_vcf_file, vardict_vcf_file, sample_id) {
  ### read in
  mutect <- vcfR2tidy(read.vcfR(mutect2_vcf_file, verbose = FALSE),
    single_frame = TRUE,
    info_types = TRUE,
    format_types = TRUE
  )
  mutect.info <- mutect$meta[mutect$meta$Tag == "INFO", ]$ID
  mutect.info <- mutect.info[!(mutect.info %in% c("SAMPLE", "PON_RefDepth", "PON_AltDepth", "PON_FISHER", "CSQ"))]
  mutect.df <- mutect[["dat"]]

  vardict <- vcfR2tidy(read.vcfR(vardict_vcf_file, verbose = FALSE),
    single_frame = TRUE,
    info_types = TRUE,
    format_types = TRUE
  )
  vardict.info <- vardict$meta[vardict$meta$Tag == "INFO", ]$ID
  vardict.info <- vardict.info[!(vardict.info %in% c("OLD_MULTIALLELIC", "SAMPLE", "CSQ"))]
  vardict.df <- vardict[["dat"]]
  vardict.df <- vardict.df[, !colnames(vardict.df) %in% c("OLD_MULTIALLELIC")]

  mutect.df <- mutect.df %>% distinct()
  vardict.df <- vardict.df %>% distinct()
  
  # Fix different CSQ length, some vardict file only has 85 fields
  # mutectCSQ <- mutect$meta$Description[mutect$meta$ID == "CSQ"] %>% 
  #   str_split(": ") %>% .[[1]] %>% .[2]
  # vardictCSQ <- vardict$meta$Description[vardict$meta$ID == "CSQ"] %>% 
  #   str_split(": ") %>% .[[1]] %>% .[2]
  # # If 2 string are differnt, just set vardict CSQ to NA
  # if (mutectCSQ != vardictCSQ) {
  #   vardict.df$CSQ <- NA
  # }

  # Prefix all FORMAT/ (gt_) columns with caller such as Vardict_gt_AD_ref & Vardict_gt_AD_alt
  colnames(mutect.df)[which(grepl("^gt_", colnames(mutect.df)))] <-
    paste0("Mutect2_", colnames(mutect.df)[which(grepl("^gt_", colnames(mutect.df)))])
  # then move all FILTER & FORMAT columns to front after ALT/FILTER column
  mutect2.pfx <- which(grepl("Mutect2_", colnames(mutect.df)))
  mutect.new.col.ord <- c(1:7, mutect2.pfx, 8:(min(mutect2.pfx) - 1))
  length(mutect.new.col.ord) == dim(mutect.df)[2]
  mutect2.final <- mutect.df[, mutect.new.col.ord]
  colnames(mutect2.final)[colnames(mutect2.final) == "FILTER"] <- "Mutect2_FILTER"
  mutect2.final$Mutect2_PASS <- ifelse(grepl("PASS", mutect2.final$Mutect2_FILTER), 1, 0) ## 56 cols

  ## prefix FMT data with Caller
  colnames(vardict.df)[which(grepl("^gt_", colnames(vardict.df)))] <-
    paste0("Vardict_", colnames(vardict.df)[which(grepl("^gt_", colnames(vardict.df)))])
  vardict.pfx <- which(grepl("Vardict_", colnames(vardict.df)))
  vardict.new.col.ord <- c(1:7, vardict.pfx, 8:(min(vardict.pfx) - 1))
  length(vardict.new.col.ord) == dim(vardict.df)[2]
  vardict.final <- vardict.df[, vardict.new.col.ord]
  colnames(vardict.final)[colnames(vardict.final) == "FILTER"] <- "Vardict_FILTER"
  vardict.final$Vardict_PASS <- ifelse(grepl("PASS", vardict.final$Vardict_FILTER), 1, 0) ## 52 cols

  ## first fill PON_2AT2 percent column, need fix for batch 12
  if (!"PON_2AT2_percent" %in% colnames(mutect2.final)) {
    mapKey <- unite(mutect2.final, "CHROM", "POS", "REF", "ALT")[, 1]
    mutect2.final$PON_2AT2_percent <- as.numeric(mapKey %in% mutectPonKey)
  }
  if (!"PON_2AT2_percent" %in% colnames(vardict.final)) {
    mapKey <- unite(vardict.final, "CHROM", "POS", "REF", "ALT")[, 1]
    vardict.final$PON_2AT2_percent <- as.numeric(mapKey %in% vardictPonKey)
  }

  mutect2.final$PON_2AT2_percent <- tidyr::replace_na(mutect2.final$PON_2AT2_percent, 0)
  colnames(mutect2.final)[colnames(mutect2.final) == "PON_2AT2_percent"] <- "Mutect2_PON_2AT2_percent"

  vardict.final$PON_2AT2_percent <- tidyr::replace_na(vardict.final$PON_2AT2_percent, 0)
  colnames(vardict.final)[colnames(vardict.final) == "PON_2AT2_percent"] <- "Vardict_PON_2AT2_percent"

  ## remove columns where all NA
  na.mutect.cols <- names(mutect2.final[, apply(mutect2.final, 2, function(x) all(is.na(x)))])
  na.vardict.cols <- names(vardict.final[, apply(vardict.final, 2, function(x) all(is.na(x)))])

  length(na.mutect.cols)
  length(na.vardict.cols)
  ## remove the columns that are from Alex's script
  na.vardict.cols <- union(na.vardict.cols, c("calpos", "calref", "calalt", "calpos_best", "calref_best", "calalt_best"))
  length(na.vardict.cols)
  mutect2.final <- mutect2.final[, -which(names(mutect2.final) %in% c(na.mutect.cols))]
  vardict.final <- vardict.final[, -which(names(vardict.final) %in% c(na.vardict.cols))]

  ## prefix INFO data with Caller
  colnames(mutect2.final)[colnames(mutect2.final) %in% mutect.info] <-
    paste0("Mutect2_", colnames(mutect2.final)[colnames(mutect2.final) %in% mutect.info])

  colnames(vardict.final)[colnames(vardict.final) %in% vardict.info] <-
    paste0("Vardict_", colnames(vardict.final)[colnames(vardict.final) %in% vardict.info])

  ## Way 3 works
  # https://stackoverflow.com/questions/16042380/merge-data-frames-and-overwrite-values
  colnames(mutect2.final)[which(colnames(mutect2.final) == "CALLER")] <- "Mutect2_CALLER"
  colnames(vardict.final)[which(colnames(vardict.final) == "CALLER")] <- "Vardict_CALLER"
  mutect2.final$Mutect2_CALLER <- 1
  vardict.final$Vardict_CALLER <- 1
  mutect2.final$SAMPLE <- sample_id
  vardict.final$SAMPLE <- sample_id

  mutect2.final <- as.data.frame(mutect2.final)
  vardict.final <- as.data.frame(vardict.final)
  

  ## vardict with mutect2
  intersection <- intersect(colnames(vardict.final), colnames(mutect2.final))
  intersection <- intersection[6:length(intersection)] ## first 5 are the grouped by columns
  intersection.cols.x <- paste0(intersection, ".x")
  intersection.cols.y <- paste0(intersection, ".y")
  final <- full_join(vardict.final, mutect2.final,
    by = c(
      "CHROM" = "CHROM", "POS" = "POS",
      "REF" = "REF", "ALT" = "ALT", "SAMPLE" = "SAMPLE"
    )
  )
  final[, intersection.cols.x][is.na(final[, intersection.cols.x])] <- final[, intersection.cols.y][is.na(final[, intersection.cols.x])]

  ## remove *.x from column names
  colnames(final)[colnames(final) %in% intersection.cols.x] <- intersection
  ## remove *.y columns completely
  final <- final[, !(colnames(final) %in% intersection.cols.y)]

  # ## remove blank CSQ columns because they are non-coding (off-target) and separate
  # final <- final[!is.na(final$CSQ), ]
  # # VEP CSQ
  # ## new VEP with gnomADg has 106 fields
  # CSQnames <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
  # CSQnames <- str_split(CSQnames, "\\|")[[1]]
  # 
  # # CSQ has 91 columns, maybe we should suffix cols with VEP?
  # final <- final %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep = "\\|", extra = "merge", fill = "right")
  
  # Remove CSQ column
  final <- final %>% dplyr::select(-CSQ) 

  final$PON_RefDepth <- as.numeric(final$PON_RefDepth)
  final$PON_AltDepth <- as.numeric(final$PON_AltDepth)
  final$PON_FISHER <- as.numeric(final$PON_FISHER)
  ## fill NA with 0
  for (col in c("Mutect2_CALLER", "Vardict_CALLER", "Mutect2_PASS", "Vardict_PASS")) {
    final[, col] <- tidyr::replace_na(final[, col], 0)
  }
  final <- final %>%
    # Keep var called by both
    dplyr::filter(Mutect2_CALLER == 1 & Vardict_CALLER == 1)

  final
}

annotateGnomad <- function(df) {
  message("annotating variants with gnomAD...")
  
  # Fix gnomADe problem => not necessary for new VEP
  # for (col in grep("^gnomADe_.*_VEP", colnames(df), value = TRUE)) {
  #   tmp <- str_split(df[, col], '&|,', simplify = TRUE) 
  #   mode(tmp) <- 'numeric'
  #   tmp <-  tmp %>% 
  #     tidyr::replace_na(0) %>%
  #     matrixStats::rowMaxs(na.rm = TRUE)
  #   df[, col] <- tmp
  # }

  for (col in colnames(df)[grepl("^gnomAD.*_.*_VEP", colnames(df))]) {
    df[, col] <- tidyr::replace_na(as.numeric(df[, col]), 0)
  }

  final.coding.gnomad.sorted <- df[with(df, order(CHROM, POS)), ]
  # gnomAD, new VEP has gnomADe and g
  gnomADe.col <- grep("gnomADe_.*_VEP", colnames(final.coding.gnomad.sorted))
  final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP <- matrixStats::rowMaxs(final.coding.gnomad.sorted[, gnomADe.col] %>% 
                                                                          as.matrix(), na.rm = TRUE)

  gnomADg.col <- grep("gnomADg_.*_VEP", colnames(final.coding.gnomad.sorted))
  final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP <- matrixStats::rowMaxs(final.coding.gnomad.sorted[, gnomADg.col] %>% 
                                                                          as.matrix(), na.rm = TRUE)

  final.coding.gnomad.sorted
}

annotateGenomicRegion <- function(df, supportData) {
  df$key <- with(df, paste(CHROM, POS, REF, ALT, sep = ":"))
  final <- df %>% distinct(key, .keep_all = TRUE)

  blacklist <- supportData[["blacklist"]]
  dups <- supportData[["dups"]]
  repeats <- supportData[["repeats"]]
  repeatMasker <- supportData[["repeatMasker"]]

  blacklist$blacklist <- TRUE
  dups$dups <- TRUE
  repeats$repeats <- TRUE
  repeatMasker$repeatMasker <- TRUE
  sample.gr <- GRanges(
    seqnames = final$CHROM,
    ranges = IRanges(
      start = final$POS,
      end = final$POS
    ),
    mcols = final[, 3:ncol(final)]
  )
  join <- plyranges::join_overlap_left(sample.gr, blacklist)
  join2 <- plyranges::join_overlap_left(join, dups)
  join3 <- plyranges::join_overlap_left(join2, repeats)
  join4 <- plyranges::join_overlap_left(join3, repeatMasker)
  columns <- colnames(final)
  final <- as.data.frame(join4)[, c(1, 2, 6:ncol(as.data.frame(join4)))]
  colnames(final) <- c(
    columns, "blacklist", "dup_score", "dups",
    "trf", "trf_score", "simpleTandemRepeat",
    "RMName", "RMScore", "repeatMasker"
  )
  final <- final %>% distinct(CHROM, POS, REF, ALT, SAMPLE, .keep_all = TRUE)

  df <- df %>%
    left_join(final %>% dplyr::select(key, blacklist, dups, simpleTandemRepeat, repeatMasker),
      by = "key"
    ) %>%
    dplyr::select(-key)

  df$blacklist[is.na(df$blacklist)] <- FALSE
  df$dups[is.na(df$dups)] <- FALSE
  df$simpleTandemRepeat[is.na(df$simpleTandemRepeat)] <- FALSE
  df$repeatMasker[is.na(df$repeatMasker)] <- FALSE

  df
}

annotateComplexity <- function(df) {
  final.passed <- df %>% distinct(key, .keep_all = TRUE)
  final.passed$context_5 <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(
    seqnames = final.passed$CHROM,
    ranges = IRanges(start = final.passed$POS - 10, end = final.passed$POS - 1)
  ))
  final.passed$context_3 <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(
    seqnames = final.passed$CHROM,
    ranges = IRanges(start = final.passed$POS + nchar(final.passed$REF), end = final.passed$POS + nchar(final.passed$REF) + 9)
  ))
  final.passed$context_5_3 <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(
    seqnames = final.passed$CHROM,
    ranges = IRanges(start = final.passed$POS - 5, end = final.passed$POS + nchar(final.passed$REF) + 4)
  ))

  pdf(NULL) # Avoid creating pdf file
  final.passed$dust_score_5 <- R453Plus1Toolbox::complexity.dust(final.passed$context_5)
  final.passed$dust_score_3 <- R453Plus1Toolbox::complexity.dust(final.passed$context_3)
  final.passed$dust_score_5_3 <- R453Plus1Toolbox::complexity.dust(final.passed$context_5_3)
  dev.off()
  final.passed$dust_score <- apply(final.passed[, c("dust_score_5", "dust_score_3", "dust_score_5_3")], 1, max)
  
  final.passed$context_5 <- as.character(final.passed$context_5)
  final.passed$context_3 <- as.character(final.passed$context_3)
  final.passed$context_5_3 <- as.character(final.passed$context_5_3)

  final.passed <- final.passed %>%
    dplyr::mutate(
      # Determines if the Downstream sequence has 3 of the same nucleotides in a row
      case_NXXX = str_detect(context_3, regex("^[A]{3,}|^[G]{3,}|^[T]{3,}|^[C]{3,}")) &
        # Checks if the last character of the mutation completes the homopolymer
        substr(ALT, nchar(ALT), nchar(ALT)) == substr(context_3, 1, 1),
      # Sandwich Mutation
      case_XNXX = ifelse(nchar(ALT) == 1, TRUE, FALSE) & # Can only occur with Single Nucleotides
        # Check if the last character of the upstream + the mutation + the first two characters of downstream
        # completes the homopolymer
        paste0(
          substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
          ALT,
          substr(context_3, start = 1, stop = 2)
        ) %>%
          str_detect(regex("[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}")),
      # Second Sandwich Case
      case_XXNX = ifelse(nchar(ALT) == 1, TRUE, FALSE) &
        paste0(
          substr(context_5, start = nchar(context_5) - 1, stop = nchar(context_5)),
          ALT,
          substr(context_3, start = 1, stop = 1)
        ) %>%
          str_detect(regex("[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}")),
      # Determines if the Upstream sequence ends with 3 of the same nucleotides in a row
      case_XXXN = sapply(substr(context_5, start = nchar(context_5) - 2, stop = nchar(context_5)), sameLetters) &
        # Checks if the first character of the mutation completes the homopolymer
        substr(ALT, 1, 1) == substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
      # Same as Case NXXX but for dinucleotides
      case_NNXX = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & str_detect(context_3, regex("^[A]{2,}|^[G]{2,}|^[T]{2,}|^[C]{2,}")) &
        substr(ALT, nchar(ALT) - 1, nchar(ALT)) == substr(context_3, 1, 2),
      # Sandiwch Case
      case_XNNX = ifelse(nchar(ALT) == 2, TRUE, FALSE) & sapply(ALT, sameLetters) &
        paste0(
          substr(context_5, start = nchar(context_5) - 1, stop = nchar(context_5)),
          ALT,
          substr(context_3, start = 1, stop = 2)
        ) %>%
          str_detect(regex("[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}")),
      # Same as Case XXXN but for dinucleotides
      case_XXNN = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & sapply(substr(context_5, start = nchar(context_5) - 1, stop = nchar(context_5)), sameLetters) &
        substr(ALT, 1, 2) == substr(context_5, start = nchar(context_5) - 1, stop = nchar(context_5))
    )

  caseList <- c("case_NXXX", "case_XNXX", "case_XXNX", "case_XXXN", "case_NNXX", "case_XNNX", "case_XXNN")
  names(caseList) <- sapply(caseList, function(x) str_split(x, "_")[[1]][2])
  final.passed$homopolymerCase <- apply(final.passed[, caseList], 1, function(x) {
    paste0(names(caseList)[x], collapse = ";")
  })

  df <- df %>%
    left_join(final.passed %>% dplyr::select(key, context_5, context_3, context_5_3, dust_score, homopolymerCase),
      by = "key"
    )

  df
}

prepareAnnotatePdData <- function(df, supportData) {
  final <- df

  AminoAcids <- c(
    "Cys" = "C", "Asp" = "D", "Ser" = "S", "Gln" = "Q", "Lys" = "K",
    "Ile" = "I", "Pro" = "P", "Thr" = "T", "Phe" = "F", "Asn" = "N",
    "Gly" = "G", "His" = "H", "Leu" = "L", "Arg" = "R", "Trp" = "W",
    "Ala" = "A", "Val" = "V", "Glu" = "E", "Tyr" = "Y", "Met" = "M",
    "%3D" = "=", "=" = "="
  )

  final$AAchange <- gsub("(.*p\\.)(.*)", "\\2", final$HGVSp_VEP)
  for (i in 1:length(AminoAcids)) {
    final$AAchange <- gsub(names(AminoAcids)[i], AminoAcids[i], final$AAchange)
  }
  final$gene_loci_p <- paste(final$SYMBOL_VEP,
    paste0(
      sapply(final$AAchange, function(x) str_split(x, "[0-9]+", n = 2)[[1]][1]),
      as.numeric(str_extract(final$AAchange, "\\d+"))
    ),
    sep = "_"
  )
  final$gene_loci_c <- paste(final$SYMBOL_VEP,
    gsub(".*:", "", final$HGVSc),
    sep = "_"
  )
  final$gene_loci_vep <- ifelse(is.na(final$gene_loci_p), final$gene_loci_c, final$gene_loci_p)
  final$key <- with(final, paste(CHROM, POS, REF, ALT, sep = ":"))
  final$gene_aachange <- with(final, paste(SYMBOL_VEP, AAchange, sep = "_"))
  final$gene_cDNAchange <- paste(final$SYMBOL_VEP, gsub(".*:", "", final$HGVSc_VEP), sep = "_")

  ## ch_pd stuff
  vars <- supportData[["vars"]]
  vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
  vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:", "", vars$HGVSc_VEP), sep = "_")

  dims <- dim(final)[[1]]
  finalTmp <- final %>%
    dplyr::filter(key %in% vars$key |
      gene_loci_vep %in% vars$gene_loci_vep |
      gene_aachange %in% vars$gene_aachange |
      gene_cDNAchange %in% vars$gene_cDNAchange) %>%
    dplyr::select(key, gene_loci_vep, gene_aachange, gene_cDNAchange)

  finalTmp$truncating <- "not"
  finalTmp <- finalTmp %>%
    mutate(
      truncating = ifelse(str_detect(gene_aachange, "Ter"), "truncating", truncating)
    )
  
  tmp <- vars[vars$key %in% finalTmp$key |
    vars$gene_loci_vep %in% finalTmp$gene_loci_vep, ]
  finalTmp <- sqldf("SELECT l.*, r.`n.loci.vep`, r.`source.totals.loci`
            FROM `finalTmp` as l
            LEFT JOIN `tmp` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
  finalTmp <- finalTmp %>% distinct()

  finalTmp <- sqldf("SELECT l.*, r.`n.loci.truncating.vep`, r.`source.totals.loci.truncating`
            FROM `finalTmp` as l
            LEFT JOIN `tmp` as r
            on (l.key = r.key AND l.truncating = r.truncating) OR (l.gene_loci_vep = r.gene_loci_vep AND l.truncating = r.truncating)")
  finalTmp <- finalTmp %>% distinct()
  finalTmp <- finalTmp %>% dplyr::select(-truncating)

  ## make sure aachange exists as in doesn't end with an '_'; example: DNMT3A_ for splice
  vars <- vars[!(grepl("_$", vars$gene_aachange) | grepl("_$", vars$gene_cDNAchange)), ]
  tmp <- vars[vars$key %in% finalTmp$key |
    vars$gene_aachange %in% finalTmp$gene_aachange, ]
  finalTmp <- sqldf("SELECT l.*, r.`n.HGVSp`, r.`source.totals.p`
            FROM `finalTmp` as l
            LEFT JOIN `tmp` as r
            on l.key = r.key OR (l.gene_aachange = r.gene_aachange)")
  finalTmp <- finalTmp %>% distinct()

  tmp <- vars[vars$key %in% finalTmp$key |
    vars$gene_cDNAchange %in% finalTmp$gene_cDNAchange, ]
  finalTmp <- sqldf("SELECT l.*, r.`n.HGVSc`, r.`source.totals.c`
            FROM `finalTmp` as l
            LEFT JOIN `tmp` as r
            on l.key = r.key OR (l.gene_cDNAchange = r.gene_cDNAchange)")
  finalTmp <- finalTmp %>% distinct()

  final <- final %>%
    left_join(finalTmp %>% dplyr::select(-gene_loci_vep, -gene_aachange, -gene_cDNAchange),
      by = "key"
    )

  paste0("dims match after sqldf: ", dim(final)[[1]] == dims)

  final$Mutect2_PON_2AT2_percent <- tidyr::replace_na(final$Mutect2_PON_2AT2_percent, 0)
  final$Vardict_PON_2AT2_percent <- tidyr::replace_na(final$Vardict_PON_2AT2_percent, 0)

  final
}

annotateCosmicNoSplit <- function(df, var_key_column, supportData) {
  dfTmp <- df

  dfTmp$HGVSp_VEP <- gsub(".*:", "", dfTmp$HGVSp_VEP)
  dfTmp$Gene_HGVSp_VEP <- with(dfTmp, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
  dfTmp$skey <- paste(dfTmp[[var_key_column]], sep = ":")

  cosmic <- supportData[["fullCosmic"]]

  # Reduce run time by optimizing sql calling
  tmp <- cosmic[cosmic$var_key %in% dfTmp$var_key | cosmic$Gene_HGVSp_VEP %in% dfTmp$Gene_HGVSp_VEP, ]
  tmp <- sqldf("SELECT l.*, r.COSMIC_ID, r.CosmicCount,
              r.heme_cosmic_count, r.myeloid_cosmic_count
              FROM `dfTmp` as l
              LEFT JOIN `tmp` as r
              on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
  tmp <- tmp[, c(colnames(dfTmp), "COSMIC_ID", "CosmicCount", "heme_cosmic_count", "myeloid_cosmic_count")]
  tmp$CosmicCount <- tidyr::replace_na(tmp$CosmicCount, 0)
  tmp$heme_cosmic_count <- tidyr::replace_na(tmp$heme_cosmic_count, 0)
  tmp$myeloid_cosmic_count <- tidyr::replace_na(tmp$myeloid_cosmic_count, 0)
  tmp <- tmp %>%
    group_by(var_key) %>%
    slice_max(order_by = heme_cosmic_count, n = 1, with_ties = F)
  tmp <- data.frame(tmp)
  cosmic.test <- tmp

  cosmic.test$skey <- paste(cosmic.test[[var_key_column]], sep = ":")
  row.names(cosmic.test) <- cosmic.test$skey

  ## put back into row order of dfTmp
  cosmic.test <- cosmic.test[dfTmp$skey, ]

  if (all(cosmic.test[, c("CHROM", "POS", "REF", "ALT")] == df[, c("CHROM", "POS", "REF", "ALT")])) {
    MUTS <- cbind(df, cosmic.test %>% dplyr::select(COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count))
    return(MUTS)
  } else {
    print("COSMIC bad #################################################")
  }
}

annotateCosmicJoin <- function(df, supportData) {
  dfTmp <- df

  dfTmp$HGVSp_VEP <- gsub(".*:", "", dfTmp$HGVSp_VEP)
  dfTmp$Gene_HGVSp_VEP <- with(dfTmp, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))

  cosmic <- supportData[["fullCosmic"]]

  # Reduce run time by optimizing sql calling
  tmp <- cosmic[cosmic$var_key %in% dfTmp$var_key | cosmic$Gene_HGVSp_VEP %in% dfTmp$Gene_HGVSp_VEP, ]

  matchKey <- dfTmp %>%
    dplyr::select(var_key) %>%
    right_join(tmp, by = "var_key") %>%
    dplyr::select(var_key, COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count) %>%
    drop_na(COSMIC_ID)

  matchHgvsp <- dfTmp %>%
    dplyr::select(Gene_HGVSp_VEP, var_key) %>%
    right_join(tmp %>% dplyr::select(-var_key), by = "Gene_HGVSp_VEP") %>%
    dplyr::select(var_key, COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count) %>%
    drop_na(COSMIC_ID)

  tmp <- rbind(matchKey, matchHgvsp) %>%
    mutate(
      CosmicCount = tidyr::replace_na(CosmicCount, 0),
      heme_cosmic_count = tidyr::replace_na(heme_cosmic_count, 0),
      myeloid_cosmic_count = tidyr::replace_na(myeloid_cosmic_count, 0)
    )
  tmp <- tmp %>%
    group_by(var_key) %>%
    arrange(desc(heme_cosmic_count), desc(CosmicCount), desc(myeloid_cosmic_count), .by_group = TRUE) %>%
    slice_head(n = 1)

  df <- df %>%
    left_join(tmp, by = "var_key") %>%
    mutate(
      CosmicCount = tidyr::replace_na(CosmicCount, 0),
      heme_cosmic_count = tidyr::replace_na(heme_cosmic_count, 0),
      myeloid_cosmic_count = tidyr::replace_na(myeloid_cosmic_count, 0)
    )
  df
}

annotateOncoKb <- function(MUTS, supportData) {

  # Rename and format
  {
    ## has 9
    ## missing 6
    # additional we can have "upstream_gene_variant", "non_coding_transcript_variant", "non_coding_transcript_exon_variant"
    # total 18 unique(unlist(sapply(MUTS$Consequence_VEP, function(x) { str_split(x,"&",simplify = TRUE)})))
    translate_consequence <- c(
      "frameshift_variant" = "frameshift_variant", # 1

      "inframe_deletion" = "inframe_deletion", # 2
      "inframe_insertion" = "inframe_insertion", # 3
      "synonymous_variant" = "synonymous_variant", # 4
      "missense_variant" = "missense_variant", # 5
      "intron_variant" = "intron_variant", # 6

      "splice_region_variant" = "splice_region_variant", # 7
      "splice_donor_variant" = "splice_region_variant", # missing
      "splice_acceptor_variant" = "splice_region_variant", # missing
      "splice_polypyrimidine_tract_variant" = "splice_region_variant", # new vep 109
      "splice_donor_region_variant" = "splice_region_variant",  # new vep 109
      "splice_donor_5th_base_variant" = "splice_region_variant",  # new vep 109
      "feature_truncation" = "feature_truncation",
      "3_prime_UTR_variant" = "3_prime_UTR_variant", # missing
      "5_prime_UTR_variant" = "5_prime_UTR_variant", # missing

      "stop_gained" = "stop_gained", # 8
      "start_lost" = "start_lost", # 9
      "stop_lost" = "missense_variant", # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
      "stop_retained_variant" = "synonymous_variant" # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
    )

    MUTS <- MUTS %>%
      mutate(VariantClass = str_split(Consequence_VEP, "&|,", simplify = TRUE)[, 1]) %>%
      mutate(VariantClass = case_when(
        (VariantClass == "coding_sequence_variant" | VariantClass == "protein_altering_variant") & (nchar(REF) - nchar(ALT)) %% 3 == 0 & nchar(REF) - nchar(ALT) > 0 ~ "inframe_deletion",
        (VariantClass == "coding_sequence_variant" | VariantClass == "protein_altering_variant") & (nchar(REF) - nchar(ALT)) %% 3 == 0 & nchar(REF) - nchar(ALT) < 0 ~ "inframe_insertion",
        (VariantClass == "coding_sequence_variant" | VariantClass == "protein_altering_variant") & (nchar(REF) - nchar(ALT)) %% 3 != 0 ~ "frameshift_variant",
        (VariantClass == "coding_sequence_variant" | VariantClass == "protein_altering_variant") & nchar(REF) == nchar(ALT) ~ "missense_variant",
        TRUE ~ VariantClass
      )) %>%
      mutate(VariantClass2 = translate_consequence[VariantClass])

    MUTS$aa_ref <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n = 2)[[1]][1])
    MUTS$aa_alt <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n = 2)[[1]][2])
    MUTS$aa_pos <- as.numeric(str_extract(MUTS$AAchange, "\\d+"))
    MUTS$var_key <- paste(MUTS$CHROM, MUTS$POS, MUTS$REF, MUTS$ALT, sep = ":")
  }

  allJsonQuery <- apply(MUTS, 1, function(x) {
    list(
      alteration = x["AAchange"],
      consequence = x["VariantClass2"],
      gene = list(
        hugoSymbol = x["SYMBOL_VEP"]
      )
    )
  })

  jsonGene <- sapply(allJsonQuery, function(x) x$gene$hugoSymbol)
  validJson <- which(jsonGene %in% supportData[["oncoKbAvailGene"]]) # Gene not in oncoKB will return Unknown
  allJsonQuery <- allJsonQuery[validJson]

  # Try redis cache first
  redisHost <- paste0("compute1-exec-", supportData[["redisHost"]], ".ris.wustl.edu")
  redisConnect <- NULL
  try({
    redisConnect <- redux::hiredis(url = redisHost, port = 8101)
  })
  if (!is.null(redisConnect)) {
    jsonKey <- sapply(allJsonQuery, function(x) paste0(unlist(x[sort(names(x))]), collapse = "_")) # reusable key even if the json structure changes a bit
    cmds <- lapply(jsonKey, function(k) redux::redis$GET(k))
    names(cmds) <- jsonKey
    redisResp <- redisConnect$pipeline(.commands = cmds)
    allJsonQuery <- allJsonQuery[sapply(redisResp, is.null)]
  }

  # If there are unavailable keys
  if (length(allJsonQuery) > 0) {
    # Query the rest from oncokb, Retry if there is no response
    resp <- NULL
    attempt <- 1
    while ((is.null(resp) || httr::status_code(resp) != 200) && attempt <= 10) {
      attempt <- attempt + 1
      try(
        resp <- httr::POST(
          url = "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange",
          httr::add_headers(c("Content-Type" = "application/json", "Authorization" = paste0("Bearer ", supportData[["oncoKbApiKey"]]))),
          body = jsonlite::toJSON(unname(allJsonQuery), auto_unbox = TRUE)
        )
      )
      if (is.null(resp) || httr::status_code(resp) != 200) {
        Sys.sleep(runif(1, min = 0.5, max = 1))
      }
    }
    respContent <- httr::content(resp)
    oncogenic <- sapply(respContent, "[[", "oncogenic")
    reviewed <- sapply(respContent, "[[", "variantSummary")

    # update redis cache
    if (!is.null(redisConnect)) {
      jsonKeyToUpdate <- jsonKey[sapply(redisResp, is.null)]
      redisResp[sapply(redisResp, is.null)] <- respContent
      names(redisResp) <- jsonKey

      cmds <- lapply(jsonKeyToUpdate, function(k) redux::redis$SET(k, redisResp[[k]]))
      redisConnect$pipeline(.commands = cmds)

      respContent <- unlist(redisResp)
    }
  } else { # if not then just return redis result
    respContent <- unlist(redisResp)
  }

  MUTS$oncoKB <- "Unknown" # set as Unknown
  MUTS$oncoKB[validJson] <- oncogenic
  MUTS$oncoKB_reviewed <- !grepl("mutation has not specifically been reviewed", reviewed)

  MUTS
}


annotatePD <- function(df, supportData) {
  # only say ch_pd if ch_pd and in gene_list
  gene_list <- supportData[["gene_list"]]

  df <- df %>% distinct(key, .keep_all = TRUE)

  MUTS <- df[, c(
    "CHROM", "POS", "REF", "ALT", "SYMBOL_VEP",
    "HGVSp_VEP", "n.HGVSp", "HGVSc_VEP", "n.HGVSc",
    "Consequence_VEP", "EXON_VEP", "AAchange"
  )]


  message("annotating variants with oncoKB...")

  MUTS <- annotateOncoKb(MUTS, supportData)
  print(table(MUTS$oncoKB))

  message("oncoKB done\n")

  ########### annotate with COSMIC########
  message("annotating variants with COSMIC...")
  # MUTS <- annotateCosmicNoSplit(MUTS, "var_key", supportData)
  MUTS <- annotateCosmicJoin(MUTS, supportData)
  MUTS$CosmicCount <- tidyr::replace_na(MUTS$CosmicCount, 0)
  MUTS$heme_cosmic_count <- tidyr::replace_na(MUTS$heme_cosmic_count, 0)
  MUTS$myeloid_cosmic_count <- tidyr::replace_na(MUTS$myeloid_cosmic_count, 0)
  message("COSMIC done\n")


  # Variant Annotation
  # Variants were annotated according to evidence for functional relevance in cancer (putative driver or CH-PD).
  MUTS$aa_pos <- as.character(MUTS$aa_pos)
  MUTS$Exon <- sapply(MUTS$EXON_VEP, function(x) str_split(x, "/")[[1]][1])
  colnames(MUTS)[colnames(MUTS) == "SYMBOL_VEP"] <- c("Gene")

  ## anything truncating in PD_table
  truncating <- c("frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "stop_gained")
  MUTS$n.HGVSp <- tidyr::replace_na(MUTS$n.HGVSp, 0)
  MUTS$n.HGVSc <- tidyr::replace_na(MUTS$n.HGVSc, 0)

  ## initial ch_pd
  MUTS$ch_pd <- 0
  MUTS$WHY_CH <- ""

  # 1. Any variant occurring in the COSMIC “haematopoietic and lymphoid” category greater than or equal to 10 times or myeloid >= 5
  MUTS <- MUTS %>%
    mutate(
      ch_pd = ifelse((heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5) &
        Gene %in% gene_list$Gene, 1, ch_pd),
      WHY_CH = ifelse((heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5) &
        Gene %in% gene_list$Gene,
      paste0(WHY_CH, "Cosmic_heme;"),
      WHY_CH
      )
    )

  # 2. Any variant noted as oncogenic or likely oncogenic in OncoKB
  MUTS <- MUTS %>%
    mutate(
      ch_pd = ifelse((oncoKB == "Oncogenic" | oncoKB == "Likely Oncogenic" | oncoKB == "Predicted Oncogenic") &
        Gene %in% gene_list$Gene,
      1,
      ch_pd
      ),
      WHY_CH = ifelse((oncoKB == "Oncogenic" | oncoKB == "Likely Oncogenic" | oncoKB == "Predicted Oncogenic") &
        Gene %in% gene_list$Gene,
      paste0(WHY_CH, " OncoKB;"),
      WHY_CH
      )
    )

  MUTS$isOncogenic <- MUTS$oncoKB == "Oncogenic" | MUTS$oncoKB == "Likely Oncogenic" | MUTS$oncoKB == "Predicted Oncogenic"
  MUTS$isTSG <- MUTS$Gene %in% gene_list$Gene[gene_list$isTSG == 1]

  # 3. Any truncating mutations (nonsense, essential splice site or frameshift indel) in known tumor suppressor genes as per the Cancer Gene Census or OncoKB.
  MUTS <- MUTS %>%
    mutate(
      ch_pd = ifelse(isTSG & VariantClass %in% truncating,
        1,
        ch_pd
      ),
      WHY_CH = ifelse(isTSG & VariantClass %in% truncating,
        paste0(WHY_CH, " TSG & Truncating;"),
        WHY_CH
      )
    )

  vars.truncating <- supportData[["vars.truncating"]]
  truncating_genes <- unique(vars.truncating[vars.truncating$n >= 10, "SYMBOL_VEP"])
  MUTS$isTruncatingHotSpot <- MUTS$Gene %in% truncating_genes & MUTS$VariantClass %in% truncating

  # 4. Any variant reported as somatic at least 20 times in COSMIC
  MUTS$CosmicCount <- as.numeric(MUTS$CosmicCount)
  MUTS <- MUTS %>%
    mutate(
      ch_pd = ifelse(CosmicCount >= 20 & Gene %in% gene_list$Gene,
        1,
        ch_pd
      ),
      WHY_CH = ifelse(CosmicCount >= 20 & Gene %in% gene_list$Gene,
        paste0(WHY_CH, " Cosmic;"),
        WHY_CH
      )
    )

  MUTS$ch_pd <- tidyr::replace_na(MUTS$ch_pd, 0)
  annotate_PD.R <- MUTS %>% distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)

  if (all(annotate_PD.R[, c("CHROM", "POS", "REF", "ALT")] == df[, c("CHROM", "POS", "REF", "ALT")])) {
    message("good same dim as orig")
    topmed.mutation.2 <- supportData[["topmed.mutation.2"]]
    kelly.mutation.2 <- supportData[["kelly.mutation.2"]]
    ## both Bick/Bolton
    matches.2.c.p <- supportData[["matches.2.c.p"]]
    vars <- supportData[["vars"]]
    annotate_PD.R <- annotate_PD.R %>%
      dplyr::mutate(
        gene_loci_p = paste(Gene, paste0(aa_ref, aa_pos), sep = "_"),
        cDNAchange = gsub(".*:", "", HGVSc_VEP),
        gene_loci_c = paste(Gene, cDNAchange, sep = "_"),
        gene_aachange = paste(Gene, AAchange, sep = "_"),
        gene_loci_vep = ifelse(is.na(AAchange) | AAchange == "", gene_loci_c, gene_loci_p),
        ch_pd = case_when(
          ch_pd == 1 ~ 1,
          (gene_loci_vep %in% vars$gene_loci_vep |
            gene_aachange %in% matches.2.c.p$exact_match |
            gene_loci_c %in% matches.2.c.p$exact_match) ~ 1,
          TRUE ~ 0
        ),
        WHY_CH = ifelse(gene_loci_vep %in% vars$gene_loci_vep,
          paste0(WHY_CH, " Gene_loci>=5;"),
          WHY_CH
        ),
        WHY_CH = ifelse(gene_aachange %in% matches.2.c.p$exact_match,
          paste0(WHY_CH, " gene_aachange>=2 BB;"),
          WHY_CH
        ),
        WHY_CH = ifelse(gene_loci_c %in% matches.2.c.p$exact_match,
          paste0(WHY_CH, " gene_loci_c>=2 BB;"),
          WHY_CH
        ),
        WHY_CH = ifelse(gene_aachange %in% topmed.mutation.2$exact_match,
          paste0(WHY_CH, " gene_aachange>=2 Bick;"),
          WHY_CH
        ),
        WHY_CH = ifelse(gene_loci_c %in% topmed.mutation.2$exact_match,
          paste0(WHY_CH, " gene_loci_c>=2 Bick;"),
          WHY_CH
        ),
        WHY_CH = ifelse(gene_aachange %in% kelly.mutation.2$exact_match,
          paste0(WHY_CH, " gene_aachange>=2 Bolton;"),
          WHY_CH
        ),
        WHY_CH = ifelse(gene_loci_c %in% kelly.mutation.2$exact_match,
          paste0(WHY_CH, " gene_loci_c>=2 Bolton;"),
          WHY_CH
        )
      )

    return(
      annotate_PD.R %>% #
        dplyr::select(
          "CHROM", "POS", "REF", "ALT",
          "COSMIC_ID", "CosmicCount", "heme_cosmic_count", "myeloid_cosmic_count",
          "VariantClass", "Gene", "oncoKB", "oncoKB_reviewed",
          "isOncogenic", "isTSG", "isTruncatingHotSpot", "ch_pd", "WHY_CH"
        )
    )
  } else {
    message("Annotate PD failed: something wrong with dims")
  }
}
