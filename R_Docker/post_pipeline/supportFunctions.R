prepare_vars_file <- function(vars) {
  vars$aa.pos <- as.numeric(str_extract(vars$loci.vep, "\\d+"))  # e.g. DNMT3A_R882 --> 882
  vars$CHROM.pos <- paste0(vars$CHROM, "_", vars$POS)  # e.g. chr2_25234373
  vars$GENE.AA.POS <- ifelse(!is.na(vars$aa.pos), paste0(vars$SYMBOL_VEP, "_", as.integer(vars$aa.pos)), NA)  # e.g. DNMT3A_882
  vars$GENE.AA.POS[vars$GENE.AA.POS == "NA_NA"] <- NA     # e.g. When VEP annotates like DNMT3A_NA or ASXL1_NA
  vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
  vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:", "", vars$HGVSc_VEP), sep = "_")
  vars
}

prepare_cosmic_file <- function(ct) {
  ct$CHROM.POS <- paste0(ct$CHROM, "_", ct$POS)
  ct$GENE.AA.POS <- paste0(ct$gene, "_", ct$aa.pos)
  ct$gene_cDNAchange <- paste0(ct$gene, "_", ct$cDNAchange)
  ct$gene_aachange <- paste0(ct$gene, "_", ct$AAchange)
  ct
}

near_BB_loci_HS <- function(df_row, vars) {
  p = c(-3:0,1:3)       # Establishes 3 AA downstream and upstream
  n = c(-9:0,1:9)       # Establishes the 9 nucleotides downstream and upstream
  
  prot = p + as.integer(df_row[["aa.pos"]])                   # e.g. DNMT3A_R882 would result in [879, 880, 881, 882, 883, 884, 885]
  vector.p = paste(df_row[["SYMBOL_VEP"]], prot ,sep = "_")   # which would result in [DNMT3A_879, DNMT3A_880, ...]
  
  nuc = n + as.integer(df_row[["POS"]])                       # e.g. DNMT3A_R882 would result in [25234364, 25234365, ...]
  vector.n = paste(df_row[["CHROM"]], nuc ,sep = "_")         # which would result in [chr2_25234364, chr2_25234365, ...]
  
  # We can store these and check if any of these AA or nucleotide positions are inside vars, it means this mutation is NEAR a BB Hotspot
  if (any(vector.p %in% vars$GENE.AA.POS)) {
    res <- unique(vars[vars$GENE.AA.POS %in% vector.p, c("gene_loci_vep", "truncating", "n.loci.truncating.vep")])
    ter_present <- grepl("Ter", df_row[["gene_aachange"]])
    res <- res %>% filter(truncating == ifelse(ter_present, "truncating", "not"))     # Only compare truncating variants to truncating hotspots
    if (nrow(res) != 0){
      # We only need to return truncating because in bick.bolton.vars3.txt the truncating is grouped by "truncating" so this column has the correct
      # value for the loci count depending on whether or not it is a truncating mutation or not. n.loci.vep will have the TOTAL number
      # e.g. If truncating = 9 and missense = 3, n.loci.vep = 12, and n.loci.truncating.vep = 9 or = 3 depending on which one it is
      return_string <- paste(res$gene_loci_vep, ":", res$n.loci.truncating.vep)
      return(paste(return_string, collapse = " | "))
    } else {
      return("")
    }
  } else if (any(vector.n %in% vars$CHROM.POS)) {
    res <- unique(vars[vars$CHROM.POS %in% vector.n, c("CHROM.POS", "n.HGVSc", "gene_cDNAchange", "gene_loci_vep")])
    del_ins_dup_present <- grepl("del|ins|dup", df_row[["gene_cDNAchange"]])
    res <- res %>% filter(grepl("del|ins|dup",gene_cDNAchange) == del_ins_dup_present)  # Like with the truncating, only compare del, ins, and dup to itself
    res <- res %>% group_by(gene_loci_vep) %>% summarise(n = sum(n.HGVSc))
    if (nrow(res) != 0){
      return_string <- paste(res$gene_loci_vep, ":", res$n)
      return(paste(return_string, collapse = " | "))
    } else {
      return("")
    }
  } else {
    return("")
  }
}

near_COSMIC_loci_HS <- function(df_row, ct) {
  p = c(-3:0,1:3)       # Establishes 3 AA downstream and upstream
  n = c(-9:0,1:9)       # Establishes the 9 nucleotides downstream and upstream
  
  prot = p + as.integer(df_row[["aa.pos"]])                   # e.g. DNMT3A_R882 would result in [879, 880, 881, 882, 883, 884, 885]
  vector.p = paste(df_row[["SYMBOL_VEP"]], prot ,sep = "_")   # which would result in [DNMT3A_879, DNMT3A_880, ...]
  
  nuc = n + as.integer(df_row[["POS"]])                       # e.g. DNMT3A_R882 would result in [25234364, 25234365, ...]
  vector.n = paste(df_row[["CHROM"]], nuc ,sep = "_")         # which would result in [chr2_25234364, chr2_25234365, ...]
  
  # We can store these and check if any of these AA or nucleotide positions are inside vars, it means this mutation is NEAR a BB Hotspot
  if (any(vector.p %in% ct$GENE.AA.POS)) {
    res <- unique(ct[ct$GENE.AA.POS %in% vector.p & (ct$cosmic_count.loci.truncating >= 25 | ct$heme_count.loci.truncating >= 10 | ct$myeloid_count.loci.truncating >= 5), c('gene_loci_vep', 'truncating', 'cosmic_count.loci.truncating', 'heme_count.loci.truncating', 'myeloid_count.loci.truncating')])
    ter_present <- grepl("Ter", df_row[["gene_aachange"]])
    res <- res %>% filter(truncating == ter_present)     # Only compare truncating variants to truncating hotspots
    if (nrow(res) != 0){
      return_string <- paste(res$gene_loci_vep, ":", res$cosmic_count.loci.truncating, ",", res$heme_count.loci.truncating, ",", res$myeloid_count.loci.truncating)
      return(paste(return_string, collapse = " | "))
    } else {
      return("")
    }
  } else if (any(vector.n %in% ct$CHROM.POS)) {
    res <- unique(ct[ct$CHROM.POS %in% vector.n, c('gene_loci_vep', 'cosmic_count.totals.c', 'heme_count.totals.c', 'myeloid_count.totals.c', 'gene_cDNAchange')])
    del_ins_dup_present <- grepl("del|ins|dup", df_row[["gene_cDNAchange"]])
    res <- res %>% filter(grepl("del|ins|dup",gene_cDNAchange) == del_ins_dup_present)  # Like with the truncating, only compare del, ins, and dup to itself
    res <- res %>% group_by(gene_loci_vep) %>% summarise(n_cosmic = sum(cosmic_count.totals.c), n_heme = sum(heme_count.totals.c), n_myeloid = sum(myeloid_count.totals.c))
    res <- res[(res$n_cosmic >= 25) | (res$n_heme >= 10) | (res$n_myeloid >= 5), ]
    if (nrow(res) != 0){
      return_string <- paste(res$gene_loci_vep, ":", res$n_cosmic, ",", res$n_heme, ",", res$n_myeloid)
      return(paste(return_string, collapse = " | "))
    } else {
      return("")
    }
  } else {
    return("")
  }
}
