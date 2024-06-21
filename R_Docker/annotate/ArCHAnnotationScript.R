#!/usr/local/bin/Rscript
options(gsubfn.engine = "R")

library(tidyselect)
library(tidyverse)
library(data.table)
library(vcfR)
library(httr)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(R453Plus1Toolbox)
library(sqldf)
library(jsonlite)
library(optparse)
source("/opt/bin/annotate/supportFunctions.R", local = TRUE)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input VCF File e.g. mutect.sample_name.final.annotated.vcf.gz", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="output", 
              help="Output Filename e.g. mutect.sample_name.final.annotated.Rscript", metavar="character"),
  make_option(c("-c", "--caller"), type="character", default=NULL,
              help="The VCF from which caller e.g. mutect | lofreq | vardict", metavar="character"),
  make_option("--bolton_bick_vars", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/bick.bolton.vars3.txt",
              help="The directory path where the 'bick.bolton.vars' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--mut2_bick", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/topmed.n2.mutation.c.p.txt",
              help="The directory path where the 'bick_topmed' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--mut2_kelly", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/kelly.n2.mutation.c.p.txt",
              help="The directory path where the 'kelly_impact' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--matches2", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/matches.2.c.p.txt",
              help="The directory path where the 'bick_kelly' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--truncating", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/BB.truncating.more.than.1.tsv",
              help="The directory path where the 'truncating' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--gene_list", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/oncoKB_CGC_pd_table_disparity_KB_BW.csv",
              help="The directory path where the 'putative driver annotations for specific genes' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--oncokb_genes", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/oncoKbCancerGeneList.tsv",
              help="The directory path where the 'oncoKB available genes' data for annotate_PD is stored. [default_path = %default]", metavar="character"),
  make_option("--cosmic_dir", type="character", default="/storage1/fs1/bolton/Active/Protected/Annotation_Files/cosmic/",
              help="The directory path where the 'cosmic' files are stored for cosmic function.", metavar="character"),
  make_option("--api_key", type="character", default="https://www.oncokb.org/account/settings",
              help="OncoKB API Key [default = %default]", metavar="double"),
  make_option("--p_value", type="double", default=2.114164905e-6,
              help="The Bonferroni Corrected P-Value [default = %default]", metavar="double"),
  make_option("--csq_string", type="character", default="Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomADe_AF|gnomADe_AFR_AF|gnomADe_AMR_AF|gnomADe_ASJ_AF|gnomADe_EAS_AF|gnomADe_FIN_AF|gnomADe_NFE_AF|gnomADe_OTH_AF|gnomADe_SAS_AF|gnomADg_AF|gnomADg_AFR_AF|gnomADg_AMI_AF|gnomADg_AMR_AF|gnomADg_ASJ_AF|gnomADg_EAS_AF|gnomADg_FIN_AF|gnomADg_MID_AF|gnomADg_NFE_AF|gnomADg_OTH_AF|gnomADg_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|CADD_PHRED|CADD_RAW|REVEL|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|pLI_gene_value|clinvar|clinvar_ID|clinvar_AF_ESP|clinvar_AF_EXAC|clinvar_AF_TGP|clinvar_CLNSIG|clinvar_CLNSIGCONF|clinvar_CLNDN|clinvar_ORIGIN",
              help="The CSQ string for the VEP annotation [default = %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt$caller <- tolower(opt$caller)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Must supply <input_file_name>", call.=FALSE)
} else if (is.null(opt$out)) {
  print_help(opt_parser)
  stop("Must supply <output_file_name>", call.=FALSE)
} else if (is.null(opt$caller)) {
  print_help(opt_parser)
  stop("Must supply <caller>", call.=FALSE)
}

startTime <- Sys.time()

# Loading in the VCF Files
message("Reading in the VCF file...")
if (opt$caller == "mutect") {
  vcf <- vcfR2tidy(read.vcfR(opt$input), single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat
} else if (opt$caller == "lofreq") {
  vcf <- vcfR2tidy(read.vcfR(opt$input), info_only = TRUE, info_types = TRUE, format_types = TRUE)$fix
} else if (opt$caller == "vardict") {
  vcf <- vcfR2tidy(read.vcfR(opt$input), single_frame = TRUE, info_types = TRUE, format_types = TRUE)$dat
  vcf <- vcf %>% dplyr::rename(QUAL = QUAL...6, ReadQual = QUAL...20)
}
vcf <- as.data.frame(vcf)
message("Finished reading in VCF file.")

message("Processing DF into desired format...")
# Remove Off-Target & Introns
if (nrow(vcf %>% dplyr::filter(CSQ != "")) != 0){
  vcf <- vcf %>% dplyr::filter(CSQ != "")
} else {
  vcf <- vcf %>% mutate(CSQ = "Off-Target | Intron")
}

vcf <- vcf %>%
  dplyr::filter(PON_FISHER <= opt$p_value)

# Manipulation of VCF Dataframe
vcf <- vcf %>% mutate(
  key = paste0(CHROM, ":", POS, ":", REF, ":", ALT),
  FP_filter = str_replace_all(FP_filter, ",", ";"),
  FILTER = ifelse(
    FP_filter == "PASS",
    FILTER,
    ifelse(
      FILTER == "PASS",
      FP_filter,
      paste0(FILTER, ";", FP_filter)
    )
  ),
  subject = str_extract(SAMPLE, "^[^_]+")
)
vcf <- vcf %>%
  rename_with(~paste0(., "_", opt$caller), 
              .cols = c(which(grepl("gt_", colnames(vcf))),
                        which(colnames(vcf) == "QUAL"),
                        which(colnames(vcf) == "FILTER")))

# Reformating of REF and ALT count columns
# Reformatting StrandBias: SB= RefFor, RefRev, AltFor, AltRev
if (opt$caller == "mutect") {
  vcf <- vcf %>%
    separate(gt_AD_mutect, c("gt_AD_ref_mutect", "gt_AD_alt_mutect"), sep = ",", extra = "merge", fill = "right") %>%
    separate(gt_SB_mutect, c("RDF_mutect", "RDR_mutect", "ADF_mutect", "ADR_mutect"), sep = ",", extra = "merge", fill = "right") %>%
    mutate(gt_AF_mutect = as.numeric(as.character(gt_AF_mutect)),
           gt_AD_ref_mutect = as.numeric(as.character(gt_AD_ref_mutect)),
           gt_AD_alt_mutect = as.numeric(as.character(gt_AD_alt_mutect)),
           RDF_mutect = as.numeric(as.character(RDF_mutect)),
           RDR_mutect = as.numeric(as.character(RDR_mutect)),
           ADF_mutect = as.numeric(as.character(ADF_mutect)),
           ADR_mutect = as.numeric(as.character(ADR_mutect)))
} else if (opt$caller == "lofreq") {
  vcf <- vcf %>%
    separate(DP4, c("RDF_lofreq", "RDR_lofreq", "ADF_lofreq", "ADR_lofreq"), sep = ",", extra = "merge", fill = "right") %>%
    mutate(gt_AD_alt_lofreq = as.numeric(as.character(ADF_lofreq)) + as.numeric(as.character(ADR_lofreq)),
           gt_AD_ref_lofreq = as.numeric(as.character(RDF_lofreq)) + as.numeric(as.character(RDR_lofreq)),
           gt_AF_lofreq = as.numeric(as.character(AF)))
} else if (opt$caller == "vardict") {
  vcf <- vcf %>%
    separate(gt_AD_vardict, c("gt_AD_ref_vardict", "gt_AD_alt_vardict"), sep = ",", extra = "merge", fill = "right") %>%
    separate(gt_RD_vardict, c("RDF_vardict", "RDR_vardict"), sep = ",", extra = "merge", fill = "right") %>%
    separate(gt_ALD_vardict, c("ADF_vardict", "ADR_vardict"), sep = ",", extra = "merge", fill = "right") %>%
    mutate(gt_AF_vardict = as.numeric(as.character(gt_AF_vardict)),
           gt_AD_ref_vardict = as.numeric(as.character(gt_AD_ref_vardict)),
           gt_AD_alt_vardict = as.numeric(as.character(gt_AD_alt_vardict)),
           RDF_vardict = as.numeric(as.character(RDF_vardict)),
           RDR_vardict = as.numeric(as.character(RDR_vardict)),
           ADF_vardict = as.numeric(as.character(ADF_vardict)),
           ADR_vardict = as.numeric(as.character(ADR_vardict)))
}

# Split VEP String into individual colums
csq_string <- str_split(opt$csq_string, "\\|")[[1]]
vcf <- vcf %>% tidyr::separate(CSQ, paste0(csq_string, "_VEP"), sep="\\|", extra = "merge", fill = "right")
message("Finished processing DF into desired format.")

message("Processing gnomAD...")
# Annotate gnomAD
vcf <- annotateGnomad(vcf) %>%
  dplyr::rename(
    max_gnomADe_AF_VEP = MAX_gnomADe_AF_VEP,
    max_gnomADg_AF_VEP = MAX_gnomADg_AF_VEP
  ) %>%
  mutate(
    pass_max_sub_gnomAD_AF = (max_gnomADe_AF_VEP < 0.005 & max_gnomADg_AF_VEP < 0.005),
    max_sub_gnomAD_AF = pmax(max_gnomADe_AF_VEP, max_gnomADg_AF_VEP, na.rm = TRUE),
    max_pop_gnomAD_AF = pmax(gnomADe_AF_VEP, gnomADg_AF_VEP, na.rm = TRUE)
  )
message("Finished processing gnomAD.")

message("Calculating Complexity...")
vcf <- annotateComplexity(vcf)
message("Done Complexity")

# Load all support data
message("Loading support files...")
{
  data.table::setDTthreads(2)
  supportData <- list()
  
  # Cosmic
  chrList <- paste0("chr", c(1:22, "X", "Y"))
  supportData[["cosmic"]] <- lapply(chrList, function(x) {
    cosmic <- fread(paste0(opt$cosmic_dir, "CosmicMutantExport.final.minimal.", x, ".tsv"), header = T, sep = "\t", quote = "", drop = c(2, 4))
    colnames(cosmic) <- c("COSMIC_ID", "var_key", "CosmicCount", "heme_cosmic_count", "myeloid_cosmic_count", "Gene_HGVSp_VEP")
    cosmic
  })
  names(supportData[["cosmic"]]) <- chrList
  supportData[["fullCosmic"]] <- supportData[["cosmic"]] %>% do.call(what = rbind)
  supportData[["cosmic"]] <- NULL
  
  supportData[["topmed.mutation.2"]] <- fread(opt$mut2_bick, sep = "\t", header = T, data.table = FALSE)
  supportData[["kelly.mutation.2"]] <- fread(opt$mut2_kelly, sep = "\t", header = T, data.table = FALSE)
  supportData[["matches.2.c.p"]] <- fread(opt$matches2, sep = "\t", header = T, data.table = FALSE)
  supportData[["vars"]] <- fread(opt$bolton_bick_vars, sep = "\t", header = T, data.table = FALSE)
  
  supportData[["vars.truncating"]] <- fread(opt$truncating, sep = "\t", header = T, data.table = FALSE)
  supportData[["gene_list"]] <- fread(opt$gene_list, sep = ",", header = T, data.table = FALSE)
  supportData[["gene_list"]] <- supportData[["gene_list"]] %>% dplyr::filter(remove != 1)
  
  # Read updated list of genes available in oncoKB, including aliases
  # The file can be directly read from the API link but download once at the beginning for consistency
  tmp <- fread(opt$oncokb_genes, data.table = FALSE) %>%
    dplyr::filter(`OncoKB Annotated` == "Yes")
  tmp <- c(tmp$`Hugo Symbol`, str_squish(unlist(str_split(tmp$`Gene Aliases`, ",")))) %>%
    unique() %>%
    .[. != ""]
  supportData[["oncoKbAvailGene"]] <- tmp
  supportData[["oncoKbApiKey"]] <- opt$api_key
}
message("Done loading support files.")

message("Preparing AnnotatePD Data...")
vcf <- prepareAnnotatePdData(vcf, supportData)
message("Finished preparing AnnotatePD Data.")

message("Run AnnotatePD")
annotation <- annotatePD(vcf, supportData)
vcf <- left_join(vcf,
                    annotation,
                    by = c("CHROM", "POS", "REF", "ALT")
)
message("Finished AnnotatePD")

print(Sys.time() - startTime)

write.table(vcf, paste0(opt$out, ".tsv"), row.names = FALSE, sep="\t")
message(paste0("Finished: Writing ", opt$caller, " TSV file: ", opt$out, ".tsv"))
