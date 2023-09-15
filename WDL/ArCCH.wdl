version 1.0

# ArCH WDL Pipeline for TERRAbio
# -------------------------------
# This pipeline is designed to process mutant/wildtype H.sapiens sequencing data from ArcherDX for low VAF variants.
# It features four variant callers (Mutect, Vardict, Lofreq, Pindel) for variant detection and
# performs various false positive filters and detection methods (fp_filter, PoN FishersTest, XGB Model).
# This pipeline also generates VEP style annotations for all called variants as well as additional putative driver
# annotations generated from various database sources (TOPMed, MSK-IMPACT, COSMIC, OncoKB, etc.)

# Created by: Irenaeus Chan
# Contact: chani@wustl.edu
# Date: 03/17/2022

# A file with a label. E.g. A bed file that requires a specific label to identify it.
struct LabelledFile {
    File file
    String label
}

# VEP can utilize other files for custom annotations outside of the normal available plugins
# This structure format is ported over from MGI's CWL Pipelines so it could potentially be better integrated
struct VepCustomAnnotation {
    Boolean check_existing
    File custom_file
    String name
    String data_format  # enum, ['bed', 'gff', 'gtf', 'vcf', 'bigwig']
    String method  # enum, ['exact', 'overlap']
    Boolean force_report_coordinates
    Array[String]? vcf_fields
    Array[File]? secondary_files
}

# Main Workflow
workflow ArCH {
    input {
        # Input Files
        File input_file                     # This can be a FASTQ, BAM, or CRAM
        File input_file_two                 # This will be R2, BAI, or CRAI
        String tumor_sample_name

        # If Normal Samples Exists
        Boolean tumor_only = false
        File? normal_bam
        File? normal_bai
        String? normal_sample_name

        # Input File Format
        String input_type = "FASTQ"         # Options: "FASTQ", "BAM", "CRAM" (Default: "FASTQ")
        Boolean aligned = false             # If input is an already aligned BAM file then set this flag
    
        # Sequence Information
        String platform = "Illumina"
        String platform_unit = "ArcherDX"
        String library = "LIBRARY"

        File target_intervals               # Interval List

        # Reference
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa

        # UMI Processing
        Boolean has_umi = true
        Boolean? umi_paired = false         # If the UMI is paired (R1 and R2) then set this flag
        String where_is_umi = "T"           # Three options "N = Name", "R = Read", or "T = Tag"
        Array[String] read_structure        # Used for the UMI processing see: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures

        # Consensus Building Parameters
        Array[Int] min_reads = [1]          # The minimum number of reads that constitutes a "read family"
        Float? max_read_error_rate = 0.05   # If 5% of the reads within a "read family" does not match, removes the family
        Float? max_base_error_rate = 0.1    # If 10% of the bases are different than the majority, masks the base with N
        Int min_base_quality = 1            # Any base less than QUAL 1 will be masked with N
        Float max_no_call_fraction = 0.5    # Maximum fraction of no-calls (N) in the read after filtering

        # BQSR
        # https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf
        Boolean apply_bqsr = false
        Array[File] bqsr_known_sites
        Array[File] bqsr_known_sites_tbi

        # Variants within gnomAD that have a VAF higher than 0.5%
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi

        # Somalier
        File chrom_sizes
        File af_only_snp_only_vcf

        # Variant Calling
        Array[Pair[File, File]] pon_bams            # List of BAMs within the Panel of Normals (PoN)
        Float? af_threshold = 0.0001                # Minimum VAF Cut-Off
        String? pon_pvalue = "2.114164905e-6"       # Bonferroni Corrected P-Value for Significance - (Default: Archer Panel)

        # See: https://github.com/bcbio/bcbio_validations/blob/master/somatic-lowfreq/README.md
        # Parameters MQ, NM, DP, and QUAL are calculated using a small subset then identifying the cut-off for 2% of the left side samples
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))"

        # PoN2
        # If a variant exists inside two or more of our Panel of Normal Samples (PoN) at 2% VAF or greater, then this variant is assumed to be
        # sequencing noise because our PoN should not have any variants (that are not germline).
        # Instructions:
        # PoN Samples are run through Tumor Only mode and Filtered for each caller.
        # Variants below 2% VAF are removed from PoN Samples
        # Variants that appear within the PoN Samples in two or more samples are kept
        File mutect_pon2_file
        File mutect_pon2_file_tbi
        File lofreq_pon2_file
        File lofreq_pon2_file_tbi
        File vardict_pon2_file
        File vardict_pon2_file_tbi

        # R Files for Putative Driver Annotations
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File truncating
        File gene_list
        File oncokb_genes
        File cosmic_dir_zip

        # VEP Parameters
        File vep_cache_dir_zip                      # WDL does not have a Directory Variable, so the entire cache needs to be ZIP
        Array[String] vep_plugins = ["Frameshift", "Wildtype"]
        File clinvar_vcf
        File clinvar_vcf_tbi
        File? synonyms_file
        Boolean? annotate_coding_only = true
        #Array[VepCustomAnnotation] vep_custom_annotations
    }

    # If the BAM file is already aligned and consensus sequencing was done, then alignment can be skipped
    if (!aligned) {
        # Archer UMIs are not typical, they have a 13 bp Adapter after their 8 bp UMI
        # There needs to be some pruning involved before we can extract the UMIs
        call filterArcherUMILengthAndbbmapRepair as repair {
            input:
            input_file = input_file,
            input_file_two = input_file_two,
            input_type = input_type,
            sample_name = tumor_sample_name,
            library_name = library,
            platform_unit = platform_unit,
            platform = platform
        }

        call processUMIs {
            input:
            bam = repair.bam,
            umi_paired = umi_paired,
            where_is_umi = where_is_umi,
            read_structure = read_structure
        }

        call sortAndMarkIlluminaAdapters {
            input:
            bam = processUMIs.umi_extracted_bam
        }

        if (has_umi) {
            # First Alignment
            call umiAlign as align {
                input:
                bam = sortAndMarkIlluminaAdapters.marked_bam,
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                reference_amb = reference_amb,
                reference_ann = reference_ann,
                reference_bwt = reference_bwt,
                reference_pac = reference_pac,
                reference_sa = reference_sa
            }

            # Create Read Families and Perform Consensus Calling
            call groupReadsAndConsensus {
                input:
                bam = align.aligned_bam,
                umi_paired = umi_paired,
            }
        }

        # Realign the Consensus Called Reads
        call umiAlign as realign {
            input:
            bam = select_first([groupReadsAndConsensus.consensus_bam, sortAndMarkIlluminaAdapters.marked_bam]),
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            reference_amb = reference_amb,
            reference_ann = reference_ann,
            reference_bwt = reference_bwt,
            reference_pac = reference_pac,
            reference_sa = reference_sa,
            realign = true
        }

        # Filter and Clip
        call filterAndClip {
            input:
            bam = realign.aligned_bam,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            sample_name = tumor_sample_name,
            min_reads = min_reads,
            max_read_error_rate = max_read_error_rate,
            max_base_error_rate = max_base_error_rate,
            min_base_quality = min_base_quality,
            max_no_call_fraction = max_no_call_fraction,
            has_umi = has_umi
        }
    } # END OF ALIGNMENT STEP

    # This will handle BQSR, CRAM to BAM, and whatever processing is necessary to produce the BAM that will be used for variant calling
    call processAlignedBAM as initialBAM {
        input:
        bam = select_first([filterAndClip.clipped_bam, input_file]),
        bai = select_first([filterAndClip.clipped_bai, input_file_two]),
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        interval_list = target_intervals,
        known_sites = bqsr_known_sites,
        known_sites_tbi = bqsr_known_sites_tbi,
        sample_name = tumor_sample_name,
        apply_bqsr = apply_bqsr,
        input_type = input_type
    }

    call fastQC {
        input:
        bam = initialBAM.initial_bam,
        bai = initialBAM.initial_bai
    }

    # Some of our callers use BED file instead of interval list
    call intervalsToBed as interval_to_bed {
        input: interval_list = target_intervals
    }

    # Perform Somalier
    call createSomalierVcf {
        input:
        interval_bed = interval_to_bed.interval_bed,
        chrom_sizes = chrom_sizes,
        af_only_snp_only_vcf = af_only_snp_only_vcf,
        reference = reference
    }

    call somalier {
        input:
        somalier_vcf = createSomalierVcf.somalier_vcf,
        reference = reference,
        bam = initialBAM.initial_bam,
        bai = initialBAM.initial_bai,
        sample_name = tumor_sample_name
    }

    # In order to parallelize as much as the workflow as possible, we analyze by chromosome
    call splitBedToChr {
        input:
            interval_bed = interval_to_bed.interval_bed
    }

    scatter (bed_chr in splitBedToChr.split_chr) {
        # Mutect
        call mutect {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            gnomad = normalized_gnomad_exclude,
            gnomad_tbi = normalized_gnomad_exclude_tbi,
            tumor_bam = initialBAM.initial_bam,
            tumor_bai = initialBAM.initial_bai,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            interval_list = bed_chr,
            tumor_only = tumor_only
        }

        # Vardict
        call vardict {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = initialBAM.initial_bam,
            tumor_bai = initialBAM.initial_bai,
            tumor_sample_name = tumor_sample_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            normal_sample_name = normal_sample_name,
            interval_bed = bed_chr,
            min_var_freq = af_threshold,
            tumor_only = tumor_only
        }

        call lofreq {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = initialBAM.initial_bam,
            tumor_bai = initialBAM.initial_bai,
            tumor_sample_name = tumor_sample_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            interval_bed = bed_chr,
            tumor_only = tumor_only
        }

        call pindel {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = initialBAM.initial_bam,
            tumor_bai = initialBAM.initial_bai,
            tumor_sample_name = tumor_sample_name,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            normal_sample_name = normal_sample_name,
            region_file = bed_chr,
            tumor_only = tumor_only
        }

        call removeEndTags {
            input: pindel_vcf = pindel.vcf
        }
    } # END OF VARIANT CALLING

    call mergeVcf as merge_mutect {
        input:
            vcfs = mutect.vcf,
            vcf_tbis = mutect.vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name
    }
    call mergeVcf as merge_vardict {
        input:
            vcfs = vardict.vcf,
            vcf_tbis = vardict.vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name
    }
    call mergeVcf as merge_lofreq {
        input:
            vcfs = lofreq.vcf,
            vcf_tbis = lofreq.vcf_tbi,
            merged_vcf_basename = "lofreq." + tumor_sample_name
    }
    call mergeVcf as merge_pindel {
        input:
            vcfs = removeEndTags.vcf,
            vcf_tbis = removeEndTags.vcf_tbi,
            merged_vcf_basename = "pindel." + tumor_sample_name
    }

    # Cleans the VCF output that don't match the expected VCF Format
    # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
    # Removes any germine variant that is reported in gnomAD
    # Removes any variants that appeared in two of our PoN Samples at 2% VAF or higher
    call sanitizeNormalizeFilter as mutect_filter {
        input:
            vcf = merge_mutect.merged_vcf,
            vcf_tbi = merge_mutect.merged_vcf_tbi,
            reference = reference,
            reference_fai = reference_fai,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            vcf2PON = mutect_pon2_file,
            vcf2PON_tbi = mutect_pon2_file_tbi,
            caller = "mutect",
            sample_name = tumor_sample_name,
            filter_string = bcbio_filter_string
    }

    call sanitizeNormalizeFilter as vardict_filter {
        input:
            vcf = merge_vardict.merged_vcf,
            vcf_tbi = merge_vardict.merged_vcf_tbi,
            reference = reference,
            reference_fai = reference_fai,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            vcf2PON = vardict_pon2_file,
            vcf2PON_tbi = vardict_pon2_file_tbi,
            caller = "vardict",
            sample_name = tumor_sample_name,
            filter_string = bcbio_filter_string
    }

    call sanitizeNormalizeFilter as lofreq_filter {
        input:
            vcf = merge_lofreq.merged_vcf,
            vcf_tbi = merge_lofreq.merged_vcf_tbi,
            reference = reference,
            reference_fai = reference_fai,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            vcf2PON = lofreq_pon2_file,
            vcf2PON_tbi = lofreq_pon2_file_tbi,
            caller = "lofreq",
            sample_name = tumor_sample_name,
            filter_string = bcbio_filter_string
    }

    call sanitizeNormalizeFilter as pindel_filter {
        input:
            vcf = merge_pindel.merged_vcf,
            vcf_tbi = merge_pindel.merged_vcf_tbi,
            reference = reference,
            reference_fai = reference_fai,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            vcf2PON = merge_pindel.merged_vcf,
            vcf2PON_tbi = merge_pindel.merged_vcf_tbi,
            caller = "pindel",
            sample_name = tumor_sample_name,
            filter_string = bcbio_filter_string
    }

    # In order to be efficient, we run all of the annotation and filtering ONCE. In order to do this, we need to merge
    # all of the callers together into one giant VCF that will be used as the main VCF for VEP, PoN, etc...
    scatter (caller_vcf in [mutect_filter.annotated_vcf, vardict_filter.annotated_vcf, lofreq_filter.annotated_vcf]){
        call createFakeVcf as fake_vcf {
            input:
            vcf = caller_vcf,
            tumor_sample_name = tumor_sample_name
        }
    }

    call mergeVcf as mergeCallers {
        input:
        vcfs = fake_vcf.fake_vcf,
        vcf_tbis = fake_vcf.fake_vcf_tbi,
        merged_vcf_basename = "all_callers." + tumor_sample_name
    }

    # This is VarScan2's FP Filter. https://github.com/ucscCancer/fpfilter-tool/blob/master/fpfilter.pl
    call fpFilter {
        input:
        reference=reference,
        reference_fai=reference_fai,
        bam = initialBAM.initial_bam,
        vcf=mergeCallers.merged_vcf,
        sample_name=tumor_sample_name,
        min_var_freq=af_threshold,
        output_vcf_basename = "all_callers." + tumor_sample_name + ".fpfilter"
    }

    # Using the fake VCF of all the variants calls, perform a "pileup" on the PoN BAMs`
    scatter (pon_bam in pon_bams) {
        call mskGetBaseCounts {
            input:
            reference = reference,
            reference_fai = reference_fai,
            normal_bam = pon_bam,
            pon_final_name = "all_callers." + tumor_sample_name + ".pon.pileup",
            vcf = mergeCallers.merged_vcf
        }
    }
    call bcftoolsMerge as pileup_merge {
        input:
            vcfs = mskGetBaseCounts.pileup,
            vcf_tbis = mskGetBaseCounts.pileup_tbi,
            merged_vcf_basename = tumor_sample_name + ".pon.total.counts",
            RD_AD = true
    }

    # Using the fake VCF of all the variants calls, perform VEP on all the variants
    call vep {
        input:
            vcf = mergeCallers.merged_vcf,
            cache_dir_zip = vep_cache_dir_zip,
            reference = reference,
            reference_fai = reference_fai,
            plugins = vep_plugins,
            synonyms_file = synonyms_file,
            coding_only = annotate_coding_only,
            clinvar = clinvar_vcf,
            clinvar_tbi = clinvar_vcf_tbi
    }

    # Using the "pileup", perform a Fisher's Exact Test with the Variants in each Caller
    call normalFisher as mutect_call_R_fisher {
        input:
        vcf = mutect_filter.annotated_vcf,
        pon = select_first([pileup_merge.pileup_vcf,pileup_merge.merged_vcf]),
        pon_tbi = select_first([pileup_merge.pileup_vcf_tbi, pileup_merge.merged_vcf_tbi]),
        p_value = pon_pvalue,
        caller = "mutect"
    }
    call annotateVcf as mutect_annotate_vcf {
        input:
        vcf = mutect_call_R_fisher.pon_filtered_vcf,
        vcf_tbi = mutect_call_R_fisher.pon_filtered_vcf_tbi,
        fp_filter = fpFilter.filtered_vcf,
        fp_filter_tbi = fpFilter.filtered_vcf_tbi,
        vep = vep.annotated_vcf,
        vep_tbi = vep.annotated_vcf_tbi,
        caller_prefix = "mutect",
        sample_name = tumor_sample_name
    }

    call normalFisher as vardict_call_R_fisher {
        input:
        vcf = vardict_filter.annotated_vcf,
        pon = select_first([pileup_merge.pileup_vcf,pileup_merge.merged_vcf]),
        pon_tbi = select_first([pileup_merge.pileup_vcf_tbi, pileup_merge.merged_vcf_tbi]),
        p_value = pon_pvalue,
        caller = "vardict"
    }
    call annotateVcf as vardict_annotate_vcf {
        input:
        vcf = vardict_call_R_fisher.pon_filtered_vcf,
        vcf_tbi = vardict_call_R_fisher.pon_filtered_vcf_tbi,
        fp_filter = fpFilter.filtered_vcf,
        fp_filter_tbi = fpFilter.filtered_vcf_tbi,
        vep = vep.annotated_vcf,
        vep_tbi = vep.annotated_vcf_tbi,
        caller_prefix = "vardict",
        sample_name = tumor_sample_name
    }

    call normalFisher as lofreq_call_R_fisher {
        input:
        vcf = lofreq_filter.annotated_vcf,
        pon = select_first([pileup_merge.pileup_vcf,pileup_merge.merged_vcf]),
        pon_tbi = select_first([pileup_merge.pileup_vcf_tbi, pileup_merge.merged_vcf_tbi]),
        p_value = pon_pvalue,
        caller = "lofreq"
    }
    call annotateVcf as lofreq_annotate_vcf {
        input:
        vcf = lofreq_call_R_fisher.pon_filtered_vcf,
        vcf_tbi = lofreq_call_R_fisher.pon_filtered_vcf_tbi,
        fp_filter = fpFilter.filtered_vcf,
        fp_filter_tbi = fpFilter.filtered_vcf_tbi,
        vep = vep.annotated_vcf,
        vep_tbi = vep.annotated_vcf_tbi,
        caller_prefix = "lofreq",
        sample_name = tumor_sample_name
    }

    call annotatePD {
        input:
            mutect_vcf = mutect_annotate_vcf.final_annotated_vcf,
            lofreq_vcf = lofreq_annotate_vcf.final_annotated_vcf,
            vardict_vcf = vardict_annotate_vcf.final_annotated_vcf,
            bolton_bick_vars = bolton_bick_vars,
            mut2_bick = mut2_bick,
            mut2_kelly = mut2_kelly,
            matches2 = matches2,
            truncating = truncating,
            gene_list = gene_list,
            oncokb_genes = oncokb_genes,
            cosmic_dir_zip = cosmic_dir_zip,
            pon_pvalue = pon_pvalue
    }

    call combine_all {
        input:
            mutect_tsv = annotatePD.mutect_vcf_annotate_pd,
            lofreq_tsv = annotatePD.lofreq_vcf_annotate_pd,
            vardict_tsv = annotatePD.vardict_vcf_annotate_pd,
            pindel_vcf = pindel_filter.annotated_vcf,
            pon = pileup_merge.merged_vcf,
            pon_pvalue = pon_pvalue,
            model = false,
            tumor_sample_name = tumor_sample_name
    }

    output {
        # Alignments
        File aligned_bam = initialBAM.initial_bam

        # Population Level Information
        File fastqc_html = fastQC.fastqc_html
        File fastqc = fastQC.fastqc
        File somalier_out = somalier.somalier_out

        # Mutect
        File mutect_vcf = mutect_filter.annotated_vcf                                   # Raw Mutect Ouput (gnomAD Filtered + PoN2 Annotated)
        File mutect_pon_annotated_vcf = mutect_call_R_fisher.pon_vcf                    # PoN Pileup Annotated
        File mutect_pon_filtered_vcf = mutect_call_R_fisher.pon_filtered_vcf            # PoN Pileup Filtered
        File mutect_fpfilter_vcf = mutect_annotate_vcf.fpfilter_annotated_vcf           # PoN Filtered + FPFilter
        File mutect_vep_annotated_vcf = mutect_annotate_vcf.final_annotated_vcf         # PoN Filtered + FPFilter + VEP Annotation

        # Lofreq
        File lofreq_vcf = lofreq_filter.annotated_vcf                                   # Raw Mutect Ouput (gnomAD Filtered + PoN2 Annotated)
        File lofreq_pon_annotated_vcf = lofreq_call_R_fisher.pon_vcf                    # PoN Pileup Annotated
        File lofreq_pon_filtered_vcf = lofreq_call_R_fisher.pon_filtered_vcf            # PoN Pileup Filtered
        File lofreq_fpfilter_vcf = lofreq_annotate_vcf.fpfilter_annotated_vcf           # PoN Filtered + FPFilter
        File lofreq_vep_annotated_vcf = lofreq_annotate_vcf.final_annotated_vcf         # PoN Filtered + FPFilter + VEP Annotation

        # Vardict
        File vardict_vcf = vardict_filter.annotated_vcf                                 # Raw Mutect Ouput (gnomAD Filtered + PoN2 Annotated)
        File vardict_pon_annotated_vcf = vardict_call_R_fisher.pon_vcf                  # PoN Pileup Annotated
        File vardict_pon_filtered_vcf = vardict_call_R_fisher.pon_filtered_vcf          # PoN Pileup Filtered
        File vardict_fpfilter_vcf = vardict_annotate_vcf.fpfilter_annotated_vcf         # PoN Filtered + FPFilter
        File vardict_vep_annotated_vcf = vardict_annotate_vcf.final_annotated_vcf       # PoN Filtered + FPFilter + VEP Annotation

        # Pindel
        File pindel_vcf = pindel_filter.annotated_vcf                                   # Raw Pindel Ouput

        File? pon_total_counts = pileup_merge.pileup_vcf                                 # PoN Pileup Results
        File fpfilter_results = fpFilter.filtered_vcf                                   # FPFilters Results
        File vep_results = vep.annotated_vcf                                            # VEP Results

        # R Things
        File mutect_annotate_pd = annotatePD.mutect_vcf_annotate_pd
        File lofreq_annotate_pd = annotatePD.lofreq_vcf_annotate_pd
        File vardict_annotate_pd = annotatePD.vardict_vcf_annotate_pd

        # Final
        File final = combine_all.final_output
        File? model_output = combine_all.model_output
    }
}

task filterArcherUMILengthAndbbmapRepair {
    input {
        File input_file                     # This can be an unaligned BAM or FASTQ
        File? input_file_two                # If unaligned BAM... we don't need this. If FASTQ, this is R2
        Int? umi_length = 8
        String input_type = "FASTQ"         # FASTQ, BAM, or CRAM
        String sample_name
        String library_name
        String platform_unit
        String platform
        Float? mem_limit_override
    }

    Float data_size = size(input_file, "GB")
    Int space_needed_gb = ceil(6 * data_size)                   # We need at least 2*3*Data because if it's FASTQ, R1 and R2, then 3 for two steps
    Float memory = select_first([mem_limit_override, 8])
    Int cores = 1
    Int java_mem = floor(memory)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/bbmap:39.01--h92535d8_1"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        # If it's ArcherDX Sequencing, we can either have a BAM input or FASTQ input. Either way, we need to filter out the 13 bp Adapter Issue
        if [[ "~{platform}" == "ArcherDX" ]]; then
            if [ "~{input_type}" == "BAM" ]; then
                samtools view --no-header ~{input_file} | awk -v regex="AACCGCCAGGAGT" -v umi_length="~{umi_length}" -F'\t' '{if ($2 == "77") { split($10,a,regex); if(length(a[1]) == umi_length){print$0}} else {print $0}}' > filtered.sam
                repair.sh -Xmx~{java_mem}g \
                repair=t \
                overwrite=true \
                interleaved=false \
                in=stdin \
                out1=R1.fixed.fastq.gz \
                out2=R2.fixed.fastq.gz
            elif [ "~{input_type}" == "FASTQ" ]; then
                zcat ~{input_file} | awk -v regex="AACCGCCAGGAGT" -v umi_length="~{umi_length}" 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; split(seq,a,regex); if (length(a[1]) == umi_length) {print header, seq, qheader, qseq}}' | \
                repair.sh -Xmx~{java_mem}g \
                repair=t \
                overwrite=true \
                interleaved=false \
                in1=stdin \
                in2=~{input_file_two} \
                out1=R1.fixed.fastq.gz \
                out2=R2.fixed.fastq.gz
            fi
            samtools import -@ ~{cores} \
            -1 R1.fixed.fastq.gz \
            -2 R2.fixed.fastq.gz \
            -O BAM \
            -o unaligned.bam \
            -r ID:A -r SM:~{sample_name} -r LB:~{library_name} -r PL:~{platform} -r PU:~{platform_unit}
        else
            # For any other type of sequencing, we either leave it as a BAM input or we convert FASTQ to BAM
            if [ "~{input_type}" == "BAM" ]; then
                cp ~{input_file} unaligned.bam
            else
                samtools import -@ ~{cores} \
                -1 ~{input_file} \
                -2 ~{input_file_two} \
                -O BAM \
                -o unaligned.bam \
                -r ID:A -r SM:~{sample_name} -r LB:~{library_name} -r PL:~{platform} -r PU:~{platform_unit}
            fi
        fi
        
    >>>

    output {
        File bam = "unaligned.bam"
    }
}

task processUMIs {
    input {
        File bam
        Boolean? umi_paired = true
        String where_is_umi = "T"       # T = Tag, R = Read Name, N = Name
        Array[String] read_structure
    }

    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(3 * data_size)
    Float memory = 2
    Int cores = 1
    Int java_mem = floor(memory)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        if [ "~{where_is_umi}" == "R" ]; then       # If UMI is in the READ itself
            if [ "~{umi_paired}" == true ]; then
                fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp ExtractUmisFromBam \
                --molecular-index-tags ZA ZB \
                --single-tag RX \
                --input ~{bam} \
                --read-structure ~{sep=" " read_structure} \
                --output umi_extracted.bam
            else
                fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp ExtractUmisFromBam \
                --molecular-index-tags ZA \
                --single-tag RX \
                --input ~{bam} \
                --read-structure ~{sep=" " read_structure} \
                --output umi_extracted.bam
            fi
        elif [ "~{where_is_umi}" == "N" ]; then     # If UMI is in the READ NAME
            fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp CopyUmiFromReadName \
            --input ~{bam} \
            --output umi_extracted.bam
        elif [ "~{where_is_umi}" == "T" ]; then     # If UMI is already tagged in the BAM
            cp ~{bam} umi_extracted.bam
        fi
    >>>

    output {
        File umi_extracted_bam = "umi_extracted.bam"
    }
}

task sortAndMarkIlluminaAdapters {
    input {
        File bam
    }

    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(4 * data_size)
    Float memory = 2
    Int cores = 1
    Int java_mem = floor(memory)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "broadinstitute/picard:2.27.5"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        # Sort and Mark Illumina Adapters
        /usr/bin/java -Xmx~{java_mem}g -Djava.io.tmpdir=`pwd`/tmp \
        -jar /usr/picard/picard.jar SortSam \
        INPUT=~{bam} \
        OUTPUT=/dev/stdout \
        SORT_ORDER=queryname | \
        /usr/bin/java -Xmx~{java_mem}g -Djava.io.tmpdir=`pwd`/tmp \
        -jar /usr/picard/picard.jar MarkIlluminaAdapters \
        INPUT=/dev/stdin \
        OUTPUT=marked.bam \
        METRICS=adapter_metrics.txt
    >>>

    output {
        File marked_bam = "marked.bam"
    }
}

task umiAlign {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa
        Boolean realign = false
        Int? mem_limit_override
    }

    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")
    Int space_needed_gb = ceil(8 * data_size + reference_size)      # 4 Cores means multiple streams going from /dev/stdin to /dev/stdout
    Float memory = select_first([mem_limit_override, 3])                # No Memory required, everything is streamed straight to output
    Int cores = 4                   # CPU only affects reading in the SAMtoFASTQ for BWA
    Int preemptible = 1
    Int maxRetries = 2

    runtime {
      docker: "mgibio/dna-alignment:1.0.0"
      memory: cores * memory + "GB"
      cpu: cores
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: 10
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String outfile = if realign then "realigned.sorted.bam" else "aligned.bam"

    command <<<
        # 40 min on 4 cores
        if ~{if realign then "true" else "false" }; then
            /bin/bash /usr/bin/umi_realignment.sh ~{bam} ~{reference} ~{cores}
            /usr/bin/java -jar /opt/picard-2.18.1/picard.jar SortSam \
            INPUT=realigned.bam \
            OUTPUT=realigned.sorted.bam \
            SORT_ORDER=queryname
        else
            /bin/bash /usr/bin/umi_alignment.sh ~{bam} ~{reference} ~{cores}
        fi
    >>>

    output {
        File aligned_bam = "~{outfile}"
        #select_first(["aligned.bam", "realigned.sorted.bam"])
        #File aligned_bam_bai = select_first(["aligned.bai", "realigned.bai"])
    }
}

task groupReadsAndConsensus {
    input {
        File bam
        Boolean umi_paired = true
        Int? reads_per_umi_group = 1
    }

    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(8 * data_size)
    Float memory = 4
    Int cores = 1
    Int java_mem = floor(memory)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        if ~{if umi_paired then "true" else "false"}; then
            fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp GroupReadsByUmi --allow-inter-contig false --strategy paired --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input ~{bam} --output umi_grouped.bam
        else
            fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp GroupReadsByUmi --allow-inter-contig false --strategy adjacency --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input ~{bam} --output umi_grouped.bam
        fi
        fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp CallMolecularConsensusReads --sort-order Queryname --input umi_grouped.bam --error-rate-pre-umi 45 --error-rate-post-umi 30 --min-input-base-quality 30 --min-reads ~{reads_per_umi_group} --output consensus_unaligned.bam
    >>>

    output {
        File consensus_bam = "consensus_unaligned.bam"
    }
}

task filterAndClip {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        String sample_name
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5
        Boolean has_umi = true
    }

    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int space_needed_gb = ceil(8 * data_size + reference_size)
    Float memory = reference_size + 4.0                # FilterConsesusReads reads the entire Reference into Memory
    Int cores = 1
    Int java_mem = floor(memory)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        if ~{if has_umi then "true" else "false"}; then
            fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp FilterConsensusReads \
            --input ~{bam} \
            --output consensus_filtered.bam \
            --ref ~{reference} \
            --min-reads ~{sep=" " min_reads} \
            --max-read-error-rate ~{max_read_error_rate} \
            --max-base-error-rate ~{max_base_error_rate} \
            --min-base-quality ~{min_base_quality} \
            --max-no-call-fraction ~{max_no_call_fraction}
            
            fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp ClipBam \
            --input consensus_filtered.bam \
            --ref ~{reference} \
            --clipping-mode Hard \
            --clip-overlapping-reads true \
            --output ~{sample_name}.bam \
            --sort-order Coordinate
        else
            fgbio -Xmx~{java_mem}g --tmp-dir=`pwd`/large_tmp ClipBam \
            --input ~{bam} \
            --ref ~{reference} \
            --clipping-mode Hard \
            --clip-overlapping-reads true \
            --output ~{sample_name}.bam \
            --sort-order Coordinate
        fi
    >>>

    output {
        File clipped_bam = "~{sample_name}.bam"
        File clipped_bai = "~{sample_name}.bai"
    }
}

task processAlignedBAM {
    input {
        File bam
        File bai
        File reference
        File reference_fai
        File reference_dict
        File interval_list
        Array[File] known_sites
        Array[File] known_sites_tbi
        String sample_name
        Boolean apply_bqsr = false
        String input_type = "FASTQ"
    }

    Float data_size = size(bam, "GB")
    Float vcf_size = size(known_sites, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int space_needed_gb = ceil(6 * data_size + reference_size + vcf_size)
    Float memory = 2
    Int cores = 2
    Int java_mem = floor(memory)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.1.8.1"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        # Due to: https://gatk.broadinstitute.org/hc/en-us/community/posts/360075246531-Is-BQSR-accurate-on-Novaseq-6000-
        # Sometimes it is worth skipping this step
        # Applies BQSR on specific intervals defined by the User, if aligned BAM is provided, starts here
        if ~{if apply_bqsr then "true" else "false"}; then
            /gatk/gatk --java-options -Xmx~{java_mem}g BaseRecalibrator -O bqsr.table -L ~{interval_list} -R ~{reference} -I ~{bam} ~{sep=" " prefix("--known-sites ", known_sites)}
            /gatk/gatk --java-options -Xmx~{java_mem}g ApplyBQSR -O ~{sample_name}.bam ~{sep=" " prefix("--static-quantized-quals ", [10, 20, 30])} -R ~{reference} -I ~{bam} -bqsr bqsr.table
            mv ~{sample_name}.bai ~{sample_name}.bam.bai
        else
            if [ "~{input_type}" == "CRAM" ]; then
                /usr/bin/samtools view -@ ~{cores} --fast -b -T ~{reference} -o ~{sample_name}.bam ~{bam}
                /usr/bin/samtools index ~{sample_name}.bam
            else
                cp ~{bam} ~{sample_name}.bam
                cp ~{bai} ~{sample_name}.bam.bai
            fi
        fi
    >>>

    output {
        File initial_bam = "~{sample_name}.bam"
        File initial_bai = "~{sample_name}.bam.bai"
    }

}

task fastQC {
    input {
        File bam
        File bai
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([bam, bai], "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 1
    Int cores = 1

    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bamroot = basename(bam, ".bam")

    command <<<
        /usr/local/bin/fastqc ~{bam} -outdir $PWD
    >>>

    output {
        File fastqc_html = "~{bamroot}_fastqc.html"
        File fastqc = "~{bamroot}_fastqc.zip"
    }
}

task intervalsToBed {
    input {
        File interval_list
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(interval_list, "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 1
    Int cores = 1

    runtime {
        docker: "ubuntu:bionic"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/perl -e '
        use feature qw(say);

        for my $line (<>) {
            chomp $line;
            next if substr($line,0,1) eq q(@); #skip header lines
            my ($chrom, $start, $stop) = split(/\t/, $line);
            say(join("\t", $chrom, $start-1, $stop));
        }' ~{interval_list} > interval_list.bed
    >>>

    output {
        File interval_bed = "interval_list.bed"
    }
}

task createSomalierVcf {
    input {
        File interval_bed
        File chrom_sizes
        File af_only_snp_only_vcf
        File reference
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([interval_bed, chrom_sizes, af_only_snp_only_vcf, reference], "GB")
    Int space_needed_gb = ceil(2 * data_size)
    Float memory = 1
    Int cores = 1

    runtime {
        docker: "kboltonlab/bst:1.0"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        bed_size=$(cat ~{interval_bed} | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
        if (( ${bed_size} > 2*65535 )); then
            cp ~{interval_bed} somalier.bed
        else
            bedtools slop -i ~{interval_bed} -g ~{chrom_sizes} -b 100 > somalier.bed
        fi
        bedtools intersect -a ~{af_only_snp_only_vcf} -b somalier.bed -header > somalier.vcf
        bcftools norm --multiallelics -any -Oz -o somalier.norm.vcf.gz -f ~{reference} somalier.vcf
    >>>

    output {
        File somalier_vcf = "somalier.norm.vcf.gz"
    }
}

task somalier {
    input {
        File somalier_vcf
        File reference
        File bam
        File bai
        String sample_name = "tumor"
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([somalier_vcf, reference, bam, bai], "GB")
    Int space_needed_gb = ceil(3 * data_size)
    Float memory = 2
    Int cores = 1

    runtime {
        docker: "brentp/somalier:latest"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        somalier extract --sites ~{somalier_vcf} -f ~{reference} ~{bam}
    >>>

    output {
        File somalier_out = "~{sample_name}.somalier"
    }
}

task splitBedToChr {
    input {
        File interval_bed
    }

    Float data_size = size(interval_bed, "GB")
    Int space_needed_gb = ceil(data_size)
    Int memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "ubuntu:bionic"
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        intervals=$(awk '{print $1}' ~{interval_bed} | uniq)
        for chr in ${intervals}; do
            grep -w $chr ~{interval_bed} > ~{basename(interval_bed, ".bed")}_${chr}.bed
        done
    >>>

    output {
        Array[File] split_chr = glob(basename(interval_bed, ".bed")+"_*.bed")
    }
}

task mutect {
    input {
        File tumor_bam
        File tumor_bai
        File interval_list
        Boolean tumor_only = true

        File reference
        File reference_fai
        File reference_dict
        
        File? gnomad
        File? gnomad_tbi
        
        File? normal_bam
        File? normal_bai

        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, reference_dict, interval_list], "GB")
    Float data_size = size([tumor_bam, tumor_bai, normal_bam, normal_bai], "GB")
    Int space_needed_gb = ceil(2 * data_size + reference_size)
    Int memory = select_first([mem_limit_override, 6])
    Int cores = select_first([cpu_override, 1])
    Int java_mem = floor(memory)
    Int preemptible = 3
    Int maxRetries = 3

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: cores * memory + "GB"
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = "mutect.filtered.vcf.gz"

    command <<<
        set -e -x -o pipefail
        if ~{if tumor_only then "true" else "false"}; then
            /gatk/gatk Mutect2 --java-options "-Xmx~{java_mem}g" \
            --native-pair-hmm-threads ~{cores} \
                -O mutect.vcf.gz \
                -R ~{reference} \
                -L ~{interval_list} \
                -I ~{tumor_bam} \
                ~{"--germline-resource " + gnomad} \
                --read-index ~{tumor_bai} \
                --f1r2-tar-gz mutect.f1r2.tar.gz \
                --max-reads-per-alignment-start 0
        else
            NORMAL=`samtools view -H ~{normal_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
            TUMOR=`samtools view -H ~{tumor_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`

            /gatk/gatk Mutect2 --java-options "-Xmx~{java_mem}g" \
            --native-pair-hmm-threads ~{cores} \
                -O mutect.vcf.gz \
                -R ~{reference} \
                -L ~{interval_list} \
                -I ~{tumor_bam} \
                --read-index ~{tumor_bai} \
                -tumor "$TUMOR" \
                -I ~{normal_bam} \
                --read-index ~{normal_bai} \
                -normal "$NORMAL" \
                ~{"--germline-resource " + gnomad} \
                --f1r2-tar-gz mutect.f1r2.tar.gz \
                --max-reads-per-alignment-start 0
        fi
        /gatk/gatk LearnReadOrientationModel \
            -I mutect.f1r2.tar.gz \
            -O mutect.read-orientation-model.tar.gz

        /gatk/gatk FilterMutectCalls \
            -R ~{reference} \
            -V mutect.vcf.gz \
            --ob-priors mutect.read-orientation-model.tar.gz \
            -O ~{output_vcf} #Running FilterMutectCalls on the output vcf.
    >>>

    output {
    File vcf = output_vcf
    File vcf_tbi = output_vcf + ".tbi"
    }
}

task vardict {
    input {
        File tumor_bam
        File tumor_bai
        File interval_bed
        Boolean tumor_only = true

        File? normal_bam
        File? normal_bai
        String? normal_sample_name = "NORMAL"

        File reference
        File reference_fai
        
        String tumor_sample_name = "TUMOR"
        Float? min_var_freq = 0.005
        
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, interval_bed], "GB")
    Float data_size = size([tumor_bam, tumor_bai, normal_bam, normal_bai], "GB")
    Int space_needed_gb = ceil(3 * data_size + reference_size)
    Int preemptible = 3
    Int maxRetries = 3
    Int memory = select_first([mem_limit_override, 2])
    Int cores = select_first([cpu_override, 2])

    runtime {
        docker: "kboltonlab/vardictjava:bedtools"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -e -x -o pipefail

        bedtools makewindows -b ~{interval_bed} -w 20250 -s 20000 > ~{basename(interval_bed, ".bed")}_windows.bed
        split -d --additional-suffix .bed -n l/16 ~{basename(interval_bed, ".bed")}_windows.bed splitBed.

        nProcs=~{cores}
        nJobs="\j"

        for fName in splitBed.*.bed; do
            # Wait until nJobs < nProcs, only start nProcs jobs at most
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done

            part=$(echo $fName | cut -d'.' -f2)

            echo ${fName}
            echo ${part}

            /opt/VarDictJava/build/install/VarDict/bin/VarDict \
                -U -G ~{reference} \
                -X 1 \
                -f ~{min_var_freq} \
                -N ~{tumor_sample_name} \
                -b ~{if defined(normal_bam) then "~{tumor_bam}|~{normal_bam}" else "~{tumor_bam}"} \
                -c 1 -S 2 -E 3 -g 4 ${fName} \
                -th ~{cores} \
                --deldupvar -Q 10 -F 0x700 --fisher > result.${part}.txt &
        done;
        # Wait for all running jobs to finish
        wait

        for fName in result.*.txt; do
            cat ${fName} >> resultCombine.txt
        done;

        cat resultCombine.txt | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
            -N ~{if defined(normal_bam) then "~{tumor_sample_name}|~{normal_sample_name}" else "~{tumor_sample_name} -E"} \
            -f ~{min_var_freq} > vardict.vcf

        /usr/bin/bgzip vardict.vcf && /usr/bin/tabix vardict.vcf.gz

        if ~{if tumor_only then "false" else "true"}; then
            # Extract tumor variants before bcbio filter or else both normal and tumor variants will be canidates for the filter
            /usr/bin/bcftools view -h \
                -s ~{tumor_sample_name} \
                --threads ~{cores} vardict.vcf.gz \
                -Oz -o vardict.vcf.gz
            /usr/bin/tabix vardict.vcf.gz
        fi
    >>>

    output {
        File vcf = "vardict.vcf.gz"
        File vcf_tbi = "vardict.vcf.gz.tbi"
    }
}

task lofreq {
    input {
        File tumor_bam
        File tumor_bai
        File interval_bed
        Boolean tumor_only = true

        File? normal_bam
        File? normal_bai

        File reference
        File reference_fai
        
        String tumor_sample_name = "TUMOR"
        Float? min_var_freq = 0.005
        
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bai, normal_bam, normal_bai], "GB")
    Int space_needed_gb = ceil(3 * bam_size + reference_size)
    Int preemptible = 3
    Int maxRetries = 3
    Int memory = select_first([mem_limit_override, 4])
    Int cores = select_first([cpu_override, 1])

    runtime {
        docker: "kboltonlab/lofreq:latest"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        # Lofreq
        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o lofreq.indel.bam ~{tumor_bam}
        samtools index lofreq.indel.bam

        if ~{if tumor_only then "true" else "false"}; then
            /opt/lofreq/bin/lofreq call-parallel --pp-threads ~{cores} -A -B -f ~{reference} --call-indels --bed ~{interval_bed} -o unsorted.lofreq.vcf lofreq.indel.bam
            cat unsorted.lofreq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > lofreq.vcf
        else
            /opt/lofreq/bin/lofreq somatic --call-indels -n ~{normal_bam} -t lofreq.indel.bam -f ~{reference} -l ~{interval_bed} -o lofreq_ --threads ~{cores}
            tabix lofreq_somatic_final.snvs.vcf.gz
            tabix lofreq_somatic_final.indels.vcf.gz
            bcftools concat -a lofreq_somatic_final.snvs.vcf.gz lofreq_somatic_final.indels.vcf.gz > unsorted.lofreq.vcf
            cat unsorted.lofreq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > lofreq.vcf
        fi

        # Reformat
        # Lofreq does not output FORMAT and SAMPLE columns, so we need to reformat the VCF with these columns
        # since many downstream tools require the VCF to be formatted in this way
        grep "##" lofreq.vcf > lofreq.reformat.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">"  >> lofreq.reformat.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> lofreq.reformat.vcf;
        grep -v '#' lofreq.vcf | \
        awk '{
            n=split($8, semi, /;/); 
            sample=""; format=""; 
            for(i in semi){
                split(semi[i], equ, /=/);
                if(equ[1]=="DP" || equ[1]=="AF" || equ[1]=="SB"){
                    format=format equ[1] ":";
                    sample=sample equ[2] ":";
                }
            }
            format=substr(format, 1, length(format)-1);
            sample=substr(sample, 1, length(sample)-1);
            {print $0, "GT:"format, "0/1:"sample}
        }' OFS='\t' >> lofreq.reformat.vcf;
        bgzip lofreq.reformat.vcf && tabix lofreq.reformat.vcf.gz
    >>>

    output {
        File vcf = "lofreq.reformat.vcf.gz"
        File vcf_tbi = "lofreq.reformat.vcf.gz.tbi"
    }
}


task pindel {
    input {
        File tumor_bam
        File tumor_bai
        File region_file
        String? chromosome
        Boolean tumor_only = true

        File? normal_bam
        File? normal_bai
        String? normal_sample_name = "NORMAL"

        File reference
        File reference_fai
        
        String tumor_sample_name = "TUMOR"
        Float? min_var_freq = 0.005
        
        Int? mem_limit_override
        Int? cpu_override

        # Pindel Specific Parameters
        Int insert_size = 400
        String ref_name = "GRCh38DH"
        String ref_date = "20161216"
        Int min_supporting_reads = 3
    }

    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bai, normal_bam, normal_bai], "GB")
    Int space_needed_gb = ceil(3 * bam_size + reference_size)
    Int preemptible = 3
    Int maxRetries = 3
    Int memory = select_first([mem_limit_override, 5])
    Int cores = select_first([cpu_override, 4])

    runtime {
        memory: cores * memory + "GB"
        docker: "mgibio/cle:v1.4.2"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -e -x -o pipefail
        
        if ~{if tumor_only then "true" else "false"}; then
            echo -e "~{tumor_bam}\t~{insert_size}\t~{tumor_sample_name}" > pindel.config
        else
            echo -e "~{normal_bam}\t~{insert_size}\t~{normal_sample_name}" > pindel.config
            echo -e "~{tumor_bam}\t~{insert_size}\t~{tumor_sample_name}" >> pindel.config
        fi

        #bedtools makewindows -b ~{region_file} -w 20250 -s 20000 > ~{basename(region_file, ".bed")}_windows.bed
        split -d --additional-suffix .bed -n l/6 ~{region_file} splitBed.

        time {
            nProcs=~{cores}
            #nJobs="\j"
            nJobs=$(jobs -p | wc -l)

            for fName in splitBed.*.bed; do
                # Wait until nJobs < nProcs, only start nProcs jobs at most
                echo ${nJobs}
                while (( ${nJobs} >= nProcs )); do
                    wait -n
                done

                part=$(echo $fName | cut -d'.' -f2)

                echo ${fName}
                echo ${part}

                # This part takes 1 Gb of memory per core
                /usr/bin/pindel \
                    -i pindel.config \
                    -w 30 \
                    -T ~{cores} \
                    -o ${part} \
                    -f ~{reference} \
                    -j ${fName} \
                    ~{if defined(chromosome) then "-c ~{chromosome}" else ""} &
            done;
            # Wait for all running jobs to finish
            wait
        }

        /bin/cat *_D *_SI *_TD *_LI *_INV | /bin/grep "ChrID" /dev/stdin > pindel.head

        time {
            if ~{if tumor_only then "true" else "false"}; then
                /usr/bin/pindel2vcf -G -p pindel.head -r ~{reference} -R ~{ref_name} -e ~{min_supporting_reads} -d ~{ref_date} -v pindel.vcf
            else
                /usr/bin/perl /usr/bin/write_pindel_filter_config.pl pindel.head ~{reference} $PWD
                /usr/bin/perl /usr/bin/somatic_indelfilter.pl filter.config
                mv pindel.out.vcf pindel.vcf
            fi
        }

        # If pindel returns empty pindel.head file, need to account for empty file.
        is_empty=$(grep "~{tumor_sample_name}" pindel.vcf)
        if [[ ${is_empty} == "" ]]; then
            grep "##" pindel.vcf > temp.vcf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> temp.vcf
            mv temp.vcf pindel.vcf
        fi
    >>>

    output {
        File vcf = "pindel.vcf"
    }
}

task removeEndTags {
    input {
        File pindel_vcf
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(pindel_vcf, "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 1
    Int cores = 1

    runtime {
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst:latest"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/bgzip -c ~{pindel_vcf} >> pindel.vcf.gz
        /usr/local/bin/tabix -p vcf pindel.vcf.gz
        /usr/local/bin/bcftools annotate -x INFO/END -Oz -o pindel.noend.vcf.gz pindel.vcf.gz
        /usr/local/bin/tabix -p vcf pindel.noend.vcf.gz
    >>>

    output {
        File vcf = "pindel.noend.vcf.gz"
        File vcf_tbi = "pindel.noend.vcf.gz.tbi"
    }
}

task sanitizeNormalizeFilter {
    input {
        File vcf
        File vcf_tbi
        File reference
        File reference_fai

        #gnomAD
        File exclude_vcf
        File exclude_vcf_tbi

        #PoN2
        File vcf2PON
        File vcf2PON_tbi
        String caller = "caller"
        String sample_name = "tumor"

        #bcbio
        String filter_string
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([vcf, vcf_tbi, exclude_vcf, exclude_vcf_tbi, vcf2PON, vcf2PON_tbi], "GB")
    Float reference_size = size([reference, reference_fai], "GB")
    Int space_needed_gb = ceil(10 + data_size + reference_size)
    Int memory = 1
    Int cores = 1

    runtime {
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        # BCBIO First for Vardict
        if [[ ~{caller} == "vardict" ]]; then
            bcftools filter -e "~{filter_string}" ~{vcf} --output-type z --output bcbio.vcf.gz -s "BCBIO" -m+
        else
            cp ~{vcf} bcbio.vcf.gz
        fi

        # Santize Step
        # 1) removes lines containing non ACTGN bases, as they conflict with the VCF spec and cause GATK to choke
        # 2) removes mutect-specific format tags containing underscores, which are likewise illegal in the vcf spec
        gunzip -c bcbio.vcf.gz | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" > sanitized.vcf
        bgzip sanitized.vcf && tabix sanitized.vcf.gz
        
        # Normalize Step
        bcftools norm --check-ref w --multiallelics -any --output-type z --output norm.vcf.gz sanitized.vcf.gz -f ~{reference}
        tabix norm.vcf.gz

        if [[ ~{caller} == "pindel" ]]; then
            mv norm.vcf.gz ~{caller}.~{sample_name}.vcf.gz
            tabix ~{caller}.~{sample_name}.vcf.gz
        else
            # gnomAD Intersection
            bcftools isec -C -w1 norm.vcf.gz ~{exclude_vcf} --output-type z --output isec.vcf.gz

            # PoN2
            export name=~{caller}.~{sample_name}.vcf.gz
            printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
            printf "##INFO=<ID=PON_NAT2_percent,Number=1,Type=Integer,Description=\"Number of samples with variant at >=2 percent\">\n" >> pon2.header;
            printf "##INFO=<ID=PON_MAX_VAF,Number=1,Type=Float,Description=\"The maximum VAF found in the PoN Samples\">\n" >> pon2.header;
            bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\t%INFO/NS\t%INFO/max_VAF\n" ~{vcf2PON} > normal2.txt
            bgzip -f normal2.txt
            tabix -f -s1 -b2 -e2 normal2.txt.gz
            bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent,PON_NAT2_percent,PON_MAX_VAF isec.vcf.gz -Oz -o $name
            tabix $name
        fi
    >>>

    output {
        File annotated_vcf = "~{caller}.~{sample_name}.vcf.gz"
        File annotated_vcf_tbi = "~{caller}.~{sample_name}.vcf.gz.tbi"
    }
}

# The fake VCF only really needs the basic information, FORMAT ans SAMPLE columns can be empty since
# it will be overwritten by the annotations and filtering
task createFakeVcf {
    input {
        File vcf
        String tumor_sample_name
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(vcf, "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 1
    Int cores = 1

    runtime {
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -o errexit
        set -o nounset

        echo -e "##fileformat=VCFv4.2" > fake.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> fake.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> fake.vcf;
        zcat ~{vcf} | grep -v '#' | awk '{print $1, $2, $3, $4, $5, $6, "PASS\t.\tGT\t0/1"}' OFS='\t' >> fake.vcf;
        bgzip fake.vcf && tabix fake.vcf.gz
    >>>

    output {
        File fake_vcf = "fake.vcf.gz"
        File fake_vcf_tbi = "fake.vcf.gz.tbi"
    }
}

task mergeVcf {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(vcfs, "GB")
    Int space_needed_gb = ceil(2 * data_size)
    Float memory = 1
    Int cores = 1

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_file = merged_vcf_basename + ".vcf.gz"

    command <<<
        /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}
        /usr/local/bin/tabix ~{output_file}
    >>>

    output {
        File merged_vcf = output_file
        File merged_vcf_tbi = "~{output_file}.tbi"
    }
}

task fpFilter {
    input {
        File reference
        File reference_fai
        File bam
        File vcf
        String output_vcf_basename = "fpfilter"
        String sample_name = "TUMOR"
        Float? min_var_freq = 0.05
        Int? mem_limit_override
    }

    Float data_size = size(vcf, "GB")
    Float reference_size = size([reference, reference_fai, bam], "GB")
    Int space_needed_gb = ceil(2 * (data_size + reference_size))
    Float memory = 1
    Int cores = 4
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: cores * memory + "GB"
        docker: "kboltonlab/fp_filter-wdl"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = output_vcf_basename + ".vcf"

    command <<<
        zgrep -v '#' ~{vcf} | split -d --additional-suffix .vcf -l 2000 - splitVCF.
        for vcf in splitVCF.*.vcf; do
            zgrep '#' ~{vcf} > tmp_file
            cat "$vcf" >> tmp_file
            mv -f tmp_file "$vcf"
        done

        nProcs=~{cores}
        #nJobs="\j"
        nJobs=$(jobs -p | wc -l)

        for vcf in splitVCF.*.vcf; do
            # Wait until nJobs < nProcs, only start nProcs jobs at most
            echo ${nJobs}
            while (( ${nJobs} >= nProcs )); do
                wait -n
            done
            fBase=$(echo $vcf | cut -d '.' -f2)
            echo ${fBase}
            /usr/bin/perl /usr/bin/fpfilter.pl --bam-readcount /usr/bin/bam-readcount \
                --samtools /opt/samtools/bin/samtools --output fpfilter.${fBase}.vcf --reference ~{reference} \
                --bam-file ~{bam} --vcf-file ${vcf} --sample ~{sample_name} --min-var-freq ~{min_var_freq} &
        done;
        wait

        zgrep '#' ~{vcf} > ~{output_vcf}
        for vcf in fpfilter.*.vcf; do
            zgrep -v '#' ${vcf} >> ~{output_vcf}
        done
        /usr/bin/bgzip ~{output_vcf} && /usr/bin/tabix ~{output_vcf}.gz
    >>>

    output {
        File filtered_vcf = "~{output_vcf}.gz"
        File filtered_vcf_tbi = "~{output_vcf}.gz.tbi"
    }
}

task mskGetBaseCounts {
    input {
        File reference
        File reference_fai
        Pair[File, File] normal_bam
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
        Int? mem_limit_override
    }

    Float reference_size = size([reference, reference_fai], "GB")
    Float data_size = size([normal_bam.left, normal_bam.right, vcf], "GB")
    Int space_needed_gb = ceil(2 * (data_size + reference_size))
    Float memory = select_first([mem_limit_override, 4])
    Int cores = 4
    Int preemptible = 1
    Int maxRetries = 2              # This task is prone to failing due to memory issues, so we allow for retries

    runtime {
      docker: "duct/getbasecount:latest"
      cpu: cores
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} HDD"
      bootDiskSizeGb: 10
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{data_size}"
        echo "SPACE_NEEDED: ~{space_needed_gb}"

        sample_name=$(samtools view -H ~{normal_bam.left} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
        if [[ $(zgrep -v '#' ~{vcf} | wc -l) -lt 1 ]]; then
            echo "PRINT: ~{vcf}"
            printf "##fileformat=VCFv4.2\n" > ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth matching reference (REF) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth matching alternate (ALT) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant frequence (AD/DP)\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPP,Number=1,Type=Integer,Description=\"Depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPN,Number=1,Type=Integer,Description=\"Depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDP,Number=1,Type=Integer,Description=\"Reference depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDN,Number=1,Type=Integer,Description=\"Reference depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADP,Number=1,Type=Integer,Description=\"Alternate depth on postitive strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADN,Number=1,Type=Integer,Description=\"Alternate depth on negative strand\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Total fragment depth\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=RDF,Number=1,Type=Float,Description=\"Fragment depth matching reference (REF) allele\">\n" >> ~{pon_final_name}.vcf
            printf "##FORMAT=<ID=ADF,Number=1,Type=Float,Description=\"Fragment depth matching alternate (ALT) allele\">\n" >> ~{pon_final_name}.vcf
            printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${sample_name}\n" >> ~{pon_final_name}.vcf
        else
            echo "SCRIPT: ~{normal_bam.left}"
            if [[ ~{vcf} == *.vcf.gz ]]; then
                bgzip -d ~{vcf}
                vcf_file=~{vcf}
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.left} --vcf "${vcf_file%.*}" --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread ~{cores} --max_block_dist 10000
            else
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.left} --vcf ~{vcf} --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread ~{cores} --max_block_dist 10000
            fi
        fi
        bgzip ~{pon_final_name}.vcf && tabix ~{pon_final_name}.vcf.gz
    >>>

    output {
        File pileup = "~{pon_final_name}.vcf.gz"
        File pileup_tbi = "~{pon_final_name}.vcf.gz.tbi"
    }
}

task bcftoolsMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
        Boolean RD_AD = false
    }

    Float data_size = size(vcfs, "GB")
    Int space_needed_gb = ceil(2 * data_size)
    Float memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_file = merged_vcf_basename + ".vcf.gz"

    command <<<
        /usr/local/bin/bcftools merge --output-type z -o ~{output_file} ~{sep=" " vcfs}
        /usr/local/bin/tabix ~{output_file}

        if ~{if RD_AD then "true" else "false"}; then
            bcftools +fill-tags -Ov ~{output_file} -- -t "PON_RefDepth=sum(RD)" | \
            bcftools +fill-tags -Oz -o pileup.vcf.gz - -- -t "PON_AltDepth=sum(AD)" && tabix pileup.vcf.gz
        fi
    >>>

    output {
        File merged_vcf = output_file
        File merged_vcf_tbi = "~{output_file}.tbi"
        File? pileup_vcf = "pileup.vcf.gz"
        File? pileup_vcf_tbi = "pileup.vcf.gz.tbi"
    }
}

task vep {
    input {
        File vcf
        File reference
        File reference_fai
        
        # Vep Stuff
        File cache_dir_zip
        Array[String] plugins = ["Frameshift", "Wildtype"]
        Boolean coding_only = false
        File clinvar
        File clinvar_tbi
        File? synonyms_file
        Int? mem_limit_override
    }

    Float cache_size = 3*size(cache_dir_zip, "GB")  # doubled to unzip
    Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
    Float reference_size = size([reference, reference_fai], "GB")
    Int space_needed_gb = ceil(cache_size + vcf_size + reference_size + size(synonyms_file, "GB"))
    Float memory = select_first([mem_limit_override, 1]) # We want the base to be around 6
    Int cores = 4
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: cores * memory + "GB"
        cpu: cores
        docker: "ensemblorg/ensembl-vep:release_109.3"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String outfile = basename(basename(vcf, ".gz"), ".vcf") + ".VEP_annotated.vcf"
    String cache_dir = basename(cache_dir_zip, ".zip")

    command <<<
        #mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}
        unzip -qq ~{cache_dir_zip}

         /opt/vep/src/ensembl-vep/vep -i ~{vcf} --vcf -o ~{outfile}  \
            --cache --offline --dir_cache ~{cache_dir}/VepData/ --merged --assembly GRCh38 --use_given_ref --species homo_sapiens \
            --symbol --transcript_version --everything --check_existing \
            ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
            ~{if coding_only then "--coding_only" else ""} \
            --buffer_size 1000  --fork ~{cores} \
            --fasta ~{reference} \
            --pick --pick_order canonical,rank,mane_select,mane_plus_clinical,ccds,appris,tsl,biotype,length \
            --force_overwrite  --no_stats \
            --custom ~{clinvar},clinvar,vcf,exact,0,ID,AF_ESP,AF_EXAC,AF_TGP,CLNSIG,CLNSIGCONF,CLNDN,ORIGIN \
            --dir_plugins ~{cache_dir}/plugin/ \
            ~{sep=" " prefix("--plugin ", plugins)} \
            --plugin CADD,~{cache_dir}/CADD/whole_genome_SNVs.tsv.gz,~{cache_dir}/CADD/gnomad.genomes.r3.0.indel.tsv.gz \
            --plugin REVEL,~{cache_dir}/REVEL/new_tabbed_revel_grch38.tsv.gz \
            --plugin SpliceAI,snv=~{cache_dir}/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=~{cache_dir}/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
            --plugin pLI,~{cache_dir}/pLI/plI_gene.txt

            bgzip ~{outfile} && tabix ~{outfile}.gz
    >>>

    output {
        File annotated_vcf = "~{outfile}.gz"
        File annotated_vcf_tbi = "~{outfile}.gz.tbi"
    }
}

task normalFisher {
    input {
        File vcf
        File pon
        File pon_tbi
        String caller = "caller"
        String? p_value = "0.05"
    }


    Float data_size = size([vcf, pon, pon_tbi], "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 2
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} HDD"
      bootDiskSizeGb: 10
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        if [[ "~{vcf}" == *.gz ]]; then
            name=$(basename ~{vcf} .vcf.gz)
        else
            name=$(basename ~{vcf} .vcf)
        fi

        printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
        printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

        sample=`bcftools query -l ~{vcf}`
        bcftools view -H ~{vcf} | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
        bgzip $sample.name;
        tabix $sample.name.gz -s1 -b2 -e2;
        bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE ~{vcf} -Oz -o $name.sample.vcf.gz && tabix $name.sample.vcf.gz;

        ## Lofreq has DP4 which splits into RefFwd, RefRev, AltFwd, AltRev
        if [ "~{caller}" == "lofreq" ]; then
            echo '
            #!/usr/bin/env Rscript
            args = commandArgs(trailingOnly=TRUE)

            df <- read.delim(args[1], header = FALSE, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            # Remember, Lofreq splits DP4 into RefFwd, RefRev and AltFwd, AltRev so technically ref = x[3] + x[4] and alt = x[5] + x[6]
            ref = x[3] + x[4]
            alt = x[5] + x[6]
            if ((x[1]+x[2]==0) | (ref+alt==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) & (ref==0 & alt!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= alt/(ref+alt)) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], ref, alt), ncol=2))$p.value)
            }
            })
            df$V2 <- sprintf("%1.0f", df$V2)
            write.table(df[, -c(9:10)], file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a ~{pon} -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        else
            echo '
            #!/usr/bin/env Rscript
            args = commandArgs(trailingOnly=TRUE)
            
            df <- read.delim(args[1], header = FALSE, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]
            
            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            if ((x[1]+x[2]==0) | (x[3]+x[4]==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) & (x[3]==0 & x[4]!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+x[4])) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
            }
            })
            df$V2 <- sprintf("%1.0f", df$V2)
            write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a ~{pon} -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        fi
        chmod u+x fisherTestInput.R

        # Depending on how we split, we might have caller_vcf that doesn't have any variants called
        if [ -s $name.fisher.input ]; then
            LC_ALL=C.UTF-8 Rscript --vanilla ./fisherTestInput.R $name.fisher.input $name.fisher.output
            bgzip -f $name.fisher.output
            tabix -f -s1 -b2 -e2 $name.fisher.output.gz
            bcftools annotate -a $name.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $name.sample.pileup.vcf.gz -Oz -o $name.pileup.fisherPON.vcf.gz && tabix $name.pileup.fisherPON.vcf.gz
            bcftools filter -i "INFO/PON_FISHER<=~{p_value}" $name.pileup.fisherPON.vcf.gz -Oz -o $name.filtered.pileup.fisherPON.vcf.gz && tabix $name.filtered.pileup.fisherPON.vcf.gz
        else
            bcftools annotate -h fisher.header $name.sample.pileup.vcf.gz -Oz -o $name.pileup.fisherPON.vcf.gz && tabix $name.pileup.fisherPON.vcf.gz
            bcftools filter -i "INFO/PON_FISHER<=~{p_value}" $name.pileup.fisherPON.vcf.gz -Oz -o $name.filtered.pileup.fisherPON.vcf.gz && tabix $name.filtered.pileup.fisherPON.vcf.gz
        fi
    >>>

    output {
        File pon_vcf = select_first([basename(vcf, ".vcf.gz") + ".pileup.fisherPON.vcf.gz",basename(vcf, ".vcf") + ".pileup.fisherPON.vcf.gz"])
        File pon_vcf_tbi = select_first([basename(vcf, ".vcf.gz") + ".pileup.fisherPON.vcf.gz.tbi",basename(vcf, ".vcf") + ".pileup.fisherPON.vcf.gz.tbi"])
        File pon_filtered_vcf = select_first([basename(vcf, ".vcf.gz") + ".filtered.pileup.fisherPON.vcf.gz",basename(vcf, ".vcf") + ".filtered.pileup.fisherPON.vcf.gz"])
        File pon_filtered_vcf_tbi = select_first([basename(vcf, ".vcf.gz") + ".filtered.pileup.fisherPON.vcf.gz.tbi",basename(vcf, ".vcf") + ".filtered.pileup.fisherPON.vcf.gz.tbi"])
    }
}

task annotateVcf {
    input {
        File vcf
        File vcf_tbi
        File fp_filter
        File fp_filter_tbi
        File vep
        File vep_tbi
        String caller_prefix
        String sample_name
    }

    Float data_size = size([vcf, vcf_tbi, fp_filter, fp_filter_tbi, vep], "GB")
    Int space_needed_gb = ceil(3 * data_size)
    Float memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      memory: cores * memory + "GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} HDD"
      bootDiskSizeGb: 10
      cpu: cores
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        #zcat ~{fp_filter} | grep '##' | tail -n +4  > fp_filter.header;
        zcat ~{fp_filter} | grep '##FILTER' > fp_filter.header;
        zcat ~{fp_filter} | grep -v '#' > fp_filter.results;
        zcat ~{vep} | grep '##' | tail -n +3 > vep.header;

        printf "##INFO=<ID=FP_filter,Number=.,Type=String,Description=\"Result from FP Filter\">" >> fp_filter.header;
        sed -i 's/;/,/g' fp_filter.results

        bgzip -f fp_filter.results
        tabix -f -s1 -b2 -e2 fp_filter.results.gz

        bcftools annotate -a fp_filter.results.gz -h fp_filter.header -c CHROM,POS,ID,REF,ALT,-,FP_filter ~{vcf} -Oz -o ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz
        tabix ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz
        bcftools annotate -a ~{vep} -h vep.header -c CSQ ~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz -Oz -o ~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz
        tabix ~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz
    >>>

    output {
        File final_annotated_vcf = "~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz"
        File final_annotated_vcf_tbi = "~{caller_prefix}.~{sample_name}.final.annotated.vcf.gz.tbi"
        File fpfilter_annotated_vcf = "~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz"
        File fpfilter_annotated_vcf_tbi = "~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz.tbi"
    }
}

task annotatePD {
    input {
        File mutect_vcf
        File lofreq_vcf
        File vardict_vcf
        File bolton_bick_vars
        File mut2_bick
        File mut2_kelly
        File matches2
        File truncating
        File gene_list
        File oncokb_genes
        File cosmic_dir_zip
        String? pon_pvalue = "2.114164905e-6"
    }

    Float caller_size = size([mutect_vcf, lofreq_vcf, vardict_vcf], "GB")
    Float file_size = size([bolton_bick_vars, mut2_bick, mut2_kelly, matches2, truncating, gene_list, oncokb_genes], "GB")
    Float cosmic_size = 6*size(cosmic_dir_zip, "GB")        #cosmic_zip.zip is 194MB, but unzipped is 895MB
    Int space_needed_gb = ceil(10 + caller_size + file_size + cosmic_size)
    Float memory = 6
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/r_docker_ichan:latest"
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} HDD"
      bootDiskSizeGb: 10
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String cosmic_dir = basename(cosmic_dir_zip, ".zip") + "/"

    command <<<
        set -eou pipefail

        unzip -qq ~{cosmic_dir_zip}

        CSQ_string=$(zgrep 'Ensembl VEP. Format:' ~{mutect_vcf} | cut -d":" -f2)
        CSQ_string="${CSQ_string%??}"
        CSQ_string="${CSQ_string// /}"

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/annotate/ArCHAnnotationScript.R --input ~{mutect_vcf} --out $(basename ~{mutect_vcf} .vcf.gz) --caller mutect \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --truncating ~{truncating} \
        --gene_list ~{gene_list} \
        --oncokb_genes ~{oncokb_genes} \
        --cosmic_dir ~{cosmic_dir} \
        --p_value ~{pon_pvalue} \
        --csq_string ${CSQ_string}
        echo "Mutect AnnotatePD Finished..."

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/annotate/ArCHAnnotationScript.R --input ~{lofreq_vcf} --out $(basename ~{lofreq_vcf} .vcf.gz) --caller lofreq \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --truncating ~{truncating} \
        --gene_list ~{gene_list} \
        --oncokb_genes ~{oncokb_genes} \
        --cosmic_dir ~{cosmic_dir} \
        --p_value ~{pon_pvalue} \
        --csq_string ${CSQ_string}
        echo "Lofreq AnnotatePD Finished..."

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/annotate/ArCHAnnotationScript.R --input ~{vardict_vcf} --out $(basename ~{vardict_vcf} .vcf.gz) --caller vardict \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --truncating ~{truncating} \
        --gene_list ~{gene_list} \
        --oncokb_genes ~{oncokb_genes} \
        --cosmic_dir ~{cosmic_dir} \
        --p_value ~{pon_pvalue} \
        --csq_string ${CSQ_string}
        echo "Vardict AnnotatePD Finished..."
    >>>

    output {
        File mutect_vcf_annotate_pd = basename(mutect_vcf, ".vcf.gz") + ".tsv"
        File lofreq_vcf_annotate_pd = basename(lofreq_vcf, ".vcf.gz") + ".tsv"
        File vardict_vcf_annotate_pd = basename(vardict_vcf, ".vcf.gz") + ".tsv"
    }
}

task combine_all {
    input {
        File mutect_tsv
        File lofreq_tsv
        File vardict_tsv
        File pindel_vcf
        File pon
        String? pon_pvalue = "2.114164905e-6"
        String tumor_sample_name
        Boolean model = false
    }

    Float data_size = size([mutect_tsv, lofreq_tsv, vardict_tsv, pindel_vcf, pon], "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "kboltonlab/r_docker_ichan:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/combine/ArCHCombine.R \
        --mutect ~{mutect_tsv} \
        --lofreq ~{lofreq_tsv} \
        --vardict ~{vardict_tsv} \
        --pindel ~{pindel_vcf} \
        --pon ~{pon} \
        --p_value ~{pon_pvalue} \
        --sample ~{tumor_sample_name} \
        --model ~{model}
    >>>

    output {
        File? model_output = "~{tumor_sample_name}.model.tsv"
        File final_output = "~{tumor_sample_name}.final.annotated.tsv"
    }
}

# Generates an output file that combines all repeated columns (features that are not unique to each caller)
# Adds additional features and runs select features through a trained (on ArcherDX Samples) XGBoost Model.
# False positive filters are run to tag additional artifacts from the final list,

# Additionally performs complex variant matching by looking for complex variants called by Vardict and matching possible
# components found in Mutect or Lofreq. If complex variants are identified they are matched with Pindel as a final stop
# to determine the validity of real complex variants.

# Additional features are added such as the Z-score for the number of reference reads for all variants in a caller.
# The Z-score tells us whether the alignment step correctly worked for all positions as one would except the reference depth to be fairly consistent for all variants.
# An outlier indicates that this region may be over or under aligned. Another feature is the Fisher’s test with the forward and reverse reference / alternate reads.
# A significant value here would indicate that the forward and reverse strands are not consistent and strand bias is occurring. We were not able to calculate a cutoff p-value so this feature is not used in our analysis.

# Finally a predicted probabiliy from the XGBoost Model is calculated for each variant
task xgb_model {
    input {
        File mutect_tsv
        File lofreq_tsv
        File vardict_tsv
        File pindel_full_vcf
        File pon
        String? pon_pvalue = "2.114164905e-6"
        Boolean model = false
        String tumor_sample_name
    }

    Float data_size = size([mutect_tsv, lofreq_tsv, vardict_tsv, pindel_full_vcf, pon], "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/xgb:ic_patch"
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} HDD"
      bootDiskSizeGb: 10
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        if [ ~{model} == true ]; then
            /opt/bin/xgbappcompiled/bin/xgbapp ~{lofreq_tsv} ~{mutect_tsv} ~{vardict_tsv} ~{pindel_full_vcf} ~{pon} --pvalue ~{pon_pvalue}
        else
            /opt/bin/xgbappcompiled/bin/xgbapp ~{lofreq_tsv} ~{mutect_tsv} ~{vardict_tsv} ~{pindel_full_vcf} ~{pon} --pvalue ~{pon_pvalue} --nomodel
        fi
        echo "Model Finished..."
    >>>


    output {
        File model_output = "output_~{tumor_sample_name}.tsv.gz"
        File model_raw_output = "output_~{tumor_sample_name}.raw.tsv.gz"
        File mutect_complex = "output_mutect_complex_~{tumor_sample_name}.tsv.gz"
        File pindel_complex = "output_pindel_complex_~{tumor_sample_name}.tsv.gz"
        File lofreq_complex = "output_lofreq_complex_~{tumor_sample_name}.tsv.gz"
        File caller_filters = "Caller_Filters.raw.tsv.gz"
    }
}
