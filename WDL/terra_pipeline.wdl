version 1.0

# Final WDL Pipeline for TERRAbio
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

# The SpliceAI Plugin requires two files be provided, so rather than having 4 files passed individually
# it makes things cleaner to have a single structure hold all four.
struct VepSpliceAIPlugin {
    File? spliceAI_snv
    File? spliceAI_snv_tbi
    File? spliceAI_indel
    File? spliceAI_indel_tbi
}

# Main Workflow
workflow boltonlab_CH {
    input {
        # Sequence Information
        String platform = "ArcherDX"
        String platform_unit = "Illumina"
        String library = "LIBRARY"
        File fastq_one
        File fastq_two
        Boolean bam_input = false           # Default input is FASTQ files, but unaligned BAM is preferred
        File unaligned_bam = ""
        Boolean? aligned = false            # If input is an already aligned BAM file then set this flag
        File? aligned_bam_file
        File? aligned_bam_file_bai
        Array[String] read_structure        # Used for the UMI processing see: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
        String tumor_sample_name
        File target_intervals               # Interval List
        Int? mem_limit_override = 6         # Some applications will require more memory depending on BAM size and BED size... (in GB)
                                            # Need to account for these types of errors

        # Reference
        File reference
        File reference_fai
        File reference_dict
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_pac
        File reference_sa

        # FASTQ Preprocessing
        Boolean has_umi = true
        Boolean? umi_paired = true          # If the UMI is paired (R1 and R2) then set this flag
        Int umi_length = 8
        Array[Int] min_reads = [1]          # The minimum number of reads that constitutes a "read family"
        Float? max_read_error_rate = 0.05   # If 5% of the reads within a "read family" does not match, removes the family
        Float? max_base_error_rate = 0.1    # If 10% of the bases are different than the majority, masks the base with N
        Int min_base_quality = 1            # Any base less than QUAL 1 will be masked with N
        Float max_no_call_fraction = 0.5    # Maximum fraction of no-calls (N) in the read after filtering

        # BQSR
        Array[File] bqsr_known_sites
        Array[File] bqsr_known_sites_tbi
        Array[String] bqsr_intervals

        # QC
        # These are essentially coordinates for regions you were able to design probes for in the reagent.
        # Typically the reagent provider has this information available in bed format and it can be converted to an interval_list with Picard BedToIntervalList.
        # Astrazeneca also maintains a repo of baits for common sequencing reagents available at https://github.com/AstraZeneca-NGS/reference_data
        File bait_intervals
        Array[LabelledFile] per_base_intervals      # If QC needs to be done at a per-base resolution
        Array[LabelledFile] per_target_intervals    # If QC needs to be done at a per-target resolution
        Array[LabelledFile] summary_intervals       # If QC needs to be done for specific intervals
        File omni_vcf                               # The omni VCF is a list of sites used by verifyBamId for identifying contamination (really mixing of multiple samples)
        File omni_vcf_tbi
        String picard_metric_accumulation_level
        Int? qc_minimum_mapping_quality = 0
        Int? qc_minimum_base_quality = 0
        File chrom_sizes
        File af_only_snp_only_vcf

        # Variant Calling
        Array[Pair[File, File]] pon_bams            # List of BAMs within the Panel of Normals (PoN)
        Float? af_threshold = 0.0001                # Minimum VAF Cut-Off
        String? pon_pvalue = "2.114164905e-6"       # Bonferroni Corrected P-Value for Significance

        # Normal BAMs
        Boolean tumor_only = true                   # Defines if Normal BAMs should be used for Variant Calling
        File normal_bam                             # TODO: Implement Normal
        File normal_bam_bai

        # Pindel
        Int pindel_insert_size = 400
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? pindel_min_supporting_reads = 3

        # See: http://bcb.io/2016/04/04/vardict-filtering/
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
        File TSG_file
        File oncoKB_curated
        File pd_annotation_file
        File truncating
        File cosmic_dir_zip

        # VEP Parameters
        File vep_cache_dir_zip                      # WDL does not have a Directory Variable, so the entire cache needs to be ZIP
        String vep_ensembl_assembly
        String vep_ensembl_version
        String vep_ensembl_species
        Array[String] vep_plugins = ["Frameshift", "Wildtype"]
        VepSpliceAIPlugin vep_plugin_spliceAI_files = {}
        File? synonyms_file
        Boolean? annotate_coding_only = true
        Array[VepCustomAnnotation] vep_custom_annotations
        String vep_pick = "pick"
        Boolean everything = true

        # Variants within gnomAD that have a VAF higher than 0.5%
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi

        # Used for Pilot Data. Leave FALSE
        Boolean? pilot = false
    }

    # If the BAM file is already aligned and consensus sequencing was done, then alignment can be skipped
    if (!aligned) {
        # If the input is BAM, we need to check prune reads that have bad UMI information
        if (platform == 'ArcherDX') {
            if (bam_input) {
                call bamToFastq {
                    input: unaligned_bam = unaligned_bam
                }
            }

            # Archer UMIs are not typical, they have a 13 bp Adapter after their 8 bp UMI
            # There needs to be some pruning involved before we can extract the UMIs
            call filterArcherUmiLength as filterUMI {
                input:
                fastq1 = select_first([bamToFastq.fastq_one, fastq_one]),
                fastq2 = select_first([bamToFastq.fastq_two, fastq_two]),
                umi_length = umi_length
            }
            call bbmapRepair as repair {
                input:
                fastq1 = filterUMI.fastq1_filtered,
                fastq2 = filterUMI.fastq2_filtered
            }
            call fastqToBam as archer_fastq_to_bam {
                input:
                fastq1 = repair.fastq1_repair,
                fastq2 = repair.fastq2_repair,
                sample_name = tumor_sample_name,
                library_name = library,
                platform_unit = platform_unit,
                platform = platform
            }
        }

        # Since we only need to apply the above for Archer
        if (!bam_input) {
            call fastqToBam as fastq_to_bam {
                input:
                fastq1 = fastq_one,
                fastq2 = fastq_two,
                sample_name = tumor_sample_name,
                library_name = library,
                platform_unit = platform_unit,
                platform = platform
            }
        }

        if (has_umi) {
            # Removes UMIs from Reads and adds them as RX Tag
            if (platform != 'MGI') {
                call extractUmis {
                    input:
                    bam = select_first([archer_fastq_to_bam.bam, fastq_to_bam.bam, unaligned_bam]),
                    read_structure = read_structure,
                    umi_paired = umi_paired
                }
                # Mark Adapters
                call markIlluminaAdapters as markAdapters{
                    input:
                    bam = extractUmis.umi_extracted_bam
                }
            }

            if (platform == 'MGI') {
                call markIlluminaAdapters as markAdapters_MGI{
                    input:
                    bam = select_first([fastq_to_bam.bam, unaligned_bam])
                }
            }

            # First Alignment
            call umiAlign as align {
                input:
                bam = select_first([markAdapters.marked_bam, markAdapters_MGI.marked_bam]),
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

        if (!has_umi) {
            # Mark Adapters
            call markIlluminaAdapters as markAdapters_NoUMI{
                input:
                bam = select_first([archer_fastq_to_bam.bam, fastq_to_bam.bam, unaligned_bam]),
            }
        }

        # Realign the Consensus Called Reads
        call realign {
            input:
            bam = select_first([groupReadsAndConsensus.consensus_bam, markAdapters_NoUMI.marked_bam]),
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            reference_amb = reference_amb,
            reference_ann = reference_ann,
            reference_bwt = reference_bwt,
            reference_pac = reference_pac,
            reference_sa = reference_sa
        }

        if (has_umi) {
            # Filter, Clip, and Collect QC Metrics
            call filterClipAndCollectMetrics {
                input:
                bam = realign.consensus_aligned_bam,
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                min_reads = min_reads,
                max_read_error_rate = max_read_error_rate,
                max_base_error_rate = max_base_error_rate,
                min_base_quality = min_base_quality,
                max_no_call_fraction = max_no_call_fraction,
                target_intervals = target_intervals,
                description = tumor_sample_name,
                umi_paired = umi_paired
            }
        }

        if (!has_umi) {
            # Filter, Clip, and Collect QC Metrics
            call clipAndCollectMetrics {
                input:
                bam = realign.consensus_aligned_bam,
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                target_intervals = target_intervals,
                description = tumor_sample_name,
                umi_paired = umi_paired
            }
        }
    }

    # Applies BQSR on specific intervals defined by the User, if aligned BAM is provided, starts here
    call bqsrApply as bqsr {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bam = select_first([filterClipAndCollectMetrics.clipped_bam, aligned_bam_file, clipAndCollectMetrics.clipped_bam]),
        bam_bai = select_first([filterClipAndCollectMetrics.clipped_bam_bai, aligned_bam_file_bai, clipAndCollectMetrics.clipped_bam_bai]),
        intervals = bqsr_intervals,
        known_sites = bqsr_known_sites,
        known_sites_tbi = bqsr_known_sites_tbi,
        output_name = tumor_sample_name
    }

    # Obtains Alignment Metrics and Insert Size Metrics
    call Metrics as collectAllMetrics {
        input:
        bam = bqsr.bqsr_bam,
        bam_bai = bqsr.bqsr_bam_bai,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        metric_accumulation_level = picard_metric_accumulation_level
    }

    # Collects QC for various levels defined by the User
    call collectHsMetrics as collectRoiHsMetrics {
        input:
        bam = bqsr.bqsr_bam,
        bam_bai = bqsr.bqsr_bam_bai,
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        metric_accumulation_level = "ALL_READS",
        bait_intervals = bait_intervals,
        target_intervals = target_intervals,
        per_target_coverage = false,
        per_base_coverage = false,
        output_prefix = "roi",
        minimum_mapping_quality=qc_minimum_mapping_quality,
        minimum_base_quality=qc_minimum_base_quality
    }

    scatter(interval in summary_intervals) {
        call collectHsMetrics as collectSummaryHsMetrics{
            input:
            bam = bqsr.bqsr_bam,
            bam_bai = bqsr.bqsr_bam_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_intervals = interval.file,
            target_intervals = interval.file,
            output_prefix = "summary-~{interval.label}",
            minimum_mapping_quality = qc_minimum_mapping_quality,
            minimum_base_quality = qc_minimum_base_quality,
            metric_accumulation_level = "ALL_READS",
            per_target_coverage = false,
            per_base_coverage = false
        }
    }

    scatter(interval in per_base_intervals) {
        call collectHsMetrics as collectPerBaseHsMetrics {
            input:
            bam = bqsr.bqsr_bam,
            bam_bai = bqsr.bqsr_bam_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_intervals = interval.file,
            target_intervals = interval.file,
            output_prefix = "base-~{interval.label}",
            minimum_mapping_quality = qc_minimum_mapping_quality,
            minimum_base_quality = qc_minimum_base_quality,
            metric_accumulation_level = "ALL_READS",
            per_target_coverage = false,
            per_base_coverage = true
        }
    }

    scatter(interval in per_target_intervals) {
        call collectHsMetrics as collectPerTargetHsMetrics{
            input:
            bam = bqsr.bqsr_bam,
            bam_bai = bqsr.bqsr_bam_bai,
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bait_intervals = interval.file,
            target_intervals = interval.file,
            output_prefix = "target-~{interval.label}",
            minimum_mapping_quality = qc_minimum_mapping_quality,
            minimum_base_quality = qc_minimum_base_quality,
            metric_accumulation_level = "ALL_READS",
            per_target_coverage = true,
            per_base_coverage = false
        }
    }

    call samtoolsFlagstat {
        input:
        bam = bqsr.bqsr_bam,
        bam_bai = bqsr.bqsr_bam_bai,
    }

    # Uses the OMNI vcf to calculate possible contamination within the Sample
    call selectVariants {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        vcf = omni_vcf,
        vcf_tbi = omni_vcf_tbi,
        interval_list = target_intervals
    }

    call verifyBamId {
        input:
        bam = bqsr.bqsr_bam,
        bam_bai = bqsr.bqsr_bam_bai,
        vcf=selectVariants.filtered_vcf
    }

    call fastQC {
        input:
        bam = bqsr.bqsr_bam,
        bam_bai = bqsr.bqsr_bam_bai
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
        bam = bqsr.bqsr_bam,
        bam_bai = bqsr.bqsr_bam_bai,
        sample_name = tumor_sample_name
    }

    # In order to parallelize as much as the workflow as possible, we analyze by chromosome
    call splitBedToChr as split_bed_to_chr {
        input:
        interval_bed = interval_to_bed.interval_bed
    }

    scatter (chr_bed in split_bed_to_chr.split_chr) {
        # Mutect
        if (tumor_only) {
            call mutectTumorOnly as mutectTumorTask {
              input:
              reference = reference,
              reference_fai = reference_fai,
              reference_dict = reference_dict,
              gnomad = normalized_gnomad_exclude,
              gnomad_tbi = normalized_gnomad_exclude_tbi,
              tumor_bam = bqsr.bqsr_bam,
              tumor_bam_bai = bqsr.bqsr_bam_bai,
              interval_list = chr_bed
            }
        }

        if (!tumor_only) {
            call mutectNormal as mutectNormalTask {
              input:
              reference = reference,
              reference_fai = reference_fai,
              reference_dict = reference_dict,
              gnomad = normalized_gnomad_exclude,
              gnomad_tbi = normalized_gnomad_exclude_tbi,
              tumor_bam = bqsr.bqsr_bam,
              tumor_bam_bai = bqsr.bqsr_bam_bai,
              normal_bam = normal_bam,
              normal_bam_bai = normal_bam_bai,
              interval_list = chr_bed
            }
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as mutectSanitizeVcf {
            input: vcf = select_first([mutectTumorTask.vcf, mutectNormalTask.vcf])
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as mutectNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = mutectSanitizeVcf.sanitized_vcf,
            vcf_tbi = mutectSanitizeVcf.sanitized_vcf_tbi
        }

        if (!pilot) {
        # Removes any germline variant that is reported in gnomAD
            call bcftoolsIsecComplement as mutect_isec_complement_gnomAD {
                input:
                vcf = mutectNormalize.normalized_vcf,
                vcf_tbi = mutectNormalize.normalized_vcf_tbi,
                exclude_vcf = normalized_gnomad_exclude,
                exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
                output_vcf_name = "mutect." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
                output_type = "z"
            }
        }
        # Removes any variants that appeared in two of our PoN Samples at 2% VAF or higher
        call pon2Percent as mutect_pon2 {
            input:
            vcf = select_first([mutect_isec_complement_gnomAD.complement_vcf, mutectNormalize.normalized_vcf]),
            vcf2PON = mutect_pon2_file,
            vcf2PON_tbi = mutect_pon2_file_tbi,
            caller = "mutect",
            sample_name = tumor_sample_name
        }

        # Vardict
        if (tumor_only) {
            call vardictTumorOnly as vardictTumorTask {
                input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = bqsr.bqsr_bam,
                tumor_bam_bai = bqsr.bqsr_bam_bai,
                interval_bed = chr_bed,
                min_var_freq = af_threshold,
                tumor_sample_name = tumor_sample_name
            }
        }

        if (!tumor_only) {
            call vardictNormal as vardictNormalTask {
                input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = bqsr.bqsr_bam,
                tumor_bam_bai = bqsr.bqsr_bam_bai,
                tumor_sample_name = tumor_sample_name,
                normal_bam = normal_bam,
                normal_bam = normal_bam_bai,
                interval_bed = chr_bed,
                min_var_freq = af_threshold
            }
        }

        # Performs the BCBIO Filtering: http://bcb.io/2016/04/04/vardict-filtering/
        call bcftoolsFilterBcbio as bcbio_filter {
            input:
            vcf = select_first([vardictTumorTask.vcf, vardictNormalTask.vcf_tumor]),
            vcf_tbi = select_first([vardictTumorTask.vcf_tbi, vardictNormalTask.vcf_tumor_tbi]),
            filter_string = bcbio_filter_string,
            filter_flag = "exclude",
            output_type = "z",
            output_vcf_prefix = "vardict.bcbiofilter"
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as vardictSanitizeVcf {
            input: vcf = bcbio_filter.filtered_vcf
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as vardictNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = vardictSanitizeVcf.sanitized_vcf,
            vcf_tbi = vardictSanitizeVcf.sanitized_vcf_tbi
        }

        if (!pilot) {
            # Removes any germline variant that is reported in gnomAD
            call bcftoolsIsecComplement as vardict_isec_complement_gnomAD {
                input:
                vcf = vardictNormalize.normalized_vcf,
                vcf_tbi = vardictNormalize.normalized_vcf_tbi,
                exclude_vcf = normalized_gnomad_exclude,
                exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
                output_vcf_name = "vardict." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
                output_type = "z"
            }
        }
        # Removes any variants that appeared in two of our PoN Samples at 2% VAF or higher
        call pon2Percent as vardict_pon2 {
            input:
            vcf = select_first([vardict_isec_complement_gnomAD.complement_vcf, vardictNormalize.normalized_vcf]),
            vcf2PON = vardict_pon2_file,
            vcf2PON_tbi = vardict_pon2_file_tbi,
            caller = "vardict",
            sample_name = tumor_sample_name
        }

        # Lofreq
        if (tumor_only) {
            call lofreqTumorOnly as lofreqTumorTask {
                input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = bqsr.bqsr_bam,
                tumor_bam_bai = bqsr.bqsr_bam_bai,
                interval_bed = chr_bed
            }
        }

        if (!tumor_only) {
            call lofreqNormal as lofreqNormalTask {
                input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = bqsr.bqsr_bam,
                tumor_bam_bai = bqsr.bqsr_bam_bai,
                normal_bam = normal_bam,
                normal_bam_bai = normal_bam_bai,
                interval_bed = chr_bed
            }
        }

        # Lofreq does not output FORMAT and SAMPLE columns, so we need to reformat the VCF with these columns
        # since many downstream tools require the VCF to be formatted in this way
        call lofreqReformat as reformat {
            input:
                vcf = select_first([lofreqTumorTask.vcf, lofreqNormalTask.vcf]),
                tumor_sample_name = tumor_sample_name
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as lofreqSanitizeVcf {
            input: vcf = reformat.reformat_vcf
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as lofreqNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = lofreqSanitizeVcf.sanitized_vcf,
            vcf_tbi = lofreqSanitizeVcf.sanitized_vcf_tbi
        }

        if (!pilot) {
            # Removes any germline variant that is reported in gnomAD
            call bcftoolsIsecComplement as lofreq_isec_complement_gnomAD {
                input:
                vcf = lofreqNormalize.normalized_vcf,
                vcf_tbi = lofreqNormalize.normalized_vcf_tbi,
                exclude_vcf = normalized_gnomad_exclude,
                exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
                output_vcf_name = "lofreq." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
                output_type = "z"
            }
        }

        # Removes any variants that appeared in two of our PoN Samples at 2% VAF or higher
        call pon2Percent as lofreq_pon2 {
            input:
            vcf = select_first([lofreq_isec_complement_gnomAD.complement_vcf, lofreqNormalize.normalized_vcf]),
            vcf2PON = lofreq_pon2_file,
            vcf2PON_tbi = lofreq_pon2_file_tbi,
            caller = "lofreq",
            sample_name = tumor_sample_name
        }

        # Pindel
        if (tumor_only) {
            call pindelTumorOnly as pindelTumorTask {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                tumor_bam = bqsr.bqsr_bam,
                tumor_bam_bai = bqsr.bqsr_bam_bai,
                region_file = chr_bed,
                insert_size = pindel_insert_size,
                tumor_sample_name = tumor_sample_name
            }
            call catOut as pindelTumorCat {
                input: pindel_outs = [pindelTumorTask.deletions, pindelTumorTask.insertions, pindelTumorTask.tandems, pindelTumorTask.long_insertions, pindelTumorTask.inversions]
            }
            call pindelToVcf {
                input:
                reference=reference,
                reference_fai=reference_fai,
                reference_dict=reference_dict,
                pindel_output_summary=pindelTumorCat.pindel_out,
                ref_name = ref_name,
                ref_date = ref_date,
                min_supporting_reads = pindel_min_supporting_reads,
                tumor_sample_name = tumor_sample_name
            }
        }

        if (!tumor_only) {
            call pindelNormal as PindelNormalTask {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                tumor_bam = bqsr.bqsr_bam,
                tumor_bam_bai = bqsr.bqsr_bam_bai,
                normal_bam = normal_bam,
                normal_bam_bai = normal_bam_bai,
                region_file = chr_bed,
                tumor_sample_name = tumor_sample_name,
                normal_sample_name = "NORMAL",
                insert_size = pindel_insert_size
            }
            call catOut as pindelNormalCat {
                input: pindel_outs = [PindelNormalTask.deletions, PindelNormalTask.insertions, PindelNormalTask.tandems, PindelNormalTask.long_insertions, PindelNormalTask.inversions]
            }
            call pindelSomaticFilter {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                pindel_output_summary = pindelNormalCat.pindel_out
            }
        }

        call removeEndTags {
            input: vcf = select_first([pindelToVcf.pindel_vcf, pindelSomaticFilter.pindel_vcf])
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as pindelSanitizeVcf {
            input: vcf = removeEndTags.processed_vcf
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as pindelNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = pindelSanitizeVcf.sanitized_vcf,
            vcf_tbi = pindelSanitizeVcf.sanitized_vcf_tbi
        }

        # In order to be efficient, we run all of the annotation and filtering ONCE. In order to do this, we need to merge
        # all of the callers together into one giant VCF that will be used as the main VCF for VEP, PoN, etc...
        scatter (caller_vcf in [mutect_pon2.annotated_vcf, vardict_pon2.annotated_vcf, lofreq_pon2.annotated_vcf]){
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
            bam = bqsr.bqsr_bam,
            vcf=mergeCallers.merged_vcf,
            sample_name=tumor_sample_name,
            min_var_freq=af_threshold,
            output_vcf_basename = "all_callers_full"
        }

        # Using the fake VCF of all the variants calls, perform a "pileup" on the PoN BAMs`
        scatter (pon_bam in pon_bams) {
            call mskGetBaseCounts {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                normal_bam = pon_bam,
                pon_final_name = "all_callers." + tumor_sample_name + ".pon.pileup",
                vcf = mergeCallers.merged_vcf
            }
        }
        call bcftoolsMerge as pileup_merge {
            input:
                vcfs = mskGetBaseCounts.pileup,
                vcf_tbis = mskGetBaseCounts.pileup_tbi
        }

        # Using the fake VCF of all the variants calls, perform VEP on all the variants
        call vep {
          input:
              vcf = mergeCallers.merged_vcf,
              cache_dir_zip = vep_cache_dir_zip,
              reference = reference,
              reference_fai = reference_fai,
              reference_dict = reference_dict,
              plugins = vep_plugins,
              spliceAI_files = vep_plugin_spliceAI_files,
              ensembl_assembly = vep_ensembl_assembly,
              ensembl_version = vep_ensembl_version,
              ensembl_species = vep_ensembl_species,
              synonyms_file = synonyms_file,
              custom_annotations = vep_custom_annotations,
              coding_only = annotate_coding_only,
              everything = everything,
              pick = vep_pick
        }

        # Using the "pileup", perform a Fisher's Exact Test with the Variants in each Caller
        call normalFisher as mutect_call_R_fisher {
            input:
            vcf = mutect_pon2.annotated_vcf,
            pon = pileup_merge.merged_vcf,
            pon_tbi = pileup_merge.merged_vcf_tbi,
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
            vcf = vardict_pon2.annotated_vcf,
            pon = pileup_merge.merged_vcf,
            pon_tbi = pileup_merge.merged_vcf_tbi,
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
            vcf = lofreq_pon2.annotated_vcf,
            pon = pileup_merge.merged_vcf,
            pon_tbi = pileup_merge.merged_vcf_tbi,
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
    }

    call mergeVcf as merge_mutect_full {
        input:
            vcfs = mutectNormalize.normalized_vcf,
            vcf_tbis = mutectNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "mutect_full." + tumor_sample_name
    }
    call mergeVcf as merge_mutect_pon {
        input:
            vcfs = mutect_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = mutect_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mergeVcf as merge_mutect_final {
        input:
            vcfs = mutect_annotate_vcf.final_annotated_vcf,
            vcf_tbis = mutect_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name + ".final.annotated"
    }

    call mergeVcf as merge_vardict_full {
        input:
            vcfs = vardictNormalize.normalized_vcf,
            vcf_tbis = vardictNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "vardict_full." + tumor_sample_name
    }
    call mergeVcf as merge_vardict_pon {
        input:
            vcfs = vardict_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = vardict_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mergeVcf as merge_vardict_final {
        input:
            vcfs = vardict_annotate_vcf.final_annotated_vcf,
            vcf_tbis = vardict_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name + ".final.annotated"
    }

    call mergeVcf as merge_lofreq_full {
        input:
            vcfs = lofreqNormalize.normalized_vcf,
            vcf_tbis = lofreqNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "lofreq_full." + tumor_sample_name
    }
    call mergeVcf as merge_lofreq_pon {
        input:
            vcfs = lofreq_annotate_vcf.pon_annotated_vcf,
            vcf_tbis = lofreq_annotate_vcf.pon_annotated_vcf_tbi,
            merged_vcf_basename = "lofreq." + tumor_sample_name + ".pileup.fisherPON"
    }
    call mergeVcf as merge_lofreq_final {
        input:
            vcfs = lofreq_annotate_vcf.final_annotated_vcf,
            vcf_tbis = lofreq_annotate_vcf.final_annotated_vcf_tbi,
            merged_vcf_basename = "lofreq." + tumor_sample_name + ".final.annotated"
    }

    call mergeVcf as merge_pindel_full {
        input:
            vcfs = pindelNormalize.normalized_vcf,
            vcf_tbis = pindelNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "pindel_full." + tumor_sample_name
    }

    call mergeVcf as merge_pon {
        input:
            vcfs = pileup_merge.merged_vcf,
            vcf_tbis = pileup_merge.merged_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".pon.total.counts"
    }

    call mergeVcf as merge_fp_filter {
        input:
            vcfs = fpFilter.filtered_vcf,
            vcf_tbis = fpFilter.filtered_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".fpfilter"
    }

    call mergeVcf as merge_vep {
        input:
            vcfs = vep.annotated_vcf,
            vcf_tbis = vep.annotated_vcf_tbi,
            merged_vcf_basename = tumor_sample_name + ".vep"
    }

    call annotatePD {
        input:
            mutect_vcf = merge_mutect_final.merged_vcf,
            lofreq_vcf = merge_lofreq_final.merged_vcf,
            vardict_vcf = merge_vardict_final.merged_vcf,
            bolton_bick_vars = bolton_bick_vars,
            mut2_bick = mut2_bick,
            mut2_kelly = mut2_kelly,
            matches2 = matches2,
            TSG_file = TSG_file,
            oncoKB_curated = oncoKB_curated,
            pd_annotation_file = pd_annotation_file,
            truncating = truncating,
            cosmic_dir_zip = cosmic_dir_zip,
            pon_pvalue = pon_pvalue
    }

    if (platform == 'ArcherDX') {
        call xgb_model {
            input:
            mutect_tsv = annotatePD.mutect_vcf_annotate_pd,
            lofreq_tsv = annotatePD.lofreq_vcf_annotate_pd,
            vardict_tsv = annotatePD.vardict_vcf_annotate_pd,
            pindel_full_vcf = merge_pindel_full.merged_vcf,
            pon = merge_pon.merged_vcf,
            pon_pvalue = pon_pvalue,
            model = true,
            tumor_sample_name = tumor_sample_name
        }
    }

    if (platform != 'ArcherDX') {
        call xgb_model as no_xgb_model {
            input:
            mutect_tsv = annotatePD.mutect_vcf_annotate_pd,
            lofreq_tsv = annotatePD.lofreq_vcf_annotate_pd,
            vardict_tsv = annotatePD.vardict_vcf_annotate_pd,
            pindel_full_vcf = merge_pindel_full.merged_vcf,
            pon = merge_pon.merged_vcf,
            pon_pvalue = pon_pvalue,
            model = false,
            tumor_sample_name = tumor_sample_name
        }
    }

    output {
        # Alignments
        File? aligned_bam = select_first([filterClipAndCollectMetrics.clipped_bam, clipAndCollectMetrics.clipped_bam])
        File bqsr_bam = bqsr.bqsr_bam

        # Tumor QC
        File tumor_insert_size_metrics = collectAllMetrics.insert_size_metrics
        File tumor_alignment_summary_metrics = collectAllMetrics.alignment_summary_metrics
        File tumor_hs_metrics = collectRoiHsMetrics.hs_metrics
        Array[File] tumor_per_target_coverage_metrics = select_all(collectPerTargetHsMetrics.per_target_coverage_metrics)
        Array[File] tumor_per_target_hs_metrics = collectPerTargetHsMetrics.hs_metrics
        Array[File] tumor_per_base_coverage_metrics = select_all(collectPerBaseHsMetrics.per_base_coverage_metrics)
        Array[File] tumor_per_base_hs_metrics = collectPerBaseHsMetrics.hs_metrics
        Array[File] tumor_summary_hs_metrics = collectSummaryHsMetrics.hs_metrics
        File tumor_flagstats = samtoolsFlagstat.flagstats
        File tumor_verify_bam_id_metrics = verifyBamId.verify_bam_id_metrics
        File tumor_verify_bam_id_depth = verifyBamId.verify_bam_id_depth
        File fastqc_html = fastQC.fastqc_html
        File fastqc = fastQC.fastqc
        File somalier_out = somalier.somalier_out


        # Mutect
        File mutect_full =  merge_mutect_full.merged_vcf                                # Raw Mutect Ouput
        File mutect_pon_annotated_vcf = merge_mutect_pon.merged_vcf                     # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        File mutect_vep_annotated_vcf = merge_mutect_final.merged_vcf                   # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Lofreq
        File lofreq_full = merge_lofreq_full.merged_vcf                                 # Raw Lofreq Ouput
        File lofreq_pon_annotated_vcf = merge_lofreq_pon.merged_vcf                     # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        File lofreq_vep_annotated_vcf = merge_lofreq_final.merged_vcf                   # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Vardict
        File vardict_full = merge_vardict_full.merged_vcf                               # Raw Vardict Ouput
        File vardict_pon_annotated_vcf = merge_vardict_pon.merged_vcf                   # gnomAD Filtered + PoN Filtered + PoN2 Annotated
        File vardict_vep_annotated_vcf = merge_vardict_final.merged_vcf                 # gnomAD Filtered + PoN Filtered + PoN2 Annotated w/ VEP Annotation

        # Pindel
        File pindel_full = merge_pindel_full.merged_vcf                                 # Raw Pindel Ouput

        File pon_total_counts = merge_pon.merged_vcf                                    # PoN Pileup Results
        File fpfilter_results = merge_fp_filter.merged_vcf
        File vep_results = merge_vep.merged_vcf

        #File gnomAD_exclude = get_gnomad_exclude.normalized_gnomad_exclude

        # R Things
        File mutect_annotate_pd = annotatePD.mutect_vcf_annotate_pd
        File lofreq_annotate_pd = annotatePD.lofreq_vcf_annotate_pd
        File vardict_annotate_pd = annotatePD.vardict_vcf_annotate_pd

        # Model
        File model_output = select_first([xgb_model.model_output, no_xgb_model.model_output])
        File model_raw_output = select_first([xgb_model.model_raw_output, no_xgb_model.model_raw_output])
        File mutect_complex = select_first([xgb_model.mutect_complex, no_xgb_model.mutect_complex])
        File pindel_complex = select_first([xgb_model.pindel_complex, no_xgb_model.pindel_complex])
        File lofreq_complex = select_first([xgb_model.lofreq_complex, no_xgb_model.lofreq_complex])
        File caller_filters = select_first([xgb_model.caller_filters, no_xgb_model.caller_filters])
    }
}

task filterArcherUmiLength {
    input {
        File fastq1
        File fastq2
        Int umi_length
        Float? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([fastq1, fastq2], "GB")
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0]) #2 GB or 1 GB
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])


    runtime {
        docker: "ubuntu:bionic"
        memory: cores*memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        zcat ~{fastq1} | awk -v regex="AACCGCCAGGAGT" -v umi_length="~{umi_length}" 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; split(seq,a,regex); if (length(a[1]) == umi_length) {print header, seq, qheader, qseq}}' > R1_filtered.fastq
        gzip R1_filtered.fastq
        cp ~{fastq2} R2_filtered.fastq.gz
    >>>

    output {
        File fastq1_filtered = "R1_filtered.fastq.gz"
        File fastq2_filtered = "R2_filtered.fastq.gz"
    }
}

task bamToFastq {
    input {
        File unaligned_bam
        Float? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(unaligned_bam, "GB")
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0]) #2 GB or 1 GB
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "mgibio/rnaseq:1.0.0"
        memory: cores*memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I=~{unaligned_bam} F=read1.fastq F2=read2.fastq
    >>>

    output {
        File fastq_one = "read1.fastq"
        File fastq_two = "read2.fastq"
    }
}

task bbmapRepair {
    input {
        File fastq1
        File fastq2
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([fastq1, fastq2], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Float memory = select_first([mem_limit_override, if 4.0 * data_size > 12.0 then ceil(4.0 * data_size) else 12.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "quay.io/biocontainers/bbmap:38.92--he522d1c_0"
        memory: cores*memory + "GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        repair.sh -Xmx10g repair=t overwrite=true interleaved=false outs=singletons.fq out1=R1.fixed.fastq.gz out2=R2.fixed.fastq.gz in1=~{fastq1} in2=~{fastq2}
    >>>

    output {
        Array[File] fastqs = ["R1.fixed.fastq.gz", "R2.fixed.fastq.gz"]
        File fastq1_repair = "R1.fixed.fastq.gz"
        File fastq2_repair = "R2.fixed.fastq.gz"
    }
}

task fastqToBam {
    input {
        File fastq1
        File fastq2
        String sample_name
        String library_name
        String platform_unit
        String platform
        Float? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([fastq1, fastq2], "GB")
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory:  cores*memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar FastqToSam FASTQ=~{fastq1} FASTQ2=~{fastq2} SAMPLE_NAME=~{sample_name} LIBRARY_NAME=~{library_name} PLATFORM_UNIT=~{platform_unit} PLATFORM=~{platform} OUTPUT=unaligned.bam
    >>>

    output {
        File bam = "unaligned.bam"
    }
}

task extractUmis {
    input {
        File bam
        Array[String] read_structure
        Boolean? umi_paired = true
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Int preemptible = 1
    Int maxRetries = 0
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: cores*memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        if [ "~{umi_paired}" == true ]; then
            /usr/local/bin/fgbio ExtractUmisFromBam --molecular-index-tags ZA ZB --single-tag RX --input ~{bam} --read-structure ~{sep=" " read_structure} --output umi_extracted.bam
        else
            /usr/local/bin/fgbio ExtractUmisFromBam --molecular-index-tags ZA --single-tag RX --input ~{bam} --read-structure ~{sep=" " read_structure} --output umi_extracted.bam
        fi
    >>>

    output {
        File umi_extracted_bam = "umi_extracted.bam"
    }
}

task markIlluminaAdapters {
    input {
        File bam
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Int preemptible = 1
    Int maxRetries = 0
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar MarkIlluminaAdapters INPUT=~{bam} OUTPUT=marked.bam METRICS=adapter_metrics.txt
    >>>

    output {
        File marked_bam = "marked.bam"
        File metrics = "adapter_metrics.txt"
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
        Float? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")
    Int space_needed_gb = ceil(10 + 10 * data_size + reference_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 48.0 then floor(memory / 24)*8 else 8])

    runtime {
      docker: "mgibio/dna-alignment:1.0.0"
      memory: cores * memory + "GB"
      cpu: cores
      # 1 + just for a buffer
      # data_size*10 because bam uncompresses and streams to /dev/stdout and /dev/stdin, could have a couple flying at once
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /bin/bash /usr/bin/umi_alignment.sh ~{bam} ~{reference} ~{cores}
    >>>

    output {
        File aligned_bam = "aligned.bam"
        File aligned_bam_bai = "aligned.bai"
    }
}

task groupReadsAndConsensus {
    input {
        File bam
        Boolean umi_paired = true
        Float? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(4.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        PAIRED=~{umi_paired}
        BAM=~{bam}

        if [ "$PAIRED" == true ]; then
            /usr/local/bin/fgbio GroupReadsByUmi --strategy paired --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        else
            /usr/local/bin/fgbio GroupReadsByUmi --strategy adjacency --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        fi
        /usr/local/bin/fgbio CallMolecularConsensusReads --input umi_grouped.bam --error-rate-pre-umi 45 --error-rate-post-umi 30 --min-input-base-quality 30 --min-reads 1 --output consensus_unaligned.bam
    >>>

    output {
        File consensus_bam = "consensus_unaligned.bam"
    }
}


task realign {
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
        Float? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")
    Int space_needed_gb = ceil(10 + 10*data_size + reference_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 48.0 then floor(memory / 24)*8 else 8])

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /bin/bash /usr/bin/umi_realignment.sh ~{bam} ~{reference} ~{cores}
    >>>

    output {
        File consensus_aligned_bam = "realigned.bam"
        File consensus_aligned_bam_bai = "realigned.bai"
    }
}

task filterClipAndCollectMetrics {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        Array[Int] min_reads = [1]
        Float? max_read_error_rate = 0.05
        Float? max_base_error_rate = 0.1
        Int min_base_quality = 1
        Float max_no_call_fraction = 0.5
        File? target_intervals
        String description
        Boolean umi_paired = true
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = ceil(10 + 2 * data_size + reference_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        /usr/local/bin/fgbio FilterConsensusReads --input ~{bam} --output consensus_filtered.bam --ref ~{reference} --min-reads ~{sep=" " min_reads} --max-read-error-rate ~{max_read_error_rate} --max-base-error-rate ~{max_base_error_rate} --min-base-quality ~{min_base_quality} --max-no-call-fraction ~{max_no_call_fraction}
        /usr/local/bin/fgbio ClipBam --input consensus_filtered.bam --ref ~{reference} --clipping-mode Hard --clip-overlapping-reads true --output clipped.bam

        PAIRED=~{umi_paired}
        DESCRIPTION=~{description}
        INTERVALS=~{target_intervals}

        if [ "$PAIRED" == true ]; then
            if [[ -z "$DESCRIPTION" ]]; then
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            else
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            fi
        else
            echo "Sample not UMI Paired" > duplex_seq.metrics.txt
        fi
    >>>

    output {
        File clipped_bam = "clipped.bam"
        File clipped_bam_bai = "clipped.bai"
        Array[File] duplex_seq_metrics = glob("duplex_seq.metrics.*")
    }
}

task clipAndCollectMetrics {
    input {
        File bam
        File reference
        File reference_fai
        File reference_dict
        File? target_intervals
        String description
        Boolean umi_paired = true
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = ceil(10 + 2 * data_size + reference_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eo pipefail

        /usr/local/bin/fgbio ClipBam --input ~{bam} --ref ~{reference} --clipping-mode Hard --clip-overlapping-reads true --output clipped.bam

        PAIRED=~{umi_paired}
        DESCRIPTION=~{description}
        INTERVALS=~{target_intervals}

        if [ "$PAIRED" == true ]; then
            if [[ -z "$DESCRIPTION" ]]; then
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            else
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            fi
        else
            echo "Sample not UMI Paired" > duplex_seq.metrics.txt
        fi
    >>>

    output {
        File clipped_bam = "clipped.bam"
        File clipped_bam_bai = "clipped.bai"
        Array[File] duplex_seq_metrics = glob("duplex_seq.metrics.*")
    }
}

task bqsrApply {
    input {
        File reference
        File reference_fai
        File reference_dict
        File bam
        File bam_bai
        String output_name = "final"
        Array[File] known_sites
        Array[File] known_sites_tbi  # secondaryFiles...
        Array[String] intervals = ["chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([bam, bam_bai], "GB")
    Float reference_size = size(known_sites, "GB") + size(known_sites_tbi, "GB") + size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = ceil(10 + 2 * data_size + reference_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 18.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.1.8.1"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /gatk/gatk --java-options -Xmx16g BaseRecalibrator -O bqsr.table ~{sep=" " prefix("-L ", intervals)} -R ~{reference} -I ~{bam} ~{sep=" " prefix("--known-sites ", known_sites)}
        /gatk/gatk --java-options -Xmx16g ApplyBQSR -O ~{output_name}.bam ~{sep=" " prefix("--static-quantized-quals ", [10, 20, 30])} -R ~{reference} -I ~{bam} -bqsr bqsr.table
    >>>

    output {
        File bqsr_bam = "~{output_name}.bam"
        File bqsr_bam_bai = "~{output_name}.bai"
    }
}

task Metrics {
    input {
        File bam
        File bam_bai
        File reference
        File reference_fai
        File reference_dict
        String metric_accumulation_level
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([bam, bam_bai], "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = ceil(10 + data_size + reference_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        memory: cores * memory + "GB"
        docker: "broadinstitute/picard:2.23.6"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bamroot = basename(bam, ".bam")
    String summary_metrics = "~{bamroot}.AlignmentSummaryMetrics.txt"
    String size_metrics = "~{bamroot}.InsertSizeMetrics.txt"
    String size_histogram = "~{bamroot}.InsertSizeHistogram.pdf"

    command <<<
        /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics INPUT=~{bam} OUTPUT=~{summary_metrics} REFERENCE_SEQUENCE=~{reference} METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level}
        /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar CollectInsertSizeMetrics O=~{size_metrics} H=~{size_histogram} I=~{bam} REFERENCE_SEQUENCE=~{reference} METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level}
    >>>

    output {
    File alignment_summary_metrics = summary_metrics
    File insert_size_histogram = size_histogram
    File insert_size_metrics = size_metrics
    }
}

task collectHsMetrics {
    input {
        File bam
        File bam_bai
        File reference
        File reference_fai
        File reference_dict
        String metric_accumulation_level

        File bait_intervals
        File target_intervals
        Boolean per_target_coverage = false
        Boolean per_base_coverage = false
        Int? minimum_base_quality
        Int? minimum_mapping_quality

        String output_prefix = "out"

        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([bam, bam_bai], "GB")
    Float reference_size = size([reference, reference_fai, reference_dict, bait_intervals, target_intervals], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = ceil(10 + data_size + reference_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        memory: cores * memory + "GB"
        docker: "broadinstitute/picard:2.23.6"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bamroot = basename(bam, ".bam")
    String hs_txt = "~{bamroot}.~{output_prefix}-HsMetrics.txt"
    String per_target_txt = "~{bamroot}.~{output_prefix}-PerTargetCoverage.txt"
    String per_base_txt = "~{bamroot}.~{output_prefix}-PerBaseCoverage.txt"

    command <<<
        /usr/bin/java -Xmx6g -jar /usr/picard/picard.jar CollectHsMetrics \
        O=~{hs_txt} \
        I=~{bam} \
        R=~{reference} \
        TARGET_INTERVALS=~{target_intervals} \
        METRIC_ACCUMULATION_LEVEL=~{metric_accumulation_level} \
        BAIT_INTERVALS=~{bait_intervals} \
        ~{if per_target_coverage then "PER_TARGET_COVERAGE=~{per_target_txt}" else ""} \
        ~{if per_base_coverage then "PER_BASE_COVERAGE=~{per_base_txt}" else ""} \
        ~{if defined(minimum_mapping_quality) then "MINIMUM_MAPPING_QUALITY=~{minimum_mapping_quality}" else ""} \
        ~{if defined(minimum_base_quality) then "MINIMUM_BASE_QUALITY=~{minimum_base_quality}" else ""}
    >>>

    output {
        File hs_metrics = hs_txt
        File? per_target_coverage_metrics = per_target_txt
        File? per_base_coverage_metrics = per_base_txt
    }
}

task samtoolsFlagstat {
    input {
        File bam
        File bam_bai
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([bam, bam_bai], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = ceil(5 + data_size)
    Float memory = select_first([mem_limit_override, if 2.0 * data_size > 6.0 then ceil(2.0 * data_size) else 6.0])
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 32) else 1])

    runtime {
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String outfile = basename(bam) + ".flagstat"

    command <<<
        /usr/local/bin/samtools flagstat ~{bam} > ~{outfile}
    >>>

    output {
        File flagstats = outfile
    }
}

task selectVariants {
    input {
        File reference
        File reference_fai
        File reference_dict

        File vcf
        File vcf_tbi

        File? interval_list
        Boolean exclude_filtered = false
        String output_vcf_basename = "select_variants"
        Array[String]? samples_to_include  # include genotypes from this sample

        # ENUM: one of ["INDEL", "SNP", "MIXED", "MNP", "SYMBOLIC", "NO_VARIATION"]
        String? select_type
    }

    Int space_needed_gb = 5 + round(size([vcf, vcf_tbi], "GB")*3 + size([reference, reference_fai, reference_dict, interval_list], "GB"))
    Int preemptible = 1
    Int maxRetries = 0
    Int cores = 1

    runtime {
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "6GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String outfile = "~{output_vcf_basename}.vcf.gz"
    Array[String] samples = if defined(samples_to_include) then prefix("--sample-name ", select_first([samples_to_include])) else []

    command <<<
        /gatk/gatk --java-options -Xmx4g SelectVariants -O ~{outfile} \
        -R ~{reference} \
        --variant ~{vcf} \
        ~{if defined(interval_list) then "-L ~{interval_list}" else ""} \
        ~{if exclude_filtered then "--exclude-filtered ~{select_first([exclude_filtered])}" else ""} \
        ~{sep=" " samples} \
        ~{if defined(select_type) then "-select-type ~{select_type}" else ""}
    >>>

    output {
        File filtered_vcf = outfile
        File filtered_vcf_tbi = "~{outfile}.tbi"
    }
}

task verifyBamId {
    input {
        File vcf
        File bam
        File bam_bai
    }

    Int space_needed_gb = 5 + round(size([bam, bam_bai, vcf], "GB"))
    Int preemptible = 1
    Int maxRetries = 0
    Int cores = 1

    runtime {
        docker: "mgibio/verify_bam_id-cwl:1.1.3"
        memory: "4GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bamroot = basename(bam, ".bam")
    String outroot = "~{bamroot}.VerifyBamId"

    command <<<
        /usr/local/bin/verifyBamID --out ~{outroot} --vcf ~{vcf} --bam ~{bam} --bai ~{bam_bai}
    >>>

    output {
        File verify_bam_id_metrics = "~{outroot}.selfSM"
        File verify_bam_id_depth = "~{outroot}.depthSM"
    }
}

task fastQC {
    input {
        File bam
        File bam_bai
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([bam, bam_bai], "GB")
    Int space_needed_gb = 10 + round(data_size)

    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
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

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(interval_list, "GB")
    Int space_needed_gb = 10 + round(data_size)

    runtime {
        docker: "ubuntu:bionic"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
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

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([interval_bed, chrom_sizes, af_only_snp_only_vcf, reference], "GB")
    Int space_needed_gb = 10 + round(data_size)

    runtime {
        docker: "kboltonlab/bst:1.0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        bedtools slop -i ~{interval_bed} -g ~{chrom_sizes} -b 100 > somalier.bed
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
        File bam_bai
        String sample_name = "tumor"
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([somalier_vcf, reference, bam, bam_bai], "GB")
    Int space_needed_gb = 10 + round(data_size)

    runtime {
        docker: "brentp/somalier:latest"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
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

    Int cores = 1
    Int space_needed_gb = 10 + round(size(interval_bed, "GB")*2)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: "6GB"
        docker: "kboltonlab/bst:latest"
        disks: "local-disk ~{space_needed_gb} SSD"
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

task mutectNormal {
    input {
        File reference
        File reference_fai
        File reference_dict
        File? pon
        File? pon_tbi
        File? gnomad
        File? gnomad_tbi

        File tumor_bam
        File tumor_bam_bai

        File? normal_bam
        File? normal_bam_bai

        File interval_list
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_list, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "32GB"
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = "mutect.filtered.vcf.gz"

    command <<<
        set -o pipefail
        set -o errexit

        NORMAL=`samtools view -H ~{normal_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
        TUMOR=`samtools view -H ~{tumor_bam} | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`


        /gatk/gatk Mutect2 --java-options "-Xmx20g" -O mutect.vcf.gz -R ~{reference} -L ~{interval_list} \
            -I ~{tumor_bam} --read-index ~{tumor_bam_bai} -tumor "$TUMOR" \
            -I ~{normal_bam} --read-index ~{normal_bam_bai} -normal "$NORMAL" \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            --f1r2-tar-gz mutect.f1r2.tar.gz --max-reads-per-alignment-start 0

        /gatk/gatk LearnReadOrientationModel -I mutect.f1r2.tar.gz -O mutect.read-orientation-model.tar.gz
        /gatk/gatk FilterMutectCalls -R ~{reference} -V mutect.vcf.gz --ob-priors mutect.read-orientation-model.tar.gz -O ~{output_vcf} #Running FilterMutectCalls on the output vcf.
    >>>

    output {
    File vcf = output_vcf
    File vcf_tbi = output_vcf + ".tbi"
    }
}

task mutectTumorOnly {
    input {
        File reference
        File reference_fai
        File reference_dict
        File? pon
        File? pon_tbi
        File? gnomad
        File? gnomad_tbi
        File tumor_bam
        File tumor_bam_bai

        File interval_list
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_list, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "32GB"
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = "mutect.filtered.vcf.gz"

    command <<<
        set -o pipefail
        set -o errexit

        /gatk/gatk Mutect2 --java-options "-Xmx20g" \
            -O mutect.vcf.gz \
            -R ~{reference} \
            -L ~{interval_list} \
            -I ~{tumor_bam} \
            ~{"--germline-resource " + gnomad} \
            ~{"-pon " + pon} \
            --read-index ~{tumor_bam_bai} \
            --f1r2-tar-gz mutect.f1r2.tar.gz \
            --max-reads-per-alignment-start 0

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

task vcfSanitize {
    input {
        File vcf
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: "6GB"
        cpu: cores
        docker: "mgibio/samtools-cwl:1.0.0"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    # outbase should match in script but I don't want to risk changing it yet
    String outbase = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        set -eou pipefail

        # 1) removes lines containing non ACTGN bases, as they conflict with the VCF spec and cause GATK to choke
        # 2) removes mutect-specific format tags containing underscores, which are likewise illegal in the vcf spec
        base=`basename ~{vcf}`
        outbase=`echo $base | perl -pe 's/.vcf(.gz)?$//g'`
        echo "~{vcf}   $base    $outbase"
        if [[ "~{vcf}" =~ ".gz" ]];then
        #gzipped input
        gunzip -c "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
        else
        #non-gzipped input
        cat "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
        fi
        /opt/htslib/bin/bgzip $outbase.sanitized.vcf
        /usr/bin/tabix -p vcf $outbase.sanitized.vcf.gz
    >>>

    output {
        File sanitized_vcf = outbase + ".sanitized.vcf.gz"
        File sanitized_vcf_tbi = outbase + ".sanitized.vcf.gz.tbi"
    }
}

task bcftoolsNorm {
    input {
        File reference
        File reference_fai

        File vcf
        File vcf_tbi
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB") * 2 + size([reference, reference_fai], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: "6GB"
        docker: "kboltonlab/bst"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/bcftools norm --check-ref w --multiallelics -any --output-type z --output bcftools_norm.vcf.gz ~{vcf} -f ~{reference}
        /usr/local/bin/tabix bcftools_norm.vcf.gz
    >>>

    output {
        File normalized_vcf = "bcftools_norm.vcf.gz"
        File normalized_vcf_tbi = "bcftools_norm.vcf.gz.tbi"
    }
}

task bcftoolsIsecComplement {
    input {
        File vcf
        File vcf_tbi
        File exclude_vcf
        File exclude_vcf_tbi
        String output_type = "z"
        String? output_vcf_name = "bcftools_isec.vcf"
    }

    Int space_needed_gb = 10 + 2*round(size([vcf, vcf_tbi, exclude_vcf, exclude_vcf_tbi], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "kboltonlab/bst:latest"
        memory: "6GB"
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/bcftools isec -C -w1 ~{vcf} ~{exclude_vcf} --output-type ~{output_type} --output ~{output_vcf_name}.gz && /usr/local/bin/tabix ~{output_vcf_name}.gz
    >>>

    output {
        File complement_vcf = "~{output_vcf_name}.gz"
        File complement_vcf_tbi = "~{output_vcf_name}.gz.tbi"
    }
}

task pon2Percent {
    input {
        File vcf
        File vcf2PON
        File vcf2PON_tbi
        String caller = "caller"
        String sample_name = "tumor"
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf2PON, vcf2PON_tbi],"GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        export name=~{caller}.~{sample_name}.pon2.annotated.vcf.gz

        printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
        printf "##INFO=<ID=PON_NAT2_percent,Number=1,Type=Integer,Description=\"Number of samples with variant at >=2 percent\">\n" >> pon2.header;
        printf "##INFO=<ID=PON_MAX_VAF,Number=1,Type=Float,Description=\"The maximum VAF found in the PoN Samples\">\n" >> pon2.header;
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\t%INFO/NS\t%INFO/max_VAF\n" ~{vcf2PON} > normal2.txt
        bgzip -f normal2.txt
        tabix -f -s1 -b2 -e2 normal2.txt.gz
        bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent,PON_NAT2_percent,PON_MAX_VAF ~{vcf} -Oz -o $name
        tabix $name
    >>>

    output {
        File annotated_vcf = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz"
        File annotated_vcf_tbi = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz.tbi"
    }
}

task vardictTumorOnly {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File interval_bed
        Float? min_var_freq = 0.005
    }

    Int cores = 16
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/vardictjava:1.0"
        memory: "96GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_name = "vardict"

    command <<<
        set -o pipefail
        set -o errexit

        /opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G ~{reference} \
            -X 1 \
            -f ~{min_var_freq} \
            -N ~{tumor_sample_name} \
            -b ~{tumor_bam} \
            -c 1 -S 2 -E 3 -g 4 ~{interval_bed} \
            -th ~{cores} | \
        /opt/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
        /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
            -N "~{tumor_sample_name}" \
            -E \
            -f ~{min_var_freq} > ~{output_name}.vcf

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

task vardictNormal {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File? normal_bam
        File? normal_bam_bai
        String? normal_sample_name = "NORMAL"
        File interval_bed
        Float? min_var_freq = 0.005
    }

    Int cores = 16
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/vardictjava:1.0"
        memory: "96GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_name = "vardict"
    command <<<
        set -o pipefail
        set -o errexit

        /opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G ~{reference} \
            -X 1 \
            -f ~{min_var_freq} \
            -N ~{tumor_sample_name} \
            -b "~{tumor_bam}|~{normal_bam}" \
            -c 1 -S 2 -E 3 -g 4 ~{interval_bed} \
            -th ~{cores} | \
        /opt/VarDictJava/build/install/VarDict/bin/testsomatic.R | \
        /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl \
            -N "~{tumor_sample_name}|~{normal_sample_name}" \
            -f ~{min_var_freq} > ~{output_name}.vcf

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz

        # Extract tumor variants before bcbio filter or else both normal and tumor variants will be canidates for the filter
        /usr/bin/bcftools view -h \
            -s ~{tumor_sample_name} \
            --threads ~{cores} ~{output_name}.vcf.gz \
            -Oz -o ~{output_name}.tumor.vcf.gz
        /usr/bin/tabix ~{output_name}.tumor.vcf.gz

        # Extract normal variants in case we need to check the variants found in the normal sample
        /usr/bin/bcftools view -h \
            -s ~{normal_sample_name} \
            --threads ~{cores} ~{output_name}.vcf.gz \
            -Oz -o ~{output_name}.normal.vcf.gz
        /usr/bin/tabix ~{output_name}.normal.vcf.gz
    >>>

    output {
        File vcf_tumor = "~{output_name}.tumor.vcf.gz"
        File vcf_tumor_tbi = "~{output_name}.tumor.vcf.gz.tbi"
        File vcf_normal = "~{output_name}.normal.vcf.gz"
        File vcf_normal_tbi = "~{output_name}.normal.vcf.gz.tbi"
    }
}

task bcftoolsFilterBcbio {
    input {
        File vcf
        File vcf_tbi
        String filter_flag = "exclude"
        String filter_string
        String? output_vcf_prefix = "bcftools_filter"
        String output_type = "z"
    }

    Int space_needed_gb = 10 + 2*round(size(vcf, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String ff = if filter_flag == "include" then "-i" else "-e"
    command <<<
        /usr/local/bin/bcftools filter ~{ff} "~{filter_string}" ~{vcf} --output-type ~{output_type} --output ~{output_vcf_prefix}.vcf.gz -s "BCBIO" -m+

        /usr/local/bin/tabix ~{output_vcf_prefix}.vcf.gz

    >>>

    output {
        File filtered_vcf = "~{output_vcf_prefix}.vcf.gz"
        File filtered_vcf_tbi = "~{output_vcf_prefix}.vcf.gz.tbi"
    }
}

task lofreqTumorOnly {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        File interval_bed
        String? output_name = "lofreq.vcf"
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/lofreq:latest"
        memory: "24GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o output.indel.bam ~{tumor_bam}
        samtools index output.indel.bam
        /opt/lofreq/bin/lofreq call-parallel --pp-threads ~{cores} -A -B -f ~{reference} --call-indels --bed ~{interval_bed} -o ~{output_name} output.indel.bam --force-overwrite
        bgzip ~{output_name} && tabix ~{output_name}.gz
    >>>

    output {
        File vcf = "~{output_name}.gz"
        File vcf_tbi = "~{output_name}.gz.tbi"
    }
}

# TODO: Test Lofreq Normal Mode
task lofreqNormal {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File interval_bed
        String? output_name = "lofreq"
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/lofreq:latest"
        memory: "24GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o output.indel.bam ~{tumor_bam}
        /opt/lofreq/bin/lofreq somatic --call-indels -n ~{normal_bam} -t output.indel.bam -f ~{reference} -l ~{interval_bed} -o lofreq_ --threads ~{cores}
        tabix lofreq_somatic_final.snvs.vcf.gz
        tabix lofreq_somatic_final.indels.vcf.gz
        bcftools concat -a lofreq_somatic_final.snvs.vcf.gz lofreq_somatic_final.indels.vcf.gz > ~{output_name}.vcf
        bgzip ~{output_name}.vcf && tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

task lofreqReformat {
    input {
        File vcf
        String tumor_sample_name
    }

    Int space_needed_gb = 10 + 2*round(size(vcf, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        zcat ~{vcf} | grep "##" > lofreq.reformat.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">"  >> lofreq.reformat.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> lofreq.reformat.vcf;
        zcat ~{vcf} | grep -v '#' | awk '{ n=split($8, semi, /;/); sample=""; format=""; for(i in semi){ split(semi[i], equ, /=/); if(i<=3){ if(i+1==4) sample=sample equ[2]; else sample=sample equ[2] ":"; if(i+1==4) format=format equ[1]; else format=format equ[1] ":";}}{print $0, "GT:"format, "0/1:"sample}}' OFS='\t' >> lofreq.reformat.vcf;
        bgzip lofreq.reformat.vcf && tabix lofreq.reformat.vcf.gz
    >>>

    output {
        File reformat_vcf = "lofreq.reformat.vcf.gz"
        File reformat_vcf_tbi = "lofreq.reformat.vcf.gz.tbi"
    }
}

task pindelNormal {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File? normal_bam
        File? normal_bam_bai
        File region_file
        String tumor_sample_name
        String? normal_sample_name
        String? chromosome
        Int insert_size = 400
    }

    Int cores = 4
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai, region_file], "GB"))
    Int preemptible = 1
    Int maxRetries = 5

    runtime {
        bootDiskSizeGb: 100
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        docker: "mgibio/cle:v1.4.2"
        memory: "24GB"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        echo -e "~{normal_bam}\t~{insert_size}\t~{normal_sample_name}" > pindel.config
        echo -e "~{tumor_bam}\t~{insert_size}\t~{tumor_sample_name}" >> pindel.config

        /usr/bin/pindel -i pindel.config -w 30 -T ~{cores} -o all -f ~{reference} \
            ~{if defined(chromosome) then "-c ~{chromosome}" else ""} \
            ~{if defined(region_file) then "-j ~{region_file}" else ""}
    >>>

    output {
        File deletions = "all_D"
        File insertions = "all_SI"
        File tandems = "all_TD"
        File long_insertions = "all_LI"
        File inversions = "all_INV"
    }
}

task pindelTumorOnly {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File region_file
        String tumor_sample_name
        String? chromosome
        Int insert_size = 400
    }

    Int cores = 4
    Int preemptible = 1
    Int maxRetries = 5
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, tumor_bam, tumor_bam_bai, region_file], "GB"))

    runtime {
        bootDiskSizeGb: 100
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        docker: "mgibio/cle:v1.4.2"
        memory: "24GB"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        echo -e "~{tumor_bam}\t~{insert_size}\t~{tumor_sample_name}" >> pindel.config

        /usr/bin/pindel -i pindel.config -w 30 -T ~{cores} -o all -f ~{reference} \
        ~{if defined(chromosome) then "-c ~{chromosome}" else ""} \
        ~{if defined(region_file) then "-j ~{region_file}" else ""}
    >>>

    output {
        File deletions = "all_D"
        File insertions = "all_SI"
        File tandems = "all_TD"
        File long_insertions = "all_LI"
        File inversions = "all_INV"
    }
}

task catOut {
  input {
    Array[File] pindel_outs
  }

  Int cores = 1
  Int preemptible = 1
  Int maxRetries = 0
  Int space_needed_gb = 10 + round(size(pindel_outs, "GB")*2)
  runtime {
    memory: "6GB"
    docker: "ubuntu:bionic"
    disks: "local-disk ~{space_needed_gb} SSD"
    cpu: cores
    preemptible: preemptible
    maxRetries: maxRetries
    continueOnReturnCode: [0,1]
  }

  command <<<
    /bin/cat ~{sep=" " pindel_outs} | /bin/grep "ChrID" /dev/stdin > pindel.head
  >>>

  output {
    File pindel_out = "pindel.head"
  }
}

task pindelToVcf {
    input {
        File pindel_output_summary
        File reference
        File reference_fai
        File reference_dict
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? min_supporting_reads = 3
        String? output_name = "pindel.vcf"
        String tumor_sample_name
    }
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, pindel_output_summary], "GB"))
    runtime {
      memory: "6GB"
      docker: "mgibio/cle:v1.3.1"
      disks: "local-disk ~{space_needed_gb} SSD"
      cpu: cores
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        /usr/bin/pindel2vcf -G -p ~{pindel_output_summary} -r ~{reference} -R ~{ref_name} -e ~{min_supporting_reads} -d ~{ref_date} -v ~{output_name}
        # If pindel returns empty pindel.head file, need to account for empty file.
        is_empty=$(grep "~{tumor_sample_name}" ~{output_name})
        if [[ ${is_empty} == "" ]]; then
            grep "##" ~{output_name} > temp.vcf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> temp.vcf
            mv temp.vcf ~{output_name}
        fi
    >>>

    output {
        File pindel_vcf = "~{output_name}"
    }
}

task pindelSomaticFilter {
    input {
        File reference
        File reference_fai
        File reference_dict
        File pindel_output_summary
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = 10 + round(size([reference, reference_fai, reference_dict, pindel_output_summary], "GB"))

    runtime {
        memory: "6GB"
        docker: "mgibio/cle:v1.3.1"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/perl /usr/bin/write_pindel_filter_config.pl ~{pindel_output_summary} ~{reference} $PWD
        /usr/bin/perl /usr/bin/somatic_indelfilter.pl filter.config
    >>>

    output {
        File pindel_vcf = "pindel.out.vcf"
    }
}

task removeEndTags {
    input {
        File vcf
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = 10 + round(size(vcf, "GB")*2)

    runtime {
        memory: "6GB"
        docker: "kboltonlab/bst:latest"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String outfile = "pindel.noend.vcf.gz"

    command <<<
        /usr/local/bin/bgzip -c ~{vcf} > pindel.vcf.gz
        /usr/local/bin/tabix -p vcf pindel.vcf.gz
        /usr/local/bin/bcftools annotate -x INFO/END -Oz -o ~{outfile} pindel.vcf.gz
        /usr/local/bin/tabix -p vcf ~{outfile}
    >>>

    output {
        File processed_vcf = outfile
        File processed_vcf_tbi = "~{outfile}.tbi"
    }
}

# The fake VCF only really needs the basic information, FORMAT ans SAMPLE columns can be empty since
# it will be overwritten by the annotations and filtering
task createFakeVcf {
    input {
        File vcf
        String tumor_sample_name
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(2*size(vcf, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: "6GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
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

    Int space_needed_gb = 10 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: "6GB"
        disks: "local-disk ~{space_needed_gb} SSD"
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
        File? reference_dict
        File bam
        File vcf
        String output_vcf_basename = "fpfilter"
        String sample_name = "TUMOR"
        Float? min_var_freq = 0.05
    }

    Int space_needed_gb = 10 + round(size(vcf, "GB")*2 + size([reference, reference_fai, reference_dict, bam], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: "6GB"
        bootDiskSizeGb: 25
        docker: "kboltonlab/fp_filter-wdl"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = output_vcf_basename + ".vcf"

    command <<<
        /usr/bin/perl /usr/bin/fpfilter.pl --bam-readcount /usr/bin/bam-readcount \
            --samtools /opt/samtools/bin/samtools --output ~{output_vcf} --reference ~{reference} \
            --bam-file ~{bam} --vcf-file ~{vcf} --sample ~{sample_name} --min-var-freq ~{min_var_freq}
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
        File reference_dict
        Pair[File, File] normal_bam
        String? pon_final_name = "pon.pileup"
        File vcf
        Int? mapq = 5
        Int? baseq = 5
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size([normal_bam.left, normal_bam.right], "GB")
    Float vcf_size = size(vcf, "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + vcf_size)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: "24GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{bam_size}"
        echo "VCF: ~{vcf_size}"
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
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.left} --vcf "${vcf_file%.*}" --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
            else
                /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.left} --vcf ~{vcf} --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16
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
    }

    Int space_needed_gb = 10 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: "6GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_file = merged_vcf_basename + ".vcf.gz"

    command <<<
        /usr/local/bin/bcftools merge --output-type z -o ~{output_file} ~{sep=" " vcfs}
        /usr/local/bin/tabix ~{output_file}
    >>>

    output {
        File merged_vcf = output_file
        File merged_vcf_tbi = "~{output_file}.tbi"
    }
}

task vep {
    input {
        File vcf
        File cache_dir_zip
        File reference
        File reference_fai
        File reference_dict
        String ensembl_assembly
        String ensembl_version
        String ensembl_species
        Array[String] plugins
        VepSpliceAIPlugin spliceAI_files = {}
        Boolean coding_only = false
        Array[VepCustomAnnotation] custom_annotations = []
        Boolean everything = true
        # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
        String pick = "flag_pick"
        String additional_args = "--pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --merged --buffer_size 1000 --af_gnomad"
        File? synonyms_file
    }

    Float cache_size = 3*size(cache_dir_zip, "GB")  # doubled to unzip
    Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float splice_AI_size = size([spliceAI_files.spliceAI_indel, spliceAI_files.spliceAI_snv], "GB")
    Int space_needed_gb = 50 + round(reference_size + vcf_size + cache_size + splice_AI_size +size(synonyms_file, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: "64GB"
        bootDiskSizeGb: 30
        cpu: 4
        docker: "kboltonlab/ic_vep"
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String annotated_path = basename(basename(vcf, ".gz"), ".vcf") + "_annotated.vcf"
    String cache_dir = basename(cache_dir_zip, ".zip")
    Int annotation_len = length(custom_annotations)
    File spliceAI_snv = spliceAI_files.spliceAI_snv
    File spliceAI_indel = spliceAI_files.spliceAI_indel

    command <<<
        if [[ ~{annotation_len} -ge 1 ]]; then
            custom_annotation=$(/usr/bin/python3 /opt/bin/jsonToVepString.py ~{write_json(custom_annotations)})
        else
            custom_annotation=""
        fi
        echo $custom_annotation

        echo ~{spliceAI_snv}
        echo ~{spliceAI_indel}

        #mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}
        unzip -qq ~{cache_dir_zip}

        /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
        --format vcf \
        --vcf \
        --fork 4 \
        --terms SO \
        --transcript_version \
        --offline \
        --cache \
        --symbol \
        -o ~{annotated_path} \
        -i ~{vcf} \
        ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
        --sift p \
        --polyphen p \
        ~{if coding_only then "--coding_only" else ""} \
        --~{pick} \
        --dir ~{cache_dir} \
        --fasta ~{reference} \
        ~{sep=" " prefix("--plugin ", plugins)}  \
        ~{if defined(spliceAI_snv) && defined(spliceAI_indel) then "--plugin SpliceAI,snv=~{spliceAI_snv},indel=~{spliceAI_indel}" else ""} \
        ~{if everything then "--everything" else ""} \
        --assembly ~{ensembl_assembly} \
        --cache_version ~{ensembl_version} \
        --species ~{ensembl_species} \
        ~{additional_args} \
        ${custom_annotation}

        bgzip ~{annotated_path} && tabix ~{annotated_path}.gz
    >>>

    output {
        File annotated_vcf = "~{annotated_path}.gz"
        File annotated_vcf_tbi = "~{annotated_path}.gz.tbi"
        File vep_summary = annotated_path + "_summary.html"
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


    Int space_needed_gb = 10 + round(size([vcf, pon, pon_tbi], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/sam_bcftools_tabix_bgzip:1.0"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        if [[ "~{vcf}" == *.gz ]]; then
            name=$(basename ~{vcf} .vcf.gz)
        else
            name=$(basename ~{vcf} .vcf)
        fi

        bcftools +fill-tags -Oz -o RD.vcf.gz ~{pon} -- -t "PON_RefDepth=sum(RD)"
        bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix RD_AD.vcf.gz

        printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
        printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

        sample=`bcftools query -l ~{vcf}`
        bcftools view -H ~{vcf} | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
        bgzip $sample.name;
        tabix $sample.name.gz -s1 -b2 -e2;
        bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE ~{vcf} -Oz -o $name.sample.vcf.gz && tabix $name.sample.vcf.gz;

        ## Varscan has AD and RD instead of comma sep AD field
        ## you can't double quote for string match in bash like you can in zsh so need to make it a variable
        pat="[Vv]arscan"

        ## Lofreq has DP4 which splits into RefFwd, RefRev, AltFwd, AltRev
        patt="[Ll]ofreq"
        if [[ ~{caller} =~ $pat ]]
        then
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "numeric", "numeric"))

            if (length(colnames(df)) != 8) {
            stop("Must supply file with 8 columns: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]", call.=FALSE)
            }

            df$fisher.exact.pval <- apply(df, 1, function(x) {
            x <- as.numeric(x[-c(1,2,3,4)])
            if ((x[1]+x[2]==0) | (x[3]+x[4]==0)){
                return(0)
            } else if (x[2]==0 & x[1]!=0) {
                return(0)
            } else if ((x[1]==0 & x[2]!=0) & (x[3]==0 & x[4]!=0)) {
                return(1)
            } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+ x[4])) {
                return(1)
            } else {
                return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
            }
            })
            df[,2]=sprintf("%1.0f", df[,2])
            write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%RD]\t[%AD]\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        elif [[ ~{caller} =~ $patt ]]
        then
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]

            if (length(colnames(df)) < 8) {
            stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4", call.=FALSE)
            }

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
            df[,2]=sprintf("%1.0f", df[,2])
            write.table(df[, -c(9:10)], file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t%INFO/DP4\n' $name.sample.pileup.vcf.gz > $name.fisher.input;
        else
            echo '
            #!/usr/bin/env Rscript

            args = commandArgs(trailingOnly=TRUE)

            if (length(args)==0) {
            stop("At least one argument must be supplied (input file).n", call.=FALSE)
            } else if (length(args)==1) {
            # default output file
            stop("Must supply (output file).n", call.=FALSE)
            }

            df = read.table(args[1], header=F, colClasses = c("character", "numeric", "character", "character", "numeric", "numeric", "character"))

            #https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
            df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]
            if (length(colnames(df)) != 8) {
            stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]", call.=FALSE)
            }

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
            df[,2]=sprintf("%1.0f", df[,2])
            write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")
            ' > fisherTestInput.R
            bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $name.sample.vcf.gz -Oz -o $name.sample.pileup.vcf.gz;
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

    Int space_needed_gb = 10 + 2*round(size([vcf, vcf_tbi, fp_filter, fp_filter_tbi, vep], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      memory: "6GB"
      docker: "kboltonlab/bst"
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: space_needed_gb
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
        File pon_annotated_vcf = "~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz"
        File pon_annotated_vcf_tbi = "~{caller_prefix}.~{sample_name}.fp_filter.annotated.vcf.gz.tbi"
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
        File TSG_file
        File oncoKB_curated
        File pd_annotation_file
        File truncating
        File cosmic_dir_zip
        String? pon_pvalue = "2.114164905e-6"
    }

    Float caller_size = size([mutect_vcf, lofreq_vcf, vardict_vcf], "GB")
    Float file_size = size([bolton_bick_vars, mut2_bick, mut2_kelly, matches2, TSG_file, oncoKB_curated, pd_annotation_file, truncating], "GB")
    Float cosmic_size = 3*size(cosmic_dir_zip, "GB")
    Int space_needed_gb = 20 + round(caller_size + file_size + cosmic_size)
    Int cores = 2
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/r_docker_ichan:latest"
      memory: "12GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String cosmic_dir = basename(cosmic_dir_zip, ".zip")

    command <<<
        set -eou pipefail

        unzip -qq ~{cosmic_dir_zip}

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{mutect_vcf} --out $(basename ~{mutect_vcf} .vcf.gz) --caller mutect \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --cosmic_dir ~{cosmic_dir} \
        --truncating ~{truncating} \
        --p_value ~{pon_pvalue}
        echo "Mutect AnnotatePD Finished..."

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{lofreq_vcf} --out $(basename ~{lofreq_vcf} .vcf.gz) --caller lofreq \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --cosmic_dir ~{cosmic_dir} \
        --truncating ~{truncating} \
        --p_value ~{pon_pvalue}
        echo "Lofreq AnnotatePD Finished..."

        LC_ALL=C.UTF-8 Rscript --vanilla /opt/bin/ArcherAnnotationScript.R --input ~{vardict_vcf} --out $(basename ~{vardict_vcf} .vcf.gz) --caller vardict \
        --bolton_bick_vars ~{bolton_bick_vars} \
        --mut2_bick ~{mut2_bick} \
        --mut2_kelly ~{mut2_kelly} \
        --matches2 ~{matches2} \
        --TSG_file ~{TSG_file} \
        --oncoKB_curated ~{oncoKB_curated} \
        --pd_annotation_file ~{pd_annotation_file} \
        --cosmic_dir ~{cosmic_dir} \
        --truncating ~{truncating} \
        --p_value ~{pon_pvalue}
        echo "Vardict AnnotatePD Finished..."
    >>>

    output {
        File mutect_vcf_annotate_pd = basename(mutect_vcf, ".vcf.gz") + ".tsv"
        File lofreq_vcf_annotate_pd = basename(lofreq_vcf, ".vcf.gz") + ".tsv"
        File vardict_vcf_annotate_pd = basename(vardict_vcf, ".vcf.gz") + ".tsv"
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
        Boolean model = true
        String tumor_sample_name
    }

    Float caller_size = size([mutect_tsv, lofreq_tsv, vardict_tsv, pindel_full_vcf, pon], "GB")
    Int space_needed_gb = 10 + round(caller_size)
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/xgb:ic_patch"
      memory: "6GB"
      disks: "local-disk ~{space_needed_gb} SSD"
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
