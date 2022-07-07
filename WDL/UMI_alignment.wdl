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
        String platform = "Illumina"
        String platform_unit = "Illumina"
        String library = "LIBRARY"
        File? fastq_one
        File? fastq_two
        Boolean bam_input = true           # Default input BAM files, can also be FASTQ files
        File? unaligned_bam
        Boolean is_umi_concensus_unaligned = true       #using consensus bam (unaligned) as input (DEFAULT)
        File? consensus_unaligned_bam 
        Boolean? aligned = false            # If input is an already aligned BAM file then set this flag
        File? aligned_bam_file
        File? aligned_bam_file_bai
        Array[String] read_structure        # Used for the UMI processing see: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
        String tumor_sample_name
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

        # FASTQ Preprocessing
        Boolean? umi_paired = true          # If the UMI is paired (R1 and R2) then set this flag
        Int umi_length = 9
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

        # Variant Calling
        Float? af_threshold = 0.0001                # Minimum VAF Cut-Off

        # Normal BAMs
        Boolean tumor_only = true                   # Defines if Normal BAMs should be used for Variant Calling
        
        # Pindel
        Int pindel_insert_size = 400
        String? ref_name = "GRCh38DH"
        String? ref_date = "20161216"
        Int? pindel_min_supporting_reads = 3

        # See: http://bcb.io/2016/04/04/vardict-filtering/
        # Parameters MQ, NM, DP, and QUAL are calculated using a small subset then identifying the cut-off for 2% of the left side samples
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))"

        # Variants within gnomAD that have a VAF higher than 0.5%
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi
    }

    if (!is_umi_concensus_unaligned) {
    # If the BAM file is already aligned and consensus sequencing was done, then alignment can be skipped
        if (!aligned) {
            # If the input is BAM, we need to check prune reads that have bad UMI information
            if (platform == 'ArcherDX') {
                if (bam_input) {
                    call bamToFastq {
                        input: unaligned_bam = select_first([unaligned_bam, ""])
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
                fastq1 = select_first([fastq_one, ""]),
                fastq2 = select_first([fastq_two, ""]),
                sample_name = tumor_sample_name,
                library_name = library,
                platform_unit = platform_unit,
                platform = platform
            }
        }

        
        # Removes UMIs from Reads and adds them as RX Tag
        call extractUmis {
            input:
                bam = select_first([archer_fastq_to_bam.bam,  unaligned_bam, fastq_to_bam.bam]),
                read_structure = read_structure,
                umi_paired = umi_paired
        }

        # Mark Adapters
        call markIlluminaAdapters {
            input:
                bam = extractUmis.umi_extracted_bam
        }
        
        
        # First Alignment
        call umiAlign as align {
            input:
                bam = markIlluminaAdapters.marked_bam,
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
    }
    
    # Realign the Consensus Called Reads
    call realign {
        input:
            bam = select_first([groupReadsAndConsensus.consensus_bam, consensus_unaligned_bam]),
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            reference_amb = reference_amb,
            reference_ann = reference_ann,
            reference_bwt = reference_bwt,
            reference_pac = reference_pac,
            reference_sa = reference_sa
    }

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
    

    # Applies BQSR on specific intervals defined by the User, if aligned BAM is provided, starts here
    call bqsrApply as bqsr {
        input:
        reference = reference,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        bam = select_first([filterClipAndCollectMetrics.clipped_bam, aligned_bam_file]),
        bam_bai = select_first([filterClipAndCollectMetrics.clipped_bam_bai, aligned_bam_file_bai]),
        intervals = bqsr_intervals,
        known_sites = bqsr_known_sites,
        known_sites_tbi = bqsr_known_sites_tbi
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


        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as mutectSanitizeVcf {
            input: vcf = select_first([mutectTumorTask.vcf,""])
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as mutectNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = mutectSanitizeVcf.sanitized_vcf,
            vcf_tbi = mutectSanitizeVcf.sanitized_vcf_tbi
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


        # Performs the BCBIO Filtering: http://bcb.io/2016/04/04/vardict-filtering/
        call bcftoolsFilterBcbio as bcbio_filter {
            input:
            vcf = select_first([vardictTumorTask.vcf, ""]),
            vcf_tbi = select_first([vardictTumorTask.vcf_tbi, ""]),
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


        # Lofreq does not output FORMAT and SAMPLE columns, so we need to reformat the VCF with these columns
        # since many downstream tools require the VCF to be formatted in this way
        call lofreqReformat as reformat {
            input:
                vcf = select_first([lofreqTumorTask.vcf, ""]),
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

        call removeEndTags {
            input: vcf = select_first([pindelToVcf.pindel_vcf, ""])
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

    }    

    call mergeVcf as merge_mutect_full {
        input:
            vcfs = mutectNormalize.normalized_vcf,
            vcf_tbis = mutectNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "mutect_full." + tumor_sample_name
    }

    call mergeVcf as merge_vardict_full {
        input:
            vcfs = vardictNormalize.normalized_vcf,
            vcf_tbis = vardictNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "vardict_full." + tumor_sample_name
    }

    call mergeVcf as merge_lofreq_full {
        input:
            vcfs = lofreqNormalize.normalized_vcf,
            vcf_tbis = lofreqNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "lofreq_full." + tumor_sample_name
    }

    call mergeVcf as merge_pindel_full {
        input:
            vcfs = pindelNormalize.normalized_vcf,
            vcf_tbis = pindelNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "pindel_full." + tumor_sample_name
    }

    output {
        # Alignments
        File? aligned_bam = filterClipAndCollectMetrics.clipped_bam
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
        File fastqc = fastQC.fastqc

        # variant calling
        File mutect_vcf =   merge_mutect_full.merged_vcf
        File mutect_vcf_tbi = merge_mutect_full.merged_vcf_tbi
        File vardict_vcf = merge_vardict_full.merged_vcf
        File vardict_vcf_tbi = merge_vardict_full.merged_vcf_tbi
        File lofreq_vcf = merge_lofreq_full.merged_vcf
        File lofreq_vcf_tbi = merge_lofreq_full.merged_vcf_tbi
        File pindel_vcf =  merge_pindel_full.merged_vcf
        File pindel_vcf_tbi =  merge_pindel_full.merged_vcf_tbi
    }
}

task filterArcherUmiLength {
    input {
        File fastq1
        File fastq2
        Int umi_length
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([fastq1, fastq2], "GB")
    Int space_needed_gb = 20 + round(3*data_size)

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
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(unaligned_bam, "GB")
    Int space_needed_gb = 5 + round(3*data_size)

    runtime {
        docker: "mgibio/rnaseq:1.0.0"
        memory: "6GB"
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
    }

    Int cores = 1
    Float data_size = size([fastq1, fastq2], "GB")
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/bbmap:38.92--he522d1c_0"
        memory: "12GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(2*data_size)
        disks: "local-disk ~{10 + round(2*data_size)} SSD"
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
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = 10 + round(2*size([fastq1, fastq2], "GB"))

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "6GB"
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
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(2*size(bam, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
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
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(2*size(bam, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "6GB"
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
    }

    Int cores = 8
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")

    runtime {
      docker: "mgibio/dna-alignment:1.0.0"
      memory: "48GB"
      cpu: cores
      # 1 + just for a buffer
      # data_size*10 because bam uncompresses and streams to /dev/stdout and /dev/stdin, could have a couple flying at once
      bootDiskSizeGb: 10 + round(10*data_size + reference_size)
      disks: "local-disk ~{10 + round(10*data_size + reference_size)} SSD"
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
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(3*size(bam, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
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
    }

    Int cores = 8
    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_amb, reference_ann, reference_bwt, reference_pac, reference_sa], "GB")

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory: "48GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(10*data_size + reference_size)
        disks: "local-disk ~{10 + round(10*data_size + reference_size)} SSD"
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
    }

    Int cores = 1
    Float data_size = size(bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "quay.io/biocontainers/fgbio:1.3.0--0"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: 10 + round(3*data_size + reference_size)
        disks: "local-disk ~{10 + round(5*data_size + reference_size)} SSD"
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
    }

    Int cores = 1
    Int space_needed_gb = 10 + round(size(known_sites, "GB") + size(known_sites_tbi, "GB") + size([reference, reference_fai, reference_dict], "GB") + size([bam, bam_bai], "GB") * 2)
    Int preemptible = 1
    Int maxRetries = 0
    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.1.8.1"
        memory: "18GB"
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
    }

    Int space_needed_gb = 10 + round(size([bam, bam_bai, reference, reference_fai, reference_dict],"GB"))
    Int preemptible = 1
    Int maxRetries = 0
    Int cores = 1

    runtime {
        memory: "6GB"
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
    }

    Int space_needed_gb = 10 + round(size([bam, bam_bai, reference, reference_fai, reference_dict, bait_intervals, target_intervals], "GB"))
    Int preemptible = 1
    Int maxRetries = 0
    Int cores = 1

    runtime {
        memory: "6GB"
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
    }

    Int space_needed_gb = 5 + round(size([bam, bam_bai], "GB")*2)
    Int preemptible = 1
    Int maxRetries = 0
    Int cores = 1

    runtime {
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: "4GB"
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
        File fastqc = "~{bamroot}_fastqc.html"
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
        /usr/local/bin/bcftools norm --multiallelics -any --output-type z --output bcftools_norm.vcf.gz ~{vcf} -f ~{reference}
        /usr/local/bin/tabix bcftools_norm.vcf.gz
    >>>

    output {
        File normalized_vcf = "bcftools_norm.vcf.gz"
        File normalized_vcf_tbi = "bcftools_norm.vcf.gz.tbi"
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