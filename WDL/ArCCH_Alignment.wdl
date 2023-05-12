version 1.0

# ArCCH BQSR Aligned BAM Creation
# -------------------------------
# Created by: Irenaeus Chan
# Contact: chani@wustl.edu
# Date: 11/10/2022

workflow ArCCH_Alignment {
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
        String sample_name
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
        Boolean has_umi = true
        Boolean? umi_paired = true          # If the UMI is paired (R1 and R2) then set this flag
        String where_is_umi = "T"           # Three options "N = Name", "R = Read", or "T = Tag"
        Int umi_length = 8
        Array[Int] min_reads = [1]          # The minimum number of reads that constitutes a "read family"
        Float? max_read_error_rate = 0.05   # If 5% of the reads within a "read family" does not match, removes the family
        Float? max_base_error_rate = 0.1    # If 10% of the bases are different than the majority, masks the base with N
        Int min_base_quality = 1            # Any base less than QUAL 1 will be masked with N
        Float max_no_call_fraction = 0.5    # Maximum fraction of no-calls (N) in the read after filtering

        # BQSR
        Boolean apply_bqsr = false
        Array[File] bqsr_known_sites
        Array[File] bqsr_known_sites_tbi
        Array[String] bqsr_intervals
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
                sample_name = sample_name,
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
                sample_name = sample_name,
                library_name = library,
                platform_unit = platform_unit,
                platform = platform
            }
        }

        if (has_umi) {
            # Removes UMIs from Reads and adds them as RX Tag
            if (where_is_umi == 'R') {
                call extractUmis {
                    input:
                    bam = select_first([archer_fastq_to_bam.bam, fastq_to_bam.bam, unaligned_bam]),
                    read_structure = read_structure,
                    umi_paired = umi_paired
                }
            }
            if (where_is_umi == 'N') {
                call copyUMIFromReadName {
                    input:
                    bam = select_first([fastq_to_bam.bam, unaligned_bam])
                }
            }

            call markIlluminaAdapters as markAdapters {
                input:
                bam = select_first([extractUmis.umi_extracted_bam, copyUMIFromReadName.umi_extracted_bam, fastq_to_bam.bam, unaligned_bam])
            }

            # First Alignment
            call umiAlign as align {
                input:
                bam = markAdapters.marked_bam,
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
                description = sample_name,
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
                description = sample_name,
                umi_paired = umi_paired
            }
        }
    }

    call indexBam {
        input:
        input_bam = select_first([filterClipAndCollectMetrics.clipped_bam, aligned_bam_file, clipAndCollectMetrics.clipped_bam]),
        sample_name = sample_name
    }

    if (apply_bqsr) {
        # Due to: https://gatk.broadinstitute.org/hc/en-us/community/posts/360075246531-Is-BQSR-accurate-on-Novaseq-6000-
        # Sometimes it is worth skipping this step
        # Applies BQSR on specific intervals defined by the User, if aligned BAM is provided, starts here
        call bqsrApply as bqsr {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bam = indexBam.bam,
            bam_bai = indexBam.bai,
            interval_list = target_intervals,
            known_sites = bqsr_known_sites,
            known_sites_tbi = bqsr_known_sites_tbi,
            output_name = sample_name
        }
    }

    output {
        # Alignments
        File? aligned_bam = align.aligned_bam
        File? aligned_bam_bai = align.aligned_bam_bai
        File consensus_bam = indexBam.bam
        File consensus_bam_bai = indexBam.bai
        File? bqsr_bam = bqsr.bqsr_bam
        File? bqsr_bam_bai = bqsr.bqsr_bam_bai
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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])


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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

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
        /usr/bin/java -Xmx~{memory_total}g -jar /opt/picard/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT I=~{unaligned_bam} F=read1.fastq F2=read2.fastq
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
    Float memory = select_first([mem_limit_override, ceil(data_size/3 + 10)]) # 12
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

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
        repair.sh -Xmx~{memory_total}g repair=t overwrite=true interleaved=false outs=singletons.fq out1=R1.fixed.fastq.gz out2=R2.fixed.fastq.gz in1=~{fastq1} in2=~{fastq2}
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
    Int space_needed_gb = ceil(10 + 4 * data_size)
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

    runtime {
        docker: "mgibio/dna-alignment:1.0.0"
        memory:  cores*memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/java -Xmx~{memory_total}g -jar /opt/picard/picard.jar FastqToSam FASTQ=~{fastq1} FASTQ2=~{fastq2} SAMPLE_NAME=~{sample_name} LIBRARY_NAME=~{library_name} PLATFORM_UNIT=~{platform_unit} PLATFORM=~{platform} OUTPUT=unaligned.bam
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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

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
            /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp ExtractUmisFromBam --molecular-index-tags ZA ZB --single-tag RX --input ~{bam} --read-structure ~{sep=" " read_structure} --output umi_extracted.bam
        else
            /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp ExtractUmisFromBam --molecular-index-tags ZA --single-tag RX --input ~{bam} --read-structure ~{sep=" " read_structure} --output umi_extracted.bam
        fi
    >>>

    output {
        File umi_extracted_bam = "umi_extracted.bam"
    }
}

task copyUMIFromReadName {
    input {
        File bam
        Float? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(10 + 2 * data_size)
    Int preemptible = 1
    Int maxRetries = 0
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

    runtime {
        docker: "quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0"
        memory: cores*memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CopyUmiFromReadName -i ~{bam} -o umi_extracted.bam
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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

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
        /usr/bin/java -Xmx~{memory_total}g -Djava.io.tmpdir=`pwd`/tmp -jar /opt/picard/picard.jar MarkIlluminaAdapters INPUT=~{bam} OUTPUT=marked.bam METRICS=adapter_metrics.txt
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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory/18)*8 else 8])

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
        Int? reads_per_umi_group = 1
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(bam, "GB")
    Int space_needed_gb = ceil(10 + 4 * data_size)
    Float memory = select_first([mem_limit_override, ceil(data_size/3 + 10)]) # 12
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

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
            /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp GroupReadsByUmi --strategy paired --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        else
            /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp GroupReadsByUmi --strategy adjacency --assign-tag MI --raw-tag RX --min-map-q 1 --edits 1 --input $BAM --output umi_grouped.bam
        fi
        /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CallMolecularConsensusReads --input umi_grouped.bam --error-rate-pre-umi 45 --error-rate-post-umi 30 --min-input-base-quality 30 --min-reads ~{reads_per_umi_group} --output consensus_unaligned.bam
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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18)*8 else 8])

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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

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

        /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp FilterConsensusReads --input ~{bam} --output consensus_filtered.bam --ref ~{reference} --min-reads ~{sep=" " min_reads} --max-read-error-rate ~{max_read_error_rate} --max-base-error-rate ~{max_base_error_rate} --min-base-quality ~{min_base_quality} --max-no-call-fraction ~{max_no_call_fraction}
        /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp ClipBam --input consensus_filtered.bam --ref ~{reference} --clipping-mode Hard --clip-overlapping-reads true --output clipped.bam

        PAIRED=~{umi_paired}
        DESCRIPTION=~{description}
        INTERVALS=~{target_intervals}

        if [ "$PAIRED" == true ]; then
            if [[ -z "$DESCRIPTION" ]]; then
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            else
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics --intervals ${INTERVALS}
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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

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

        /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp ClipBam --input ~{bam} --ref ~{reference} --clipping-mode Hard --clip-overlapping-reads true --output clipped.bam

        PAIRED=~{umi_paired}
        DESCRIPTION=~{description}
        INTERVALS=~{target_intervals}

        if [ "$PAIRED" == true ]; then
            if [[ -z "$DESCRIPTION" ]]; then
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --output duplex_seq.metrics --intervals ${INTERVALS}
                fi
            else
                if [[ -z "$INTERVALS" ]]; then
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics
                else
                    /usr/local/bin/fgbio -Xmx~{memory_total}g --tmp-dir=`pwd`/large_tmp CollectDuplexSeqMetrics --input clipped.bam --description ${DESCRIPTION} --output duplex_seq.metrics --intervals ${INTERVALS}
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

task indexBam {
    input {
        File input_bam
        String sample_name
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(input_bam, "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])

    runtime {
        cpu: cores
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bam_link = sub(basename(input_bam), basename(basename(input_bam, ".bam"), ".cram"), sample_name)

    command <<<
        cp ~{input_bam} ~{bam_link}
        /usr/local/bin/samtools index ~{bam_link}
    >>>

    output {
        File bam = bam_link
        File bai = sub(sub(bam_link, "bam$", "bam.bai"), "cram$", "cram.crai")
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
        #Array[String] intervals = ["chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        File interval_list
        Float? mem_limit_override
        Int? cpu_override
        Int? disk_size_override
    }

    Float data_size = size([bam, bam_bai], "GB")
    Float reference_size = size(known_sites, "GB") + size(known_sites_tbi, "GB") + size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 1
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 4 * data_size + reference_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/3 + 10)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int memory_total = floor(memory)-2

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.1.8.1"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /gatk/gatk --java-options -Xmx~{memory_total}g BaseRecalibrator -O bqsr.table -L ~{interval_list} -R ~{reference} -I ~{bam} ~{sep=" " prefix("--known-sites ", known_sites)}
        /gatk/gatk --java-options -Xmx~{memory_total}g ApplyBQSR -O ~{output_name}.bam ~{sep=" " prefix("--static-quantized-quals ", [10, 20, 30])} -R ~{reference} -I ~{bam} -bqsr bqsr.table
    >>>

    output {
        File bqsr_bam = "~{output_name}.bam"
        File bqsr_bam_bai = "~{output_name}.bai"
    }
}
