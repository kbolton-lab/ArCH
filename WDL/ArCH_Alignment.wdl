version 1.0

# ArCH Consensus Building Workflow
# -------------------------------
# Created by: Irenaeus Chan
# Contact: chani@wustl.edu
# Date: 11/10/2022

workflow ArCH_Alignment {
    input {
        # Input Files
        File input_file                     # This can be a FASTQ, BAM, or CRAM
        File input_file_two                 # This will be R2, BAI, or CRAI
        String sample_name

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
    }

    call filterArcherUMILengthAndbbmapRepair as repair {
        input:
        input_file = input_file,
        input_file_two = input_file_two,
        input_type = input_type,
        sample_name = sample_name,
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
        sample_name = sample_name,
        min_reads = min_reads,
        max_read_error_rate = max_read_error_rate,
        max_base_error_rate = max_base_error_rate,
        min_base_quality = min_base_quality,
        max_no_call_fraction = max_no_call_fraction,
        has_umi = has_umi
    }

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
        sample_name = sample_name,
        apply_bqsr = apply_bqsr,
        input_type = input_type
    }

    output {
        # Alignments
        File aligned_bam = initialBAM.initial_bam
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