version 1.0

meta {
    author: "Your Name"
    email: "tran.n@wustl.edu"
    description: "Variant calling workflow with U2AF1 region fixed"
    version: "1.0.0"
    last_modified: "2024-11-14"
}

# WORKFLOW DEFINITION 
# This workflow performs variant calling with U2AF1 region fixed
# Dependencies:
#  - Docker containers: 
#     - kboltonlab/bst:latest
#     - kboltonlab/bwa_sam_gatk:1.0 
#     - broadinstitute/gatk:4.2.0.0
#     - staphb/bcftools:latest
#     - kboltonlab/vardictjava:bedtools
#     - ubuntu:bionic

workflow variant_calling {
    input {
        File aligned_bam_file
        File aligned_bai_file
        String tumor_sample_name
        File target_intervals               # Interval List

        # Reference
        File reference
        File reference_fai
        File reference_dict

        # Reference with U2AF1 fixed
        File reference_u2af1
        File reference_u2af1_fai
        File reference_u2af1_dict
        File exclusion_intervals
        # BWA Index
        File bwa_index_tar

        Float af_threshold = 0.0001                # Minimum VAF Cut-Off

        # See: https://github.com/bcbio/bcbio_validations/blob/master/somatic-lowfreq/README.md
        # Parameters MQ, NM, DP, and QUAL are calculated using a small subset then identifying the cut-off for 2% of the left side samples
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 3) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < 10) || (INFO/QUAL < 30)))"

        # PoN2
        # If a variant exists inside two or more of our Panel of Normal Samples (PoN) at 2% VAF or greater, then this variant is assumed to be
        # sequencing noise because our PoN should not have any variants (that are not germline).
        # Instructions:
        # PoN Samples are run through Tumor Only mode and Filtered for each caller.
        # Variants below 2% VAF are removed from PoN Samples
        # Variants that appear within the PoN Samples in two or more samples are kept
        File? mutect_pon2_file
        File? mutect_pon2_file_tbi
        File? vardict_pon2_file
        File? vardict_pon2_file_tbi

        # Variants within gnomAD that have a VAF higher than 0.5%
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi


    }

    # Some of our callers use BED file instead of interval list
    call intervalsToBed as interval_to_bed {
        input: interval_list = target_intervals
    }

    # In order to parallelize as much as the workflow as possible, we analyze by chromosome
    call splitBedToChr {
        input:
            interval_bed = interval_to_bed.interval_bed
    }

    call bamIndex {
        input:
            bam_file = aligned_bam_file
    }

    call read_extract {
        input:
            aligned_bam_file = aligned_bam_file,
            aligned_bai_file = bamIndex.bam_index,
            reference = reference,
            reference_fai = reference_fai,
            exclusion_intervals = exclusion_intervals
    }

    call bwa_align {
        input:
            extracted_fastq = read_extract.fastq,
            rg = read_extract.rg,
            reference_u2af1 = reference_u2af1,
            reference_u2af1_fai = reference_u2af1_fai,
            bwa_index_tar = bwa_index_tar
    }

    call merge_bam {
        input:
            realigned_bam = bwa_align.realigned_bam,
            aligned_bam_file = aligned_bam_file,
            aligned_bai_file = aligned_bai_file,
            reference = reference,
            reference_fai = reference_fai,
            interval_bed = interval_to_bed.interval_bed
    }

    # scatter (bed_chr in splitBedToChr.split_chr) {
    #     # Mutect
    #     call mutect {
    #         input:
    #         reference = reference,
    #         reference_fai = reference_fai,
    #         reference_dict = reference_dict,
    #         gnomad = normalized_gnomad_exclude,
    #         gnomad_tbi = normalized_gnomad_exclude_tbi,
    #         tumor_bam = aligned_bam_file, 
    #         tumor_bam_bai = bamIndex.bam_index, 
    #         interval_list = bed_chr
    #     }

    #     call mutect_pass {
    #         input:
    #         mutect_vcf = mutect.vcf,
    #         mutect_vcf_tbi = mutect.vcf_tbi
    #     }

    #     # Vardict
    #     call vardict {
    #         input:
    #         reference = reference,
    #         reference_fai = reference_fai,
    #         tumor_bam = aligned_bam_file, 
    #         tumor_bam_bai = bamIndex.bam_index,
    #         interval_bed = bed_chr,
    #         min_var_freq = af_threshold,
    #         tumor_sample_name = tumor_sample_name,
    #         mutect_vcf = mutect_pass.vcf
    #     }
    # }

    # call mergeVcf as merge_mutect {
    #     input:
    #         vcfs = mutect.vcf,
    #         vcf_tbis = mutect.vcf_tbi,
    #         merged_vcf_basename = "mutect." + tumor_sample_name
    # }

    # call mergeVcf as merge_vardict {
    #     input:
    #         vcfs = vardict.vcf,
    #         vcf_tbis = vardict.vcf_tbi,
    #         merged_vcf_basename = "vardict." + tumor_sample_name
    # }

    # call sanitizeNormalizeFilter as mutect_filter {
    #     input:
    #         vcf = merge_mutect.merged_vcf,
    #         vcf_tbi = merge_mutect.merged_vcf_tbi,
    #         reference = reference,
    #         reference_fai = reference_fai,
    #         exclude_vcf = normalized_gnomad_exclude,
    #         exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
    #         vcf2PON = mutect_pon2_file,
    #         vcf2PON_tbi = mutect_pon2_file_tbi,
    #         caller = "mutect",
    #         sample_name = tumor_sample_name,
    #         filter_string = bcbio_filter_string
    # }

    # call sanitizeNormalizeFilter as vardict_filter {
    #     input:
    #         vcf = merge_vardict.merged_vcf,
    #         vcf_tbi = merge_vardict.merged_vcf_tbi,
    #         reference = reference,
    #         reference_fai = reference_fai,
    #         exclude_vcf = normalized_gnomad_exclude,
    #         exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
    #         vcf2PON = vardict_pon2_file,
    #         vcf2PON_tbi = vardict_pon2_file_tbi,
    #         caller = "vardict",
    #         sample_name = tumor_sample_name,
    #         filter_string = bcbio_filter_string
    # }

    output {
        # File mutect_vcf = mutect_filter.annotated_vcf
        # File vardict_vcf = vardict_filter.annotated_vcf
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

task bamIndex {
    input {
        File bam_file
    }

    Float data_size = size([bam_file], "GB")
    Int space_needed_gb = ceil(10 + data_size)
    Int memory = 4
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst:latest"
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bam_link = basename(bam_file)

    command <<<
        set -e -x -o pipefail

        ln -s ~{bam_file} ~{bam_link}
        samtools index ~{bam_link}
    >>>

    output {
        File bam_index = sub(sub(bam_link, "bam$", "bam.bai"), "cram$", "cram.crai")
    }
}

task read_extract {
    input {
        File aligned_bam_file 
        File aligned_bai_file
        File reference
        File reference_fai
        File exclusion_intervals
    }

    Float reference_size = size([reference, reference_fai], "GB")
    Float data_size = size([aligned_bam_file, aligned_bai_file], "GB")
    Int space_needed_gb = ceil(10 + reference_size + data_size)
    Int memory = 4
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -e -x -o pipefail

        exclusionRegions=$(grep -v "#" ~{exclusion_intervals} | awk '{print $1":"$2"-"$3}' | tr '\n' ' ')

        # Extract the BAM corresponds to the input bed regions and excludes reads unmapped, mate unmapped, not primary alignment, read fails platform/vendor quality checks, read is PCR or optical duplicate, and supllementary alignment (-F 3852)
        samtools view -hb -T ~{reference}  -t ~{reference_fai} -@ ~{cores} -F 3852 -X ~{aligned_bam_file} ~{aligned_bai_file} -o extracted_reads.bam $exclusionRegions

        # Read names unique in file and extract BAM corresponds to paired-end reads
        samtools view -hb -@ ~{cores} extracted_reads.bam | \
            samtools sort -@ ~{cores} -n -  > extracted_reads_sorted_by_read_names.bam 
        # Extract the read names that are not paired
        samtools view -@ ~{cores} extracted_reads_sorted_by_read_names.bam | cut -f1 | sort | uniq -c | awk '$1!=2 {print $2}' > nonpairs_rnames.txt
        # Extract the paired-end reads to a new BAM file
        samtools view -h -@ ~{cores} extracted_reads_sorted_by_read_names.bam | \
            grep -F -vf nonpairs_rnames.txt | \
            samtools view -hb -@ ~{cores} - > pairs_only.bam 
        # Sort the BAM file by read names
        samtools sort -m ~{memory}G -@ ~{cores} -n -o pairs_only_sorted_by_read_names.bam pairs_only.bam

        # Extract the paired-end reads
        samtools view -H pairs_only_sorted_by_read_names.bam | grep "^@RG" | sed 's/	/\\t/g'  > extracted_reads_RG.txt
        samtools fastq -@ ~{cores} pairs_only_sorted_by_read_names.bam -1 extract_1.fastq -2 extract_2.fastq -0 /dev/null -s /dev/null -n
    >>>

    output {
        Array[File] fastq = ["extract_1.fastq", "extract_2.fastq"]
        File rg = "extracted_reads_RG.txt"
    }
}

task bwa_align {
    input {
        Array[File] extracted_fastq
        File rg 
        File reference_u2af1
        File reference_u2af1_fai
        File bwa_index_tar
    }

    Int preemptible = 1
    Int maxRetries = 0
    Int memory = 4
    Int cores = 2
    Float reference_size = size([reference_u2af1, bwa_index_tar], "GB")
    Int space_needed_gb = ceil(20 + reference_size)

    runtime {
        docker: "kboltonlab/bwa_sam_gatk:1.0"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -e -x -o pipefail

        # Ln the reference
        ln -s ~{reference_u2af1} .
        ln -s ~{reference_u2af1_fai} .
        # cp ~{reference_u2af1} .
        # cp ~{reference_u2af1_fai} .

        # Unpack the bwa index
        tar -xvzf ~{bwa_index_tar}
        
        # Get reference name
        reference_u2af1_name=$(basename ~{reference_u2af1})

        # Align the reads
        # Keep only primary alignments, i.e., remove flag 0x100 (-F)
        #! 4GB is not enough for this step
        bwa mem -t ~{cores} -R $(head -n 1 ~{rg}) -M ${reference_u2af1_name} ~{extracted_fastq[0]} ~{extracted_fastq[1]} | \
            samtools view -hb -F 0x100 - > realigned.bam
        samtools sort -m ~{memory}G -o realigned.sorted.bam realigned.bam
        samtools index realigned.sorted.bam
        # Keep read from U2AF1 region (chr21:43082956-43117578) plus minus ~50kb, chr21:43,050,000-43,150,000
        # Minimize the change to the original BAM file
        samtools view -hb realigned.sorted.bam chr21:43,050,000-43,150,000 > realigned_u2af1.bam
        # Get all non-paired reads
        samtools view -@ ~{cores} realigned_u2af1.bam | cut -f1 | sort | uniq -c | awk '$1!=2 {print $2}' > nonpairs_rnames.txt
        # Extract the paired-end reads to a new BAM file
        samtools view -h -@ ~{cores} realigned_u2af1.bam | \
            grep -F -vf nonpairs_rnames.txt | \
            samtools view -hb -@ ~{cores} - > realigned_u2af1_paired.bam
    >>>

    output {
        File realigned_bam = "realigned_u2af1_paired.bam"
    }
}

task merge_bam {
    input {
        File realigned_bam
        File aligned_bam_file 
        File aligned_bai_file
        File reference
        File reference_fai
        File interval_bed
    }

    Int preemptible = 1
    Int maxRetries = 0
    Int memory = 6
    Int cores = 4
    Float data_size = size([aligned_bam_file, aligned_bai_file], "GB")
    Float reference_size = size([reference], "GB")
    Int space_needed_gb = ceil(20 + data_size*3 + reference_size)

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -e -x -o pipefail
 
        # Extract read name from the realigned bam
        samtools view ~{realigned_bam} | cut -f1 | sort | uniq > realigned_read_names.txt
        # Remove the realigned reads from the original BAM file
        # Extract the reads that are not realigned
        # Then merge the realigned reads with the original BAM file
        #! Some read names have additional non-primary alignments in the original BAM file, so after removing the realigned reads by read names, the number of reads in the merged BAM file is less than the original BAM file. But it should not affect the downstream analysis.
        samtools view -h -@ ~{cores} -T ~{reference} -t ~{reference_fai} ~{aligned_bam_file} | \
            grep -F -vf realigned_read_names.txt | \
            samtools merge -u -@ ~{cores} - - ~{realigned_bam} | \
            samtools sort -u -m 4G -@ ~{cores} - | \
            samtools view -hC -1 -@ ~{cores} -T ~{reference} -t ~{reference_fai} - > fixed.sorted.cram 
        samtools index fixed.sorted.cram
    >>>

    output {
        File merged_bam = "fixed.sorted.cram"
        File merged_bai = "fixed.sorted.cram.bai"
    }
}




task mutect {
    input {
        File reference
        File reference_fai
        File reference_dict
        File? gnomad
        File? gnomad_tbi
        File tumor_bam
        File tumor_bam_bai
        File interval_list
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, reference_dict, interval_list], "GB")
    Float data_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = ceil(10 + data_size + reference_size)
    Int memory = select_first([mem_limit_override, 4])
    Int cores = select_first([cpu_override, 1])
    Int preemptible = 3
    Int maxRetries = 3

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: cores * memory + "GB"
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = "mutect.filtered.vcf.gz"

    command <<<
        set -e -x -o pipefail

        /gatk/gatk Mutect2 --java-options "-Xmx~{memory}g" \
        --native-pair-hmm-threads ~{cores} \
            -O mutect.vcf.gz \
            -R ~{reference} \
            -L ~{interval_list} \
            -I ~{tumor_bam} \
            ~{"--germline-resource " + gnomad} \
            --read-index ~{tumor_bam_bai} \
            --f1r2-tar-gz mutect.f1r2.tar.gz

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
        File? vcf2PON
        File? vcf2PON_tbi
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

        # gnomAD Intersection
        bcftools isec -C -w1 norm.vcf.gz ~{exclude_vcf} --output-type z --output isec.vcf.gz

        # PoN2
        export name=~{caller}.~{sample_name}.vcf.gz
        if [[ -n ~{vcf2PON} ]]; then
            printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
            printf "##INFO=<ID=PON_NAT2_percent,Number=1,Type=Integer,Description=\"Number of samples with variant at >=2 percent\">\n" >> pon2.header;
            printf "##INFO=<ID=PON_MAX_VAF,Number=1,Type=Float,Description=\"The maximum VAF found in the PoN Samples\">\n" >> pon2.header;
            bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\t%INFO/NS\t%INFO/max_VAF\n" ~{vcf2PON} > normal2.txt
            bgzip -f normal2.txt
            tabix -f -s1 -b2 -e2 normal2.txt.gz
            bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent,PON_NAT2_percent,PON_MAX_VAF isec.vcf.gz -Oz -o $name
        else
            cp isec.vcf.gz $name
        fi
        tabix $name
    >>>

    output {
        File annotated_vcf = "~{caller}.~{sample_name}.vcf.gz"
        File annotated_vcf_tbi = "~{caller}.~{sample_name}.vcf.gz.tbi"
    }
}

task mutect_pass {
    input {
        File mutect_vcf
        File mutect_vcf_tbi
    }

    Float data_size = size([mutect_vcf, mutect_vcf_tbi], "GB")
    Int space_needed_gb = ceil(10 + data_size)
    Int memory = 2
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "staphb/bcftools:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        bcftools filter -i 'FILTER="PASS" || FILTER="weak_evidence" || FILTER="strand_bias" || FILTER="weak_evidence;strand_bias" || (CHROM="chr20" && POS=32434638)' ~{mutect_vcf} > mutect_passed.vcf
        bcftools sort -m ~{memory}G mutect_passed.vcf -Oz -o mutect_passed.sorted.vcf.gz 
    >>>

    output {
        File vcf = "mutect_passed.sorted.vcf.gz"
    }
}

task vardict {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File interval_bed
        File mutect_vcf
        Float min_var_freq = 0.005
        Int JavaXmx = 24
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, interval_bed], "GB")
    Float data_size = size([tumor_bam, tumor_bam_bai, mutect_vcf], "GB")
    Int space_needed_gb = ceil(10 + 2.5 * data_size + reference_size)
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

    String output_name = "vardict"

    command <<<
        set -e -x -o pipefail

        echo ~{space_needed_gb}

        # TODO: Account for when Mutect File is "Empty"..

        # Make windows and intersect with vcf file
        bedtools makewindows -b ~{interval_bed} -w 1150 -s 1000 > ~{basename(interval_bed, ".bed")}_windows.bed
        bedtools intersect -u -wa -a ~{basename(interval_bed, ".bed")}_windows.bed -b ~{mutect_vcf} > interval_list_mutect.bed

        # Merge small intervals
        bedtools merge -i interval_list_mutect.bed > interval_list_mutect_merged.bed
        # intersect with bed file
        samtools view -u -b -M -L interval_list_mutect_merged.bed -T ~{reference} -t ~{reference_fai} \
            -o intersected.bam -X ~{tumor_bam} ~{tumor_bam_bai}
        samtools index intersected.bam

        # Make windows
        bedtools makewindows -b interval_list_mutect_merged.bed -w 10150 -s 10000 > interval_list_mutect_merged_windows.bed

        # Split bed file into 16 equal parts
        split -d --additional-suffix .bed -n l/16 interval_list_mutect_merged_windows.bed splitBed.

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
                -b intersected.bam \
                -c 1 -S 2 -E 3 -g 4 ${fName} \
                --deldupvar -Q 10 -F 0x700 --fisher > result.${part}.txt &
        done;
        # Wait for all running jobs to finish
        wait

        for fName in result.*.txt; do
            cat ${fName} >> resultCombine.txt
        done;

        cat resultCombine.txt | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
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

task mergeVcf {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Float data_size = size(vcfs, "GB")
    Float data_size_tbi = size(vcf_tbis, "GB")
    Int space_needed_gb = ceil(10 + 2 * (data_size + data_size_tbi))
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
        /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}
        /usr/local/bin/tabix ~{output_file}
    >>>

    output {
        File merged_vcf = output_file
        File merged_vcf_tbi = "~{output_file}.tbi"
    }
}

