version 1.0

workflow PoN_Pileup {
    input {
        Array[File] pon_bam
        Array[File] pon_bai
        Array[String] sample_name
        Array[File] input_vcfs
        Array[File] input_tbi
        File reference
        File reference_fai
        File reference_dict
    }

    Array[Pair[String, Pair[String, String]]] pon_bam_bai = zip(sample_name, zip(pon_bam, pon_bai))
    Array[Pair[String, String]] vcf_tbi = zip(input_vcfs, input_tbi)

    # For each PoN BAM, get the pileup
    scatter (bam in pon_bam_bai) {
        # Separate the job for each chromosome since the tool reads the entire VCF/BAM into memory
        scatter(vcf in vcf_tbi) {
            call mskGetBaseCounts {
                input:
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict,
                normal_bam = bam.right,
                pon_final_name = bam.left + ".pon.pileup",
                input_vcf = vcf
            }
        }

        # Need to merge the pileups for each chromosome
        call bcftoolsConcat {
            input:
            vcfs = mskGetBaseCounts.pileup,
            vcf_tbis = mskGetBaseCounts.pileup_tbi
        }
    }

    call bcftoolsMergePileup as pileup_merge {
        input:
            vcfs = bcftoolsConcat.merged_vcf,
            vcf_tbis = bcftoolsConcat.merged_vcf_tbi
    }

    output {
        File pileup = pileup_merge.merged_vcf
    }
}

task mskGetBaseCounts {
    input {
        File reference
        File reference_fai
        File reference_dict
        Pair[File, File] normal_bam
        String? pon_final_name = "pon.pileup"
        Pair[File, File] input_vcf
        Int? mapq = 5
        Int? baseq = 5
        Int? mem_limit_override
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float vcf_size = size([input_vcf.left, input_vcf.right], "GB")
    Float data_size = size([normal_bam.left, normal_bam.right], "GB")
    Int space_needed_gb = ceil(2 * data_size + vcf_size + reference_size)
    Float memory = select_first([mem_limit_override, 32])
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} HDD"
      bootDiskSizeGb: 10
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        sample_name=$(samtools view -H ~{normal_bam.left} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)

        # We need to pull the chromosome name from the VCF
        chr=$(zgrep -v '#' ~{input_vcf.left} | awk '{print $1}' | sort | uniq)

        # Split the BAM into the specific chromosome
        samtools view -@ ~{cores} --fast -b -T ~{reference} -o ${sample_name}_${chr}.bam ~{normal_bam.left} ${chr}
        samtools index ${sample_name}_${chr}.bam

        # Run the tool
        bgzip -c -d ~{input_vcf.left} > ~{basename(input_vcf.left, ".gz")}
        /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} \
            --bam ${sample_name}:${sample_name}_${chr}.bam \
            --vcf ~{basename(input_vcf.left, ".gz")} \
            --output ~{pon_final_name}.vcf \
            --maq ~{mapq} --baq ~{baseq}
        bgzip ~{pon_final_name}.vcf && tabix ~{pon_final_name}.vcf.gz
    >>>

    output {
        File pileup = "~{pon_final_name}.vcf.gz"
        File pileup_tbi = "~{pon_final_name}.vcf.gz.tbi"
    }
}

task bcftoolsConcat {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Float data_size = size(vcfs, "GB")
    Int space_needed_gb = ceil(2*data_size)
    Float memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        cpu: cores
        bootDiskSizeGb: 10
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

task bcftoolsMergePileup {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Float data_size = size(vcfs, "GB")
    Int space_needed_gb = ceil(2*data_size)
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
        /usr/local/bin/bcftools merge --output-type z -o merged.vcf.gz ~{sep=" " vcfs}
        /usr/local/bin/tabix merged.vcf.gz
        /usr/local/bin/bcftools +fill-tags -Oz -o RD.vcf.gz ~{output_file} -- -t "PON_RefDepth=sum(RD)"
        /usr/local/bin/bcftools +fill-tags -Oz -o pon_pileup.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix pon_pileup.vcf.gz
    >>>

    output {
        File merged_vcf = "pon_pileup.vcf.gz"
        File merged_vcf_tbi = "pon_pileup.vcf.gz.tbi"
    }
}
