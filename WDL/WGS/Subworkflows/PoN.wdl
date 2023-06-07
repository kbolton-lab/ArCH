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

    scatter (bam in pon_bam_bai) {
        call mskGetBaseCounts {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            normal_bam = bam.right,
            pon_final_name = bam.left + ".pon.pileup",
            vcfs = input_vcfs,
            tbis = input_tbi
        }
    }

    call bcftoolsMergePileup as pileup_merge {
        input:
            vcfs = mskGetBaseCounts.pileup,
            vcf_tbis = mskGetBaseCounts.pileup_tbi
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
        Array[File] vcfs
        Array[File] tbis
        Int? mapq = 5
        Int? baseq = 5
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float vcf_size = size(vcfs, "GB")
    Float data_size = size([normal_bam.left, normal_bam.right], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size + vcf_size + reference_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/4 + 5)])
    Int cores = 4
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      docker: "kboltonlab/msk_getbasecounts:3.0"
      cpu: cores
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: space_needed_gb
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        echo "REFERENCE: ~{reference_size}"
        echo "BAM: ~{data_size}"
        echo "SPACE_NEEDED: ~{space_needed_gb}"

        /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates -o combined.vcf ~{sep=" " vcfs}
        sample_name=$(samtools view -H ~{normal_bam.left} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)

        /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} --bam ${sample_name}:~{normal_bam.left} --vcf combined.vcf --output ~{pon_final_name}.vcf --maq ~{mapq} --baq ~{baseq} --thread 16

        bgzip ~{pon_final_name}.vcf && tabix ~{pon_final_name}.vcf.gz
    >>>

    output {
        File pileup = "~{pon_final_name}.vcf.gz"
        File pileup_tbi = "~{pon_final_name}.vcf.gz.tbi"
    }
}

task bcftoolsMergePileup {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(vcfs, "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)])
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
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
