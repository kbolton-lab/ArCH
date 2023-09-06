version 1.0

task splitBAMToChr {
    input {
        File bam_file
        File bai_file
        File interval_bed
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([interval_bed, bam_file], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst:latest"
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        intervals=$(awk '{print $1}' ~{interval_bed} | uniq)
        for chr in ${intervals}; do
            samtools view -b ~{bam_file} $chr > ~{basename(bam_file, ".bam")}_${chr}.bam
            samtools index ~{basename(bam_file, ".bam")}_${chr}.bam
        done
    >>>

    output {
        Array[Pair[File, File]] split_bam_chr = zip(glob(basename(bam_file, ".bam")+"_*.bam"), glob(basename(bam_file, ".bam")+"_*.bam.bai"))
    }
}

task splitBedToChr {
    input {
        File interval_bed
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(interval_bed, "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst:latest"
        bootDiskSizeGb: space_needed_gb
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

task echoTest {
    input {
        File interval_bed
        File bam_file
        File bai_file
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([interval_bed, bam_file, bai_file], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst:latest"
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        ls ~{interval_bed}
        ls ~{bam_file}
        ls ~{bai_file}
        echo "It works"
    >>>

    output {
        String result = read_string(stdout())
    }
}


workflow wf {
    input {
        File interval_bed
        File bam_file
        File bai_file
    }

    call splitBedToChr {
        input:
            interval_bed = interval_bed
    }

    call splitBAMToChr {
        input:
            bam_file = bam_file,
            bai_file = bai_file,
            interval_bed = interval_bed
    }

    Array[Pair[File, Pair[File, File]]] splitBedandBAM = zip(splitBedToChr.split_chr, splitBAMToChr.split_bam_chr)

    scatter (bed_bam_chr in splitBedandBAM) {
        call echoTest {
            input:
                interval_bed = bed_bam_chr.left,
                bam_file = bed_bam_chr.right.left,
                bai_file = bed_bam_chr.right.right
        }
    }

    output {
        Array[String] file_output = echoTest.result
    }
}
