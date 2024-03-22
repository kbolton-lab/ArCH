version 1.0

workflow test_block{
    input {
        Array[File] mutect_vcfs
        Array[File] vardict_vcfs
        Int batch_number
    }

    call create_sample_blocks {
        input: 
            mutect_vcfs = mutect_vcfs, 
            vardict_vcfs = vardict_vcfs, 
            batch_number = batch_number
    }

    call register_sample_variants_ as mutect_variants {
        input:
            input_vcfs = create_sample_blocks.mutect_blocks,
            batch_number = batch_number
    }

    call register_sample_variants_ as vardict_variants {
        input:
            input_vcfs = create_sample_blocks.vardict_blocks,
            batch_number = batch_number
    }

    call import_sample_vcf_ as mutect_sample_vcfs {
        input:
            input_vcfs = create_sample_blocks.mutect_blocks,
            caller = "mutect",
            batch_number = batch_number
    }

    call import_sample_vcf_ as vardict_sample_vcfs {
        input:
            input_vcfs = create_sample_blocks.vardict_blocks,
            caller = "vardict",
            batch_number = batch_number
    }

    output {
        Array[File] mutect_variants_db = mutect_variants.sample_vdbs
        Array[File] vardict_variants_db = vardict_variants.sample_vdbs
        Array[File] mutect_sample_vcfs_db = mutect_sample_vcfs.sample_cdbs
        Array[File] vardict_sample_vcfs_db = vardict_sample_vcfs.sample_cdbs
    }
}

task create_sample_blocks {
    input {
        Array[File] mutect_vcfs
        Array[File] vardict_vcfs
        Int batch_number
        Int cores = 16
    }

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: "128GB"
        cpu: cores
    }

    command <<<
        # For Mutect we need to create a SAMPLE variable
        zcat ~{select_first(mutect_vcfs)} | grep "##" > mutect_sample.header
        echo -e "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" >> mutect_sample.header
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> mutect_sample.header
        
        cp /storage1/fs1/bolton/Active/Users/IrenaeusChan/WGS/add_SAMPLE_to_mutect_info_tag.sh add_SAMPLE_to_mutect_info_tag.sh
        cp /storage1/fs1/bolton/Active/Users/IrenaeusChan/WGS/change_sample_name_to_SAMPLE_vardict.sh change_sample_name_to_SAMPLE_vardict.sh

        nProcs=~{cores}
        nJobs="\j"
        # Parallelize the creation of the Mutect VCFs
        for vcf in ~{sep=" " mutect_vcfs}; do
            echo ${vcf}
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done
            bash add_SAMPLE_to_mutect_info_tag.sh ${vcf} &
        done;

        # Do the same for Vardict
        for vcf in ~{sep=" " vardict_vcfs}; do
            echo ${vcf}
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done
            bash change_sample_name_to_SAMPLE_vardict.sh ${vcf} &
        done;
        wait
        
        # Here we have SAMPLE in both Vardict and Mutect, so we can just concatenate them into 1Gb blocks
        group=(); group_size=0; group_num=1
        for vcf in mutect.*.wSAMPLE.vcf.gz; do
            echo ${vcf}
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done
            file_size=$(du -b "$vcf" | cut -f1)
            if (( group_size + file_size > 1073741824 )); then
                /usr/local/bin/bcftools concat --allow-overlaps -Oz -o mutect.block_${group_num}.vcf.gz "${group[@]}" &
                group=("$vcf")
                group_size=$file_size
                group_num=$((group_num + 1))
            else
                group+=("$vcf")
                group_size=$((group_size + file_size))
            fi
        done;
        # Concatenate the last group
        /usr/local/bin/bcftools concat --allow-overlaps -Oz -o mutect.block_${group_num}.vcf.gz "${group[@]}" &

        group=()
        group_size=0
        group_num=1
        for vcf in vardict.*.wSAMPLE.vcf.gz; do
            echo ${vcf}
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done
            file_size=$(du -b "$vcf" | cut -f1)
            if (( group_size + file_size > 1073741824 )); then
                /usr/local/bin/bcftools concat --allow-overlaps -Oz -o vardict.block_${group_num}.vcf.gz "${group[@]}" &
                group=("$vcf")
                group_size=$file_size
                group_num=$((group_num + 1))
            else
                group+=("$vcf")
                group_size=$((group_size + file_size))
            fi
        done;
        # Concatenate the last group
        /usr/local/bin/bcftools concat --allow-overlaps -Oz -o vardict.block_${group_num}.vcf.gz "${group[@]}" &
        wait

        rm *.wSAMPLE.vcf.gz*
    >>>

    output {
        Array[File] mutect_blocks = glob("mutect.block_*.vcf.gz")
        Array[File] vardict_blocks = glob("vardict.block_*.vcf.gz")
    }
}

task register_sample_variants_ {
    input {
        Array[File] input_vcfs
        Int batch_number
        Int cores = 32
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.4.2"
        memory: "256GB"
        cpu: cores
    }

    command <<<
        nProcs=$((~{cores}-1))
        nJobs="\j"
        #nJobs=$(jobs -p | wc -l)

        for vcf in ~{sep=" " input_vcfs}; do
            echo ${vcf}
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done
            sample_name=$(basename ${vcf} .vcf.gz)
            echo ${sample_name}
            /storage1/fs1/bolton/Active/Users/IrenaeusChan/ch-toolkit/venv/bin/ch-toolkit import-sample-variants --input-vcf ${vcf} --vdb ${sample_name}.db --batch-number ~{batch_number} &
        done;
        wait
    >>>

    output {
        Array[File] sample_vdbs = glob("*.db")
    }
}

task import_sample_vcf_ {
    input {
        Array[File] input_vcfs
        String caller
        Int batch_number
        Int cores = 32
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.4.2"
        memory: "172GB"
        cpu: cores
    }

    command <<<
        nProcs=~{cores}
        nJobs="\j"
        #nJobs=$(jobs -p | wc -l)

        for vcf in ~{sep=" " input_vcfs}; do
            echo ${vcf}
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done
            sample_name=$(basename ${vcf} .vcf.gz)
            echo ${sample_name}
            /storage1/fs1/bolton/Active/Users/IrenaeusChan/ch-toolkit/venv/bin/ch-toolkit import-sample-vcf --caller ~{caller} --input-vcf ${vcf} --cdb ${sample_name}.db --batch-number ~{batch_number} &
        done;
        wait
    >>>

    output {
        Array[File] sample_cdbs = glob("*.db")
    }
}