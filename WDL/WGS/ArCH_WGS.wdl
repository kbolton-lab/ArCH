version 1.0

workflow WGS {
    input {
        #General
        String db_path                         # We need a centralized path for the DBs
        Int batch_number

        # Information pertaining to Samples
        File samples_csv
        String sdb_name = "samples.db"

        # Information pertaining to Variants
        String vdb_name = "variants.db"

        # Information pertaining to Callers
        String mutect_cdb_name = "mutect.db"
        String vardict_cdb_name = "vardict.db"
        #Array[Pair[File, File]] vcfs
        Array[File] mutect_vcfs
        Array[File] vardict_vcfs

        # Information pertaining to Annotations
        String adb_name = "annotations.db"
        File? vep_input

        # Reference
        File reference
        File reference_fai
        File reference_dict

        #PoN
        File hiseq_bam_fof
        File novaseq_bam_fof
        String hiseq_pdb_name = "hiseq_pileup.db"
        String novaseq_pdb_name = "novaseq_pileup.db"
    }

    # Array[Pair[File, File]] vcfs = zip(mutect_vcfs, vardict_vcfs)

    call import_samples as samples {
        input:
            db_path = db_path,
            samples_csv = samples_csv,
            samples_db = sdb_name,
            batch_number = batch_number
    }

    call create_sample_blocks {
        input:
            mutect_vcfs = mutect_vcfs,
            vardict_vcfs = vardict_vcfs,
            batch_number = batch_number
    }

    call register_sample_variants as mutect_variants {
        input:
            input_vcfs = create_sample_blocks.mutect_blocks,
            batch_number = batch_number
    }

    call register_sample_variants as vardict_variants {
        input:
            input_vcfs = create_sample_blocks.vardict_blocks,
            batch_number = batch_number
    }

    call merge_batch_variants as variants {
        input:
            sample_variants = flatten([mutect_variants.sample_vdbs, vardict_variants.sample_vdbs]),
            variants_db = vdb_name,
            batch_number = batch_number,
            db_path = db_path
    }

    call dump_variants {
        input:
            db_path = db_path,
            variants_db = vdb_name,
            batch_number = batch_number,
            status_variants = variants.status
    }

    call import_sample_vcf as mutect_sample_vcfs {
        input:
            input_vcfs = create_sample_blocks.mutect_blocks,
            caller = "mutect",
            batch_number = batch_number,
            cores = 4
    }

    call import_sample_vcf as vardict_sample_vcfs {
        input:
            input_vcfs = create_sample_blocks.vardict_blocks,
            caller = "vardict",
            batch_number = batch_number
    }

    call merge_batch_vcfs as merge_mutect_vcfs {
        input:
            vcfs = mutect_sample_vcfs.sample_cdbs,
            db_path = db_path,
            caller_db = mutect_cdb_name,
            caller = "mutect",
            variants_db = vdb_name,
            samples_db = sdb_name,
            batch_number = batch_number,
            status_variants = variants.status
    }

    call merge_batch_vcfs as merge_vardict_vcfs {
        input:
            vcfs = vardict_sample_vcfs.sample_cdbs,
            db_path = db_path,
            caller_db = vardict_cdb_name,
            caller = "vardict",
            variants_db = vdb_name,
            samples_db = sdb_name,
            batch_number = batch_number,
            status_variants = variants.status
    }

    scatter(vcf in dump_variants.chr_vcf) {
        call run_vep {
            input:
                fake_vcf = vcf,
                reference = reference,
                reference_fai = reference_fai,
                reference_dict = reference_dict
        }

        call mskGetBaseCounts as hiseq_pileup{
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bam_fof = hiseq_bam_fof,
            input_vcf = vcf
        }

        call mskGetBaseCounts as novaseq_pileup {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            bam_fof = novaseq_bam_fof,
            input_vcf = vcf
        }
    }

    # Need to merge the pileups for each chromosome
    call bcftoolsConcat as merge_hiseq_pileup {
        input:
        vcfs = hiseq_pileup.pileup,
        vcf_tbis = hiseq_pileup.pileup_tbi
    }

    call bcftoolsConcat as merge_novaseq_pileup {
        input:
        vcfs = novaseq_pileup.pileup,
        vcf_tbis = novaseq_pileup.pileup_tbi
    }

    call merge_vep {
        input:
            vep_tsv = run_vep.vep
    }

    call import_vep {
        input:
            db_path = db_path,
            variants_db = vdb_name,
            annotations_db = adb_name,
            vep = select_first([vep_input, merge_vep.vep]),
            batch_number = batch_number
    }

    call import_pon_pileup as hiseq_pon {
        input:
            db_path = db_path,
            variants_db = vdb_name,
            pileup_db = hiseq_pdb_name,
            pileup = merge_hiseq_pileup.merged_vcf,
            batch_number = batch_number
    }

    call import_pon_pileup as novaseq_pon {
        input:
            db_path = db_path,
            variants_db = vdb_name,
            pileup_db = novaseq_pdb_name,
            pileup = merge_novaseq_pileup.merged_vcf,
            batch_number = batch_number
    }

    call dump_annotations {
        input:
            db_path = db_path,
            annotations_db = adb_name,
            batch_number = batch_number,
            status_vep = import_vep.status
    }

    call run_annotatePD {
        input:
            csv_input = dump_annotations.csv_file
    }

    call import_annotate_pd {
        input:
            db_path = db_path,
            annotations_db = adb_name,
            annotate_pd = run_annotatePD.annotate_pd,
            batch_number = batch_number
    }

    output {
        String done = import_annotate_pd.status
        String variant_database = db_path + "/" + vdb_name
        String sample_database = db_path + "/" + sdb_name
        String annotation_database = db_path + "/" + adb_name
        String mutect_database = db_path + "/" + mutect_cdb_name
        String vardict_database =  db_path + "/" + vardict_cdb_name
        Array[File] chr_vcf = dump_variants.chr_vcf
    }
}

task import_samples {
    input {
        String db_path
        File samples_csv
        String samples_db
        Int batch_number
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "1GB"
        cpu: 1
    }

    command <<<
        ch-toolkit import-samples --samples ~{samples_csv} --sdb ~{db_path}/~{samples_db} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}

task create_sample_blocks {
    input {
        Array[File] mutect_vcfs
        Array[File] vardict_vcfs
        Int batch_number
        Int cores = 32
    }

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: "256GB"
        cpu: cores
    }

    command <<<
        # For Mutect we need to create a SAMPLE variable
        zcat ~{select_first(mutect_vcfs)} | grep "##" > mutect_sample.header
        echo -e "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" >> mutect_sample.header
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> mutect_sample.header
        
        cp /storage1/fs1/bolton/Active/Users/IrenaeusChan/WGS/add_SAMPLE_to_mutect_info_tag.sh add_SAMPLE_to_mutect_info_tag.sh
        cp /storage1/fs1/bolton/Active/Users/IrenaeusChan/WGS/change_sample_name_to_SAMPLE_vardict.sh change_sample_name_to_SAMPLE_vardict.sh

        nProcs=$((~{cores}-1))
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

        find . -type f -name "*.wSAMPLE.vcf.gz*" -delete
    >>>

    output {
        Array[File] mutect_blocks = glob("mutect.block_*.vcf.gz")
        Array[File] vardict_blocks = glob("vardict.block_*.vcf.gz")
    }
}

task register_sample_variants {
    input {
        Array[File] input_vcfs
        Int batch_number
        Int cores = 32
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
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
            ch-toolkit import-sample-variants --input-vcf ${vcf} --vdb ${sample_name}.db --batch-number ~{batch_number} &
        done;
        wait
    >>>

    output {
        Array[File] sample_vdbs = glob("*.db")
    }
}

task merge_batch_variants {
    input {
        Array[File] sample_variants
        Int batch_number
        String db_path
        String variants_db
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "256GB"
        cpu: 1
    }

    command <<<
        mkdir variants
        for variant_db in ~{sep=" " sample_variants}; do
            cp ${variant_db} variants
        done;
        ch-toolkit merge-batch-variants --db-path variants --vdb ~{db_path}/~{variants_db} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}

task dump_variants {
    input {
        String db_path
        String variants_db
        Int batch_number
        String status_variants
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "256GB"
        cpu: 1
    }

    command <<<
        echo ~{status_variants}
        #ch-toolkit dump-variants --vdb ~{db_path}/~{variants_db} --header-type dummy --batch-number ~{batch_number}
        for chr in {1..22} X Y; do
            ch-toolkit dump-variants --vdb ~{db_path}/~{variants_db} --header-type dummy --batch-number ~{batch_number} --chromosome chr${chr}
        done;
    >>>

    output {
        Array[File] chr_vcf = glob("chr*.vcf.gz")
        String status = "Finished."
    }
}

task import_sample_vcf {
    input {
        Array[File] input_vcfs
        String caller
        Int batch_number
        Int cores = 6
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
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
            ch-toolkit import-sample-vcf --caller ~{caller} --input-vcf ${vcf} --cdb ${sample_name}.db --batch-number ~{batch_number} &
        done;
        wait
    >>>

    output {
        Array[File] sample_cdbs = glob("*.db")
    }
}

task merge_batch_vcfs {
    input {
        Array[File] vcfs
        String db_path
        String caller_db
        String caller
        String samples_db
        String variants_db
        Int batch_number
        String status_variants
        Int cores = 8
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "256GB"
        cpu: cores
    }

    command <<<
        echo "Dump Variants is: ~{status_variants}"

        mkdir ~{caller}
        cp ~{sep=" " vcfs} ~{caller}
        ch-toolkit merge-batch-vcf --db-path ~{caller} --cdb ~{db_path}/~{caller_db} --caller ~{caller} --vdb ~{db_path}/~{variants_db} --sdb ~{db_path}/~{samples_db} --batch-number ~{batch_number} --threads ~{cores}
    >>>

    output {
        String status = "Finished."
    }
}

task run_vep {
    input {
        File fake_vcf
        File? reference
        File? reference_fai
        File? reference_dict
        #String vep_dir = "/storage1/fs1/bolton/Active/Users/DucTran/VEP/"
        String vep_dir = "/storage1/fs1/bolton/Active/Users/IrenaeusChan/VEP_cache/"
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float data_size = size(fake_vcf, "GB")
    Int space_needed_gb = ceil(10 + 4 * data_size + reference_size)
    Int cores = 16
    runtime {
        docker: "ensemblorg/ensembl-vep:release_109.3"
        memory: "64GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
    }

    String chrom = basename(fake_vcf, ".vcf.gz")

    command <<<
        outfile=~{chrom}_VEP_annotated.tsv
        
        /opt/vep/src/ensembl-vep/vep -i ~{fake_vcf} --tab -o ${outfile}  \
            --cache --offline --dir_cache ~{vep_dir}/VepData/ --merged --assembly GRCh38 --use_given_ref --species homo_sapiens \
            --symbol --transcript_version --everything --check_existing \
            --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
            --buffer_size 1000  --fork ~{cores} \
            --fasta ~{reference} \
            --pick --pick_order canonical,rank,mane_select,mane_plus_clinical,ccds,appris,tsl,biotype,length \
            --force_overwrite  --no_stats \
            --custom ~{vep_dir}/Clinvar/clinvar_20230430.vcf.gz,clinvar,vcf,exact,0,ID,AF_ESP,AF_EXAC,AF_TGP,CLNSIG,CLNSIGCONF,CLNDN,ORIGIN \
            --dir_plugins ~{vep_dir}/plugin/ \
            --plugin Frameshift --plugin Wildtype \
            --plugin CADD,~{vep_dir}/CADD/whole_genome_SNVs.tsv.gz,~{vep_dir}/CADD/gnomad.genomes.r3.0.indel.tsv.gz \
            --plugin REVEL,~{vep_dir}/REVEL/new_tabbed_revel_grch38.tsv.gz \
            --plugin SpliceAI,snv=~{vep_dir}/spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=~{vep_dir}/spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
            --plugin pLI,~{vep_dir}/pLI/plI_gene.txt
    >>>
            

    output {
        File vep = "~{chrom}_VEP_annotated.tsv"
    }
}

task merge_vep {
    input {
        Array[File] vep_tsv
    }

    runtime {
        docker: "ubuntu:latest"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        grep '#Uploaded_variation' ~{select_first(vep_tsv)} >> VEP_annotated.tsv
        for tsv in ~{sep=" " vep_tsv}; do
            grep -v '^#' ${tsv} >> VEP_annotated.tsv
        done
    >>>

    output {
        File vep = "VEP_annotated.tsv"
    }
}

task import_vep {
    input {
        String db_path
        String variants_db
        String annotations_db
        File vep
        Int batch_number
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "256GB"
        cpu: 1
    }

    command <<<
        ch-toolkit import-vep --vdb ~{db_path}/~{variants_db} --adb ~{db_path}/~{annotations_db} --vep ~{vep} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}

task dump_annotations {
    input {
        String db_path
        String annotations_db
        Int batch_number
        String status_vep
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "256GB"
        cpu: 1
    }

    command <<<
        echo "Importing of VEP is: ~{status_vep}"
        ch-toolkit dump-annotations --adb ~{db_path}/~{annotations_db} --batch-number ~{batch_number}
    >>>

    output {
        File csv_file = "batch-~{batch_number}-forAnnotatePD.csv"
        String status = "Finished."
    }
}

task run_annotatePD {
    input {
        File csv_input
        Int cores = 4
    }

    runtime {
        docker: "indraniel/chip-pipeline-annotation:v1"
        memory: "256GB"
        cpu: cores
    }

    command <<<
        #Ensure each CSV is only around 2 million lines (?)
        tail -n +2 ~{csv_input} | split -d --additional-suffix .csv -l 2000000 - forAnnotatePD.
        for file in forAnnotatePD.*.csv; do
            head -n 1 ~{csv_input} > tmp_file
            cat "$file" >> tmp_file
            mv -f tmp_file "$file"
        done

        nProcs=$((~{cores}-1))
        nJobs="\j"

        for file in forAnnotatePD.*; do
            # Wait until nJobs < nProcs, only start nProcs jobs at most
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done

            /opt/bolton-lab/R-4.2.3/bin/Rscript /storage1/fs1/bolton/Active/Users/IrenaeusChan/ch-toolkit/ch/resources/annotate_pd/run_annotePD.R ${file} $(basename $file .csv).results.csv &
        done;
        wait

        head -n 1 forAnnotatePD.00.results.csv > annotatePD_results.csv

        for file in *.results.csv; do
            tail -n +2 ${file} >> annotatePD_results.csv
        done;
    >>>

    output {
        File annotate_pd = "annotatePD_results.csv"
    }
}


task import_annotate_pd {
    input {
        String db_path
        String annotations_db
        File annotate_pd
        Int batch_number
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "256GB"
        cpu: 1
    }

    command <<<
        ch-toolkit import-annotate-pd --adb ~{db_path}/~{annotations_db} --pd ~{annotate_pd} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}

task mskGetBaseCounts {
    input {
        File reference
        File reference_fai
        File reference_dict
        File bam_fof
        File input_vcf
        Int? mapq = 5
        Int? baseq = 5
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float vcf_size = size([input_vcf], "GB")
    Float data_size = size([bam_fof], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(2 * data_size + vcf_size + reference_size)])
    Float memory = select_first([mem_limit_override, 32])
    Int cores = select_first([cpu_override, 8])
    Int preemptible = 1
    Int maxRetries = 0

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
        # Run the tool
        bgzip -c -d ~{input_vcf} > ~{basename(input_vcf, ".gz")}
        /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta ~{reference} \
            --bam_fof ~{bam_fof} \
            --vcf ~{basename(input_vcf, ".gz")} \
            --output pon.pileup.vcf \
            --maq ~{mapq} --baq ~{baseq}
        bgzip pon.pileup.vcf && tabix pon.pileup.vcf.gz
    >>>

    output {
        File pileup = "pon.pileup.vcf.gz"
        File pileup_tbi = "pon.pileup.vcf.gz.tbi"
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
        /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o merged.vcf.gz ~{sep=" " vcfs}
        /usr/local/bin/tabix merged.vcf.gz
        /usr/local/bin/bcftools +fill-tags -Ov merged.vcf.gz -- -t "PON_RefDepth=sum(RD)" | \
        /usr/local/bin/bcftools +fill-tags -Oz -o pon_pileup.vcf.gz -- -t "PON_AltDepth=sum(AD)"
        /usr/local/bin/tabix pon_pileup.vcf.gz
    >>>

    output {
        File merged_vcf = "pon_pileup.vcf.gz"
        File merged_vcf_tbi = "pon_pileup.vcf.gz.tbi"
    }
}

task import_pon_pileup {
    input {
        String db_path
        String variants_db
        String pileup_db
        File pileup
        Int batch_number
    }

        runtime {
        docker: "kboltonlab/ch-toolkit:v2.6.1"
        memory: "256GB"
        cpu: 1
    }

    command <<<
        ch-toolkit import-pon-pileup --vdb ~{db_path}/~{variants_db} --pdb ~{db_path}/~{pileup_db} --pon-pileup ~{pileup} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}