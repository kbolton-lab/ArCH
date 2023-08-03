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
        Array[Pair[File, File]] vcfs

        # Information pertaining to Annotations
        String adb_name = "annotations.db"
        Boolean vep_done = false
        File? vep_input

        # Reference
        File reference
        File reference_fai
        File reference_dict
    }

    call import_samples as samples {
        input:
            db_path = db_path,
            samples_csv = samples_csv,
            samples_db = sdb_name,
            batch_number = batch_number
    }

    scatter (vcf in vcfs) {
        call register_sample_variants as mutect_variants {
            input:
                input_vcf = vcf.left,
                batch_number = batch_number
        }

        call register_sample_variants as vardict_variants {
            input:
                input_vcf = vcf.right,
                batch_number = batch_number
        }
    }

    call merge_batch_variants as variants {
        input:
            sample_variants = flatten([mutect_variants.sample_vdb, vardict_variants.sample_vdb]),
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

    scatter(vcf in vcfs) {
        call import_sample_vcf as mutect_vcfs {
            input:
                input_vcf = vcf.left,
                caller = "mutect",
                batch_number = batch_number,
                status_samples = samples.status,
                status_variants = variants.status
        }

        call import_sample_vcf as vardict_vcfs {
            input:
                input_vcf = vcf.right,
                caller = "vardict",
                batch_number = batch_number,
                status_samples = samples.status,
                status_variants = variants.status
        }
    }

    call merge_batch_vcfs {
        input:
            mutect_vcfs = mutect_vcfs.sample_cdb,
            vardict_vcfs = vardict_vcfs.sample_cdb,
            db_path = db_path,
            mutect_db = mutect_cdb_name,
            vardict_db = vardict_cdb_name,
            variants_db = vdb_name,
            samples_db = sdb_name,
            batch_number = batch_number,
            status_dump = dump_variants.status
    }

    if (!vep_done) {
        scatter(vcf in dump_variants.chr_vcf) {
            call run_vep {
                input:
                    fake_vcf = vcf,
                    reference = reference,
                    reference_fai = reference_fai,
                    reference_dict = reference_dict
            }
        }

        call merge_vep {
            input:
                vep_tsv = run_vep.vep
        }
    }

    call import_vep {
        input:
            db_path = db_path,
            variants_db = vdb_name,
            annotations_db = adb_name,
            vep = select_first([vep_input, merge_vep.vep]),
            batch_number = batch_number,
            status_import = merge_batch_vcfs.status
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
        docker: "kboltonlab/ch-toolkit:v2.2.2"
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

task register_sample_variants {
    input {
        File input_vcf
        Int batch_number
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "4GB"
        cpu: 1
    }

    command <<<
        sample_name=~{basename(input_vcf, ".vcf.gz")}.db
        ch-toolkit import-sample-variants --input-vcf ~{input_vcf} --vdb ${sample_name} --batch-number ~{batch_number}
    >>>

    output {
        File sample_vdb = basename(input_vcf, ".vcf.gz") + ".db"
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
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        mkdir variants
        cp ~{sep=" " sample_variants} variants
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
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "32GB"
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
        File input_vcf
        String caller
        Int batch_number
        String status_samples
        String status_variants
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "16GB"
        cpu: 1
    }

    String sample_name = basename(input_vcf, ".vcf.gz")
    command <<<
        echo "Samples is: ~{status_samples}"
        echo "Variants is: ~{status_variants}" 
        ch-toolkit import-sample-vcf --caller ~{caller} --input-vcf ~{input_vcf} --cdb ~{sample_name}.db --batch-number ~{batch_number}
    >>>

    output {
        File sample_cdb = sample_name + ".db"
    }
}

task merge_batch_vcfs {
    input {
        Array[File] mutect_vcfs
        Array[File] vardict_vcfs
        String db_path
        String mutect_db
        String vardict_db
        String samples_db
        String variants_db
        Int batch_number
        String status_dump
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "86GB"
        cpu: 1
    }

    command <<<
        echo "Dump Variants is: ~{status_dump}"
        mkdir mutect
        mkdir vardict
        cp ~{sep=" " mutect_vcfs} mutect
        cp ~{sep=" " vardict_vcfs} vardict
        ch-toolkit merge-batch-vcf --db-path mutect --cdb ~{db_path}/~{mutect_db} --caller mutect --vdb ~{db_path}/~{variants_db} --sdb ~{db_path}/~{samples_db} --batch-number ~{batch_number}
        ch-toolkit merge-batch-vcf --db-path vardict --cdb ~{db_path}/~{vardict_db} --caller vardict --vdb ~{db_path}/~{variants_db} --sdb ~{db_path}/~{samples_db} --batch-number ~{batch_number}
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
        String vep_dir = "/storage1/fs1/bolton/Active/Users/DucTran/VEP/"
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float data_size = size(fake_vcf, "GB")
    Int space_needed_gb = ceil(10 + 4 * data_size + reference_size)
    Int cores = 4
    runtime {
        docker: "ensemblorg/ensembl-vep:release_109.3"
        memory: "32GB"
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
            --plugin SpliceAI,snv=/storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/spliceai_scores.raw.indel.hg38.vcf.gz \
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
        memory: "8GB"
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
        String status_import
    }

    runtime {
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "86GB"
        cpu: 1
    }

    command <<<
        echo "VCF Importing is: ~{status_import}"
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
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "12GB"
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
    }

    runtime {
        docker: "indraniel/chip-pipeline-annotation:v1"
        memory: "32GB"
        cpu: 2
    }

    command <<<
        #Ensure each CSV is only around 2 million lines (?)
        tail -n +2 ~{csv_input} | split -d --additional-suffix .csv -l 2000000 - forAnnotatePD.
        for file in forAnnotatePD.*.csv; do
            head -n 1 ~{csv_input} > tmp_file
            cat "$file" >> tmp_file
            mv -f tmp_file "$file"
        done

        nProcs=2
        nJobs="\j"

        for file in forAnnotatePD.*; do
            # Wait until nJobs < nProcs, only start nProcs jobs at most
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done

            /opt/bolton-lab/R-4.2.3/bin/Rscript /storage1/fs1/bolton/Active/Projects/chip-toolkit/scripts/run_annotePD.R ${file} $(basename $file .csv).results.csv &
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
        docker: "kboltonlab/ch-toolkit:v2.2.2"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        ch-toolkit import-annotate-pd --adb ~{db_path}/~{annotations_db} --pd ~{annotate_pd} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}

