version 1.0

workflow WGS {
    input {
        #General
        String db_path                         # We need a centralized path for the DBs
        String chip_toolkit
        Int batch_number

        # Information pertaining to Samples
        File samples_csv
        String sdb_name = "samples.db"

        # Information pertaining to Variants
        String vdb_name = "variants.db"
        File? pileup_input
        File? pileup_input_tbi

        # Information pertaining to Callers
        String mutect_cdb_name = "mutect.db"
        String vardict_cdb_name = "vardict.db"
        Array[Pair[File, File]] vcfs

        # Information pertaining to Annotations
        String adb_name = "annotations.db"
        File? vep_input
        File? annotate_pd

    }

    call import_samples as samples {
        input:
            db_path = db_path,
            chip_toolkit = chip_toolkit,
            samples_csv = samples_csv,
            samples_db = sdb_name
    }

    scatter (vcf in vcfs) {
        call register_sample_variants {
            input:
                chip_toolkit = chip_toolkit,
                input_vcf = vcf,
                batch_number = batch_number
        }
    }

    call merge_batch_variants as variants {
        input:
            chip_toolkit = chip_toolkit,
            sample_variants = flatten(register_sample_variants.sample_vdb),
            variants_db = vdb_name,
            batch_number = batch_number,
            db_path = db_path
    }

    call dump_variants {
        input:
            chip_toolkit = chip_toolkit,
            db_path = db_path,
            variants_db = vdb_name,
            batch_number = batch_number,
            status_variants = variants.status
    }

    scatter(vcf in vcfs) {
        File mutect = vcf.left
        File vardict = vcf.right
    }

    # This step can't be parallized because only ONE file can access the DB at a time..
    call import_sample_vcf {
        input:
            chip_toolkit = chip_toolkit,
            mutect_vcfs = mutect,
            vardict_vcfs = vardict,
            db_path = db_path,
            variants_db = vdb_name,
            samples_db = sdb_name,
            batch_number = batch_number,
            status_samples = samples.status,
            status_variants = variants.status,
            status_dump = dump_variants.status
    }

    call merge_batch_vcfs {
        input:
            chip_toolkit = chip_toolkit,
            mutect_vcfs = import_sample_vcf.mutect_cdb,
            vardict_vcfs = import_sample_vcf.vardict_cdb,
            db_path = db_path,
            mutect_caller = mutect_cdb_name,
            vardict_caller = vardict_cdb_name,
            batch_number = batch_number
    }

    scatter(vcf in dump_variants.chr_vcf) {
        call run_vep {
            input:
                fake_vcf = vcf
        }
    }

    call merge_vep {
        input:
            vep_tsv = run_vep.vep,
            vep_input = vep_input
    }

    call import_vep {
        input:
            chip_toolkit = chip_toolkit,
            db_path = db_path,
            variants_db = vdb_name,
            annotations_db = adb_name,
            vep = merge_vep.vep,
            batch_number = batch_number,
            status_import = merge_batch_vcfs.status
    }

    call dump_annotations {
        input:
            chip_toolkit = chip_toolkit,
            db_path = db_path,
            annotations_db = adb_name,
            batch_number = batch_number,
            status_vep = import_vep.status
    }

    call run_annotatePD {
        input:
            csv_input = dump_annotations.csv_file,
            annotate_pd_input = annotate_pd
    }

    call import_annotate_pd {
        input:
            chip_toolkit = chip_toolkit,
            db_path = db_path,
            annotations_db = adb_name,
            annotate_pd = run_annotatePD.annotate_pd,
            batch_number = batch_number
    }

    output {
        String variant_database = db_path + "/" + vdb_name
        String sample_database = db_path + "/" + sdb_name
        String annotation_database = db_path + "/" + adb_name
        String mutect_database = db_path + "/" + mutect_cdb_name
        String vardict_database =  db_path + "/" + vardict_cdb_name
        File variants_vcf = dump_variants.fake_vcf
    }
}

task import_samples {
    input {
        String db_path
        String chip_toolkit
        File samples_csv
        String samples_db
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "1GB"
        cpu: 1
    }

    command <<<
        ~{chip_toolkit} import-samples --samples ~{samples_csv} --sdb ~{db_path}/~{samples_db}
    >>>

    output {
        String status = "Finished."
    }
}

task register_sample_variants {
    input {
        String chip_toolkit
        Pair[File, File] input_vcf
        Int batch_number
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "4GB"
        cpu: 1
    }

    command <<<
        # Mutect
        sample_name=~{basename(input_vcf.left, ".vcf.gz")}.db
        ~{chip_toolkit} import-sample-variants --input-vcf ~{input_vcf.left} --vdb ${sample_name} --batch-number ~{batch_number}

        # Vardict
        sample_name=~{basename(input_vcf.right, ".vcf.gz")}.db
        ~{chip_toolkit} import-sample-variants --input-vcf ~{input_vcf.right} --vdb ${sample_name} --batch-number ~{batch_number}
    >>>

    output {
        Array[File] sample_vdb = glob("*.db")
    }
}

task merge_batch_variants {
    input {
        String chip_toolkit
        Array[File] sample_variants
        Int batch_number
        String db_path
        String variants_db
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        mkdir variants
        cp ~{sep=" " sample_variants} variants
        ~{chip_toolkit} merge-batch-variants --db-path variants --vdb ~{db_path}/~{variants_db} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}

task dump_variants {
    input {
        String chip_toolkit
        String db_path
        String variants_db
        Int batch_number
        String status_variants
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        ~{chip_toolkit} dump-variants --vdb ~{db_path}/~{variants_db} --header-type dummy --batch-number ~{batch_number}
        CHROMS=($(seq 1 22))
        CHROMS+=( "X" "Y" )
        for chr in ${CHROMS[@]}; do
            ~{chip_toolkit} dump-variants --vdb ~{db_path}/~{variants_db} --header-type dummy --batch-number ~{batch_number} --chromosome chr${chr}
        done;
    >>>

    output {
        File fake_vcf = "batch-~{batch_number}.vcf.gz"
        Array[File] chr_vcf = glob("chr*.vcf.gz")
        String status = "Finished."
    }
}

task import_sample_vcf {
    input {
        String chip_toolkit
        #Array[Pair[File, File]] input_vcfs
        #Pair[Array[File], Array[File]] input_vcfs
        Array[File] mutect_vcfs
        Array[File] vardict_vcfs
        String db_path
        String variants_db
        String samples_db
        Int batch_number
        String status_samples
        String status_variants
        String status_dump
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "16GB"
        cpu: 1
    }

    command <<<
        # Mutect
        for vcf in ~{sep=" " mutect_vcfs}; do
            sample_name=$(basename ${vcf} ".vcf.gz").db
            ~{chip_toolkit} import-sample-vcf --input-vcf ${vcf} \
                --cdb ${sample_name} --vdb ~{db_path}/~{variants_db} --sdb ~{db_path}/~{samples_db} \
                --batch-number ~{batch_number} --caller mutect
        done
        #Vardict
        for vcf in ~{sep=" " vardict_vcfs}; do
            sample_name=$(basename ${vcf} ".vcf.gz").db
            ~{chip_toolkit} import-sample-vcf --input-vcf ${vcf} \
                --cdb ${sample_name} --vdb ~{db_path}/~{variants_db} --sdb ~{db_path}/~{samples_db} \
                --batch-number ~{batch_number} --caller vardict
        done
    >>>

    output {
        Array[File] mutect_cdb = glob("mutect*.db")
        Array[File] vardict_cdb = glob("vardict*.db")
    }
}

task merge_batch_vcfs {
    input {
        String chip_toolkit
        Array[File] mutect_vcfs
        Array[File] vardict_vcfs
        String db_path
        String mutect_caller
        String vardict_caller
        Int batch_number
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        mkdir mutect
        mkdir vardict
        cp ~{sep=" " mutect_vcfs} mutect
        cp ~{sep=" " vardict_vcfs} vardict
        ~{chip_toolkit} merge-batch-vcf --db-path mutect --cdb ~{db_path}/~{mutect_caller} --caller mutect --batch-number ~{batch_number}
        ~{chip_toolkit} merge-batch-vcf --db-path vardict --cdb ~{db_path}/~{vardict_caller} --caller vardict --batch-number ~{batch_number}
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
    }

    runtime {
        docker: "kboltonlab/ic_vep:latest"
        memory: "32GB"
        cpu: 4
    }

    String chrom = basename(fake_vcf, ".vcf.gz")

    command <<<
        outfile=~{chrom}_VEP_annotated.tsv
        /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
            --fork 4 -i ~{fake_vcf} --tab -o ${outfile} \
            --offline --cache --buffer_size 1000 \
            --symbol --transcript_version --assembly GRCh38 --cache_version 104 --species homo_sapiens --merged --use_given_ref \
            --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
            --pick --pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length \
            --dir /storage1/fs1/bolton/Active/Projects/mocha/UKBB/ukbb_calls/pvcf/vep_zip \
            --fasta /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/all_sequences.fa \
            --af_gnomad \
            --plugin Frameshift --plugin Wildtype \
            --plugin SpliceAI,snv=/storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/spliceai_scores.raw.indel.hg38.vcf.gz \
            --everything \
            --check_existing --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
            --check_existing --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20181028/custom.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
            --check_existing --custom /storage1/fs1/bolton/Active/Protected/Data/hg38/vcf/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.vcf.gz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
            --force_overwrite --no_stats
    >>>

    output {
        File vep = "~{chrom}_VEP_annotated.tsv"
    }
}

task merge_vep {
    input {
        Array[File] vep_tsv
        File? vep_input
    }

    runtime {
        docker: "ubuntu:latest"
        memory: "8GB"
        cpu: 1
    }

    command <<<
        header=$(grep '#Uploaded_variation' ~{select_first(vep_tsv)})
        echo ${header} >> VEP_annotated.tsv
        for tsv in ~{sep=" " vep_tsv}; do
            grep -v '^#' ${tsv} >> VEP_annotated.tsv
        done
    >>>

    output {
        File vep = select_first([vep_input, "VEP_annotated.tsv"])
    }
}

task import_vep {
    input {
        String chip_toolkit
        String db_path
        String variants_db
        String annotations_db
        File vep
        Int batch_number
        String status_import
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "86GB"
        cpu: 1
    }

    command <<<
        ~{chip_toolkit} import-vep --vdb ~{db_path}/~{variants_db} --adb ~{db_path}/~{annotations_db} --vep ~{vep} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}

task dump_annotations {
    input {
        String chip_toolkit
        String db_path
        String annotations_db
        Int batch_number
        String status_vep
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "4GB"
        cpu: 1
    }

    command <<<
        ~{chip_toolkit} dump-annotations --adb ~{db_path}/~{annotations_db} --batch-number ~{batch_number}
    >>>

    output {
        File csv_file = "batch-~{batch_number}-forAnnotatePD.csv"
        String status = "Finished."
    }
}

task run_annotatePD {
    input {
        File csv_input
        File? annotate_pd_input
    }

    runtime {
        docker: "kboltonlab/r_docker_ichan:latest"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        Rscript /storage1/fs1/bolton/Active/Projects/chip-toolkit/scripts/run_annotePD.R ~{csv_input}
    >>>

    output {
        File annotate_pd = select_first([annotate_pd_input, "annotatePD_results.csv"])
    }
}


task import_annotate_pd {
    input {
        String chip_toolkit
        String db_path
        String annotations_db
        File annotate_pd
        Int batch_number
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "4GB"
        cpu: 1
    }

    command <<<
        ~{chip_toolkit} import-annotate-pd --adb ~{db_path}/~{annotations_db} --pd ~{annotate_pd} --batch-number ~{batch_number}
    >>>

    output {
        String status = "Finished."
    }
}
