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
        File vep_input

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
            sample_variants = register_sample_variants.sample_vdb,
            variants_db = vdb_name
    }

    output {
        File variant_database = db_path + vdb_name
        String sample_database = db_path + sdb_name
        #File annotation_database = .adb
        #File mutect_database = mutect.cdb
        #File vardict_database =  vardict.cdb
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
        ~{chip_toolkit} import-samples --samples ~{samples_csv} --sdb ~{db_path}/~{samples_db} --clobber
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
        memory: "8GB"
        cpu: 1
    }

    command <<<
        # Mutect
        mutect_sample_name=~{basename(input_vcf.left, ".vcf.gz")}.db
        ~{chip_toolkit} register-sample-variants --input-vcf ~{input_vcf.left} --db ${mutect_sample_name} --batch-number ~{batch_number} --clobber

        # Vardict
        vardict_sample_name=~{basename(input_vcf.right, ".vcf.gz")}.db
        ~{chip_toolkit} register-sample-variants --input-vcf ~{input_vcf.right} --db ${vardict_sample_name} --batch-number ~{batch_number} --clobber
    >>>

    output {
        Pair[File, File] sample_vdb = ("~{basename(input_vcf.left, '.vcf.gz')}.db", "~{basename(input_vcf.right, '.vcf.gz')}.db")
    }
}

task merge_batch_variants {
    input {
        String chip_toolkit
        Array[Pair[File, File]] sample_variants
        String variants_db
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "32GB"
        cpu: 1
    }

    command <<<
        mkdir mutect
        for sample in ~{sample_variants}; do
            echo ${sample}
        done
    >>>

    output {
        String status = "Finished."
    }
}
