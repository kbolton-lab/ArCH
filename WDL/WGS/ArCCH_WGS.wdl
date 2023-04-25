version 1.0

workflow WGS {
    input {
        #General
        String db_path                         # We need a centralized path for the DBs
        Int batch_number

        # Information pertaining to Samples
        File samples_csv
        String sdb_name = "samples.db"
        #File sdb = "samples.db"

        # Information pertaining to Variants
        String vdb_name = "variants.db"
        #File vdb = "variants.db"
        File? pileup_input
        File? pileup_input_tbi

        # Information pertaining to Callers
        String mutect_cdb_name = "mutect.db"
        String vardict_cdb_name = "vardict.db"
        #File mutect_cdb = "mutect.db"
        #File vardict_cdb = "vardict.db"
        Array[Pair[File, File]] vcfs

        # Information pertaining to Annotations
        String adb_name = "annotations.db"
        #File adb = "annotations.db"
        File vep_input

    }

    call import_samples as samples {
        input:
            db_path = db_path,
            samples_csv = samples_csv,
            samples_db = sdb_name
    }

    output {
        #File variant_database = variants.vdb
        File sample_database = samples.sdb
        #File annotation_database = .adb
        #File mutect_database = mutect.cdb
        #File vardict_database =  vardict.cdb
    }
}

task import_samples {
    input {
        String db_path
        File samples_csv
        String samples_db
    }

    runtime {
        docker: "indraniel/bolton-db-toolkit:v1"
        memory: "1GB"
        cpu: 1
        disks: "local-disk 10GB SSD"
        bootDiskSizeGb: "10GB"
    }

    command <<<
        SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
        VIRTUALENV=${SCRIPT_DIR}/../venv_ic
        BASE_CMD=${VIRTUALENV}/bin/chip-variant-db

        echo ${BASE_CMD}

        ${BASE_CMD} import-samples --samples samples_csv --sdb ~{db_path}/samples_db
    >>>

    output {
        File sdb = "samples.db"
    }
}
