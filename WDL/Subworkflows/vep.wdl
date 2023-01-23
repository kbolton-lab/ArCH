# VEP can utilize other files for custom annotations outside of the normal available plugins
# This structure format is ported over from MGI's CWL Pipelines so it could potentially be better integrated
struct VepCustomAnnotation {
    Boolean check_existing
    File custom_file
    String name
    String data_format  # enum, ['bed', 'gff', 'gtf', 'vcf', 'bigwig']
    String method  # enum, ['exact', 'overlap']
    Boolean force_report_coordinates
    Array[String]? vcf_fields
    Array[File]? secondary_files
}

# The SpliceAI Plugin requires two files be provided, so rather than having 4 files passed individually
# it makes things cleaner to have a single structure hold all four.
struct VepSpliceAIPlugin {
    File? spliceAI_snv
    File? spliceAI_snv_tbi
    File? spliceAI_indel
    File? spliceAI_indel_tbi
}


workflow vep {
    input {
        File vcf
        #File sql_db

        # Reference File
        File reference
        File reference_fai
        File reference_dict

        # VEP Parameters
        File vep_cache_dir_zip                      # WDL does not have a Directory Variable, so the entire cache needs to be ZIP
        String vep_ensembl_assembly
        String vep_ensembl_version
        String vep_ensembl_species
        Array[String] vep_plugins = ["Frameshift", "Wildtype"]
        VepSpliceAIPlugin vep_plugin_spliceAI_files = {}
        File? synonyms_file
        Boolean? annotate_coding_only = true
        Array[VepCustomAnnotation] vep_custom_annotations
        String vep_pick = "pick"
        Boolean everything = true
    }

    call dumpVCF {
        input:
            sql_db = sql_db,
            options = options
    }

    call vep {
      input:
          vcf = dumpVCF.vcf,
          cache_dir_zip = vep_cache_dir_zip,
          reference = reference,
          reference_fai = reference_fai,
          reference_dict = reference_dict,
          plugins = vep_plugins,
          spliceAI_files = vep_plugin_spliceAI_files,
          ensembl_assembly = vep_ensembl_assembly,
          ensembl_version = vep_ensembl_version,
          ensembl_species = vep_ensembl_species,
          synonyms_file = synonyms_file,
          custom_annotations = vep_custom_annotations,
          coding_only = annotate_coding_only,
          everything = everything,
          pick = vep_pick
    }

    call annotate_DB {
        input:
        vcf = vep.annotated_vcf,
        sql_db = sql_db,
        which_task = "vep"
    }

    output {
        File sql_db = annotate_DB.sql_db
    }
}

task vep {
    input {
        File vcf
        File cache_dir_zip
        File reference
        File reference_fai
        File reference_dict
        String ensembl_assembly
        String ensembl_version
        String ensembl_species
        Array[String] plugins
        VepSpliceAIPlugin spliceAI_files = {}
        Boolean coding_only = false
        Array[VepCustomAnnotation] custom_annotations = []
        Boolean everything = true
        # one of [pick, flag_pick, pick-allele, per_gene, pick_allele_gene, flag_pick_allele, flag_pick_allele_gene]
        String pick = "flag_pick"
        String additional_args = "--pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --merged --buffer_size 1000 --af_gnomad"
        File? synonyms_file
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float cache_size = 3*size(cache_dir_zip, "GB")  # doubled to unzip
    Float vcf_size = 2*size(vcf, "GB")  # doubled for output vcf
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float splice_AI_size = size([spliceAI_files.spliceAI_indel, spliceAI_files.spliceAI_snv], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + cache_size + vcf_size + reference_size + splice_AI_size + size(synonyms_file, "GB"))])
    Float memory = select_first([mem_limit_override, ceil(vcf_size/2 + 8)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18)*4 else 4])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: cores * memory + "GB"
        cpu: cores
        docker: "kboltonlab/ic_vep"
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String annotated_path = basename(basename(vcf, ".gz"), ".vcf") + "_annotated.vcf"
    String cache_dir = basename(cache_dir_zip, ".zip")
    Int annotation_len = length(custom_annotations)
    File spliceAI_snv = spliceAI_files.spliceAI_snv
    File spliceAI_indel = spliceAI_files.spliceAI_indel

    command <<<
        if [[ ~{annotation_len} -ge 1 ]]; then
            custom_annotation=$(/usr/bin/python3 /opt/bin/jsonToVepString.py ~{write_json(custom_annotations)})
        else
            custom_annotation=""
        fi
        echo $custom_annotation

        echo ~{spliceAI_snv}
        echo ~{spliceAI_indel}

        #mkdir ~{cache_dir} && unzip -qq ~{cache_dir_zip} -d ~{cache_dir}
        unzip -qq ~{cache_dir_zip}

        /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
        --format vcf \
        --vcf \
        --fork 4 \
        --terms SO \
        --transcript_version \
        --offline \
        --cache \
        --symbol \
        -o ~{annotated_path} \
        -i ~{vcf} \
        ~{if defined(synonyms_file) then "--synonyms ~{synonyms_file}" else ""} \
        --sift p \
        --polyphen p \
        ~{if coding_only then "--coding_only" else ""} \
        --~{pick} \
        --dir ~{cache_dir} \
        --fasta ~{reference} \
        ~{sep=" " prefix("--plugin ", plugins)}  \
        ~{if defined(spliceAI_snv) && defined(spliceAI_indel) then "--plugin SpliceAI,snv=~{spliceAI_snv},indel=~{spliceAI_indel}" else ""} \
        ~{if everything then "--everything" else ""} \
        --assembly ~{ensembl_assembly} \
        --cache_version ~{ensembl_version} \
        --species ~{ensembl_species} \
        ~{additional_args} \
        ${custom_annotation}

        bgzip ~{annotated_path} && tabix ~{annotated_path}.gz
    >>>

    output {
        File annotated_vcf = "~{annotated_path}.gz"
        File annotated_vcf_tbi = "~{annotated_path}.gz.tbi"
        File vep_summary = annotated_path + "_summary.html"
    }
}

task dump_VCF {
    input {
        File sql_db
        String options
    }

    runtime {
        memory: cores * memory + "GB"
        cpu: cores
        docker: "kboltonlab/chip-toolkit"
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        chip-variant-db dump-vcf --output-vcf=dump.vcf.gz /path/to/sample.db
    >>>

    output {
        File vcf =
    }
}

task annotate_DB {
    input {
        File vcf
        File sql_db

    }

    runtime {
        memory: cores * memory + "GB"
        cpu: cores
        docker: "kboltonlab/chip-toolkit"
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        chip-variant-db update-db --input-vcf=dump.vcf.gz /path/to/sample.db
    >>>

    output {
        File sql_db =
    }
}
