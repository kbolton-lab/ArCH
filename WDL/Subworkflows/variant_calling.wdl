version 1.0

workflow variant_calling {
    input {
        File aligned_bam_file
        File aligned_bai_file
        String tumor_sample_name
        File target_intervals               # Interval List
        Int? mem_limit_override = 6         # Some applications will require more memory depending on BAM size and BED size... (in GB)
                                            # Need to account for these types of errors
        # Reference
        File reference
        File reference_fai
        File reference_dict

        Float? af_threshold = 0.0001                # Minimum VAF Cut-Off

        # See: http://bcb.io/2016/04/04/vardict-filtering/
        # Parameters MQ, NM, DP, and QUAL are calculated using a small subset then identifying the cut-off for 2% of the left side samples
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 6) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 3.0) || (FMT/DP < 6500) || (INFO/QUAL < 27)))"

        # PoN2
        # If a variant exists inside two or more of our Panel of Normal Samples (PoN) at 2% VAF or greater, then this variant is assumed to be
        # sequencing noise because our PoN should not have any variants (that are not germline).
        # Instructions:
        # PoN Samples are run through Tumor Only mode and Filtered for each caller.
        # Variants below 2% VAF are removed from PoN Samples
        # Variants that appear within the PoN Samples in two or more samples are kept
        File mutect_pon2_file
        File mutect_pon2_file_tbi
        File vardict_pon2_file
        File vardict_pon2_file_tbi

        # Variants within gnomAD that have a VAF higher than 0.5%
        File normalized_gnomad_exclude
        File normalized_gnomad_exclude_tbi

    }

    # Some of our callers use BED file instead of interval list
    call intervalsToBed as interval_to_bed {
        input: interval_list = target_intervals
    }

    # In order to parallelize as much as the workflow as possible, we analyze by chromosome
    call splitBedToChr {
        input:
            interval_bed = interval_to_bed.interval_bed
    }

    call splitBAMToChr {
        input:
            bam_file = aligned_bam_file,
            bai_file = aligned_bai_file,
            interval_bed = interval_to_bed.interval_bed
    }

    Array[Pair[File, Pair[File, File]]] splitBedandBAM = zip(splitBedToChr.split_chr, splitBAMToChr.split_bam_chr)

    scatter (bed_bam_chr in splitBedandBAM) {
        # Mutect
        call mutect {
          input:
          reference = reference,
          reference_fai = reference_fai,
          reference_dict = reference_dict,
          gnomad = normalized_gnomad_exclude,
          gnomad_tbi = normalized_gnomad_exclude_tbi,
          tumor_bam = bed_bam_chr.right.left,
          tumor_bam_bai = bed_bam_chr.right.right,
          interval_list = bed_bam_chr.left
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as mutectSanitizeVcf {
            input: vcf = mutect.vcf
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as mutectNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = mutectSanitizeVcf.sanitized_vcf,
            vcf_tbi = mutectSanitizeVcf.sanitized_vcf_tbi
        }

        call bcftoolsIsecComplement as mutect_isec_complement_gnomAD {
            input:
            vcf = mutectNormalize.normalized_vcf,
            vcf_tbi = mutectNormalize.normalized_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "mutect." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }

        call pon2Percent as mutect_pon2 {
            input:
            vcf = mutect_isec_complement_gnomAD.complement_vcf,
            vcf2PON = mutect_pon2_file,
            vcf2PON_tbi = mutect_pon2_file_tbi,
            caller = "mutect",
            sample_name = tumor_sample_name
        }

        # Vardict
        call vardict {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = bed_bam_chr.right.left,
            tumor_bam_bai = bed_bam_chr.right.right,
            interval_bed = bed_bam_chr.left,
            min_var_freq = af_threshold,
            tumor_sample_name = tumor_sample_name
        }

        # Performs the BCBIO Filtering: http://bcb.io/2016/04/04/vardict-filtering/
        call bcftoolsFilterBcbio as bcbio_filter {
            input:
            vcf = vardict.vcf,
            vcf_tbi = vardict.vcf_tbi,
            filter_string = bcbio_filter_string,
            filter_flag = "exclude",
            output_type = "z",
            output_vcf_prefix = "vardict.bcbiofilter"
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as vardictSanitizeVcf {
            input: vcf = bcbio_filter.filtered_vcf
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as vardictNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = vardictSanitizeVcf.sanitized_vcf,
            vcf_tbi = vardictSanitizeVcf.sanitized_vcf_tbi
        }

        call bcftoolsIsecComplement as vardict_isec_complement_gnomAD {
            input:
            vcf = vardictNormalize.normalized_vcf,
            vcf_tbi = vardictNormalize.normalized_vcf_tbi,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            output_vcf_name = "vardict." + tumor_sample_name + ".gnomAD_AF_filter.vcf",
            output_type = "z"
        }

        # Removes any variants that appeared in two of our PoN Samples at 2% VAF or higher
        call pon2Percent as vardict_pon2 {
            input:
            vcf = vardict_isec_complement_gnomAD.complement_vcf,
            vcf2PON = vardict_pon2_file,
            vcf2PON_tbi = vardict_pon2_file_tbi,
            caller = "vardict",
            sample_name = tumor_sample_name
        }
    }

    call mergeVcf as merge_mutect {
        input:
            vcfs = mutect_pon2.annotated_vcf,
            vcf_tbis = mutect_pon2.annotated_vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name
    }

    call mergeVcf as merge_vardict {
        input:
            vcfs = vardict_pon2.annotated_vcf,
            vcf_tbis = vardict_pon2.annotated_vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name
    }

    output {
        File mutect_vcf = merge_mutect.merged_vcf
        File vardict_vcf = merge_vardict.merged_vcf
    }
}

task intervalsToBed {
    input {
        File interval_list
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(interval_list, "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])

    runtime {
        docker: "ubuntu:bionic"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} SSD"
        bootDiskSizeGb: space_needed_gb
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/bin/perl -e '
        use feature qw(say);

        for my $line (<>) {
            chomp $line;
            next if substr($line,0,1) eq q(@); #skip header lines
            my ($chrom, $start, $stop) = split(/\t/, $line);
            say(join("\t", $chrom, $start-1, $stop));
        }' ~{interval_list} > interval_list.bed
    >>>

    output {
        File interval_bed = "interval_list.bed"
    }
}

task indexBam {
    input {
        File input_bam
        String sample_name
        File reference
        File reference_fai
        File reference_dict
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size(input_bam, "GB")
    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Int preemptible = 1
    Int maxRetries = 0
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 4 * data_size + reference_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])

    runtime {
        cpu: cores
        docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bam_link = sub(basename(input_bam), basename(basename(input_bam, ".bam"), ".cram"), sample_name)

    command <<<
        ln -s ~{input_bam} ~{bam_link}
        if [[ ~{bam_link} == *.cram ]]; then
            /usr/local/bin/samtools view -b -T ~{reference} -o ~{sample_name}.bam ~{bam_link}
        fi
        /usr/local/bin/samtools index ~{sample_name}.bam
    >>>

    output {
        File bam = "~{sample_name}.bam"
        File bai = "~{sample_name}.bam.bai"
        #File bai = sub(sub(bam_link, "bam$", "bam.bai"), "cram$", "cram.crai")
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

task mutect {
    input {
        File reference
        File reference_fai
        File reference_dict
        File? gnomad
        File? gnomad_tbi
        File tumor_bam
        File tumor_bam_bai
        File interval_list
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, reference_dict, interval_list], "GB")
    Float data_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size + reference_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/2 + 20)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: cores * memory + "GB"
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = "mutect.filtered.vcf.gz"

    command <<<
        set -o pipefail
        set -o errexit

        THREADS=$((~{cores}*4))

        /gatk/gatk Mutect2 --java-options "-Xmx20g" \
        --native-pair-hmm-threads ${THREADS} \
            -O mutect.vcf.gz \
            -R ~{reference} \
            -L ~{interval_list} \
            -I ~{tumor_bam} \
            ~{"--germline-resource " + gnomad} \
            --read-index ~{tumor_bam_bai} \
            --f1r2-tar-gz mutect.f1r2.tar.gz

        /gatk/gatk LearnReadOrientationModel \
            -I mutect.f1r2.tar.gz \
            -O mutect.read-orientation-model.tar.gz

        /gatk/gatk FilterMutectCalls \
            -R ~{reference} \
            -V mutect.vcf.gz \
            --ob-priors mutect.read-orientation-model.tar.gz \
            -O ~{output_vcf} #Running FilterMutectCalls on the output vcf.
    >>>

    output {
    File vcf = output_vcf
    File vcf_tbi = output_vcf + ".tbi"
    }
}

task vcfSanitize {
    input {
        File vcf
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(vcf, "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])

    runtime {
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        docker: "mgibio/samtools-cwl:1.0.0"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    # outbase should match in script but I don't want to risk changing it yet
    String outbase = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        set -eou pipefail

        # 1) removes lines containing non ACTGN bases, as they conflict with the VCF spec and cause GATK to choke
        # 2) removes mutect-specific format tags containing underscores, which are likewise illegal in the vcf spec
        base=`basename ~{vcf}`
        outbase=`echo $base | perl -pe 's/.vcf(.gz)?$//g'`
        echo "~{vcf}   $base    $outbase"
        if [[ "~{vcf}" =~ ".gz" ]];then
        #gzipped input
        gunzip -c "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
        else
        #non-gzipped input
        cat "~{vcf}" | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" >$outbase.sanitized.vcf
        fi
        /opt/htslib/bin/bgzip $outbase.sanitized.vcf
        /usr/bin/tabix -p vcf $outbase.sanitized.vcf.gz
    >>>

    output {
        File sanitized_vcf = outbase + ".sanitized.vcf.gz"
        File sanitized_vcf_tbi = outbase + ".sanitized.vcf.gz.tbi"
    }
}

task bcftoolsNorm {
    input {
        File reference
        File reference_fai
        File vcf
        File vcf_tbi
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([vcf, vcf_tbi], "GB")
    Float reference_size = size([reference, reference_fai], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size + reference_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/bcftools norm --check-ref w --multiallelics -any --output-type z --output bcftools_norm.vcf.gz ~{vcf} -f ~{reference}
        /usr/local/bin/tabix bcftools_norm.vcf.gz
    >>>

    output {
        File normalized_vcf = "bcftools_norm.vcf.gz"
        File normalized_vcf_tbi = "bcftools_norm.vcf.gz.tbi"
    }
}

task bcftoolsIsecComplement {
    input {
        File vcf
        File vcf_tbi
        File exclude_vcf
        File exclude_vcf_tbi
        String output_type = "z"
        String? output_vcf_name = "bcftools_isec.vcf"
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([vcf, vcf_tbi, exclude_vcf, exclude_vcf_tbi], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/bcftools isec -C -w1 ~{vcf} ~{exclude_vcf} --output-type ~{output_type} --output ~{output_vcf_name}.gz && /usr/local/bin/tabix ~{output_vcf_name}.gz
    >>>

    output {
        File complement_vcf = "~{output_vcf_name}.gz"
        File complement_vcf_tbi = "~{output_vcf_name}.gz.tbi"
    }
}

task pon2Percent {
    input {
        File vcf
        File vcf2PON
        File vcf2PON_tbi
        String caller = "caller"
        String sample_name = "tumor"
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([vcf, vcf2PON, vcf2PON_tbi],"GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      bootDiskSizeGb: space_needed_gb
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        export name=~{caller}.~{sample_name}.pon2.annotated.vcf.gz

        printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
        printf "##INFO=<ID=PON_NAT2_percent,Number=1,Type=Integer,Description=\"Number of samples with variant at >=2 percent\">\n" >> pon2.header;
        printf "##INFO=<ID=PON_MAX_VAF,Number=1,Type=Float,Description=\"The maximum VAF found in the PoN Samples\">\n" >> pon2.header;
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\t%INFO/NS\t%INFO/max_VAF\n" ~{vcf2PON} > normal2.txt
        bgzip -f normal2.txt
        tabix -f -s1 -b2 -e2 normal2.txt.gz
        bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent,PON_NAT2_percent,PON_MAX_VAF ~{vcf} -Oz -o $name
        tabix $name
    >>>

    output {
        File annotated_vcf = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz"
        File annotated_vcf_tbi = "~{caller}.~{sample_name}.pon2.annotated.vcf.gz.tbi"
    }
}

task vardict {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File interval_bed
        Float? min_var_freq = 0.005
        Int? mem_limit_override
        Int? cpu_override
        Int? JavaXmx = 24
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, interval_bed], "GB")
    Float data_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 4 * data_size + reference_size)])
    Int preemptible = 1
    Int maxRetries = 0
    Float memory = select_first([mem_limit_override, ceil(data_size/2 + 6)]) # We want the base to be around 8
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18)*16 else 16])

    runtime {
        docker: "kboltonlab/vardictjava:bedtools"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_name = "vardict"

    command <<<
        set -o pipefail
        set -o errexit

        export VAR_DICT_OPTS='"-Xms256m" "-Xmx~{JavaXmx}g"'
        echo ${VAR_DICT_OPTS}
        echo ~{space_needed_gb}

        samtools index ~{tumor_bam}
        bedtools makewindows -b interval_bed -w 50150 -s 5000 > ~{basename(interval_bed, ".bed")}_windows.bed

        /opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G ~{reference} \
            -X 1 \
            -f ~{min_var_freq} \
            -N ~{tumor_sample_name} \
            -b ~{tumor_bam} \
            -c 1 -S 2 -E 3 -g 4 ~{basename(interval_bed, ".bed")}_windows.bed \
            -th ~{cores} | \
        /opt/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
        /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
            -N "~{tumor_sample_name}" \
            -E \
            -f ~{min_var_freq} > ~{output_name}.vcf

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

task bcftoolsFilterBcbio {
    input {
        File vcf
        File vcf_tbi
        String filter_flag = "exclude"
        String filter_string
        String? output_vcf_prefix = "bcftools_filter"
        String output_type = "z"
        Int? disk_size_override
        Int? mem_limit_override
        Int? cpu_override
    }

    Float data_size = size([vcf, vcf_tbi], "GB")
    Int space_needed_gb = select_first([disk_size_override, ceil(10 + 2 * data_size)])
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "kboltonlab/bst:latest"
      memory: cores * memory + "GB"
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    String ff = if filter_flag == "include" then "-i" else "-e"
    command <<<
        /usr/local/bin/bcftools filter ~{ff} "~{filter_string}" ~{vcf} --output-type ~{output_type} --output ~{output_vcf_prefix}.vcf.gz -s "BCBIO" -m+
        /usr/local/bin/tabix ~{output_vcf_prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{output_vcf_prefix}.vcf.gz"
        File filtered_vcf_tbi = "~{output_vcf_prefix}.vcf.gz.tbi"
    }
}

task mergeVcf {
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
    Float memory = select_first([mem_limit_override, ceil(data_size/6 + 5)]) # We want the base to be around 6
    Int cores = select_first([cpu_override, if memory > 36.0 then floor(memory / 18) else 1])
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
        /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}
        /usr/local/bin/tabix ~{output_file}
    >>>

    output {
        File merged_vcf = output_file
        File merged_vcf_tbi = "~{output_file}.tbi"
    }
}
