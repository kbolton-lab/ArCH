version 1.0

# ArCCH BQSR Aligned BAM Creation
# -------------------------------
# Created by: Irenaeus Chan
# Contact: chani@wustl.edu
# Date: 11/10/2022

workflow ArCCH_PoN2 {
    input {
        Array[File] aligned_bam
        Array[File] aligned_bai
        Array[String] sample_name
        File interval_list

        File reference
        File reference_fai
        File reference_dict
    }

    Array[Pair[String, Pair[String, String]]] aligned_bam_bai = zip(sample_name, zip(aligned_bam, aligned_bai))

    # Some of our callers use BED file instead of interval list
    call intervalsToBed as interval_to_bed {
        input: interval_list = interval_list
    }

    scatter (bam in aligned_bam_bai) {
        # Mutect
        call mutectTumorOnly as mutectTask {
          input:
          reference = reference,
          reference_fai = reference_fai,
          reference_dict = reference_dict,
          tumor_bam = bam.right.left,
          tumor_bam_bai = bam.right.right,
          interval_list = interval_list
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as mutectSanitizeVcf {
            input: vcf = mutectTask.vcf
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as mutectNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = mutectSanitizeVcf.sanitized_vcf,
            vcf_tbi = mutectSanitizeVcf.sanitized_vcf_tbi
        }

        call bcftoolsFilter as mutectFilter {
            input:
            vcf = mutectNormalize.normalized_vcf,
            vcf_tbi = mutectNormalize.normalized_vcf_tbi
        }

        call vardictTumorOnly as vardictTask {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = bam.right.left,
            tumor_bam_bai = bam.right.right,
            interval_bed = interval_to_bed.interval_bed,
            tumor_sample_name = bam.left
        }

        # Cleans the VCF output that don't match the expected VCF Format
        call vcfSanitize as vardictSanitizeVcf {
            input: vcf = vardictTask.vcf
        }

        # Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        call bcftoolsNorm as vardictNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = vardictSanitizeVcf.sanitized_vcf,
            vcf_tbi = vardictSanitizeVcf.sanitized_vcf_tbi
        }

        call lofreq_indelqual {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = bam.right.left,
            tumor_bam_bai = bam.right.right
        }

        call lofreqTumorOnly as lofreqTask {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = lofreq_indelqual.output_indel_qual_bam,
            tumor_bam_bai = lofreq_indelqual.output_indel_qual_bai,
            interval_bed = interval_to_bed.interval_bed
        }

        call lofreqReformat as reformat {
            input:
            vcf = lofreqTask.vcf,
            tumor_sample_name = bam.left
        }

        call vcfSanitize as lofreqSanitizeVcf {
            input: vcf = reformat.reformat_vcf
        }

        call bcftoolsNorm as lofreqNormalize {
            input:
            reference = reference,
            reference_fai = reference_fai,
            vcf = lofreqSanitizeVcf.sanitized_vcf,
            vcf_tbi = lofreqSanitizeVcf.sanitized_vcf_tbi
        }

        call bcftoolsFilter as lofreqFilter {
            input:
            vcf = lofreqNormalize.normalized_vcf,
            vcf_tbi = lofreqNormalize.normalized_vcf_tbi
        }
    }

    call bcftoolsMerge as merge_mutect_pon2 {
        input:
            vcfs = mutectFilter.filtered_vcf,
            vcf_tbis = mutectFilter.filtered_vcf_tbi,
            merged_vcf_basename = "mutect_pon2"
    }

    call bcftoolsMerge as merge_vardict_pon2 {
        input:
            vcfs = vardictNormalize.normalized_vcf,
            vcf_tbis = vardictNormalize.normalized_vcf_tbi,
            merged_vcf_basename = "vardict_pon2"
    }

    call bcftoolsMerge as merge_lofreq_pon2 {
        input:
            vcfs = lofreqFilter.filtered_vcf,
            vcf_tbis = lofreqFilter.filtered_vcf_tbi,
            merged_vcf_basename = "lofreq_pon2"
    }

    call bcftoolsPoN2 as mutect_PoN2 {
        input:
            vcf = merge_mutect_pon2.merged_vcf,
            caller = "mutect"
    }

    call bcftoolsPoN2 as vardict_PoN2 {
        input:
            vcf = merge_vardict_pon2.merged_vcf,
            caller = "vardict"
    }

    call bcftoolsPoN2 as lofreq_PoN2 {
        input:
            vcf = merge_lofreq_pon2.merged_vcf,
            caller = "lofreq"
    }


    output {
        File mutect_pon2_file = mutect_PoN2.pon2_vcf
        File mutect_pon2_file_tbi = mutect_PoN2.pon2_vcf_tbi
        File lofreq_pon2_file = lofreq_PoN2.pon2_vcf
        File lofreq_pon2_file_tbi = lofreq_PoN2.pon2_vcf_tbi
        File vardict_pon2_file = vardict_PoN2.pon2_vcf
        File vardict_pon2_file_tbi = vardict_PoN2.pon2_vcf_tbi
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

task mutectTumorOnly {
    input {
        File reference
        File reference_fai
        File reference_dict
        File tumor_bam
        File tumor_bam_bai
        File interval_list
    }

    Float reference_size = size([reference, reference_fai, reference_dict], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_list, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: "32GB"
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
            --read-index ~{tumor_bam_bai} \
            --f1r2-tar-gz mutect.f1r2.tar.gz \
            --max-reads-per-alignment-start 0

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
    }

    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: "6GB"
        cpu: cores
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
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB") * 2 + size([reference, reference_fai], "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: "6GB"
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

task bcftoolsFilter {
    input {
        File vcf
        File vcf_tbi
    }

    Int space_needed_gb = 10 + round(size([vcf, vcf_tbi], "GB") * 2)
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        memory: "6GB"
        docker: "kboltonlab/bst"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
    /usr/local/bin/bcftools filter -i 'FILTER~"PASS" && FMT/AF >= 0.02' -Oz -o bcftools_filtered.vcf.gz ~{vcf}
    /usr/local/bin/tabix bcftools_filtered.vcf.gz
    >>>

    output {
        File filtered_vcf = "bcftools_filtered.vcf.gz"
        File filtered_vcf_tbi = "bcftools_filtered.vcf.gz.tbi"
    }
}

task vardictTumorOnly {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name = "TUMOR"
        File interval_bed
        Int? mem_limit_override
        Int? cpu_override
        Int? JavaXmx = 24
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 4*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0
    Int memory = select_first([mem_limit_override, ceil(bam_size/6 + 5)]) # We want the base to be around 6
    Int memory_total = cores * memory

    runtime {
        docker: "kboltonlab/vardictjava:bedtools"
        memory: memory_total + "GB"
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

        how_many_lines=$(wc -l ~{interval_bed} | cut -d' ' -f1)
        if [[ how_many_lines -lt 25 ]]; then
            bedtools makewindows -b ~{interval_bed} -w 50150 -s 50000 > ~{basename(interval_bed, ".bed")}_windows.bed
        else
            cp ~{interval_bed} ~{basename(interval_bed, ".bed")}_windows.bed
        fi

        /opt/VarDictJava/build/install/VarDict/bin/VarDict \
            -U -G ~{reference} \
            -X 1 \
            -f 0.02 \
            -N ~{tumor_sample_name} \
            -b ~{tumor_bam} \
            -c 1 -S 2 -E 3 -g 4 ~{basename(interval_bed, ".bed")}_windows.bed \
            -th ~{cores} | \
        /opt/VarDictJava/build/install/VarDict/bin/teststrandbias.R | \
        /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
            -N "~{tumor_sample_name}" \
            -E \
            -f 0.02 > ~{output_name}.vcf

        /usr/bin/bgzip ~{output_name}.vcf && /usr/bin/tabix ~{output_name}.vcf.gz
    >>>

    output {
        File vcf = "~{output_name}.vcf.gz"
        File vcf_tbi = "~{output_name}.vcf.gz.tbi"
    }
}

task lofreq_indelqual {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
    }

    Int cores = 1
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 4*bam_size)
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/lofreq:latest"
        memory: "6GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o output.indel.bam ~{tumor_bam}
        samtools index output.indel.bam
    >>>

    output {
        File output_indel_qual_bam = "output.indel.bam"
        File output_indel_qual_bai = "output.indel.bam.bai"
    }
}

task lofreqTumorOnly {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        File interval_bed
        String? output_name = "lofreq.vcf"
    }

    Int cores = 4
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = 10 + round(reference_size + 2*bam_size + size(interval_bed, "GB"))
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/lofreq:latest"
        memory: "24GB"
        cpu: cores
        bootDiskSizeGb: space_needed_gb
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /opt/lofreq/bin/lofreq call-parallel --pp-threads ~{cores} -A -B -f ~{reference} --call-indels --bed ~{interval_bed} -o ~{output_name} ~{tumor_bam} --force-overwrite
        bgzip ~{output_name} && tabix ~{output_name}.gz
    >>>

    output {
        File vcf = "~{output_name}.gz"
        File vcf_tbi = "~{output_name}.gz.tbi"
    }
}

task lofreqReformat {
    input {
        File vcf
        String tumor_sample_name
    }

    Int space_needed_gb = 10 + 2*round(size(vcf, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
      cpu: cores
      docker: "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
      memory: "6GB"
      bootDiskSizeGb: space_needed_gb
      disks: "local-disk ~{space_needed_gb} SSD"
      preemptible: preemptible
      maxRetries: maxRetries
    }

    command <<<
        zcat ~{vcf} | grep "##" > lofreq.reformat.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">"  >> lofreq.reformat.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> lofreq.reformat.vcf;
        zcat ~{vcf} | grep -v '#' | awk '{ n=split($8, semi, /;/); sample=""; format=""; for(i in semi){ split(semi[i], equ, /=/); if(i<=3){ if(i+1==4) sample=sample equ[2]; else sample=sample equ[2] ":"; if(i+1==4) format=format equ[1]; else format=format equ[1] ":";}}{print $0, "GT:"format, "0/1:"sample}}' OFS='\t' >> lofreq.reformat.vcf;
        bgzip lofreq.reformat.vcf && tabix lofreq.reformat.vcf.gz
    >>>

    output {
        File reformat_vcf = "lofreq.reformat.vcf.gz"
        File reformat_vcf_tbi = "lofreq.reformat.vcf.gz.tbi"
    }
}

task bcftoolsMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Int space_needed_gb = 10 + round(2*(size(vcfs, "GB") + size(vcf_tbis, "GB")))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: "6GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_file = merged_vcf_basename + ".vcf.gz"

    command <<<
        /usr/local/bin/bcftools merge -m none --output-type z -o ~{output_file} --threads ~{cores} ~{sep=" " vcfs}
        /usr/local/bin/tabix ~{output_file}
    >>>

    output {
        File merged_vcf = output_file
        File merged_vcf_tbi = "~{output_file}.tbi"
    }
}

task bcftoolsPoN2 {
    input {
        File vcf
        String caller
    }

    Int space_needed_gb = 10 + round(2*size(vcf, "GB"))
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: "6GB"
        disks: "local-disk ~{space_needed_gb} SSD"
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        if [[ "~{caller}" =~ "mutect" ]];then
            /usr/local/bin/bcftools +fill-tags -- ~{vcf} -t NS > NS.vcf
            /usr/local/bin/bcftools filter -i 'INFO/NS >= 2' NS.vcf > 2N.vcf
            /usr/local/bin/bcftools +fill-tags -- 2N.vcf -t 'max_VAF=max(AF)' > ~{caller}.2N.maxVAF.vcf
            bgzip ~{caller}.2N.maxVAF.vcf
            tabix ~{caller}.2N.maxVAF.vcf.gz
        else
            /usr/local/bin/bcftools +fill-tags -- ~{vcf} -t NS > NS.vcf
            /usr/local/bin/bcftools filter -i 'INFO/NS >= 2' NS.vcf > 2N.vcf
            /usr/local/bin/bcftools +fill-tags -- 2N.vcf -t 'max_VAF=max(FORMAT/AF)' > ~{caller}.2N.maxVAF.vcf
            bgzip 2N.maxVAF.vcf
            tabix 2N.maxVAF.vcf.gz
        fi
    >>>

    output {
        File pon2_vcf = "~{caller}.2N.maxVAF.vcf.gz"
        File pon2_vcf_tbi = "~{caller}.2N.maxVAF.vcf.gz.tbi"
    }
}
