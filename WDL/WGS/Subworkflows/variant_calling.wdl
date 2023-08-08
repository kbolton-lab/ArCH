version 1.0

workflow variant_calling {
    input {
        File aligned_bam_file
        File aligned_bai_file
        String tumor_sample_name
        File target_intervals               # Interval List

        # Reference
        File reference
        File reference_fai
        File reference_dict

        Float af_threshold = 0.0001                # Minimum VAF Cut-Off

        # See: https://github.com/bcbio/bcbio_validations/blob/master/somatic-lowfreq/README.md
        # Parameters MQ, NM, DP, and QUAL are calculated using a small subset then identifying the cut-off for 2% of the left side samples
        String bcbio_filter_string = "((FMT/AF * FMT/DP < 3) && ((INFO/MQ < 55.0 && INFO/NM > 1.0) || (INFO/MQ < 60.0 && INFO/NM > 2.0) || (FMT/DP < 10) || (INFO/QUAL < 30)))"

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

    call bamIndex {
        input:
            bam_file = aligned_bam_file
    }

    scatter (bed_chr in splitBedToChr.split_chr) {
        # Mutect
        call mutect {
            input:
            reference = reference,
            reference_fai = reference_fai,
            reference_dict = reference_dict,
            gnomad = normalized_gnomad_exclude,
            gnomad_tbi = normalized_gnomad_exclude_tbi,
            tumor_bam = aligned_bam_file, 
            tumor_bam_bai = bamIndex.bam_index, 
            interval_list = bed_chr
        }

        call mutect_pass {
            input:
            mutect_vcf = mutect.vcf,
            mutect_vcf_tbi = mutect.vcf_tbi
        }

        # Vardict
        call vardict {
            input:
            reference = reference,
            reference_fai = reference_fai,
            tumor_bam = aligned_bam_file, 
            tumor_bam_bai = bamIndex.bam_index,
            interval_bed = bed_chr,
            min_var_freq = af_threshold,
            tumor_sample_name = tumor_sample_name,
            mutect_vcf = mutect_pass.vcf
        }
    }

    call mergeVcf as merge_mutect {
        input:
            vcfs = mutect.vcf,
            vcf_tbis = mutect.vcf_tbi,
            merged_vcf_basename = "mutect." + tumor_sample_name
    }

    call mergeVcf as merge_vardict {
        input:
            vcfs = vardict.vcf,
            vcf_tbis = vardict.vcf_tbi,
            merged_vcf_basename = "vardict." + tumor_sample_name
    }

    call sanitizeNormalizeFilter as mutect_filter {
        input:
            vcf = merge_mutect.merged_vcf,
            vcf_tbi = merge_mutect.merged_vcf_tbi,
            reference = reference,
            reference_fai = reference_fai,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            vcf2PON = mutect_pon2_file,
            vcf2PON_tbi = mutect_pon2_file_tbi,
            caller = "mutect",
            sample_name = tumor_sample_name,
            filter_string = bcbio_filter_string
    }

    call sanitizeNormalizeFilter as vardict_filter {
        input:
            vcf = merge_vardict.merged_vcf,
            vcf_tbi = merge_vardict.merged_vcf_tbi,
            reference = reference,
            reference_fai = reference_fai,
            exclude_vcf = normalized_gnomad_exclude,
            exclude_vcf_tbi = normalized_gnomad_exclude_tbi,
            vcf2PON = vardict_pon2_file,
            vcf2PON_tbi = vardict_pon2_file_tbi,
            caller = "vardict",
            sample_name = tumor_sample_name,
            filter_string = bcbio_filter_string
    }

    output {
        File mutect_vcf = mutect_filter.annotated_vcf
        File vardict_vcf = vardict_filter.annotated_vcf
    }
}

task intervalsToBed {
    input {
        File interval_list
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size(interval_list, "GB")
    Int space_needed_gb = ceil(data_size)
    Float memory = 1
    Int cores = 1

    runtime {
        docker: "ubuntu:bionic"
        memory: cores * memory + "GB"
        cpu: cores
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
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

task splitBedToChr {
    input {
        File interval_bed
    }

    Float data_size = size(interval_bed, "GB")
    Int space_needed_gb = ceil(data_size)
    Int memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "ubuntu:bionic"
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
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
        String sample_name
        File reference
        File reference_fai
        File reference_dict
    }

    Float data_size = size([interval_bed, bam_file, bai_file], "GB")
    Int space_needed_gb = ceil(6 * data_size)
    Int memory = 2
    Int cores = 4
    Int preemptible = 2
    Int maxRetries = 2

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst:latest"
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bam_link = sub(basename(bam_file), basename(basename(bam_file, ".bam"), ".cram"), sample_name)

    command <<<
        ln -s ~{bam_file} ~{bam_link}
        if [[ ~{bam_link} == *.cram ]]; then
            /usr/local/bin/samtools index ~{bam_link}
        fi
        intervals=$(awk '{print $1}' ~{interval_bed} | uniq)
        for chr in ${intervals}; do
            samtools view -@ ~{cores} --fast -b -T ~{reference} -o ~{sample_name}_${chr}.bam ~{bam_link} $chr
            samtools index ~{sample_name}_${chr}.bam
        done
    >>>

    output {
        Array[Pair[File, File]] split_bam_chr = zip(glob(sample_name+"_*.bam"), glob(sample_name+"_*.bam.bai"))
    }
}

task bamIndex {
        input {
        File bam_file
    }

    Float data_size = size([bam_file], "GB")
    Int space_needed_gb = ceil(10 + data_size)
    Int memory = 4
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst:latest"
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String bam_link = basename(bam_file)

    command <<<
        set -e -x -o pipefail

        ln -s ~{bam_file} ~{bam_link}
        samtools index ~{bam_link}
    >>>

    output {
        File bam_index = basename(bam_file)+".crai"
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
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, reference_dict, interval_list], "GB")
    Float data_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = ceil(10 + data_size + reference_size)
    Int memory = select_first([mem_limit_override, 4])
    Int cores = select_first([cpu_override, 1])
    Int preemptible = 3
    Int maxRetries = 3

    runtime {
        cpu: cores
        docker: "broadinstitute/gatk:4.2.0.0"
        memory: cores * memory + "GB"
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_vcf = "mutect.filtered.vcf.gz"

    command <<<
        set -e -x -o pipefail

        /gatk/gatk Mutect2 --java-options "-Xmx~{memory}g" \
        --native-pair-hmm-threads ~{cores} \
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

task sanitizeNormalizeFilter {
    input {
        File vcf
        File vcf_tbi
        File reference
        File reference_fai

        #gnomAD
        File exclude_vcf
        File exclude_vcf_tbi

        #PoN2
        File vcf2PON
        File vcf2PON_tbi
        String caller = "caller"
        String sample_name = "tumor"

        #bcbio
        String filter_string
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([vcf, vcf_tbi, exclude_vcf, exclude_vcf_tbi, vcf2PON, vcf2PON_tbi], "GB")
    Float reference_size = size([reference, reference_fai], "GB")
    Int space_needed_gb = ceil(10 + data_size + reference_size)
    Int memory = 1
    Int cores = 1

    runtime {
        memory: cores * memory + "GB"
        docker: "kboltonlab/bst"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        cpu: cores
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        set -eou pipefail

        # BCBIO First for Vardict
        if [[ ~{caller} == "vardict" ]]; then
            bcftools filter -e "~{filter_string}" ~{vcf} --output-type z --output bcbio.vcf.gz -s "BCBIO" -m+
        else
            cp ~{vcf} bcbio.vcf.gz
        fi

        # Santize Step
        # 1) removes lines containing non ACTGN bases, as they conflict with the VCF spec and cause GATK to choke
        # 2) removes mutect-specific format tags containing underscores, which are likewise illegal in the vcf spec
        gunzip -c bcbio.vcf.gz | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" > sanitized.vcf
        bgzip sanitized.vcf && tabix sanitized.vcf.gz
        
        # Normalize Step
        bcftools norm --check-ref w --multiallelics -any --output-type z --output norm.vcf.gz sanitized.vcf.gz -f ~{reference}
        tabix norm.vcf.gz

        # gnomAD Intersection
        bcftools isec -C -w1 norm.vcf.gz ~{exclude_vcf} --output-type z --output isec.vcf.gz

        # PoN2
        export name=~{caller}.~{sample_name}.vcf.gz
        printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
        printf "##INFO=<ID=PON_NAT2_percent,Number=1,Type=Integer,Description=\"Number of samples with variant at >=2 percent\">\n" >> pon2.header;
        printf "##INFO=<ID=PON_MAX_VAF,Number=1,Type=Float,Description=\"The maximum VAF found in the PoN Samples\">\n" >> pon2.header;
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\t%INFO/NS\t%INFO/max_VAF\n" ~{vcf2PON} > normal2.txt
        bgzip -f normal2.txt
        tabix -f -s1 -b2 -e2 normal2.txt.gz
        bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent,PON_NAT2_percent,PON_MAX_VAF isec.vcf.gz -Oz -o $name
        tabix $name
    >>>

    output {
        File annotated_vcf = "~{caller}.~{sample_name}.vcf.gz"
        File annotated_vcf_tbi = "~{caller}.~{sample_name}.vcf.gz.tbi"
    }
}

task mutect_pass {
    input {
        File mutect_vcf
        File mutect_vcf_tbi
    }

    Float data_size = size([mutect_vcf, mutect_vcf_tbi], "GB")
    Int space_needed_gb = ceil(10 + data_size)
    Int memory = 2
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        cpu: cores
        docker: "staphb/bcftools:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        bcftools filter -i 'FILTER="PASS" || FILTER="weak_evidence" || FILTER="strand_bias" || FILTER="weak_evidence;strand_bias"' ~{mutect_vcf} > mutect_passed.vcf
        bcftools view -r 'chr20:32434638' --no-header ~{mutect_vcf} >> mutect_passed.vcf
        bcftools sort -m ~{memory}G mutect_passed.vcf -Oz -o mutect_passed.sorted.vcf.gz 
    >>>

    output {
        File vcf = "mutect_passed.sorted.vcf.gz"
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
        File mutect_vcf
        Float min_var_freq = 0.005
        Int JavaXmx = 24
        Int? mem_limit_override
        Int? cpu_override
    }

    Float reference_size = size([reference, reference_fai, interval_bed], "GB")
    Float data_size = size([tumor_bam, tumor_bam_bai, mutect_vcf], "GB")
    Int space_needed_gb = ceil(10 + 2.5 * data_size + reference_size)
    Int preemptible = 3
    Int maxRetries = 3
    Int memory = select_first([mem_limit_override, 2])
    Int cores = select_first([cpu_override, 2])

    runtime {
        docker: "kboltonlab/vardictjava:bedtools"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    String output_name = "vardict"

    command <<<
        set -e -x -o pipefail

        echo ~{space_needed_gb}

        # TODO: Account for when Mutect File is "Empty"..

        # Make windows and intersect with vcf file
        bedtools makewindows -b ~{interval_bed} -w 1150 -s 1000 > ~{basename(interval_bed, ".bed")}_windows.bed
        bedtools intersect -u -wa -a ~{basename(interval_bed, ".bed")}_windows.bed -b ~{mutect_vcf} > interval_list_mutect.bed

        # Merge small intervals
        bedtools merge -i interval_list_mutect.bed > interval_list_mutect_merged.bed
        # intersect with bed file
        samtools view -u -b -M -L interval_list_mutect_merged.bed -T ~{reference} -t ~{reference_fai} \
            -o intersected.bam -X ~{tumor_bam} ~{tumor_bam_bai}
        samtools index intersected.bam

        # Make windows
        bedtools makewindows -b interval_list_mutect_merged.bed -w 10150 -s 10000 > interval_list_mutect_merged_windows.bed

        # Split bed file into 16 equal parts
        split -d --additional-suffix .bed -n l/16 interval_list_mutect_merged_windows.bed splitBed.

        nProcs=~{cores}
        nJobs="\j"

        for fName in splitBed.*.bed; do
            # Wait until nJobs < nProcs, only start nProcs jobs at most
            echo ${nJobs@P}
            while (( ${nJobs@P} >= nProcs )); do
                wait -n
            done

            part=$(echo $fName | cut -d'.' -f2)

            echo ${fName}
            echo ${part}

            /opt/VarDictJava/build/install/VarDict/bin/VarDict \
                -U -G ~{reference} \
                -X 1 \
                -f ~{min_var_freq} \
                -N ~{tumor_sample_name} \
                -b intersected.bam \
                -c 1 -S 2 -E 3 -g 4 ${fName} \
                --deldupvar -Q 10 -F 0x700 --fisher > result.${part}.txt &
        done;
        # Wait for all running jobs to finish
        wait

        for fName in result.*.txt; do
            cat ${fName} >> resultCombine.txt
        done;

        cat resultCombine.txt | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl \
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

task mergeVcf {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Float data_size = size(vcfs, "GB")
    Float data_size_tbi = size(vcf_tbis, "GB")
    Int space_needed_gb = ceil(10 + 2 * (data_size + data_size_tbi))
    Float memory = 1
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 0

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        bootDiskSizeGb: 10
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

