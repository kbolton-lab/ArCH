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

    # In order to parallelize as much as the workflow as possible, we analyze by chromosome
    call splitBedToChr {
        input:
            interval_bed = interval_to_bed.interval_bed
    }

    scatter (bam in aligned_bam_bai) {

        scatter (bed_chr in splitBedToChr.split_chr) {
            # Mutect
            call mutect {
              input:
              reference = reference,
              reference_fai = reference_fai,
              reference_dict = reference_dict,
              tumor_bam = bam.right.left,
              tumor_bam_bai = bam.right.right,
              interval_list = bed_chr
            }

            call sanitizeNormalizeFilter as mutectFilter {
                input:
                    vcf = mutect.vcf,
                    vcf_tbi = mutect.vcf_tbi,
                    reference = reference,
                    reference_fai = reference_fai
            }

            call vardict {
                input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = bam.right.left,
                tumor_bam_bai = bam.right.right,
                interval_bed = bed_chr,
                tumor_sample_name = bam.left
            }

            call sanitizeNormalizeFilter as vardictFilter {
                input:
                    vcf = vardict.vcf,
                    vcf_tbi = vardict.vcf_tbi,
                    reference = reference,
                    reference_fai = reference_fai
            }

            call lofreq {
                input:
                reference = reference,
                reference_fai = reference_fai,
                tumor_bam = bam.right.left,
                tumor_bam_bai = bam.right.right,
                tumor_sample_name = bam.left,
                interval_bed = bed_chr
            }

            call sanitizeNormalizeFilter as lofreqFilter {
                input:
                    vcf = lofreq.vcf,
                    vcf_tbi = lofreq.vcf_tbi,
                    reference = reference,
                    reference_fai = reference_fai
            }
        }

        call bcftoolsConcat as merge_mutect {
            input:
                vcfs = mutectFilter.filtered_vcf,
                vcf_tbis = mutectFilter.filtered_vcf_tbi,
                merged_vcf_basename = "mutect." + bam.left
        }

        call bcftoolsConcat as merge_vardict {
            input:
                vcfs = vardictFilter.filtered_vcf,
                vcf_tbis = vardictFilter.filtered_vcf_tbi,
                merged_vcf_basename = "vardict." + bam.left
        }

        call bcftoolsConcat as merge_lofreq {
            input:
                vcfs = lofreqFilter.filtered_vcf,
                vcf_tbis = lofreqFilter.filtered_vcf_tbi,
                merged_vcf_basename = "lofreq." + bam.left
        }
    }

    call bcftoolsMerge as merge_mutect_pon2 {
        input:
            vcfs = merge_mutect.merged_vcf,
            vcf_tbis = merge_mutect.merged_vcf_tbi,
            merged_vcf_basename = "mutect_pon2"
    }

    call bcftoolsMerge as merge_vardict_pon2 {
        input:
            vcfs = merge_vardict.merged_vcf,
            vcf_tbis = merge_vardict.merged_vcf_tbi,
            merged_vcf_basename = "vardict_pon2"
    }

    call bcftoolsMerge as merge_lofreq_pon2 {
        input:
            vcfs = merge_lofreq.merged_vcf,
            vcf_tbis = merge_lofreq.merged_vcf_tbi,
            merged_vcf_basename = "lofreq_pon2"
    }

    call bcftoolsPoN2 as mutect_PoN2 {
        input:
            vcf = merge_mutect_pon2.merged_vcf,
            vcf_tbi = merge_mutect_pon2.merged_vcf_tbi,
            caller = "mutect"
    }

    call bcftoolsPoN2 as vardict_PoN2 {
        input:
            vcf = merge_vardict_pon2.merged_vcf,
            vcf_tbi = merge_vardict_pon2.merged_vcf_tbi,
            caller = "vardict"
    }

    call bcftoolsPoN2 as lofreq_PoN2 {
        input:
            vcf = merge_lofreq_pon2.merged_vcf,
            vcf_tbi = merge_lofreq_pon2.merged_vcf_tbi,
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
    }

    Int preemptible = 1
    Int maxRetries = 0
    Float data_size = size([vcf, vcf_tbi], "GB")
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

        # Santize Step
        # 1) removes lines containing non ACTGN bases, as they conflict with the VCF spec and cause GATK to choke
        # 2) removes mutect-specific format tags containing underscores, which are likewise illegal in the vcf spec
        gunzip -c ~{vcf} | perl -a -F'\t' -ne 'print $_ if $_ =~ /^#/ || $F[3] !~ /[^ACTGNactgn]/' | sed -e "s/ALT_F1R2/ALTF1R2/g;s/ALT_F2R1/ALTF2R1/g;s/REF_F1R2/REFF1R2/g;s/REF_F2R1/REFF2R1/g" > sanitized.vcf
        bgzip sanitized.vcf && tabix sanitized.vcf.gz

        # Normalize Step
        bcftools norm --check-ref w --multiallelics -any --output-type z --output norm.vcf.gz sanitized.vcf.gz -f ~{reference}
        tabix norm.vcf.gz

        bcftools filter -i 'FILTER~"PASS" && FMT/AF >= 0.02' -Oz -o filtered.vcf.gz norm.vcf.gz
        tabix filtered.vcf.gz
    >>>

    output {
        File filtered_vcf = "filtered.vcf.gz"
        File filtered_vcf_tbi = "filtered.vcf.gz.tbi"
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
    }

    Float reference_size = size([reference, reference_fai, interval_bed], "GB")
    Float data_size = size([tumor_bam, tumor_bam_bai], "GB")
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

        bedtools makewindows -b ~{interval_bed} -w 20250 -s 20000 > ~{basename(interval_bed, ".bed")}_windows.bed

        # Split bed file into 16 equal parts
        split -d --additional-suffix .bed -n l/16 ~{basename(interval_bed, ".bed")}_windows.bed splitBed.

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
                -b ~{tumor_bam} \
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

task lofreq {
    input {
        File reference
        File reference_fai
        File tumor_bam
        File tumor_bam_bai
        String tumor_sample_name
        File interval_bed
        String? output_name = "lofreq"
        Int? mem_limit_override
        Int? cpu_override
    }

    Int memory = select_first([mem_limit_override, 4])
    Int cores = select_first([cpu_override, 1])
    Float reference_size = size([reference, reference_fai], "GB")
    Float bam_size = size([tumor_bam, tumor_bam_bai], "GB")
    Int space_needed_gb = ceil(3 * bam_size + reference_size)
    Int preemptible = 3
    Int maxRetries = 3

    runtime {
        docker: "kboltonlab/lofreq:latest"
        memory: cores * memory + "GB"
        cpu: cores
        bootDiskSizeGb: 10
        disks: "local-disk ~{space_needed_gb} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        # Lofreq
        /opt/lofreq/bin/lofreq indelqual --dindel -f ~{reference} -o lofreq.indel.bam ~{tumor_bam}
        samtools index lofreq.indel.bam
        /opt/lofreq/bin/lofreq call-parallel --pp-threads ~{cores} -A -B -f ~{reference} --call-indels --bed ~{interval_bed} -o unsorted.lofreq.vcf lofreq.indel.bam
        cat unsorted.lofreq.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > lofreq.vcf

        # Reformat
        grep "##" lofreq.vcf > lofreq.reformat.vcf
        echo -e "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw Depth\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">"  >> lofreq.reformat.vcf;
        echo -e "##FORMAT=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias at this position\">"  >> lofreq.reformat.vcf;
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{tumor_sample_name}" >> lofreq.reformat.vcf;
        grep -v '#' lofreq.vcf | \
        awk '{
            n=split($8, semi, /;/); 
            sample=""; format=""; 
            for(i in semi){
                split(semi[i], equ, /=/);
                if(equ[1]=="DP" || equ[1]=="AF" || equ[1]=="SB"){
                    format=format equ[1] ":";
                    sample=sample equ[2] ":";
                }
            }
            format=substr(format, 1, length(format)-1);
            sample=substr(sample, 1, length(sample)-1);
            {print $0, "GT:"format, "0/1:"sample}
        }' OFS='\t' >> lofreq.reformat.vcf;
        bgzip lofreq.reformat.vcf && tabix lofreq.reformat.vcf.gz
    >>>

    output {
        File vcf = "lofreq.reformat.vcf.gz"
        File vcf_tbi = "lofreq.reformat.vcf.gz.tbi"
    }
}

task bcftoolsMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Int space_needed_gb = 10 + round((size(vcfs, "GB") + size(vcf_tbis, "GB")))
    Int memory = 6
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 2

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
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

task bcftoolsConcat {
    input {
        Array[File] vcfs
        Array[File] vcf_tbis
        String merged_vcf_basename = "merged"
    }

    Float data_size = size(vcfs, "GB")
    Int space_needed_gb = ceil(10 + data_size)
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
        /usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o ~{output_file} ~{sep=" " vcfs}
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
        File vcf_tbi
        String caller
    }

    Int space_needed_gb = 10 + round(size(vcf, "GB"))
    Int memory = 12
    Int cores = 1
    Int preemptible = 1
    Int maxRetries = 2

    runtime {
        docker: "kboltonlab/bst:latest"
        memory: cores * memory + "GB"
        disks: "local-disk ~{space_needed_gb} HDD"
        cpu: cores
        bootDiskSizeGb: 10
        preemptible: preemptible
        maxRetries: maxRetries
    }

    command <<<
        /usr/local/bin/bcftools +fill-tags -Oz -o NS.vcf.gz -- ~{vcf} -t NS
        /usr/local/bin/bcftools filter -i 'INFO/NS >= 2' -Oz -o 2N.vcf.gz NS.vcf.gz
        if [[ "~{caller}" =~ "mutect" ]];then
            /usr/local/bin/bcftools +fill-tags -Oz -o ~{caller}.2N.maxVAF.vcf.gz -- 2N.vcf.gz -t 'max_VAF=max(AF)'
            tabix ~{caller}.2N.maxVAF.vcf.gz
        else
            /usr/local/bin/bcftools +fill-tags -Oz -o ~{caller}.2N.maxVAF.vcf.gz -- 2N.vcf.gz -t 'max_VAF=max(FORMAT/AF)'
            tabix ~{caller}.2N.maxVAF.vcf.gz
        fi
    >>>

    output {
        File pon2_vcf = "~{caller}.2N.maxVAF.vcf.gz"
        File pon2_vcf_tbi = "~{caller}.2N.maxVAF.vcf.gz.tbi"
    }
}