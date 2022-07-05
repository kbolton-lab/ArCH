version 1.0

# WDL Pipeline to merge VCFs
# -------------------------------
# This pipeline takes a list of VCFs and merges them into a single VCF.
# The bcftools merge workflow is the same used in the terra_pipeline. 

# Created by: Wendy Wong
# Date: 06/29/2022

workflow merge_vcfs {
    input {
        
		#list of VCFs to merge
		Array[File] vcf_files
        Array[File] vcf_files_tbi

        # Reference
        File reference
        File reference_fai

		#output vcf basename
		String? merged_vcf_basename = "merged"
        
	}

	Array[Pair[File,File]] vcf_pairs = zip(vcf_files, vcf_files_tbi)

	# Normalize the VCF by left aligning and trimming indels. Also splits multiallelics into multiple rows
        scatter (vcf_pair in vcf_pairs) {
		
			call bcftoolsNorm  {
            	input:
            	reference = reference,
            	reference_fai = reference_fai,
            	vcf = vcf_pair.left,
            	vcf_tbi = vcf_pair.right
        	}
		}	


        call bcftoolsMerge {
            input:
            vcfs = bcftoolsNorm.normalized_vcf,
            vcf_tbis = bcftoolsNorm.normalized_vcf_tbi,
            merged_vcf_basename = merged_vcf_basename
        }

		output {
			# merged VCF
			File merged_vcf = bcftoolsMerge.merged_vcf
			File merged_vcf_tbi = bcftoolsMerge.merged_vcf_tbi
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
        /usr/local/bin/bcftools norm --multiallelics -any --output-type z --output bcftools_norm.vcf.gz ~{vcf} -f ~{reference}
        /usr/local/bin/tabix bcftools_norm.vcf.gz
    >>>

    output {
        File normalized_vcf = "bcftools_norm.vcf.gz"
        File normalized_vcf_tbi = "bcftools_norm.vcf.gz.tbi"
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
        /usr/local/bin/bcftools merge --output-type z -o ~{output_file} ~{sep=" " vcfs}
        /usr/local/bin/tabix ~{output_file}
    >>>

    output {
        File merged_vcf = output_file
        File merged_vcf_tbi = "~{output_file}.tbi"
    }
}