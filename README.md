# BoltonLab (Ar)tifact Filtering (C)lonal (H)ematapoiesis Variant Calling Pipeline
ArCH is a somatic variant calling pipeline designed to detect low variant allele fraction (VAF) clonal hematopoiesjsonsis (CH) variants. Starting from either unaligned FASTQ/BAM/CRAM files or aligned BAM/CRAM files, ArCH utilizes four variant callers (Mutect2, VarDictJava, LoFreq2, and Pindel) to detect somatic variants. These variants are then filtered using a variety of false positive filters and detection methods (false positive filters, panel of normal, etc.). The pipeline also generates VEP style annotations for all called variants as well as additional putative driver annotations generated from various database sources (TOPMed, MSK-IMPACT, COSMIC, OncoKB, etc.).

If you end up using this tool in your publication, please cite this paper:
```
Irenaeus C C Chan, Alex Panchot, Evelyn Schmidt, et al. ArCH: improving the performance of clonal hematopoiesis variant calling and interpretation, Bioinformatics, Volume 40, Issue 4, April 2024, btae121, https://doi.org/10.1093/bioinformatics/btae121
```

## Installation
This pipeline requires several files to be downloaded and configured prior to running. The following files are required for the pipeline to run:
1. Reference Genome
2. Gene Panel Interval List
3. Panel of Normals
4. gnomAD VCF
5. VEP Cache & Plugins
6. Somalier VCF
7. COSMIC VCF

Step-by-step instructions to prepare each file will be provided below.

### Panel of Normals
One of the most important pieces of this pipeline is the Panel of Normals (PoN). The PoN is a collection of BAM files from young individuals that are used to filter out false positives from the tumor samples by representing a base threshold of noise. Typically, 10 to 20 samples yields the best performance for the PoN.

Two separate filters are generated from the PoN:
1. Threshold of Noise: This is a Bonferroni corrected Fisher's Exact Test that is used to determine the threshold of noise for the PoN. This is used to filter out any variants that are found in the tumor samples at a lower frequency than the threshold of noise.
2. Potential Germline Variants: These are non-hotspot variants that are found in 2 or more of the PoN samples at a 2% VAF or higher. These variants are then used to filter out any variants found in the tumor samples as possible germline variants.

To generate the PoN, the following steps are required:
1. Run the (ArCH Alignment WDL Workflow)[https://github.com/kbolton-lab/ArCH/blob/main/WDL/ArCH_Alignment.wdl] - This will create UMI consensus aligned BAM files for all the PoN samples.
2. Check the PoN Files for potential CH hotspots and remove them from the PoN samples.
```sh
# Using this Docker: duct/getbasecount:latest
docker run -v /path/to/PoN:/mnt duct/getbasecount:latest /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $REF --bam ${sample_name}:/mnt/PoN.bam --vcf AnnotatePD_Files/bick_kelly_HGVSp5_pileup.vcf --output ${sample_name}.pileup.vcf --maq 5 --baq 5
bgzip ${sample_name}.pileup.vcf && tabix ${sample_name}.pileup.vcf.gz

# Check the resulting pileup files for potential CH hotspots
bcftools view -i 'FORMAT/VF>0.02' ${sample_name}.pileup.vcf.gz
```
NOTE: For the variant: `chr20:32434638:A:AG`. Only remove the PoN Sample if this variant is found above 5% VAF in the PoN sample.

3. Run the UMI Consensus Aligned BAM files through the [pon2_creation.wdl](https://github.com/kbolton-lab/ArCH/blob/main/WDL/pon2_creation.wdl)

This pipeline will generate 3 files:
- mutect.2N.maxVAF.vcf.gz
- lofreq.2N.maxVAF.vcf.gz
- vardict.2N.maxVAF.vcf.gz

### gnomAD resource:
```sh
for chr in {1..22} X Y; do
  bcftools view -f PASS -i 'INFO/AF>=0.005' -Ou gnomad.exomes.v4.1.sites.chr${chr}.vcf.bgz | bcftools annotate -x ^INFO/AC,INFO/AF -Ou - | bcftools norm --multiallelics -any -Oz -o gnomad.exomes.v4.1.sites.chr${chr}.AF_only.exclude_0.005.normalized.vcf.gz -  
done
bcftools concat -Oz -o gnomad.exomes.v4.1.AF_only.exclude_0.005.normalized.vcf.gz gnomad.exomes.v4.1.sites.chr*.AF_only.exclude_0.005.normalized.vcf.gz
```

### VEP Cache:
This Pipeline's annotation step has been configured to use VEP cache files, which can be downloaded from the [Ensembl FTP - Homo sapiens v109](ftp://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_merged_vep_109_GRCh38.tar.gz). The cache files should be downloaded along with all necessary plugin files and zipped into a single file for the pipeline to use.

Create the VepData which will contain the VEP v109 Cache
```sh
mkdir VEP_cache && mkdir VEP_cache/VepData;
cd VEP_cache;
curl -O ftp://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_merged_vep_109_GRCh38.tar.gz
tar xzf homo_sapiens_merged_vep_109_GRCh38.tar.gz -C VepData
rm homo_sapiens_merged_vep_109_GRCh38.tar.gz
```

Create the plugin directory that will contain all the VEP plugins used in this pipeline
```sh
mkdir plugin
curl -o plugin/CADD.pm https://github.com/Ensembl/VEP_plugins/blob/431b1516431fcb4ee6431120c749769e6516a23e/CADD.pm
curl -o plugin/REVEL.pm https://github.com/Ensembl/VEP_plugins/blob/431b1516431fcb4ee6431120c749769e6516a23e/REVEL.pm
curl -o plugin/SpliceAI.pm https://github.com/Ensembl/VEP_plugins/blob/431b1516431fcb4ee6431120c749769e6516a23e/SpliceAI.pm
curl -o plugin/pLI.pm https://github.com/Ensembl/VEP_plugins/blob/431b1516431fcb4ee6431120c749769e6516a23e/pLI.pm
curl -o plugin/Frameshift.pm https://raw.githubusercontent.com/griffithlab/pVACtools/v2.0.0/tools/pvacseq/VEP_plugins/Frameshift.pm
curl -o plugin/Wildtype.pm https://raw.githubusercontent.com/griffithlab/pVACtools/v2.0.0/tools/pvacseq/VEP_plugins/Wildtype.pm
```

Create the individual raw resources for each of the VEP Plugins
```sh
# Synonyms File
mkdir Synonyms
curl -o Synonyms/chromAlias.txt https://hgdownload.soe.ucsc.edu/hubs/GCF/000/001/405/GCF_000001405.39/GCF_000001405.39.chromAlias.txt

# CADD
mkdir CADD
curl -o CADD/whole_genome_SNVs.tsv.gz https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz
curl -o CADD/gnomad.genomes.r4.0.indel.tsv.gz https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz

# Clinvar
mkdir Clinvar
curl -o Clinvar/clinvar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
tabix Clinvar/clinvar.vcf.gz

# pLI
mkdir pLI
curl -o pLI/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt https://ftp.broadinstitute.org/pub/ExAC_release/release0.3/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt
awk '{print $2, $20}' pLI/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt > pLI/pLI_gene.txt

# REVEL
mkdir REVEL
curl -O https://www.google.com/url?q=https%3A%2F%2Frothsj06.dmz.hpc.mssm.edu%2Frevel-v1.3_all_chromosomes.zip&sa=D&sntz=1&usg=AOvVaw2DS2TWUYl__0vqijzzxp5M

# SpliceAI
mkdir spliceAI
# Files for spliceAI can be downloaded from: https://basespace.illumina.com/s/otSPW8hnhaZR
# You will need:
# - spliceai_scores.raw.indel.hg38.vcf.gz
# - spliceai_scores.raw.indel.hg38.vcf.gz.tbi
# - spliceai_scores.raw.snv.hg38.vcf.gz
# - spliceai_scores.raw.snv.hg38.vcf.gz.tbi
```

ZIP all the files into a single file for the pipeline to use
```sh
zip -r VEP_cache.zip VEP_cache
```

### Somalier
[Somalier (v0.2.15) written by Brentp](https://github.com/brentp/somalier/releases/tag/v0.2.15) is run as a part of this pipeline to ensure that the samples are correctly identified. 

The VCF file containing all the sites used for this step can be prepared as follows
```sh
curl -O https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
```

### COSMIC Counts Files
To generate the necessary file inputs for these files requires a lot of work and we have not optimized a pipeline to automatically perform this task. Please download the necessary zip file which has been configured for [COSMIC v94](https://arch-example-files.s3.us-east-2.amazonaws.com/cosmic_zip.zip)

## Usage
For ease of use, a [WDL pipeline](https://github.com/kbolton-lab/ArCH/blob/main/WDL/ArCH.wdl) is available to run the entire ArCH pipeline from either unaligned FASTQ/BAM/CRAM or aligned BAM/CRAM files.

An example sample has been provided which can be accessed through the following links:
```
https://arch-example-files.s3.us-east-2.amazonaws.com/ArCH_S1_R1.fastq.gz
https://arch-example-files.s3.us-east-2.amazonaws.com/ArCH_S1_R2.fastq.gz
```

This Pipeline has been tested and configured to run using Cromwell-70 as well as on TERRAbio.

### Inputs
For basic usage, please use the following [JSON](https://github.com/kbolton-lab/ArCH/tree/main/Example/ArCH_pipeline.json) as a base template for your input file.

| Variable | Type | Definition |
| --- | --- | --- |
|input_file|File|This will be the first input file. It can be a R1 FASTQ, BAM, or CRAM|
|input_file_two|File|This will be the second input file. It can be R2 FASTQ, BAI, or CRAI|
|tumor_sample_name|String|Name of the tumor sample|
|normal_bam|File|Optional matched normal sample BAM for variant calling|
|normal_bai|File|Optional matched normal sample BAM index for variant calling|
|normal_sample_name|String|Name of the normal sample|
|input_type|String|Three options: "BAM", "CRAM", or "FASTQ" (Default: "FASTQ")|
|aligned|Boolean|Set TRUE if UMI consensus sequencing building ArCH_Alignment.wdl was done prior to pipeline, FALSE if unaligned|
|target_intervals|File|Interval list for the sequencing panel|

### Sequence and UMI Information
| Variable | Type | Definition |
| --- | --- | --- |
|platform|String|PL Tag for BAM Metadata e.g. Illumina, PacBio, ElementBio, etc...|
|platform_unit|String|PU Tag for BAM Metadata e.g. the specific sequencing machine used|
|library|String|LB Tag for BAM Metadata e.g. ArcherDX VariantPlex, MGI, IlluminaWES|
|has_umi|Boolean|Set TRUE if the sequencing data has UMIs, FALSE if it does not|
|umi_paired|Boolean|Set TRUE if the sequencing data has paired UMIs, FALSE if it does not|
|where_is_umi|String|Three options: "N = Name", "R = Read", or "T = Tag"|
|read_structure|Array[String]|https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures|
|min_reads|Array[Int]|Minimum number of reads that constitutes a "read family" (Default: 1)|

### Consensus Building
|min_base_quality|Integer|During consensus building, any base with a QUAL less than this value is masked with an N (Default: 1)|
|max_base_error_rate|Float|During consensus building, if this percent of the bases within a "read family" do not match, the base is masked with an N (Default: 0.1)|
|max_read_error_rate|Float|During consensus building, if this percent of the reads within a "read family" do not match, the entire family is removed (Default: 0.05)|
|max_no_call_fraction|Float|During consensus building, the maximum fraction of no-calls (N) within the read after filtering allowed (Default: 0.5)|

### Reference
| Variable | Type | Definition |
| --- | --- | --- |
|reference|File|Reference|
|reference_fai|File|Reference FAI|
|reference_dict|File|Reference Dictionary|
|reference_amb|File|Reference AMB|
|reference_ann|File|Reference ANN|
|reference_bwt|File|Reference BWT|
|reference_pac|File|Reference PAC|
|reference_sa|File|Reference SA|

### Quality Control
|Variable|Type|Definition|
|---|---|---|
|apply_bqsr|Boolean|Set TRUE if Base Quality Score Recalibration should be applied, FALSE if it should not (Default: False)|
|bqsr_known_sites|Array[File]|A series of VCF denoted sites in which have known variation to avoid confusing real variation with errors. Can be downloaded from: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/|
|bqsr_known_sites_tbi|Array[File]|The index for the VCFs within bqsr_known_sites|
|af_only_snp_only_vcf|File|A VCF file that contains specific SNPs sites of interest, used for Somalier. Can be from https://github.com/brentp/somalier/releases/tag/v0.2.15|

### Variant Callers
|Variable|Type|Definition|
|---|---|---|
|tumor_only|Boolean|Set TRUE if the analysis will be done using only Tumor samples, FALSE if there is a matched normal available|
|af_threshold|Float|Minimum VAF cut-off (Default: 0.0001)|
|bcbio_filter_string|String|https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/vardict.py#L251|

### Filtering Parameters
|Variable|Type|Definition|
|---|---|---|
|pon_bams|Array[Pair[File, File]]|The Panel of Normal BAMs and their associated index files|
|pon_pvalue|Float|Minimum Bonferroni corrected p-value for Fisher's Exact Test of the Panel of Normals (Default: 2.114164905e-6)|
|normalized_gnomad_exclude|File|Filtered gnomAD VCF with VAFs higher than 0.5%
|normalized_gnomad_exclude_tbi|File|Filtered gnomAD VCF index|
|mutect_pon2_file|File|Mutect2 called variants from PoN BAMs that are found in two or more samples above 2% VAF|
|mutect_pon2_file_tbi|File|Mutect2 called variants from PoN BAMs index|
|lofreq_pon2_file|File|Lofreq2 called variants from PoN BAMs that are found in two or more samples above 2% VAF|
|lofreq_pon2_file_tbi|File|Lofreq2 called variants from PoN BAMs index|
|vardict_pon2_file|File|VarDictJava called variants from PoN BAMs that are found in two or more samples above 2% VAF|
|vardict_pon2_file_tbi|File|VarDictJava called variants from PoN BAMs index|

### Annotation Parameters
|Variable|Type|Definition|
|---|---|---|
|vep_cache_dir_zip|File|The VEP cache directory in ZIP format|
|vep_plugins|Array[String]|List of plugins to be used in VEP (Default: "Frameshift", "Wildtype")|
|synonyms_file|File|File of chromosome synonyms|
|annotate_coding_only|Boolean|Set TRUE if VEP should return consequences that fall within the coding only regions of the transcript, FALSE if all consequences should be returned|

### BCBio Filter Parameters
According to BCBIO, VarDict has multiple false positive calls at regions of low depth and allelic fractions. These are the [default](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/vardict.py#L251) parameters recommended by BCBIO. However, we have found that these parameters are too stringent for our purposes and have modified them to the following:
```
- Low mapping quality and multiple mismatches in a read (NM)
  For bwa only: MQ < 55.0 and NM > 1.0 or MQ < 60.0 and NM > 3.0
- Low depth (DP < n) where n is calculated as 0.25 of the average read depth e.g. If the average read depth is 20,000 bp, then the n is 5000
- Low QUAL (QUAL < 27)
```

## Output
|File|Description|
|---|---|
|sample.bam|Aligned BAM file|
|sample.bam.bai|BAM index file|
|sample_fastqc.html|FASTQC report for the aligned BAM file|
|sample_fastqc.zip|FASTQC zip file for the aligned BAM file|
|sample.somalier|Somalier output file|
|variant_caller.sample.vcf.gz|Base variant calls from the callers without any filtering|
|variant_caller.sample.pileup.fisherPON.fp_filter.VEP.vcf.gz|Variant calls after being annotated with the PoN, FP filters, and VEP|
|variant_caller.sample.pileup.fisherPON.filtered.fp_filter.VEP.vcf.gz|Variant calls after being annotated with the PoN, FP filters, VEP and filtered by the Bonferroni corrected p-value|
|all_callers.sample.fpfilter.vcf.gz|Results from Varscan's FP Filter on ALL of the variants found for every variant caller in VCF Format|
|sample.pon.total.counts.vcf.gz|Results from the PoN Pileup. Contains information regarding the reference depth, alternate depth, depth for strand, etc.|
|all_callers.sample.VEP_annotated.vcf.gz|
|variant_caller.sample.pileup.fisherPON.filtered.fp_filter.VEP.tsv|Variant calls after being annotated with our putative driver annotation script|
|sample.final.annotated.tsv|The FINAL output file that contains all merged information from the pipeline|

## Post Pipeline Steps
After the pipeline has finished running there should be several files that are generated. The final output file would be `${sample_name}.final.annotated.tsv`. All of these files generated for all the samples should be combined together into a final file using the following command
```sh
# Grab the header from one of the output files
cat ${sample_name}.final.annotated.tsv | head -n1 > final.combined.tsv

# Merge all output files together
for dir in $(ls -d */); do
  zcat $dir/${sample_name}.final.annotated.tsv | tail -n+2 >> final.combined.tsv;
done

# If the resulting file is too large, we can pre-filter the results prior to running our post filtering script
# Find the relative index position for the "all_fp_pass" filter
head -n1 final.combined.tsv | awk -F'\t' -vs='all_fp_pass' '{for (i=1;i<=NF;i++)if($i~"^"s"$"){print i;exit;}}'

# For all of the output files, only keep the variants that had passed all our applied filters via checking the "all_fp_pass" column
# In this case, our column is index position 170
for dir in $(ls -d */); do 
  zcat $dir/${sample_name}.final.annotated.tsv | tail -n+2 | awk -F'\t' '{if($170=="TRUE")print $0}' >> final.combined.FPpass.tsv; 
done

# Now the resulting final.combined.FPpass.tsv can be used as an input into our ArCHPostPipeline.R 
LC_ALL=C.UTF-8 Rscript --vanilla ArCHPostPipeline.R --tsv final.combined.FPpass.tsv --bolton_bick_vars AnnotatePD_Files/bick.bolton.vars3.txt --gene_list AnnotatePD_Files/oncoKB_CGC_pd_table_disparity_KB_BW.csv --cosmic AnnotatePD_Files/COSMIC.heme.myeloid.hotspot.w_truncating_counts.tsv --pd_table AnnotatePD_Files/pd_table_kbreview_bick_trunc4_oncoKB_SAFE.filtered_genes_oncoKB_CGC.tsv
```
The output from ArCHPostPipeline.R will produce three output files:
|File|Description|
|---|---|
| final.all.csv | All variants that passed all post pipeline filters with additional annotations |
| final.pass.csv | Variants that passed all post pipeline filters and have been identified as potential CH variants that are putative drivers |
| final.review.csv | Variants that passed all post pipeline filters and could potentially be CH variants that are putative drivers, but need to be manually reviewed by an expert |

## Limitations | TODO
- VEP Annotations are hardcoded to work with the v109 cache be consistent so additional plugins or lack of plugins will cause issues with downstream annotation

## Contact Information
Created by: Irenaeus Chan <br />
Email: chani@wustl.edu <br />
