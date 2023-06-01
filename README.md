# BoltonLab (Ar)tifact Filtering (C)lonal (H)ematapoiesis Variant Calling Pipeline

This pipeline is designed to process mutant/wildtype H.sapiens sequencing data from Illumina based sequencing for low VAF CH variants. It features four variant callers (Mutect2, VarDictJava, Lofreq2, and Pindel) for variant detection and performs various false positive filters and detection methods (fp_filter, PoN Fisher’s Exact Test, XGB model, etc.). This pipeline also generates VEP style annotations for all called variants as well as additional putative driver annotations generated from various database sources (TOPMed, MSK-IMPACT, COSMIC, OncoKB, etc.) <br />

## Installation

## Usage
The default input for this pipeline is Illumina FASTQ files with UMI tags within the individual reads. However, unaligned BAM files with the same format structure is allowed. If consensus sequencing was performed prior to using this pipeline, an aligned consensus BAM can also be provided. 

This Pipeline has been tested and configured to run using Cromwell-70 as well as on TERRAbio.

### Inputs
| Variable | Type | Definition |
| --- | --- | ---|
|bam_input|Boolean|Set TRUE for BAM input, FALSE for FASTQ input|
|unaligned_bam|File|Unaligned BAM Input, leave empty for no unaligned BAM input|
|aligned_bam_file|File|Aligned BAM Input, leave empty for no aligned BAM input|
|aligned_bam_file_bai|File|Aligned BAM Index, leave empty for no aligned BAM input|
|aligned|Boolean|Set TRUE if UMI consensus sequencing was done prior to pipeline, FALSE if unaligned|
|fastq_one|File|FastQ R1 Input|
|fastq_two|File|FastQ R2 Input|
|normal_bam|File|Optional matched normal sample BAM for variant calling|
|normal_bam_bai|File|Optional matched normal sample BAM index for variant calling|
|tumor_sample_name|String|Name of the tumor sample|

### Sequence Information
| Variable | Type | Definition |
| --- | --- | ---|
|platform|String|PL Tag for BAM Metadata|
|platform_unit|String|PU Tag for BAM Metadata|
|library|String|LB Tag for BAM Metadata|
|read_structure|Array[String]|https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures|
|min_reads|Array[Int]|Minimum number of reads that constitutes a "read family"|
|where_is_umi|String|Three options: "N = Name", "R = Read", or "T = Tag"|
|umi_paired|Boolean|Set TRUE if using UMI pairs, FALSE for single read UMIs|
|umi_length|Integer|ArcherDX specifies that all UMIs are a specific length, anything shorter or longer is thrown out|
|min_base_quality|Integer|During consensus building, any base with a QUAL less than this value is masked with an N|
|max_base_error_rate|Float|During consensus building, if this percent of the bases within a "read family" do not match, the base is masked with an N|
|max_read_error_rate|Float|During consensus building, if this percent of the reads within a "read family" do not match, the entire family is removed|
|max_no_call_fraction|Float|During consensus building, the maximum fraction of no-calls (N) within the read after filtering allowed|

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
|bait_intervals|File|Coordinates for regions you were able to design probes for in the reagent. Typically the reagent provider has this information available in bed format|
|per_base_intervals|Array[File]|If QC needs to be done at a per-base resolution, provide a per base interval list|
|per_target_intervals|Array[File]|If QC needs to be done at a per-target resolution, provide a per target interval list|
|summary_intervals|Array[File]|If QC needs to be done for specific intervals, provide an interval list|
|omni_vcf|File|A VCF of previously identified sites used by verifyBamId for identifying contamination|
|omni_vcf_tbi|File|The index of the OMNI VCF|
|picard_metric_accumulation_level|String|The level at which you want the quality control to accumulate metrics|
|qc_minimum_base_quality|Integer|Minimum base quality for a base to contribute coverage|
|qc_minimum_mapping_quality|Integer|Minimum mapping quality for a read to contribute coverage|
|bqsr_known_sites|Array[File]|A series of VCF denoted sites in which have known variation to avoid confusing real variation with errors|
|bqsr_known_sites_tbi|Array[File]|The index for the VCFs within bqsr_known_sites|
|bqsr_intervals|Array[String]|A list of chromosomes which to apply the base quality score recalibration|
|chrom_sizes|File|A file containing the chromosome and the total size (bp) of the chromosome|
|af_only_snp_only_vcf|File|A VCF file that contains specific SNPs sites of interest, used for Somalier. Can be from https://github.com/brentp/somalier/releases/tag/v0.2.15|

### Variant Callers
|Variable|Type|Definition|
|---|---|---|
|tumor_only|Boolean|Set TRUE if the analysis will be done using only Tumor samples, FALSE if there is a matched normal available|
|normal_bam|File|Matched normal sample BAM, leave empty for no matched normal BAM|
|normal_bam_bai|File|Matched normal sample BAM index, leave empty for no matched normal BAM|
|target_intervals|File|Interval list for the panel|
|af_threshold|Float|Minimum VAF cut-off|
|bcbio_filter_string|String|http://bcb.io/2016/04/04/vardict-filtering/|
|pindel_insert_size|Integer|Average size of the inserts in the sequencing|
|pindel_min_supporting_reads|Integer|Minimum amount of reads supporting the INDEL|
|ref_date|String|The date of the reference|
|ref_name|String|The name of the reference|

### Filtering Parameters
|Variable|Type|Definition|
|---|---|---|
|pon_bams|Array[Pair[File, File]]|The Panel of Normal BAMs and their associated index files|
|pon_pvalue|Float|Minimum Bonferroni corrected p-value for Fisher's Exact Test of the Panel of Normals|
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
|vep_ensembl_assembly|String|When using VEP cache, denote the assembly|
|vep_ensembl_version|String|When using VEP cache, denote which version of the assembly|
|vep_ensembl_species|String|When using VEP cache, specify which species for the assembly|
|vep_plugins|Array[String]|List of plugins to be used in VEP|
|vep_plugin_spliceAI_files|CustomStruct[Files]|The four SpliceAI scores obtained from Illumina: https://github.com/Illumina/SpliceAI|
|synonyms_file|File|File of chromosome synonyms|
|annotate_coding_only|Boolean|Set TRUE if VEP should return consequences that fall within the coding only regions of the transcript, FALSE if all consequences should be returned|
|vep_custom_annotations|CustomStruct|Extra files that VEP should use for additional annotations e.g. gnomAD, clinvar, etc|
|vep_pick|String|Pick one line of block of consequence data per variant rather than returning all consequences|
|everything|Boolean| Shortcut to turn on all following VEP flags: <br/>--sift b<br/>--polyphen b<br/>--ccds<br/>--hgvs<br/>--symbol<br/>--numbers<br/>--domains<br/>--regulatory<br/>--canonical<br/>--protein<br/>--biotype<br/>--uniprot<br/>--tsl<br/>--appris<br/>--gene_phenotype<br/>--af<br/>--af_1kg<br/>--af_esp<br/>--af_gnomad<br/>--max_af<br/>--pubmed<br/>--var_synonyms<br/>--variant_class<br/>--mane|

### Required Files
Most of the files can be found inside our Files directory, however, the following list of files will need to be generated based on each project's specific requirements and availability.
1. Reference Genome
2.	Gene Panel Interval List
3.	Panel of Normals 2%
4.	SpliceAI Scores
5.	Panel of Normal Aligned Consensus BAMs
6.	BCBio Filter Parameters

## Example
An example input can be found inside pipeline.json

## Output
|File|Description|
|---|---|
|sample.fpfilter.vcf.gz|Results from Varscan's FP Filter on ALL of the variants found for every variant caller in VCF Format|
|sample.pon.total.counts.vcf.gz|Results from the PoN Pileup. Contains information regarding the reference depth, alternate depth, depth for strand, etc.|
|CALLER_FILTERS.tsv.gz|Original untouched filters from the individual callers for each specific variant|
|caller.sample.final.annotated.tsv|Variant calls after being annotated with our putative driver annotation script|
|caller.sample.final.annotated.vcf.gz|Filtered variant calls that have been annotated with VEP|
|caller.sample.pileup.fisherPON.vcf.gz|Variant calls after being filtered with our PoN Fisher's Exact Test|
|caller_full.sample.vcf.gz|Base variant calls from the callers without any post variant caller filtering|
|output_sample.tsv.gz|The FINAL output file that contains all merged information from the pipeline|
|output_caller_complex_sample.tsv.gz|Individual files containing variant calls annotated as complex with their assumed support|

## Post Pipeline Steps
After the pipeline has finished running there should be several files that are generated. The final output file would be output_sample.tsv.gz. All of these files generated for all the samples should be combined together into a final file using the following command
```sh
# Grab the header from one of the output files
zcat output_sample.tsv.gz | head -n1 > final.combined.tsv

# Merge all output files together
for dir in $(ls -d */); do
  sample_name=$(basename $dir);
  zcat $dir/output_$sample_name.tsv.gz | tail -n+2 >> final.combined.tsv;
done

# If the resulting file is too large, we can pre-filter the results prior to running our post filtering script
# Find the relative index position for the "all_fp_pass_XGB" filter
head -n1 final.combined.tsv | awk -F'\t' -vs='all_fp_pass_XGB' '{for (i=1;i<=NF;i++)if($i~"^"s"$"){print i;exit;}}'

# For all of the output files, only keep the variants that had passed all our applied filters via checking the "all_fp_pass_XGB" column
# In this case, our column is index position 337
for dir in $(ls -d */); do 
  sample_name=$(basename $dir); 
  zcat $dir/output_$sample_name.tsv.gz | tail -n+2 | awk -F'\t' '{if($337=="true")print $0}' >> final.combined.FPpass.tsv; 
done
```
Now the resulting final.combined.FPpass.tsv can be read into our PostPipelineFilters.R and the resulting file will be the final combined list of all variants that have been filtered, annotated, and checked.

## Limitations | TODO
- VEP Annotations are hardcoded to be consistent so additional plugins or lack of plugins will cause issues with downstream annotation

## Contact Information
Created by: Irenaeus Chan <br />
Email: chani@wustl.edu <br />
