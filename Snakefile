# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


report: "report/workflow.rst"


"""
Snakefile to compute a VCF file per sample based on a reference genome and paired-end Illumina DNA-Seq reads.
"""
from glob import glob
from subprocess import check_output
import pandas as pd
import os
from helpers import get_fastq_file_name


#############################
# Load pipeline configuration
#############################
configfile: "config/config.yaml"

TEMP_DIR = config["tempdir"]
RESULT_DIR = config["resultdir"]

wildcard_constraints:
  sample="[A-Za-z0-9]"
wildcard_constraints:
  unit="L[0-9]"

# get samples and units
samples = pd.read_table(config["units"], index_col="sample")
SAMPLES = list(set(samples.index.values))

units = pd.read_table(config["units"], index_col=["sample", "unit"], dtype=str)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
UNITS = units.index.get_level_values('unit').unique().tolist()

wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "L[0-9]+"

##################
# Helper functions
##################

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastq(wildcards):
     """Get fastq files of given sample-unit."""
     return units.loc[(wildcards.sample, wildcards.unit), ["fq1","fq2"]].dropna()


#def merge_bams(wildcards):
#    "collects all bams files corresponding to the same library"
#    bam_files = glob(TEMP_DIR + "mapped/" + wildcards.sample + "_" + "L[0-9]+\.sorted.dedup.bam")
#    return bam_files


#################
# Desired output
#################
QC = RESULT_DIR + "multiqc_report.html"

BAMS = expand(TEMP_DIR + "mapped/{sample}_{unit}.bam",
    sample=SAMPLES,
    unit=UNITS)

STATS = [
  expand(RESULT_DIR + "stats/{sample}.stats.txt", sample = SAMPLES),
  expand(RESULT_DIR + "stats/{sample}.bigwig", sample = SAMPLES),
  RESULT_DIR + "stats/genome_coverage_depth.tsv"
]

GLOBAL_VCF  = [
  RESULT_DIR + "freebayes_variants.vcf.gz",
  RESULT_DIR + "samtools_variants.vcf.gz"
  ]

if config["remove_workdir"]:
    rule all:
        input:
            QC,
            BAMS,
            GLOBAL_VCF,
            STATS
        message:"all done! Cleaning working directory"
        shell:
            "rm -r {TEMP_DIR}"
else: 
    rule all:
        input:
            QC,
            BAMS,
            GLOBAL_VCF,
            STATS 
        message:"All done! Keeping temporary directory"


#########################
# Call SNPs with samtools 
#########################
rule samtool_pileup:
    input:
        expand(RESULT_DIR + "mapped/{sample}.bam", sample = SAMPLES)
    output:
        RESULT_DIR + "samtools_variants.vcf.gz"
    message:
        "Calling variants with samtools"
    params:
        genome = config["refs"]["genome"]
    threads: 8
    shell:
        "bcftools mpileup -a AD,DP,SP --output-type u --fasta-ref {params.genome} {input} |"
        "bcftools call --threads {threads} --multiallelic-caller --variants-only --output-type z > {output}"  


##########################
# Call SNPs with freebayes
##########################

if len(SAMPLES) == 1: # only one sample
    rule merge_variants:
        input:
            vcf = expand(RESULT_DIR + "vcf/{sample}.vcf.gz", sample = SAMPLES),
            index = expand(RESULT_DIR + "vcf/{sample}.vcf.gz.csi", sample = SAMPLES)
        output:
            RESULT_DIR + "freebayes_variants.vcf.gz"
        message: 
            "Copying {input.vcf} to {output}"
        threads: 1
        shell:
            "mv {input.vcf} {output}"
else:
    rule merge_variants:
        input:
            vcf = expand(RESULT_DIR + "vcf/{sample}.vcf.gz", sample = SAMPLES),
            index = expand(RESULT_DIR + "vcf/{sample}.vcf.gz.csi", sample = SAMPLES)
        output:
            RESULT_DIR + "freebayes_variants.vcf.gz"
        message:
            "Merging all VCF files"
        threads: 1
        shell:
            "bcftools merge {input.vcf} -Oz --output {output}"

rule index_variants:
    input:
        RESULT_DIR + "vcf/{sample}.vcf.gz"
    output:
        RESULT_DIR + "vcf/{sample}.vcf.gz.csi"
    message:
        "Compressing and indexing {wildcards.sample} VCF file."
    threads: 4
    shell:
        "bcftools index {input}"  

rule compress_variants:
    input:
        RESULT_DIR + "vcf/{sample}.vcf"
    output:
        RESULT_DIR + "vcf/{sample}.vcf.gz"
    message:
        "Compressing and indexing {wildcards.sample} VCF file."
    threads: 4
    shell:
        "bgzip {input}"   

rule call_variants:
    input:
        ref = config["refs"]["genome"],
        bam = RESULT_DIR + "mapped/{sample}.bam"
    output:
        RESULT_DIR + "vcf/{sample}.vcf"
    message:
        "calling variants for {wildcards.sample}"
    threads: 4
    shell:
        "freebayes -f {input.ref} {input.bam} > {output}"   

#################################
# Compute read mapping statistics
#################################
rule read_depth_per_genome_position:
    input:
        expand(RESULT_DIR + "mapped/{sample}.bam", sample = SAMPLES)
    output:
        RESULT_DIR + "stats/genome_coverage_depth.tsv"
    message:
        "Computing genome coverage depth for all samples"
    params:
    shell:
        "samtools depth {input} > {output}"

rule genome_coverage_bigwig:
    input:
        RESULT_DIR + "mapped/{sample}.bam"
    output:
        RESULT_DIR + "stats/{sample}.bigwig"
    message:
        "Computing genome coverage for {wildcards.sample}"
    params:
        bin_size = config["deeptools"]["binsize"]
    shell:
        "bamCoverage -b {input} "
        "--outFileFormat 'bigwig' "
        "--binSize {params.bin_size} "
        "-o {output}"

rule samtools_stats:
    input: 
        RESULT_DIR + "mapped/{sample}.bam"
    output:
        RESULT_DIR + "stats/{sample}.stats.txt"
    message:
        "Computing statistics from {wildcards.sample} BAM file"
    shell:
        "samtools stats {input} > {output}"

##############################
# Merge BAMs from same library
##############################

rule merge_bams:
    input:
        expand(TEMP_DIR + "mapped/{{sample}}_{unit}.sorted.fixed.sorted.dedup.bam", unit=UNITS)
    output:
        bam = RESULT_DIR + "mapped/{sample}.bam",
        bai = RESULT_DIR + "mapped/{sample}.bam.bai",
    message:
        "merging all BAM files for {wildcards.sample}"
    shell:
        "samtools merge {output.bam} {input};"
        "samtools index {output.bam}"


##########################################
# Align to genome, sort and mark duplicate
##########################################
rule mark_duplicate:
    input:
        TEMP_DIR + "mapped/{sample}_{unit}.sorted.fixed.sorted.bam"
    output:
        temp(TEMP_DIR + "mapped/{sample}_{unit}.sorted.fixed.sorted.dedup.bam")
    message:
        "marking duplicates in {wildcards.sample} {wildcards.unit} bam file"
    threads: 4
    shell:
        "samtools markdup -@ {threads} {input} {output}"


rule samtools_sort_by_coordinates:
    input:
        TEMP_DIR + "mapped/{sample}_{unit}.sorted.fixed.bam"
    output:
       temp(TEMP_DIR + "mapped/{sample}_{unit}.sorted.fixed.sorted.bam")
    message:
        "sorting {wildcards.sample} {wildcards.unit} bam file by coordinate"
    threads: 4
    shell:
        "samtools sort -@ {threads} {input} > {output}"

rule samtools_fixmate:
    input:
        TEMP_DIR + "mapped/{sample}_{unit}.sorted.bam"
    output:
        temp(TEMP_DIR + "mapped/{sample}_{unit}.sorted.fixed.bam")
    message:
        "Fixing mate in {wildcards.sample} {wildcards.unit} bam file"
    threads: 4
    shell:
        "samtools fixmate -m -@ {threads} {input} {output}"


rule samtools_sort_by_qname:
    input:
        TEMP_DIR + "mapped/{sample}_{unit}.bam"
    output:
        temp(TEMP_DIR + "mapped/{sample}_{unit}.sorted.bam")
    message:
        "sorting {wildcards.sample} {wildcards.unit} bam file by read name (QNAME field)"
    threads: 4
    shell:
        "samtools sort -n -@ {threads} {input} > {output}"


rule bwa_align:
    input:
        index = [TEMP_DIR + "index/genome." + ext for ext in ["amb","ann","pac","sa","bwt"]],
        forward = TEMP_DIR + "trimmed/{sample}_{unit}_forward.fastq",
        reverse = TEMP_DIR + "trimmed/{sample}_{unit}_reverse.fastq"
    output:
        TEMP_DIR + "mapped/{sample}_{unit}.bam"
    message:"mapping {wildcards.sample} {wildcards.unit} reads to genomic reference"
    params:
        db_prefix = TEMP_DIR + "index/genome"
    threads: 4
    run:
        # Building the read group id (sequencer_id + flowcell_name + lane_number + barcode)
        SEQUENCER_ID=check_output("head -n 1 " + input.forward + " |cut -d: -f1",shell=True).decode().strip()
        FLOWCELL_NAME=check_output("head -n 1 " + input.forward + " |cut -d: -f3",shell=True).decode().strip()
        FLOWCELL_LANE=check_output("head -n 1 " + input.forward + " |cut -d: -f4",shell=True).decode().strip()
        BARCODE=check_output("head -n 1 " + input.forward + " |cut -d' ' -f2 |cut -d: -f4",shell=True).decode().strip()
        # Feeding the READ_GROUP_ID to bwa
        READ_GROUP = SEQUENCER_ID + "." + FLOWCELL_NAME + "." + FLOWCELL_LANE + "." + BARCODE
        # If sample is single end, feeding only one fastq file (other outputs an empty BAM file)
        if is_single_end(wildcards.sample, wildcards.unit):
            shell("bwa mem -v 1 -t {threads} -R '@RG\\tID:{READ_GROUP}\\tPL:ILLUMINA\\tLB:{wildcards.unit}\\tSM:{wildcards.sample}' {params.db_prefix} {input.forward} >{output}")
        else:
            shell("bwa mem -v 1 -t {threads} -R '@RG\\tID:{READ_GROUP}\\tPL:ILLUMINA\\tLB:{wildcards.unit}\\tSM:{wildcards.sample}' {params.db_prefix} {input.forward} {input.reverse} >{output}")


rule uncompress:
    input:
        forward = TEMP_DIR + "trimmed/" + "{sample}_{unit}_R1_trimmed.fq.gz",
        reverse = TEMP_DIR + "trimmed/" + "{sample}_{unit}_R2_trimmed.fq.gz"
    output:
        forward = temp(TEMP_DIR + "trimmed/{sample}_{unit}_forward.fastq"),
        reverse = temp(TEMP_DIR + "trimmed/{sample}_{unit}_reverse.fastq")
    message:"uncompressing {wildcards.sample} {wildcards.unit} reads"
    run:
        if is_single_end(wildcards.sample, wildcards.unit):
            shell("gzip -cd {input.forward} > {output.forward};touch {output.reverse}")
        else:
            shell("gzip -cd {input.forward} > {output.forward}")
            shell("gzip -cd {input.reverse} > {output.reverse}")

rule bwa_index:
    input:
        genome = config["refs"]["genome"]
    output:
        TEMP_DIR + "index/genome.amb",
        TEMP_DIR + "index/genome.ann",
        TEMP_DIR + "index/genome.pac",
        TEMP_DIR + "index/genome.sa",
        TEMP_DIR + "index/genome.bwt" 
    message:"building BWA index for the genomic reference"
    params:
        db_prefix = TEMP_DIR + "index/genome"
    shell:
        "bwa index -p {params.db_prefix} {input}"


####################
# QC + trimming
####################
rule multiqc:
    input:
        expand(TEMP_DIR + "fastp/{sample}_{unit}_fastp.json", sample = SAMPLES, unit = UNITS)
    output:
        RESULT_DIR + "multiqc_report.html"
    params:
        fastp_directory = TEMP_DIR + "fastp/",
        outdir = RESULT_DIR
    message: "Summarising fastp reports with multiqc"
    shell:
        "multiqc --force "
        "--outdir {params.outdir} "
        "{params.fastp_directory} "


rule fastp:
    input:
        get_fastq
    output:
        fq1  = temp(TEMP_DIR + "trimmed/" + "{sample}_{unit}_R1_trimmed.fq.gz"),
        fq2  = temp(TEMP_DIR + "trimmed/" + "{sample}_{unit}_R2_trimmed.fq.gz"),
        html = TEMP_DIR + "fastp/{sample}_{unit}_fastp.html",
        json = TEMP_DIR + "fastp/{sample}_{unit}_fastp.json"
    message:
        "trimming {wildcards.sample} reads from {wildcards.unit}"
    threads: 4
    log:
        RESULT_DIR + "fastp/{sample}_{unit}.log.txt"
    params:
        sample_name = "{sample}",
        unit_name = "{unit}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"],
        maximum_read_length = config["fastp"]["maximum_length"]
    run:
        if is_single_end(wildcards.sample, wildcards.unit):
            shell("fastp --thread {threads}  --html {output.html} --json {output.json} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --max_len1 {params.maximum_read_length} \
            --in1 {input} --out1 {output} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} --json {output.json} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --max_len1 {params.maximum_read_length} \
            --max_len2 {params.maximum_read_length} \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")
