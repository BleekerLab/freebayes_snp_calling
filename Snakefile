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


def merge_bams(wildcards):
    "collects all bams files corresponding to the same library"
    bam_files = glob(TEMP_DIR + "mapped/" + wildcards.sample + "_" + "L[0-9]+\.sorted.dedup.bam")
    return bam_files



#################
# Desired output
#################
QC = expand(RESULT_DIR + "fastp/{sample}_{unit}.html", 
    sample=SAMPLES, 
    unit=UNITS)

BAMS = expand(TEMP_DIR + "mapped/{sample}_{unit}.bam",
    sample=SAMPLES,
    unit=UNITS)

VCFs = expand(RESULT_DIR + "vcf/{sample}.vcf",
    sample=SAMPLES)

if config["remove_workdir"]:
    rule all:
        input:
            QC,
            BAMS,
            VCFs
        message:"all done! Cleaning working directory"
        shell:
            "rm -r {TEMP_DIR}"
else: 
    rule all:
        input:
            QC,
            BAMS,
            VCFs
        message:"All done! Keeping temporary directory"
  
###################
# Save master files
###################
rule copy_master_files:
    input:
        "Snakefile",
        "config.yaml",
        "units.tsv",
        "environment.yaml"
    output:
        RESULT_DIR + "Snakefile",
        RESULT_DIR + "config.yaml",
        RESULT_DIR + "units.tsv",
        RESULT_DIR + "environment.yaml"
    message:
        "copying master files"
    shell:
         "copy {input} {RESULT_DIR}"

##############################
# Call SNPs with freebayes
##############################
rule call_variants:
    input:
        ref = config["refs"]["genome"],
        bam = RESULT_DIR + "mapped/{sample}.bam"
    output:
        RESULT_DIR + "vcf/{sample}.vcf"
    message:
        "calling variants for {wildcards.sample}"
    threads: 10
    shell:
        "freebayes -f {input.ref} {input.bam} > {output}"   

##############################
# Merge BAMs from same library
##############################

rule merge_bams:
    input:
        expand(TEMP_DIR + "mapped/{{sample}}_{unit}.sorted.dedup.bam",unit=UNITS)
    output:
        RESULT_DIR + "mapped/{sample}.bam"
    message:
        "merging all BAM files for {wildcards.sample}"
    shell:
        "samtools merge {output} {input}"


##########################################
# Align to genome, sort and mark duplicate
##########################################
rule mark_duplicate:
    input:
        TEMP_DIR + "mapped/{sample}_{unit}.sorted.bam"
    output:
        TEMP_DIR + "mapped/{sample}_{unit}.sorted.dedup.bam"
    message:
        "marking duplicates in {wildcards.sample} {wildcards.unit} bam file"
    log:
        RESULT_DIR + "logs/picard/{sample}.{unit}.metrics.txt"
    shell:
        "picard MarkDuplicates "
        "I={input} "
        "O={output} "
        "M={log} "
        "REMOVE_DUPLICATES=false"


rule samtools_sort:
    input:
        TEMP_DIR + "mapped/{sample}_{unit}.bam"
    output:
        TEMP_DIR + "mapped/{sample}_{unit}.sorted.bam"
    message:"sorting {wildcards.sample} {wildcards.unit} bam file"
    threads: 5
    shell:
        "samtools sort -@ {threads} {input} > {output}"
 

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
    threads: 10
    run:
        # Building the read group id (sequencer_id + flowcell_name + lane_number + barcode)
        SEQUENCER_ID=check_output("head -n 1 " + input.forward + " |cut -d: -f1",shell=True).decode().strip()
        FLOWCELL_NAME=check_output("head -n 1 " + input.forward + " |cut -d: -f3",shell=True).decode().strip()
        FLOWCELL_LANE=check_output("head -n 1 " + input.forward + " |cut -d: -f4",shell=True).decode().strip()
        BARCODE=check_output("head -n 1 " + input.forward + " |cut -d' ' -f2 |cut -d: -f4",shell=True).decode().strip()
        # Feeding the READ_GROUP_ID to bwa
        READ_GROUP = SEQUENCER_ID + "." + FLOWCELL_NAME + "." + FLOWCELL_LANE + "." + BARCODE
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
rule fastp:
    input:
        get_fastq
    output:
        fq1  = TEMP_DIR + "trimmed/" + "{sample}_{unit}_R1_trimmed.fq.gz",
        fq2  = TEMP_DIR + "trimmed/" + "{sample}_{unit}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}_{unit}.html",
        json = RESULT_DIR + "fastp/{sample}_{unit}.json"
    message:"trimming {wildcards.sample} reads from {wildcards.unit}"
    threads: 10
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
