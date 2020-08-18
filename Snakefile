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

WORKING_DIR = config["workdir"]
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

def sample_is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    if "fq2" not in samples.columns:
        return True
    else:
        return pd.isnull(samples.loc[(sample), "fq2"])

def get_fastq(wildcards):
    """ This function checks if the sample has paired end or single end reads
    and returns 1 or 2 names of the fastq files """
    if sample_is_single_end(wildcards.sample):
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna()
    else:
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1","fq2"]].dropna()


def merge_bams(wildcards):
    "collects all bams files corresponding to the same library"
    bam_files = glob(WORKING_DIR + "mapped/" + wildcards.sample + "_" + "L[0-9]+\.sorted.dedup.bam")
    return bam_files



#################
# Desired output
#################
QC = expand(RESULT_DIR + "fastp/{sample}_{unit}.html", 
    sample=SAMPLES, 
    unit=UNITS)

BAMS = expand(WORKING_DIR + "mapped/{sample}_{unit}.bam",
    sample=SAMPLES,
    unit=UNITS)

VCFs = expand(RESULT_DIR + "vcf/{sample}.vcf",
    sample=SAMPLES)

rule all:
    input:
        QC,
        BAMS,
        VCFs
    message:"all done! Cleaning working directory"
    shell:
        "rm -r {WORKING_DIR}"
  
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
    message:"copying master files"
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
    message:"calling variants for {wildcards.sample}"
    threads: 10
    shell:
        "freebayes -f {input.ref} {input.bam} > {output}"   

##############################
# Merge BAMs from same library
##############################

rule merge_bams:
    input:
        #merge_bams
        expand(WORKING_DIR + "mapped/{{sample}}_{unit}.sorted.dedup.bam",unit=UNITS)
    output:
        RESULT_DIR + "mapped/{sample}.bam"
    message:"merging all BAM files for {wildcards.sample}"
    shell:
        "samtools merge {output} {input}"


##########################################
# Align to genome, sort and mark duplicate
##########################################
rule mark_duplicate:
    input:
        WORKING_DIR + "mapped/{sample}_{unit}.sorted.bam"
    output:
        WORKING_DIR + "mapped/{sample}_{unit}.sorted.dedup.bam"
    message:"marking duplicates in {wildcards.sample} {wildcards.unit} bam file"
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
        WORKING_DIR + "mapped/{sample}_{unit}.bam"
    output:
        WORKING_DIR + "mapped/{sample}_{unit}.sorted.bam"
    message:"sorting {wildcards.sample} {wildcards.unit} bam file"
    threads: 5
    shell:
        "samtools sort -@ {threads} {input} > {output}"

# TODO: check if read_group assignment works with BWA MEM
# TODO: if it does not work, try to replace it with Picard
#rule add_read_groups:
#    input:
#       WORKING_DIR + "mapped/{sample}_{unit}.bam"
#    output:
#       WORKING_DIR + "mapped/{sample}_{unit}.rg.bam"
#    message:"adding read group to {input}"
#    shell:
#        "picard AddOrReplaceReadGroups "
#        "I={input} "
#        "O={output}
#
 

rule bwa_align:
    input:
        index = [WORKING_DIR + "index/genome." + ext for ext in ["amb","ann","pac","sa","bwt"]],
        forward = WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq",
        reverse = WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq"
    output:
        WORKING_DIR + "mapped/{sample}_{unit}.bam"
    message:"mapping {wildcards.sample} {wildcards.unit} reads to genomic reference"
    params:
        db_prefix = WORKING_DIR + "index/genome"
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
        forward = WORKING_DIR + "trimmed/" + "{sample}_{unit}_R1_trimmed.fq.gz",
        reverse = WORKING_DIR + "trimmed/" + "{sample}_{unit}_R2_trimmed.fq.gz"
    output:
        forward = WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq",
        reverse = WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq"
    message:"uncompressing {wildcards.sample} {wildcards.unit} reads"
    shell:
        "gzip -cd {input.forward} > {output.forward};"
        "gzip -cd {input.reverse} > {output.reverse}"

rule bwa_index:
    input:
        genome = config["refs"]["genome"]
    output:
        WORKING_DIR + "index/genome.amb",
        WORKING_DIR + "index/genome.ann",
        WORKING_DIR + "index/genome.pac",
        WORKING_DIR + "index/genome.sa",
        WORKING_DIR + "index/genome.bwt" 
    message:"building BWA index for the genomic reference"
    params:
        db_prefix = WORKING_DIR + "index/genome"
    shell:
        "bwa index -p {params.db_prefix} {input}"


####################
# QC + trimming
####################
rule fastp:
    input:
        get_fastq
    output:
        fq1  = WORKING_DIR + "trimmed/" + "{sample}_{unit}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/" + "{sample}_{unit}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}_{unit}.html",
        json = RESULT_DIR + "fastp/{sample}_{unit}.json"
    message:"trimming {wildcards.sample} reads from {wildcards.unit}"
    threads: 10
    log:
        RESULT_DIR + "fastp/{sample}_{unit}.log.txt"
    params:
        sampleName = "{sample}",
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    run:
        if sample_is_single_end(params.sampleName):
            shell("fastp --thread {threads}  --html {output.html} --json {output.json} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --in1 {input} --out1 {output} \
            2> {log}; \
            touch {output.fq2}")
        else:
            shell("fastp --thread {threads}  --html {output.html} --json {output.json} \
            --qualified_quality_phred {params.qualified_quality_phred} \
            --detect_adapter_for_pe \
            --in1 {input[0]} --in2 {input[1]} --out1 {output.fq1} --out2 {output.fq2}; \
            2> {log}")
