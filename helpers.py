#!/usr/bin/env python3

import os

def get_fastq_file_name(fastq_file):
	if fastq_file.endswith(".fastq.gz"):
		file_name = os.path.basename(path_to_file).split(".fastq.gz")[0]
		return file_name
	if fastq_file.endswith(".fq.gz"):
		file_name = os.path.basename(path_to_file).split(".fq.gz")[0]
		return file_name
	if fastq_file.endswith(".fastq"):
		file_name = os.path.basename(path_to_file).split(".fastq")[0]
		return file_name
	if fastq_file.endswith(".fq"):
		file_name = os.path.basename(path_to_file).split(".fq")[0]
		return file_name
