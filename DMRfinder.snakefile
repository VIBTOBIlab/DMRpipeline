# Parallelisation options
import multiprocessing
import sys
import itertools
import os
import collections
import json
import glob
import pandas as pd
cpuCount = (multiprocessing.cpu_count() - 2)

configfile: "config.yaml"
nameList: config["nameList"]
report: "report/workflow.rst"

rule all:
    input:
        "multiqc_report.html"

# select and name samples
names = pd.read_csv(nameList)

rule name_samples:
    input:


#Extract methylation counts from alignment files
rule sort_samtools:
    input:
        bam = "rm_optical_dups/{sample}_R1_001_val_1_bismark_bt2_pe.bam"
    output:
        bam = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam",
        bai = "sorted_reads/{sample}_R1_001_val_1_bismark_bt2_pe.bam.bai"
    threads: 4
    shell:
        "ml purge && ml SAMtools/1.9-intel-2018b; "
        "samtools sort -@ {threads} {input.bam} -o {output.bam}; "
        "samtools index {output.bam} {output.bai} "
