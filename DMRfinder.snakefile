# Parallelisation options
import multiprocessing
import sys
import itertools
import os
import collections
import json
import glob
import pandas as pd
pd.options.mode.chained_assignment = None

# input files
configfile: "config.yaml"
nameList = config["nameList"]
names = pd.read_csv(nameList)
sampleFolder = config["sampleFolder"].rstrip("/")
samples = glob.glob(sampleFolder + '/*.cov')

# scripts
scriptFolder = config["scriptFolder"].rstrip("/")

# output files
outputFolder = config["outputFolder"].rstrip("/")
report: "report/workflow.rst"

# check if samples in the nameList are also in the sampleFolder:
# if yes, replace ID in nameList by sample-path. If not, remove.
for sample in samples:
    file_name = os.path.splitext(os.path.basename(sample))[0]
    smapleID = file_name.split('_')[0]
    if smapleID in names.values:
         names.replace([smapleID], sample, regex=True, inplace=True)
path_nameDF = names[names.smapleID.str.contains(sampleFolder)]
missing_samples = names[~names.smapleID.str.contains(sampleFolder)]
print(f"There are {len(missing_samples)} missing samples: \n{missing_samples}")

# make all possible combos with the available tumor types
path_nameDF['name'] = path_nameDF['name'].map(lambda x: x.split(' ')[0])
tumor_types = path_nameDF['name'].unique()
combos = list(itertools.combinations(tumor_types, 2))

# group samples together according to tumor type
grouped_samples = dict(path_nameDF.groupby('name')['smapleID'].apply(list))

# link combos to sample paths
combo_dict = {}
for combo in combos:
    T1, T2 = combo
    G1 = grouped_samples[T1]
    G2 = grouped_samples[T2]
    combo = '-'.join(combo)
    combo_dict[combo] = G1+G2

if config["method"] == DMRfinder:

    rule all:
        input:
            expand(outputFolder + "/results{COMBOS}.csv", COMBOS = list(combo_dict.keys()))

    # combine CpGs into regions
    rule combine_CpGs:
        input:
            lambda wildcards: combo_dict[wildcards.COMBOS]
        output:
            outputFolder + "/combined{COMBOS}.csv"
        params:
            s = scriptFolder + "/combine_CpG_sites.py"
        shell:
            "python3 {params.s} -v -o {output} {input}"

    # Test regions for differential methylation
    ## In params section the tumor groups (T1 and T2) and samples paths of those tumor groups (G1 and G2) are specified:
    ## in G1 you cannot simply refer to T1 because when T1 is re-evaluated (next combo in lambda function), G1 is not automatically re-evaluated as well.
    ## thats why the lambda function in T1 needs to be repeated in G1
    rule find_DMRs:
        input:
            outputFolder + "/combined{COMBOS}.csv"
        output:
            outputFolder + "/results{COMBOS}.csv"
        params:
            s = scriptFolder + "/findDMRs.r",
            T1 = lambda wildcards: wildcards.COMBOS.split('-')[0],
            T2 = lambda wildcards: wildcards.COMBOS.split('-')[1],
            G1 = lambda wildcards: str([os.path.splitext(os.path.basename(sample))[0].split('.')[0] for sample in grouped_samples[wildcards.COMBOS.split('-')[0]]]).strip("[]").replace(" ",""),
            G2 = lambda wildcards: str([os.path.splitext(os.path.basename(sample))[0].split('.')[0] for sample in grouped_samples[wildcards.COMBOS.split('-')[1]]]).strip("[]").replace(" ","")
        shell:
            "ml purge && ml R/4.1.2-foss-2021b; "
            "Rscript {params.s} -i {input} -o {output} -v {params.G1} {params.G2} -n {params.T1},{params.T2}"

if config["method"] == DSS:

    rule all:
        input:
            expand(outputFolder + "/results{COMBOS}.csv", COMBOS = list(combo_dict.keys()))

    # Test regions for differential methylation
    ## In params section the tumor groups (T1 and T2) and samples paths of those tumor groups (G1 and G2) are specified:
    ## in G1 you cannot simply refer to T1 because when T1 is re-evaluated (next combo in lambda function), G1 is not automatically re-evaluated as well.
    ## thats why the lambda function in T1 needs to be repeated in G1
    rule find_DMRs:
        input:
            G1 = lambda wildcards: grouped_samples[wildcards.COMBOS.split('-')[0]],
            G2 = lambda wildcards: grouped_samples[wildcards.COMBOS.split('-')[1]]
        output:
            outputFolder + "/results{COMBOS}.csv"
        params:
            s = scriptFolder + "/DMRs_DSS.R",
            T1 = lambda wildcards: wildcards.COMBOS.split('-')[0],
            T2 = lambda wildcards: wildcards.COMBOS.split('-')[1]

        shell:
            "ml purge && ml R/4.1.2-foss-2021b; "
            "Rscript {params.s} -g1 {input.G1} -g2 {input.G2} -n1 {params.T1} -n2 {params.T2} -o {output}"

else:

    sys.exit("specify valid method")
