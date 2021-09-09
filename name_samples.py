### Name samples according to tumor type and put in seperate folders
# import libraries
import glob
import os

# Read arguments from command line
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', help='path to samples files')
parser.add_argument('--names', '-n', help='path to list with new names')
parser.add_argument('--output', '-o', help='output directory')
args = parser.parse_args()

# Input files
input = args.input
if array == True:
    samples = glob.glob(input + '*.txt')
if NGS == True:
    samples = glob.glob(input + '*.cov')
names = args.manifest

# Output files
output = args.output
sample_output = output + "samples_clustered.csv"
metadata_output = output + "metadata_clusters.csv"

# replace sample name with folder name for easy annotation in reference dataframe
cfRRBSfolders = glob.glob("./reference/cfRRBSreferences/*/")
for folder in cfRRBSfolders:
    samples = glob.glob(os.path.join(folder, "*.cov"))
    label = os.path.basename(os.path.dirname(folder))
    i = 1
    for sample in samples:
        dst = os.path.dirname(sample) + "/" + label + "_" + str(i) + ".cov"
        os.rename(sample, dst)
        i += 1

names_df = pd.read_csv(names, sep=';',header=None,index_col=[0],names=['names'])
samples_df = samples_df.transpose()
samples_df = pd.merge(samples_df,names_df,how='left',left_index=True,right_index=True)
samples_df = samples_df.set_index('names')
samples_df = samples_df.transpose()
samples_df['clusterID'] = samples_df['clusterID'].astype(int)

## using cfRRBS samples as references: group samples per tumor type and take median
samples_df.columns = samples_df.columns.str.split(' ').str[0]
samples_df = samples_df.groupby(by=samples_df.columns, axis=1).median()
