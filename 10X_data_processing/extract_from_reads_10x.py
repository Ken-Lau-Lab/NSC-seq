# Copy-pasted from extract_from_reads.ipynb for easy running
# Inputs: Header and read data
# Outputs: .csv table of combined values
# Copy-pasted as of 11/4/2021

import argparse # TODO: maybe use arguments?
import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# from Bio.Seq import Seq # TODO

filename = sys.argv[1]

# These are the 10x Grep Sequences
Grep_Seqs = ["CAACACATAA", "ACAATCCTAG", "AAAATCTAGA", "TACCCTATAG", "TGAGGGTCTT", "AAATGGGCCT",
             "CTCAAACCTA", "ATAATCTTAT", "TGAACCGAAT", "CCAATGCTAA", "AATCACATAA", "GGGGAATAGG"]

# Read text files as dataframes
Seq_dfs = {}
for seq in Grep_Seqs:
    R2 = pd.read_table('./Data/'+filename+'_'+seq+'.txt', header=None, sep='\n')
    # Collapse header and read lines
    R2 = pd.DataFrame(np.reshape(R2.values,(int(R2.shape[0] / 2),2)), columns=['Header', 0])
    # Set grep sequence
    R2['GrepSeq'] = seq

    Seq_dfs[seq] = R2


## Helper functions ##
def extract_seq(df, name, pattern):
    """
    Extracts sequences that follow pattern and adds them as new row into dataframe

    Parameters:
        df (pd.DataFrame): Dataframe with individual reads. df[0] must contain reads.
        name (string): Name of Dataframe column to be added
        pattern (re.Pattern): Regular Expression pattern to be applied to df[0]

    Returns:
        df (pd.DataFrame): Inputted dataframe with added column. Will modify inputted dataframe.
    """
    df[name] = df[0].apply(pattern.findall)
    return df

def extract_from_list(x):
    """ TODO: Document this. I feel like this shouldn't be necessary but we need it rn """
    if len(x) > 0:
        return x[0]
    else:
        return ""

## TODO: From cody <-- credit him
def reverse_comp(sequence, nucleic_acid):
    """
    Generate reverse compliment to given sequence, s,
    in desired nucleic acid format (i.e. DNA or RNA).
    Output is reported 3' to 5'.
​
    Parameters:
        sequence (str): sequence to test.
        nucleic_acid (str): "DNA" or "RNA" defining desired output
    """
    # define base complements
    bases = {
        "A": {"DNA": "T", "RNA": "U"},
        "T": {"DNA": "A", "RNA": "A"},
        "G": {"DNA": "C", "RNA": "C"},
        "C": {"DNA": "G", "RNA": "G"},
        "N": {"DNA": "N", "RNA": "N"},
        "a": {"DNA": "T", "RNA": "U"},
        "t": {"DNA": "A", "RNA": "A"},
        "g": {"DNA": "C", "RNA": "C"},
        "c": {"DNA": "G", "RNA": "G"},
        "n": {"DNA": "n", "RNA": "n"},
        " ": {"DNA": " ", "RNA": " "},
        "\t": {"DNA": "\t", "RNA": "\t"},
        "\n": {"DNA": "\n", "RNA": "\n"},
        ",": {"DNA": ",", "RNA": ","},
    }
    # return reverse complimentary string (3' - 5')
    return "".join([bases[x][nucleic_acid] for x in sequence])[::-1]

## Extract Information from reads ##

## Set up regex patterns
# the key will be the column name, and the value will be applied to the reads

# Patterns
patterns = {}
for seq in Grep_Seqs:
    patterns[seq] = re.compile('(.*)'+seq)

## Apply regex patterns

# TODO: there's probably a better way to iterate through these items
# Iterate through patterns and add columns
# Also extract all from the 1 element list
for pattern in patterns.items():
    extract_seq(Seq_dfs[pattern[0]], 'hBC', pattern[1])
    Seq_dfs[pattern[0]]['hBC'] = Seq_dfs[pattern[0]]['hBC'].apply(extract_from_list)


## Filtering and post processing
# Loop through grep sequence dataframes and merge
Combined = pd.DataFrame()
for df in Seq_dfs.items():
    Combined = Combined.append(df[1], ignore_index=True)

# Filter out pts where BC1 or BC2 are missing (Cannot determine what cell it came from)
# Combined['BC1_R1'].replace('', np.nan, inplace=True)
# Combined.dropna(subset=['BC1_R1', 'BC2_R2'], inplace=True)

print(Combined)


## Identify cells and umi counts per cell

# # All unique Cells
# cells = np.unique(Combined["BC1_BC2"])

# umi_counts = {}
# for cell in cells:
#     count = len(np.unique(Combined[Combined["BC1_BC2"] == cell]["UMI_R1"]))
#     umi_counts.update({cell: count})

# # Make a umi/cell plot with cell rank and stuff
# umi_counts_df = pd.DataFrame(list(umi_counts.items()),columns = ['Cell','UMI Count']) 
# # Rank each cell
# umi_counts_df = umi_counts_df.sort_values('UMI Count', ascending=False)
# umi_counts_df['Cell Rank'] = range(1, len(umi_counts_df) + 1)

# Using outer merge now
# Combined = Combined.merge(umi_counts_df, how='outer', left_on='BC1_BC2', right_on='Cell')

## This is the final combined dataframe
Combined.to_csv("./Sequences/"+filename+"_10x.csv")
#Combined.to_csv("./Sequences/"+filename+'-TATTTTAACTTG- CAAGTTAAAATA.csv')
