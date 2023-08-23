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
mito = sys.argv[2]

# Read text files as dataframes
R1 = pd.read_table('./Data/'+filename+'_R1.txt', header=None)
R2 = pd.read_table('./Data/'+filename+'_R2.txt', header=None)

# Collapse header and read lines
R1 = pd.DataFrame(np.reshape(R1.values,(int(R1.shape[0] / 2),2)), columns=['Header', 0])
R2 = pd.DataFrame(np.reshape(R2.values,(int(R2.shape[0] / 2),2)), columns=['Header', 0])

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
# R1 patterns
R1_patterns = {
    # hBC is barcode that falls between these 2 sequences
    'hBC': re.compile('AACTTGCTAT(.*)ACCCC'),
    # BC1 is from the start to the start of the W1 sequence
    'BC1': re.compile('^(.*)GAGTGATTGC'),
    # BC2 and UMI are a 14 bp sequence right after the W1 sequence
    'BC2_UMI': re.compile('GTGACGCCTT(.{1,14})'), # TODO: 1,14 means 1 to 14 characters
    # The TSO UMI is at the very end of the read
    'TSO_UMI': re.compile('ACCCCC(.{6})'),
}

# R2 patterns
R2_patterns = {
    # hBC is barcode that falls between these 2 sequences
    'hBC': re.compile('GGGGT(.*)ATAGCAAGTT'),
    # BC1 is right after the end of the W1 sequence
    'BC1': re.compile('GCAATCACTC(.*)'),
    # BC2 and UMI are a 14 bp sequence right before the W1 sequence
    # 'BC2_UMI': re.compile('(.{14})GCAATCACTC'), OLD PATTERN THAT WAS WRONG
    # 'BC2_UMI': re.compile('(.{14})AAGGCGTCAC'), # This cut-off for larger hBC's
    'BC2_UMI': re.compile('TCACTCCTGA(.{1,14})'), # 1,14 means 1 to 14 characters
    # The TSO UMI is 6bp at the very start of the read
    'TSO_UMI': re.compile('^(.{6})'),
}

## Apply regex patterns

# TODO: there's probably a better way to iterate through these items
# Iterate through patterns and add columns
# Also extract all from the 1 element list
for pattern in R1_patterns.items():
    extract_seq(R1, pattern[0], pattern[1])
    if(pattern[0] == 'TSO_UMI'):
        # Because the first one we don't want; we want the last pattern match b/c it will be at the end
        R1['TSO_UMI'] = R1['TSO_UMI'].apply(lambda x: x[-1:])
    R1[pattern[0]] = R1[pattern[0]].apply(extract_from_list)

for pattern in R2_patterns.items():
    extract_seq(R2, pattern[0], pattern[1])
    R2[pattern[0]] = R2[pattern[0]].apply(extract_from_list)

## Specific modifications
## These include separating UMI and BC2, Adding sequences to hBC, etc.

# Separate UMI and BC2
# R1: 8bp BC2 + 6bp UMI
R1["BC2"] = R1["BC2_UMI"].apply(lambda x: x[0:8])
R1["UMI"] = R1["BC2_UMI"].apply(lambda x: x[8:])
# R2: 6bp UMI + 8bp BC2
R2["UMI"] = R2["BC2_UMI"].apply(lambda x: x[0:6])
R2["BC2"] = R2["BC2_UMI"].apply(lambda x: x[6:]) # TODO: CHECK THIS AND MAKE SURE NUMBERS ARE GOOD!!

# Add ACC to end or GGT to start of hBC for R1 and R2
R1["hBC"] = R1["hBC"].apply(lambda x: (x + "ACC") if(len(x)>0) else x)
R2["hBC"] = R2["hBC"].apply(lambda x: ("GGT" + x) if(len(x)>0) else x) # TODO: Make sure these are correct
# TODO: TODO: Should this not include the start and end? It is often there even if there is no actual sequence???

# Create combined BC1+BC2 column in R1
R1["BC1_BC2"] = R1["BC1"] + R1["BC2"]

## Combining R1 and R2 based on read headers ##
## Process R2 for reverse complement so we can align

# We neeed to get reverse complement of R2 data to align this stuff
R2['BC2_revcomp'] = R2['BC2'].apply(lambda x: reverse_comp(x, 'DNA'))
R2['UMI_revcomp'] = R2['UMI'].apply(lambda x: reverse_comp(x, 'DNA'))

# We will use the BC2 and UMI combo column to align
R2['BC2_UMI_revcomp'] = R2['BC2_UMI'].apply(lambda x: reverse_comp(x, 'DNA'))

## Modify headers so R1 and R2 match
R2['Header'] = R2['Header'].str.replace('2:N:0','1:N:0')

## Combine using pandas merge
#Combined = pd.merge(R1, R2, how='inner', on='Header', suffixes=('_R1', '_R2'))
Combined = pd.merge(R1, R2, how='outer', on='Header', suffixes=('_R1', '_R2'))

## Filtering and post processing

# Filter out pts where BC1 or BC2 are missing (Cannot determine what cell it came from)
Combined['BC1_R1'].replace('', np.nan, inplace=True)
Combined['BC2_R1'].replace('', np.nan, inplace=True)
Combined.dropna(subset=['BC1_R1', 'BC2_R2'], inplace=True)


## Identify cells and umi counts per cell

# All unique Cells
cells = np.unique(Combined["BC1_BC2"])

umi_counts = {}
for cell in cells:
    count = len(np.unique(Combined[Combined["BC1_BC2"] == cell]["UMI_R1"]))
    umi_counts.update({cell: count})

# Make a umi/cell plot with cell rank and stuff
umi_counts_df = pd.DataFrame(list(umi_counts.items()),columns = ['Cell','UMI Count']) 
# Rank each cell
umi_counts_df = umi_counts_df.sort_values('UMI Count', ascending=False)
umi_counts_df['Cell Rank'] = range(1, len(umi_counts_df) + 1)

# TODO: Using outer merge now
Combined = Combined.merge(umi_counts_df, how='outer', left_on='BC1_BC2', right_on='Cell')

#### TODO TODO TODO: Make sure this still works--I rearranged some code

## This is the final combined dataframe
Combined.to_csv("./output_mito/"+filename+'_'+mito+".csv")
#Combined.to_csv("./Sequences/"+filename+'-TATTTTAACTTG- CAAGTTAAAATA.csv')
