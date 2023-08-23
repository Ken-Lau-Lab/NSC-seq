# Will isolate sequences and input them into MatchR1R2.py (DNA)
# Output is a table of matched reads
# Input:
# $1 Fastq file for R2

# TODO: Adjust arguments

declare -a Grep_Seqs=("CAACACATAA" "ACAATCCTAG" "AAAATCTAGA" "TACCCTATAG" "TGAGGGTCTT" "AAATGGGCCT" "CTCAAACCTA" "ATAATCTTAT" "TGAACCGAAT" "CCAATGCTAA" "AATCACATAA" "GGGGAATAGG")

for grep_seq in ${Grep_Seqs[@]}
do
    echo running $grep_seq
#    zgrep --no-group-separator -A 0 -B 1 $3 $1 > ./Data/$(basename $1 _S1_L005_R1_001.fastq.gz)_R1.txt
    zgrep --no-group-separator -A 0 -B 1 $grep_seq $1 > ./Data/$(basename $1 _S1_L001_R2_001.fastq.gz)_$grep_seq.txt

done

python3 ./extract_from_reads_10x.py $(basename $1 _S1_L001_R2_001.fastq.gz)