# Runs loop to extract sequences on all files in folder
# Output is tables of matched reads
# Input:
# $1 is sequence for R1
# $2 is sequence for R2
# e.g. ./loop.sh ACGGACTAGCCT AGGCTAGTCCGT

R1=($(find ../* | grep fastq | grep R1))
R2=($(find ../* | grep fastq | grep R2))
#R1=($(find ~/Documents/BarcodeIndelseek/0-raw-data/* | grep fastq | grep R1))
#R2=($(find ~/Documents/BarcodeIndelseek/0-raw-data/* | grep fastq | grep R2))

LEN=${#R1[@]}
wd=$(pwd)
echo Running Extract Seq on files with sequences $1 and $2

for i in $(seq 0 $(($LEN-1)))
do
	echo ${R1[i]}
	#echo $i
	zgrep --no-group-separator -A 0 -B 1 $1 ${R1[i]} > ./Data/$(basename ${R1[i]} _data_R1.fastq.gz)_R1.txt
	zgrep --no-group-separator -A 0 -B 1 $2 ${R2[i]} > ./Data/$(basename ${R2[i]} _data_R2.fastq.gz)_R2.txt
	python3 ./extract_from_reads_mito.py $(basename ${R1[i]} _data_R1.fastq.gz) $3
	echo $(basename ${R1[i]} _data_R1.fastq.gz) completed
done

echo Finished Loop

