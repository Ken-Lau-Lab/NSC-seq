#!/bin/bash
len=$(wc -l < mito.txt)
for i  in $(seq 1 $len) 
do
	content=($(sed -n "${i} p" mito.txt))
	bash loop_mito.sh ${content[1]} ${content[2]::-1} ${content[0]}
done






