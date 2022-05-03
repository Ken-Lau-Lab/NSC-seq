#!/bin/bash
for file in *R2*.gz
do
        output=$(echo $file | sed 's/_S1_L005_R2_001.fastq.gz.*//')
        output+='_R2.txt'
        echo $output
        zgrep -A 0 -B 0 'CTAGAAATAG' $file | sed '/--/d' >$output
        #new line
done


for file in *R2.txt
do
        output=$(echo $file | sed 's/R2.txt.*//')
        output+='R2_trimmed.txt'
        echo $output
        sed 's/GTTAACCTAAGGCTA..*//' $file | cut -c 9-|sed 's/^GGGT/GGT/g'| sed 's/^GGGGT/GGT/g'|grep ^GGT >$output
done


for file in *R2_trimmed.txt
do
        Rscript mutation_barcodes_test.R $file
done


for file in *R2_trimmed_mutation.txt
do
        Rscript mutation_calling.R $file
done

