#!/bin/bash
len=$(wc -l < mito.txt)
for i  in $(seq 1 $len) 
do
	#echo $i
	content=($(sed -n "${i} p" mito.txt))
	bash loop_mito.sh ${content[1]} ${content[2]::-1} ${content[0]}
done

#content=($(sed -n '1 p' mito.txt))

#echo ${#content[0]}
#echo ${#content[1]}
#echo ${content[2]::-1}

#content=($(sed -n "2p" mito.txt))
#bash loop_mito.sh ${content[1]} ${content[2]::-1} ${content[0]}





