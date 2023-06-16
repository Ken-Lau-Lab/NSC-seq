#!/bin/bash
for file in output/*.csv
do
	echo $file
        Rscript Processing_bulk_DNA_data.R $file
done
