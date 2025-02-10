#!/usr/bin/env bash

file=$1
sample=${file%.*}

for ct in alpha beta gamma delta acinar ductal stellate immune EC; do
    samtools view -H $file > SAM_header
    samtools view $file | LC_ALL=C grep -F -f /media/leon/Masha/ATAC/sample_ct_barcodes/${sample}_${ct}.txt > filtered_SAM_body
    cat SAM_header filtered_SAM_body | samtools view -b > "${sample}_${ct}.bam"
    rm SAM_header filtered_SAM_body
done
