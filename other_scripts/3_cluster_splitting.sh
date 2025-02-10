#!/bin/bash

for cluster in $(seq 15 -1 0); do
  for sample in $(seq 61 94); do
    cd /media/leon/Masha/ATAC/CR_output/sample${sample}/outs
    samtools view -H filtered_${sample}.bam > SAM_header
    samtools view filtered_${sample}.bam | LC_ALL=C grep -F -f ${cluster}.txt > filtered_SAM_body
    cat SAM_header filtered_SAM_body > ${sample}_${cluster}.sam
    samtools view -b ${i}_${cluster}.sam > ${sample}_${cluster}.bam
    rm SAM_header filtered_SAM_body ${sample}_${cluster}.sam
  done
  samtools merge - /media/leon/Masha/ATAC/CR_output/sample*/outs/*${cluster}.bam | samtools sort -o /media/leon/Masha/ATAC/clusters/${cluster}_merged.bam -
done
