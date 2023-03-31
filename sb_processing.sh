#!/bin/bash

# Memory request for 12G
#$ -l h_vmem=30g

# Cores
#$ -pe smp 1
#$ -binding linear:1

# Runtime request
#$ -l h_rt=4:00:00
#$ -l os=RedHat7
#$ -R y

# use these
source /broad/software/scripts/useuse
reuse UGER
reuse Glimpse
reuse Google-Cloud-SDK
use .python-3.8.3
use Anaconda3
reuse -q Anaconda



######## Change these parameters ########

sample=samplename

slidetags_directory=/path/to/slidetags/directory

R1_path=/path/to/R1.fastq.gz
R2_path=/path/to/R2.fastq.gz



######## grep for UP site ########

reads=25000000

## get line numbers to extract
zgrep -n "TCTTCAGCGTTCCCGAGA" ${R2_path} | cut -f1 -d: > ${slidetags_directory}/slidetags_${sample}/up_line_numbers_${sample}.txt


## extract line numbers from file
zcat ${R1_path} | awk 'NR==FNR{for(i=($1-1);i<=($1+2);i++)a[i];next}FNR in a' ${slidetags_directory}/slidetags_${sample}/up_line_numbers_${sample}.txt - > ${slidetags_directory}/slidetags_${sample}/up_only_${sample}_R1.fastq

zcat ${R2_path} | awk 'NR==FNR{for(i=($1-1);i<=($1+2);i++)a[i];next}FNR in a' ${slidetags_directory}/slidetags_${sample}/up_line_numbers_${sample}.txt - > ${slidetags_directory}/slidetags_${sample}/up_only_${sample}_R2.fastq


## output metrics
echo "output stats for" ${sample} ${samplename}
echo "number of lines in fastq"
zcat ${R2_path} | wc -l

echo "number of lines in fastqs with exact UP"
wc -l  ${slidetags_directory}/slidetags_${sample}/up_only_${sample}_R1.fastq ${slidetags_directory}/slidetags_${sample}/up_only_${sample}_R2.fastq
echo "number of reads with approximate UP"
zcat ${R2_path} | agrep -1 -n  "TCTTCAGCGTTCCCGAGA" | cut -f1 -d: > ${slidetags_directory}/slidetags_${sample}/up_hamming_1_line_numbers_${sample}.txt

wc -l ${slidetags_directory}/slidetags_${sample}/up_hamming_1_line_numbers_${sample}.txt



######## downsample to specified number of reads ########

# downsample
seqtk sample -s 100 ${slidetags_directory}/slidetags_${sample}/up_only_${sample}_R1.fastq ${reads} > ${slidetags_directory}/slidetags_${sample}/sub_${sample}_${reads}_R1.fastq
seqtk sample -s 100 ${slidetags_directory}/slidetags_${sample}/up_only_${sample}_R2.fastq ${reads} > ${slidetags_directory}/slidetags_${sample}/sub_${sample}_${reads}_R2.fastq







