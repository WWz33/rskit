#!/bin/bash

#SBATCH --qos=qos_hejianbo
#SBATCH -o %x.%j.out
#SBATCH -N 1
#SBATCH -J salmon
#SBATCH -p sri-com
#SBATCH -c 96

rskit all -S coldata.csv -g Gmax_508_Wm82.a4.v1.fa -gtf Gmax_508_Wm82.a4.v1.gene_exons.gff3.gtf -gf gffread_transcript.fa -o rnaseq-tools/data/ -p 90 --trim  
