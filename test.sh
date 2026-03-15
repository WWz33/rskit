#!/bin/bash

#SBATCH --qos=qos_hejianbo
#SBATCH -o %x.%j.out
#SBATCH -N 1
#SBATCH -J salmon
#SBATCH -p sri-com
#SBATCH -c 96
for i in lhy-D-rep2 lhy-D-rep1
do
python rskit/cli.py quant -s ${i} -1 ../${i}_1 -2 ../${i}_2 --index-dir ../STAR_Gmax_index/  -t 96 -g ../Gmax_508_Wm82.a4.v1.fa -gtf ../Gmax_508_Wm82.a4.v1.gene_exons.gff3.gtf -gf ../gffread_transcript.fa  -o data/
done
