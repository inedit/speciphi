#!/bin/bash
#PBS -N 144g_pile
#PBS -e 144g_pile.err
#PBS -o 144g_pile.out
#PBS -l select=1:ncpus=12:mem=23Gb
#PBS -l cput=100:00:00
#PBS -l walltime=200:00:00
#PBS -l place=pack:shared
#PBS -q standard
#PBS -W group_list=mbsulli
#PBS -m bea
#PBS -M ssolonen@email.arizona.edu

# loading modules and adjusting environmental variables
. /usr/share/Modules/init/bash
module load samtools/0.1.18

/usr/bin/time -v samtools mpileup -f /rsgrps1/mbsulli/sergei/ann_genomes/144genomes.fa /rsgrps1/mbsulli/sergei/ann_genomes/bowtie_runs/very_sensitive_local_concatenated_144genomes.sorted.bam >very_sensitive_local_concatenated_144genomes.mpileup 2>very_sensitive_local_concatenated_144genomes.mpileup.err
