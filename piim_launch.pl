open (OUT,">$ARGV[0].pbs");

print OUT
"#!/bin/bash
#PBS -N $ARGV[0]_pm
#PBS -e $ARGV[0]_pm.err
#PBS -o $ARGV[0]_pm.out
#PBS -l select=1:ncpus=1:mem=2Gb
#PBS -l cput=400:00:00
#PBS -l walltime=800:00:00
#PBS -l place=pack:shared
#PBS -q standard
#PBS -W group_list=mbsulli
#PBS -m bea
#PBS -M ssolonen\@email.arizona.edu

# loading modules and adjusting environmental variables
. /usr/share/Modules/init/bash

cd /rsgrps1/mbsulli/sergei/ann_genomes/bowtie_runs/cat_all_vs/xmls/

/usr/bin/time -v ./piim_r -i $ARGV[0].xml -lk lk_n50_t0.01_allsmooth -joint-est >$ARGV[0].xml.piim.out 2>$ARGV[0].xml.piim.err\n";
