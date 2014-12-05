#!/bin/bash

DIR=$1

cd $DIR

for file in *.out; do
	fil=${file%.ncm.out}
	fil2=${file%.out}
        inf1=$fil2.pi
        inf2=$fil2.frq
#	echo "Rscript /rsgrps1/mbsulli/sergei/ann_genomes/bin/scripts/hist.R $DIR/$inf1 $DIR/$inf2 $fil $fil.pdf"
	RHIST=`Rscript /rsgrps1/mbsulli/sergei/ann_genomes/bin/scripts/hist.R $DIR/$inf1 $DIR/$inf2 $fil $fil.pdf`
done

