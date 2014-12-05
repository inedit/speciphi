#!/bin/bash

DIR=$1

cd $DIR

COUNT=0
for file in `ls *.gc`; do
	ARRAY[$COUNT]=$file
	COUNT=`expr $COUNT + 1`
done


for file in `ls *.gc`; do
	fil=${file%.gc}
	for file2 in "${ARRAY[@]}"; do
		fil2=${file2%.gc}
#		echo "Rscript $DIR/scripts/hist.R $DIR/$file $DIR/$file2 $fil.vs.$fil2 $fil.vs.$fil2.pdf"
		RHIST=`Rscript $DIR/hist.R $DIR/$file $DIR/$file2 $fil.vs.$fil2 $fil.vs.$fil2.pdf`
	done
done

