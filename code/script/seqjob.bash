#!/bin/bash

JOBS="jobs"
JOBSTMP="$JOBS.tmp"
CMD="cmd"

if [[ -e $CMD ]]; then
	echo "cmd file $CMD already exists!"
	exit -1
fi

if [[ -e $JOBSTMP ]]; then
	echo "jobstmp file $JOBSTMP already exists"
	exit -1
fi

while [[ -s $JOBS ]] ; do
    head -1 $JOBS > $CMD
    bash $CMD 
    line_left=$(($(wc -l <$JOBS)-1))
    tail -$line_left $JOBS > $JOBS.tmp
    mv $JOBS.tmp $JOBS
done

rm -rf $CMD
rm -rf $JOBSTMP
