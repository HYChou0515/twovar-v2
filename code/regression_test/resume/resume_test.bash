#!/bin/bash
DATA=../../../data/heart_scale
TRAIN=../../train
CUT_ITER=300
TTL_ITER=500
SMRY=12
ID=$(echo $@ | sed 's/ //g' | sed 's/-//g')
LOG=".log$ID"
RES=".res$ID"
CMP=".cmp$ID"
LOGCMP=".logcmp$ID"

if [ -f 1$RES ]
then
	rm 1$RES
fi
$TRAIN $@ -m $CUT_ITER $DATA 1$LOG 1$RES > 1$CMP
exit
$TRAIN $@ -m $TTL_ITER $DATA 2$LOG 1$RES > 2$CMP
$TRAIN $@ -m $TTL_ITER $DATA 3$LOG > 3$CMP

if [ "$(wc -l < 2$LOG)" -lt $(($SMRY+1)) ]; then
	>&2 echo "$ID 2$LOG empty"
	exit 1
fi
LOG2_ITER_START="$(head -1 2$LOG | cut -d' ' -f 2)"
head -n$(($LOG2_ITER_START-1)) 1$LOG > 1$LOG.tmp
cat 2$LOG >> 1$LOG.tmp
viewlog.bash -f 1,3- 1$LOG.tmp > 1$LOGCMP
viewlog.bash -f 1,3- 3$LOG > 3$LOGCMP
if [ -n "$(diff 1$LOGCMP 3$LOGCMP)" ]
then
	>&2 echo $ID
fi
