#!/bin/bash
GRID=../../grid
DATA=../../../data/heart_scale

rm 1* 2* 3* 4* 5*
$GRID $DATA < test1.job
$GRID $DATA < test2.job
$GRID $DATA < test3.job
parallel --jobs 8 < test.job 2> err
RED='\033[0;31m'
NC='\033[0m'
while read p; do
	printf "${RED}${p}${NC}\n"
done <err
