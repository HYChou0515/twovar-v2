#!/bin/bash

NEW_HASH=$(git branch | grep \* | cut -d ' ' -f2)
if [ "detech from" == *"$NEW_HASH"* ]; then
	NEW_HASH=$(git rev-parse HEAD)
fi
if [ -z $1 ]; then
	OLD_HASH=$(git rev-parse HEAD^)
else
	OLD_HASH=$(git rev-parse $1)
fi

rm -ri new
rm -ri old

mkdir -p new
mkdir -p old

WORK_DIR=../../

git checkout $OLD_HASH
make -C $WORK_DIR clean
cp -r $WORK_DIR/* old/

git checkout $NEW_HASH
make -C $WORK_DIR clean
cp -r $WORK_DIR/* new/