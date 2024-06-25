#! /bin/sh
ID=`git log -1 --format='%cd %h' --date=short`
VERFILE=version.h
echo "        character*32 VERSIONID" > $VERFILE
echo "        parameter (VERSIONID='$ID')" >> $VERFILE
