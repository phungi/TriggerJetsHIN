#!/bin/sh --login                                                                                                                                

#BSUB -q 1nh                                                                                                                                     

WorkDir=$1
File=$2
Output=$3
ID=$4




source $WorkDir/Setup_CMSSW.sh

echo Input files are: $File
echo WorkDir is: $WorkDir
echo output dir is: $Output

root $WorkDir/Molly_script.cpp\(\"$File\"\)
 
cp out.root ${Output}/${ID}.root

rm out.root