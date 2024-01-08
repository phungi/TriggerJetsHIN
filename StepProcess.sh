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




root $WorkDir/Plot_Final_datacalo_JEC.C\(\"$File\"\)
 
cp out.root ${Output}/OutHist_PbPb22_data_corr_030_HFcut${ID}.root

rm out.root
