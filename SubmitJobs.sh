#!/bin/sh                                                                                                                                                                                    

source /afs/cern.ch/user/l/lcunquei/TestTrigger23/Setup_FileLocation.sh



mkdir -p Log
SubmissionFile=Step1.condor

echo "Universe   = vanilla" > $SubmissionFile
echo "Executable = `pwd`/StepProcess.sh" >> $SubmissionFile
echo "+JobFlavour = microcentury" >> $SubmissionFile
echo "should_transfer_files = NO" >> $SubmissionFile
# echo "transfer_output_files = DONE.txt" >> $SubmissionFile                                                                                                                                 
echo >> $SubmissionFile
Count=0
Input=$MyInputFiles
for i in $Input
do
   echo "Arguments = `pwd` $i /eos/home-l/lcunquei/OutputTriggerCorrect $Count" >> $SubmissionFile
   echo 'Output    = Log/Step1NoPU.out.$(Process)' >> $SubmissionFile
   echo 'Error     = Log/Step1NoPU.err.$(Process)' >> $SubmissionFile
   echo 'Log       = Log/Step1NoPU.log.$(Process)' >> $SubmissionFile
   echo 'Queue' >> $SubmissionFile
#   Count=`echo $Count | ./AddConst 1`
    Count=$((Count+1))
done



#condor_submit $SubmissionFile

