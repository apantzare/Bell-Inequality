#!/bin/bash
## small script to run the analysis
analysis="main_HWWAnalysis"

##OPTION
echo Which option should I run? 
echo Options are:
echo 0 = run all data and MC one after another
echo 10 = run the big data file
echo 1= run the merged simulation file
echo 30= run the merged simulation file
echo 31= run the merged background file
echo 32= run the merged Higgs file
echo 2,3,4,5,6 = run MC samples only \(can be run in parallel\)
read varname
echo Option is $varname
option=$varname

#echo Should I use PROOF? \(will make things faster\)
#echo Options are:
#echo 0 = NO
#echo 1 = YES
#read proofvarname
#echo PROOF option is $proofvarname
parallel=0 #$proofvarname


## execute and run ROOT
echo "starting ROOT"
##
root -l -b << EOF
.L $analysis.C+
$analysis($parallel,$option)
EOF
##
echo "end of ROOT execution"
