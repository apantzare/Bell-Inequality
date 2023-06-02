#!/bin/bash
########################################
#
# General plotting script
#
########################################

# names of 12 analyses
analysisCollection=("HWWAnalysis")

# location of their outputs
outputpath=("../Analysis/HWWAnalysis/Output_HWWAnalysis")

# begin
    analysisPath=${outputpath[${choice}]}

    # run main command
    root -l -b << EOF
	.L Plotting.cxx+
    	Plotting *m=new Plotting();
	m->SetLumi(10064);
	m->SetOption("HWWAnalysis");
    	m->SetInputLocation("$analysisPath")
    	m->run()
    	.q
EOF
########################################
