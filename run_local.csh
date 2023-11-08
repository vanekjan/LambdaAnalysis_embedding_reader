#!/bin/csh

starver SL21d

rm local.log

root4star -l -b -q -x 'run_StMcAnalysisMaker.C("./input_local/pythia8.event.root","./output_local/output.root",0)' 
