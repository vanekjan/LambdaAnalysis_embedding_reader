#!/bin/csh

starver SL21d

rm local.log

#root4star -l -b -q -x 'run_StMcAnalysisMaker.C("./input_local/test_false/pythia8.event.root","./output_local/output.root",0)' 
root4star -l -b -q -x 'run_StMcAnalysisMaker.C("./input_local/test_true/pythia8.event.root","./output_local/output.root",0)' 
#root4star -l -b -q -x 'run_StMcAnalysisMaker.C("/star/u/vanekjan/pwg/vanekjan/ppEmbedding/Lambda_PYTHIA_pp/production/bfc/Run12_L_Lbar_1k_full/out_6A64DE0C542F06FF80341054720D98C2_0.event.root","./output_local/output.root",0)' 
