###Code Authors:  
[Mustafa Mustafa](http://github.com/MustafaMustafa)  - mmustafa@lbl.gov  

- - -
###How to build this code:  
```bash
mkdir myMcAnalysis
cd myMcAnalysis

# Replace address below with your own fork if you have one
git clone git@github.com:MustafaMustafa/auau200GeVRun14Ana.git

# Link all needed code under one StRoot directory:
mkdir StRoot
ln -s `pwd`/auau200GeVRun14Ana/StRoot/StMcAnalysisMaker StRoot

# Compile
starver SL15c
cons
```

###How to run this code:  
```bash
# For testing we can run the code on one file:
ln -s `pwd`/auau200GeVRun14Ana/StRoot/macros/run_StMcAnalysisMaker.C
root4star -l -b -q -x 'run_StMcAnalysisMaker.C("input.event.root","test_out")'
```
