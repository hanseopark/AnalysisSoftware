#! /bin/bash

   echo "Please tell me the data path.";
   read datapath

   echo "Please tell me the output folder.";
   read outputFolder

   echo "Please tell me the QA fileName.";
   read fileNameQA

   echo "Please tell me the Material fileName.";
   read fileNameMaterial

   echo "Please tell me the cutnumber."
   read cutnumber

   echo $fileNameQA
   echo $fileNameMaterial
   
   kFirstRun=1;
  ls $datapath/*/$fileNameQA > listOfRuns.txt
#    fileNames=`ls $datapath/`
   fileNames=`cat listOfRuns.txt`
   nMaxfiles=`cat listOfRuns.txt | wc -l`
   nCurrentFile=1

   mergedOutputNR=1
   echo "maximum number of file " $nMaxfiles
   
   until [ $nMaxfiles -lt 0 ]; do
      echo $nMaxfiles                   
      Run1=$(head -n 1 listOfRuns.txt)
      sed -i 1d listOfRuns.txt
      Run2=$(head -n 1 listOfRuns.txt)
      sed -i 1d listOfRuns.txt
      Run3=$(head -n 1 listOfRuns.txt)
      sed -i 1d listOfRuns.txt
      Run4=$(head -n 1 listOfRuns.txt)
      sed -i 1d listOfRuns.txt
      Run5=$(head -n 1 listOfRuns.txt)
      sed -i 1d listOfRuns.txt
      
      echo $Run1 $Run2 $Run3 $Run4 $Run5
      mkdir -p $datapath\_merged/$outputFolder$mergedOutputNR
      hadd -f $datapath\_merged/$outputFolder$mergedOutputNR/$fileNameQA $Run1 $Run2 $Run3 $Run4 $Run5
      
      nMaxfiles=`expr $nMaxfiles - 5`
      mergedOutputNR=`expr $mergedOutputNR + 1`
   done


#    for i in `seq 1 10`; do
#      
#      echo $i
#    done    
   
#    for fileName in $fileNames; do
#      
#    
#    done;

   