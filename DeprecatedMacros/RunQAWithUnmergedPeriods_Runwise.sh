! /bin/bash

  echo "Please tell me the data path.";
  read datapath

   echo "Please tell me the QA fileName.";
  read fileNameQA

  echo "Please tell me the cutnumber."
  read cutnumber

  echo "Please tell me the QA cutnumber."
  read cutnumberQA
 
  echo "Please tell me the OutputName"
  read outputName
   
  echo "Is it MC? 0/1"
  read kMC

  echo "With which digits do the runs start?"
  read digits	
   
  echo $datapath
  echo $fileNameQA
  echo $cutnumber




   kFirstRun=1;

   runNumbers=`ls $datapath/ | grep $digits`
   for runNumber in $runNumbers; do   
	echo $runNumber
        inputFiles=`find $datapath/$runNumber -name $fileNameQA`
        nMaxfiles=`find $datapath/$runNumber -name $fileNameQA | wc -l`
        nCurrentFile=1
        echo "maximum number of file " $nMaxfiles
        for inputFile in $inputFiles; do
           echo "**********************************************************"
           echo $nCurrentFile/$nMaxfiles
           echo $inputFile
           echo "**********************************************************"
           if [ -f $inputFile ]; then
              nCurrentFile=`expr $nCurrentFile + 1`
              if [ $kFirstRun = 1 ]; then
                 root -l -b -q -x TaskV1/BuildHistogramsForGammaQAAdvV3.C\+\+\(\"$inputFile\"\,\"$cutnumber\"\,\"$cutnumberQA\"\,$kMC\,0\,0\,\"$datapath/$runNumber/$outputName\"\)
                 kFirstRun=0;
              else    
                 root -l -b -q -x TaskV1/BuildHistogramsForGammaQAAdvV3.C\+\+\(\"$inputFile\"\,\"$cutnumber\"\,\"$cutnumberQA\"\,$kMC\,1\,0\,\"$datapath/$runNumber/$outputName\"\)
                 kFirstRun=0;
              fi
           fi
        done;
   done;
   


