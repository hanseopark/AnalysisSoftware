#! /bin/bash

# variable: $1 = username, $2 = AOD/ESD, $3 = filename

# usernames: fbock, leardini, passfeld, amarin, mwilde, pgonzales
# specific location is added the same username will have a different output directory
# i.e fbockGSI or leardiniALICESERV1, passfeldMAF

if [ $1 = "fbock" ]; then 
   BASEDIR=/run/media/fredi/ff316d3c-0ff4-40e6-8b73-e06f46108e1b/Grid/OutputLegoTrains/pPb
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pPb
elif [ $1 = "bachelor" ]; then 
   BASEDIR=/home/admin1/leardini/GridOutput/pp  #/alidata50/alice_u/leardini/GridOutput/pp
elif [ $1 = "leardini" ]; then 
   BASEDIR=/home/admin1/leardini/GridOutput/PbPb
elif [ $1 = "leardiniALICESERV1" ]; then 
   BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
elif [ $1 = "leardiniGSI" ]; then 
   BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/pPb
elif [ $1 = "passfeld" ]; then 
   BASEDIR=~/work/Gridoutput/pPb
elif [ $1 = "passfeldMAF" ]; then 
   BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then  
   BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pPb
elif [ $1 = "amarin" ]; then     
   BASEDIR=/Users/marin/
elif [ $1 = "amarinGSI" ]; then     
   BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/pPb 
elif [ $1 = "amarinALICESERV1" ]; then     
   BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/   
elif [ $1 = "mwilde" ]; then        
   BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
   BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pPb 
elif [ $1 = "pgonzales" ]; then     
   BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
   BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pPb
fi

mkdir -p $BASEDIR

 
TRAINDIR=Legotrain-vAN-20141104
   LHC11hData=77_20141105-0956;  #78_20141105-1028;  #
   LHC14a1a=266_20141105-1002;   #267_20141105-1053;   #
   LHC14a1b=268_20141105-1003;   #269_20141105-1004;   # 
  
OUTPUTDIR=$BASEDIR/$TRAINDIR;
# if [ $3 = "AODdata" ]; then
    TRAINPATHData=GA_PbPb_AOD
# elif [ $4 = "AODmc" ]; then
    TRAINPATHMC=GA_PbPb_MC_AOD
# fi   
OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData
OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a
OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b
  
mkdir -p $OUTPUTDIR_LHC11h
mkdir -p $OUTPUTDIR_LHC14a1a
mkdir -p $OUTPUTDIR_LHC14a1b

fileName=$2  #GammaConvV1_QA_5010001_022092970003190000000.root
counter=0;

#if [ $3 = "AODdata" ]; then
#   echo "copying LHC11h data" 
#   #rm runNumbersLHC11hFailedRunMerge.txt
#   #category=`echo /home/admin1/leardini/photonconv/AnalysisSoftware/pastel.txt | cut -d "/" -f 7 | cut -d "." -f1`
#   #runNumbers=`cat /home/admin1/leardini/photonconv/AnalysisSoftware/pastel.txt`  #LHC11h_AOD145_full.txt` #LHC11h.txt`
#   runNumbers=`cat /home/admin1/leardini/photonconv/AnalysisSoftware/LHC11h_AOD145_goodruns.txt` #LHC11h.txt`
#   for runNumber in $runNumbers; do
#      echo $runNumber
#      mkdir -p $OUTPUTDIR_LHC11h/$runNumber
#      if [ -f $OUTPUTDIR_LHC11h/$runNumber/$fileName ]; then
#           echo "file " $fileName  " has already been copied for run " $runNumber
#      else 
#         alien_cp alien:/alice/data/2011/LHC11h_2/000$runNumber/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/$fileName file:$OUTPUTDIR_LHC11h/$runNumber/
#      fi
#      if [ -f $OUTPUTDIR_LHC11h/$runNumber/$fileName ]; then
#           echo "file " $fileName  " has already been copied for sucessfully for run " $runNumber
#      else 
#         echo $runNumber >> runNumbersLHC11hFailedRunMerge.txt
#      fi    
#   done;

#   runNumbersBroken=`cat runNumbersLHC11hFailedRunMerge.txt`
#   for runNumber in $runNumbersBroken; do
#      echo "copying stage_1 output for " $runNumber
#      mkdir -p $OUTPUTDIR_LHC11h/$runNumber
#      stageOutputs=`alien_ls /alice/data/2011/LHC11h_2/000$runNumber/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/Stage_1/`
#      for stageOutput in $stageOutputs; do
#         mkdir -p $OUTPUTDIR_LHC11h/$runNumber/Stage_1/$stageOutput
#         if [ -f $OUTPUTDIR_LHC11h/$runNumber/Stage_1/$stageOutput/$fileName ]; then
#           echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
#         else 
#          alien_cp alien:/alice/data/2011/LHC11h_2/000$runNumber/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC11h/$runNumber/Stage_1/$stageOutput/
#         fi
#      done;

#   done;



##################################################
#  for runNumber in $runNumbers; do
##   ls $OUTPUTDIR_LHC11h/$runNumber/GammaConvV1_*.root > fileLHC11hChangeStructure.txt
##   fileCutNumbers=`cat fileLHC11hChangeStructure.txt`
##   for fileName in $fileNumbers; do
#      echo $fileName
#      number=`echo $fileName | cut -d "_" -f 2 | cut -d "." -f1` 
#      echo $number
#      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/$runNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/$runNumber/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/$runNumber/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/$runNumber/CutSelection_LHC11h_$number.log\"\)
##   done;
#  done;
###################################################
#  for runNumber in $runNumbers; do
#      echo $runNumber
#      mkdir -p $OUTPUTDIR_LHC11h/$category
#      if [ "$(ls $OUTPUTDIR_LHC11h/$runNumber/GammaConv*)" ]; then
#        if [ $counter = 0 ]; then
#	   cp $OUTPUTDIR_LHC11h/$runNumber/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root $OUTPUTDIR_LHC11h/intermediate_$category.root 
#	   ls $OUTPUTDIR_LHC11h
#	   counter=$(($counter+1));
#	   echo $counter;
#        else
#	   hadd -T -f $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_110_$category.root $OUTPUTDIR_LHC11h/intermediate_$category.root $OUTPUTDIR_LHC11h/$runNumber/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root
#	   #root -l -b -x -q CheckBadRuns.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_110_$category.root\"\)
#	   mv $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_110_$category.root $OUTPUTDIR_LHC11h/intermediate_$category.root 
#        fi
#     else
#	echo "file not there"
#     fi
#  done;
#  mv $OUTPUTDIR_LHC11h/intermediate_$category.root  $OUTPUTDIR_LHC11h/$category/GammaConvV1_GA_PbPb_LHC11h-pass2_110_$category.root
#  echo "Final merge category " $category 
##  root -l -b -x -q CheckBadRuns.C\(\"$OUTPUTDIR_LHC11h/intermediate_$category.root\"\)	   
##  root -l -b -x -q CheckBadRuns.C\(\"$OUTPUTDIR_LHC11h/$category/GammaConvV1_GA_PbPb_LHC11h-pass2_110_$category.root\"\)


#fi


# if [ $4 = "AODmc" ]; then
#   echo "copying LHC14a1a data" 
#   #runNumbers=`alien_ls /alice/data/2011/LHC11h_2` 
#   rm runNumbersLHC11hFailedRunMergeLHC14a1a.txt
#   runNumbers=`cat /home/admin1/leardini/photonconv/AnalysisSoftware/LHC11h_AOD145_goodruns.txt`  #LHC11h_AOD145_full.txt` #LHC11h.txt`
#   for runNumber in $runNumbers; do
#      echo $runNumber
#      mkdir -p $OUTPUTDIR_LHC14a1a/$runNumber
#      if [ -f $OUTPUTDIR_LHC14a1a/$runNumber/$fileName ]; then
#           echo "file " $fileName  " has already been copied for run " $runNumber
#      else 
#	 alien_cp alien:/alice/sim/2014/LHC14a1a/$runNumber/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1a/$fileName file:$OUTPUTDIR_LHC14a1a/$runNumber
#      fi
#      if [ -f $OUTPUTDIR_LHC14a1a/$runNumber/$fileName ]; then
#           echo "file " $fileName  " has already been copied for sucessfully for run " $runNumber
#      else 
#         echo $runNumber >> runNumbersLHC11hFailedRunMergeLHC14a1a.txt
#      fi     
#   done;

#   runNumbersBroken=`cat runNumbersLHC11hFailedRunMergeLHC14a1a.txt`
#   for runNumber in $runNumbersBroken; do
#      echo "copying stage_1 output for " $runNumber
#      mkdir -p $OUTPUTDIR_LHC14a1a/$runNumber
#      stageOutputs=`alien_ls /alice/sim/2014/LHC14a1a/$runNumber/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1a/Stage_1/`
#      for stageOutput in $stageOutputs; do
#         mkdir -p $OUTPUTDIR_LHC14a1a/$runNumber/Stage_1/$stageOutput
#         if [ -f $OUTPUTDIR_LHC14a1a/$runNumber/Stage_1/$stageOutput/$fileName ]; then
#           echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
#         else 
#	  alien_cp alien:/alice/sim/2014/LHC14a1a/$runNumber/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1a/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC14a1a/$runNumber/Stage_1/$stageOutput/
#         fi
#      done;

#   done;

   echo "copying LHC14a1b data" 
   #runNumbers=`alien_ls /alice/data/2011/LHC11h_2` 
   rm runNumbersLHC11hFailedRunMergeLHC14a1b.txt
   runNumbers=`cat /home/admin1/leardini/photonconv/AnalysisSoftware/LHC11h_AOD145_goodruns.txt`  #LHC11h_AOD145_full.txt` #LHC11h.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      mkdir -p $OUTPUTDIR_LHC14a1b/$runNumber
      if [ -f $OUTPUTDIR_LHC14a1b/$runNumber/$fileName ]; then
           echo "file " $fileName  " has already been copied for run " $runNumber
      else 
	 alien_cp alien:/alice/sim/2014/LHC14a1b/$runNumber/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1b/$fileName file:$OUTPUTDIR_LHC14a1b/$runNumber
      fi
      if [ -f $OUTPUTDIR_LHC14a1b/$runNumber/$fileName ]; then
           echo "file " $fileName  " has already been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC11hFailedRunMergeLHC14a1b.txt
      fi     
   done;

   runNumbersBroken=`cat runNumbersLHC11hFailedRunMergeLHC14a1b.txt`
   for runNumber in $runNumbersBroken; do
      echo "copying stage_1 output for " $runNumber
      mkdir -p $OUTPUTDIR_LHC14a1b/$runNumber
      stageOutputs=`alien_ls /alice/sim/2014/LHC14a1b/$runNumber/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1b/Stage_1/`
      for stageOutput in $stageOutputs; do
         mkdir -p $OUTPUTDIR_LHC14a1b/$runNumber/Stage_1/$stageOutput
         if [ -f $OUTPUTDIR_LHC14a1b/$runNumber/Stage_1/$stageOutput/$fileName ]; then
           echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
         else 
	  alien_cp alien:/alice/sim/2014/LHC14a1b/$runNumber/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1b/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC14a1b/$runNumber/Stage_1/$stageOutput/
         fi
      done;

   done;

# fi

