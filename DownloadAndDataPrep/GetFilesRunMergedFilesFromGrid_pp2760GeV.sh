#! /bin/bash

# variable: $1 = username, $2 = AOD/ESD, $3 = filename

# usernames: fbock, leardini, passfeld, amarin, mwilde, pgonzales
# specific location is added the same username will have a different output directory
# i.e fbockGSI or leardiniALICESERV1, passfeldMAF

if [ $1 = "fbock" ]; then 
   BASEDIR=/home/fbock/Photon/Grid/OutputLegoTrains/pp/
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp
elif [ $1 = "leardini" ]; then 
   BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then 
   BASEDIR=/alidata50/alice_u/leardini/GridOutput/pp/
elif [ $1 = "leardiniGSI" ]; then 
   BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/pp
elif [ $1 = "passfeld" ]; then 
   BASEDIR=~/work/Gridoutput/pp
elif [ $1 = "passfeldMAF" ]; then 
   BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then  
   BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pp
elif [ $1 = "amarin" ]; then     
   BASEDIR=/Users/marin/
elif [ $1 = "amarinGSI" ]; then     
   BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/pp 
elif [ $1 = "amarinALICESERV1" ]; then     
   BASEDIR=/alidata50/alice_u/amarin/GridOutput/pp/   
elif [ $1 = "mwilde" ]; then        
   BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
   BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pp 
elif [ $1 = "pgonzales" ]; then     
   BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
   BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pp
fi

mkdir -p $BASEDIR
fileName=$3 
 

TRAINDIR=Legotrain-vAN-20150409-QA-ConvCalo
	LHC11a="496_20150411-1435"; #ESD
	LHC12f1a="452_20150410-1320" ; #ESD 
	LHC12f1b="467_20150411-2319" ; #ESD 
	LHC12i3="464_20150410-1335" ; #ESD 

OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $2 = "AOD" ]; then
   TRAINPATHData=GA_pp_AOD
else    
   TRAINPATHData=GA_pp
fi

if [ $2 = "AOD" ]; then
   TRAINPATHMC=GA_pp_MC_AOD
else    
   TRAINPATHMC=GA_pp_MC
fi

if [ $LHC12f1a != "" ]; then
   OUTPUTDIR_LHC12f1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12f1a
   echo $TRAINPATHMC-$LHC12f1a;
   mkdir -p $OUTPUTDIR_LHC12f1a
   rm runNumbersLHC12f1aFailedRunMerge.txt
   echo "copying LHC12f1a MC" 
   runNumbers=`cat runNumbersLHC11a_pass4_wSDD.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC12f1a/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC12f1a/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1a/$fileName"
         alien_cp alien:/alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1a/$fileName file:$OUTPUTDIR_LHC12f1a/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC12f1a/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC12f1aFailedRunMerge.txt
      fi     
   done;
fi

if [ $LHC12f1b != "" ]; then
   OUTPUTDIR_LHC12f1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12f1b
   echo $TRAINPATHMC-$LHC12f1b;
   mkdir -p $OUTPUTDIR_LHC12f1b
   rm runNumbersLHC12f1bFailedRunMerge.txt
   echo "copying LHC12f1b MC" 
   runNumbers=`cat runNumbersLHC11a_pass4_wSDD.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC12f1b/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC12f1b/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1b/$fileName"
         alien_cp alien:/alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1b/$fileName file:$OUTPUTDIR_LHC12f1b/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC12f1b/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC12f1bFailedRunMerge.txt
      fi     
   done;
fi

if [ $LHC12i3 != "" ]; then
   OUTPUTDIR_LHC12i3=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12i3
   echo $TRAINPATHMC-$LHC12i3;
   mkdir -p $OUTPUTDIR_LHC12i3
   rm runNumbersLHC12i3FailedRunMerge.txt
   echo "copying LHC12i3 MC" 
   runNumbers=`cat runNumbersLHC11a_pass4_wSDD.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC12i3/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC12i3/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2012/LHC12i3/$runNumber/PWGGA/$TRAINPATHMC/$LHC12i3/$fileName"
         alien_cp alien:/alice/sim/2012/LHC12i3/$runNumber/PWGGA/$TRAINPATHMC/$LHC12i3/$fileName file:$OUTPUTDIR_LHC12i3/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC12i3/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC12i3FailedRunMerge.txt
      fi     
   done;
fi


if [ $LHC11a != "" ]; then
   OUTPUTDIR_LHC11a=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11a
   echo $TRAINPATHData-$LHC11a;
   mkdir -p $OUTPUTDIR_LHC11a
   rm runNumbersLHC11aFailedRunMerge.txt
   echo "copying LHC11a" 
   runNumbers=`cat runNumbersLHC11a_pass4_wSDD.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      mkdir -p $OUTPUTDIR_LHC11a/$runNumber
      if [ -f $OUTPUTDIR_LHC11a/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then 
            echo "copying /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11a/$fileName"
            alien_cp alien:/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11a/$fileName file:$OUTPUTDIR_LHC11a/$runNumber
         else 
            echo "copying /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11a/$fileName"
            alien_cp alien:/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11a/$fileName file:$OUTPUTDIR_LHC11a/$runNumber
         fi
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC11a/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC11aFailedRunMerge.txt
      fi     
   done;
fi