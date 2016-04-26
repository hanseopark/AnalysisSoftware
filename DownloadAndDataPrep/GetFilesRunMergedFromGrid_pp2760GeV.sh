#! /bin/bash

# variable: $1 = username, $2 = AOD/ESD, $3 = filename

# usernames: fbock, leardini, passfeld, amarin, mwilde, pgonzales
# specific location is added the same username will have a different output directory
# i.e fbockGSI or leardiniALICESERV1, passfeldMAF

if [ $1 = "fbock" ]; then 
   BASEDIR=/home/fbock/Photon/Grid/OutputLegoTrains/pp
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp
elif [ $1 = "leardini" ]; then 
   BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then 
   BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
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
   BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/   
elif [ $1 = "mwilde" ]; then        
   BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
   BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pp 
elif [ $1 = "pgonzales" ]; then     
   BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
   BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pp
elif [ $1 = "dmuhlheim" ]; then 
   BASEDIR=/home/daniel/Desktop/Grid
fi

mkdir -p $BASEDIR
fileName=$3 
 
TRAINDIR=Legotrain-vAN-20140918-ConvCalo
LHC11aData="310_20140919-1322"; #ESD
LHC12f1aMC="174_20140919-1327";
LHC12f1bMC="175_20140919-1327";  

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

if [ $LHC12f1aMC != "" ]; then
   OUTPUTDIR_LHC12f1aMC=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12f1aMC
   echo $TRAINPATHMC-$LHC12f1aMC;
   mkdir -p $OUTPUTDIR_LHC12f1aMC
   rm runNumbersLHC12f1aFailedRunMerge.txt
   echo "copying LHC12f1a MC" 
   runNumbers=`cat runNumbersLHC12f1a.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC12f1aMC/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC12f1aMC/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1aMC/$fileName"
         alien_cp alien:/alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1aMC/$fileName file:$OUTPUTDIR_LHC12f1aMC/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC12f1a/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC12f1aFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC12f1a.txt`
   runNumbersBroken=`cat runNumbersLHC12f1aFailedRunMerge.txt`
   for runNumber in $runNumbersBroken; do
      echo "copying stage_1 output for " $runNumber
      mkdir -p $OUTPUTDIR_LHC12f1aMC/$runNumber
      stageOutputs=`alien_ls /alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1aMC/Stage_2/`
      for stageOutput in $stageOutputs; do
		 testOut=`alien_ls /alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1aMC/Stage_2/$stageOutput/$fileName`
		 if [ $testOut = $fileName ]; then 
			mkdir -p $OUTPUTDIR_LHC12f1aMC/$runNumber/Stage_2/$stageOutput
			if [ -f $OUTPUTDIR_LHC12f1aMC/$runNumber/Stage_2/$stageOutput/$fileName ]; then
				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
			else 
				echo "copying /alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1aMC/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/sim/2012/LHC12f1a/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1aMC/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC12f1aMC/$runNumber/Stage_2/$stageOutput/
			fi
		 fi	
      done;
   done;      
fi

if [ $LHC12f1bMC != "" ]; then
   OUTPUTDIR_LHC12f1bMC=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12f1bMC
   echo $TRAINPATHMC-$LHC12f1bMC;
   mkdir -p $OUTPUTDIR_LHC12f1bMC
   rm runNumbersLHC12f1bFailedRunMerge.txt
   echo "copying LHC12f1b MC" 
   runNumbers=`cat runNumbersLHC12f1b.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC12f1bMC/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC12f1bMC/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1bMC/$fileName"
         alien_cp alien:/alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1bMC/$fileName file:$OUTPUTDIR_LHC12f1bMC/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC12f1b/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC12f1bFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC12f1b.txt`
   runNumbersBroken=`cat runNumbersLHC12f1bFailedRunMerge.txt`
   for runNumber in $runNumbersBroken; do
      echo "copying stage_1 output for " $runNumber
      mkdir -p $OUTPUTDIR_LHC12f1bMC/$runNumber
      stageOutputs=`alien_ls /alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1bMC/Stage_2/`
      for stageOutput in $stageOutputs; do
		 testOut=`alien_ls /alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1bMC/Stage_2/$stageOutput/$fileName`
		 if [ $testOut = $fileName ]; then 
			mkdir -p $OUTPUTDIR_LHC12f1bMC/$runNumber/Stage_2/$stageOutput
			if [ -f $OUTPUTDIR_LHC12f1bMC/$runNumber/Stage_2/$stageOutput/$fileName ]; then
				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
			else 
				echo "copying /alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1bMC/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/sim/2012/LHC12f1b/$runNumber/PWGGA/$TRAINPATHMC/$LHC12f1bMC/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC12f1bMC/$runNumber/Stage_2/$stageOutput/
			fi
		 fi	
      done;
   done;      
fi

# 
# 
if [ $LHC11aData != "" ]; then
   OUTPUTDIR_LHC11a=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11aData
   echo $TRAINPATHData-$LHC11aData;
   mkdir -p $OUTPUTDIR_LHC11a
   rm runNumbersLHC11aFailedRunMerge.txt
   echo "copying LHC11a MC" 
   runNumbers=`cat runNumbersLHC11a.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      mkdir -p $OUTPUTDIR_LHC11a/$runNumber
      if [ -f $OUTPUTDIR_LHC11a/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then 
            echo "copying /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/AOD139/PWGGA/$TRAINPATHData/$LHC11aData/$fileName"
            alien_cp alien:/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/AOD139/PWGGA/$TRAINPATHData/$LHC11aData/$fileName file:$OUTPUTDIR_LHC11a/$runNumber
         else 
            echo "copying /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11aData/$fileName"
            alien_cp alien:/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11aData/$fileName file:$OUTPUTDIR_LHC11a/$runNumber
         fi
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC11a/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC11aFailedRunMerge.txt
      fi     
   done;

   runNumbersBroken=`cat runNumbersLHC11aFailedRunMerge.txt`
   for runNumber in $runNumbersBroken; do
      echo "copying stage_1 output for " $runNumber
      mkdir -p $OUTPUTDIR_LHC11a/$runNumber
      if [ $2 = "AOD" ]; then
         stageOutputs=`alien_ls /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/AOD139/PWGGA/$TRAINPATHData/$LHC11aData/Stage_2/`
      else 
         stageOutputs=`alien_ls /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11aData/Stage_2/`
      fi
      for stageOutput in $stageOutputs; do
		 testOut=`alien_ls /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11aData/Stage_2/$stageOutput/$fileName`
		 if [ $testOut = $fileName ]; then 
			mkdir -p $OUTPUTDIR_LHC11a/$runNumber/Stage_2/$stageOutput
			if [ -f $OUTPUTDIR_LHC11a/$runNumber/Stage_2/$stageOutput/$fileName ]; then
				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
			else 
				if [ $2 = "AOD" ]; then
				echo "copying /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/AOD139/PWGGA/$TRAINPATHData/$LHC11aData/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/AOD139/PWGGA/$TRAINPATHData/$LHC11aData/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC11a/$runNumber/Stage_2/$stageOutput/
				else
				echo "copying /alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11aData/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/$TRAINPATHData/$LHC11aData/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC11a/$runNumber/Stage_2/$stageOutput/
				fi
			fi
         fi
      done;
   done;      
fi

# if [ $LHC13c != "" ]; then
#    OUTPUTDIR_LHC13c=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC13c
#    echo $TRAINPATHData-$LHC13c;
#    mkdir -p $OUTPUTDIR_LHC13c
#    rm runNumbersLHC13cFailedRunMerge.txt
#    echo "copying LHC13c MC" 
#    runNumbers=`cat runNumbersLHC13c.txt`
#    for runNumber in $runNumbers; do
#       echo $runNumber
#       
#       
#       mkdir -p $OUTPUTDIR_LHC13c/$runNumber
#       if [ -f $OUTPUTDIR_LHC13c/$runNumber/$fileName ]; then
#          echo "file " $fileName  " has already been copied for run " $runNumber
#       else 
#          if [ $2 = "AOD" ]; then
#             echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/$fileName"
#             alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/$fileName file:$OUTPUTDIR_LHC13c/$runNumber
#          else
#             echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/$fileName"
#             alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/$fileName file:$OUTPUTDIR_LHC13c/$runNumber
#          fi   
#       fi
#       if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13c/$runNumber/$fileName ]; then
#          echo "file " $fileName  " been copied for sucessfully for run " $runNumber
#       else 
#          echo $runNumber >> runNumbersLHC13cFailedRunMerge.txt
#       fi     
#    done;
# 
#    runNumbersBroken=`cat runNumbersLHC13cFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC13c/$runNumber
#       if [ $2 = "AOD" ]; then
#          stageOutputs=`alien_ls /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/`
#       else 
#          stageOutputs=`alien_ls /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/`
#       fi
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				if [ $2 = "AOD" ]; then
# 				echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput/
# 				else 
# 				echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput/
# 				fi
# 			fi
# 		 fi	
#       done;
#    done;      
