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
elif [ $1 = "nschmidt" ]; then        
   BASEDIR=/media/nschmidt/Daten/GridOutput/pp
fi

mkdir -p $BASEDIR
fileName=$3 
 

TRAINDIR=Legotrain-vAN-20150331-QA-pass4
	LHC10b="471_20150406-0640"; #ESD
	LHC10c="472_20150406-0642"; #ESD
	LHC10d="473_20150406-0644"; #ESD
	LHC10e="474_20150406-0642"; #ESD
	LHC10f="475_20150406-0642"; #ESD
	LHC10g="476_20150406-0641"; #ESD
	LHC14j4b="429_20150401-1155" ; #ESD 
	LHC14j4c="430_20150401-1156" ; #ESD 
	LHC14j4d="431_20150401-1158" ; #ESD 
	LHC14j4e="432_20150401-1157" ; #ESD 
	LHC14j4f="433_20150401-1157" ; #ESD 
 

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

if [ $LHC14j4b != "" ]; then
   OUTPUTDIR_LHC14j4b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14j4b
   echo $TRAINPATHMC-$LHC14j4b;
   mkdir -p $OUTPUTDIR_LHC14j4b
   rm runNumbersLHC14j4bFailedRunMerge.txt
   echo "copying LHC14j4b MC" 
   runNumbers=`cat runNumbersLHC14j4b.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC14j4b/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC14j4b/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2014/LHC14j4b/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4b/$fileName"
         alien_cp alien:/alice/sim/2014/LHC14j4b/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4b/$fileName file:$OUTPUTDIR_LHC14j4b/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC14j4b/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC14j4bFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC14j4b.txt`
#    runNumbersBroken=`cat runNumbersLHC14j4bFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC14j4b/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2014/LHC14j4b/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4b/Stage_2/`
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/sim/2014/LHC14j4b/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4b/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC14j4b/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC14j4b/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				echo "copying /alice/sim/2014/LHC14j4b/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4b/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/sim/2014/LHC14j4b/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4b/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC14j4b/$runNumber/Stage_2/$stageOutput/
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC14j4c != "" ]; then
   OUTPUTDIR_LHC14j4c=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14j4c
   echo $TRAINPATHMC-$LHC14j4c;
   mkdir -p $OUTPUTDIR_LHC14j4c
   rm runNumbersLHC14j4cFailedRunMerge.txt
   echo "copying LHC14j4c MC" 
   runNumbers=`cat runNumbersLHC14j4c.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC14j4c/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC14j4c/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2014/LHC14j4c/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4c/$fileName"
         alien_cp alien:/alice/sim/2014/LHC14j4c/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4c/$fileName file:$OUTPUTDIR_LHC14j4c/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC14j4c/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC14j4cFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC14j4c.txt`
#    runNumbersBroken=`cat runNumbersLHC14j4cFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC14j4c/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2014/LHC14j4c/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4c/Stage_2/`
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/sim/2014/LHC14j4c/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4c/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC14j4c/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC14j4c/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				echo "copying /alice/sim/2014/LHC14j4c/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4c/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/sim/2014/LHC14j4c/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4c/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC14j4c/$runNumber/Stage_2/$stageOutput/
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC14j4d != "" ]; then
   OUTPUTDIR_LHC14j4d=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14j4d
   echo $TRAINPATHMC-$LHC14j4d;
   mkdir -p $OUTPUTDIR_LHC14j4d
   rm runNumbersLHC14j4dFailedRunMerge.txt
   echo "copying LHC14j4d MC" 
   runNumbers=`cat runNumbersLHC14j4d.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC14j4d/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC14j4d/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2014/LHC14j4d/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4d/$fileName"
         alien_cp alien:/alice/sim/2014/LHC14j4d/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4d/$fileName file:$OUTPUTDIR_LHC14j4d/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC14j4d/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC14j4dFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC14j4d.txt`
#    runNumbersBroken=`cat runNumbersLHC14j4dFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC14j4d/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2014/LHC14j4d/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4d/Stage_2/`
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/sim/2014/LHC14j4d/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4d/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC14j4d/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC14j4d/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				echo "copying /alice/sim/2014/LHC14j4d/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4d/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/sim/2014/LHC14j4d/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4d/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC14j4d/$runNumber/Stage_2/$stageOutput/
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC14j4e != "" ]; then
   OUTPUTDIR_LHC14j4e=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14j4e
   echo $TRAINPATHMC-$LHC14j4e;
   mkdir -p $OUTPUTDIR_LHC14j4e
   rm runNumbersLHC14j4eFailedRunMerge.txt
   echo "copying LHC14j4e MC" 
   runNumbers=`cat runNumbersLHC14j4e.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC14j4e/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC14j4e/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2014/LHC14j4e/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4e/$fileName"
         alien_cp alien:/alice/sim/2014/LHC14j4e/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4e/$fileName file:$OUTPUTDIR_LHC14j4e/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC14j4e/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC14j4eFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC14j4e.txt`
#    runNumbersBroken=`cat runNumbersLHC14j4eFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC14j4e/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2014/LHC14j4e/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4e/Stage_2/`
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/sim/2014/LHC14j4e/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4e/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC14j4e/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC14j4e/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				echo "copying /alice/sim/2014/LHC14j4e/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4e/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/sim/2014/LHC14j4e/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4e/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC14j4e/$runNumber/Stage_2/$stageOutput/
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC14j4f != "" ]; then
   OUTPUTDIR_LHC14j4f=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14j4f
   echo $TRAINPATHMC-$LHC14j4f;
   mkdir -p $OUTPUTDIR_LHC14j4f
   rm runNumbersLHC14j4fFailedRunMerge.txt
   echo "copying LHC14j4f MC" 
   runNumbers=`cat runNumbersLHC14j4f.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC14j4f/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC14j4f/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2014/LHC14j4f/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4f/$fileName"
         alien_cp alien:/alice/sim/2014/LHC14j4f/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4f/$fileName file:$OUTPUTDIR_LHC14j4f/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC14j4f/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC14j4fFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC14j4f.txt`
#    runNumbersBroken=`cat runNumbersLHC14j4fFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC14j4f/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2014/LHC14j4f/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4f/Stage_2/`
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/sim/2014/LHC14j4f/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4f/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC14j4f/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC14j4f/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				echo "copying /alice/sim/2014/LHC14j4f/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4f/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/sim/2014/LHC14j4f/$runNumber/PWGGA/$TRAINPATHMC/$LHC14j4f/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC14j4f/$runNumber/Stage_2/$stageOutput/
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC10b != "" ]; then
   OUTPUTDIR_LHC10b=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10b
   echo $TRAINPATHData-$LHC10b;
   mkdir -p $OUTPUTDIR_LHC10b
   rm runNumbersLHC10bFailedRunMerge.txt
   echo "copying LHC10b" 
   runNumbers=`cat runNumbersLHC10b_pass4.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      mkdir -p $OUTPUTDIR_LHC10b/$runNumber
      if [ -f $OUTPUTDIR_LHC10b/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then 
            echo "copying /alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/$fileName"
            alien_cp alien:/alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/$fileName file:$OUTPUTDIR_LHC10b/$runNumber
         else 
            echo "copying /alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/$fileName"
            alien_cp alien:/alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/$fileName file:$OUTPUTDIR_LHC10b/$runNumber
         fi
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC10b/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC10bFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC10bFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC10b/$runNumber
#       if [ $2 = "AOD" ]; then
#          stageOutputs=`alien_ls /alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/Stage_2/`
#       else 
#          stageOutputs=`alien_ls /alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/Stage_2/`
#       fi
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC10b/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC10b/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				if [ $2 = "AOD" ]; then
# 				echo "copying /alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10b/$runNumber/Stage_2/$stageOutput/
# 				else
# 				echo "copying /alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10b/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10b/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10b/$runNumber/Stage_2/$stageOutput/
# 				fi
# 			fi
#          fi
#       done;
#    done;      
fi

if [ $LHC10c != "" ]; then
   OUTPUTDIR_LHC10c=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10c
   echo $TRAINPATHData-$LHC10c;
   mkdir -p $OUTPUTDIR_LHC10c
   rm runNumbersLHC10cFailedRunMerge.txt
   echo "copying LHC10c" 
   runNumbers=`cat runNumbersLHC10c_pass4.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      
      mkdir -p $OUTPUTDIR_LHC10c/$runNumber
      if [ -f $OUTPUTDIR_LHC10c/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then
            echo "copying /alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/$fileName"
            alien_cp alien:/alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/$fileName file:$OUTPUTDIR_LHC10c/$runNumber
         else
            echo "copying /alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/$fileName"
            alien_cp alien:/alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/$fileName file:$OUTPUTDIR_LHC10c/$runNumber
         fi   
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC10c/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC10cFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC10cFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC10c/$runNumber
#       if [ $2 = "AOD" ]; then
#          stageOutputs=`alien_ls /alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/Stage_2/`
#       else 
#          stageOutputs=`alien_ls /alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/Stage_2/`
#       fi
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC10c/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC10c/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				if [ $2 = "AOD" ]; then
# 				echo "copying /alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10c/$runNumber/Stage_2/$stageOutput/
# 				else 
# 				echo "copying /alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10c/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10c/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10c/$runNumber/Stage_2/$stageOutput/
# 				fi
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC10d != "" ]; then
   OUTPUTDIR_LHC10d=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10d
   echo $TRAINPATHData-$LHC10d;
   mkdir -p $OUTPUTDIR_LHC10d
   rm runNumbersLHC10dFailedRunMerge.txt
   echo "copying LHC10d" 
   runNumbers=`cat runNumbersLHC10d_pass4.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      
      mkdir -p $OUTPUTDIR_LHC10d/$runNumber
      if [ -f $OUTPUTDIR_LHC10d/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then
            echo "copying /alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/$fileName"
            alien_cp alien:/alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/$fileName file:$OUTPUTDIR_LHC10d/$runNumber
         else
            echo "copying /alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/$fileName"
            alien_cp alien:/alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/$fileName file:$OUTPUTDIR_LHC10d/$runNumber
         fi   
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC10d/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC10dFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC10dFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC10d/$runNumber
#       if [ $2 = "AOD" ]; then
#          stageOutputs=`alien_ls /alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/Stage_2/`
#       else 
#          stageOutputs=`alien_ls /alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/Stage_2/`
#       fi
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC10d/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC10d/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				if [ $2 = "AOD" ]; then
# 				echo "copying /alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10d/$runNumber/Stage_2/$stageOutput/
# 				else 
# 				echo "copying /alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10d/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10d/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10d/$runNumber/Stage_2/$stageOutput/
# 				fi
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC10e != "" ]; then
   OUTPUTDIR_LHC10e=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10e
   echo $TRAINPATHData-$LHC10e;
   mkdir -p $OUTPUTDIR_LHC10e
   rm runNumbersLHC10eFailedRunMerge.txt
   echo "copying LHC10e" 
   runNumbers=`cat runNumbersLHC10e_pass4.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      
      mkdir -p $OUTPUTDIR_LHC10e/$runNumber
      if [ -f $OUTPUTDIR_LHC10e/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then
            echo "copying /alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/$fileName"
            alien_cp alien:/alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/$fileName file:$OUTPUTDIR_LHC10e/$runNumber
         else
            echo "copying /alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/$fileName"
            alien_cp alien:/alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/$fileName file:$OUTPUTDIR_LHC10e/$runNumber
         fi   
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC10e/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC10eFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC10eFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC10e/$runNumber
#       if [ $2 = "AOD" ]; then
#          stageOutputs=`alien_ls /alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/Stage_2/`
#       else 
#          stageOutputs=`alien_ls /alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/Stage_2/`
#       fi
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC10e/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC10e/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				if [ $2 = "AOD" ]; then
# 				echo "copying /alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10e/$runNumber/Stage_2/$stageOutput/
# 				else 
# 				echo "copying /alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10e/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10e/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10e/$runNumber/Stage_2/$stageOutput/
# 				fi
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC10f != "" ]; then
   OUTPUTDIR_LHC10f=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10f
   echo $TRAINPATHData-$LHC10f;
   mkdir -p $OUTPUTDIR_LHC10f
   rm runNumbersLHC10fFailedRunMerge.txt
   echo "copying LHC10f" 
   runNumbers=`cat runNumbersLHC10f_pass4.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      
      mkdir -p $OUTPUTDIR_LHC10f/$runNumber
      if [ -f $OUTPUTDIR_LHC10f/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then
            echo "copying /alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/$fileName"
            alien_cp alien:/alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/$fileName file:$OUTPUTDIR_LHC10f/$runNumber
         else
            echo "copying /alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/$fileName"
            alien_cp alien:/alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/$fileName file:$OUTPUTDIR_LHC10f/$runNumber
         fi   
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC10f/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC10fFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC10fFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC10f/$runNumber
#       if [ $2 = "AOD" ]; then
#          stageOutputs=`alien_ls /alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/Stage_2/`
#       else 
#          stageOutputs=`alien_ls /alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/Stage_2/`
#       fi
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC10f/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC10f/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				if [ $2 = "AOD" ]; then
# 				echo "copying /alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10f/$runNumber/Stage_2/$stageOutput/
# 				else 
# 				echo "copying /alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10f/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10f/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10f/$runNumber/Stage_2/$stageOutput/
# 				fi
# 			fi
# 		 fi	
#       done;
#    done;      
fi

if [ $LHC10g != "" ]; then
   OUTPUTDIR_LHC10g=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10g
   echo $TRAINPATHData-$LHC10g;
   mkdir -p $OUTPUTDIR_LHC10g
   rm runNumbersLHC10gFailedRunMerge.txt
   echo "copying LHC10g" 
   runNumbers=`cat runNumbersLHC10g_pass4.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      
      mkdir -p $OUTPUTDIR_LHC10g/$runNumber
      if [ -f $OUTPUTDIR_LHC10g/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then
            echo "copying /alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/$fileName"
            alien_cp alien:/alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/$fileName file:$OUTPUTDIR_LHC10g/$runNumber
         else
            echo "copying /alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/$fileName"
            alien_cp alien:/alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/$fileName file:$OUTPUTDIR_LHC10g/$runNumber
         fi   
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC10g/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC10gFailedRunMerge.txt
      fi     
   done;

#    runNumbersBroken=`cat runNumbersLHC10gFailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC10g/$runNumber
#       if [ $2 = "AOD" ]; then
#          stageOutputs=`alien_ls /alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/Stage_2/`
#       else 
#          stageOutputs=`alien_ls /alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/Stage_2/`
#       fi
#       for stageOutput in $stageOutputs; do
# 		 testOut=`alien_ls /alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/Stage_2/$stageOutput/$fileName`
# 		 if [ $testOut = $fileName ]; then 
# 			mkdir -p $OUTPUTDIR_LHC10g/$runNumber/Stage_2/$stageOutput
# 			if [ -f $OUTPUTDIR_LHC10g/$runNumber/Stage_2/$stageOutput/$fileName ]; then
# 				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
# 			else 
# 				if [ $2 = "AOD" ]; then
# 				echo "copying /alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10g/$runNumber/Stage_2/$stageOutput/
# 				else 
# 				echo "copying /alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/Stage_2/$stageOutput/$fileName"
# 				alien_cp alien:/alice/data/2010/LHC10g/000$runNumber/pass4/PWGGA/$TRAINPATHData/$LHC10g/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC10g/$runNumber/Stage_2/$stageOutput/
# 				fi
# 			fi
# 		 fi	
#       done;
#    done;      
fi
