#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# $1 = option dataset selection

BASEDIR=/mnt/External2GBHarddrive/Grid/OutputLegoTrains/PbPb

TRAINDIR=Legotrain-v5-05-49-AN;
TRAINDIR=Legotrain-v5-05-49-AN;
# LHC11hData=33_20131217-0038;
# LHC11hData=34_20131217-1839;
# LHC11hData=35_20131217-2316;
# LHC11hData=36_20131217-2234;
# LHC11hData=37_20131217-2235;

LHC11hMC1=92_20131218-1816; #LHC12a17a_fix
# LHC11hMC1=101_20131219-2057; #LHC12a17a_fix

# LHC11hMC2=102_20131223-1132; #LHC12a17b_fix
# LHC11hMC2=106_20131223-1153; #LHC12a17b_fix
LHC11hMC2=110_20131223-1204; #LHC12a17b_fix
# LHC11hMC2=114_20131223-1234; #LHC12a17b_fix

LHC11hMC3=118_20131223-1257; #LHC12a17c_fix
# LHC11hMC3=122_20131223-1259; #LHC12a17c_fix

# LHC11hMC4=93_20131218-1803; #LHC12a17d_fix
# LHC11hMC4=100_20131219-2129; #LHC12a17d_fix
 
# LHC11hMC5=103_20131223-1135; #LHC12a17e_fix
# LHC11hMC5=107_20131223-1155; #LHC12a17e_fix
LHC11hMC5=111_20131223-1205; #LHC12a17e_fix
# LHC11hMC5=115_20131223-1235; #LHC12a17e_fix

LHC11hMC6=119_20131223-1257; #LHC12a17f_fix
# LHC11hMC6=123_20131223-1300; #LHC12a17f_fix

LHC11hMC7=94_20131218-1809; #LHC12a17g_fix
# LHC11hMC7=98_20131219-2057; #LHC12a17g_fix

# LHC11hMC8=104_20131223-1146; #LHC12a17h_fix
# LHC11hMC8=108_20131223-1156; #LHC12a17h_fix
LHC11hMC8=112_20131223-1206; #LHC12a17h_fix
# LHC11hMC8=116_20131223-1236; #LHC12a17h_fix

LHC11hMC9=120_20131223-1306; #LHC12a17i_fix
# LHC11hMC9=124_20131223-1306; #LHC12a17i_fix

LHC11hMC10=95_20131218-1810; #LHC13c7a_fix
# LHC11hMC10=99_20131219-2057; #LHC13c7a_fix

# LHC11hMC11=105_20131223-1147; #LHC13c7b_fix
# LHC11hMC11=109_20131223-1157; #LHC13c7b_fix
LHC11hMC11=113_20131223-1216; #LHC13c7b_fix
# LHC11hMC11=117_20131223-1236; #LHC13c7b_fix

LHC11hMC12=121_20131223-1258; #LHC13c7c_fix
# LHC11hMC12=125_20131223-1301; #LHC13c7c_fix


OUTPUTDIR=$BASEDIR/$TRAINDIR
#/alidata50/alice_u/leardini/OutputFiles/Legotrain-v5-05-42-AN/GA_PbPb-19_20131202-1122
OUTPUTDIR_Data=$BASEDIR/$TRAINDIR/GA_PbPb_AOD-$LHC11hData
OUTPUTDIR_MC1=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC1
OUTPUTDIR_MC2=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC2
OUTPUTDIR_MC3=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC3
OUTPUTDIR_MC4=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC4
OUTPUTDIR_MC5=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC5
OUTPUTDIR_MC6=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC6
OUTPUTDIR_MC7=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC7
OUTPUTDIR_MC8=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC8
OUTPUTDIR_MC9=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC9
OUTPUTDIR_MC10=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC10
OUTPUTDIR_MC11=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC11
OUTPUTDIR_MC12=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC11hMC12

if [ $1 = "data" ]; then
	mkdir -p $OUTPUTDIR_Data
	alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/$LHC11hData/merge/Gam* > filesToCopy.txt
	files=`cat filesToCopy.txt`
	for fileToCopy in $files; do
      if [ -f $OUTPUTDIR_Data/$fileToCopy ]; then
         echo "already copied $OUTPUTDIR_Data/$fileToCopy ";
      else 
         alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/$LHC11hData/merge/$fileToCopy file:$OUTPUTDIR_Data/
      fi
	done
	
	ls $OUTPUTDIR_Data/GammaConvV1_*.root > fileData.txt
        fileNumbers=`cat fileData.txt`
        for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_Data/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_Data_$number.root\"\,\"GammaConvV1_$number\"\)
      	root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_Data_$LHC11hData_$number.root\"\,\"$OUTPUTDIR/CutSelection_Data_AOD_$number.log\"\)
	done
	
elif [ $1 = "dataRuns" ]; then
   mkdir -p $OUTPUTDIR_Data/
   runNumbers=`cat listOfRunsLHC11h.txt`
   for runNumber in $runNumbers; do
      echo $runNumber 
      alien_ls /alice/data/2011/LHC11h_2/000$runNumber/ESDs/pass2/AOD142/PWGGA/GA_PbPb_AOD/$LHC11hData/Gam* > filesToCopy.txt
      mkdir -p $OUTPUTDIR_Data/$runNumber
      files=`cat filesToCopy.txt`
      for fileToCopy in $files; do
         if [ -f $OUTPUTDIR_Data/$runNumber/$fileToCopy ]; then
            echo "already copied $OUTPUTDIR_Data/$runNumber/$fileToCopy ";
         else 
            alien_cp alien:/alice/data/2011/LHC11h_2/000$runNumber/ESDs/pass2/AOD142/PWGGA/GA_PbPb_AOD/$LHC11hData/$fileToCopy file:$OUTPUTDIR_Data/$runNumber/ 
         fi
      done
   done

elif [ $1 = "MC" ]; then
	echo "MC no runwise";
	if [ $2 = "runwise" ]; then
      mkdir -p $OUTPUTDIR_MC1
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17a_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC1/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC1/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC1/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC1/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17a_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC1/$fileToCopy file:$OUTPUTDIR_MC1/$runNumber/ 
            fi
         done
      done

      mkdir -p $OUTPUTDIR_MC2
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17b_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC2/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC2/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC2/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC2/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17b_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC2/$fileToCopy file:$OUTPUTDIR_MC2/$runNumber/ 
            fi
         done
      done

      mkdir -p $OUTPUTDIR_MC3
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17c_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC3/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC3/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC3/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC3/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17c_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC3/$fileToCopy file:$OUTPUTDIR_MC3/$runNumber/ 
            fi
         done
      done

      mkdir -p $OUTPUTDIR_MC4
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17d_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC4/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC4/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC4/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC4/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17d_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC4/$fileToCopy file:$OUTPUTDIR_MC4/$runNumber/ 
            fi
         done
      done
      
      mkdir -p $OUTPUTDIR_MC5
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17e_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC5/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC5/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC5/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC5/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17e_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC5/$fileToCopy file:$OUTPUTDIR_MC5/$runNumber/ 
            fi
         done
      done
      
      mkdir -p $OUTPUTDIR_MC6
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17f_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC6/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC6/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC6/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC6/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17f_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC6/$fileToCopy file:$OUTPUTDIR_MC6/$runNumber/ 
            fi
         done
      done
      
      mkdir -p $OUTPUTDIR_MC7
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17g_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC7/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC7/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC7/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC7/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17g_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC7/$fileToCopy file:$OUTPUTDIR_MC7/$runNumber/ 
            fi
         done
      done
      
      mkdir -p $OUTPUTDIR_MC8
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17h_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC8/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC8/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC8/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC8/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17h_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC8/$fileToCopy file:$OUTPUTDIR_MC8/$runNumber/ 
            fi
         done
      done
      
      mkdir -p $OUTPUTDIR_MC9
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2012/LHC12a17i_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC9/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC9/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC9/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC9/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2012/LHC12a17i_fix/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC9/$fileToCopy file:$OUTPUTDIR_MC9/$runNumber/ 
            fi
         done
      done

      mkdir -p $OUTPUTDIR_MC10
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2013/LHC13c7a/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC10/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC10/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC10/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC10/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2013/LHC13c7a/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC10/$fileToCopy file:$OUTPUTDIR_MC10/$runNumber/ 
            fi
         done
      done
      
      mkdir -p $OUTPUTDIR_MC11
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2013/LHC13c7b/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC11/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC11/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC11/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC11/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2013/LHC13c7b/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC11/$fileToCopy file:$OUTPUTDIR_MC11/$runNumber/ 
            fi
         done
      done

      mkdir -p $OUTPUTDIR_MC12
      runNumbers=`cat listOfRunsLHC11h.txt`
      for runNumber in $runNumbers; do
         echo $runNumber 
         alien_ls /alice/sim/2013/LHC13c7c/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC12/Gam* > filesToCopy.txt
         mkdir -p $OUTPUTDIR_MC12/$runNumber
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC12/$runNumber/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC12/$runNumber/$fileToCopy ";
            else 
               alien_cp alien:/alice/sim/2013/LHC13c7c/$runNumber/PWGGA/GA_PbPb_MC/$LHC11hMC12/$fileToCopy file:$OUTPUTDIR_MC12/$runNumber/ 
            fi
         done
      done
    
   else

      if [ -z $LHC11hMC1 ]; then
         echo "no output defined for LHC12a17a_fix";
      else
         mkdir -p $OUTPUTDIR_MC1
         
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC1/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC1/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC1/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC1/merge/$fileToCopy file:$OUTPUTDIR_MC1/
            fi
         done
            
         ls $OUTPUTDIR_MC1/GammaConvV1_*.root > fileMC_$LHC11hMC1.txt
            fileNumbers=`cat fileMC_$LHC11hMC1.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC1/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17a_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17a_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17a_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17a_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi 
      
      if [ -z $LHC11hMC2 ]; then
         echo "no output defined for LHC12a17b_fix";
      else
         mkdir -p $OUTPUTDIR_MC2
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC2/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC2/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC2/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC2/merge/$fileToCopy file:$OUTPUTDIR_MC2/
            fi
         done
         
         ls $OUTPUTDIR_MC2/GammaConvV1_*.root > fileMC_$LHC11hMC2.txt
            fileNumbers=`cat fileMC_$LHC11hMC2.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC2/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17b_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17b_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17b_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17b_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi   

      if [ -z $LHC11hMC3 ]; then
         echo "no output defined for LHC12a17c_fix";
      else
         mkdir -p $OUTPUTDIR_MC3
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC3/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC3/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC3/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC3/merge/$fileToCopy file:$OUTPUTDIR_MC3/
            fi
         done
         
         
         ls $OUTPUTDIR_MC3/GammaConvV1_*.root > fileMC_$LHC11hMC3.txt
            fileNumbers=`cat fileMC_$LHC11hMC3.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC3/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17c_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17c_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17c_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17c_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi 
      
      if [ -z $LHC11hMC4 ]; then
         echo "no output defined for LHC12a17d_fix";
      else
         mkdir -p $OUTPUTDIR_MC4
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC4/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC4/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC4/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC4/merge/$fileToCopy file:$OUTPUTDIR_MC4/
            fi
         done
         
         
         ls $OUTPUTDIR_MC4/GammaConvV1_*.root > fileMC_$LHC11hMC4.txt
            fileNumbers=`cat fileMC_$LHC11hMC4.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC4/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17d_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17d_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17d_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17d_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
       fi

       if [ -z $LHC11hMC5 ]; then
         echo "no output defined for LHC12a17e_fix";
      else
         mkdir -p $OUTPUTDIR_MC5
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC5/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC5/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC5/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC5/merge/$fileToCopy file:$OUTPUTDIR_MC5/
            fi
         done
         
         
         ls $OUTPUTDIR_MC5/GammaConvV1_*.root > fileMC_$LHC11hMC5.txt
            fileNumbers=`cat fileMC_$LHC11hMC5.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC5/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17e_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17e_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17e_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17e_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
       fi  
      
       if [ -z $LHC11hMC6 ]; then
         echo "no output defined for LHC12a17f_fix";
      else
         mkdir -p $OUTPUTDIR_MC6
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC6/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC6/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC6/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC6/merge/$fileToCopy file:$OUTPUTDIR_MC6/
            fi
         done

         
         ls $OUTPUTDIR_MC6/GammaConvV1_*.root > fileMC_$LHC11hMC6.txt
            fileNumbers=`cat fileMC_$LHC11hMC6.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC6/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17f_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17f_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17f_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17f_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi
      
      if [ -z $LHC11hMC7 ]; then
         echo "no output defined for LHC12a17g_fix";
      else
         mkdir -p $OUTPUTDIR_MC7
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC7/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC7/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC7/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC7/merge/$fileToCopy file:$OUTPUTDIR_MC7/
            fi
         done

         
         ls $OUTPUTDIR_MC7/GammaConvV1_*.root > fileMC_$LHC11hMC7.txt
            fileNumbers=`cat fileMC_$LHC11hMC7.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC7/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17g_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17g_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17g_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17g_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
       fi
       
       if [ -z $LHC11hMC8 ]; then
         echo "no output defined for LHC12a17h_fix";
      else
         mkdir -p $OUTPUTDIR_MC8
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC8/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC8/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC8/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC8/merge/$fileToCopy file:$OUTPUTDIR_MC8/
            fi
         done

         
         ls $OUTPUTDIR_MC8/GammaConvV1_*.root > fileMC_$LHC11hMC8.txt
            fileNumbers=`cat fileMC_$LHC11hMC8.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC8/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17h_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17h_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17h_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17h_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi
      
      if [ -z $LHC11hMC9 ]; then
         echo "no output defined for LHC12a17i_fix";
      else
         mkdir -p $OUTPUTDIR_MC9
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC9/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC9/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC9/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC9/merge/$fileToCopy file:$OUTPUTDIR_MC9/
            fi
         done

         
         
         ls $OUTPUTDIR_MC9/GammaConvV1_*.root > fileMC_$LHC11hMC9.txt
            fileNumbers=`cat fileMC_$LHC11hMC9.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC9/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17i_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC12a17i_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12a17i_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC12a17i_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi
       
      if [ -z $LHC11hMC10  ]; then
         echo "no output defined for LHC13c7a_fix";
      else
         mkdir -p $OUTPUTDIR_MC10
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC10/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC10/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC10/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC10/merge/$fileToCopy file:$OUTPUTDIR_MC10/
            fi
         done

         
         ls $OUTPUTDIR_MC10/GammaConvV1_*.root > fileMC_$LHC11hMC10.txt
            fileNumbers=`cat fileMC_$LHC11hMC10.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC10/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13c7a_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13c7a_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13c7a_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC13c7a_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi
      
      if [ -z $LHC11hMC11 ]; then
          echo "no output defined for LHC13c7b_fix";
      else
     
         mkdir -p $OUTPUTDIR_MC11
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC11/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC11/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC11/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC11/merge/$fileToCopy file:$OUTPUTDIR_MC11/
            fi
         done

         
         ls $OUTPUTDIR_MC11/GammaConvV1_*.root > fileMC_$LHC11hMC11.txt
            fileNumbers=`cat fileMC_$LHC11hMC11.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC11/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13c7b_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13c7b_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13c7b_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC13c7b_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done
      fi

      if [ -z $LHC11hMC12 ]; then
         echo "no output defined for LHC13c7c_fix";
      else
         mkdir -p $OUTPUTDIR_MC12
         alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC12/merge/Gam* > filesToCopy.txt
         files=`cat filesToCopy.txt`
         for fileToCopy in $files; do
            if [ -f $OUTPUTDIR_MC12/$fileToCopy ]; then
               echo "already copied $OUTPUTDIR_MC12/$fileToCopy ";
            else 
               alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC11hMC12/merge/$fileToCopy file:$OUTPUTDIR_MC12/
            fi
         done

         
         ls $OUTPUTDIR_MC12/GammaConvV1_*.root > fileMC_$LHC11hMC12.txt
            fileNumbers=`cat fileMC_$LHC11hMC12.txt`
            for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC12/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13c7c_fix_$number.root\"\,\"GammaConvV1_$number\"\)
               root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13c7c_fix_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13c7c_fix_$number.log\"\)
         #cp $OUTPUTDIR/CutSelection_MC_LHC13c7c_fix_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
         done	
      fi
	fi

# elif [ $1 = "MC" ]; then
# 	mkdir -p $OUTPUTDIR_MC
# 	alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/$LHC10hMC/merge/Gam* file:$OUTPUTDIR_MC/
# 	
# 	ls $OUTPUTDIR_MC/GammaConvV1_*.root > fileMC.txt
#        fileNumbers=`cat fileMC.txt`
#        for fileName in $fileNumbers; do
#        echo $fileName
#        number=`echo $fileName | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
#        echo $number
#        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC10h_$number.root\"\,\"GammaConvV1_$number\"\)
#      	root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC10h_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_AOD_LHC10h_$number.log\"\)
# 	cp $OUTPUTDIR/CutSelection_MC_AOD_LHC10h_$number.log  /alidata50/alice_u/leardini/photonconv/AnalysisSoftware/.
# 	done

fi