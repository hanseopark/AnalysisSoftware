#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

if [ $1 = "fbock" ]; then 
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/PbPb
    NSlashes=8
    NSlashes2=7
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/PbPb
elif [ $1 = "leardini" ]; then 
   BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then 
   BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
elif [ $1 = "leardiniGSI" ]; then 
   BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/PbPb
elif [ $1 = "passfeld" ]; then 
   BASEDIR=~/work/Gridoutput/PbPb
elif [ $1 = "passfeldMAF" ]; then 
   BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then  
   BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/PbPb
elif [ $1 = "amarin" ]; then     
   BASEDIR=/Users/marin/
elif [ $1 = "amarinGSI" ]; then     
   BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/PbPb 
elif [ $1 = "amarinALICESERV1" ]; then     
   BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/   
elif [ $1 = "mwilde" ]; then        
   BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
   BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/PbPb 
elif [ $1 = "pgonzales" ]; then     
   BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
   BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/PbPb
elif [ $1 = "dmuhlheim" ]; then 
   BASEDIR=/home/daniel/Desktop/Grid
fi
mkdir -p $BASEDIR
  
# 	TRAINDIR=Legotrain-vAN-20141110-Calo
# 	LHC11hData=144_20141111-1840; #ESD
# 	LHC14a1a=176_20141111-1928; #ESD
# 	LHC14a1b=177_20141111-2033; #ESD 
# 	LHC14a1c=; 

	TRAINDIR=Legotrain-vAN-20141116-Calo
	LHC11hData=147_20141117-1517; #ESD
	LHC14a1a=186_20141124-1555;
	LHC14a1b=187_20141124-1556; #ESD 

	
	
OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $2 = "AODdata" ]; then
   TRAINPATHData=GA_PbPb_AOD
else
   TRAINPATHData=GA_PbPb
fi   
OUTPUTDIR_LHC11h_A=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData-allRuns
OUTPUTDIR_LHC11h_B=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData-PCMgood
OUTPUTDIR_LHC11h_C=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData-TPCOROCC8Problems
OUTPUTDIR_LHC11h_D=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData-TPCIROCC13Problems

if [ $3 = "AODmc" ]; then
   TRAINPATHMC=GA_PbPb_MC_AOD
else
   TRAINPATHMC=GA_PbPb_MC
fi   
OUTPUTDIR_LHC14a1a_A=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a-allRuns
OUTPUTDIR_LHC14a1b_A=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b-allRuns
# OUTPUTDIR_LHC14a1c_A=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1c-allRuns
OUTPUTDIR_LHC14a1a_B=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a-PCMgood
OUTPUTDIR_LHC14a1b_B=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b-PCMgood
# OUTPUTDIR_LHC14a1c_B=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1c-PCMgood
OUTPUTDIR_LHC14a1a_C=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a-TPCOROCC8Problems
OUTPUTDIR_LHC14a1b_C=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b-TPCOROCC8Problems
# OUTPUTDIR_LHC14a1c_C=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1c-TPCOROCC8Problems
OUTPUTDIR_LHC14a1a_D=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a-TPCIROCC13Problems
OUTPUTDIR_LHC14a1b_D=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b-TPCIROCC13Problems
# OUTPUTDIR_LHC14a1c_D=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1c-TPCIROCC13Problems


mkdir -p $OUTPUTDIR_LHC11h_A
mkdir -p $OUTPUTDIR_LHC14a1a_A
mkdir -p $OUTPUTDIR_LHC14a1b_A
# mkdir -p $OUTPUTDIR_LHC14a1c_A
mkdir -p $OUTPUTDIR_LHC11h_B
mkdir -p $OUTPUTDIR_LHC14a1a_B
mkdir -p $OUTPUTDIR_LHC14a1b_B
# mkdir -p $OUTPUTDIR_LHC14a1c_B
mkdir -p $OUTPUTDIR_LHC11h_C
mkdir -p $OUTPUTDIR_LHC14a1a_C
mkdir -p $OUTPUTDIR_LHC14a1b_C
# mkdir -p $OUTPUTDIR_LHC14a1c_C
mkdir -p $OUTPUTDIR_LHC11h_D
mkdir -p $OUTPUTDIR_LHC14a1a_D
mkdir -p $OUTPUTDIR_LHC14a1b_D
# mkdir -p $OUTPUTDIR_LHC14a1c_D

# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge/GammaConvCalo* file:$OUTPUTDIR_LHC14a1a_A/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge/GammaConvCalo* file:$OUTPUTDIR_LHC14a1b_A/

alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_1/GammaConvCalo* file:$OUTPUTDIR_LHC14a1a_A/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/GammaConvCalo* file:$OUTPUTDIR_LHC14a1b_A/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1c/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC14a1c_A/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_1/GammaConvCalo* file:$OUTPUTDIR_LHC11h_A/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_2/GammaConvCalo* file:$OUTPUTDIR_LHC14a1a_B/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_2/GammaConvCalo* file:$OUTPUTDIR_LHC14a1b_B/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1c/merge_runlist_2/GammaConvCalo* file:$OUTPUTDIR_LHC14a1c_B/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_2/GammaConvCalo* file:$OUTPUTDIR_LHC11h_B/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_3/GammaConvCalo* file:$OUTPUTDIR_LHC14a1a_C/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_3/GammaConvCalo* file:$OUTPUTDIR_LHC14a1b_C/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1c/merge_runlist_3/GammaConvCalo* file:$OUTPUTDIR_LHC14a1c_C/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_3/GammaConvCalo* file:$OUTPUTDIR_LHC11h_C/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_4/GammaConvCalo* file:$OUTPUTDIR_LHC14a1a_D/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_4/GammaConvCalo* file:$OUTPUTDIR_LHC14a1b_D/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1c/merge_runlist_4/GammaConvCalo* file:$OUTPUTDIR_LHC14a1c_D/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_4/GammaConvCalo* file:$OUTPUTDIR_LHC11h_D/



if [ $3 = "AODmc" ]; then
   ls $OUTPUTDIR_LHC14a1a_A/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC14a1a_AOD_$number.log\"\)

   done;
   
   ls $OUTPUTDIR_LHC14a1b_A/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC14a1b_AOD_$number.log\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_A/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC14a1c_AOD_$number.log\"\)
   done;

   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_AllRCTGoodRuns_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_AllRCTGoodRuns_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 8 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_AllRCTGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_AllRCTGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_AllRCTGoodRuns_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_AllRCTGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_AllRCTGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_AllRCTGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_AllRCTGoodRuns_$number.root
      fi
   done

   # runlist PCM good runs
   ls $OUTPUTDIR_LHC14a1a_B/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1b_B/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_B/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_PCMGoodRuns_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_PCMGoodRuns_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 8 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_PCMGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_PCMGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_PCMGoodRuns_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_PCMGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_PCMGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_PCMGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_PCMGoodRuns_$number.root
      fi
   done

   # runlist TPC OROC C8 problems
   ls $OUTPUTDIR_LHC14a1a_C/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1b_C/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_C/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_TPCOROCC8Problems_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCOROCC8Problems_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 8 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCOROCC8Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCOROCC8Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_TPCOROCC8Problems_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_TPCOROCC8Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCOROCC8Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCOROCC8Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCOROCC8Problems_$number.root
      fi
   done

    # runlist TPC IROC C13 problems
   ls $OUTPUTDIR_LHC14a1a_D/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1b_D/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_D/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_TPCIROCC13Problems_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCIROCC13Problems_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 8 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCIROCC13Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCIROCC13Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1c_TPCIROCC13Problems_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_LHC14a1c_TPCIROCC13Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1a_TPCIROCC13Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCIROCC13Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_AOD_LHC14a1b_TPCIROCC13Problems_$number.root
      fi
   done
 
   
else 
   # all good runs according to RCT
   ls $OUTPUTDIR_LHC14a1a_A/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC14a1a_$number.log\"\)

   done;
   
   ls $OUTPUTDIR_LHC14a1b_A/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC14a1b_$number.log\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_A/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC14a1c_$number.log\"\)
   done;
   
   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_AllRCTGoodRuns_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_AllRCTGoodRuns_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 7 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_AllRCTGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_AllRCTGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_AllRCTGoodRuns_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_AllRCTGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_AllRCTGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_AllRCTGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_AllRCTGoodRuns_$number.root
      fi
   done
   
   # all PCM good runs
   ls $OUTPUTDIR_LHC14a1a_B/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1b_B/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_B/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_PCMGoodRuns_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_PCMGoodRuns_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 7 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_PCMGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_PCMGoodRuns_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_PCMGoodRuns_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_PCMGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_PCMGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_PCMGoodRuns_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_PCMGoodRuns_$number.root
      fi
   done
   
   # all TPC OROC C8 problems
   ls $OUTPUTDIR_LHC14a1a_C/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1b_C/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_C/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_TPCOROCC8Problems_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCOROCC8Problems_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 7 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCOROCC8Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_TPCOROCC8Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_TPCOROCC8Problems_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_TPCOROCC8Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCOROCC8Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_TPCOROCC8Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_TPCOROCC8Problems_$number.root
      fi
   done
   
   # all TPC IROC C13 problems
   ls $OUTPUTDIR_LHC14a1a_D/GammaConvCalo_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1a_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1b_D/GammaConvCalo_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1b_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   ls $OUTPUTDIR_LHC14a1c_D/GammaConvCalo_*.root > fileLHC14a1c.txt
   fileNumbersc=`cat fileLHC14a1c.txt`
   for fileName in $fileNumbersc; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14a1c_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
   rm $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_TPCIROCC13Problems_*.root
   ls $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCIROCC13Problems_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 7 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCIROCC13Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_TPCIROCC13Problems_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_TPCIROCC13Problems_$number.root ]; then
         hadd -f $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_LHC14a1b_LHC14a1c_TPCIROCC13Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1a_TPCIROCC13Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1b_TPCIROCC13Problems_$number.root $OUTPUTDIR/GammaConvCalo_GA_PbPb_MC_LHC14a1c_TPCIROCC13Problems_$number.root
      fi
   done

fi


if [ $2 = "AODdata" ]; then
   ls $OUTPUTDIR_LHC11h_A/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC11h_$number.log\"\)
   done;
   
   ls $OUTPUTDIR_LHC11h_B/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

   ls $OUTPUTDIR_LHC11h_C/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

   ls $OUTPUTDIR_LHC11h_D/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;
   
else 
   ls $OUTPUTDIR_LHC11h_A/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_A/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_AllRCTGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_AllRCTGoodRuns_$number.root\"\,\"$OUTPUTDIR/CutSelectionConvCalo_LHC11h_$number.log\"\)
   done;

   ls $OUTPUTDIR_LHC11h_B/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_B/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_PCMGoodRuns_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

   ls $OUTPUTDIR_LHC11h_C/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_C/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_TPCOROCC8Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

   ls $OUTPUTDIR_LHC11h_D/GammaConvCalo_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC11h_D/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_GA_PbPb_LHC11h-pass2_TPCIROCC13Problems_$number.root\"\,\"GammaConvCalo_$number\"\)
   done;

fi   
