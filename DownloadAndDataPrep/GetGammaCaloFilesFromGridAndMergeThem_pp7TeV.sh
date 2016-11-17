#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

DEPTH=12;
NSlashes=9
NSlashes2=8
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
elif [ $1 = "dmuhlhei" ]; then
   BASEDIR=~/data/work/Grid
elif [ $1 = "bsahlmul" ]; then
   BASEDIR=/Users/sahlmul/Documents/ALICE/mycode/photonconv/AnalysisSoftware/grid
  NSlashes=12
  NSlashes2=11
fi

#TRAINDIR=Legotrain-vAN-20150922-LHC10_p4-QA
#LHC10bData=997_20150922-0953;
#LHC10cData=998_20150922-0953;
#LHC10dData=999_20150922-0954;
#LHC10eData=1000_20150922-0954;
#LHC10fData=1001_20150922-0954;
#LHC14j4bMC=1080_20150922-0945;
#LHC14j4cMC=1081_20150922-0946;
#LHC14j4dMC=1082_20150922-0946;
#LHC14j4eMC=1093_20150923-1117;
#LHC14j4fMC=1084_20150922-0947;

TRAINDIR=Legotrain-vAN-20151123-LHC10_NL_iter_1
LHC10bData=1097_20151123-1136;
LHC10cData=1098_20151123-1137;
LHC10dData=1099_20151123-1137;
LHC10eData=1100_20151123-1137;
LHC10fData=1101_20151123-1138;
LHC14j4bMC=1321_20151123-1140;
LHC14j4cMC=1322_20151123-1140;
LHC14j4dMC=1323_20151123-1141;
LHC14j4eMC=1324_20151123-1141;
LHC14j4fMC=1325_20151123-1141;

OUTPUTDIR=$BASEDIR/$TRAINDIR
mkdir -p $OUTPUTDIR
OUTPUTDIR_LHC10b=$BASEDIR/$TRAINDIR/GA_pp-$LHC10bData
OUTPUTDIR_LHC10c=$BASEDIR/$TRAINDIR/GA_pp-$LHC10cData
OUTPUTDIR_LHC10d=$BASEDIR/$TRAINDIR/GA_pp-$LHC10dData
OUTPUTDIR_LHC10e=$BASEDIR/$TRAINDIR/GA_pp-$LHC10eData
OUTPUTDIR_LHC10f=$BASEDIR/$TRAINDIR/GA_pp-$LHC10fData
OUTPUTDIR_LHC14j4b=$BASEDIR/$TRAINDIR/GA_pp-$LHC14j4bMC
OUTPUTDIR_LHC14j4c=$BASEDIR/$TRAINDIR/GA_pp-$LHC14j4cMC
OUTPUTDIR_LHC14j4d=$BASEDIR/$TRAINDIR/GA_pp-$LHC14j4dMC
OUTPUTDIR_LHC14j4e=$BASEDIR/$TRAINDIR/GA_pp-$LHC14j4eMC
OUTPUTDIR_LHC14j4f=$BASEDIR/$TRAINDIR/GA_pp-$LHC14j4fMC
mkdir -p $OUTPUTDIR_LHC10b
mkdir -p $OUTPUTDIR_LHC10c
mkdir -p $OUTPUTDIR_LHC10d
mkdir -p $OUTPUTDIR_LHC10e
mkdir -p $OUTPUTDIR_LHC10f
mkdir -p $OUTPUTDIR_LHC14j4b
mkdir -p $OUTPUTDIR_LHC14j4c
mkdir -p $OUTPUTDIR_LHC14j4d
mkdir -p $OUTPUTDIR_LHC14j4e
mkdir -p $OUTPUTDIR_LHC14j4f

# get files from grid
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10bData/merge_runlist_3/GammaCalo* file:$OUTPUTDIR_LHC10b/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10cData/merge_runlist_3/GammaCalo* file:$OUTPUTDIR_LHC10c/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10dData/merge_runlist_3/GammaCalo* file:$OUTPUTDIR_LHC10d/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10eData/merge_runlist_3/GammaCalo* file:$OUTPUTDIR_LHC10e/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10fData/merge_runlist_3/GammaCalo* file:$OUTPUTDIR_LHC10f/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4bMC/merge/GammaCalo* file:$OUTPUTDIR_LHC14j4b/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4cMC/merge/GammaCalo* file:$OUTPUTDIR_LHC14j4c/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4dMC/merge/GammaCalo* file:$OUTPUTDIR_LHC14j4d/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4eMC/merge/GammaCalo* file:$OUTPUTDIR_LHC14j4e/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4fMC/merge/GammaCalo* file:$OUTPUTDIR_LHC14j4f/

# Change structure of data files
ls $OUTPUTDIR_LHC10b/GammaCalo_*.root > fileLHC10b.txt
fileNumbers=`cat fileLHC10b.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC10b/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC10b-pass4_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_LHC10b-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC10b_$number.log\"\,4\)
done;

ls $OUTPUTDIR_LHC10c/GammaCalo_*.root > fileLHC10c.txt
fileNumbers=`cat fileLHC10c.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC10c/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC10c-pass4_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_LHC10c-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC10c_$number.log\"\,4\)
done;

ls $OUTPUTDIR_LHC10d/GammaCalo_*.root > fileLHC10d.txt
fileNumbers=`cat fileLHC10d.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC10d/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC10d-pass4_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_LHC10d-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC10d_$number.log\"\,4\)
done;

ls $OUTPUTDIR_LHC10e/GammaCalo_*.root > fileLHC10e.txt
fileNumbers=`cat fileLHC10e.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC10e/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC10e-pass4_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_LHC10e-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC10e_$number.log\"\,4\)
done;

ls $OUTPUTDIR_LHC10f/GammaCalo_*.root > fileLHC10f.txt
fileNumbers=`cat fileLHC10f.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC10f/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC10f-pass4_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_LHC10f-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC10f_$number.log\"\,4\)
done;

# Change structure of MC files
ls $OUTPUTDIR_LHC14j4b/GammaCalo_*.root > fileLHC14j4b.txt
fileNumbers=`cat fileLHC14j4b.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC14j4b/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC14j4b_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14j4b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14j4b_$number.log\"\,4\)
done;

ls $OUTPUTDIR_LHC14j4c/GammaCalo_*.root > fileLHC14j4c.txt
fileNumbers=`cat fileLHC14j4c.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC14j4c/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC14j4c_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14j4c_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14j4c_$number.log\"\,4\)
done;

ls $OUTPUTDIR_LHC14j4d/GammaCalo_*.root > fileLHC14j4d.txt
fileNumbers=`cat fileLHC14j4d.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC14j4d/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC14j4d_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14j4d_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14j4d_$number.log\"\,4\)
done;

ls $OUTPUTDIR_LHC14j4e/GammaCalo_*.root > fileLHC14j4e.txt
fileNumbers=`cat fileLHC14j4e.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC14j4e/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC14j4e_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14j4e_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14j4e_$number.log\"\,4\)
done;


ls $OUTPUTDIR_LHC14j4f/GammaCalo_*.root > fileLHC14j4f.txt
fileNumbers=`cat fileLHC14j4f.txt`
for fileName in $fileNumbers; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
echo $number
root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC14j4f/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC14j4f_$number.root\"\,\"GammaCalo_$number\"\)
root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14j4f_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14j4f_$number.log\"\,4\)
done;

# Merge data files
rm $OUTPUTDIR/GammaCalo_LHC10bc*-pass4_*.root
ls $OUTPUTDIR/GammaCalo_LHC10b-pass4_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_LHC10b-pass4_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC10c-pass4_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_LHC10bc-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10b-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10c-pass4_$number.root
fi
done

ls $OUTPUTDIR/GammaCalo_LHC10bc-pass4_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_LHC10bc-pass4_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC10d-pass4_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_LHC10bcd-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10bc-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10d-pass4_$number.root
fi
done

ls $OUTPUTDIR/GammaCalo_LHC10bcd-pass4_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_LHC10bcd-pass4_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC10e-pass4_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_LHC10bcde-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10bcd-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10e-pass4_$number.root
fi
done

ls $OUTPUTDIR/GammaCalo_LHC10bcde-pass4_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_LHC10bcde-pass4_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC10f-pass4_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_LHC10bcdef-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10bcde-pass4_$number.root $OUTPUTDIR/GammaCalo_LHC10f-pass4_$number.root
fi
done

# merge MC files
rm $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_*.root
ls $OUTPUTDIR/GammaCalo_MC_LHC14j4b_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4c_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4b_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4c_$number.root
fi
done

ls $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4d_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4d_$number.root
fi
done

ls $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 6 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4e_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_LHC14j4e_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4e_$number.root
fi
done

ls $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_LHC14j4e_*.root > filesForMerging.txt
filesForMerging=`cat filesForMerging.txt`
for fileName in $filesForMerging; do
echo $fileName
number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 7 | cut -d "." -f1`
echo $number
if [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_LHC14j4e_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC14j4f_$number.root ] ; then
hadd -f $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_LHC14j4e_LHC14j4f_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4b_LHC14j4c_LHC14j4d_LHC14j4e_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14j4f_$number.root
fi
done


