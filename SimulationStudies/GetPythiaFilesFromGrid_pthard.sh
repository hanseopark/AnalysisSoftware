#! /bin/bash

# copies files from grid
# creates directory

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
elif [ $1 = "bsahlmul" ]; then
  BASEDIR=/Users/sahlmul/Documents/ALICE/mycode/AnalysisSoftware/SimulationStudies/grid/
elif [ $1 = "bsahlmul2" ]; then
    BASEDIR=/Users/sahlmul/alice/photonconv/AnalysisSoftware/SimulationStudies/grid/
fi

#TRAINDIR=Legotrain_vAN-20150823-Pythia8
#Bin5Monash="237_20150824-1144"
#Bin4Monash="236_20150824-1143"
#Bin3Monash="235_20150824-1142"
#Bin2Monash="234_20150824-1142"
#Bin1Monash="233_20150824-1142"
#MBMonash="232_20150824-1132"
#Bin5Tune4C="229_20150824-1121"
#Bin4Tune4C="228_20150824-1120"
#Bin3Tune4C="227_20150824-1120"
#Bin2Tune4C="226_20150824-1122"
#Bin1Tune4C="225_20150824-1114"
#MBTune4C="230_20150824-1125"

#TRAINDIR=Legotrain_vAN-20150825-Pythia6
#Bin5Perugia2011="245_20150825-2114"
#Bin4Perugia2011="244_20150825-2113"
#Bin3Perugia2011="243_20150825-2112"
#Bin2Perugia2011="242_20150825-2112"
#Bin1Perugia2011="241_20150825-2112"

TRAINDIR=Legotrain_vAN-20150825-Pythia8
Bin1Monash="247_20150826-1233"
Bin2Monash="248_20150826-1234"
Bin3Monash="249_20150826-1234"
Bin4Monash="250_20150826-1235"
Bin5Monash="246_20150826-1220"
Bin6Monash="251_20150826-1512"
Bin7Monash="252_20150826-1513"
Bin8Monash="253_20150826-1514"
Bin9Monash="254_20150826-1513"
MBMonash="264_20150828-0950"

#TRAINDIR=Legotrain_vAN-20150825-Pythia8-Wide
#Bin1Monash="255_20150826-1531"
#Bin2Monash="256_20150826-1531"
#Bin3Monash="257_20150826-1531"
#Bin4Monash="258_20150826-1531"
#MBMonash="264_20150828-0950"

#TRAINDIR=Legotrain_vAN-20150825-Pythia8-WideFull
#Bin1Monash="255_20150826-1531"
#Bin2Monash="256_20150826-1531"
#Bin3Monash="257_20150826-1531"
#Bin4Monash="258_20150826-1531"
#Bin1MoreMonash="265_20150828-1012"
#Bin2MoreMonash="266_20150828-1012"
#Bin3MoreMonash="267_20150828-1013"
#Bin4MoreMonash="268_20150828-1013"
#MBMonash="264_20150828-0950"

OUTPUTDIR=$BASEDIR/$TRAINDIR

OUTPUTDIRBin9Monash=$BASEDIR/$TRAINDIR/$Bin9Monash
OUTPUTDIRBin8Monash=$BASEDIR/$TRAINDIR/$Bin8Monash
OUTPUTDIRBin7Monash=$BASEDIR/$TRAINDIR/$Bin7Monash
OUTPUTDIRBin6Monash=$BASEDIR/$TRAINDIR/$Bin6Monash
OUTPUTDIRBin5Monash=$BASEDIR/$TRAINDIR/$Bin5Monash
OUTPUTDIRBin4Monash=$BASEDIR/$TRAINDIR/$Bin4Monash
OUTPUTDIRBin3Monash=$BASEDIR/$TRAINDIR/$Bin3Monash
OUTPUTDIRBin2Monash=$BASEDIR/$TRAINDIR/$Bin2Monash
OUTPUTDIRBin1Monash=$BASEDIR/$TRAINDIR/$Bin1Monash
OUTPUTDIRMBMonash=$BASEDIR/$TRAINDIR/$MBMonash
#OUTPUTDIRBin4MoreMonash=$BASEDIR/$TRAINDIR/$Bin4MoreMonash
#OUTPUTDIRBin3MoreMonash=$BASEDIR/$TRAINDIR/$Bin3MoreMonash
#OUTPUTDIRBin2MoreMonash=$BASEDIR/$TRAINDIR/$Bin2MoreMonash
#OUTPUTDIRBin1MoreMonash=$BASEDIR/$TRAINDIR/$Bin1MoreMonash
#OUTPUTDIRBin5Tune4C=$BASEDIR/$TRAINDIR/$Bin5Tune4C
#OUTPUTDIRBin4Tune4C=$BASEDIR/$TRAINDIR/$Bin4Tune4C
#OUTPUTDIRBin3Tune4C=$BASEDIR/$TRAINDIR/$Bin3Tune4C
#OUTPUTDIRBin2Tune4C=$BASEDIR/$TRAINDIR/$Bin2Tune4C
#OUTPUTDIRBin1Tune4C=$BASEDIR/$TRAINDIR/$Bin1Tune4C
#OUTPUTDIRMBTune4C=$BASEDIR/$TRAINDIR/$MBTune4C
#OUTPUTDIRBin5Perugia2011=$BASEDIR/$TRAINDIR/$Bin5Perugia2011
#OUTPUTDIRBin4Perugia2011=$BASEDIR/$TRAINDIR/$Bin4Perugia2011
#OUTPUTDIRBin3Perugia2011=$BASEDIR/$TRAINDIR/$Bin3Perugia2011
#OUTPUTDIRBin2Perugia2011=$BASEDIR/$TRAINDIR/$Bin2Perugia2011
#OUTPUTDIRBin1Perugia2011=$BASEDIR/$TRAINDIR/$Bin1Perugia2011

mkdir -p $OUTPUTDIRBin9Monash
mkdir -p $OUTPUTDIRBin8Monash
mkdir -p $OUTPUTDIRBin7Monash
mkdir -p $OUTPUTDIRBin6Monash
mkdir -p $OUTPUTDIRBin5Monash
mkdir -p $OUTPUTDIRBin4Monash
mkdir -p $OUTPUTDIRBin3Monash
mkdir -p $OUTPUTDIRBin2Monash
mkdir -p $OUTPUTDIRBin1Monash
mkdir -p $OUTPUTDIRMBMonash
#mkdir -p $OUTPUTDIRBin4MoreMonash
#mkdir -p $OUTPUTDIRBin3MoreMonash
#mkdir -p $OUTPUTDIRBin2MoreMonash
#mkdir -p $OUTPUTDIRBin1MoreMonash
#mkdir -p $OUTPUTDIRBin5Tune4C
#mkdir -p $OUTPUTDIRBin4Tune4C
#mkdir -p $OUTPUTDIRBin3Tune4C
#mkdir -p $OUTPUTDIRBin2Tune4C
#mkdir -p $OUTPUTDIRBin1Tune4C
#mkdir -p $OUTPUTDIRMBTune4C
#mkdir -p $OUTPUTDIRBin5Perugia2011
#mkdir -p $OUTPUTDIRBin4Perugia2011
#mkdir -p $OUTPUTDIRBin3Perugia2011
#mkdir -p $OUTPUTDIRBin2Perugia2011
#mkdir -p $OUTPUTDIRBin1Perugia2011

alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin9Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin9Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin8Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin8Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin7Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin7Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin6Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin6Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin5Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin5Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin4Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin4Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin3Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin3Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin2Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin2Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin1Monash/merge/AnalysisResults.root file:$OUTPUTDIRBin1Monash
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin4MoreMonash/merge/AnalysisResults.root file:$OUTPUTDIRBin4MoreMonash
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin3MoreMonash/merge/AnalysisResults.root file:$OUTPUTDIRBin3MoreMonash
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin2MoreMonash/merge/AnalysisResults.root file:$OUTPUTDIRBin2MoreMonash
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin1MoreMonash/merge/AnalysisResults.root file:$OUTPUTDIRBin1MoreMonash
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin5Tune4C/merge/AnalysisResults.root file:$OUTPUTDIRBin5Tune4C
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin4Tune4C/merge/AnalysisResults.root file:$OUTPUTDIRBin4Tune4C
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin3Tune4C/merge/AnalysisResults.root file:$OUTPUTDIRBin3Tune4C
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin2Tune4C/merge/AnalysisResults.root file:$OUTPUTDIRBin2Tune4C
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin1Tune4C/merge/AnalysisResults.root file:$OUTPUTDIRBin1Tune4C
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBTune4C/merge/AnalysisResults.root file:$OUTPUTDIRMBTune4C
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin5Perugia2011/merge/AnalysisResults.root file:$OUTPUTDIRBin5Perugia2011
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin4Perugia2011/merge/AnalysisResults.root file:$OUTPUTDIRBin4Perugia2011
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin3Perugia2011/merge/AnalysisResults.root file:$OUTPUTDIRBin3Perugia2011
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin2Perugia2011/merge/AnalysisResults.root file:$OUTPUTDIRBin2Perugia2011
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$Bin1Perugia2011/merge/AnalysisResults.root file:$OUTPUTDIRBin1Perugia2011

### let's move these files to some better naming scheme
echo "copying files to common directory"

cp $OUTPUTDIRBin9Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin9Monash.root
cp $OUTPUTDIRBin8Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin8Monash.root
cp $OUTPUTDIRBin7Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin7Monash.root
cp $OUTPUTDIRBin6Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin6Monash.root
cp $OUTPUTDIRBin5Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin5Monash.root
cp $OUTPUTDIRBin4Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin4Monash.root
cp $OUTPUTDIRBin3Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin3Monash.root
cp $OUTPUTDIRBin2Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin2Monash.root
cp $OUTPUTDIRBin1Monash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin1Monash.root
#cp $OUTPUTDIRBin4MoreMonash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin4MoreMonash.root
#cp $OUTPUTDIRBin3MoreMonash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin3MoreMonash.root
#cp $OUTPUTDIRBin2MoreMonash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin2MoreMonash.root
#cp $OUTPUTDIRBin1MoreMonash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin1MoreMonash.root
cp $OUTPUTDIRMBMonash/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash.root
#cp $OUTPUTDIRBin5Tune4C/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin5Tune4C.root
#cp $OUTPUTDIRBin4Tune4C/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin4Tune4C.root
#cp $OUTPUTDIRBin3Tune4C/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin3Tune4C.root
#cp $OUTPUTDIRBin2Tune4C/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin2Tune4C.root
#cp $OUTPUTDIRBin1Tune4C/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin1Tune4C.root
#cp $OUTPUTDIRMBTune4C/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBTune4C.root
#cp $OUTPUTDIRBin5Perugia2011/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin5Perugia2011.root
#cp $OUTPUTDIRBin4Perugia2011/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin4Perugia2011.root
#cp $OUTPUTDIRBin3Perugia2011/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin3Perugia2011.root
#cp $OUTPUTDIRBin2Perugia2011/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin2Perugia2011.root
#cp $OUTPUTDIRBin1Perugia2011/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin1Perugia2011.root

# dirty hack
#rm $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin4Monash.root
#hadd $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin4Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin4Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin4MoreMonash.root
#rm $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin3Monash.root
#hadd $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin3Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin3Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin3MoreMonash.root
#rm $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin2Monash.root
#hadd $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin2Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin2Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin2MoreMonash.root
#rm $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin1Monash.root
#hadd $BASEDIR/$TRAINDIR/PythiaAnalysisResultsFullBin1Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin1Monash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsBin1MoreMonash.root

