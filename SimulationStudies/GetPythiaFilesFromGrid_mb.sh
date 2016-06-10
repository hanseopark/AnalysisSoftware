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
  BASEDIR=/Users/sahlmul/Documents/ALICE/mycode/PCMEMC_Software/AnalysisSoftware/SimulationStudies
elif [ $1 = "bsahlmul2" ]; then
    BASEDIR=/Users/sahlmul/alice/photonconv/AnalysisSoftware/SimulationStudies/grid/
fi


TRAINDIR=Legotrain_vAN-20150825-Pythia8

MBMonash1="272_20150903-1229"
MBMonash2="273_20150903-2049"
MBMonash3="274_20150903-2142"
MBMonash4="275_20150904-2202"
MBMonash5="276_20150905-2115"
MBMonash6="277_20150906-1130"
MBMonash7="278_20150909-2123"
MBMonash8="279_20150910-0942"
MBMonash9="280_20150911-0922"
MBMonash10="281_20150912-1143"

OUTPUTDIR=$BASEDIR/$TRAINDIR

OUTPUTDIRMBMonash1=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash2=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash3=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash4=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash5=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash6=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash7=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash8=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash9=$BASEDIR/$TRAINDIR/$MBMonash
OUTPUTDIRMBMonash10=$BASEDIR/$TRAINDIR/$MBMonash

mkdir -p $OUTPUTDIRMBMonash1
mkdir -p $OUTPUTDIRMBMonash2
mkdir -p $OUTPUTDIRMBMonash3
mkdir -p $OUTPUTDIRMBMonash4
mkdir -p $OUTPUTDIRMBMonash5
mkdir -p $OUTPUTDIRMBMonash6
mkdir -p $OUTPUTDIRMBMonash7
mkdir -p $OUTPUTDIRMBMonash8
mkdir -p $OUTPUTDIRMBMonash9
mkdir -p $OUTPUTDIRMBMonash10

alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash1/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash1
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash2/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash2
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash3/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash3
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash4/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash4
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash5/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash5
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash6/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash6
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash7/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash7
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash8/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash8
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash9/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash9
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGZZ/MCGen_pp/$MBMonash10/merge/AnalysisResults.root file:$OUTPUTDIRMBMonash10

### let's move these files to some better naming scheme
echo "copying files to common directory"

cp $OUTPUTDIRMBMonash1/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash1.root
cp $OUTPUTDIRMBMonash2/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash2.root
cp $OUTPUTDIRMBMonash3/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash3.root
cp $OUTPUTDIRMBMonash4/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash4.root
cp $OUTPUTDIRMBMonash5/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash5.root
cp $OUTPUTDIRMBMonash6/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash6.root
cp $OUTPUTDIRMBMonash7/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash7.root
cp $OUTPUTDIRMBMonash8/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash8.root
cp $OUTPUTDIRMBMonash9/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash9.root
cp $OUTPUTDIRMBMonash10/AnalysisResults.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash10.root

### let's hadd them
hadd $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash?.root $BASEDIR/$TRAINDIR/PythiaAnalysisResultsMBMonash??.root