#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

BASEDIR=~/Photon/Results/pPbAnalysisJune/pPbRootFiles

TRAINDIR=Legotrain-v5-04-68-AN
LHC13bData=49_20130618-1451; 
LHC13cData=46_20130618-1134; 
LHC13b2_efix_p1MC=27_20130618-1830; 

OUTPUTDIR=$BASEDIR/$TRAINDIR
OUTPUTDIR_LHC13b=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13bData
OUTPUTDIR_LHC13c=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13cData
OUTPUTDIR_LHC13b2_efix_p1=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p1MC
mkdir -p $OUTPUTDIR_LHC13b
mkdir -p $OUTPUTDIR_LHC13c
mkdir -p $OUTPUTDIR_LHC13b2_efix_p1

alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13bData/merge/Gam* file:$OUTPUTDIR_LHC13b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge/Gam* file:$OUTPUTDIR_LHC13c/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/merge/Gam* file:$OUTPUTDIR_LHC13b2_efix_p1/

root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b/GammaConvV1_1.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b-pass3_1.root\"\,\"GammaConvV1_1\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b/GammaConvV1_2.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b-pass3_2.root\"\,\"GammaConvV1_2\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b/GammaConvV1_3.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b-pass3_3.root\"\,\"GammaConvV1_3\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b/GammaConvV1_4.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b-pass3_4.root\"\,\"GammaConvV1_4\"\)

root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13c/GammaConvV1_1.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13c-pass2_1.root\"\,\"GammaConvV1_1\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13c/GammaConvV1_2.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13c-pass2_2.root\"\,\"GammaConvV1_2\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13c/GammaConvV1_3.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13c-pass2_3.root\"\,\"GammaConvV1_3\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13c/GammaConvV1_4.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13c-pass2_4.root\"\,\"GammaConvV1_4\"\)

root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_1.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b2_efix_part1_1.root\"\,\"GammaConvV1_1\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_2.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b2_efix_part1_2.root\"\,\"GammaConvV1_2\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_3.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b2_efix_part1_3.root\"\,\"GammaConvV1_3\"\)
root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_4.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b2_efix_part1_4.root\"\,\"GammaConvV1_4\"\)

hadd $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass2_1.root $OUTPUTDIR/GammaConvV1_LHC13b-pass3_1.root $OUTPUTDIR/GammaConvV1_LHC13c-pass2_1.root
hadd $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass2_2.root $OUTPUTDIR/GammaConvV1_LHC13b-pass3_2.root $OUTPUTDIR/GammaConvV1_LHC13c-pass2_2.root
hadd $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass2_3.root $OUTPUTDIR/GammaConvV1_LHC13b-pass3_3.root $OUTPUTDIR/GammaConvV1_LHC13c-pass2_3.root
hadd $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass2_4.root $OUTPUTDIR/GammaConvV1_LHC13b-pass3_4.root $OUTPUTDIR/GammaConvV1_LHC13c-pass2_4.root