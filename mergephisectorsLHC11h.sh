#! /bin/bash

BASEDIR=/home/admin1/leardini/GridOutput/PbPb
TRAINDIR=Legotrain-vAN-20160510-1;
OUTPUTDIR=$BASEDIR/$TRAINDIR;
OUTPUTDIRmerged=$BASEDIR/$TRAINDIR/mergeFilesSTDesd;
mkdir -p $OUTPUTDIRmerged

standardData=234;
standardPhiData=236;
standard=234;
standardPhi=236;
added=235;
addedPhi=237;

# standardData=226;
# standardPhiData=228;
# standard=226;
# standardPhi=228;
# added=227;
# addedPhi=229;

#std on ESD with tag of 10th may
LHC11hData=222_20160510-1826; #---> list 7
LHC11hDataPhi=222_20160510-1826; #---> list 8

LHC14a1awithEta=271_20160510-1833;
LHC14a1bwithEta=273_20160510-1834;
LHC14a1awithPi0=272_20160510-1844;
LHC14a1bwithPi0=274_20160510-1836;

#std RP after TOF fix
# LHC11hData=306_20160429-1406;
# LHC11hDataPhi=307_20160429-1328;
#
# LHC14a1awithEta=731_20160429-1109;
# LHC14a1bwithEta=733_20160429-1110;
# LHC14a1awithPi0=732_20160502-1855;
# LHC14a1bwithPi0=734_20160429-1518;

#RP backgorund - pre TOF fix
# LHC11hData=298_20160418-1349;
# LHC11hDataPhi=299_20160418-1350;
#
# LHC14a1awithEta=698_20160418-1611;
# LHC14a1bwithEta=700_20160418-1352;
# LHC14a1awithPi0=699_20160418-1352;
# LHC14a1bwithPi0=701_20160418-1353;


echo "For traingconfig data:" $standardData
# cat $OUTPUTDIR/GA_PbPb_AOD-$LHC11hData/merge_runlist_7/CutSelection_LHC11h_$standardData.log
#
# ###################################### LHC14a1a
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standardPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$addedPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPi0/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$addedPhi.root
#
#
# ####################################### LHC14a1b
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standardPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithEta/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$addedPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPi0/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$addedPhi.root
#
#
# ##########################################   DATA   ###########################################
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardData.root $OUTPUTDIR/GA_PbPb_AOD-$LHC11hData/merge_runlist_7/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardData.root $OUTPUTDIR/GA_PbPb_AOD-$LHC11hDataPhi/merge/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardPhiData.root
#



cat $OUTPUTDIR/GA_PbPb-$LHC11hData/merge_runlist_7/CutSelection_LHC11h_$standardData.log

###################################### esdLHC14a1a
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithPi0/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1a_$standardPhi.root
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1a_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithEta/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1a_$addedPhi.root
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1a_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithPi0/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1a_$addedPhi.root


####################################### esdLHC14a1b
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithPi0/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1b_$standardPhi.root
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1b_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithEta/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithEta/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1b_$addedPhi.root
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1b_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithPi0/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithPi0/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1b_$addedPhi.root


##########################################   esdDATA   ###########################################
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardData.root $OUTPUTDIR/GA_PbPb-$LHC11hData/merge_runlist_7/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardData.root $OUTPUTDIR/GA_PbPb-$LHC11hDataPhi/merge_runlist_8/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardPhiData.root

