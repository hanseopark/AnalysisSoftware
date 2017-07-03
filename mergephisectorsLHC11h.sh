#! /bin/bash

BASEDIR=/home/admin1/leardini/GridOutput/PbPb
TRAINDIR=Legotrain-vAN-20170626-1;
OUTPUTDIR=$BASEDIR/$TRAINDIR;
OUTPUTDIRmerged=$BASEDIR/$TRAINDIR;
mkdir -p $OUTPUTDIRmerged

standard=226;
standardPhi=227;
added=228;
addedPhi=229;

LHC11hData=299_20170626-1844; 
LHC11hDataPhi=299_20170626-1844; 

LHC14a1a=554_20170626-1859;
LHC14a1b=555_20170627-1029;
LHC14a1awithPhi=554_20170626-1859;
LHC14a1bwithPhi=555_20170627-1029;


echo "For traingconfig data:" $standard
cat $OUTPUTDIR/GA_PbPb-$LHC11hData/merge_runlist_7/CutSelection_LHC11h_$standard.log

##########################################   esdDATA   ###########################################
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_LHC11h-pass2_$standard.root $OUTPUTDIR/GA_PbPb-$LHC11hData/merge_runlist_7/GammaConvV1_GA_PbPb_LHC11h-pass2_$standard.root $OUTPUTDIR/GA_PbPb-$LHC11hDataPhi/merge_runlist_8/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardPhi.root

###################################### esdLHC14a1a
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1a/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithPhi/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1a_$standardPhi.root
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1a_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1a/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1a/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1a_$addedPhi.root
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1a_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithPhi/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1awithPhi/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1a_$addedPhi.root

####################################### esdLHC14a1b
hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_LHC11h-pass2_$added.root $OUTPUTDIR/GA_PbPb-$LHC11hData/merge_runlist_7/GammaConvV1_GA_PbPb_LHC11h-pass2_$added.root $OUTPUTDIR/GA_PbPb-$LHC11hDataPhi/merge_runlist_8/GammaConvV1_GA_PbPb_LHC11h-pass2_$addedPhi.root

hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1b/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithPhi/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1b_$addedPhi.root
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1b_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1b/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1b/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1b_$addedPhi.root
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_LHC14a1b_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithPhi/merge_runlist_7/GammaConvV1_GA_PbPb_MC_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC-$LHC14a1bwithPhi/merge_runlist_8/GammaConvV1_GA_PbPb_MC_LHC14a1b_$addedPhi.root



# cat $OUTPUTDIR/GA_PbPb_AOD-$LHC11hData/merge_runlist_7/CutSelection_LHC11h_$standard.log
#
# ###################################### LHC14a1a
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1a/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPhi/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$standardPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1a/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1a/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$addedPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPhi/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1awithPhi/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$addedPhi.root
#
#
# ####################################### LHC14a1b
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1b/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standard.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPhi/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$standardPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_Eta_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1b/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1b/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$addedPhi.root
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_Pi0_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPhi/merge_runlist_6/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$added.root $OUTPUTDIR/GA_PbPb_MC_AOD-$LHC14a1bwithPhi/merge_runlist_7/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$addedPhi.root
#
#
# ##########################################   DATA   ###########################################
#
# hadd -f $OUTPUTDIRmerged/GammaConvV1_GA_PbPb_LHC11h-pass2_$standard.root $OUTPUTDIR/GA_PbPb_AOD-$LHC11hData/merge_runlist_7/GammaConvV1_GA_PbPb_LHC11h-pass2_$standard.root $OUTPUTDIR/GA_PbPb_AOD-$LHC11hDataPhi/merge/GammaConvV1_GA_PbPb_LHC11h-pass2_$standardPhi.root
#



