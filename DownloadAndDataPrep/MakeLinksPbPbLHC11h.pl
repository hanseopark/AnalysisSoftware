#! /usr/local/bin/perl
#
# $TRAINDIRMain="Legotrain-vAN-20180731-1_withppMBW/";
# $TRAINDIR="Legotrain-vAN-20180731-1_withppMBW/AddSignal_withppMBW";
# $trainConfigMBA="162";
# $trainConfigAddSigPi0A="Pi0_163";
# $trainConfigAddSigEtaA="Eta_163";
#
# # system("echo $PWD");
# $outputMergedBase="/home/admin1/leardini/GridOutput/PbPb/$TRAINDIR";
# $outputMergedBaseOrg="/home/admin1/leardini/GridOutput/PbPb/$TRAINDIRMain";
# system("echo $outputMergedBase");
# system("mkdir -p $outputMergedBase");
# system("mkdir $outputMergedBase/mergedAddSignal");
# system("mkdir $outputMergedBase/mergedMinBias");
#
# ### Mesons linking
# system("cp $outputMergedBaseOrg/*_$trainConfigMBA.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigAddSigPi0A.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigAddSigEtaA.* $outputMergedBase/");
# ### AOD
# # system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC11h-pass2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
# # system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root");
# # system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1b_MC_A.root");
# # system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Pi0_A.root");
# # system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Pi0_A.root");
# # system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Eta_A.root");
# # system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Eta_A.root");
# # #system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_LHC14a1b_MC_Eta_A.root");
#
# ### ESD
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC11h-pass2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1b_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Pi0_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Pi0_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Eta_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Eta_A.root");


### gamma linking
$TRAINDIRMain="Legotrain-vAN-20180824-1_DPwithppMBW/";
$TRAINDIR="Legotrain-vAN-20180824-1_DPwithppMBW/AddSignalForGamma_STDwithppMBW";
$trainConfigMBA="186";
$trainConfigAddSigPi0A="313";
$trainConfigMBB="188";
$trainConfigAddSigPi0B="315";

$outputMergedBase="/home/admin1/leardini/GridOutput/PbPb/$TRAINDIR";
$outputMergedBaseOrg="/home/admin1/leardini/GridOutput/PbPb/$TRAINDIRMain";
system("echo $outputMergedBase");
system("mkdir -p $outputMergedBase");
system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

system("mv $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigAddSigPi0A.root $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_Pi0_$trainConfigAddSigPi0A.root");
system("mv $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigAddSigPi0B.root $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_Pi0_$trainConfigAddSigPi0B.root");

system("cp $outputMergedBaseOrg/*_$trainConfigMBA.* $outputMergedBase/");
system("cp $outputMergedBaseOrg/*_$trainConfigMBB.* $outputMergedBase/");
system("cp $outputMergedBaseOrg/*_$trainConfigAddSigPi0A.* $outputMergedBase/");
system("cp $outputMergedBaseOrg/*_$trainConfigAddSigPi0B.* $outputMergedBase/");


system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC11h-pass2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC11h-pass2_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1Data_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1b_MC_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Pi0_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigAddSigPi0B.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Pi0_B.root");
