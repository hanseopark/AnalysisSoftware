#! /usr/local/bin/perl

$TRAINDIRMain="Legotrain-vAN-20160510-1/mergeFilesSTDesd";
$TRAINDIR="Legotrain-vAN-20160510-1/AddedSigFiles_mergedSTDesd";
$trainConfigMBA="234";
$trainConfigAddSigPi0A="Pi0_235";
$trainConfigAddSigEtaA="Eta_235";

# system("echo $PWD");
$outputMergedBase="/home/admin1/leardini/GridOutput/PbPb/$TRAINDIR";
$outputMergedBaseOrg="/home/admin1/leardini/GridOutput/PbPb/$TRAINDIRMain";
system("echo $outputMergedBase");
system("mkdir -p $outputMergedBase");
# system("cp $outputMergedBaseOrg/*_$trainConfigMBData.* $outputMergedBase/");
 system("cp $outputMergedBaseOrg/*_$trainConfigMBA.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigMBB.* $outputMergedBase/");
# # system("cp $outputMergedBaseOrg/*_$trainConfigMBC.* $outputMergedBase/");
 system("cp $outputMergedBaseOrg/*_$trainConfigAddSigPi0A.* $outputMergedBase/");
 system("cp $outputMergedBaseOrg/*_$trainConfigAddSigEtaA.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigAddSigB.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigAddSigC.* $outputMergedBase/");

system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC11h-pass2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1b_MC_A.root");
# #system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1a_LHC14a1b_MC_A.root");
#
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Pi0_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Pi0_A.root");
# #system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_LHC14a1b_MC_Pi0_A.root");
#
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Eta_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Eta_A.root");
# #system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_LHC14a1b_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_LHC14a1b_MC_Eta_A.root");


system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC11h-pass2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC14a1b_MC_A.root");

system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Pi0_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigAddSigPi0A.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Pi0_A.root");

system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1a_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Eta_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC14a1b_$trainConfigAddSigEtaA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Eta_A.root");
