#!/usr/local/bin/perl


# $TRAINDIRMain="Legotrain-v5-05-60-AN-20140123";
# $TRAINDIR="Legotrain-v5-05-60-AN-20140123_1";
# $trainConfigMBA="9";
# $trainConfigAddSigA="10"; 
# $trainConfigMBB="11";
# $trainConfigAddSigB="12"; 
# $trainConfigMBC="13";
# $trainConfigAddSigC="14"; 

# $TRAINDIR="Legotrain-v5-05-60-AN-20140123_2";
# $trainConfigMBA="15";
# $trainConfigAddSigA="16"; 
# $trainConfigMBB="17";
# $trainConfigAddSigB="18"; 
# $trainConfigMBC="19";
# $trainConfigAddSigC="20"; 

# $TRAINDIRMain="Legotrain-v5-05-63-AN-SystematicErr";
# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_1";
# $trainConfigMBA="1";
# $trainConfigAddSigA="2"; 
# $trainConfigMBB="3";
# $trainConfigAddSigB="4"; 
# $trainConfigMBC="5";
# $trainConfigAddSigC="6"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_2";
# $trainConfigMBA="7";
# $trainConfigAddSigA="8"; 
# $trainConfigMBB="9";
# $trainConfigAddSigB="10"; 
# $trainConfigMBC="11";
# $trainConfigAddSigC="12"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_3";
# $trainConfigMBA="13";
# $trainConfigAddSigA="14"; 
# $trainConfigMBB="15";
# $trainConfigAddSigB="16"; 
# $trainConfigMBC="17";
# $trainConfigAddSigC="18"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_4";
# $trainConfigMBA="19";
# $trainConfigAddSigA="20"; 
# $trainConfigMBB="21";
# $trainConfigAddSigB="22"; 
# $trainConfigMBC="17";
# $trainConfigAddSigC="18"; 

# $TRAINDIRMain="Legotrain-v5-05-63-AN-SystematicErrDiffEta";
# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErrDiffEta";
# $trainConfigMBData="1_3";
# $trainConfigMBA="1";
# $trainConfigAddSigA="2"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_5";
# $trainConfigMBData="1_3";
# $trainConfigMBA="1";
# $trainConfigAddSigA="2"; 
# $trainConfigMBB="3";
# $trainConfigAddSigB="4"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_6";
# $trainConfigMBData="5_7";
# $trainConfigMBA="5";
# $trainConfigAddSigA="6"; 
# $trainConfigMBB="7";
# $trainConfigAddSigB="8"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_7";
# $trainConfigMBData="9_11";
# $trainConfigMBA="9";
# $trainConfigAddSigA="10"; 
# $trainConfigMBB="11";
# $trainConfigAddSigB="12"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_8";
# $trainConfigMBData="13_15";
# $trainConfigMBA="13";
# $trainConfigAddSigA="14"; 
# $trainConfigMBB="15";
# $trainConfigAddSigB="16"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_9";
# $trainConfigMBData="17_19";
# $trainConfigMBA="17";
# $trainConfigAddSigA="18"; 
# $trainConfigMBB="19";
# $trainConfigAddSigB="20"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErr_10";
# $trainConfigMBData="21_43";
# $trainConfigMBA="21";
# $trainConfigAddSigA="22"; 
# $trainConfigMBB="43";
# $trainConfigAddSigB="44"; 

# $TRAINDIRMain="Legotrain-v5-05-63-AN-SystematicErrEtaMB";
# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErrEta_1";
# $trainConfigMBData="133_135";
# $trainConfigMBA="133";
# $trainConfigAddSigA="134"; 
# $trainConfigMBB="135";
# $trainConfigAddSigB="136"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErrEta_2";
# $trainConfigMBData="137_139";
# $trainConfigMBA="137";
# $trainConfigAddSigA="138"; 
# $trainConfigMBB="139";
# $trainConfigAddSigB="140"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErrEta_3";
# $trainConfigMBData="141_143";
# $trainConfigMBA="141";
# $trainConfigAddSigA="142"; 
# $trainConfigMBB="143";
# $trainConfigAddSigB="144"; 
# 
# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErrEta_4";
# $trainConfigMBData="145_147";
# $trainConfigMBA="145";
# $trainConfigAddSigA="146"; 
# $trainConfigMBB="147";
# $trainConfigAddSigB="148"; 

# $TRAINDIR="Legotrain-v5-05-63-AN-SystematicErrEta_5";
# $trainConfigMBData="149_151";
# $trainConfigMBA="149";
# $trainConfigAddSigA="150"; 
# $trainConfigMBB="151";
# $trainConfigAddSigB="152"; 

# $TRAINDIR="Legotrain-vAN-20140410";
# $trainConfigMBData="173";
# $trainConfigMBA="173";
# $trainConfigAddSigA="174"; 
# # $trainConfigMBB="151";
# # $trainConfigAddSigB="152"; 

# $TRAINDIR="Legotrain-vAN-20140420";
# $trainConfigMBData="175";
# $trainConfigMBA="175";
# $trainConfigAddSigA="176"; 
# $trainConfigMBB="177";
# $trainConfigAddSigB="178"; 


$TRAINDIR="Legotrain-vAN-20141124-ConvV1";
$trainConfigMBData="5";
$trainConfigMBA="5";
$trainConfigAddSigA="6"; 
$trainConfigMBB="177";
$trainConfigAddSigB="178"; 


# system("echo $PWD");
$outputMergedBase="/home/fbock/Photon/Grid/OutputLegoTrains/pPb/$TRAINDIR";
$outputMergedBaseOrg="/home/fbock/Photon/Grid/OutputLegoTrains/pPb/$TRAINDIRMain";
system("echo $outputMergedBase");
# system("mkdir -p $outputMergedBase");
# system("cp $outputMergedBaseOrg/*_$trainConfigMBData.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigMBA.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigMBB.* $outputMergedBase/");
# # system("cp $outputMergedBaseOrg/*_$trainConfigMBC.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigAddSigA.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigAddSigB.* $outputMergedBase/");
# # system("cp $outputMergedBaseOrg/*_$trainConfigAddSigC.* $outputMergedBase/");

system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

# system("ln -sf $outputMergedBase/GammaConvV1_LHC13b-pass3_LHC13c-pass3_$trainConfigMBData.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_LHC13b-pass3_LHC13c-pass3_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_LHC13b-pass3_LHC13c-pass3_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1Data_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_LHC13b-pass3_LHC13c-pass3_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1Data_C.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13e7_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_HIJING_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13e7_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_HIJING_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13e7_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_HIJING_MC_C.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_DPMJET_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_DPMJET_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_DPMJET_MC_C.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13e7_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConvV1_HIJING_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13e7_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConvV1_HIJING_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC13e7_$trainConfigAddSigC.root $outputMergedBase/mergedAddSignal/GammaConvV1_HIJING_MC_C.root");
