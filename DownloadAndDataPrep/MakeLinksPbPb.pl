#!/usr/local/bin/perl

#private train: weigting first iteration V0 Mult cuts pi0
# $TRAINDIR="Legotrain-v5-05-11-AN-20130826";
# $trainConfigMBA="1";
# $trainConfigAddSigA="5"; 
# $trainConfigMBB="2";
# $trainConfigAddSigB="6"; 

#private train: weigting first iteration TPC Mult cuts pi0
# $TRAINDIR="Legotrain-v5-05-13-AN-20130902";
# $trainConfigMBA="12";
# $trainConfigAddSigA="16"; 
# $trainConfigMBB="13";
# $trainConfigAddSigB="17"; 

#private train: weigting first iteration TPC Mult cuts pi0
# $TRAINDIR="Legotrain-v5-05-13-AN-20130904";
# $trainConfigMBA="14";
# $trainConfigAddSigA="20"; 
# $trainConfigMBB="15";
# $trainConfigAddSigB="21"; 

#private train: weigting first iteration TPC Mult cuts pi0
# $TRAINDIR="Legotrain-v5-05-13-AN-20130903";
# $trainConfigMBA="3";
# $trainConfigAddSigA="10"; 
# $trainConfigMBB="4";
# $trainConfigAddSigB="11"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130907";
# $trainConfigMBA="7";
# $trainConfigAddSigA="10"; 
# $trainConfigMBB="8";
# $trainConfigAddSigB="11"; 
# $trainConfigMBC="9";
# $trainConfigAddSigC="12"; 

# $TRAINDIRMain="Legotrain-v5-05-15-AN-20130909";
# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_1";
# $trainConfigMBA="13";
# $trainConfigAddSigA="16"; 
# $trainConfigMBB="14";
# $trainConfigAddSigB="17"; 
# $trainConfigMBC="15";
# $trainConfigAddSigC="18"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_2";
# $trainConfigMBA="19";
# $trainConfigAddSigA="22"; 
# $trainConfigMBB="20";
# $trainConfigAddSigB="23"; 
# $trainConfigMBC="21";
# $trainConfigAddSigC="24"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_3";
# $trainConfigMBA="25";
# $trainConfigAddSigA="28"; 
# $trainConfigMBB="26";
# $trainConfigAddSigB="29"; 
# $trainConfigMBC="27";
# $trainConfigAddSigC="30"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_4";
# $trainConfigMBA="31";
# $trainConfigAddSigA="34"; 
# $trainConfigMBB="32";
# $trainConfigAddSigB="35"; 
# $trainConfigMBC="33";
# $trainConfigAddSigC="36"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_5";
# $trainConfigMBA="37";
# $trainConfigAddSigA="40"; 
# $trainConfigMBB="38";
# $trainConfigAddSigB="41"; 
# $trainConfigMBC="39";
# $trainConfigAddSigC="42"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_6";
# $trainConfigMBA="43";
# $trainConfigAddSigA="46"; 
# $trainConfigMBB="44";
# $trainConfigAddSigB="47"; 
# $trainConfigMBC="45";
# $trainConfigAddSigC="48"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_7";
# $trainConfigMBA="49";
# $trainConfigAddSigA="52"; 
# $trainConfigMBB="50";
# $trainConfigAddSigB="53"; 
# $trainConfigMBC="51";
# $trainConfigAddSigC="54"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_8";
# $trainConfigMBA="55";
# $trainConfigAddSigA="58"; 
# $trainConfigMBB="56";
# $trainConfigAddSigB="59"; 
# $trainConfigMBC="57";
# $trainConfigAddSigC="60"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_9";
# $trainConfigMBA="61";
# $trainConfigAddSigA="64"; 
# $trainConfigMBB="62";
# $trainConfigAddSigB="65"; 
# $trainConfigMBC="63";
# $trainConfigAddSigC="66"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_10";
# $trainConfigMBA="67";
# $trainConfigAddSigA="70"; 
# $trainConfigMBB="68";
# $trainConfigAddSigB="71"; 
# $trainConfigMBC="69";
# $trainConfigAddSigC="72"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_11";
# $trainConfigMBA="73";
# $trainConfigAddSigA="76"; 
# $trainConfigMBB="74";
# $trainConfigAddSigB="77"; 
# $trainConfigMBC="75";
# $trainConfigAddSigC="78"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_12";
# $trainConfigMBA="79";
# $trainConfigAddSigA="82"; 
# $trainConfigMBB="80";
# $trainConfigAddSigB="83"; 
# $trainConfigMBC="81";
# $trainConfigAddSigC="84"; 

# $TRAINDIR="Legotrain-v5-05-15-AN-20130909_13";
# $trainConfigMBA="85";
# $trainConfigAddSigA="88"; 
# $trainConfigMBB="86";
# $trainConfigAddSigB="89"; 
# $trainConfigMBC="87";
# $trainConfigAddSigC="90"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_1";
# $trainConfigMBA="1";
# $trainConfigAddSigA="4"; 
# $trainConfigMBB="2";
# $trainConfigAddSigB="5"; 
# $trainConfigMBC="3";
# $trainConfigAddSigC="6"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_2";
# $trainConfigMBA="7";
# $trainConfigAddSigA="10"; 
# $trainConfigMBB="8";
# $trainConfigAddSigB="11"; 
# $trainConfigMBC="9";
# $trainConfigAddSigC="12"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_3";
# $trainConfigMBA="13";
# $trainConfigAddSigA="16"; 
# $trainConfigMBB="14";
# $trainConfigAddSigB="17"; 
# $trainConfigMBC="15";
# $trainConfigAddSigC="18"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_4";
# $trainConfigMBA="19";
# $trainConfigAddSigA="22"; 
# $trainConfigMBB="20";
# $trainConfigAddSigB="23"; 
# $trainConfigMBC="21";
# $trainConfigAddSigC="24"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_5";
# $trainConfigMBA="25";
# $trainConfigAddSigA="28"; 
# $trainConfigMBB="26";
# $trainConfigAddSigB="29"; 
# $trainConfigMBC="27";
# $trainConfigAddSigC="30"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_6";
# $trainConfigMBA="31";
# $trainConfigAddSigA="34"; 
# $trainConfigMBB="32";
# $trainConfigAddSigB="35"; 
# $trainConfigMBC="33";
# $trainConfigAddSigC="36"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_7";
# $trainConfigMBA="37";
# $trainConfigAddSigA="40"; 
# $trainConfigMBB="38";
# $trainConfigAddSigB="41"; 
# $trainConfigMBC="39";
# $trainConfigAddSigC="42"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_8";
# $trainConfigMBA="43";
# $trainConfigAddSigA="46"; 
# $trainConfigMBB="44";
# $trainConfigAddSigB="47"; 
# $trainConfigMBC="45";
# $trainConfigAddSigC="48"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_9";
# $trainConfigMBA="49";
# $trainConfigAddSigA="52"; 
# $trainConfigMBB="50";
# $trainConfigAddSigB="53"; 
# $trainConfigMBC="51";
# $trainConfigAddSigC="54"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_10";
# $trainConfigMBA="55";
# $trainConfigAddSigA="58"; 
# $trainConfigMBB="56";
# $trainConfigAddSigB="59"; 
# $trainConfigMBC="57";
# $trainConfigAddSigC="60"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_11";
# $trainConfigMBA="61";
# $trainConfigAddSigA="64"; 
# $trainConfigMBB="62";
# $trainConfigAddSigB="65"; 
# $trainConfigMBC="63";
# $trainConfigAddSigC="66"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_12";
# $trainConfigMBA="67";
# $trainConfigAddSigA="70"; 
# $trainConfigMBB="68";
# $trainConfigAddSigB="71"; 
# $trainConfigMBC="69";
# $trainConfigAddSigC="72"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_13";
# $trainConfigMBA="73";
# $trainConfigAddSigA="76"; 
# $trainConfigMBB="74";
# $trainConfigAddSigB="77"; 
# $trainConfigMBC="75";
# $trainConfigAddSigC="78"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_14";
# $trainConfigMBA="79";
# $trainConfigAddSigA="82"; 
# $trainConfigMBB="80";
# $trainConfigAddSigB="83"; 
# $trainConfigMBC="81";
# $trainConfigAddSigC="84"; 

# $TRAINDIR="Legotrain-v5-05-16-AN-20130910_15";
# $trainConfigMBA="85";
# $trainConfigAddSigA="88"; 
# $trainConfigMBB="86";
# $trainConfigAddSigB="89"; 
# $trainConfigMBC="87";
# $trainConfigAddSigC="90"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_1";
# $trainConfigMBA="1";
# $trainConfigAddSigA="4"; 
# $trainConfigMBB="2";
# $trainConfigAddSigB="5"; 
# $trainConfigMBC="3";
# $trainConfigAddSigC="6"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_2";
# $trainConfigMBA="7";
# $trainConfigAddSigA="10"; 
# $trainConfigMBB="8";
# $trainConfigAddSigB="11"; 
# $trainConfigMBC="9";
# $trainConfigAddSigC="12"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_3";
# $trainConfigMBA="13";
# $trainConfigAddSigA="16"; 
# $trainConfigMBB="14";
# $trainConfigAddSigB="17"; 
# $trainConfigMBC="15";
# $trainConfigAddSigC="18"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_4";
# $trainConfigMBA="19";
# $trainConfigAddSigA="22"; 
# $trainConfigMBB="20";
# $trainConfigAddSigB="23"; 
# $trainConfigMBC="21";
# $trainConfigAddSigC="24"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_5";
# $trainConfigMBA="25";
# $trainConfigAddSigA="28"; 
# $trainConfigMBB="26";
# $trainConfigAddSigB="29"; 
# $trainConfigMBC="27";
# $trainConfigAddSigC="30"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_6";
# $trainConfigMBA="31";
# $trainConfigAddSigA="34"; 
# $trainConfigMBB="32";
# $trainConfigAddSigB="35"; 
# $trainConfigMBC="33";
# $trainConfigAddSigC="36"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_7";
# $trainConfigMBA="37";
# $trainConfigAddSigA="40"; 
# $trainConfigMBB="38";
# $trainConfigAddSigB="41"; 
# $trainConfigMBC="39";
# $trainConfigAddSigC="42"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_8";
# $trainConfigMBA="43";
# $trainConfigAddSigA="46"; 
# $trainConfigMBB="44";
# $trainConfigAddSigB="47"; 
# $trainConfigMBC="45";
# $trainConfigAddSigC="48"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_9";
# $trainConfigMBA="49";
# $trainConfigAddSigA="52"; 
# $trainConfigMBB="50";
# $trainConfigAddSigB="53"; 
# $trainConfigMBC="51";
# $trainConfigAddSigC="54"; 

# $TRAINDIR="Legotrain-v5-05-18-AN-20130919_10";
# $trainConfigMBA="55";
# $trainConfigAddSigA="58"; 
# $trainConfigMBB="56";
# $trainConfigAddSigB="59"; 
# $trainConfigMBC="57";
# $trainConfigAddSigC="60"; 

# $TRAINDIRMain="Legotrain-v5-05-21-AN-20130928";
# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_1";
# $trainConfigMBData="1_2";
# $trainConfigMBA="1";
# $trainConfigAddSigA="3"; 
# $trainConfigMBB="2";
# $trainConfigAddSigB="4"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_2";
# $trainConfigMBData="5_6";
# $trainConfigMBA="5";
# $trainConfigAddSigA="7"; 
# $trainConfigMBB="6";
# $trainConfigAddSigB="8"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_3";
# $trainConfigMBData="9_10";
# $trainConfigMBA="9";
# $trainConfigAddSigA="11"; 
# $trainConfigMBB="10";
# $trainConfigAddSigB="12"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_4";
# $trainConfigMBData="13_14";
# $trainConfigMBA="13";
# $trainConfigAddSigA="15"; 
# $trainConfigMBB="14";
# $trainConfigAddSigB="16"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_5";
# $trainConfigMBData="17_18";
# $trainConfigMBA="17";
# $trainConfigAddSigA="19"; 
# $trainConfigMBB="18";
# $trainConfigAddSigB="20"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_6";
# $trainConfigMBData="21_22";
# $trainConfigMBA="21";
# $trainConfigAddSigA="23"; 
# $trainConfigMBB="22";
# $trainConfigAddSigB="24"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_7";
# $trainConfigMBData="25_26";
# $trainConfigMBA="25";
# $trainConfigAddSigA="27"; 
# $trainConfigMBB="26";
# $trainConfigAddSigB="28"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_8";
# $trainConfigMBData="29_30";
# $trainConfigMBA="29";
# $trainConfigAddSigA="31"; 
# $trainConfigMBB="30";
# $trainConfigAddSigB="32"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_9";
# $trainConfigMBData="33_34";
# $trainConfigMBA="33";
# $trainConfigAddSigA="35"; 
# $trainConfigMBB="34";
# $trainConfigAddSigB="36"; 

# $TRAINDIR="Legotrain-v5-05-21-AN-20130928_10";
# $trainConfigMBData="37_38";
# $trainConfigMBA="37";
# $trainConfigAddSigA="39"; 
# $trainConfigMBB="38";
# $trainConfigAddSigB="40"; 

# $TRAINDIR="newStandardCuts";
# $trainConfigMBData="13";
# $trainConfigMBA="13";
# $trainConfigAddSigA="15"; 
# $trainConfigMBB="14";
# $trainConfigAddSigB="16"; 

# $TRAINDIR="Legotrain-v5-05-29-AN-20131024_1";
# # $trainConfigMBData="13";
# $trainConfigMBA="17";
# $trainConfigAddSigA="19"; 
# $trainConfigMBB="18";
# $trainConfigAddSigB="20"; 

# $TRAINDIR="Legotrain-v5-05-29-AN-20131024_2";
# # $trainConfigMBData="13";
# $trainConfigMBA="21";
# $trainConfigAddSigA="23"; 
# $trainConfigMBB="22";
# $trainConfigAddSigB="24"; 

# $TRAINDIRMain="PhotonQualityTest";
# $TRAINDIR="PhotonQualityTest_1";
# $trainConfigMBData="17_18";
# $trainConfigMBA="17";
# $trainConfigAddSigA="19"; 
# $trainConfigMBB="18";
# $trainConfigAddSigB="20"; 

# $TRAINDIRMain="PhotonQualityTest";
# $TRAINDIR="PhotonQualityTest_2";
# $trainConfigMBData="21_22";
# $trainConfigMBA="21";
# $trainConfigAddSigA="23"; 
# $trainConfigMBB="22";
# $trainConfigAddSigB="24"; 

# $TRAINDIRMain="PhotonQualityTest";
# $TRAINDIR="PhotonQualityTest_3";
# $trainConfigMBData="25_26";
# $trainConfigMBA="25";
# $trainConfigAddSigA="27"; 
# $trainConfigMBB="26";
# $trainConfigAddSigB="28"; 

# $TRAINDIRMain="PhotonQualityTest";
# $TRAINDIR="PhotonQualityTest_4";
# $trainConfigMBData="13_14";
# $trainConfigMBA="13";
# $trainConfigAddSigA="15"; 
# $trainConfigMBB="14";
# $trainConfigAddSigB="16"; 

# $TRAINDIRMain="PhotonQualityTest";
# $TRAINDIR="PhotonQualityTest_5";
# $trainConfigMBData="29_30";
# $trainConfigMBA="29";
# $trainConfigAddSigA="31"; 
# $trainConfigMBB="30";
# $trainConfigAddSigB="32"; 

$TRAINDIRMain="PhotonQualityTest";
$TRAINDIR="PhotonQualityTest_6";
$trainConfigMBData="33_34";
$trainConfigMBA="33";
$trainConfigAddSigA="35"; 
$trainConfigMBB="34";
$trainConfigAddSigB="36"; 


# system("echo $PWD");
$outputMergedBase="/home/fbock/Photon/Grid/OutputLegoTrains/PbPb/$TRAINDIR";
$outputMergedBaseOrg="/home/fbock/Photon/Grid/OutputLegoTrains/PbPb/$TRAINDIRMain";
system("echo $outputMergedBase");
system("mkdir -p $outputMergedBase");
system("cp $outputMergedBaseOrg/*_$trainConfigMBData.* $outputMergedBase/");
system("cp $outputMergedBaseOrg/*_$trainConfigMBA.* $outputMergedBase/");
system("cp $outputMergedBaseOrg/*_$trainConfigMBB.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigMBC.* $outputMergedBase/");
system("cp $outputMergedBaseOrg/*_$trainConfigAddSigA.* $outputMergedBase/");
system("cp $outputMergedBaseOrg/*_$trainConfigAddSigB.* $outputMergedBase/");
# system("cp $outputMergedBaseOrg/*_$trainConfigAddSigC.* $outputMergedBase/");

system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC10h-pass2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC10h-pass2_$trainConfigMBData.root $outputMergedBase/mergedMinBias/GammaConvV1Data_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2b_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2b_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_A.root");

# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC10h-pass2_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1Data_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC10h-pass2_$trainConfigMBData.root $outputMergedBase/mergedMinBias/GammaConvV1Data_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_MC_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2b_MC_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2b_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2b_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_B.root");

# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_LHC10h-pass2_$trainConfigMBData.root $outputMergedBase/mergedMinBias/GammaConvV1Data_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2b_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2b_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2b_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_$trainConfigMBC.root $outputMergedBase/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_C.root");

system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2b_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2b_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2b_MC_A.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_A.root");

system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_MC_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2b_MC_B.root");
system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2b_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2b_MC_B.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_B.root");

# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$trainConfigAddSigC.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$trainConfigAddSigC.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2b_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_$trainConfigAddSigC.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_$trainConfigAddSigC.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2b_$trainConfigAddSigC.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2b_MC_C.root");
# system("ln -sf $outputMergedBase/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_$trainConfigAddSigC.root $outputMergedBase/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_C.root");

