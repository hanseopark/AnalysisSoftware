 #!/usr/local/bin/perl


# $TRAINDIR="Legotrain-vAN-20150101-ConvWeight1stIte";
# $trainConfigMBA="12";
# $trainConfigAddSigA="13"; 
# $trainConfigMBB="60";
# $trainConfigAddSigB="60"; 

# $TRAINDIR="Legotrain-vAN-20150101-ConvWeight2ndIte";
# $trainConfigMBA="12";
# $trainConfigAddSigA="13"; 
# $trainConfigMBB="60";
# $trainConfigAddSigB="60"; 

# $TRAINDIR="Legotrain-vAN-20150101-ConvWeight3rdIte";
# $trainConfigMBA="12";
# $trainConfigAddSigA="13"; 
# $trainConfigMBB="60";
# $trainConfigAddSigB="60"; 

$TRAINDIR="Legotrain-vAN-20150101-ConvWeight4thIte";
$trainConfigMBA="12";
$trainConfigAddSigA="13"; 
$trainConfigMBB="60";
$trainConfigAddSigB="60"; 

$outputMergedBase="/home/fbock/Photon/Grid/OutputLegoTrains/pp8TeV/$TRAINDIR";

system("echo $outputMergedBase");

system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

system("ln -sf $outputMergedBase/GammaConvV1_LHC12abcdfghi-pass1_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_Data_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_LHC12abcdfghi-pass1_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_Data_B.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC14e2a_LHC14e2b_LHC14e2c_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC14e2a_LHC14e2b_LHC14e2c_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_MC_B.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC14e2b_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConv_MC_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC14e2b_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConv_MC_B.root");
