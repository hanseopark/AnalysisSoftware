 #!/usr/local/bin/perl

# $TRAINDIR="Legotrain-vAN-20150101-ConvWeighting1stIte";
# $trainConfigMBA="12";
# $trainConfigAddSigA="13"; 
# $trainConfigMBB="60";
# $trainConfigAddSigB="60"; 
# $AddName="WOSDD";

# $TRAINDIR="Legotrain-vAN-20150101-ConvWeighting2ndIte";
# $trainConfigMBA="12";
# $trainConfigAddSigA="13"; 
# $trainConfigMBB="60";
# $trainConfigAddSigB="60"; 
# $AddName="WOSDD";

# $TRAINDIR="Legotrain-vAN-20150101-ConvWeighting3rdIte";
# $trainConfigMBA="12";
# $trainConfigAddSigA="13"; 
# $trainConfigMBB="60";
# $trainConfigAddSigB="60"; 
# $AddName="WOSDD";

$TRAINDIR="Legotrain-vAN-20150101-ConvWeighting4thIte";
$trainConfigMBA="12";
$trainConfigAddSigA="13"; 
$trainConfigMBB="60";
$trainConfigAddSigB="60"; 
$AddName="WOSDD";

$outputMergedBase="/home/fbock/Photon/Grid/OutputLegoTrains/pp/$TRAINDIR";
system("echo $outputMergedBase");

system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

system("ln -sf $outputMergedBase/GammaConvV1_LHC11a-pass4-$AddName\_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_LHC11a-pass4-$AddName\_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_B.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12f1a_LHC12f1b_LHC12i3-$AddName\_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_MC-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12f1a_LHC12f1b_LHC12i3-$AddName\_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_MC-$AddName\_B.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12i3-$AddName\_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConv_MC-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12i3-$AddName\_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConv_MC-$AddName\_B.root");

