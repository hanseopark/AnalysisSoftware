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

#$TRAINDIR="Legotrain-vAN-20150101-ConvWeighting4thIte";
#$trainConfigMBA="12";
#$trainConfigAddSigA="13"; 
#$trainConfigMBB="60";
#$trainConfigAddSigB="60"; 
#$AddName="WOSDD";


#Efficiency study for pp 2.76TeV LHC11a by Hikari  
#$TRAINDIR="Legotrain-vAN-20161115-ConvTesting";# 1th Iteration 
#$trainConfigMBA="139";#trainconfig 139 = 00003113
#$trainConfigAddSigA="140";#trainconfig 140 = 00003123 
#$trainConfigMBB="139";
#$trainConfigAddSigB="140"; 
#$AddName="WSDD";

#$TRAINDIR="Legotrain-vAN-20161119-ConvTesting";# 2th Iteration
#$trainConfigMBA="139";
#$trainConfigAddSigA="140"; 
#$trainConfigMBB="139";
#$trainConfigAddSigB="140"; 
#$AddName="WSDD";

$TRAINDIR="Legotrain-vAN-20161128-ConvTesting";# 3th Iteration
$trainConfigMBA="139";
$trainConfigAddSigA="140"; 
$trainConfigMBB="139";
$trainConfigAddSigB="140"; 
$AddName="WSDD";

#$outputMergedBase="~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161115-ConvTesting";
#$outputMergedBase="~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161119-ConvTesting";
$outputMergedBase="~/Hikari/Photon/Grid/OutputLegoTrains/pp/Legotrain-vAN-20161128-ConvTesting";
system("echo $outputMergedBase");

system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

system("ln -sf $outputMergedBase/GammaConvV1_LHC11a-pass4-$AddName\_7.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_LHC11a-pass4-$AddName\_7.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_B.root");
#system("ln -sf $outputMergedBase/GammaConvV1_LHC11a-pass4-$AddName\_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_A.root");
#system("ln -sf $outputMergedBase/GammaConvV1_LHC11a-pass4-$AddName\_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_B.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12f1a_LHC12f1b_LHC12i3-$AddName\_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_MC-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12f1a_LHC12f1b_LHC12i3-$AddName\_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_MC-$AddName\_B.root");

system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12i3-$AddName\_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConv_MC-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvV1_MC_LHC12i3-$AddName\_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConv_MC-$AddName\_B.root");

