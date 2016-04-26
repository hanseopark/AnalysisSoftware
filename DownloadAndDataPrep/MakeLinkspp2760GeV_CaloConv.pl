 #!/usr/local/bin/perl

$TRAINDIR="Legotrain-vAN-20150319-ConvCalo-1stIterationFix2";
$trainConfigMBA="12";
$trainConfigAddSigA="12"; 
$trainConfigMBB="12";
$trainConfigAddSigB="12"; 
$AddName="WSDD";

$outputMergedBase="/home/fbock/Photon/Grid/OutputLegoTrains/pp/$TRAINDIR";
system("echo $outputMergedBase");

system("mkdir $outputMergedBase/mergedAddSignal");
system("mkdir $outputMergedBase/mergedMinBias");

system("ln -sf $outputMergedBase/GammaConvCalo_LHC11a-pass3\_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvCalo_LHC11a-pass3\_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_Data-$AddName\_B.root");

system("ln -sf $outputMergedBase/GammaConvCalo_MC_LHC12f1a_LHC12f1b_LHC12i3\_$trainConfigMBA.root $outputMergedBase/mergedMinBias/GammaConv_MC-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvCalo_MC_LHC12f1a_LHC12f1b_LHC12i3\_$trainConfigMBB.root $outputMergedBase/mergedMinBias/GammaConv_MC-$AddName\_B.root");

system("ln -sf $outputMergedBase/GammaConvCalo_MC_LHC12i3\_$trainConfigAddSigA.root $outputMergedBase/mergedAddSignal/GammaConv_MC-$AddName\_A.root");
system("ln -sf $outputMergedBase/GammaConvCalo_MC_LHC12i3\_$trainConfigAddSigB.root $outputMergedBase/mergedAddSignal/GammaConv_MC-$AddName\_B.root");

