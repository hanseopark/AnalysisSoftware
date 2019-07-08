#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=1
MERGEONFASTAndWOSDD=1
SINGLERUN=1
SEPARATEON=0
MERGEONSINGLEData=1
MERGEONSINGLEMC=1
CLEANUP=1
CLEANUPMAYOR=$2
number=""
FAST="_FAST"
SEPARATEONLYConv=1
# check if train configuration has actually been given
HAVELHC13b=1
HAVELHC13c=1
HAVELHC13d=1
HAVELHC13e=1
HAVELHC13f=1
HAVELHC18j51=1
HAVELHC18j52=1
HAVELHC18j53=1
HAVELHC16q=1
HAVELHC16t=1
HAVELHC16qF=1
HAVELHC16tF=1
HAVETOBUILDData=0
HAVELHC17f2b=1
HAVELHC17f2bF=1
HAVETOBUILDLHC17f2b=0
HAVELHC17f2afix=1
HAVELHC17f2afixF=1
HAVETOBUILDLHC17f2afix=0
HAVELHC17g8a=1
HAVELHC17g8aF=1
HAVETOBUILDLHC17g8a=0

# default trainconfigurations
LHC13bData="";
LHC13cData="";
LHC13dData="";
LHC13eData="";
LHC13fData="";
LHC13beData="";
LHC18j5MC="";
LHC18j5_1MC="";
LHC18j5_2MC="";
LHC18j5_3MC="";
LHC16qData="";
LHC16tData="";
LHC16qDataFast="";
LHC16tDataFast="";
LHC16qtData="";
LHC17f2bMCMoth="";
LHC17f2bMC="";
LHC17f2bMCFast="";
LHC17f2a_fixMCMoth="";
LHC17f2a_fixMC="";
LHC17f2a_fixMCFast="";
LHC17g8aMCMoth=""
LHC17g8aMC=""
LHC17g8aMCFast=""

passNr="1";

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorageExternal/OutputLegoTrains/pPb
elif [ $1 = "fbockExt" ]; then
    BASEDIR=/media/fbock/Elements/OutputLegoTrains/pPb
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
fi

if [ $3 = "AOD" ]; then
    baseLegoData=GA_pPb_AOD
    baseLegoMC=GA_pPb_MC_AOD
    pathDataR1=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathMCR1=AOD214/PWGGA/GA_pPb_MC_AOD
    pathDataR2=pass1/AOD210/PWGGA/GA_pPb_AOD
    pathMCR2=AOD214/PWGGA/GA_pPb_MC_AOD
else
    baseLegoData=GA_pPb
    baseLegoMC=GA_pPb_MC
    pathDataR1=pass4/PWGGA/GA_pPb
    pathMCR1=PWGGA/GA_pPb_MC
    pathDataR2=pass1/PWGGA/GA_pPb
    pathMCR2=PWGGA/GA_pPb_MC
fi


# Definitition of number of slashes in your path to different depths
NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
NSlashes=`expr $NSlashesBASE + 4`
NSlashes2=`expr $NSlashes - 1`
NSlashes3=`expr $NSlashes + 1`
NSlashes4=`expr $NSlashes + 2`
echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"

# TRAINDIR=Legotrain-vAN20171005-dirGammaRun2
# woSDD (CENT)
# LHC16qtData="679"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1090";
# LHC17f2a_fixMC="1088";
#cell QA
# LHC16qtData="674"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1075";

# FAST
# LHC16qtData="681"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1089";
# LHC17f2a_fixMC="1087";
#cell QA
# LHC16qtData="673"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1076";

# WSDD (CENT)
# LHC16qtData="680"; #pass 1
# LHC16qData="child_1"; #pass 1
# LHC16tData="child_2"; #pass 1
# LHC17f2bMC="1092";
# LHC17f2a_fixMC="1091";

# TRAINDIR=Legotrain-vAN20180122-CellQAPHOSandEMC
# # woSDD (CENT) EMC
# LHC16qtData="707"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1150";
# LHC17f2a_fixMC="1148_20180124";

# FAST EMC
# LHC16qtData="708"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1151";
# LHC17f2a_fixMC="1149_20180124";

# # woSDD (CENT) PHOS
# LHC16qtData="710"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1146";
# LHC17f2a_fixMC="1144";

# FAST PHOS
# LHC16qtData="709_20180125"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1147";
# LHC17f2a_fixMC="1145";

# TRAINDIR=Legotrain-vAN20180122-CellQAPHOSandEMC
# woSDD (CENT) PHOS
# LHC16qtData="712"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1162";
# LHC17f2a_fixMC="1160";

# FAST PHOS
# LHC16qtData="711"; #pass 2
# LHC16qData="child_1"; #pass 3
# LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1163";
# LHC17f2a_fixMC="1161";

# TRAINDIR=Legotrain-vAN20180206-EMC
# woSDD (CENT) PHOS
# LHC16qtData="718"; #pass 2
# LHC16qDataFast="child_1"; #pass 3
# LHC16tDataFast="child_2"; #pass 2
# LHC16qData="child_3"; #pass 3
# LHC16tData="child_4"; #pass 2
# LHC17f2bMCMoth="1168";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";
# LHC17f2a_fixMCMoth="1167";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";

# LHC17f2bMCMoth="1171";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";
# LHC17f2a_fixMCMoth="1170";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";
# LHC17g8aMCMoth="1173"
# LHC17g8aMC="child_2"
# LHC17g8aMCFast="child_1"

# TRAINDIR=Legotrain-vAN20180220-EMCNonLin
# woSDD EMC
# LHC16qtData="720"; #pass 2
# LHC16qDataFast="child_1"; #pass 3
# LHC16tDataFast="child_2"; #pass 2
# LHC16qData="child_3"; #pass 3
# LHC16tData="child_4"; #pass 2
# LHC17f2bMCMoth="1205";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";
# LHC17f2a_fixMCMoth="1207";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";

# TRAINDIR=Legotrain-vAN20180220-EMCNonLin2
# LHC17f2bMCMoth="1211";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";
# LHC17f2a_fixMCMoth="1210";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";
# LHC17g8aMCMoth="1216"
# LHC17g8aMC="child_2"
# LHC17g8aMCFast="child_1"

# TRAINDIR=Legotrain-vAN20180322-RerunAll
# LHC16qtData="724"; #pass 2
# LHC16qtData="725"; #pass 2
# LHC16qDataFast="child_1"; #pass 3
# LHC16tDataFast="child_2"; #pass 2
# LHC16qData="child_3"; #pass 3
# LHC16tData="child_4"; #pass 2
# LHC17f2bMCMoth="1218";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";
# LHC17f2a_fixMCMoth="1217";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";
# LHC17g8aMCMoth="1216"
# LHC17g8aMC="child_2"
# LHC17g8aMCFast="child_1"

# LHC17f2a_fixMCMoth="1231";
# LHC17f2bMCMoth="1232";
# LHC17f2a_fixMCMoth="1230";
# LHC17f2bMCMoth="1218";
# LHC17f2a_fixMCMoth="1227";
# LHC17f2bMCMoth="1229";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";

# TRAINDIR=Legotrain-vAN20180410-RerunAll2
# LHC17f2a_fixMCMoth="1239";
# LHC17f2bMCMoth="1240";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";

# TRAINDIR=Legotrain-vAN20180427-AODValidationESD
# LHC16qtData="739"; #pass 2
# LHC16qDataFast="child_1"; #pass 3
# LHC16tDataFast="child_2"; #pass 2
# LHC16qData="child_3"; #pass 3
# LHC16tData="child_4"; #pass 2
# LHC17f2a_fixMCMoth="1251";
# LHC17f2bMCMoth="1252";
# LHC17f2a_fixMC="child_2";
# LHC17f2a_fixMCFast="child_1";
# LHC17f2bMC="child_2";
# LHC17f2bMCFast="child_1";

# new reco run1 pPb
TRAINDIR=20190612-QAPass4PCMEMCPHOS
LHC13beData="534"
LHC13bData="child_1"
LHC13cData="child_2"
LHC13dData="child_3"
LHC13eData="child_4"
LHC13fData="535"
LHC18j5MC="699"
LHC18j5_2MC="child_2"
LHC18j5_1MC="child_1"
LHC18j5_3MC="child_3"

OUTPUTDIR=$BASEDIR/$TRAINDIR
ALIENDIRData="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoData/"
OUTPUTDIRData=$BASEDIR/$TRAINDIR/$baseLegoData
ALIENDIRMC="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoMC/"
OUTPUTDIRMC=$BASEDIR/$TRAINDIR/$baseLegoMC
mkdir -p $OUTPUTDIR/CutSelections

# finding run1 data paths
FindCorrectTrainDirectory $LHC13bData $OUTPUTDIRData $ALIENDIRData $LHC13beData
HAVELHC13b=$tempBool
LHC13bData=$tempDir
OUTPUTDIR_LHC13b=$tempPath
echo "13b: $HAVELHC13b $LHC13bData $OUTPUTDIR_LHC13b"
FindCorrectTrainDirectory $LHC13cData $OUTPUTDIRData $ALIENDIRData $LHC13beData
HAVELHC13c=$tempBool
LHC13cData=$tempDir
OUTPUTDIR_LHC13c=$tempPath
echo "13c: $HAVELHC13c $LHC13cData $OUTPUTDIR_LHC13c"
FindCorrectTrainDirectory $LHC13dData $OUTPUTDIRData $ALIENDIRData $LHC13beData
HAVELHC13d=$tempBool
LHC13dData=$tempDir
OUTPUTDIR_LHC13d=$tempPath
echo "13d: $HAVELHC13d $LHC13dData $OUTPUTDIR_LHC13d"
FindCorrectTrainDirectory $LHC13eData $OUTPUTDIRData $ALIENDIRData $LHC13beData
HAVELHC13e=$tempBool
LHC13eData=$tempDir
OUTPUTDIR_LHC13e=$tempPath
echo "13e: $HAVELHC13e $LHC13eData $OUTPUTDIR_LHC13e"
FindCorrectTrainDirectory $LHC13fData $OUTPUTDIRData $ALIENDIRData
HAVELHC13f=$tempBool
LHC13fData=$tempDir
OUTPUTDIR_LHC13f=$tempPath
echo "13f: $HAVELHC13f $LHC13fData $OUTPUTDIR_LHC13f"

# finding run2 data paths
FindCorrectTrainDirectory $LHC16qData $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16q=$tempBool
LHC16qData=$tempDir
OUTPUTDIR_LHC16q=$tempPath
echo "16q: $HAVELHC16q $LHC16qData $OUTPUTDIR_LHC16q"
FindCorrectTrainDirectory $LHC16tData $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16t=$tempBool
LHC16tData=$tempDir
OUTPUTDIR_LHC16t=$tempPath
echo "16t: $HAVELHC16t $LHC16tData $OUTPUTDIR_LHC16t"
FindCorrectTrainDirectory $LHC16qDataFast $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16qF=$tempBool
LHC16qData=$tempDir
OUTPUTDIR_LHC16qF=$tempPath
echo "16q: $HAVELHC16qF $LHC16qDataFast $OUTPUTDIR_LHC16qF"
FindCorrectTrainDirectory $LHC16tDataFast $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16tF=$tempBool
LHC16tData=$tempDir
OUTPUTDIR_LHC16tF=$tempPath
echo "16t: $HAVELHC16tF $LHC16tDataFast $OUTPUTDIR_LHC16tF"

# finding run1 MC path
FindCorrectTrainDirectory $LHC18j5_1MC $OUTPUTDIRMC $ALIENDIRMC $LHC18j5MC
HAVELHC18j51=$tempBool
LHC18j5_1MC=$tempDir
OUTPUTDIR_LHC18j51=$tempPath
echo "18j5_1 anchored to 13bf: $HAVELHC18j51 $LHC18j5_1MC $OUTPUTDIR_LHC18j51"
FindCorrectTrainDirectory $LHC18j5_2MC $OUTPUTDIRMC $ALIENDIRMC $LHC18j5MC
HAVELHC18j52=$tempBool
LHC18j5_2MC=$tempDir
OUTPUTDIR_LHC18j52=$tempPath
echo "18j5_2 anchored to 13bf: $HAVELHC18j52 $LHC18j5_2MC $OUTPUTDIR_LHC18j52"
FindCorrectTrainDirectory $LHC18j5_3MC $OUTPUTDIRMC $ALIENDIRMC $LHC18j5MC
HAVELHC18j53=$tempBool
LHC18j5_3MC=$tempDir
OUTPUTDIR_LHC18j53=$tempPath
echo "18j5_3 anchored to 13bf: $HAVELHC18j53 $LHC18j5_3MC $OUTPUTDIR_LHC18j53"

# finding run2 MC path
FindCorrectTrainDirectory $LHC17f2bMC $OUTPUTDIRMC $ALIENDIRMC $LHC17f2bMCMoth
HAVELHC17f2b=$tempBool
LHC17f2bMC=$tempDir
OUTPUTDIR_LHC17f2b=$tempPath
echo "17f2b anchored to 16qt: $HAVELHC17f2b $LHC17f2bMC $OUTPUTDIR_LHC17f2b"
FindCorrectTrainDirectory $LHC17f2bMCFast $OUTPUTDIRMC $ALIENDIRMC $LHC17f2bMCMoth
HAVELHC17f2bF=$tempBool
LHC17f2bMCFast=$tempDir
OUTPUTDIR_LHC17f2bF=$tempPath
echo "17f2b_fast anchored to 16qt: $HAVELHC17f2b $LHC17f2bMCFast $OUTPUTDIR_LHC17f2bF"
FindCorrectTrainDirectory $LHC17f2a_fixMC $OUTPUTDIRMC $ALIENDIRMC $LHC17f2a_fixMCMoth
HAVELHC17f2afix=$tempBool
LHC17f2a_fixMC=$tempDir
OUTPUTDIR_LHC17f2a_fix=$tempPath
echo "17f2a_fix anchored to 16qt: $HAVELHC17f2afix $LHC17f2a_fixMC $OUTPUTDIR_LHC17f2a_fix"
FindCorrectTrainDirectory $LHC17f2a_fixMCFast $OUTPUTDIRMC $ALIENDIRMC $LHC17f2a_fixMCMoth
HAVELHC17f2afixF=$tempBool
LHC17f2a_fixMCFast=$tempDir
OUTPUTDIR_LHC17f2a_fixF=$tempPath
echo "17f2a_fix_fast anchored to 16qt: $HAVELHC17f2afixF $LHC17f2a_fixMCFast $OUTPUTDIR_LHC17f2a_fixF"
FindCorrectTrainDirectory $LHC17g8aMC $OUTPUTDIRMC $ALIENDIRMC $LHC17g8aMCMoth
HAVELHC17g8a=$tempBool
LHC17g8aMC=$tempDir
OUTPUTDIR_LHC17g8a=$tempPath
echo "17g8a anchored to 16qt: $HAVELHC17g8a $LHC17g8aMC $OUTPUTDIR_LHC17g8a"
FindCorrectTrainDirectory $LHC17g8aMCFast $OUTPUTDIRMC $ALIENDIRMC $LHC17g8aMCMoth
HAVELHC17g8aF=$tempBool
LHC17g8aMCFast=$tempDir
OUTPUTDIR_LHC17g8aF=$tempPath
echo "17g8a_fast anchored to 16qt: $HAVELHC17g8aF $LHC17g8aMCFast $OUTPUTDIR_LHC17g8aF"


currentDir=$PWD
if [ $CLEANUPMAYOR == 0 ]; then
        #echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackIncAccAndEMC" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13b" $HAVELHC13b $OUTPUTDIR_LHC13b $LHC13bData $pathDataR1 $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13c" $HAVELHC13c $OUTPUTDIR_LHC13c $LHC13cData $pathDataR1 $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13d" $HAVELHC13d $OUTPUTDIR_LHC13d $LHC13dData $pathDataR1 $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13e" $HAVELHC13e $OUTPUTDIR_LHC13e $LHC13eData $pathDataR1 $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13f" $HAVELHC13f $OUTPUTDIR_LHC13f $LHC13fData $pathDataR1 $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo

    cd $currentDir
    echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j5_1" $HAVELHC18j51 $OUTPUTDIR_LHC18j51 $LHC18j5_1MC $pathMCR1 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j5_2" $HAVELHC18j52 $OUTPUTDIR_LHC18j52 $LHC18j5_2MC $pathMCR1 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j5_3" $HAVELHC18j53 $OUTPUTDIR_LHC18j53 $LHC18j5_3MC $pathMCR1 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    CopyRunwiseAndMergeAccordingToRunlistData "LHC16q" $HAVELHC16q $OUTPUTDIR_LHC16q $LHC16qData $pathDataR2 $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    # "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData/"
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16qF" $HAVELHC16qF $OUTPUTDIR_LHC16qF $LHC16qDataFast $pathDataR2 $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$FAST" GammaConvCalo
    #   "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16qDataFast" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16qDataFast/"
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16t" $HAVELHC16t $OUTPUTDIR_LHC16t $LHC16tData $pathDataR2 $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    # "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData/"
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16tF" $HAVELHC16tF $OUTPUTDIR_LHC16tF $LHC16tDataFast $pathDataR2 $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$FAST" GammaConvCalo
    #  "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16tDataFast" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16tDataFast/"


    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f2b" $HAVELHC17f2b $OUTPUTDIR_LHC17f2b $LHC17f2bMC $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17f2b$5/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMC" $NSlashes3 "/alice/sim/2017/LHC17f2b$5/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMC/"
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f2b_fast" $HAVELHC17f2bF $OUTPUTDIR_LHC17f2bF $LHC17f2bMCFast $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17f2b_fast/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMCFast" $NSlashes3 "/alice/sim/2017/LHC17f2b_fast/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMCFast/"
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f2a_fix" $HAVELHC17f2afix $OUTPUTDIR_LHC17f2a_fix $LHC17f2a_fixMC $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17f2a$5_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC" $NSlashes3 "/alice/sim/2017/LHC17f2a$5_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/"
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f2a_fix_fast" $HAVELHC17f2afixF $OUTPUTDIR_LHC17f2a_fixF $LHC17f2a_fixMCFast $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17f2a_fast_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast" $NSlashes3 "/alice/sim/2017/LHC17f2a_fast_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast/"
    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC17g8a" $HAVELHC17g8a $OUTPUTDIR_LHC17g8a $LHC17g8aMC $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17g8a$5/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMC" $NSlashes3 "/alice/sim/2017/LHC17g8a$5/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMC/"
    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC17g8a_fast" $HAVELHC17g8aF $OUTPUTDIR_LHC17g8aF $LHC17g8aMCFast $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17g8a_fast/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMCFast" $NSlashes3 "/alice/sim/2017/LHC17g8a_fast/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMCFast/"

    listsToMerge=`cat runlistsToMerge.txt`
    for runListName in $listsToMerge; do
        if [ $HAVELHC13b == 1 ]; then
            ls $OUTPUTDIR_LHC13b/GammaConvCalo-$runListName\_*.root > fileLHC13b.txt
            fileNumbers=`cat fileLHC13b.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13c == 1 ]; then
            ls $OUTPUTDIR_LHC13c/GammaConvCalo-$runListName\_*.root > fileLHC13c.txt
            fileNumbers=`cat fileLHC13c.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13d == 1 ]; then
            ls $OUTPUTDIR_LHC13d/GammaConvCalo-$runListName\_*.root > fileLHC13d.txt
            fileNumbers=`cat fileLHC13d.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13d $NSlashes "LHC13d-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13e == 1 ]; then
            ls $OUTPUTDIR_LHC13e/GammaConvCalo-$runListName\_*.root > fileLHC13e.txt
            fileNumbers=`cat fileLHC13e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13e $NSlashes "LHC13e-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13f == 1 ]; then
            ls $OUTPUTDIR_LHC13f/GammaConvCalo-$runListName\_*.root > fileLHC13f.txt
            fileNumbers=`cat fileLHC13f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13f $NSlashes "LHC13f-pass4-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC18j51 == 1 ]; then
            ls $OUTPUTDIR_LHC18j51/GammaConvCalo-$runListName\_*.root > fileLHC18j5.txt
            fileNumbers=`cat fileLHC18j5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j51 $NSlashes "MC_LHC18j5_1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j52 == 1 ]; then
            ls $OUTPUTDIR_LHC18j52/GammaConvCalo-$runListName\_*.root > fileLHC18j5.txt
            fileNumbers=`cat fileLHC18j5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j52 $NSlashes "MC_LHC18j5_2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j53 == 1 ]; then
            ls $OUTPUTDIR_LHC18j53/GammaConvCalo-$runListName\_*.root > fileLHC18j5.txt
            fileNumbers=`cat fileLHC18j5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j53 $NSlashes "MC_LHC18j5_3-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC16q == 1 ]; then
            ls $OUTPUTDIR_LHC16q/GammaConvCalo-$runListName\_*.root > fileLHC16q.txt
            fileNumbers=`cat fileLHC16q.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16q $NSlashes "LHC16q-pass1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16qF == 1 ]; then
            ls $OUTPUTDIR_LHC16qF/GammaConvCalo-$runListName\_*.root > fileLHC16qF.txt
            fileNumbers=`cat fileLHC16qF.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16qF $NSlashes "LHC16q_fast-pass1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16t == 1 ]; then
            ls $OUTPUTDIR_LHC16t/GammaConvCalo-$runListName\_*.root > fileLHC16t.txt
            fileNumbers=`cat fileLHC16t.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16t $NSlashes "LHC16t-pass1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16tF == 1 ]; then
            ls $OUTPUTDIR_LHC16tF/GammaConvCalo-$runListName\_*.root > fileLHC16tF.txt
            fileNumbers=`cat fileLHC16tF.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16tF $NSlashes "LHC16t_fast-pass1-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC17f2b == 1 ]; then
            ls $OUTPUTDIR_LHC17f2b/GammaConvCalo-$runListName\_*.root > fileLHC17f2b.txt
            fileNumbers=`cat fileLHC17f2b.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f2b $NSlashes "MC_LHC17f2b-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17f2bF == 1 ]; then
            ls $OUTPUTDIR_LHC17f2bF/GammaConvCalo-$runListName\_*.root > fileLHC17f2bF.txt
            fileNumbers=`cat fileLHC17f2bF.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f2bF $NSlashes "MC_LHC17f2b_fast-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17f2afix == 1 ]; then
            ls $OUTPUTDIR_LHC17f2a_fix/GammaConvCalo-$runListName\_*.root > fileLHC17f2a_fix.txt
            fileNumbers=`cat fileLHC17f2a_fix.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f2a_fix $NSlashes "MC_LHC17f2a_fix-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17f2afixF == 1 ]; then
            ls $OUTPUTDIR_LHC17f2a_fixF/GammaConvCalo-$runListName\_*.root > fileLHC17f2a_fixF.txt
            fileNumbers=`cat fileLHC17f2a_fixF.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f2a_fixF $NSlashes "MC_LHC17f2a_fix_fast-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17g8a == 1 ]; then
            ls $OUTPUTDIR_LHC17g8a/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
            fileNumbers=`cat fileLHC17g8a.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17g8a $NSlashes "MC_LHC17g8a-$runListName" "-$runListName"
            done;
#             for binNumber in $binNumbersJJ; do
#                 echo $binNumber
#                 ls $OUTPUTDIR_LHC17g8a/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
#                 fileNumbers=`cat fileLHC17g8a.txt`
#                 for fileName in $fileNumbers; do
#                     echo $fileName
#                     GetFileNumberMerging $fileName $((NSlashes)) 2
#                     cp $OUTPUTDIR_LHC17g8a/$binNumber/GammaConvCalo-$runListName\_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber\_$3-$runListName\_$number.root
#                 done
#             done;
        fi
        if [ $HAVELHC17g8aF == 1 ]; then
            ls $OUTPUTDIR_LHC17g8aF/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
            fileNumbers=`cat fileLHC17g8a.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17g8aF $NSlashes "MC_LHC17g8a_fast-$runListName" "-$runListName"
            done;
#             for binNumber in $binNumbersJJ; do
#                 echo $binNumber
#                 ls $OUTPUTDIR_LHC17g8aF/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
#                 fileNumbers=`cat fileLHC17g8a.txt`
#                 for fileName in $fileNumbers; do
#                     echo $fileName
#                     GetFileNumberMerging $fileName $((NSlashes)) 2
#                     cp $OUTPUTDIR_LHC17g8aF/$binNumber/GammaConvCalo-$runListName\_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a_fast-$binNumber\_$3-$runListName\_$number.root
#                 done
#             done;
        fi
    done


    if [ $MERGEON == 1 ]; then
        echo "Starting Merging"

        echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
        listsToMerge=`cat runlistsToMerge.txt`
        for runListName in $listsToMerge; do
            ls $OUTPUTDIR/GammaConvCalo_LHC13f-pass4-$runListName\_*.root > filesForMerging.txt
            filesForMerging=`cat filesForMerging.txt`
            #13er daten
            periodList=`echo -e "b-pass4\nc-pass4\nd-pass4\ne-pass4\nf-pass4"`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                echo $number
                nameOut=""
                rm listCurrMerge.txt
                echo $fileName
                for periodID in $periodList; do
                    echo $periodID
                    currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-$runListName\_$number.root
                    if [ -f $currFile ]; then
                        outAdd=`echo $periodID  | cut -d "-" -f 1 `
                        nameOut+=$outAdd
                        echo -e "$currFile\n" >> listCurrMerge.txt
                    else
                        echo $currFile " does not exist"
                    fi
                done
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass4-$runListName\_$number.root
            done

            ls $OUTPUTDIR/GammaConvCalo_LHC13e-pass4-$runListName\_*.root > filesForMerging.txt
            periodList=`echo -e "b-pass4\nc-pass4\nd-pass4\ne-pass4"`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                echo $number
                nameOut=""
                rm listCurrMerge.txt
                echo $fileName
                for periodID in $periodList; do
                    echo $periodID
                    currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-$runListName\_$number.root
                    if [ -f $currFile ]; then
                        outAdd=`echo $periodID  | cut -d "-" -f 1 `
                        nameOut+=$outAdd
                        echo -e "$currFile\n" >> listCurrMerge.txt
                    else
                        echo $currFile " does not exist"
                    fi
                done
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass4-$runListName\_$number.root
            done

            ls $OUTPUTDIR/GammaConvCalo_MC_LHC18j5_1-$runListName\_*.root > filesForMerging.txt
            periodList=`echo -e "j5_1\nj5_2\nj5_3"`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                echo $number
                nameOut=""
                rm listCurrMerge.txt
                echo $fileName
                for periodID in $periodList; do
                    echo $periodID
                    currFile=$OUTPUTDIR/GammaConvCalo_MC_LHC18$periodID-$runListName\_$number.root
                    if [ -f $currFile ]; then
                        outAdd=`echo $periodID  | cut -d "-" -f 1 `
                        nameOut+=$outAdd
                        echo -e "$currFile\n" >> listCurrMerge.txt
                    else
                        echo $currFile " does not exist"
                    fi
                done
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC18j5x-$runListName\_$number.root
            done
        done

        echo -e "$3\nfast" > listReconstruction.txt
        listReconstruction=`cat listReconstruction.txt`
        for reco in $listReconstruction; do
            ls $OUTPUTDIR/GammaConvCalo_LHC16q_$reco-pass$passNr-DPGTrack\_*.root > filesForMerging.txt
            echo -e "DPGTrack" > runlistsToMerge.txt
            filesForMerging=`cat filesForMerging.txt`
            listsToMerge=`cat runlistsToMerge.txt`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 4
                echo $number
                for runListName in $listsToMerge; do
                    rm listCurrMerge.txt
                    fileQ="$OUTPUTDIR/GammaConvCalo_LHC16q_$reco-pass$passNr-$runListName""_$number.root"
                    fileT="$OUTPUTDIR/GammaConvCalo_LHC16t_$reco-pass$passNr-$runListName""_$number.root"
                    echo -e "$fileQ\n$fileT" > listCurrMerge.txt
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC16qt_$reco-pass$passNr-$runListName\_$number.root
                done
            done
        done
    fi

    if [ $MERGEONFASTAndWOSDD == 1 ]; then
        ls $OUTPUTDIR/GammaConvCalo_LHC16qt_fast-pass$passNr-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        echo -e "DPGTrack" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        listsToMerge=`cat runlistsToMerge.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_LHC16qt_fast-pass$passNr-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_LHC16qt_woSDD-pass$passNr-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC16qt_fast-woSDD-pass$passNr-$runListName\_$number.root
            done
        done

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 6
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-woSDD-$runListName\_$number.root
            done
        done

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC17f2b_fast-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2b_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2b_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17f2b_fast-woSDD-$runListName\_$number.root
            done
        done

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC17f2b_fast-woSDD-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2b_fast-woSDD-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_LHC17f2b_fast-woSDD-$runListName\_$number.root
            done
        done


        ls $OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_fast-DPGTrack\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_fast-woSDD-$runListName\_$number.root
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    rm listCurrMerge.txt
                    fileF="$OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber""_fast-$runListName""_$number.root"
                    fileW="$OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber""_woSDD-$runListName""_$number.root"
                    echo -e "$fileF\n$fileW" > listCurrMerge.txt
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber\_fast-woSDD-$runListName\_$number.root
                done
            done
        done
    fi
else
    if [ $HAVELHC16q == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16q";
#         rm $OUTPUTDIR_LHC16q/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16q/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16t";
#         rm $OUTPUTDIR_LHC16t/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16t/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC16qF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16q";
#         rm $OUTPUTDIR_LHC16q/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16qF/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC16tF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16t";
#         rm $OUTPUTDIR_LHC16t/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16tF/*/*/*GammaConvCalo_*.root
    fi

    if [ $HAVELHC17f2b == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17f2b";
#         rm $OUTPUTDIR_LHC17f2b/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17f2b/*/Stage*
    fi
    if [ $HAVELHC17f2afix == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17f2a_fix";
#         rm $OUTPUTDIR_LHC17f2a_fix/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17f2a_fix/*/Stage*
    fi
    if [ $HAVELHC17f2bF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17f2b";
#         rm $OUTPUTDIR_LHC17f2b/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17f2bF/*/Stage*
    fi
    if [ $HAVELHC17f2afixF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17f2a_fix";
#         rm $OUTPUTDIR_LHC17f2a_fix/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17f2a_fixF/*/Stage*
    fi
    if [ $HAVELHC17g8a == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17g8a";
#         rm $OUTPUTDIR_LHC17g8a/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17g8a/*/*/Stage*
    fi
    if [ $HAVELHC17g8aF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17g8a";
#         rm $OUTPUTDIR_LHC17g8a/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17g8aF/*/*/Stage*
    fi
fi
