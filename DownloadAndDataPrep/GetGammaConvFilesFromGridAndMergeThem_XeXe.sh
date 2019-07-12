#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=1
SINGLERUN=1
SEPARATEON=0
MERGEONSINGLEData=1
MERGEONSINGLEMC=1
CLEANUP=1
CLEANUPMAYOR=$2
number=""

# check if train configuration has actually been given
HAVELHC17n=1
HAVETOBUILDLHC18d2=0
HAVELHC18d2a=1
HAVELHC18d2b=1
HAVELHC18d2c=1
HAVELHC18d2d=1

# default trainconfigurations
LHC17nData="";
LHC18d2MC="";
LHC18d2MCa="";
LHC18d2MCb="";
LHC18d2MCc="";
LHC18d2MCd="";
passNr="1";

if [ $1 = "fbock" ]; then
    BASEDIR=/media/fbock/Elements/OutputLegoTrains/XeXe
elif [ $1 = "fbockExt" ]; then
    BASEDIR=/media/fbock/Elements/OutputLegoTrains/XeXe
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
fi

# Definitition of number of slashes in your path to different depths
NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
NSlashes=`expr $NSlashesBASE + 4`
NSlashes2=`expr $NSlashes - 1`
NSlashes3=`expr $NSlashes + 1`
NSlashes4=`expr $NSlashes + 2`
echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"

# TRAINDIR=Legotrain-vAN20180531-EMCAndPHOSAndPCM
# ISAOD=0
# LHC17nData="406"; #pass 1
# LHC18d2MC="943"; # PCM
# LHC18d2MCa="child_1";
# LHC18d2MCb="child_2";
# LHC18d2MCc="child_3";
# LHC18d2MCd="child_4";

# TRAINDIR=Legotrain-vAN20180612-EMCAndPHOSAndPCM
# ISAOD=0
# LHC17nData="412"; #pass 1
# LHC17nData="407"; #pass 1
# LHC18d2MC="950"; # PHOS nl
# LHC18d2MC="944"; # PHOS nl
# LHC18d2MCa="child_1";
# LHC18d2MCb="child_2";
# LHC18d2MCc="child_3";
# LHC18d2MCd="child_4";

# TRAINDIR=20190520-PCMtrees-QA
# ISAOD=1
# # LHC17nData="510";
# LHC18d2MC="1317";
# LHC18d2MCa="child_1";
# LHC18d2MCb="child_2";
# LHC18d2MCc="child_3";
# LHC18d2MCd="child_4";

# TRAINDIR=20190618-PCMtrees-QA
# ISAOD=1
# LHC17nData="518";
# LHC18d2MC="1339";
# LHC18d2MCa="child_1";
# LHC18d2MCb="child_2";
# LHC18d2MCc="child_3";
# LHC18d2MCd="child_4";

TRAINDIR=20190711-PCMtrees-QA
ISAOD=1
LHC17nData="522";

OUTPUTDIR=$BASEDIR/$TRAINDIR
mkdir -p $OUTPUTDIR/CutSelections
if [ $ISAOD -eq 1 ]; then
    baseLegoData=GA_PbPb_AOD
    baseLegoMC=GA_PbPb_MC_AOD
    pathData=pass1/PWGGA/GA_PbPb_AOD
    pathMC=PWGGA/GA_PbPb_MC_AOD
else
    baseLegoData=GA_PbPb
    baseLegoMC=GA_PbPb_MC
    pathData=pass1/PWGGA/GA_PbPb_AOD
    pathMC=PWGGA/GA_PbPb_MC_AOD
fi
mkdir -p $OUTPUTDIR/CutSelections
ALIENDIRData="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoData/"
OUTPUTDIRData=$BASEDIR/$TRAINDIR/$baseLegoData
ALIENDIRMC="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoMC/"
OUTPUTDIRMC=$BASEDIR/$TRAINDIR/$baseLegoMC

# finding run1 data paths
FindCorrectTrainDirectory $LHC17nData $OUTPUTDIRData $ALIENDIRData
HAVELHC17n=$tempBool
LHC17nData=$tempDir
OUTPUTDIR_LHC17n=$tempPath
echo "17n: $HAVELHC17n $LHC17nData $OUTPUTDIR_LHC17n"
FindCorrectTrainDirectory $LHC18d2MCa $OUTPUTDIRMC $ALIENDIRMC $LHC18d2MC
HAVELHC18d2a=$tempBool
LHC18d2MCa=$tempDir
OUTPUTDIR_LHC18d2a=$tempPath
echo "18d2a: $HAVELHC18d2a $LHC18d2MCa $OUTPUTDIR_LHC18d2a"
FindCorrectTrainDirectory $LHC18d2MCb $OUTPUTDIRMC $ALIENDIRMC $LHC18d2MC
HAVELHC18d2b=$tempBool
LHC18d2MCb=$tempDir
OUTPUTDIR_LHC18d2b=$tempPath
echo "18d2b: $HAVELHC18d2b $LHC18d2MCb $OUTPUTDIR_LHC18d2b"
FindCorrectTrainDirectory $LHC18d2MCc $OUTPUTDIRMC $ALIENDIRMC $LHC18d2MC
HAVELHC18d2c=$tempBool
LHC18d2MCc=$tempDir
OUTPUTDIR_LHC18d2c=$tempPath
echo "18d2c: $HAVELHC18d2c $LHC18d2MCc $OUTPUTDIR_LHC18d2c"
FindCorrectTrainDirectory $LHC18d2MCd $OUTPUTDIRMC $ALIENDIRMC $LHC18d2MC
HAVELHC18d2d=$tempBool
LHC18d2MCd=$tempDir
OUTPUTDIR_LHC18d2d=$tempPath
echo "18d2a: $HAVELHC18d2d $LHC18d2MCd $OUTPUTDIR_LHC18d2d"



if [ $CLEANUPMAYOR == 0 ]; then

    currentDir=$PWD
    echo -e "all" > runlistsToMerge.txt
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17n" $HAVELHC17n $OUTPUTDIR_LHC17n $LHC17nData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvV1
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d2_1" $HAVELHC18d2a $OUTPUTDIR_LHC18d2a $LHC18d2MCa $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvV1
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d2_2" $HAVELHC18d2b $OUTPUTDIR_LHC18d2b $LHC18d2MCb $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvV1
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d2_3" $HAVELHC18d2c $OUTPUTDIR_LHC18d2c $LHC18d2MCc $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvV1
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d2_4" $HAVELHC18d2d $OUTPUTDIR_LHC18d2d $LHC18d2MCd $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvV1

    if [ $HAVELHC17n == 1 ]; then
        ls $OUTPUTDIR_LHC17n/GammaConvV1-all_*.root > fileLHC17n.txt
        fileNumbers=`cat fileLHC17n.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17n $NSlashes "LHC17n-pass$passNr-all" "-all"
        done;
    fi

    if [ $HAVELHC18d2a == 1 ]; then
        ls $OUTPUTDIR_LHC18d2a/GammaConvV1-all_*.root > fileLHC18d2a.txt
        fileNumbers=`cat fileLHC18d2a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2a $NSlashes "MC_LHC18d2-4-all" "-all"
        done;
    fi
    if [ $HAVELHC18d2b == 1 ]; then
        ls $OUTPUTDIR_LHC18d2b/GammaConvV1-all_*.root > fileLHC18d2b.txt
        fileNumbers=`cat fileLHC18d2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2b $NSlashes "MC_LHC18d2-1-all" "-all"
        done;
    fi
    if [ $HAVELHC18d2c == 1 ]; then
        ls $OUTPUTDIR_LHC18d2c/GammaConvV1-all_*.root > fileLHC18d2c.txt
        fileNumbers=`cat fileLHC18d2c.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2c $NSlashes "MC_LHC18d2-2-all" "-all"
        done;
    fi
    if [ $HAVELHC18d2d == 1 ]; then
        ls $OUTPUTDIR_LHC18d2d/GammaConvV1-all_*.root > fileLHC18d2d.txt
        fileNumbers=`cat fileLHC18d2d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2d $NSlashes "MC_LHC18d2-3-all" "-all"
        done;
    fi

    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaConvV1_MC_LHC18d2-4-all_*.root > filesForMerging.txt
        echo -e "\nall" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            rm listCurrMerge.txt
            fileA="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-4-all_$number.root"
            fileB="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-1-all_$number.root"
            fileC="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-2-all_$number.root"
            fileD="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-3-all_$number.root"
            echo -e "$fileA\n$fileB\n$fileC\n$fileD" > listCurrMerge.txt
            MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvV1_MC_LHC18d2x-all_$number.root
        done
    fi


else
    if [ $HAVELHC17n == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC17n";
        rm $OUTPUTDIR_LHC17n/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC17n/*/*/*GammaConvV1_*.root
    fi
    if [ $HAVELHC18d2a == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC18d2_4";
        rm $OUTPUTDIR_LHC18d2a/*/GammaConvV1_*.root
    fi
    if [ $HAVELHC18d2b == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC18d2_1";
        rm $OUTPUTDIR_LHC18d2b/*/GammaConvV1_*.root
    fi
    if [ $HAVELHC18d2c == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC18d2_2";
        rm $OUTPUTDIR_LHC18d2c/*/GammaConvV1_*.root
    fi
    if [ $HAVELHC18d2d == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC18d2_3";
        rm $OUTPUTDIR_LHC18d2d/*/GammaConvV1_*.root
    fi

fi
