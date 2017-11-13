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
HAVELHC17j7=1
HAVELHC17f2afix=1

# default trainconfigurations
LHC17nData="";
LHC17j7MC="";

passNr="1";
NSlashes=10;

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/XeXe
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
    NSlashes4=10
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
    NSlashes=9;
fi

# TRAINDIR=Legotrain-vAN20171029-XeXeQA
# LHC17nData="340"; #pass 3
# LHC17j7MC="1090";

# TRAINDIR=Legotrain-vAN20171105-XeXeQA
# LHC17nData="341"; #pass 1
# LHC17j7MC="664"
# LHC17j7MC="665"

TRAINDIR=Legotrain-vAN20171109-XeXeQA
# LHC17nData="345"; #pass 1
# LHC17nData="344"; #pass 1
# LHC17j7MC="673";
LHC17j7MC="672";
# LHC17j7MC="674";


OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ "$LHC17nData" == "" ]; then
    HAVELHC17n=0;
fi

if [ "$LHC17j7MC" == "" ]; then
    HAVELHC17j7=0;
    echo $LHC17j7MC
fi

mkdir -p $OUTPUTDIR/CutSelections


if [ $HAVELHC17n == 1 ]; then
    LHC17nData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/ | grep $LHC17nData\_`
    if [ "$LHC17nData" == "" ]; then
        HAVELHC17n=0;
    else
        OUTPUTDIR_LHC17n=$BASEDIR/$TRAINDIR/GA_PbPb-$LHC17nData
    fi
    echo $OUTPUTDIR_LHC17n
fi

if [ $HAVELHC17j7 == 1 ]; then
    LHC17j7MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/ | grep $LHC17j7MC\_`
    echo $LHC17j7MC
    if [ "$LHC17j7MC" == "" ]; then
        HAVELHC17j7=0;
    else
        OUTPUTDIR_LHC17j7=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC17j7MC
        echo $OUTPUTDIR_LHC17j7
    fi
fi

if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC17n == 1 ]; then
        echo "downloading LHC17n"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17n_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17n/$runNumber "/alice/data/2017/LHC17n/000$runNumber/pass$passNr/PWGGA/GA_PbPb/$LHC17nData" $NSlashes3 "/alice/data/2017/LHC17n/000$runNumber/pass$passNr/PWGGA/GA_PbPb/$LHC17nData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC17n/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC17n/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17n_all.txt`
                ls $OUTPUTDIR_LHC17n/$firstrunNumber/GammaConvV1_*.root > fileLHC17n.txt
                MergeAccordingToSpecificRunlist fileLHC17n.txt $OUTPUTDIR_LHC17n $NSlashes3 GammaConvV1 All runlists/runNumbersLHC17n_all.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17n "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/$LHC17nData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi

    currentDir=$PWD
    if [ $HAVELHC17j7 == 1 ]; then
        echo "downloading LHC17j7"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17j7_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17j7/$runNumber "/alice/sim/2017/LHC17j7/$runNumber/PWGGA/GA_PbPb_MC/$LHC17j7MC" $NSlashes3 "/alice/sim/2017/LHC17j7/$runNumber/PWGGA/GA_PbPb_MC/$LHC17j7MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17j7/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17j7/GammaConvV1*.root*
                echo runlists/runNumbersLHC17j7_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17j7_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17j7/$firstrunNumber/GammaConvV1_*.root > fileLHC17j7.txt
                MergeAccordingToSpecificRunlist fileLHC17j7.txt $OUTPUTDIR_LHC17j7 $NSlashes3 GammaConvV1 All runlists/runNumbersLHC17j7_all.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17j7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC17j7MC/merge" $NSlashes "" kTRUE
        fi
    fi

    if [ $HAVELHC17n == 1 ]; then
        ls $OUTPUTDIR_LHC17n/GammaConvV1-All_*.root > fileLHC17n.txt
        fileNumbers=`cat fileLHC17n.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17n $NSlashes "LHC17n-pass$passNr-All" "-All"
        done;
    fi

    if [ $HAVELHC17j7 == 1 ]; then
        ls $OUTPUTDIR_LHC17j7/GammaConvV1-All_*.root > fileLHC17j7.txt
        fileNumbers=`cat fileLHC17j7.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17j7 $NSlashes "MC_LHC17j7-All" "-All"
        done;
    fi

else
    if [ $HAVELHC17n == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17n";
        rm $OUTPUTDIR_LHC17n/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC17n/*/*/*GammaConvV1_*.root
    fi
    if [ $HAVELHC17j7 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17j7";
        rm $OUTPUTDIR_LHC17j7/*/GammaConvV1_*.root
    fi
fi
