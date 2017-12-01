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
HAVELHC16k=1
HAVELHC16l=1
HAVETOBUILDData=0
HAVELHC17d20a1=1
HAVELHC17d20a1Ex=1
HAVELHC17d20a2=1
HAVELHC17d20a2Ex=1
HAVETOBUILDMC=0

# default trainconfigurations
LHC16kData="";
LHC16lData="";
LHC16Data="";
LHC17MCPythia="";
LHC17MCEPOS="";
LHC17d20a1MC="";
LHC17d20a1ExMC="";
LHC17d20a2MC="";
LHC17d20a2ExMC="";

passNr="1";
NSlashes=10;

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorageExternal/OutputLegoTrains/pp13TeV
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
    NSlashes4=10
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
    NSlashes=9;
fi

TRAINDIR=Legotrain-vAN20171105-QA
LHC16Data="2254"; #pass 1
#LHC16kData="child_6"; #pass 1
LHC16lData="child_7"; #pass 1
LHC17MCPythia="3176"; #pass 1
#LHC17d20a1MC="child_7";
#LHC17d20a1ExMC="child_8";
LHC17d20a2MC="child_9";
#LHC17d20a2ExMC="child_10";

OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ "$LHC16kData" == "" ]; then
    HAVELHC16k=0;
fi
if [ "$LHC16lData" = "" ]; then
    HAVELHC16l=0;
fi
if [ "$LHC16Data" != "" ]; then
    HAVETOBUILDData=1;
fi

if [ "$LHC17d20a1MC" == "" ]; then
    HAVELHC17d20a1=0;
    echo $LHC17d20a1MC
fi
if [ "$LHC17d20a1ExMC" == "" ]; then
    HAVELHC17d20a1Ex=0;
fi
if [ "$LHC17d20a2MC" == "" ]; then
    HAVELHC17d20a2=0;
    echo $LHC17d20a2MC
fi
if [ "$LHC17d20a2ExMC" == "" ]; then
    HAVELHC17d20a2Ex=0;
fi
if [ "$LHC17MCPythia" != "" ]; then
    HAVETOBUILDMC=1;
fi
if [ "$LHC17MCEPOS" != "" ]; then
    HAVETOBUILDMC=1;
fi


mkdir -p $OUTPUTDIR/CutSelections

# Get data directory for 16k period
if [ $HAVELHC16k == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16kData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16kData`
    else
        LHC16kData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16kData\_`
    fi
    if [ "$LHC16kData" == "" ]; then
        HAVELHC16k=0;
    else
        OUTPUTDIR_LHC16k=$BASEDIR/$TRAINDIR/GA_pp-$LHC16kData
    fi
    echo $OUTPUTDIR_LHC16k
fi
# Get data directory for 16l period
if [ $HAVELHC16l == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16lData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16lData`
    else
        LHC16lData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16lData\_`
    fi
    if [ "$LHC16lData" == "" ]; then
        HAVELHC16l=0;
    else
        OUTPUTDIR_LHC16l=$BASEDIR/$TRAINDIR/GA_pp-$LHC16lData
    fi
    echo $OUTPUTDIR_LHC16l
fi

# Get MC directory for LHC17d20a1 MC anchored to LHC16k
if [ $HAVELHC17d20a1 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a1MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a1MC`
    else
        LHC17d20a1MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a1MC\_`
    fi
    if [ "$LHC17d20a1MC" == "" ]; then
        HAVELHC17d20a1=0;
    else
        OUTPUTDIR_LHC17d20a1=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a1MC
    fi
    echo $OUTPUTDIR_LHC17d20a1
fi

# Get MC directory for LHC17d20a1_extra MC anchored to LHC16k
if [ $HAVELHC17d20a1Ex == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a1ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a1ExMC`
    else
        LHC17d20a1ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a1ExMC\_`
    fi
    if [ "$LHC17d20a1ExMC" == "" ]; then
        HAVELHC17d20a1Ex=0;
    else
        OUTPUTDIR_LHC17d20a1Ex=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a1ExMC
    fi
    echo $OUTPUTDIR_LHC17d20a1Ex
fi

# Get MC directory for LHC17d20a2 MC anchored to LHC16k
if [ $HAVELHC17d20a2 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a2MC`
    else
        LHC17d20a2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a2MC\_`
    fi
    if [ "$LHC17d20a2MC" == "" ]; then
        HAVELHC17d20a2=0;
    else
        OUTPUTDIR_LHC17d20a2=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a2MC
    fi
    echo $OUTPUTDIR_LHC17d20a2
fi

# Get MC directory for LHC17d20a2_extra MC anchored to LHC16k
if [ $HAVELHC17d20a2Ex == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a2ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a2ExMC`
    else
        LHC17d20a2ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a2ExMC\_`
    fi
    if [ "$LHC17d20a2ExMC" == "" ]; then
        HAVELHC17d20a2Ex=0;
    else
        OUTPUTDIR_LHC17d20a2Ex=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a2ExMC
    fi
    echo $OUTPUTDIR_LHC17d20a2Ex
fi



if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC16k == 1 ]; then
        echo "downloading LHC16k"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16k.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16k/$runNumber "/alice/data/2016/LHC16k/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16kData" $NSlashes3 "/alice/data/2016/LHC16k/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16kData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16k/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16k/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16k_pass1.txt`
                ls $OUTPUTDIR_LHC16k/$firstrunNumber/GammaConvCalo_*.root > fileLHC16k.txt

                MergeAccordingToSpecificRunlist fileLHC16k.txt $OUTPUTDIR_LHC16k $NSlashes3 GammaConvCalo All runlists/runNumbersLHC16k_pass1.txt
                MergeAccordingToSpecificRunlist fileLHC16k.txt $OUTPUTDIR_LHC16k $NSlashes3 GammaConvCalo DPGTrack runlists/runNumbersLHC16k_pass1_DPG.txt
                MergeAccordingToSpecificRunlist fileLHC16k.txt $OUTPUTDIR_LHC16k $NSlashes3 GammaConvCalo DPGTrackAndCalo runlists/runNumbersLHC16k_pass1_DPGTrackAndCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16k "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16kData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16l == 1 ]; then
        echo "downloading LHC16l"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16l_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16l/$runNumber "/alice/data/2016/LHC16l/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16lData" $NSlashes3 "/alice/data/2016/LHC16l/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16lData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16l/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16l/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16l_pass1.txt`
                ls $OUTPUTDIR_LHC16l/$firstrunNumber/GammaConvCalo_*.root > fileLHC16l.txt
                fileNumbers=`cat fileLHC16l.txt`
                MergeAccordingToSpecificRunlist fileLHC16l.txt $OUTPUTDIR_LHC16l $NSlashes3 GammaConvCalo All runlists/runNumbersLHC16l_pass1.txt
                MergeAccordingToSpecificRunlist fileLHC16l.txt $OUTPUTDIR_LHC16l $NSlashes3 GammaConvCalo DPGTrack runlists/runNumbersLHC16l_pass1_DPG.txt
                MergeAccordingToSpecificRunlist fileLHC16l.txt $OUTPUTDIR_LHC16l $NSlashes3 GammaConvCalo DPGTrackAndCalo runlists/runNumbersLHC16l_pass1_DPGTrackAndCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16l "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16lData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi

    currentDir=$PWD
    if [ $HAVELHC17d20a1 == 1 ]; then
        echo "downloading LHC17d20a1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1/$runNumber "/alice/sim/2017/LHC17d20a1/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1MC" $NSlashes3 "/alice/sim/2017/LHC17d20a1/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d20a1/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a1/GammaConvCalo*.root*
                echo runlists/runNumbersLHC17d20a1.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a1.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d20a1/$firstrunNumber/GammaConvCalo_*.root > fileLHC17d20a1.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a1.txt $OUTPUTDIR_LHC17d20a1 $NSlashes3 GammaConvCalo All runlists/runNumbersLHC17d20a1.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a1.txt $OUTPUTDIR_LHC17d20a1 $NSlashes3 GammaConvCalo DPGTrack runlists/runNumbersLHC17d20a1_DPG.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a1.txt $OUTPUTDIR_LHC17d20a1 $NSlashes3 GammaConvCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a1_DPGTrackAndCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a1MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d20a1Ex == 1 ]; then
        echo "downloading LHC17d20a1_extra"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a1_extra.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1Ex/$runNumber "/alice/sim/2017/LHC17d20a1_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1ExMC" $NSlashes3 "/alice/sim/2017/LHC17d20a1_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1ExMC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17d20a1Ex/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a1Ex/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a1_extra.txt`
                ls $OUTPUTDIR_LHC17d20a1Ex/$firstrunNumber/GammaConvCalo_*.root > fileLHC17d20a1_extra.txt
                fileNumbers=`cat fileLHC17d20a1_extra.txt`
                MergeAccordingToSpecificRunlist fileLHC17d20a1_extra.txt $OUTPUTDIR_LHC17d20a1Ex $NSlashes3 GammaConvCalo All runlists/runNumbersLHC17d20a1_extra.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a1_extra.txt $OUTPUTDIR_LHC17d20a1Ex $NSlashes3 GammaConvCalo DPGTrack runlists/runNumbersLHC17d20a1_extra_DPG.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a1_extra.txt $OUTPUTDIR_LHC17d20a1Ex $NSlashes3 GammaConvCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a1_extra_DPGTrackAndCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1Ex "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a1ExMC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d20a2 == 1 ]; then
        echo "downloading LHC17d20a2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a2.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2/$runNumber "/alice/sim/2017/LHC17d20a2/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2MC" $NSlashes3 "/alice/sim/2017/LHC17d20a2/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d20a2/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a2/GammaConvCalo*.root*
                echo runlists/runNumbersLHC17d20a2.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a2.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d20a2/$firstrunNumber/GammaConvCalo_*.root > fileLHC17d20a2.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a2.txt $OUTPUTDIR_LHC17d20a2 $NSlashes3 GammaConvCalo All runlists/runNumbersLHC17d20a2.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a2.txt $OUTPUTDIR_LHC17d20a2 $NSlashes3 GammaConvCalo DPGTrack runlists/runNumbersLHC17d20a2_DPG.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a2.txt $OUTPUTDIR_LHC17d20a2 $NSlashes3 GammaConvCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a2_DPG.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a2MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d20a2Ex == 1 ]; then
        echo "downloading LHC17d20a2_extra"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a2_extra.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2Ex/$runNumber "/alice/sim/2017/LHC17d20a2_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2ExMC" $NSlashes3 "/alice/sim/2017/LHC17d20a2_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2ExMC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17d20a2Ex/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a2Ex/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a2_extra.txt`
                ls $OUTPUTDIR_LHC17d20a2Ex/$firstrunNumber/GammaConvCalo_*.root > fileLHC17d20a2_extra.txt
                fileNumbers=`cat fileLHC17d20a2_extra.txt`
                MergeAccordingToSpecificRunlist fileLHC17d20a2_extra.txt $OUTPUTDIR_LHC17d20a2Ex $NSlashes3 GammaConvCalo All runlists/runNumbersLHC17d20a2_extra.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a2_extra.txt $OUTPUTDIR_LHC17d20a2Ex $NSlashes3 GammaConvCalo DPGTrack runlists/runNumbersLHC17d20a2_extra_DPG.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a2_extra.txt $OUTPUTDIR_LHC17d20a2Ex $NSlashes3 GammaConvCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a2_extra_DPGTrackAndCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2Ex "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a2ExMC/merge" $NSlashes "" kTRUE
        fi
    fi


    if [ $HAVELHC16k == 1 ]; then
        ls $OUTPUTDIR_LHC16k/GammaConvCalo-All_*.root > fileLHC16k.txt
        fileNumbers=`cat fileLHC16k.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16k/GammaConvCalo-DPGTrack_*.root > fileLHC16k.txt
        fileNumbers=`cat fileLHC16k.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16k/GammaConvCalo-DPGTrackAndCalo_*.root > fileLHC16k.txt
        fileNumbers=`cat fileLHC16k.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16l == 1 ]; then
        ls $OUTPUTDIR_LHC16l/GammaConvCalo-All_*.root > fileLHC16l.txt
        fileNumbers=`cat fileLHC16l.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16l/GammaConvCalo-DPGTrack_*.root > fileLHC16l.txt
        fileNumbers=`cat fileLHC16l.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16l/GammaConvCalo-DPGTrackAndCalo_*.root > fileLHC16l.txt
        fileNumbers=`cat fileLHC16l.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d20a1 == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a1/GammaConvCalo-All_*.root > fileLHC17d20a1.txt
        fileNumbers=`cat fileLHC17d20a1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a1/GammaConvCalo-DPGTrack_*.root > fileLHC17d20a1.txt
        fileNumbers=`cat fileLHC17d20a1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a1/GammaConvCalo-DPGTrackAndCalo_*.root > fileLHC17d20a1.txt
        fileNumbers=`cat fileLHC17d20a1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d20a1Ex == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a1Ex/GammaConvCalo-All_*.root > fileLHC17d20a1_extra.txt
        fileNumbers=`cat fileLHC17d20a1_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a1Ex/GammaConvCalo-DPGTrack_*.root > fileLHC17d20a1_extra.txt
        fileNumbers=`cat fileLHC17d20a1_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a1Ex/GammaConvCalo-DPGTrackAndCalo*.root > fileLHC17d20a1_extra.txt
        fileNumbers=`cat fileLHC17d20a1_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17d20a2 == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a2/GammaConvCalo-All_*.root > fileLHC17d20a2.txt
        fileNumbers=`cat fileLHC17d20a2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a2/GammaConvCalo-DPGTrack_*.root > fileLHC17d20a2.txt
        fileNumbers=`cat fileLHC17d20a2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a2/GammaConvCalo-DPGTrackAndCalo_*.root > fileLHC17d20a2.txt
        fileNumbers=`cat fileLHC17d20a2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d20a2Ex == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a2Ex/GammaConvCalo-All_*.root > fileLHC17d20a2_extra.txt
        fileNumbers=`cat fileLHC17d20a2_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2Ex $NSlashes "MC_LHC17d20a2_extra-anchor16l-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a2Ex/GammaConvCalo-DPGTrack_*.root > fileLHC17d20a2_extra.txt
        fileNumbers=`cat fileLHC17d20a2_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2Ex $NSlashes "MC_LHC17d20a2_extra-anchor16l-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a2Ex/GammaConvCalo-DPGTrackAndCalo*.root > fileLHC17d20a2_extra.txt
        fileNumbers=`cat fileLHC17d20a2_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2Ex $NSlashes "MC_LHC17d20a2_extra-anchor16l-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrackAndCalo\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root
            ls $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root
            if [ -f $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_LHC16kl-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root
            fi
        done
        ls $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrack\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrack\_$number.root
            ls $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-DPGTrack\_$number.root
            if [ -f $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrack\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-DPGTrack\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_LHC16kl-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-DPGTrack\_$number.root
            fi
        done
        ls $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-All\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-All\_$number.root
            ls $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-All\_$number.root
            if [ -f $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-All\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-All\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_LHC16kl-pass$passNr-All\_$number.root $OUTPUTDIR/GammaConvCalo_LHC16k-pass$passNr-All\_$number.root $OUTPUTDIR/GammaConvCalo_LHC16l-pass$passNr-All\_$number.root
            fi
        done



    fi
else
    if [ $HAVELHC16k == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16k";
        rm $OUTPUTDIR_LHC16k/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16k/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC16l == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16l";
        rm $OUTPUTDIR_LHC16l/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16l/*/*/*GammaConvCalo_*.root
    fi

    if [ $HAVELHC17d20a1 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a1";
        rm $OUTPUTDIR_LHC17d20a1/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC17d20a1/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC17d20a1Ex == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a1Ex";
        rm $OUTPUTDIR_LHC17d20a1Ex/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC17d20a1Ex/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC17d20a2 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a2";
        rm $OUTPUTDIR_LHC17d20a2/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC17d20a2/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC17d20a2Ex == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a2Ex";
        rm $OUTPUTDIR_LHC17d20a2Ex/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC17d20a2Ex/*/*/*GammaConvCalo_*.root
    fi
fi
