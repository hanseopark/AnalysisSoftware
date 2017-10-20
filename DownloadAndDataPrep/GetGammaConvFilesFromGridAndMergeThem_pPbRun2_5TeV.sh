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
MERGEONSINGLEMC=0
CLEANUP=1
CLEANUPMAYOR=$2
number=""

# check if train configuration has actually been given
HAVELHC16q=1
HAVELHC16t=1
HAVETOBUILDData=0
HAVELHC17f2b=1
HAVELHC17f2afix=1

# default trainconfigurations
LHC16qData="";
LHC16tData="";
LHC16qtData="";
LHC17f2bMC="";
LHC17f2a_fixMC="";

passNr="1";
NSlashes=10;

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pPb
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
    NSlashes4=10
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
    NSlashes=9;
fi

TRAINDIR=Legotrain-vAN20171005-dirGammaRun2
LHC16qtData="679"; #pass 2
LHC16qData="child_1"; #pass 3
LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1090";
# LHC17f2a_fixMC="1088";

OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ "$LHC16qData" == "" ]; then
    HAVELHC16q=0;
fi
if [ "$LHC16tData" = "" ]; then
    HAVELHC16t=0;
fi
if [ "$LHC16qtData" != "" ]; then
    HAVETOBUILDData=1;
fi

if [ "$LHC17f2bMC" == "" ]; then
    HAVELHC17f2b=0;
    echo $LHC17f2bMC
fi
if [ "$LHC17f2a_fixMC" == "" ]; then
    HAVELHC17f2afix=0;
fi
if [ "$LHC13e7MC" = "" ]; then
    HAVELHC13e7=0;
fi

mkdir -p $OUTPUTDIR/CutSelections


if [ $HAVELHC16q == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\_ | grep $LHC16qData`
    else
        LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qData\_`
    fi
    if [ "$LHC16qData" == "" ]; then
        HAVELHC16q=0;
    else
        OUTPUTDIR_LHC16q=$BASEDIR/$TRAINDIR/GA_pPb-$LHC16qData
    fi
    echo $OUTPUTDIR_LHC16q
fi
if [ $HAVELHC16t == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16tData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\_ | grep $LHC16tData`
    else
        LHC16tData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16tData\_`
    fi
    if [ "$LHC16tData" == "" ]; then
        HAVELHC16t=0;
    else
        OUTPUTDIR_LHC16t=$BASEDIR/$TRAINDIR/GA_pPb-$LHC16tData
    fi
    echo $OUTPUTDIR_LHC16t
fi

if [ $HAVELHC17f2b == 1 ]; then
    LHC17f2bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMC\_`
    echo $LHC17f2bMC
    if [ "$LHC17f2bMC" == "" ]; then
        HAVELHC17f2b=0;
    else
        OUTPUTDIR_LHC17f2b=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2bMC
        echo $OUTPUTDIR_LHC17f2b
    fi
fi
if [ $HAVELHC17f2afix == 1 ]; then
    LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMC\_`
    if [ "$LHC17f2a_fixMC" == "" ]; then
        HAVELHC17f2afix=0;
    else
        OUTPUTDIR_LHC17f2a_fix=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2a_fixMC
    fi
fi

if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC16q == 1 ]; then
        echo "downloading LHC16q"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16q_all_$3.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16q/$runNumber "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData/Stage_1/"
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16q/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16q/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16q_all_$3.txt`
                ls $OUTPUTDIR_LHC16q/$firstrunNumber/GammaConvV1_*.root > fileLHC16q.txt

                MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16q $NSlashes3 GammaConvV1 All runlists/runNumbersLHC16q_all_$3.txt
                MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16q $NSlashes3 GammaConvV1 DPGTrack runlists/runNumbersLHC16q_dpgTracks_$3.txt
                MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16q $NSlashes3 GammaConvV1 DPGTrackAndCalo runlists/runNumbersLHC16q_dpgTracksAndCalo_$3.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16q "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16qData/merge_runlist_1" $NSlashes
        fi
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "downloading LHC16t"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16t_all_$3.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16t/$runNumber "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData/Stage_1/"
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16t/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16t/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16t_all_$3.txt`
                ls $OUTPUTDIR_LHC16t/$firstrunNumber/GammaConvV1_*.root > fileLHC16t.txt
                fileNumbers=`cat fileLHC16t.txt`
                MergeAccordingToSpecificRunlist fileLHC16t.txt $OUTPUTDIR_LHC16t $NSlashes3 GammaConvV1 All runlists/runNumbersLHC16t_all_$3.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16t "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16tData/merge_runlist_1" $NSlashes
        fi
    fi

    currentDir=$PWD
    if [ $HAVELHC17f2b == 1 ]; then
        echo "downloading LHC17f2b"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f2b.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2b/$runNumber "/alice/sim/2017/LHC17f2b$5/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMC" $NSlashes3 "/alice/sim/2017/LHC17f2b$5/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMC/Stage_1/"
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17f2b/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f2b/GammaConvV1*.root*
                echo runlists/runNumbersLHC17f2b_all_$3.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2b_all_$3.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17f2b/$firstrunNumber/GammaConvV1_*.root > fileLHC17f2b.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConvV1 All runlists/runNumbersLHC16q_all_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConvV1 DPGTrack runlists/runNumbersLHC16q_dpgTracks_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConvV1 DPGTrackAndCalo runlists/runNumbersLHC16q_dpgTracksAndCalo_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConvV1 All runlists/runNumbersLHC16t_all_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConvV1 DPGTrack runlists/runNumbersLHC16t_all_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConvV1 DPGTrackAndCalo runlists/runNumbersLHC16t_all_$3.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMC/merge" $NSlashes
        fi
    fi
    if [ $HAVELHC17f2afix == 1 ]; then
        echo "downloading LHC17f2a_fix"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f2a_fix_all_$3.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2a_fix/$runNumber "/alice/sim/2017/LHC17f2a$5_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC" $NSlashes3 "/alice/sim/2017/LHC17f2a$5_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/Stage_1/"
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17f2a_fix/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f2a_fix/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2a_fix_all_$3.txt`
                ls $OUTPUTDIR_LHC17f2a_fix/$firstrunNumber/GammaConvV1_*.root > fileLHC17f2a_fix.txt
                fileNumbers=`cat fileLHC17f2a_fix.txt`
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConvV1 All runlists/runNumbersLHC16q_all_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConvV1 DPGTrack runlists/runNumbersLHC16q_dpgTracks_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConvV1 DPGTrackAndCalo runlists/runNumbersLHC16q_dpgTracksAndCalo_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConvV1 All runlists/runNumbersLHC16t_all_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConvV1 DPGTrack runlists/runNumbersLHC16t_all_$3.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConvV1 DPGTrackAndCalo runlists/runNumbersLHC16t_all_$3.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2a_fix "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/merge" $NSlashes
        fi
    fi

    if [ $HAVELHC16q == 1 ]; then
        ls $OUTPUTDIR_LHC16q/GammaConvV1_*.root > fileLHC16q.txt
        fileNumbers=`cat fileLHC16q.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16q $NSlashes "LHC16q-pass$passNr"
        done;
    fi

    if [ $HAVELHC16t == 1 ]; then
        ls $OUTPUTDIR_LHC16t/GammaConvV1_*.root > fileLHC16t.txt
        fileNumbers=`cat fileLHC16t.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16t $NSlashes "LHC16t-pass$passNr"
        done;
    fi

    if [ $HAVELHC17f2b == 1 ]; then
        ls $OUTPUTDIR_LHC17f2b/GammaConvV1_*.root > fileLHC17f2b.txt
        fileNumbers=`cat fileLHC17f2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f2b $NSlashes "MC_LHC17f2b"
        done;
    fi

    if [ $HAVELHC17f2afix == 1 ]; then
        ls $OUTPUTDIR_LHC17f2a_fix/GammaConvV1_*.root > fileLHC17f2a_fix.txt
        fileNumbers=`cat fileLHC17f2a_fix.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f2a_fix $NSlashes "MC_LHC17f2a_fix"
        done;
    fi

    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaConvV1_LHC16q-pass$passNr\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 3
            echo $number
            ls $OUTPUTDIR/GammaConvV1_LHC16q-pass$passNr\_$number.root
            ls $OUTPUTDIR/GammaConvV1_LHC16t-pass$passNr\_$number.root
            if [ -f $OUTPUTDIR/GammaConvV1_LHC16q-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16t-pass$passNr\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16qt-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC16q-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC16t-pass$passNr\_$number.root
            fi
        done
    fi
else
    if [ $HAVELHC16q == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16q";
        rm $OUTPUTDIR_LHC16q/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC16q/*/*/*GammaConvV1_*.root
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16t";
        rm $OUTPUTDIR_LHC16t/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC16t/*/*/*GammaConvV1_*.root
    fi

    if [ $HAVELHC17f2b == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17f2b";
        rm $OUTPUTDIR_LHC17f2b/*/GammaConvV1_*.root
    fi
    if [ $HAVELHC17f2afix == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13b2efix2";
        rm $OUTPUTDIR_LHC17f2a_fix/*/GammaConvV1_*.root
    fi
fi
