#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=0
MERGEONFASTAndWOSDD=1
SINGLERUN=1
SEPARATEON=0
MERGEONSINGLEData=0
MERGEONSINGLEMC=1
CLEANUP=1
CLEANUPMAYOR=$2
number=""

# check if train configuration has actually been given
HAVELHC16q=1
HAVELHC16t=1
HAVETOBUILDData=0
HAVELHC17f2b=1
HAVETOBUILDLHC17f2b=0
HAVELHC17f2afix=1
HAVETOBUILDLHC17f2afix=0
# default trainconfigurations
LHC16qData="";
LHC16tData="";
LHC16qtData="";
LHC17f2bMC="";
LHC17f2bMCchild="";
LHC17f2a_fixMC="";
LHC17f2a_fixMCchild="";
passNr="1";

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorageExternal/OutputLegoTrains/pPb
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

TRAINDIR=Legotrain-vAN20180206-EMC
# woSDD (CENT) PHOS
# LHC16qtData="718"; #pass 2
# LHC16qData="child_3"; #pass 3
# LHC16tData="child_4"; #pass 2
# LHC17f2bMC="1168";
# LHC17f2bMCchild="child_2";
# LHC17f2a_fixMC="1167";
# LHC17f2a_fixMCchild="child_2";

# LHC17f2bMC="1171";
# LHC17f2bMCchild="child_2";
# LHC17f2a_fixMC="1170";
# LHC17f2a_fixMCchild="child_2";

# FAST PHOS
LHC16qtData="718"; #pass 2
LHC16qData="child_1"; #pass 3
LHC16tData="child_2"; #pass 2
# LHC17f2bMC="1168";
# LHC17f2bMCchild="child_1";
# LHC17f2a_fixMC="1167";
# LHC17f2a_fixMCchild="child_1";

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
if [ "$LHC17f2bMCchild" != "" ]; then
    HAVETOBUILDLHC17f2b=1;
fi
if [ "$LHC17f2a_fixMC" == "" ]; then
    HAVELHC17f2afix=0;
fi
if [ "$LHC17f2a_fixMCchild" != "" ]; then
    HAVETOBUILDLHC17f2afix=1;
fi

mkdir -p $OUTPUTDIR/CutSelections


if [ $HAVELHC16q == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\_ | grep $LHC16qData`
#         LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\- | grep $LHC16qData`
    else
        LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qData\_`
#         LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qData\-`
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
#         LHC16tData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\- | grep $LHC16tData`
    else
        LHC16tData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16tData\_`
#         LHC16tData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16tData\-`
    fi
    if [ "$LHC16tData" == "" ]; then
        HAVELHC16t=0;
    else
        OUTPUTDIR_LHC16t=$BASEDIR/$TRAINDIR/GA_pPb-$LHC16tData
    fi
    echo $OUTPUTDIR_LHC16t
fi

if [ $HAVELHC17f2b == 1 ]; then
    if [ $HAVETOBUILDLHC17f2b == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMC\_
        LHC17f2bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMC\_ | grep $LHC17f2bMCchild`
    else
        LHC17f2bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMC\_`
    fi
    echo $LHC17f2bMC
    if [ "$LHC17f2bMC" == "" ]; then
        HAVELHC17f2b=0;
    else
        OUTPUTDIR_LHC17f2b=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2bMC
        echo $OUTPUTDIR_LHC17f2b
    fi
fi
if [ $HAVELHC17f2afix == 1 ]; then
    if [ $HAVETOBUILDLHC17f2afix == 1 ]; then
        LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMC\_ | grep $LHC17f2a_fixMCchild`
    else
        LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMC\_`
#     LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMC\-`
    fi
    echo $LHC17f2a_fixMC
    if [ "$LHC17f2a_fixMC" == "" ]; then
        HAVELHC17f2afix=0;
    else
        OUTPUTDIR_LHC17f2a_fix=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2a_fixMC
    fi
    echo $LHC17f2a_fixMC
    echo $OUTPUTDIR_LHC17f2a_fix
fi

if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC16q == 1 ]; then
        echo "downloading LHC16q"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16q_$3_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16q/$runNumber "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16q/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16q/GammaCaloMerged*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16q_$3_dpgTracks.txt`
                ls $OUTPUTDIR_LHC16q/$firstrunNumber/GammaCaloMerged_*.root > fileLHC16q.txt

#                 MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16q $NSlashes3 GammaCaloMerged All runlists/runNumbersLHC16q_$3_all.txt
                MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16q $NSlashes3 GammaCaloMerged DPGTrack runlists/runNumbersLHC16q_$3_dpgTracks.txt
                MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16q $NSlashes3 GammaCaloMerged DPGTrackAndCalo runlists/runNumbersLHC16q_$3_dpgTracksAndCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16q "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16qData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "downloading LHC16t"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16t_$3_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16t/$runNumber "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16t/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16t/GammaCaloMerged*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16t_$3_dpgTracks.txt`
                ls $OUTPUTDIR_LHC16t/$firstrunNumber/GammaCaloMerged_*.root > fileLHC16t.txt
                fileNumbers=`cat fileLHC16t.txt`
#                 MergeAccordingToSpecificRunlist fileLHC16t.txt $OUTPUTDIR_LHC16t $NSlashes3 GammaCaloMerged All runlists/runNumbersLHC16t_$3_all.txt
                MergeAccordingToSpecificRunlist fileLHC16t.txt $OUTPUTDIR_LHC16t $NSlashes3 GammaCaloMerged DPGTrack runlists/runNumbersLHC16t_$3_dpgTracks.txt
                MergeAccordingToSpecificRunlist fileLHC16t.txt $OUTPUTDIR_LHC16t $NSlashes3 GammaCaloMerged DPGTrackAndCalo runlists/runNumbersLHC16t_$3_dpgTracksAndCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16t "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16tData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi

    currentDir=$PWD
    if [ $HAVELHC17f2b == 1 ]; then
        echo "downloading LHC17f2b"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f2b.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2b/$runNumber "/alice/sim/2017/LHC17f2b$5/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMC" $NSlashes3 "/alice/sim/2017/LHC17f2b$5/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMC/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17f2b/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f2b/GammaCaloMerged*.root*
                echo runlists/runNumbersLHC17f2b_$3_dpgTracks.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2b_$3_dpgTracks.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17f2b/$firstrunNumber/GammaCaloMerged_*.root > fileLHC17f2b.txt
#                 MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged All runlists/runNumbersLHC17f2b_$3_all.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged DPGTrack runlists/runNumbersLHC17f2b_$3_dpgTracks.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged DPGTrackAndCalo runlists/runNumbersLHC17f2b_$3_dpgTracks.txt
#                 MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged All-LHC16q runlists/runNumbersLHC16q_woSDD_all.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged DPGTrack-LHC16q runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16q.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged DPGTrackAndCalo-LHC16q runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16q.txt
#                 MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged All-LHC16t runlists/runNumbersLHC16t_woSDD_all.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged DPGTrack-LHC16t runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16t.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaCaloMerged DPGTrackAndCalo-LHC16t runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16t.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17f2afix == 1 ]; then
        echo "downloading LHC17f2a_fix"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f2a_fix_$3_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2a_fix/$runNumber "/alice/sim/2017/LHC17f2a$5_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC" $NSlashes3 "/alice/sim/2017/LHC17f2a$5_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17f2a_fix/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f2a_fix/GammaCaloMerged*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2a_fix_$3_dpgTracks.txt`
                ls $OUTPUTDIR_LHC17f2a_fix/$firstrunNumber/GammaCaloMerged_*.root > fileLHC17f2a_fix.txt
                fileNumbers=`cat fileLHC17f2a_fix.txt`
#                 MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged All runlists/runNumbersLHC17f2a_fix_$3_all.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged DPGTrack runlists/runNumbersLHC17f2a_fix_$3_dpgTracks.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged DPGTrackAndCalo runlists/runNumbersLHC17f2a_fix_$3_dpgTracksAndCalo.txt
#                 MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged All-LHC16q runlists/runNumbersLHC16q_woSDD_all.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged DPGTrack-LHC16q runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracks-LHC16q.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged DPGTrackAndCalo-LHC16q runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracksAndCalo-LHC16q.txt
#                 MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged All-LHC16t runlists/runNumbersLHC16t_woSDD_all.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged DPGTrack-LHC16t runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracks-LHC16t.txt
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaCaloMerged DPGTrackAndCalo-LHC16t runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracksAndCalo-LHC16t.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2a_fix "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/merge" $NSlashes "" kTRUE
        fi
    fi

    if [ $HAVELHC16q == 1 ]; then
        ls $OUTPUTDIR_LHC16q/GammaCaloMerged-All_*.root > fileLHC16q.txt
        fileNumbers=`cat fileLHC16q.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC16q $NSlashes "LHC16q_$3-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16q/GammaCaloMerged-DPGTrack_*.root > fileLHC16q.txt
        fileNumbers=`cat fileLHC16q.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC16q $NSlashes "LHC16q_$3-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16q/GammaCaloMerged-DPGTrackAndCalo_*.root > fileLHC16q.txt
        fileNumbers=`cat fileLHC16q.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC16q $NSlashes "LHC16q_$3-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16t == 1 ]; then
        ls $OUTPUTDIR_LHC16t/GammaCaloMerged-All_*.root > fileLHC16t.txt
        fileNumbers=`cat fileLHC16t.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC16t $NSlashes "LHC16t_$3-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16t/GammaCaloMerged-DPGTrack_*.root > fileLHC16t.txt
        fileNumbers=`cat fileLHC16t.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC16t $NSlashes "LHC16t_$3-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16t/GammaCaloMerged-DPGTrackAndCalo_*.root > fileLHC16t.txt
        fileNumbers=`cat fileLHC16t.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC16t $NSlashes "LHC16t_$3-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17f2b == 1 ]; then
        ls $OUTPUTDIR_LHC17f2b/GammaCaloMerged-All_*.root > fileLHC17f2b.txt
        fileNumbers=`cat fileLHC17f2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC17f2b $NSlashes "MC_LHC17f2b_$3-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17f2b/GammaCaloMerged-DPGTrack_*.root > fileLHC17f2b.txt
        fileNumbers=`cat fileLHC17f2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC17f2b $NSlashes "MC_LHC17f2b_$3-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17f2b/GammaCaloMerged-DPGTrackAndCalo_*.root > fileLHC17f2b.txt
        fileNumbers=`cat fileLHC17f2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC17f2b $NSlashes "MC_LHC17f2b_$3-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17f2afix == 1 ]; then
        ls $OUTPUTDIR_LHC17f2a_fix/GammaCaloMerged-All_*.root > fileLHC17f2a_fix.txt
        fileNumbers=`cat fileLHC17f2a_fix.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC17f2a_fix $NSlashes "MC_LHC17f2a_fix_$3-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17f2a_fix/GammaCaloMerged-DPGTrack_*.root > fileLHC17f2a_fix.txt
        fileNumbers=`cat fileLHC17f2a_fix.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC17f2a_fix $NSlashes "MC_LHC17f2a_fix_$3-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17f2a_fix/GammaCaloMerged-DPGTrackAndCalo*.root > fileLHC17f2a_fix.txt
        fileNumbers=`cat fileLHC17f2a_fix.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCaloMerged $fileName $OUTPUTDIR_LHC17f2a_fix $NSlashes "MC_LHC17f2a_fix_$3-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaCaloMerged_LHC16q_$3-pass$passNr-DPGTrackAndCalo\_*.root > filesForMerging.txt
        echo -e "DPGTrackAndCalo\nDPGTrack" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        listsToMerge=`cat runlistsToMerge.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileQ="$OUTPUTDIR/GammaCaloMerged_LHC16q_$3-pass$passNr-$runListName""_$number.root"
                fileT="$OUTPUTDIR/GammaCaloMerged_LHC16t_$3-pass$passNr-$runListName""_$number.root"
                echo -e "$fileQ\n$fileT" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCaloMerged_LHC16qt_$3-pass$passNr-$runListName\_$number.root
            done
        done
    fi

    if [ $MERGEONFASTAndWOSDD == 1 ]; then
        ls $OUTPUTDIR/GammaCaloMerged_LHC16q_fast-pass$passNr-DPGTrackAndCalo\_*.root > filesForMerging.txt
        echo -e "DPGTrackAndCalo\nDPGTrack" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        listsToMerge=`cat runlistsToMerge.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaCaloMerged_LHC16qt_fast-pass$passNr-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaCaloMerged_LHC16qt_woSDD-pass$passNr-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCaloMerged_LHC16qt_fast-woSDD-pass$passNr-$runListName\_$number.root
            done
        done
    fi

else
    if [ $HAVELHC16q == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16q";
#         rm $OUTPUTDIR_LHC16q/*/GammaCaloMerged_*.root
        rm $OUTPUTDIR_LHC16q/*/*/*GammaCaloMerged_*.root
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16t";
#         rm $OUTPUTDIR_LHC16t/*/GammaCaloMerged_*.root
        rm $OUTPUTDIR_LHC16t/*/*/*GammaCaloMerged_*.root
    fi

    if [ $HAVELHC17f2b == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17f2b";
#         rm $OUTPUTDIR_LHC17f2b/*/GammaCaloMerged_*.root
        rm -rf $OUTPUTDIR_LHC17f2b/*/Stage*
    fi
    if [ $HAVELHC17f2afix == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13b2efix2";
#         rm $OUTPUTDIR_LHC17f2a_fix/*/GammaCaloMerged_*.root
        rm -rf $OUTPUTDIR_LHC17f2a_fix/*/Stage*
    fi
fi
