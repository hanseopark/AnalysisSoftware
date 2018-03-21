#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=1
MERGEONFASTAndWOSDD=1
SINGLERUN=1
SEPARATEON=1
MERGEONSINGLEData=1
MERGEONSINGLEMC=1
CLEANUP=1
CLEANUPMAYOR=$2
number=""
FAST="_FAST"
# check if train configuration has actually been given
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
LHC16qtData="718"; #pass 2
LHC16qDataFast="child_1"; #pass 3
LHC16tDataFast="child_2"; #pass 2
LHC16qData="child_3"; #pass 3
LHC16tData="child_4"; #pass 2


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
# LHC17g8aMCMoth="1216"
# LHC17g8aMC="child_2"
# LHC17g8aMCFast="child_1"


OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ "$LHC16qData" == "" ]; then
    HAVELHC16q=0;
fi
if [ "$LHC16tData" = "" ]; then
    HAVELHC16t=0;
fi
if [ "$LHC16qDataFast" == "" ]; then
    HAVELHC16qF=0;
fi
if [ "$LHC16tDataFast" = "" ]; then
    HAVELHC16tF=0;
fi
if [ "$LHC16qtData" != "" ]; then
    HAVETOBUILDData=1;
fi

if [ "$LHC17f2bMC" == "" ]; then
    HAVELHC17f2b=0;
fi
if [ "$LHC17f2bMCFast" == "" ]; then
    HAVELHC17f2bF=0;
fi
if [ "$LHC17f2bMCMoth" != "" ]; then
    HAVETOBUILDLHC17f2b=1;
fi
if [ "$LHC17f2a_fixMC" == "" ]; then
    HAVELHC17f2afix=0;
fi
if [ "$LHC17f2a_fixMCFast" == "" ]; then
    HAVELHC17f2afixF=0;
fi
if [ "$LHC17f2a_fixMCMoth" != "" ]; then
    HAVETOBUILDLHC17f2afix=1;
fi
if [ "$LHC17g8aMC" == "" ]; then
    HAVELHC17g8a=0;
fi
if [ "$LHC17g8aMCFast" == "" ]; then
    HAVELHC17g8aF=0;
fi
if [ "$LHC17g8aMCMoth" != "" ]; then
    HAVETOBUILDLHC17g8a=1;
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
if [ $HAVELHC16qF == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16qDataFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\_ | grep $LHC16qDataFast`
#         LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\- | grep $LHC16qData`
    else
        LHC16qDataFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qDataFast\_`
#         LHC16qData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qData\-`
    fi
    if [ "$LHC16qDataFast" == "" ]; then
        HAVELHC16qF=0;
    else
        OUTPUTDIR_LHC16qF=$BASEDIR/$TRAINDIR/GA_pPb-$LHC16qDataFast
    fi
    echo $OUTPUTDIR_LHC16qF
fi
if [ $HAVELHC16tF == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16tDataFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\_ | grep $LHC16tDataFast`
#         LHC16tData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16qtData\- | grep $LHC16tData`
    else
        LHC16tDataFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16tDataFast\_`
#         LHC16tData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC16tData\-`
    fi
    if [ "$LHC16tDataFast" == "" ]; then
        HAVELHC16tF=0;
    else
        OUTPUTDIR_LHC16tF=$BASEDIR/$TRAINDIR/GA_pPb-$LHC16tDataFast
    fi
    echo $OUTPUTDIR_LHC16tF
fi


if [ $HAVELHC17f2b == 1 ]; then
    if [ $HAVETOBUILDLHC17f2b == 1 ]; then
        LHC17f2bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMCMoth\_ | grep $LHC17f2bMC`
    else
        LHC17f2bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMC\_`
    fi
    if [ "$LHC17f2bMC" == "" ]; then
        HAVELHC17f2b=0;
    else
        OUTPUTDIR_LHC17f2b=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2bMC
    fi
fi
if [ $HAVELHC17f2bF == 1 ]; then
    if [ $HAVETOBUILDLHC17f2b == 1 ]; then
        LHC17f2bMCFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMCMoth\_ | grep $LHC17f2bMCFast`
    else
        LHC17f2bMCFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2bMCFast\_`
    fi
    if [ "$LHC17f2bMCFast" == "" ]; then
        HAVELHC17f2bF=0;
    else
        OUTPUTDIR_LHC17f2bF=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2bMCFast
    fi
fi
if [ $HAVELHC17f2afix == 1 ]; then
    if [ $HAVETOBUILDLHC17f2afix == 1 ]; then
        LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMCMoth\_ | grep $LHC17f2a_fixMC`
    else
        LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMC\_`
#     LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMC\-`
    fi
    if [ "$LHC17f2a_fixMC" == "" ]; then
        HAVELHC17f2afix=0;
    else
        OUTPUTDIR_LHC17f2a_fix=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2a_fixMC
    fi
fi
if [ $HAVELHC17f2afixF == 1 ]; then
    if [ $HAVETOBUILDLHC17f2afix == 1 ]; then
        LHC17f2a_fixMCFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMCMoth\_ | grep $LHC17f2a_fixMCFast`
    else
        LHC17f2a_fixMCFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMCFast\_`
#     LHC17f2a_fixMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17f2a_fixMC\-`
    fi
    if [ "$LHC17f2a_fixMCFast" == "" ]; then
        HAVELHC17f2afixF=0;
    else
        OUTPUTDIR_LHC17f2a_fixF=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17f2a_fixMCFast
    fi
fi

if [ $HAVELHC17g8a == 1 ]; then
    if [ $HAVETOBUILDLHC17g8a == 1 ]; then
        LHC17g8aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17g8aMCMoth\_ | grep $LHC17g8aMC`
    else
        LHC17g8aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17g8aMC\_`
#     LHC17g8aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17g8aMC\-`
    fi
    if [ "$LHC17g8aMC" == "" ]; then
        HAVELHC17g8a=0;
    else
        OUTPUTDIR_LHC17g8a=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17g8aMC
    fi
fi
if [ $HAVELHC17g8aF == 1 ]; then
    if [ $HAVETOBUILDLHC17g8a == 1 ]; then
        LHC17g8aMCFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17g8aMCMoth\_ | grep $LHC17g8aMCFast`
    else
        LHC17g8aMCFast=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17g8aMCFast\_`
#     LHC17g8aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC17g8aMC\-`
    fi
    if [ "$LHC17g8aMCFast" == "" ]; then
        HAVELHC17g8aF=0;
    else
        OUTPUTDIR_LHC17g8aF=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC17g8aMCFast
    fi
fi

if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC16q == 1 ]; then
        echo "downloading LHC16q"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16q_$3_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC16q/$runNumber "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData/Stage_1/" kTRUE
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16q/$runNumber "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16qData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16q/mergedAllConv.txt ]; then
#                 rm $OUTPUTDIR_LHC16q/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16q_$3_dpgTracks.txt`
                ls $OUTPUTDIR_LHC16q/$firstrunNumber/GammaConv_*.root > fileLHC16q.txt

                MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16q $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC16q_$3_dpgTracks.txt "no"
                rm fileNumbers4.txt
                GetFileNumberListPCM $OUTPUTDIR_LHC16q $NSlashes fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC16q/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;

            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16q "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16qData/merge_runlist_1" All $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16qF == 1 ]; then
        echo "downloading LHC16q fast"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16q_fast_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC16qF/$runNumber "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16qDataFast" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16qDataFast/Stage_1/" kTRUE
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16qF/$runNumber "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16qDataFast" $NSlashes3 "/alice/data/2016/LHC16q/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16qDataFast/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16qF/mergedAllConv.txt ]; then
#                 rm $OUTPUTDIR_LHC16qF/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16q_fast_dpgTracks.txt`
                ls $OUTPUTDIR_LHC16qF/$firstrunNumber/GammaConv_*.root > fileLHC16q.txt
                MergeAccordingToSpecificRunlist fileLHC16q.txt $OUTPUTDIR_LHC16qF $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC16q_fast_dpgTracks.txt "no"
                rm fileNumbers4.txt
                GetFileNumberListPCM $OUTPUTDIR_LHC16qF $NSlashes fileNumbers4.txt
                cat fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC16qF/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16qF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16qDataFast/merge_runlist_1" All $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "downloading LHC16t"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16t_$3_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC16t/$runNumber "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData/Stage_1/" kTRUE
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16t/$runNumber "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$4/PWGGA/GA_pPb/$LHC16tData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16t/mergedAllConv.txt ]; then
#                 rm $OUTPUTDIR_LHC16t/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16t_$3_dpgTracks.txt`
                ls $OUTPUTDIR_LHC16t/$firstrunNumber/GammaConv_*.root > fileLHC16t.txt
                fileNumbers=`cat fileLHC16t.txt`
                MergeAccordingToSpecificRunlist fileLHC16t.txt $OUTPUTDIR_LHC16t $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC16t_$3_dpgTracks.txt "no"
                rm fileNumbers4.txt
                GetFileNumberListPCM $OUTPUTDIR_LHC16t $NSlashes fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC16t/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16t "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16tData/merge_runlist_1" All $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16tF == 1 ]; then
        echo "downloading LHC16t fast"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16t_fast_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC16tF/$runNumber "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16tDataFast" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16tDataFast/Stage_1/" kTRUE
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16tF/$runNumber "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16tDataFast" $NSlashes3 "/alice/data/2016/LHC16t/000$runNumber/pass$passNr$FAST/PWGGA/GA_pPb/$LHC16tDataFast/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16tF/mergedAllConv.txt ]; then
#                 rm $OUTPUTDIR_LHC16tF/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16t_fast_dpgTracks.txt`
                ls $OUTPUTDIR_LHC16tF/$firstrunNumber/GammaConv_*.root > fileLHC16t.txt
                fileNumbers=`cat fileLHC16t.txt`
                MergeAccordingToSpecificRunlist fileLHC16t.txt $OUTPUTDIR_LHC16tF $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC16t_fast_dpgTracks.txt "no"
                rm fileNumbers4.txt
                GetFileNumberListPCM $OUTPUTDIR_LHC16tF $NSlashes fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC16tF/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;

            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16tF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC16tDataFast/merge_runlist_1" All $NSlashes "" kTRUE
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
                rm $OUTPUTDIR_LHC17f2b/GammaConv*.root*
                echo runlists/runNumbersLHC17f2b_$3_dpgTracks.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2b_$3_dpgTracks.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17f2b/$firstrunNumber/GammaConv_*.root > fileLHC17f2b.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC17f2b_$3_dpgTracks.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConv DPGTrack-LHC16q runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16q.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2b $NSlashes3 GammaConv DPGTrack-LHC16t runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16t.txt "no"
                GetFileNumberListPCM $OUTPUTDIR_LHC17f2b $NSlashes fileNumbers4.txt
                cat fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC17f2b/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMC/merge_runlist_2" DPGTrack $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMC/merge_runlist_4" DPGTrack-LHC16q $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMC/merge_runlist_6" DPGTrack-LHC16t $NSlashes3 "" kTRUE
        fi
    fi
    if [ $HAVELHC17f2bF == 1 ]; then
        echo "downloading LHC17f2b fast"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f2b.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2bF/$runNumber "/alice/sim/2017/LHC17f2b_fast/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMCFast" $NSlashes3 "/alice/sim/2017/LHC17f2b_fast/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2bMCFast/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17f2bF/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f2bF/GammaConv*.root*
                echo runlists/runNumbersLHC17f2b_fast_dpgTracks.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2b_fast_dpgTracks.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17f2bF/$firstrunNumber/GammaConv_*.root > fileLHC17f2b.txt
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2bF $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC17f2b_fast_dpgTracks.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2bF $NSlashes3 GammaConv DPGTrack-LHC16q runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16q.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2b.txt $OUTPUTDIR_LHC17f2bF $NSlashes3 GammaConv DPGTrack-LHC16t runlists/runNumbersLHC17f2b_woSDD_dpgTracks-LHC16t.txt "no"
                GetFileNumberListPCM $OUTPUTDIR_LHC17f2bF $NSlashes fileNumbers4.txt
                cat fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC17f2bF/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2bF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMCFast/merge_runlist_2" DPGTrack $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2bF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMCFast/merge_runlist_4" DPGTrack-LHC16q $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2bF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2bMCFast/merge_runlist_6" DPGTrack-LHC16t $NSlashes3 "" kTRUE
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
                rm $OUTPUTDIR_LHC17f2a_fix/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2a_fix_$3_dpgTracks.txt`
                ls $OUTPUTDIR_LHC17f2a_fix/$firstrunNumber/GammaConv_*.root > fileLHC17f2a_fix.txt
                fileNumbers=`cat fileLHC17f2a_fix.txt`
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC17f2a_fix_$3_dpgTracks.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConv DPGTrack-LHC16q runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracks-LHC16q.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fix $NSlashes3 GammaConv DPGTrack-LHC16t runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracks-LHC16t.txt "no"
                GetFileNumberListPCM $OUTPUTDIR_LHC17f2a_fix $NSlashes fileNumbers4.txt
                cat fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC17f2a_fix/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2a_fix "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/merge_runlist_2" DPGTrack $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2a_fix "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/merge_runlist_4" DPGTrack-LHC16q $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2a_fix "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/merge_runlist_6" DPGTrack-LHC16t $NSlashes3 "" kTRUE
        fi
    fi
    if [ $HAVELHC17f2afixF == 1 ]; then
        echo "downloading LHC17f2a_fix fast"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f2a_fix_fast_dpgTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f2a_fixF/$runNumber "/alice/sim/2017/LHC17f2a_fast_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast" $NSlashes3 "/alice/sim/2017/LHC17f2a_fast_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17f2a_fixF/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f2a_fixF/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f2a_fix_fast_dpgTracks.txt`
                ls $OUTPUTDIR_LHC17f2a_fixF/$firstrunNumber/GammaConv_*.root > fileLHC17f2a_fix.txt
                fileNumbers=`cat fileLHC17f2a_fix.txt`
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fixF $NSlashes3 GammaConv DPGTrack runlists/runNumbersLHC17f2a_fix_fast_dpgTracks.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fixF $NSlashes3 GammaConv DPGTrack-LHC16q runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracks-LHC16q.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f2a_fix.txt $OUTPUTDIR_LHC17f2a_fixF $NSlashes3 GammaConv DPGTrack-LHC16t runlists/runNumbersLHC17f2a_fix_woSDD_dpgTracks-LHC16t.txt "no"
                GetFileNumberListPCM $OUTPUTDIR_LHC17f2a_fixF $NSlashes fileNumbers4.txt
                cat fileNumbers4.txt
                fileNumbers=`cat fileNumbers4.txt`
                for fileNumber in $fileNumbers; do
                    echo $fileNumber
                    SeparateCutsIfNeeded $OUTPUTDIR_LHC17f2a_fixF/GammaConv-DPGTrack_$fileNumber 0 kTRUE
                done;
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2a_fixF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast/merge_runlist_2" DPGTrack $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2a_fixF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast/merge_runlist_4" DPGTrack-LHC16q $NSlashes3 "" kTRUE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f2a_fixF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast/merge_runlist_6" DPGTrack-LHC16t $NSlashes3 "" kTRUE
        fi
    fi

    if [ $HAVELHC17g8a == 1 ]; then
        echo "LHC17g8a"
            if [ $SINGLERUN == 1 ]; then
                runNumbers=`cat runlists/runNumbersLHC17g8a_cent_woSDD.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binsJetJetLHC17g8a_cent_woSDD.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC17g8a/$binNumber/$runNumber "/alice/sim/2017/LHC17g8a$5/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMC" $NSlashes3 "/alice/sim/2017/LHC17g8a$5/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMC/" kTRUE
                done;
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17g8a/mergedAllConv.txt ]; then
                echo "HERERRE"
                cd $currentDir
                rm $OUTPUTDIR_LHC17g8a/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17g8a_cent_woSDD.txt`
                firstbinNumber=`head -n1 runlists/binsJetJetLHC17g8a_cent_woSDD.txt`
                ls $OUTPUTDIR_LHC17g8a/$firstbinNumber/$firstrunNumber/GammaConv_*.root > fileLHC17g8a.txt
                fileNumbers=`cat fileLHC17g8a.txt`
                MergeAccordingToSpecificRunlist fileLHC17g8a.txt $OUTPUTDIR_LHC17g8a $NSlashes4 GammaConv DPGTrack runlists/runNumbersLHC16qt_dpgTracks.txt runlists/binsJetJetLHC17g8a_cent_woSDD.txt
#                 MergeAccordingToSpecificRunlist fileLHC17g8a.txt $OUTPUTDIR_LHC17g8a $NSlashes4 GammaConv DPGTrack-LHC16q runlists/runNumbersLHC16q_dpgTracks.txt runlists/binsJetJetLHC17g8a_cent_woSDD.txt
#                 MergeAccordingToSpecificRunlist fileLHC17g8a.txt $OUTPUTDIR_LHC17g8a $NSlashes4 GammaConv DPGTrack-LHC16t runlists/runNumbersLHC16t_dpgTracks.txt runlists/binsJetJetLHC17g8a_cent_woSDD.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17g8a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17g8aMC/merge" All $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17g8aF == 1 ]; then
        echo "LHC17g8a fast"
            if [ $SINGLERUN == 1 ]; then
                runNumbers=`cat runlists/runNumbersLHC17g8a_fast.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binsJetJetLHC17g8a_fast.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC17g8aF/$binNumber/$runNumber "/alice/sim/2017/LHC17g8a_fast/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMCFast" $NSlashes3 "/alice/sim/2017/LHC17g8a_fast/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMCFast/" kTRUE
                done;
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17g8aF/mergedAllConv.txt ]; then
                echo "HERERRE"
                cd $currentDir
                rm $OUTPUTDIR_LHC17g8aF/GammaConv*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17g8a_cent_woSDD.txt`
                firstbinNumber=`head -n1 runlists/binsJetJetLHC17g8a_cent_woSDD.txt`
                ls $OUTPUTDIR_LHC17g8aF/$firstbinNumber/$firstrunNumber/GammaConv_*.root > fileLHC17g8a.txt
                fileNumbers=`cat fileLHC17g8a.txt`
                MergeAccordingToSpecificRunlist fileLHC17g8a.txt $OUTPUTDIR_LHC17g8aF $NSlashes4 GammaConv DPGTrack runlists/runNumbersLHC16qt_dpgTracks.txt runlists/binsJetJetLHC17g8a_fast.txt
#                 MergeAccordingToSpecificRunlist fileLHC17g8a.txt $OUTPUTDIR_LHC17g8aF $NSlashes4 GammaConv DPGTrack-LHC16q runlists/runNumbersLHC16q_dpgTracks.txt runlists/binsJetJetLHC17g8a_fast.txt
#                 MergeAccordingToSpecificRunlist fileLHC17g8a.txt $OUTPUTDIR_LHC17g8aF $NSlashes4 GammaConv DPGTrack-LHC16t runlists/runNumbersLHC16t_dpgTracks.txt runlists/binsJetJetLHC17g8a_fast.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17g8aF "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC17g8aMCFast/merge" All $NSlashes "" kTRUE
        fi
    fi


    if [ $HAVELHC16q == 1 ]; then
        ls $OUTPUTDIR_LHC16q/GammaConv-DPGTrack_*.root > fileLHC16q.txt
        fileNumbers=`cat fileLHC16q.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16q $NSlashes "LHC16q_$3-pass$passNr-DPGTrack" "-DPGTrack"
        done;
    fi
    if [ $HAVELHC16qF == 1 ]; then
        ls $OUTPUTDIR_LHC16qF/GammaConv-DPGTrack_*.root > fileLHC16q.txt
        fileNumbers=`cat fileLHC16q.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16qF $NSlashes "LHC16q_fast-pass$passNr-DPGTrack" "-DPGTrack"
        done;
    fi

    if [ $HAVELHC16t == 1 ]; then
        ls $OUTPUTDIR_LHC16t/GammaConv-DPGTrack_*.root > fileLHC16t.txt
        fileNumbers=`cat fileLHC16t.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16t $NSlashes "LHC16t_$3-pass$passNr-DPGTrack" "-DPGTrack"
        done;
    fi
    if [ $HAVELHC16tF == 1 ]; then
        ls $OUTPUTDIR_LHC16tF/GammaConv-DPGTrack_*.root > fileLHC16t.txt
        fileNumbers=`cat fileLHC16t.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16tF $NSlashes "LHC16t_fast-pass$passNr-DPGTrack" "-DPGTrack"
        done;
    fi


    if [ $HAVELHC17f2b == 1 ]; then
        ls $OUTPUTDIR_LHC17f2b/GammaConv-DPGTrack_*.root > fileLHC17f2b.txt
        fileNumbers=`cat fileLHC17f2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f2b $NSlashes "MC_LHC17f2b_$3-DPGTrack" "-DPGTrack"
        done;
    fi
    if [ $HAVELHC17f2bF == 1 ]; then
        ls $OUTPUTDIR_LHC17f2bF/GammaConv-DPGTrack_*.root > fileLHC17f2b.txt
        fileNumbers=`cat fileLHC17f2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f2bF $NSlashes "MC_LHC17f2b_fast-DPGTrack" "-DPGTrack"
        done;
    fi

    if [ $HAVELHC17f2afix == 1 ]; then
        ls $OUTPUTDIR_LHC17f2a_fix/GammaConv-DPGTrack_*.root > fileLHC17f2a_fix.txt
        fileNumbers=`cat fileLHC17f2a_fix.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f2a_fix $NSlashes "MC_LHC17f2a_fix_$3-DPGTrack" "-DPGTrack"
        done;
    fi
    if [ $HAVELHC17f2afixF == 1 ]; then
        ls $OUTPUTDIR_LHC17f2a_fixF/GammaConv-DPGTrack_*.root > fileLHC17f2a_fix.txt
        fileNumbers=`cat fileLHC17f2a_fix.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f2a_fixF $NSlashes "MC_LHC17f2a_fix_fast-DPGTrack" "-DPGTrack"
        done;
    fi
    if [ $HAVELHC17g8a == 1 ]; then
        ls $OUTPUTDIR_LHC17g8a/GammaConv-DPGTrack_*.root > fileLHC17g8a.txt
        fileNumbers=`cat fileLHC17g8a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17g8a $NSlashes "MC_LHC17g8a_$3-DPGTrack" "-DPGTrack"
        done;
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            ls $OUTPUTDIR_LHC17g8a/GammaConv-DPGTrack_*.root > fileLHC17g8a.txt
            fileNumbers=`cat fileLHC17g8a.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes)) 2
                cp $OUTPUTDIR_LHC17g8a/$binNumber/GammaConv-DPGTrack_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConv_MC_LHC17g8a-$binNumber\_$3-DPGTrack_$number.root
            done
        done;
    fi
    if [ $HAVELHC17g8aF == 1 ]; then
        ls $OUTPUTDIR_LHC17g8aF/GammaConv-DPGTrack_*.root > fileLHC17g8a.txt
        fileNumbers=`cat fileLHC17g8a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17g8aF $NSlashes "MC_LHC17g8a_fast-DPGTrack" "-DPGTrack"
        done;
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            ls $OUTPUTDIR_LHC17g8aF/GammaConv-DPGTrack_*.root > fileLHC17g8a.txt
            fileNumbers=`cat fileLHC17g8a.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes)) 2
                cp $OUTPUTDIR_LHC17g8aF/$binNumber/GammaConv-DPGTrack_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConv_MC_LHC17g8a-$binNumber\_fast-DPGTrack_$number.root
            done
        done;

    fi

    if [ $MERGEON == 1 ]; then
        echo -e "$3\nfast" > listReconstruction.txt
        listReconstruction=`cat listReconstruction.txt`
        for reco in $listReconstruction; do
            ls $OUTPUTDIR/GammaConv_LHC16q_$reco-pass$passNr-DPGTrack\_*.root > filesForMerging.txt
            echo -e "DPGTrack" > runlistsToMerge.txt
            filesForMerging=`cat filesForMerging.txt`
            listsToMerge=`cat runlistsToMerge.txt`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 4
                echo $number
                for runListName in $listsToMerge; do
                    rm listCurrMerge.txt
                    fileQ="$OUTPUTDIR/GammaConv_LHC16q_$reco-pass$passNr-$runListName""_$number.root"
                    fileT="$OUTPUTDIR/GammaConv_LHC16t_$reco-pass$passNr-$runListName""_$number.root"
                    echo -e "$fileQ\n$fileT" > listCurrMerge.txt
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConv_LHC16qt_$reco-pass$passNr-$runListName\_$number.root
                done
            done
        done
    fi

    if [ $MERGEONFASTAndWOSDD == 1 ]; then
        ls $OUTPUTDIR/GammaConv_LHC16qt_fast-pass$passNr-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        echo -e "DPGTrack" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        listsToMerge=`cat runlistsToMerge.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConv_LHC16qt_fast-pass$passNr-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConv_LHC16qt_woSDD-pass$passNr-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConv_LHC16qt_fast-woSDD-pass$passNr-$runListName\_$number.root
            done
        done

        ls $OUTPUTDIR/GammaConv_MC_LHC17f2a_fix_fast-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 6
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConv_MC_LHC17f2a_fix_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConv_MC_LHC17f2a_fix_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConv_MC_LHC17f2a_fix_fast-woSDD-$runListName\_$number.root
            done
        done

        ls $OUTPUTDIR/GammaConv_MC_LHC17f2b_fast-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConv_MC_LHC17f2b_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConv_MC_LHC17f2b_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConv_MC_LHC17f2b_fast-woSDD-$runListName\_$number.root
            done
        done

        ls $OUTPUTDIR/GammaConv_MC_LHC17f2b_fast-woSDD-DPGTrack\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConv_MC_LHC17f2b_fast-woSDD-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConv_MC_LHC17f2a_fix_fast-woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConv_MC_LHC17f2a_fix_LHC17f2b_fast-woSDD-$runListName\_$number.root
            done
        done


        ls $OUTPUTDIR/GammaConv_MC_LHC17g8a_fast-DPGTrack\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConv_MC_LHC17g8a_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConv_MC_LHC17g8a_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConv_MC_LHC17g8a_fast-woSDD-$runListName\_$number.root
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    rm listCurrMerge.txt
                    fileF="$OUTPUTDIR/JJMCSingleBins/GammaConv_MC_LHC17g8a-$binNumber""_fast-$runListName""_$number.root"
                    fileW="$OUTPUTDIR/JJMCSingleBins/GammaConv_MC_LHC17g8a-$binNumber""_woSDD-$runListName""_$number.root"
                    echo -e "$fileF\n$fileW" > listCurrMerge.txt
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/JJMCSingleBins/GammaConv_MC_LHC17g8a-$binNumber\_fast-woSDD-$runListName\_$number.root
                done
            done
        done
    fi


else
    if [ $HAVELHC16q == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16q";
#         rm $OUTPUTDIR_LHC16q/*/GammaConv_*.root
        rm $OUTPUTDIR_LHC16q/*/*/*GammaConv_*.root
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16t";
#         rm $OUTPUTDIR_LHC16t/*/GammaConv_*.root
        rm $OUTPUTDIR_LHC16t/*/*/*GammaConv_*.root
    fi
    if [ $HAVELHC16qF == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16q";
#         rm $OUTPUTDIR_LHC16q/*/GammaConv_*.root
        rm $OUTPUTDIR_LHC16qF/*/*/*GammaConv_*.root
    fi
    if [ $HAVELHC16tF == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16t";
#         rm $OUTPUTDIR_LHC16t/*/GammaConv_*.root
        rm $OUTPUTDIR_LHC16tF/*/*/*GammaConv_*.root
    fi

    if [ $HAVELHC17f2b == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17f2b";
#         rm $OUTPUTDIR_LHC17f2b/*/GammaConv_*.root
        rm -rf $OUTPUTDIR_LHC17f2b/*/Stage*
    fi
    if [ $HAVELHC17f2afix == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17f2a_fix";
#         rm $OUTPUTDIR_LHC17f2a_fix/*/GammaConv_*.root
        rm -rf $OUTPUTDIR_LHC17f2a_fix/*/Stage*
    fi
    if [ $HAVELHC17f2bF == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17f2b";
#         rm $OUTPUTDIR_LHC17f2b/*/GammaConv_*.root
        rm -rf $OUTPUTDIR_LHC17f2bF/*/Stage*
    fi
    if [ $HAVELHC17f2afixF == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17f2a_fix";
#         rm $OUTPUTDIR_LHC17f2a_fix/*/GammaConv_*.root
        rm -rf $OUTPUTDIR_LHC17f2a_fixF/*/Stage*
    fi
    if [ $HAVELHC17g8a == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17g8a";
#         rm $OUTPUTDIR_LHC17g8a/*/GammaConv_*.root
        rm -rf $OUTPUTDIR_LHC17g8a/*/*/Stage*
    fi
    if [ $HAVELHC17g8aF == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17g8a";
#         rm $OUTPUTDIR_LHC17g8a/*/GammaConv_*.root
        rm -rf $OUTPUTDIR_LHC17g8aF/*/*/Stage*
    fi
fi
