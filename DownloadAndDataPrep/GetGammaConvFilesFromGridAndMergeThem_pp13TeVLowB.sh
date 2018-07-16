

#! /bin/bash
source basicFunction.sh

# Switches to enable/disable certain procedures
DOWNLOADON=1
MERGEONDATA=0
# Switch on download per run
SINGLERUN=1

MERGEONMC=0
MERGEONMCJJ=0
MERGEONBINSSingle=1
MERGEONBINS=1
passNr=1

# Check if train configuration has actually been given
HAVELHC17g=0
HAVETOBUILDData=0
HAVETOBUILD18c=0
HAVELHC18c_woSDD=0
HAVELHC18c_wSDD=0

# MonteCarlo files
HAVETOBUILDMC=0
HAVELHC17h3MC=0

# Default trainconfigurations
LHClowB="";
LHC17gData="";
LHClowBMC="";
LHC17h3MC="";
LHC18clowB="";
LHC18c_woSDD="";
LHC18c_wSDD="";

# Copies files from grid
# Creates directory
# Changes internal structure
# Merges files according to the needs
if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp
elif [ $1 = "jlietave" ]; then
    BASEDIR=/home/jaklie/Gridoutput
fi
echo $BASEDIR

# Definitition of number of slashes in your path to different depths
    NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
    NSlashes=`expr $NSlashesBASE + 4`
    NSlashes2=`expr $NSlashes - 1`
    NSlashes3=`expr $NSlashes + 1`
    NSlashes4=`expr $NSlashes + 2`
    echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"

TRAINDIR=Legotrain-vAN-20170612-LowBtest
LHClowB="2376"
# LHC17gData="child_2";
LHClowBMC="3340";
# LHC17h3MC="child_2";
LHC18clowB="2396";
LHC18c_woSDD="child_3";
LHC18c_wSDD="child_5";


OUTPUTDIR=$BASEDIR/$TRAINDIR

# Decides what to download
    if [ "$LHClowB" != "" ]; then
        HAVETOBUILDData=1;
    fi

    if [ "$LHC17gData" != "" ]; then
        HAVELHC17g=1;
    fi

    if [ "$LHClowBMC" != "" ]; then
        HAVETOBUILDMC=1;
    fi

    if [ "$LHC17h3MC" != "" ]; then
        HAVELHC17h3MC=1;
    fi

    if [ "$LHC18clowB" != "" ]; then
        HAVETOBUILD18c=1;
    fi

    if [ "$LHC18c_woSDD" != "" ]; then
        HAVELHC18c_woSDD=1;
    fi

    if [ "$LHC18c_wSDD" != "" ]; then
        HAVELHC18c_wSDD=1;
    fi

# Creates directory
    if [ $HAVELHC17g == 1 ]; then
	    echo "outputdir LHC17g"
        if [ $HAVETOBUILDData == 1 ]; then
            LHC17gData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp | grep $LHClowB\_ | grep $LHC17gData`
        else
            LHC17gData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp | grep $LHC17gData\_`
        fi
        OUTPUTDIR_LHC17g=$BASEDIR/$TRAINDIR/GA_pp-$LHC17gData
        mkdir -p $OUTPUTDIR_LHC17g
    fi

    if [ $HAVELHC17h3MC == 1 ]; then
	    echo "outputdir LHC17h3"
        if [ $HAVETOBUILDMC == 1 ]; then
            LHC17h3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC | grep $LHClowBMC\_ | grep $LHC17h3MC`
        else
            LHC17h3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC | grep $LHC17h3MC\_`
        fi
        OUTPUTDIR_LHC17h3MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17h3MC
        mkdir -p $OUTPUTDIR_LHC17h3MC
    fi

    if [ $HAVELHC18c_woSDD == 1 ]; then
    echo "outputdir LHC18c_woSDD"
        if [ $HAVETOBUILD18c == 1 ]; then
            LHC18c_woSDD=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp | grep $LHC18clowB\_ | grep $LHC18c_woSDD`
        else
            LHC18c_woSDD=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp | grep $LHC18c_woSDD\_`
        fi
        OUTPUTDIR_LHC18c_woSDD=$BASEDIR/$TRAINDIR/GA_pp-$LHC18c_woSDD
        mkdir -p $OUTPUTDIR_LHC18c_woSDD
    fi

    if [ $HAVELHC18c_wSDD == 1 ]; then
	    echo "outputdir LHC18c_wSDD"
        if [ $HAVETOBUILD18c == 1 ]; then
            LHC18c_wSDD=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp | grep $LHC18clowB\_ | grep $LHC18c_wSDD`
        else
            LHC18c_wSDD=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp | grep $LHC18c_wSDD\_`
        fi
        OUTPUTDIR_LHC18c_wSDD=$BASEDIR/$TRAINDIR/GA_pp-$LHC18c_wSDD
        mkdir -p $OUTPUTDIR_LHC18c_wSDD
    fi

# Download:
    if [ $DOWNLOADON == 1 ]; then

    # LHC17g
        if [ $HAVELHC17g == 1 ]; then
            echo "downloading LHC17g"
            if [ $SINGLERUN == 1 ]; then
                runNumbers=`cat runlists/runNumbersLHC17g_dpgTracks.txt`
                echo $runNumbers
                for runNumber in $runNumbers; do
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC17g/$runNumber "/alice/data/2017/LHC17g/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC17gData" $NSlashes3 "/alice/data/2017/LHC17g/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC17gData/" kTRUE
                done;
                if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC17g/mergedAllConv.txt ]; then
                    #rm $OUTPUTDIR_LHC17g/GammaConvV1*.root*
                    firstrunNumber=`head -n1 runlists/runNumbersLHC17g_dpgTracks.txt`
                    ls $OUTPUTDIR_LHC17g/$firstrunNumber/GammaConvV1_*.root > fileLHC17g.txt

                    MergeAccordingToSpecificRunlist fileLHC17g.txt $OUTPUTDIR_LHC17g $NSlashes3 GammaConvV1 DPGTrack runlists/runNumbersLHC17g_dpgTracks.txt "no"
                    rm fileNumbers4.txt
                    GetFileNumberListPCM $OUTPUTDIR_LHC17g $NSlashes fileNumbers4.txt
                    fileNumbers=`cat fileNumbers4.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        SeparateCutsIfNeeded $OUTPUTDIR_LHC17g/GammaConvV1-DPGTrack_$fileNumber 0 kTRUE
                    done;

                fi
            else
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC17gData/merge_runlist_1" DPGTrack $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC17gData/merge_runlist_2" DPGTrackAndCalo $NSlashes3 "" kFALSE
            fi
        fi

    # LHC18c_woSDD
        if [ $HAVELHC18c_woSDD == 1 ]; then
            echo "downloading LHC18c_woSDD"
            if [ $SINGLERUN == 1 ]; then
                runNumbers=`cat runlists/runNumbersLHC18c_woSDD_all.txt`
                echo $runNumbers
                for runNumber in $runNumbers; do
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC18c_woSDD/$runNumber "/alice/data/2018/LHC18c/000$runNumber/pass"$passNr"_CENT_woSDD/PWGGA/GA_pp/$LHC18c_woSDD" $NSlashes3 "/alice/data/2018/LHC18c/000$runNumber/pass"$passNr"_CENT_woSDD/PWGGA/GA_pp/$LHC18c_woSDD/" kTRUE
                done;
                if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC18c_woSDD/mergedAllConv.txt ]; then
                    #rm $OUTPUTDIR_LHC18c_woSDD/GammaConvV1*.root*
                    firstrunNumber=`head -n1 runlists/runNumbersLHC18c_woSDD_all.txt`
                    ls $OUTPUTDIR_LHC18c_woSDD/$firstrunNumber/GammaConvV1_*.root > fileLHC18c.txt

                    MergeAccordingToSpecificRunlist fileLHC18c.txt $OUTPUTDIR_LHC18c_woSDD $NSlashes3 GammaConvV1 all runlists/runNumbersLHC18c_woSDD_all.txt "no"
                    rm fileNumbers4.txt
                    GetFileNumberListPCM $OUTPUTDIR_LHC18c_woSDD $NSlashes fileNumbers4.txt
                    fileNumbers=`cat fileNumbers4.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        SeparateCutsIfNeeded $OUTPUTDIR_LHC18c_woSDD/GammaConvV1-all_$fileNumber 0 kTRUE
                    done;

                fi
            else
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18c_woSDD "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC18c_woSDD/merge" all $NSlashes3 "" kFALSE
            fi
        fi

    #LHC18c_wSDD
        if [ $HAVELHC18c_wSDD == 1 ]; then
            echo "downloading LHC18c_wSDD"
            if [ $SINGLERUN == 1 ]; then
                runNumbers=`cat runlists/runNumbersLHC18c_wSDD_all.txt`
                echo $runNumbers
                for runNumber in $runNumbers; do
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC18c_wSDD/$runNumber "/alice/data/2018/LHC18c/000$runNumber/pass"$passNr"_CENT/PWGGA/GA_pp/$LHC18c_wSDD" $NSlashes3 "/alice/data/2018/LHC18c/000$runNumber/pass"$passNr"_CENT/PWGGA/GA_pp/$LHC18c_wSDD/" kTRUE
                done;
                if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC18c_wSDD/mergedAllConv.txt ]; then
                    #rm $OUTPUTDIR_LHC18c_wSDD/GammaConvV1*.root*
                    firstrunNumber=`head -n1 runlists/runNumbersLHC18c_wSDD_all.txt`
                    ls $OUTPUTDIR_LHC18c_wSDD/$firstrunNumber/GammaConvV1_*.root > fileLHC18c.txt

                    MergeAccordingToSpecificRunlist fileLHC18c.txt $OUTPUTDIR_LHC18c_wSDD $NSlashes3 GammaConvV1 all runlists/runNumbersLHC18c_wSDD_all.txt "no"
                    rm fileNumbers4.txt
                    GetFileNumberListPCM $OUTPUTDIR_LHC18c_wSDD $NSlashes fileNumbers4.txt
                    fileNumbers=`cat fileNumbers4.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        SeparateCutsIfNeeded $OUTPUTDIR_LHC18c_wSDD/GammaConvV1-all_$fileNumber 0 kTRUE
                    done;

                fi
            else
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18c_wSDD "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC18c_wSDD/merge" all $NSlashes3 "" kFALSE
            fi
        fi

    # MC outputs
        if [ $HAVELHC17h3MC == 1 ]; then
            echo "downloading LHC17h3MC"
            if [ $SINGLERUN == 1 ]; then
                runNumbers=`cat runlists/runNumbersLHC17h3_DPGTrack.txt`
                echo $runNumbers
                for runNumber in $runNumbers; do
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC17h3/$runNumber 
                    "/alice/sim/2017/LHC17h3$5/$runNumber/PWGGA/GA_pp_MC/$LHC17h3MC" $NSlashes3 "/alice/sim/2017/LHC17h3$5/$runNumber/PWGGA/GA_pp_MC/$LHC17h3MC/" kTRUE
                done;
                if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17h3/mergedAllConv.txt ]; then
                    cd $currentDir
                    rm $OUTPUTDIR_LHC17h3/GammaConvV1*.root*
                    echo runlists/runNumbersLHC17h3_$3_dpgTracks.txt
                    firstrunNumber=`head -n1 runlists/runNumbersLHC17h3_$3_dpgTracks.txt`
                    echo $firstrunNumber
                    ls $OUTPUTDIR_LHC17h3/$firstrunNumber/GammaConvV1_*.root > fileLHC17h3.txt
                    MergeAccordingToSpecificRunlist fileLHC17h3.txt $OUTPUTDIR_LHC17h3 $NSlashes3 GammaConvV1 DPGTrack runlists/runNumbersLHC17h3_$3_dpgTracks.txt "no"
                    GetFileNumberListPCM $OUTPUTDIR_LHC17h3 $NSlashes fileNumbers4.txt
                    cat fileNumbers4.txt
                    fileNumbers=`cat fileNumbers4.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        SeparateCutsIfNeeded $OUTPUTDIR_LHC17h3/GammaConvV1-DPGTrack_$fileNumber 0 kTRUE
                    done;
                fi
            else
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17h3MC "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17h3MC/merge_runlist_1" DPGTrack $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17h3MC "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17h3MC/merge_runlist_2" DPGTrackAndCalo $NSlashes3 "" kFALSE
            fi
        fi

    fi

# Changes internal structure
    if [ $HAVELHC17g == 1 ]; then
        ls $OUTPUTDIR_LHC17g/GammaConvV1-DPGTrack_*.root > fileLHC17g.txt
        fileNumbers=`cat fileLHC17g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17g $NSlashes "LHC17g-pass$passNr-DPGTrack" "-DPGTrack"
        done;
    fi

    if [ $HAVELHC17h3MC == 1 ]; then
        ls $OUTPUTDIR_LHC17h3MC/GammaConvV1-DPGTrack_*.root > fileLHC17h3MC.txt
        fileNumbers=`cat fileLHC17h3MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17h3MC $NSlashes "LHC17h3MC-pass$passNr-DPGTrack" "-DPGTrack"
        done;
    fi

    if [ $HAVELHC18c_woSDD == 1 ]; then
        ls $OUTPUTDIR_LHC18c_woSDD/GammaConvV1-all_*.root > fileLHC18c_woSDD.txt
        fileNumbers=`cat fileLHC18c_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18c_woSDD $NSlashes "LHC18c_woSDD-pass$passNr-all" "-all"
        done;
    fi

    if [ $HAVELHC18c_wSDD == 1 ]; then
        ls $OUTPUTDIR_LHC18c_wSDD/GammaConvV1-all_*.root > fileLHC18c_wSDD.txt
        fileNumbers=`cat fileLHC18c_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18c_wSDD $NSlashes "LHC18c_wSDD-pass$passNr-all" "-all"
        done;
    fi

# ChangeStructureToStandard
   if [ $HAVELHC17h3MC == 1 ]; then
        ls $OUTPUTDIR_LHC17h3MC/GammaConvV1_*.root > fileLHC17h3MC.txt
        fileNumbers=`cat fileLHC17h3MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17h3MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17h3_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17h3_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17h3_$number.log\"\)
        done;
    fi

# Merging files
