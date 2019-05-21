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
HAVETOBUILDLHC17j7=0
HAVELHC17j7a=1
HAVELHC17j7b=1
HAVELHC17j7c=1
HAVETOBUILDLHC18d2=0
HAVELHC18d2a=1
HAVELHC18d2b=1
HAVELHC18d2c=1
HAVELHC18d2d=1

# default trainconfigurations
LHC17nData="";
LHC17j7MC="";
LHC17j7MCa="";
LHC17j7MCb="";
LHC17j7MCc="";
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
# LHC17j7MC="936";
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";
# LHC18d2MC="943"; # PCM
# LHC18d2MCa="child_1";
# LHC18d2MCb="child_2";
# LHC18d2MCc="child_3";
# LHC18d2MCd="child_4";

# TRAINDIR=Legotrain-vAN20180612-EMCAndPHOSAndPCM
# ISAOD=0
# LHC17nData="412"; #pass 1
# LHC17nData="407"; #pass 1
# LHC17j7MC="935"; # EMC nl
# LHC17j7MC="940"; # PHOS nl
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";
# LHC18d2MC="950"; # PHOS nl
# LHC18d2MC="944"; # PHOS nl
# LHC18d2MCa="child_1";
# LHC18d2MCb="child_2";
# LHC18d2MCc="child_3";
# LHC18d2MCd="child_4";

TRAINDIR=20190520-PCMtrees-QA
ISAOD=1
# LHC17nData="510";
LHC18d2MC="1317";
LHC18d2MCa="child_1";
LHC18d2MCb="child_2";
LHC18d2MCc="child_3";
LHC18d2MCd="child_4";


OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ "$LHC17nData" == "" ]; then
    HAVELHC17n=0;
fi

if [ "$LHC17j7MC" != "" ]; then
    HAVETOBUILDLHC17j7=1;
fi

if [ "$LHC17j7MCa" == "" ]; then
    HAVELHC17j7a=0;
    echo $LHC17j7MCa
fi
if [ "$LHC17j7MCb" == "" ]; then
    HAVELHC17j7b=0;
    echo $LHC17j7MCb
fi
if [ "$LHC17j7MCc" == "" ]; then
    HAVELHC17j7c=0;
    echo $LHC17j7MCc
fi

if [ "$LHC18d2MC" != "" ]; then
    HAVETOBUILDLHC18d2=1;
fi

if [ "$LHC18d2MCa" == "" ]; then
    HAVELHC18d2a=0;
    echo $LHC18d2MCa
fi
if [ "$LHC18d2MCb" == "" ]; then
    HAVELHC18d2b=0;
    echo $LHC18d2MCb
fi
if [ "$LHC18d2MCc" == "" ]; then
    HAVELHC18d2c=0;
    echo $LHC18d2MCc
fi
if [ "$LHC18d2MCd" == "" ]; then
    HAVELHC18d2d=0;
    echo $LHC18d2MCd
fi


mkdir -p $OUTPUTDIR/CutSelections
if [ $ISAOD -eq 1 ]; then
    pathTrainData="GA_PbPb_AOD";
    pathTrainMC="GA_PbPb_MC_AOD";
    addName="AOD"
else
    pathTrainData="GA_PbPb";
    pathTrainMC="GA_PbPb_MC";
fi

if [ $HAVELHC17n == 1 ]; then
    LHC17nData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainData/ | grep $LHC17nData\_`
    if [ "$LHC17nData" == "" ]; then
        HAVELHC17n=0;
    else
        OUTPUTDIR_LHC17n=$BASEDIR/$TRAINDIR/$pathTrainData-$LHC17nData
    fi
    echo $OUTPUTDIR_LHC17n
fi


if [ $HAVELHC17j7a == 1 ]; then
    if [ $HAVETOBUILDLHC17j7 == 1 ]; then
        LHC17j7MCa=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC17j7MC\_ | grep $LHC17j7MCa`
    else
        LHC17j7MCa=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC17j7MCa\_`
    fi
    echo $LHC17j7MCa
    if [ "$LHC17j7MCa" == "" ]; then
        HAVELHC17j7a=0;
    else
        OUTPUTDIR_LHC17j7a=$BASEDIR/$TRAINDIR/$pathTrainMC-$LHC17j7MCa
        echo $OUTPUTDIR_LHC17j7a
    fi
fi

if [ $HAVELHC17j7b == 1 ]; then
    if [ $HAVETOBUILDLHC17j7 == 1 ]; then
        LHC17j7MCb=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC17j7MC\_ | grep $LHC17j7MCb`
    else
        LHC17j7MCb=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC17j7MCb\_`
    fi
    echo $LHC17j7MCb
    if [ "$LHC17j7MCb" == "" ]; then
        HAVELHC17j7b=0;
    else
        OUTPUTDIR_LHC17j7b=$BASEDIR/$TRAINDIR/$pathTrainMC-$LHC17j7MCb
        echo $OUTPUTDIR_LHC17j7b
    fi
fi
if [ $HAVELHC17j7c == 1 ]; then
    if [ $HAVETOBUILDLHC17j7 == 1 ]; then
        LHC17j7MCc=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC17j7MC\_ | grep $LHC17j7MCc`
    else
        LHC17j7MCc=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC17j7MCc\_`
    fi
    echo $LHC17j7MCc
    if [ "$LHC17j7MCc" == "" ]; then
        HAVELHC17j7c=0;
    else
        OUTPUTDIR_LHC17j7c=$BASEDIR/$TRAINDIR/$pathTrainMC-$LHC17j7MCc
        echo $OUTPUTDIR_LHC17j7c
    fi
fi

if [ $HAVELHC18d2a == 1 ]; then
    if [ $HAVETOBUILDLHC18d2 == 1 ]; then
        LHC18d2MCa=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MC\_ | grep $LHC18d2MCa`
    else
        LHC18d2MCa=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MCa\_`
    fi
    echo $LHC18d2MCa
    if [ "$LHC18d2MCa" == "" ]; then
        HAVELHC18d2a=0;
    else
        OUTPUTDIR_LHC18d2a=$BASEDIR/$TRAINDIR/$pathTrainMC-$LHC18d2MCa
        echo $OUTPUTDIR_LHC18d2a
    fi
fi

if [ $HAVELHC18d2b == 1 ]; then
    if [ $HAVETOBUILDLHC18d2 == 1 ]; then
        LHC18d2MCb=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MC\_ | grep $LHC18d2MCb`
    else
        LHC18d2MCb=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MCb\_`
    fi
    echo $LHC18d2MCb
    if [ "$LHC18d2MCb" == "" ]; then
        HAVELHC18d2b=0;
    else
        OUTPUTDIR_LHC18d2b=$BASEDIR/$TRAINDIR/$pathTrainMC-$LHC18d2MCb
        echo $OUTPUTDIR_LHC18d2b
    fi
fi
if [ $HAVELHC18d2c == 1 ]; then
    if [ $HAVETOBUILDLHC18d2 == 1 ]; then
        LHC18d2MCc=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MC\_ | grep $LHC18d2MCc`
    else
        LHC18d2MCc=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MCc\_`
    fi
    echo $LHC18d2MCc
    if [ "$LHC18d2MCc" == "" ]; then
        HAVELHC18d2c=0;
    else
        OUTPUTDIR_LHC18d2c=$BASEDIR/$TRAINDIR/$pathTrainMC-$LHC18d2MCc
        echo $OUTPUTDIR_LHC18d2c
    fi
fi
if [ $HAVELHC18d2d == 1 ]; then
    if [ $HAVETOBUILDLHC18d2 == 1 ]; then
        LHC18d2MCd=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MC\_ | grep $LHC18d2MCc`
    else
        LHC18d2MCd=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/ | grep $LHC18d2MCd\_`
    fi
    echo $LHC18d2MCd
    if [ "$LHC18d2MCd" == "" ]; then
        HAVELHC18d2d=0;
    else
        OUTPUTDIR_LHC18d2d=$BASEDIR/$TRAINDIR/$pathTrainMC-$LHC18d2MCd
        echo $OUTPUTDIR_LHC18d2d
    fi
fi


if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC17n == 1 ]; then
        echo "downloading LHC17n"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17n_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17n/$runNumber "/alice/data/2017/LHC17n/000$runNumber/pass$passNr/PWGGA/$pathTrainData/$LHC17nData" $NSlashes3 "/alice/data/2017/LHC17n/000$runNumber/pass$passNr/PWGGA/$pathTrainData/$LHC17nData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC17n/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC17n/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17n_all.txt`
                ls $OUTPUTDIR_LHC17n/$firstrunNumber/GammaConvV1_*.root > fileLHC17n.txt
                MergeAccordingToSpecificRunlist fileLHC17n.txt $OUTPUTDIR_LHC17n $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC17n_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17n "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainData/$LHC17nData/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi

    currentDir=$PWD
    if [ $HAVELHC17j7a == 1 ]; then
        echo "downloading LHC17j7"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17j7_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17j7a/$runNumber "/alice/sim/2017/LHC17j7/$runNumber/PWGGA/$pathTrainMC/$LHC17j7MCa" $NSlashes3 "/alice/sim/2017/LHC17j7/$runNumber/PWGGA/$pathTrainMC/$LHC17j7MCa/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17j7a/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17j7a/GammaConvV1*.root*
                echo runlists/runNumbersLHC17j7_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17j7_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17j7a/$firstrunNumber/GammaConvV1_*.root > fileLHC17j7.txt
                MergeAccordingToSpecificRunlist fileLHC17j7.txt $OUTPUTDIR_LHC17j7a $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC17j7_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17j7a "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC17j7MCa/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi
    if [ $HAVELHC17j7b == 1 ]; then
        echo "downloading LHC17j7"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17j7_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17j7b/$runNumber "/alice/sim/2017/LHC17j7_ZDCfix/$runNumber/PWGGA/$pathTrainMC/$LHC17j7MCb" $NSlashes3 "/alice/sim/2017/LHC17j7_ZDCfix/$runNumber/PWGGA/$pathTrainMC/$LHC17j7MCb/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17j7b/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17j7b/GammaConvV1*.root*
                echo runlists/runNumbersLHC17j7_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17j7_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17j7b/$firstrunNumber/GammaConvV1_*.root > fileLHC17j7.txt
                MergeAccordingToSpecificRunlist fileLHC17j7.txt $OUTPUTDIR_LHC17j7b $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC17j7_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17j7b "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC17j7MCb/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi
    if [ $HAVELHC17j7c == 1 ]; then
        echo "downloading LHC17j7"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17j7_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17j7c/$runNumber "/alice/sim/2017/LHC17j7_ZDCfix_extra/$runNumber/PWGGA/$pathTrainMC/$LHC17j7MCc" $NSlashes3 "/alice/sim/2017/LHC17j7_ZDCfix_extra/$runNumber/PWGGA/$pathTrainMC/$LHC17j7MCc/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17j7c/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17j7c/GammaConvV1*.root*
                echo runlists/runNumbersLHC17j7_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17j7_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17j7c/$firstrunNumber/GammaConvV1_*.root > fileLHC17j7.txt
                MergeAccordingToSpecificRunlist fileLHC17j7.txt $OUTPUTDIR_LHC17j7c $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC17j7_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17j7c "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC17j7MCc/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi

      if [ $HAVELHC18d2a == 1 ]; then
        echo "downloading LHC18d2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC18d2_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC18d2a/$runNumber "/alice/sim/2018/LHC18d2_1/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCa" $NSlashes3 "/alice/sim/2018/LHC18d2_1/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCa/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC18d2a/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC18d2a/GammaConvV1*.root*
                echo runlists/runNumbersLHC18d2_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC18d2_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC18d2a/$firstrunNumber/GammaConvV1_*.root > fileLHC18d2.txt
                MergeAccordingToSpecificRunlist fileLHC18d2.txt $OUTPUTDIR_LHC18d2a $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC18d2_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18d2a "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC18d2MCa/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi
    if [ $HAVELHC18d2b == 1 ]; then
        echo "downloading LHC18d2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC18d2_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC18d2b/$runNumber "/alice/sim/2018/LHC18d2_2/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCb" $NSlashes3 "/alice/sim/2018/LHC18d2_2/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCb/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC18d2b/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC18d2b/GammaConvV1*.root*
                echo runlists/runNumbersLHC18d2_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC18d2_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC18d2b/$firstrunNumber/GammaConvV1_*.root > fileLHC18d2.txt
                MergeAccordingToSpecificRunlist fileLHC18d2.txt $OUTPUTDIR_LHC18d2b $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC18d2_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18d2b "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC18d2MCb/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi
    if [ $HAVELHC18d2c == 1 ]; then
        echo "downloading LHC18d2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC18d2_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC18d2c/$runNumber "/alice/sim/2018/LHC18d2_3/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCc" $NSlashes3 "/alice/sim/2018/LHC18d2_3/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCc/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC18d2c/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC18d2c/GammaConvV1*.root*
                echo runlists/runNumbersLHC18d2_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC18d2_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC18d2c/$firstrunNumber/GammaConvV1_*.root > fileLHC18d2.txt
                MergeAccordingToSpecificRunlist fileLHC18d2.txt $OUTPUTDIR_LHC18d2c $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC18d2_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18d2c "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC18d2MCc/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi
    if [ $HAVELHC18d2d == 1 ]; then
        echo "downloading LHC18d2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC18d2_all.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC18d2d/$runNumber "/alice/sim/2018/LHC18d2_4/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCd" $NSlashes3 "/alice/sim/2018/LHC18d2_4/$runNumber/PWGGA/$pathTrainMC/$LHC18d2MCd/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC18d2d/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC18d2d/GammaConvV1*.root*
                echo runlists/runNumbersLHC18d2_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC18d2_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC18d2d/$firstrunNumber/GammaConvV1_*.root > fileLHC18d2.txt
                MergeAccordingToSpecificRunlist fileLHC18d2.txt $OUTPUTDIR_LHC18d2d $NSlashes3 GammaConvV1 All$addName runlists/runNumbersLHC18d2_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18d2d "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC18d2MCd/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi


    if [ $HAVELHC17n == 1 ]; then
        ls $OUTPUTDIR_LHC17n/GammaConvV1-All$addName_*.root > fileLHC17n.txt
        fileNumbers=`cat fileLHC17n.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17n $NSlashes "LHC17n-pass$passNr-All$addName" "-All$addName"
        done;
    fi

    if [ $HAVELHC17j7a == 1 ]; then
        ls $OUTPUTDIR_LHC17j7a/GammaConvV1-All$addName_*.root > fileLHC17j7a.txt
        fileNumbers=`cat fileLHC17j7a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17j7a $NSlashes "MC_LHC17j7-All$addName" "-All$addName"
        done;
    fi
    if [ $HAVELHC17j7b == 1 ]; then
        ls $OUTPUTDIR_LHC17j7b/GammaConvV1-All$addName_*.root > fileLHC17j7b.txt
        fileNumbers=`cat fileLHC17j7b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17j7b $NSlashes "MC_LHC17j7ZDCfix-All$addName" "-All$addName"
        done;
    fi
    if [ $HAVELHC17j7c == 1 ]; then
        ls $OUTPUTDIR_LHC17j7c/GammaConvV1-All$addName_*.root > fileLHC17j7c.txt
        fileNumbers=`cat fileLHC17j7c.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17j7c $NSlashes "MC_LHC17j7ZDCfixExtra-All$addName" "-All$addName"
        done;
    fi

    if [ $HAVELHC18d2a == 1 ]; then
        ls $OUTPUTDIR_LHC18d2a/GammaConvV1-All$addName_*.root > fileLHC18d2a.txt
        fileNumbers=`cat fileLHC18d2a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2a $NSlashes "MC_LHC18d2-4-All$addName" "-All$addName"
        done;
    fi
    if [ $HAVELHC18d2b == 1 ]; then
        ls $OUTPUTDIR_LHC18d2b/GammaConvV1-All$addName_*.root > fileLHC18d2b.txt
        fileNumbers=`cat fileLHC18d2b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2b $NSlashes "MC_LHC18d2-1-All$addName" "-All$addName"
        done;
    fi
    if [ $HAVELHC18d2c == 1 ]; then
        ls $OUTPUTDIR_LHC18d2c/GammaConvV1-All$addName_*.root > fileLHC18d2c.txt
        fileNumbers=`cat fileLHC18d2c.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2c $NSlashes "MC_LHC18d2-2-All$addName" "-All$addName"
        done;
    fi
    if [ $HAVELHC18d2d == 1 ]; then
        ls $OUTPUTDIR_LHC18d2d/GammaConvV1-All$addName_*.root > fileLHC18d2d.txt
        fileNumbers=`cat fileLHC18d2d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC18d2d $NSlashes "MC_LHC18d2-3-All$addName" "-All$addName"
        done;
    fi

    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaConvV1_MC_LHC17j7-All$addName\_*.root > filesForMerging.txt
        echo -e "\nAll" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            rm listCurrMerge.txt
            fileA="$OUTPUTDIR/GammaConvV1_MC_LHC17j7-All$addName""_$number.root"
            fileB="$OUTPUTDIR/GammaConvV1_MC_LHC17j7ZDCfix-All$addName""_$number.root"
            fileC="$OUTPUTDIR/GammaConvV1_MC_LHC17j7ZDCfixExtra-All$addName""_$number.root"
            echo -e "$fileA\n$fileB\n$fileC" > listCurrMerge.txt
            MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvV1_MC_LHC17j7x-All$addName\_$number.root
        done
        ls $OUTPUTDIR/GammaConvV1_MC_LHC18d2-4-All$addName\_*.root > filesForMerging.txt
        echo -e "\nAll" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            rm listCurrMerge.txt
            fileA="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-4-All$addName""_$number.root"
            fileB="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-1-All$addName""_$number.root"
            fileC="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-2-All$addName""_$number.root"
            fileD="$OUTPUTDIR/GammaConvV1_MC_LHC18d2-3-All$addName""_$number.root"
            echo -e "$fileA\n$fileB\n$fileC\n$fileD" > listCurrMerge.txt
            MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvV1_MC_LHC18d2x-All$addName\_$number.root
        done
    fi


else
    if [ $HAVELHC17n == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC17n";
        rm $OUTPUTDIR_LHC17n/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC17n/*/*/*GammaConvV1_*.root
    fi
    if [ $HAVELHC17j7a == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC17j7";
        rm $OUTPUTDIR_LHC17j7a/*/GammaConvV1_*.root
    fi
    if [ $HAVELHC17j7b == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC17j7ZDCfix";
        rm $OUTPUTDIR_LHC17j7b/*/GammaConvV1_*.root
    fi
    if [ $HAVELHC17j7c == 1 ]; then
        echo "removing all GammaConvV1 files in runFolders for LHC17j7ZDCfixExtra";
        rm $OUTPUTDIR_LHC17j7c/*/GammaConvV1_*.root
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
