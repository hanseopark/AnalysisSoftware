#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=1
SINGLERUN=0
SEPARATEON=0
MERGEONSINGLEData=0
MERGEONSINGLEMC=0
CLEANUP=1
CLEANUPMAYOR=$2
number=""

# check if train configuration has actually been given
HAVELHC17n=1
HAVETOBUILDLHC17j7=0
HAVELHC17j7a=1
HAVELHC17j7b=1
HAVELHC17j7c=1

# default trainconfigurations
LHC17nData="";
LHC17j7MC="";
LHC17j7MCa="";
LHC17j7MCb="";
LHC17j7MCc="";
passNr="1";

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/XeXe
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


# TRAINDIR=Legotrain-vAN20171115-XeXeQA
# LHC17nData="348"; #pass 1
#  LHC17nData="349"; #pass 1
# LHC17j7MC="687";
# LHC17j7MC="688";
# LHC17j7MC="689";
# LHC17j7MC="690";

# TRAINDIR=Legotrain-vAN20180201-XeXeQA
# LHC17nData="382"; #pass 1
# LHC17j7MC="811";

# TRAINDIR=Legotrain-vAN20180207-XeXeQA
# LHC17nData="387"; #pass 1
# LHC17j7MC="816";
# LHC17j7MC="821";
# LHC17j7MC="822";

# TRAINDIR=Legotrain-vAN20180226-MCmorestat
# LHC17nData="392"; #pass 1
# LHC17j7MC="835";
# LHC17j7MC="836";
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";
# ISAOD=1
# LHC17nData="448"; #pass 1
# LHC17j7MC="874";
# LHC17j7MC="873";
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";

# TRAINDIR=Legotrain-vAN20180228-nonlinV1
# ISAOD=0
# LHC17nData="392"; #pass 1
# LHC17j7MC="836";
# LHC17j7MC="837";
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";
# ISAOD=1
# LHC17nData="449"; #pass 1
# LHC17j7MC="875";
# LHC17j7MC="877";
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";


TRAINDIR=Legotrain-vAN20180228-EMClowth
# ISAOD=0
# LHC17nData="392"; #pass 1
# LHC17j7MC="836";
# LHC17j7MC="837";
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";
ISAOD=1
LHC17nData="450"; #pass 1
LHC17j7MC="876";
LHC17j7MCa="child_1";
LHC17j7MCb="child_2";
LHC17j7MCc="child_3";

TRAINDIR=Legotrain-vAN20180228-EMClowth
ISAOD=0
LHC17nData="395"; #pass 1
LHC17j7MC="906";
LHC17j7MCa="child_1";
LHC17j7MCb="child_2";
LHC17j7MCc="child_3";
# ISAOD=1
# LHC17nData="450"; #pass 1
# LHC17j7MC="876";
# LHC17j7MCa="child_1";
# LHC17j7MCb="child_2";
# LHC17j7MCc="child_3";


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
                rm $OUTPUTDIR_LHC17n/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17n_all.txt`
                ls $OUTPUTDIR_LHC17n/$firstrunNumber/GammaCalo_*.root > fileLHC17n.txt
                MergeAccordingToSpecificRunlist fileLHC17n.txt $OUTPUTDIR_LHC17n $NSlashes3 GammaCalo All$addName runlists/runNumbersLHC17n_all.txt "no"
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
                rm $OUTPUTDIR_LHC17j7a/GammaCalo*.root*
                echo runlists/runNumbersLHC17j7_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17j7_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17j7a/$firstrunNumber/GammaCalo_*.root > fileLHC17j7.txt
                MergeAccordingToSpecificRunlist fileLHC17j7.txt $OUTPUTDIR_LHC17j7a $NSlashes3 GammaCalo All$addName runlists/runNumbersLHC17j7_all.txt "no"
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
                rm $OUTPUTDIR_LHC17j7b/GammaCalo*.root*
                echo runlists/runNumbersLHC17j7_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17j7_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17j7b/$firstrunNumber/GammaCalo_*.root > fileLHC17j7.txt
                MergeAccordingToSpecificRunlist fileLHC17j7.txt $OUTPUTDIR_LHC17j7b $NSlashes3 GammaCalo All$addName runlists/runNumbersLHC17j7_all.txt "no"
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
                rm $OUTPUTDIR_LHC17j7c/GammaCalo*.root*
                echo runlists/runNumbersLHC17j7_all.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17j7_all.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17j7c/$firstrunNumber/GammaCalo_*.root > fileLHC17j7.txt
                MergeAccordingToSpecificRunlist fileLHC17j7.txt $OUTPUTDIR_LHC17j7c $NSlashes3 GammaCalo All$addName runlists/runNumbersLHC17j7_all.txt "no"
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17j7c "/alice/cern.ch/user/a/alitrain/PWGGA/$pathTrainMC/$LHC17j7MCc/merge" All$addName $NSlashes3 "" kTRUE
        fi
    fi

    if [ $HAVELHC17n == 1 ]; then
        ls $OUTPUTDIR_LHC17n/GammaCalo-All$addName_*.root > fileLHC17n.txt
        fileNumbers=`cat fileLHC17n.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17n $NSlashes "LHC17n-pass$passNr-All$addName" "-All$addName"
        done;
    fi

    if [ $HAVELHC17j7a == 1 ]; then
        ls $OUTPUTDIR_LHC17j7a/GammaCalo-All$addName_*.root > fileLHC17j7a.txt
        fileNumbers=`cat fileLHC17j7a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17j7a $NSlashes "MC_LHC17j7-All$addName" "-All$addName"
        done;
    fi
    if [ $HAVELHC17j7b == 1 ]; then
        ls $OUTPUTDIR_LHC17j7b/GammaCalo-All$addName_*.root > fileLHC17j7b.txt
        fileNumbers=`cat fileLHC17j7b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17j7b $NSlashes "MC_LHC17j7ZDCfix-All$addName" "-All$addName"
        done;
    fi
    if [ $HAVELHC17j7c == 1 ]; then
        ls $OUTPUTDIR_LHC17j7c/GammaCalo-All$addName_*.root > fileLHC17j7c.txt
        fileNumbers=`cat fileLHC17j7c.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17j7c $NSlashes "MC_LHC17j7ZDCfixExtra-All$addName" "-All$addName"
        done;
    fi

    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaCalo_MC_LHC17j7-All$addName\_*.root > filesForMerging.txt
        echo -e "\nAll" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            rm listCurrMerge.txt
            fileA="$OUTPUTDIR/GammaCalo_MC_LHC17j7-All$addName""_$number.root"
            fileB="$OUTPUTDIR/GammaCalo_MC_LHC17j7ZDCfix-All$addName""_$number.root"
            fileC="$OUTPUTDIR/GammaCalo_MC_LHC17j7ZDCfixExtra-All$addName""_$number.root"
            echo -e "$fileA\n$fileB\n$fileC" > listCurrMerge.txt
            MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCalo_MC_LHC17j7x-All$addName\_$number.root
        done
    fi


else
    if [ $HAVELHC17n == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC17n";
        rm $OUTPUTDIR_LHC17n/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC17n/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC17j7a == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC17j7";
        rm $OUTPUTDIR_LHC17j7a/*/GammaCalo_*.root
    fi
    if [ $HAVELHC17j7b == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC17j7ZDCfix";
        rm $OUTPUTDIR_LHC17j7b/*/GammaCalo_*.root
    fi
    if [ $HAVELHC17j7c == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC17j7ZDCfixExtra";
        rm $OUTPUTDIR_LHC17j7c/*/GammaCalo_*.root
    fi
fi
