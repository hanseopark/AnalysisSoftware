#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

source basicFunction.sh

DOWNLOADON=0
REMERGE=0
MERGEON=1
MERGEONData=1
MERGEONMC=0
SINGLERUN=0
SEPARATEON=0
MERGEONSINGLEData=0
MERGEONSINGLEMC=0
CLEANUP=1
CLEANUPMAYOR=$2
number=""


echo $1
echo $2
# echo $PATH

# check if train configuration has actually been given
HAVELHC16f=1
HAVELHC17g=1
HAVELHC18cwoSDD=1
HAVELHC18cfast=1

HAVELHC17d1=1;
HAVELHC17h3=1;
HAVELHC18h1woSDD=1;
HAVELHC18h1fast=1;

# default trainconfigurations
LHC16fData="";
LHC17gData="";
LHC18Data=""
    LHC18cwoSDDData="";
    LHC18cfastData="";

LHC17d1MC="";
LHC17h3MC="";

LHC18xMCPHY="";
    LHC18h1woSDDMC=""
    LHC18h1fastMC=""

passNr="1";
if [ $1 = "fbock" ]; then
    BASEDIR=/media/fbock/Elements/OutputLegoTrains/pp13TeV
elif [ $1 = "hannahbossi" ]; then
    BASEDIR=/Volumes/external_memory/CERN_data/QA
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
elif [ $1 = "jlueh" ]; then
    BASEDIR=~/Daten/GridDownload
fi

if [ $3 = "AOD" ]; then
    baseLegoData=GA_pp_AOD
    baseLegoMC=GA_pp_MC_AOD
    pathData=pass$passNr/AOD212/PWGGA/GA_pp_AOD
    pathData2=pass$passNr\_CENT\_woSDD/AOD212/PWGGA/GA_pp_AOD
    pathData3=pass$passNr\_FAST/AOD212/PWGGA/GA_pp_AOD
    pathMC=AOD213/PWGGA/GA_pp_MC_AOD
    pathMC=AOD213/PWGGA/GA_pp_MC_AOD
    pathMC=AOD213/PWGGA/GA_pp_MC_AOD
else
    baseLegoData=GA_pp
    baseLegoMC=GA_pp_MC
    pathData=pass$passNr/PWGGA/GA_pp
    pathMC=PWGGA/GA_pp_MC
fi

# Definitition of number of slashes in your path to different depths
NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
NSlashes=`expr $NSlashesBASE + 4`
NSlashes2=`expr $NSlashes - 1`
NSlashes3=`expr $NSlashes + 1`
NSlashes4=`expr $NSlashes + 2`
echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"

TRAINDIR=Legotrain-QA2019LowB

LHC16fData="755";
LHC17gData="756";
LHC18Data="754"
    LHC18cwoSDDData="child_1";
    LHC18cfastData="child_2";

LHC17d1MC="1416";
LHC17h3MC="1417";
LHC18xMCPHY="1418";
    LHC18h1woSDDMC="child_1"
    LHC18h1fastMC="child_2"

OUTPUTDIR=$BASEDIR/$TRAINDIR

ALIENDIRData="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoData/"
OUTPUTDIRData=$BASEDIR/$TRAINDIR/$baseLegoData
ALIENDIRMC="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoMC/"
OUTPUTDIRMC=$BASEDIR/$TRAINDIR/$baseLegoMC
mkdir -p $OUTPUTDIR/CutSelections

FindCorrectTrainDirectory $LHC16fData $OUTPUTDIRData $ALIENDIRData
HAVELHC16f=$tempBool
LHC16fData=$tempDir
OUTPUTDIR_LHC16f=$tempPath
echo "16f: $HAVELHC16f $LHC16fData $OUTPUTDIR_LHC16f"

FindCorrectTrainDirectory $LHC17gData $OUTPUTDIRData $ALIENDIRData
HAVELHC17g=$tempBool
LHC17gData=$tempDir
OUTPUTDIR_LHC17g=$tempPath
echo "17g: $HAVELHC17g $LHC17gData $OUTPUTDIR_LHC17g"

FindCorrectTrainDirectory $LHC18cwoSDDData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18cwoSDD=$tempBool
LHC18cwoSDDData=$tempDir
OUTPUTDIR_LHC18cwoSDD=$tempPath
echo "18c_woSDD: $HAVELHC18cwoSDD $LHC18cwoSDDData $OUTPUTDIR_LHC18cwoSDD"

FindCorrectTrainDirectory $LHC18cfastData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18cfast=$tempBool
LHC18cfastData=$tempDir
OUTPUTDIR_LHC18cfast=$tempPath
echo "18c_fast: $HAVELHC18cfast $LHC18cfastData $OUTPUTDIR_LHC18cfast"


# start with finding MC directories
FindCorrectTrainDirectory $LHC17d1MC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC17d1=$tempBool
LHC17d1MC=$tempDir
OUTPUTDIR_LHC17d1=$tempPath
echo "17d1 anchored to 16f: $HAVELHC17d1 $LHC17d1MC $OUTPUTDIR_LHC17d1"

FindCorrectTrainDirectory $LHC17h3MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17h3=$tempBool
LHC17h3MC=$tempDir
OUTPUTDIR_LHC17h3=$tempPath
echo "17h3 anchored to 17g: $HAVELHC17h3 $LHC17h3MC $OUTPUTDIR_LHC17h3"

FindCorrectTrainDirectory $LHC18h1woSDDMC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18h1woSDD=$tempBool
LHC18h1woSDDMC=$tempDir
OUTPUTDIR_LHC18h1woSDD=$tempPath
echo "18h1_woSDD anchored to 18c_woSDD: $HAVELHC18h1woSDD $LHC18h1woSDDMC $OUTPUTDIR_LHC18h1woSDD"

FindCorrectTrainDirectory $LHC18h1fastMC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18h1fast=$tempBool
LHC18h1fastMC=$tempDir
OUTPUTDIR_LHC18h1fast=$tempPath
echo "18h1_fast anchored to 18c_fast: $HAVELHC18h1fast $LHC18h1fastMC $OUTPUTDIR_LHC18h1fast"


# exit


if [ $CLEANUPMAYOR == 0 ]; then
    if [ $REMERGE == 1 ]; then
        echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMergeLowB.txt
        CopyRunwiseAndMergeAccordingToRunlistData "LHC16f" $HAVELHC16f $OUTPUTDIR_LHC16f $LHC16fData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMergeLowB.txt "pass1_lowB" GammaCalo
        CopyRunwiseAndMergeAccordingToRunlistData "LHC17g" $HAVELHC17g $OUTPUTDIR_LHC17g $LHC17gData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMergeLowB.txt "pass1" GammaCalo
        CopyRunwiseAndMergeAccordingToRunlistData "LHC18c" $HAVELHC18cwoSDD $OUTPUTDIR_LHC18cwoSDD $LHC18cwoSDDData $pathData2 $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMergeLowB.txt "pass1_CENT_woSDD" GammaCalo
        CopyRunwiseAndMergeAccordingToRunlistData "LHC18c" $HAVELHC18cfast $OUTPUTDIR_LHC18cfast $LHC18cfastData $pathData3 $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMergeLowB.txt "pass1_FAST" GammaCalo

        currentDir=$PWD
        echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMergeLowB.txt
        CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d1" $HAVELHC17d1 $OUTPUTDIR_LHC17d1 $LHC17d1MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMergeLowB.txt GammaCalo "_lowB"
        CopyRunwiseAndMergeAccordingToRunlistMC "LHC17h3" $HAVELHC17h3 $OUTPUTDIR_LHC17h3 $LHC17h3MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMergeLowB.txt GammaCalo
        CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h1_cent_woSDD" $HAVELHC18h1woSDD $OUTPUTDIR_LHC18h1woSDD $LHC18h1woSDDMC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMergeLowB.txt GammaCalo
        CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h1_fast" $HAVELHC18h1fast $OUTPUTDIR_LHC18h1fast $LHC18h1fastMC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMergeLowB.txt GammaCalo
    fi
    echo "Change Structure If Needed"

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMergeLowB.txt
    listsToMerge=`cat runlistsToMergeLowB.txt`
    for runListName in $listsToMerge; do
        if [ $HAVELHC16f == 1 ]; then
            ls $OUTPUTDIR_LHC16f/GammaCalo-$runListName\_*.root > fileLHC16f.txt
            fileNumbers=`cat fileLHC16f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16f $NSlashes "LHC16f-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17g == 1 ]; then
            ls $OUTPUTDIR_LHC17g/GammaCalo-$runListName\_*.root > fileLHC17g.txt
            fileNumbers=`cat fileLHC17g.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17g $NSlashes "LHC17g-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18cwoSDD == 1 ]; then
            ls $OUTPUTDIR_LHC18cwoSDD/GammaCalo-$runListName\_*.root > fileLHC18cwoSDD.txt
            fileNumbers=`cat fileLHC18cwoSDD.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18cwoSDD $NSlashes "LHC18c_woSDD-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18cfast == 1 ]; then
            ls $OUTPUTDIR_LHC18cfast/GammaCalo-$runListName\_*.root > fileLHC18cfast.txt
            fileNumbers=`cat fileLHC18cfast.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18cfast $NSlashes "LHC18c_fast-pass$passNr-$runListName" "-$runListName"
            done;
        fi
    done

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMergeLowB.txt
#         echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMergeLowB.txt
    listsToMerge=`cat runlistsToMergeLowB.txt`
    for runListName in $listsToMerge; do
        # MC for LHC16f
        if [ $HAVELHC17d1 == 1 ]; then
            ls $OUTPUTDIR_LHC17d1/GammaCalo-$runListName\_*.root > fileLHC17d1.txt
            fileNumbers=`cat fileLHC17d1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d1 $NSlashes "MC_LHC17d1-anchor16f-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC17g
        if [ $HAVELHC17h3 == 1 ]; then
            ls $OUTPUTDIR_LHC17h3/GammaCalo-$runListName\_*.root > fileLHC17h3.txt
            fileNumbers=`cat fileLHC17h3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17h3 $NSlashes "MC_LHC17h3-anchor17g-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC18cwoSDD
        if [ $HAVELHC18h1woSDD == 1 ]; then
            ls $OUTPUTDIR_LHC18h1woSDD/GammaCalo-$runListName\_*.root > fileLHC18h1woSDD.txt
            fileNumbers=`cat fileLHC18h1woSDD.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18h1woSDD $NSlashes "MC_LHC18h1_woSDD-anchor18c-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC18cfast
        if [ $HAVELHC18h1fast == 1 ]; then
            ls $OUTPUTDIR_LHC18h1fast/GammaCalo-$runListName\_*.root > fileLHC18h1fast.txt
            fileNumbers=`cat fileLHC18h1fast.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18h1fast $NSlashes "MC_LHC18h1_fast-anchor18c-$runListName" "-$runListName"
            done;
        fi
    done

    echo "Download Done"

    if [ $MERGEON == 1 ]; then
        echo "Starting Merging"

        if [ $MERGEONData == 1 ]; then
            echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMergeLowB.txt
            listsToMerge=`cat runlistsToMergeLowB.txt`
            for runListName in $listsToMerge; do
                ls $OUTPUTDIR/GammaCalo_LHC17g-pass$passNr-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                periodList=`echo -e "16f-pass1\n17g-pass1\n18c_woSDD-pass1\n18c_fast-pass1"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaCalo_LHC$periodID-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            outAdd=`echo $periodID  | cut -d "-" -f 1 `
                            nameOut+=$outAdd"-"
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        else
                            echo $currFile " does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCalo_LHC$nameOut\pass$passNr-$runListName\_$number.root
                done
                periodList=`echo -e "16f-pass1\n17g-pass1"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaCalo_LHC$periodID-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            outAdd=`echo $periodID  | cut -d "-" -f 1 `
                            nameOut+=$outAdd"-"
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        else
                            echo $currFile " does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCalo_LHC$nameOut\pass$passNr-$runListName\_$number.root
                done
                periodList=`echo -e "18c_woSDD-pass1\n18c_fast-pass1"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaCalo_LHC$periodID-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            outAdd=`echo $periodID  | cut -d "-" -f 1 `
                            nameOut+=$outAdd"-"
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        else
                            echo $currFile " does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCalo_LHC$nameOut\pass$passNr-$runListName\_$number.root
                done
            done
        fi
        if [ $MERGEONMC == 1 ]; then
            echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMergeLowB.txt
            listsToMerge=`cat runlistsToMergeLowB.txt`
            for runListName in $listsToMerge; do
                ls $OUTPUTDIR/GammaCalo_MC_LHC17h3-anchor17g-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                periodList=`echo -e "17d1-anchor16f\n17h3-anchor17g\n18h1_woSDD-anchor18c\n18h1_fast-anchor18c"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 4 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaCalo_MC_LHC$periodID-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            outAdd=`echo $periodID  | cut -d "-" -f 1 `
                            nameOut+=$outAdd"-"
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        else
                            echo $currFile " does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCalo_MC_LHC$nameOut$runListName\_$number.root
                done
          done
        fi

        echo "Merging Done"
    fi
else
    if [ $HAVELHC16f == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16f";
        rm $OUTPUTDIR_LHC16f/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16f/*/*/*GammaCalo_*.root
    fi
fi
