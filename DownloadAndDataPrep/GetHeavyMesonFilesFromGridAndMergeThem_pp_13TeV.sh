# WARNING: this script uses regular expressions (e.g. ^...$). Checkout if your terminal uses the same regex version!

#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=0
SEPARATEON=0 # used in basisFunctions.sh?
CLEANUP=1 # used in basisFunctions.sh?
CLEANUPMAYOR=$2
number=""
AliTrainDirDt="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp" # data train path on the GRID
AliTrainDirMC="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC" # MC train path on the GRID
TRAINDIR="vAN-20180123-1" # local main directory to which the files will be downloaded

# Print input arguments
    echo "User name: $1"
    echo "Clean up mayor: $2"
    # echo "Path: $PATH"
#

# Set number of slashes
    if [ $1 = "fbock" ]; then
        BASEDIR=/mnt/additionalStorageExternal/OutputLegoTrains/pp13TeV
    elif [ $1 = "hannahbossi" ]; then
        BASEDIR=/Volumes/external_memory/CERN_data/QA
    elif [ $1 = "dmuhlhei" ]; then
        BASEDIR=~/data/work/Grid
    elif [ $1 = "jlueh" ]; then
        BASEDIR=~/Daten/GridDownload
    elif [ $1 = "redeboer" ]; then
        BASEDIR=~/alice/Gridoutput
    else
        echo "NO VALID USER NAME GIVEN"
        exit 1
    fi
    NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
    NSlashes=`expr $NSlashesBASE + 4`
    NSlashes2=`expr $NSlashes - 1`
    NSlashes3=`expr $NSlashes + 1`
    NSlashes4=`expr $NSlashes + 2`
    OUTPUTDIR=$BASEDIR/$TRAINDIR
#

# Train configurations
    passNr="1";
    # DEFAULT LHC data files (for if you comment below input)
        LHC16Data="";
        LHC16dData="";
        LHC16gData="";
        LHC16hData="";
        LHC16iData="";
        LHC16jData="";
        LHC16kData="";
        LHC16lData="";
        LHC16oData="";
        LHC16pData="";
        LHC16eData="";
    # DEFAULT Monte Carlo data files (for if you comment below input)
        LHC17MC="";
        LHC17f6MC="";
        LHC17f9MC="";
        LHC17d17MC="";
        LHC17f5MC="";
        LHC17d3MC="";
        LHC17e5MC="";
        LHC17d20a1MC="";
        LHC17d20a2MC="";
        LHC17d16MC="";
        LHC17d18MC="";
    # LHC data files
        LHC16Data="2374"
        LHC16dData="child_1"
        # LHC16gData="child_2"
        # LHC16hData="child_3"
        # LHC16iData="child_4"
        # LHC16jData="child_5"
        # LHC16kData="child_6"
        # LHC16lData="child_7"
        # LHC16oData="child_8"
        # LHC16pData="child_9"
        # LHC16eData="child_10"
    # Monte Carlo data files
        LHC17MC="3225"
        # LHC17f6MC="child_1"
        # LHC17f9MC="child_2"
        # LHC17d17MC="child_3"
        # LHC17f5MC="child_4"
        # LHC17d3MC="child_5"
        # LHC17e5MC="child_6"
        # LHC17d20a1MC="child_7"
        # LHC17d20a2MC="child_8"
        # LHC17d16MC="child_9"
        # LHC17d18MC="child_10"
#

# Set local and GRID directory names
    mkdir -p $OUTPUTDIR/CutSelections # this folder will contain txt files with cutnumbers
    # LHC data files
        if [ -n "$LHC16dData" ]; then # child 1
            LHC16dData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16dData$"`
            OUTPUTDIR_LHC16d=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16dData
        fi
        if [ -n "$LHC16gData" ]; then # child 2
            LHC16gData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16gData$"`
            OUTPUTDIR_LHC16g=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16gData
        fi
        if [ -n "$LHC16hData" ]; then # child 3
            LHC16hData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16hData$"`
            OUTPUTDIR_LHC16h=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16hData
        fi
        if [ -n "$LHC16iData" ]; then # child 4
            LHC16iData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16iData$"`
            OUTPUTDIR_LHC16i=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16iData
        fi
        if [ -n "$LHC16jData" ]; then # child 5
            LHC16jData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16jData$"`
            OUTPUTDIR_LHC16j=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16jData
        fi
        if [ -n "$LHC16kData" ]; then # child 6
            LHC16kData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16kData$"`
            OUTPUTDIR_LHC16k=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16kData
        fi
        if [ -n "$LHC16lData" ]; then # child 7
            LHC16lData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16lData$"`
            OUTPUTDIR_LHC16l=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16lData
        fi
        if [ -n "$LHC16oData" ]; then # child 8
            LHC16oData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16oData$"`
            OUTPUTDIR_LHC16o=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16oData
        fi
        if [ -n "$LHC16pData" ]; then # child 9
            LHC16pData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16pData$"`
            OUTPUTDIR_LHC16p=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16pData
        fi
        if [ -n "$LHC16eData" ]; then # child 10
            LHC16eData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16eData$"`
            OUTPUTDIR_LHC16e=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16eData
        fi
    # Monte Carlo files
        if [ -n "$LHC17f6Data" ]; then # child 1
            LHC17f6Data=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC17f6Data$"`
            OUTPUTDIR_LHC17f6=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f6Data
        fi
        if [ -n "$LHC17f9Data" ]; then # child 2
            LHC17f9Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17f9Data$"`
            OUTPUTDIR_LHC17f9=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f9Data
        fi
        if [ -n "$LHC17d17Data" ]; then # child 3
            LHC17d17Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17d17Data$"`
            OUTPUTDIR_LHC17d17=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d17Data
        fi
        if [ -n "$LHC17f5Data" ]; then # child 4
            LHC17f5Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17f5Data$"`
            OUTPUTDIR_LHC17f5=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f5Data
        fi
        if [ -n "$LHC17d3Data" ]; then # child 5
            LHC17d3Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17d3Data$"`
            OUTPUTDIR_LHC17d3=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d3Data
        fi
        if [ -n "$LHC17e5Data" ]; then # child 6
            LHC17e5Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17e5Data$"`
            OUTPUTDIR_LHC17e5=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17e5Data
        fi
        if [ -n "$LHC17d20a1Data" ]; then # child 7
            LHC17d20a1Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17d20a1Data$"`
            OUTPUTDIR_LHC17d20a1=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a1Data
        fi
        if [ -n "$LHC17d20a2Data" ]; then # child 8
            LHC17d20a2Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17d20a2Data$"`
            OUTPUTDIR_LHC17d20a2=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a2Data
        fi
        if [ -n "$LHC17d16Data" ]; then # child 9
            LHC17d16Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17d16Data$"`
            OUTPUTDIR_LHC17d16=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d16Data
        fi
        if [ -n "$LHC17d18Data" ]; then # child 10
            LHC17d18Data=`alien_ls $AliTrainDirMC/ | grep "^$LHC16Data.*$LHC17d18Data$"`
            OUTPUTDIR_LHC17d18=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d18Data
        fi
#

if [ $CLEANUPMAYOR == 0 ]; then
# Download files from the GRID!
    # LHC data files
        if [ -n "$LHC16dData" ]; then # child 1
            echo ""
            echo "--> Downloading LHC16d"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16dData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16dData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16gData" ]; then # child 2
            echo ""
            echo "--> Downloading LHC16g"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16gData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16gData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16hData" ]; then # child 3
            echo ""
            echo "--> Downloading LHC16h"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16hData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16hData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16iData" ]; then # child 4
            echo ""
            echo "--> Downloading LHC16i"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16iData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16iData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16jData" ]; then # child 5
            echo ""
            echo "--> Downloading LHC16j"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16jData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16jData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16kData" ]; then # child 6
            echo ""
            echo "--> Downloading LHC16k"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16kData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16kData/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16lData" ]; then # child 7
            echo ""
            echo "--> Downloading LHC16l"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16lData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16lData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16oData" ]; then # child 8
            echo ""
            echo "--> Downloading LHC16o"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16oData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16oData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16pData" ]; then # child 9
            echo ""
            echo "--> Downloading LHC16p"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16eData" ]; then # child 10
            echo ""
            echo "--> Downloading LHC16e"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
    # Monte Carlo files
        if [ -n "$LHC17f6Data" ]; then # child 1
            echo ""
            echo "--> Downloading LHC17f6"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f9Data" ]; then # child 2
            echo ""
            echo "--> Downloading LHC17f9"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d17Data" ]; then # child 3
            echo ""
            echo "--> Downloading LHC17d17"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f5Data" ]; then # child 4
            echo ""
            echo "--> Downloading LHC17f5"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d3Data" ]; then # child 5
            echo ""
            echo "--> Downloading LHC17d3"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17e5Data" ]; then # child 6
            echo ""
            echo "--> Downloading LHC17e5"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a1Data" ]; then # child 7
            echo ""
            echo "--> Downloading LHC17d20a1"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a2Data" ]; then # child 8
            echo ""
            echo "--> Downloading LHC17d20a2"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d16Data" ]; then # child 9
            echo ""
            echo "--> Downloading LHC17d16"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d18Data" ]; then # child 10
            echo ""
            echo "--> Downloading LHC17d18"
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
#

# Copy files with proper name
echo "Change Structure If Needed"
    # LHC data files
        if [ -n "$LHC16dData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16d/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16d/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16d/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
            echo "16d done"
        fi
        if [ -n "$LHC16gData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16g/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16g/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16g/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16hData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16iData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16jData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16kData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16lData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16oData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16pData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC16eData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17f6Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-IncAcc_*.root` $fileNumbers; do
            ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
    # Monte Carlo files
        if [ -n "$LHC17f9Data" ]; then
            echo "Testen"
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                echo fileName $fileName
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                echo fileName $fileName
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                echo fileName $fileName
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d17Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17f5Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d3Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17e5Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17f6Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17f9Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d17Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17f5Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d3Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17e5Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17d3-anchor16j-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d20a1Data" ]; then
            echo $OUTPUTDIR_LHC17d20a1
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d20a1ExData" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d20a1Ex/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d20a1Ex/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d20a1Ex/HeavyNeutralMesonToGG-DPGTrackAndCalo*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d20a2Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d16Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi
        if [ -n "$LHC17d18Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi

        if [ -n "$LHC17d16Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi

        if [ -n "$LHC17d18Data" ]; then
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-IncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-IncAcc" "-IncAcc"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTrack_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrack" "-DPGTrack"
            done;
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTrackAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrackAndCalo" "-DPGTrackAndCalo"
            done;
        fi

    echo "--> All downloads done!"
    echo ""
# Merge ROOT files
    if [ $MERGEON == 1 ]; then
        echo ""
        echo "--> Merging ROOT files"
        filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-DPGTrack\_*.root`
        for fileName in $filesForMerging; do
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-DPGTrack\_$number.root
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-DPGTrack\_$number.root
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-DPGTrack\_$number.root
            if [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-DPGTrack\_$number.root ] && [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-DPGTrack\_$number.root ] && [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-DPGTrack\_$number.root ] ; then
                hadd -f $OUTPUTDIR/HeavyNeutralMeson_LHC16klo-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-DPGTrack\_$number.root
            fi
        done
        filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMeson_LHC16*-pass$passNr-IncAcc\_*.root`
        for fileName in $filesForMerging; do
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-IncAcc\_$number.root
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-IncAcc\_$number.root
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-IncAcc\_$number.root
            if [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-IncAcc\_$number.root ] && [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-IncAcc\_$number.root ] && [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-IncAcc\_$number.root ]; then
                hadd -f $OUTPUTDIR/HeavyNeutralMeson_LHC16klo-pass$passNr-IncAcc\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-IncAcc\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-IncAcc\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-IncAcc\_$number.root
            fi
        done
        filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMeson_LHC16*-pass$passNr-DPGTrackAndCalo\_*.root`
        for fileName in $filesForMerging; do
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-DPGTrackAndCalo\_$number.root
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root
            ls $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root
            if [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root ] && [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root ]  && [ -f $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-DPGTrackAndCalo\_$number.root ] ; then
                hadd -f $OUTPUTDIR/HeavyNeutralMeson_LHC16klo-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/HeavyNeutralMeson_LHC16o-pass$passNr-DPGTrackAndCalo\_$number.root
            fi
        done


    echo "--> Merging done"
    fi
else
    echo ""
    echo "--> Cleaning up ROOT files"
    # LHC data files
        if [ -n "$LHC16dData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16d";
            rm $OUTPUTDIR_LHC16d/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16d/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16gData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16g";
            rm $OUTPUTDIR_LHC16g/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16g/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16hData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16h";
            rm $OUTPUTDIR_LHC16h/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16h/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16iData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16i";
            rm $OUTPUTDIR_LHC16i/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16i/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16jData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16j";
            rm $OUTPUTDIR_LHC16j/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16j/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16kData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16k";
            rm $OUTPUTDIR_LHC16k/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16k/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16lData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16l";
            rm $OUTPUTDIR_LHC16l/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16l/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16oData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16o";
            rm $OUTPUTDIR_LHC16o/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16o/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16pData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16p";
            rm $OUTPUTDIR_LHC16p/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16p/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16eData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16e";
            rm $OUTPUTDIR_LHC16e/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16e/*/*/*HeavyNeutralMeson_*.root
        fi
    # Monte Carlo files
        if [ -n "$LHC16eData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16e";
            rm $OUTPUTDIR_LHC16e/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16e/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16eData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16e";
            rm $OUTPUTDIR_LHC16e/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16e/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16eData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16e";
            rm $OUTPUTDIR_LHC16e/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16e/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16eData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16e";
            rm $OUTPUTDIR_LHC16e/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16e/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC16eData" ]; then
            echo "removing all GammaConv files in runFolders for LHC16e";
            rm $OUTPUTDIR_LHC16e/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC16e/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC17d20a1Data" ]; then
            echo "removing all GammaConv files in runFolders for LHC17d20a1";
            rm $OUTPUTDIR_LHC17d20a1/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC17d20a1/*/*/*HeavyNeutralMeson_*.root
        fi
        if [ -n "$LHC17d20a2Data" ]; then
            echo "removing all GammaConv files in runFolders for LHC17d20a2";
            rm $OUTPUTDIR_LHC17d20a2/*/HeavyNeutralMeson_*.root
            rm $OUTPUTDIR_LHC17d20a2/*/*/*HeavyNeutralMeson_*.root
        fi
fi