# WARNING: this script uses regular expressions (e.g. ^...$). Checkout if your terminal uses the same regex version!

# Run as:
# bash GetHeavyMesonFilesFromGridAndMergeThem_pp_13TeV.sh <your_username>

#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=0
MERGEON=1
SEPARATEON=0 # used in basisFunctions.sh?
CLEANUP=1 # used in basisFunctions.sh?
CLEANUPMAYOR=0
number=""
AliTrainDirDt="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp" # data train path on the GRID
AliTrainDirMC="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC" # MC train path on the GRID

# Print input arguments
    echo "User name: $1"
    # echo "Path: $PATH"
#

# Set number of slashes
    if [ $1 = "fbock" ]; then
        BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp13TeV
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
#

# Train configurations
    passNr="1"
    # DEFAULT LHC data files (for if you comment below input)
        LHC16Data=""
        LHC16dData=""
        LHC16gData=""
        LHC16hData=""
        LHC16iData=""
        LHC16jData=""
        LHC16kData=""
        LHC16lData=""
        LHC16oData=""
        LHC16pData=""
        LHC16eData=""
    # DEFAULT Monte Carlo data files (for if you comment below input)
        LHC17MC=""
        LHC17f6MC=""
        LHC17f9MC=""
        LHC17d17MC=""
        LHC17f5MC=""
        LHC17d3MC=""
        LHC17e5MC=""
        LHC17d20a1MC=""
        LHC17d20a2MC=""
        LHC17d16MC=""
        LHC17d18MC=""
    # LHC data files
        LHC16Data="2374" # comment: "eta' test"
        LHC16dData="child_1"
        LHC16gData="child_2"
        LHC16hData="child_3"
        LHC16iData="child_4"
        LHC16jData="child_5"
        LHC16kData="child_6"
        LHC16lData="child_7"
        LHC16oData="child_8"
        LHC16pData="child_9"
        LHC16eData="child_10"
    # Monte Carlo data files
        LHC17MC="3341" # comment: "omega + eta'"
        LHC17f6MC="child_1"
        LHC17f9MC="child_2"
        LHC17d17MC="child_3"
        LHC17f5MC="child_4"
        LHC17d3MC="child_5"
        LHC17e5MC="child_6"
        LHC17d20a1MC="child_7"
        LHC17d20a2MC="child_8"
        LHC17d16MC="child_9"
        LHC17d18MC="child_10"
#

#Set output directory
    TRAINDIR="vAN-20180123-1" # local main directory to which the files will be downloaded
    OUTPUTDIR=$BASEDIR/$TRAINDIR
#

# Set local and GRID directory names
    echo ""
    echo "--> Setting GRID and local output directories..."
    mkdir -p $OUTPUTDIR/CutSelections # this folder will contain txt files with cutnumbers
    # LHC data files
        if [ -n "$LHC16dData" ]; then # child 1
            LHC16dData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16dData$"`
            OUTPUTDIR_LHC16d=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16dData
            echo "Output dir LHC16d: $OUTPUTDIR_LHC16d"
        fi
        if [ -n "$LHC16gData" ]; then # child 2
            LHC16gData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16gData$"`
            OUTPUTDIR_LHC16g=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16gData
            echo "Output dir LHC16g: $OUTPUTDIR_LHC16g"
        fi
        if [ -n "$LHC16hData" ]; then # child 3
            LHC16hData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16hData$"`
            OUTPUTDIR_LHC16h=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16hData
            echo "Output dir LHC16h: $OUTPUTDIR_LHC16h"
        fi
        if [ -n "$LHC16iData" ]; then # child 4
            LHC16iData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16iData$"`
            OUTPUTDIR_LHC16i=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16iData
            echo "Output dir LHC16i: $OUTPUTDIR_LHC16i"
        fi
        if [ -n "$LHC16jData" ]; then # child 5
            LHC16jData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16jData$"`
            OUTPUTDIR_LHC16j=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16jData
            echo "Output dir LHC16j: $OUTPUTDIR_LHC16j"
        fi
        if [ -n "$LHC16kData" ]; then # child 6
            LHC16kData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16kData$"`
            OUTPUTDIR_LHC16k=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16kData
            echo "Output dir LHC16k: $OUTPUTDIR_LHC16k"
        fi
        if [ -n "$LHC16lData" ]; then # child 7
            LHC16lData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16lData$"`
            OUTPUTDIR_LHC16l=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16lData
            echo "Output dir LHC16l: $OUTPUTDIR_LHC16l"
        fi
        if [ -n "$LHC16oData" ]; then # child 8
            LHC16oData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16oData$"`
            OUTPUTDIR_LHC16o=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16oData
            echo "Output dir LHC16o: $OUTPUTDIR_LHC16o"
        fi
        if [ -n "$LHC16pData" ]; then # child 9
            LHC16pData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16pData$"`
            OUTPUTDIR_LHC16p=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16pData
            echo "Output dir LHC16p: $OUTPUTDIR_LHC16p"
        fi
        if [ -n "$LHC16eData" ]; then # child 10
            LHC16eData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16eData$"`
            OUTPUTDIR_LHC16e=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16eData
            echo "Output dir LHC16e: $OUTPUTDIR_LHC16e"
        fi
    # Monte Carlo files
        if [ -n "$LHC17f6MC" ]; then # child 1
            LHC17f6MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f6MC$"`
            OUTPUTDIR_LHC17f6=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f6MC
            echo "Output dir LHC17f6: $OUTPUTDIR_LHC17f6"
        fi
        if [ -n "$LHC17f9MC" ]; then # child 2
            LHC17f9MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f9MC$"`
            OUTPUTDIR_LHC17f9=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f9MC
            echo "Output dir LHC17f9: $OUTPUTDIR_LHC17f9"
        fi
        if [ -n "$LHC17d17MC" ]; then # child 3
            LHC17d17MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d17MC$"`
            OUTPUTDIR_LHC17d17=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d17MC
            echo "Output dir LHC17d17: $OUTPUTDIR_LHC17d17"
        fi
        if [ -n "$LHC17f5MC" ]; then # child 4
            LHC17f5MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f5MC$"`
            OUTPUTDIR_LHC17f5=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f5MC
            echo "Output dir LHC17f5: $OUTPUTDIR_LHC17f5"
        fi
        if [ -n "$LHC17d3MC" ]; then # child 5
            LHC17d3MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d3MC$"`
            OUTPUTDIR_LHC17d3=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d3MC
            echo "Output dir LHC17d3: $OUTPUTDIR_LHC17d3"
        fi
        if [ -n "$LHC17e5MC" ]; then # child 6
            LHC17e5MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17e5MC$"`
            OUTPUTDIR_LHC17e5=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17e5MC
            echo "Output dir LHC17e5: $OUTPUTDIR_LHC17e5"
        fi
        if [ -n "$LHC17d20a1MC" ]; then # child 7
            LHC17d20a1MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d20a1MC$"`
            OUTPUTDIR_LHC17d20a1=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a1MC
            echo "Output dir LHC17d20a1: $OUTPUTDIR_LHC17d20a1"
        fi
        if [ -n "$LHC17d20a2MC" ]; then # child 8
            LHC17d20a2MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d20a2MC$"`
            OUTPUTDIR_LHC17d20a2=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a2MC
            echo "Output dir LHC17d20a2: $OUTPUTDIR_LHC17d20a2"
        fi
        if [ -n "$LHC17d16MC" ]; then # child 9
            LHC17d16MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d16MC$"`
            OUTPUTDIR_LHC17d16=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d16MC
            echo "Output dir LHC17d16: $OUTPUTDIR_LHC17d16"
        fi
        if [ -n "$LHC17d18MC" ]; then # child 10
            LHC17d18MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d18MC$"`
            OUTPUTDIR_LHC17d18=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d18MC
            echo "Output dir LHC17d18: $OUTPUTDIR_LHC17d18"
        fi
#

# Download files from the GRID!
    if [ $DOWNLOADON == 1 ]; then
        echo ""
        echo ""
        echo "--> Downloading files from the GRID..."
    # LHC data files
        if [ -n "$LHC16dData" ]; then # child 1
            echo ""
            echo "--> Downloading LHC16d..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16dData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16dData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16gData" ]; then # child 2
            echo ""
            echo "--> Downloading LHC16g..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16gData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16gData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16hData" ]; then # child 3
            echo ""
            echo "--> Downloading LHC16h..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16hData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16hData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16iData" ]; then # child 4
            echo ""
            echo "--> Downloading LHC16i..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16iData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16iData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16jData" ]; then # child 5
            echo ""
            echo "--> Downloading LHC16j..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16jData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16jData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16kData" ]; then # child 6
            echo ""
            echo "--> Downloading LHC16k..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16kData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16kData/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16lData" ]; then # child 7
            echo ""
            echo "--> Downloading LHC16l..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16lData/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16lData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16oData" ]; then # child 8
            echo ""
            echo "--> Downloading LHC16o..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16oData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16oData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16pData" ]; then # child 9
            echo ""
            echo "--> Downloading LHC16p..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16eData" ]; then # child 10
            echo ""
            echo "--> Downloading LHC16e..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16e "$AliTrainDirDt/$LHC16eData/merge" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16e "$AliTrainDirDt/$LHC16eData/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
    # Monte Carlo files
        if [ -n "$LHC17f6MC" ]; then # child 1
            echo ""
            echo "--> Downloading LHC17f6..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f9MC" ]; then # child 2
            echo ""
            echo "--> Downloading LHC17f9..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d17MC" ]; then # child 3 (has no IncAcc)
            echo ""
            echo "--> Downloading LHC17d17..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_2" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f5MC" ]; then # child 4
            echo ""
            echo "--> Downloading LHC17f5..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d3MC" ]; then # child 5
            echo ""
            echo "--> Downloading LHC17d3..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17e5MC" ]; then # child 6
            echo ""
            echo "--> Downloading LHC17e5..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a1MC" ]; then # child 7
            echo ""
            echo "--> Downloading LHC17d20a1..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a2MC" ]; then # child 8
            echo ""
            echo "--> Downloading LHC17d20a2..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d16MC" ]; then # child 9
            echo ""
            echo "--> Downloading LHC17d16..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d18MC" ]; then # child 10
            echo ""
            echo "--> Downloading LHC17d18..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_3" DPGTracksIncAcc $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
#

# Copy files with proper name
        echo ""
        echo ""
        echo "--> Changing structure if needed..."
    # LHC data files
        if [ -n "$LHC16dData" ]; then # child 1
            for fileName in `ls $OUTPUTDIR_LHC16d/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16d/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16d/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16gData" ]; then # child 2
            for fileName in `ls $OUTPUTDIR_LHC16g/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16g/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16g/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16hData" ]; then # child 3 (has no IncAcc)
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16iData" ]; then # child 4
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16jData" ]; then # child 5
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16kData" ]; then # child 6
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16lData" ]; then # child 7
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16oData" ]; then # child 8
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16pData" ]; then # child 9
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16eData" ]; then # child 10
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
    # Monte Carlo files
        if [ -n "$LHC17f6MC" ]; then    # child 1
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root` $fileNumbers; do
            ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17f9MC" ]; then    # child 2
            echo "Testen"
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                echo fileName $fileName
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                echo fileName $fileName
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                echo fileName $fileName
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d17MC" ]; then   # child 3
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17f5MC" ]; then    # child 4
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d3MC" ]; then    # child 5
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17e5MC" ]; then    # child 6
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d20a1MC" ]; then # child 7
            echo $OUTPUTDIR_LHC17d20a1
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d20a2MC" ]; then # child 8
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d16MC" ]; then   # child 9
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d18MC" ]; then   # child 10
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
    fi # end of [ $DOWNLOADON == 1 ]

# Merge ROOT files
    if [ $MERGEON == 1 ]; then
        echo ""
        echo ""
        echo "--> Merging ROOT files"
        filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMesonToGG_LHC16d-pass$passNr-DPGTracks\_*.root` # will be used to get the number, so period doesn't matter
        # LHC data files
            echo ""
            echo "--> Merging LHC data files"
            periodList=(d e g h i j k l o p)
            rm -f runlistsToMerge.txt
            echo -e "DPGTracks" >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo "$fileName --> number: $number"
                for runListName in $listsToMerge; do
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for periodID in ${periodList[@]}; do
                        echo "Period ID: $periodID"
                        currFile=$OUTPUTDIR/HeavyNeutralMesonToGG_LHC16$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile" >> listCurrMerge.txt
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/HeavyNeutralMesonToGG_LHC16$nameOut\-pass$passNr-$runListName\_$number.root
                done
            done
        # Monte Carlo files
            echo ""
            echo "--> Merging Monte Carlo files"
            periodListMC=(f6 f9 d17 f5 d3 e5 d20a1 d20a2 d16 d18)
            rm -f runlistsToMerge.txt
            echo -e "DPGTracks" >> runlistsToMerge.txt
            echo -e "DPGTracksIncAcc" >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo ""
                echo "--> Merging $number"
                for runListName in $listsToMerge; do
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for ((i = 0; i < ${#periodListMC[@]}; i++)); do
                        echo "Period ID: 17${periodListMC[i]}-anchor16${periodList[i]}"
                        currFile=$OUTPUTDIR/HeavyNeutralMesonToGG_MC_LHC17${periodListMC[i]}-anchor16${periodList[i]}-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=${periodListMC[i]}
                            echo -e "$currFile" >> listCurrMerge.txt # maybe "$currFile\n"?
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/HeavyNeutralMesonToGG_MC_LHC17$nameOut\-$runListName\_$number.root
                done
            done
        echo "--> Merging done"
    fi # end of [ $MERGEON == 1 ]
#

# Clean ROOT files
    if [ $CLEANUPMAYOR == 1 ]; then
        echo ""
        echo "--> Cleaning up ROOT files"
        # LHC data files
            if [ -n "$LHC16dData" ]; then # child 1
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16d"
                rm $OUTPUTDIR_LHC16d/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16d/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16gData" ]; then # child 2
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16g"
                rm $OUTPUTDIR_LHC16g/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16g/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16hData" ]; then # child 3
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16h"
                rm $OUTPUTDIR_LHC16h/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16h/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16iData" ]; then # child 4
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16i"
                rm $OUTPUTDIR_LHC16i/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16i/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16jData" ]; then # child 5
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16j"
                rm $OUTPUTDIR_LHC16j/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16j/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16kData" ]; then # child 6
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16k"
                rm $OUTPUTDIR_LHC16k/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16k/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16lData" ]; then # child 7
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16l"
                rm $OUTPUTDIR_LHC16l/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16l/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16oData" ]; then # child 8
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16o"
                rm $OUTPUTDIR_LHC16o/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16o/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16pData" ]; then # child 9
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16p"
                rm $OUTPUTDIR_LHC16p/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16p/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC16eData" ]; then # child 10
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC16e"
                rm $OUTPUTDIR_LHC16e/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC16e/*/*/*HeavyNeutralMesonToGG_*.root
            fi
        # Monte Carlo files
            if [ -n "$LHC17f6MC" ]; then    # child 1
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17f6MC"
                rm $OUTPUTDIR_LHC17f6/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17f6/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17f9MC" ]; then    # child 2
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17f9MC"
                rm $OUTPUTDIR_LHC17f9/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17f9/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17d17MC" ]; then   # child 3
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d17MC"
                rm $OUTPUTDIR_LHC17d17/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17d17/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17f5MC" ]; then    # child 4
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17f5MC"
                rm $OUTPUTDIR_LHC17f5/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17f5/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17d3MC" ]; then    # child 5
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d3MC"
                rm $OUTPUTDIR_LHC17d3/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17d3/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17e5MC" ]; then    # child 6
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17e5MC"
                rm $OUTPUTDIR_LHC17e5/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17e5/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17d20a1MC" ]; then # child 7
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d20a1MC"
                rm $OUTPUTDIR_LHC17d20a1/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17d20a1/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17d20a2MC" ]; then # child 8
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d20a2MC"
                rm $OUTPUTDIR_LHC17d20a2/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17d20a2/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17d16MC" ]; then   # child 9
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d16MC"
                rm $OUTPUTDIR_LHC17d16/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17d16/*/*/*HeavyNeutralMesonToGG_*.root
            fi
            if [ -n "$LHC17d18MC" ]; then   # child 10
                echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d18MC"
                rm $OUTPUTDIR_LHC17d18/*/HeavyNeutralMesonToGG_*.root
                rm $OUTPUTDIR_LHC17d18/*/*/*HeavyNeutralMesonToGG_*.root
            fi
    fi # end of [ $CLEANUPMAYOR == 1 ]
