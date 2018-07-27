# WARNING: this script uses regular expressions (e.g. ^...$). Checkout if your terminal uses the same regex version!

# Run as:
# bash GetGammaConvFilesFromGridAndMergeThem_pp_13TeV.sh <your_username>

#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
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
        BASEDIR=~/alice/GridOutput
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
        LHC16=""
        LHC16d=""
        LHC16g=""
        LHC16h=""
        LHC16i=""
        LHC16j=""
        LHC16k=""
        LHC16l=""
        LHC16o=""
        LHC16p=""
        LHC16e=""
    # DEFAULT Monte Carlo data files (for if you comment below input)
        LHC17MC=""
        LHC17MCem=""
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
    # DEFAULT Monte Carlo data files EMcal (for if you comment below input)
        LHC17MCem=""
        LHC17f6MCem=""
        LHC17f9MCem=""
        LHC17d17MCem=""
        LHC17f5MCem=""
        LHC17d3MCem=""
        LHC17e5MCem=""
        LHC17d20a1MCem=""
        LHC17d20a2MCem=""
        LHC17d16MCem=""
        LHC17d18MCem=""
    # LHC data files: LHC16_13TeV_pass1
        LHC16="2367" # comment: "Hannah, Joshua, Ana, Jens..." date: 2018 Jun 04
        LHC16d="child_1"
        LHC16g="child_2"
        # LHC16h="child_3"
        LHC16i="child_4"
        LHC16j="child_5"
        LHC16k="child_6"
        LHC16l="child_7"
        LHC16o="child_8"
        LHC16p="child_9"
        LHC16e="child_10"
    # Monte Carlo data files: LHC17_PYT8_13TeV_anchLHC16 PCM
        LHC17MC="3320" # comment: "Request Jens, Joshua" date: 2018 Jun 04 (Calo+ConvV1)
        LHC17f6MC="child_1"    # anchored to LHC16d
        LHC17f9MC="child_2"    # anchored to LHC16e
        LHC17d17MC="child_3"   # anchored to LHC16g
        # LHC17f5MC="child_4"    # anchored to LHC16h
        LHC17d3MC="child_5"    # anchored to LHC16i
        LHC17e5MC="child_6"    # anchored to LHC16j
        LHC17d20a1MC="child_7" # anchored to LHC16k
        LHC17d20a2MC="child_8" # anchored to LHC16l
        LHC17d16MC="child_9"   # anchored to LHC16o
        LHC17d18MC="child_10"  # anchored to LHC16p
    # Monte Carlo data files: LHC17_PYT8_13TeV_anchLHC16 EMcal only (ConvCalo)
        LHC17MCem="3321" # comment: "Hannah request" date: 2018 Jun 04
        LHC17f6MCem="child_1"    # anchored to LHC16d
        LHC17f9MCem="child_2"    # anchored to LHC16e
        LHC17d17MCem="child_3"   # anchored to LHC16g
        # LHC17f5MCem="child_4"    # anchored to LHC16h
        LHC17d3MCem="child_5"    # anchored to LHC16i
        LHC17e5MCem="child_6"    # anchored to LHC16j
        LHC17d20a1MCem="child_7" # anchored to LHC16k
        LHC17d20a2MCem="child_8" # anchored to LHC16l
        LHC17d16MCem="child_9"   # anchored to LHC16o
        LHC17d18MCem="child_10"  # anchored to LHC16p
#

#Set output directory
    TRAINDIR="vAN-20180603-1" # local main directory to which the files will be downloaded
    OUTPUTDIR=$BASEDIR/$TRAINDIR
#

# Set local and GRID directory names
    echo -e "\n\n--== SETTING GRID DIRECTORY NAMES AND OUTPUT DIRECTORIES ==--\n"
    mkdir -p $OUTPUTDIR/CutSelections # this folder will contain txt files with cutnumbers
    # LHC data files
        if [ -n "$LHC16" ]; then
            echo -e "___________________"
            echo -e "LHC 2016 data files"
        fi
        if [ -n "$LHC16d" ]; then # child 1
            LHC16d=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16d$"`
            OUTPUTDIR_LHC16d=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16d
            echo "Directories for LHC16d:"
            echo "  local: \"$OUTPUTDIR_LHC16d\""
            echo "  alien: \"$AliTrainDirDt/$LHC16d\""
        fi
        if [ -n "$LHC16g" ]; then # child 2
            LHC16g=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16g$"`
            OUTPUTDIR_LHC16g=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16g
            echo "Directories for LHC16g:"
            echo "  local: \"$OUTPUTDIR_LHC16g\""
            echo "  alien: \"$AliTrainDirDt/$LHC16g\""
        fi
        if [ -n "$LHC16h" ]; then # child 2
            LHC16h=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16h$"`
            OUTPUTDIR_LHC16h=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16h
            echo "Directories for LHC16h:"
            echo "  local: \"$OUTPUTDIR_LHC16h\""
            echo "  alien: \"$AliTrainDirDt/$LHC16h\""
        fi
        if [ -n "$LHC16i" ]; then # child 4
            LHC16i=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16i$"`
            OUTPUTDIR_LHC16i=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16i
            echo "Directories for LHC16i:"
            echo "  local: \"$OUTPUTDIR_LHC16i\""
            echo "  alien: \"$AliTrainDirDt/$LHC16i\""
        fi
        if [ -n "$LHC16j" ]; then # child 5
            LHC16j=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16j$"`
            OUTPUTDIR_LHC16j=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16j
            echo "Directories for LHC16j:"
            echo "  local: \"$OUTPUTDIR_LHC16j\""
            echo "  alien: \"$AliTrainDirDt/$LHC16j\""
        fi
        if [ -n "$LHC16k" ]; then # child 6
            LHC16k=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16k$"`
            OUTPUTDIR_LHC16k=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16k
            echo "Directories for LHC16k:"
            echo "  local: \"$OUTPUTDIR_LHC16k\""
            echo "  alien: \"$AliTrainDirDt/$LHC16k\""
        fi
        if [ -n "$LHC16l" ]; then # child 7
            LHC16l=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16l$"`
            OUTPUTDIR_LHC16l=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16l
            echo "Directories for LHC16l:"
            echo "  local: \"$OUTPUTDIR_LHC16l\""
            echo "  alien: \"$AliTrainDirDt/$LHC16l\""
        fi
        if [ -n "$LHC16o" ]; then # child 8
            LHC16o=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16o$"`
            OUTPUTDIR_LHC16o=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16o
            echo "Directories for LHC16o:"
            echo "  local: \"$OUTPUTDIR_LHC16o\""
            echo "  alien: \"$AliTrainDirDt/$LHC16o\""
        fi
        if [ -n "$LHC16p" ]; then # child 9
            LHC16p=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16p$"`
            OUTPUTDIR_LHC16p=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16p
            echo "Directories for LHC16p:"
            echo "  local: \"$OUTPUTDIR_LHC16p\""
            echo "  alien: \"$AliTrainDirDt/$LHC16p\""
        fi
        if [ -n "$LHC16e" ]; then # child 10
            LHC16e=`alien_ls $AliTrainDirDt/ | grep "^$LHC16.*$LHC16e$"`
            OUTPUTDIR_LHC16e=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16e
            echo "Directories for LHC16e:"
            echo "  local: \"$OUTPUTDIR_LHC16e\""
            echo "  alien: \"$AliTrainDirDt/$LHC16e\""
        fi
    # Monte Carlo files 3320 (Calo+ConvV1)
        if [ -n "$LHC17MC" ]; then
            echo -e "\n______________________"
            echo -e   "Monte Carlo files $LHC17MC"
        fi
        if [ -n "$LHC17f6MC" ]; then # child 1
            LHC17f6MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f6MC$"`
            OUTPUTDIR_LHC17f6="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f6MC"
            echo "Directories for LHC17f6:"
            echo "  local: \"$OUTPUTDIR_LHC17f6\""
            echo "  alien: \"$AliTrainDirMC/$LHC17f6MC\""
        fi
        if [ -n "$LHC17f9MC" ]; then # child 2
            LHC17f9MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f9MC$"`
            OUTPUTDIR_LHC17f9="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f9MC"
            echo "Directories for LHC17f9:"
            echo "  local: \"$OUTPUTDIR_LHC17f9\""
            echo "  alien: \"$AliTrainDirMC/$LHC17f9MC\""
        fi
        if [ -n "$LHC17d17MC" ]; then # child 3
            LHC17d17MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d17MC$"`
            OUTPUTDIR_LHC17d17="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d17MC"
            echo "Directories for LHC17d17:"
            echo "  local: \"$OUTPUTDIR_LHC17d17\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d17MC\""
        fi
        if [ -n "$LHC17f5MC" ]; then # child 4
            LHC17f5MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f5MC$"`
            OUTPUTDIR_LHC17f5="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f5MC"
            echo "Directories for LHC17f5:"
            echo "  local: \"$OUTPUTDIR_LHC17f5\""
            echo "  alien: \"$AliTrainDirMC/$LHC17f5MC\""
        fi
        if [ -n "$LHC17d3MC" ]; then # child 5
            LHC17d3MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d3MC$"`
            OUTPUTDIR_LHC17d3="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d3MC"
            echo "Directories for LHC17d3:"
            echo "  local: \"$OUTPUTDIR_LHC17d3\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d3MC\""
        fi
        if [ -n "$LHC17e5MC" ]; then # child 6
            LHC17e5MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17e5MC$"`
            OUTPUTDIR_LHC17e5="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17e5MC"
            echo "Directories for LHC17e5:"
            echo "  local: \"$OUTPUTDIR_LHC17e5\""
            echo "  alien: \"$AliTrainDirMC/$LHC17e5MC\""
        fi
        if [ -n "$LHC17d20a1MC" ]; then # child 7
            LHC17d20a1MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d20a1MC$"`
            OUTPUTDIR_LHC17d20a1="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a1MC"
            echo "Directories for LHC17d20a1:"
            echo "  local: \"$OUTPUTDIR_LHC17d20a1\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d20a1MC\""
        fi
        if [ -n "$LHC17d20a2MC" ]; then # child 8
            LHC17d20a2MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d20a2MC$"`
            OUTPUTDIR_LHC17d20a2="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a2MC"
            echo "Directories for LHC17d20a2:"
            echo "  local: \"$OUTPUTDIR_LHC17d20a2\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d20a2MC\""
        fi
        if [ -n "$LHC17d16MC" ]; then # child 9
            LHC17d16MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d16MC$"`
            OUTPUTDIR_LHC17d16="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d16MC"
            echo "Directories for LHC17d16:"
            echo "  local: \"$OUTPUTDIR_LHC17d16\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d16MC\""
        fi
        if [ -n "$LHC17d18MC" ]; then # child 10
            LHC17d18MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d18MC$"`
            OUTPUTDIR_LHC17d18="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d18MC"
            echo "Directories for LHC17d18:"
            echo "  local: \"$OUTPUTDIR_LHC17d18\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d18MC\""
        fi
    # Monte Carlo files 3321 (ConvCalo)
        if [ -n "$LHC17MC" ]; then
            echo -e "\n______________________"
            echo -e   "Monte Carlo files $LHC17MCem"
        fi
        if [ -n "$LHC17f6MCem" ]; then # child 1
            LHC17f6MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17f6MCem$"`
            OUTPUTDIR_LHC17f6em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f6MCem"
            echo "Directories for LHC17f6em:"
            echo "  local: \"$OUTPUTDIR_LHC17f6em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17f6MCem\""
        fi
        if [ -n "$LHC17f9MCem" ]; then # child 2
            LHC17f9MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17f9MCem$"`
            OUTPUTDIR_LHC17f9em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f9MCem"
            echo "Directories for LHC17f9em:"
            echo "  local: \"$OUTPUTDIR_LHC17f9em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17f9MCem\""
        fi
        if [ -n "$LHC17d17MCem" ]; then # child 3
            LHC17d17MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17d17MCem$"`
            OUTPUTDIR_LHC17d17em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d17MCem"
            echo "Directories for LHC17d17em:"
            echo "  local: \"$OUTPUTDIR_LHC17d17em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d17MCem\""
        fi
        if [ -n "$LHC17f5MCem" ]; then # child 4
            LHC17f5MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17f5MCem$"`
            OUTPUTDIR_LHC17f5em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f5MCem"
            echo "Directories for LHC17f5em:"
            echo "  local: \"$OUTPUTDIR_LHC17f5em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17f5MCem\""
        fi
        if [ -n "$LHC17d3MCem" ]; then # child 5
            LHC17d3MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17d3MCem$"`
            OUTPUTDIR_LHC17d3em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d3MCem"
            echo "Directories for LHC17d3em:"
            echo "  local: \"$OUTPUTDIR_LHC17d3em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d3MCem\""
        fi
        if [ -n "$LHC17e5MCem" ]; then # child 6
            LHC17e5MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17e5MCem$"`
            OUTPUTDIR_LHC17e5em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17e5MCem"
            echo "Directories for LHC17e5em:"
            echo "  local: \"$OUTPUTDIR_LHC17e5em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17e5MCem\""
        fi
        if [ -n "$LHC17d20a1MCem" ]; then # child 7
            LHC17d20a1MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17d20a1MCem$"`
            OUTPUTDIR_LHC17d20a1em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a1MCem"
            echo "Directories for LHC17d20a1em:"
            echo "  local: \"$OUTPUTDIR_LHC17d20a1em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d20a1MCem\""
        fi
        if [ -n "$LHC17d20a2MCem" ]; then # child 8
            LHC17d20a2MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17d20a2MCem$"`
            OUTPUTDIR_LHC17d20a2em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a2MCem"
            echo "Directories for LHC17d20a2em:"
            echo "  local: \"$OUTPUTDIR_LHC17d20a2em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d20a2MCem\""
        fi
        if [ -n "$LHC17d16MCem" ]; then # child 9
            LHC17d16MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17d16MCem$"`
            OUTPUTDIR_LHC17d16em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d16MCem"
            echo "Directories for LHC17d16em:"
            echo "  local: \"$OUTPUTDIR_LHC17d16em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d16MCem\""
        fi
        if [ -n "$LHC17d18MCem" ]; then # child 10
            LHC17d18MCem=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCem.*$LHC17d18MCem$"`
            OUTPUTDIR_LHC17d18em="$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d18MCem"
            echo "Directories for LHC17d18em:"
            echo "  local: \"$OUTPUTDIR_LHC17d18em\""
            echo "  alien: \"$AliTrainDirMC/$LHC17d18MCem\""
        fi
#

# Download files from the GRID!
    if [ $DOWNLOADON == 1 ]; then
        echo -e "\n\n--== DOWNLOADING FILES FROM THE GRID ==--\n"
        # LHC data files
            if [ -n "$LHC16d" ]; then # child 1
                echo -e "\n\n--> Downloading LHC16d ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16d/merge" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16d/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16g" ]; then # child 2
                echo -e "\n\n--> Downloading LHC16g ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16g/merge" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16g/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16h" ]; then # child 2
                echo -e "\n\n--> Downloading LHC16h ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16h/merge" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16h/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16i" ]; then # child 4
                echo -e "\n\n--> Downloading LHC16i ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16i/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16i/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16j" ]; then # child 5
                echo -e "\n\n--> Downloading LHC16j ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16j/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16j/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16k" ]; then # child 6
                echo -e "\n\n--> Downloading LHC16k ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16k/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16k/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16l" ]; then # child 7
                echo -e "\n\n--> Downloading LHC16l ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16l/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16l/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16o" ]; then # child 8
                echo -e "\n\n--> Downloading LHC16o ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16o/merge" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16o/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16p" ]; then # child 9
                echo -e "\n\n--> Downloading LHC16p ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16p/merge" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16p/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC16e" ]; then # child 10
                echo -e "\n\n--> Downloading LHC16e ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16e "$AliTrainDirDt/$LHC16e/merge" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16e "$AliTrainDirDt/$LHC16e/merge" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
        # Monte Carlo files 3320 (Calo+ConvV1)
            if [ -n "$LHC17f6MC" ]; then # child 1
                echo -e "\n\n--> Downloading LHC17f6 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17f9MC" ]; then # child 2
                echo -e "\n\n--> Downloading LHC17f9 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d17MC" ]; then # child 3 (has no IncAcc)
                echo -e "\n\n--> Downloading LHC17d17 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17f5MC" ]; then # child 4
                echo -e "\n\n--> Downloading LHC17f5 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d3MC" ]; then # child 5
                echo -e "\n\n--> Downloading LHC17d3 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17e5MC" ]; then # child 6
                echo -e "\n\n--> Downloading LHC17e5 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d20a1MC" ]; then # child 7
                echo -e "\n\n--> Downloading LHC17d20a1 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d20a2MC" ]; then # child 8
                echo -e "\n\n--> Downloading LHC17d20a2 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d16MC" ]; then # child 9
                echo -e "\n\n--> Downloading LHC17d16 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d18MC" ]; then # child 10
                echo -e "\n\n--> Downloading LHC17d18 ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
        # Monte Carlo files 3321 (ConvCalo)
            if [ -n "$LHC17f6MCem" ]; then # child 1
                echo -e "\n\n--> Downloading LHC17f6em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6em "$AliTrainDirMC/$LHC17f6MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6em "$AliTrainDirMC/$LHC17f6MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17f9MCem" ]; then # child 2
                echo -e "\n\n--> Downloading LHC17f9em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9em "$AliTrainDirMC/$LHC17f9MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9em "$AliTrainDirMC/$LHC17f9MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d17MCem" ]; then # child 3 (has no IncAcc)
                echo -e "\n\n--> Downloading LHC17d17em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17em "$AliTrainDirMC/$LHC17d17MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17em "$AliTrainDirMC/$LHC17d17MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17f5MCem" ]; then # child 4
                echo -e "\n\n--> Downloading LHC17f5em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5em "$AliTrainDirMC/$LHC17f5MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5em "$AliTrainDirMC/$LHC17f5MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d3MCem" ]; then # child 5
                echo -e "\n\n--> Downloading LHC17d3em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3em "$AliTrainDirMC/$LHC17d3MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3em "$AliTrainDirMC/$LHC17d3MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17e5MCem" ]; then # child 6
                echo -e "\n\n--> Downloading LHC17e5em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5em "$AliTrainDirMC/$LHC17e5MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5em "$AliTrainDirMC/$LHC17e5MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d20a1MCem" ]; then # child 7
                echo -e "\n\n--> Downloading LHC17d20a1em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1em "$AliTrainDirMC/$LHC17d20a1MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1em "$AliTrainDirMC/$LHC17d20a1MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d20a2MCem" ]; then # child 8
                echo -e "\n\n--> Downloading LHC17d20a2em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2em "$AliTrainDirMC/$LHC17d20a2MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2em "$AliTrainDirMC/$LHC17d20a2MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d16MCem" ]; then # child 9
                echo -e "\n\n--> Downloading LHC17d16em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16em "$AliTrainDirMC/$LHC17d16MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16em "$AliTrainDirMC/$LHC17d16MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
            if [ -n "$LHC17d18MCem" ]; then # child 10
                echo -e "\n\n--> Downloading LHC17d18em ..."
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18em "$AliTrainDirMC/$LHC17d18MCem/merge_runlist_2" DPGTracks $NSlashes3 "" kFALSE
                CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18em "$AliTrainDirMC/$LHC17d18MCem/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            fi
#

# Copy files with proper name
    echo -e "\n\n--== CHANGING STRUCTURE IF NEEDED ==--\n"
    # LHC data files
        if [ -n "$LHC16d" ]; then # child 1
            for fileName in `ls $OUTPUTDIR_LHC16d/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16d/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16d/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16d/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16d/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16d/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16g" ]; then # child 2
            for fileName in `ls $OUTPUTDIR_LHC16g/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16g/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16g/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16g/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16g/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16g/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16h" ]; then # child 3 (has no IncAcc)
            for fileName in `ls $OUTPUTDIR_LHC16h/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16i" ]; then # child 4
            for fileName in `ls $OUTPUTDIR_LHC16i/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16j" ]; then # child 5
            for fileName in `ls $OUTPUTDIR_LHC16j/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16k" ]; then # child 6
            for fileName in `ls $OUTPUTDIR_LHC16k/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16l" ]; then # child 7
            for fileName in `ls $OUTPUTDIR_LHC16l/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16o" ]; then # child 8
            for fileName in `ls $OUTPUTDIR_LHC16o/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16p" ]; then # child 9
            for fileName in `ls $OUTPUTDIR_LHC16p/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC16e" ]; then # child 10
            for fileName in `ls $OUTPUTDIR_LHC16e/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
    # Monte Carlo files 3320 (Calo+ConvV1)
        if [ -n "$LHC17f6MC" ]; then    # child 1
            for fileName in `ls $OUTPUTDIR_LHC17f6/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17f9MC" ]; then    # child 2
            for fileName in `ls $OUTPUTDIR_LHC17f9/GammaConvV1-DPGTracks_*.root`;do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/GammaConvV1-DPGTracksAndCalo_*.root`;do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/GammaCalo-DPGTracks_*.root`;do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/GammaCalo-DPGTracksAndCalo_*.root`;do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d17MC" ]; then   # child 3
            for fileName in `ls $OUTPUTDIR_LHC17d17/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17f5MC" ]; then    # child 4
            for fileName in `ls $OUTPUTDIR_LHC17f5/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d3MC" ]; then    # child 5
            for fileName in `ls $OUTPUTDIR_LHC17d3/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17e5MC" ]; then    # child 6
            for fileName in `ls $OUTPUTDIR_LHC17e5/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d20a1MC" ]; then # child 7
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d20a2MC" ]; then # child 8
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d16MC" ]; then   # child 9
            for fileName in `ls $OUTPUTDIR_LHC17d16/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d18MC" ]; then   # child 10
            for fileName in `ls $OUTPUTDIR_LHC17d18/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/GammaConvV1-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/GammaCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/GammaCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
    # Monte Carlo files 3321 (ConvCalo)
        if [ -n "$LHC17f6MCem" ]; then    # child 1
            for fileName in `ls $OUTPUTDIR_LHC17f6em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6em $NSlashes "MC_LHC17f6-anchor16d-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6em $NSlashes "MC_LHC17f6-anchor16d-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17f9MCem" ]; then    # child 2
            for fileName in `ls $OUTPUTDIR_LHC17f9em/GammaConvCalo-DPGTracks_*.root`;do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9em $NSlashes "MC_LHC17f9-anchor16e-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9em/GammaConvCalo-DPGTracksAndCalo_*.root`;do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9em $NSlashes "MC_LHC17f9-anchor16e-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d17MCem" ]; then   # child 3
            for fileName in `ls $OUTPUTDIR_LHC17d17em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17em $NSlashes "MC_LHC17d17-anchor16g-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17em $NSlashes "MC_LHC17d17-anchor16g-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17f5MCem" ]; then    # child 4
            for fileName in `ls $OUTPUTDIR_LHC17f5em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5em $NSlashes "MC_LHC17f5-anchor16h-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5em $NSlashes "MC_LHC17f5-anchor16h-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d3MCem" ]; then    # child 5
            for fileName in `ls $OUTPUTDIR_LHC17d3em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3em $NSlashes "MC_LHC17d3-anchor16i-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3em $NSlashes "MC_LHC17d3-anchor16i-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17e5MCem" ]; then    # child 6
            for fileName in `ls $OUTPUTDIR_LHC17e5em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5em $NSlashes "MC_LHC17e5-anchor16j-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5em $NSlashes "MC_LHC17e5-anchor16j-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d20a1MCem" ]; then # child 7
            for fileName in `lsem $OUTPUTDIR_LHC17d20a1em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1em $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1em $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `lsem $OUTPUTDIR_LHC17d20a1/GammaConvV1-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCM $fileName $OUTPUTDIR_LHC17d20a1em $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
            for fileName in `lsem $OUTPUTDIR_LHC17d20a1/GammaColo-DPGTracks_*.root`; do
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d20a1em $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
        fi
        if [ -n "$LHC17d20a2MCem" ]; then # child 8
            for fileName in `ls $OUTPUTDIR_LHC17d20a2em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2em $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2em $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d16MCem" ]; then   # child 9
            for fileName in `ls $OUTPUTDIR_LHC17d16em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16em $NSlashes "MC_LHC17d16-anchor16o-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16em $NSlashes "MC_LHC17d16-anchor16o-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
        if [ -n "$LHC17d18MCem" ]; then   # child 10
            for fileName in `ls $OUTPUTDIR_LHC17d18em/GammaConvCalo-DPGTracks_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18em $NSlashes "MC_LHC17d18-anchor16p-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18em/GammaConvCalo-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18em $NSlashes "MC_LHC17d18-anchor16p-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
        fi
    fi # end of [ $DOWNLOADON == 1 ]
#

# Merge ROOT files
    if [ $MERGEON == 1 ]; then
        filesForMerging=`ls $OUTPUTDIR/GammaConvCalo_LHC16d-pass$passNr-DPGTracks\_*.root` # will be used to get the number, so period doesn't matter
        # LHC data files
            echo -e "\n\n--== MERGING LHC DATA FILES ==-\n"
            periodList=(
                d # child 1
                e # child 10
                g # child 2
                h # child 3
                i # child 4
                j # child 5
                k # child 6
                l # child 7
                o # child 8
                p # child 9
            )
            rm -f runlistsToMerge.txt
            echo -e "DPGTracks" >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo "$fileName --> number: $number"
                for runListName in $listsToMerge; do
                # ConvCalo
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for periodID in ${periodList[@]}; do
                        # echo "Period ID: $periodID"
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC16$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile" >> listCurrMerge.txt
                        else
                            echo "--> WARNING: File \"$currFile\" does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvV1_LHC16$nameOut\-pass$passNr-$runListName\_$number.root
                # ConvV1
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for periodID in ${periodList[@]}; do
                        # echo "Period ID: $periodID"
                        currFile=$OUTPUTDIR/GammaConvV1_LHC16$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile" >> listCurrMerge.txt
                        else
                            echo "--> WARNING: File \"$currFile\" does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC16$nameOut\-pass$passNr-$runListName\_$number.root
                # Calo
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for periodID in ${periodList[@]}; do
                        # echo "Period ID: $periodID"
                        currFile=$OUTPUTDIR/GammaCalo_LHC16$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile" >> listCurrMerge.txt
                        else
                            echo "--> WARNING: File \"$currFile\" does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCalo_LHC16$nameOut\-pass$passNr-$runListName\_$number.root
                done
            done
        # Monte Carlo files
            echo -e "\n\n--== MERGING MONTE CARLO FILES ==--\n"
            periodListMC=(
                f6    # child 1  --> anchored to LHC16d
                f9    # child 2  --> anchored to LHC16e
                d17   # child 3  --> anchored to LHC16g
                f5    # child 4  --> anchored to LHC16h
                d3    # child 5  --> anchored to LHC16i
                e5    # child 6  --> anchored to LHC16j
                d20a1 # child 7  --> anchored to LHC16k
                d20a2 # child 8  --> anchored to LHC16l
                d16   # child 9  --> anchored to LHC16o
                d18   # child 10 --> anchored to LHC16p
            )
            rm -f runlistsToMerge.txt
            echo -e "DPGTracks" >> runlistsToMerge.txt
            # echo -e "DPGTracksIncAcc" >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo -e "\n--> Merging $number"
                for runListName in $listsToMerge; do
                # ConvV1
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for ((i = 0; i < ${#periodListMC[@]}; i++)); do
                        # echo "Period ID: 17${periodListMC[i]}-anchor16${periodList[i]}"
                        currFile=$OUTPUTDIR/GammaConvCalo_MC_LHC17${periodListMC[i]}-anchor16${periodList[i]}-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=${periodListMC[i]}
                            echo -e "$currFile" >> listCurrMerge.txt # maybe "$currFile\n"?
                        else
                            echo "--> WARNING: File \"$currFile\" does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17$nameOut\-$runListName\_$number.root
                # Calo
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for ((i = 0; i < ${#periodListMC[@]}; i++)); do
                        # echo "Period ID: 17${periodListMC[i]}-anchor16${periodList[i]}"
                        currFile=$OUTPUTDIR/GammaCalo_MC_LHC17${periodListMC[i]}-anchor16${periodList[i]}-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=${periodListMC[i]}
                            echo -e "$currFile" >> listCurrMerge.txt # maybe "$currFile\n"?
                        else
                            echo "--> WARNING: File \"$currFile\" does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaCalo_MC_LHC17$nameOut\-$runListName\_$number.root
                # ConvCalo
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for ((i = 0; i < ${#periodListMC[@]}; i++)); do
                        # echo "Period ID: 17${periodListMC[i]}-anchor16${periodList[i]}"
                        currFile=$OUTPUTDIR/GammaConvCalo_MC_LHC17${periodListMC[i]}-anchor16${periodList[i]}-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=${periodListMC[i]}
                            echo -e "$currFile" >> listCurrMerge.txt # maybe "$currFile\n"?
                        else
                            echo "--> WARNING: File \"$currFile\" does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17$nameOut\-$runListName\_$number.root
                done
            done
    fi # end of [ $MERGEON == 1 ]
#

# Clean ROOT files
    if [ $CLEANUPMAYOR == 1 ]; then
        echo -e "\n\n--== CLEANING UP ROOT FILES ==--\n"
        # LHC data files
            if [ -n "$LHC16d" ]; then # child 1
                echo "removing all Gamma* files in runFolders for LHC16d"
                rm $OUTPUTDIR_LHC16d/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16d/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16g" ]; then # child 2
                echo "removing all Gamma* files in runFolders for LHC16g"
                rm $OUTPUTDIR_LHC16g/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16g/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16h" ]; then # child 3
                echo "removing all Gamma* files in runFolders for LHC16h"
                rm $OUTPUTDIR_LHC16h/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16h/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16i" ]; then # child 4
                echo "removing all Gamma* files in runFolders for LHC16i"
                rm $OUTPUTDIR_LHC16i/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16i/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16j" ]; then # child 5
                echo "removing all Gamma* files in runFolders for LHC16j"
                rm $OUTPUTDIR_LHC16j/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16j/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16k" ]; then # child 6
                echo "removing all Gamma* files in runFolders for LHC16k"
                rm $OUTPUTDIR_LHC16k/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16k/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16l" ]; then # child 7
                echo "removing all Gamma* files in runFolders for LHC16l"
                rm $OUTPUTDIR_LHC16l/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16l/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16o" ]; then # child 8
                echo "removing all Gamma* files in runFolders for LHC16o"
                rm $OUTPUTDIR_LHC16o/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16o/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16p" ]; then # child 9
                echo "removing all Gamma* files in runFolders for LHC16p"
                rm $OUTPUTDIR_LHC16p/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16p/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC16e" ]; then # child 10
                echo "removing all Gamma* files in runFolders for LHC16e"
                rm $OUTPUTDIR_LHC16e/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC16e/*/*/*Gamma*_*.root
            fi
        # Monte Carlo files PCM
            if [ -n "$LHC17f6MC" ]; then    # child 1
                echo "removing all Gamma* files in runFolders for LHC17f6MC"
                rm $OUTPUTDIR_LHC17f6/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17f6/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17f9MC" ]; then    # child 2
                echo "removing all Gamma* files in runFolders for LHC17f9MC"
                rm $OUTPUTDIR_LHC17f9/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17f9/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d17MC" ]; then   # child 3
                echo "removing all Gamma* files in runFolders for LHC17d17MC"
                rm $OUTPUTDIR_LHC17d17/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d17/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17f5MC" ]; then    # child 4
                echo "removing all Gamma* files in runFolders for LHC17f5MC"
                rm $OUTPUTDIR_LHC17f5/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17f5/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d3MC" ]; then    # child 5
                echo "removing all Gamma* files in runFolders for LHC17d3MC"
                rm $OUTPUTDIR_LHC17d3/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d3/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17e5MC" ]; then    # child 6
                echo "removing all Gamma* files in runFolders for LHC17e5MC"
                rm $OUTPUTDIR_LHC17e5/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17e5/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d20a1MC" ]; then # child 7
                echo "removing all Gamma* files in runFolders for LHC17d20a1MC"
                rm $OUTPUTDIR_LHC17d20a1/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d20a1/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d20a2MC" ]; then # child 8
                echo "removing all Gamma* files in runFolders for LHC17d20a2MC"
                rm $OUTPUTDIR_LHC17d20a2/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d20a2/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d16MC" ]; then   # child 9
                echo "removing all Gamma* files in runFolders for LHC17d16MC"
                rm $OUTPUTDIR_LHC17d16/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d16/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d18MC" ]; then   # child 10
                echo "removing all Gamma* files in runFolders for LHC17d18MC"
                rm $OUTPUTDIR_LHC17d18/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d18/*/*/*Gamma*_*.root
            fi
        # Monte Carlo files EMcal only
            if [ -n "$LHC17f6MCem" ]; then    # child 1
                echo "removing all Gamma* files in runFolders for LHC17f6MC"
                rm $OUTPUTDIR_LHC17f6em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17f6em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17f9MCem" ]; then    # child 2
                echo "removing all Gamma* files in runFolders for LHC17f9MC"
                rm $OUTPUTDIR_LHC17f9em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17f9em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d17MCem" ]; then   # child 3
                echo "removing all Gamma* files in runFolders for LHC17d17MC"
                rm $OUTPUTDIR_LHC17d17em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d17em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17f5MCem" ]; then    # child 4
                echo "removing all Gamma* files in runFolders for LHC17f5MC"
                rm $OUTPUTDIR_LHC17f5em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17f5em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d3MCem" ]; then    # child 5
                echo "removing all Gamma* files in runFolders for LHC17d3MC"
                rm $OUTPUTDIR_LHC17d3em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d3em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17e5MCem" ]; then    # child 6
                echo "removing all Gamma* files in runFolders for LHC17e5MC"
                rm $OUTPUTDIR_LHC17e5em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17e5em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d20a1MCem" ]; then # child 7
                echo "removing all Gamma* files in runFolders for LHC17d20a1MC"
                rm $OUTPUTDIR_LHC17d20a1em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d20a1em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d20a2MCem" ]; then # child 8
                echo "removing all Gamma* files in runFolders for LHC17d20a2MC"
                rm $OUTPUTDIR_LHC17d20a2em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d20a2em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d16MCem" ]; then   # child 9
                echo "removing all Gamma* files in runFolders for LHC17d16MC"
                rm $OUTPUTDIR_LHC17d16em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d16em/*/*/*Gamma*_*.root
            fi
            if [ -n "$LHC17d18MCem" ]; then   # child 10
                echo "removing all Gamma* files in runFolders for LHC17d18MC"
                rm $OUTPUTDIR_LHC17d18em/*/Gamma*_*.root
                rm $OUTPUTDIR_LHC17d18em/*/*/*Gamma*_*.root
            fi
    fi # end of [ $CLEANUPMAYOR == 1 ]
