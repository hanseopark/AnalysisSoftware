# WARNING: this script uses regular expressions (e.g. ^...$). Check if your terminal uses the same regex version!

# Run as:
# bash GetHeavyMesonFilesFromGridAndMergeThem_pp_13TeV.sh <your_username>

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=0
SEPARATEON=0
CLEANUP=1
CLEANUPMAYOR=0
MERGEACCROSSYEARS=0 # ! use with caution: read what the function does
number=""
AliTrainDirDt="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp" # data train path on the GRID
AliTrainDirMC="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC" # MC train path on the GRID

# * Print input arguments
    echo "User name: $1"
    # echo "Path: $PATH"
#

# * Set user name and number of slashes
    if [ "$1" = "fbock" ]; then
        BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp13TeV
    elif [ "$1" = "hannahbossi" ]; then
        BASEDIR=/Volumes/external_memory/CERN_data/QA
    elif [ "$1" = "dmuhlhei" ]; then
        BASEDIR=~/data/work/Grid
    elif [ "$1" = "jlueh" ]; then
        BASEDIR=~/Daten/GridDownload
    elif [ "$1" = "redeboer" ]; then
        BASEDIR=~/alice/GridOutput
    else
        echo "ERROR: User name \"$1\" not valid"
        exit 1
    fi
    NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
    NSlashes=`expr $NSlashesBASE + 4`
    NSlashes2=`expr $NSlashes - 1`
    NSlashes3=`expr $NSlashes + 1`
    NSlashes4=`expr $NSlashes + 2`
#

# NOTE: Train configurations
    passNr="1"
    # DEFAULT (for if you comment below input)
        # Data files
            # LHC16_13TeV_pass1 data files
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
            # LHC17_13TeV_pass1 data files
                LHC17Data=""
                LHC17cData=""
                LHC17eData=""
                LHC17fData=""
                LHC17iData=""
                LHC17jData=""
                LHC17kData=""
                LHC17lData=""
                LHC17mData=""
                LHC17oData=""
                LHC17rData=""
                LHC17hData=""
            # LHC18_13TeV_pass1 data files
                LHC18Data=""
                LHC18bData=""
                LHC18dData=""
                LHC18eData=""
                LHC18fData=""
                LHC18gData=""
                LHC18hData=""
                LHC18iData=""
                LHC18jData=""
        # Monte Carlo files
            # LHC17_PYT8_13TeV_anchLHC16 files
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
            # LHC17_PYT8_13TeV_anchLHC16_extra files
                LHC17MCextra=""
                LHC17f6MCextra=""
                LHC17f9MCextra=""
                LHC17d17MCextra=""
                LHC17f5MCextra=""
                LHC17d3MCextra=""
                LHC17e5MCextra=""
                LHC17d20a1MCextra=""
                LHC17d20a2MCextra=""
                LHC17d16MCextra=""
                LHC17d18MCextra=""
    # NOTE: LHC16_13TeV_pass1 DATA files
        # * vAN-20180614
            # LHC16Data="2374" # PCM + PCMEMC + EMC (2018/07/18) # ! test
        # * vAN-20180727
            # LHC16Data="2422" # PCMPHOS + PHOSPHOS (2018/07/27)
            # LHC16Data="2423" # PCMEMC + EMC (2018/07/30)
        # * vAN-20180814
            # LHC16Data="2462" # PCM (2018/08/20)
        # LHC16dData="child_1"
        # LHC16gData="child_2"
        # LHC16hData="child_3"
        # LHC16iData="child_4"
        # LHC16jData="child_5"
        # LHC16kData="child_6"
        # LHC16lData="child_7"
        # LHC16oData="child_8"
        # LHC16pData="child_9"
        # LHC16eData="child_10"
    # NOTE: LHC17_13TeV_pass1 DATA files
        # * vAN-20180814
            LHC17Data="2444" # PCMEMC + EMC (2018/08/15)
            # LHC17Data="2446" # PCM + PCMPHOS + PHOSPHOS (2018/08/15) # ! misses children 6, 7, 8, 9, 10
        LHC17cData="child_1"
        LHC17eData="child_2"
        LHC17fData="child_3"
        LHC17iData="child_4"
        LHC17jData="child_5"
        LHC17kData="child_6"
        LHC17lData="child_7"
        LHC17mData="child_8"
        LHC17oData="child_9"
        LHC17rData="child_10"
        LHC17hData="child_11"
    # NOTE: LHC18_13TeV_pass1 DATA files
        # * vAN-20180814
            # LHC18Data="2445" # PCMEMC + EMC (2018/08/16)
            # LHC18Data="2449" # PCM + PCMPHOS + PHOSPHOS (2018/08/16)
        # LHC18bData="child_1"
        # LHC18dData="child_2"
        # LHC18eData="child_3"
        # LHC18fData="child_4"
        # LHC18gData="child_5"
        # LHC18hData="child_6"
        # LHC18iData="child_7"
        # LHC18jData="child_8"
    # NOTE: Monte Carlo LHC17_PYT8_13TeV_anchLHC16 files
        # * vAN-20180727
            # LHC17MC="3341" # PCM + EMC  (2018/06/18) # ! test
        # * vAN-20180727
            # LHC17MC="3398" # PCM + PHOS (2018/07/27)
        # * vAN-20180727
            # LHC17MC="3409" # PCM + EMC  (2018/08/02)
        # LHC17f6MC="child_1"    # anchored to LHC16d (child 1)
        # LHC17f9MC="child_2"    # anchored to LHC16e (child 10)
        # LHC17d17MC="child_3"   # anchored to LHC16g (child 2)
        # LHC17f5MC="child_4"    # anchored to LHC16h (child 3)
        # LHC17d3MC="child_5"    # anchored to LHC16i (child 4)
        # LHC17e5MC="child_6"    # anchored to LHC16j (child 5)
        # LHC17d20a1MC="child_7" # anchored to LHC16k (child 6)
        # LHC17d20a2MC="child_8" # anchored to LHC16l (child 7)
        # LHC17d16MC="child_9"   # anchored to LHC16o (child 8)
        # LHC17d18MC="child_10"  # anchored to LHC16p (child 9)
    # NOTE: Monte Carlo LHC17_PYT8_13TeV_anchLHC16_extra files EXTRA
        # * vAN-20180727
            # LHC17MCextra="3399" # PCM + PHOS (2018/07/30)
        # * vAN-20180731
            # LHC17MCextra="3410" # PCM + EMC (2018/08/02)
        # LHC17f6MCextra=$LHC17f6MC       # anchored to LHC16d (child 1)
        # LHC17f9MCextra=$LHC17f9MC       # anchored to LHC16e (child 10)
        # LHC17d17MCextra=$LHC17d17MC     # anchored to LHC16g (child 2)
        # LHC17f5MCextra=$LHC17f5MC       # anchored to LHC16h (child 3)
        # LHC17d3MCextra=$LHC17d3MC       # anchored to LHC16i (child 4)
        # LHC17e5MCextra=$LHC17e5MC       # anchored to LHC16j (child 5)
        # LHC17d20a1MCextra=$LHC17d20a1MC # anchored to LHC16k (child 6)
        # LHC17d20a2MCextra=$LHC17d20a2MC # anchored to LHC16l (child 7)
        # LHC17d16MCextra=$LHC17d16MC     # anchored to LHC16o (child 8)
        # LHC17d18MCextra=$LHC17d18MC     # anchored to LHC16p (child 9)
    # NOTE: Which periods do you want to merge?
        periodList16=(
            d # child 1
            g # child 2
            h # child 3
            i # child 4
            j # child 5
            k # child 6
            l # child 7
            o # child 8
            p # child 9
            e # child 10
        )
        periodList17=(
            c # child 1
            e # child 2
            f # child 3
            i # child 4
            j # child 5
            k # child 6
            l # child 7
            m # child 8
            o # child 9
            r # child 10
            h # child 11
        )
        periodList18=(
            b # child 1
            d # child 2
            e # child 3
            f # child 4
            g # child 5
            h # child 6
            i # child 7
            j # child 8
        )
        periodList17MC=(
            f6    # child 1  --> anchored to LHC16d
            d17   # child 3  --> anchored to LHC16g
            f5    # child 4  --> anchored to LHC16h
            d3    # child 5  --> anchored to LHC16i
            e5    # child 6  --> anchored to LHC16j
            d20a1 # child 7  --> anchored to LHC16k
            d20a2 # child 8  --> anchored to LHC16l
            d16   # child 9  --> anchored to LHC16o
            d18   # child 10 --> anchored to LHC16p
            f9    # child 2  --> anchored to LHC16e
        )
#

# NOTE: Set output directory (offline and remote)
    # TRAINDIR="vAN-20180617-1" # 1st EtaPrime analysis # ! test run
    # TRAINDIR="vAN-20180727-1" # 2nd EtaPrime analysis run: includes 2016 data and 2017 MC only
    TRAINDIR="vAN-20180814-1" # 3rd EtaPrime analysis run: includes all DATA from 2016, 2017, and 2018
    OUTPUTDIR=$BASEDIR/$TRAINDIR
#

# * Set local and GRID directory names *
    echo -e "\n\n--== SETTING REMOTE AND LOCAL OUTPUT DIRECTORIES ==--\n"
    mkdir -p $OUTPUTDIR/CutSelections # this folder will contain txt files with cutnumbers
    # LHC16_13TeV_pass1 data files (10 children)
        if [ -n "$LHC16Data" ]; then
            echo -e "\nChildren for \"LHC16_13TeV_pass1\""
            if [ -n "$LHC16dData" ]; then # child 1
                LHC16dData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16dData$"`
                OUTPUTDIR_LHC16d=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16dData
                echo "Directories for LHC16d:"
                echo "  remote \"$LHC16dData\""
                echo "  local  \"$OUTPUTDIR_LHC16d\""
            fi
            if [ -n "$LHC16gData" ]; then # child 2
                LHC16gData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16gData$"`
                OUTPUTDIR_LHC16g=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16gData
                echo "Directories for LHC16g:"
                echo "  remote \"$LHC16gData\""
                echo "  local  \"$OUTPUTDIR_LHC16g\""
            fi
            if [ -n "$LHC16hData" ]; then # child 3
                LHC16hData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16hData$"`
                OUTPUTDIR_LHC16h=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16hData
                echo "Directories for LHC16h:"
                echo "  remote \"$LHC16hData\""
                echo "  local  \"$OUTPUTDIR_LHC16h\""
            fi
            if [ -n "$LHC16iData" ]; then # child 4
                LHC16iData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16iData$"`
                OUTPUTDIR_LHC16i=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16iData
                echo "Directories for LHC16i:"
                echo "  remote \"$LHC16iData\""
                echo "  local  \"$OUTPUTDIR_LHC16i\""
            fi
            if [ -n "$LHC16jData" ]; then # child 5
                LHC16jData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16jData$"`
                OUTPUTDIR_LHC16j=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16jData
                echo "Directories for LHC16j:"
                echo "  remote \"$LHC16jData\""
                echo "  local  \"$OUTPUTDIR_LHC16j\""
            fi
            if [ -n "$LHC16kData" ]; then # child 6
                LHC16kData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16kData$"`
                OUTPUTDIR_LHC16k=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16kData
                echo "Directories for LHC16k:"
                echo "  remote \"$LHC16kData\""
                echo "  local  \"$OUTPUTDIR_LHC16k\""
            fi
            if [ -n "$LHC16lData" ]; then # child 7
                LHC16lData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16lData$"`
                OUTPUTDIR_LHC16l=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16lData
                echo "Directories for LHC16l:"
                echo "  remote \"$LHC16lData\""
                echo "  local  \"$OUTPUTDIR_LHC16l\""
            fi
            if [ -n "$LHC16oData" ]; then # child 8
                LHC16oData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16oData$"`
                OUTPUTDIR_LHC16o=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16oData
                echo "Directories for LHC16o:"
                echo "  remote \"$LHC16oData\""
                echo "  local  \"$OUTPUTDIR_LHC16o\""
            fi
            if [ -n "$LHC16pData" ]; then # child 9
                LHC16pData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16pData$"`
                OUTPUTDIR_LHC16p=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16pData
                echo "Directories for LHC16p:"
                echo "  remote \"$LHC16pData\""
                echo "  local  \"$OUTPUTDIR_LHC16p\""
            fi
            if [ -n "$LHC16eData" ]; then # child 10
                LHC16eData=`alien_ls $AliTrainDirDt/ | grep "^$LHC16Data.*$LHC16eData$"`
                OUTPUTDIR_LHC16e=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC16eData
                echo "Directories for LHC16e:"
                echo "  remote \"$LHC16eData\""
                echo "  local  \"$OUTPUTDIR_LHC16e\""
            fi
        fi
    # LHC17_13TeV_pass1 data files (11 children)
        if [ -n "$LHC17Data" ]; then
            echo -e "\nChildren for \"LHC17_13TeV_pass1\""
            if [ -n "$LHC17cData" ]; then # child 1
                LHC17cData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17cData$"`
                OUTPUTDIR_LHC17c=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17cData
                echo "Directories for LHC17c:"
                echo "  remote \"$LHC17cData\""
                echo "  local  \"$OUTPUTDIR_LHC17c\""
            fi
            if [ -n "$LHC17eData" ]; then # child 2
                LHC17eData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17eData$"`
                OUTPUTDIR_LHC17e=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17eData
                echo "Directories for LHC17e:"
                echo "  remote \"$LHC17eData\""
                echo "  local  \"$OUTPUTDIR_LHC17e\""
            fi
            if [ -n "$LHC17fData" ]; then # child 3
                LHC17fData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17fData$"`
                OUTPUTDIR_LHC17f=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17fData
                echo "Directories for LHC17f:"
                echo "  remote \"$LHC17fData\""
                echo "  local  \"$OUTPUTDIR_LHC17f\""
            fi
            if [ -n "$LHC17iData" ]; then # child 4
                LHC17iData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17iData$"`
                OUTPUTDIR_LHC17i=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17iData
                echo "Directories for LHC17i:"
                echo "  remote \"$LHC17iData\""
                echo "  local  \"$OUTPUTDIR_LHC17i\""
            fi
            if [ -n "$LHC17jData" ]; then # child 5
                LHC17jData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17jData$"`
                OUTPUTDIR_LHC17j=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17jData
                echo "Directories for LHC17j:"
                echo "  remote \"$LHC17jData\""
                echo "  local  \"$OUTPUTDIR_LHC17j\""
            fi
            if [ -n "$LHC17kData" ]; then # child 6
                LHC17kData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17kData$"`
                OUTPUTDIR_LHC17k=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17kData
                echo "Directories for LHC17k:"
                echo "  remote \"$LHC17kData\""
                echo "  local  \"$OUTPUTDIR_LHC17k\""
            fi
            if [ -n "$LHC17lData" ]; then # child 7
                LHC17lData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17lData$"`
                OUTPUTDIR_LHC17l=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17lData
                echo "Directories for LHC17l:"
                echo "  remote \"$LHC17lData\""
                echo "  local  \"$OUTPUTDIR_LHC17l\""
            fi
            if [ -n "$LHC17mData" ]; then # child 8
                LHC17mData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17mData$"`
                OUTPUTDIR_LHC17m=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17mData
                echo "Directories for LHC17m:"
                echo "  remote \"$LHC17mData\""
                echo "  local  \"$OUTPUTDIR_LHC17m\""
            fi
            if [ -n "$LHC17oData" ]; then # child 9
                LHC17oData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17oData$"`
                OUTPUTDIR_LHC17o=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17oData
                echo "Directories for LHC17o:"
                echo "  remote \"$LHC17oData\""
                echo "  local  \"$OUTPUTDIR_LHC17o\""
            fi
            if [ -n "$LHC17rData" ]; then # child 10
                LHC17rData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17rData$"`
                OUTPUTDIR_LHC17r=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17rData
                echo "Directories for LHC17r:"
                echo "  remote \"$LHC17rData\""
                echo "  local  \"$OUTPUTDIR_LHC17r\""
            fi
            if [ -n "$LHC17hData" ]; then # child 11
                LHC17hData=`alien_ls $AliTrainDirDt/ | grep "^$LHC17Data.*$LHC17hData$"`
                OUTPUTDIR_LHC17h=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC17hData
                echo "Directories for LHC17h:"
                echo "  remote \"$LHC17hData\""
                echo "  local  \"$OUTPUTDIR_LHC17h\""
            fi
        fi
    # LHC18_13TeV_pass1 data files (8 children)
        if [ -n "$LHC18Data" ]; then
            echo -e "\nChildren for \"LHC18_13TeV_pass1\""
            if [ -n "$LHC18bData" ]; then # child 1
                LHC18bData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18bData$"`
                OUTPUTDIR_LHC18b=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18bData
                echo "Directories for LHC18b:"
                echo "  remote \"$LHC18bData\""
                echo "  local  \"$OUTPUTDIR_LHC18b\""
            fi
            if [ -n "$LHC18dData" ]; then # child 2
                LHC18dData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18dData$"`
                OUTPUTDIR_LHC18d=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18dData
                echo "Directories for LHC18d:"
                echo "  remote \"$LHC18dData\""
                echo "  local  \"$OUTPUTDIR_LHC18d\""
            fi
            if [ -n "$LHC18eData" ]; then # child 3
                LHC18eData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18eData$"`
                OUTPUTDIR_LHC18e=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18eData
                echo "Directories for LHC18e:"
                echo "  remote \"$LHC18eData\""
                echo "  local  \"$OUTPUTDIR_LHC18e\""
            fi
            if [ -n "$LHC18fData" ]; then # child 4
                LHC18fData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18fData$"`
                OUTPUTDIR_LHC18f=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18fData
                echo "Directories for LHC18f:"
                echo "  remote \"$LHC18fData\""
                echo "  local  \"$OUTPUTDIR_LHC18f\""
            fi
            if [ -n "$LHC18gData" ]; then # child 5
                LHC18gData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18gData$"`
                OUTPUTDIR_LHC18g=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18gData
                echo "Directories for LHC18g:"
                echo "  remote \"$LHC18gData\""
                echo "  local  \"$OUTPUTDIR_LHC18g\""
            fi
            if [ -n "$LHC18hData" ]; then # child 6
                LHC18hData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18hData$"`
                OUTPUTDIR_LHC18h=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18hData
                echo "Directories for LHC18h:"
                echo "  remote \"$LHC18hData\""
                echo "  local  \"$OUTPUTDIR_LHC18h\""
            fi
            if [ -n "$LHC18iData" ]; then # child 7
                LHC18iData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18iData$"`
                OUTPUTDIR_LHC18i=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18iData
                echo "Directories for LHC18i:"
                echo "  remote \"$LHC18iData\""
                echo "  local  \"$OUTPUTDIR_LHC18i\""
            fi
            if [ -n "$LHC18jData" ]; then # child 8
                LHC18jData=`alien_ls $AliTrainDirDt/ | grep "^$LHC18Data.*$LHC18jData$"`
                OUTPUTDIR_LHC18j=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirDt)-$LHC18jData
                echo "Directories for LHC18j:"
                echo "  remote \"$LHC18jData\""
                echo "  local  \"$OUTPUTDIR_LHC18j\""
            fi
        fi
    # Monte Carlo LHC17_PYT8_13TeV_anchLHC16 files (10 children)
        if [ -n "$LHC17MC" ]; then
            echo -e "\nChildren for \"LHC17_PYT8_13TeV_anchLHC16\""
            if [ -n "$LHC17f6MC" ]; then # child 1
                LHC17f6MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f6MC$"`
                OUTPUTDIR_LHC17f6=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f6MC
                echo "Directories for LHC17f6 (MC):"
                echo "  remote \"$LHC17f6MC\""
                echo "  local  \"$OUTPUTDIR_LHC17f6\""
            fi
            if [ -n "$LHC17f9MC" ]; then # child 2
                LHC17f9MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f9MC$"`
                OUTPUTDIR_LHC17f9=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f9MC
                echo "Directories for LHC17f9 (MC):"
                echo "  remote \"$LHC17f9MC\""
                echo "  local  \"$OUTPUTDIR_LHC17f9\""
            fi
            if [ -n "$LHC17d17MC" ]; then # child 3
                LHC17d17MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d17MC$"`
                OUTPUTDIR_LHC17d17=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d17MC
                echo "Directories for LHC17d17 (MC):"
                echo "  remote \"$LHC17d17MC\""
                echo "  local  \"$OUTPUTDIR_LHC17d17\""
            fi
            if [ -n "$LHC17f5MC" ]; then # child 4
                LHC17f5MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17f5MC$"`
                OUTPUTDIR_LHC17f5=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f5MC
                echo "Directories for LHC17f5 (MC):"
                echo "  remote \"$LHC17f5MC\""
                echo "  local  \"$OUTPUTDIR_LHC17f5\""
            fi
            if [ -n "$LHC17d3MC" ]; then # child 5
                LHC17d3MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d3MC$"`
                OUTPUTDIR_LHC17d3=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d3MC
                echo "Directories for LHC17d3 (MC):"
                echo "  remote \"$LHC17d3MC\""
                echo "  local  \"$OUTPUTDIR_LHC17d3\""
            fi
            if [ -n "$LHC17e5MC" ]; then # child 6
                LHC17e5MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17e5MC$"`
                OUTPUTDIR_LHC17e5=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17e5MC
                echo "Directories for LHC17e5 (MC):"
                echo "  remote \"$LHC17e5MC\""
                echo "  local  \"$OUTPUTDIR_LHC17e5\""
            fi
            if [ -n "$LHC17d20a1MC" ]; then # child 7
                LHC17d20a1MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d20a1MC$"`
                OUTPUTDIR_LHC17d20a1=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a1MC
                echo "Directories for LHC17d20a1 (MC):"
                echo "  remote \"$LHC17d20a1MC\""
                echo "  local  \"$OUTPUTDIR_LHC17d20a1\""
            fi
            if [ -n "$LHC17d20a2MC" ]; then # child 8
                LHC17d20a2MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d20a2MC$"`
                OUTPUTDIR_LHC17d20a2=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a2MC
                echo "Directories for LHC17d20a2 (MC):"
                echo "  remote \"$LHC17d20a2MC\""
                echo "  local  \"$OUTPUTDIR_LHC17d20a2\""
            fi
            if [ -n "$LHC17d16MC" ]; then # child 9
                LHC17d16MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d16MC$"`
                OUTPUTDIR_LHC17d16=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d16MC
                echo "Directories for LHC17d16 (MC):"
                echo "  remote \"$LHC17d16MC\""
                echo "  local  \"$OUTPUTDIR_LHC17d16\""
            fi
            if [ -n "$LHC17d18MC" ]; then # child 10
                LHC17d18MC=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MC.*$LHC17d18MC$"`
                OUTPUTDIR_LHC17d18=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d18MC
                echo "Directories for LHC17d18 (MC):"
                echo "  remote \"$LHC17d18MC\""
                echo "  local  \"$OUTPUTDIR_LHC17d18\""
            fi
        fi
    # Monte Carlo LHC17_PYT8_13TeV_anchLHC16_extra files (10 children)
        if [ -n "$LHC17MCextra" ]; then
            echo -e "\nChildren for \"LHC17_PYT8_13TeV_anchLHC16\""
            if [ -n "$LHC17f6MCextra" ]; then # child 1
                LHC17f6MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17f6MCextra$"`
                OUTPUTDIR_LHC17f6extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f6MCextra
                echo "Directories for LHC17f6 (MC extra):"
                echo "  remote \"$LHC17f6MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17f6extra\""
            fi
            if [ -n "$LHC17f9MCextra" ]; then # child 2
                LHC17f9MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17f9MCextra$"`
                OUTPUTDIR_LHC17f9extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f9MCextra
                echo "Directories for LHC17f9 (MC extra):"
                echo "  remote \"$LHC17f9MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17f9extra\""
            fi
            if [ -n "$LHC17d17MCextra" ]; then # child 3
                LHC17d17MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17d17MCextra$"`
                OUTPUTDIR_LHC17d17extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d17MCextra
                echo "Directories for LHC17d17 (MC extra):"
                echo "  remote \"$LHC17d17MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17d17extra\""
            fi
            if [ -n "$LHC17f5MCextra" ]; then # child 4
                LHC17f5MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17f5MCextra$"`
                OUTPUTDIR_LHC17f5extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17f5MCextra
                echo "Directories for LHC17f5 (MC extra):"
                echo "  remote \"$LHC17f5MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17f5extra\""
            fi
            if [ -n "$LHC17d3MCextra" ]; then # child 5
                LHC17d3MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17d3MCextra$"`
                OUTPUTDIR_LHC17d3extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d3MCextra
                echo "Directories for LHC17d3 (MC extra):"
                echo "  remote \"$LHC17d3MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17d3extra\""
            fi
            if [ -n "$LHC17e5MCextra" ]; then # child 6
                LHC17e5MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17e5MCextra$"`
                OUTPUTDIR_LHC17e5extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17e5MCextra
                echo "Directories for LHC17e5 (MC extra):"
                echo "  remote \"$LHC17e5MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17e5extra\""
            fi
            if [ -n "$LHC17d20a1MCextra" ]; then # child 7
                LHC17d20a1MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17d20a1MCextra$"`
                OUTPUTDIR_LHC17d20a1extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a1MCextra
                echo "Directories for LHC17d20a1 (MC extra):"
                echo "  remote \"$LHC17d20a1MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17d20a1extra\""
            fi
            if [ -n "$LHC17d20a2MCextra" ]; then # child 8
                LHC17d20a2MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17d20a2MCextra$"`
                OUTPUTDIR_LHC17d20a2extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d20a2MCextra
                echo "Directories for LHC17d20a2 (MC extra):"
                echo "  remote \"$LHC17d20a2MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17d20a2extra\""
            fi
            if [ -n "$LHC17d16MCextra" ]; then # child 9
                LHC17d16MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17d16MCextra$"`
                OUTPUTDIR_LHC17d16extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d16MCextra
                echo "Directories for LHC17d16 (MC extra):"
                echo "  remote \"$LHC17d16MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17d16extra\""
            fi
            if [ -n "$LHC17d18MCextra" ]; then # child 10
                LHC17d18MCextra=`alien_ls $AliTrainDirMC/ | grep "^$LHC17MCextra.*$LHC17d18MCextra$"`
                OUTPUTDIR_LHC17d18extra=$BASEDIR/$TRAINDIR/$(basename $AliTrainDirMC)-$LHC17d18MCextra
                echo "Directories for LHC17d18 (MC extra):"
                echo "  remote \"$LHC17d18MCextra\""
                echo "  local  \"$OUTPUTDIR_LHC17d18extra\""
            fi
        fi
#

# * Download files from the GRID *
    if [ $DOWNLOADON == 1 ]; then
        echo -e "\n\n--== DOWNLOADING FILES FROM THE GRID ==--\n"
    # * LHC16_13TeV_pass1 data files # ! misses DPGTracksIncAcc # ! EMCgood == DPGTracksAndCalo
        if [ -n "$LHC16dData" ]; then # child 1  # ! misses EMCgood
            echo -e "\n--> Downloading LHC16d..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16dData/merge" DPGTracks       $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16d "$AliTrainDirDt/$LHC16dData/merge" DPGTracksIncAcc $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16gData" ]; then # child 2  # ! misses EMCgood
            echo -e "\n--> Downloading LHC16g..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16gData/merge" DPGTracks       $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16g "$AliTrainDirDt/$LHC16gData/merge" DPGTracksIncAcc $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16hData" ]; then # child 3
            echo -e "\n--> Downloading LHC16h..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16hData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16hData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16h "$AliTrainDirDt/$LHC16hData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16iData" ]; then # child 4
            echo -e "\n--> Downloading LHC16i..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16iData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16iData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16i "$AliTrainDirDt/$LHC16iData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16jData" ]; then # child 5
            echo -e "\n--> Downloading LHC16j..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16jData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16jData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16j "$AliTrainDirDt/$LHC16jData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16kData" ]; then # child 6
            echo -e "\n--> Downloading LHC16k..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16kData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16kData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16k "$AliTrainDirDt/$LHC16kData/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16lData" ]; then # child 7
            echo -e "\n--> Downloading LHC16l..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16lData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16lData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16l "$AliTrainDirDt/$LHC16lData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16oData" ]; then # child 8  # ! misses EMCgood
            echo -e "\n--> Downloading LHC16o..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16oData/merge" DPGTracks       $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16o "$AliTrainDirDt/$LHC16oData/merge" DPGTracksIncAcc $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16pData" ]; then # child 9  # ! misses EMCgood
            echo -e "\n--> Downloading LHC16p..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracks       $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16p "$AliTrainDirDt/$LHC16pData/merge" DPGTracksIncAcc $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC16eData" ]; then # child 10 # ! misses EMCgood
            echo -e "\n--> Downloading LHC16e..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16e "$AliTrainDirDt/$LHC16eData/merge" DPGTracks       $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16e "$AliTrainDirDt/$LHC16eData/merge" DPGTracksIncAcc $NSlashes3 "" kFALSE
        fi
    # * LHC17_13TeV_pass1 data files
        if [ -n "$LHC17cData" ]; then # child 1  # ! misses IncAcc
            echo -e "\n--> Downloading LHC17c..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17c "$AliTrainDirDt/$LHC17cData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17c "$AliTrainDirDt/$LHC17cData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17c "$AliTrainDirDt/$LHC17cData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17eData" ]; then # child 2  # ! misses IncAcc
            echo -e "\n--> Downloading LHC17e..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e "$AliTrainDirDt/$LHC17eData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e "$AliTrainDirDt/$LHC17eData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e "$AliTrainDirDt/$LHC17eData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17fData" ]; then # child 3
            echo -e "\n--> Downloading LHC17f..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f "$AliTrainDirDt/$LHC17fData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f "$AliTrainDirDt/$LHC17fData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f "$AliTrainDirDt/$LHC17fData/merge_runlist_4" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17iData" ]; then # child 4
            echo -e "\n--> Downloading LHC17i..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17i "$AliTrainDirDt/$LHC17iData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17i "$AliTrainDirDt/$LHC17iData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17i "$AliTrainDirDt/$LHC17iData/merge_runlist_4" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17jData" ]; then # child 5  # ! misses IncAcc
            echo -e "\n--> Downloading LHC17j..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17j "$AliTrainDirDt/$LHC17jData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17j "$AliTrainDirDt/$LHC17jData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17j "$AliTrainDirDt/$LHC17jData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17kData" ]; then # child 6
            echo -e "\n--> Downloading LHC17k..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17k "$AliTrainDirDt/$LHC17kData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17k "$AliTrainDirDt/$LHC17kData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17k "$AliTrainDirDt/$LHC17kData/merge_runlist_4" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17lData" ]; then # child 7
            echo -e "\n--> Downloading LHC17l..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17l "$AliTrainDirDt/$LHC17lData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17l "$AliTrainDirDt/$LHC17lData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17l "$AliTrainDirDt/$LHC17lData/merge_runlist_4" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17mData" ]; then # child 8  # ! misses IncAcc
            echo -e "\n--> Downloading LHC17m..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17m "$AliTrainDirDt/$LHC17mData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17m "$AliTrainDirDt/$LHC17mData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17m "$AliTrainDirDt/$LHC17mData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17oData" ]; then # child 9
            echo -e "\n--> Downloading LHC17o..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17o "$AliTrainDirDt/$LHC17oData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17o "$AliTrainDirDt/$LHC17oData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17o "$AliTrainDirDt/$LHC17oData/merge_runlist_4" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17rData" ]; then # child 10 # ! misses IncAcc
            echo -e "\n--> Downloading LHC17r..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17r "$AliTrainDirDt/$LHC17rData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17r "$AliTrainDirDt/$LHC17rData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17r "$AliTrainDirDt/$LHC17rData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17hData" ]; then # child 11
            echo -e "\n--> Downloading LHC17h..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17h "$AliTrainDirDt/$LHC17hData/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17h "$AliTrainDirDt/$LHC17hData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17h "$AliTrainDirDt/$LHC17hData/merge_runlist_4" DPGTracksIncAcc  $NSlashes3 "" kFALSE
        fi
    # * LHC18_13TeV_pass1 data files # ! DPGTracksIncTPC == DPGTracksIncAcc
        if [ -n "$LHC18bData" ]; then # child 1
            echo -e "\n--> Downloading LHC18b..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18b "$AliTrainDirDt/$LHC18bData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18b "$AliTrainDirDt/$LHC18bData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18b "$AliTrainDirDt/$LHC18bData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC18dData" ]; then # child 2
            echo -e "\n--> Downloading LHC18d..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18d "$AliTrainDirDt/$LHC18dData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18d "$AliTrainDirDt/$LHC18dData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18d "$AliTrainDirDt/$LHC18dData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC18eData" ]; then # child 3
            echo -e "\n--> Downloading LHC18e..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18e "$AliTrainDirDt/$LHC18eData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18e "$AliTrainDirDt/$LHC18eData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18e "$AliTrainDirDt/$LHC18eData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC18fData" ]; then # child 4
            echo -e "\n--> Downloading LHC18f..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18f "$AliTrainDirDt/$LHC18fData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18f "$AliTrainDirDt/$LHC18fData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18f "$AliTrainDirDt/$LHC18fData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC18gData" ]; then # child 5
            echo -e "\n--> Downloading LHC18g..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18g "$AliTrainDirDt/$LHC18gData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18g "$AliTrainDirDt/$LHC18gData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18g "$AliTrainDirDt/$LHC18gData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC18hData" ]; then # child 6
            echo -e "\n--> Downloading LHC18h..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18h "$AliTrainDirDt/$LHC18hData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18h "$AliTrainDirDt/$LHC18hData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18h "$AliTrainDirDt/$LHC18hData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC18iData" ]; then # child 7
            echo -e "\n--> Downloading LHC18i..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18i "$AliTrainDirDt/$LHC18iData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18i "$AliTrainDirDt/$LHC18iData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18i "$AliTrainDirDt/$LHC18iData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC18jData" ]; then # child 8
            echo -e "\n--> Downloading LHC18j..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18j "$AliTrainDirDt/$LHC18jData/merge_runlist_1" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18j "$AliTrainDirDt/$LHC18jData/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC18j "$AliTrainDirDt/$LHC18jData/merge_runlist_3" DPGTracksAndCalo $NSlashes3 "" kFALSE
        fi
    # * Monte Carlo LHC17_PYT8_13TeV_anchLHC16 files
        if [ -n "$LHC17f6MC" ]; then    # child 1  # ! misses EMCgood
            echo -e "\n--> Downloading LHC17f6..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6 "$AliTrainDirMC/$LHC17f6MC/merge_runlist_4" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f9MC" ]; then    # child 2  # ! misses EMCgood
            echo -e "\n--> Downloading LHC17f9..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9 "$AliTrainDirMC/$LHC17f9MC/merge_runlist_4" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d17MC" ]; then   # child 3  # ! misses IncAcc and EMCgood
            echo -e "\n--> Downloading LHC17d17..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_2" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17 "$AliTrainDirMC/$LHC17d17MC/merge_runlist_4" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f5MC" ]; then    # child 4
            echo -e "\n--> Downloading LHC17f5..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5 "$AliTrainDirMC/$LHC17f5MC/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d3MC" ]; then    # child 5
            echo -e "\n--> Downloading LHC17d3..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3 "$AliTrainDirMC/$LHC17d3MC/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17e5MC" ]; then    # child 6
            echo -e "\n--> Downloading LHC17e5..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5 "$AliTrainDirMC/$LHC17e5MC/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a1MC" ]; then # child 7
            echo -e "\n--> Downloading LHC17d20a1..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1 "$AliTrainDirMC/$LHC17d20a1MC/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a2MC" ]; then # child 8
            echo -e "\n--> Downloading LHC17d20a2..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2 "$AliTrainDirMC/$LHC17d20a2MC/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d16MC" ]; then   # child 9  # ! misses EMCgood
            echo -e "\n--> Downloading LHC17d16..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16 "$AliTrainDirMC/$LHC17d16MC/merge_runlist_4" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d18MC" ]; then   # child 10 # ! misses EMCgood
            echo -e "\n--> Downloading LHC17d18..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18 "$AliTrainDirMC/$LHC17d18MC/merge_runlist_4" EMCgood          $NSlashes3 "" kFALSE
        fi
    # * Monte Carlo LHC17_PYT8_13TeV_anchLHC16_extra files EXTRA
        if [ -n "$LHC17f6MCextra" ]; then    # child 1 EXTRA # ! misses EMCgood
            echo -e "\n--> Downloading LHC17f6 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6extra "$AliTrainDirMC/$LHC17f6MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6extra "$AliTrainDirMC/$LHC17f6MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6extra "$AliTrainDirMC/$LHC17f6MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6extra "$AliTrainDirMC/$LHC17f6MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f6extra "$AliTrainDirMC/$LHC17f6MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f9MCextra" ]; then    # child 2 EXTRA # ! misses EMCgood
            echo -e "\n--> Downloading LHC17f9 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9extra "$AliTrainDirMC/$LHC17f9MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9extra "$AliTrainDirMC/$LHC17f9MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9extra "$AliTrainDirMC/$LHC17f9MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9extra "$AliTrainDirMC/$LHC17f9MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f9extra "$AliTrainDirMC/$LHC17f9MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d17MCextra" ]; then   # child 3 EXTRA # ! misses IncAcc and EMCgood
            echo -e "\n--> Downloading LHC17d17 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17extra "$AliTrainDirMC/$LHC17d17MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17extra "$AliTrainDirMC/$LHC17d17MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17extra "$AliTrainDirMC/$LHC17d17MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d17extra "$AliTrainDirMC/$LHC17d17MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17f5MCextra" ]; then    # child 4 EXTRA
            echo -e "\n--> Downloading LHC17f5 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5extra "$AliTrainDirMC/$LHC17f5MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5extra "$AliTrainDirMC/$LHC17f5MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5extra "$AliTrainDirMC/$LHC17f5MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5extra "$AliTrainDirMC/$LHC17f5MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5extra "$AliTrainDirMC/$LHC17f5MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17f5extra "$AliTrainDirMC/$LHC17f5MCextra/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d3MCextra" ]; then    # child 5 EXTRA
            echo -e "\n--> Downloading LHC17d3 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3extra "$AliTrainDirMC/$LHC17d3MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3extra "$AliTrainDirMC/$LHC17d3MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3extra "$AliTrainDirMC/$LHC17d3MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3extra "$AliTrainDirMC/$LHC17d3MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3extra "$AliTrainDirMC/$LHC17d3MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d3extra "$AliTrainDirMC/$LHC17d3MCextra/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17e5MCextra" ]; then    # child 6 EXTRA
            echo -e "\n--> Downloading LHC17e5 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5extra "$AliTrainDirMC/$LHC17e5MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5extra "$AliTrainDirMC/$LHC17e5MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5extra "$AliTrainDirMC/$LHC17e5MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5extra "$AliTrainDirMC/$LHC17e5MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5extra "$AliTrainDirMC/$LHC17e5MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17e5extra "$AliTrainDirMC/$LHC17e5MCextra/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a1MCextra" ]; then # child 7 EXTRA
            echo -e "\n--> Downloading LHC17d20a1 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1extra "$AliTrainDirMC/$LHC17d20a1MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1extra "$AliTrainDirMC/$LHC17d20a1MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1extra "$AliTrainDirMC/$LHC17d20a1MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1extra "$AliTrainDirMC/$LHC17d20a1MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1extra "$AliTrainDirMC/$LHC17d20a1MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a1extra "$AliTrainDirMC/$LHC17d20a1MCextra/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d20a2MCextra" ]; then # child 8 EXTRA
            echo -e "\n--> Downloading LHC17d20a2 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2extra "$AliTrainDirMC/$LHC17d20a2MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2extra "$AliTrainDirMC/$LHC17d20a2MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2extra "$AliTrainDirMC/$LHC17d20a2MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2extra "$AliTrainDirMC/$LHC17d20a2MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2extra "$AliTrainDirMC/$LHC17d20a2MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d20a2extra "$AliTrainDirMC/$LHC17d20a2MCextra/merge_runlist_6" EMCgood          $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d16MCextra" ]; then   # child 9 EXTRA  # ! misses EMCgood
            echo -e "\n--> Downloading LHC17d16 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16extra "$AliTrainDirMC/$LHC17d16MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16extra "$AliTrainDirMC/$LHC17d16MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16extra "$AliTrainDirMC/$LHC17d16MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16extra "$AliTrainDirMC/$LHC17d16MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d16extra "$AliTrainDirMC/$LHC17d16MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
        fi
        if [ -n "$LHC17d18MCextra" ]; then   # child 10 EXTRA # ! misses EMCgood
            echo -e "\n--> Downloading LHC17d18 extra..."
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18extra "$AliTrainDirMC/$LHC17d18MCextra/merge_runlist_1" ListAll          $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18extra "$AliTrainDirMC/$LHC17d18MCextra/merge_runlist_2" DPGTracks        $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18extra "$AliTrainDirMC/$LHC17d18MCextra/merge_runlist_3" DPGTracksIncAcc  $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18extra "$AliTrainDirMC/$LHC17d18MCextra/merge_runlist_4" DPGTracksAndCalo $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC17d18extra "$AliTrainDirMC/$LHC17d18MCextra/merge_runlist_5" HadronPID        $NSlashes3 "" kFALSE
        fi
#

# * Copy files with proper name *
        echo -e "\n\n--== CHANGING STRUCTURE IF NEEDED ==--\n"
    # LHC16_13TeV_pass1 data files
        if [ -n "$LHC16dData" ]; then # child 1 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC16d/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16gData" ]; then # child 2 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC16g/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16hData" ]; then # child 3 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16h/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16iData" ]; then # child 4
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16i/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16jData" ]; then # child 5
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16j/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16kData" ]; then # child 6 # !! PASS 2, check expr addition below
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass`expr $passNr + 1`-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass`expr $passNr + 1`-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16k/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass`expr $passNr + 1`-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16lData" ]; then # child 7 # !! PASS 2, check expr addition below
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass`expr $passNr + 1`-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass`expr $passNr + 1`-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16l/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass`expr $passNr + 1`-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16oData" ]; then # child 8 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16o/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16pData" ]; then # child 9 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16p/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC16eData" ]; then # child 10 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC16e/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
    # LHC17_13TeV_pass1 data files
        if [ -n "$LHC17cData" ]; then # child 1 # ! misses IncAcc
            for fileName in `ls $OUTPUTDIR_LHC17c/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17c $NSlashes "LHC17c-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17c/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17c $NSlashes "LHC17c-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17c/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17c $NSlashes "LHC17c-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17eData" ]; then # child 2 # ! misses IncAcc
            for fileName in `ls $OUTPUTDIR_LHC17e/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e $NSlashes "LHC17e-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e $NSlashes "LHC17e-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e $NSlashes "LHC17e-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17fData" ]; then # child 3
            for fileName in `ls $OUTPUTDIR_LHC17f/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f $NSlashes "LHC17f-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f $NSlashes "LHC17f-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f $NSlashes "LHC17f-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17iData" ]; then # child 4
            for fileName in `ls $OUTPUTDIR_LHC17i/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17i $NSlashes "LHC17i-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17i/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17i $NSlashes "LHC17i-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17i/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17i $NSlashes "LHC17i-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17jData" ]; then # child 5 # ! misses IncAcc
            for fileName in `ls $OUTPUTDIR_LHC17j/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17j $NSlashes "LHC17j-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17j/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17j $NSlashes "LHC17j-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17j/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17j $NSlashes "LHC17j-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17kData" ]; then # child 6
            for fileName in `ls $OUTPUTDIR_LHC17k/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17k $NSlashes "LHC17k-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17k/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17k $NSlashes "LHC17k-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17k/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17k $NSlashes "LHC17k-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17lData" ]; then # child 7
            for fileName in `ls $OUTPUTDIR_LHC17l/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17l $NSlashes "LHC17l-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17l/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17l $NSlashes "LHC17l-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17l/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17l $NSlashes "LHC17l-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17mData" ]; then # child 8 # ! misses IncAcc
            for fileName in `ls $OUTPUTDIR_LHC17m/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17m $NSlashes "LHC17m-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17m/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17m $NSlashes "LHC17m-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17m/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17m $NSlashes "LHC17m-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17oData" ]; then # child 9
            for fileName in `ls $OUTPUTDIR_LHC17o/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17o $NSlashes "LHC17o-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17o/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17o $NSlashes "LHC17o-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17o/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17o $NSlashes "LHC17o-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17rData" ]; then # child 10 # ! misses IncAcc
            for fileName in `ls $OUTPUTDIR_LHC17r/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17r $NSlashes "LHC17r-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17r/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17r $NSlashes "LHC17r-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17r/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17r $NSlashes "LHC17r-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC17hData" ]; then # child 11
            for fileName in `ls $OUTPUTDIR_LHC17h/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17h $NSlashes "LHC17h-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17h/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17h $NSlashes "LHC17h-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17h/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17h $NSlashes "LHC17h-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
    # LHC18_13TeV_pass1 data files # ! DPGTracksIncTPC == DPGTracksIncAcc
        if [ -n "$LHC18bData" ]; then # child 1
            for fileName in `ls $OUTPUTDIR_LHC18b/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18b $NSlashes "LHC18b-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18b/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18b $NSlashes "LHC18b-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18b/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18b $NSlashes "LHC18b-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC18dData" ]; then # child 2
            for fileName in `ls $OUTPUTDIR_LHC18d/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18d $NSlashes "LHC18d-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18d/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18d $NSlashes "LHC18d-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18d/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18d $NSlashes "LHC18d-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC18eData" ]; then # child 3
            for fileName in `ls $OUTPUTDIR_LHC18e/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18e $NSlashes "LHC18e-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18e/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18e $NSlashes "LHC18e-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18e/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18e $NSlashes "LHC18e-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC18fData" ]; then # child 4
            for fileName in `ls $OUTPUTDIR_LHC18f/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18f $NSlashes "LHC18f-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18f/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18f $NSlashes "LHC18f-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18f/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18f $NSlashes "LHC18f-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC18gData" ]; then # child 5
            for fileName in `ls $OUTPUTDIR_LHC18g/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18g $NSlashes "LHC18g-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18g/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18g $NSlashes "LHC18g-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18g/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18g $NSlashes "LHC18g-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC18hData" ]; then # child 6
            for fileName in `ls $OUTPUTDIR_LHC18h/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18h $NSlashes "LHC18h-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18h/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18h $NSlashes "LHC18h-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18h/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18h $NSlashes "LHC18h-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC18iData" ]; then # child 7
            for fileName in `ls $OUTPUTDIR_LHC18i/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18i $NSlashes "LHC18i-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18i/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18i $NSlashes "LHC18i-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18i/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18i $NSlashes "LHC18i-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
        if [ -n "$LHC18jData" ]; then # child 8
            for fileName in `ls $OUTPUTDIR_LHC18j/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18j $NSlashes "LHC18j-pass$passNr-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC18j/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18j $NSlashes "LHC18j-pass$passNr-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC18j/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC18j $NSlashes "LHC18j-pass$passNr-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
        fi
    # Monte Carlo LHC17_PYT8_13TeV_anchLHC16 files
        if [ -n "$LHC17f6MC" ]; then    # child 1 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17f9MC" ]; then    # child 2 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17d17MC" ]; then   # child 3 # ! misses IncAcc and EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17f5MC" ]; then    # child 4
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d3MC" ]; then    # child 5
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17e5MC" ]; then    # child 6
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d20a1MC" ]; then # child 7
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d20a2MC" ]; then # child 8
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d16MC" ]; then   # child 9 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17d18MC" ]; then   # child 10 # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-HadronPID" "-HadronPID"
            done
        fi
    # Monte Carlo LHC17_PYT8_13TeV_anchLHC16_extra files EXTRA
        if [ -n "$LHC17f6MCextra" ]; then    # child 1 EXTRA # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17f6extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6extra $NSlashes "MCextra_LHC17f6-anchor16d-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6extra $NSlashes "MCextra_LHC17f6-anchor16d-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6extra $NSlashes "MCextra_LHC17f6-anchor16d-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6extra $NSlashes "MCextra_LHC17f6-anchor16d-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f6extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f6extra $NSlashes "MCextra_LHC17f6-anchor16d-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17f9MCextra" ]; then    # child 2 EXTRA # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17f9extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9extra $NSlashes "MCextra_LHC17f9-anchor16e-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9extra $NSlashes "MCextra_LHC17f9-anchor16e-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9extra $NSlashes "MCextra_LHC17f9-anchor16e-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9extra $NSlashes "MCextra_LHC17f9-anchor16e-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f9extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f9extra $NSlashes "MCextra_LHC17f9-anchor16e-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17d17MCextra" ]; then   # child 3 EXTRA # ! misses IncAcc and EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17d17extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17extra $NSlashes "MCextra_LHC17d17-anchor16g-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17extra $NSlashes "MCextra_LHC17d17-anchor16g-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17extra $NSlashes "MCextra_LHC17d17-anchor16g-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17extra $NSlashes "MCextra_LHC17d17-anchor16g-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d17extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d17extra $NSlashes "MCextra_LHC17d17-anchor16g-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17f5MCextra" ]; then    # child 4 EXTRA
            for fileName in `ls $OUTPUTDIR_LHC17f5extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5extra $NSlashes "MCextra_LHC17f5-anchor16h-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5extra $NSlashes "MCextra_LHC17f5-anchor16h-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5extra $NSlashes "MCextra_LHC17f5-anchor16h-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5extra $NSlashes "MCextra_LHC17f5-anchor16h-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5extra $NSlashes "MCextra_LHC17f5-anchor16h-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17f5extra/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17f5extra $NSlashes "MCextra_LHC17f5-anchor16h-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d3MCextra" ]; then    # child 5 EXTRA
            for fileName in `ls $OUTPUTDIR_LHC17d3extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3extra $NSlashes "MCextra_LHC17d3-anchor16i-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3extra $NSlashes "MCextra_LHC17d3-anchor16i-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3extra $NSlashes "MCextra_LHC17d3-anchor16i-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3extra $NSlashes "MCextra_LHC17d3-anchor16i-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3extra $NSlashes "MCextra_LHC17d3-anchor16i-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d3extra/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d3extra $NSlashes "MCextra_LHC17d3-anchor16i-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17e5MCextra" ]; then    # child 6 EXTRA
            for fileName in `ls $OUTPUTDIR_LHC17e5extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5extra $NSlashes "MCextra_LHC17e5-anchor16j-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5extra $NSlashes "MCextra_LHC17e5-anchor16j-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5extra $NSlashes "MCextra_LHC17e5-anchor16j-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5extra $NSlashes "MCextra_LHC17e5-anchor16j-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5extra $NSlashes "MCextra_LHC17e5-anchor16j-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17e5extra/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17e5extra $NSlashes "MCextra_LHC17e5-anchor16j-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d20a1MCextra" ]; then # child 7 EXTRA
            for fileName in `ls $OUTPUTDIR_LHC17d20a1extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1extra $NSlashes "MCextra_LHC17d20a1-anchor16k-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1extra $NSlashes "MCextra_LHC17d20a1-anchor16k-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1extra $NSlashes "MCextra_LHC17d20a1-anchor16k-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1extra $NSlashes "MCextra_LHC17d20a1-anchor16k-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1extra $NSlashes "MCextra_LHC17d20a1-anchor16k-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a1extra/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a1extra $NSlashes "MCextra_LHC17d20a1-anchor16k-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d20a2MCextra" ]; then # child 8 EXTRA
            for fileName in `ls $OUTPUTDIR_LHC17d20a2extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2extra $NSlashes "MCextra_LHC17d20a2-anchor16l-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2extra $NSlashes "MCextra_LHC17d20a2-anchor16l-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2extra $NSlashes "MCextra_LHC17d20a2-anchor16l-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2extra $NSlashes "MCextra_LHC17d20a2-anchor16l-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2extra $NSlashes "MCextra_LHC17d20a2-anchor16l-HadronPID" "-HadronPID"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d20a2extra/HeavyNeutralMesonToGG-EMCgood_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d20a2extra $NSlashes "MCextra_LHC17d20a2-anchor16l-EMCgood" "-EMCgood"
            done
        fi
        if [ -n "$LHC17d16MCextra" ]; then   # child 9 EXTRA # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17d16extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16extra $NSlashes "MCextra_LHC17d16-anchor16o-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16extra $NSlashes "MCextra_LHC17d16-anchor16o-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16extra $NSlashes "MCextra_LHC17d16-anchor16o-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16extra $NSlashes "MCextra_LHC17d16-anchor16o-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d16extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d16extra $NSlashes "MCextra_LHC17d16-anchor16o-HadronPID" "-HadronPID"
            done
        fi
        if [ -n "$LHC17d18MCextra" ]; then   # child 10 EXTRA # ! misses EMCgood
            for fileName in `ls $OUTPUTDIR_LHC17d18extra/HeavyNeutralMesonToGG-ListAll_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18extra $NSlashes "MCextra_LHC17d18-anchor16p-ListAll" "-ListAll"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18extra/HeavyNeutralMesonToGG-DPGTracks_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18extra $NSlashes "MCextra_LHC17d18-anchor16p-DPGTracks" "-DPGTracks"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18extra/HeavyNeutralMesonToGG-DPGTracksIncAcc_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18extra $NSlashes "MCextra_LHC17d18-anchor16p-DPGTracksIncAcc" "-DPGTracksIncAcc"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18extra/HeavyNeutralMesonToGG-DPGTracksAndCalo_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18extra $NSlashes "MCextra_LHC17d18-anchor16p-DPGTracksAndCalo" "-DPGTracksAndCalo"
            done
            for fileName in `ls $OUTPUTDIR_LHC17d18extra/HeavyNeutralMesonToGG-HadronPID_*.root`; do
                ChangeStructureIfNeededHeavy $fileName $OUTPUTDIR_LHC17d18extra $NSlashes "MCextra_LHC17d18-anchor16p-HadronPID" "-HadronPID"
            done
        fi
    fi # end of [ $DOWNLOADON == 1 ]
#

# * Merge ROOT files *
# NOTE: Change period name in ls command: choose one of the periods you have downloaded
    if [ $MERGEON == 1 ]; then
        echo -e "\n\n--== MERGING ROOT FILES ==--\n"
        # LHC16_13TeV_pass1 data files
        if [ ${#periodList16[@]} -gt 0 ]; then
            echo -e "\n--> Merging LHC16_13TeV_pass1 data files"
            rm -f runlistsToMerge.txt
            echo -e "DPGTracks"        >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            echo -e "DPGTracksIncAcc"  >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMesonToGG_LHC16h-pass$passNr-DPGTracks\_*.root`
            # ! Be careful with pass nr here: periods k and l are pass2
            echo "  \"filesForMerging\" contains ${#filesForMerging[@]} files"
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo "$(basename $fileName) --> number: $number"
                if [ -n "$LHC16Data" ]; then
                    for runListName in $listsToMerge; do
                        echo -e "\nMerging runlist \"$runListName\""
                        rm -f listCurrMerge.txt
                        nameOut=""
                        for periodID in ${periodList16[@]}; do
                            currFile=`ls $OUTPUTDIR/HeavyNeutralMesonToGG_LHC16$periodID-pass?-$runListName\_$number.root`
                                # ? currFile constructed using ls because of different pass numbers
                            if [ -f $currFile ]; then
                                echo "Period ID: $periodID"
                                nameOut+=$periodID
                                echo -e "$currFile" >> listCurrMerge.txt
                            fi
                        done
                        MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/HeavyNeutralMesonToGG_LHC16$nameOut\-pass$passNr-$runListName\_$number.root
                    done
                fi
            done
        fi
        # LHC17_13TeV_pass1 data files
        if [ ${#periodList17[@]} -gt 0 ]; then
            echo -e "\n--> Merging LHC17_13TeV_pass1 data files"
            rm -f runlistsToMerge.txt
            echo -e "DPGTracks"        >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            echo -e "DPGTracksIncAcc"  >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMesonToGG_LHC17c-pass$passNr-DPGTracks\_*.root`
            echo "  \"filesForMerging\" contains ${#filesForMerging[@]} files"
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo "$(basename $fileName) --> number: $number"
                if [ -n "$LHC17Data" ]; then
                    for runListName in $listsToMerge; do
                        echo -e "\nMerging runlist \"$runListName\""
                        rm -f listCurrMerge.txt
                        nameOut=""
                        for periodID in ${periodList17[@]}; do
                            currFile=$OUTPUTDIR/HeavyNeutralMesonToGG_LHC17$periodID-pass$passNr-$runListName\_$number.root
                            if [ -f $currFile ]; then
                                echo "Period ID: $periodID"
                                nameOut+=$periodID
                                echo -e "$currFile" >> listCurrMerge.txt
                            fi
                        done
                        MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/HeavyNeutralMesonToGG_LHC17$nameOut\-pass$passNr-$runListName\_$number.root
                    done
                fi
            done
        fi
        # LHC18_13TeV_pass1 data files
        if [ ${#periodList18[@]} -gt 0 ]; then
            echo -e "\n--> Merging LHC18_13TeV_pass1 data files"
            rm -f runlistsToMerge.txt
            echo -e "DPGTracks"        >> runlistsToMerge.txt
            echo -e "DPGTracksIncAcc"  >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMesonToGG_LHC18b-pass$passNr-DPGTracks\_*.root`
            echo "  \"filesForMerging\" contains ${#filesForMerging[@]} files"
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo "$(basename $fileName) --> number: $number"
                if [ -n "$LHC18Data" ]; then
                    for runListName in $listsToMerge; do
                        echo -e "\nMerging runlist \"$runListName\""
                        rm -f listCurrMerge.txt
                        nameOut=""
                        for periodID in ${periodList18[@]}; do
                            currFile=$OUTPUTDIR/HeavyNeutralMesonToGG_LHC18$periodID-pass$passNr-$runListName\_$number.root
                            if [ -f $currFile ]; then
                                echo "Period ID: $periodID"
                                nameOut+=$periodID
                                echo -e "$currFile" >> listCurrMerge.txt
                            fi
                        done
                        MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/HeavyNeutralMesonToGG_LHC18$nameOut\-pass$passNr-$runListName\_$number.root
                    done
                fi
            done
        fi
        # Monte Carlo LHC17_PYT8_13TeV_anchLHC16 and LHC17_PYT8_13TeV_anchLHC16_extra files
        if [ ${#periodList17MC[@]} -gt 0 ]; then
            echo -e "\n--> Merging Monte Carlo files"
            rm -f runlistsToMerge.txt
            echo -e "ListAll"          >> runlistsToMerge.txt
            echo -e "DPGTracks"        >> runlistsToMerge.txt
            echo -e "DPGTracksIncAcc"  >> runlistsToMerge.txt
            echo -e "DPGTracksAndCalo" >> runlistsToMerge.txt
            echo -e "HadronPID"        >> runlistsToMerge.txt
            echo -e "EMCgood"          >> runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            filesForMerging=`ls $OUTPUTDIR/HeavyNeutralMesonToGG_MC_LHC17d20a2-anchor16*-DPGTracks\_*.root`
            echo "  \"filesForMerging\" contains ${#filesForMerging[@]} files"
            for fileName in $filesForMerging; do
                GetFileNumberMerging $fileName $((NSlashes-1)) 4
                echo -e "\n--> Merging $number"
                for runListName in $listsToMerge; do
                    rm -f listCurrMerge.txt
                    nameOut=""
                    for ((i = 0; i < ${#periodListMC[@]}; i++)); do
                        # 'Normal' Monte Carlo
                        if [ -n "$LHC17MC" ]; then
                            currFile=$OUTPUTDIR/HeavyNeutralMesonToGG_MC_LHC17${periodList17MC[i]}-anchor16${periodList16[i]}-$runListName\_$number.root
                            if [ -f $currFile ]; then
                                echo "Period ID: 17${periodList17MC[i]}-anchor16${periodList16[i]}"
                                nameOut+=${periodListMC[i]}
                                echo -e "$currFile" >> listCurrMerge.txt
                            fi
                        fi
                        # Extra Monte Carlo files
                        if [ -n "$LHC17MCextra" ]; then
                            currFile=$OUTPUTDIR/HeavyNeutralMesonToGG_MCextra_LHC17${periodList17MC[i]}-anchor16${periodList16[i]}-$runListName\_$number.root
                            if [ -f $currFile ]; then
                                echo "Period ID: 17${periodList17MC[i]}-anchor16${periodList16[i]}"
                                # nameOut+=${periodListMC[i]}
                                echo -e "$currFile" >> listCurrMerge.txt
                            fi
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/HeavyNeutralMesonToGG_MC_LHC17$nameOut\-$runListName\_$number.root
                done
            done
            echo "--> Merging done"
        fi
    fi # end of [ $MERGEON == 1 ]
#

# * Clean ROOT files *
    if [ $CLEANUPMAYOR == 1 ]; then
        echo -e "\n--> Cleaning up ROOT files"
    # LHC16_13TeV_pass1 data files
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
    # LHC17_13TeV_pass1 data files
        if [ -n "$LHC17cData" ]; then # child 1
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17c"
            rm $OUTPUTDIR_LHC17c/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17c/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17eData" ]; then # child 2
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17e"
            rm $OUTPUTDIR_LHC17e/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17e/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17fData" ]; then # child 3
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17f"
            rm $OUTPUTDIR_LHC17f/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17f/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17iData" ]; then # child 4
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17i"
            rm $OUTPUTDIR_LHC17i/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17i/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17jData" ]; then # child 5
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17j"
            rm $OUTPUTDIR_LHC17j/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17j/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17kData" ]; then # child 6
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17k"
            rm $OUTPUTDIR_LHC17k/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17k/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17lData" ]; then # child 7
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17l"
            rm $OUTPUTDIR_LHC17l/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17l/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17mData" ]; then # child 8
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17m"
            rm $OUTPUTDIR_LHC17m/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17m/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17oData" ]; then # child 9
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17o"
            rm $OUTPUTDIR_LHC17o/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17o/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17rData" ]; then # child 10
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17r"
            rm $OUTPUTDIR_LHC17r/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17r/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17hData" ]; then # child 11
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17h"
            rm $OUTPUTDIR_LHC17h/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17h/*/*/*HeavyNeutralMesonToGG_*.root
        fi
    # LHC18_13TeV_pass1 data files
        if [ -n "$LHC18dData" ]; then # child 1
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18d"
            rm $OUTPUTDIR_LHC18d/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18d/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC18gData" ]; then # child 2
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18g"
            rm $OUTPUTDIR_LHC18g/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18g/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC18hData" ]; then # child 3
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18h"
            rm $OUTPUTDIR_LHC18h/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18h/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC18iData" ]; then # child 4
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18i"
            rm $OUTPUTDIR_LHC18i/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18i/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC18jData" ]; then # child 5
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18j"
            rm $OUTPUTDIR_LHC18j/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18j/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC18kData" ]; then # child 6
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18k"
            rm $OUTPUTDIR_LHC18k/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18k/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC18lData" ]; then # child 7
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18l"
            rm $OUTPUTDIR_LHC18l/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18l/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC18oData" ]; then # child 8
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC18o"
            rm $OUTPUTDIR_LHC18o/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC18o/*/*/*HeavyNeutralMesonToGG_*.root
        fi
    # Monte Carlo LHC17_PYT8_13TeV_anchLHC16 files
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
    # Monte Carlo LHC17_PYT8_13TeV_anchLHC16_extra files EXTRA
        if [ -n "$LHC17f6MCextra" ]; then    # child 1 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17f6MCextra"
            rm $OUTPUTDIR_LHC17f6extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17f6extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17f9MCextra" ]; then    # child 2 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17f9MCextra"
            rm $OUTPUTDIR_LHC17f9extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17f9extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17d17MCextra" ]; then   # child 3 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d17MCextra"
            rm $OUTPUTDIR_LHC17d17extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17d17extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17f5MCextra" ]; then    # child 4 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17f5MCextra"
            rm $OUTPUTDIR_LHC17f5extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17f5extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17d3MCextra" ]; then    # child 5 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d3MCextra"
            rm $OUTPUTDIR_LHC17d3extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17d3extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17e5MCextra" ]; then    # child 6 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17e5MCextra"
            rm $OUTPUTDIR_LHC17e5extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17e5extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17d20a1MCextra" ]; then # child 7 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d20a1MCextra"
            rm $OUTPUTDIR_LHC17d20a1extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17d20a1extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17d20a2MCextra" ]; then # child 8 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d20a2MCextra"
            rm $OUTPUTDIR_LHC17d20a2extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17d20a2extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17d16MCextra" ]; then   # child 9 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d16MCextra"
            rm $OUTPUTDIR_LHC17d16extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17d16extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
        if [ -n "$LHC17d18MCextra" ]; then   # child 10 EXTRA
            echo "removing all HeavyNeutralMesonToGG files in runFolders for LHC17d18MCextra"
            rm $OUTPUTDIR_LHC17d18extra/*/HeavyNeutralMesonToGG_*.root
            rm $OUTPUTDIR_LHC17d18extra/*/*/*HeavyNeutralMesonToGG_*.root
        fi
    fi # end of [ $CLEANUPMAYOR == 1 ]
#

# * Merge periods accross years *
# NOTE: 2nd argument should be a file name that exists in the download output you generate
    if [ $MERGEACCROSSYEARS == 1 ]; then
        MergeFilesAccrossYears "$BASEDIR/$TRAINDIR" "HeavyNeutralMesonToGG_LHC17c-pass1-DPGTracks_1_400.root"
    fi
#