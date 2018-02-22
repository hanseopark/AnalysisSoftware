#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=0
SINGLERUN=1
SEPARATEON=0
MERGEONSINGLEData=1
MERGEONSINGLEMC=1
CLEANUP=1
CLEANUPMAYOR=$2
number=""


echo $1
echo $2
echo $PATH

# check if train configuration has actually been given

HAVELHC16d=1
HAVELHC16g=1
HAVELHC16h=1
HAVELHC16i=1
HAVELHC16j=1
HAVELHC16k=1
HAVELHC16l=1
HAVELHC16o=1
HAVELHC16p=1
HAVELHC16e=1
HAVETOBUILDData=1
HAVELHC17f6=1;
HAVELHC17f9=1;
HAVELHC17d17=1;
HAVELHC17f5=1;
HAVELHC17d3=1;
HAVELHC17e5=1;
HAVELHC17d20a1=1
HAVELHC17d20a1Ex=1
HAVELHC17d20a2=1
HAVELHC17d20a2Ex=1
HAVETOBUILDMC=1
HAVELHC17d16=1;
HAVELHC17d18=1;

# default trainconfigurations
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
LHC16Data="";
LHC17MCPythia="";
LHC17MCEPOS="";
LHC17f6MC="";
LHC17f9MC="";
LHC17d17MC="";
LHC17f5MC="";
LHC17d3MC="";
LHC17e5MC="";
LHC17d20a1MC="";
LHC17d20a1ExMC="";
LHC17d20a2MC="";
LHC17d20a2ExMC="";
LHC17d16MC="";
LHC17d18MC="";
passNr="1";
NSlashes=10;

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorageExternal/OutputLegoTrains/pp13TeV
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
    NSlashes4=10
elif [ $1 = "hannahbossi" ]; then
    BASEDIR=/Volumes/external_memory/CERN_data/QA
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
    NSlashes4=10
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
    NSlashes=9;
elif [ $1 = "jlueh" ]; then
    BASEDIR=~/Daten/GridDownload
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
    NSlashes4=10
fi

#TRAINDIR=Legotrain-vAN-20171126-1-QA # FIRST TRAIN RUN
TRAINDIR=Legotrain-vAN-20180123-1-QA # SECOND TRAIN RUN

# LHC 16 data
LHC16Data="2296"; #pass 1 SECOND TRAIN RUN
# LHC16dData="child_1"; #pass 1
# LHC16gData="child_2"; #pass 1
# LHC16hData="child_3"; #pass 1
# LHC16iData="child_4"; #pass 1
# LHC16jData="child_5"; #pass 1
# LHC16kData="child_6"; #pass 1
# LHC16lData="child_7"; #pass 1
# LHC16oData="child_8"; #pass 1
# LHC16pData="child_9"; #pass 1
# LHC16eData="child_10"; #pass 1

# LHC17MCPythia="3184"; #pass 1 FIRST TRAIN RUN
LHC17MCPythia="3225"; #pass 1 FIRST TRAIN RUN
LHC17f6MC="child_1";
LHC17f9MC="child_2";
LHC17d17MC="child_3";
LHC17f5MC="child_4";
LHC17d3MC="child_5";
LHC17e5MC="child_6";
LHC17d20a1MC="child_7";
LHC17d20a1ExMC="child_8";
LHC17d20a2MC="child_9";
LHC17d20a2ExMC="child_10";
LHC17d16MC="child_11";
LHC17d18MC="child_12";

OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ "$LHC16dData" == "" ]; then
    HAVELHC16d=0;
fi
if [ "$LHC16gData" == "" ]; then
    HAVELHC16g=0;
fi
if [ "$LHC16hData" == "" ]; then
    HAVELHC16h=0;
fi
if [ "$LHC16iData" == "" ]; then
    HAVELHC16i=0;
fi
if [ "$LHC16jData" == "" ]; then
    HAVELHC16j=0;
fi
if [ "$LHC16kData" == "" ]; then
    HAVELHC16k=0;
fi
if [ "$LHC16lData" = "" ]; then
    HAVELHC16l=0;
fi
if [ "$LHC16oData" = "" ]; then
    HAVELHC16o=0;
fi
if [ "$LHC16pData" = "" ]; then
    HAVELHC16p=0;
fi
if [ "$LHC16eData" = "" ]; then
    HAVELHC16e=0;
fi
if [ "$LHC16Data" != "" ]; then
    HAVETOBUILDData=1;
fi
if [ "$LHC17f6MC" == "" ]; then
    HAVELHC17f6=0;
fi
if [ "$LHC17f9MC" == "" ]; then
    HAVELHC17f9=0;
fi
if [ "$LHC17d17MC" == "" ]; then
    HAVELHC17d17=0;
fi
if [ "$LHC17f5MC" == "" ]; then
    HAVELHC17f5=0;
fi
if [ "$LHC17d3MC" == "" ]; then
    HAVELHC17d3=0;
fi
if [ "$LHC17e5MC" == "" ]; then
    HAVELHC17e5=0;
fi
if [ "$LHC17d20a1MC" == "" ]; then
    HAVELHC17d20a1=0;
    echo $LHC17d20a1MC
fi
if [ "$LHC17d20a1ExMC" == "" ]; then
    HAVELHC17d20a1Ex=0;
fi
if [ "$LHC17d20a2MC" == "" ]; then
    HAVELHC17d20a2=0;
    echo $LHC17d20a2MC
fi
if [ "$LHC17d20a2ExMC" == "" ]; then
    HAVELHC17d20a2Ex=0;
fi
if [ "$LHC17d16MC" == "" ]; then
    HAVELHC17d16=0;
fi
if [ "$LHC17d18MC" == "" ]; then
    HAVELHC17d18=0;
fi
if [ "$LHC17MCPythia" != "" ]; then
    HAVETOBUILDMC=1;
fi
if [ "$LHC17MCEPOS" != "" ]; then
    HAVETOBUILDMC=1;
fi


mkdir -p $OUTPUTDIR/CutSelections
# Get data directory for 16d period
if [ $HAVELHC16d == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC16dData
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC16dData="$InterMediate"
    else
        LHC16dData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16dData\_`
fi
    #if [ "$LHC16dData" == "" ]; then
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC16d=0;
    else
        OUTPUTDIR_LHC16d=$BASEDIR/$TRAINDIR/GA_pp-$LHC16dData
    fi
    rm listGrid.txt
    echo OUTPUTDIR_LHC16d $OUTPUTDIR_LHC16d
fi
# Get data directory for 16g period
if [ $HAVELHC16g == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16gData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16gData`
    else
        LHC16gData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16gData\_`
    fi
    if [ "$LHC16gData" == "" ]; then
        HAVELHC16g=0;
    else
        OUTPUTDIR_LHC16g=$BASEDIR/$TRAINDIR/GA_pp-$LHC16gData
    fi
    echo OUTPUTDIR_LHC16g $OUTPUTDIR_LHC16g
fi
# Get data directory for 16h period
if [ $HAVELHC16h == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16hData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16hData`
    else
        LHC16hData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16hData\_`
    fi
    if [ "$LHC16hData" == "" ]; then
        HAVELHC16h=0;
    else
        OUTPUTDIR_LHC16h=$BASEDIR/$TRAINDIR/GA_pp-$LHC16hData
    fi
    echo OUTPUTDIR_LHC16h $OUTPUTDIR_LHC16h
fi
# Get data directory for 16i period
if [ $HAVELHC16i == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16iData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16iData`
    else
        LHC16iData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16iData\_`
    fi
    if [ "$LHC16iData" == "" ]; then
        HAVELHC16i=0;
    else
        OUTPUTDIR_LHC16i=$BASEDIR/$TRAINDIR/GA_pp-$LHC16iData
    fi
    echo OUTPUTDIR_LHC16i $OUTPUTDIR_LHC16i
fi
# Get data directory for 16j period
if [ $HAVELHC16j == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16jData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16jData`
    else
        LHC16jData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16jData\_`
    fi
    if [ "$LHC16jData" == "" ]; then
        HAVELHC16j=0;
    else
        OUTPUTDIR_LHC16j=$BASEDIR/$TRAINDIR/GA_pp-$LHC16jData
    fi
    echo OUTPUTDIR_LHC16j $OUTPUTDIR_LHC16j
fi
# Get data directory for 16k period
if [ $HAVELHC16k == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16kData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16kData`
    else
        LHC16kData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16kData\_`
    fi
    if [ "$LHC16kData" == "" ]; then
        HAVELHC16k=0;
    else
        OUTPUTDIR_LHC16k=$BASEDIR/$TRAINDIR/GA_pp-$LHC16kData
    fi
    echo OUTPUTDIR_LHC16k $OUTPUTDIR_LHC16k
fi
# Get data directory for 16l period
if [ $HAVELHC16l == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16lData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16lData`
    else
        LHC16lData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16lData\_`
    fi
    if [ "$LHC16lData" == "" ]; then
        HAVELHC16l=0;
    else
        OUTPUTDIR_LHC16l=$BASEDIR/$TRAINDIR/GA_pp-$LHC16lData
    fi
    echo OUTPUTDIR_LHC16l $OUTPUTDIR_LHC16l
fi
# Get data directory for 16o period
if [ $HAVELHC16o == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16oData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16oData`
    else
        LHC16oData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16oData\_`
    fi
    if [ "$LHC16oData" == "" ]; then
        HAVELHC16o=0;
    else
        OUTPUTDIR_LHC16o=$BASEDIR/$TRAINDIR/GA_pp-$LHC16oData
    fi
    echo OUTPUTDIR_LHC16o $OUTPUTDIR_LHC16o
fi
# Get data directory for 16p period
if [ $HAVELHC16p == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC16pData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ | grep $LHC16pData`
    else
        LHC16pData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16pData\_`
    fi
    if [ "$LHC16pData" == "" ]; then
        HAVELHC16p=0;
    else
        OUTPUTDIR_LHC16p=$BASEDIR/$TRAINDIR/GA_pp-$LHC16pData
    fi
    echo OUTPUTDIR_LHC16p $OUTPUTDIR_LHC16p
fi
# Get data directory for 16e period
if [ $HAVELHC16e == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16Data\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC16eData
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC16eData="$InterMediate"
    else
        LHC16eData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC16eData\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC16e=0;
    else
        OUTPUTDIR_LHC16e=$BASEDIR/$TRAINDIR/GA_pp-$LHC16eData
    fi
    rm listGrid.txt
    echo OUTPUTDIR_LHC16e $OUTPUTDIR_LHC16e
fi

if [ $HAVELHC17f6 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17f6MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17f6MC="$InterMediate"
    else
        LHC17f6MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17f6MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17f6=0;
    else
       OUTPUTDIR_LHC17f6=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f6MC
    fi
    rm listGrid.txt
    echo OUTPUTDIR_LHC17f6 $OUTPUTDIR_LHC17f6
fi

if [ $HAVELHC17f9 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17f9MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17f9MC="$InterMediate"
    else
        LHC17f9MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17f9MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17f9=0;
    else
        OUTPUTDIR_LHC17f9=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f9MC
    fi
rm listGrid.txt
echo OUTPUTDIR_LHC17f9 $OUTPUTDIR_LHC17f9
fi

if [ $HAVELHC17d17 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17d17MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17d17MC="$InterMediate"
    else
        LHC17d17MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d17MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17d17=0;
    else
        OUTPUTDIR_LHC17d17=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d17MC
    fi
rm listGrid.txt
echo OUTPUTDIR_LHC17d17 $OUTPUTDIR_LHC17d17
fi

if [ $HAVELHC17f5 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17f5MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17f5MC="$InterMediate"
    else
        LHC17f5MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17f5MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17f5=0;
    else
        OUTPUTDIR_LHC17f5=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f5MC
    fi
rm listGrid.txt
echo OUTPUTDIR_LHC17f5 $OUTPUTDIR_LHC17f5
fi

if [ $HAVELHC17d3 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17d3MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17d3MC="$InterMediate"
    else
        LHC17d3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d3MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17d3=0;
    else
        OUTPUTDIR_LHC17d3=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d3MC
    fi
rm listGrid.txt
echo OUTPUTDIR_LHC17d3 $OUTPUTDIR_LHC17d3
fi

if [ $HAVELHC17e5 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17e5MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17e5MC="$InterMediate"
    else
        LHC17MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17e5MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17e5=0;
    else
        OUTPUTDIR_LHC17e5=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17e5MC
    fi
rm listGrid.txt
echo OUTPUTDIR_LHC17e5 $OUTPUTDIR_LHC17e5
fi

# Get MC directory for LHC17d20a1 MC anchored to LHC16k
if [ $HAVELHC17d20a1 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a1MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a1MC`
    else
        LHC17d20a1MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a1MC\_`
    fi
    if [ "$LHC17d20a1MC" == "" ]; then
        HAVELHC17d20a1=0;
    else
        OUTPUTDIR_LHC17d20a1=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a1MC
    fi
    echo OUTPUTDIR_LHC17d20a1 $OUTPUTDIR_LHC17d20a1
fi

# Get MC directory for LHC17d20a1_extra MC anchored to LHC16k
if [ $HAVELHC17d20a1Ex == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a1ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a1ExMC`
    else
        LHC17d20a1ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a1ExMC\_`
    fi
    if [ "$LHC17d20a1ExMC" == "" ]; then
        HAVELHC17d20a1Ex=0;
    else
        OUTPUTDIR_LHC17d20a1Ex=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a1ExMC
    fi
    echo OUTPUTDIR_LHC17d20a1Ex $OUTPUTDIR_LHC17d20a1Ex
fi

# Get MC directory for LHC17d20a2 MC anchored to LHC16k
if [ $HAVELHC17d20a2 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a2MC`
    else
        LHC17d20a2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a2MC\_`
    fi
    if [ "$LHC17d20a2MC" == "" ]; then
        HAVELHC17d20a2=0;
    else
        OUTPUTDIR_LHC17d20a2=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a2MC
    fi
    echo OUTPUTDIR_LHC17d20a2 $OUTPUTDIR_LHC17d20a2
fi

# Get MC directory for LHC17d20a2_extra MC anchored to LHC16k
if [ $HAVELHC17d20a2Ex == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC17d20a2ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ | grep $LHC17d20a2ExMC`
    else
        LHC17d20a2ExMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d20a2ExMC\_`
    fi
    if [ "$LHC17d20a2ExMC" == "" ]; then
        HAVELHC17d20a2Ex=0;
    else
        OUTPUTDIR_LHC17d20a2Ex=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a2ExMC
    fi
    echo OUTPUTDIR_LHC17d20a2Ex $OUTPUTDIR_LHC17d20a2Ex
fi

if [ $HAVELHC17d16 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17d16MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17d16MC="$InterMediate"
    else
        LHC17d16MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d16MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17d16=0;
    else
        OUTPUTDIR_LHC17d16=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d16MC
    fi
rm listGrid.txt
echo OUTPUTDIR_LHC17d16 $OUTPUTDIR_LHC17d16
fi

if [ $HAVELHC17d18 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17MCPythia\_ > listGrid.txt
        sort listGrid.txt -o listGrid.txt
        InterMediate=`head -n1 listGrid.txt`"_"$LHC17d18MC
        InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
        LHC17d18MC="$InterMediate"
    else
        LHC17d18MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC17d18MC\_`
    fi
    if [ "$InterMediateExists" == "" ]; then
        HAVELHC17d18=0;
    else
        OUTPUTDIR_LHC17d18=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d18MC
    fi
rm listGrid.txt
echo OUTPUTDIR_LHC17d18 $OUTPUTDIR_LHC17d18
fi

if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC16d == 1 ]; then
        echo "downloading LHC16d"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16d_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16d/$runNumber "/alice/data/2016/LHC16d/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16dData" $NSlashes3 "/alice/data/2016/LHC16d/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16dData" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16d/mergedAllCalo.txt ]; then
                rm $OUTPUTDIR_LHC16d/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16d_pass1.txt`
                ls $OUTPUTDIR_LHC16d/$firstrunNumber/GammaCalo_*.root > fileLHC16d.txt

                MergeAccordingToSpecificRunlist fileLHC16d.txt $OUTPUTDIR_LHC16d $NSlashes3 GammaCalo All runlists/runNumbersLHC16d_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16d.txt $OUTPUTDIR_LHC16d $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16d_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16d.txt $OUTPUTDIR_LHC16d $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16d_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16d "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16dData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16g == 1 ]; then
        echo "downloading LHC16g"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16g_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16g/$runNumber "/alice/data/2016/LHC16g/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16gData" $NSlashes3 "/alice/data/2016/LHC16g/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16gData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16g/mergedAllCalo.txt ]; then
                rm $OUTPUTDIR_LHC16g/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16g_pass1.txt`
                ls $OUTPUTDIR_LHC16g/$firstrunNumber/GammaCalo_*.root > fileLHC16g.txt

                MergeAccordingToSpecificRunlist fileLHC16g.txt $OUTPUTDIR_LHC16g $NSlashes3 GammaCalo All runlists/runNumbersLHC16g_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16g.txt $OUTPUTDIR_LHC16g $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16g_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16g.txt $OUTPUTDIR_LHC16g $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16g_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16gData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16h == 1 ]; then
        echo "downloading LHC16h"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16h_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16h/$runNumber "/alice/data/2016/LHC16h/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16hData" $NSlashes3 "/alice/data/2016/LHC16h/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16hData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16h/mergedAllCalo.txt ]; then
                rm $OUTPUTDIR_LHC16h/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16h_pass1.txt`
                ls $OUTPUTDIR_LHC16h/$firstrunNumber/GammaCalo_*.root > fileLHC16h.txt

                MergeAccordingToSpecificRunlist fileLHC16h.txt $OUTPUTDIR_LHC16h $NSlashes3 GammaCalo All runlists/runNumbersLHC16h_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16h.txt $OUTPUTDIR_LHC16h $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16h_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16h.txt $OUTPUTDIR_LHC16h $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16h_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16h "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16hData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16i == 1 ]; then
        echo "downloading LHC16i"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16i_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16i/$runNumber "/alice/data/2016/LHC16i/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16iData" $NSlashes3 "/alice/data/2016/LHC16i/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16iData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16i/mergedAllCalo.txt ]; then
                rm $OUTPUTDIR_LHC16i/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16i_pass1.txt`
                ls $OUTPUTDIR_LHC16i/$firstrunNumber/GammaCalo_*.root > fileLHC16i.txt

                MergeAccordingToSpecificRunlist fileLHC16i.txt $OUTPUTDIR_LHC16i $NSlashes3 GammaCalo All runlists/runNumbersLHC16i_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16i.txt $OUTPUTDIR_LHC16i $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16i_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16i.txt $OUTPUTDIR_LHC16i $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16i_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16i "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16iData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16j == 1 ]; then
        echo "downloading LHC16j"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16j_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16j/$runNumber "/alice/data/2016/LHC16j/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16jData" $NSlashes3 "/alice/data/2016/LHC16j/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16jData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16j/mergedAllCalo.txt ]; then
                rm $OUTPUTDIR_LHC16j/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16j_pass1.txt`
                ls $OUTPUTDIR_LHC16j/$firstrunNumber/GammaCalo_*.root > fileLHC16j.txt

                MergeAccordingToSpecificRunlist fileLHC16j.txt $OUTPUTDIR_LHC16j $NSlashes3 GammaCalo All runlists/runNumbersLHC16j_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16j.txt $OUTPUTDIR_LHC16j $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16j_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16j.txt $OUTPUTDIR_LHC16j $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16j_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16j "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16jData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16k == 1 ]; then
        echo "downloading LHC16k"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16k_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16k/$runNumber "/alice/data/2016/LHC16k/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16kData" $NSlashes3 "/alice/data/2016/LHC16k/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16kData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16k/mergedAllCalo.txt ]; then
                rm $OUTPUTDIR_LHC16k/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16k_pass1.txt`
                ls $OUTPUTDIR_LHC16k/$firstrunNumber/GammaCalo_*.root > fileLHC16k.txt

                MergeAccordingToSpecificRunlist fileLHC16k.txt $OUTPUTDIR_LHC16k $NSlashes3 GammaCalo All runlists/runNumbersLHC16k_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16k.txt $OUTPUTDIR_LHC16k $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16k_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16k.txt $OUTPUTDIR_LHC16k $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16k_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16k "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16kData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16l == 1 ]; then
        echo "downloading LHC16l"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16l_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16l/$runNumber "/alice/data/2016/LHC16l/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16lData" $NSlashes3 "/alice/data/2016/LHC16l/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16lData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16l/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16l/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16l_pass1.txt`
                ls $OUTPUTDIR_LHC16l/$firstrunNumber/GammaCalo_*.root > fileLHC16l.txt
                fileNumbers=`cat fileLHC16l.txt`
                MergeAccordingToSpecificRunlist fileLHC16l.txt $OUTPUTDIR_LHC16l $NSlashes3 GammaCalo All runlists/runNumbersLHC16l_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16l.txt $OUTPUTDIR_LHC16l $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16l_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16l.txt $OUTPUTDIR_LHC16l $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16l_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16l "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16lData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16o == 1 ]; then
        echo "downloading LHC16o"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16o_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16o/$runNumber "/alice/data/2016/LHC16o/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16oData" $NSlashes3 "/alice/data/2016/LHC16o/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16oData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16o/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16o/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16o_pass1.txt`
                ls $OUTPUTDIR_LHC16o/$firstrunNumber/GammaCalo_*.root > fileLHC16o.txt
                fileNumbers=`cat fileLHC16o.txt`
                MergeAccordingToSpecificRunlist fileLHC16o.txt $OUTPUTDIR_LHC16o $NSlashes3 GammaCalo All runlists/runNumbersLHC16o_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16o.txt $OUTPUTDIR_LHC16o $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16o_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16o.txt $OUTPUTDIR_LHC16o $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16o_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16o "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16oData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16p == 1 ]; then
        echo "downloading LHC16p"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16p_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16p/$runNumber "/alice/data/2016/LHC16p/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16pData" $NSlashes3 "/alice/data/2016/LHC16p/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16pData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16p/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16p/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16p_pass1.txt`
                ls $OUTPUTDIR_LHC16p/$firstrunNumber/GammaCalo_*.root > fileLHC16p.txt
                fileNumbers=`cat fileLHC16p.txt`
                MergeAccordingToSpecificRunlist fileLHC16p.txt $OUTPUTDIR_LHC16p $NSlashes3 GammaCalo All runlists/runNumbersLHC16p_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16p.txt $OUTPUTDIR_LHC16p $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16p_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16p.txt $OUTPUTDIR_LHC16p $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16p_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16p "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16pData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC16e == 1 ]; then
        echo "downloading LHC16e"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16e_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC16e/$runNumber "/alice/data/2016/LHC16e/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16eData" $NSlashes3 "/alice/data/2016/LHC16e/000$runNumber/pass$passNr/PWGGA/GA_pp/$LHC16eData/" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC16e/mergedAllConv.txt ]; then
                rm $OUTPUTDIR_LHC16e/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16e_pass1.txt`
                ls $OUTPUTDIR_LHC16e/$firstrunNumber/GammaCalo_*.root > fileLHC16e.txt
                fileNumbers=`cat fileLHC16p.txt`
                MergeAccordingToSpecificRunlist fileLHC16e.txt $OUTPUTDIR_LHC16e $NSlashes3 GammaCalo All runlists/runNumbersLHC16e_pass1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16e.txt $OUTPUTDIR_LHC16e $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC16e_pass1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC16e.txt $OUTPUTDIR_LHC16e $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC16e_pass1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC16p "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16pData/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi

    currentDir=$PWD
    if [ $HAVELHC17f6 == 1 ]; then
        echo "downloading LHC17f6"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f6.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f6/$runNumber "/alice/sim/2017/LHC17f6/$runNumber/PWGGA/GA_pp_MC/$LHC17f6MC" $NSlashes3 "/alice/sim/2017/LHC17f6/$runNumber/PWGGA/GA_pp_MC/$LHC17f6MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17f6/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f6/GammaCalo*.root*
                echo runlists/runNumbersLHC17f6.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f6.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17f6/$firstrunNumber/GammaCalo_*.root > fileLHC17f6.txt
                MergeAccordingToSpecificRunlist fileLHC17f6.txt $OUTPUTDIR_LHC17f6 $NSlashes3 GammaCalo All runlists/runNumbersLHC17f6.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f6.txt $OUTPUTDIR_LHC17f6 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17f6_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f6.txt $OUTPUTDIR_LHC17f6 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17f6_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17f6 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f6MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17f9 == 1 ]; then
        echo "downloading LHC17f9"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f9.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f9/$runNumber "/alice/sim/2017/LHC17f9/$runNumber/PWGGA/GA_pp_MC/$LHC17f9MC" $NSlashes3 "/alice/sim/2017/LHC17f9/$runNumber/PWGGA/GA_pp_MC/$LHC17f9MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17f9/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f9/GammaCalo*.root*
                echo runlists/runNumbersLHC17f9.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f9.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17f9/$firstrunNumber/GammaCalo_*.root > fileLHC17f9.txt
                MergeAccordingToSpecificRunlist fileLHC17f9.txt $OUTPUTDIR_LHC17f9 $NSlashes3 GammaCalo All runlists/runNumbersLHC17f9.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f9.txt $OUTPUTDIR_LHC17f9 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17f9_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f9.txt $OUTPUTDIR_LHC17f9 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17f9_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17f9 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f9MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d17 == 1 ]; then
        echo "downloading LHC17d17"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d17.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d17/$runNumber "/alice/sim/2017/LHC17d17/$runNumber/PWGGA/GA_pp_MC/$LHC17d17MC" $NSlashes3 "/alice/sim/2017/LHC17d17/$runNumber/PWGGA/GA_pp_MC/$LHC17d17MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d17/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d17/GammaCalo*.root*
                echo runlists/runNumbersLHC17d17.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d17.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d17/$firstrunNumber/GammaCalo_*.root > fileLHC17d17.txt
                MergeAccordingToSpecificRunlist fileLHC17d17.txt $OUTPUTDIR_LHC17d17 $NSlashes3 GammaCalo All runlists/runNumbersLHC17d17.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d17.txt $OUTPUTDIR_LHC17d17 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d17_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d17.txt $OUTPUTDIR_LHC17d17 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d17_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d17 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d17MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17f5 == 1 ]; then
        echo "downloading LHC17f5"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17f5.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17f5/$runNumber "/alice/sim/2017/LHC17f5/$runNumber/PWGGA/GA_pp_MC/$LHC17f5MC" $NSlashes3 "/alice/sim/2017/LHC17f5/$runNumber/PWGGA/GA_pp_MC/$LHC17f5MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17f5/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17f5/GammaCalo*.root*
                echo runlists/runNumbersLHC17f5.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17f5.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17f5/$firstrunNumber/GammaCalo_*.root > fileLHC17f5.txt
                MergeAccordingToSpecificRunlist fileLHC17f5.txt $OUTPUTDIR_LHC17f5 $NSlashes3 GammaCalo All runlists/runNumbersLHC17f5.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f5.txt $OUTPUTDIR_LHC17f5 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17f5_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17f5.txt $OUTPUTDIR_LHC17f5 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17f5_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17f5 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f5MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d3 == 1 ]; then
        echo "downloading LHC17d3"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d3.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d3/$runNumber "/alice/sim/2017/LHC17d3/$runNumber/PWGGA/GA_pp_MC/$LHC17d3MC" $NSlashes3 "/alice/sim/2017/LHC17d3/$runNumber/PWGGA/GA_pp_MC/$LHC17d3MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d3/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d3/GammaCalo*.root*
                echo runlists/runNumbersLHC17d3.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d3.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d3/$firstrunNumber/GammaCalo_*.root > fileLHC17d3.txt
                MergeAccordingToSpecificRunlist fileLHC17d3.txt $OUTPUTDIR_LHC17d3 $NSlashes3 GammaCalo All runlists/runNumbersLHC17d3.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d3.txt $OUTPUTDIR_LHC17d3 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d3_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d3.txt $OUTPUTDIR_LHC17d3 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d3_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d3MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17e5 == 1 ]; then
        echo "downloading LHC17e5"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17e5.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17e5/$runNumber "/alice/sim/2017/LHC17e5/$runNumber/PWGGA/GA_pp_MC/$LHC17e5MC" $NSlashes3 "/alice/sim/2017/LHC17e5/$runNumber/PWGGA/GA_pp_MC/$LHC17e5MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17e5/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17e5/GammaCalo*.root*
                echo runlists/runNumbersLHC17e5.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17e5.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17e5/$firstrunNumber/GammaCalo_*.root > fileLHC17e5.txt
                MergeAccordingToSpecificRunlist fileLHC17e5.txt $OUTPUTDIR_LHC17e5 $NSlashes3 GammaCalo All runlists/runNumbersLHC17e5.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17e5.txt $OUTPUTDIR_LHC17e5 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17e5_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17e5.txt $OUTPUTDIR_LHC17e5 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17e5_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17e5 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17e5MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d20a1 == 1 ]; then
        echo "downloading LHC17d20a1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1/$runNumber "/alice/sim/2017/LHC17d20a1/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1MC" $NSlashes3 "/alice/sim/2017/LHC17d20a1/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d20a1/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a1/GammaCalo*.root*
                echo runlists/runNumbersLHC17d20a1.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a1.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d20a1/$firstrunNumber/GammaCalo_*.root > fileLHC17d20a1.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a1.txt $OUTPUTDIR_LHC17d20a1 $NSlashes3 GammaCalo All runlists/runNumbersLHC17d20a1.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a1.txt $OUTPUTDIR_LHC17d20a1 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d20a1_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a1.txt $OUTPUTDIR_LHC17d20a1 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a1_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a1MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d20a1Ex == 1 ]; then
        echo "downloading LHC17d20a1_extra"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a1_extra.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1Ex/$runNumber "/alice/sim/2017/LHC17d20a1_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1ExMC" $NSlashes3 "/alice/sim/2017/LHC17d20a1_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a1ExMC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17d20a1Ex/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a1Ex/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a1_extra.txt`
                ls $OUTPUTDIR_LHC17d20a1Ex/$firstrunNumber/GammaCalo_*.root > fileLHC17d20a1_extra.txt
                fileNumbers=`cat fileLHC17d20a1_extra.txt`
                MergeAccordingToSpecificRunlist fileLHC17d20a1_extra.txt $OUTPUTDIR_LHC17d20a1Ex $NSlashes3 GammaCalo All runlists/runNumbersLHC17d20a1_extra.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a1_extra.txt $OUTPUTDIR_LHC17d20a1Ex $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d20a1_extra_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a1_extra.txt $OUTPUTDIR_LHC17d20a1Ex $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a1_extra_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a1Ex "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a1ExMC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d20a2 == 1 ]; then
        echo "downloading LHC17d20a2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a2.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2/$runNumber "/alice/sim/2017/LHC17d20a2/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2MC" $NSlashes3 "/alice/sim/2017/LHC17d20a2/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d20a2/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a2/GammaCalo*.root*
                echo runlists/runNumbersLHC17d20a2.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a2.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d20a2/$firstrunNumber/GammaCalo_*.root > fileLHC17d20a2.txt
                MergeAccordingToSpecificRunlist fileLHC17d20a2.txt $OUTPUTDIR_LHC17d20a2 $NSlashes3 GammaCalo All runlists/runNumbersLHC17d20a2.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a2.txt $OUTPUTDIR_LHC17d20a2 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d20a2_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a2.txt $OUTPUTDIR_LHC17d20a2 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a2_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a2MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d20a2Ex == 1 ]; then
        echo "downloading LHC17d20a2_extra"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d20a2_extra.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2Ex/$runNumber "/alice/sim/2017/LHC17d20a2_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2ExMC" $NSlashes3 "/alice/sim/2017/LHC17d20a2_extra/$runNumber/PWGGA/GA_pp_MC/$LHC17d20a2ExMC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC17d20a2Ex/mergedAllConv.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d20a2Ex/GammaCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d20a2_extra.txt`
                ls $OUTPUTDIR_LHC17d20a2Ex/$firstrunNumber/GammaCalo_*.root > fileLHC17d20a2_extra.txt
                fileNumbers=`cat fileLHC17d20a2_extra.txt`
                MergeAccordingToSpecificRunlist fileLHC17d20a2_extra.txt $OUTPUTDIR_LHC17d20a2Ex $NSlashes3 GammaCalo All runlists/runNumbersLHC17d20a2_extra.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a2_extra.txt $OUTPUTDIR_LHC17d20a2Ex $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d20a2_extra_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d20a2_extra.txt $OUTPUTDIR_LHC17d20a2Ex $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d20a2_extra_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d20a2Ex "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a2ExMC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d16 == 1 ]; then
        echo "downloading LHC17d16"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d16.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d16/$runNumber "/alice/sim/2017/LHC17d16/$runNumber/PWGGA/GA_pp_MC/$LHC17d16MC" $NSlashes3 "/alice/sim/2017/LHC17d16/$runNumber/PWGGA/GA_pp_MC/$LHC17d16MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d16/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d16/GammaCalo*.root*
                echo runlists/runNumbersLHC17d16.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d16.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d16/$firstrunNumber/GammaCalo_*.root > fileLHC17d16.txt
                MergeAccordingToSpecificRunlist fileLHC17d16.txt $OUTPUTDIR_LHC17d16 $NSlashes3 GammaCalo All runlists/runNumbersLHC17d16.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d16.txt $OUTPUTDIR_LHC17d16 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d16_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d16.txt $OUTPUTDIR_LHC17d16 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d16_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d16 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d16MC/merge" $NSlashes "" kTRUE
        fi
    fi
    if [ $HAVELHC17d18 == 1 ]; then
        echo "downloading LHC17d18"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC17d18.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC17d18/$runNumber "/alice/sim/2017/LHC17d18/$runNumber/PWGGA/GA_pp_MC/$LHC17d18MC" $NSlashes3 "/alice/sim/2017/LHC17d18/$runNumber/PWGGA/GA_pp_MC/$LHC17d18MC/Stage_1/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC17d18/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $OUTPUTDIR_LHC17d18/GammaCalo*.root*
                echo runlists/runNumbersLHC17d18.txt
                firstrunNumber=`head -n1 runlists/runNumbersLHC17d18.txt`
                echo $firstrunNumber
                ls $OUTPUTDIR_LHC17d18/$firstrunNumber/GammaCalo_*.root > fileLHC17d18.txt
                MergeAccordingToSpecificRunlist fileLHC17d18.txt $OUTPUTDIR_LHC17d18 $NSlashes3 GammaCalo All runlists/runNumbersLHC17d18.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d18.txt $OUTPUTDIR_LHC17d18 $NSlashes3 GammaCalo DPGTrack runlists/runNumbersLHC17d18_DPG.txt "no"
                MergeAccordingToSpecificRunlist fileLHC17d18.txt $OUTPUTDIR_LHC17d18 $NSlashes3 GammaCalo DPGTrackAndCalo runlists/runNumbersLHC17d18_DPGTrackAndCalo.txt "no"
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC17d18 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d18MC/merge" $NSlashes "" kTRUE
        fi
    fi
echo "Change Structure If Needed"
    if [ $HAVELHC16d == 1 ]; then
        ls $OUTPUTDIR_LHC16d/GammaCalo-All_*.root > fileLHC16d.txt
        fileNumbers=`cat fileLHC16d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16d/GammaCalo-DPGTrack_*.root > fileLHC16d.txt
        fileNumbers=`cat fileLHC16d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16d/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16d.txt
        fileNumbers=`cat fileLHC16d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
        echo "16d done"
    fi

    if [ $HAVELHC16g == 1 ]; then
        ls $OUTPUTDIR_LHC16g/GammaCalo-All_*.root > fileLHC16g.txt
        fileNumbers=`cat fileLHC16g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16g/GammaCalo-DPGTrack_*.root > fileLHC16g.txt
        fileNumbers=`cat fileLHC16g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16g/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16g.txt
        fileNumbers=`cat fileLHC16g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16h == 1 ]; then
        ls $OUTPUTDIR_LHC16h/GammaCalo-All_*.root > fileLHC16h.txt
        fileNumbers=`cat fileLHC16h.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16h/GammaCalo-DPGTrack_*.root > fileLHC16h.txt
        fileNumbers=`cat fileLHC16h.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16h/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16h.txt
        fileNumbers=`cat fileLHC16h.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16i == 1 ]; then
        ls $OUTPUTDIR_LHC16i/GammaCalo-All_*.root > fileLHC16i.txt
        fileNumbers=`cat fileLHC16i.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16i/GammaCalo-DPGTrack_*.root > fileLHC16i.txt
        fileNumbers=`cat fileLHC16i.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16i/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16i.txt
        fileNumbers=`cat fileLHC16i.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16j == 1 ]; then
        ls $OUTPUTDIR_LHC16j/GammaCalo-All_*.root > fileLHC16j.txt
        fileNumbers=`cat fileLHC16j.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16j/GammaCalo-DPGTrack_*.root > fileLHC16j.txt
        fileNumbers=`cat fileLHC16j.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16j/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16j.txt
        fileNumbers=`cat fileLHC16j.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16k == 1 ]; then
        ls $OUTPUTDIR_LHC16k/GammaCalo-All_*.root > fileLHC16k.txt
        fileNumbers=`cat fileLHC16k.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16k/GammaCalo-DPGTrack_*.root > fileLHC16k.txt
        fileNumbers=`cat fileLHC16k.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16k/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16k.txt
        fileNumbers=`cat fileLHC16k.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16k $NSlashes "LHC16k-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16l == 1 ]; then
        ls $OUTPUTDIR_LHC16l/GammaCalo-All_*.root > fileLHC16l.txt
        fileNumbers=`cat fileLHC16l.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16l/GammaCalo-DPGTrack_*.root > fileLHC16l.txt
        fileNumbers=`cat fileLHC16l.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16l/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16l.txt
        fileNumbers=`cat fileLHC16l.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16l $NSlashes "LHC16l-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16o == 1 ]; then
        ls $OUTPUTDIR_LHC16o/GammaCalo-All_*.root > fileLHC16o.txt
        fileNumbers=`cat fileLHC16o.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16o/GammaCalo-DPGTrack_*.root > fileLHC16o.txt
        fileNumbers=`cat fileLHC16o.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16o/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16o.txt
        fileNumbers=`cat fileLHC16o.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16p == 1 ]; then
        ls $OUTPUTDIR_LHC16p/GammaCalo-All_*.root > fileLHC16p.txt
        fileNumbers=`cat fileLHC16p.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16p/GammaCalo-DPGTrack_*.root > fileLHC16p.txt
        fileNumbers=`cat fileLHC16p.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16p/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16p.txt
        fileNumbers=`cat fileLHC16p.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC16e == 1 ]; then
        ls $OUTPUTDIR_LHC16e/GammaCalo-All_*.root > fileLHC16e.txt
        fileNumbers=`cat fileLHC16e.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-All" "-All"
        done;
        ls $OUTPUTDIR_LHC16e/GammaCalo-DPGTrack_*.root > fileLHC16e.txt
        fileNumbers=`cat fileLHC16e.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC16e/GammaCalo-DPGTrackAndCalo_*.root > fileLHC16e.txt
        fileNumbers=`cat fileLHC16e.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17f6 == 1 ]; then
        ls $OUTPUTDIR_LHC17f6/GammaCalo-All_*.root > fileLHC17f6.txt
            fileNumbers=`cat fileLHC17f6.txt`
            for fileName in $fileNumbers; do
            echo $fileName
        ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-All" "-All"
    done;
    ls $OUTPUTDIR_LHC17f6/GammaCalo-DPGTrack_*.root > fileLHC17f6.txt
    fileNumbers=`cat fileLHC17f6.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrack" "-DPGTrack"
    done;
    ls $OUTPUTDIR_LHC17f6/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17f6.txt
    fileNumbers=`cat fileLHC17f6.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrackAndCalo" "-DPGTrackAndCalo"
    done;
fi

    if [ $HAVELHC17f9 == 1 ]; then
        echo "Testen"
        ls $OUTPUTDIR_LHC17f9/GammaCalo-All_*.root > fileLHC17f9.txt
        fileNumbers=`cat fileLHC17f9.txt`
        for fileName in $fileNumbers; do
            echo fileName $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17f9/GammaCalo-DPGTrack_*.root > fileLHC17f9.txt
        fileNumbers=`cat fileLHC17f9.txt`
        for fileName in $fileNumbers; do
            echo fileName $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17f9/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17f9.txt
        fileNumbers=`cat fileLHC17f9.txt`
        for fileName in $fileNumbers; do
            echo fileName $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d17 == 1 ]; then
        ls $OUTPUTDIR_LHC17d17/GammaCalo-All_*.root > fileLHC17d17.txt
        fileNumbers=`cat fileLHC17d17.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d17/GammaCalo-DPGTrack_*.root > fileLHC17d17.txt
        fileNumbers=`cat fileLHC17d17.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d17/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d17.txt
        fileNumbers=`cat fileLHC17d17.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17f5 == 1 ]; then
        ls $OUTPUTDIR_LHC17f5/GammaCalo-All_*.root > fileLHC17f5.txt
        fileNumbers=`cat fileLHC17f5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17f5/GammaCalo-DPGTrack_*.root > fileLHC17f5.txt
        fileNumbers=`cat fileLHC17f5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17f5/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17f5.txt
        fileNumbers=`cat fileLHC17f5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d3 == 1 ]; then
        ls $OUTPUTDIR_LHC17d3/GammaCalo-All_*.root > fileLHC17d3.txt
        fileNumbers=`cat fileLHC17d3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d3/GammaCalo-DPGTrack_*.root > fileLHC17d3.txt
        fileNumbers=`cat fileLHC17d3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d3/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d3.txt
        fileNumbers=`cat fileLHC17d3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17e5 == 1 ]; then
        ls $OUTPUTDIR_LHC17e5/GammaCalo-All_*.root > fileLHC17e5.txt
        fileNumbers=`cat fileLHC17e5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17e5/GammaCalo-DPGTrack_*.root > fileLHC17e5.txt
        fileNumbers=`cat fileLHC17e5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17e5/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17e5.txt
        fileNumbers=`cat fileLHC17e5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17f6 == 1 ]; then
        ls $OUTPUTDIR_LHC17f6/GammaCalo-All_*.root > fileLHC17f6.txt
        fileNumbers=`cat fileLHC17f6.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17f6/GammaCalo-DPGTrack_*.root > fileLHC17f6.txt
        fileNumbers=`cat fileLHC17f6.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17f6/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17f6.txt
        fileNumbers=`cat fileLHC17f6.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17f9 == 1 ]; then
        ls $OUTPUTDIR_LHC17f9/GammaCalo-All_*.root > fileLHC17f9.txt
        fileNumbers=`cat fileLHC17f9.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17f9/GammaCalo-DPGTrack_*.root > fileLHC17f9.txt
        fileNumbers=`cat fileLHC17f9.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17f9/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17f9.txt
        fileNumbers=`cat fileLHC17f9.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17d17 == 1 ]; then
        ls $OUTPUTDIR_LHC17d17/GammaCalo-All_*.root > fileLHC17d17.txt
        fileNumbers=`cat fileLHC17d17.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d17/GammaCalo-DPGTrack_*.root > fileLHC17d17.txt
        fileNumbers=`cat fileLHC17d17.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d17/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d17.txt
        fileNumbers=`cat fileLHC17d17.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17f5 == 1 ]; then
        ls $OUTPUTDIR_LHC17f5/GammaCalo-All_*.root > fileLHC17f5.txt
        fileNumbers=`cat fileLHC17f5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17f5/GammaCalo-DPGTrack_*.root > fileLHC17f5.txt
        fileNumbers=`cat fileLHC17f5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17f5/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17f5.txt
        fileNumbers=`cat fileLHC17f5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17d3 == 1 ]; then
        ls $OUTPUTDIR_LHC17d3/GammaCalo-All_*.root > fileLHC17d3.txt
        fileNumbers=`cat fileLHC17d3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d3/GammaCalo-DPGTrack_*.root > fileLHC17d3.txt
        fileNumbers=`cat fileLHC17d3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d3/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d3.txt
        fileNumbers=`cat fileLHC17d3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17e5 == 1 ]; then
        ls $OUTPUTDIR_LHC17e5/GammaCalo-All_*.root > fileLHC17e5.txt
        fileNumbers=`cat fileLHC17e5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17e5/GammaCalo-DPGTrack_*.root > fileLHC17e5.txt
        fileNumbers=`cat fileLHC17e5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17d3-anchor16j-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17e5/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17e5.txt
        fileNumbers=`cat fileLHC17e5.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17d20a1 == 1 ]; then
    echo $OUTPUTDIR_LHC17d20a1
        ls $OUTPUTDIR_LHC17d20a1/GammaCalo-All_*.root > fileLHC17d20a1.txt
        fileNumbers=`cat fileLHC17d20a1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a1/GammaCalo-DPGTrack_*.root > fileLHC17d20a1.txt
        fileNumbers=`cat fileLHC17d20a1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a1/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d20a1.txt
        fileNumbers=`cat fileLHC17d20a1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1 $NSlashes "MC_LHC17d20a1-anchor16k-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d20a1Ex == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a1Ex/GammaCalo-All_*.root > fileLHC17d20a1_extra.txt
        fileNumbers=`cat fileLHC17d20a1_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a1Ex/GammaCalo-DPGTrack_*.root > fileLHC17d20a1_extra.txt
        fileNumbers=`cat fileLHC17d20a1_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a1Ex/GammaCalo-DPGTrackAndCalo*.root > fileLHC17d20a1_extra.txt
        fileNumbers=`cat fileLHC17d20a1_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a1Ex $NSlashes "MC_LHC17d20a1_extra-anchor16k-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17d20a2 == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a2/GammaCalo-All_*.root > fileLHC17d20a2.txt
        fileNumbers=`cat fileLHC17d20a2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a2/GammaCalo-DPGTrack_*.root > fileLHC17d20a2.txt
        fileNumbers=`cat fileLHC17d20a2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a2/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d20a2.txt
        fileNumbers=`cat fileLHC17d20a2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2 $NSlashes "MC_LHC17d20a2-anchor16l-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d20a2Ex == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a2Ex/GammaCalo-All_*.root > fileLHC17d20a2_extra.txt
        fileNumbers=`cat fileLHC17d20a2_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2Ex $NSlashes "MC_LHC17d20a2_extra-anchor16l-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d20a2Ex/GammaCalo-DPGTrack_*.root > fileLHC17d20a2_extra.txt
        fileNumbers=`cat fileLHC17d20a2_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2Ex $NSlashes "MC_LHC17d20a2_extra-anchor16l-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d20a2Ex/GammaCalo-DPGTrackAndCalo*.root > fileLHC17d20a2_extra.txt
        fileNumbers=`cat fileLHC17d20a2_extra.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d20a2Ex $NSlashes "MC_LHC17d20a2_extra-anchor16l-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17d16 == 1 ]; then
        ls $OUTPUTDIR_LHC17d16/GammaCalo-All_*.root > fileLHC17d16.txt
        fileNumbers=`cat fileLHC17d16.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d16/GammaCalo-DPGTrack_*.root > fileLHC17d16.txt
        fileNumbers=`cat fileLHC17d16.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d16/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d16.txt
        fileNumbers=`cat fileLHC17d16.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi
    if [ $HAVELHC17d18 == 1 ]; then
        ls $OUTPUTDIR_LHC17d18/GammaCalo-All_*.root > fileLHC17d18.txt
        fileNumbers=`cat fileLHC17d18.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d18/GammaCalo-DPGTrack_*.root > fileLHC17d18.txt
        fileNumbers=`cat fileLHC17d18.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d18/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d18.txt
        fileNumbers=`cat fileLHC17d18.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d16 == 1 ]; then
        ls $OUTPUTDIR_LHC17d16/GammaCalo-All_*.root > fileLHC17d16.txt
        fileNumbers=`cat fileLHC17d16.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d16/GammaCalo-DPGTrack_*.root > fileLHC17d16.txt
        fileNumbers=`cat fileLHC17d16.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d16/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d16.txt
        fileNumbers=`cat fileLHC17d16.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    if [ $HAVELHC17d18 == 1 ]; then
        ls $OUTPUTDIR_LHC17d18/GammaCalo-All_*.root > fileLHC17d18.txt
        fileNumbers=`cat fileLHC17d18.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-All" "-All"
        done;
        ls $OUTPUTDIR_LHC17d18/GammaCalo-DPGTrack_*.root > fileLHC17d18.txt
        fileNumbers=`cat fileLHC17d18.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrack" "-DPGTrack"
        done;
        ls $OUTPUTDIR_LHC17d18/GammaCalo-DPGTrackAndCalo_*.root > fileLHC17d18.txt
        fileNumbers=`cat fileLHC17d18.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-DPGTrackAndCalo" "-DPGTrackAndCalo"
        done;
    fi

    echo "Download Done"

    if [ $MERGEON == 1 ]; then
        echo "Starting Merging"
        ls $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrackAndCalo\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrackAndCalo\_$number.root
            ls $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root
            ls $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root
            if [ -f $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root ]  && [ -f $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrackAndCalo\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCalo_LHC16klo-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-DPGTrackAndCalo\_$number.root $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrackAndCalo\_$number.root
            fi
        done
        ls $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrack\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrack\_$number.root
            ls $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-DPGTrack\_$number.root
            ls $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-DPGTrack\_$number.root
            if [ -f $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-DPGTrack\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-DPGTrack\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrack\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCalo_LHC16klo-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-DPGTrack\_$number.root $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-DPGTrack\_$number.root
            fi
        done
        ls $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-All\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            ls $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-All\_$number.root
            ls $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-All\_$number.root
            ls $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-All\_$number.root
            if [ -f $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-All\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-All\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-All\_$number.root ]; then
                hadd -f $OUTPUTDIR/GammaCalo_LHC16klo-pass$passNr-All\_$number.root $OUTPUTDIR/GammaCalo_LHC16k-pass$passNr-All\_$number.root $OUTPUTDIR/GammaCalo_LHC16l-pass$passNr-All\_$number.root $OUTPUTDIR/GammaCalo_LHC16o-pass$passNr-All\_$number.root
            fi
        done


    echo "Merging Done"
    fi
else
    if [ $HAVELHC16d == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16d";
        rm $OUTPUTDIR_LHC16d/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16d/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16g == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16g";
        rm $OUTPUTDIR_LHC16g/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16g/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16h == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16h";
        rm $OUTPUTDIR_LHC16h/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16h/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16i == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16i";
        rm $OUTPUTDIR_LHC16i/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16i/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16j == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16j";
        rm $OUTPUTDIR_LHC16j/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16j/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16k == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16k";
        rm $OUTPUTDIR_LHC16k/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16k/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16l == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16l";
        rm $OUTPUTDIR_LHC16l/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16l/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16o == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16o";
        rm $OUTPUTDIR_LHC16o/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16o/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16p == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16p";
        rm $OUTPUTDIR_LHC16p/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16p/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC16e == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16e";
        rm $OUTPUTDIR_LHC16e/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC16e/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC17d20a1 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a1";
        rm $OUTPUTDIR_LHC17d20a1/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC17d20a1/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC17d20a1Ex == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a1Ex";
        rm $OUTPUTDIR_LHC17d20a1Ex/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC17d20a1Ex/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC17d20a2 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a2";
        rm $OUTPUTDIR_LHC17d20a2/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC17d20a2/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC17d20a2Ex == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC17d20a2Ex";
        rm $OUTPUTDIR_LHC17d20a2Ex/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC17d20a2Ex/*/*/*GammaCalo_*.root
    fi
fi
