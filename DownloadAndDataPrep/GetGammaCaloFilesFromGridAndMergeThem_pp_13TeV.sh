#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

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
# echo $PATH

# check if train configuration has actually been given
HAVELHC16d=1
HAVELHC16e=1
HAVELHC16f=1
HAVELHC16g=1
HAVELHC16h=1
HAVELHC16i=1
HAVELHC16j=1
HAVELHC16k=1
HAVELHC16l=1
HAVELHC16o=1
HAVELHC16p=1

HAVELHC17c=1
HAVELHC17e=1
HAVELHC17f=1
HAVELHC17h=1
HAVELHC17i=1
HAVELHC17j=1
HAVELHC17k=1
HAVELHC17l=1
HAVELHC17m=1
HAVELHC17o=1
HAVELHC17r=1

HAVELHC18b=1
HAVELHC18d=1
HAVELHC18e=1
HAVELHC18f=1
HAVELHC18g=1
HAVELHC18h=1
HAVELHC18i=1
HAVELHC18j=1
HAVELHC18k=1
HAVELHC18l=1
HAVELHC18m=1
HAVELHC18n=1
HAVELHC18o=1
HAVELHC18p=1

HAVELHC17f6=1;
HAVELHC17f9=1;
HAVELHC17d17=1;
HAVELHC17f5=1;
HAVELHC17d3=1;
HAVELHC17e5=1;
HAVELHC18f1=1
HAVELHC18d8=1
HAVELHC17d16=1;
HAVELHC17d18=1;
HAVELHC17d1=1;

HAVELHC18g4=1
HAVELHC18g5=1
HAVELHC18g6=1
HAVELHC18h2=1
HAVELHC18h4=1
HAVELHC18j1=1
HAVELHC18j4=1
HAVELHC18k1=1
HAVELHC18k2=1
HAVELHC18k3=1

# default trainconfigurations
LHC16Data="";
    LHC16dData="";
    LHC16eData="";
    LHC16fData="";
    LHC16gData="";
    LHC16hData="";
    LHC16iData="";
    LHC16jData="";
    LHC16kData="";
    LHC16lData="";
    LHC16oData="";
    LHC16pData="";

LHC17Data="";
    LHC17cData="";
    LHC17eData="";
    LHC17fData="";
    LHC17hData="";
    LHC17iData="";
    LHC17jData="";
    LHC17kData="";
    LHC17lData="";
    LHC17mData="";
    LHC17oData="";
    LHC17rData="";

LHC18Data="";
    LHC18bData="";
    LHC18dData="";
    LHC18eData="";
    LHC18fData="";
    LHC18gData="";
    LHC18hData="";
    LHC18iData="";
    LHC18jData="";
    LHC18kData="";
    LHC18lData="";
    LHC18mData="";
    LHC18nData="";
    LHC18oData="";
    LHC18pData="";
LHC17MCEPOS="";

LHC16xMCPHY="";
    LHC17f6MC="";
    LHC17f9MC="";
    LHC17d17MC="";
    LHC17f5MC="";
    LHC17d3MC="";
    LHC17e5MC="";
    LHC18f1MC="";
    LHC18d8MC="";
    LHC17d16MC="";
    LHC17d18MC="";
    LHC17d1MC="";

LHC18xMCPHY="";
    LHC18g4MC=""
    LHC18g5MC=""
    LHC18g6MC=""
    LHC18h2MC=""
    LHC18h4MC=""
    LHC18j1MC=""
    LHC18j4MC=""
    LHC18k1MC=""
    LHC18k2MC=""
    LHC18k3MC=""

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
    pathData=pass$passNr/AOD208/PWGGA/GA_pp_AOD
    pathData3=pass$passNr\_withTRDtracking/AOD208/PWGGA/GA_pp_AOD
    pathData2=pass2/AOD208/PWGGA/GA_pp_AOD
    pathMC=AOD209/PWGGA/GA_pp_MC_AOD
else
    baseLegoData=GA_pp
    baseLegoMC=GA_pp_MC
    pathData=pass$passNr/PWGGA/GA_pp
    pathData3=pass$passNr\_withTRDtracking/PWGGA/GA_pp
    pathData2=pass2/PWGGA/GA_pp
    pathMC=PWGGA/GA_pp_MC
fi

# Definitition of number of slashes in your path to different depths
NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
NSlashes=`expr $NSlashesBASE + 4`
NSlashes2=`expr $NSlashes - 1`
NSlashes3=`expr $NSlashes + 1`
NSlashes4=`expr $NSlashes + 2`
echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"

TRAINDIR=Legotrain-vAN-QA2019 # SECOND TRAIN RUN

# LHC16 data
# LHC16Data="679"; #pass 1 SECOND TRAIN RUN
# LHC16dData="child_1"; #pass 1
# LHC16eData="child_2"; #pass 1
# LHC16fData="child_3"; #pass 1
# LHC16gData="child_4"; #pass 1
# LHC16hData="child_5"; #pass 1
# LHC16iData="child_6"; #pass 1
# LHC16jData="child_7"; #pass 1
# LHC16kData="child_8"; #pass 2
# LHC16lData="child_9"; #pass 2
# LHC16oData="child_10"; #pass 1
# LHC16pData="child_11"; #pass 1

# LHC17 data
# LHC17Data="635"; #pass 1 SECOND TRAIN RUN
# LHC17cData="child_1"; #pass 1
# LHC17eData="child_2"; #pass 1
# LHC17fData="child_3"; #pass 1
# LHC17hData="child_4"; #pass 1
# LHC17iData="child_5"; #pass 1
# LHC17jData="child_6"; #pass 1
# LHC17kData="child_7"; #pass 1
# LHC17lData="child_8"; #pass 1
# LHC17mData="child_9"; #pass 1
# LHC17oData="child_10"; #pass 1
# LHC17rData="child_11"; #pass 1


# LHC18 data
LHC18Data="705"; #pass 1 SECOND TRAIN RUN
LHC18bData="child_1";
LHC18dData="child_2";
LHC18eData="child_3";
LHC18fData="child_4";
LHC18gData="child_5";
LHC18hData="child_6";
LHC18iData="child_7";
LHC18jData="child_8";
LHC18kData="child_9";
LHC18lData="child_10";
LHC18mData="child_11";
LHC18nData="child_12";
LHC18oData="child_13";
LHC18pData="child_14";

# LHC16 MC
# LHC16xMCPHY="1276"; #pass 1
# LHC17f6MC="child_1";
# LHC17f9MC="child_2";
# LHC17d1MC="child_11";
# LHC17d17MC="child_3";
# LHC17f5MC="child_4";
# LHC17d3MC="child_5";
# LHC17e5MC="child_6";
# LHC18f1MC="child_9";
# LHC18d8MC="child_10";
# LHC17d16MC="child_7";
# LHC17d18MC="child_8";


# LHC18 MC
LHC18xMCPHY="1325"; #pass 1 FIRST TRAIN RUN
LHC18g4MC="child_1"
LHC18g5MC="child_2"
LHC18g6MC="child_3"
LHC18h2MC="child_4"
LHC18h4MC="child_5"
LHC18j1MC="child_6"
LHC18j4MC="child_7"
LHC18k1MC="child_8"
LHC18k2MC="child_9"
LHC18k3MC="child_10"

OUTPUTDIR=$BASEDIR/$TRAINDIR

ALIENDIRData="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoData/"
OUTPUTDIRData=$BASEDIR/$TRAINDIR/$baseLegoData
ALIENDIRMC="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoMC/"
OUTPUTDIRMC=$BASEDIR/$TRAINDIR/$baseLegoMC
mkdir -p $OUTPUTDIR/CutSelections

FindCorrectTrainDirectory $LHC16dData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16d=$tempBool
LHC16dData=$tempDir
OUTPUTDIR_LHC16d=$tempPath
echo "16d: $HAVELHC16d $LHC16dData $OUTPUTDIR_LHC16d"

FindCorrectTrainDirectory $LHC16eData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16e=$tempBool
LHC16eData=$tempDir
OUTPUTDIR_LHC16e=$tempPath
echo "16e: $HAVELHC16e $LHC16eData $OUTPUTDIR_LHC16e"

FindCorrectTrainDirectory $LHC16fData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16f=$tempBool
LHC16fData=$tempDir
OUTPUTDIR_LHC16f=$tempPath
echo "16f: $HAVELHC16f $LHC16fData $OUTPUTDIR_LHC16f"

FindCorrectTrainDirectory $LHC16gData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16g=$tempBool
LHC16gData=$tempDir
OUTPUTDIR_LHC16g=$tempPath
echo "16g: $HAVELHC16g $LHC16gData $OUTPUTDIR_LHC16g"

FindCorrectTrainDirectory $LHC16hData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16h=$tempBool
LHC16hData=$tempDir
OUTPUTDIR_LHC16h=$tempPath
echo "16h: $HAVELHC16h $LHC16hData $OUTPUTDIR_LHC16h"

FindCorrectTrainDirectory $LHC16iData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16i=$tempBool
LHC16iData=$tempDir
OUTPUTDIR_LHC16i=$tempPath
echo "16i: $HAVELHC16i $LHC16iData $OUTPUTDIR_LHC16i"

FindCorrectTrainDirectory $LHC16jData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16j=$tempBool
LHC16jData=$tempDir
OUTPUTDIR_LHC16j=$tempPath
echo "16j: $HAVELHC16j $LHC16jData $OUTPUTDIR_LHC16j"

FindCorrectTrainDirectory $LHC16kData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16k=$tempBool
LHC16kData=$tempDir
OUTPUTDIR_LHC16k=$tempPath
echo "16k: $HAVELHC16k $LHC16kData $OUTPUTDIR_LHC16k"

FindCorrectTrainDirectory $LHC16lData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16l=$tempBool
LHC16lData=$tempDir
OUTPUTDIR_LHC16l=$tempPath
echo "16l: $HAVELHC16l $LHC16lData $OUTPUTDIR_LHC16l"

FindCorrectTrainDirectory $LHC16oData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16o=$tempBool
LHC16oData=$tempDir
OUTPUTDIR_LHC16o=$tempPath
echo "16o: $HAVELHC16o $LHC16oData $OUTPUTDIR_LHC16o"

FindCorrectTrainDirectory $LHC16pData $OUTPUTDIRData $ALIENDIRData $LHC16Data
HAVELHC16p=$tempBool
LHC16pData=$tempDir
OUTPUTDIR_LHC16p=$tempPath
echo "16p: $HAVELHC16p $LHC16pData $OUTPUTDIR_LHC16p"

FindCorrectTrainDirectory $LHC17cData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17c=$tempBool
LHC17cData=$tempDir
OUTPUTDIR_LHC17c=$tempPath
echo "17c: $HAVELHC17c $LHC17cData $OUTPUTDIR_LHC17c"

FindCorrectTrainDirectory $LHC17eData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17e=$tempBool
LHC17eData=$tempDir
OUTPUTDIR_LHC17e=$tempPath
echo "17e: $HAVELHC17e $LHC17eData $OUTPUTDIR_LHC17e"

FindCorrectTrainDirectory $LHC17fData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17f=$tempBool
LHC17fData=$tempDir
OUTPUTDIR_LHC17f=$tempPath
echo "17f: $HAVELHC17f $LHC17fData $OUTPUTDIR_LHC17f"

FindCorrectTrainDirectory $LHC17hData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17h=$tempBool
LHC17hData=$tempDir
OUTPUTDIR_LHC17h=$tempPath
echo "17h: $HAVELHC17h $LHC17hData $OUTPUTDIR_LHC17h"

FindCorrectTrainDirectory $LHC17iData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17i=$tempBool
LHC17iData=$tempDir
OUTPUTDIR_LHC17i=$tempPath
echo "17i: $HAVELHC17i $LHC17iData $OUTPUTDIR_LHC17i"

FindCorrectTrainDirectory $LHC17jData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17j=$tempBool
LHC17jData=$tempDir
OUTPUTDIR_LHC17j=$tempPath
echo "17j: $HAVELHC17j $LHC17jData $OUTPUTDIR_LHC17j"

FindCorrectTrainDirectory $LHC17kData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17k=$tempBool
LHC17kData=$tempDir
OUTPUTDIR_LHC17k=$tempPath
echo "17k: $HAVELHC17k $LHC17kData $OUTPUTDIR_LHC17k"

FindCorrectTrainDirectory $LHC17lData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17l=$tempBool
LHC17lData=$tempDir
OUTPUTDIR_LHC17l=$tempPath
echo "17l: $HAVELHC17l $LHC17lData $OUTPUTDIR_LHC17l"

FindCorrectTrainDirectory $LHC17mData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17m=$tempBool
LHC17mData=$tempDir
OUTPUTDIR_LHC17m=$tempPath
echo "17m: $HAVELHC17m $LHC17mData $OUTPUTDIR_LHC17m"

FindCorrectTrainDirectory $LHC17oData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17o=$tempBool
LHC17oData=$tempDir
OUTPUTDIR_LHC17o=$tempPath
echo "17o: $HAVELHC17o $LHC17oData $OUTPUTDIR_LHC17o"

FindCorrectTrainDirectory $LHC17rData $OUTPUTDIRData $ALIENDIRData $LHC17Data
HAVELHC17r=$tempBool
LHC17rData=$tempDir
OUTPUTDIR_LHC17r=$tempPath
echo "17r: $HAVELHC17r $LHC17rData $OUTPUTDIR_LHC17r"


FindCorrectTrainDirectory $LHC18bData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18b=$tempBool
LHC18bData=$tempDir
OUTPUTDIR_LHC18b=$tempPath
echo "18b: $HAVELHC18b $LHC18bData $OUTPUTDIR_LHC18b"

FindCorrectTrainDirectory $LHC18dData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18d=$tempBool
LHC18dData=$tempDir
OUTPUTDIR_LHC18d=$tempPath
echo "18d: $HAVELHC18d $LHC18dData $OUTPUTDIR_LHC18d"

FindCorrectTrainDirectory $LHC18eData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18e=$tempBool
LHC18eData=$tempDir
OUTPUTDIR_LHC18e=$tempPath
echo "18e: $HAVELHC18e $LHC18eData $OUTPUTDIR_LHC18e"

FindCorrectTrainDirectory $LHC18fData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18f=$tempBool
LHC18fData=$tempDir
OUTPUTDIR_LHC18f=$tempPath
echo "18f: $HAVELHC18f $LHC18fData $OUTPUTDIR_LHC18f"

FindCorrectTrainDirectory $LHC18gData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18g=$tempBool
LHC18gData=$tempDir
OUTPUTDIR_LHC18g=$tempPath
echo "18g: $HAVELHC18g $LHC18gData $OUTPUTDIR_LHC18g"

FindCorrectTrainDirectory $LHC18hData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18h=$tempBool
LHC18hData=$tempDir
OUTPUTDIR_LHC18h=$tempPath
echo "18h: $HAVELHC18h $LHC18hData $OUTPUTDIR_LHC18h"

FindCorrectTrainDirectory $LHC18iData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18i=$tempBool
LHC18iData=$tempDir
OUTPUTDIR_LHC18i=$tempPath
echo "18i: $HAVELHC18i $LHC18iData $OUTPUTDIR_LHC18i"

FindCorrectTrainDirectory $LHC18jData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18j=$tempBool
LHC18jData=$tempDir
OUTPUTDIR_LHC18j=$tempPath
echo "18j: $HAVELHC18j $LHC18jData $OUTPUTDIR_LHC18j"

FindCorrectTrainDirectory $LHC18kData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18k=$tempBool
LHC18kData=$tempDir
OUTPUTDIR_LHC18k=$tempPath
echo "18k: $HAVELHC18k $LHC18kData $OUTPUTDIR_LHC18k"

FindCorrectTrainDirectory $LHC18lData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18l=$tempBool
LHC18lData=$tempDir
OUTPUTDIR_LHC18l=$tempPath
echo "18l: $HAVELHC18l $LHC18lData $OUTPUTDIR_LHC18l"

FindCorrectTrainDirectory $LHC18mData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18m=$tempBool
LHC18mData=$tempDir
OUTPUTDIR_LHC18m=$tempPath
echo "18m: $HAVELHC18m $LHC18mData $OUTPUTDIR_LHC18m"

FindCorrectTrainDirectory $LHC18nData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18n=$tempBool
LHC18nData=$tempDir
OUTPUTDIR_LHC18n=$tempPath
echo "18n: $HAVELHC18n $LHC18nData $OUTPUTDIR_LHC18n"

FindCorrectTrainDirectory $LHC18oData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18o=$tempBool
LHC18oData=$tempDir
OUTPUTDIR_LHC18o=$tempPath
echo "18o: $HAVELHC18o $LHC18oData $OUTPUTDIR_LHC18o"

FindCorrectTrainDirectory $LHC18pData $OUTPUTDIRData $ALIENDIRData $LHC18Data
HAVELHC18p=$tempBool
LHC18pData=$tempDir
OUTPUTDIR_LHC18p=$tempPath
echo "18p: $HAVELHC18p $LHC18pData $OUTPUTDIR_LHC18p"

# start with finding MC directories
FindCorrectTrainDirectory $LHC17f6MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f6=$tempBool
LHC17f6MC=$tempDir
OUTPUTDIR_LHC17f6=$tempPath
echo "17f6 anchored to 16d: $HAVELHC17f6 $LHC17f6MC $OUTPUTDIR_LHC17f6"

FindCorrectTrainDirectory $LHC17f9MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f9=$tempBool
LHC17f9MC=$tempDir
OUTPUTDIR_LHC17f9=$tempPath
echo "17f9 anchored to 16e: $HAVELHC17f9 $LHC17f9MC $OUTPUTDIR_LHC17f9"

FindCorrectTrainDirectory $LHC17d1MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d1=$tempBool
LHC17d1MC=$tempDir
OUTPUTDIR_LHC17d1=$tempPath
echo "17d1 anchored to 16f: $HAVELHC17d1 $LHC17d1MC $OUTPUTDIR_LHC17d1"

FindCorrectTrainDirectory $LHC17d17MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d17=$tempBool
LHC17d17MC=$tempDir
OUTPUTDIR_LHC17d17=$tempPath
echo "17d17 anchored to 16g: $HAVELHC17d17 $LHC17d17MC $OUTPUTDIR_LHC17d17"

FindCorrectTrainDirectory $LHC17f5MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f5=$tempBool
LHC17f5MC=$tempDir
OUTPUTDIR_LHC17f5=$tempPath
echo "17f5 anchored to 16h: $HAVELHC17f5 $LHC17f5MC $OUTPUTDIR_LHC17f5"

FindCorrectTrainDirectory $LHC17d3MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d3=$tempBool
LHC17d3MC=$tempDir
OUTPUTDIR_LHC17d3=$tempPath
echo "17d3 anchored to 16i: $HAVELHC17d3 $LHC17d3MC $OUTPUTDIR_LHC17d3"

FindCorrectTrainDirectory $LHC17e5MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17e5=$tempBool
LHC17e5MC=$tempDir
OUTPUTDIR_LHC17e5=$tempPath
echo "17e5 anchored to 16j: $HAVELHC17e5 $LHC17e5MC $OUTPUTDIR_LHC17e5"

FindCorrectTrainDirectory $LHC18f1MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC18f1=$tempBool
LHC18f1MC=$tempDir
OUTPUTDIR_LHC18f1=$tempPath
echo "18f1 anchored to 16k: $HAVELHC18f1 $LHC18f1MC $OUTPUTDIR_LHC18f1"

FindCorrectTrainDirectory $LHC18d8MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC18d8=$tempBool
LHC18d8MC=$tempDir
OUTPUTDIR_LHC18d8=$tempPath
echo "18d8 anchored to 16l: $HAVELHC18d8 $LHC18d8MC $OUTPUTDIR_LHC18d8"

FindCorrectTrainDirectory $LHC17d16MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d16=$tempBool
LHC17d16MC=$tempDir
OUTPUTDIR_LHC17d16=$tempPath
echo "17d16 anchored to 16o: $HAVELHC17d16 $LHC17d16MC $OUTPUTDIR_LHC17d16"

FindCorrectTrainDirectory $LHC17d18MC $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d18=$tempBool
LHC17d18MC=$tempDir
OUTPUTDIR_LHC17d18=$tempPath
echo "17d18 anchored to 16p: $HAVELHC17d18 $LHC17d18MC $OUTPUTDIR_LHC17d18"

FindCorrectTrainDirectory $LHC18g4MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18g4=$tempBool
LHC18g4MC=$tempDir
OUTPUTDIR_LHC18g4=$tempPath
echo "18g4 anchored to 18b: $HAVELHC18g4 $LHC18g4MC $OUTPUTDIR_LHC18g4"

FindCorrectTrainDirectory $LHC18g5MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18g5=$tempBool
LHC18g5MC=$tempDir
OUTPUTDIR_LHC18g5=$tempPath
echo "18g5 anchored to 18d: $HAVELHC18g5 $LHC18g5MC $OUTPUTDIR_LHC18g5"

FindCorrectTrainDirectory $LHC18g6MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18g6=$tempBool
LHC18g6MC=$tempDir
OUTPUTDIR_LHC18g6=$tempPath
echo "18g6 anchored to 18e: $HAVELHC18g6 $LHC18g6MC $OUTPUTDIR_LHC18g6"

FindCorrectTrainDirectory $LHC18h2MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18h2=$tempBool
LHC18h2MC=$tempDir
OUTPUTDIR_LHC18h2=$tempPath
echo "18h2 anchored to 18f: $HAVELHC18h2 $LHC18h2MC $OUTPUTDIR_LHC18h2"

FindCorrectTrainDirectory $LHC18h4MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18h4=$tempBool
LHC18h4MC=$tempDir
OUTPUTDIR_LHC18h4=$tempPath
echo "18h4 anchored to 18g,h,i,j,k: $HAVELHC18h4 $LHC18h4MC $OUTPUTDIR_LHC18h4"

FindCorrectTrainDirectory $LHC18j1MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18j1=$tempBool
LHC18j1MC=$tempDir
OUTPUTDIR_LHC18j1=$tempPath
echo "18j1 anchored to 18l: $HAVELHC18j1 $LHC18j1MC $OUTPUTDIR_LHC18j1"

FindCorrectTrainDirectory $LHC18j4MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18j4=$tempBool
LHC18j4MC=$tempDir
OUTPUTDIR_LHC18j4=$tempPath
echo "18j4 anchored to 18m: $HAVELHC18j4 $LHC18j4MC $OUTPUTDIR_LHC18j4"

FindCorrectTrainDirectory $LHC18k1MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18k1=$tempBool
LHC18k1MC=$tempDir
OUTPUTDIR_LHC18k1=$tempPath
echo "18k1 anchored to 18n: $HAVELHC18k1 $LHC18k1MC $OUTPUTDIR_LHC18k1"

FindCorrectTrainDirectory $LHC18k2MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18k2=$tempBool
LHC18k2MC=$tempDir
OUTPUTDIR_LHC18k2=$tempPath
echo "18k2 anchored to 18o: $HAVELHC18k2 $LHC18k2MC $OUTPUTDIR_LHC18k2"

FindCorrectTrainDirectory $LHC18k3MC $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY
HAVELHC18k3=$tempBool
LHC18k3MC=$tempDir
OUTPUTDIR_LHC18k3=$tempPath
echo "18k3 anchored to 18p: $HAVELHC18k3 $LHC18k3MC $OUTPUTDIR_LHC18k3"

# exit


if [ $CLEANUPMAYOR == 0 ]; then
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt

    CopyRunwiseAndMergeAccordingToRunlistData "LHC16d" $HAVELHC16d $OUTPUTDIR_LHC16d $LHC16dData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16e" $HAVELHC16e $OUTPUTDIR_LHC16e $LHC16eData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16f" $HAVELHC16f $OUTPUTDIR_LHC16f $LHC16fData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16g" $HAVELHC16g $OUTPUTDIR_LHC16g $LHC16gData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16h" $HAVELHC16h $OUTPUTDIR_LHC16h $LHC16hData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16i" $HAVELHC16i $OUTPUTDIR_LHC16i $LHC16iData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16j" $HAVELHC16j $OUTPUTDIR_LHC16j $LHC16jData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16k" $HAVELHC16k $OUTPUTDIR_LHC16k $LHC16kData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass2" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16l" $HAVELHC16l $OUTPUTDIR_LHC16l $LHC16lData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass2" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16o" $HAVELHC16o $OUTPUTDIR_LHC16o $LHC16oData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16p" $HAVELHC16p $OUTPUTDIR_LHC16p $LHC16pData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo

    CopyRunwiseAndMergeAccordingToRunlistData "LHC17c" $HAVELHC17c $OUTPUTDIR_LHC17c $LHC17cData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17e" $HAVELHC17e $OUTPUTDIR_LHC17e $LHC17eData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17f" $HAVELHC17f $OUTPUTDIR_LHC17f $LHC17fData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17h" $HAVELHC17h $OUTPUTDIR_LHC17h $LHC17hData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17i" $HAVELHC17i $OUTPUTDIR_LHC17i $LHC17iData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17j" $HAVELHC17j $OUTPUTDIR_LHC17j $LHC17jData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17k" $HAVELHC17k $OUTPUTDIR_LHC17k $LHC17kData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17l" $HAVELHC17l $OUTPUTDIR_LHC17l $LHC17lData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17m" $HAVELHC17m $OUTPUTDIR_LHC17m $LHC17mData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17o" $HAVELHC17o $OUTPUTDIR_LHC17o $LHC17oData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17r" $HAVELHC17r $OUTPUTDIR_LHC17r $LHC17rData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo

    CopyRunwiseAndMergeAccordingToRunlistData "LHC18b" $HAVELHC18b $OUTPUTDIR_LHC18b $LHC18bData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18d" $HAVELHC18d $OUTPUTDIR_LHC18d $LHC18dData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18e" $HAVELHC18e $OUTPUTDIR_LHC18e $LHC18eData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18f" $HAVELHC18f $OUTPUTDIR_LHC18f $LHC18fData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18g" $HAVELHC18g $OUTPUTDIR_LHC18g $LHC18gData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18h" $HAVELHC18h $OUTPUTDIR_LHC18h $LHC18hData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18i" $HAVELHC18i $OUTPUTDIR_LHC18i $LHC18iData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18j" $HAVELHC18j $OUTPUTDIR_LHC18j $LHC18jData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18k" $HAVELHC18k $OUTPUTDIR_LHC18k $LHC18kData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18l" $HAVELHC18l $OUTPUTDIR_LHC18l $LHC18lData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18m" $HAVELHC18m $OUTPUTDIR_LHC18m $LHC18mData $pathData3 $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1_withTRDtracking" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18n" $HAVELHC18n $OUTPUTDIR_LHC18n $LHC18nData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18o" $HAVELHC18o $OUTPUTDIR_LHC18o $LHC18oData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18p" $HAVELHC18p $OUTPUTDIR_LHC18p $LHC18pData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaCalo

    currentDir=$PWD
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f6" $HAVELHC17f6 $OUTPUTDIR_LHC17f6 $LHC17f6MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f9" $HAVELHC17f9 $OUTPUTDIR_LHC17f9 $LHC17f9MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d1" $HAVELHC17d1 $OUTPUTDIR_LHC17d1 $LHC17d1MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d17" $HAVELHC17d17 $OUTPUTDIR_LHC17d17 $LHC17d17MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f5" $HAVELHC17f5 $OUTPUTDIR_LHC17f5 $LHC17f5MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d3" $HAVELHC17d3 $OUTPUTDIR_LHC17d3 $LHC17d3MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17e5" $HAVELHC17e5 $OUTPUTDIR_LHC17e5 $LHC17e5MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18f1" $HAVELHC18f1 $OUTPUTDIR_LHC18f1 $LHC18f1MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d8" $HAVELHC18d8 $OUTPUTDIR_LHC18d8 $LHC18d8MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d16" $HAVELHC17d16 $OUTPUTDIR_LHC17d16 $LHC17d16MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d18" $HAVELHC17d18 $OUTPUTDIR_LHC17d18 $LHC17d18MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaCalo

    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g4" $HAVELHC18g4 $OUTPUTDIR_LHC18g4 $LHC18g4MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g5" $HAVELHC18g5 $OUTPUTDIR_LHC18g5 $LHC18g5MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g6" $HAVELHC18g6 $OUTPUTDIR_LHC18g6 $LHC18g6MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h2" $HAVELHC18h2 $OUTPUTDIR_LHC18h2 $LHC18h2MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h4" $HAVELHC18h4 $OUTPUTDIR_LHC18h4 $LHC18h4MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j1" $HAVELHC18j1 $OUTPUTDIR_LHC18j1 $LHC18j1MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j4" $HAVELHC18j4 $OUTPUTDIR_LHC18j4 $LHC18j4MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k1" $HAVELHC18k1 $OUTPUTDIR_LHC18k1 $LHC18k1MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k2" $HAVELHC18k2 $OUTPUTDIR_LHC18k2 $LHC18k2MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k3" $HAVELHC18k3 $OUTPUTDIR_LHC18k3 $LHC18k3MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaCalo

    echo "Change Structure If Needed"

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    listsToMerge=`cat runlistsToMerge.txt`
    for runListName in $listsToMerge; do
        if [ $HAVELHC16d == 1 ]; then
            ls $OUTPUTDIR_LHC16d/GammaCalo-$runListName\_*.root > fileLHC16d.txt
            fileNumbers=`cat fileLHC16d.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16e == 1 ]; then
            ls $OUTPUTDIR_LHC16e/GammaCalo-$runListName\_*.root > fileLHC16e.txt
            fileNumbers=`cat fileLHC16e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16f == 1 ]; then
            ls $OUTPUTDIR_LHC16f/GammaCalo-$runListName\_*.root > fileLHC16f.txt
            fileNumbers=`cat fileLHC16f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16f $NSlashes "LHC16f-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16g == 1 ]; then
            ls $OUTPUTDIR_LHC16g/GammaCalo-$runListName\_*.root > fileLHC16g.txt
            fileNumbers=`cat fileLHC16g.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16h == 1 ]; then
            ls $OUTPUTDIR_LHC16h/GammaCalo-$runListName\_*.root > fileLHC16h.txt
            fileNumbers=`cat fileLHC16h.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16i == 1 ]; then
            ls $OUTPUTDIR_LHC16i/GammaCalo-$runListName\_*.root > fileLHC16i.txt
            fileNumbers=`cat fileLHC16i.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16j == 1 ]; then
            ls $OUTPUTDIR_LHC16j/GammaCalo-$runListName\_*.root > fileLHC16j.txt
            fileNumbers=`cat fileLHC16j.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16k == 1 ]; then
            ls $OUTPUTDIR_LHC16k/GammaCalo-$runListName\_*.root > fileLHC16k.txt
            fileNumbers=`cat fileLHC16k.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16k-pass2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16l == 1 ]; then
            ls $OUTPUTDIR_LHC16l/GammaCalo-$runListName\_*.root > fileLHC16l.txt
            fileNumbers=`cat fileLHC16l.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16l-pass2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16o == 1 ]; then
            ls $OUTPUTDIR_LHC16o/GammaCalo-$runListName\_*.root > fileLHC16o.txt
            fileNumbers=`cat fileLHC16o.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16p == 1 ]; then
            ls $OUTPUTDIR_LHC16p/GammaCalo-$runListName\_*.root > fileLHC16p.txt
            fileNumbers=`cat fileLHC16p.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-$runListName" "-$runListName"
            done;
        fi


        if [ $HAVELHC17c == 1 ]; then
            ls $OUTPUTDIR_LHC17c/GammaCalo-$runListName\_*.root > fileLHC17c.txt
            fileNumbers=`cat fileLHC17c.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17c $NSlashes "LHC17c-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17e == 1 ]; then
            ls $OUTPUTDIR_LHC17e/GammaCalo-$runListName\_*.root > fileLHC17e.txt
            fileNumbers=`cat fileLHC17e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17e $NSlashes "LHC17e-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17f == 1 ]; then
            ls $OUTPUTDIR_LHC17f/GammaCalo-$runListName\_*.root > fileLHC17f.txt
            fileNumbers=`cat fileLHC17f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f $NSlashes "LHC17f-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17h == 1 ]; then
            ls $OUTPUTDIR_LHC17h/GammaCalo-$runListName\_*.root > fileLHC17h.txt
            fileNumbers=`cat fileLHC17h.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17h $NSlashes "LHC17h-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17i == 1 ]; then
            ls $OUTPUTDIR_LHC17i/GammaCalo-$runListName\_*.root > fileLHC17i.txt
            fileNumbers=`cat fileLHC17i.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17i $NSlashes "LHC17i-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17j == 1 ]; then
            ls $OUTPUTDIR_LHC17j/GammaCalo-$runListName\_*.root > fileLHC17j.txt
            fileNumbers=`cat fileLHC17j.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17j $NSlashes "LHC17j-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17k == 1 ]; then
            ls $OUTPUTDIR_LHC17k/GammaCalo-$runListName\_*.root > fileLHC17k.txt
            fileNumbers=`cat fileLHC17k.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17k $NSlashes "LHC17k-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17l == 1 ]; then
            ls $OUTPUTDIR_LHC17l/GammaCalo-$runListName\_*.root > fileLHC17l.txt
            fileNumbers=`cat fileLHC17l.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17l $NSlashes "LHC17l-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17m == 1 ]; then
            ls $OUTPUTDIR_LHC17m/GammaCalo-$runListName\_*.root > fileLHC17m.txt
            fileNumbers=`cat fileLHC17m.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17m $NSlashes "LHC17m-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17o == 1 ]; then
            ls $OUTPUTDIR_LHC17o/GammaCalo-$runListName\_*.root > fileLHC17o.txt
            fileNumbers=`cat fileLHC17o.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17o $NSlashes "LHC17o-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17r == 1 ]; then
            ls $OUTPUTDIR_LHC17r/GammaCalo-$runListName\_*.root > fileLHC17r.txt
            fileNumbers=`cat fileLHC17r.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17r $NSlashes "LHC17r-pass$passNr-$runListName" "-$runListName"
            done;
        fi



        if [ $HAVELHC18b == 1 ]; then
            ls $OUTPUTDIR_LHC18b/GammaCalo-$runListName\_*.root > fileLHC18b.txt
            fileNumbers=`cat fileLHC18b.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18b $NSlashes "LHC18b-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18d == 1 ]; then
            ls $OUTPUTDIR_LHC18d/GammaCalo-$runListName\_*.root > fileLHC18d.txt
            fileNumbers=`cat fileLHC18d.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18d $NSlashes "LHC18d-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18e == 1 ]; then
            ls $OUTPUTDIR_LHC18e/GammaCalo-$runListName\_*.root > fileLHC18e.txt
            fileNumbers=`cat fileLHC18e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18e $NSlashes "LHC18e-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18f == 1 ]; then
            ls $OUTPUTDIR_LHC18f/GammaCalo-$runListName\_*.root > fileLHC18f.txt
            fileNumbers=`cat fileLHC18f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18f $NSlashes "LHC18f-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18g == 1 ]; then
            ls $OUTPUTDIR_LHC18g/GammaCalo-$runListName\_*.root > fileLHC18g.txt
            fileNumbers=`cat fileLHC18g.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18g $NSlashes "LHC18g-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h == 1 ]; then
            ls $OUTPUTDIR_LHC18h/GammaCalo-$runListName\_*.root > fileLHC18h.txt
            fileNumbers=`cat fileLHC18h.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18h $NSlashes "LHC18h-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18i == 1 ]; then
            ls $OUTPUTDIR_LHC18i/GammaCalo-$runListName\_*.root > fileLHC18i.txt
            fileNumbers=`cat fileLHC18i.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18i $NSlashes "LHC18i-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j == 1 ]; then
            ls $OUTPUTDIR_LHC18j/GammaCalo-$runListName\_*.root > fileLHC18j.txt
            fileNumbers=`cat fileLHC18j.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18j $NSlashes "LHC18j-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k == 1 ]; then
            ls $OUTPUTDIR_LHC18k/GammaCalo-$runListName\_*.root > fileLHC18k.txt
            fileNumbers=`cat fileLHC18k.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18k $NSlashes "LHC18k-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18l == 1 ]; then
            ls $OUTPUTDIR_LHC18l/GammaCalo-$runListName\_*.root > fileLHC18l.txt
            fileNumbers=`cat fileLHC18l.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18l $NSlashes "LHC18l-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18m == 1 ]; then
            ls $OUTPUTDIR_LHC18m/GammaCalo-$runListName\_*.root > fileLHC18m.txt
            fileNumbers=`cat fileLHC18m.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18m $NSlashes "LHC18m-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18n == 1 ]; then
            ls $OUTPUTDIR_LHC18n/GammaCalo-$runListName\_*.root > fileLHC18n.txt
            fileNumbers=`cat fileLHC18n.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18n $NSlashes "LHC18n-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18o == 1 ]; then
            ls $OUTPUTDIR_LHC18o/GammaCalo-$runListName\_*.root > fileLHC18o.txt
            fileNumbers=`cat fileLHC18o.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18o $NSlashes "LHC18o-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18p == 1 ]; then
            ls $OUTPUTDIR_LHC18p/GammaCalo-$runListName\_*.root > fileLHC18p.txt
            fileNumbers=`cat fileLHC18p.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18p $NSlashes "LHC18p-pass$passNr-$runListName" "-$runListName"
            done;
        fi
    done

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
#         echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    listsToMerge=`cat runlistsToMerge.txt`
    for runListName in $listsToMerge; do
        # MC for LHC16d
        if [ $HAVELHC17f6 == 1 ]; then
            ls $OUTPUTDIR_LHC17f6/GammaCalo-$runListName\_*.root > fileLHC17f6.txt
            fileNumbers=`cat fileLHC17f6.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16e
        if [ $HAVELHC17f9 == 1 ]; then
            ls $OUTPUTDIR_LHC17f9/GammaCalo-$runListName\_*.root > fileLHC17f9.txt
            fileNumbers=`cat fileLHC17f9.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16f
        if [ $HAVELHC17d1 == 1 ]; then
            ls $OUTPUTDIR_LHC17d1/GammaCalo-$runListName\_*.root > fileLHC17d1.txt
            fileNumbers=`cat fileLHC17d1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d1 $NSlashes "MC_LHC17d1-anchor16f-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16g
        if [ $HAVELHC17d17 == 1 ]; then
            ls $OUTPUTDIR_LHC17d17/GammaCalo-$runListName\_*.root > fileLHC17d17.txt
            fileNumbers=`cat fileLHC17d17.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16h
        if [ $HAVELHC17f5 == 1 ]; then
            ls $OUTPUTDIR_LHC17f5/GammaCalo-$runListName\_*.root > fileLHC17f5.txt
            fileNumbers=`cat fileLHC17f5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16i
        if [ $HAVELHC17d3 == 1 ]; then
            ls $OUTPUTDIR_LHC17d3/GammaCalo-$runListName\_*.root > fileLHC17d3.txt
            fileNumbers=`cat fileLHC17d3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16j
        if [ $HAVELHC17e5 == 1 ]; then
            ls $OUTPUTDIR_LHC17e5/GammaCalo-$runListName\_*.root > fileLHC17e5.txt
            fileNumbers=`cat fileLHC17e5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16k
        if [ $HAVELHC18f1 == 1 ]; then
            ls $OUTPUTDIR_LHC18f1/GammaCalo-$runListName\_*.root > fileLHC18f1.txt
            fileNumbers=`cat fileLHC18f1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18f1 $NSlashes "MC_LHC18f1-anchor16k-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16l
        if [ $HAVELHC18d8 == 1 ]; then
            ls $OUTPUTDIR_LHC18d8/GammaCalo-$runListName\_*.root > fileLHC18d8.txt
            fileNumbers=`cat fileLHC18d8.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18d8 $NSlashes "MC_LHC18d8-anchor16l-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16o
        if [ $HAVELHC17d16 == 1 ]; then
            ls $OUTPUTDIR_LHC17d16/GammaCalo-$runListName\_*.root > fileLHC17d16.txt
            fileNumbers=`cat fileLHC17d16.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16p
        if [ $HAVELHC17d18 == 1 ]; then
            ls $OUTPUTDIR_LHC17d18/GammaCalo-$runListName\_*.root > fileLHC17d18.txt
            fileNumbers=`cat fileLHC17d18.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC18b
        if [ $HAVELHC18g4 == 1 ]; then
            ls $OUTPUTDIR_LHC18g4/GammaCalo-$runListName\_*.root > fileLHC18g4.txt
            fileNumbers=`cat fileLHC18g4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18g4 $NSlashes "MC_LHC18g4-anchor18b-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC18d
        if [ $HAVELHC18g5 == 1 ]; then
            ls $OUTPUTDIR_LHC18g5/GammaCalo-$runListName\_*.root > fileLHC18g5.txt
            fileNumbers=`cat fileLHC18g5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18g5 $NSlashes "MC_LHC18g5-anchor18d-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18g6 == 1 ]; then
            ls $OUTPUTDIR_LHC18g6/GammaCalo-$runListName\_*.root > fileLHC18g6.txt
            fileNumbers=`cat fileLHC18g6.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18g6 $NSlashes "MC_LHC18g6-anchor18e-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h2 == 1 ]; then
            ls $OUTPUTDIR_LHC18h2/GammaCalo-$runListName\_*.root > fileLHC18h2.txt
            fileNumbers=`cat fileLHC18h2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18h2 $NSlashes "MC_LHC18h2-anchor18f-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h4 == 1 ]; then
            ls $OUTPUTDIR_LHC18h4/GammaCalo-$runListName\_*.root > fileLHC18h4.txt
            fileNumbers=`cat fileLHC18h4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18h4 $NSlashes "MC_LHC18h4-anchor18ghijk-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j1 == 1 ]; then
            ls $OUTPUTDIR_LHC18j1/GammaCalo-$runListName\_*.root > fileLHC18j1.txt
            fileNumbers=`cat fileLHC18j1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18j1 $NSlashes "MC_LHC18j1-anchor18l-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j4 == 1 ]; then
            ls $OUTPUTDIR_LHC18j4/GammaCalo-$runListName\_*.root > fileLHC18j4.txt
            fileNumbers=`cat fileLHC18j4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18j4 $NSlashes "MC_LHC18j4-anchor18m-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k1 == 1 ]; then
            ls $OUTPUTDIR_LHC18k1/GammaCalo-$runListName\_*.root > fileLHC18k1.txt
            fileNumbers=`cat fileLHC18k1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18k1 $NSlashes "MC_LHC18k1-anchor18n-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k2 == 1 ]; then
            ls $OUTPUTDIR_LHC18k2/GammaCalo-$runListName\_*.root > fileLHC18k2.txt
            fileNumbers=`cat fileLHC18k2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18k2 $NSlashes "MC_LHC18k2-anchor18o-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k3 == 1 ]; then
            ls $OUTPUTDIR_LHC18k3/GammaCalo-$runListName\_*.root > fileLHC18k3.txt
            fileNumbers=`cat fileLHC18k3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC18k3 $NSlashes "MC_LHC18k3-anchor18o-$runListName" "-$runListName"
            done;
        fi
    done

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
    if [ $HAVELHC18f1 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC18f1";
        rm $OUTPUTDIR_LHC18f1/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC18f1/*/*/*GammaCalo_*.root
    fi
    if [ $HAVELHC18d8 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC18d8";
        rm $OUTPUTDIR_LHC18d8/*/GammaCalo_*.root
        rm $OUTPUTDIR_LHC18d8/*/*/*GammaCalo_*.root
    fi
fi
