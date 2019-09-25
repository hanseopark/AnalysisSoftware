#! /bin/bash
source basicFunction.sh

DOWNLOADON=1
MERGEON=1
MERGEONFASTAndWOSDD=1
SINGLERUN=1
SEPARATEON=1
MERGEONSINGLEData=1
MERGEONSINGLEMC=1
CLEANUP=1
CLEANUPMAYOR=$2
number=""
FAST="_FAST"
SEPARATEONLYConv=1
# check if train configuration has actually been given
HAVELHC13b=1
HAVELHC13c=1
HAVELHC13d=1
HAVELHC13e=1
HAVELHC13f=1
HAVELHC18j51=1
HAVELHC18j52=1
HAVELHC18j53=1
HAVELHC16q=1
HAVELHC16t=1
HAVELHC16qF=1
HAVELHC16tF=1
HAVETOBUILDData=0
HAVELHC18f31=1
HAVELHC18f32=1
HAVELHC18f3F1=1
HAVELHC18f3F2=1
HAVETOBUILDLHC18f3=0
HAVELHC17f2afix=1
HAVELHC17f2afixF=1
HAVETOBUILDLHC17f2afix=0
HAVELHC17g8a=1
HAVELHC17g8aF=1
HAVETOBUILDLHC17g8a=0
HAVELHC19a41=1
HAVELHC19a42=1

# default trainconfigurations
LHC13bData="";
LHC13cData="";
LHC13dData="";
LHC13eData="";
LHC13fData="";
LHC13beData="";
LHC18j5MC="";
LHC18j5_1MC="";
LHC18j5_2MC="";
LHC18j5_3MC="";
LHC16qData="";
LHC16tData="";
LHC16qDataFast="";
LHC16tDataFast="";
LHC16qtData="";
LHC18f3MCMoth="";
LHC18f3MC1="";
LHC18f3MC2="";
LHC18f3MCFast1="";
LHC18f3MCFast2="";
LHC17f2a_fixMCMoth="";
LHC17f2a_fixMC="";
LHC17f2a_fixMCFast="";
LHC17g8aMCMoth=""
LHC17g8aMC=""
LHC17g8aMCFast=""
LHC19a4MCMother=""
LHC19a4MC1=""
LHC19a4MC2=""
passNr="1";

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorageExternal/OutputLegoTrains/pPb
elif [ $1 = "fbockExt" ]; then
    BASEDIR=/media/fbock/BackupSeagate/OutputLegoTrains/pPb
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
fi

if [ $3 = "AOD" ]; then
    baseLegoData=GA_pPb_AOD
    baseLegoMC=GA_pPb_MC_AOD
    pathDataR1b=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathDataR1c=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathDataR1d=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathDataR1e=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathDataR1f=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathMCR1=AOD214/PWGGA/GA_pPb_MC_AOD
    pathDataR2WOSDD=pass1_CENT_woSDD/AOD190/PWGGA/GA_pPb_AOD
    pathDataR2FAST=pass1_FAST/AOD190/PWGGA/GA_pPb_AOD
    pathMCR2=AOD202/PWGGA/GA_pPb_MC_AOD
    FAST="_FAST"
    WOSDD="_CENT_woSDD"
elif [ $3 = "AODSKIMMB" ]; then
    baseLegoData=GA_pPb_AOD
    baseLegoMC=GA_pPb_MC_AOD
    pathDataR1b=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathDataR1c=pass4/AOD210/PWGGA/GA_pPb_AOD
    pathDataR1d=pass4/AOD210/PWGGA/GA_pPb_AOD/547_20190603-1141/PWGGA/GA_pPb_AOD
    pathDataR1e=pass4/AOD210/PWGGA/GA_pPb_AOD/548_20190603-1142/PWGGA/GA_pPb_AOD
    pathDataR1f=pass4/AOD210/PWGGA/GA_pPb_AOD/540_20190524-1221/PWGGA/GA_pPb_AOD
    pathMCR1=AOD214/PWGGA/GA_pPb_MC_AOD
    pathDataR2WOSDD=pass1_CENT_woSDD/AOD190/PWGGA/GA_pPb_AOD
    pathDataR2FAST=pass1_FAST/AOD190/PWGGA/GA_pPb_AOD
    pathMCR2=AOD214/PWGGA/GA_pPb_MC_AOD
    FAST="_FAST"
    WOSDD="_CENT_woSDD"
elif [ $3 = "AODSKIMEMC7" ]; then
    baseLegoData=GA_pPb_AOD
    baseLegoMC=GA_pPb_MC_AOD
    pathDataR1b=pass4/AOD210/PWGGA/GA_pPb_AOD/549_20190603-1142/PWGGA/GA_pPb_AOD
    pathDataR1c=pass4/AOD210/PWGGA/GA_pPb_AOD/550_20190603-1143/PWGGA/GA_pPb_AOD
    pathDataR1d=pass4/AOD210/PWGGA/GA_pPb_AOD/551_20190603-1143/PWGGA/GA_pPb_AOD
    pathDataR1e=pass4/AOD210/PWGGA/GA_pPb_AOD/552_20190603-1143/PWGGA/GA_pPb_AOD
    pathDataR1f=pass4/AOD210/PWGGA/GA_pPb_AOD/539_20190524-1221/PWGGA/GA_pPb_AOD
    pathMCR1=AOD214/PWGGA/GA_pPb_MC_AOD
    pathDataR2WOSDD=pass1_CENT_woSDD/AOD190/PWGGA/GA_pPb_AOD
    pathDataR2FAST=pass1_FAST/AOD190/PWGGA/GA_pPb_AOD
    pathMCR2=AOD214/PWGGA/GA_pPb_MC_AOD
    FAST="_FAST"
    WOSDD="_CENT_woSDD"
elif [ $3 = "AODSKIMEGAJE" ]; then
    baseLegoData=GA_pPb_AOD
    baseLegoMC=GA_pPb_MC_AOD
    pathDataR1b=pass4/AOD210/PWGGA/GA_pPb_AOD/544_20190524-1413/PWGGA/GA_pPb_AOD
    pathDataR1c=pass4/AOD210/PWGGA/GA_pPb_AOD/543_20190524-1412/PWGGA/GA_pPb_AOD
    pathDataR1d=pass4/AOD210/PWGGA/GA_pPb_AOD/542_20190524-1412/PWGGA/GA_pPb_AOD
    pathDataR1e=pass4/AOD210/PWGGA/GA_pPb_AOD/541_20190524-1411/PWGGA/GA_pPb_AOD
    pathDataR1f=pass4/AOD210/PWGGA/GA_pPb_AOD/537_20190524-1218/PWGGA/GA_pPb_AOD
    pathMCR1=AOD214/PWGGA/GA_pPb_MC_AOD
    pathDataR2WOSDD=pass1_CENT_woSDD/AOD190/PWGGA/GA_pPb_AOD
    pathDataR2FAST=pass1_FAST/AOD190/PWGGA/GA_pPb_AOD
    pathMCR2=AOD214/PWGGA/GA_pPb_MC_AOD
    FAST="_FAST"
    WOSDD="_CENT_woSDD"
elif [ $3 = "AODSKIMPHI7" ]; then
    baseLegoData=GA_pPb_AOD
    baseLegoMC=GA_pPb_MC_AOD
    pathDataR1b=pass4/AOD210/PWGGA/GA_pPb_AOD/558_20190604-1320/PWGGA/GA_pPb_AOD
    pathDataR1c=pass4/AOD210/PWGGA/GA_pPb_AOD/559_20190604-1320/PWGGA/GA_pPb_AOD
    pathDataR1d=pass4/AOD210/PWGGA/GA_pPb_AOD/560_20190604-1320/PWGGA/GA_pPb_AOD
    pathDataR1e=pass4/AOD210/PWGGA/GA_pPb_AOD/561_20190604-1321/PWGGA/GA_pPb_AOD
    pathDataR1f=pass4/AOD210/PWGGA/GA_pPb_AOD/562_20190604-1321/PWGGA/GA_pPb_AOD
    pathMCR1=AOD214/PWGGA/GA_pPb_MC_AOD
    pathDataR2WOSDD=pass1_CENT_woSDD/AOD190/PWGGA/GA_pPb_AOD
    pathDataR2FAST=pass1_FAST/AOD190/PWGGA/GA_pPb_AOD
    pathMCR2=AOD214/PWGGA/GA_pPb_MC_AOD
    FAST="_FAST"
    WOSDD="_CENT_woSDD"
else
    baseLegoData=GA_pPb
    baseLegoMC=GA_pPb_MC
    pathDataR1b=pass4/PWGGA/GA_pPb
    pathDataR1c=pass4/PWGGA/GA_pPb
    pathDataR1d=pass4/PWGGA/GA_pPb
    pathDataR1e=pass4/PWGGA/GA_pPb
    pathDataR1f=pass4/PWGGA/GA_pPb
    pathMCR1=PWGGA/GA_pPb_MC
    pathDataR2WOSDD=pass1/PWGGA/GA_pPb
    pathDataR2FAST=pass1/PWGGA/GA_pPb
    pathMCR2=PWGGA/GA_pPb_MC
fi


# Definitition of number of slashes in your path to different depths
NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
NSlashes=`expr $NSlashesBASE + 4`
NSlashes2=`expr $NSlashes - 1`
NSlashes3=`expr $NSlashes + 1`
NSlashes4=`expr $NSlashes + 2`
echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"

if [ $3 = "AODSKIMMB" ]; then
#     TRAINDIR=20190831-EMCNonLin
#     LHC13beData="";
#     LHC13bcData="589"
#     LHC13bData="child_1"
#     LHC13cData="child_2"
#     LHC13deData="590" #skim MB
#     LHC13dData="child_1"
#     LHC13eData="child_2"
#     LHC13fData="591" #skim MB
#     TRAINDIR=20190903-EMCtriggerStat
#     LHC13beData="";
#     LHC13bcData="602"
#     LHC13bData="child_1"
#     LHC13cData="child_2"
#     LHC13deData="603" #skim MB
#     LHC13dData="child_1"
#     LHC13eData="child_2"
#     LHC13fData="604" #skim MB
    TRAINDIR=20190916-PCMEMCSys
    LHC13beData="";
#     LHC13bcData="619"
    LHC13bcData="620"
#     LHC13bcData="621"
    LHC13bData="child_1"
    LHC13cData="child_2"
#     LHC13deData="612" #skim MB
#     LHC13dData="child_1"
#     LHC13eData="child_2"
#     LHC13fData="613" #skim MB
#     LHC13fData="617" #skim MB


elif [ $3 = "AODSKIMEMC7" ]; then
    TRAINDIR=20190831-EMCNonLin2
    # LHC13beData="600" #skim EMC7
    # LHC13bData="child_1"
    # LHC13cData="child_2"
    # LHC13dData="child_3"
    # LHC13eData="child_4"
    # LHC13fData="596" #skim EMC7
elif [ $3 = "AODSKIMEGAJE" ]; then
    # TRAINDIR=20190831-EMCNonLin
    # LHC13beData="593" #skim EGA
    # LHC13bData="child_1"
    # LHC13cData="child_2"
    # LHC13dData="child_3"
    # LHC13eData="child_4"
    # LHC13fData="594" #skim EGA

#     TRAINDIR=20190903-EMCtriggerStat
#     LHC13beData="605" #skim EGA
#     LHC13bData="child_1"
#     LHC13cData="child_2"
#     LHC13dData="child_3"
#     LHC13eData="child_4"
#     LHC13fData="606" #skim EGA

    TRAINDIR=20190916-PCMEMCSys
    LHC13beData="627" #skim EGA
    LHC13bData="child_1"
    LHC13cData="child_2"
    LHC13dData="child_3"
    LHC13eData="child_4"
    LHC13fData="628" #skim EGA
elif [ $3 = "AODSKIMPHI7" ]; then
#     TRAINDIR=20190831-EMCNonLin
#     LHC13beData="597" #skim PHI7
#     LHC13bData="child_1"
#     LHC13cData="child_2"
#     LHC13dData="child_3"
#     LHC13eData="child_4"
#     LHC13fData="598" #skim PHI7
    TRAINDIR=20190903-PHOStriggerStat
    LHC13beData="609" #skim PHI7
    LHC13bData="child_1"
    LHC13cData="child_2"
    LHC13dData="child_3"
    LHC13eData="child_4"
    LHC13fData="610" #skim PHI7

else
#     TRAINDIR=20190831-EMCNonLin
    # LHC18j5MC="730"
    # LHC18j5_1MC="child_1"
    # LHC18j5_2MC="child_2"
    # LHC18j5_3MC="child_3"
#     TRAINDIR=20190903-EMCtriggerStat
#     LHC16qtData="601";
# #     LHC16qDataFast="child_1";
#     LHC16tDataFast="child_2";
#     LHC16qData="child_3";
#     LHC16tData="child_4";
    TRAINDIR=20190916-PCMEMCSys
#     LHC16qtData="622";
#     LHC16qtData="623";
#     LHC16qtData="637";
#     LHC16qtData="638";
#     LHC16qDataFast="child_1";
#     LHC16tDataFast="child_2";
#     LHC16qData="child_3";
#     LHC16tData="child_4";
#     LHC18j5MC="744"
#     LHC18j5_1MC="child_1"
#     LHC18j5_2MC="child_2"
#     LHC18j5_3MC="child_3"
#     LHC18f3MCMoth="745";
#     LHC18f3MC1="child_2";
#     LHC18f3MC2="child_4";
#     LHC18f3MCFast1="child_1";
#     LHC18f3MCFast2="child_3";
    # LHC17f2a_fixMCMoth="1217";
    # LHC17f2a_fixMC="child_2";
    # LHC17f2a_fixMCFast="child_1";

    LHC19a4MCMother="748";
    LHC19a4MC1="child_1";
    LHC19a4MC2="child_2";
fi


OUTPUTDIR=$BASEDIR/$TRAINDIR
ALIENDIRData="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoData/"
OUTPUTDIRData=$BASEDIR/$TRAINDIR/$baseLegoData
ALIENDIRMC="/alice/cern.ch/user/a/alitrain/PWGGA/$baseLegoMC/"
OUTPUTDIRMC=$BASEDIR/$TRAINDIR/$baseLegoMC
mkdir -p $OUTPUTDIR/CutSelections
mkdir -p $OUTPUTDIR/SinglePeriods

# finding run1 data paths
if [ "$LHC13beData" == "" ]; then
  FindCorrectTrainDirectory $LHC13bData $OUTPUTDIRData $ALIENDIRData $LHC13bcData
  HAVELHC13b=$tempBool
  LHC13bData=$tempDir
  OUTPUTDIR_LHC13b=$tempPath
  echo "13b: $HAVELHC13b $LHC13bData $OUTPUTDIR_LHC13b"
  FindCorrectTrainDirectory $LHC13cData $OUTPUTDIRData $ALIENDIRData $LHC13bcData
  HAVELHC13c=$tempBool
  LHC13cData=$tempDir
  OUTPUTDIR_LHC13c=$tempPath
  echo "13c: $HAVELHC13c $LHC13cData $OUTPUTDIR_LHC13c"
  FindCorrectTrainDirectory $LHC13dData $OUTPUTDIRData $ALIENDIRData $LHC13deData
  HAVELHC13d=$tempBool
  LHC13dData=$tempDir
  OUTPUTDIR_LHC13d=$tempPath
  echo "13d: $HAVELHC13d $LHC13dData $OUTPUTDIR_LHC13d"
  FindCorrectTrainDirectory $LHC13eData $OUTPUTDIRData $ALIENDIRData $LHC13deData
  HAVELHC13e=$tempBool
  LHC13eData=$tempDir
  OUTPUTDIR_LHC13e=$tempPath
  echo "13e: $HAVELHC13e $LHC13eData $OUTPUTDIR_LHC13e"
else
  FindCorrectTrainDirectory $LHC13bData $OUTPUTDIRData $ALIENDIRData $LHC13beData
  HAVELHC13b=$tempBool
  LHC13bData=$tempDir
  OUTPUTDIR_LHC13b=$tempPath
  echo "13b: $HAVELHC13b $LHC13bData $OUTPUTDIR_LHC13b"
  FindCorrectTrainDirectory $LHC13cData $OUTPUTDIRData $ALIENDIRData $LHC13beData
  HAVELHC13c=$tempBool
  LHC13cData=$tempDir
  OUTPUTDIR_LHC13c=$tempPath
  echo "13c: $HAVELHC13c $LHC13cData $OUTPUTDIR_LHC13c"
  FindCorrectTrainDirectory $LHC13dData $OUTPUTDIRData $ALIENDIRData $LHC13beData
  HAVELHC13d=$tempBool
  LHC13dData=$tempDir
  OUTPUTDIR_LHC13d=$tempPath
  echo "13d: $HAVELHC13d $LHC13dData $OUTPUTDIR_LHC13d"
  FindCorrectTrainDirectory $LHC13eData $OUTPUTDIRData $ALIENDIRData $LHC13beData
  HAVELHC13e=$tempBool
  LHC13eData=$tempDir
  OUTPUTDIR_LHC13e=$tempPath
  echo "13e: $HAVELHC13e $LHC13eData $OUTPUTDIR_LHC13e"
fi

FindCorrectTrainDirectory $LHC13fData $OUTPUTDIRData $ALIENDIRData
HAVELHC13f=$tempBool
LHC13fData=$tempDir
OUTPUTDIR_LHC13f=$tempPath
echo "13f: $HAVELHC13f $LHC13fData $OUTPUTDIR_LHC13f"

# finding run2 data paths
FindCorrectTrainDirectory $LHC16qData $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16q=$tempBool
LHC16qData=$tempDir
OUTPUTDIR_LHC16q=$tempPath
echo "16q: $HAVELHC16q $LHC16qData $OUTPUTDIR_LHC16q"
FindCorrectTrainDirectory $LHC16tData $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16t=$tempBool
LHC16tData=$tempDir
OUTPUTDIR_LHC16t=$tempPath
echo "16t: $HAVELHC16t $LHC16tData $OUTPUTDIR_LHC16t"
FindCorrectTrainDirectory $LHC16qDataFast $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16qF=$tempBool
LHC16qDataFast=$tempDir
OUTPUTDIR_LHC16qF=$tempPath
echo "16q: $HAVELHC16qF $LHC16qDataFast $OUTPUTDIR_LHC16qF"
FindCorrectTrainDirectory $LHC16tDataFast $OUTPUTDIRData $ALIENDIRData $LHC16qtData
HAVELHC16tF=$tempBool
LHC16tDataFast=$tempDir
OUTPUTDIR_LHC16tF=$tempPath
echo "16t: $HAVELHC16tF $LHC16tDataFast $OUTPUTDIR_LHC16tF"

# finding run1 MC path
FindCorrectTrainDirectory $LHC18j5_1MC $OUTPUTDIRMC $ALIENDIRMC $LHC18j5MC
HAVELHC18j51=$tempBool
LHC18j5_1MC=$tempDir
OUTPUTDIR_LHC18j51=$tempPath
echo "18j5_1 anchored to 13bf: $HAVELHC18j51 $LHC18j5_1MC $OUTPUTDIR_LHC18j51"
FindCorrectTrainDirectory $LHC18j5_2MC $OUTPUTDIRMC $ALIENDIRMC $LHC18j5MC
HAVELHC18j52=$tempBool
LHC18j5_2MC=$tempDir
OUTPUTDIR_LHC18j52=$tempPath
echo "18j5_2 anchored to 13bf: $HAVELHC18j52 $LHC18j5_2MC $OUTPUTDIR_LHC18j52"
FindCorrectTrainDirectory $LHC18j5_3MC $OUTPUTDIRMC $ALIENDIRMC $LHC18j5MC
HAVELHC18j53=$tempBool
LHC18j5_3MC=$tempDir
OUTPUTDIR_LHC18j53=$tempPath
echo "18j5_3 anchored to 13bf: $HAVELHC18j53 $LHC18j5_3MC $OUTPUTDIR_LHC18j53"

FindCorrectTrainDirectory $LHC19a4MC1 $OUTPUTDIRMC $ALIENDIRMC $LHC19a4MCMother
HAVELHC19a41=$tempBool
LHC19a4MC1=$tempDir
OUTPUTDIR_LHC19a41=$tempPath
echo "19a4_1 JJ anchored to 13bf: $HAVELHC19a41 $LHC19a4MC1 $OUTPUTDIR_LHC19a41"
FindCorrectTrainDirectory $LHC19a4MC2 $OUTPUTDIRMC $ALIENDIRMC $LHC19a4MCMother
HAVELHC19a42=$tempBool
LHC19a4MC2=$tempDir
OUTPUTDIR_LHC19a42=$tempPath
echo "19a4_2 JJ anchored to 13bf: $HAVELHC19a42 $LHC19a4MC2 $OUTPUTDIR_LHC19a42"

# finding run2 MC path
FindCorrectTrainDirectory $LHC18f3MC1 $OUTPUTDIRMC $ALIENDIRMC $LHC18f3MCMoth
HAVELHC18f31=$tempBool
LHC18f3MC1=$tempDir
OUTPUTDIR_LHC18f31=$tempPath
echo "18f3_1 anchored to 16qt: $HAVELHC18f31 $LHC18f3MC1 $OUTPUTDIR_LHC18f31"
FindCorrectTrainDirectory $LHC18f3MC2 $OUTPUTDIRMC $ALIENDIRMC $LHC18f3MCMoth
HAVELHC18f32=$tempBool
LHC18f3MC2=$tempDir
OUTPUTDIR_LHC18f32=$tempPath
echo "18f3_2 anchored to 16qt: $HAVELHC18f32 $LHC18f3MC2 $OUTPUTDIR_LHC18f32"
FindCorrectTrainDirectory $LHC18f3MCFast1 $OUTPUTDIRMC $ALIENDIRMC $LHC18f3MCMoth
HAVELHC18f3F1=$tempBool
LHC18f3MCFast1=$tempDir
OUTPUTDIR_LHC18f3F1=$tempPath
echo "18f3_fast_1 anchored to 16qt: $HAVELHC18f3F1 $LHC18f3MCFast1 $OUTPUTDIR_LHC18f3F1"
FindCorrectTrainDirectory $LHC18f3MCFast2 $OUTPUTDIRMC $ALIENDIRMC $LHC18f3MCMoth
HAVELHC18f3F2=$tempBool
LHC18f3MCFast2=$tempDir
OUTPUTDIR_LHC18f3F2=$tempPath
echo "18f3_fast_2 anchored to 16qt: $HAVELHC18f3F2 $LHC18f3MCFast2 $OUTPUTDIR_LHC18f3F2"
FindCorrectTrainDirectory $LHC17f2a_fixMC $OUTPUTDIRMC $ALIENDIRMC $LHC17f2a_fixMCMoth
HAVELHC17f2afix=$tempBool
LHC17f2a_fixMC=$tempDir
OUTPUTDIR_LHC17f2a_fix=$tempPath
echo "17f2a_fix anchored to 16qt: $HAVELHC17f2afix $LHC17f2a_fixMC $OUTPUTDIR_LHC17f2a_fix"
FindCorrectTrainDirectory $LHC17f2a_fixMCFast $OUTPUTDIRMC $ALIENDIRMC $LHC17f2a_fixMCMoth
HAVELHC17f2afixF=$tempBool
LHC17f2a_fixMCFast=$tempDir
OUTPUTDIR_LHC17f2a_fixF=$tempPath
echo "17f2a_fix_fast anchored to 16qt: $HAVELHC17f2afixF $LHC17f2a_fixMCFast $OUTPUTDIR_LHC17f2a_fixF"
FindCorrectTrainDirectory $LHC17g8aMC $OUTPUTDIRMC $ALIENDIRMC $LHC17g8aMCMoth
HAVELHC17g8a=$tempBool
LHC17g8aMC=$tempDir
OUTPUTDIR_LHC17g8a=$tempPath
echo "17g8a anchored to 16qt: $HAVELHC17g8a $LHC17g8aMC $OUTPUTDIR_LHC17g8a"
FindCorrectTrainDirectory $LHC17g8aMCFast $OUTPUTDIRMC $ALIENDIRMC $LHC17g8aMCMoth
HAVELHC17g8aF=$tempBool
LHC17g8aMCFast=$tempDir
OUTPUTDIR_LHC17g8aF=$tempPath
echo "17g8a_fast anchored to 16qt: $HAVELHC17g8aF $LHC17g8aMCFast $OUTPUTDIR_LHC17g8aF"


currentDir=$PWD
if [ $CLEANUPMAYOR == 0 ]; then
        #echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackIncAccAndEMC" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13b" $HAVELHC13b $OUTPUTDIR_LHC13b $LHC13bData $pathDataR1b $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13c" $HAVELHC13c $OUTPUTDIR_LHC13c $LHC13cData $pathDataR1c $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13d" $HAVELHC13d $OUTPUTDIR_LHC13d $LHC13dData $pathDataR1d $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13e" $HAVELHC13e $OUTPUTDIR_LHC13e $LHC13eData $pathDataR1e $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC13f" $HAVELHC13f $OUTPUTDIR_LHC13f $LHC13fData $pathDataR1f $baseLegoData "/alice/data/2013" $NSlashes3 runlistsToMerge.txt "pass4" GammaConvCalo

    cd $currentDir
    echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j5_1" $HAVELHC18j51 $OUTPUTDIR_LHC18j51 $LHC18j5_1MC $pathMCR1 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j5_2" $HAVELHC18j52 $OUTPUTDIR_LHC18j52 $LHC18j5_2MC $pathMCR1 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j5_3" $HAVELHC18j53 $OUTPUTDIR_LHC18j53 $LHC18j5_3MC $pathMCR1 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC19a4_1" $HAVELHC19a41 $OUTPUTDIR_LHC19a41 $LHC19a4MC1 $pathMCR1 $baseLegoMC "/alice/sim/2019" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC19a4_2" $HAVELHC19a42 $OUTPUTDIR_LHC19a42 $LHC19a4MC2 $pathMCR1 $baseLegoMC "/alice/sim/2019" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    
    cd $currentDir
    echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
    echo "LHC16q" $HAVELHC16q $OUTPUTDIR_LHC16q $LHC16qData $pathDataR2WOSDD $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$WOSDD" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16q" $HAVELHC16q $OUTPUTDIR_LHC16q $LHC16qData $pathDataR2WOSDD $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$WOSDD" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16q" $HAVELHC16qF $OUTPUTDIR_LHC16qF $LHC16qDataFast $pathDataR2FAST $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$FAST" GammaConvCalo
    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16t" $HAVELHC16t $OUTPUTDIR_LHC16t $LHC16tData $pathDataR2WOSDD $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$WOSDD" GammaConvCalo
    cd $currentDir
    echo "LHC16t" $HAVELHC16tF $OUTPUTDIR_LHC16tF $LHC16tDataFast $pathDataR2FAST $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$FAST" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16t" $HAVELHC16tF $OUTPUTDIR_LHC16tF $LHC16tDataFast $pathDataR2FAST $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1$FAST" GammaConvCalo


    cd $currentDir
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18f3_cent_woSDD_1" $HAVELHC18f31 $OUTPUTDIR_LHC18f31 $LHC18f3MC1 $pathMCR2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18f3_cent_woSDD_2" $HAVELHC18f32 $OUTPUTDIR_LHC18f32 $LHC18f3MC2 $pathMCR2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC18f3$5/$runNumber/PWGGA/GA_pPb_MC/$LHC18f3MC" $NSlashes3 "/alice/sim/2017/LHC18f3$5/$runNumber/PWGGA/GA_pPb_MC/$LHC18f3MC/"
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18f3_fast_1" $HAVELHC18f3F1 $OUTPUTDIR_LHC18f3F1 $LHC18f3MCFast1 $pathMCR2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18f3_fast_2" $HAVELHC18f3F2 $OUTPUTDIR_LHC18f3F2 $LHC18f3MCFast2 $pathMCR2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC18f3_fast/$runNumber/PWGGA/GA_pPb_MC/$LHC18f3MCFast" $NSlashes3 "/alice/sim/2017/LHC18f3_fast/$runNumber/PWGGA/GA_pPb_MC/$LHC18f3MCFast/"
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f2a_fix" $HAVELHC17f2afix $OUTPUTDIR_LHC17f2a_fix $LHC17f2a_fixMC $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17f2a$5_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMC" $NSlashes3 "/alice/sim/2017/LHC17f2a$5_fix/$runNumber /PWGGA/GA_pPb_MC/$LHC17f2a_fixMC/"
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f2a_fix_fast" $HAVELHC17f2afixF $OUTPUTDIR_LHC17f2a_fixF $LHC17f2a_fixMCFast $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17f2a_fast_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast" $NSlashes3 "/alice/sim/2017/LHC17f2a_fast_fix/$runNumber/PWGGA/GA_pPb_MC/$LHC17f2a_fixMCFast/"
    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC17g8a" $HAVELHC17g8a $OUTPUTDIR_LHC17g8a $LHC17g8aMC $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17g8a$5/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMC" $NSlashes3 "/alice/sim/2017/LHC17g8a$5/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMC/"
    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC17g8a_fast" $HAVELHC17g8aF $OUTPUTDIR_LHC17g8aF $LHC17g8aMCFast $pathMCR2 $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    # "/alice/sim/2017/LHC17g8a_fast/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMCFast" $NSlashes3 "/alice/sim/2017/LHC17g8a_fast/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC17g8aMCFast/"

    listsToMerge=`cat runlistsToMerge.txt`
    for runListName in $listsToMerge; do
        if [ $HAVELHC13b == 1 ]; then
            ls $OUTPUTDIR_LHC13b/GammaConvCalo-$runListName\_*.root > fileLHC13b.txt
            fileNumbers=`cat fileLHC13b.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13c == 1 ]; then
            ls $OUTPUTDIR_LHC13c/GammaConvCalo-$runListName\_*.root > fileLHC13c.txt
            fileNumbers=`cat fileLHC13c.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13d == 1 ]; then
            ls $OUTPUTDIR_LHC13d/GammaConvCalo-$runListName\_*.root > fileLHC13d.txt
            fileNumbers=`cat fileLHC13d.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13d $NSlashes "LHC13d-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13e == 1 ]; then
            ls $OUTPUTDIR_LHC13e/GammaConvCalo-$runListName\_*.root > fileLHC13e.txt
            fileNumbers=`cat fileLHC13e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13e $NSlashes "LHC13e-pass4-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC13f == 1 ]; then
            ls $OUTPUTDIR_LHC13f/GammaConvCalo-$runListName\_*.root > fileLHC13f.txt
            fileNumbers=`cat fileLHC13f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13f $NSlashes "LHC13f-pass4-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC18j51 == 1 ]; then
            ls $OUTPUTDIR_LHC18j51/GammaConvCalo-$runListName\_*.root > fileLHC18j5.txt
            fileNumbers=`cat fileLHC18j5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j51 $NSlashes "MC_LHC18j5_1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j52 == 1 ]; then
            ls $OUTPUTDIR_LHC18j52/GammaConvCalo-$runListName\_*.root > fileLHC18j5.txt
            fileNumbers=`cat fileLHC18j5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j52 $NSlashes "MC_LHC18j5_2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j53 == 1 ]; then
            ls $OUTPUTDIR_LHC18j53/GammaConvCalo-$runListName\_*.root > fileLHC18j5.txt
            fileNumbers=`cat fileLHC18j5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j53 $NSlashes "MC_LHC18j5_3-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC16q == 1 ]; then
            ls $OUTPUTDIR_LHC16q/GammaConvCalo-$runListName\_*.root > fileLHC16q.txt
            fileNumbers=`cat fileLHC16q.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16q $NSlashes "LHC16q_woSDD-pass1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16qF == 1 ]; then
            ls $OUTPUTDIR_LHC16qF/GammaConvCalo-$runListName\_*.root > fileLHC16qF.txt
            fileNumbers=`cat fileLHC16qF.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16qF $NSlashes "LHC16q_fast-pass1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16t == 1 ]; then
            ls $OUTPUTDIR_LHC16t/GammaConvCalo-$runListName\_*.root > fileLHC16t.txt
            fileNumbers=`cat fileLHC16t.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16t $NSlashes "LHC16t_woSDD-pass1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16tF == 1 ]; then
            ls $OUTPUTDIR_LHC16tF/GammaConvCalo-$runListName\_*.root > fileLHC16tF.txt
            fileNumbers=`cat fileLHC16tF.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16tF $NSlashes "LHC16t_fast-pass1-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC18f31 == 1 ]; then
            ls $OUTPUTDIR_LHC18f31/GammaConvCalo-$runListName\_*.root > fileLHC18f3.txt
            fileNumbers=`cat fileLHC18f3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f31 $NSlashes "MC_LHC18f3_1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18f32 == 1 ]; then
            ls $OUTPUTDIR_LHC18f32/GammaConvCalo-$runListName\_*.root > fileLHC18f3.txt
            fileNumbers=`cat fileLHC18f3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f32 $NSlashes "MC_LHC18f3_2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18f3F1 == 1 ]; then
            ls $OUTPUTDIR_LHC18f3F1/GammaConvCalo-$runListName\_*.root > fileLHC18f3F1.txt
            fileNumbers=`cat fileLHC18f3F1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f3F1 $NSlashes "MC_LHC18f3_fast_1-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18f3F2 == 1 ]; then
            ls $OUTPUTDIR_LHC18f3F2/GammaConvCalo-$runListName\_*.root > fileLHC18f3F2.txt
            fileNumbers=`cat fileLHC18f3F2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f3F2 $NSlashes "MC_LHC18f3_fast_2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17f2afix == 1 ]; then
            ls $OUTPUTDIR_LHC17f2a_fix/GammaConvCalo-$runListName\_*.root > fileLHC17f2a_fix.txt
            fileNumbers=`cat fileLHC17f2a_fix.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f2a_fix $NSlashes "MC_LHC17f2a_fix-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17f2afixF == 1 ]; then
            ls $OUTPUTDIR_LHC17f2a_fixF/GammaConvCalo-$runListName\_*.root > fileLHC17f2a_fixF.txt
            fileNumbers=`cat fileLHC17f2a_fixF.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f2a_fixF $NSlashes "MC_LHC17f2a_fix_fast-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17g8a == 1 ]; then
            ls $OUTPUTDIR_LHC17g8a/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
            fileNumbers=`cat fileLHC17g8a.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17g8a $NSlashes "MC_LHC17g8a-$runListName" "-$runListName"
            done;
#             for binNumber in $binNumbersJJ; do
#                 echo $binNumber
#                 ls $OUTPUTDIR_LHC17g8a/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
#                 fileNumbers=`cat fileLHC17g8a.txt`
#                 for fileName in $fileNumbers; do
#                     echo $fileName
#                     GetFileNumberMerging $fileName $((NSlashes)) 2
#                     cp $OUTPUTDIR_LHC17g8a/$binNumber/GammaConvCalo-$runListName\_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber\_$3-$runListName\_$number.root
#                 done
#             done;
        fi
        if [ $HAVELHC17g8aF == 1 ]; then
            ls $OUTPUTDIR_LHC17g8aF/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
            fileNumbers=`cat fileLHC17g8a.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17g8aF $NSlashes "MC_LHC17g8a_fast-$runListName" "-$runListName"
            done;
#             for binNumber in $binNumbersJJ; do
#                 echo $binNumber
#                 ls $OUTPUTDIR_LHC17g8aF/GammaConvCalo-$runListName\_*.root > fileLHC17g8a.txt
#                 fileNumbers=`cat fileLHC17g8a.txt`
#                 for fileName in $fileNumbers; do
#                     echo $fileName
#                     GetFileNumberMerging $fileName $((NSlashes)) 2
#                     cp $OUTPUTDIR_LHC17g8aF/$binNumber/GammaConvCalo-$runListName\_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a_fast-$binNumber\_$3-$runListName\_$number.root
#                 done
#             done;
        fi
    done


    if [ $MERGEON == 1 ]; then
        echo "Starting Merging"

        echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
        listsToMerge=`cat runlistsToMerge.txt`
        for runListName in $listsToMerge; do
            ls $OUTPUTDIR/GammaConvCalo_LHC13f-pass4-$runListName\_*.root > filesForMerging.txt
            filesForMerging=`cat filesForMerging.txt`
            #13er daten
            periodList=`echo -e "b-pass4\nc-pass4\nd-pass4\ne-pass4\nf-pass4"`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                echo $number
                nameOut=""
                rm listCurrMerge.txt
                echo $fileName
                for periodID in $periodList; do
                    echo $periodID
                    currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-$runListName\_$number.root
                    if [ -f $currFile ]; then
                        outAdd=`echo $periodID  | cut -d "-" -f 1 `
                        nameOut+=$outAdd
                        echo -e "$currFile\n" >> listCurrMerge.txt
                    else
                        echo $currFile " does not exist"
                    fi
                done
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass4-$runListName\_$number.root
                for periodID in $periodList; do
                    mv $OUTPUTDIR/GammaConvCalo_LHC13$periodID-$runListName\_$number.root $OUTPUTDIR/SinglePeriods/
                done
            done

#             ls $OUTPUTDIR/GammaConvCalo_LHC13e-pass4-$runListName\_*.root > filesForMerging.txt
#             periodList=`echo -e "b-pass4\nc-pass4\nd-pass4\ne-pass4"`
#             for fileName in $filesForMerging; do
#                 echo $fileName
#                 GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
#                 echo $number
#                 nameOut=""
#                 rm listCurrMerge.txt
#                 echo $fileName
#                 for periodID in $periodList; do
#                     echo $periodID
#                     currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-$runListName\_$number.root
#                     if [ -f $currFile ]; then
#                         outAdd=`echo $periodID  | cut -d "-" -f 1 `
#                         nameOut+=$outAdd
#                         echo -e "$currFile\n" >> listCurrMerge.txt
#                     else
#                         echo $currFile " does not exist"
#                     fi
#                 done
#                 MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass4-$runListName\_$number.root
#             done

            ls $OUTPUTDIR/GammaConvCalo_MC_LHC18j5_1-$runListName\_*.root > filesForMerging.txt
            filesForMerging=`cat filesForMerging.txt`
            periodList=`echo -e "j5_1\nj5_2\nj5_3"`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 5 "bla" 1
                echo $number
                nameOut=""
                rm listCurrMerge.txt
                echo $fileName
                for periodID in $periodList; do
                    echo $periodID
                    currFile=$OUTPUTDIR/GammaConvCalo_MC_LHC18$periodID-$runListName\_$number.root
                    if [ -f $currFile ]; then
                        outAdd=`echo $periodID  | cut -d "-" -f 1 `
                        nameOut+=$outAdd
                        echo -e "$currFile\n" >> listCurrMerge.txt
                    else
                        echo $currFile " does not exist"
                    fi
                done
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC18j5x-$runListName\_$number.root
                for periodID in $periodList; do
                    mv $OUTPUTDIR/GammaConvCalo_MC_LHC18$periodID-$runListName\_$number.root $OUTPUTDIR/SinglePeriods/
                done
            done
        done

        echo -e "_woSDD\n_fast" > listReconstruction.txt
        listReconstruction=`cat listReconstruction.txt`
        for reco in $listReconstruction; do
            ls $OUTPUTDIR/GammaConvCalo_LHC16q$reco-pass$passNr-DPGTrackIncAccAndEMC\_*.root > filesForMerging.txt
            echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
            filesForMerging=`cat filesForMerging.txt`
            listsToMerge=`cat runlistsToMerge.txt`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 4
                echo $number
                for runListName in $listsToMerge; do
                    rm listCurrMerge.txt
                    fileQ="$OUTPUTDIR/GammaConvCalo_LHC16q$reco-pass$passNr-$runListName""_$number.root"
                    fileT="$OUTPUTDIR/GammaConvCalo_LHC16t$reco-pass$passNr-$runListName""_$number.root"
                    echo -e "$fileQ\n$fileT" > listCurrMerge.txt
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC16qt$reco-pass$passNr-$runListName\_$number.root
                    mv $fileQ $OUTPUTDIR/SinglePeriods/
                    mv $fileT $OUTPUTDIR/SinglePeriods/
                done
            done
        done
    fi

    if [ $MERGEONFASTAndWOSDD == 1 ]; then
        ls $OUTPUTDIR/GammaConvCalo_LHC16qt_fast-pass$passNr-DPGTrackIncAccAndEMC\_*.root | grep -v "WTree" > filesForMerging.txt
        echo -e "DPGTrackIncAccAndEMC" > runlistsToMerge.txt
        filesForMerging=`cat filesForMerging.txt`
        listsToMerge=`cat runlistsToMerge.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_LHC16qt_fast-pass$passNr-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_LHC16qt_woSDD-pass$passNr-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC16qt_fast-woSDD-pass$passNr-$runListName\_$number.root
            done
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_LHC16qt_fast-pass$passNr-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_LHC16qt_woSDD-pass$passNr-$runListName""_$number.root"
                fileR1="$OUTPUTDIR/GammaConvCalo_LHC13bcdef-pass4-$runListName""_$number.root"
                echo -e "$fileF\n$fileW\n$fileR1" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13bcdef-pass4_LHC16qt_fast-woSDD-pass$passNr-$runListName\_$number.root
            done
        done
        

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-DPGTrackIncAccAndEMC\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 6
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-woSDD-$runListName\_$number.root
                mv $fileW $OUTPUTDIR/SinglePeriods/
                mv $fileF $OUTPUTDIR/SinglePeriods/
            done
        done

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC18f3_fast_1-DPGTrackIncAccAndEMC\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 6
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF1="$OUTPUTDIR/GammaConvCalo_MC_LHC18f3_fast_1-$runListName""_$number.root"
                fileF2="$OUTPUTDIR/GammaConvCalo_MC_LHC18f3_fast_2-$runListName""_$number.root"
                fileW1="$OUTPUTDIR/GammaConvCalo_MC_LHC18f3_1-$runListName""_$number.root"
                fileW2="$OUTPUTDIR/GammaConvCalo_MC_LHC18f3_2-$runListName""_$number.root"
                echo -e "$fileF1\n$fileF2\n$fileW1\n$fileW2" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC18f3x-$runListName\_$number.root
                mv $fileF1 $OUTPUTDIR/SinglePeriods/
                mv $fileF2 $OUTPUTDIR/SinglePeriods/
                mv $fileW1 $OUTPUTDIR/SinglePeriods/
                mv $fileW2 $OUTPUTDIR/SinglePeriods/
            done
        done

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC18f3x-DPGTrackIncAccAndEMC\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF1="$OUTPUTDIR/GammaConvCalo_MC_LHC18f3x-$runListName""_$number.root"
                fileW2="$OUTPUTDIR/GammaConvCalo_MC_LHC18j5x-$runListName""_$number.root"
                echo -e "$fileF1\n$fileW2" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC18j5x_LHC18f3x-$runListName\_$number.root
            done
        done

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC18f3x-DPGTrackIncAccAndEMC\_*.root | grep -v "WTree" > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_MC_LHC18f3_fast-woSDD-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_fast-woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17f2a_fix_LHC18f3_fast-woSDD-$runListName\_$number.root
            done
        done


        ls $OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_fast-DPGTrackIncAccAndEMC\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 5
            echo $number
            for runListName in $listsToMerge; do
                rm listCurrMerge.txt
                fileF="$OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_fast-$runListName""_$number.root"
                fileW="$OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_woSDD-$runListName""_$number.root"
                echo -e "$fileF\n$fileW" > listCurrMerge.txt
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC17g8a_fast-woSDD-$runListName\_$number.root
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    rm listCurrMerge.txt
                    fileF="$OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber""_fast-$runListName""_$number.root"
                    fileW="$OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber""_woSDD-$runListName""_$number.root"
                    echo -e "$fileF\n$fileW" > listCurrMerge.txt
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC17g8a-$binNumber\_fast-woSDD-$runListName\_$number.root
                done
            done
        done
    fi
else
    if [ $HAVELHC16q == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16q";
#         rm $OUTPUTDIR_LHC16q/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16q/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC16t == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16t";
#         rm $OUTPUTDIR_LHC16t/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16t/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC16qF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16q";
#         rm $OUTPUTDIR_LHC16q/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16qF/*/*/*GammaConvCalo_*.root
    fi
    if [ $HAVELHC16tF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC16t";
#         rm $OUTPUTDIR_LHC16t/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16tF/*/*/*GammaConvCalo_*.root
    fi

    if [ $HAVELHC18f3 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC18f3";
#         rm $OUTPUTDIR_LHC18f3/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC18f3/*/Stage*
    fi
    if [ $HAVELHC17f2afix == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17f2a_fix";
#         rm $OUTPUTDIR_LHC17f2a_fix/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17f2a_fix/*/Stage*
    fi
    if [ $HAVELHC18f3F1 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC18f3";
#         rm $OUTPUTDIR_LHC18f3/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC18f3F1/*/Stage*
    fi
    if [ $HAVELHC18f3F2 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC18f3";
#         rm $OUTPUTDIR_LHC18f3/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC18f3F2/*/Stage*
    fi
    if [ $HAVELHC17f2afixF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17f2a_fix";
#         rm $OUTPUTDIR_LHC17f2a_fix/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17f2a_fixF/*/Stage*
    fi
    if [ $HAVELHC17g8a == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17g8a";
#         rm $OUTPUTDIR_LHC17g8a/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17g8a/*/*/Stage*
    fi
    if [ $HAVELHC17g8aF == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC17g8a";
#         rm $OUTPUTDIR_LHC17g8a/*/GammaConvCalo_*.root
        rm -rf $OUTPUTDIR_LHC17g8aF/*/*/Stage*
    fi
fi
