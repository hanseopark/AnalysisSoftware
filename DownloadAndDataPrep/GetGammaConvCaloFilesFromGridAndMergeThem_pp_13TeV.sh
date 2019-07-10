#! /bin/bash


thisuser=`echo ${USER}`
d=`date +%Y_%m_%d`
if [[ $thisuser = "adrian" || $thisuser = "amechler" ]]
then

    rm Error.log
    #######################################################
    ######################   Daten   ######################
    #######################################################

    # ### LHC16Data="679"; #pass 1 SECOND TRAIN RUN
    # bash DownScript.sh 679 -Name_LHC16Data GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### LHC17Data="635"; #pass 1 SECOND TRAIN RUN
    # bash DownScript.sh 635 -Name_LHC17Data GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### LHC18Data="705"; #pass 1 SECOND TRAIN RUN
    # bash DownScript.sh 705 -Name_LHC18Data_old GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### LHC18Data="753"; #pass 1 SECOND TRAIN RUN
    # bash DownScript.sh 753 -Name_LHC18Data GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### Daten LHC16Data="773"
    # bash DownScript.sh 773 -Name_Data_pp13TeV_16 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### Daten LHC17Data="774"
    # bash DownScript.sh 774 -Name_Data_pp13TeV_17 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### Daten LHC18Data="753"
    # bash DownScript.sh 753 -Name_Data_pp13TeV_18 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC -RL_listDPGEDC -runwise

    ### merging 16 + 17
    # bash DownScript.sh 773 774 -Name_Data_pp13TeV_16_17_v0 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC -noDown #-runwise

    # ### merging 16 + 17 + 18
    # bash DownScript.sh 773 774 753 -Name_Data_pp13TeV_16_17_18 GA_pp_AOD ?_GammaConvCalo_ -mergechilds -childsareperiods  -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC -RL_listDPGEDC -noDown #-runwise

    # ### Daten LHC16Data_triggered="775"
    # bash DownScript.sh 775 -Name_Skim_Data_pp13TeV_16 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDCtrigger -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### Daten LHC17Data_triggered="776"
    # bash DownScript.sh 776 -Name_Skim_Data_pp13TeV_17 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDCtrigger -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### merging 16 + 17 triggered
    # bash DownScript.sh 775 776 -Name_Skim_Data_pp13TeV_16_17 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDCtrigger -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC -noDown #-runwise


##### new Trains
    ### LHC16 Data 829
    bash DownScript.sh 829 -Name_LHC16Data_v1 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -runwise -totalLog

    ### LHC17 Data 830
    bash DownScript.sh 830 -Name_LHC17Data_v1 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -runwise -totalLog

    ### LHC18 Data 734
    # bash DownScript.sh 834 -Name_LHC18Data_v1 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -runwise -totalLog

    ### merging 16 + 17
    bash DownScript.sh 829 830 -Name_Data_pp13TeV_16_17_v1 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -noDown -runwise -totalLog
    ### merging 16 + 17 + 18
    # bash DownScript.sh 829 830 834 -Name_Data_pp13TeV_16_17_18_v1 GA_pp_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -noDown -runwise -totalLog

    #######################################################
    ######################   MC   ######################
    #######################################################


    # ### PYT8_13TeV_anchLHC16_AOD209   LHC16xMCPHY="1276"; #pass 1
    # bash DownScript.sh 1276 -Name_LHC16xMCPHY GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### PYT8_13TeV_anchLHC16_AOD209_extra   LHC16xMCPHY_extra="1273"; #pass 1
    # bash DownScript.sh 1273 -Name_LHC16xMCPHY_extra GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### PYT8_13TeV_anchLHC17_AOD209_extra   LHC17xMCPHY_extra="1408";
    # bash DownScript.sh 1408 -Name_LHC17xMCPHY_extra GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### PYT8_13TeV_anchLHC17_AOD209  LHC17xJJMC="1393";
    # bash DownScript.sh 1393 -Name_LHC17xJJMC GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### PYT8_13TeV_anchLHC17_AOD209  LHC18xMCPHY="1325"; #pass 1 FIRST TRAIN RUN
    # bash DownScript.sh 1325 -Name_LHC18xMCPHY GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### PYT8_13TeV_anchLHC17_AOD209  LHC18xMCPHY="1415"; #pass 1 FIRST TRAIN RUN
    # bash DownScript.sh 1415 -Name_LHC18xMCPHY GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC #-runwise

    # ### PYT8_13TeV_anchLHC16_AOD209
    # bash DownScript.sh 1470 -Name_MC_pp13TeV_16_17 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### PYT8_13TeV_anchLHC16_AOD209_extra
    # bash DownScript.sh 1471 -Name_MC_pp13TeV_16_17 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### PYT8_13TeV_anchLHC16_AOD209_extra2
    # bash DownScript.sh 1472 -Name_MC_pp13TeV_16_17 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### PYT8_13TeV_anchLHC17_AOD209_extra
    # bash DownScript.sh 1473 -Name_MC_pp13TeV_16_17 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### PYT8_13TeV_anchLHC17_AOD209
    # bash DownScript.sh 1474 -Name_MC_pp13TeV_16_17 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

    # ### merging 16 MCs
    # bash DownScript.sh 1470 1471 1472 -Name_MC_pp13TeV_16_v0 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC
    #
    # ### merging 17 MC
    # bash DownScript.sh 1473 1474 -Name_MC_pp13TeV_17_v0 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC
    #
    # ### merging 16 + 17 MCs
    # bash DownScript.sh 1470 1471 1472 1473 1474 -Name_MC_pp13TeV_16_17_v0 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -RL_listDPGEDC -RL_listDPGCalo -RL_listDPGalo -RL_listDPCIncAccTPC -RL_listDPGIncAccTPC -RL_listDPGIncTPC -RL_listDPGEMC #-runwise

##### new Trains
    ### LHC16 Data 1546 + 1547 + 1548
    bash DownScript.sh 1546 1547 1548 -Name_LHC16MC_v1 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -runwise -totalLog

    ### LHC17 Data 1549 + 1550
    bash DownScript.sh 1549 1550 -Name_LHC17MC_v1 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -runwise -totalLog

    ### LHC18 Data 1556 1557
    # bash DownScript.sh 1556 1557 -Name_LHC18MC_v1 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods -uSR_PCMEDC -runwise -totalLog

    ### merging 16 + 17
    bash DownScript.sh 1546 1547 1548 1549 1550 -Name_MC_pp13TeV_16_17_v1 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods -uSR_PCMEDC -noDown -totalLog -runwise
    ### merging 16 + 17 + 18
    # bash DownScript.sh 1546 1547 1548 1549 1550 1556 1557 -Name_MC_pp13TeV_16_17_18_v1 GA_pp_MC_AOD ?_GammaConvCalo -mergechilds -mergetrains -childsareperiods  -uSR_PCMEDC -noDown -totalLog -runwise

    cat Error.log
    exit
fi



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

LHC16xMCPHY_extra="";
    LHC17f6MC_extra="";
    LHC17f9MC_extra="";
    LHC17d17MC_extra="";
    LHC17f5MC_extra="";
    LHC17d3MC_extra="";
    LHC17e5MC_extra="";
    LHC18f1MC_extra="";
    LHC18d8MC_extra="";
    LHC17d16MC_extra="";
    LHC17d18MC_extra="";
    LHC17d1MC_extra="";

LHC16xMCPHY_extra2="";
    LHC17f6MC_extra2="";
    LHC17f9MC_extra2="";
    LHC17d17MC_extra2="";
    LHC17f5MC_extra2="";
    LHC17d3MC_extra2="";
    LHC17e5MC_extra2="";
    LHC17d16MC_extra2="";
    LHC17d18MC_extra2="";
    LHC17d1MC_extra2="";

LHC17xMCPHY="";
    LHC18d3=""
    LHC17k4=""
    LHC17k4=""
    LHC17h11=""
    LHC18c13=""
    LHC18a8=""
    LHC17l5=""
    LHC18a9=""
    LHC18a1=""

LHC17xMCPHY_extra="";
    LHC18d3_extra=""
    LHC18c12_extra=""
    LHC17k4_extra=""
    LHC17h11_extra=""
    LHC18c13_extra=""
    LHC18a8_extra=""
    LHC17l5_extra=""
    LHC18a9_extra=""
    LHC18a1_extra=""

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
LHC18xMCPHY_extra="";
    LHC18g4MC_extra=""
    LHC18g5MC_extra=""
    LHC18g6MC_extra=""
    LHC18h2MC_extra=""
    LHC18h4MC_extra=""
    LHC18j1MC_extra=""
    LHC18j4MC_extra=""
    LHC18k1MC_extra=""
    LHC18k2MC_extra=""
    LHC18k3MC_extra=""

LHC17xJJMC="";
  LHC18f5_1="";
  LHC18f5_2="";

passNr="1";
if [ $1 = "fbock" ]; then
    BASEDIR=/media/fbock/Elements/OutputLegoTrains/pp13TeV
elif [ $1 = "hannahbossi" ]; then
    BASEDIR=/Volumes/external_memory/CERN_data/QA
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
elif [ $1 = "jlueh" ]; then
    BASEDIR=~/Daten/GridDownload
elif [ $1 = "jokonig" ]; then
    BASEDIR=/media/gustav/external_drive/Data
fi

if [ $3 = "AOD" ]; then
    baseLegoData=GA_pp_AOD
    baseLegoMC=GA_pp_MC_AOD
    pathData=pass$passNr/AOD208/PWGGA/GA_pp_AOD
    pathData3=pass$passNr\_withTRDtracking/AOD208/PWGGA/GA_pp_AOD
    pathData2=pass2/AOD208/PWGGA/GA_pp_AOD
    pathMC=AOD209/PWGGA/GA_pp_MC_AOD
    pathMC2=PWGGA/GA_pp_MC_AOD
elif [ $3 = "AOD_JJ" ]; then
    baseLegoData=GA_pp_AOD
    baseLegoMC=GA_pp_MC_AOD
    pathData=pass$passNr/AOD208/PWGGA/GA_pp_AOD
    pathData3=pass$passNr\_withTRDtracking/AOD208/PWGGA/GA_pp_AOD
    pathData2=pass2/AOD208/PWGGA/GA_pp_AOD
    pathMC=PWGGA/GA_pp_MC_AOD
    pathMC2=PWGGA/GA_pp_MC_AOD
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

# TRAINDIR=pp_13TeV_AOD_LHC17x # SECOND TRAIN RUN
# TRAINDIR=Legotrain-vAN-QA2019 # SECOND TRAIN RUN
#
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
# LHC18Data="705"; #pass 1 SECOND TRAIN RUN
# LHC18bData="child_1";
# LHC18dData="child_2";
# LHC18eData="child_3";
# LHC18fData="child_4";
# LHC18gData="child_5";
# LHC18hData="child_6";
# LHC18iData="child_7";
# LHC18jData="child_8";
# LHC18kData="child_9";
# LHC18lData="child_10";
# LHC18mData="child_11";
# LHC18nData="child_12";
# LHC18oData="child_13";
# LHC18pData="child_14";

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

# LHC16xMCPHY_extra="1276"; #pass 1
# LHC17f6MC_extra="child_1";
# LHC17f9MC_extra="child_2";
# LHC17d1MC_extra="child_11";
# LHC17d17MC_extra="child_3";
# LHC17f5MC_extra="child_4";
# LHC17d3MC_extra="child_5";
# LHC17e5MC_extra="child_6";
# LHC18f1MC_extra="child_9";
# LHC18d8MC_extra="child_10";
# LHC17d16MC_extra="child_7";
# LHC17d18MC_extra="child_8";

# LHC16xMCPHY_extra2="1276"; #pass 1
# LHC17f6MC_extra2="child_1";
# LHC17f9MC_extra2="child_2";
# LHC17d1MC_extra2="child_11";
# LHC17d17MC_extra2="child_3";
# LHC17f5MC_extra2="child_4";
# LHC17d3MC_extra2="child_5";
# LHC17e5MC_extra2="child_6";
# LHC17d16MC_extra2="child_7";
# LHC17d18MC_extra2="child_8";

# LHC17xMCPHY="";
# LHC18d3="child_2"
# LHC18c12="child_3"
# LHC18k4="child_4"
# LHC17h11="child_5"
# LHC18c13="child_6"
# LHC18a8="child_7"
# LHC17l5="child_8"
# LHC18a9="child_9"
# LHC18a1="child_10"
#
# LHC17xMCPHY_extra="1408";
# LHC18d3_extra="child_2"
# LHC18c12_extra="child_3"
# LHC18k4_extra="child_4"
# LHC17h11_extra="child_5"
# LHC18c13_extra="child_6"
# LHC18a8_extra="child_7"
# LHC17l5_extra="child_8"
# LHC18a9_extra="child_9"
# LHC18a1_extra="child_10"

# LHC17xJJMC="1393";
# LHC18f5_1="child_1"
# LHC18f5_2="child_2"

# LHC18 MC
# LHC18xMCPHY="1325"; #pass 1 FIRST TRAIN RUN
# LHC18g4MC="child_1"
# LHC18g5MC="child_2"
# LHC18g6MC="child_3"
# LHC18h2MC="child_4"
# LHC18h4MC="child_5"
# LHC18j1MC="child_6"
# LHC18j4MC="child_7"
# LHC18k1MC="child_8"
# LHC18k2MC="child_9"
# LHC18k3MC="child_10"

# TRAINDIR=Legotrain-QA20190403 # SECOND TRAIN RUN
# LHC18 data
# LHC18Data="753"; #pass 1 SECOND TRAIN RUN
# LHC18bData="child_1";
# LHC18dData="child_2";
# LHC18eData="child_3";
# LHC18fData="child_4";
# LHC18gData="child_5";
# LHC18hData="child_6";
# LHC18iData="child_7";
# LHC18jData="child_8";
# LHC18kData="child_9";
# LHC18lData="child_10";
# LHC18mData="child_11";
# LHC18nData="child_12";
# LHC18oData="child_13";
# LHC18pData="child_14";

# LHC18 MC
# LHC18xMCPHY="1415"; #pass 1 FIRST TRAIN RUN
# LHC18g4MC="child_1"
# LHC18g5MC="child_2"
# LHC18g6MC="child_3"
# LHC18h2MC="child_4"
# LHC18h4MC="child_5"
# LHC18j1MC="child_6"
# LHC18j4MC="child_7"
# LHC18k1MC="child_8"
# LHC18k2MC="child_9"
# LHC18k3MC="child_10"

TRAINDIR=Legotrain-QA20190610 # THIRD TRAIN RUN
# LHC18 data
# LHC18Data="846";
# LHC18bData="child_1";
# LHC18dData="child_2";
# LHC18eData="child_3";
# LHC18fData="child_4";
# LHC18gData="child_5";
# LHC18hData="child_6";
# LHC18iData="child_7";
# LHC18jData="child_8";
# LHC18kData="child_9";
# LHC18lData="child_10";
# LHC18mData="child_11";
# LHC18nData="child_12";
# LHC18oData="child_13";
# LHC18pData="child_14";

# LHC18 MC
LHC18xMCPHY="1574";
# LHC18g4MC="child_1"
# LHC18g5MC="child_2"
# LHC18g6MC="child_3"
# LHC18h2MC="child_4"
# LHC18h4MC="child_5"
# LHC18j1MC="child_6"
# LHC18j4MC="child_7"
LHC18k1MC="child_8"
# LHC18k2MC="child_9"
# LHC18k3MC="child_10"

LHC18xMCPHY_extra="1575";
# LHC18g4MC_extra="child_1"
# LHC18g5MC_extra="child_2"
# LHC18g6MC_extra="child_3"
# LHC18h2MC_extra="child_4"
LHC18h4MC_extra="child_5"
# LHC18j1MC_extra="child_6"
# LHC18j4MC_extra="child_7"
# LHC18k1MC_extra="child_8"
LHC18k2MC_extra="child_9"
# LHC18k3MC_extra="child_10"

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

######### 16 extra MCs

FindCorrectTrainDirectory $LHC17f6MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f6_extra=$tempBool
LHC17f6MC_extra=$tempDir
OUTPUTDIR_LHC17f6_extra=$tempPath
echo "17f6_extra anchored to 16d: $HAVELHC17f6_extra $LHC17f6MC_extra $OUTPUTDIR_LHC17f6_extra"

FindCorrectTrainDirectory $LHC17f9MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f9_extra=$tempBool
LHC17f9MC_extra=$tempDir
OUTPUTDIR_LHC17f9_extra=$tempPath
echo "17f9_extra anchored to 16e: $HAVELHC17f9_extra $LHC17f9MC_extra $OUTPUTDIR_LHC17f9_extra"

FindCorrectTrainDirectory $LHC17d1MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d1_extra=$tempBool
LHC17d1MC_extra=$tempDir
OUTPUTDIR_LHC17d1_extra=$tempPath
echo "17d1_extra anchored to 16f: $HAVELHC17d1_extra $LHC17d1MC_extra $OUTPUTDIR_LHC17d1_extra"

FindCorrectTrainDirectory $LHC17d17MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d17_extra=$tempBool
LHC17d17MC_extra=$tempDir
OUTPUTDIR_LHC17d17_extra=$tempPath
echo "17d17_extra anchored to 16g: $HAVELHC17d17_extra $LHC17d17MC_extra $OUTPUTDIR_LHC17d17_extra"

FindCorrectTrainDirectory $LHC17f5MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f5_extra=$tempBool
LHC17f5MC_extra=$tempDir
OUTPUTDIR_LHC17f5_extra=$tempPath
echo "17f5_extra anchored to 16h: $HAVELHC17f5_extra $LHC17f5MC_extra $OUTPUTDIR_LHC17f5_extra"

FindCorrectTrainDirectory $LHC17d3MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d3_extra=$tempBool
LHC17d3MC_extra=$tempDir
OUTPUTDIR_LHC17d3_extra=$tempPath
echo "17d3_extra anchored to 16i: $HAVELHC17d3_extra $LHC17d3MC_extra $OUTPUTDIR_LHC17d3_extra"

FindCorrectTrainDirectory $LHC17e5MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17e5_extra=$tempBool
LHC17e5MC_extra=$tempDir
OUTPUTDIR_LHC17e5_extra=$tempPath
echo "17e5_extra anchored to 16j: $HAVELHC17e5_extra $LHC17e5MC_extra $OUTPUTDIR_LHC17e5_extra"

FindCorrectTrainDirectory $LHC18f1MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC18f1_extra=$tempBool
LHC18f1MC_extra=$tempDir
OUTPUTDIR_LHC18f1_extra=$tempPath
echo "18f1_extra anchored to 16k: $HAVELHC18f1_extra $LHC18f1MC_extra $OUTPUTDIR_LHC18f1_extra"

FindCorrectTrainDirectory $LHC18d8MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC18d8_extra=$tempBool
LHC18d8MC_extra=$tempDir
OUTPUTDIR_LHC18d8_extra=$tempPath
echo "18d8_extra anchored to 16l: $HAVELHC18d8_extra $LHC18d8MC_extra $OUTPUTDIR_LHC18d8_extra"

FindCorrectTrainDirectory $LHC17d16MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d16_extra=$tempBool
LHC17d16MC_extra=$tempDir
OUTPUTDIR_LHC17d16_extra=$tempPath
echo "17d16_extra anchored to 16o: $HAVELHC17d16_extra $LHC17d16MC_extra $OUTPUTDIR_LHC17d16_extra"

FindCorrectTrainDirectory $LHC17d18MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d18_extra=$tempBool
LHC17d18MC_extra=$tempDir
OUTPUTDIR_LHC17d18_extra=$tempPath
echo "17d18_extra anchored to 16p: $HAVELHC17d18_extra $LHC17d18MC_extra $OUTPUTDIR_LHC17d18_extra"


############ 16 extra2 MCs


FindCorrectTrainDirectory $LHC17f6MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f6_extra2=$tempBool
LHC17f6MC_extra2=$tempDir
OUTPUTDIR_LHC17f6_extra2=$tempPath
echo "17f6_extra2 anchored to 16d: $HAVELHC17f6_extra2 $LHC17f6MC_extra2 $OUTPUTDIR_LHC17f6_extra2"

FindCorrectTrainDirectory $LHC17f9MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f9_extra2=$tempBool
LHC17f9MC_extra2=$tempDir
OUTPUTDIR_LHC17f9_extra2=$tempPath
echo "17f9_extra2 anchored to 16e: $HAVELHC17f9_extra2 $LHC17f9MC_extra2 $OUTPUTDIR_LHC17f9_extra2"

FindCorrectTrainDirectory $LHC17d1MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d1_extra2=$tempBool
LHC17d1MC_extra2=$tempDir
OUTPUTDIR_LHC17d1_extra2=$tempPath
echo "17d1_extra2 anchored to 16f: $HAVELHC17d1_extra2 $LHC17d1MC_extra2 $OUTPUTDIR_LHC17d1_extra2"

FindCorrectTrainDirectory $LHC17d17MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d17_extra2=$tempBool
LHC17d17MC_extra2=$tempDir
OUTPUTDIR_LHC17d17_extra2=$tempPath
echo "17d17_extra2 anchored to 16g: $HAVELHC17d17_extra2 $LHC17d17MC_extra2 $OUTPUTDIR_LHC17d17_extra2"

FindCorrectTrainDirectory $LHC17f5MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17f5_extra2=$tempBool
LHC17f5MC_extra2=$tempDir
OUTPUTDIR_LHC17f5_extra2=$tempPath
echo "17f5_extra2 anchored to 16h: $HAVELHC17f5_extra2 $LHC17f5MC_extra2 $OUTPUTDIR_LHC17f5_extra2"

FindCorrectTrainDirectory $LHC17d3MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d3_extra2=$tempBool
LHC17d3MC_extra2=$tempDir
OUTPUTDIR_LHC17d3_extra2=$tempPath
echo "17d3_extra2 anchored to 16i: $HAVELHC17d3_extra2 $LHC17d3MC_extra2 $OUTPUTDIR_LHC17d3_extra2"

FindCorrectTrainDirectory $LHC17e5MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17e5_extra2=$tempBool
LHC17e5MC_extra2=$tempDir
OUTPUTDIR_LHC17e5_extra2=$tempPath
echo "17e5_extra2 anchored to 16j: $HAVELHC17e5_extra2 $LHC17e5MC_extra2 $OUTPUTDIR_LHC17e5_extra2"

FindCorrectTrainDirectory $LHC17d16MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d16_extra2=$tempBool
LHC17d16MC_extra2=$tempDir
OUTPUTDIR_LHC17d16_extra2=$tempPath
echo "17d16_extra2 anchored to 16o: $HAVELHC17d16_extra2 $LHC17d16MC_extra2 $OUTPUTDIR_LHC17d16_extra2"

FindCorrectTrainDirectory $LHC17d18MC_extra2 $OUTPUTDIRMC $ALIENDIRMC $LHC16xMCPHY
HAVELHC17d18_extra2=$tempBool
LHC17d18MC_extra2=$tempDir
OUTPUTDIR_LHC17d18_extra2=$tempPath
echo "17d18_extra2 anchored to 16p: $HAVELHC17d18_extra2 $LHC17d18MC_extra2 $OUTPUTDIR_LHC17d18_extra2"




############# 17 MCs

FindCorrectTrainDirectory $LHC18d3 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC18d3=$tempBool
LHC18d3=$tempDir
OUTPUTDIR_LHC18d3=$tempPath
echo "LHC18d3 anchored to 17c,e,f: $HAVELHC18d3 $LHC18d3MC $OUTPUTDIR_LHC18d3"

FindCorrectTrainDirectory $LHC18c12 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC18c12=$tempBool
LHC18c12=$tempDir
OUTPUTDIR_LHC18c12=$tempPath
echo "LHC18c12 anchored to 17h: $HAVELHC18c12 $LHC18c12MC $OUTPUTDIR_LHC18c12"

FindCorrectTrainDirectory $LHC17k4 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC17k4=$tempBool
LHC17k4=$tempDir
OUTPUTDIR_LHC17k4=$tempPath
echo "LHC17k4 anchored to 17i: $HAVELHC17k4 $LHC17k4MC $OUTPUTDIR_LHC17k4"

FindCorrectTrainDirectory $LHC17h11 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC17h11=$tempBool
LHC17h11=$tempDir
OUTPUTDIR_LHC17h11=$tempPath
echo "LHC17h11 anchored to 17j: $HAVELHC17h11 $LHC17h11MC $OUTPUTDIR_LHC17h11"

FindCorrectTrainDirectory $LHC18c13 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC18c13=$tempBool
LHC18c13=$tempDir
OUTPUTDIR_LHC18c13=$tempPath
echo "LHC18c13 anchored to 17k: $HAVELHC18c13 $LHC18c13MC $OUTPUTDIR_LHC18c13"

FindCorrectTrainDirectory $LHC18a8 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC18a8=$tempBool
LHC18a8=$tempDir
OUTPUTDIR_LHC18a8=$tempPath
echo "LHC18a8 anchored to 17l: $HAVELHC18a8 $LHC18a8MC $OUTPUTDIR_LHC18a8"

FindCorrectTrainDirectory $LHC17l5 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC17l5=$tempBool
LHC17l5=$tempDir
OUTPUTDIR_LHC17l5=$tempPath
echo "LHC17l5 anchored to 17m: $HAVELHC17l5 $LHC17l5MC $OUTPUTDIR_LHC17l5"

FindCorrectTrainDirectory $LHC18a9 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC18a9=$tempBool
LHC18a9=$tempDir
OUTPUTDIR_LHC18a9=$tempPath
echo "LHC18a9 anchored to 17o: $HAVELHC18a9 $LHC18a9MC $OUTPUTDIR_LHC18a9"

FindCorrectTrainDirectory $LHC18a1 $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC18a1=$tempBool
LHC18a1=$tempDir
OUTPUTDIR_LHC18a1=$tempPath
echo "LHC18a1 anchored to 17r: $HAVELHC18a1 $LHC18a1MC $OUTPUTDIR_LHC18a1"


##### extra MCs

FindCorrectTrainDirectory $LHC18d3_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC18d3_extra=$tempBool
LHC18d3_extra=$tempDir
OUTPUTDIR_LHC18d3_extra=$tempPath
echo "LHC18d3_extra anchored to 17c,e,f: $HAVELHC18d3_extra $LHC18d3MC_extra $OUTPUTDIR_LHC18d3_extra"

FindCorrectTrainDirectory $LHC18c12_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC18c12_extra=$tempBool
LHC18c12_extra=$tempDir
OUTPUTDIR_LHC18c12_extra=$tempPath
echo "LHC18c12_extra anchored to 17h: $HAVELHC18c12_extra $LHC18c12MC_extra $OUTPUTDIR_LHC18c12_extra"

FindCorrectTrainDirectory $LHC17k4_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC17k4_extra=$tempBool
LHC17k4_extra=$tempDir
OUTPUTDIR_LHC17k4_extra=$tempPath
echo "LHC17k4_extra anchored to 17i: $HAVELHC17k4_extra $LHC17k4MC_extra $OUTPUTDIR_LHC17k4_extra"

FindCorrectTrainDirectory $LHC17h11_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC17h11_extra=$tempBool
LHC17h11_extra=$tempDir
OUTPUTDIR_LHC17h11_extra=$tempPath
echo "LHC17h11_extra anchored to 17j: $HAVELHC17h11_extra $LHC17h11MC_extra $OUTPUTDIR_LHC17h11_extra"

FindCorrectTrainDirectory $LHC18c13_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC18c13_extra=$tempBool
LHC18c13_extra=$tempDir
OUTPUTDIR_LHC18c13_extra=$tempPath
echo "LHC18c13_extra anchored to 17k: $HAVELHC18c13_extra $LHC18c13MC_extra $OUTPUTDIR_LHC18c13_extra"

FindCorrectTrainDirectory $LHC18a8_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC18a8_extra=$tempBool
LHC18a8_extra=$tempDir
OUTPUTDIR_LHC18a8_extra=$tempPath
echo "LHC18a8_extra anchored to 17l: $HAVELHC18a8_extra $LHC18a8MC_extra $OUTPUTDIR_LHC18a8_extra"

FindCorrectTrainDirectory $LHC17l5_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY
HAVELHC17l5_extra=$tempBool
LHC17l5_extra=$tempDir
OUTPUTDIR_LHC17l5_extra=$tempPath
echo "LHC17l5_extra anchored to 17m: $HAVELHC17l5_extra $LHC17l5MC_extra $OUTPUTDIR_LHC17l5_extra"

FindCorrectTrainDirectory $LHC18a9_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC18a9_extra=$tempBool
LHC18a9_extra=$tempDir
OUTPUTDIR_LHC18a9_extra=$tempPath
echo "LHC18a9_extra anchored to 17o: $HAVELHC18a9_extra $LHC18a9MC_extra $OUTPUTDIR_LHC18a9_extra"

FindCorrectTrainDirectory $LHC18a1_extra $OUTPUTDIRMC $ALIENDIRMC $LHC17xMCPHY_extra
HAVELHC18a1_extra=$tempBool
LHC18a1_extra=$tempDir
OUTPUTDIR_LHC18a1_extra=$tempPath
echo "LHC18a1_extra anchored to 17r: $HAVELHC18a1_extra $LHC18a1MC_extra $OUTPUTDIR_LHC18a1_extra"


# JJ MC anchored to LHC17x
FindCorrectTrainDirectory $LHC18f5_1 $OUTPUTDIRMC $ALIENDIRMC $LHC17xJJMC
HAVELHC18f5_1=$tempBool
LHC18f5_1=$tempDir
OUTPUTDIR_LHC18f5_1=$tempPath
echo "LHC18f5_1 anchored to 17x: $HAVELHC18f5_1 $LHC18f5_1 $OUTPUTDIR_LHC18f5_1"

FindCorrectTrainDirectory $LHC18f5_2 $OUTPUTDIRMC $ALIENDIRMC $LHC17xJJMC
HAVELHC18f5_2=$tempBool
LHC18f5_2=$tempDir
OUTPUTDIR_LHC18f5_2=$tempPath
echo "LHC18f5_2 anchored to 17x: $HAVELHC18f5_2 $LHC18f5_2 $OUTPUTDIR_LHC18f5_2"



# LHC18x MCs
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

# LHC18x extra MCs
FindCorrectTrainDirectory $LHC18g4MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18g4_extra=$tempBool
LHC18g4MC_extra=$tempDir
OUTPUTDIR_LHC18g4_extra=$tempPath
echo "18g4_extra anchored to 18b: $HAVELHC18g4_extra $LHC18g4MC_extra $OUTPUTDIR_LHC18g4_extra"

FindCorrectTrainDirectory $LHC18g5MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18g5_extra=$tempBool
LHC18g5MC_extra=$tempDir
OUTPUTDIR_LHC18g5_extra=$tempPath
echo "18g5_extra anchored to 18d: $HAVELHC18g5_extra $LHC18g5MC_extra $OUTPUTDIR_LHC18g5_extra"

FindCorrectTrainDirectory $LHC18g6MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18g6_extra=$tempBool
LHC18g6MC_extra=$tempDir
OUTPUTDIR_LHC18g6_extra=$tempPath
echo "18g6_extra anchored to 18e: $HAVELHC18g6_extra $LHC18g6MC_extra $OUTPUTDIR_LHC18g6_extra"

FindCorrectTrainDirectory $LHC18h2MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18h2_extra=$tempBool
LHC18h2MC_extra=$tempDir
OUTPUTDIR_LHC18h2_extra=$tempPath
echo "18h2_extra anchored to 18f: $HAVELHC18h2_extra $LHC18h2MC_extra $OUTPUTDIR_LHC18h2_extra"

FindCorrectTrainDirectory $LHC18h4MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18h4_extra=$tempBool
LHC18h4MC_extra=$tempDir
OUTPUTDIR_LHC18h4_extra=$tempPath
echo "18h4_extra anchored to 18g,h,i,j,k: $HAVELHC18h4_extra $LHC18h4MC_extra $OUTPUTDIR_LHC18h4_extra"

FindCorrectTrainDirectory $LHC18j1MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18j1_extra=$tempBool
LHC18j1MC_extra=$tempDir
OUTPUTDIR_LHC18j1_extra=$tempPath
echo "18j1_extra anchored to 18l: $HAVELHC18j1_extra $LHC18j1MC_extra $OUTPUTDIR_LHC18j1_extra"

FindCorrectTrainDirectory $LHC18j4MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18j4_extra=$tempBool
LHC18j4MC_extra=$tempDir
OUTPUTDIR_LHC18j4_extra=$tempPath
echo "18j4_extra anchored to 18m: $HAVELHC18j4_extra $LHC18j4MC_extra $OUTPUTDIR_LHC18j4_extra"

FindCorrectTrainDirectory $LHC18k1MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18k1_extra=$tempBool
LHC18k1MC_extra=$tempDir
OUTPUTDIR_LHC18k1_extra=$tempPath
echo "18k1_extra anchored to 18n: $HAVELHC18k1_extra $LHC18k1MC_extra $OUTPUTDIR_LHC18k1_extra"

FindCorrectTrainDirectory $LHC18k2MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18k2_extra=$tempBool
LHC18k2MC_extra=$tempDir
OUTPUTDIR_LHC18k2_extra=$tempPath
echo "18k2_extra anchored to 18o: $HAVELHC18k2_extra $LHC18k2MC_extra $OUTPUTDIR_LHC18k2_extra"

FindCorrectTrainDirectory $LHC18k3MC_extra $OUTPUTDIRMC $ALIENDIRMC $LHC18xMCPHY_extra
HAVELHC18k3_extra=$tempBool
LHC18k3MC_extra=$tempDir
OUTPUTDIR_LHC18k3_extra=$tempPath
echo "18k3_extra anchored to 18p: $HAVELHC18k3_extra $LHC18k3MC_extra $OUTPUTDIR_LHC18k3_extra"

# exit


if [ $CLEANUPMAYOR == 0 ]; then
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt

    CopyRunwiseAndMergeAccordingToRunlistData "LHC16d" $HAVELHC16d $OUTPUTDIR_LHC16d $LHC16dData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16e" $HAVELHC16e $OUTPUTDIR_LHC16e $LHC16eData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16f" $HAVELHC16f $OUTPUTDIR_LHC16f $LHC16fData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16g" $HAVELHC16g $OUTPUTDIR_LHC16g $LHC16gData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16h" $HAVELHC16h $OUTPUTDIR_LHC16h $LHC16hData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16i" $HAVELHC16i $OUTPUTDIR_LHC16i $LHC16iData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16j" $HAVELHC16j $OUTPUTDIR_LHC16j $LHC16jData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16k" $HAVELHC16k $OUTPUTDIR_LHC16k $LHC16kData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass2" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16l" $HAVELHC16l $OUTPUTDIR_LHC16l $LHC16lData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass2" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16o" $HAVELHC16o $OUTPUTDIR_LHC16o $LHC16oData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC16p" $HAVELHC16p $OUTPUTDIR_LHC16p $LHC16pData $pathData $baseLegoData "/alice/data/2016" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo

    CopyRunwiseAndMergeAccordingToRunlistData "LHC17c" $HAVELHC17c $OUTPUTDIR_LHC17c $LHC17cData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17e" $HAVELHC17e $OUTPUTDIR_LHC17e $LHC17eData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17f" $HAVELHC17f $OUTPUTDIR_LHC17f $LHC17fData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17h" $HAVELHC17h $OUTPUTDIR_LHC17h $LHC17hData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17i" $HAVELHC17i $OUTPUTDIR_LHC17i $LHC17iData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17j" $HAVELHC17j $OUTPUTDIR_LHC17j $LHC17jData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17k" $HAVELHC17k $OUTPUTDIR_LHC17k $LHC17kData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17l" $HAVELHC17l $OUTPUTDIR_LHC17l $LHC17lData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17m" $HAVELHC17m $OUTPUTDIR_LHC17m $LHC17mData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17o" $HAVELHC17o $OUTPUTDIR_LHC17o $LHC17oData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistData "LHC17r" $HAVELHC17r $OUTPUTDIR_LHC17r $LHC17rData $pathData $baseLegoData "/alice/data/2017" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18b" $HAVELHC18b $OUTPUTDIR_LHC18b $LHC18bData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18d" $HAVELHC18d $OUTPUTDIR_LHC18d $LHC18dData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18e" $HAVELHC18e $OUTPUTDIR_LHC18e $LHC18eData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18f" $HAVELHC18f $OUTPUTDIR_LHC18f $LHC18fData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18g" $HAVELHC18g $OUTPUTDIR_LHC18g $LHC18gData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18h" $HAVELHC18h $OUTPUTDIR_LHC18h $LHC18hData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18i" $HAVELHC18i $OUTPUTDIR_LHC18i $LHC18iData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18j" $HAVELHC18j $OUTPUTDIR_LHC18j $LHC18jData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18k" $HAVELHC18k $OUTPUTDIR_LHC18k $LHC18kData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2" > runlistsToMerge.txt
#     echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4\nDPGTrackIncAccTPCandEMC_5" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_5" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18l" $HAVELHC18l $OUTPUTDIR_LHC18l $LHC18lData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPC_2\nDPGTrackIncAccTPC_3\nDPGTrackIncAccTPC_4\nDPGTrackIncAccTPC_5\nDPGTrackIncAccTPC_6\nDPGTrackIncAccTPC_7\nDPGTrackAndCalo_2\nDPGTrackAndCalo_3\nDPGTrackAndCalo_4\nDPGTrackAndCalo_5\nDPGTrackAndCalo_6\nDPGTrackAndCalo_7" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4\nDPGTrackIncAccTPCandEMC_5\nDPGTrackIncAccTPCandEMC_6\nDPGTrackIncAccTPCandEMC_7" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18m" $HAVELHC18m $OUTPUTDIR_LHC18m $LHC18mData $pathData3 $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1_withTRDtracking" GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18n" $HAVELHC18n $OUTPUTDIR_LHC18n $LHC18nData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18o" $HAVELHC18o $OUTPUTDIR_LHC18o $LHC18oData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo
    echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistData "LHC18p" $HAVELHC18p $OUTPUTDIR_LHC18p $LHC18pData $pathData $baseLegoData "/alice/data/2018" $NSlashes3 runlistsToMerge.txt "pass1" GammaConvCalo

    currentDir=$PWD
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f6" $HAVELHC17f6 $OUTPUTDIR_LHC17f6 $LHC17f6MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f9" $HAVELHC17f9 $OUTPUTDIR_LHC17f9 $LHC17f9MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d1" $HAVELHC17d1 $OUTPUTDIR_LHC17d1 $LHC17d1MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d17" $HAVELHC17d17 $OUTPUTDIR_LHC17d17 $LHC17d17MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f5" $HAVELHC17f5 $OUTPUTDIR_LHC17f5 $LHC17f5MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d3" $HAVELHC17d3 $OUTPUTDIR_LHC17d3 $LHC17d3MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17e5" $HAVELHC17e5 $OUTPUTDIR_LHC17e5 $LHC17e5MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18f1" $HAVELHC18f1 $OUTPUTDIR_LHC18f1 $LHC18f1MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d8" $HAVELHC18d8 $OUTPUTDIR_LHC18d8 $LHC18d8MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d16" $HAVELHC17d16 $OUTPUTDIR_LHC17d16 $LHC17d16MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d18" $HAVELHC17d18 $OUTPUTDIR_LHC17d18 $LHC17d18MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f6_extra" $HAVELHC17f6_extra $OUTPUTDIR_LHC17f6_extra $LHC17f6MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d1_extra" $HAVELHC17d1_extra $OUTPUTDIR_LHC17d1_extra $LHC17d1MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f9_extra" $HAVELHC17f9_extra $OUTPUTDIR_LHC17f9_extra $LHC17f9MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d17_extra" $HAVELHC17d17_extra $OUTPUTDIR_LHC17d17_extra $LHC17d17MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f5_extra" $HAVELHC17f5_extra $OUTPUTDIR_LHC17f5_extra $LHC17f5MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17e5_extra" $HAVELHC17e5_extra $OUTPUTDIR_LHC17e5_extra $LHC17e5MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d3_extra" $HAVELHC17d3_extra $OUTPUTDIR_LHC17d3_extra $LHC17d3MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18f1_extra" $HAVELHC18f1_extra $OUTPUTDIR_LHC18f1_extra $LHC18f1MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d8_extra" $HAVELHC18d8_extra $OUTPUTDIR_LHC18d8_extra $LHC18d8MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d16_extra" $HAVELHC17d16_extra $OUTPUTDIR_LHC17d16_extra $LHC17d16MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d18_extra" $HAVELHC17d18_extra $OUTPUTDIR_LHC17d18_extra $LHC17d18MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f6_extra2" $HAVELHC17f6_extra2 $OUTPUTDIR_LHC17f6_extra2 $LHC17f6MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d1_extra2" $HAVELHC17d1_extra2 $OUTPUTDIR_LHC17d1_extra2 $LHC17d1MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f9_extra2" $HAVELHC17f9_extra2 $OUTPUTDIR_LHC17f9_extra2 $LHC17f9MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d17_extra2" $HAVELHC17d17_extra2 $OUTPUTDIR_LHC17d17_extra2 $LHC17d17MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17f5_extra2" $HAVELHC17f5_extra2 $OUTPUTDIR_LHC17f5_extra2 $LHC17f5MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17e5_extra2" $HAVELHC17e5_extra2 $OUTPUTDIR_LHC17e5_extra2 $LHC17e5MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d3_extra2" $HAVELHC17d3_extra2 $OUTPUTDIR_LHC17d3_extra2 $LHC17d3MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d16_extra2" $HAVELHC17d16_extra2 $OUTPUTDIR_LHC17d16_extra2 $LHC17d16MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17d18_extra2" $HAVELHC17d18_extra2 $OUTPUTDIR_LHC17d18_extra2 $LHC17d18MC_extra2 $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC18f5_1" $HAVELHC18f5_1 $OUTPUTDIR_LHC18f5_1 $LHC18f5_1 $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistJJMC "LHC18f5_2" $HAVELHC18f5_2 $OUTPUTDIR_LHC18f5_2 $LHC18f5_2 $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d3" $HAVELHC18d3 $OUTPUTDIR_LHC18d3 $LHC18d3MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18c12" $HAVELHC18c12 $OUTPUTDIR_LHC18c12 $LHC18c12MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17k4" $HAVELHC17k4 $OUTPUTDIR_LHC17k4 $LHC17k4MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17h11" $HAVELHC17h11 $OUTPUTDIR_LHC17h11 $LHC17h11MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18c13" $HAVELHC18c13 $OUTPUTDIR_LHC18c13 $LHC18c13MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18a8" $HAVELHC18a8 $OUTPUTDIR_LHC18a8 $LHC18a8MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17l5" $HAVELHC17l5 $OUTPUTDIR_LHC17l5 $LHC17l5MC $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18a9" $HAVELHC18a9 $OUTPUTDIR_LHC18a9 $LHC18a9MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18a1" $HAVELHC18a1 $OUTPUTDIR_LHC18a1 $LHC18a1MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18d3_extra" $HAVELHC18d3_extra $OUTPUTDIR_LHC18d3_extra $LHC18d3MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18c12_extra" $HAVELHC18c12_extra $OUTPUTDIR_LHC18c12_extra $LHC18c12MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17k4_extra" $HAVELHC17k4_extra $OUTPUTDIR_LHC17k4_extra $LHC17k4MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17h11_extra" $HAVELHC17h11_extra $OUTPUTDIR_LHC17h11_extra $LHC17h11MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18c13_extra" $HAVELHC18c13_extra $OUTPUTDIR_LHC18c13_extra $LHC18c13MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18a8_extra" $HAVELHC18a8_extra $OUTPUTDIR_LHC18a8_extra $LHC18a8MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC17l5_extra" $HAVELHC17l5_extra $OUTPUTDIR_LHC17l5_extra $LHC17l5MC_extra $pathMC $baseLegoMC "/alice/sim/2017" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18a9_extra" $HAVELHC18a9_extra $OUTPUTDIR_LHC18a9_extra $LHC18a9MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18a1_extra" $HAVELHC18a1_extra $OUTPUTDIR_LHC18a1_extra $LHC18a1MC_extra $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g4" $HAVELHC18g4 $OUTPUTDIR_LHC18g4 $LHC18g4MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g4_extra" $HAVELHC18g4_extra $OUTPUTDIR_LHC18g4_extra $LHC18g4MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g5" $HAVELHC18g5 $OUTPUTDIR_LHC18g5 $LHC18g5MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g5_extra" $HAVELHC18g5_extra $OUTPUTDIR_LHC18g5_extra $LHC18g5MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g6" $HAVELHC18g6 $OUTPUTDIR_LHC18g6 $LHC18g6MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18g6_extra" $HAVELHC18g6_extra $OUTPUTDIR_LHC18g6_extra $LHC18g6MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h2" $HAVELHC18h2 $OUTPUTDIR_LHC18h2 $LHC18h2MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h2_extra" $HAVELHC18h2_extra $OUTPUTDIR_LHC18h2_extra $LHC18h2MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h4" $HAVELHC18h4 $OUTPUTDIR_LHC18h4 $LHC18h4MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18h4_extra" $HAVELHC18h4_extra $OUTPUTDIR_LHC18h4_extra $LHC18h4MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
#     echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4\nDPGTrackIncAccTPCandEMC_5" > runlistsToMerge.txt
    echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_5" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j1" $HAVELHC18j1 $OUTPUTDIR_LHC18j1 $LHC18j1MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j1_extra" $HAVELHC18j1_extra $OUTPUTDIR_LHC18j1_extra $LHC18j1MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPC_2\nDPGTrackIncAccTPC_3\nDPGTrackIncAccTPC_4\nDPGTrackIncAccTPC_5\nDPGTrackIncAccTPC_6\nDPGTrackIncAccTPC_7\nDPGTrackAndCalo_2\nDPGTrackAndCalo_3\nDPGTrackAndCalo_4\nDPGTrackAndCalo_5\nDPGTrackAndCalo_6\nDPGTrackAndCalo_7" > runlistsToMerge.txt
#     echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4\nDPGTrackIncAccTPCandEMC_5\nDPGTrackIncAccTPCandEMC_6\nDPGTrackIncAccTPCandEMC_7" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j4" $HAVELHC18j4 $OUTPUTDIR_LHC18j4 $LHC18j4MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18j4_extra" $HAVELHC18j4_extra $OUTPUTDIR_LHC18j4_extra $LHC18j4MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k1" $HAVELHC18k1 $OUTPUTDIR_LHC18k1 $LHC18k1MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k1_extra" $HAVELHC18k1_extra $OUTPUTDIR_LHC18k1_extra $LHC18k1MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k2" $HAVELHC18k2 $OUTPUTDIR_LHC18k2 $LHC18k2MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k2_extra" $HAVELHC18k2_extra $OUTPUTDIR_LHC18k2_extra $LHC18k2MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    echo -e "DPGTrackIncAccTPCandEMC\nDPGTrackIncAccTPCandEMC_1\nDPGTrackIncAccTPCandEMC_2\nDPGTrackIncAccTPCandEMC_3\nDPGTrackIncAccTPCandEMC_4" > runlistsToMerge.txt
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k3" $HAVELHC18k3 $OUTPUTDIR_LHC18k3 $LHC18k3MC $pathMC $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo
    CopyRunwiseAndMergeAccordingToRunlistMC "LHC18k3_extra" $HAVELHC18k3_extra $OUTPUTDIR_LHC18k3_extra $LHC18k3MC_extra $pathMC2 $baseLegoMC "/alice/sim/2018" $NSlashes3 runlistsToMerge.txt GammaConvCalo

    echo "Change Structure If Needed"

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
    listsToMerge=`cat runlistsToMerge.txt`
    for runListName in $listsToMerge; do
        if [ $HAVELHC16d == 1 ]; then
            ls $OUTPUTDIR_LHC16d/GammaConvCalo-$runListName\_*.root > fileLHC16d.txt
            fileNumbers=`cat fileLHC16d.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16d $NSlashes "LHC16d-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16e == 1 ]; then
            ls $OUTPUTDIR_LHC16e/GammaConvCalo-$runListName\_*.root > fileLHC16e.txt
            fileNumbers=`cat fileLHC16e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16e $NSlashes "LHC16e-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16f == 1 ]; then
            ls $OUTPUTDIR_LHC16f/GammaConvCalo-$runListName\_*.root > fileLHC16f.txt
            fileNumbers=`cat fileLHC16f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16f $NSlashes "LHC16f-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16g == 1 ]; then
            ls $OUTPUTDIR_LHC16g/GammaConvCalo-$runListName\_*.root > fileLHC16g.txt
            fileNumbers=`cat fileLHC16g.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16g $NSlashes "LHC16g-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16h == 1 ]; then
            ls $OUTPUTDIR_LHC16h/GammaConvCalo-$runListName\_*.root > fileLHC16h.txt
            fileNumbers=`cat fileLHC16h.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16h $NSlashes "LHC16h-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16i == 1 ]; then
            ls $OUTPUTDIR_LHC16i/GammaConvCalo-$runListName\_*.root > fileLHC16i.txt
            fileNumbers=`cat fileLHC16i.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16i $NSlashes "LHC16i-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16j == 1 ]; then
            ls $OUTPUTDIR_LHC16j/GammaConvCalo-$runListName\_*.root > fileLHC16j.txt
            fileNumbers=`cat fileLHC16j.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16j-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16k == 1 ]; then
            ls $OUTPUTDIR_LHC16k/GammaConvCalo-$runListName\_*.root > fileLHC16k.txt
            fileNumbers=`cat fileLHC16k.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16k-pass2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16l == 1 ]; then
            ls $OUTPUTDIR_LHC16l/GammaConvCalo-$runListName\_*.root > fileLHC16l.txt
            fileNumbers=`cat fileLHC16l.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16j $NSlashes "LHC16l-pass2-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16o == 1 ]; then
            ls $OUTPUTDIR_LHC16o/GammaConvCalo-$runListName\_*.root > fileLHC16o.txt
            fileNumbers=`cat fileLHC16o.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16o $NSlashes "LHC16o-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC16p == 1 ]; then
            ls $OUTPUTDIR_LHC16p/GammaConvCalo-$runListName\_*.root > fileLHC16p.txt
            fileNumbers=`cat fileLHC16p.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16p $NSlashes "LHC16p-pass$passNr-$runListName" "-$runListName"
            done;
        fi


        if [ $HAVELHC17c == 1 ]; then
            ls $OUTPUTDIR_LHC17c/GammaConvCalo-$runListName\_*.root > fileLHC17c.txt
            fileNumbers=`cat fileLHC17c.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17c $NSlashes "LHC17c-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17e == 1 ]; then
            ls $OUTPUTDIR_LHC17e/GammaConvCalo-$runListName\_*.root > fileLHC17e.txt
            fileNumbers=`cat fileLHC17e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e $NSlashes "LHC17e-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17f == 1 ]; then
            ls $OUTPUTDIR_LHC17f/GammaConvCalo-$runListName\_*.root > fileLHC17f.txt
            fileNumbers=`cat fileLHC17f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f $NSlashes "LHC17f-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17h == 1 ]; then
            ls $OUTPUTDIR_LHC17h/GammaConvCalo-$runListName\_*.root > fileLHC17h.txt
            fileNumbers=`cat fileLHC17h.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17h $NSlashes "LHC17h-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17i == 1 ]; then
            ls $OUTPUTDIR_LHC17i/GammaConvCalo-$runListName\_*.root > fileLHC17i.txt
            fileNumbers=`cat fileLHC17i.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17i $NSlashes "LHC17i-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17j == 1 ]; then
            ls $OUTPUTDIR_LHC17j/GammaConvCalo-$runListName\_*.root > fileLHC17j.txt
            fileNumbers=`cat fileLHC17j.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17j $NSlashes "LHC17j-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17k == 1 ]; then
            ls $OUTPUTDIR_LHC17k/GammaConvCalo-$runListName\_*.root > fileLHC17k.txt
            fileNumbers=`cat fileLHC17k.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17k $NSlashes "LHC17k-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17l == 1 ]; then
            ls $OUTPUTDIR_LHC17l/GammaConvCalo-$runListName\_*.root > fileLHC17l.txt
            fileNumbers=`cat fileLHC17l.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17l $NSlashes "LHC17l-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17m == 1 ]; then
            ls $OUTPUTDIR_LHC17m/GammaConvCalo-$runListName\_*.root > fileLHC17m.txt
            fileNumbers=`cat fileLHC17m.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17m $NSlashes "LHC17m-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17o == 1 ]; then
            ls $OUTPUTDIR_LHC17o/GammaConvCalo-$runListName\_*.root > fileLHC17o.txt
            fileNumbers=`cat fileLHC17o.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17o $NSlashes "LHC17o-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC17r == 1 ]; then
            ls $OUTPUTDIR_LHC17r/GammaConvCalo-$runListName\_*.root > fileLHC17r.txt
            fileNumbers=`cat fileLHC17r.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17r $NSlashes "LHC17r-pass$passNr-$runListName" "-$runListName"
            done;
        fi



        if [ $HAVELHC18b == 1 ]; then
            ls $OUTPUTDIR_LHC18b/GammaConvCalo-$runListName\_*.root > fileLHC18b.txt
            fileNumbers=`cat fileLHC18b.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18b $NSlashes "LHC18b-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18d == 1 ]; then
            ls $OUTPUTDIR_LHC18d/GammaConvCalo-$runListName\_*.root > fileLHC18d.txt
            fileNumbers=`cat fileLHC18d.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18d $NSlashes "LHC18d-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18e == 1 ]; then
            ls $OUTPUTDIR_LHC18e/GammaConvCalo-$runListName\_*.root > fileLHC18e.txt
            fileNumbers=`cat fileLHC18e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18e $NSlashes "LHC18e-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18f == 1 ]; then
            ls $OUTPUTDIR_LHC18f/GammaConvCalo-$runListName\_*.root > fileLHC18f.txt
            fileNumbers=`cat fileLHC18f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f $NSlashes "LHC18f-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18g == 1 ]; then
            ls $OUTPUTDIR_LHC18g/GammaConvCalo-$runListName\_*.root > fileLHC18g.txt
            fileNumbers=`cat fileLHC18g.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18g $NSlashes "LHC18g-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h == 1 ]; then
            ls $OUTPUTDIR_LHC18h/GammaConvCalo-$runListName\_*.root > fileLHC18h.txt
            fileNumbers=`cat fileLHC18h.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18h $NSlashes "LHC18h-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18i == 1 ]; then
            ls $OUTPUTDIR_LHC18i/GammaConvCalo-$runListName\_*.root > fileLHC18i.txt
            fileNumbers=`cat fileLHC18i.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18i $NSlashes "LHC18i-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j == 1 ]; then
            ls $OUTPUTDIR_LHC18j/GammaConvCalo-$runListName\_*.root > fileLHC18j.txt
            fileNumbers=`cat fileLHC18j.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j $NSlashes "LHC18j-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k == 1 ]; then
            ls $OUTPUTDIR_LHC18k/GammaConvCalo-$runListName\_*.root > fileLHC18k.txt
            fileNumbers=`cat fileLHC18k.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18k $NSlashes "LHC18k-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18l == 1 ]; then
            ls $OUTPUTDIR_LHC18l/GammaConvCalo-$runListName\_*.root > fileLHC18l.txt
            fileNumbers=`cat fileLHC18l.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18l $NSlashes "LHC18l-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18m == 1 ]; then
            ls $OUTPUTDIR_LHC18m/GammaConvCalo-$runListName\_*.root > fileLHC18m.txt
            fileNumbers=`cat fileLHC18m.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18m $NSlashes "LHC18m-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18n == 1 ]; then
            ls $OUTPUTDIR_LHC18n/GammaConvCalo-$runListName\_*.root > fileLHC18n.txt
            fileNumbers=`cat fileLHC18n.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18n $NSlashes "LHC18n-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18o == 1 ]; then
            ls $OUTPUTDIR_LHC18o/GammaConvCalo-$runListName\_*.root > fileLHC18o.txt
            fileNumbers=`cat fileLHC18o.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18o $NSlashes "LHC18o-pass$passNr-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18p == 1 ]; then
            ls $OUTPUTDIR_LHC18p/GammaConvCalo-$runListName\_*.root > fileLHC18p.txt
            fileNumbers=`cat fileLHC18p.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18p $NSlashes "LHC18p-pass$passNr-$runListName" "-$runListName"
            done;
        fi
    done

    echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
#         echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo" > runlistsToMerge.txt
    listsToMerge=`cat runlistsToMerge.txt`
    for runListName in $listsToMerge; do
        # MC for LHC16d
        if [ $HAVELHC17f6 == 1 ]; then
            ls $OUTPUTDIR_LHC17f6/GammaConvCalo-$runListName\_*.root > fileLHC17f6.txt
            fileNumbers=`cat fileLHC17f6.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6 $NSlashes "MC_LHC17f6-anchor16d-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16e
        if [ $HAVELHC17f9 == 1 ]; then
            ls $OUTPUTDIR_LHC17f9/GammaConvCalo-$runListName\_*.root > fileLHC17f9.txt
            fileNumbers=`cat fileLHC17f9.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9 $NSlashes "MC_LHC17f9-anchor16e-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16f
        if [ $HAVELHC17d1 == 1 ]; then
            ls $OUTPUTDIR_LHC17d1/GammaConvCalo-$runListName\_*.root > fileLHC17d1.txt
            fileNumbers=`cat fileLHC17d1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d1 $NSlashes "MC_LHC17d1-anchor16f-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16g
        if [ $HAVELHC17d17 == 1 ]; then
            ls $OUTPUTDIR_LHC17d17/GammaConvCalo-$runListName\_*.root > fileLHC17d17.txt
            fileNumbers=`cat fileLHC17d17.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17 $NSlashes "MC_LHC17d17-anchor16g-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16h
        if [ $HAVELHC17f5 == 1 ]; then
            ls $OUTPUTDIR_LHC17f5/GammaConvCalo-$runListName\_*.root > fileLHC17f5.txt
            fileNumbers=`cat fileLHC17f5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5-anchor16h-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16i
        if [ $HAVELHC17d3 == 1 ]; then
            ls $OUTPUTDIR_LHC17d3/GammaConvCalo-$runListName\_*.root > fileLHC17d3.txt
            fileNumbers=`cat fileLHC17d3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3 $NSlashes "MC_LHC17d3-anchor16i-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16j
        if [ $HAVELHC17e5 == 1 ]; then
            ls $OUTPUTDIR_LHC17e5/GammaConvCalo-$runListName\_*.root > fileLHC17e5.txt
            fileNumbers=`cat fileLHC17e5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5 $NSlashes "MC_LHC17e5-anchor16j-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16k
        if [ $HAVELHC18f1 == 1 ]; then
            ls $OUTPUTDIR_LHC18f1/GammaConvCalo-$runListName\_*.root > fileLHC18f1.txt
            fileNumbers=`cat fileLHC18f1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f1 $NSlashes "MC_LHC18f1-anchor16k-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16l
        if [ $HAVELHC18d8 == 1 ]; then
            ls $OUTPUTDIR_LHC18d8/GammaConvCalo-$runListName\_*.root > fileLHC18d8.txt
            fileNumbers=`cat fileLHC18d8.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18d8 $NSlashes "MC_LHC18d8-anchor16l-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16o
        if [ $HAVELHC17d16 == 1 ]; then
            ls $OUTPUTDIR_LHC17d16/GammaConvCalo-$runListName\_*.root > fileLHC17d16.txt
            fileNumbers=`cat fileLHC17d16.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16 $NSlashes "MC_LHC17d16-anchor16o-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16p
        if [ $HAVELHC17d18 == 1 ]; then
            ls $OUTPUTDIR_LHC17d18/GammaConvCalo-$runListName\_*.root > fileLHC17d18.txt
            fileNumbers=`cat fileLHC17d18.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18 $NSlashes "MC_LHC17d18-anchor16p-$runListName" "-$runListName"
            done;
        fi

        #extra MC for LHC16d
        if [ $HAVELHC17f6_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17f6_extra/GammaConvCalo-$runListName\_*.root > fileLHC17f6_extra.txt
            fileNumbers=`cat fileLHC17f6_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6_extra $NSlashes "MC_LHC17f6_extra-anchor16d-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16e
        if [ $HAVELHC17f9_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17f9_extra/GammaConvCalo-$runListName\_*.root > fileLHC17f9_extra.txt
            fileNumbers=`cat fileLHC17f9_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9_extra $NSlashes "MC_LHC17f9_extra-anchor16e-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16f
        if [ $HAVELHC17d1_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17d1_extra/GammaConvCalo-$runListName\_*.root > fileLHC17d1_extra.txt
            fileNumbers=`cat fileLHC17d1_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d1_extra $NSlashes "MC_LHC17d1_extra-anchor16f-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16g
        if [ $HAVELHC17d17_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17d17_extra/GammaConvCalo-$runListName\_*.root > fileLHC17d17_extra.txt
            fileNumbers=`cat fileLHC17d17_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17_extra $NSlashes "MC_LHC17d17-anchor16g_extra-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16h
        if [ $HAVELHC17f5_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17f5_extra/GammaConvCalo-$runListName\_*.root > fileLHC17f5_extra.txt
            fileNumbers=`cat fileLHC17f5_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5_extra-anchor16h-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16i
        if [ $HAVELHC17d3_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17d3_extra/GammaConvCalo-$runListName\_*.root > fileLHC17d3_extra.txt
            fileNumbers=`cat fileLHC17d3_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3_extra $NSlashes "MC_LHC17d3_extra-anchor16i-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16j
        if [ $HAVELHC17e5_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17e5_extra/GammaConvCalo-$runListName\_*.root > fileLHC17e5_extra.txt
            fileNumbers=`cat fileLHC17e5_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5_extra $NSlashes "MC_LHC17e5_extra-anchor16j-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16k
        if [ $HAVELHC18f1_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18f1_extra/GammaConvCalo-$runListName\_*.root > fileLHC18f1_extra.txt
            fileNumbers=`cat fileLHC18f1_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f1_extra $NSlashes "MC_LHC18f1-anchor16k-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16l
        if [ $HAVELHC18d8_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18d8_extra/GammaConvCalo-$runListName\_*.root > fileLHC18d8_extra.txt
            fileNumbers=`cat fileLHC18d8_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18d8_extra $NSlashes "MC_LHC18d8_extra-anchor16l-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16o
        if [ $HAVELHC17d16_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17d16_extra/GammaConvCalo-$runListName\_*.root > fileLHC17d16_extra.txt
            fileNumbers=`cat fileLHC17d16_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16_extra $NSlashes "MC_LHC17d16_extra-anchor16o-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16p
        if [ $HAVELHC17d18_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17d18_extra/GammaConvCalo-$runListName\_*.root > fileLHC17d18_extra.txt
            fileNumbers=`cat fileLHC17d18_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18_extra $NSlashes "MC_LHC17d18_extra-anchor16p-$runListName" "-$runListName"
            done;
        fi

        #extra2 MC for LHC16d
        if [ $HAVELHC17f6_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17f6_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17f6_extra2.txt
            fileNumbers=`cat fileLHC17f6_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f6_extra2 $NSlashes "MC_LHC17f6_extra2-anchor16d-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16e
        if [ $HAVELHC17f9_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17f9_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17f9_extra2.txt
            fileNumbers=`cat fileLHC17f9_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f9_extra2 $NSlashes "MC_LHC17f9_extra2-anchor16e-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16f
        if [ $HAVELHC17d1_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17d1_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17d1_extra2.txt
            fileNumbers=`cat fileLHC17d1_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d1_extra2 $NSlashes "MC_LHC17d1_extra2-anchor16f-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16g
        if [ $HAVELHC17d17_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17d17_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17d17_extra2.txt
            fileNumbers=`cat fileLHC17d17_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d17_extra2 $NSlashes "MC_LHC17d17-anchor16g_extra2-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16h
        if [ $HAVELHC17f5_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17f5_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17f5_extra2.txt
            fileNumbers=`cat fileLHC17f5_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17f5 $NSlashes "MC_LHC17f5_extra2-anchor16h-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16i
        if [ $HAVELHC17d3_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17d3_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17d3_extra2.txt
            fileNumbers=`cat fileLHC17d3_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d3_extra2 $NSlashes "MC_LHC17d3_extra2-anchor16i-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16j
        if [ $HAVELHC17e5_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17e5_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17e5_extra2.txt
            fileNumbers=`cat fileLHC17e5_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17e5_extra2 $NSlashes "MC_LHC17e5_extra2-anchor16j-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16o
        if [ $HAVELHC17d16_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17d16_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17d16_extra2.txt
            fileNumbers=`cat fileLHC17d16_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d16_extra2 $NSlashes "MC_LHC17d16_extra2-anchor16o-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC16p
        if [ $HAVELHC17d18_extra2 == 1 ]; then
            ls $OUTPUTDIR_LHC17d18_extra2/GammaConvCalo-$runListName\_*.root > fileLHC17d18_extra2.txt
            fileNumbers=`cat fileLHC17d18_extra2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17d18_extra2 $NSlashes "MC_LHC17d18_extra2-anchor16p-$runListName" "-$runListName"
            done;
        fi

        # MC for LHC17c,e,f
        if [ $HAVELHC18d3 == 1 ]; then
            ls $OUTPUTDIR_LHC18d3/GammaConvCalo-$runListName\_*.root > fileLHC18d3.txt
            fileNumbers=`cat fileLHC18d3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18d3 $NSlashes "MC_LHC18d3-anchor17c-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC17h
        if [ $HAVELHC18c12 == 1 ]; then
            ls $OUTPUTDIR_LHC18c12/GammaConvCalo-$runListName\_*.root > fileLHC18c12.txt
            fileNumbers=`cat fileLHC18c12.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18c12 $NSlashes "MC_LHC18c12-anchor17h-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17i
        if [ $HAVELHC17k4 == 1 ]; then
            ls $OUTPUTDIR_LHC17k4/GammaConvCalo-$runListName\_*.root > fileLHC17k4.txt
            fileNumbers=`cat fileLHC17k4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17k4 $NSlashes "MC_LHC17k4-anchor17i-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17j
        if [ $HAVELHC17h11 == 1 ]; then
            ls $OUTPUTDIR_LHC17h11/GammaConvCalo-$runListName\_*.root > fileLHC17h11.txt
            fileNumbers=`cat fileLHC17h11.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17h11 $NSlashes "MC_LHC17h11-anchor17j-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17k
        if [ $HAVELHC18c13 == 1 ]; then
            ls $OUTPUTDIR_LHC18c13/GammaConvCalo-$runListName\_*.root > fileLHC18c13.txt
            fileNumbers=`cat fileLHC18c13.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18c13 $NSlashes "MC_LHC18c13-anchor17k-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17l
        if [ $HAVELHC18a8 == 1 ]; then
            ls $OUTPUTDIR_LHC18a8/GammaConvCalo-$runListName\_*.root > fileLHC18a8.txt
            fileNumbers=`cat fileLHC18a8.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18a8 $NSlashes "MC_LHC18a8-anchor17l-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17m
        if [ $HAVELHC17l5 == 1 ]; then
            ls $OUTPUTDIR_LHC17l5/GammaConvCalo-$runListName\_*.root > fileLHC17l5.txt
            fileNumbers=`cat fileLHC17l5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17l5 $NSlashes "MC_LHC17l5-anchor17m-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17o
        if [ $HAVELHC18a9 == 1 ]; then
            ls $OUTPUTDIR_LHC18a9/GammaConvCalo-$runListName\_*.root > fileLHC18a9.txt
            fileNumbers=`cat fileLHC18a9.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18a9 $NSlashes "MC_LHC18a9-anchor17r-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17r
        if [ $HAVELHC18a1 == 1 ]; then
            ls $OUTPUTDIR_LHC18a1/GammaConvCalo-$runListName\_*.root > fileLHC18a1.txt
            fileNumbers=`cat fileLHC18a1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18a1 $NSlashes "MC_LHC18a1-anchor17r-$runListName" "-$runListName"
            done;
        fi



        # MC for LHC17c,e,f
        if [ $HAVELHC18d3_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18d3_extra/GammaConvCalo-$runListName\_*.root > fileLHC18d3_extra.txt
            fileNumbers=`cat fileLHC18d3_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18d3 $NSlashes "MC_LHC18d3_extra-anchor17c-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC17h
        if [ $HAVELHC18c12_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18c12_extra/GammaConvCalo-$runListName\_*.root > fileLHC18c12_extra.txt
            fileNumbers=`cat fileLHC18c12_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18c12_extra $NSlashes "MC_LHC18c12_extra-anchor17h-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17i
        if [ $HAVELHC17k4_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17k4_extra/GammaConvCalo-$runListName\_*.root > fileLHC17k4_extra.txt
            fileNumbers=`cat fileLHC17k4_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17k4_extra $NSlashes "MC_LHC17k4_extra-anchor17i-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17j
        if [ $HAVELHC17h11_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17h11_extra/GammaConvCalo-$runListName\_*.root > fileLHC17h11_extra.txt
            fileNumbers=`cat fileLHC17h11_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17h11_extra $NSlashes "MC_LHC17h11_extra-anchor17j-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17k
        if [ $HAVELHC18c13_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18c13_extra/GammaConvCalo-$runListName\_*.root > fileLHC18c13_extra.txt
            fileNumbers=`cat fileLHC18c13_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18c13_extra $NSlashes "MC_LHC18c13_extra-anchor17k-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17l
        if [ $HAVELHC18a8_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18a8_extra/GammaConvCalo-$runListName\_*.root > fileLHC18a8_extra.txt
            fileNumbers=`cat fileLHC18a8_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18a8_extra $NSlashes "MC_LHC18a8_extra-anchor17l-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17m
        if [ $HAVELHC17l5_extra == 1 ]; then
            ls $OUTPUTDIR_LHC17l5_extra/GammaConvCalo-$runListName\_*.root > fileLHC17l5_extra.txt
            fileNumbers=`cat fileLHC17l5_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC17l5_extra $NSlashes "MC_LHC17l5_extra-anchor17m-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17o
        if [ $HAVELHC18a9_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18a9_extra/GammaConvCalo-$runListName\_*.root > fileLHC18a9_extra.txt
            fileNumbers=`cat fileLHC18a9_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18a9_extra $NSlashes "MC_LHC18a9_extra-anchor17r-$runListName" "-$runListName"
            done;
        fi
        #MC for LHC17r
        if [ $HAVELHC18a1_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18a1_extra/GammaConvCalo-$runListName\_*.root > fileLHC18a1_extra.txt
            fileNumbers=`cat fileLHC18a1_extra.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18a1_extra $NSlashes "MC_LHC18a1_extra-anchor17r-$runListName" "-$runListName"
            done;
        fi

        # JJ MC anch to LHC17x
        if [ $HAVELHC18f5_1 == 1 ]; then
            ls $OUTPUTDIR_LHC18f5_1/GammaConvCalo-$runListName\_*.root > fileLHC18f5_1.txt
            fileNumbers=`cat fileLHC18f5_1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f5_1 $NSlashes "MC_LHC18f5_1-anchor17r-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18f5_2 == 1 ]; then
            ls $OUTPUTDIR_LHC18f5_2/GammaConvCalo-$runListName\_*.root > fileLHC18f5_2.txt
            fileNumbers=`cat fileLHC18f5_2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18f5_2 $NSlashes "MC_LHC18f5_2-anchor17r-$runListName" "-$runListName"
            done;
        fi


        # MC for LHC18b
        if [ $HAVELHC18g4 == 1 ]; then
            ls $OUTPUTDIR_LHC18g4/GammaConvCalo-$runListName\_*.root > fileLHC18g4.txt
            fileNumbers=`cat fileLHC18g4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18g4 $NSlashes "MC_LHC18g4-anchor18b-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18g4_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18g4_extra/GammaConvCalo-$runListName\_*.root > fileLHC18g4.txt
            fileNumbers=`cat fileLHC18g4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18g4_extra $NSlashes "MC_LHC18g4_extra-anchor18b-$runListName" "-$runListName"
            done;
        fi
        # MC for LHC18d
        if [ $HAVELHC18g5 == 1 ]; then
            ls $OUTPUTDIR_LHC18g5/GammaConvCalo-$runListName\_*.root > fileLHC18g5.txt
            fileNumbers=`cat fileLHC18g5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18g5 $NSlashes "MC_LHC18g5-anchor18d-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18g5_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18g5_extra/GammaConvCalo-$runListName\_*.root > fileLHC18g5.txt
            fileNumbers=`cat fileLHC18g5.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18g5_extra $NSlashes "MC_LHC18g5_extra-anchor18d-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18g6 == 1 ]; then
            ls $OUTPUTDIR_LHC18g6/GammaConvCalo-$runListName\_*.root > fileLHC18g6.txt
            fileNumbers=`cat fileLHC18g6.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18g6 $NSlashes "MC_LHC18g6-anchor18e-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18g6_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18g6_extra/GammaConvCalo-$runListName\_*.root > fileLHC18g6.txt
            fileNumbers=`cat fileLHC18g6.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18g6_extra $NSlashes "MC_LHC18g6_extra-anchor18e-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h2 == 1 ]; then
            ls $OUTPUTDIR_LHC18h2/GammaConvCalo-$runListName\_*.root > fileLHC18h2.txt
            fileNumbers=`cat fileLHC18h2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18h2 $NSlashes "MC_LHC18h2-anchor18f-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h2_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18h2_extra/GammaConvCalo-$runListName\_*.root > fileLHC18h2.txt
            fileNumbers=`cat fileLHC18h2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18h2_extra $NSlashes "MC_LHC18h2_extra-anchor18f-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h4 == 1 ]; then
            ls $OUTPUTDIR_LHC18h4/GammaConvCalo-$runListName\_*.root > fileLHC18h4.txt
            fileNumbers=`cat fileLHC18h4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18h4 $NSlashes "MC_LHC18h4-anchor18ghijk-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18h4_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18h4_extra/GammaConvCalo-$runListName\_*.root > fileLHC18h4.txt
            fileNumbers=`cat fileLHC18h4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18h4_extra $NSlashes "MC_LHC18h4_extra-anchor18ghijk-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j1 == 1 ]; then
            ls $OUTPUTDIR_LHC18j1/GammaConvCalo-$runListName\_*.root > fileLHC18j1.txt
            fileNumbers=`cat fileLHC18j1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j1 $NSlashes "MC_LHC18j1-anchor18l-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j1_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18j1_extra/GammaConvCalo-$runListName\_*.root > fileLHC18j1.txt
            fileNumbers=`cat fileLHC18j1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j1_extra $NSlashes "MC_LHC18j1_extra-anchor18l-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j4 == 1 ]; then
            ls $OUTPUTDIR_LHC18j4/GammaConvCalo-$runListName\_*.root > fileLHC18j4.txt
            fileNumbers=`cat fileLHC18j4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j4 $NSlashes "MC_LHC18j4-anchor18m-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18j4_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18j4_extra/GammaConvCalo-$runListName\_*.root > fileLHC18j4.txt
            fileNumbers=`cat fileLHC18j4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18j4_extra $NSlashes "MC_LHC18j4_extra-anchor18m-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k1 == 1 ]; then
            ls $OUTPUTDIR_LHC18k1/GammaConvCalo-$runListName\_*.root > fileLHC18k1.txt
            fileNumbers=`cat fileLHC18k1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18k1 $NSlashes "MC_LHC18k1-anchor18n-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k1_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18k1_extra/GammaConvCalo-$runListName\_*.root > fileLHC18k1.txt
            fileNumbers=`cat fileLHC18k1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18k1_extra $NSlashes "MC_LHC18k1_extra-anchor18n-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k2 == 1 ]; then
            ls $OUTPUTDIR_LHC18k2/GammaConvCalo-$runListName\_*.root > fileLHC18k2.txt
            fileNumbers=`cat fileLHC18k2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18k2_extra $NSlashes "MC_LHC18k2_extra-anchor18o-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k2_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18k2_extra/GammaConvCalo-$runListName\_*.root > fileLHC18k2.txt
            fileNumbers=`cat fileLHC18k2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18k2 $NSlashes "MC_LHC18k2-anchor18o-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k3 == 1 ]; then
            ls $OUTPUTDIR_LHC18k3/GammaConvCalo-$runListName\_*.root > fileLHC18k3.txt
            fileNumbers=`cat fileLHC18k3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18k3 $NSlashes "MC_LHC18k3-anchor18p-$runListName" "-$runListName"
            done;
        fi
        if [ $HAVELHC18k3_extra == 1 ]; then
            ls $OUTPUTDIR_LHC18k3_extra/GammaConvCalo-$runListName\_*.root > fileLHC18k3.txt
            fileNumbers=`cat fileLHC18k3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC18k3_extra $NSlashes "MC_LHC18k3_extra-anchor18p-$runListName" "-$runListName"
            done;
        fi
    done

    echo "Download Done"

    if [ $MERGEON == 1 ]; then
        echo "Starting Merging"

        if [ $MERGEONData == 1 ]; then
            echo -e "DPGTrack\nDPGTrackIncAccTPC\nDPGTrackAndCalo\nDPGTrackIncAccTPCandEMC" > runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            for runListName in $listsToMerge; do
                ls $OUTPUTDIR/GammaConvCalo_LHC16o-pass$passNr-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                #16er daten
                periodList=`echo -e "d-pass1\ne-pass1\nf-pass1\ng-pass1\nh-pass1\ni-pass1\nj-pass1\nk-pass2\nl-pass2\no-pass1\np-pass1"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC16$periodID-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            outAdd=`echo $periodID  | cut -d "-" -f 1 `
                            nameOut+=$outAdd
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        else
                            echo $currFile " does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC16$nameOut-pass$passNr-2-$runListName\_$number.root
                done
                #17er daten
                ls $OUTPUTDIR/GammaConvCalo_LHC17o-pass$passNr-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                periodList=`echo -e "c-pass1\ne-pass1\nf-pass1\nh-pass1\ni-pass1\nj-pass1\nk-pass1\nl-pass1\nm-pass1\no-pass1\nr-pass1"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC17$periodID-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            outAdd=`echo $periodID  | cut -d "-" -f 1 `
                            nameOut+=$outAdd
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        else
                            echo $currFile " does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC17$nameOut-pass$passNr-2-$runListName\_$number.root
                done
                #18er daten
                ls $OUTPUTDIR/GammaConvCalo_LHC18m-pass$passNr-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                periodList=`echo -e "b-pass1\nd-pass1\ne-pass1\nf-pass1\ng-pass1\nh-pass1\ni-pass1\nj-pass1\nk-pass1\nl-pass1\nm-pass1\nm-pass1\no-pass1\np-pass1"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC18$periodID-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            outAdd=`echo $periodID  | cut -d "-" -f 1 `
                            nameOut+=$outAdd
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        else
                            echo $currFile " does not exist"
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC18$nameOut-pass$passNr-2-$runListName\_$number.root
                done
            done
        fi

        echo "Merging Done"
    fi
else
    if [ $HAVELHC16d == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC16d";
        rm $OUTPUTDIR_LHC16d/*/GammaConvCalo_*.root
        rm $OUTPUTDIR_LHC16d/*/*/*GammaConvCalo_*.root
    fi
fi
