#! /bin/bash

# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEONDATA=0
MERGEONMC=0
MERGEONMCJJ=0
MERGEONBINSSingle=1
MERGEONBINS=1

# check if train configuration has actually been given


HAVELHC15n=1
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

HAVELHC15n=0
HAVELHC16d=0
HAVELHC16e=0
HAVELHC16f=0
HAVELHC16g=0
HAVELHC16h=0
HAVELHC16i=0
HAVELHC16j=0
HAVELHC16k=0
HAVELHC16l=0
HAVELHC16o=0
HAVELHC16p=0

#########################-----MonteCarlo files----------------################
#15n
HAVELHC17e2MC=1
#16d
HAVELHC17f6MC=1
#16e
HAVELHC17f9MC=1
#16f-lowB
HAVELHC17d12MC=1
#16f
HAVELHC17d1MC=1
#16g
HAVELHC17d17MC=1
#16h
HAVELHC17f5MC=1
#16i
HAVELHC17d3MC=1
#16j
HAVELHC17e5MC=1
#16k
HAVELHC17d20a1MC=1
HAVELHC17d20a1_extraMC=1
#16l?
HAVELHC17d20a2MC=1
HAVELHC17d20a2_extraMC=1
#16o
HAVELHC17d16MC=1
#16p
HAVELHC17d18MC=1

HAVELHC17f8bJJMC=1
HAVELHC17f8cJJMC=1
HAVELHC17f8dJJMC=1
HAVELHC17f8aJJMC=1
HAVELHC17f8eJJMC=1


#15n
HAVELHC17e2MC=0
#16d
HAVELHC17f6MC=0
#16e
HAVELHC17f9MC=0
#16f-lowB
HAVELHC17d12MC=0
#16f
HAVELHC17d1MC=0
#16g
HAVELHC17d17MC=0
#16h
HAVELHC17f5MC=0
#16i
HAVELHC17d3MC=0
#16j
HAVELHC17e5MC=0
#16k
HAVELHC17d20a1MC=0
HAVELHC17d20a1_extraMC=0
#16l?
HAVELHC17d20a2MC=0
HAVELHC17d20a2_extraMC=0
#16o
HAVELHC17d16MC=0
#16p
HAVELHC17d18MC=0

#HAVELHC17f8bJJMC=0
#HAVELHC17f8cJJMC=0
HAVELHC17f8dJJMC=0
HAVELHC17f8aJJMC=0
HAVELHC17f8eJJMC=0

# default trainconfigurations


LHC15nData="";
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


LHC17f6MC="";
LHC17f9MC="";
LHC17d12MC="";
LHC17d1MC="";
LHC17d17MC="";
LHC17f5MC="";
LHC17d3MC="";
LHC17e5MC="";
LHC17d20a1MC="";
LHC17d20a1_extraMC="";
LHC17d20a2MC="";
LHC17d20a2_extraMC="";
LHC17d16MC="";
LHC17d18MC="";
LHC17f8bJJMC="";
LHC17f8cJJMC="";
LHC17f8dJJMC="";
LHC17f8aJJMC="";
LHC17f8eJJMC="";



# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pp2760GeV needs
#NSlashes=9
#NSlashes2=11

NSlashes=11
NSlashes2=10
NSlashes3=11

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp
    NSlashes=8
    NSlashes2=7
elif [ $1 = "fbockGSI" ]; then
	BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp
elif [ $1 = "leardini" ]; then
	BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then
	BASEDIR=/alidata50/alice_u/leardini/GridOutput/pp/
elif [ $1 = "leardiniGSI" ]; then
	BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/pp
elif [ $1 = "passfeld" ]; then
	BASEDIR=~/work/Gridoutput/pp
elif [ $1 = "passfeldMAF" ]; then
	BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then
	BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pp
elif [ $1 = "amarin" ]; then
	BASEDIR=/Users/marin/analysis/2017/Grid/OutputLegoTrains/pp
elif [ $1 = "amarinGSI" ]; then
	BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/pp
elif [ $1 = "amarinALICESERV1" ]; then
	BASEDIR=/alidata50/alice_u/amarin/GridOutput/pp/
elif [ $1 = "mwilde" ]; then
	BASEDIR=~/work/GridOutput
elif [ $1 = "mwildeGSI" ]; then
	BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pp
elif [ $1 = "pgonzales" ]; then
	BASEDIR=~/work/GridOutput
elif [ $1 = "pgonzalesGSI" ]; then
	BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pp
fi

# TRAINDIR=Legotrain-vAN-20151025-MultiplicityDependence

#TRAINDIR=Legotrain-mCalo-20160813_SecEffiAndTMStudiesRerun
#LHC11aData="1782_20160813-2128";
#LHC12f1aMC="2426_20160813-2128";
#LHC12f1bMC="2427_20160813-2129";

#TRAINDIR=Legotrain-vAN-20170426-1-FirstLook
#LHC16fData="2060_20170426-1817";
#LHC16iData="2059_20170426-1819";
#LHC16kData="2058_20170426-1819";
#LHC16lData="2057_20170426-1819";
#LHC16oData="2056_20170426-1819";


TRAINDIR=Legotrain-vAN-20170710-13TeV
#LHC16fData="2060_20170426-1817";
#LHC16iData="2059_20170426-1819";
#LHC16kData="2058_20170426-1819";
#LHC16lData="2057_20170426-1819";
#LHC16oData="2056_20170426-1819";

LHC15nData="2184_20170711-2100";
#LHC15nDatapass3="2185_20170711-2100";
LHC16dData="2183_20170711-2100";
LHC16eData="2182_20170711-2059";
LHC16fData="2181_20170711-2059";
LHC16gData="2180_20170711-2059";
LHC16hData="2179_20170711-2059";
LHC16iData="2178_20170711-2059";
LHC16jData="2176_20170711-2058";
LHC16kData="2175_20170711-2058";
LHC16lData="2174_20170711-2058";
LHC16oData="2173_20170711-2058";
LHC16pData="2172_20170711-2100";

#15n-pass4
LHC17e2MC="3057_20170714-1443";
LHC17f6MC="3054_20170714-1442";
#LHC17f9MC="";
#lowBfield
#LHC17d12MC="";
LHC17d1MC="3063_20170714-1444";
LHC17d17MC="3064_20170714-1444";
LHC17f5MC="3055_20170714-1443";
LHC17d3MC="3058_20170714-1444";
LHC17e5MC="3056_20170714-1443";
LHC17d20a1MC="3061_20170714-1444";
LHC17d20a1_extraMC="3062_20170714-1444";
LHC17d20a2MC="3059_20170714-1444";
LHC17d20a2_extraMC="3060_20170714-1444";
LHC17d16MC="3065_20170714-1445";
LHC17d18MC="3066_20170714-1445";


#LHC17f8bJJMC="3052_20170714-1643";
#LHC17f8cJJMC="3051_20170714-1643";
#LHC17f8dJJMC="3050_20170714-1643";
#LHC17f8aJJMC="3053_20170714-1642";
#LHC17f8eJJMC="3049_20170714-1643";
#-----JetJet Slow train for ntrials
LHC17f8bJJMC="3079_20170725-1138";
LHC17f8cJJMC="3078_20170725-1138";
LHC17f8dJJMC="3077_20170725-1138";
LHC17f8aJJMC="3080_20170725-1138";
LHC17f8eJJMC="3076_20170725-1138";



#DATAADD=""

OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ $2 == "LHC16f" ]; then

#    if [ "$LHC16fData" == "" ]; then
#        HAVELHC16f=0;
#    fi

    if [ $HAVELHC15n == 1 ]; then
	echo "outputdir LHC15n"
        OUTPUTDIR_LHC15n=$BASEDIR/$TRAINDIR/GA_pp-$LHC15nData
        mkdir -p $OUTPUTDIR_LHC15n
    fi


    if [ $HAVELHC16d == 1 ]; then
	echo "outputdir LHC16d"
        OUTPUTDIR_LHC16d=$BASEDIR/$TRAINDIR/GA_pp-$LHC16dData
        mkdir -p $OUTPUTDIR_LHC16d
    fi

    if [ $HAVELHC16e == 1 ]; then
	echo "outputdir LHC16e"
        OUTPUTDIR_LHC16e=$BASEDIR/$TRAINDIR/GA_pp-$LHC16eData
        mkdir -p $OUTPUTDIR_LHC16e
    fi

    if [ $HAVELHC16f == 1 ]; then
	echo "outputdir LHC16f"
        OUTPUTDIR_LHC16f=$BASEDIR/$TRAINDIR/GA_pp-$LHC16fData
        mkdir -p $OUTPUTDIR_LHC16f
    fi

    if [ $HAVELHC16g == 1 ]; then
	echo "outputdir LHC16g"
        OUTPUTDIR_LHC16g=$BASEDIR/$TRAINDIR/GA_pp-$LHC16gData
        mkdir -p $OUTPUTDIR_LHC16g
    fi

    if [ $HAVELHC16h == 1 ]; then
	echo "outputdir LHC16h"
        OUTPUTDIR_LHC16h=$BASEDIR/$TRAINDIR/GA_pp-$LHC16hData
        mkdir -p $OUTPUTDIR_LHC16h
    fi


    if [ $HAVELHC16i == 1 ]; then
	echo "outputdir LHC16i"
        OUTPUTDIR_LHC16i=$BASEDIR/$TRAINDIR/GA_pp-$LHC16iData
        mkdir -p $OUTPUTDIR_LHC16i
    fi

    if [ $HAVELHC16j == 1 ]; then
	echo "outputdir LHC16j"
        OUTPUTDIR_LHC16j=$BASEDIR/$TRAINDIR/GA_pp-$LHC16jData
        mkdir -p $OUTPUTDIR_LHC16j
    fi

    if [ $HAVELHC16k == 1 ]; then
	echo "outputdir LHC16k"
        OUTPUTDIR_LHC16k=$BASEDIR/$TRAINDIR/GA_pp-$LHC16kData
        mkdir -p $OUTPUTDIR_LHC16k
    fi

    if [ $HAVELHC16l == 1 ]; then
	echo "outputdir LHC16l"
        OUTPUTDIR_LHC16l=$BASEDIR/$TRAINDIR/GA_pp-$LHC16lData
        mkdir -p $OUTPUTDIR_LHC16l
    fi

    if [ $HAVELHC16o == 1 ]; then
	echo "outputdir LHC16o"
        OUTPUTDIR_LHC16o=$BASEDIR/$TRAINDIR/GA_pp-$LHC16oData
        mkdir -p $OUTPUTDIR_LHC16o
    fi

    if [ $HAVELHC16p == 1 ]; then
	echo "outputdir LHC16p"
        OUTPUTDIR_LHC16p=$BASEDIR/$TRAINDIR/GA_pp-$LHC16pData
        mkdir -p $OUTPUTDIR_LHC16p
    fi


########################################################
############going to create directories for MC##########
########################################################

    if [ $HAVELHC17e2MC == 1 ]; then
	echo "outputdir LHC17e2MC"
        OUTPUTDIR_LHC17e2MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17e2MC
	mkdir -p $OUTPUTDIR_LHC17e2MC
    fi


    if [ $HAVELHC17f6MC == 1 ]; then
	echo "outputdir LHC17f6MC"
        OUTPUTDIR_LHC17f6MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f6MC
	mkdir -p $OUTPUTDIR_LHC17f6MC
    fi

    if [ $HAVELHC17d1MC == 1 ]; then
	echo "outputdir LHC17d1MC"
        OUTPUTDIR_LHC17d1MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d1MC
	mkdir -p $OUTPUTDIR_LHC17d1MC
    fi

    if [ $HAVELHC17d17MC == 1 ]; then
	echo "outputdir LHC17d17MC"
        OUTPUTDIR_LHC17d17MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d17MC
	mkdir -p $OUTPUTDIR_LHC17d17MC
    fi

    if [ $HAVELHC17f5MC == 1 ]; then
	echo "outputdir LHC17f5MC"
        OUTPUTDIR_LHC17f5MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f5MC
	mkdir -p $OUTPUTDIR_LHC17f5MC
    fi

    if [ $HAVELHC17d3MC == 1 ]; then
	echo "outputdir LHC17d3MC"
        OUTPUTDIR_LHC17d3MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d3MC
	mkdir -p $OUTPUTDIR_LHC17d3MC
    fi


    if [ $HAVELHC17e5MC == 1 ]; then
	echo "outputdir LHC17e5MC"
        OUTPUTDIR_LHC17e5MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17e5MC
	mkdir -p $OUTPUTDIR_LHC17e5MC
    fi

    if [ $HAVELHC17d20a1MC == 1 ]; then
	echo "outputdir LHC17d20a1MC"
        OUTPUTDIR_LHC17d20a1MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a1MC
	mkdir -p $OUTPUTDIR_LHC17d20a1MC
    fi

    if [ $HAVELHC17d20a2MC == 1 ]; then
	echo "outputdir LHC17d20a2MC"
        OUTPUTDIR_LHC17d20a2MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a2MC
	mkdir -p $OUTPUTDIR_LHC17d20a2MC
    fi

    if [ $HAVELHC17d20a1_extraMC == 1 ]; then
	echo "outputdir LHC17d20a1_extraMC"
        OUTPUTDIR_LHC17d20a1_extraMC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a1_extraMC
	mkdir -p $OUTPUTDIR_LHC17d20a1_extraMC
    fi

    if [ $HAVELHC17d20a2_extraMC == 1 ]; then
	echo "outputdir LHC17d20a2_extraMC"
        OUTPUTDIR_LHC17d20a2_extraMC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d20a2_extraMC
	mkdir -p $OUTPUTDIR_LHC17d20a2_extraMC
    fi


    if [ $HAVELHC17d16MC == 1 ]; then
	echo "outputdir LHC17d16MC"
        OUTPUTDIR_LHC17d16MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d16MC
	mkdir -p $OUTPUTDIR_LHC17d16MC
    fi

    if [ $HAVELHC17d18MC == 1 ]; then
	echo "outputdir LHC17d18MC"
        OUTPUTDIR_LHC17d18MC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17d18MC
	mkdir -p $OUTPUTDIR_LHC17d18MC
    fi

    if [ $HAVELHC17f8bJJMC == 1 ]; then
	echo "outputdir LHC17f8bJJMC"
        OUTPUTDIR_LHC17f8bJJMC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f8bJJMC
	mkdir -p $OUTPUTDIR_LHC17f8bJJMC
    fi

    if [ $HAVELHC17f8cJJMC == 1 ]; then
	echo "outputdir LHC17f8cJJMC"
        OUTPUTDIR_LHC17f8cJJMC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f8cJJMC
	mkdir -p $OUTPUTDIR_LHC17f8cJJMC
    fi

    if [ $HAVELHC17f8dJJMC == 1 ]; then
	echo "outputdir LHC17f8dJJMC"
        OUTPUTDIR_LHC17f8dJJMC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f8dJJMC
	mkdir -p $OUTPUTDIR_LHC17f8dJJMC
    fi

    if [ $HAVELHC17f8aJJMC == 1 ]; then
	echo "outputdir LHC17f8aJJMC"
        OUTPUTDIR_LHC17f8aJJMC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f8aJJMC
	mkdir -p $OUTPUTDIR_LHC17f8aJJMC
    fi

    if [ $HAVELHC17f8eJJMC == 1 ]; then
	echo "outputdir LHC17f8eJJMC"
        OUTPUTDIR_LHC17f8eJJMC=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC17f8eJJMC
	mkdir -p $OUTPUTDIR_LHC17f8eJJMC
    fi


    if [ $DOWNLOADON == 1 ]; then

        if [ $HAVELHC15n == 1 ]; then
            echo "downloading LHC15n"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC15nData/merge_runlist_3/root_archive.zip file:$OUTPUTDIR_LHC15n/
            unzip -u $OUTPUTDIR_LHC15n/root_archive.zip -d $OUTPUTDIR_LHC15n/
        fi

        if [ $HAVELHC16d == 1 ]; then
            echo "downloading LHC16d"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16dData/merge/root_archive.zip file:$OUTPUTDIR_LHC16d/
            unzip -u $OUTPUTDIR_LHC16d/root_archive.zip -d $OUTPUTDIR_LHC16d/
        fi

        if [ $HAVELHC16e == 1 ]; then
            echo "downloading LHC16e"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16eData/merge/root_archive.zip file:$OUTPUTDIR_LHC16e/
            unzip -u $OUTPUTDIR_LHC16e/root_archive.zip -d $OUTPUTDIR_LHC16e/
        fi


        if [ $HAVELHC16f == 1 ]; then
            echo "downloading LHC16f"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16fData/merge/root_archive.zip file:$OUTPUTDIR_LHC16f/
            unzip -u $OUTPUTDIR_LHC16f/root_archive.zip -d $OUTPUTDIR_LHC16f/
        fi

        if [ $HAVELHC16g == 1 ]; then
            echo "downloading LHC16g"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16gData/merge/root_archive.zip file:$OUTPUTDIR_LHC16g/
            unzip -u $OUTPUTDIR_LHC16g/root_archive.zip -d $OUTPUTDIR_LHC16g/
        fi

        if [ $HAVELHC16h == 1 ]; then
            echo "downloading LHC16h"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16hData/merge/root_archive.zip file:$OUTPUTDIR_LHC16h/
            unzip -u $OUTPUTDIR_LHC16h/root_archive.zip -d $OUTPUTDIR_LHC16h/
        fi

        if [ $HAVELHC16i == 1 ]; then
            echo "downloading LHC16i"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16iData/merge/root_archive.zip file:$OUTPUTDIR_LHC16i/
            unzip -u $OUTPUTDIR_LHC16i/root_archive.zip -d $OUTPUTDIR_LHC16i/
        fi

        if [ $HAVELHC16j == 1 ]; then
            echo "downloading LHC16j"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16jData/merge/root_archive.zip file:$OUTPUTDIR_LHC16j/
            unzip -u $OUTPUTDIR_LHC16j/root_archive.zip -d $OUTPUTDIR_LHC16j/
        fi



	if [ $HAVELHC16k == 1 ]; then
            echo "downloading LHC16k"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16kData/merge/root_archive.zip file:$OUTPUTDIR_LHC16k/
            unzip -u $OUTPUTDIR_LHC16k/root_archive.zip -d $OUTPUTDIR_LHC16k/
        fi

	if [ $HAVELHC16l == 1 ]; then
            echo "downloading LHC16l"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16lData/merge_runlist_3/root_archive.zip file:$OUTPUTDIR_LHC16l/
            unzip -u $OUTPUTDIR_LHC16l/root_archive.zip -d $OUTPUTDIR_LHC16l/
        fi

	if [ $HAVELHC16o == 1 ]; then
            echo "downloading LHC16o"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16oData/merge/root_archive.zip file:$OUTPUTDIR_LHC16o/
            unzip -u $OUTPUTDIR_LHC16o/root_archive.zip -d $OUTPUTDIR_LHC16o/
        fi

        if [ $HAVELHC16p == 1 ]; then
            echo "downloading LHC16p"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC16pData/merge/root_archive.zip file:$OUTPUTDIR_LHC16p/
            unzip -u $OUTPUTDIR_LHC16p/root_archive.zip -d $OUTPUTDIR_LHC16p/
        fi
        ########--------MC outputs--------#

	if [ $HAVELHC17e2MC == 1 ]; then
            echo "downloading LHC17e2"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17e2MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17e2MC/
            unzip -u $OUTPUTDIR_LHC17e2MC/root_archive.zip -d $OUTPUTDIR_LHC17e2MC/
        fi


	if [ $HAVELHC17f6MC == 1 ]; then
            echo "downloading LHC17f6"
#            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f6MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17f6MC/
	    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f6MC/merge/GammaConvV1_300.root file:$OUTPUTDIR_LHC17f6MC/
	    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f6MC/merge/GammaConvV1_301.root file:$OUTPUTDIR_LHC17f6MC/
	    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f6MC/merge/GammaConvV1_302.root file:$OUTPUTDIR_LHC17f6MC/
	    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f6MC/merge/GammaConv_Material_21.root file:$OUTPUTDIR_LHC17f6MC/
	    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f6MC/merge/GammaConv_Material_121.root file:$OUTPUTDIR_LHC17f6MC/
#            unzip -u $OUTPUTDIR_LHC17f6MC/root_archive.zip -d $OUTPUTDIR_LHC17f6MC/
        fi

	if [ $HAVELHC17d1MC == 1 ]; then
            echo "downloading LHC17d1"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d1MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d1MC/
            unzip -u $OUTPUTDIR_LHC17d1MC/root_archive.zip -d $OUTPUTDIR_LHC17d1MC/
        fi

	if [ $HAVELHC17d17MC == 1 ]; then
            echo "downloading LHC17d17"
 #           alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d17MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d17MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d17MC/merge/GammaConvV1_300.root file:$OUTPUTDIR_LHC17d17MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d17MC/merge/GammaConvV1_301.root file:$OUTPUTDIR_LHC17d17MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d17MC/merge/GammaConvV1_302.root file:$OUTPUTDIR_LHC17d17MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d17MC/merge/GammaConv_Material_21.root file:$OUTPUTDIR_LHC17d17MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d17MC/merge/GammaConv_Material_121.root file:$OUTPUTDIR_LHC17d17MC/
 #           unzip -u $OUTPUTDIR_LHC17d17MC/root_archive.zip -d $OUTPUTDIR_LHC17d17MC/
        fi

	if [ $HAVELHC17f5MC == 1 ]; then
            echo "downloading LHC17f5"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f5MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17f5MC/
            unzip -u $OUTPUTDIR_LHC17f5MC/root_archive.zip -d $OUTPUTDIR_LHC17f5MC/
        fi

	if [ $HAVELHC17d3MC == 1 ]; then
            echo "downloading LHC17d3"
 #           alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d3MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d3MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d3MC/merge/GammaConvV1_300.root file:$OUTPUTDIR_LHC17d3MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d3MC/merge/GammaConvV1_301.root file:$OUTPUTDIR_LHC17d3MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d3MC/merge/GammaConvV1_302.root file:$OUTPUTDIR_LHC17d3MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d3MC/merge/GammaConv_Material_21.root file:$OUTPUTDIR_LHC17d3MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d3MC/merge/GammaConv_Material_121.root file:$OUTPUTDIR_LHC17d3MC/
 #           unzip -u $OUTPUTDIR_LHC17d3MC/root_archive.zip -d $OUTPUTDIR_LHC17d3MC/
        fi

	if [ $HAVELHC17e5MC == 1 ]; then
            echo "downloading LHC17e5"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17e5MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17e5MC/
            unzip -u $OUTPUTDIR_LHC17e5MC/root_archive.zip -d $OUTPUTDIR_LHC17e5MC/
        fi

	if [ $HAVELHC17d20a1MC == 1 ]; then
            echo "downloading LHC17d20a1"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a1MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d20a1MC/
            unzip -u $OUTPUTDIR_LHC17d20a1MC/root_archive.zip -d $OUTPUTDIR_LHC17d20a1MC/
        fi

	if [ $HAVELHC17d20a1_extraMC == 1 ]; then
            echo "downloading LHC17d20a1_extra"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a1_extraMC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d20a1_extraMC/
            unzip -u $OUTPUTDIR_LHC17d20a1_extraMC/root_archive.zip -d $OUTPUTDIR_LHC17d20a1_extraMC/
        fi

	if [ $HAVELHC17d20a2MC == 1 ]; then
            echo "downloading LHC17d20a2"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a2MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d20a2MC/
            unzip -u $OUTPUTDIR_LHC17d20a2MC/root_archive.zip -d $OUTPUTDIR_LHC17d20a2MC/
        fi

	if [ $HAVELHC17d20a2_extraMC == 1 ]; then
            echo "downloading LHC17d20a2_extra"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d20a2_extraMC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d20a2_extraMC/
            unzip -u $OUTPUTDIR_LHC17d20a2_extraMC/root_archive.zip -d $OUTPUTDIR_LHC17d20a2_extraMC/
        fi

	if [ $HAVELHC17d16MC == 1 ]; then
            echo "downloading LHC17d16"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d16MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d16MC/
            unzip -u $OUTPUTDIR_LHC17d16MC/root_archive.zip -d $OUTPUTDIR_LHC17d16MC/
        fi

	if [ $HAVELHC17d18MC == 1 ]; then
            echo "downloading LHC17d18"
 #           alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d18MC/merge/root_archive.zip file:$OUTPUTDIR_LHC17d18MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d18MC/merge/GammaConvV1_300.root file:$OUTPUTDIR_LHC17d18MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d18MC/merge/GammaConvV1_301.root file:$OUTPUTDIR_LHC17d18MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d18MC/merge/GammaConvV1_302.root file:$OUTPUTDIR_LHC17d18MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d18MC/merge/GammaConv_Material_21.root file:$OUTPUTDIR_LHC17d18MC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17d18MC/merge/GammaConv_Material_121.root file:$OUTPUTDIR_LHC17d18MC/
 #           unzip -u $OUTPUTDIR_LHC17d18MC/root_archive.zip -d $OUTPUTDIR_LHC17d18MC/
        fi

	if [ $HAVELHC17f8bJJMC == 1 ]; then
            echo "downloading LHC17f8bJJ"
 #           alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8bJJMC/merge/root_archive.zip file:$OUTPUTDIR_LHC17f8bJJMC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8bJJMC/merge/GammaConvV1_300.root  file:$OUTPUTDIR_LHC17f8bJJMC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8bJJMC/merge/GammaConvV1_301.root  file:$OUTPUTDIR_LHC17f8bJJMC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8bJJMC/merge/GammaConvV1_302.root  file:$OUTPUTDIR_LHC17f8bJJMC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8bJJMC/merge/GammaConv_Material_21.root  file:$OUTPUTDIR_LHC17f8bJJMC/
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8bJJMC/merge/GammaConv_Material_121.root  file:$OUTPUTDIR_LHC17f8bJJMC/
 #           unzip -u $OUTPUTDIR_LHC17f8bJJMC/root_archive.zip -d $OUTPUTDIR_LHC17f8bJJMC/
        fi

	if [ $HAVELHC17f8cJJMC == 1 ]; then
            echo "downloading LHC17f8cJJ"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8cJJMC/merge/root_archive.zip file:$OUTPUTDIR_LHC17f8cJJMC/
            unzip -u $OUTPUTDIR_LHC17f8cJJMC/root_archive.zip -d $OUTPUTDIR_LHC17f8cJJMC/
        fi

	if [ $HAVELHC17f8dJJMC == 1 ]; then
            echo "downloading LHC17f8dJJ"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8dJJMC/merge/root_archive.zip file:$OUTPUTDIR_LHC17f8dJJMC/
            unzip -u $OUTPUTDIR_LHC17f8dJJMC/root_archive.zip -d $OUTPUTDIR_LHC17f8dJJMC/
        fi

	if [ $HAVELHC17f8aJJMC == 1 ]; then
            echo "downloading LHC17f8aJJ"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8aJJMC/merge/root_archive.zip file:$OUTPUTDIR_LHC17f8aJJMC/
            unzip -u $OUTPUTDIR_LHC17f8aJJMC/root_archive.zip -d $OUTPUTDIR_LHC17f8aJJMC/
        fi

	if [ $HAVELHC17f8eJJMC == 1 ]; then
            echo "downloading LHC17f8eJJ"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC17f8eJJMC/merge/root_archive.zip file:$OUTPUTDIR_LHC17f8eJJMC/
            unzip -u $OUTPUTDIR_LHC17f8eJJMC/root_archive.zip -d $OUTPUTDIR_LHC17f8eJJMC/
        fi

        echo "copying LHC17f8cJJ MC bins (16f)"
        runNumbers=`cat runlists/runNumbersLHC17f8cJJ.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat runlists/binNumbersJJ_13TeV.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC17f8cJJMC == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC17f8cJJMC/$binNumber/$runNumber
                        alien_cp alien:/alice/sim/2017/LHC17f8c/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC17f8cJJMC/root_archive.zip file:$OUTPUTDIR_LHC17f8cJJMC/$binNumber/$runNumber/
                        unzip -u $OUTPUTDIR_LHC17f8cJJMC/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC17f8cJJMC/$binNumber/$runNumber/
                    fi
                done;
            fi
        done;


        echo "copying LHC17f8bJJ MC bins (16f)"
        runNumbers=`cat runlists/runNumbersLHC17f8bJJ.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat runlists/binNumbersJJ_13TeV.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC17f8bJJMC == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC17f8bJJMC/$binNumber/$runNumber
                        alien_cp alien:/alice/sim/2017/LHC17f8b/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC17f8bJJMC/root_archive.zip file:$OUTPUTDIR_LHC17f8bJJMC/$binNumber/$runNumber/
                        unzip -u $OUTPUTDIR_LHC17f8bJJMC/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC17f8bJJMC/$binNumber/$runNumber/
                    fi
                done;
            fi
        done;


       echo "copying LHC17f8dJJ MC bins (16)"
        runNumbers=`cat runlists/runNumbersLHC17f8dJJ.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat runlists/binNumbersJJ_13TeV.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC17f8dJJMC == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC17f8dJJMC/$binNumber/$runNumber
                        alien_cp alien:/alice/sim/2017/LHC17f8d/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC17f8dJJMC/root_archive.zip file:$OUTPUTDIR_LHC17f8dJJMC/$binNumber/$runNumber/
                        unzip -u $OUTPUTDIR_LHC17f8dJJMC/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC17f8dJJMC/$binNumber/$runNumber/
                    fi
                done;
            fi
        done;


      echo "copying LHC17f8aJJ MC bins (16)"
        runNumbers=`cat runlists/runNumbersLHC17f8aJJ.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat runlists/binNumbersJJ_13TeV.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC17f8aJJMC == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC17f8aJJMC/$binNumber/$runNumber
                        alien_cp alien:/alice/sim/2017/LHC17f8d/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC17f8aJJMC/root_archive.zip file:$OUTPUTDIR_LHC17f8aJJMC/$binNumber/$runNumber/
                        unzip -u $OUTPUTDIR_LHC17f8aJJMC/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC17f8aJJMC/$binNumber/$runNumber/
                    fi
                done;
            fi
        done;

      echo "copying LHC17f8eJJ MC bins (16)"
        runNumbers=`cat runlists/runNumbersLHC17f8eJJ.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat runlists/binNumbersJJ_13TeV.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC17f8eJJMC == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC17f8eJJMC/$binNumber/$runNumber
                        alien_cp alien:/alice/sim/2017/LHC17f8d/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC17f8eJJMC/root_archive.zip file:$OUTPUTDIR_LHC17f8eJJMC/$binNumber/$runNumber/
                        unzip -u $OUTPUTDIR_LHC17f8eJJMC/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC17f8eJJMC/$binNumber/$runNumber/
                    fi
                done;
            fi
        done;



    fi
#################
    if [ $HAVELHC15n == 1 ]; then
        ls $OUTPUTDIR_LHC15n/GammaConvV1_*.root > fileLHC15n.txt
        fileNumbers=`cat fileLHC15n.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC15n/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC15n\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC15n\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC15n_$number.log\"\)
        done;
    fi

    if [ $HAVELHC16d == 1 ]; then
        ls $OUTPUTDIR_LHC16d/GammaConvV1_*.root > fileLHC16d.txt
        fileNumbers=`cat fileLHC16d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16d/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16d\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16d\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16d_$number.log\"\)
        done;
    fi

    if [ $HAVELHC16e == 1 ]; then
        ls $OUTPUTDIR_LHC16e/GammaConvV1_*.root > fileLHC16e.txt
        fileNumbers=`cat fileLHC16e.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16e/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16e\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16e\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16e_$number.log\"\)
        done;
    fi


    if [ $HAVELHC16f == 1 ]; then
        ls $OUTPUTDIR_LHC16f/GammaConvV1_*.root > fileLHC16f.txt
        fileNumbers=`cat fileLHC16f.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16f/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16f\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16f\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16f_$number.log\"\)
        done;
    fi


    if [ $HAVELHC16g == 1 ]; then
        ls $OUTPUTDIR_LHC16f/GammaConvV1_*.root > fileLHC16g.txt
        fileNumbers=`cat fileLHC16g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16f/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16g\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16g\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16g_$number.log\"\)
        done;
    fi

    if [ $HAVELHC16h == 1 ]; then
        ls $OUTPUTDIR_LHC16h/GammaConvV1_*.root > fileLHC16h.txt
        fileNumbers=`cat fileLHC16h.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16h/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16h\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16h\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16h_$number.log\"\)
        done;
    fi


    if [ $HAVELHC16i == 1 ]; then
        ls $OUTPUTDIR_LHC16i/GammaConvV1_*.root > fileLHC16i.txt
        fileNumbers=`cat fileLHC16i.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16i/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16i\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16i\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16i_$number.log\"\)
        done;
    fi

    if [ $HAVELHC16j == 1 ]; then
        ls $OUTPUTDIR_LHC16j/GammaConvV1_*.root > fileLHC16j.txt
        fileNumbers=`cat fileLHC16j.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16j/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16j\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16j\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16j_$number.log\"\)
        done;
    fi




    if [ $HAVELHC16k == 1 ]; then
        ls $OUTPUTDIR_LHC16k/GammaConvV1_*.root > fileLHC16k.txt
        fileNumbers=`cat fileLHC16k.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16k/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16k\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16k\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16k_$number.log\"\)
        done;
    fi




    if [ $HAVELHC16l == 1 ]; then
        ls $OUTPUTDIR_LHC16l/GammaConvV1_*.root > fileLHC16l.txt
        fileNumbers=`cat fileLHC16l.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16l/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16l\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16l\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16l_$number.log\"\)
        done;
    fi

#    if [ $HAVELHC16m == 1 ]; then
#        ls $OUTPUTDIR_LHC16m/GammaConvV1_*.root > fileLHC16m.txt
#        fileNumbers=`cat fileLHC16m.txt`
#        for fileName in $fileNumbers; do
#            echo $fileName
#            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#            echo $number
#            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16m/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16m\_$number.root\"\,\"GammaConvV1_$number\"\)
#            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16m\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16m_$number.log\"\)
#        done;
#    fi



    if [ $HAVELHC16o == 1 ]; then
        ls $OUTPUTDIR_LHC16o/GammaConvV1_*.root > fileLHC16o.txt
        fileNumbers=`cat fileLHC16o.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16o/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16o\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16o\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16o_$number.log\"\)
        done;
    fi

    if [ $HAVELHC16p == 1 ]; then
        ls $OUTPUTDIR_LHC16p/GammaConvV1_*.root > fileLHC16p.txt
        fileNumbers=`cat fileLHC16i.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16p/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC16p\_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC16p\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC16p_$number.log\"\)
        done;
    fi

#############################################

   if [ $HAVELHC17e2MC == 1 ]; then
        ls $OUTPUTDIR_LHC17e2MC/GammaConvV1_*.root > fileLHC17e2MC.txt
        fileNumbers=`cat fileLHC17e2MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17e2MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17e2_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17e2_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17e2_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17f6MC == 1 ]; then
        ls $OUTPUTDIR_LHC17f6MC/GammaConvV1_*.root > fileLHC17f6MC.txt
        fileNumbers=`cat fileLHC17f6MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f6MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f6_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17f6_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17f6_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17d1MC == 1 ]; then
        ls $OUTPUTDIR_LHC17d1MC/GammaConvV1_*.root > fileLHC17d1MC.txt
        fileNumbers=`cat fileLHC17d1MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d1MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d1_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17d17MC == 1 ]; then
        ls $OUTPUTDIR_LHC17d17MC/GammaConvV1_*.root > fileLHC17d17MC.txt
        fileNumbers=`cat fileLHC17d17MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d17MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d17_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d17_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d17_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17f5MC == 1 ]; then
        ls $OUTPUTDIR_LHC17f5MC/GammaConvV1_*.root > fileLHC17f5MC.txt
        fileNumbers=`cat fileLHC17f5MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f5MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f5_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17f5_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17f5_$number.log\"\)
        done;
    fi


    if [ $HAVELHC17d3MC == 1 ]; then
        ls $OUTPUTDIR_LHC17d3MC/GammaConvV1_*.root > fileLHC17d3MC.txt
        fileNumbers=`cat fileLHC17d3MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d3MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d3_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d3_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d3_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17e5MC == 1 ]; then
        ls $OUTPUTDIR_LHC17e5MC/GammaConvV1_*.root > fileLHC17e5MC.txt
        fileNumbers=`cat fileLHC17e5MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17e5MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17e5_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17e5_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17e5_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17d20a1MC == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a1MC/GammaConvV1_*.root > fileLHC17d20a1MC.txt
        fileNumbers=`cat fileLHC17d20a1MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d20a1MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d20a1_$number.log\"\)
        done;
    fi
    if [ $HAVELHC17d20a1_extraMC == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a1_extraMC/GammaConvV1_*.root > fileLHC17d20a1_extraMC.txt
        fileNumbers=`cat fileLHC17d20a1_extraMC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d20a1_extraMC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_extra_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_extra_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d20a1_extra_$number.log\"\)
        done;
    fi

   if [ $HAVELHC17d20a2MC == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a2MC/GammaConvV1_*.root > fileLHC17d20a2MC.txt
        fileNumbers=`cat fileLHC17d20a2MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d20a2MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d20a2_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17d20a2_extraMC == 1 ]; then
        ls $OUTPUTDIR_LHC17d20a2_extraMC/GammaConvV1_*.root > fileLHC17d20a2_extraMC.txt
        fileNumbers=`cat fileLHC17d20a2_extraMC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d20a2_extraMC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_extra_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_extra_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d20a2_extra_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17d16MC == 1 ]; then
        ls $OUTPUTDIR_LHC17d16MC/GammaConvV1_*.root > fileLHC17d16MC.txt
        fileNumbers=`cat fileLHC17d16MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d16MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d16_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d16_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d16_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17d18MC == 1 ]; then
        ls $OUTPUTDIR_LHC17d18MC/GammaConvV1_*.root > fileLHC17d18MC.txt
        fileNumbers=`cat fileLHC17d18MC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17d18MC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17d18_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17d18_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17d18_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17f8bJJMC == 1 ]; then
        ls $OUTPUTDIR_LHC17f8bJJMC/GammaConvV1_*.root > fileLHC17f8bJJMC.txt
        fileNumbers=`cat fileLHC17f8bJJMC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8bJJMC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17f8bJJ_$number.log\"\)
        done;
    fi

    if [ $HAVELHC17f8cJJMC == 1 ]; then
        ls $OUTPUTDIR_LHC17f8cJJMC/GammaConvV1_*.root > fileLHC17f8cJJMC.txt
        fileNumbers=`cat fileLHC17f8cJJMC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8cJJMC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17f8cJJ_$number.log\"\)
        done;
    fi


   if [ $HAVELHC17f8dJJMC == 1 ]; then
        ls $OUTPUTDIR_LHC17f8dJJMC/GammaConvV1_*.root > fileLHC17f8dJJMC.txt
        fileNumbers=`cat fileLHC17f8dJJMC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8dJJMC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17f8dJJ_$number.log\"\)
        done;
    fi


   if [ $HAVELHC17f8aJJMC == 1 ]; then
        ls $OUTPUTDIR_LHC17f8aJJMC/GammaConvV1_*.root > fileLHC17f8aJJMC.txt
        fileNumbers=`cat fileLHC17f8aJJMC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8aJJMC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17f8aJJ_$number.log\"\)
        done;
    fi


   if [ $HAVELHC17f8eJJMC == 1 ]; then
        ls $OUTPUTDIR_LHC17f8eJJMC/GammaConvV1_*.root > fileLHC17f8eJJMC.txt
        fileNumbers=`cat fileLHC17f8cJJMC.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8eJJMC/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17f8eJJ_$number.log\"\)
        done;
    fi

######---------JetJet MonteCarlos----------------######

    if [ $MERGEONBINSSingle = 1 ]; then
        binNumbersJJ=`cat runlists/binNumbersJJ_13TeV.txt`
        echo $binNumbersJJ
        ls $OUTPUTDIR_LHC17f8bJJMC/GammaConvV1_*.root > filetemp.txt
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC17f8bJJMC == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC17f8bJJMC/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC17f8bJJMC/$binNumber/*/GammaConvV1_$number.root
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8bJJMC/$binNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ$binNumber\_$number.root\"\,\"GammaConvV1_$number\"\)
                fi
            done;
        done;


	ls $OUTPUTDIR_LHC17f8cJJMC/GammaConvV1_*.root > filetemp.txt
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC17f8cJJMC == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC17f8cJJMC/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC17f8cJJMC/$binNumber/*/GammaConvV1_$number.root
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8cJJMC/$binNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ$binNumber\_$number.root\"\,\"GammaConvV1_$number\"\)
                fi
            done;
        done;


       ls $OUTPUTDIR_LHC17f8dJJMC/GammaConvV1_*.root > filetemp.txt
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC17f8dJJMC == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC17f8dJJMC/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC17f8dJJMC/$binNumber/*/GammaConvV1_$number.root
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8dJJMC/$binNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ$binNumber\_$number.root\"\,\"GammaConvV1_$number\"\)
                fi
            done;
        done;



       ls $OUTPUTDIR_LHC17f8aJJMC/GammaConvV1_*.root > filetemp.txt
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC17f8aJJMC == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC17f8aJJMC/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC17f8aJJMC/$binNumber/*/GammaConvV1_$number.root
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8aJJMC/$binNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ$binNumber\_$number.root\"\,\"GammaConvV1_$number\"\)
                fi
            done;
        done;



       ls $OUTPUTDIR_LHC17f8eJJMC/GammaConvV1_*.root > filetemp.txt
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC17f8eJJMC == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC17f8eJJMC/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC17f8eJJMC/$binNumber/*/GammaConvV1_$number.root
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17f8eJJMC/$binNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ$binNumber\_$number.root\"\,\"GammaConvV1_$number\"\)
                fi
            done;
        done;
    fi



    ################################# Merging files###############

    if [ $MERGEONDATA == 1 ]; then
	echo " Merging files"
#	rm $OUTPUTDIR/GammaConvV1_LHC16defghijklop_*.root
	ls $OUTPUTDIR/GammaConvV1_LHC16d_*.root  > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_LHC16d_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16e_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16de_$number.root $OUTPUTDIR/GammaConvV1_LHC16d_$number.root $OUTPUTDIR/GammaConvV1_LHC16e_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_LHC16f_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16g_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16fg_$number.root $OUTPUTDIR/GammaConvV1_LHC16f_$number.root $OUTPUTDIR/GammaConvV1_LHC16g_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_LHC16h_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16i_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16hi_$number.root $OUTPUTDIR/GammaConvV1_LHC16h_$number.root $OUTPUTDIR/GammaConvV1_LHC16i_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_LHC16j_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16k_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16jk_$number.root $OUTPUTDIR/GammaConvV1_LHC16j_$number.root $OUTPUTDIR/GammaConvV1_LHC16k_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_LHC16l_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16o_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16lo_$number.root $OUTPUTDIR/GammaConvV1_LHC16l_$number.root $OUTPUTDIR/GammaConvV1_LHC16o_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_LHC16de_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16fg_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16defg_$number.root $OUTPUTDIR/GammaConvV1_LHC16de_$number.root $OUTPUTDIR/GammaConvV1_LHC16fg_$number.root
            fi

	    if [ -f $OUTPUTDIR/GammaConvV1_LHC16hi_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16jk_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16hijk_$number.root $OUTPUTDIR/GammaConvV1_LHC16hi_$number.root $OUTPUTDIR/GammaConvV1_LHC16jk_$number.root
            fi

	    if [ -f $OUTPUTDIR/GammaConvV1_LHC16lo_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16p_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16lop_$number.root $OUTPUTDIR/GammaConvV1_LHC16lo_$number.root $OUTPUTDIR/GammaConvV1_LHC16p_$number.root
            fi

	    if [ -f $OUTPUTDIR/GammaConvV1_LHC16defg_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16hijk_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16defghijk_$number.root $OUTPUTDIR/GammaConvV1_LHC16defg_$number.root $OUTPUTDIR/GammaConvV1_LHC16hijk_$number.root
            fi

	    if [ -f $OUTPUTDIR/GammaConvV1_LHC16defghijk_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC16lop_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC16defghijklop_$number.root $OUTPUTDIR/GammaConvV1_LHC16defghijk_$number.root $OUTPUTDIR/GammaConvV1_LHC16lop_$number.root
            fi
        done
    fi
###################---------MonteCarlos------------############

    if [ $MERGEONMC == 1 ]; then
#	rm $OUTPUTDIR/GammaConvV1_MC_LHC17f6_LHC17d1_d17_f5_d3_e5_d20a1_d20a1ex_d20a2_d20a2ex_d16_d18_*.root

        ls $OUTPUTDIR/GammaConvV1_MC_LHC17f6_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d1_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f6_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d1_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17e5_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17d3_e5_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d3_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17e5_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_extra_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_a1ex_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_extra_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_extra_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_a2ex_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_extra_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d16_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d18_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17d16_d18_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d16_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d18_$number.root
            fi
##-------------------------------------------#
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d3_e5_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_d3_e5_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d3_e5_$number.root
            fi
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_a1ex_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_a2ex_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_a1ex_a2_a2ex_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_a1ex_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d20a2_a2ex_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_d3_e5_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_a1ex_a2_a2ex_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_d3_e5_a1_a1ex_a2_a2ex_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_d3_e5_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d20a1_a1ex_a2_a2ex_$number.root
            fi
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_d3_e5_a1_a1ex_a2_a2ex_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17d16_d18_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_d3_e5_a1_a1ex_a2_a2ex_d16_d18_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f6_d1_d3_e5_a1_a1ex_a2_a2ex_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17d16_d18_$number.root
            fi


       done


    fi


   if [ $MERGEONMCJJ == 1 ]; then
#	rm $OUTPUTDIR/GammaConvV1_MC_LHC17f6_LHC17d1_d17_f5_d3_e5_d20a1_d20a1ex_d20a2_d20a2ex_d16_d18_*.root

        ls $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_f8cJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_f8aJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_f8cJJ_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_f8aJJ_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_f8cJJ_f8dJJ_f8aJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_f8cJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_f8aJJ_$number.root
            fi

            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_f8cJJ_f8dJJ_f8aJJ_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_f8cJJ_f8dJJ_f8aJJ_f8eJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_f8cJJ_f8dJJ_f8aJJ_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ_$number.root
            fi
	done
   fi


    if [ $MERGEONBINS == 1 ]; then
        echo "entered single bin merging LHC17f8bJJ"
        ls $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ$bin\_$number.root
		if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJ$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8bJJFinerPtHardBins_$number.root $TOMERGE
        done;


        echo "entered single bin merging LHC17f8c"
        ls $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ$bin\_$number.root
		if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJ$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8cJJFinerPtHardBins_$number.root $TOMERGE
        done;



        echo "entered single bin merging LHC17f8dJJ"
        ls $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ$bin\_$number.root
		if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJ$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8dJJFinerPtHardBins_$number.root $TOMERGE
        done;




        echo "entered single bin merging LHC17f8aJJ"
        ls $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ$bin\_$number.root
		if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJ$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8aJJFinerPtHardBins_$number.root $TOMERGE
        done;



        echo "entered single bin merging LHC17f8eJJ"
        ls $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ$bin\_$number.root
		if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJ$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17f8eJJFinerPtHardBins_$number.root $TOMERGE
        done;
    fi

fi
