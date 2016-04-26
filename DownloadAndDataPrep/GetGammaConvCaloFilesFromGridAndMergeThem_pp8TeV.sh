#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pp8TeV needs
NSlashes=9
NSlashes2=8
if [ $1 = "fbock" ]; then 
	BASEDIR=/home/fbock/Photon/Grid/OutputLegoTrains/pp8TeV
	NSlashes=10
	NSlashes2=9
elif [ $1 = "fbockGSI" ]; then 
	BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp8TeV
elif [ $1 = "leardini" ]; then 
	BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then 
	BASEDIR=/alidata50/alice_u/leardini/GridOutput/pp8TeV/
elif [ $1 = "leardiniGSI" ]; then 
	BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/pp8TeV
elif [ $1 = "passfeld" ]; then 
	BASEDIR=~/work/Gridoutput/pp8TeV
elif [ $1 = "passfeldMAF" ]; then 
	BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then  
	BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pp8TeV
elif [ $1 = "amarin" ]; then     
	BASEDIR=/Users/marin/
elif [ $1 = "amarinGSI" ]; then     
	BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/pp8TeV 
elif [ $1 = "amarinALICESERV1" ]; then     
	BASEDIR=/alidata50/alice_u/amarin/GridOutput/pp8TeV/   
elif [ $1 = "mwilde" ]; then        
	BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
	BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pp8TeV 
elif [ $1 = "pgonzales" ]; then     
	BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
	BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pp8TeV
elif [ $1 = "dmuhlhei" ]; then
	BASEDIR=~/data/work/Grid
fi

# EMCal
#TRAINDIR=Legotrain-vAN-20150603
#LHC12aData=638_20150604-1804;
#LHC12bData=639_20150604-1813;
#LHC12cData=640_20150604-1830;
#LHC12dData=641_20150604-1815;
#LHC12eData=;
#LHC12fData=642_20150604-1832;
#LHC12gData=643_20150604-1828;
#LHC12hData=644_20150604-1823;
#LHC12iData=645_20150604-1833;
#LHC14e2aMC=606_20150604-1752;
#LHC14e2bMC=607_20150604-2153;
#LHC14e2cMC=608_20150604-1755;

#TRAINDIR=Legotrain-vAN-20150616
#LHC12aData=680_20150616-1911;
#LHC12bData=681_20150616-1911;
#LHC12cData=682_20150616-1911;
#LHC12dData=683_20150616-1912;
#LHC12eData=;
#LHC12fData=690_20150617-1612;
#LHC12gData=686_20150616-1920;
#LHC12hData=679_20150616-1828;
#LHC12iData=678_20150616-1828;
#LHC14e2aMC=;
#LHC14e2bMC=;
#LHC14e2cMC=;

#TRAINDIR=Legotrain-vAN-20150619
#LHC12aData=691_20150620-0025;
#LHC12bData=692_20150620-0025;
#LHC12cData=693_20150620-0026;
#LHC12dData=694_20150620-0026;
#LHC12eData=;
#LHC12fData=695_20150620-0029;
#LHC12gData=696_20150620-0029;
#LHC12hData=697_20150620-0029;
#LHC12iData=698_20150620-0030;
#LHC14e2aMC=;
#LHC14e2bMC=;
#LHC14e2cMC=;

#TRAINDIR=Legotrain-vAN-20150624
#LHC12aData=700_20150624-1719;
#LHC12bData=701_20150624-1720;
#LHC12cData=702_20150624-1720;
#LHC12dData=703_20150624-1721;
#LHC12eData=;
#LHC12fData=704_20150624-1721;
#LHC12gData=705_20150624-1722;
#LHC12hData=706_20150624-1722;
#LHC12iData=707_20150624-1723;
#LHC14e2aMC=663_20150624-1411;
#LHC14e2bMC=664_20150624-1411;
#LHC14e2cMC=665_20150624-1416;

#TRAINDIR=Legotrain-vAN-20150625
#LHC12aData=708_20150625-0919;
#LHC12bData=709_20150625-0920;
#LHC12cData=710_20150625-0922;
#LHC12dData=711_20150625-0922;
#LHC12eData=712_20150625-0923;
#LHC12fData=713_20150625-0924;
#LHC12gData=714_20150625-0925;
#LHC12hData=715_20150625-0925;
#LHC12iData=716_20150625-0927;
##LHC14e2aMC=;
##LHC14e2bMC=;
##LHC14e2cMC=;

#TRAINDIR=Legotrain-vAN-20150728-LHC12_QA-2nd
#LHC12aData=813_20150728-0757;
#LHC12bData=814_20150728-0758;
#LHC12cData=815_20150728-0758;
#LHC12dData=816_20150728-0759;
#LHC12eData=;
#LHC12fData=817_20150728-0800;
#LHC12gData=818_20150728-0801;
#LHC12hData=819_20150728-0802;
#LHC12iData=820_20150728-0802;
#LHC14e2aMC=798_20150728-0851;
#LHC14e2bMC=799_20150728-0852;
#LHC14e2cMC=800_20150728-0854;

#TRAINDIR=Legotrain-vAN-20150728-LHC12_p2_QA
#LHC12aData=821_20150728-0843;
#LHC12bData=822_20150728-0845;
#LHC12cData=823_20150728-0846;
#LHC12dData=824_20150728-0846;
#LHC12eData=;
#LHC12fData=825_20150728-0847;
#LHC12gData=826_20150728-0849;
#LHC12hData=827_20150728-0850;
#LHC12iData=828_20150728-0850;

#LHC14e2aMC=;
#LHC14e2bMC=;
#LHC14e2cMC=;

#TRAINDIR=Legotrain-vAN-20150825-LHC12_NonLinearity
#LHC12aData=813_20150728-0757;
#LHC12bData=908_20150826-2136;
#LHC12cData=909_20150826-2136;
#LHC12dData=816_20150728-0759;
#LHC12eData=;
#LHC12fData=817_20150728-0800;
#LHC12gData=818_20150728-0801;
#LHC12hData=819_20150728-0802;
#LHC12iData=820_20150728-0802;
#LHC14e2aMC=916_20150824-1151;
#LHC14e2bMC=917_20150824-1152;
#LHC14e2cMC=918_20150824-1152;

#TRAINDIR=Legotrain-vAN-20150728-LHC12_pass2
#LHC12aData=821_20150728-0843;
#LHC12bData=822_20150728-0845;
#LHC12cData=823_20150728-0846;
#LHC12dData=824_20150728-0846;
#LHC12eData=;
#LHC12fData=825_20150728-0847;
#LHC12gData=826_20150728-0849;
#LHC12hData=827_20150728-0850;
#LHC12iData=828_20150728-0850;
#LHC14e2aMC=916_20150824-1151;
#LHC14e2bMC=917_20150824-1152;
#LHC14e2cMC=918_20150824-1152;

#TRAINDIR=Legotrain-vAN-20150903-LHC12_pass2
#LHC12aData=927_20150903-2031;
#LHC12bData=928_20150903-2031;
#LHC12cData=929_20150903-2031;
#LHC12dData=931_20150903-2033;
#LHC12eData=;
#LHC12fData=930_20150903-2032;
#LHC12gData=932_20150903-2033;
#LHC12hData=934_20150903-2034;
#LHC12iData=936_20150903-2037;
#LHC14e2aMC=916_20150824-1151;
#LHC14e2bMC=917_20150824-1152;
#LHC14e2cMC=918_20150824-1152;

#TRAINDIR=Legotrain-vAN-20150903-LHC12_pass2-oldXrootD
#LHC12aData=977_20150909-1023;
#LHC12bData=978_20150909-1021;
#LHC12cData=979_20150909-1022;
#LHC12dData=980_20150909-1025;
#LHC12eData=;
#LHC12fData=981_20150909-1024;
#LHC12gData=982_20150909-1029;
#LHC12hData=983_20150909-1029;
#LHC12iData=984_20150909-1026;
#LHC14e2aMC=916_20150824-1151;
#LHC14e2bMC=917_20150824-1152;
#LHC14e2cMC=918_20150824-1152;

#TRAINDIR=Legotrain-vAN-20151002-LHC12_pass2-QA
#LHC12aData=1011_20151002-1943;
#LHC12bData=1012_20151002-1943;
#LHC12cData=1013_20151002-1946;
#LHC12dData=1014_20151002-1946;
#LHC12eData=;
#LHC12fData=1015_20151002-1947;
#LHC12gData=1016_20151002-1947;
#LHC12hData=1018_20151002-1948;
#LHC12iData=1019_20151002-1948;
#LHC14e2aMC=1145_20151002-1951;
#LHC14e2bMC=1146_20151002-1952;
#LHC14e2cMC=1147_20151002-1952;

TRAINDIR=Legotrain-vAN-20151025-LHC12_pass2-TriggerQA
LHC12aData=1046_20151026-1818;
LHC12bData=1047_20151026-1819;
LHC12cData=1048_20151026-1819;
LHC12dData=1049_20151026-1819;
#LHC12eData=;
LHC12fData=1050_20151026-1819;
LHC12gData=1051_20151026-1820;
LHC12hData=1052_20151026-1821;
LHC12iData=1053_20151026-1821;
#LHC14e2aMC=1145_20151002-1951;
#LHC14e2bMC=1146_20151002-1952;
#LHC14e2cMC=1147_20151002-1952;


OUTPUTDIR=$BASEDIR/$TRAINDIR
OUTPUTDIR_LHC12a=$BASEDIR/$TRAINDIR/GA_pp-$LHC12aData
OUTPUTDIR_LHC12b=$BASEDIR/$TRAINDIR/GA_pp-$LHC12bData
OUTPUTDIR_LHC12c=$BASEDIR/$TRAINDIR/GA_pp-$LHC12cData
OUTPUTDIR_LHC12d=$BASEDIR/$TRAINDIR/GA_pp-$LHC12dData
OUTPUTDIR_LHC12e=$BASEDIR/$TRAINDIR/GA_pp-$LHC12eData
OUTPUTDIR_LHC12f=$BASEDIR/$TRAINDIR/GA_pp-$LHC12fData
OUTPUTDIR_LHC12g=$BASEDIR/$TRAINDIR/GA_pp-$LHC12gData
OUTPUTDIR_LHC12h=$BASEDIR/$TRAINDIR/GA_pp-$LHC12hData
OUTPUTDIR_LHC12i=$BASEDIR/$TRAINDIR/GA_pp-$LHC12iData
OUTPUTDIR_LHC14e2a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC
OUTPUTDIR_LHC14e2aa=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/a
OUTPUTDIR_LHC14e2ab=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/b
OUTPUTDIR_LHC14e2ac=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/c
OUTPUTDIR_LHC14e2ad=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/d
OUTPUTDIR_LHC14e2af=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/f
OUTPUTDIR_LHC14e2ag=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/g
OUTPUTDIR_LHC14e2ah=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/h
OUTPUTDIR_LHC14e2ai=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC/i
OUTPUTDIR_LHC14e2b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC
OUTPUTDIR_LHC14e2ba=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/a
OUTPUTDIR_LHC14e2bb=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/b
OUTPUTDIR_LHC14e2bc=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/c
OUTPUTDIR_LHC14e2bd=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/d
OUTPUTDIR_LHC14e2bf=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/f
OUTPUTDIR_LHC14e2bg=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/g
OUTPUTDIR_LHC14e2bh=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/h
OUTPUTDIR_LHC14e2bi=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC/i
OUTPUTDIR_LHC14e2c=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC
OUTPUTDIR_LHC14e2ca=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/a
OUTPUTDIR_LHC14e2cb=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/b
OUTPUTDIR_LHC14e2cc=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/c
OUTPUTDIR_LHC14e2cd=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/d
OUTPUTDIR_LHC14e2cf=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/f
OUTPUTDIR_LHC14e2cg=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/g
OUTPUTDIR_LHC14e2ch=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/h
OUTPUTDIR_LHC14e2ci=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC/i

mkdir -p $OUTPUTDIR_LHC12a
mkdir -p $OUTPUTDIR_LHC12b
mkdir -p $OUTPUTDIR_LHC12c
mkdir -p $OUTPUTDIR_LHC12d
mkdir -p $OUTPUTDIR_LHC12e
mkdir -p $OUTPUTDIR_LHC12f
mkdir -p $OUTPUTDIR_LHC12g
mkdir -p $OUTPUTDIR_LHC12h
mkdir -p $OUTPUTDIR_LHC12i
mkdir -p $OUTPUTDIR_LHC14e2a
mkdir -p $OUTPUTDIR_LHC14e2aa
mkdir -p $OUTPUTDIR_LHC14e2ab
mkdir -p $OUTPUTDIR_LHC14e2ac
mkdir -p $OUTPUTDIR_LHC14e2ad
mkdir -p $OUTPUTDIR_LHC14e2af
mkdir -p $OUTPUTDIR_LHC14e2ag
mkdir -p $OUTPUTDIR_LHC14e2ah
mkdir -p $OUTPUTDIR_LHC14e2ai
mkdir -p $OUTPUTDIR_LHC14e2b
mkdir -p $OUTPUTDIR_LHC14e2ba
mkdir -p $OUTPUTDIR_LHC14e2bb
mkdir -p $OUTPUTDIR_LHC14e2bc
mkdir -p $OUTPUTDIR_LHC14e2bd
mkdir -p $OUTPUTDIR_LHC14e2bf
mkdir -p $OUTPUTDIR_LHC14e2bg
mkdir -p $OUTPUTDIR_LHC14e2bh
mkdir -p $OUTPUTDIR_LHC14e2bi
mkdir -p $OUTPUTDIR_LHC14e2c
mkdir -p $OUTPUTDIR_LHC14e2ca
mkdir -p $OUTPUTDIR_LHC14e2cb
mkdir -p $OUTPUTDIR_LHC14e2cc
mkdir -p $OUTPUTDIR_LHC14e2cd
mkdir -p $OUTPUTDIR_LHC14e2cf
mkdir -p $OUTPUTDIR_LHC14e2cg
mkdir -p $OUTPUTDIR_LHC14e2ch
mkdir -p $OUTPUTDIR_LHC14e2ci

# get files from grid
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12aData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12a/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12bData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12cData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12c/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12dData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12d/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12eData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12e/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12fData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12f/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12gData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12g/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12hData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12h/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12iData/merge/GammaConvCalo* file:$OUTPUTDIR_LHC12i/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_14/GammaConvCalo* file:$OUTPUTDIR_LHC14e2a/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_6/GammaConvCalo* file:$OUTPUTDIR_LHC14e2aa/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_7/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ab/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_8/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ac/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_9/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ad/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_10/GammaConvCalo* file:$OUTPUTDIR_LHC14e2af/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_11/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ag/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_12/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ah/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_16/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ai/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_14/GammaConvCalo* file:$OUTPUTDIR_LHC14e2b/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_6/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ba/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_7/GammaConvCalo* file:$OUTPUTDIR_LHC14e2bb/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_8/GammaConvCalo* file:$OUTPUTDIR_LHC14e2bc/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_9/GammaConvCalo* file:$OUTPUTDIR_LHC14e2bd/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_10/GammaConvCalo* file:$OUTPUTDIR_LHC14e2bf/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_11/GammaConvCalo* file:$OUTPUTDIR_LHC14e2bg/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_12/GammaConvCalo* file:$OUTPUTDIR_LHC14e2bh/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_16/GammaConvCalo* file:$OUTPUTDIR_LHC14e2bi/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_14/GammaConvCalo* file:$OUTPUTDIR_LHC14e2c/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_6/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ca/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_7/GammaConvCalo* file:$OUTPUTDIR_LHC14e2cb/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_8/GammaConvCalo* file:$OUTPUTDIR_LHC14e2cc/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_9/GammaConvCalo* file:$OUTPUTDIR_LHC14e2cd/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_10/GammaConvCalo* file:$OUTPUTDIR_LHC14e2cf/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_11/GammaConvCalo* file:$OUTPUTDIR_LHC14e2cg/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_12/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ch/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_16/GammaConvCalo* file:$OUTPUTDIR_LHC14e2ci/

	# Change structure of data files
	ls $OUTPUTDIR_LHC12a/GammaConvCalo_*.root > fileLHC12a.txt
	fileNumbers=`cat fileLHC12a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12a/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12a-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12a-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12a_$number.log\"\)

	done;

	ls $OUTPUTDIR_LHC12b/GammaConvCalo_*.root > fileLHC12b.txt
	fileNumbers=`cat fileLHC12b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12b/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12b-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12b-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12b_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC12c/GammaConvCalo_*.root > fileLHC12c.txt
	fileNumbers=`cat fileLHC12c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12c/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12c-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12c-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12c_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC12d/GammaConvCalo_*.root > fileLHC12d.txt
	fileNumbers=`cat fileLHC12d.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12d/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12d-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12d-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12d_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12e/GammaConvCalo_*.root > fileLHC12e.txt
	fileNumbers=`cat fileLHC12e.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12e/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12e-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12e-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12e_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC12f/GammaConvCalo_*.root > fileLHC12f.txt
	fileNumbers=`cat fileLHC12f.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12f/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12f-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12f-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12f_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12g/GammaConvCalo_*.root > fileLHC12g.txt
	fileNumbers=`cat fileLHC12g.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12g/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12g-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12g-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12g_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12h/GammaConvCalo_*.root > fileLHC12h.txt
	fileNumbers=`cat fileLHC12h.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12h/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12h-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12h-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12h_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12i/GammaConvCalo_*.root > fileLHC12i.txt
	fileNumbers=`cat fileLHC12i.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC12i/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_LHC12i-pass1_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC12i-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12i_$number.log\"\)
	done;

	# Change structure of MC files
	ls $OUTPUTDIR_LHC14e2a/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2a/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2aa/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2aa/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12a_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12a_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ab/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ab/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12b_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12b_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ac/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ac/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12c_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12c_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12c_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ad/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ad/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12d_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12d_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12d_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2af/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2af/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12f_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12f_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12f_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ag/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ag/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12g_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12g_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12g_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ah/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ah/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12h_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12h_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12h_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ai/GammaConvCalo_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ai/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12i_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2a-LHC12i_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a-LHC12i_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC14e2b/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2b/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ba/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ba/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12a_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12a_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2bb/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2bb/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12b_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12b_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2bc/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2bc/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12c_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12c_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12c_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2bd/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2bd/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12d_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12d_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12d_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2bf/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2bf/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12f_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12f_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12f_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2bg/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2bg/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12g_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12g_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12g_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2bh/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2bh/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12h_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12h_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12h_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2bi/GammaConvCalo_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2bi/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12i_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2b-LHC12i_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b-LHC12i_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2c/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2c/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ca/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ca/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12a_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12a_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2cb/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2cb/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12b_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12b_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2cc/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2cc/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12c_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12c_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12c_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2cd/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2cd/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12d_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12d_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12d_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2cf/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2cf/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12f_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12f_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12f_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2cg/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2cg/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12g_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12g_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12g_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ch/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ch/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12h_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12h_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12h_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2ci/GammaConvCalo_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $((NSlashes+1)) | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$OUTPUTDIR_LHC14e2ci/GammaConvCalo_$number.root\"\,\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12i_$number.root\"\,\"GammaConvCalo_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC14e2c-LHC12i_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c-LHC12i_$number.log\"\)
	done;
	
	# Merge data files
	rm $OUTPUTDIR/GammaConvCalo_LHC12ab*-pass1_*.root
	ls $OUTPUTDIR/GammaConvCalo_LHC12a-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvCalo_LHC12a-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12b-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvCalo_LHC12ab-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12a-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12b-pass1_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvCalo_LHC12ab-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvCalo_LHC12ab-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12c-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abc-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12ab-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12c-pass1_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvCalo_LHC12abc-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abc-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12d-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcd-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abc-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12d-pass1_$number.root
		fi
	done

# 	if [ "$LHC12eData" != "" ]; then 
# 		ls $OUTPUTDIR/GammaConvCalo_LHC12abcd-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcd-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12e-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcde-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcd-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12e-pass1_$number.root
# 			fi
# 		done
# 		
# 		ls $OUTPUTDIR/GammaConvCalo_LHC12abcde-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcde-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12f-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdef-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcde-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12f-pass1_$number.root
# 			fi
# 		done
# 
# 		ls $OUTPUTDIR/GammaConvCalo_LHC12abcdef-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcdef-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12g-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdefg-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcdef-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12g-pass1_$number.root
# 			fi
# 		done
# 		
# 		ls $OUTPUTDIR/GammaConvCalo_LHC12abcdefg-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcdefg-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12h-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdefgh-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcdefg-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12h-pass1_$number.root
# 			fi
# 		done
# 
# 		ls $OUTPUTDIR/GammaConvCalo_LHC12abcdefgh-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcdefgh-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12i-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdefghi-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcdefgh-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12i-pass1_$number.root
# 			fi
# 		done
# 	elif 
		ls $OUTPUTDIR/GammaConvCalo_LHC12abcd-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcd-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12f-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdf-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcd-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12f-pass1_$number.root
			fi
		done
		
		ls $OUTPUTDIR/GammaConvCalo_LHC12abcdf-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcdf-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12g-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdfg-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcdf-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12g-pass1_$number.root
			fi
		done

		ls $OUTPUTDIR/GammaConvCalo_LHC12abcdfg-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcdfg-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12h-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdfgh-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcdfg-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12h-pass1_$number.root
			fi
		done
		
		ls $OUTPUTDIR/GammaConvCalo_LHC12abcdfgh-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvCalo_LHC12abcdfgh-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC12i-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvCalo_LHC12abcdfghi-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12abcdfgh-pass1_$number.root $OUTPUTDIR/GammaConvCalo_LHC12i-pass1_$number.root
			fi
		done	
# 	fi
	
	# merge MC files
	rm $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_LHC14e2b_*.root
	ls $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC14e2b_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_LHC14e2b_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC14e2b_$number.root
		fi
	done

	ls $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_LHC14e2b_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_LHC14e2b_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC14e2c_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_LHC14e2b_LHC14e2c_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC14e2a_LHC14e2b_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC14e2c_$number.root
		fi
	done
	
	
