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

# GammaConvNeutralMeson Conv Mode
#TRAINDIR=Legotrain-vAN-20150529-NeutralMeson
#N=0;
#LHC12aData=600_20150601-0953;
#LHC12bData=601_20150601-0955;
#LHC12cData=602_20150601-0956;
#LHC12dData=603_20150601-0957;
# LHC12eData=;
#LHC12fData=605_20150601-1002;
#LHC12gData=606_20150601-1003;
#LHC12hData=607_20150601-1012;
#LHC12iData=608_20150601-1014;
#LHC14e2aMC=582_20150601-1540;
#LHC14e2bMC=583_20150601-1540;
#LHC14e2cMC=584_20150601-1540;

# GammaConvNeutralMeson Calo+Mixed Mode
TRAINDIR=Legotrain-vAN-20150529-NeutralMeson
#N=1;
N=2;
LHC12aData=613_20150601-1149;
LHC12bData=614_20150601-1230;
LHC12cData=615_20150601-1232;
LHC12dData=616_20150601-1236;
#LHC12eData=;
LHC12fData=617_20150601-1237;
LHC12gData=618_20150601-1240;
LHC12hData=619_20150601-1241;
LHC12iData=620_20150601-1242;
LHC14e2aMC=579_20150601-1309;
LHC14e2bMC=580_20150601-1309;
LHC14e2cMC=581_20150601-1310;


OUTPUTDIR=$BASEDIR/$TRAINDIR
OUTPUTDIR_LHC12a=$BASEDIR/$TRAINDIR/GA_pp-$LHC12aData
OUTPUTDIR_LHC12b=$BASEDIR/$TRAINDIR/GA_pp-$LHC12bData
OUTPUTDIR_LHC12c=$BASEDIR/$TRAINDIR/GA_pp-$LHC12cData
OUTPUTDIR_LHC12d=$BASEDIR/$TRAINDIR/GA_pp-$LHC12dData
#OUTPUTDIR_LHC12e=$BASEDIR/$TRAINDIR/GA_pp-$LHC12eData
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

if [ $2 = "DL" ]; then
mkdir -p $OUTPUTDIR_LHC12a
mkdir -p $OUTPUTDIR_LHC12b
mkdir -p $OUTPUTDIR_LHC12c
mkdir -p $OUTPUTDIR_LHC12d
#mkdir -p $OUTPUTDIR_LHC12e
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
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12aData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12a/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12bData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12cData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12c/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12dData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12d/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12eData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12e/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12fData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12f/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12gData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12g/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12hData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12h/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12iData/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC12i/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_14/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2a/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_6/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2aa/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_7/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ab/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_8/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ac/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_9/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ad/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_10/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2af/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_11/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ag/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_12/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ah/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_16/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ai/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_14/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_6/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ba/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_7/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2bb/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_8/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2bc/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_9/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2bd/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_10/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2bf/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_11/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2bg/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_12/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2bh/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_16/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2bi/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_14/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2c/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_6/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ca/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_7/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2cb/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_8/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2cc/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_9/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2cd/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_10/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2cf/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_11/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2cg/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_12/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ch/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_16/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14e2ci/
return
fi
if [ $2 = "PROCESS" ]; then
	# Merge data files
	rm $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_*.root
	ls $OUTPUTDIR_LHC12a/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR_LHC12a/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] && [ -f $OUTPUTDIR_LHC12b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12ab-pass1_$number.root $OUTPUTDIR_LHC12a/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root $OUTPUTDIR_LHC12b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12ab-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12ab-pass1_$number.root ] && [ -f $OUTPUTDIR_LHC12c/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abc-pass1_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12ab-pass1_$number.root $OUTPUTDIR_LHC12c/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abc-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abc-pass1_$number.root ] && [ -f $OUTPUTDIR_LHC12d/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcd-pass1_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abc-pass1_$number.root $OUTPUTDIR_LHC12d/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done

	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcd-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcd-pass1_$number.root ] && [ -f $OUTPUTDIR_LHC12f/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdf-pass1_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcd-pass1_$number.root $OUTPUTDIR_LHC12f/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done

	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdf-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdf-pass1_$number.root ] && [ -f $OUTPUTDIR_LHC12g/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfg-pass1_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdf-pass1_$number.root $OUTPUTDIR_LHC12g/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done

	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfg-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfg-pass1_$number.root ] && [ -f $OUTPUTDIR_LHC12h/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfgh-pass1_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfg-pass1_$number.root $OUTPUTDIR_LHC12h/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done
		
	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfgh-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfgh-pass1_$number.root ] && [ -f $OUTPUTDIR_LHC12i/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfghi-pass1_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC12abcdfgh-pass1_$number.root $OUTPUTDIR_LHC12i/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done
# 	fi
	
#	# merge MC files
#	rm $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_LHC14e2b_*.root
#	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_*.root > filesForMerging.txt
#	filesForMerging=`cat filesForMerging.txt`
#	for fileName in $filesForMerging; do
#		echo $fileName
#		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
#		echo $number
#		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_$number.root ] && [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2b_$number.root ] ; then
#			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_LHC14e2b_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2b_$number.root
#		fi
#	done
#
#	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_LHC14e2b_*.root > filesForMerging.txt
#	filesForMerging=`cat filesForMerging.txt`
#	for fileName in $filesForMerging; do
#		echo $fileName
#		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
#		echo $number
#		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_LHC14e2b_$number.root ] && [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2c_$number.root ] ; then
#			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_LHC14e2b_LHC14e2c_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2a_LHC14e2b_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_MC_LHC14e2c_$number.root
#		fi
#	done
	
return
fi
