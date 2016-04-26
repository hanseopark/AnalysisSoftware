#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pp8TeV needs
NSlashes=9
NSlashes2=8
if [ $1 = "fbock" ]; then 
	BASEDIR=/home/fbock/Photon/Grid/OutputLegoTrains/pp7TeV
	NSlashes=10
	NSlashes2=9
elif [ $1 = "dmuhlhei" ]; then
	BASEDIR=~/data/work/Grid
fi

#GammaConvNeutralMeson Conv Mode
#TRAINDIR=Legotrain-vAN-20151102-NeutralMeson7TeV
#N=0;
#LHC10bData=1054_20151027-1449;
#LHC10cData=1055_20151027-1451;
#LHC10dData=1056_20151027-1452;
#LHC10eData=1057_20151027-1452;
#LHC10fData=1058_20151027-1453;
#LHC14j4b=1265_20151027-1529;
#LHC14j4c=1266_20151027-1530;
#LHC14j4d=1267_20151027-1531;
#LHC14j4e=1268_20151027-1531;
#LHC14j4f=1269_20151027-1531;

TRAINDIR=Legotrain-vAN-20151123-NeutralMeson7TeV_Calo
N=2;
LHC10bData=1102_20151123-1746;
LHC10cData=1103_20151123-1745;
LHC10dData=1104_20151123-1745;
LHC10eData=1105_20151123-1745;
LHC10fData=1106_20151123-1745;
LHC14j4b=1326_20151123-1744;
LHC14j4c=1327_20151123-1745;
LHC14j4d=1328_20151123-1745;
LHC14j4e=1329_20151123-1745;
LHC14j4f=1330_20151123-1745;


OUTPUTDIR=$BASEDIR/$TRAINDIR
OUTPUTDIR_LHC10b=$BASEDIR/$TRAINDIR/GA_pp-$LHC10bData
OUTPUTDIR_LHC10c=$BASEDIR/$TRAINDIR/GA_pp-$LHC10cData
OUTPUTDIR_LHC10d=$BASEDIR/$TRAINDIR/GA_pp-$LHC10dData
OUTPUTDIR_LHC10e=$BASEDIR/$TRAINDIR/GA_pp-$LHC10eData
OUTPUTDIR_LHC10f=$BASEDIR/$TRAINDIR/GA_pp-$LHC10fData
OUTPUTDIR_LHC14j4b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14j4b
OUTPUTDIR_LHC14j4c=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14j4c
OUTPUTDIR_LHC14j4d=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14j4d
OUTPUTDIR_LHC14j4e=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14j4e
OUTPUTDIR_LHC14j4f=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14j4f

mkdir -p $OUTPUTDIR_LHC10b
mkdir -p $OUTPUTDIR_LHC10c
mkdir -p $OUTPUTDIR_LHC10d
mkdir -p $OUTPUTDIR_LHC10e
mkdir -p $OUTPUTDIR_LHC10f
mkdir -p $OUTPUTDIR_LHC14j4b
mkdir -p $OUTPUTDIR_LHC14j4c
mkdir -p $OUTPUTDIR_LHC14j4d
mkdir -p $OUTPUTDIR_LHC14j4e
mkdir -p $OUTPUTDIR_LHC14j4f

# get files from grid
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10bData/merge_runlist_1/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC10b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10cData/merge_runlist_1/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC10c/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10dData/merge_runlist_1/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC10d/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10eData/merge_runlist_1/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC10e/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC10fData/merge_runlist_1/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC10f/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4b/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14j4b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4c/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14j4c/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4d/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14j4d/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4e/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14j4e/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14j4f/merge/GammaConvNeutralMesonPiPlPiMiPiZero* file:$OUTPUTDIR_LHC14j4f/
	
	# Merge data files
	rm $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_*.root
	ls $OUTPUTDIR_LHC10b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR_LHC10b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] && [ -f $OUTPUTDIR_LHC10c/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bc-pass4_$number.root $OUTPUTDIR_LHC10b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root $OUTPUTDIR_LHC10c/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bc-pass4_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bc-pass4_$number.root ] && [ -f $OUTPUTDIR_LHC10d/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcd-pass4_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bc-pass4_$number.root $OUTPUTDIR_LHC10d/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcd-pass4_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcd-pass4_$number.root ] && [ -f $OUTPUTDIR_LHC10e/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcde-pass4_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcd-pass4_$number.root $OUTPUTDIR_LHC10e/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done

	ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcde-pass4_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcde-pass4_$number.root ] && [ -f $OUTPUTDIR_LHC10f/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcdef-pass4_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_LHC10bcde-pass4_$number.root $OUTPUTDIR_LHC10f/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
		fi
	done
		
#	# merge MC files
        rm $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_*.root
        ls $OUTPUTDIR_LHC14j4b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
                echo $number
                if [ -f $OUTPUTDIR_LHC14j4b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] && [ -f $OUTPUTDIR_LHC14j4c/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
                        hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bc-pass4_$number.root $OUTPUTDIR_LHC14j4b/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root $OUTPUTDIR_LHC14j4c/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
                fi
        done

        ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bc-pass4_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
                echo $number
                if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bc-pass4_$number.root ] && [ -f $OUTPUTDIR_LHC14j4d/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
                        hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcd-pass4_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bc-pass4_$number.root $OUTPUTDIR_LHC14j4d/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
                fi
        done

        ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcd-pass4_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
                echo $number
                if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcd-pass4_$number.root ] && [ -f $OUTPUTDIR_LHC14j4e/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
                        hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcde-pass4_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcd-pass4_$number.root $OUTPUTDIR_LHC14j4e/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
                fi
        done

        ls $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcde-pass4_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
                echo $number
                if [ -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcde-pass4_$number.root ] && [ -f $OUTPUTDIR_LHC14j4f/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root ] ; then
                        hadd -f $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcdef-pass4_$number.root $OUTPUTDIR/GammaConvNeutralMesonPiPlPiMiPiZero_MC_${N}_LHC14j4bcde-pass4_$number.root $OUTPUTDIR_LHC14j4f/GammaConvNeutralMesonPiPlPiMiPiZero_${N}_$number.root
                fi
        done
