#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pp8TeV needs
NSlashes=9
NSlashes2=8
if [ $1 = "fbock" ]; then 
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp8TeV
    NSlashes=8
    NSlashes2=7
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
elif [ $1 = "nschmidt" ]; then        
    BASEDIR=/media/nschmidt/Daten/GridOutput/pp
fi


# 
# TRAINDIR=Legotrain-vAN-20150412
# # LHC12aData=498_20150414-1631-LHC12a;
# # LHC12bData=501_20150414-1634-LHC12b;
# # LHC12cData=504_20150414-1641-LHC12c; 
# # LHC12dData=508_20150414-1643-LHC12d; 
# # LHC12eData=510_20150414-1644-LHC12e; 
# # LHC12fData=513_20150414-1646-LHC12f; 
# # LHC12gData=516_20150414-1647-LHC12g; 
# # LHC12hData=519_20150414-1649-LHC12h; 
# # LHC12iData=522_20150414-1650-LHC12i; 
# 
# LHC14e2aMC=470_20150414-2134; 
# LHC14e2bMC=476_20150414-2137; 
# LHC14e2cMC=473_20150414-2136; 

TRAINDIR=Legotrain-vAN-20150825-LHC12_NonLinearity
#LHC12aData=813_20150728-0757;
LHC12bData=908_20150826-2136;
LHC12cData=909_20150826-2136;
#LHC12dData=816_20150728-0759;
#LHC12eData=;
#LHC12fData=817_20150728-0800;
#LHC12gData=818_20150728-0801;
#LHC12hData=819_20150728-0802;
#LHC12iData=820_20150728-0802;
LHC14e2aMC=916_20150824-1151;
LHC14e2bMC=917_20150824-1152;
LHC14e2cMC=918_20150824-1152;

OUTPUTDIR=$BASEDIR/$TRAINDIR
 #OUTPUTDIR_LHC12a=$BASEDIR/$TRAINDIR/GA_pp-$LHC12aData
 OUTPUTDIR_LHC12b=$BASEDIR/$TRAINDIR/GA_pp-$LHC12bData
 OUTPUTDIR_LHC12c=$BASEDIR/$TRAINDIR/GA_pp-$LHC12cData
# OUTPUTDIR_LHC12d=$BASEDIR/$TRAINDIR/GA_pp-$LHC12dData
 #OUTPUTDIR_LHC12e=$BASEDIR/$TRAINDIR/GA_pp-$LHC12eData
 #OUTPUTDIR_LHC12f=$BASEDIR/$TRAINDIR/GA_pp-$LHC12fData
 #OUTPUTDIR_LHC12g=$BASEDIR/$TRAINDIR/GA_pp-$LHC12gData
 #OUTPUTDIR_LHC12h=$BASEDIR/$TRAINDIR/GA_pp-$LHC12hData
 #OUTPUTDIR_LHC12i=$BASEDIR/$TRAINDIR/GA_pp-$LHC12iData
OUTPUTDIR_LHC14e2a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2aMC
OUTPUTDIR_LHC14e2b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2bMC
OUTPUTDIR_LHC14e2c=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC14e2cMC
# mkdir -p $OUTPUTDIR_LHC12a
# mkdir -p $OUTPUTDIR_LHC12b
mkdir -p $OUTPUTDIR_LHC12c
mkdir -p $OUTPUTDIR_LHC12d
# mkdir -p $OUTPUTDIR_LHC12e
# mkdir -p $OUTPUTDIR_LHC12f
# mkdir -p $OUTPUTDIR_LHC12g
# mkdir -p $OUTPUTDIR_LHC12h
# mkdir -p $OUTPUTDIR_LHC12i
mkdir -p $OUTPUTDIR_LHC14e2a
mkdir -p $OUTPUTDIR_LHC14e2b
mkdir -p $OUTPUTDIR_LHC14e2c

# get files from grid
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12aData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12a/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12bData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12cData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12c/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12dData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12d/
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12eData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12e/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12fData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12f/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12gData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12g/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12hData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12h/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12iData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC12i/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2aMC/merge_runlist_14/GammaConvV1* file:$OUTPUTDIR_LHC14e2a/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2bMC/merge_runlist_14/GammaConvV1* file:$OUTPUTDIR_LHC14e2b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC14e2cMC/merge_runlist_14/GammaConvV1* file:$OUTPUTDIR_LHC14e2c/


	# Change structure of data files
	ls $OUTPUTDIR_LHC12a/GammaConvV1_*.root > fileLHC12a.txt
	fileNumbers=`cat fileLHC12a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12a-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12a-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12a_$number.log\"\)

	done;

	ls $OUTPUTDIR_LHC12b/GammaConvV1_*.root > fileLHC12b.txt
	fileNumbers=`cat fileLHC12b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12b-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12b-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12b_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC12c/GammaConvV1_*.root > fileLHC12c.txt
	fileNumbers=`cat fileLHC12c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12c/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12c-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12c-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12c_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC12d/GammaConvV1_*.root > fileLHC12d.txt
	fileNumbers=`cat fileLHC12d.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12d/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12d-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12d-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12d_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12e/GammaConvV1_*.root > fileLHC12e.txt
	fileNumbers=`cat fileLHC12e.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12e/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12e-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12e-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12e_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC12f/GammaConvV1_*.root > fileLHC12f.txt
	fileNumbers=`cat fileLHC12f.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12f/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12f-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12f-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12f_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12g/GammaConvV1_*.root > fileLHC12g.txt
	fileNumbers=`cat fileLHC12g.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12g/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12g-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12g-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12g_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12h/GammaConvV1_*.root > fileLHC12h.txt
	fileNumbers=`cat fileLHC12h.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12h/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12h-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12h-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12h_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC12i/GammaConvV1_*.root > fileLHC12i.txt
	fileNumbers=`cat fileLHC12i.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12i/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC12i-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC12i-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC12i_$number.log\"\)
	done;

	# Change structure of MC files
	ls $OUTPUTDIR_LHC14e2a/GammaConvV1_*.root > fileLHC14e2a.txt
	fileNumbers=`cat fileLHC14e2a.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14e2a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC14e2a_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC14e2a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2a_$number.log\"\)
	done;
	
	ls $OUTPUTDIR_LHC14e2b/GammaConvV1_*.root > fileLHC14e2b.txt
	fileNumbers=`cat fileLHC14e2b.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14e2b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC14e2b_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC14e2b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2b_$number.log\"\)
	done;

	ls $OUTPUTDIR_LHC14e2c/GammaConvV1_*.root > fileLHC14e2c.txt
	fileNumbers=`cat fileLHC14e2c.txt`
	for fileName in $fileNumbers; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
		echo $number
		root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14e2c/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC14e2c_$number.root\"\,\"GammaConvV1_$number\"\)
		root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC14e2c_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC14e2c_$number.log\"\)
	done;
	
	# Merge data files
	rm $OUTPUTDIR/GammaConvV1_LHC12ab*-pass1_*.root
	ls $OUTPUTDIR/GammaConvV1_LHC12a-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC12a-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12b-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC12ab-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12a-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12b-pass1_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvV1_LHC12ab-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC12ab-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12c-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC12abc-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12ab-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12c-pass1_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvV1_LHC12abc-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC12abc-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12d-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcd-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abc-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12d-pass1_$number.root
		fi
	done

# 	if [ "$LHC12eData" != "" ]; then 
# 		ls $OUTPUTDIR/GammaConvV1_LHC12abcd-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcd-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12e-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcde-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcd-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12e-pass1_$number.root
# 			fi
# 		done
# 		
# 		ls $OUTPUTDIR/GammaConvV1_LHC12abcde-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcde-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12f-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdef-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcde-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12f-pass1_$number.root
# 			fi
# 		done
# 
# 		ls $OUTPUTDIR/GammaConvV1_LHC12abcdef-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcdef-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12g-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdefg-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcdef-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12g-pass1_$number.root
# 			fi
# 		done
# 		
# 		ls $OUTPUTDIR/GammaConvV1_LHC12abcdefg-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcdefg-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12h-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdefgh-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcdefg-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12h-pass1_$number.root
# 			fi
# 		done
# 
# 		ls $OUTPUTDIR/GammaConvV1_LHC12abcdefgh-pass1_*.root > filesForMerging.txt
# 		filesForMerging=`cat filesForMerging.txt`
# 		for fileName in $filesForMerging; do
# 			echo $fileName
# 			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
# 			echo $number
# 			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcdefgh-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12i-pass1_$number.root ] ; then
# 				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdefghi-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcdefgh-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12i-pass1_$number.root
# 			fi
# 		done
# 	elif 
		ls $OUTPUTDIR/GammaConvV1_LHC12abcd-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcd-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12f-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdf-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcd-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12f-pass1_$number.root
			fi
		done
		
		ls $OUTPUTDIR/GammaConvV1_LHC12abcdf-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcdf-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12g-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdfg-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcdf-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12g-pass1_$number.root
			fi
		done

		ls $OUTPUTDIR/GammaConvV1_LHC12abcdfg-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcdfg-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12h-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdfgh-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcdfg-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12h-pass1_$number.root
			fi
		done
		
		ls $OUTPUTDIR/GammaConvV1_LHC12abcdfgh-pass1_*.root > filesForMerging.txt
		filesForMerging=`cat filesForMerging.txt`
		for fileName in $filesForMerging; do
			echo $fileName
			number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
			echo $number
			if [ -f $OUTPUTDIR/GammaConvV1_LHC12abcdfgh-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC12i-pass1_$number.root ] ; then
				hadd -f $OUTPUTDIR/GammaConvV1_LHC12abcdfghi-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12abcdfgh-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC12i-pass1_$number.root
			fi
		done	
# 	fi
	
	# merge MC files
	rm $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_LHC14e2b_*.root
	ls $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC14e2b_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_LHC14e2b_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC14e2b_$number.root
		fi
	done
	
	ls $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_LHC14e2b_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_LHC14e2b_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC14e2c_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_LHC14e2b_LHC14e2c_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC14e2a_LHC14e2b_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC14e2c_$number.root
		fi
	done
	
	
