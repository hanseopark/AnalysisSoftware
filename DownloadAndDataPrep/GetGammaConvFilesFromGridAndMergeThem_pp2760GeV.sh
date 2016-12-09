#! /bin/bash

# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=0
MERGEONBINSSingle=0
MERGEONBINS=0

# check if train configuration has actually been given
HAVELHC11a=1
HAVELHC12f1a=1
HAVELHC12f1b=1
HAVELHC12i3=1
HAVELHC15g1a=1
HAVELHC13g=1
HAVELHC15g2=1
HAVELHC15a3a=1
HAVELHC15a3aplus=1

# default trainconfigurations
LHC11aData="";
LHC12f1aMC="";
LHC12f1bMC="";
LHC12i3MC="";
LHC15g1aMC="";
LHC13gData="";
LHC15a3aMC=""; 
LHC15a3aplusMC=""; 
LHC15g2MC=""; 


# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pp2760GeV needs
NSlashes=9
NSlashes2=8
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
	BASEDIR=/Users/marin/
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
# LHC11aData="";
# LHC12f1aMC="1256_20151026-1603";
# LHC12f1bMC="1257_20151026-1604";
# LHC12i3MC="";
# LHC15g1aMC="1262_20151026-1610";
# 
# LHC13gData="";
# LHC15a3aMC="1260_20151026-1609"; et
# LHC15a3aplusMC="1261_20151026-1609"; 
# LHC15g2MC="1263_20151026-1626"; 
# LHC15g2MC="1263_20151026-1626"; 

# TRAINDIR=Legotrain-vAN-20151103-MultiplicityDependence2
# # LHC11aData="371_20141220-0947";
# LHC12f1aMC="1276_20151103-1737"; 
# LHC12f1bMC="1277_20151103-1737"; 
# LHC15g1aMC="1275_20151103-1735"; 

# TRAINDIR=Legotrain-vAN-20151106-MultiplicityDependence3
# # LHC11aData="371_20141220-0947";
# LHC12f1aMC="1280_20151106-2332"; 
# LHC12f1bMC="1281_20151106-2332"; 
# # LHC15g1aMC="1275_20151103-1735"; 

# TRAINDIR=Legotrain-vAN-20151106-MultiplicityDependence3_woWeight
# # LHC11aData="371_20141220-0947";
# LHC12f1aMC="1282_20151106-2336"; 
# LHC12f1bMC="1283_20151106-2336"; 
# # LHC15g1aMC="1275_20151103-1735"; 

# TRAINDIR=Legotrain-vAN-20151218-TriggersWithoutPileupAndPCMTriggers
# # LHC13gData="1136_20151219-1548"
# # LHC15g2MC="1373_20151218-2155";
# # LHC15a3aMC="1377_20151219-2110"; 
# # LHC15a3aplusMC="1375_20151219-2110"; 
# LHC13gData="1137_20151219-1549"
# LHC15g2MC="1374_20151218-2157";
# LHC15a3aMC="1378_20151219-1302"; 
# LHC15a3aplusMC="1376_20151219-1300"; 
# LHC11aData="1135_20151219-1547";
# LHC12f1aMC="1371_20151218-2151"; 
# LHC12f1bMC="1372_20151218-2153"; 
# LHC15g1aMC="1379_20151219-2111"; 

TRAINDIR=Legotrain-vAN20160826_ConvReweigting
LHC11aData="1782_20160813-2128";
# LHC12f1aMC="2467_20160826-1001"; 
# LHC12f1bMC="2468_20160826-1001"; 
# LHC12i3MC="2570_20161024-1234"

DATAADD="WOSDD"

OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ $2 == "LHC11a" ]; then 

    if [ "$LHC11aData" == "" ]; then 
        HAVELHC11a=0;
    fi
    if [ "$LHC12f1aMC" = "" ]; then 
        HAVELHC12f1a=0; 
    fi
    if [ "$LHC12f1bMC" = "" ]; then 
        HAVELHC12f1b=0; 
    fi
    if [ "$LHC12i3MC" = "" ]; then 
        HAVELHC12i3=0; 
    fi
    if [ "$LHC15g1aMC" = "" ]; then 
        HAVELHC15g1a=0; 
    fi
    
    if [ $HAVELHC11a == 1 ]; then
        OUTPUTDIR_LHC11a=$BASEDIR/$TRAINDIR/GA_pp-$LHC11aData
        mkdir -p $OUTPUTDIR_LHC11a
    fi
    if [ $HAVELHC12f1a == 1 ]; then
        OUTPUTDIR_LHC12f1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1aMC
        mkdir -p $OUTPUTDIR_LHC12f1a
    fi  
    
    if [ $HAVELHC12f1b == 1 ]; then
        OUTPUTDIR_LHC12f1b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1bMC
        mkdir -p $OUTPUTDIR_LHC12f1b
    fi
    
    if [ $HAVELHC12i3 == 1 ]; then
       OUTPUTDIR_LHC12i3=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12i3MC
       mkdir -p $OUTPUTDIR_LHC12i3
    fi
    if [ $HAVELHC15g1a == 1 ]; then
       OUTPUTDIR_LHC15g1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g1aMC
       mkdir -p $OUTPUTDIR_LHC15g1a
    fi
       
    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC11a == 1 ]; then
            echo "downloading LHC11a"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC11aData/merge_runlist_1/root_archive.zip file:$OUTPUTDIR_LHC11a/
            unzip -u $OUTPUTDIR_LHC11a/root_archive.zip -d $OUTPUTDIR_LHC11a/
        fi  
        if [ $HAVELHC12f1a == 1 ]; then
            echo "downloading LHC12f1a"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1aMC/merge_runlist_1/root_archive.zip file:$OUTPUTDIR_LHC12f1a/
            unzip -u $OUTPUTDIR_LHC12f1a/root_archive.zip -d $OUTPUTDIR_LHC12f1a/
        fi    
        if [ $HAVELHC12f1b == 1 ]; then
            echo "downloading LHC12f1b"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1bMC/merge_runlist_1/root_archive.zip file:$OUTPUTDIR_LHC12f1b/
            unzip -u $OUTPUTDIR_LHC12f1b/root_archive.zip -d $OUTPUTDIR_LHC12f1b/
        fi
        if [ $HAVELHC12i3 == 1 ]; then
            echo "downloading LHC12i3"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12i3MC/merge_runlist_1/root_archive.zip file:$OUTPUTDIR_LHC12i3/
            unzip -u $OUTPUTDIR_LHC12i3/root_archive.zip -d $OUTPUTDIR_LHC12i3/
        fi
        if [ $HAVELHC15g1a == 1 ]; then
            echo "downloading LHC15g1a"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g1aMC/merge/root_archive.zip file:$OUTPUTDIR_LHC15g1a/
            unzip -u $OUTPUTDIR_LHC15g1a/root_archive.zip -d $OUTPUTDIR_LHC15g1a/
            
            echo "copying LHC15g1a in single bins" 
            runNumbers=`cat runNumbersLHC11aJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat binNumbersJJ.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    mkdir -p $OUTPUTDIR_LHC15g1a/$binNumber/$runNumber
                    alien_cp alien:/alice/sim/2015/LHC15g1a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15g1aMC/root_archive.zip file:$OUTPUTDIR_LHC15g1a/$binNumber/$runNumber/
                    unzip -u $OUTPUTDIR_LHC15g1a/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC15g1a/$binNumber/$runNumber/
                done;   
            done;
        fi
    fi 

    if [ $HAVELHC11a == 1 ]; then
        ls $OUTPUTDIR_LHC11a/GammaConvV1_*.root > fileLHC11a.txt
        fileNumbers=`cat fileLHC11a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            cp $OUTPUTDIR_LHC11a/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_LHC11a-pass4-$DATAADD\_$number.root
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC11a-pass4-$DATAADD\_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC11a_$number.log\"\,0\)
        done;
    fi
    
    if [ $HAVELHC12f1a == 1 ]; then
        ls $OUTPUTDIR_LHC12f1a/GammaConvV1_*.root > fileLHC12f1a.txt
        fileNumbers=`cat fileLHC12f1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            cp $OUTPUTDIR_LHC12f1a/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_$number.root
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12f1a_$number.log\"\,0\)
        done;
    fi

    if [ $HAVELHC12f1b == 1 ]; then
        ls $OUTPUTDIR_LHC12f1b/GammaConvV1_*.root > fileLHC12f1b.txt
        fileNumbers=`cat fileLHC12f1b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            cp $OUTPUTDIR_LHC12f1b/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12f1b-$DATAADD\_$number.root
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC12f1b-$DATAADD\_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12f1b_$number.log\"\,0\)
        done;
    fi
    
    if [ $HAVELHC12i3 == 1 ]; then
        ls $OUTPUTDIR_LHC12i3/GammaConvV1_*.root > fileLHC12i3.txt
        fileNumbers=`cat fileLHC12i3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            cp $OUTPUTDIR_LHC12i3/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12i3-$DATAADD\_$number.root
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC12i3-$DATAADD\_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC12i3_$number.log\"\,0\)
        done;
    fi

    if [ $HAVELHC15g1a == 1 ]; then
        ls $OUTPUTDIR_LHC15g1a/GammaConvV1_*.root > fileLHC15g1a.txt
        fileNumbers=`cat fileLHC15g1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            cp $OUTPUTDIR_LHC15g1a/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC15g1a_$number.root
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC15g1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15g1a_$number.log\"\,0\)
        done;
    fi
    
    if [ $MERGEONBINSSingle == 1 ]; then
        binNumbersJJ=`cat binNumbersJJ.txt`
        echo $binNumbersJJ
        ls $OUTPUTDIR_LHC15g1a/GammaConvV1_*.root > filetemp.txt
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number            
                hadd -f $OUTPUTDIR_LHC15g1a/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC15g1a/$binNumber/*/GammaConvV1_$number.root
                cp $OUTPUTDIR_LHC15g1a/$binNumber/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC15g1a$binNumber\_$number.root
            done;
        done;
    fi
    
    if [ $MERGEON == 1 ]; then
        rm $OUTPUTDIR/GammaConvV1_MC_LHC12f1a_LHC12f1b-$DATAADD\_*.root
        ls $OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC12f1b-$DATAADD\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC12f1a_LHC12f1b-$DATAADD\_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12f1b-$DATAADD\_$number.root
            fi
        done
        
        ls $OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC12i3-$DATAADD\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC12f1a_LHC12i3-$DATAADD\_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12i3-$DATAADD\_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12f1a-$DATAADD\_$number.root
            fi
        done
        
        
        ls $OUTPUTDIR/GammaConvV1_MC_LHC12f1a_LHC12f1b-$DATAADD\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC12f1a_LHC12f1b-$DATAADD\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC12i3-$DATAADD\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC12f1a_LHC12f1b_LHC12i3-$DATAADD\_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12f1a_LHC12f1b-$DATAADD\_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC12i3-$DATAADD\_$number.root
            fi
        done
    fi    
    
    if [ $MERGEONBINS == 1 ]; then    
        echo "entered single bin merging"
        ls $OUTPUTDIR/GammaConvV1_MC_LHC15g1a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/GammaConvV1_MC_LHC15g1a$bin\_$number.root
                if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC15g1a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC15g1a$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/GammaConvV1_MC_LHC15g1a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC15g1aFinerPtHardBins_$number.root $TOMERGE
        done;
    fi
    
elif [ $2 == "LHC13g" ]; then 

    if [ "$LHC13gData" == "" ]; then 
        HAVELHC13g=0;
    fi
    
    if [ "$LHC15a3aMC" == "" ]; then 
        HAVELHC15a3a=0;
    fi

    if [ "$LHC15a3aplusMC" == "" ]; then 
        HAVELHC15a3aplus=0;
    fi
    
    if [ "$LHC15g2MC" == "" ]; then 
        HAVELHC15g2=0;
    fi
    

    if [ $HAVELHC13g == 1 ]; then
        OUTPUTDIR_LHC13g=$BASEDIR/$TRAINDIR/GA_pp-$LHC13gData
        mkdir -p $OUTPUTDIR_LHC13g
    fi

    if [ $HAVELHC15g2 == 1 ]; then
        OUTPUTDIR_LHC15g2=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g2MC
        mkdir -p $OUTPUTDIR_LHC15g2
    fi
    
    if [ $HAVELHC15a3a == 1 ]; then
        OUTPUTDIR_LHC15a3a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aMC
        mkdir -p $OUTPUTDIR_LHC15a3a
    fi

    if [ $HAVELHC15a3aplus == 1 ]; then
        OUTPUTDIR_LHC15a3aplus=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aplusMC
        mkdir -p $OUTPUTDIR_LHC15a3aplus
    fi
 
    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC13g == 1 ]; then
            echo "downloading LHC13g"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC13gData/merge/root_archive.zip file:$OUTPUTDIR_LHC13g/
            unzip -u $OUTPUTDIR_LHC13g/root_archive.zip -d $OUTPUTDIR_LHC13g/
        fi    
        if [ $HAVELHC15g2 == 1 ]; then
            echo "downloading LHC15g2"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g2MC/merge/root_archive.zip file:$OUTPUTDIR_LHC15g2/
            unzip -u $OUTPUTDIR_LHC15g2/root_archive.zip -d $OUTPUTDIR_LHC15g2/
        fi    

        if [ $HAVELHC15a3a == 1 ]; then
            echo "downloading LHC15a3a"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aMC/merge/root_archive.zip file:$OUTPUTDIR_LHC15a3a/
            unzip -u $OUTPUTDIR_LHC15a3a/root_archive.zip -d $OUTPUTDIR_LHC15a3a/
        fi     
        
        if [ $HAVELHC15a3aplus == 1 ]; then
            echo "downloading LHC15a3a_plus"
            alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aplusMC/merge/root_archive.zip file:$OUTPUTDIR_LHC15a3aplus/
            unzip -u $OUTPUTDIR_LHC15a3aplus/root_archive.zip -d $OUTPUTDIR_LHC15a3aplus/
        fi
    
        echo "copying LHC13gJetJet bins" 
        runNumbers=`cat runNumbersLHC13gJetJet.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat binNumbersJJ.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC15a3a == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC15a3a/$binNumber/$runNumber
                        alien_cp alien:/alice/sim/2015/LHC15a3a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aMC/root_archive.zip file:$OUTPUTDIR_LHC15a3a/$binNumber/$runNumber/
                        unzip -u $OUTPUTDIR_LHC15a3a/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC15a3a/$binNumber/$runNumber/
                    fi
                    if [ $HAVELHC15a3aplus == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber
                        alien_cp alien:/alice/sim/2015/LHC15a3a_plus/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aplusMC/root_archive.zip file:$OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber/
                        unzip -u $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber/root_archive.zip -d $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber/
                    fi    
                done;   
            fi  
        done;
    fi
    
    if [ $HAVELHC13g == 1 ]; then
        ls $OUTPUTDIR_LHC13g/GammaConvV1_*.root > fileLHC13g.txt
        fileNumbers=`cat fileLHC13g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13g/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13g-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13g-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13g_$number.log\"\)
        done;
    fi

    if [ $HAVELHC15a3a == 1 ]; then  
        ls $OUTPUTDIR_LHC15a3a/GammaConvV1_*.root > fileLHC15a3a.txt
        fileNumbers=`cat fileLHC15a3a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC15a3a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC15a3a_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC15a3a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15a3a_$number.log\"\)
        done;
    fi    
    
    if [ $HAVELHC15a3aplus == 1 ]; then  
        ls $OUTPUTDIR_LHC15a3aplus/GammaConvV1_*.root > fileLHC15a3aplus.txt
        fileNumbers=`cat fileLHC15a3aplus.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC15a3aplus/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC15a3aplus_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC15a3aplus_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15a3aplus_$number.log\"\)
        done;
    fi
        
    if [ $MERGEONBINSSingle = 1 ]; then
        binNumbersJJ=`cat binNumbersJJ.txt`
        echo $binNumbersJJ
        ls $OUTPUTDIR_LHC15a3a/GammaConvV1_*.root > filetemp.txt
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC15a3a == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC15a3a/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC15a3a/$binNumber/*/GammaConvV1_$number.root
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC15a3a/$binNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC15a3a$binNumber\_$number.root\"\,\"GammaConvV1_$number\"\)
                fi    
                if [ $HAVELHC15a3aplus == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC15a3aplus/$binNumber/*/GammaConvV1_$number.root
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC15a3aplus/$binNumber/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC15a3aplus$binNumber\_$number.root\"\,\"GammaConvV1_$number\"\)
                fi    
            done;
        done;
    fi
  
    if [ $HAVELHC15g2 == 1 ]; then    
        ls $OUTPUTDIR_LHC15g2/GammaConvV1_*.root > fileLHC15g2.txt
        fileNumbers=`cat fileLHC15g2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC15g2/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC15g2_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC15g2_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15g2_$number.log\"\)
        done;
    fi
    
    rm $OUTPUTDIR/GammaConvV1_MC_LHC15a3a_LHC15a3aplus_*.root
    rm $OUTPUTDIR/GammaConvV1_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_*.root

    ls $OUTPUTDIR/GammaConvV1_MC_LHC15a3a_*.root > filesForMerging.txt
    if [ $MERGEON == 1 ]; then    
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3a_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3aplus_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3a_LHC15a3aplus_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC15a3a_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC15a3aplus_$number.root
            fi
        done
    fi
    
    if [ $MERGEONBINS == 1 ]; then    
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/GammaConvV1_MC_LHC15a3a$bin\_$number.root
                if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC15a3a$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/GammaConvV1_MC_LHC15a3a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3aFinerPtHardBins_$number.root $TOMERGE

            TOMERGE="";     
            for bin in $binsForMerging; do
                if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3aplus$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/GammaConvV1_MC_LHC15a3aplus$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3aplusFinerPtHardBins_$number.root $TOMERGE

            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC15a3aFinerPtHardBins_$number.root
        done;
    fi

    
fi
