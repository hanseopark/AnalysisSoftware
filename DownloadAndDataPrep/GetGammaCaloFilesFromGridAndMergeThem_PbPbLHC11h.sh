#! /bin/bash

function CopyFileIfNonExisitent()
{
    if [ -f $1/root_archive.zip ] && [ -s $1/root_archive.zip ]; then
        echo "$1/root_archive.zip exists";
    else
        mkdir -p $1
        alien_cp alien:$2/root_archive.zip file:$1/
    fi
    unzip -u $1/root_archive.zip -d $1/
}

function ChangeStructureIfNeeded()
{
    if [ -f $2 ]; then
        echo "already changed"
    else
        mv $1 $2
    fi
}



# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1

# check if train configuration has actually been given
HAVELHC11h=1
HAVELHC14a1a=1
HAVELHC14a1b=1
HAVELHC14a1c=1

# default trainconfigurations
LHC11hData="";
LHC14a1aMC="";
LHC14a1bMC="";
LHC14a1cMC="";

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/PbPb
elif [ $1 = "leardini" ]; then
   BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then
   BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
elif [ $1 = "leardiniGSI" ]; then
   BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/PbPb
elif [ $1 = "dmuhlheim" ]; then
   BASEDIR=/home/daniel/Desktop/Grid
elif [ $1 = "loizides" ]; then
   BASEDIR=/home/loizides/Downloads
fi

# Definitition of number of slashes in your path to different depths
NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
NSlashes=`expr $NSlashesBASE + 4`
NSlashes2=`expr $NSlashes - 1`
NSlashes3=`expr $NSlashes + 1`
NSlashes4=`expr $NSlashes + 2`
echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"

#   TRAINDIR=Legotrain-vAN-20141110-Calo
#   LHC11hData=144_20141111-1840; #ESD
#   LHC14a1a=176_20141111-1928; #ESD
#   LHC14a1b=177_20141111-2033; #ESD
#   LHC14a1c=;

#   TRAINDIR=Legotrain-vAN-20141116-Calo
#   LHC11hData=147_20141117-1517; #ESD
#   LHC14a1a=186_20141124-1555;
#   LHC14a1b=187_20141124-1556; #ESD

# TRAINDIR=Legotrain-vAN20171215_EMCrerun
# LHC11hData=364;
# LHC14a1aMC=715;

TRAINDIR=Legotrain-vAN20180125_EMCtrackMult
LHC11hData=378;
LHC14a1aMC=801;

OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ "$LHC11hData" == "" ]; then
    HAVELHC11h=0;
fi
if [ "$LHC14a1aMC" = "" ]; then
    HAVELHC14a1a=0;
fi
if [ "$LHC14a1aMC" = "" ]; then
    HAVELHC14a1b=0;
fi
if [ "$LHC14a1aMC" = "" ]; then
    HAVELHC14a1c=0;
fi

if [ $HAVELHC11h == 1 ]; then
    LHC11hData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/ | grep $LHC11hData\_`
    OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/GA_PbPb-$LHC11hData
fi
if [ $HAVELHC14a1a == 1 ]; then
    LHC14a1aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/ | grep $LHC14a1aMC\_`
    OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC14a1aMC
fi
if [ $HAVELHC14a1b == 1 ]; then
    LHC14a1bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/ | grep $LHC14a1bMC\_`
    OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC14a1bMC
fi
if [ $HAVELHC14a1c == 1 ]; then
    LHC14a1cMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/ | grep $LHC14a1cMC\_`
    OUTPUTDIR_LHC14a1c=$BASEDIR/$TRAINDIR/GA_PbPb_MC-$LHC14a1cMC
fi

if [ $DOWNLOADON == 1 ]; then
    if [ $HAVELHC11h == 1 ]; then
        echo "downloading LHC11h"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC11h/runList7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/$LHC11hData/merge_runlist_7"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC11h/runList8 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/$LHC11hData/merge_runlist_8"
        ls $OUTPUTDIR_LHC11h/runList7/GammaCalo_*.root > fileLHC11h2.txt
        fileNumbers=`cat fileLHC11h2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            hadd -f $OUTPUTDIR_LHC11h/GammaCalo_$number.root $OUTPUTDIR_LHC11h/runList7/GammaCalo_$number.root $OUTPUTDIR_LHC11h/runList8/GammaCalo_$number.root
        done;
    fi
    if [ $HAVELHC14a1a == 1 ]; then
        echo "downloading LHC14a1a"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC14a1a/runList7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC14a1aMC/merge_runlist_7"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC14a1a/runList8 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC14a1aMC/merge_runlist_8"
        ls $OUTPUTDIR_LHC14a1a/runList7/GammaCalo_*.root > fileLHC11h2.txt
        fileNumbers=`cat fileLHC11h2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            hadd -f $OUTPUTDIR_LHC14a1a/GammaCalo_$number.root $OUTPUTDIR_LHC14a1a/runList7/GammaCalo_$number.root $OUTPUTDIR_LHC14a1a/runList8/GammaCalo_$number.root
        done;
    fi
#     if [ $HAVELHC14a1b == 1 ]; then
#         echo "downloading LHC14a1b"
#         CopyFileIfNonExisitent $OUTPUTDIR_LHC14a1b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC14a1bMC/merge_runlist_1"
#     fi
#     if [ $HAVELHC14a1c == 1 ]; then
#         echo "downloading LHC14a1c"
#         CopyFileIfNonExisitent $OUTPUTDIR_LHC14a1c "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/$LHC14a1cMC/merge_runlist_1"
#     fi
fi

if [ $HAVELHC11h == 1 ]; then
    ls $OUTPUTDIR_LHC11h/GammaCalo_*.root > fileLHC11h.txt
    fileNumbers=`cat fileLHC11h.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC11h/GammaCalo_$number.root $OUTPUTDIR/GammaCalo_LHC11h-pass2_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCalo_LHC11h_$number.log\"\,4\)
    done;
fi

if [ $HAVELHC14a1a == 1 ]; then
    ls $OUTPUTDIR_LHC14a1a/GammaCalo_*.root > fileLHC14a1a.txt
    fileNumbers=`cat fileLHC14a1a.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC14a1a/GammaCalo_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14a1a_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14a1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCalo_MC_LHC14a1a_$number.log\"\,4\)
    done;
fi
#
# if [ $HAVELHC14a1b == 1 ]; then
#     ls $OUTPUTDIR_LHC14a1b/GammaCalo_*.root > fileLHC14a1b.txt
#     fileNumbers=`cat fileLHC14a1b.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         ChangeStructureIfNeeded $OUTPUTDIR_LHC14a1b/GammaCalo_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14a1b_$number.root $number
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14a1b_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCalo_MC_LHC14a1b_$number.log\"\,4\)
#     done;
# fi
#
# if [ $HAVELHC14a1c == 1 ]; then
#     ls $OUTPUTDIR_LHC14a1c/GammaCalo_*.root > fileLHC14a1c.txt
#     fileNumbers=`cat fileLHC14a1c.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         ChangeStructureIfNeeded $OUTPUTDIR_LHC14a1c/GammaCalo_$number.root $OUTPUTDIR/GammaCalo_MC_LHC14a1c_$number.root $number
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC14a1c_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCalo_MC_LHC14a1c_$number.log\"\,4\)
#     done;
# fi
