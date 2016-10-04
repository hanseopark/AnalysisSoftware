#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs
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
        root -l -b -q -x ChangeStructureToStandardCaloMerged.C\(\"$1\"\,\"$2\"\,\"GammaCaloMerged_$3\"\)
    fi    
}

NSlashes=10
NSlashes2=10
NSlashes3=11

# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1
MERGEONBINSSingle=0
MERGEONBINS=0

# check if train configuration has actually been given
HAVELHC12a=1
HAVELHC12b=1
HAVELHC12c=1
HAVELHC12d=1
HAVELHC12f=1
HAVELHC12h=1
HAVELHC12i=1
HAVELHC16c2=1

# default trainconfigurations
LHC12aData="";
LHC12bData="";
LHC12cData="";
LHC12dData="";
LHC12fData="";
LHC12hData="";
LHC12iData="";
LHC16c2MC="";

if [ $1 = "fbock" ]; then 
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp8TeV
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp8TeV
elif [ $1 = "leardini" ]; then 
   BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then 
   BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
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
   BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/   
elif [ $1 = "mwilde" ]; then        
   BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
   BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pp 
elif [ $1 = "pgonzales" ]; then     
   BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
   BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pp
elif [ $1 = "dmuhlheim" ]; then 
   BASEDIR=/home/daniel/Desktop/Grid
   NSlashes=9
   NSlashes2=8
   NSlashes3=10
elif [ $1 = "bsahlmul" ]; then
  BASEDIR=/Users/sahlmul/Documents/ALICE/mycode/AnalysisSoftware/grid
  NSlashes=11
  NSlashes2=10
  NSlashes3=12
fi


# TRAINDIR=Legotrain-vAN20160627-firstFullTrial
# LHC12aData="1660"; 
# LHC12bData="1661"; 
# LHC12cData="1662"; 
# LHC12dData="1663"; 
# LHC12fData="1664"; 
# LHC12hData="1665"; 
# LHC12iData="1666"; 
# LHC16c2MC="2225";

# TRAINDIR=Legotrain-vAN20160702-withNewMCLabelAlgo
# LHC12aData="1667"; 
# LHC12bData="1668"; 
# LHC12cData="1669"; 
# LHC12dData="1670"; 
# LHC12fData="1671"; 
# LHC12hData="1672"; 
# LHC12iData="1673"; 
# LHC16c2MC="2243";

# TRAINDIR=Legotrain-mCalo-vAN20160705-sys
# LHC12aData="1680"; 
# LHC12bData="1683"; 
# LHC12cData="1686"; 
# LHC12dData="1689"; 
# LHC12fData="1692"; 
# LHC12hData="1695"; 
# LHC12iData="1698"; 
# LHC16c2MC="2263";
# LHC12aData="1681"; 
# LHC12bData="1684"; 
# LHC12cData="1687"; 
# LHC12dData="1690"; 
# LHC12fData="1693"; 
# LHC12hData="1696"; 
# LHC12iData="1699"; 
# LHC16c2MC="2264";
# LHC12aData="1682"; 
# LHC12bData="1685"; 
# LHC12cData="1688"; 
# LHC12dData="1701"; 
# LHC12fData="1694"; 
# LHC12hData="1697"; 
# LHC12iData="1700"; 
# LHC16c2MC="2265";
# LHC16c2MC="2266";


# TRAINDIR=Legotrain-mCalo-vAN20160812-V1Clus
# LHC12aData="1763"; 
# LHC12bData="1764"; 
# LHC12cData="1765"; 
# LHC12dData="1766"; 
# LHC12fData="1767"; 
# LHC12hData="1768"; 
# LHC12iData="1769"; 
# # LHC16c2MC="2403";
# LHC16c2MC="2404";
# LHC16c2MC="2405";

# TRAINDIR=Legotrain-mCalo-vAN20160812-V2ClusNewSecTreatWithSys
# LHC12aData="1770"; 
# LHC12bData="1771"; 
# LHC12cData="1772"; 
# LHC12dData="1773"; 
# LHC12fData="1774"; 
# LHC12hData="1775"; 
# LHC12iData="1776"; 
# LHC16c2MC="2402";

TRAINDIR=Legotrain-mCalo-vAN20160831-V2ClusSys
LHC12aData="1829"; 
LHC12bData="1830"; 
LHC12cData="1831"; 
LHC12dData="1832"; 
LHC12fData="1833"; 
LHC12hData="1842"; 
LHC12iData="1835"; 
# LHC16c2MC="2502";
# LHC16c2MC="2503";
# LHC16c2MC="2504";
# LHC16c2MC="2505";
LHC16c2MC="2538";


OUTPUTDIR=$BASEDIR/$TRAINDIR

    
if [ "$LHC12aData" == "" ]; then 
    HAVELHC12a=0;
fi
if [ "$LHC12bData" == "" ]; then 
    HAVELHC12b=0;
fi
if [ "$LHC12cData" == "" ]; then 
    HAVELHC12c=0;
fi
if [ "$LHC12dData" == "" ]; then 
    HAVELHC12d=0;
fi
if [ "$LHC12fData" == "" ]; then 
    HAVELHC12f=0;
fi
if [ "$LHC12hData" == "" ]; then 
    HAVELHC12h=0;
fi
if [ "$LHC12iData" == "" ]; then 
    HAVELHC12i=0;
fi
if [ "$LHC16c2MC" = "" ]; then 
    HAVELHC16c2=0; 
fi

if [ $HAVELHC12a == 1 ]; then
    LHC12aData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC12aData\_`
    if [ "$LHC12aData" == "" ]; then 
        HAVELHC12a=0;
    else     
        OUTPUTDIR_LHC12a=$BASEDIR/$TRAINDIR/GA_pp-$LHC12aData
    fi
fi
if [ $HAVELHC12b == 1 ]; then
    LHC12bData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC12bData\_`
    if [ "$LHC12bData" == "" ]; then 
        HAVELHC12b=0;
    else     
        OUTPUTDIR_LHC12b=$BASEDIR/$TRAINDIR/GA_pp-$LHC12bData
    fi
fi
if [ $HAVELHC12c == 1 ]; then
    LHC12cData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC12cData\_`
    if [ "$LHC12cData" == "" ]; then 
        HAVELHC12c=0;
    else     
        OUTPUTDIR_LHC12c=$BASEDIR/$TRAINDIR/GA_pp-$LHC12cData
    fi
fi
if [ $HAVELHC12d == 1 ]; then
    LHC12dData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC12dData\_`
    if [ "$LHC12dData" == "" ]; then 
        HAVELHC12d=0;
    else     
        OUTPUTDIR_LHC12d=$BASEDIR/$TRAINDIR/GA_pp-$LHC12dData
    fi
fi
if [ $HAVELHC12f == 1 ]; then
    LHC12fData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC12fData\_`
    if [ "$LHC12fData" == "" ]; then 
        HAVELHC12f=0;
    else     
        OUTPUTDIR_LHC12f=$BASEDIR/$TRAINDIR/GA_pp-$LHC12fData
    fi
fi
if [ $HAVELHC12h == 1 ]; then
    LHC12hData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC12hData\_`
    if [ "$LHC12hData" == "" ]; then 
        HAVELHC12h=0;
    else     
        OUTPUTDIR_LHC12h=$BASEDIR/$TRAINDIR/GA_pp-$LHC12hData
    fi
fi
if [ $HAVELHC12i == 1 ]; then
    LHC12iData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC12iData\_`
    if [ "$LHC12iData" == "" ]; then 
        HAVELHC12i=0;
    else     
        OUTPUTDIR_LHC12i=$BASEDIR/$TRAINDIR/GA_pp-$LHC12iData
    fi
fi
if [ $HAVELHC16c2 == 1 ]; then
    LHC16c2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC16c2MC\_`
    if [ "$LHC16c2MC" == "" ]; then 
        HAVELHC16c2=0;
    else     
        OUTPUTDIR_LHC16c2=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC16c2MC
    fi    
fi

if [ $DOWNLOADON == 1 ]; then
    if [ $HAVELHC12a == 1 ]; then
        echo "downloading LHC12a"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12a/INT7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12aData/merge_runlist_9"
#         CopyFileIfNonExisitent $OUTPUTDIR_LHC12a/EMC7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12aData/merge_runlist_10"
#         CopyFileIfNonExisitent $OUTPUTDIR_LHC12a/EGA "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12aData/merge_runlist_11"
#         runNumbers=`cat runNumbersLHC12a_pass2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC12a/$runNumber "/alice/data/2012/LHC12a/000$runNumber/pass2/PWGGA/GA_pp/$LHC12aData"
#         done;
    fi    
    if [ $HAVELHC12b == 1 ]; then
        echo "downloading LHC12b"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12b/INT7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12bData/merge_runlist_9"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12b/EMC7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12bData/merge_runlist_10"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12b/EGA "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12bData/merge_runlist_11"
#         runNumbers=`cat runNumbersLHC12b_pass2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC12b/$runNumber "/alice/data/2012/LHC12b/000$runNumber/pass2/PWGGA/GA_pp/$LHC12bData"
#         done;
    fi    
    if [ $HAVELHC12c == 1 ]; then
        echo "downloading LHC12c"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12c/INT7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12cData/merge_runlist_9"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12c/EMC7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12cData/merge_runlist_10"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12c/EGA "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12cData/merge_runlist_11"
#         runNumbers=`cat runNumbersLHC12c_pass2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC12c/$runNumber "/alice/data/2012/LHC12c/000$runNumber/pass2/PWGGA/GA_pp/$LHC12cData"
#         done;
    fi    
    if [ $HAVELHC12d == 1 ]; then
        echo "downloading LHC12d"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12d/INT7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12dData/merge_runlist_9"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12d/EMC7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12dData/merge_runlist_10"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12d/EGA "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12dData/merge_runlist_11"
#         runNumbers=`cat runNumbersLHC12d_pass2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC12d/$runNumber "/alice/data/2012/LHC12d/000$runNumber/pass2/PWGGA/GA_pp/$LHC12dData"
#         done;
    fi    
    if [ $HAVELHC12f == 1 ]; then
        echo "downloading LHC12f"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12f/INT7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12fData/merge_runlist_9"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12f/EMC7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12fData/merge_runlist_10"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12f/EGA "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12fData/merge_runlist_11"
#         runNumbers=`cat runNumbersLHC12f_pass2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC12f/$runNumber "/alice/data/2012/LHC12f/000$runNumber/pass2/PWGGA/GA_pp/$LHC12fData"
#         done;
    fi    
    if [ $HAVELHC12h == 1 ]; then
        echo "downloading LHC12h"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12h/INT7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12hData/merge_runlist_9"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12h/EMC7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12hData/merge_runlist_10"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12h/EGA "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12hData/merge_runlist_11"
#         runNumbers=`cat runNumbersLHC12h_pass2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC12h/$runNumber "/alice/data/2012/LHC12h/000$runNumber/pass2/PWGGA/GA_pp/$LHC12hData"
#         done;
    fi    
    if [ $HAVELHC12i == 1 ]; then
        echo "downloading LHC12i"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12i/INT7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12iData/merge_runlist_9"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12i/EMC7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12iData/merge_runlist_10"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC12i/EGA "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC12iData/merge_runlist_11"
#         runNumbers=`cat runNumbersLHC12i_pass2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC12i/$runNumber "/alice/data/2012/LHC12i/000$runNumber/pass2/PWGGA/GA_pp/$LHC12iData"
#         done;
    fi    
    if [ $HAVELHC16c2 == 1 ]; then
        echo "downloading LHC16c2"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC16c2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC16c2MC/merge"

#         echo "copying LHC12xJetJet" 
#         runNumbers=`cat runNumbersLHC16c2.txt`
#         echo $runNumbers
#         for runNumber in $runNumbers; do
#             echo $runNumber
#             binNumbersJJ=`cat binsJetJetLHC16c2.txt`
#             for binNumber in $binNumbersJJ; do
#                 echo $binNumber
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC16c2/$binNumber/$runNumber "/alice/sim/2016/LHC16c2/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC16c2MC"
#             done;   
#         done;
    fi    
fi

mkdir $OUTPUTDIR/data-SinglePeriods

if [ $HAVELHC12a == 1 ]; then
    ls $OUTPUTDIR_LHC12a/*/GammaCaloMerged_*.root > fileLHC12a.txt
    fileNumbers=`cat fileLHC12a.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12a/INT7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12a-pass2-INT7Runlist_$number.root $number
#         ChangeStructureIfNeeded $OUTPUTDIR_LHC12a/EMC7/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_LHC12a-pass2-EMC7Runlist_$number.root $number
#         ChangeStructureIfNeeded $OUTPUTDIR_LHC12a/EGA/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_LHC12a-pass2-EGARunlist_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12a-pass2-INT7Runlist_$number.root\"\,\"$OUTPUTDIR/data-SinglePeriods/CutSelection_LHC12a_$number.log\"\)
    done;
fi

if [ $HAVELHC12b == 1 ]; then
    ls $OUTPUTDIR_LHC12b/*/GammaCaloMerged_*.root > fileLHC12b.txt
    fileNumbers=`cat fileLHC12b.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12b/INT7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12b-pass2-INT7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12b/EMC7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12b-pass2-EMC7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12b/EGA/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12b-pass2-EGARunlist_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12b-pass2-INT7Runlist_$number.root\"\,\"$OUTPUTDIR/data-SinglePeriods/CutSelection_LHC12b_$number.log\"\)
    done;
fi
if [ $HAVELHC12c == 1 ]; then
    ls $OUTPUTDIR_LHC12c/*/GammaCaloMerged_*.root > fileLHC12c.txt
    fileNumbers=`cat fileLHC12c.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12c/INT7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12c-pass2-INT7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12c/EMC7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12c-pass2-EMC7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12c/EGA/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12c-pass2-EGARunlist_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12c-pass2-INT7Runlist_$number.root\"\,\"$OUTPUTDIR/data-SinglePeriods/CutSelection_LHC12c_$number.log\"\)
    done;
fi
if [ $HAVELHC12d == 1 ]; then
    ls $OUTPUTDIR_LHC12d/*/GammaCaloMerged_*.root > fileLHC12d.txt
    fileNumbers=`cat fileLHC12d.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12d/INT7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12d-pass2-INT7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12d/EMC7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12d-pass2-EMC7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12d/EGA/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12d-pass2-EGARunlist_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12d-pass2-INT7Runlist_$number.root\"\,\"$OUTPUTDIR/data-SinglePeriods/CutSelection_LHC12d_$number.log\"\)
    done;
fi
if [ $HAVELHC12f == 1 ]; then
    ls $OUTPUTDIR_LHC12f/*/GammaCaloMerged_*.root > fileLHC12f.txt
    fileNumbers=`cat fileLHC12f.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12f/INT7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12f-pass2-INT7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12f/EMC7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12f-pass2-EMC7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12f/EGA/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12f-pass2-EGARunlist_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12f-pass2-INT7Runlist_$number.root\"\,\"$OUTPUTDIR/data-SinglePeriods/CutSelection_LHC12f_$number.log\"\)
    done;
fi
if [ $HAVELHC12h == 1 ]; then
    ls $OUTPUTDIR_LHC12h/*/GammaCaloMerged_*.root > fileLHC12h.txt
    fileNumbers=`cat fileLHC12h.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12h/INT7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12h-pass2-INT7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12h/EMC7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12h-pass2-EMC7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12h/EGA/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12h-pass2-EGARunlist_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12h-pass2-INT7Runlist_$number.root\"\,\"$OUTPUTDIR/data-SinglePeriods/CutSelection_LHC12h_$number.log\"\)
    done;
fi

if [ $HAVELHC12i == 1 ]; then
    ls $OUTPUTDIR_LHC12i/*/GammaCaloMerged_*.root > fileLHC12i.txt
    fileNumbers=`cat fileLHC12i.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12i/INT7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12i-pass2-INT7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12i/EMC7/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12i-pass2-EMC7Runlist_$number.root $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC12i/EGA/GammaCaloMerged_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12i-pass2-EGARunlist_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12i-pass2-INT7Runlist_$number.root\"\,\"$OUTPUTDIR/data-SinglePeriods/CutSelection_LHC12i_$number.log\"\)
    done;
fi


if [ $HAVELHC16c2 == 1 ]; then  
    ls $OUTPUTDIR_LHC16c2/GammaCaloMerged_*.root > fileLHC16c2.txt
    fileNumbers=`cat fileLHC16c2.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC16c2/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC16c2_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC16c2_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC16c2_$number.log\"\)
    done;
    
    binNumbersJJ=`cat binsJetJetLHC16c2.txt`
    echo $binNumbersJJ
    if [ $MERGEONBINSSingle = 1 ]; then
        ls $OUTPUTDIR_LHC16c2/GammaCaloMerged_*.root > filetemp.txt
        mkdir $OUTPUTDIR/LHC16c2FineBins
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                hadd -f $OUTPUTDIR_LHC16c2/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC16c2/$binNumber/*/GammaCaloMerged_$number.root
                ChangeStructureIfNeeded $OUTPUTDIR_LHC16c2/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC16c2FineBins/GammaCaloMerged_MC_LHC16c2$binNumber\_$number.root $number
            done;
        done;
    fi
fi


if [ $MERGEON == 1 ]; then
    rm $OUTPUTDIR/GammaCaloMerged_data_*.root
    ls $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12i-pass2-INT7Runlist_*.root > configsForMerging.txt
    configsForMerging=`cat configsForMerging.txt`
    for configForMerging in $configsForMerging; do
        echo $configForMerging
        number=`echo $configForMerging  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
        echo $number
        ls $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12*INT7Runlist_$number.root > filesForMerging.txt
        sort filesForMerging.txt
        OUTPUTNAME="";
        fileNames=`cat filesForMerging.txt`
        for fileName in $fileNames; do
            periodName=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "-" -f1 | cut -d "2" -f2`
            echo $periodName
            OUTPUTNAME+="$periodName";
        done;
        hadd -f $OUTPUTDIR/GammaCaloMerged_data_LHC12$OUTPUTNAME-pass2-INT7Runlist_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12*INT7Runlist_$number.root
    done;
    ls $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12i-pass2-EMC7Runlist_*.root > configsForMerging.txt
    configsForMerging=`cat configsForMerging.txt`
    for configForMerging in $configsForMerging; do
        echo $configForMerging
        number=`echo $configForMerging  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
        echo $number        
        ls $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12*EMC7Runlist_$number.root > filesForMerging.txt
        sort filesForMerging.txt
        OUTPUTNAME="";
        fileNames=`cat filesForMerging.txt`
        for fileName in $fileNames; do
            periodName=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "-" -f1 | cut -d "2" -f2`
            echo $periodName
            OUTPUTNAME+="$periodName";
        done;
        hadd -f $OUTPUTDIR/GammaCaloMerged_data_LHC12$OUTPUTNAME-pass2-EMC7Runlist_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12*EMC7Runlist_$number.root
    done;
    ls $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12i-pass2-EGARunlist_*.root > configsForMerging.txt
    configsForMerging=`cat configsForMerging.txt`
    for configForMerging in $configsForMerging; do
        echo $configForMerging
        number=`echo $configForMerging  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
        ls $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12*EGARunlist_$number.root > filesForMerging.txt
        sort filesForMerging.txt
        OUTPUTNAME="";
        fileNames=`cat filesForMerging.txt`
        for fileName in $fileNames; do
            periodName=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "-" -f1 | cut -d "2" -f2`
            echo $periodName
            OUTPUTNAME+="$periodName";
        done;
        hadd -f $OUTPUTDIR/GammaCaloMerged_data_LHC12$OUTPUTNAME-pass2-EGARunlist_$number.root $OUTPUTDIR/data-SinglePeriods/GammaCaloMerged_LHC12*EGARunlist_$number.root
    done;
fi

if [ $MERGEONBINS == 1 ]; then
    rm $OUTPUTDIR/GammaCaloMerged_MC_LHC16c2Remerged_*.root
    ls $OUTPUTDIR/GammaCaloMerged_MC_LHC16c2_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        binsForMerging=`cat binsJetJetLHC16c2.txt`
        number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
        TOMERGE="";
        for bin in $binsForMerging; do
            echo $OUTPUTDIR/LHC16c2FineBins/GammaCaloMerged_MC_LHC16c2$bin\_$number.root
            if [ -f $OUTPUTDIR/LHC16c2FineBins/GammaCaloMerged_MC_LHC16c2$bin\_$number.root ]; then
                TOMERGE="$TOMERGE $OUTPUTDIR/LHC16c2FineBins/GammaCaloMerged_MC_LHC16c2$bin"
                TOMERGE+="_$number.root"
            else 
                echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC16c2FineBins/GammaCaloMerged_MC_LHC16c2$bin\_$number.root";
            fi
        done;
        hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC16c2Remerged_$number.root $TOMERGE
    done;
fi

