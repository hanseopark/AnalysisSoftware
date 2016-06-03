#! /bin/bash

# This script has to be run with "bash"

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
        root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$1\"\,\"$2\"\,\"GammaConvCalo_$3\"\)
    fi    
}


# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1
MERGEONBINSSingle=1
MERGEONBINS=1

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
# merges files according to the pp needs
NSlashes=10
NSlashes2=9

if [ $1 = "fbock" ]; then 
   BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp
   NSlashes=8
   NSlashes2=7
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp
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
elif [ $1 = "dmuhlhei" ]; then 
   BASEDIR=~/data/work/Grid
   NSlashes=9
   NSlashes2=8
fi


# TRAINDIR=Legotrain-ConvCalo-LHC11a-Sys-MB #
# LHC11aData="946"; #1,3,4,5
# LHC11aData="952"; #6,7,8,9 
# LHC11aData="948"; #10,11,12,14
# LHC12f1aMC="988"; #1,3,4,5
# LHC12f1aMC="989"; #6,7,8,9
# LHC12f1aMC="990"; #10,11,12,14
# LHC12f1bMC="1046"; #1,3,4,5
# LHC12f1bMC="992"; #6,7,8,9
# LHC12f1bMC="993"; #10,11,12,14
# LHC12i3MC="975"; #1,3,4,5
# LHC12i3MC="976"; #6,7,8,9
# LHC12i3MC="977"; #10,11,12,14
# LHC15g1aMC="1168";  #1,3,4 - merging jobs failed - fixed
# LHC15g1aMC="1169";  #5,6,7 - merging jobs failed - fixed
# LHC15g1aMC="1170";  #8,9,10 - fixed
# LHC15g1aMC="1171";  #11,12,14 -fixed

# TRAINDIR=Legotrain-ConvCalo-LHC11a-Sys-EMC1
# LHC11aData="951"; #16,17,18,19
# LHC11aData="1020"; #20,21,22,23 
# LHC11aData="953"; #24,25,26
# LHC11aData="955"; #27
# LHC12f1aMC="1043"; #16,17,18,19
# LHC12f1aMC="995"; #20,21,22,23
# LHC12f1aMC="1148"; #24,25,26,27 
# LHC12f1bMC="1045"; #16,17,18,19
# LHC12f1bMC="1149"; #20,21,22,23 
# LHC12f1bMC="999"; #24,25,26,27
# LHC12i3MC=""; #
# LHC15g1aMC="1191";  #16,17,18 - merging jobs failed - fixed
# LHC15g1aMC="1192";  #19,20,21 - merging jobs failed - fixed
# LHC15g1aMC="1193";  #22,23,24 - merging jobs failed
# LHC15g1aMC="1194";  #25,26,27

# TRAINDIR=Legotrain-ConvCalo-LHC11a-Sys-lowE
# LHC11aData="1021"; #2,15
# LHC12f1aMC="1000"; #2,15
# LHC12f1bMC="1001"; #2,15
# LHC12i3MC="1002"; # 2,15
# LHC15g1aMC="1181";  #2, 15

# TRAINDIR=Legotrain-ConvCalo-LHC11a-Sys-V1Cluster
# LHC11aData="1009"; #2,15
# LHC12f1aMC="1182"; #2,15
# LHC12f1bMC="1183"; #2,15
# LHC12i3MC="1184"; # 2,15
# LHC15g1aMC="1180";  #2, 15

# TRAINDIR=Legotrain-ConvCalo-LHC13g-Sys-INT7
# LHC13gData="956_20150906-1145"; #41,43,44,45
# LHC13gData="957_20150906-1146"; #46,47,48,49
# LHC13gData="958_20150906-1147"; #50,51,52,54
# LHC15a3aMC="1172_20151009-1721";  #43,44,45 
# LHC15a3aMC="1173_20151009-1721";  #46,47,48 - merging jobs failed - fixed
# LHC15a3aMC="1174_20151009-1729";  #49,50,51
# LHC15a3aMC="1175_20151009-1727";  #52,54
# LHC15a3aMC="1240_20151015-2232";  #41
# LHC15a3aplusMC="1176_20151009-1732";  #43,44,45 - merging jobs failed - fixed
# LHC15a3aplusMC="1177_20151009-1725";  #46,47,48 - merging jobs failed - fixed
# LHC15a3aplusMC="1178_20151009-1733";  #49,50,51
# LHC15a3aplusMC="1179_20151009-1731";  #52,54
# LHC15a3aplusMC="1241_20151015-2233";  #41
# LHC15g2MC="1009_20150909-1232"; #43,44,45,46
# LHC15g2MC="1152_20151006-1511"; #47,48,49,50 
# LHC15g2MC="1011_20150909-1416"; #51,52,54
# LHC15g2MC="1242_20151015-2234"; #41


# TRAINDIR=Legotrain-ConvCalo-LHC13g-Sys-EMC7
# LHC13gData="960"; # 56,57,58,59 
# LHC13gData="961"; # 60,61,62,64 
# LHC13gData="962"; # 65,66,68
# LHC13gData="1024"; #63
# LHC15a3aMC="1199"; # 56,57,58 
# LHC15a3aMC="1200"; # 59,60,61 
# LHC15a3aMC="1201"; # 62,63,64 
# LHC15a3aMC="1202"; # 65,68 
# LHC15a3aplusMC="1195";  # 56,57,58 
# LHC15a3aplusMC="1196";  # 59,60,61
# LHC15a3aplusMC="1197";  # 62,63,64 
# LHC15a3aplusMC="1198";  # 65,68 

# TRAINDIR=Legotrain-ConvCalo-LHC13g-Sys-EG2
# LHC13gData="963"; # 70,71,72,73
# LHC13gData="964"; # 74,75,76,77
# LHC13gData="965"; # 78,79,81
# LHC15a3aMC="1207"; # 70,71,72
# LHC15a3aMC="1208"; # 73,74,75 
# LHC15a3aMC="1209"; # 76,77,78 
# LHC15a3aMC="1210"; # 79, 81 
# LHC15a3aplusMC="1203";  # 70,71,72
# LHC15a3aplusMC="1204";  # 73,74,75
# LHC15a3aplusMC="1205";  # 76,77,78
# LHC15a3aplusMC="1206";  # 79,81

# TRAINDIR=Legotrain-ConvCalo-LHC13g-Sys-EG1
# LHC13gData="966"; # 83,84,85,86 
# LHC13gData="967"; # 87,88,89,90
# LHC13gData="968"; # 91,92,94
# LHC15a3aMC="1211"; # 83,84,85
# LHC15a3aMC="1212"; # 86,87,88 
# LHC15a3aMC="1213"; # 89,90,91
# LHC15a3aMC="1214"; # 92,94
# LHC15a3aplusMC="1215"; # 83,84,85
# LHC15a3aplusMC="1216"; # 86,87,88 
# LHC15a3aplusMC="1217"; # 89,90,91
# LHC15a3aplusMC="1218"; # 92,94

# TRAINDIR=Legotrain-ConvCalo-LHC13g-Sys-lowE - downloaded
# LHC13gData="959"; # 42,55,69,82
# LHC15a3aMC="1186";  # 42,55
# LHC15a3aplusMC="1185";  # 42,55,69,82
# LHC15g2MC="1187"; # 42,55,69,82

# TRAINDIR=Legotrain-ConvCalo-LHC13g-Sys-V1Cluster - downloaded
# LHC13gData="1022"; # 41
# LHC15a3aMC="1189";  # 41
# LHC15a3aplusMC="1190";  # 41
# LHC15g2MC="1188"; # 41

# TRAINDIR=Legotrain-vAN-20151019-2.76TeV-QA
#LHC13gData="1034"; # 95,96
#LHC15a3aMC="1249";  # 95,96
#LHC15a3aplusMC="1250";  # 95,96
#LHC15g2MC="1248"; # 95,96
# LHC11aData="1035"; #1
# LHC12f1aMC="1245"; #1
# LHC12f1bMC="1246"; #1
# LHC15g1aMC="1247";  #1

# TRAINDIR=Legotrain-vAN-20151025-MultiplicityDependence
# # LHC11aData="371";
# LHC12f1aMC="1258"; #calo98, convCalo98
# LHC12f1bMC="1259"; #calo98, convCalo98
# LHC15g1aMC="1262"; #calo98, convCalo98
# 
# # LHC15a3aMC="1260"; 
# # LHC15a3aplusMC="1261"; 
# # LHC15g2MC="1263"; 
# LHC15g2MC="1264"; 

# TRAINDIR=Legotrain-vAN-20151103-MultiplicityDependence2
# # LHC11aData="371";
# LHC12f1aMC="1273"; #calo98, convCalo98
# LHC12f1bMC="1274"; #calo98, convCalo98
# LHC15g1aMC="1275"; #calo98, convCalo98

# LHC15a3aMC="1260"; 
# LHC15a3aplusMC="1261"; 
# LHC15g2MC="1263"; 
# LHC15g2MC="1264"; 

# TRAINDIR=Legotrain-vAN-20151106-MultiplicityDependence3
# # LHC11aData="371";
# LHC12f1aMC="1278"; 
# LHC12f1bMC="1279"; 
# # LHC15g1aMC="1275"; 

# TRAINDIR=Legotrain-vAN-20151101-NonLinearityStudies
# LHC11aData="1059";
# LHC12f1aMC="1270"; 
# LHC12f1bMC="1271"; 
# LHC15g1aMC="1275"; 

# TRAINDIR=Legotrain-vAN-20151211-MultiplicityDependence4
# LHC12f1aMC="1364"; 
# LHC12f1bMC="1365"; 
# LHC15g2MC="1363"; 

# TRAINDIR=Legotrain-vAN-20151213-CalibDefaultTiming
# LHC11aData="1125";
# LHC12f1aMC="1360"; 
# LHC12f1bMC="1361"; 
# LHC15g1aMC="1388";
# LHC12i3MC="1385";
  
# LHC13gData="1126"
# LHC15g2MC="1362";

# TRAINDIR=Legotrain-vAN-20151213-CalibOpenTiming
# LHC11aData="1127";
# LHC12f1aMC="1360"; 
# LHC12f1bMC="1361"; 

# LHC13gData="1128"
# LHC15g2MC="1362";

# TRAINDIR=Legotrain-vAN-20151213-CalibOpenTimingDeltaT
# LHC11aData="1129";
# LHC12f1aMC="1360"; 
# LHC12f1bMC="1361"; 
# 
# LHC13gData="1130"
# LHC15g2MC="1362";

# TRAINDIR=Legotrain-vAN-20151218-TriggersWithoutPileupAndPCMTriggers
# LHC13gData="1136"
# LHC15g2MC="1373";
# LHC15a3aMC="1377"; 
# LHC15a3aplusMC="1375"; 

# LHC13gData="1137"
# LHC15g2MC="1374";
# LHC15a3aMC="1378"; 
# LHC15a3aplusMC="1376"; 
# LHC11aData="1135";
# LHC12f1aMC="1371"; 
# LHC12f1bMC="1372"; 
# LHC15g1aMC="1379"; 

# TRAINDIR=Legotrain-vAN-20160108-CalibOpenTiming
# LHC11aData="1169";
# LHC12f1aMC="1467"; 
# LHC12f1bMC="1468"; 
# LHC15g1aMC="1470";
# LHC12i3MC="1469";

# LHC13gData="1170"
# LHC15g2MC="1475";
# LHC15a3aMC="1472"; 
# LHC15a3aplusMC="1473"; 

# TRAINDIR=Legotrain-vAN-20160108-CalibOpenTimingDeltaT
# LHC11aData="1171";
# # LHC12f1aMC="1467"; 
# # LHC12f1bMC="1468"; 
# # LHC15g1aMC="1470";
# # LHC12i3MC="1469";
# 
# LHC13gData="1172"
# # LHC15g2MC="1475";
# # LHC15a3aMC="1472"; 
# # LHC15a3aplusMC="1473"; 

# TRAINDIR=Legotrain-vAN-20160124-NLCorrV1
# LHC11aData="1176";
# LHC12f1aMC="1526"; 
# LHC12f1bMC="1527"; 
# LHC12i3MC="1528";
# LHC15g1aMC="1529";
# 
# # LHC13gData="1179"
# # LHC15g2MC="1530";
# # LHC15g2MC="1531";
# # LHC15a3aMC="1532"; 
# # LHC15a3aMC="1535"; 
# LHC15a3aMC="LHC15a3a"; 
# # LHC15a3aplusMC="1533"; 
# # LHC15a3aplusMC="1534"; 
# LHC15a3aplusMC="LHC15a3aplus"; 

# TRAINDIR=Legotrain-SysAdd-ConvCalo
# # LHC11aData="1209";
# LHC11aData="1210";
# LHC12f1aMC="1575"; 
# LHC12f1bMC="1577"; 
# # LHC12f1aMC="1576"; 
# # LHC12f1bMC="1578"; 
# # LHC15g1aMC="1579";
# LHC15g1aMC="1595";
# 
# # LHC13gData="1214"
# # LHC13gData="1215"
# # LHC13gData="1216"
# # LHC15g2MC="1580";
# # LHC15g2MC="1612";
# # LHC15a3aMC="1581"; 
# # LHC15a3aplusMC="1584"; 
# # LHC15a3aMC="1593"; 
# # LHC15a3aplusMC="1594"; 
# # LHC15a3aMC="1581"; 
# # LHC15a3aplusMC="1584"; 
# LHC15a3aMC="1582";
# LHC15a3aplusMC="1585"; 
# # LHC15a3aMC="1583"; 
# # LHC15a3aplusMC="1586"; 

# TRAINDIR=Legotrain-SysAdd-Time100ns
# LHC11aData="1218";
# LHC13gData="1220"

# TRAINDIR=Legotrain-SysAdd-Time200ns
# LHC11aData="1219";
# LHC13gData="1221"

# TRAINDIR=Legotrain-SysAdd-ConvCalo2
# LHC11aData="1277";
# LHC15g1aMC="1658";
# LHC13gData="1279"
# LHC15g2MC="1648"; #weighted 40
# LHC15g2MC="1653"; #
# LHC15a3aMC="1655"; 
# LHC15a3aplusMC="1657"; 
# LHC15a3aMC="1654"; 
# LHC15a3aplusMC="1656"; 

# TRAINDIR=Legotrain-SysAdd-Calo2
# LHC13gData="1278"
# # LHC15g2MC="1648"; #weighted 40
# LHC15g2MC="1653"; #
# LHC15a3aMC="1651"; 
# LHC15a3aplusMC="1652"; 

# TRAINDIR=Legotrain-SysAdd-ConvCalo_Calo_FinalJJMC
# LHC15g1aMC="1661";
# 
# LHC15a3aMC="1659"; 
# LHC15a3aplusMC="1660"; 

# TRAINDIR=Legotrain-Rerun-FullCalo-vAN20160512
# LHC11aData="1564";
# LHC12f1aMC="2068"; 
# LHC12f1bMC="2069"; 
# LHC15g1aMC="2071";
# 
# LHC13gData="1568"
# LHC15g2MC="2070";
# LHC15a3aMC="2071"; 
# LHC15a3aplusMC="2072"; 

TRAINDIR=Legotrain-Rerun-FullCalo-vAN20160527_woSort
LHC11aData="1564";
LHC12f1aMC="2125"; 
LHC12f1bMC="2126"; 
LHC15g1aMC="2128";

LHC13gData="1568"
LHC15g2MC="2127";
LHC15a3aMC="2123"; 
LHC15a3aplusMC="2124"; 

OUTPUTDIR=$BASEDIR/$TRAINDIR

echo "************************************************************************************************";
echo "Directory to be writen into " $OUTPUTDIR
echo "************************************************************************************************";
echo ""
echo ""

if [ $2 = "LHC11a" ]; then 
echo "************************************************************************************************";
echo "********************************* Copying LHC11a ***********************************************";
echo "************************************************************************************************";

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
        LHC11aData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC11aData\_`
        OUTPUTDIR_LHC11a=$BASEDIR/$TRAINDIR/GA_pp-$LHC11aData
    fi
    
    if [ $HAVELHC12f1a == 1 ]; then
        LHC12f1aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12f1aMC\_`
        OUTPUTDIR_LHC12f1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1aMC
    fi  
    
    if [ $HAVELHC12f1b == 1 ]; then
        LHC12f1bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12f1bMC\_`
        OUTPUTDIR_LHC12f1b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1bMC
    fi
    
    if [ $HAVELHC12i3 == 1 ]; then
        LHC12i3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12i3MC\_`
        OUTPUTDIR_LHC12i3=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12i3MC
    fi
    if [ $HAVELHC15g1a == 1 ]; then
        LHC15g1aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g1aMC\_`
        OUTPUTDIR_LHC15g1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g1aMC
    fi

    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC11a == 1 ]; then
            echo "downloading LHC11a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC11a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC11aData/merge_runlist_1"
        fi    
        if [ $HAVELHC12f1a == 1 ]; then
            echo "downloading LHC12f1a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1aMC/merge_runlist_1"
        fi    
        if [ $HAVELHC12f1b == 1 ]; then
            echo "downloading LHC12f1b"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1bMC/merge_runlist_1"
        fi    
        if [ $HAVELHC12i3 == 1 ]; then
            echo "LHC12i3"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12i3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12i3MC/merge_runlist_1"
        fi
        if [ $HAVELHC15g1a == 1 ]; then
            echo "LHC15g1a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g1aMC/merge"
            
            echo "copying LHC11aJetJet" 
            runNumbers=`cat runNumbersLHC11aJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat binNumbersJJToMerge.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1a/$binNumber/$runNumber "/alice/sim/2015/LHC15g1a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15g1aMC"
                done;   
            done;
        fi    
    fi

    if [ $HAVELHC11a == 1 ]; then
        ls $OUTPUTDIR_LHC11a/GammaConvCalo_*.root > fileLHC11a.txt
        fileNumbers=`cat fileLHC11a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC11a/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_LHC11a-pass4_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC11a-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_LHC11a_$number.log\"\)
        done;
    fi

    if [ $HAVELHC12f1a == 1 ]; then
        ls $OUTPUTDIR_LHC12f1a/GammaConvCalo_*.root > fileLHC12f1a.txt
        fileNumbers=`cat fileLHC12f1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC12f1a/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC12f1a_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC12f1b == 1 ]; then
        ls $OUTPUTDIR_LHC12f1b/GammaConvCalo_*.root > fileLHC12f1b.txt
        fileNumbers=`cat fileLHC12f1b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC12f1b/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12f1b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC12f1b_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC12f1b_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC12i3 == 1 ]; then
        ls $OUTPUTDIR_LHC12i3/GammaConvCalo_*.root > fileLHC12i3.txt
        fileNumbers=`cat fileLHC12i3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC12i3/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12i3_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC12i3_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC12i3_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC15g1a == 1 ]; then
        if [ $MERGEONBINSSingle = 1 ]; then
            binNumbersJJ=`cat binNumbersJJToMerge.txt`
            echo $binNumbersJJ
            ls $OUTPUTDIR_LHC15g1a/GammaConvCalo_*.root > filetemp.txt
            mkdir -p $OUTPUTDIR/LHC15g1aFineBins/
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number 
                    
                    hadd -f $OUTPUTDIR_LHC15g1a/$binNumber/GammaConvCalo_$number.root $OUTPUTDIR_LHC15g1a/$binNumber/*/GammaConvCalo_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1a/$binNumber/GammaConvCalo_$number.root $OUTPUTDIR/LHC15g1aFineBins/GammaConvCalo_MC_LHC15g1a$binNumber\_$number.root $number
                done;
            done;
        fi
        ls $OUTPUTDIR_LHC15g1a/GammaConvCalo_*.root > fileLHC15g1a.txt
        fileNumbers=`cat fileLHC15g1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1a/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15g1a_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC15g1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC15g1a_$number.log\"\)
        done;
    fi
   
    if [ $MERGEON == 1 ]; then
        rm $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_LHC12f1b_*.root
        ls $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC12f1b_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_LHC12f1b_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12f1b_$number.root
            fi
        done

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC12i3_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_LHC12i3_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12i3_$number.root
            fi
        done

            
        ls $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_LHC12f1b_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_LHC12f1b_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC12i3_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_LHC12f1b_LHC12i3_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12f1a_LHC12f1b_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC12i3_$number.root
            fi
        done
    fi

    if [ $MERGEONBINS == 1 ]; then    
        ls $OUTPUTDIR/GammaConvCalo_MC_LHC15g1a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15g1aFineBins/GammaConvCalo_MC_LHC15g1a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15g1aFineBins/GammaConvCalo_MC_LHC15g1a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15g1aFineBins/GammaConvCalo_MC_LHC15g1a$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15g1aFineBins/GammaConvCalo_MC_LHC15g1a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC15g1aFinerPtHardBins_$number.root $TOMERGE
        done;
    fi
    
elif [ $2 = "LHC13g" ]; then   
    echo "************************************************************************************************";
    echo "********************************* Copying LHC13g ***********************************************";
    echo "************************************************************************************************";

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
        LHC13gData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC13gData\_`
        OUTPUTDIR_LHC13g=$BASEDIR/$TRAINDIR/GA_pp-$LHC13gData
    fi

    if [ $HAVELHC15g2 == 1 ]; then
        LHC15g2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g2MC\_`
        OUTPUTDIR_LHC15g2=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g2MC
    fi
    
    if [ $HAVELHC15a3a == 1 ]; then
        LHC15a3aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3aMC\_`
        OUTPUTDIR_LHC15a3a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aMC
    fi

    if [ $HAVELHC15a3aplus == 1 ]; then
        LHC15a3aplusMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3aplusMC\_`
        OUTPUTDIR_LHC15a3aplus=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aplusMC
    fi
    
    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC13g == 1 ]; then
            echo "downloading LHC13g"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC13gData/merge"
        fi
        
        if [ $HAVELHC15g2 == 1 ]; then
            echo "downloading LHC15g2"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g2MC/merge"
        fi    

        
        if [ $HAVELHC15a3a == 1 ]; then
            echo "downloading LHC15a3a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aMC/merge"
        fi
        
        if [ $HAVELHC15a3aplus == 1 ]; then
            echo "downloading LHC15a3a_plus"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aplusMC/merge"
        fi
        
        echo "copying LHC13gJetJet in bins" 
        runNumbers=`cat runNumbersLHC13gJetJet.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat binNumbersJJToMerge.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC15a3a == 1 ]; then
                        CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aMC"
                    fi
                    if [ $HAVELHC15a3aplus == 1 ]; then
                        CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a_plus/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aplusMC"
                    fi    
                done;   
            fi  
        done;
    fi 

    if [ $HAVELHC13g == 1 ]; then
        ls $OUTPUTDIR_LHC13g/GammaConvCalo_*.root > fileLHC13g.txt
        fileNumbers=`cat fileLHC13g.txt`
        for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13g/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_LHC13g-pass3_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_LHC13g-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_LHC13g_$number.log\"\)
        done;
    fi

    if [ $HAVELHC15g2 == 1 ]; then
        ls $OUTPUTDIR_LHC15g2/GammaConvCalo_*.root > fileLHC15g2.txt
        fileNumbers=`cat fileLHC15g2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g2/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15g2_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC15g2_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC15g2_$number.log\"\)
        done;
    fi

    if [ $HAVELHC15a3a == 1 ]; then
        ls $OUTPUTDIR_LHC15a3a/GammaConvCalo_*.root > fileLHC15a3a.txt
        fileNumbers=`cat fileLHC15a3a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3a/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15a3a_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC15a3a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC15a3a_$number.log\"\)
        done;
    fi    

    if [ $HAVELHC15a3aplus == 1 ]; then
        ls $OUTPUTDIR_LHC15a3aplus/GammaConvCalo_*.root > fileLHC15a3aplus.txt
        fileNumbers=`cat fileLHC15a3aplus.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3aplus/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aplus_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogConvCalo.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC15a3aplus_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC15a3aplus_$number.log\"\)
        done;
    fi
    
    if [ $MERGEONBINSSingle = 1 ]; then
        binNumbersJJ=`cat binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        ls $OUTPUTDIR_LHC15a3a/GammaConvCalo_*.root > filetemp.txt
        mkdir -p $OUTPUTDIR/LHC15a3aXFineBins/
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC15a3a == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC15a3a/$binNumber/GammaConvCalo_$number.root $OUTPUTDIR_LHC15a3a/$binNumber/*/GammaConvCalo_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3a/$binNumber/GammaConvCalo_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3a$binNumber\_$number.root $number
                fi
                if [ $HAVELHC15a3aplus == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaConvCalo_$number.root $OUTPUTDIR_LHC15a3aplus/$binNumber/*/GammaConvCalo_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaConvCalo_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3aplus$binNumber\_$number.root $number
                fi    
            done;
        done;
    fi
    
    rm $OUTPUTDIR/GammaConvCalo_MC_LHC15a3a_LHC15a3aplus_*.root
    rm $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_*.root
    
    ls $OUTPUTDIR/GammaConvCalo_MC_LHC15a3a_*.root > filesForMerging.txt
    if [ $MERGEON = 1 ]; then    
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC15a3a_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aplus_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC15a3a_LHC15a3aplus_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15a3a_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aplus_$number.root
            fi
        done
    fi

    if [ $MERGEONBINS = 1 ]; then    
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3a$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_$number.root $TOMERGE

            TOMERGE="";     
            for bin in $binsForMerging; do
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3aplus$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaConvCalo_MC_LHC15a3aplus$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aplusFinerPtHardBins_$number.root $TOMERGE

            hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC15a3aFinerPtHardBins_$number.root
        done;
    fi
fi
# 
