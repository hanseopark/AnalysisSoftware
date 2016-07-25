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
NSlashes2=9

# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1
MERGEONBINSSingle=1
MERGEONBINS=1

# check if train configuration has actually been given
HAVELHC11a=1
HAVELHC15g1a=1
HAVELHC15g1b=1
HAVELHC13g=1
HAVELHC15a3a=1
HAVELHC15a3aplus=1
HAVELHC15a3b=1
HAVELHC12f1a=1
HAVELHC12f1b=1
HAVELHC15g2=1

# default trainconfigurations
LHC11aData="";
LHC12f1aMC="";
LHC12f1bMC="";
LHC15g1aMC="";
LHC15g1bMC="";
LHC13gData="";
LHC15a3aMC=""; 
LHC15a3aplusMC=""; 
LHC15a3bMC="";
LHC15g2MC=""; 

if [ $1 = "fbock" ]; then 
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp
    NSlashes=8
    NSlashes2=7
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp
elif [ $1 = "leardini" ]; then 
   BASEDIR=/Users/lucy/
elif [ $1 = "dmuhlheim" ]; then 
   BASEDIR=/home/daniel/Desktop/Grid
   NSlashes=9
   NSlashes=8
fi

# TRAINDIR=Legotrain-vAN-20150909-MergedFirstFullScaleTest
# LHC11aData="985_20150910-2344";
# LHC15g1aMC="1056_20150910-2356";
# LHC15g1bMC="1057_20150910-2354";

# LHC13gData="949_20150904-2205";
# LHC15a3aMC="1053_20150910-2349"; 
# LHC15a3aplusMC="1054_20150910-2347"; 
# LHC15a3bMC="1055_20150910-2351";

# TRAINDIR=Legotrain-vAN-20150924-Merged2ndFullScaleTest
# LHC11aData="1009_20150925-1638";
# LHC15g1aMC="1100_20150925-1647";
# LHC15g1bMC="1101_20150925-1647";
# 
# LHC13gData="1010_20150926-2245";
# LHC15a3aMC="1053_20150910-2349"; 
# LHC15a3aplusMC="1054_20150910-2347"; 
# LHC15a3bMC="1055_20150910-2351";

# TRAINDIR=Legotrain-vAN-20151122-Merged1stAna
# LHC11aData="1107_20151124-0138";
# LHC15g1aMC="1334_20151124-0136";
# LHC15g1bMC="1335_20151124-0138";

# LHC13gData="";
# LHC15a3aMC=""; 
# LHC15a3aplusMC=""; 
# LHC15a3bMC="";

# TRAINDIR=Legotrain-vAN-20151127-Merged1stAna
# LHC11aData="1113_20151129-1913";
# LHC15g1aMC="1344_20151129-1918";
# LHC15g1bMC="1335_20151124-0138";

# LHC13gData="1114_20151129-1913";
# # LHC15a3aMC=""; 
# LHC15a3aplusMC="1342_20151129-1917"; 
# LHC15a3bMC="";

# TRAINDIR=Legotrain-vAN-20151201-MergedGustavoFix
# LHC11aData="1115_20151201-1842";
# LHC15g1aMC="1349_20151201-1836";
# LHC15g1bMC="1350_20151201-1839";
# 
# LHC13gData="1120_20151203-1111";
# LHC15a3aMC="1346_20151201-1841"; 
# LHC15a3aplusMC="1347_20151201-1833"; 
# LHC15a3bMC="1348_20151201-1834";

# TRAINDIR=Legotrain-vAN-20151211-MergedTMenabled
# LHC11aData="1122_20151211-2234";
# LHC15g1aMC="1353_20151211-2237";
# LHC15g1bMC="1354_20151211-2237";
# 
# LHC13gData="1123_20151212-2225";
# LHC15a3aMC="1355_20151212-2226"; 
# LHC15a3aplusMC="1356_20151212-2226"; 
# LHC15a3bMC="1357_20151212-2226";

# TRAINDIR=Legotrain-vAN-20151214-MergedTMenabled-fixNLM
# LHC11aData="1131_20151214-1849";
# LHC15g1aMC="1366_20151214-1851";
# LHC15g1bMC="1367_20151214-1851";
# 
# LHC13gData="1132_20151214-1850";
# LHC15a3aMC="1368_20151214-1852"; 
# LHC15a3aplusMC="1369_20151214-1852"; 
# LHC15a3bMC="1370_20151214-1852";

# TRAINDIR=Legotrain-vAN-20151218-MergedTMenabled-triggerWOPileupRejec
# LHC11aData="1133_20151219-1327";
# LHC15g1aMC="1383_20151219-1308";
# LHC15g1bMC="1384_20151219-1309";
# 
# LHC13gData="1134_20151219-1328";
# LHC15a3aMC="1380_20151219-1305"; 
# LHC15a3aplusMC="1381_20151219-1306"; 
# LHC15a3bMC="1382_20151219-1307";

# TRAINDIR=Legotrain-vAN-20160203-NLMCheckGustavo
# LHC11aData="1222_20160203-0025";
# # LHC15g1aMC="1383_20151219-1308";
# # LHC15g1bMC="1384_20151219-1309";
# 
# LHC13gData="1223_20160203-0026";
# # LHC15a3aMC="1380_20151219-1305"; 
# # LHC15a3aplusMC="1381_20151219-1306"; 
# # LHC15a3bMC="1382_20151219-1307";

# TRAINDIR=Legotrain-Sys-mergedCalo
# LHC11aData="1320";
# LHC15g1aMC="1765";

# LHC11aData="1321";
# LHC15g1aMC="1766";
# 
# LHC11aData="1333";
# LHC15g1aMC="1767";

# LHC11aData="1323";
# LHC15g1aMC="1769";

# LHC13gData="1324";
# LHC15a3aMC="1759"; 
# LHC15a3aplusMC="1762"; 

# LHC13gData="1325";
# LHC15a3aMC="1760"; 
# LHC15a3aplusMC="1763"; 

# LHC13gData="1326";
# LHC15a3aMC="1783"; 
# LHC15a3aplusMC="1768"; 

# LHC13gData="1327";
# LHC15a3aMC="1761"; 
# LHC15a3aplusMC="1787"; 

# LHC13gData="1328";
# LHC15a3aMC="1784"; 
# LHC15a3aplusMC="1788"; 

# LHC15a3aplusMC="1789"; 


# LHC13gData="1329";
# LHC15a3aMC="1785"; 
# LHC15a3aplusMC="1790"; 

# LHC13gData="1330";
# LHC15a3aMC="1786"; 
# LHC15a3aplusMC="1791"; 

# LHC13gData="1367";

# TRAINDIR=Legotrain-SysNew-mergedCalo
# LHC11aData="1369"; 
# LHC15g1aMC="1823";

# LHC11aData="1370";
# LHC15g1aMC="1824";

# LHC11aData="1371";
# LHC15g1aMC="1825";

# LHC11aData="1372";
# LHC15g1aMC="1826";

# LHC11aData="1373";
# LHC15g1aMC="1827";

# LHC13gData="1388";
# LHC15a3aMC="1828"; 
# LHC15a3aplusMC="1837"; 

# LHC13gData="1389";
# LHC15a3aMC="1829"; 
# LHC15a3aplusMC="1838"; 

# LHC13gData="1390";
# LHC15a3aMC="1830"; 
# LHC15a3aplusMC="1839"; 

# LHC13gData="1377";
# LHC15a3aMC="1831"; 
# LHC15a3aplusMC="1840"; 

# LHC13gData="1378";
# LHC15a3aMC="1832"; 
# LHC15a3aplusMC="1841"; 

# LHC13gData="1379";
# LHC15a3aMC="1833"; 
# LHC15a3aplusMC="1842"; 

# LHC13gData="1401";
# LHC15a3aMC="1834"; 
# LHC15a3aplusMC="1843"; 

# LHC13gData="1402";

# LHC13gData="1381";
# LHC15a3aMC="1835"; 
# LHC15a3aplusMC="1844"; 

# LHC13gData="1391";
# LHC15a3aMC="1836"; 
# LHC15a3aplusMC="1845"; 


# TRAINDIR=Legotrain-mergedCalo-FullStat
# LHC11aData="1461"; 
# LHC15g1aMC="1890";
# LHC15g1bMC="1891";
# 
# LHC13gData="1425";
# LHC15a3aMC="1887"; 
# LHC15a3aplusMC="1888"; 
# LHC15a3bMC="1889"; 

# TRAINDIR=Legotrain-Sys-mergedCalo-50ns
# LHC11aData="1462"; 
# LHC13gData="1429";

# TRAINDIR=Legotrain-Sys-mergedCalo-100ns
# LHC11aData="1463"; 
# LHC13gData="1430";
# 
# TRAINDIR=Legotrain-Sys-mergedCalo-200ns
# LHC11aData="1464"; 
# LHC13gData="1431";

# TRAINDIR=Legotrain-mergedCalo-20160408
# LHC11aData="1478"; 
# LHC15g1aMC="1986";
# LHC15g1bMC="1987";

# LHC13gData="1482"; #not optimul statistics
# LHC13gData="1487";
# LHC13gData="1488";
# LHC15a3aMC="1988"; 
# LHC15a3aplusMC="1989"; 
# LHC15a3bMC="1990"; 
# LHC15a3aMC="1991"; 
# LHC15a3aplusMC="1992"; 
# LHC15a3bMC="1993"; 

# TRAINDIR=Legotrain-mergedCalo-20160408_V2Clusterizer
# LHC11aData="1479"; 
# LHC15g1aMC="1994";
# LHC15g1bMC="1995";

# LHC13gData="1481"; 
# LHC15a3aMC="1996"; 
# LHC15a3aplusMC="1998"; 
# LHC15a3bMC="2000"; 
# LHC15a3aMC="1997"; 
# LHC15a3aplusMC="1999"; 
# LHC15a3bMC="2001"; 

# TRAINDIR=Legotrain-mCalo-20160412_V2ClusterizerSys
# LHC11aData="1491"; 
# LHC15g1aMC="2010";
# # LHC11aData="1492"; 
# # LHC15g1aMC="2011";
# 
# 
# LHC13gData="1498"; 
# LHC15a3aMC="2013"; 
# LHC15a3aplusMC="2016"; 
# LHC13gData="1495"; 
# LHC15a3aMC="2014"; 
# LHC15a3aplusMC="2017"; 
# LHC13gData="1496"; 
# LHC15a3aMC="2015"; 
# LHC15a3aplusMC="2018"; 

# TRAINDIR=Legotrain-mCalo-20160412_V2ClusterizerSys100ns
# LHC11aData="1500"; 
# LHC13gData="1503"; 

# TRAINDIR=Legotrain-mCalo-20160412_V2ClusterizerSys200ns
# LHC11aData="1501"; 
# LHC13gData="1504"; 
# 
# TRAINDIR=Legotrain-mCalo-20160412_V2ClusterizerSys50ns
# LHC11aData="1502"; 
# LHC13gData="1505"; 

# TRAINDIR=Legotrain-mCalo-20160412_V1ClusterizerAddInfo
# LHC11aData="1493"; 
# LHC15g1aMC="2012";

# LHC13gData="1499"; 
# LHC15a3aMC="2019"; 
# LHC15a3aplusMC="2020"; 

# TRAINDIR=Legotrain-mCalo-20160501_V2Clusterizer30MeVDiff
# LHC11aData="1534"; 
# LHC15g1aMC="2049";
# 
# LHC13gData="1535"; 
# LHC15a3aMC="2047"; 
# LHC15a3aplusMC="2048"; 

# TRAINDIR=Legotrain-mCalo-20160527_V2ClusNewMCdef
# LHC11aData="1566"; 
# LHC15g1aMC="2131";
# 
# LHC13gData="1567"; 
# LHC15a3aMC="2129"; 
# LHC15a3aplusMC="2130"; 

# TRAINDIR=Legotrain-mCalo-20160527_V2ClusNewMCdef_woSort
# LHC11aData="1566"; 
# LHC15g1aMC="2135";
# 
# LHC13gData="1567"; 
# LHC15a3aMC="2133"; 
# LHC15a3aplusMC="2134"; 

# TRAINDIR=Legotrain-mCalo-20160626-v2ClusterWithTriggMimickInData
# LHC11aData="1637"; 
# LHC13gData="1638"; 


# TRAINDIR=Legotrain-mCalo-20160629_V2Clus
# LHC11aData="1654"; 
# LHC11aData="1655"; 
# LHC11aData="1656"; 
# LHC15g1aMC="2213";
# LHC15g1aMC="2216";
# LHC15g1aMC="2219";
# LHC15g1bMC="2191";

# LHC13gData="1707";
# LHC13gData="1657"; 
# LHC13gData="1658"; 
# LHC13gData="1659"; 
# LHC15a3aMC="2211"; 
# LHC15a3aMC="2214"; 
# LHC15a3aMC="2217"; 
# LHC15a3aplusMC="2212"; 
# LHC15a3aplusMC="2215"; 
# LHC15a3aplusMC="2218"; 
# LHC15a3bMC="2189"; 



# TRAINDIR=Legotrain-mCalo-20160629_V1Clus
# # LHC11aData="1639"; 
# # LHC15g1aMC="2224";
# # 
# LHC13gData="1640"; 
# LHC15a3aMC="2222"; 
# LHC15a3aplusMC="2223"; 

# TRAINDIR=Legotrain-mCalo-20160702_V2ClusNewDefault
# LHC11aData="1674"; 
# LHC15g1aMC="2244";
# 
# LHC13gData="1675"; 
# # LHC15a3aMC="2245"; 
# LHC15a3aMC="2246"; 
# # LHC15a3aplusMC="2247"; 
# LHC15a3aplusMC="2248"; 

# TRAINDIR=Legotrain-mCalo-20160705_V2newSys
# LHC11aData="1677"; 
# LHC15g1aMC="2253";
# LHC15g1aMC="2254";
# LHC12f1aMC="2267"; 
# LHC12f1bMC="2268"; 

# LHC13gData="1705"; 
# LHC13gData="1706"; 
# LHC13gData="1679"; 
# LHC15a3aMC="2255"; 
# LHC15a3aMC="2256"; 
# LHC15a3aMC="2257"; 
# LHC15a3aMC="2258"; 
# LHC15a3aplusMC="2259"; 
# LHC15a3aplusMC="2260"
# LHC15a3aplusMC="2261"
# LHC15a3aplusMC="2262"
# LHC15g2MC="2269";

# TRAINDIR=Legotrain-mCalo-20160714_CorrNL
# LHC11aData="1674"; 
# LHC15g1aMC="2310";
# LHC12f1aMC="2308"; 
# LHC12f1bMC="2309"; 
# 
# LHC13gData="1675"; 
# LHC15a3aMC="2312"; 
# LHC15a3aplusMC="2313"; 
# LHC15g2MC="2311";

TRAINDIR=Legotrain-mCalo-20160720_SecEffiAndTMStudies
LHC11aData="1713"; 
LHC15g1aMC="2320";
LHC12f1aMC="2316"; 
LHC12f1bMC="2317"; 

LHC13gData="1715"; 
LHC15a3aMC="2318"; 
LHC15a3aplusMC="2319"; 
LHC15g2MC="2321";

OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ $2 = "LHC11a" ]; then 
    
    if [ "$LHC11aData" == "" ]; then 
        HAVELHC11a=0;
    fi
    if [ "$LHC15g1aMC" = "" ]; then 
        HAVELHC15g1a=0; 
    fi
    if [ "$LHC15g1bMC" = "" ]; then 
        HAVELHC15g1b=0; 
    fi
    if [ "$LHC12f1aMC" = "" ]; then 
        HAVELHC12f1a=0; 
    fi
    if [ "$LHC12f1bMC" = "" ]; then 
        HAVELHC12f1b=0; 
    fi

    
    if [ $HAVELHC11a == 1 ]; then
        LHC11aData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC11aData\_`
        if [ "$LHC11aData" == "" ]; then 
            HAVELHC11a=0;
        else     
            OUTPUTDIR_LHC11a=$BASEDIR/$TRAINDIR/GA_pp-$LHC11aData
        fi
    fi
    if [ $HAVELHC15g1a == 1 ]; then
        LHC15g1aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g1aMC\_`
        if [ "$LHC15g1aMC" == "" ]; then 
            HAVELHC15g1a=0;
        else     
            OUTPUTDIR_LHC15g1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g1aMC
        fi    
    fi
    if [ $HAVELHC15g1b == 1 ]; then
       LHC15g1bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g1bMC\_`
        if [ "$LHC15g1bMC" == "" ]; then 
            HAVELHC15g1b=0;
        else     
            OUTPUTDIR_LHC15g1b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g1bMC
        fi
    fi
    if [ $HAVELHC12f1a == 1 ]; then
        LHC12f1aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12f1aMC\_`
        if [ "$LHC12f1aMC" == "" ]; then 
            HAVELHC12f1a=0;
        else     
            OUTPUTDIR_LHC12f1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1aMC
        fi
    fi  
    
    if [ $HAVELHC12f1b == 1 ]; then
        LHC12f1bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12f1bMC\_`
        if [ "$LHC12f1bMC" == "" ]; then 
            HAVELHC12f1b=0;
        else     
            OUTPUTDIR_LHC12f1b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1bMC
        fi
    fi


    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC11a == 1 ]; then
            echo "downloading LHC11a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC11a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC11aData/merge_runlist_1"
#             runNumbers=`cat runNumbersLHC11a_pass4_wSDD.txt`
#             echo $runNumbers
#             for runNumber in $runNumbers; do
#                 echo $runNumber
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC11a/$runNumber "/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/GA_pp/$LHC11aData"
#             done;
        fi    
        if [ $HAVELHC12f1a == 1 ]; then
            echo "downloading LHC12f1a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1aMC/merge_runlist_1"
        fi    
        if [ $HAVELHC12f1b == 1 ]; then
            echo "downloading LHC12f1b"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1bMC/merge_runlist_1"
        fi    

        if [ $HAVELHC15g1a == 1 ]; then
            echo "downloading LHC15g1a"
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
        if [ $HAVELHC15g1b == 1 ]; then
            echo "downloading LHC15g1b"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g1bMC/merge"
            echo "copying LHC11aJetJet - gamma" 
            runNumbers=`cat runNumbersLHC11aJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJGamma=`cat binNumbersJJGamma.txt`
                for binNumber in $binNumbersJJGamma; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1b/$binNumber/$runNumber "/alice/sim/2015/LHC15g1b/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15g1bMC"
                done;   
            done;
        fi    
        
    fi
    
    if [ $HAVELHC11a == 1 ]; then
        ls $OUTPUTDIR_LHC11a/GammaCaloMerged_*.root > fileLHC11a.txt
        fileNumbers=`cat fileLHC11a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC11a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_LHC11a-pass4_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_LHC11a-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC11a_$number.log\"\)
        done;
    fi

    if [ $HAVELHC12f1a == 1 ]; then
        ls $OUTPUTDIR_LHC12f1a/GammaCaloMerged_*.root > fileLHC12f1a.txt
        fileNumbers=`cat fileLHC12f1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC12f1a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC12f1a_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC12f1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCaloMerged_MC_LHC12f1a_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC12f1b == 1 ]; then
        ls $OUTPUTDIR_LHC12f1b/GammaCaloMerged_*.root > fileLHC12f1b.txt
        fileNumbers=`cat fileLHC12f1b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC12f1b/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC12f1b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC12f1b_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCaloMerged_MC_LHC12f1b_$number.log\"\)
        done;
    fi
    
    
    if [ $HAVELHC15g1a == 1 ]; then  
        ls $OUTPUTDIR_LHC15g1a/GammaCaloMerged_*.root > fileLHC15g1a.txt
        fileNumbers=`cat fileLHC15g1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1a_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15g1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15g1a_$number.log\"\)
        done;
        
        binNumbersJJ=`cat binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        if [ $MERGEONBINSSingle = 1 ]; then
            ls $OUTPUTDIR_LHC15g1a/GammaCaloMerged_*.root > filetemp.txt
            mkdir $OUTPUTDIR/LHC15g1aFineBins
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    hadd -f $OUTPUTDIR_LHC15g1a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15g1a/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$binNumber\_$number.root $number
                done;
            done;
        fi
    fi
    
    if [ $HAVELHC15g1b == 1 ]; then
        ls $OUTPUTDIR_LHC15g1b/GammaCaloMerged_*.root > fileLHC15g1b.txt
        fileNumbers=`cat fileLHC15g1b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1b/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15g1b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15g1b_$number.log\"\)
        done;
        if [ $MERGEONBINSSingle = 1 ]; then
            binNumbersJJGamma=`cat binNumbersJJGamma.txt`
            ls $OUTPUTDIR_LHC15g1b/GammaCaloMerged_*.root > filetemp.txt
            mkdir $OUTPUTDIR/LHC15g1aFineBins
            for binNumber in $binNumbersJJGamma; do
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    hadd -f $OUTPUTDIR_LHC15g1b/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15g1b/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1b/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1b$binNumber\_$number.root $number
                done;
            done;
        fi
    fi
    
    if [ $MERGEONBINS == 1 ]; then
        rm $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1aFinerPtHardBins_*.root
        ls $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1aFinerPtHardBins_$number.root $TOMERGE
        done;
    fi
    
elif [ $2 = "LHC13g" ]; then 

    if [ "$LHC13gData" == "" ]; then 
        HAVELHC13g=0;
    fi
    if [ "$LHC15a3aMC" = "" ]; then 
        HAVELHC15a3a=0; 
    fi
    if [ "$LHC15a3aplusMC" = "" ]; then 
        HAVELHC15a3aplus=0; 
    fi
    if [ "$LHC15a3bMC" = "" ]; then 
        HAVELHC15a3b=0; 
    fi
    if [ "$LHC15g2MC" == "" ]; then 
        HAVELHC15g2=0;
    fi

    
    if [ $HAVELHC13g == 1 ]; then
        LHC13gData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC13gData\_`
        if [ "$LHC13gData" == "" ]; then 
            HAVELHC13g=0;
        else     
            OUTPUTDIR_LHC13g=$BASEDIR/$TRAINDIR/GA_pp-$LHC13gData
        fi
    fi
    
    if [ $HAVELHC15g2 == 1 ]; then
        LHC15g2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g2MC\_`
        if [ "$LHC15g2MC" == "" ]; then 
            HAVELHC15g2=0;
        else     
            OUTPUTDIR_LHC15g2=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g2MC
        fi
    fi
    
    if [ $HAVELHC15a3a == 1 ]; then
        LHC15a3aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3aMC\_`
        if [ "$LHC15a3aMC" == "" ]; then 
            HAVELHC15a3a=0;
        else     
            OUTPUTDIR_LHC15a3a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aMC
        fi    
    fi
    
    if [ $HAVELHC15a3aplus == 1 ]; then
        LHC15a3aplusMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3aplusMC\_`
        if [ "$LHC15a3aplusMC" == "" ]; then 
            HAVELHC15a3aplus=0;
        else     
            OUTPUTDIR_LHC15a3aplus=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aplusMC
        fi    
    fi
    
    if [ $HAVELHC15a3b == 1 ]; then
        LHC15a3bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3bMC\_`
        if [ "$LHC15a3bMC" == "" ]; then 
            HAVELHC15a3b=0;
        else     
            OUTPUTDIR_LHC15a3b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3bMC
        fi    
    fi
    
    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC13g == 1 ]; then
            echo "copying LHC13g" 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC13gData/merge"
            runNumbers=`cat runNumbersLHC13g_pass1.txt`
            echo $runNumbers
#             for runNumber in $runNumbers; do
#                 echo $runNumber    
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC13g/$runNumber "/alice/data/2013/LHC13g/000$runNumber/pass1/PWGGA/GA_pp/$LHC13gData"
#             done;            
        fi    
        
        if [ $HAVELHC15g2 == 1 ]; then
            echo "downloading LHC15g2"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g2MC/merge"
        fi

        if [ $HAVELHC15a3a == 1 ]; then
            echo "copying LHC15a3a" 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aMC/merge"

            echo "copying LHC15a3a single bins" 
            runNumbers=`cat runNumbersLHC13gJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat binNumbersJJToMerge.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aMC"
                done;   
            done;

        fi

        if [ $HAVELHC15a3aplus == 1 ]; then
            echo "copying LHC15a3a_plus" 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aplusMC/merge"
            echo "copying LHC15a3a_plus single bins" 
            runNumbers=`cat runNumbersLHC13gJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat binNumbersJJToMerge.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a_plus/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aplusMC"
                done;   
            done;
        fi

        if [ $HAVELHC15a3b == 1 ]; then
            echo "copying LHC15a3b" 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3bMC/merge"
            echo "copying LHC15a3b single bins" 
            runNumbers=`cat runNumbersLHC13gJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJGamma=`cat binNumbersJJGamma.txt`
                for binNumber in $binNumbersJJGamma; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3b/$binNumber/$runNumber "/alice/sim/2015/LHC15a3b/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3bMC"
                done;   
            done;
        fi
    fi

    if [ $HAVELHC13g == 1 ]; then
        ls $OUTPUTDIR_LHC13g/GammaCaloMerged_*.root > fileLHC13g.txt
        fileNumbers=`cat fileLHC13g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC13g/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_LHC13g-pass1_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_LHC13g-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13g_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC15g2 == 1 ]; then
        ls $OUTPUTDIR_LHC15g2/GammaCaloMerged_*.root > fileLHC15g2.txt
        fileNumbers=`cat fileLHC15g2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g2/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15g2_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15g2_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCaloMerged_MC_LHC15g2_$number.log\"\)
        done;
    fi

    
    if [ $HAVELHC15a3a == 1 ]; then
        binNumbersJJ=`cat binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        ls $OUTPUTDIR_LHC15a3a/GammaCaloMerged_*.root > filetemp.txt
        mkdir -p $OUTPUTDIR/LHC15a3aXFineBins
        if [ $MERGEONBINSSingle = 1 ]; then
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number

                    hadd -f $OUTPUTDIR_LHC15a3a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15a3a/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$binNumber\_$number.root $number
                done;
            done;
        fi
    
        ls $OUTPUTDIR_LHC15a3a/GammaCaloMerged_*.root > fileLHC15a3a.txt
        fileNumbers=`cat fileLHC15a3a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root  $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15a3a_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC15a3aplus == 1 ]; then
        binNumbersJJ=`cat binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        mkdir -p $OUTPUTDIR/LHC15a3aXFineBins
        ls $OUTPUTDIR_LHC15a3aplus/GammaCaloMerged_*.root > filetemp.txt
        if [ $MERGEONBINSSingle = 1 ]; then
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    hadd -f $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15a3aplus/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$binNumber\_$number.root $number
                done;
            done;
        fi

        ls $OUTPUTDIR_LHC15a3aplus/GammaCaloMerged_*.root > fileLHC15a3aplus.txt
        fileNumbers=`cat fileLHC15a3aplus.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3aplus/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15a3aplus_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC15a3b == 1 ]; then
        ls $OUTPUTDIR_LHC15a3b/GammaCaloMerged_*.root > fileLHC15a3b.txt
        fileNumbers=`cat fileLHC15a3b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3b/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLogCaloMerged.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15a3b_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC15a3b_$number.log\"\)
        done;
    fi
    
    if [ $MERGEON = 1 ]; then 
        rm $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_LHC15a3aplus_*.root
        ls $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root ] && [ -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_LHC15a3aplus_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root
            fi
        done
    fi
    
    if [ $MERGEONBINS = 1 ]; then 
        rm $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_*.root
        ls $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_$number.root $TOMERGE

            TOMERGE="";		
            for bin in $binsForMerging; do
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$bin"
                    TOMERGE+="_$number.root"
                else 
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplusFinerPtHardBins_$number.root $TOMERGE

            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_$number.root
        done;
    fi
fi