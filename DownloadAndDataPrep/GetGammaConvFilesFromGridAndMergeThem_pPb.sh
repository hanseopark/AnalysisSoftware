#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

if [ $1 = "fbock" ]; then 
   BASEDIR=/home/fbock/Photon/Grid/OutputLegoTrains/pPb
elif [ $1 = "fbockGSI" ]; then 
   BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pPb
elif [ $1 = "leardini" ]; then 
   BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then 
   BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
elif [ $1 = "leardiniGSI" ]; then 
   BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/pPb
elif [ $1 = "passfeld" ]; then 
   BASEDIR=~/work/Gridoutput/pPb
elif [ $1 = "passfeldMAF" ]; then 
   BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then  
   BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pPb
elif [ $1 = "amarin" ]; then     
   BASEDIR=/Users/marin/
elif [ $1 = "amarinGSI" ]; then     
   BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/pPb 
elif [ $1 = "amarinALICESERV1" ]; then     
   BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/   
elif [ $1 = "mwilde" ]; then        
   BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
   BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pPb 
elif [ $1 = "pgonzales" ]; then     
   BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
   BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pPb
fi

# TRAINDIR=Legotrain-v5-04-68-AN
# LHC13bData=49_20130618-1451; 
# LHC13cData=46_20130618-1134; 
# LHC13b2_efix_p1MC=27_20130618-1830; 

# TRAINDIR=Legotrain-v5-05-20-AN
# LHC13bData=61_20130924-1346; 
# LHC13cData=62_20130924-1350;
# # LHC13dData=63_20130924-1413; 
# # LHC13eData=64_20130924-1351; 
# # LHC13fData=65_20130924-1351; 
# LHC13b2_efix_p1MC=41_20130924-2114; 
# LHC13b2_efix_p2MC=42_20130924-2114; 
# LHC13b2_efix_p3MC=43_20130924-2114; 
# LHC13b2_efix_p4MC=44_20130924-2114 ; 

# TRAINDIR=Legotrain-v5-05-43-AN
# LHC13bData=10_20131203-2000; #AOD
# LHC13cData=11_20131203-2000; #AOD
# LHC13bData=87_20131203-1422; #ESD
# LHC13cData=88_20131203-1427; #ESD
# 
# LHC13b2_efix_p1MC=60_20131203-1936; 
# LHC13b2_efix_p2MC=61_20131203-1936; 
# LHC13b2_efix_p3MC=62_20131203-1937; 
# LHC13b2_efix_p4MC=63_20131203-1937; 
# LHC13e7MC=64_20131203-1938; 
# LHC13e7MC=65_20131203-1938; 

# TRAINDIR=Legotrain-v5-05-56-AN-20140112
# LHC13bData=92_20140112-1928; #ESD
# LHC13cData=97_20140113-1043; #ESD
# LHC13bData=94_20140112-1936; #ESD
# LHC13cData=98_20140113-1043; #ESD
# LHC13bData=96_20140112-1936; #ESD
# LHC13cData=99_20140113-1043; #ESD

# LHC13b2_efix_p1MC=60_20131203-1936; 
# LHC13b2_efix_p2MC=61_20131203-1936; 
# LHC13b2_efix_p3MC=62_20131203-1937; 
# LHC13b2_efix_p4MC=63_20131203-1937; 
# LHC13e7MC="77_20140117-2035" ;
# LHC13e7MC=65_20131203-1938; 

# TRAINDIR=Legotrain-v5-05-60-AN-20140123b
# LHC13bData="101_20140123-1054"; #ESD
# LHC13cData="102_20140123-1612"; #ESD
# LHC13bData="103_20140123-1055"; #ESD
# LHC13cData="104_20140123-1615"; #ESD

# LHC13e7MC="83_20140123-2017" ;
# LHC13e7MC="84_20140123-2018";
# LHC13e7MC="85_20140123-2018";
# LHC13e7MC="86_20140123-2019";
# LHC13e7MC="99_20140130-1430";
# LHC13e7MC="100_20140130-1430";


# TRAINDIR=Legotrain-v5-05-63-AN-SysErr
# LHC13bData="101_20140123-1054"; #ESD
# LHC13cData="102_20140123-1612"; #ESD
# LHC13e7MC="107_20140203-0845";
# LHC13e7MC="108_20140203-0846";
# LHC13e7MC="110_20140203-0848";
# LHC13e7MC="111_20140203-0849";

# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErr
# LHC13bData="140_20140213-1553"; #ESD
# LHC13bData="142_20140213-1558";
# LHC13bData="144_20140214-0723";
# LHC13bData="146_20140214-0723";
# LHC13cData="141_20140213-1554"; #ESD
# LHC13cData="143_20140213-2054"; #ESD
# LHC13cData="145_20140213-2054"; #ESD
# LHC13cData="147_20140213-2054"; #ESD
# LHC13b2_efix1=""; 
# LHC13b2_efix2=""; 
# LHC13b2_efix3=""; 
# LHC13b2_efix4=""; 
# LHC13e7MC="133_20140212-1256";
# LHC13e7MC="134_20140213-0840";
# LHC13e7MC="135_20140213-0840";
# LHC13e7MC="136_20140214-1718";
# LHC13e7MC="137_20140214-1718";
# LHC13e7MC="138_20140214-1719";
# LHC13e7MC="139_20140214-1719";
# LHC13e7MC="140_20140214-1720";  
# LHC13e7MC="141_20140214-1720";  
# LHC13e7MC="142_20140214-1735";  

# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErrDiffEta
# LHC13bData="153_20140214-1738"; #ESD
# LHC13cData="154_20140214-1739"; #ESD
# LHC13e7MC="146_20140214-1742";
# LHC13e7MC="157_20140228-1015";

# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErrEtaMB
# LHC13bData="153_20140214-1738"; #ESD
# LHC13cData="154_20140214-1739"; #ESD
# LHC13e7MC="146_20140214-1742";
# LHC13e7MC="102_20140203-0801";
# LHC13e7MC="157_20140228-1015";
# LHC13e7MC="104_20140203-0804";
# LHC13e7MC="105_20140203-0805";
# LHC13e7MC="106_20140203-0809";
# LHC13e7MC="109_20140203-0847";
# LHC13e7MC="119_20140206-1627";

# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErr_7
# LHC13e7MC="175_20140401-1137";
# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErr_8
# LHC13e7MC="181_20140402-1242";
# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErr_9
# LHC13e7MC="182_20140402-1242";

# TRAINDIR=Legotrain-vAN-20140408
# LHC13bData="169_20140408-2125"; #ESD
# LHC13cData="170_20140408-2126"; #ESD
# LHC13e7MC="183_20140408-1806";

# TRAINDIR=Legotrain-vAN-20140410
# LHC13bData="171_20140411-0008"; #ESD
# LHC13cData="172_20140411-0007"; #ESD
# LHC13e7MC="186_20140411-0008";

# TRAINDIR=Legotrain-vAN-20140420
# LHC13bData="175_20140419-1616"; #ESD
# LHC13cData="174_20140417-1958"; #ESD
# LHC13e7MC="190_20140417-2006";
# LHC13e7MC="191_20140417-2017";

# TRAINDIR=Legotrain-vAN-20140811-ConvV1
# LHC13bData="221_20140811-2021"; #ESD
# LHC13cData="222_20140811-2022"; #ESD
# # LHC13dData="201_20140729-0931"; #ESD
# # LHC13eData="202_20140729-0931"; #ESD
# # LHC13fData="203_20140729-0932"; #ESD
# LHC13e7MC="258_20140811-1943";
# LHC13b2_efix_p1MC="254_20140811-1942"; 
# LHC13b2_efix_p2MC="255_20140811-1942"; 
# LHC13b2_efix_p3MC="256_20140811-1942"; 
# LHC13b2_efix_p4MC="257_20140811-1943"; 

# TRAINDIR=Legotrain-vAN-20141124-ConvV1
# LHC13bData="270_20141125-1028"; #ESD
# LHC13cData="271_20141125-1029"; #ESD
# # LHC13bData="272_20141125-1050"; #ESD
# # LHC13cData="273_20141125-1050"; #ESD
# LHC13e7MC="380_20141125-1041";
# LHC13b2_efix_p1MC="376_20141125-1039"; 
# LHC13b2_efix_p2MC=""; 
# LHC13b2_efix_p3MC="378_20141125-1040"; 
# LHC13b2_efix_p4MC="379_20141125-1040"; 
# # LHC13e7MC="385_20141125-1043";
# # LHC13b2_efix_p1MC="381_20141125-1042"; 
# # LHC13b2_efix_p2MC="382_20141125-1042"; 
# # LHC13b2_efix_p3MC="383_20141125-1046"; 
# # LHC13b2_efix_p4MC="384_20141125-1043"; 

TRAINDIR=Legotrain-vAN-20150113-newPrimDefTest
   LHC13bData="298_20150115-1317"; #ESD
   LHC13cData="299_20150115-1317"; #ESD
#    LHC13e7MC="420_20150115-1307" ; #ESD 
#    LHC13b2_efix_p1MC="416_20150115-1305" ; #ESD
#    LHC13b2_efix_p2MC="417_20150115-1305" ; #ESDs
#    LHC13b2_efix_p3MC="418_20150115-1307" ; #ESD
#    LHC13b2_efix_p4MC="419_20150115-1307" ; #ESD
   LHC13e7MC="425_20150119-1718"; #ESD 
   LHC13b2_efix_p1MC="421_20150119-1710"; #ESD
   LHC13b2_efix_p2MC="422_20150119-1718"; #ESDs
   LHC13b2_efix_p3MC="423_20150119-1718"; #ESD
   LHC13b2_efix_p4MC="424_20150119-1718"; #ESD


OUTPUTDIR=$BASEDIR/$TRAINDIR
OUTPUTDIR_LHC13b=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13bData
OUTPUTDIR_LHC13c=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13cData
# OUTPUTDIR_LHC13d=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13dData
# OUTPUTDIR_LHC13e=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13eData
# OUTPUTDIR_LHC13f=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13fData
OUTPUTDIR_LHC13b2_efix_p1=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p1MC
OUTPUTDIR_LHC13b2_efix_p2=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p2MC
OUTPUTDIR_LHC13b2_efix_p3=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p3MC
OUTPUTDIR_LHC13b2_efix_p4=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p4MC
OUTPUTDIR_LHC13e7=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13e7MC
mkdir -p $OUTPUTDIR_LHC13b
mkdir -p $OUTPUTDIR_LHC13c
# mkdir -p $OUTPUTDIR_LHC13d
# mkdir -p $OUTPUTDIR_LHC13e
# mkdir -p $OUTPUTDIR_LHC13f
mkdir -p $OUTPUTDIR_LHC13b2_efix_p1
mkdir -p $OUTPUTDIR_LHC13b2_efix_p2
mkdir -p $OUTPUTDIR_LHC13b2_efix_p3
mkdir -p $OUTPUTDIR_LHC13b2_efix_p4
mkdir -p $OUTPUTDIR_LHC13e7

alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13bData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC13b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC13c/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13dData/merge/GammaConvV1* file:$OUTPUTDIR_LHC13d/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13eData/merge/GammaConvV1* file:$OUTPUTDIR_LHC13e/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13fData/merge/GammaConvV1* file:$OUTPUTDIR_LHC13f/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/merge/GammaConvV1* file:$OUTPUTDIR_LHC13b2_efix_p1/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/merge/GammaConvV1* file:$OUTPUTDIR_LHC13b2_efix_p2/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/merge/GammaConvV1* file:$OUTPUTDIR_LHC13b2_efix_p3/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/merge/GammaConvV1* file:$OUTPUTDIR_LHC13b2_efix_p4/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge/GammaConvV1* file:$OUTPUTDIR_LHC13e7/


   ls $OUTPUTDIR_LHC13b/GammaConvV1_*.root > fileLHC13b.txt
   fileNumbers=`cat fileLHC13b.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13b-pass3_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13b-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13b_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC13c/GammaConvV1_*.root > fileLHC13c.txt
   fileNumbers=`cat fileLHC13c.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13c/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13c-pass3_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13c-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13c_$number.log\"\)
   done;
#    
#    ls $OUTPUTDIR_LHC13d/GammaConvV1_*.root > fileLHC13d.txt
#    fileNumbers=`cat fileLHC13d.txt`
#    for fileName in $fileNumbers; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13d/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13d-pass3_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13d-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13d_$number.log\"\)
#    done;
#    
#    ls $OUTPUTDIR_LHC13e/GammaConvV1_*.root > fileLHC13e.txt
#    fileNumbers=`cat fileLHC13e.txt`
#    for fileName in $fileNumbers; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13e/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13e-pass3_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13e-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13e_$number.log\"\)
#    done;
#    
#    ls $OUTPUTDIR_LHC13f/GammaConvV1_*.root > fileLHC13f.txt
#    fileNumbers=`cat fileLHC13f.txt`
#    for fileName in $fileNumbers; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13f/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC13f-pass3_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13f-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13f_$number.log\"\)
#    done;
# 
   ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_*.root > fileLHC13b2_efix_p1.txt
   fileNumbers=`cat fileLHC13b2_efix_p1.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p1_$number.log\"\)
   done;
   
   ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_*.root > fileLHC13b2_efix_p2.txt
   fileNumbers=`cat fileLHC13b2_efix_p2.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p2_$number.log\"\)
   done;

   ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_*.root > fileLHC13b2_efix_p3.txt
   fileNumbers=`cat fileLHC13b2_efix_p3.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p3_$number.log\"\)
   done;

   ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_*.root > fileLHC13b2_efix_p4.txt
   fileNumbers=`cat fileLHC13b2_efix_p4.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p4_$number.log\"\)
   done;

   ls $OUTPUTDIR_LHC13e7/GammaConvV1_*.root > fileLHC13e7.txt
   fileNumbers=`cat fileLHC13e7.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC13e7/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC13e7_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13e7_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13e7_$number.log\"\)
   done;

   
   rm $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_*.root
   ls $OUTPUTDIR/GammaConvV1_LHC13b-pass3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
      echo $number
      ls $OUTPUTDIR/GammaConvV1_LHC13b-pass3_$number.root
      ls $OUTPUTDIR/GammaConvV1_LHC13c-pass3_$number.root
      if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13c-pass3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13c-pass3_$number.root
      fi
   done
   
   ls $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 4 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13d-pass3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13d-pass3_$number.root
      fi
   done
#    
#    ls $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_*.root > filesForMerging.txt
#    filesForMerging=`cat filesForMerging.txt`
#    for fileName in $filesForMerging; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 5 | cut -d "." -f1`
#       echo $number
#       if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13e-pass3_$number.root ] ; then
#          hadd -f $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_LHC13e-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13e-pass3_$number.root
#       fi
#    done
# 
#    rm $OUTPUTDIR/GammaConvV1_LHC13d-pass3_LHC13e-pass3_*.root
#     ls $OUTPUTDIR/GammaConvV1_LHC13d-pass3_*.root > filesForMerging.txt
#    filesForMerging=`cat filesForMerging.txt`
#    for fileName in $filesForMerging; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
#       echo $number
#       if [ -f $OUTPUTDIR/GammaConvV1_LHC13d-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13e-pass3_$number.root ] ; then
#          hadd -f $OUTPUTDIR/GammaConvV1_LHC13d-pass3_LHC13e-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13d-pass3_$number.root $OUTPUTDIR/GammaConvV1_LHC13e-pass3_$number.root
#       fi
#    done
# 
#    
   rm $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_*.root
   ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 6 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root
      fi
   done
#    
   ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 7 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root
      fi
   done
   
   ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 8 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root
      fi
   done
#    