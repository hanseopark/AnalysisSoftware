#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

NSlashes=10;

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pPb
    NSlashes=8
    NSlashes2=7
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
elif [ $1 = "dmuhlhei" ]; then 
   BASEDIR=~/data/work/Grid
   NSlashes=9;
fi


# TRAINDIR=Legotrain-vAN-20140813-Calo
# LHC13bData="225_20140813-1903"; #ESD
# LHC13cData="226_20140813-1904"; #ESD
# LHC13dData="201_20140729-0931"; #ESD
# LHC13eData="202_20140729-0931"; #ESD
# LHC13fData="203_20140729-0932"; #ESD
# LHC13e7MC="258_20140811-1943";
# LHC13b2_efix_p1MC="263_20140813-1914"; 
# LHC13b2_efix_p2MC="264_20140813-1915"; 
# LHC13b2_efix_p3MC="266_20140813-1916" ; 
# LHC13b2_efix_p4MC="265_20140813-1915"; 

# LHC13bData="227_20140813-1904"; #ESD
# LHC13cData="228_20140813-1904"; #ESD
# # LHC13dData="201_20140729-0931"; #ESD
# # LHC13eData="202_20140729-0931"; #ESD
# # LHC13fData="203_20140729-0932"; #ESD
# # LHC13e7MC="258_20140811-1943";
# LHC13b2_efix_p1MC="267_20140813-1917"; 
# LHC13b2_efix_p2MC="268_20140813-1917"; 
# LHC13b2_efix_p3MC="269_20140813-1918"; 
# LHC13b2_efix_p4MC="270_20140813-1918"; 

# TRAINDIR=Legotrain-vAN-20141012-Calo
# LHC13bData="258_20141012-1342"; #ESD
# LHC13cData="259_20141012-1342"; #ESD
# LHC13dData="201_20140729-0931"; #ESD
# LHC13eData="202_20140729-0931"; #ESD
# LHC13fData="203_20140729-0932"; #ESD
# LHC13e7MC="258_20140811-1943";
# LHC13b2_efix_p1MC="345_20141012-1344"; 
# LHC13b2_efix_p2MC="346_20141012-1344"; 
# LHC13b2_efix_p3MC="347_20141012-1345" ; 
# LHC13b2_efix_p4MC="348_20141012-1345"; 

# TRAINDIR=Legotrain-vAN-20141024-Calo
# LHC13bData="260_20141025-1157"; #ESD
# LHC13cData="261_20141025-1157"; #ESD
# LHC13dData="262_20141025-1243"; #ESD
# LHC13eData="263_20141025-1244"; #ESD
# LHC13fData="264_20141025-1244"; #ESD
# # LHC13e7MC="258_20140811-1943";
# LHC13b2_efix_p1MC="351_20141025-1141"; 
# LHC13b2_efix_p2MC="352_20141025-1146"; 
# LHC13b2_efix_p3MC="353_20141025-1147" ; 
# LHC13b2_efix_p4MC="354_20141025-1147"; 

# TRAINDIR=Legotrain-vAN-20141110-Calo
# LHC13bData="265_20141111-1007"; #ESD
# LHC13cData="266_20141111-0957"; #ESD
# # LHC13dData="262_20141025-1243"; #ESD
# # LHC13eData="263_20141025-1244"; #ESD
# # LHC13fData="264_20141025-1244"; #ESD
# # LHC13e7MC="258_20140811-1943";
# LHC13b2_efix_p1MC="361_20141111-1018"; 
# LHC13b2_efix_p2MC="362_20141111-1019"; 
# LHC13b2_efix_p3MC="363_20141111-1019" ; 
# LHC13b2_efix_p4MC="364_20141111-1022"; 

# TRAINDIR=Legotrain-vAN-20141112-Calo
# LHC13bData="269_20141117-1001"; #ESD
# LHC13cData="268_20141113-1023"; #ESD
# # LHC13dData="262_20141025-1243"; #ESD
# # LHC13eData="263_20141025-1244"; #ESD
# # LHC13fData="264_20141025-1244"; #ESD
# LHC13e7MC="373_20141117-1008";
# LHC13b2_efix_p1MC="365_20141113-1427"; 
# LHC13b2_efix_p2MC="366_20141113-1427"; 
# LHC13b2_efix_p3MC="367_20141113-1428" ; 
# LHC13b2_efix_p4MC="368_20141113-1429"; 
# LHC13e7MC="";
# LHC13b2_efix_p1MC="369_20141117-1003 "; 
# LHC13b2_efix_p2MC=""; 
# LHC13b2_efix_p3MC="371_20141117-1006" ; 
# LHC13b2_efix_p4MC="372_20141117-1007"; 

# TRAINDIR=Legotrain-vAN-20150507-ConvCaloRerunIncTriggerEMCAL
# LHC13bData="314_20150507-1746"; #ESD
# LHC13cData="315_20150507-1747"; #ESD
# LHC13dData="316_20150508-0858"; #ESD
# LHC13eData="317_20150508-0858"; #ESD
# LHC13fData="318_20150508-0859"; #ESD
# LHC13e7MC="";
# LHC13b2_efix_p1MC="464_20150508-1453"; 
# LHC13b2_efix_p2MC="465_20150508-1453"; 
# LHC13b2_efix_p3MC="466_20150508-1457" ; 
# LHC13b2_efix_p4MC="467_20150508-1457"; 

#TRAINDIR=Legotrain-vAN-20150507-ConvCaloRerunIncTriggerPHOS
#LHC13bData="319_20150508-1345"; #ESD
#LHC13cData="320_20150508-1346"; #ESD
#LHC13dData="321_20150508-1346"; #ESD
#LHC13eData="322_20150508-1347"; #ESD
#LHC13fData="323_20150508-1348"; #ESD
#LHC13e7MC="";
#LHC13b2_efix_p1MC="468_20150508-1456";
#LHC13b2_efix_p2MC="469_20150508-1458";
#LHC13b2_efix_p3MC="470_20150508-1459";
#LHC13b2_efix_p4MC="471_20150508-1500";

# TRAINDIR=Legotrain-vAN-20150824-pPb-Calo
# LHC13bData="380_20150824-1156"; #ESD
# LHC13cData="381_20150824-1157"; #ESD
# LHC13dData="382_20150824-1157"; #ESD
# LHC13eData="384_20150824-1158"; #ESD
# #LHC13fData="323_20150508-1348"; #ESD
# #LHC13e7MC="";
# LHC13b2_efix_p1MC="525_20150825-2130";
# LHC13b2_efix_p2MC="529_20150825-2133";
# LHC13b2_efix_p3MC="534_20150825-2135";
# LHC13b2_efix_p4MC="538_20150825-2133";

#LHC13b2_efix_p1MC="526_20150825-2131";
#LHC13b2_efix_p2MC="530_20150825-2134";
#LHC13b2_efix_p3MC="535_20150825-2135";
#LHC13b2_efix_p4MC="539_20150825-2131";

TRAINDIR=Legotrain-vAN-20151010-CaloNewNonLin
LHC13bData="452_20151011-1814"; #ESD
LHC13cData="453_20151011-1815"; #ESD
LHC13dData="456_20151011-1816"; #ESD
LHC13eData="457_20151011-1817"; #ESD
LHC13fData=""; #ESD
LHC13e7MC="";
LHC13b2_efix_p1MC=""; 
LHC13b2_efix_p2MC=""; 
LHC13b2_efix_p3MC="" ; 
LHC13b2_efix_p4MC=""; 


OUTPUTDIR=$BASEDIR/$TRAINDIR
OUTPUTDIR_LHC13b=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13bData
OUTPUTDIR_LHC13c=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13cData
OUTPUTDIR_LHC13d=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13dData
OUTPUTDIR_LHC13e=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13eData
OUTPUTDIR_LHC13f=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13fData
OUTPUTDIR_LHC13b2_efix_p1=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p1MC
OUTPUTDIR_LHC13b2_efix_p2=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p2MC
OUTPUTDIR_LHC13b2_efix_p3=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p3MC
OUTPUTDIR_LHC13b2_efix_p4=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p4MC
OUTPUTDIR_LHC13e7=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13e7MC
mkdir -p $OUTPUTDIR_LHC13b
mkdir -p $OUTPUTDIR_LHC13c
mkdir -p $OUTPUTDIR_LHC13d
mkdir -p $OUTPUTDIR_LHC13e
mkdir -p $OUTPUTDIR_LHC13f
mkdir -p $OUTPUTDIR_LHC13b2_efix_p1
mkdir -p $OUTPUTDIR_LHC13b2_efix_p2
mkdir -p $OUTPUTDIR_LHC13b2_efix_p3
mkdir -p $OUTPUTDIR_LHC13b2_efix_p4
#mkdir -p $OUTPUTDIR_LHC13e7

alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13bData/merge_runlist_1/GammaCalo* file:$OUTPUTDIR_LHC13b/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge_runlist_1/GammaCalo* file:$OUTPUTDIR_LHC13c/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13dData/merge/GammaCalo* file:$OUTPUTDIR_LHC13d/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13eData/merge_runlist_2/GammaCalo* file:$OUTPUTDIR_LHC13e/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13fData/merge/GammaCalo* file:$OUTPUTDIR_LHC13f/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/merge/GammaCalo* file:$OUTPUTDIR_LHC13b2_efix_p1/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/merge/GammaCalo* file:$OUTPUTDIR_LHC13b2_efix_p2/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/merge/GammaCalo* file:$OUTPUTDIR_LHC13b2_efix_p3/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/merge/GammaCalo* file:$OUTPUTDIR_LHC13b2_efix_p4/
#alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge/GammaCalo* file:$OUTPUTDIR_LHC13e7/


   ls $OUTPUTDIR_LHC13b/GammaCalo_*.root > fileLHC13b.txt
   fileNumbers=`cat fileLHC13b.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13b/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC13b-pass3_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_LHC13b-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13b_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC13c/GammaCalo_*.root > fileLHC13c.txt
   fileNumbers=`cat fileLHC13c.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13c/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC13c-pass3_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_LHC13c-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13c_$number.log\"\)
   done;
   
   ls $OUTPUTDIR_LHC13d/GammaCalo_*.root > fileLHC13d.txt
   fileNumbers=`cat fileLHC13d.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13d/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC13d-pass3_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_LHC13d-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13d_$number.log\"\)
   done;
   
   ls $OUTPUTDIR_LHC13e/GammaCalo_*.root > fileLHC13e.txt
   fileNumbers=`cat fileLHC13e.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13e/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC13e-pass3_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_LHC13e-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13e_$number.log\"\)
   done;
   
   ls $OUTPUTDIR_LHC13f/GammaCalo_*.root > fileLHC13f.txt
   fileNumbers=`cat fileLHC13f.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13f/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_LHC13f-pass3_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_LHC13f-pass3_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13f_$number.log\"\)
   done;

   ls $OUTPUTDIR_LHC13b2_efix_p1/GammaCalo_*.root > fileLHC13b2_efix_p1.txt
   fileNumbers=`cat fileLHC13b2_efix_p1.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13b2_efix_p1/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p1_$number.log\"\)
   done;
#    
   ls $OUTPUTDIR_LHC13b2_efix_p2/GammaCalo_*.root > fileLHC13b2_efix_p2.txt
   fileNumbers=`cat fileLHC13b2_efix_p2.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13b2_efix_p2/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p2_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p2_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p2_$number.log\"\)
   done;
# 
   ls $OUTPUTDIR_LHC13b2_efix_p3/GammaCalo_*.root > fileLHC13b2_efix_p3.txt
   fileNumbers=`cat fileLHC13b2_efix_p3.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13b2_efix_p3/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p3_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p3_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p3_$number.log\"\)
   done;

   ls $OUTPUTDIR_LHC13b2_efix_p4/GammaCalo_*.root > fileLHC13b2_efix_p4.txt
   fileNumbers=`cat fileLHC13b2_efix_p4.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13b2_efix_p4/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p4_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p4_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13b2_efix_p4_$number.log\"\)
   done;

   ls $OUTPUTDIR_LHC13e7/GammaCalo_*.root > fileLHC13e7.txt
   fileNumbers=`cat fileLHC13e7.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandardCalo.C\(\"$OUTPUTDIR_LHC13e7/GammaCalo_$number.root\"\,\"$OUTPUTDIR/GammaCalo_MC_LHC13e7_$number.root\"\,\"GammaCalo_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLogCalo.C\(\"$OUTPUTDIR/GammaCalo_MC_LHC13e7_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC13e7_$number.log\"\)
   done;

   
   rm $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_*.root
   ls $OUTPUTDIR/GammaCalo_LHC13b-pass3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
      echo $number
      ls $OUTPUTDIR/GammaCalo_LHC13b-pass3_$number.root
      ls $OUTPUTDIR/GammaCalo_LHC13c-pass3_$number.root
      if [ -f $OUTPUTDIR/GammaCalo_LHC13b-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC13c-pass3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13b-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13c-pass3_$number.root
      fi
   done
   
   ls $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 4 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC13d-pass3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13d-pass3_$number.root
      fi
   done
   
   ls $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 5 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC13e-pass3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_LHC13e-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13b-pass3_LHC13c-pass3_LHC13d-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13e-pass3_$number.root
      fi
   done

   rm $OUTPUTDIR/GammaCalo_LHC13d-pass3_LHC13e-pass3_*.root
    ls $OUTPUTDIR/GammaCalo_LHC13d-pass3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaCalo_LHC13d-pass3_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC13e-pass3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaCalo_LHC13d-pass3_LHC13e-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13d-pass3_$number.root $OUTPUTDIR/GammaCalo_LHC13e-pass3_$number.root
      fi
   done

   
   rm $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_*.root
   ls $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 6 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p2_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p2_$number.root
      fi
   done
#    
   ls $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 7 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p3_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p3_$number.root
      fi
   done
   
   ls $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_*.root > filesForMerging.txt
   filesForMerging=`cat filesForMerging.txt`
   for fileName in $filesForMerging; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 8 | cut -d "." -f1`
      echo $number
      if [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p4_$number.root ] ; then
         hadd -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p4_$number.root
      fi
   done
