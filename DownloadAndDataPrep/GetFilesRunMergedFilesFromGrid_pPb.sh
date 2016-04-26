#! /bin/bash

# variable: $1 = username, $2 = AOD/ESD, $3 = filename

# usernames: fbock, leardini, passfeld, amarin, mwilde, pgonzales
# specific location is added the same username will have a different output directory
# i.e fbockGSI or leardiniALICESERV1, passfeldMAF

if [ $1 = "fbock" ]; then 
   BASEDIR=/home/fbock/Photon/Grid/OutputLegoTrains/pPb/
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

mkdir -p $BASEDIR
fileName=$3 
 
# TRAINDIR=Legotrain-v5-05-31-AN-20131105
#    LHC13b="78_20131105-1504"
#    LHC13c="79_20131105-1510"
#    LHC13e7="55_20131106-1005" ; #ESD

# TRAINDIR=Legotrain-v5-05-35-AN-20131110
#    LHC13b="82_20131111-0951"; #
#    LHC13c="81_20131111-0951"; #
#    LHC13e7="57_20131111-1117" ; #ESD #56_20131110-1001
#    LHC13b2_efix1="51_20131105-1615" ; #ESD
#    LHC13b2_efix2="52_20131105-1616" ; #ESDs
#    LHC13b2_efix3="53_20131105-1616" ; #ESD
#    LHC13b2_efix4="54_20131105-1626" ; #ESD

# TRAINDIR=Legotrain-v5-05-56-AN-20140112
#    LHC13b="97_20140113-1043"; #
#    LHC13c="92_20140112-1928"; #
#    LHC13b="98_20140113-1043"; #
#    LHC13c="94_20140112-1936"; #
#    LHC13b="99_20140113-1043"; #
#    LHC13c="96_20140112-1936"; # 
#    LHC13e7="77_20140117-2035" ;
#    LHC13b2_efix1="51_20131105-1615" ; #ESD
#    LHC13b2_efix2="52_20131105-1616" ; #ESDs
#    LHC13b2_efix3="53_20131105-1616" ; #ESD
#    LHC13b2_efix4="54_20131105-1626" ; #ESD
   
# TRAINDIR=Legotrain-v5-05-49-AN-20121221
# LHC13b="89_20131223-1109"; #ESD
# LHC13c="90_20131223-1108"; #ESD
# 
# LHC13b2_efix1="67_20131218-1324"; 
# LHC13b2_efix2="69_20131218-1753"; 
# LHC13b2_efix3="70_20131218-1753"; 
# LHC13b2_efix4="71_20131218-1753"; 
# LHC13e7="66_20131218-1324";

# TRAINDIR=Legotrain-v5-05-63-AN-20140131-QA
# LHC13b="109_20140131-1052"; #ESD
# LHC13c="110_20140131-1055"; #ESD
# 
# LHC13b2_efix1=""; 
# LHC13b2_efix2=""; 
# LHC13b2_efix3=""; 
# LHC13b2_efix4=""; 
# LHC13e7="101_20140131-1139";

# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErr

# LHC13b="140_20140213-1553"; #ESD
# LHC13b="142_20140213-1558";
# LHC13b="144_20140214-0723";
# LHC13b="146_20140214-0723";
# LHC13c="141_20140213-1554"; #ESD
# LHC13c="143_20140213-2054"; #ESD
# LHC13c="145_20140213-2054"; #ESD
# LHC13c="147_20140213-2054"; #ESD
LHC13b2_efix1=""; 
LHC13b2_efix2=""; 
LHC13b2_efix3=""; 
LHC13b2_efix4=""; 
# LHC13e7="133_20140212-1256";
# LHC13e7="134_20140213-0840";
# LHC13e7="135_20140213-0840";
# LHC13e7="136_20140214-1718";
# LHC13e7="137_20140214-1718";
# LHC13e7="138_20140214-1719";
# LHC13e7="139_20140214-1719";
# LHC13e7="140_20140214-1720";  
# LHC13e7="141_20140214-1720";  
# LHC13e7="142_20140214-1735";  
 
# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErr0020 
# LHC13b="119_20140131-0858"; 
# LHC13b="113_20140130-1611"; 
# LHC13b="115_20140130-1612"; 
# LHC13b="117_20140130-1613"; 
# LHC13b="132_20140131-1241"; 
# LHC13c="112_20140130-1610"; 
# LHC13c="114_20140130-1611"; 
# LHC13c="116_20140130-1612"; 
# LHC13c="118_20140130-1614"; 
# LHC13c="133_20140131-1241"; 
# LHC13e7="138_20140214-1719";
# LHC13e7="139_20140214-1719";
# LHC13e7="140_20140214-1720";  
# LHC13e7="141_20140214-1720";  
# LHC13e7="142_20140214-1735";
# LHC13e7="143_20140214-1736";


# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErrDiffEta 
# LHC13b=""; 
# LHC13c=""; 
# LHC13e7="158_20140228-1015";
# LHC13e7="157_20140228-1015";

# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErrEtaMB 
# LHC13e7="158_20140228-1015";
# LHC13e7="157_20140228-1015";
# LHC13b="120_20140131-1219"; #133,135,137
# LHC13b="122_20140131-1220"; #139,141,143
# LHC13b="124_20140131-1222"; #145,147,149
# LHC13b="126_20140131-1224"; #151,153
# LHC13c="121_20140131-1219"; #133,135,137
# LHC13c="123_20140131-1221"; #139,141,143
# LHC13c="125_20140131-1223"; #145,147,149
# LHC13c="127_20140131-1224"; #151,153
# LHC13e7="102_20140203-0801"; #133,134,135,136
# LHC13e7="109_20140203-0847"; #137,138,139,140  
# LHC13e7="104_20140203-0804"; #141,142,143,144
# LHC13e7="105_20140203-0805"; #145,146,147,148
# LHC13e7="106_20140203-0809"; #149,150,151,152
# LHC13e7="119_20140206-1627"; #153,154

# TRAINDIR=Legotrain-v5-05-63-AN-SystematicErrEta0020 
# LHC13b="126_20140131-1224"; #155
# LHC13b="128_20140131-1225"; #157,159,161
# LHC13b="130_20140131-1226"; #163,165,167
# LHC13b="132_20140131-1241"; #169,171
# LHC13c="127_20140131-1224"; #155
# LHC13c="129_20140131-1225"; #157,159,161
# LHC13c="131_20140131-1227"; #163,165,167
# LHC13c="133_20140131-1241"; #169,171
# LHC13e7="119_20140206-1627"; #155,156

# LHC13e7="120_20140203-0921"; #157,158,159,160 
# LHC13e7="126_20140207-1340"; #161,162,163,164
# LHC13e7="127_20140207-1340"; #165,166,167,168 
# LHC13e7="128_20140210-2146"; #169,170,171,172 

# TRAINDIR=MaterialBudget
# LHC13b="178_20140612-1941"; 
# LHC13c="179_20140612-1942"; 
# LHC13e7="198_20140612-1155";

# TRAINDIR=MaterialBudgetWithoutRCut
# LHC13b="180_20140618-1057"; 
# LHC13c="181_20140618-1058"; 
# LHC13e7="200_20140618-1112";

# TRAINDIR=MaterialBudgetWithNTrackinTree
# LHC13b="182_20140620-0906"; 
# LHC13c="183_20140620-0906"; 
# LHC13e7="201_20140620-0907";

TRAINDIR=Legotrain-vAN-20150113-newPrimDefTest
   LHC13b="298_20150115-1317"; #ESD
   LHC13c="299_20150115-1317"; #ESD
   LHC13e7="420_20150115-1307" ; #ESD 
   LHC13b2_efix1="416_20150115-1305" ; #ESD
   LHC13b2_efix2="417_20150115-1305" ; #ESDs
   LHC13b2_efix3="418_20150115-1307" ; #ESD
   LHC13b2_efix4="419_20150115-1307" ; #ESD
 

OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $2 = "AOD" ]; then
   TRAINPATHData=GA_pPb_AOD
else    
   TRAINPATHData=GA_pPb
fi

if [ $2 = "AOD" ]; then
   TRAINPATHMC=GA_pPb_MC_AOD
else    
   TRAINPATHMC=GA_pPb_MC
fi

if [ $LHC13e7 != "" ]; then
   OUTPUTDIR_LHC13e7=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13e7
   echo $TRAINPATHMC-$LHC13e7;
   mkdir -p $OUTPUTDIR_LHC13e7
   rm runNumbersLHC13e7FailedRunMerge.txt
   echo "copying LHC13e7 MC" 
   runNumbers=`cat runNumbersLHC13e7.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      runNumberlocal=$runNumber
      if [ $2 = "AOD" ]; then
        runNumber=$runNumber/AOD158
        
      fi
      mkdir -p $OUTPUTDIR_LHC13e7/$runNumberlocal
      if [ -f $OUTPUTDIR_LHC13e7/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumberlocal
      else 
         echo "copying /alice/sim/2013/LHC13e7/$runNumber/PWGGA/$TRAINPATHMC/$LHC13e7/$fileName"
         alien_cp alien:/alice/sim/2013/LHC13e7/$runNumber/PWGGA/$TRAINPATHMC/$LHC13e7/$fileName file:$OUTPUTDIR_LHC13e7/$runNumberlocal
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13e7/$runNumberlocal/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumberlocal
      else 
         echo $runNumberlocal >> runNumbersLHC13e7FailedRunMerge.txt
      fi     
   done;

   runNumbersBroken=`cat runNumbersLHC13e7.txt`
   runNumbersBroken=`cat runNumbersLHC13e7FailedRunMerge.txt`
   for runNumber in $runNumbersBroken; do
      echo "copying stage_1 output for " $runNumber
      mkdir -p $OUTPUTDIR_LHC13e7/$runNumber
      stageOutputs=`alien_ls /alice/sim/2013/LHC13e7/$runNumber/PWGGA/$TRAINPATHMC/$LHC13e7/Stage_2/`
      for stageOutput in $stageOutputs; do
		 testOut=`alien_ls /alice/sim/2013/LHC13e7/$runNumber/PWGGA/$TRAINPATHMC/$LHC13e7/Stage_2/$stageOutput/$fileName`
		 if [ $testOut = $fileName ]; then 
			mkdir -p $OUTPUTDIR_LHC13e7/$runNumber/Stage_2/$stageOutput
			if [ -f $OUTPUTDIR_LHC13e7/$runNumber/Stage_2/$stageOutput/$fileName ]; then
				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
			else 
				echo "copying /alice/sim/2013/LHC13e7/$runNumber/PWGGA/$TRAINPATHMC/$LHC13e7/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/sim/2013/LHC13e7/$runNumber/PWGGA/$TRAINPATHMC/$LHC13e7/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13e7/$runNumber/Stage_2/$stageOutput/
			fi
		 fi	
      done;
   done;      
fi
# 
# if [ $LHC13b2_efix1 != "" ]; then
#    OUTPUTDIR_LHC13b2_efix1=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13b2_efix1
#    echo $TRAINPATHMC-$LHC13b2_efix1;
#    mkdir -p $OUTPUTDIR_LHC13b2_efix1
#    rm runNumbersLHC13b2_efix1FailedRunMerge.txt
#    echo "copying LHC13b2_efix1 MC" 
#    runNumbers=`cat runNumbersLHC13b2_efix1.txt`
#    for runNumber in $runNumbers; do
#       echo $runNumber
#       
#       
#       mkdir -p $OUTPUTDIR_LHC13b2_efix1/$runNumber
#       if [ -f $OUTPUTDIR_LHC13b2_efix1/$runNumber/$fileName ]; then
#          echo "file " $fileName  " has already been copied for run " $runNumber
#       else 
#          echo "copying /alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix1/$fileName"
#          alien_cp alien:/alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix1/$fileName file:$OUTPUTDIR_LHC13b2_efix1/$runNumber
#       fi
#       if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13b2_efix1/$runNumber/$fileName ]; then
#          echo "file " $fileName  " been copied for sucessfully for run " $runNumber
#       else 
#          echo $runNumber >> runNumbersLHC13b2_efix1FailedRunMerge.txt
#       fi     
#    done;
# 
#    runNumbersBroken=`cat runNumbersLHC13b2_efix1FailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC13b2_efix1/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix1/Stage_2/`
#       for stageOutput in $stageOutputs; do
#          mkdir -p $OUTPUTDIR_LHC13b2_efix1/$runNumber/Stage_2/$stageOutput
#          if [ -f $OUTPUTDIR_LHC13b2_efix1/$runNumber/Stage_2/$stageOutput/$fileName ]; then
#             echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
#          else 
#             echo "copying /alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix1/Stage_2/$stageOutput/$fileName"
#             alien_cp alien:/alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix1/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13b2_efix1/$runNumber/Stage_2/$stageOutput/
#          fi
#       done;
#    done;      
# fi
# 
# if [ $LHC13b2_efix2 != "" ]; then
#    OUTPUTDIR_LHC13b2_efix2=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13b2_efix2
#    echo $TRAINPATHMC-$LHC13b2_efix2;
#    mkdir -p $OUTPUTDIR_LHC13b2_efix2
#    rm runNumbersLHC13b2_efix2FailedRunMerge.txt
#    echo "copying LHC13b2_efix2 MC" 
#    runNumbers=`cat runNumbersLHC13b2_efix2.txt`
#    for runNumber in $runNumbers; do
#       echo $runNumber
#       
#       
#       mkdir -p $OUTPUTDIR_LHC13b2_efix2/$runNumber
#       if [ -f $OUTPUTDIR_LHC13b2_efix2/$runNumber/$fileName ]; then
#          echo "file " $fileName  " has already been copied for run " $runNumber
#       else 
#          echo "copying /alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix2/$fileName"
#          alien_cp alien:/alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix2/$fileName file:$OUTPUTDIR_LHC13b2_efix2/$runNumber
#       fi
#       if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13b2_efix2/$runNumber/$fileName ]; then
#          echo "file " $fileName  " been copied for sucessfully for run " $runNumber
#       else 
#          echo $runNumber >> runNumbersLHC13b2_efix2FailedRunMerge.txt
#       fi     
#    done;
# 
#    runNumbersBroken=`cat runNumbersLHC13b2_efix2FailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC13b2_efix2/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix2/Stage_2/`
#       for stageOutput in $stageOutputs; do
#          mkdir -p $OUTPUTDIR_LHC13b2_efix2/$runNumber/Stage_2/$stageOutput
#          if [ -f $OUTPUTDIR_LHC13b2_efix2/$runNumber/Stage_2/$stageOutput/$fileName ]; then
#             echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
#          else 
#             echo "copying /alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix2/Stage_2/$stageOutput/$fileName"
#             alien_cp alien:/alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix2/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13b2_efix2/$runNumber/Stage_2/$stageOutput/
#          fi
#       done;
#    done;      
# fi
# 
# if [ $LHC13b2_efix3 != "" ]; then
#    OUTPUTDIR_LHC13b2_efix3=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13b2_efix3
#    echo $TRAINPATHMC-$LHC13b2_efix3;
#    mkdir -p $OUTPUTDIR_LHC13b2_efix3
#    rm runNumbersLHC13b2_efix3FailedRunMerge.txt
#    echo "copying LHC13b2_efix3 MC" 
#    runNumbers=`cat runNumbersLHC13b2_efix3.txt`
#    for runNumber in $runNumbers; do
#       echo $runNumber
#       
#       
#       mkdir -p $OUTPUTDIR_LHC13b2_efix3/$runNumber
#       if [ -f $OUTPUTDIR_LHC13b2_efix3/$runNumber/$fileName ]; then
#          echo "file " $fileName  " has already been copied for run " $runNumber
#       else 
#          echo "copying /alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix3/$fileName"
#          alien_cp alien:/alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix3/$fileName file:$OUTPUTDIR_LHC13b2_efix3/$runNumber
#       fi
#       if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13b2_efix3/$runNumber/$fileName ]; then
#          echo "file " $fileName  " been copied for sucessfully for run " $runNumber
#       else 
#          echo $runNumber >> runNumbersLHC13b2_efix3FailedRunMerge.txt
#       fi     
#    done;
# 
#    runNumbersBroken=`cat runNumbersLHC13b2_efix3FailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC13b2_efix3/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix3/Stage_2/`
#       for stageOutput in $stageOutputs; do
#          mkdir -p $OUTPUTDIR_LHC13b2_efix3/$runNumber/Stage_2/$stageOutput
#          if [ -f $OUTPUTDIR_LHC13b2_efix3/$runNumber/Stage_2/$stageOutput/$fileName ]; then
#             echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
#          else 
#             echo "copying /alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix3/Stage_2/$stageOutput/$fileName"
#             alien_cp alien:/alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix3/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13b2_efix3/$runNumber/Stage_2/$stageOutput/
#          fi
#       done;
#    done;      
# fi
# 
# if [ $LHC13b2_efix4 != "" ]; then
#    OUTPUTDIR_LHC13b2_efix4=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13b2_efix4
#    echo $TRAINPATHMC-$LHC13b2_efix4;
#    mkdir -p $OUTPUTDIR_LHC13b2_efix4
#    rm runNumbersLHC13b2_efix4FailedRunMerge.txt
#    echo "copying LHC13b2_efix4 MC" 
#    runNumbers=`cat runNumbersLHC13b2_efix4.txt`
#    for runNumber in $runNumbers; do
#       echo $runNumber
#       
#       
#       mkdir -p $OUTPUTDIR_LHC13b2_efix4/$runNumber
#       if [ -f $OUTPUTDIR_LHC13b2_efix4/$runNumber/$fileName ]; then
#          echo "file " $fileName  " has already been copied for run " $runNumber
#       else 
#          echo "copying /alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix4/$fileName"
#          alien_cp alien:/alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix4/$fileName file:$OUTPUTDIR_LHC13b2_efix4/$runNumber
#       fi
#       if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13b2_efix4/$runNumber/$fileName ]; then
#          echo "file " $fileName  " been copied for sucessfully for run " $runNumber
#       else 
#          echo $runNumber >> runNumbersLHC13b2_efix4FailedRunMerge.txt
#       fi     
#    done;
# 
#    runNumbersBroken=`cat runNumbersLHC13b2_efix4FailedRunMerge.txt`
#    for runNumber in $runNumbersBroken; do
#       echo "copying stage_1 output for " $runNumber
#       mkdir -p $OUTPUTDIR_LHC13b2_efix4/$runNumber
#       stageOutputs=`alien_ls /alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix4/Stage_2/`
#       for stageOutput in $stageOutputs; do
#          mkdir -p $OUTPUTDIR_LHC13b2_efix4/$runNumber/Stage_2/$stageOutput
#          if [ -f $OUTPUTDIR_LHC13b2_efix4/$runNumber/Stage_2/$stageOutput/$fileName ]; then
#             echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
#          else 
#             echo "copying /alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix4/Stage_2/$stageOutput/$fileName"
#             alien_cp alien:/alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/$TRAINPATHMC/$LHC13b2_efix4/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13b2_efix4/$runNumber/Stage_2/$stageOutput/
#          fi
#       done;
#    done;      
# fi
# 
# 
if [ $LHC13b != "" ]; then
   OUTPUTDIR_LHC13b=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC13b
   echo $TRAINPATHData-$LHC13b;
   mkdir -p $OUTPUTDIR_LHC13b
   rm runNumbersLHC13bFailedRunMerge.txt
   echo "copying LHC13b MC" 
   runNumbers=`cat runNumbersLHC13b.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      
      mkdir -p $OUTPUTDIR_LHC13b/$runNumber
      if [ -f $OUTPUTDIR_LHC13b/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then 
            echo "copying /alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/AOD139/PWGGA/$TRAINPATHData/$LHC13b/$fileName"
            alien_cp alien:/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/AOD139/PWGGA/$TRAINPATHData/$LHC13b/$fileName file:$OUTPUTDIR_LHC13b/$runNumber
         else 
            echo "copying /alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/$TRAINPATHData/$LHC13b/$fileName"
            alien_cp alien:/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/$TRAINPATHData/$LHC13b/$fileName file:$OUTPUTDIR_LHC13b/$runNumber
         fi
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13b/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC13bFailedRunMerge.txt
      fi     
   done;

   runNumbersBroken=`cat runNumbersLHC13bFailedRunMerge.txt`
   for runNumber in $runNumbersBroken; do
      echo "copying stage_1 output for " $runNumber
      mkdir -p $OUTPUTDIR_LHC13b/$runNumber
      if [ $2 = "AOD" ]; then
         stageOutputs=`alien_ls /alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/AOD139/PWGGA/$TRAINPATHData/$LHC13b/Stage_2/`
      else 
         stageOutputs=`alien_ls /alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/$TRAINPATHData/$LHC13b/Stage_2/`
      fi
      for stageOutput in $stageOutputs; do
		 testOut=`alien_ls /alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/$TRAINPATHData/$LHC13b/Stage_2/$stageOutput/$fileName`
		 if [ $testOut = $fileName ]; then 
			mkdir -p $OUTPUTDIR_LHC13b/$runNumber/Stage_2/$stageOutput
			if [ -f $OUTPUTDIR_LHC13b/$runNumber/Stage_2/$stageOutput/$fileName ]; then
				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
			else 
				if [ $2 = "AOD" ]; then
				echo "copying /alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/AOD139/PWGGA/$TRAINPATHData/$LHC13b/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/AOD139/PWGGA/$TRAINPATHData/$LHC13b/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13b/$runNumber/Stage_2/$stageOutput/
				else
				echo "copying /alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/$TRAINPATHData/$LHC13b/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/$TRAINPATHData/$LHC13b/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13b/$runNumber/Stage_2/$stageOutput/
				fi
			fi
         fi
      done;
   done;      
fi

if [ $LHC13c != "" ]; then
   OUTPUTDIR_LHC13c=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC13c
   echo $TRAINPATHData-$LHC13c;
   mkdir -p $OUTPUTDIR_LHC13c
   rm runNumbersLHC13cFailedRunMerge.txt
   echo "copying LHC13c MC" 
   runNumbers=`cat runNumbersLHC13c.txt`
   for runNumber in $runNumbers; do
      echo $runNumber
      
      
      mkdir -p $OUTPUTDIR_LHC13c/$runNumber
      if [ -f $OUTPUTDIR_LHC13c/$runNumber/$fileName ]; then
         echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         if [ $2 = "AOD" ]; then
            echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/$fileName"
            alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/$fileName file:$OUTPUTDIR_LHC13c/$runNumber
         else
            echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/$fileName"
            alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/$fileName file:$OUTPUTDIR_LHC13c/$runNumber
         fi   
      fi
      if [ $? -ne 0 ] && [ -f $OUTPUTDIR_LHC13c/$runNumber/$fileName ]; then
         echo "file " $fileName  " been copied for sucessfully for run " $runNumber
      else 
         echo $runNumber >> runNumbersLHC13cFailedRunMerge.txt
      fi     
   done;

   runNumbersBroken=`cat runNumbersLHC13cFailedRunMerge.txt`
   for runNumber in $runNumbersBroken; do
      echo "copying stage_1 output for " $runNumber
      mkdir -p $OUTPUTDIR_LHC13c/$runNumber
      if [ $2 = "AOD" ]; then
         stageOutputs=`alien_ls /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/`
      else 
         stageOutputs=`alien_ls /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/`
      fi
      for stageOutput in $stageOutputs; do
		 testOut=`alien_ls /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName`
		 if [ $testOut = $fileName ]; then 
			mkdir -p $OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput
			if [ -f $OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput/$fileName ]; then
				echo "file " $fileName  " has already been copied for run " $runNumber"/"$stageOutput
			else 
				if [ $2 = "AOD" ]; then
				echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/AOD139/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput/
				else 
				echo "copying /alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName"
				alien_cp alien:/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC13c/Stage_2/$stageOutput/$fileName file:$OUTPUTDIR_LHC13c/$runNumber/Stage_2/$stageOutput/
				fi
			fi
		 fi	
      done;
   done;      
fi