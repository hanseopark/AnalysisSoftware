# /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

BASEDIR=/home/admin1/leardini/GridOutput/PbPb
mkdir -p $BASEDIR
TRAINDIR=Legotrain-vAN-20160211-1;


# #========= tight cut train =========
LHC11hData=265_20160215-1750; 

LHC14a1a=636_20160215-1752; 
LHC14a1b=639_20160215-1752;


# # #========= tight cut train =========
# LHC11hData=260_20160211-1824; 
# LHC11hDataWithPhi=261_20160211-1825; 
# 
# LHC14a1a=632_20160211-1832; 
# LHC14a1b=634_20160211-1833;
# LHC14a1aWithPhi=633_20160211-1832;
# LHC14a1bWithPhi=635_20160211-1833 ;



# #========= std (slow train for data), no MC weighting =========
#-> for cross check of weighting effect:
# LHC11hData=250_20160121-1452; 
# LHC11hDataWithPhi=251_20160121-1452; 

# LHC14a1a=627_20160121-1450; 
# LHC14a1b=628_20160121-1450;
# LHC14a1aWithPhi=626_20160121-1450;
# LHC14a1bWithPhi=625_20160121-1450;

# #========= std (slow train for data), MC multi weighting =========
# LHC11hData=241_20151209-1800; #cut 226, 230, 235, 239
# LHC11hDataWithPhi=242_20151209-1704; 

# LHC14a1a=623_20160114-1458; 
# LHC14a1b=624_20160114-1458;
# LHC14a1aWithPhi=623_20160114-1458;
# LHC14a1bWithPhi=624_20160114-1458;

#---> for added signals part:
#std slow, double rej, mult not cut
# LHC14a1a=598_20151210-0914; 
# LHC14a1b=600_20151210-1012;
# LHC14a1aWithPhi=599_20151210-0914;
# LHC14a1bWithPhi=601_20151210-1033;


# #========= std (slow train) no multi cut =========
# LHC11hData=241_20151209-1800; #cut 226, 230, 235, 239
# LHC11hDataWithPhi=242_20151209-1704; 

#std slow, double rej, mult not cut
# LHC14a1a=598_20151210-0914; 
# LHC14a1b=600_20151210-1012;
# LHC14a1aWithPhi=599_20151210-0914;
# LHC14a1bWithPhi=601_20151210-1033;


#Chi2 30. and psi pair 1D (not linked)
# LHC14a1a=602_20151209-1801; 
# LHC14a1b=604_20151209-1819;
# LHC14a1aWithPhi=603_20151209-1711;
# LHC14a1bWithPhi=605_20151209-1818;


#first syst check
# LHC14a1a=606_20151219-1059; 
# LHC14a1b=608_20151219-1059;
# LHC14a1aWithPhi=607_20151219-1059;
# LHC14a1bWithPhi=609_20151219-1059;


#second syst check
# LHC14a1a=611_20151220-1102; 
# LHC14a1b=613_20151220-1103;
# LHC14a1aWithPhi=617_20160103-1107;
# LHC14a1bWithPhi=614_20151220-1103;


#third syst check
# LHC11hData=244_20160105-1530; #cut 235 nuovo
# LHC11hDataWithPhi=245_20160105-1530; 
# 
# LHC14a1a=618_20160104-1012; 
# LHC14a1b=622_20160104-1212;
# LHC14a1aWithPhi=619_20160104-1013;
# LHC14a1bWithPhi=621_20160104-1213;



#========= only MC 0-10%, diffrent V0 status =========
# LHC14a1a=595_20151130-1818; 
# LHC14a1a=596_20151130-1917;
# LHC14a1aWithPhi=597_20151130-1826 ;

# #========= std with V0 event cut (slow train) =========
# LHC11hData=238_20151126-1238; 
# LHC11hDataWithPhi=239_20151126-1041; 
# 
# LHC14a1a=591_20151126-1347; 
# LHC14a1b=593_20151126-1042;
# LHC14a1aWithPhi=592_20151126-1042;
# LHC14a1bWithPhi=594_20151126-1042;


	
OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $1 = "AODdata" ]; then
   TRAINPATHData=GA_PbPb_AOD
else
   TRAINPATHData=GA_PbPb
fi  

OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_1
mkdir -p $OUTPUTDIR_LHC11h
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC11h/


# OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_7
# mkdir -p $OUTPUTDIR_LHC11h
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC11h/
# 
# OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge 
# mkdir -p $OUTPUTDIR_LHC11hWithPhi
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge/Gam* file:$OUTPUTDIR_LHC11hWithPhi/
 
 

if [ $2 = "AODmc" ]; then
   TRAINPATHMC=GA_PbPb_MC_AOD
else
   TRAINPATHMC=GA_PbPb_MC
 	OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge
	OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge
	mkdir -p $OUTPUTDIR_LHC14a1a
	mkdir -p $OUTPUTDIR_LHC14a1b

fi  

OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_1
OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_1

# OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_6
# OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_6
# OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_7
# OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_7

# OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_7
# OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_7
# OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_6
# OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_6

mkdir -p $OUTPUTDIR_LHC14a1a
mkdir -p $OUTPUTDIR_LHC14a1b
mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
mkdir -p $OUTPUTDIR_LHC14a1bWithPhi


alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC14a1a/
alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC14a1b/

# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1a/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1b/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1a/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1b/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/



# 
# # OUTPUTDIR_LHC14a1c=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1c
# # mkdir -p $OUTPUTDIR_LHC14a1c
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1c/merge/Gam* file:$OUTPUTDIR_LHC14a1c/
# 

if [ $2 = "AODmc" ]; then
##################################### normal selection cut ######################################################
   ls $OUTPUTDIR_LHC14a1a/GammaConvV1_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1` 
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/CutSelection_LHC14a1a_AOD_$number.log\"\)

   done;
   
   ls $OUTPUTDIR_LHC14a1b/GammaConvV1_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/CutSelection_LHC14a1b_AOD_$number.log\"\)
   done;

# #    ls $OUTPUTDIR_LHC14a1c/GammaConvV1_*.root > fileLHC14a1c.txt
# #    fileNumbersb=`cat fileLHC14a1c.txt`
# #    for fileName in $fileNumbersb; do
# #       echo $fileName
# #       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
# #       echo $number
# #       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1c/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1c/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1c_$number.root\"\,\"GammaConvV1_$number\"\)
# #       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1c/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1c_$number.root\"\,\"$OUTPUTDIR_LHC14a1c/CutSelection_LHC14a1c_AOD_$number.log\"\)
# #    done;

   
######################################### with phi cut ######################################################
#    ls $OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_*.root > fileLHC14a1aWithPhi.txt
#    fileNumbers=`cat fileLHC14a1aWithPhi.txt`
#    for fileName in $fileNumbers; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1` 
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/CutSelection_LHC14a1a_AOD_$number.log\"\)
# 
#    done;
#    
#    ls $OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_*.root > fileLHC14a1bWithPhi.txt
#    fileNumbersb=`cat fileLHC14a1bWithPhi.txt`
#    for fileName in $fileNumbersb; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/CutSelection_LHC14a1b_AOD_$number.log\"\)
#    done;


else

   ls $OUTPUTDIR_LHC14a1a/GammaConvV1_*.root > fileLHC14a1a.txt
   fileNumbers=`cat fileLHC14a1a.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1a/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/CutSelection_LHC14a1a_$number.log\"\)

   done;
   
   ls $OUTPUTDIR_LHC14a1b/GammaConvV1_*.root > fileLHC14a1b.txt
   fileNumbersb=`cat fileLHC14a1b.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/CutSelection_LHC14a1b_$number.log\"\)
   done;
   
 
# # #   rm $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1a_LHC14a1b_*.root
# # #   ls $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1a_*.root > filesForMerging.txt
# # #   filesForMerging=`cat filesForMerging.txt`
# # #   for fileName in $filesForMerging; do
# # #      echo $fileName
# # #      number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 6 | cut -d "." -f1`
# # #      echo $number
# # #      if [ -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1c_$number.root ]; then
# # #         hadd -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1a_LHC14a1b_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root
# # #      fi
# # #   done
# 
fi



if [ $1 = "AODdata" ]; then
################################## normal selection cut ######################################################

   ls $OUTPUTDIR_LHC11h/GammaConvV1_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/CutSelection_LHC11h_$number.log\"\)
   done;

######################################## with phi cut ######################################################
#    ls $OUTPUTDIR_LHC11hWithPhi/GammaConvV1_*.root > fileLHC11hWithPhi.txt
#    fileNumbersData=`cat fileLHC11hWithPhi.txt`
#    for fileName in $fileNumbersData; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/CutSelection_LHC11h_$number.log\"\)
#    done;
# 
# else
#    ls $OUTPUTDIR_LHC11h/GammaConvV1_*.root > fileLHC11h.txt
#    fileNumbersData=`cat fileLHC11h.txt`
#    for fileName in $fileNumbersData; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/CutSelection_LHC11h_$number.log\"\)
#    done;
fi 




###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

# 
# 
# ##=========================================================================================
# ##=========================================================================================
# # OUTPUTDIR_LHC12a11a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12a11a/
# # OUTPUTDIR_LHC12a11b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12a11b/
# # OUTPUTDIR_LHC12a11d=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12a11d/
# # OUTPUTDIR_LHC12a11e=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC12a11e/
# # 
# # mkdir -p $OUTPUTDIR_LHC12a11a
# # mkdir -p $OUTPUTDIR_LHC12a11b
# # mkdir -p $OUTPUTDIR_LHC12a11d
# # mkdir -p $OUTPUTDIR_LHC12a11e
# # 
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC12a11a/merge/Gam* file:$OUTPUTDIR_LHC12a11a
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC12a11b/merge/Gam* file:$OUTPUTDIR_LHC12a11b
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC12a11d/merge/Gam* file:$OUTPUTDIR_LHC12a11d
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC12a11e/merge/Gam* file:$OUTPUTDIR_LHC12a11e
# # 
# #    ls $OUTPUTDIR_LHC12a11a/GammaConvV1_*.root > fileLHC12a11a.txt
# #    fileNumbers=`cat fileLHC12a11a.txt`
# #    for fileName in $fileNumbers; do
# #       echo $fileName
# #       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
# #       echo $number
# #       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12a11a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC12a11a/GammaConvV1_GA_PbPb_MC_LHC12a11a_$number.root\"\,\"GammaConvV1_$number\"\)
# #       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC12a11a/GammaConvV1_GA_PbPb_MC_LHC12a11a_$number.root\"\,\"$OUTPUTDIR_LHC12a11a/CutSelection_LHC12a11a_$number.log\"\)
# # 
# #    done;
# #    
# #    
# #    ls $OUTPUTDIR_LHC12a11b/GammaConvV1_*.root > fileLHC12a11b.txt
# #    fileNumbers=`cat fileLHC12a11b.txt`
# #    for fileName in $fileNumbers; do
# #       echo $fileName
# #       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
# #       echo $number
# #       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12a11b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC12a11b/GammaConvV1_GA_PbPb_MC_LHC12a11b_$number.root\"\,\"GammaConvV1_$number\"\)
# #       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC12a11b/GammaConvV1_GA_PbPb_MC_LHC12a11b_$number.root\"\,\"$OUTPUTDIR_LHC12a11b/CutSelection_LHC12a11b_$number.log\"\)
# # 
# #    done;
# # 
# #    
# #    ls $OUTPUTDIR_LHC12a11d/GammaConvV1_*.root > fileLHC12a11d.txt
# #    fileNumbers=`cat fileLHC12a11d.txt`
# #    for fileName in $fileNumbers; do
# #       echo $fileName
# #       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
# #       echo $number
# #       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12a11d/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC12a11d/GammaConvV1_GA_PbPb_MC_LHC12a11d_$number.root\"\,\"GammaConvV1_$number\"\)
# #       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC12a11d/GammaConvV1_GA_PbPb_MC_LHC12a11d_$number.root\"\,\"$OUTPUTDIR_LHC12a11d/CutSelection_LHC12a11d_$number.log\"\)
# # 
# #    done;
# # 
# # 
# #    
# #    ls $OUTPUTDIR_LHC12a11e/GammaConvV1_*.root > fileLHC12a11e.txt
# #    fileNumbers=`cat fileLHC12a11e.txt`
# #    for fileName in $fileNumbers; do
# #       echo $fileName
# #       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
# #       echo $number
# #       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC12a11e/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC12a11e/GammaConvV1_GA_PbPb_MC_LHC12a11e_$number.root\"\,\"GammaConvV1_$number\"\)
# #       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC12a11e/GammaConvV1_GA_PbPb_MC_LHC12a11e_$number.root\"\,\"$OUTPUTDIR_LHC12a11e/CutSelection_LHC12a11e_$number.log\"\)
# # 
# #    done;
# ##=========================================================================================
# ##=========================================================================================
# 
#   
# 
# 


#========= std with double rejection (slow train) =========
# LHC11hData=234_20151119-1738; #226 and 234
# LHC11hDataWithPhi=237_20151119-1720; #228

# LHC14a1a=576_20151119-1742; #226, 227_eta, 229_eta
# LHC14a1b=578_20151119-1742;
# LHC14a1aWithPhi=577_20151119-1742;
# LHC14a1bWithPhi=579_20151119-1742;

# LHC14a1c=584_20151120-1018; #234

#-------------------------------- no psi pair
# LHC11hData=235_20151119-1719 #230, 235, 239
# LHC11hDataWithPhi=236_20151119-1719; #232, 237, 240
# 
# LHC14a1a=580_20151119-1745; #230, 231_eta, 233_eta
# LHC14a1b=582_20151119-1746;
# LHC14a1aWithPhi=581_20151119-1757;
# LHC14a1bWithPhi=583_20151119-1747;

#-------------------------------- no chi2
# LHC14a1a=585_20151120-1018; #235, 236_eta, 238_eta
# LHC14a1b=587_20151120-1019;
# LHC14a1aWithPhi=586_20151120-1018;
# LHC14a1bWithPhi=588_20151120-1019;


# LHC14a1a=589_20151120-1019; #239 (no phi), #240 (with phi)
# LHC14a1bWithPhi=589_20151120-1019;

#========= std with new V0cut =========
# LHC11hData=226_20151102-1458;
# LHC11hDataWithPhi=227_20151102-1458;
# 
# LHC14a1a=567_20151103-1029;
# LHC14a1b=569_20151103-1030;
# LHC14a1aWithPhi=568_20151103-1030;
# LHC14a1bWithPhi=570_20151103-1030;


#========= std with new V0cut and double rejection =========
# LHC11hData=226_20151102-1458;
# LHC11hDataWithPhi=227_20151102-1458;
# 
# LHC14a1a=563_20151102-1415;
# LHC14a1b=565_20151102-1416;
# LHC14a1aWithPhi=564_20151102-1416;
# LHC14a1bWithPhi=566_20151102-1416;



#========== vertex cross checks ========
#no qt and psi pair 224
#no qt 222 
#no psi pair 220 
# LHC14a1a=561_20151030-1051;
# LHC14a1b=562_20151030-1051;

#standard
# LHC14a1a=556_20151025-1543;
# LHC14a1b=555_20151025-1543;

#only phi cut - 218, 219
#plus phi cut - 216, 217
# LHC14a1a=548_20151020-1120;
# LHC14a1b=550_20151020-1038;

#plus R_min - 214, 215
#plus pion dE/dx - 212, 213
#plus qT - 210, 211
# LHC14a1a=547_20151020-1119;
# LHC14a1b=549_20151020-1037;

#plus cosPA - 208, 209
# LHC14a1a=551_20151020-1048;
# LHC14a1b=552_20151020-1048;

#plus single pt - 206, 207
# LHC14a1a=551_20151020-1048;
# LHC14a1b=552_20151020-1048;

#dE/dx and psi pair - 204, 205
# LHC14a1a=545_20151016-1637;
# LHC14a1b=546_20151016-1655;

#only psi pair - 202, 203
# LHC14a1a=543_20151015-1452;
# LHC14a1b=544_20151015-1457;

#only dE/dx - 200, 201
# LHC14a1a=545_20151016-1637;
# LHC14a1b=546_20151016-1655;

#open cut - 198, 199
# LHC14a1a=541_20151013-1806;
# LHC14a1b=542_20151013-1806;



#========= slow train with xrootd =========
# LHC11hData=215_20150921-1117;
# LHC11hDataWithPhi=216_20150921-1116;
# 
# LHC14a1a=527_20150921-1117;
# LHC14a1b=529_20150921-1117;
# LHC14a1aWithPhi=528_20150921-1132;
# LHC14a1bWithPhi=530_20150921-1117;


# LHC14a1a=523_20150907-1340;  #no MC weighting
# LHC14a1b=526_20150907-1305; #no MC weighting


#========= pile-up for R > 35cm study ==========
#   LHC11hData=167_20150803-1043;
#   LHC11hData=168_20150803-1044;


# #========= V0 efficiency - 198 ==========
#   LHC14a1a=225_20150808-1324;
#   LHC14a1b=224_20150808-1324;

# #========= V0 efficiency - 162 ==========
#   LHC14a1a=222_20150804-0801;
#   LHC14a1b=223_20150803-1855;


#=== AMPT simulations - 1 ====
# LHC12a11a=215_20150729-1006;
# LHC12a11b=212_20150729-0957; 
# LHC12a11d=213_20150729-1008;
# LHC12a11e=214_20150729-1008;

#=== AMPT simulations - 2 ====
# LHC12a11a=218_20150803-1828;
# LHC12a11b=219_20150803-1828; 
# LHC12a11d=220_20150803-1828;
# LHC12a11e=221_20150803-1840;


#=== open dEdx cut ====
#   LHC11hData=192_20150727-1243;
#   LHC11hDataWithPhi=193_20150727-1243;
#   LHC14a1a=492_20150727-1059;
#   LHC14a1b=494_20150727-1100;
#   LHC14a1aWithPhi=493_20150727-1059;
#   LHC14a1bWithPhi=495_20150727-1100;


#=== R > 35cm ====
#   LHC11hData=192_20150727-1243;
#   LHC11hDataWithPhi=193_20150727-1243;
#   LHC14a1a=496_20150727-1100;
#   LHC14a1b=498_20150727-1101;
#   LHC14a1aWithPhi=497_20150727-1101;
#   LHC14a1bWithPhi=499_20150727-1101;
    
    
#=== re-do cut 146 (new tag) ====
#   LHC11hData=194_20150729-1034;
#   LHC11hDataWithPhi=195_20150729-1034;
#   LHC14a1a=500_20150729-0907;
#   LHC14a1b=502_20150729-0908;
#   LHC14a1aWithPhi=501_20150729-0908;
#   LHC14a1bWithPhi=503_20150729-0908;
    
#==================================================================================
#   LHC11hDataWithPhi=187_20150630-1308; #186, 188
#   LHC14a1aWithPhi=457_20150711-1550;  #186, 187eta, 189eta
#   LHC14a1bWithPhi=459_20150711-1600; 

#   LHC14a1aWithPhi=458_20150711-1558; #188, 187pi0, 189pi0
#   LHC14a1bWithPhi=460_20150711-1601; 

#   LHC11hData=185_20150630-1307; # 154, 158    #----> done
#   LHC11hDataWithPhi=186_20150630-1307; #156, 160
    
#   LHC14a1a=449_20150702-1626;  #154, 155eta, 157eta
#   LHC14a1b=451_20150702-1627; 
#   LHC14a1aWithPhi=450_20150702-1626; # 156, 155pi, 157pi
#   LHC14a1bWithPhi=452_20150702-1627; 
# 
#   LHC14a1a=453_20150702-1627;  #158, 159eta, 161eta
#   LHC14a1b=455_20150702-1628; 
#   LHC14a1aWithPhi=454_20150702-1627; #160, 159pi, 161pi
#   LHC14a1bWithPhi=456_20150702-1628; 

#   LHC11hData=183_20150630-1258; # 138, 142, 146      #----> done
#   LHC11hDataWithPhi=184_20150630-1258; # 140, 144, 148
# 
#   LHC14a1a=480_20150721-1027; #wrong 437_20150702-1527;  #138, 139eta, 141eta
#   LHC14a1b=482_20150721-1631 #wrong 439_20150702-1528; 
#   LHC14a1aWithPhi=481_20150721-1101; #wrong 438_20150702-1527; #140, 139pi, 141pi
#   LHC14a1bWithPhi=483_20150722-2156 #wrong 440_20150702-1528; 

#   LHC14a1a=484_20150723-1042; # wrong 441_20150702-1530;  #142, 143eta, 145eta
#   LHC14a1b=486_20150723-1043; #wrong 443_20150702-1531; 
#   LHC14a1aWithPhi=485_20150723-1043; # wrong 442_20150702-1530; #144, 143pi, 145pi
#   LHC14a1bWithPhi=487_20150723-1043; #wrong 444_20150702-1531; 

#   LHC14a1a=445_20150702-1531;  #146, 147eta, 149eta ----> wrong (see above)
#   LHC14a1b=447_20150702-1625; 
#   LHC14a1aWithPhi=446_20150702-1532; #148, 147pi, 149pi
#   LHC14a1bWithPhi=448_20150702-1626; 
# 


#------------------------------------------------------------
# LHC11hData=196_20150731-0835;    #126, 130, 134  
# LHC11hDataWithPhi=197_20150731-0835; #128, 132, 136
# 
# #new qt cut!!!
#   LHC14a1a=504_20150730-2013;#130, 131eta, 133eta
#   LHC14a1b=506_20150730-2014;
#   LHC14a1aWithPhi=505_20150730-2013;#132, 131pi, 133pi
#   LHC14a1bWithPhi=507_20150730-2014; 
# 
#--------------> data redone becuase something off, qt cut changed
# #     LHC11hData=181_20150630-1032; # 126, 130, 134   #----> done
# #     LHC11hDataWithPhi=182_20150630-1032; #128, 132, 136
# 
# #     LHC14a1a=425_20150702-1632; #126, 127eta, 129eta  
# #     LHC14a1b=427_20150702-1632; 
# #     LHC14a1aWithPhi=426_20150702-1632; #128, 127pi, 129pi
# #     LHC14a1bWithPhi=428_20150702-1633; 
# 
# #     LHC14a1a=472_20150719-2149 #wrong 429_20150702-1633;  #130, 131eta, 133eta
# #     LHC14a1b=474_20150720-1024 #wrong 431_20150702-1633; 
# #     LHC14a1aWithPhi=473_20150720-1019 #wrong 429_20150702-1633; #132, 131pi, 133pi
# #     LHC14a1bWithPhi=475_20150720-1030 #wrong 432_20150702-1633; 
# 
# #     LHC14a1a=476_20150720-1740 #wrong 433_20150702-1418;  #134, 135eta, 137eta
# #     LHC14a1b=478_20150721-1026 #wrong 435_20150702-1420; 
# #     LHC14a1aWithPhi=477_20150720-1815 #wrong 434_20150702-1419; #136, 135pi, 137pi
# #     LHC14a1bWithPhi=479_20150720-1740 #wrong 436_20150702-1422; 
    

#   LHC11hData=179_20150630-1028; #118, 122, 150
#   LHC11hDataWithPhi=180_20150630-1029; # 120, 124, 152

#   LHC14a1a=468_20150719-1927; #wrong: 413_20150630-1316; #118, 119eta, 121eta
#   LHC14a1b=470_20150719-1939  #wrong: 415_20150630-1318;
#   LHC14a1aWithPhi=469_20150719-1917; #worng: 414_20150630-1316; #120, 119pi,  121pi
#   LHC14a1bWithPhi=471_20150719-2149; #wrong: 416_20150630-1319; 

#   LHC14a1a=488_20150723-1044; #wrong 417_20150630-1453;   #122, 123eta, 125eta
#   LHC14a1b=490_20150723-1042; #wrong 419_20150630-1517; 
#   LHC14a1aWithPhi=489_20150723-1044; #wrong 418_20150630-1514; #124, 123pi, 125pi
#   LHC14a1bWithPhi=491_20150721-1656; #wrong 420_20150630-1455; 

#   LHC14a1a=421_20150630-1510;  #150, 151eta, 153eta
#   LHC14a1b=423_20150630-1511; 
#   LHC14a1aWithPhi=422_20150630-1511; #152, 151pi, 153pi 
#   LHC14a1bWithPhi=424_20150630-1512; 

#   LHC11hData=177_20150619-2051; # 106, 110, 114
#   LHC11hDataWithPhi=189_20150720-1037; #better stat 178_20150619-2051; #108, 112, 116

#   LHC14a1a=401_20150626-1142;  #106, 107eta, 109eta
#   LHC14a1b=403_20150626-1144; 
#   LHC14a1aWithPhi=402_20150626-1144; #108, 107pi0, 109pi
#   LHC14a1bWithPhi=404_20150626-1145; 
    
#   LHC14a1a=405_20150630-1018;  #110, 111eta, 113eta
#   LHC14a1b=407_20150630-1019; 
#   LHC14a1aWithPhi=406_20150630-1019; #112, 11pi, 113pi
#   LHC14a1bWithPhi=408_20150630-1020; 

#   LHC14a1a=409_20150630-1313;  #114, 115ta, 117eta
#   LHC14a1b=411_20150630-1314; 
#   LHC14a1aWithPhi=410_20150630-1313; # 116, 115pi, 117pi
#   LHC14a1bWithPhi=412_20150630-1314; 

#   LHC11hData=175_20150615-1939;  #162 not flat, 102
#   LHC11hDataWithPhi=176_20150615-1940;  #164 not flat, 104

#   LHC14a1a=397_20150624-1312;  #102, 103eta, 105eta
#   LHC14a1b=399_20150624-1313; 
#   LHC14a1aWithPhi=398_20150624-1312; #104, 103pi, 105pi
#   LHC14a1bWithPhi=400_20150624-1314; 

#   LHC11hData=190_20150720-1040; #better stat 173_20150615-1331; #90, 94, 98
#   LHC11hDataWithPhi=191_20150720-1418; #better stat 174_20150615-1332; # 92, 96, 100

#   LHC14a1a=385_20150617-0811; # 90, 91eta, 93eta  
#   LHC14a1b=387_20150617-0812; 
#   LHC14a1aWithPhi=386_20150617-0811; #92, 91pi0, 93pi 
#   LHC14a1bWithPhi=388_20150617-0813; 

#   LHC14a1a=389_20150619-2053;  # 94, 95eta, 97eta
#   LHC14a1b=391_20150619-2054; 
#   LHC14a1aWithPhi=390_20150619-2054; #96 95pi, 97pi
#   LHC14a1bWithPhi=392_20150619-2056; 
    
#   LHC14a1a=393_20150624-1059;  # 98, 99eta, 101eta
#   LHC14a1b=395_20150624-1100; 
#   LHC14a1aWithPhi=394_20150624-1059; #100, 99pi, 101pi
#   LHC14a1bWithPhi=396_20150624-1100; 
        
#   LHC11hData=171_20150614-1506; #74, 82, 86
#   LHC11hDataWithPhi=172_20150614-1508; #76, 84, 88

#   LHC14a1a=373_20150614-1513;  #74, 75eta, 77eta
#   LHC14a1b=375_20150614-1514; 
#   LHC14a1aWithPhi=374_20150614-1513; #76, 75pi0, 77pi0
#   LHC14a1bWithPhi=376_20150614-1515; 
    
#   LHC14a1a=377_20150615-0837;  #82, 83eta, 85eta
#   LHC14a1b=379_20150615-0838; 
#   LHC14a1aWithPhi=378_20150615-0837; # 84, 83pi0, 85pi0
#   LHC14a1bWithPhi=380_20150615-0838; 
    
#   LHC14a1a=381_20150616-1151; #86, 87eta, 89eta  
#   LHC14a1b=383_20150616-1152; 
#   LHC14a1aWithPhi=382_20150616-1151; #88, 87pi0, 89pi0 
#   LHC14a1bWithPhi=384_20150616-1153; 
    
#   #data cut is with std analysis
#   LHC14a1a=369_20150612-1010;  #70, 71eta, 73eta
#   LHC14a1b=371_20150612-1008; 
#   LHC14a1aWithPhi=466_20150719-1912; #wrong: 370_20150612-1009; #72, 71pi0, 73pi0
#   LHC14a1bWithPhi=467_20150719-1912; #wrong: 372_20150612-1007; 


#   LHC11hData=169_20150612-1147; #std 
#   LHC11hDataWithPhi=170_20150612-1147; #std
# 
#   LHC14a1a=365_20150610-1352;  #std
#   LHC14a1b=367_20150610-1352; #std
#   LHC14a1aWithPhi=366_20150610-1408; #std
#   LHC14a1bWithPhi=368_20150610-1402; #std
    
        
#   TRAINDIR=Legotrain-vAN-20140408_1stIteration
#   LHC11hData=45_20140408-2045; #AOD
#   LHC14a1a=156_20140410-1126; #ESD
#   LHC14a1b=155_20140410-1125; #ESD
#   LHC14a1c=154_20140410-1124; #ESD
    
#   TRAINDIR=Legotrain-vAN-20140408_1stIteration_Pi0
#   LHC14a1a=151_20140410-2239; #ESD pi0
#   LHC14a1b=161_20140412-2336; #ESD pi0
#   LHC14a1c=160_20140412-2334; #ESD pi0
#   LHC14a1b=178_20140415-1847; #AOD pi0
#   LHC14a1c=180_20140415-1849; #AOD pi0

#   TRAINDIR=Legotrain-vAN-20140408_1stIteration_Eta
#   LHC14a1a=157_20140410-2240; #ESD 
#   LHC14a1b=158_20140412-2335; #ESD eta
#   LHC14a1c=159_20140412-2334; #ESD eta
#   LHC14a1b=177_20140415-1846; #AOD eta
#   LHC14a1c=179_20140415-1848; #AOD eta

#   TRAINDIR=Legotrain-vAN-20140408_2ndIteration
#   LHC11hData=50_20140415-1904; #AOD
#   LHC14a1a=181_20140416-1127; #AOD pi0
#   LHC14a1b=183_20140416-1144; #AOD pi0
#   LHC14a1c=154_20140410-1124; #AOD
#   LHC14a1a=182_20140416-1143; #AOD eta
#   LHC14a1b=184_20140416-1144; #AOD eta

#   TRAINDIR=Legotrain-vAN-20140408_3rdIteration
#   LHC11hData=50_20140415-1904; #AOD
#   LHC14a1a=185_20140417-1300; #AOD pi0
#   LHC14a1b=187_20140417-1300; #AOD pi0
#   LHC14a1c=154_20140410-1124; #AOD
#   LHC14a1a=186_20140417-1300; #AOD eta
#   LHC14a1b=188_20140417-1302; #AOD eta
