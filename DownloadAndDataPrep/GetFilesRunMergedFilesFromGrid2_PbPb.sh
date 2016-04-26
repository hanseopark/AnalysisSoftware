#! /bin/bash
BASEDIR=/run/media/fredi/ff316d3c-0ff4-40e6-8b73-e06f46108e1b/Grid/OutputLegoTrains
mkdir -p $BASEDIR
 
#  TRAINDIR=Legotrain-v5-05-18-AN-20130919
# #   LHC13d2=56_20130919-1817 ; #ESD
# #   LHC13d2=57_20130919-1819 ; #ESD
# #   LHC13d2b=58_20130919-1820 ; #ESD
# 
# #   LHC10hData=97_20130919-1800;
#   LHC10hData=98_20130919-1800;
# #   LHC10hData=99_20130919-1801;
# #   LHC10hData=100_20130919-1801;
# #   LHC10hData=101_20130919-1803;
# #   LHC10hData=102_20130919-1803 ;
# #   LHC10hData=103_20130919-1804;
# #   LHC10hData=104_20130919-1804;
# #   LHC10hData=105_20130919-1806;
# #   LHC10hData=106_20130919-1806;
# #   LHC10hData=107_20130919-1807;
# #   LHC10hData=108_20130919-1808;
# #   LHC10hData=109_20130919-1808;
# #   LHC10hData=110_20130919-1811;
# #   LHC10hData=111_20130919-1809;
# #   LHC13d2=98_20130919-1827  ; #AOD
# #   LHC13d2=99_20130919-1826; #AOD
# #   LHC13d2=100_20130919-1831; #AOD
# #   LHC13d2=101_20130919-1831; #AOD
# #   LHC13d2=102_20130919-1831; #AOD
# #   LHC13d2=103_20130919-1832; #AOD
# #   LHC13d2=104_20130919-1836; #AOD
# #   LHC13d2=105_20130919-1836 ; #AOD
# #   LHC13d2=106_20130919-1835; #AOD
# #   LHC13d2=107_20130919-1837; #AOD
# #   LHC13d2=108_20130919-1839; #AOD
# #   LHC13d2=109_20130919-1840; #AOD
# #   LHC13d2=110_20130919-1845; #AOD
# #   LHC13d2b=111_20130919-1854  ; #AOD
# #   LHC13d2b=112_20130919-1854; #AOD
# #   LHC13d2b=113_20130919-1856; #AOD
# #   LHC13d2b=114_20130919-1905; #AOD
# #   LHC13d2b=115_20130919-1906 ; #AOD
# #   LHC13d2b=116_20130919-1907; #AOD
# #   LHC13d2b=117_20130919-1906; #AOD
# #   LHC13d2b=118_20130919-1907; #AOD
# #   LHC13d2b=119_20130919-1910; #AOD
  
TRAINDIR=Legotrain-v5-05-27-AN-20131018
   LHC13d2=72_20131021-1045 ; #ESD
   LHC13d2b=73_20131021-1046  ; #ESD
   LHC10hData=120_20131018-1524;
  
OUTPUTDIR=$BASEDIR/$TRAINDIR
# if [ $1 = "AOD" ]; then
#    TRAINPATHData=GA_PbPb_AOD
# else
   TRAINPATHData=GA_PbPb
   nFilesToCopy=10;
# fi   
OUTPUTDIR_LHC10h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10hData
echo $TRAINPATHData-$LHC10hData;
# OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData

# if [ $1 = "AOD" ]; then
#    TRAINPATHMC=GA_PbPb_MC_AOD
# else
   TRAINPATHMC=GA_PbPb_MC
# fi   
OUTPUTDIR_LHC13d2=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13d2
echo $TRAINPATHMC-$LHC13d2;
OUTPUTDIR_LHC13d2b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13d2b
echo $TRAINPATHMC-$LHC13d2b;
mkdir -p $OUTPUTDIR_LHC10h
# mkdir -p $OUTPUTDIR_LHC11h
mkdir -p $OUTPUTDIR_LHC13d2
mkdir -p $OUTPUTDIR_LHC13d2b

fileName=GammaConvV1_QA_546000106009266304480300000.root

   echo "copying LHC10h data" 
   runNumbers=`cat runNumbersLHC10h.tex`
   for runNumber in $runNumbers; do
      echo $runNumber
      mkdir -p $OUTPUTDIR_LHC10h/$runNumber
      nFiles=`ls $OUTPUTDIR_LHC10h/$runNumber | wc -l`
#       if [ $nFiles = $nFilesToCopy ]; then
#          echo "files copied alread for run " $runNumber
      if [ -f $OUTPUTDIR_LHC10h/$runNumber/$fileName ]; then
           echo "file " $fileName  " has already been copied for run " $runNumber
      else 
         alien_cp alien:/alice/data/2010/LHC10h/000$runNumber/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC10hData/$fileName file:$OUTPUTDIR_LHC10h/$runNumber
      fi
   done;

   echo "copying LHC13d2 MC" 
   runNumbers=`cat runNumbersLHC10h.tex`
   for runNumber in $runNumbers; do
      echo $runNumber
      mkdir -p $OUTPUTDIR_LHC13d2/$runNumber
      nFiles=`ls $OUTPUTDIR_LHC13d2/$runNumber | wc -l`
#       if [ $nFiles = $nFilesToCopy ]; then
#          echo "files copied alread for run " $runNumber
      if [ -f $OUTPUTDIR_LHC13d2/$runNumber/$fileName ]; then
           echo "file " $fileName  " has already been copied for run " $runNumber
      else  
         alien_cp alien:/alice/sim/2013/LHC13d2/$runNumber/PWGGA/$TRAINPATHMC/$LHC13d2/$fileName file:$OUTPUTDIR_LHC13d2/$runNumber
      fi
   done;
#     
   echo "copying LHC13d2b MC"  
   runNumbers=`cat runNumbersLHC10h.tex`
   for runNumber in $runNumbers; do
      echo $runNumber
      mkdir -p $OUTPUTDIR_LHC13d2b/$runNumber
        nFiles=`ls $OUTPUTDIR_LHC13d2b/$runNumber | wc -l`
#       if [ $nFiles = $nFilesToCopy ]; then
#          echo "files copied alread for run " $runNumber
      if [ -f $OUTPUTDIR_LHC13d2b/$runNumber/$fileName ]; then
           echo "file " $fileName  " has already been copied for run " $runNumber
     
      else  
         alien_cp alien:/alice/sim/2013/LHC13d2b/$runNumber/PWGGA/$TRAINPATHMC/$LHC13d2b/$fileName file:$OUTPUTDIR_LHC13d2b/$runNumber
      fi
   done;
   
