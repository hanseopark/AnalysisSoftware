# /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

BASEDIR=/home/admin1/leardini/GridOutput/PbPb
mkdir -p $BASEDIR
TRAINDIR=Legotrain-vAN-20170211-1;


#photon cut studies
LHC11hData=263_20170216-1326; #400_20170212-1233;
LHC14a1a=350_20170216-1102; #351_20170216-1335
LHC14a1b=352_20170216-1322; #353_20170216-1103

LHC11hDataWithPhi=263_20170216-1326; #332_20160706-1144;
LHC14a1aWithPhi=350_20170216-1102; #351_20170216-1335
LHC14a1bWithPhi=352_20170216-1322; #353_20170216-1103

fileNameData=AnalysisResults_2040.root;
fileName14a1a=AnalysisResults.root;
fileName14a1b=AnalysisResults.root;

OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $1 = "AODdata" ]; then
   TRAINPATHData=GA_PbPb_AOD
   if [ $3 = "runwise" ]; then
    rm runNumbersLHC11hNotMerged.txt
    OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData;
    echo $OUTPUTDIR_LHC11h
    runs=`cat runlists/runNumbersLHC11h.txt`
    for run in $runs; do
      echo $run
      mkdir -p $OUTPUTDIR_LHC11h/$run
      if [ -f $OUTPUTDIR_LHC11h/$run/$fileName ]; then
           echo "file has already been copied for run " $runNumber
      else
        alien_cp alien:/alice/data/2011/LHC11h_2/000$run/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/AnalysisResults_2040.root file:$OUTPUTDIR_LHC11h/$run/AnalysisResults.root
      fi
      if [ -f $OUTPUTDIR_LHC11h/$run/$fileName ]; then
           echo "file has already been copied for sucessfully for run " $runNumber
      else
         echo $run >> runNumbersLHC11hNotMerged.txt
      fi
    done;

    NotMergedruns=`cat runNumbersLHC11hNotMerged.txt`
    for run in $NotMergedruns; do
        echo "copying stage_1 output for " $run
        mkdir -p $run
        stageOutputs=`alien_ls /alice/data/2011/LHC11h_2/000$run/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/Stage_1/`
        for stageOutput in $stageOutputs; do
          mkdir -p $OUTPUTDIR_LHC11h/$run/Stage_1/$stageOutput
          if [ -f $OUTPUTDIR_LHC11h/$run/Stage_1/$stageOutput/$fileNameData ]; then
            echo "file has already been copied for run " $runNumber"/"$stageOutput
          else
            alien_cp alien:/alice/data/2011/LHC11h_2/000$run/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/Stage_1/$stageOutput/AnalysisResults_2040.root file:$OUTPUTDIR_LHC11h/$run/Stage_1/$stageOutput/AnalysisResults.root
          fi
        done;
#         counter=0;
#         if [ -f $OUTPUTDIR_LHC11h/$run/$fileNameData ]; then
#           echo "file is already there for run " $run
#         else
#           for stageOutput in $stageOutputs; do
#             if [ -f $OUTPUTDIR_LHC11h/$run/Stage_1/$stageOutput/$fileNameData ]; then
#               if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC11h/$run/Stage_1/$stageOutput/$fileNameData $OUTPUTDIR_LHC11h/$run/intermediate_$run.root
#                 counter=$(($counter+1));
#                 echo $counter;
#               else
#                 hadd -f $OUTPUTDIR_LHC11h/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC11h/$run/intermediate_$run.root $OUTPUTDIR_LHC11h/$run/Stage_1/$stageOutput/$fileNameData
#                 mv $OUTPUTDIR_LHC11h/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC11h/$run/intermediate_$run.root
#               fi
#             fi
#           done;
#           mv $OUTPUTDIR_LHC11h/$run/intermediate_$run.root $OUTPUTDIR_LHC11h/$run/$fileNameData
#         fi
    done;

   else

    OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_7
    OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge_runlist_8
#     mkdir -p $OUTPUTDIR_LHC11h
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/Gamm* file:$OUTPUTDIR_LHC11h/
#     mkdir -p $OUTPUTDIR_LHC11hWithPhi
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge_runlist_8/Gamm* file:$OUTPUTDIR_LHC11hWithPhi/

#     OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_7
#     mkdir -p $OUTPUTDIR_LHC11h
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC11h/
#
#     OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge
#     mkdir -p $OUTPUTDIR_LHC11hWithPhi
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge/Gam* file:$OUTPUTDIR_LHC11hWithPhi/
#
   fi
else
   TRAINPATHData=GA_PbPb
    OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_7
    mkdir -p $OUTPUTDIR_LHC11h
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/Gamm* file:$OUTPUTDIR_LHC11h/

    OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge_runlist_8
    mkdir -p $OUTPUTDIR_LHC11hWithPhi
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge_runlist_8/Gamm* file:$OUTPUTDIR_LHC11hWithPhi/

fi

# # for the standard analysis:
# OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_1
# mkdir -p $OUTPUTDIR_LHC11h
# alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC11h/


# #for the DCA analysis (runwise or stageoutput wise):
# if [ $1 = "DCAdata" ]; then
#
#   OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_1
# #   stageOutputs=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_1/Stage_1`
# # #   for stageOutput in $stageOutputs; do
# # #     mkdir -p $OUTPUTDIR_LHC11h/Stage_1/$stageOutput
#   runs=`cat lhc11hforDCA.txt`
#   for run in $runs; do
#     mkdir -p $OUTPUTDIR_LHC11h/$run
#     alien_cp alien:/alice/data/2011/LHC11h_2/000$run/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC11hData/Gamm* file:$OUTPUTDIR_LHC11h/$run
#   done;
# fi


if [ $2 = "AODmc" ]; then
   TRAINPATHMC=GA_PbPb_MC_AOD
   if [ $3 = "runwise" ]; then
#     rm runNumbersLHC14a1aNotMerged.txt
#     OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a;
#     echo $OUTPUTDIR_LHC14a1a
#     runs=`cat runlists/runNumbersLHC14a1a.txt`
#     for run in $runs; do
#       echo $run
#       mkdir -p $OUTPUTDIR_LHC14a1a/$run
#       if [ -f $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a ]; then
#            echo "file has already been copied for run " $runNumber
#       else
#         alien_cp alien:/alice/sim/2014/LHC14a1a/$run/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1a/$fileName14a1a file:$OUTPUTDIR_LHC14a1a/$run/
#       fi
#       if [ -f $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a ]; then
#            echo "file has already been copied for sucessfully for run " $runNumber
#       else
#          echo $run >> runNumbersLHC14a1aNotMerged.txt
#       fi
#     done;
#
#     NotMergedruns=`cat runNumbersLHC14a1aNotMerged.txt`
#     for run in $NotMergedruns; do
#         echo "copying stage_1 output for " $run
#         mkdir -p $run
#         stageOutputs=`alien_ls /alice/sim/2014/LHC14a1a/$run/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1a/Stage_1/`
#         for stageOutput in $stageOutputs; do
#           mkdir -p $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput
#           if [ -f $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a ]; then
#             echo "file has already been copied for run " $runNumber"/"$stageOutput
#           else
#             alien_cp alien:/alice/sim/2014/LHC14a1a/$run/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1a/Stage_1/$stageOutput/$fileName14a1a file:$OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/
#           fi
#         done;
#         counter=0;
#         if [ -f $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a ]; then
#           echo "file is already there for run " $run
#         else
#           for stageOutput in $stageOutputs; do
#             if [ -f $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a ]; then
#               if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root
#                 counter=$(($counter+1));
#                 echo $counter;
#               else
#                 hadd -f $OUTPUTDIR_LHC14a1a/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a
#                 mv $OUTPUTDIR_LHC14a1a/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root
#               fi
#             fi
#           done;
#           mv $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a
#         fi
#     done;

    rm runNumbersLHC14a1bNotMerged.txt
    OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b;
    echo $OUTPUTDIR_LHC14a1b
    runs=`cat runlists/runNumbersLHC14a1a.txt`
    for run in $runs; do
      echo $run
      mkdir -p $OUTPUTDIR_LHC14a1b/$run
      if [ -f $OUTPUTDIR_LHC14a1b/$run/$fileName14a1b ]; then
           echo "file has already been copied for run " $runNumber
      else
        alien_cp alien:/alice/sim/2014/LHC14a1b/$run/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1b/$fileName14a1b file:$OUTPUTDIR_LHC14a1b/$run/
      fi
      if [ -f $OUTPUTDIR_LHC14a1b/$run/$fileName14a1b ]; then
           echo "file has already been copied for sucessfully for run " $runNumber
      else
         echo $run >> runNumbersLHC14a1bNotMerged.txt
      fi
    done;

    NotMergedruns=`cat runNumbersLHC14a1bNotMerged.txt`
    for run in $NotMergedruns; do
        echo "copying stage_1 output for " $run
        mkdir -p $run
        stageOutputs=`alien_ls /alice/sim/2014/LHC14a1b/$run/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1b/Stage_1/`
        for stageOutput in $stageOutputs; do
          mkdir -p $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput
          if [ -f $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b ]; then
            echo "file has already been copied for run " $runNumber"/"$stageOutput
          else
            alien_cp alien:/alice/sim/2014/LHC14a1b/$run/AOD149/PWGGA/$TRAINPATHMC/$LHC14a1b/Stage_1/$stageOutput/$fileName14a1b file:$OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/
          fi
        done;
#         counter=0;
#         if [ -f $OUTPUTDIR_LHC14a1b/$run/$fileName14a1b ]; then
#           echo "file is already there for run " $run
#         else
#           for stageOutput in $stageOutputs; do
#             if [ -f $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b ]; then
#               if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root
#                 counter=$(($counter+1));
#                 echo $counter;
#               else
#                 hadd -f $OUTPUTDIR_LHC14a1b/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b
#                 mv $OUTPUTDIR_LHC14a1b/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root
#               fi
#             fi
#           done;
#           mv $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1b/$run/$fileName14a1b
#         fi
    done;

   else

    OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_6
    OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_6
    OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_7
    OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_7

#     OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_7
#     OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_7
#     OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_6
#     OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_6

#     mkdir -p $OUTPUTDIR_LHC14a1a
#     mkdir -p $OUTPUTDIR_LHC14a1b
#     mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
#     mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1a/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1b/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1a/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1b/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/
   fi

else
   TRAINPATHMC=GA_PbPb_MC

    OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_7
    OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_7
    OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_8
    OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_8

#     OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_8
#     OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_8
#     OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_7
#     OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_7

    mkdir -p $OUTPUTDIR_LHC14a1a
    mkdir -p $OUTPUTDIR_LHC14a1b
    mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
    mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1a/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1b/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1a/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1b/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

fi

# # OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_1
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC14a1a/
# # OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_1
# # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/Gam* file:$OUTPUTDIR_LHC14a1b/




# if [ $2 = "DCAmc" ]; then
#
#   OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_1
#   stageOutputs=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/Stage_1`
#   for stageOutput in $stageOutputs; do
#     mkdir -p $OUTPUTDIR_LHC14a1b/Stage_1/$stageOutput
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_1/Stage_1/$stageOutput/Gamm* file:$OUTPUTDIR_LHC14a1b/Stage_1/$stageOutput/
#   done;
#   runs=`cat lhc11hforDCA.txt`
#   for run in $runs; do
#     mkdir -p $OUTPUTDIR_LHC14a1b/$run
#     alien_cp alien:/alice/sim/2014/LHC14a1b/$run/PWGGA/$TRAINPATHMC/$LHC14a1b/Gamm* file:$OUTPUTDIR_LHC14a1b/$run
#   done;

# fi



# # # OUTPUTDIR_LHC14a1c=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1c
# # # mkdir -p $OUTPUTDIR_LHC14a1c
# # # alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1c/merge/Gam* file:$OUTPUTDIR_LHC14a1c/
# # #


if [ $2 = "AODmc" ]; then
#################################### normal selection cut ######################################################
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

#    ls $OUTPUTDIR_LHC14a1c/GammaConvV1_*.root > fileLHC14a1c.txt
#    fileNumbersb=`cat fileLHC14a1c.txt`
#    for fileName in $fileNumbersb; do
#       echo $fileName
#       number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
#       echo $number
#       root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1c/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1c/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1c_$number.root\"\,\"GammaConvV1_$number\"\)
#       root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1c/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1c_$number.root\"\,\"$OUTPUTDIR_LHC14a1c/CutSelection_LHC14a1c_AOD_$number.log\"\)
#    done;


######################################## with phi cut ######################################################
   ls $OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_*.root > fileLHC14a1aWithPhi.txt
   fileNumbers=`cat fileLHC14a1aWithPhi.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/CutSelection_LHC14a1a_AOD_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_*.root > fileLHC14a1bWithPhi.txt
   fileNumbersb=`cat fileLHC14a1bWithPhi.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/CutSelection_LHC14a1b_AOD_$number.log\"\)
   done;


# elif [ $2 = "DCAmc" ]; then
#     echo $OUTPUTDIR_LHC14a1b
#     runs=`cat lhc11hforDCA.txt`
#     for run in $runs; do
#       ls $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_*.root > fileLHC11hDCA.txt
#       fileNumbersDCAmc=`cat fileLHC11hDCA.txt`
#       for fileName in $fileNumbersDCAmc; do
#           echo $fileName
#           number=`echo $fileName  | cut -d "/" -f 11 | cut -d "_" -f 2 | cut -d "." -f1`
#           echo $number
# #           root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
# #           root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/$run/CutSelection_LHC14a1b_$number.log\"\)
#       done;
#     done;
#
#     counter=0;
#     number=40;
#     mkdir -p $OUTPUTDIR_LHC14a1b/merged/
#     echo $number
#     for run in $runs; do
#        echo "run number ---> " $run
#        mv $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root
#         if [ -f $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root ]; then
#           if [ $counter = 0 ]; then
#             cp $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root $OUTPUTDIR_LHC14a1b/intermediate.root
#             counter=$(($counter+1));
#             echo $counter;
#           else
#             hadd -f $OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_merged_$number.root $OUTPUTDIR_LHC14a1b/intermediate.root $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root
#             mv $OUTPUTDIR_LHC14a1b/GammaConvV1_GA_PbPb_MC_LHC14a1b_merged_$number.root $OUTPUTDIR_LHC14a1b/intermediate.root
#           fi
#         else
#           echo "file not there for run " $run
#         fi
#     done;
#     mv $OUTPUTDIR_LHC14a1b/intermediate.root  $OUTPUTDIR_LHC14a1b/merged/GammaConvV1_GA_PbPb_MC_LHC14a1b_FinalMerge_$number.root


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

# ######################################## with phi cut ######################################################
   ls $OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_*.root > fileLHC14a1aWithPhi.txt
   fileNumbers=`cat fileLHC14a1aWithPhi.txt`
   for fileName in $fileNumbers; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1aWithPhi/CutSelection_LHC14a1a_$number.log\"\)

   done;

   ls $OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_*.root > fileLHC14a1bWithPhi.txt
   fileNumbersb=`cat fileLHC14a1bWithPhi.txt`
   for fileName in $fileNumbersb; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConvV1_GA_PbPb_MC_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1bWithPhi/CutSelection_LHC14a1b_$number.log\"\)
   done;

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
   ls $OUTPUTDIR_LHC11hWithPhi/GammaConvV1_*.root > fileLHC11hWithPhi.txt
   fileNumbersData=`cat fileLHC11hWithPhi.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/CutSelection_LHC11h_$number.log\"\)
   done;


# # elif [ $1 = "DCAdata" ]; then
# #     runs=`cat lhc11hforDCA.txt`
# #     for run in $runs; do
# #       ls $OUTPUTDIR_LHC11h/$run/GammaConvV1_*.root > fileLHC11hDCA.txt
# #       fileNumbersData=`cat fileLHC11hDCA.txt`
# #       for fileName in $fileNumbersData; do
# #           echo $fileName
# #           number=`echo $fileName  | cut -d "/" -f 11 | cut -d "_" -f 2 | cut -d "." -f1`
# #           echo $number
# #           root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
# #           root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/$run/CutSelection_LHC11h_$number.log\"\)
# #       done;
# #     done;
# #     counter=0;
# #     number=40;
# #     echo $number
# #     mkdir -p $OUTPUTDIR_LHC11h/merged
# #     for run in $runs; do
# #        echo "run number ---> " $run
# #        if [ $counter = 0 ]; then
# #           cp $OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root $OUTPUTDIR_LHC11h/intermediate.root
# #           counter=$(($counter+1));
# #           echo $counter;
# #        else
# #           hadd -f $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_merged_$number.root $OUTPUTDIR_LHC11h/intermediate.root $OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root
# #           mv $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_merged_$number.root $OUTPUTDIR_LHC11h/intermediate.root
# #        fi
# #     done;
# #     mv $OUTPUTDIR_LHC11h/intermediate.root  $OUTPUTDIR_LHC11h/merged/GammaConvV1_GA_PbPb_LHC11h-pass2_FinalMerge_$number.root


else

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
   ls $OUTPUTDIR_LHC11hWithPhi/GammaConvV1_*.root > fileLHC11hWithPhi.txt
   fileNumbersData=`cat fileLHC11hWithPhi.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/CutSelection_LHC11h_$number.log\"\)
   done;

fi


