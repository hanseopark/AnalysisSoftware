# /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

BASEDIR=/home/admin1/leardini/GridOutput/PbPb
mkdir -p $BASEDIR
TRAINDIR=Legotrain-vAN-20180625-1_MBcheck;


LHC11hData=422_20180626-1832 #421_20180626-1334; #
LHC14a1a=966_20180627-0956 #964_20180626-1335; #
LHC14a1b=967_20180627-0957  #965_20180626-1335; #

LHC11hDataWithPhi=422_20180626-1832;
LHC14a1aWithPhi=966_20180627-0956;
LHC14a1bWithPhi=967_20180627-0957;

fileNameData=GammaConv_Material_103.root;
fileName14a1a=GammaConv_Material_103.root;
fileName14a1b=GammaConv_Material_103.root;

OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $1 = "no" ]; then
   echo "Not dowloading data";
elif [ $1 = "AODdata" ]; then
   TRAINPATHData=GA_PbPb_AOD
   if [ $3 = "runwise" ]; then
    rm runNumbersLHC11hNotMerged.txt
    OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData;
    echo $OUTPUTDIR_LHC11h
    runs=`cat runlists/runNumbersLHC11h.txt`
    for run in $runs; do
      echo $run
      mkdir -p $OUTPUTDIR_LHC11h/$run
      if [ -f $OUTPUTDIR_LHC11h/$run/$fileNameData ]; then
           echo "file has already been copied for run " $runNumber
      else
        alien_cp alien:/alice/data/2011/LHC11h_2/000$run/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/GammaConv_Material* file:$OUTPUTDIR_LHC11h/$run/
      fi
      if [ -f $OUTPUTDIR_LHC11h/$run/$fileNameData ]; then
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
            alien_cp alien:/alice/data/2011/LHC11h_2/000$run/ESDs/pass2/AOD145/PWGGA/$TRAINPATHData/$LHC11hData/Stage_1/$stageOutput/GammaConv_Material* file:$OUTPUTDIR_LHC11h/$run/Stage_1/$stageOutput/
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

elif [ $1 = "material" ]; then

    TRAINPATHData=GA_PbPb;
    OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/Material_LHC11h
    mkdir -p $OUTPUTDIR_LHC11h
    OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/Material_LHC11h_withPhi
    mkdir -p $OUTPUTDIR_LHC11hWithPhi

    if [ $LHC11hData != "" ]; then
        echo "Downloading " $LHC11hData
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/GammaConv_Material_* file:$OUTPUTDIR_LHC11h/
        ls $OUTPUTDIR_LHC11h/GammaConv*Material*.root > fileLHC11h_mat.txt
        fileNumbers=`cat fileLHC11h_mat.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC11h/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC11h_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC11hDataWithPhi != "" ]; then
        echo "Downloading " $LHC11hDataWithPhi
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge_runlist_8/GammaConv_Material_* file:$OUTPUTDIR_LHC11hWithPhi/
        ls $OUTPUTDIR_LHC11hWithPhi/GammaConv*Material*.root > fileLHC11h_mat.txt
        fileNumbers=`cat fileLHC11h_mat.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC11h_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

else
   TRAINPATHData=GA_PbPb;
   if [ $3 = "runwise" ]; then
        rm runNumbersLHC11hNotMerged.txt
        OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData;
        echo $OUTPUTDIR_LHC11h
        runs=`cat runlists/runNumbersLHC11h.txt`
#         for run in $runs; do
#         echo $run
#         mkdir -p $OUTPUTDIR_LHC11h/$run
#         if [ -f $OUTPUTDIR_LHC11h/$run/$fileNameData ]; then
#             echo "file has already been copied for run " $runNumber
#         else
#             alien_cp alien:/alice/data/2011/LHC11h_2/000$run/ESDs/pass2/PWGGA/$TRAINPATHData/$LHC11hData/GammaConv_Material* file:$OUTPUTDIR_LHC11h/$run/
#         fi
#         if [ -f $OUTPUTDIR_LHC11h/$run/$fileNameData ]; then
#             echo "file has already been copied for sucessfully for run " $runNumber
#         else
#             echo $run >> runNumbersLHC11hNotMerged.txt
#         fi
#         done;

        for run in $runs; do
            echo $run
            ls $OUTPUTDIR_LHC11h/$run/Gamm*.root > fileLHC11h.txt
            fileNumbers=`cat fileLHC11h.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 3 | cut -d "." -f1`
                echo $OUTPUTDIR_LHC11h/$run/GammaConv_Material_$number.root
                echo $number
                root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC11h/$run/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR_LHC11h/$run/GammaConv_Material_LHC11h_$number.root\"\,\"GammaConvMaterial_$number\"\)
            done;
        done;

    else
        OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_7
        mkdir -p $OUTPUTDIR_LHC11h
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_7/Gamm* file:$OUTPUTDIR_LHC11h/

        OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge_runlist_8
        mkdir -p $OUTPUTDIR_LHC11hWithPhi
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge_runlist_8/Gamm* file:$OUTPUTDIR_LHC11hWithPhi/

        OUTPUTDIR_LHC11h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hData/merge_runlist_8
        mkdir -p $OUTPUTDIR_LHC11h
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hData/merge_runlist_8/Gamm* file:$OUTPUTDIR_LHC11h/

        OUTPUTDIR_LHC11hWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC11hDataWithPhi/merge_runlist_7
        mkdir -p $OUTPUTDIR_LHC11hWithPhi
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC11hDataWithPhi/merge_runlist_7/Gamm* file:$OUTPUTDIR_LHC11hWithPhi/
    fi
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


if [ $2 = "no" ]; then
   echo "Not dowloading MC";
elif [ $2 = "AODmc" ]; then
   TRAINPATHMC=GA_PbPb_MC_AOD
   if [ $3 = "runwise" ]; then
    rm runNumbersLHC14a1aNotMerged.txt
    OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a;
    echo $OUTPUTDIR_LHC14a1a
    runs=`cat runlists/runNumbersLHC14a1a.txt`
    for run in $runs; do
      echo $run
      mkdir -p $OUTPUTDIR_LHC14a1a/$run
      if [ -f $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a ]; then
           echo "file has already been copied for run " $runNumber
      else
        alien_cp alien:/alice/sim/2014/LHC14a1a/$run/PWGGA/$TRAINPATHMC/$LHC14a1a/$fileName14a1a file:$OUTPUTDIR_LHC14a1a/$run/
      fi
      if [ -f $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a ]; then
           echo "file has already been copied for sucessfully for run " $runNumber
      else
         echo $run >> runNumbersLHC14a1aNotMerged.txt
      fi
    done;

    NotMergedruns=`cat runNumbersLHC14a1aNotMerged.txt`
    for run in $NotMergedruns; do
        echo "copying stage_1 output for " $run
        mkdir -p $run
        stageOutputs=`alien_ls /alice/sim/2014/LHC14a1a/$run/PWGGA/$TRAINPATHMC/$LHC14a1a/Stage_1/`
        for stageOutput in $stageOutputs; do
          mkdir -p $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput
          if [ -f $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a ]; then
            echo "file has already been copied for run " $runNumber"/"$stageOutput
          else
            alien_cp alien:/alice/sim/2014/LHC14a1a/$run/PWGGA/$TRAINPATHMC/$LHC14a1a/Stage_1/$stageOutput/$fileName14a1a file:$OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/
          fi
        done;
        counter=0;
        if [ -f $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a ]; then
          echo "file is already there for run " $run
        else
          for stageOutput in $stageOutputs; do
            if [ -f $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a ]; then
              if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root
                counter=$(($counter+1));
                echo $counter;
              else
                hadd -f $OUTPUTDIR_LHC14a1a/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1a/$run/Stage_1/$stageOutput/$fileName14a1a
                mv $OUTPUTDIR_LHC14a1a/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root
              fi
            fi
          done;
          mv $OUTPUTDIR_LHC14a1a/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a
        fi
    done;

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
        alien_cp alien:/alice/sim/2014/LHC14a1b/$run/PWGGA/$TRAINPATHMC/$LHC14a1b/$fileName14a1b file:$OUTPUTDIR_LHC14a1b/$run/
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
        stageOutputs=`alien_ls /alice/sim/2014/LHC14a1b/$run/PWGGA/$TRAINPATHMC/$LHC14a1b/Stage_1/`
        for stageOutput in $stageOutputs; do
          mkdir -p $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput
          if [ -f $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b ]; then
            echo "file has already been copied for run " $runNumber"/"$stageOutput
          else
            alien_cp alien:/alice/sim/2014/LHC14a1b/$run/PWGGA/$TRAINPATHMC/$LHC14a1b/Stage_1/$stageOutput/$fileName14a1b file:$OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/
          fi
        done;
        counter=0;
        if [ -f $OUTPUTDIR_LHC14a1b/$run/$fileName14a1b ]; then
          echo "file is already there for run " $run
        else
          for stageOutput in $stageOutputs; do
            if [ -f $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b ]; then
              if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root
                counter=$(($counter+1));
                echo $counter;
              else
                hadd -f $OUTPUTDIR_LHC14a1b/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1b/$run/Stage_1/$stageOutput/$fileName14a1b
                mv $OUTPUTDIR_LHC14a1b/$run/AnalysisResults_temp.root $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root
              fi
            fi
          done;
          mv $OUTPUTDIR_LHC14a1b/$run/intermediate_$run.root $OUTPUTDIR_LHC14a1b/$run/$fileName14a1b
        fi
    done;

    for run in $runs; do
      echo $run
      ls $OUTPUTDIR_LHC14a1a/$run/GammaConvV1_*.root > fileLHC14a1a.txt
      fileNumbers=`cat fileLHC14a1a.txt`
      for fileName in $fileNumbers; do
          echo $fileName
          number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 2 | cut -d "." -f1`
          echo $number
          root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1a/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/$run/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"GammaConvV1_$number\"\)
          root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1a/$run/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1a_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/$run/CutSelection_LHC14a1a_AOD_$number.log\"\)
      done;
      ls $OUTPUTDIR_LHC14a1b/$run/GammaConvV1_*.root > fileLHC14a1b.txt
      fileNumbersb=`cat fileLHC14a1b.txt`
      for fileName in $fileNumbersb; do
          echo $fileName
          number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
          echo $number
          root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"GammaConvV1_$number\"\)
          root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC14a1b/$run/GammaConvV1_GA_PbPb_MC_AOD_LHC14a1b_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/$run/CutSelection_LHC14a1b_AOD_$number.log\"\)
      done;
    done;

   else

#     OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_6
#     OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_6
#     OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_7
#     OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_7

    OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_7
    OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_7
    OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_6
    OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_6

    mkdir -p $OUTPUTDIR_LHC14a1a
    mkdir -p $OUTPUTDIR_LHC14a1b
    mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
    mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1a/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1b/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
#     alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1a/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1b/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
    alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_6/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/
   fi


elif [ $2 = "material" ]; then
   TRAINPATHMC=GA_PbPb_MC

        OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/Material_LHC14a1a
        OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/Material_LHC14a1b
        OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/Material_LHC14a1a_withPhi
        OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/Material_LHC14a1b_withPhi

        mkdir -p $OUTPUTDIR_LHC14a1a
        mkdir -p $OUTPUTDIR_LHC14a1b
        mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
        mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

    if [ $LHC14a1a != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/GammaConv_Material_* file:$OUTPUTDIR_LHC14a1a/
        ls $OUTPUTDIR_LHC14a1a/GammaConv*Material*.root > fileLHC14a1a_mat.txt
        fileNumbers=`cat fileLHC14a1a_mat.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC14a1a/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC14a1a_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC14a1b != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/GammaConv_Material_* file:$OUTPUTDIR_LHC14a1b/
        ls $OUTPUTDIR_LHC14a1b/GammaConv*Material*.root > fileLHC14a1b_mat.txt
        fileNumbers=`cat fileLHC14a1b_mat.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC14a1b/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC14a1b_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

        if [ $LHC14a1aWithPhi != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/GammaConv_Material_* file:$OUTPUTDIR_LHC14a1aWithPhi/
        ls $OUTPUTDIR_LHC14a1aWithPhi/GammaConv*Material*.root > fileLHC14a1aWithPhi_mat.txt
        fileNumbers=`cat fileLHC14a1aWithPhi_mat.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC14a1aWithPhi/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC14a1aWithPhi_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC14a1bWithPhi != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/GammaConv_Material_* file:$OUTPUTDIR_LHC14a1bWithPhi/
        ls $OUTPUTDIR_LHC14a1bWithPhi/GammaConv*Material*.root > fileLHC14a1bWithPhi_mat.txt
        fileNumbers=`cat fileLHC14a1bWithPhi_mat.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC14a1bWithPhi/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC14a1bWithPhi_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi


else
   TRAINPATHMC=GA_PbPb_MC
   if [ $3 = "runwise" ]; then
#         rm runNumbersLHC14a1aNotMerged.txt
#         OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a;
#         echo $OUTPUTDIR_LHC14a1a
#         runs=`cat runlists/runNumbersLHC14a1a.txt`
#         for run in $runs; do
#         echo $run
#         mkdir -p $OUTPUTDIR_LHC14a1a/$run
#         if [ -f $OUTPUTDIR_LHC14a1a/$run/$fileName14a1a ]; then
#             echo "file has already been copied for run " $runNumber
#         else
#             alien_cp alien:/alice/sim/2014/LHC14a1a/$run/PWGGA/$TRAINPATHMC/$LHC14a1a/Gamm* file:$OUTPUTDIR_LHC14a1a/$run/
#         fi
#         done;

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
            alien_cp alien:/alice/sim/2014/LHC14a1b/$run/PWGGA/$TRAINPATHMC/$LHC14a1b/Gamm* file:$OUTPUTDIR_LHC14a1b/$run/
        fi
        done;

        for run in $runs; do
            echo $run
#             ls $OUTPUTDIR_LHC14a1a/$run/Gamm*.root > fileLHC14a1a.txt
#             fileNumbers=`cat fileLHC14a1a.txt`
#             for fileName in $fileNumbers; do
#                 echo $fileName
#                 number=`echo $fileName  | cut -d "/" -f 10| cut -d "_" -f 3 | cut -d "." -f1`
#                 echo $OUTPUTDIR_LHC14a1a/$run/GammaConv_Material_$number.root
#                 echo $number
#                 root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC14a1a/$run/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR_LHC14a1a/$run/GammaConv_Material_LHC14a1a_$number.root\"\,\"GammaConvMaterial_$number\"\)
#             done;
            ls $OUTPUTDIR_LHC14a1b/$run/Gamm*.root > fileLHC14a1b.txt
            fileNumbersb=`cat fileLHC14a1b.txt`
            for fileName in $fileNumbersb; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 3 | cut -d "." -f1`
                echo $number
                root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC14a1b/$run/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR_LHC14a1b/$run/GammaConv_Material_LHC14a1b_$number.root\"\,\"GammaConvMaterial_$number\"\)
            done;
        done;

    else

        OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_7
        OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_7
        OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_8
        OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_8

        mkdir -p $OUTPUTDIR_LHC14a1a
        mkdir -p $OUTPUTDIR_LHC14a1b
        mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
        mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1a/
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1b/
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

        OUTPUTDIR_LHC14a1a=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1a/merge_runlist_8
        OUTPUTDIR_LHC14a1b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1b/merge_runlist_8
        OUTPUTDIR_LHC14a1aWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1aWithPhi/merge_runlist_7
        OUTPUTDIR_LHC14a1bWithPhi=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC14a1bWithPhi/merge_runlist_7

        mkdir -p $OUTPUTDIR_LHC14a1a
        mkdir -p $OUTPUTDIR_LHC14a1b
        mkdir -p $OUTPUTDIR_LHC14a1aWithPhi
        mkdir -p $OUTPUTDIR_LHC14a1bWithPhi

        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1a/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1a/
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1b/merge_runlist_8/Gam* file:$OUTPUTDIR_LHC14a1b/
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1aWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1aWithPhi/
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC14a1bWithPhi/merge_runlist_7/Gam* file:$OUTPUTDIR_LHC14a1bWithPhi/

    fi

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
################################# normal selection cut ######################################################

   ls $OUTPUTDIR_LHC11h/GammaConvV1_*.root > fileLHC11h.txt
   fileNumbersData=`cat fileLHC11h.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/CutSelection_LHC11h_$number.log\"\)
   done;

####################################### with phi cut ######################################################
   ls $OUTPUTDIR_LHC11hWithPhi/GammaConvV1_*.root > fileLHC11hWithPhi.txt
   fileNumbersData=`cat fileLHC11hWithPhi.txt`
   for fileName in $fileNumbersData; do
      echo $fileName
      number=`echo $fileName  | cut -d "/" -f 10 | cut -d "_" -f 2 | cut -d "." -f1`
      echo $number
      root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
      root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11hWithPhi/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11hWithPhi/CutSelection_LHC11h_$number.log\"\)
   done;


# elif [ $1 = "DCAdata" ]; then
#     runs=`cat lhc11hforDCA.txt`
#     for run in $runs; do
#       ls $OUTPUTDIR_LHC11h/$run/GammaConvV1_*.root > fileLHC11hDCA.txt
#       fileNumbersData=`cat fileLHC11hDCA.txt`
#       for fileName in $fileNumbersData; do
#           echo $fileName
#           number=`echo $fileName  | cut -d "/" -f 11 | cut -d "_" -f 2 | cut -d "." -f1`
#           echo $number
#           root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"GammaConvV1_$number\"\)
#           root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root\"\,\"$OUTPUTDIR_LHC11h/$run/CutSelection_LHC11h_$number.log\"\)
#       done;
#     done;
#     counter=0;
#     number=40;
#     echo $number
#     mkdir -p $OUTPUTDIR_LHC11h/merged
#     for run in $runs; do
#        echo "run number ---> " $run
#        if [ $counter = 0 ]; then
#           cp $OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root $OUTPUTDIR_LHC11h/intermediate.root
#           counter=$(($counter+1));
#           echo $counter;
#        else
#           hadd -f $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_merged_$number.root $OUTPUTDIR_LHC11h/intermediate.root $OUTPUTDIR_LHC11h/$run/GammaConvV1_GA_PbPb_LHC11h-pass2_$number.root
#           mv $OUTPUTDIR_LHC11h/GammaConvV1_GA_PbPb_LHC11h-pass2_merged_$number.root $OUTPUTDIR_LHC11h/intermediate.root
#        fi
#     done;
#     mv $OUTPUTDIR_LHC11h/intermediate.root  $OUTPUTDIR_LHC11h/merged/GammaConvV1_GA_PbPb_LHC11h-pass2_FinalMerge_$number.root


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

####################################### with phi cut ######################################################
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


