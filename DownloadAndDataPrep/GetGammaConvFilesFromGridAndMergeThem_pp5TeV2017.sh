#! /bin/bash
source basicFunction.sh

# download script for pp 5TeV from 2017
BASEDIR=/home/admin1/leardini/GridOutput/pp
mkdir -p $BASEDIR
TRAINDIR=Legotrain-vAN-20180429-1_fullJJ

### Data
LHC17pqMETA=2340_20180420-1304
#QA pile On 337_20180220-1343
#QA pile Off 336_20180220-1253
LHC17p_fast=$LHC17pqMETA\_child_1
LHC17p_woSDD=$LHC17pqMETA\_child_2
LHC17q_fast=$LHC17pqMETA\_child_3
LHC17q_woSDD=$LHC17pqMETA\_child_4
LHC17p_wSDD=;
LHC17q_wSDD=;

# LHC17p_fast=$LHC17pqMETA\_child_1
# LHC17p_wSDD=$LHC17pqMETA\_child_2
# LHC17p_woSDD=$LHC17pqMETA\_child_3
# LHC17q_fast=$LHC17pqMETA\_child_4
# LHC17q_wSDD=$LHC17pqMETA\_child_5
# LHC17q_woSDD=$LHC17pqMETA\_child_6

### MC
LHC17lMETA=3273_20180420-1304 #with MBW 3276_20180426-0910
#3263_20180414-0206 #->ESD #791_20180417-1040 #->AOD #
#3253_20180331-2110 #wo MB weights
#3254_20180331-2111  #with MB weights
#QA pile On 644_20180220-1349
#QA pile Off 643_20180220-1232
LHC17l3b_fast=$LHC17lMETA\_child_1
LHC17l4b_fast=;
LHC17l3b_woSDD=$LHC17lMETA\_child_2
LHC17l4b_woSDD=;
LHC17l3b_cent=;
LHC17l4b_cent=;

# LHC17l3b_fast=$LHC17lMETA\_child_1
# LHC17l3b_cent=$LHC17lMETA\_child_2
# LHC17l3b_woSDD=$LHC17lMETA\_child_3
# LHC17l4b_fast=$LHC17lMETA\_child_4
# LHC17l4b_cent=$LHC17lMETA\_child_5
# LHC17l4b_woSDD=$LHC17lMETA\_child_6

### MC Jet Jet
LHC18b8META=819_20180430-1553
LHC18b8_fast=$LHC18b8META\_child_1
LHC18b8_cent=$LHC18b8META\_child_2
LHC18b8_woSDD=$LHC18b8META\_child_3

OUTPUTDIR=$BASEDIR/$TRAINDIR
TRAINPATHData=GA_pp #_AOD
TRAINPATHMC=GA_pp_MC_AOD

NSlashes=9
NSlashes2=10
NSlashes3=11

# fileName=AnalysisResults.root;
fileName=GammaConvV1_400.root;
# fileName=GammaConv_Material_121.root;

if [ $1 = "no" ]; then
   echo "Not dowloading data";
elif [ $1 = "yes" ]; then

    OUTPUTDIR_LHC17p_fast=$BASEDIR/$TRAINDIR/LHC17p_fast
    OUTPUTDIR_LHC17p_wSDD=$BASEDIR/$TRAINDIR/LHC17p_wSDD
    OUTPUTDIR_LHC17p_woSDD=$BASEDIR/$TRAINDIR/LHC17p_woSDD
    mkdir -p $OUTPUTDIR_LHC17p_fast
    mkdir -p $OUTPUTDIR_LHC17p_wSDD
    mkdir -p $OUTPUTDIR_LHC17p_woSDD

    OUTPUTDIR_LHC17q_fast=$BASEDIR/$TRAINDIR/LHC17q_fast
    OUTPUTDIR_LHC17q_wSDD=$BASEDIR/$TRAINDIR/LHC17q_wSDD
    OUTPUTDIR_LHC17q_woSDD=$BASEDIR/$TRAINDIR/LHC17q_woSDD
    mkdir -p $OUTPUTDIR_LHC17q_fast
    mkdir -p $OUTPUTDIR_LHC17q_wSDD
    mkdir -p $OUTPUTDIR_LHC17q_woSDD


    if [ $3 = "runwise" ]; then

        ################################ LHC17p_fast
        rm runNumbersLHC17pfastNotMerged.txt
        echo $OUTPUTDIR_LHC17p_fast
        runs=`cat runlists/runNumbersLHC17p_all.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17p_fast/$run
            if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/$fileName file:$OUTPUTDIR_LHC17p_fast/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17pfastNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17pfastNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17p_fast/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17p_fast/$run/$fileName $OUTPUTDIR_LHC17p_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17p_fast/$run/$fileName $OUTPUTDIR_LHC17p_fast/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17p_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_fast/$run/$fileName
            fi
        done;


        ################################ LHC17p_wSDD
        rm runNumbersLHC17pwSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17p_wSDD
        runs=`cat runlists/runNumbersLHC17p_all.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17p_wSDD/$run
            if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/$fileName file:$OUTPUTDIR_LHC17p_wSDD/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17pwSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17pwSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17p_wSDD/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17p_wSDD/$run/$fileName $OUTPUTDIR_LHC17p_wSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17p_wSDD/$run/$fileName $OUTPUTDIR_LHC17p_wSDD/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17p_wSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_wSDD/$run/$fileName
            fi
        done;


        ################################ LHC17p_woSDD
        rm runNumbersLHC17pwoSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17p_woSDD
        runs=`cat runlists/runNumbersLHC17p_all.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run
            if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/$fileName file:$OUTPUTDIR_LHC17p_woSDD/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17pwoSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17pwoSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17p_woSDD/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName $OUTPUTDIR_LHC17p_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17p_woSDD/$run/$fileName $OUTPUTDIR_LHC17p_woSDD/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17p_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_woSDD/$run/$fileName
            fi
        done;


        ################################ LHC17q_fast
        rm runNumbersLHC17qfastNotMerged.txt
        echo $OUTPUTDIR_LHC17q_fast
        runs=`cat runlists/runNumbersLHC17q_all.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17q_fast/$run
            if [ -f $OUTPUTDIR_LHC17q_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/$fileName file:$OUTPUTDIR_LHC17q_fast/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17q_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17qfastNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17qfastNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17q_fast/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17q_fast/$run/$fileName $OUTPUTDIR_LHC17q_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17q_fast/$run/$fileName $OUTPUTDIR_LHC17q_fast/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17q_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_fast/$run/$fileName
            fi
        done;


        ################################ LHC17q_wSDD
        rm runNumbersLHC17qwSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17q_wSDD
        runs=`cat runlists/runNumbersLHC17q_all.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17q_wSDD/$run
            if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/$fileName file:$OUTPUTDIR_LHC17q_wSDD/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17qwSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17qwSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17q_wSDD/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17q_wSDD/$run/$fileName $OUTPUTDIR_LHC17q_wSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17q_wSDD/$run/$fileName $OUTPUTDIR_LHC17q_wSDD/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17q_wSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_wSDD/$run/$fileName
            fi
        done;


        ################################ LHC17q_woSDD
        rm runNumbersLHC17qwoSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17q_woSDD
        runs=`cat runlists/runNumbersLHC17q_all.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17q_woSDD/$run
            if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/$fileName file:$OUTPUTDIR_LHC17q_woSDD/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17qwoSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17qwoSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17q_woSDD/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17q_woSDD/$run/$fileName $OUTPUTDIR_LHC17q_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17q_woSDD/$run/$fileName $OUTPUTDIR_LHC17q_woSDD/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17q_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_woSDD/$run/$fileName
            fi
        done;

        if [ $fileName != "AnalysisResults.root" ]; then
            runs=`cat runlists/runNumbersLHC17p_all.txt`
            for run in $runs; do
                ls $OUTPUTDIR_LHC17p_fast/$run/GammaConvV1_*.root > fileLHC17p_fast.txt
                fileNumbers=`cat fileLHC17p_fast.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_fast/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17p_fast/$run/GammaConvV1_LHC17p_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17p_fast/$run/GammaConvV1_LHC17p_fast-pass1_$number.root\"\,\"$OUTPUTDIR_LHC17p_fast/$run/CutSelection_LHC17p_fast_$number.log\"\)
                done;

                ls $OUTPUTDIR_LHC17p_wSDD/$run/GammaConvV1_*.root > fileLHC17p_wSDD.txt
                fileNumbers=`cat fileLHC17p_wSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_wSDD/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17p_wSDD/$run/GammaConvV1_LHC17p_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17p_wSDD/$run/GammaConvV1_LHC17p_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR_LHC17p_wSDD/$run/CutSelection_LHC17p_wSDD_$number.log\"\)
                done;

                ls $OUTPUTDIR_LHC17p_woSDD/$run/GammaConvV1_*.root > fileLHC17p_woSDD.txt
                fileNumbers=`cat fileLHC17p_woSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_woSDD/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17p_woSDD/$run/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17p_woSDD/$run/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR_LHC17p_woSDD/$run/CutSelection_LHC17p_woSDD_$number.log\"\)
                done;
            done;

            runs=`cat runlists/runNumbersLHC17q_all.txt`
            for run in $runs; do

                ls $OUTPUTDIR_LHC17q_fast/$run/GammaConvV1_*.root > fileLHC17q_fast.txt
                fileNumbers=`cat fileLHC17q_fast.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_fast/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17q_fast/$run/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17q_fast/$run/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"$OUTPUTDIR_LHC17q_fast/$run/CutSelection_LHC17q_fast_$number.log\"\)
                done;

                ls $OUTPUTDIR_LHC17q_wSDD/$run/GammaConvV1_*.root > fileLHC17q_wSDD.txt
                fileNumbers=`cat fileLHC17q_wSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_wSDD/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17q_wSDD/$run/GammaConvV1_LHC17q_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17q_wSDD/$run/GammaConvV1_LHC17q_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR_LHC17q_wSDD/$run/CutSelection_LHC17q_wSDD_$number.log\"\)
                done;

                ls $OUTPUTDIR_LHC17q_woSDD/$run/GammaConvV1_*.root > fileLHC17q_woSDD.txt
                fileNumbers=`cat fileLHC17q_woSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_woSDD/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17q_woSDD/$run/GammaConvV1_LHC17q_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17q_woSDD/$run/GammaConvV1_LHC17q_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR_LHC17q_woSDD/$run/CutSelection_LHC17q_woSDD_$number.log\"\)
                done;
            done;
        fi
    fi


    if [ $LHC17p_fast != "" ]; then
        echo "Downloading " $LHC17p_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_fast/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC17p_fast/
        ls $OUTPUTDIR_LHC17p_fast/GammaConvV1_*.root > fileLHC17p_fast.txt
        fileNumbers=`cat fileLHC17p_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_fast_$number.log\"\)
        done;
    fi

    if [ $LHC17p_wSDD != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_wSDD/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC17p_wSDD/
        ls $OUTPUTDIR_LHC17p_wSDD/GammaConvV1_*.root > fileLHC17p_wSDD.txt
        fileNumbers=`cat fileLHC17p_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_wSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_wSDD_$number.log\"\)
        done;
    fi

    if [ $LHC17p_woSDD != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_woSDD/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC17p_woSDD/
        ls $OUTPUTDIR_LHC17p_woSDD/GammaConvV1_*.root > fileLHC17p_woSDD.txt
        fileNumbers=`cat fileLHC17p_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_woSDD_$number.log\"\)
        done;
    fi

    if [ $LHC17q_fast != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17q_fast/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC17q_fast/
        ls $OUTPUTDIR_LHC17q_fast/GammaConvV1_*.root > fileLHC17q_fast.txt
        fileNumbers=`cat fileLHC17q_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17q_fast_$number.log\"\)
        done;
    fi

    if [ $LHC17q_wSDD != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17q_wSDD/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC17q_wSDD/
        ls $OUTPUTDIR_LHC17q_wSDD/GammaConvV1_*.root > fileLHC17q_wSDD.txt
        fileNumbers=`cat fileLHC17q_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_wSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17q_wSDD_$number.log\"\)
        done;
    fi

    if [ $LHC17q_woSDD != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17q_woSDD/merge_runlist_1/GammaConvV1* file:$OUTPUTDIR_LHC17q_woSDD/
        ls $OUTPUTDIR_LHC17q_woSDD/GammaConvV1_*.root > fileLHC17q_woSDD.txt
        fileNumbers=`cat fileLHC17q_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17q_woSDD_$number.log\"\)
        done;
    fi
fi

if [ $1 = "mergefills" ]; then

    OUTPUTDIR_LHC17p_fast=$BASEDIR/$TRAINDIR/LHC17p_fast
    OUTPUTDIR_LHC17p_wSDD=$BASEDIR/$TRAINDIR/LHC17p_wSDD
    OUTPUTDIR_LHC17p_woSDD=$BASEDIR/$TRAINDIR/LHC17p_woSDD
    mkdir -p $OUTPUTDIR_LHC17p_fast
    mkdir -p $OUTPUTDIR_LHC17p_wSDD
    mkdir -p $OUTPUTDIR_LHC17p_woSDD

    OUTPUTDIR_LHC17q_fast=$BASEDIR/$TRAINDIR/LHC17q_fast
    OUTPUTDIR_LHC17q_wSDD=$BASEDIR/$TRAINDIR/LHC17q_wSDD
    OUTPUTDIR_LHC17q_woSDD=$BASEDIR/$TRAINDIR/LHC17q_woSDD
    mkdir -p $OUTPUTDIR_LHC17q_fast
    mkdir -p $OUTPUTDIR_LHC17q_wSDD
    mkdir -p $OUTPUTDIR_LHC17q_woSDD

    fileName="GammaConvV1_400.root"
    fileName2="GammaConvV1_400"
    ######################################################################### fill 1
    runs=`cat runlists/runNumbersLHC17p_fill1.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17p_woSDD/$run/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill1.root $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill1.root $OUTPUTDIR_LHC17p_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill1.root

    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17p_fast/$run/$fileName $OUTPUTDIR_LHC17p_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17p_fast/$fileName2\_fill1.root $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17p_fast/$fileName2\_fill1.root $OUTPUTDIR_LHC17p_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$fileName2\_fill1.root

    ######################################################################### fill 2
    runs=`cat runlists/runNumbersLHC17p_fill2.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17p_woSDD/$run/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill2.root $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill2.root $OUTPUTDIR_LHC17p_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill2.root
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17p_fast/$run/$fileName $OUTPUTDIR_LHC17p_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17p_fast/$fileName2\_fill2.root $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17p_fast/$fileName2\_fill2.root $OUTPUTDIR_LHC17p_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$fileName2\_fill2.root

    ######################################################################### fill 3
    runs=`cat runlists/runNumbersLHC17p_fill3.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17p_woSDD/$run/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill3.root $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill3.root $OUTPUTDIR_LHC17p_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$fileName2\_fill3.root
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17p_fast/$run/$fileName $OUTPUTDIR_LHC17p_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17p_fast/$fileName2\_fill3.root $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17p_fast/$fileName2\_fill3.root $OUTPUTDIR_LHC17p_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$fileName2\_fill3.root

    ######################################################################### fill 4
    runs=`cat runlists/runNumbersLHC17q_fast.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17q_woSDD/$run/$fileName $OUTPUTDIR_LHC17q_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17q_woSDD/$fileName2\_fill4.root $OUTPUTDIR_LHC17q_woSDD/intermediate.root $OUTPUTDIR_LHC17q_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17q_woSDD/$fileName2\_fill4.root $OUTPUTDIR_LHC17q_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17q_woSDD/intermediate.root $OUTPUTDIR_LHC17q_woSDD/$fileName2\_fill4.root
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17q_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17q_fast/$run/$fileName $OUTPUTDIR_LHC17q_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17q_fast/$fileName2\_fill4.root $OUTPUTDIR_LHC17q_fast/intermediate.root $OUTPUTDIR_LHC17q_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17q_fast/$fileName2\_fill4.root $OUTPUTDIR_LHC17q_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17q_fast/intermediate.root $OUTPUTDIR_LHC17q_fast/$fileName2\_fill4.root
fi

if [ $1 = "material" ]; then

    OUTPUTDIR_LHC17p_fast=$BASEDIR/$TRAINDIR/Material_LHC17p_fast_LowInt
    mkdir -p $OUTPUTDIR_LHC17p_fast
    if [ $LHC17p_fast != "" ]; then
        echo "Downloading " $LHC17p_fast " in " $OUTPUTDIR_LHC17p_fast
        runs=`cat runlists/runNumbersLHC17p_LowInt.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17p_fast/$run
            if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/$fileName file:$OUTPUTDIR_LHC17p_fast/$run/
            fi
        done;

        hadd -f $OUTPUTDIR_LHC17p_fast/$fileName $OUTPUTDIR_LHC17p_fast/282008/$fileName $OUTPUTDIR_LHC17p_fast/282016/$fileName $OUTPUTDIR_LHC17p_fast/282021/$fileName $OUTPUTDIR_LHC17p_fast/282025/$fileName $OUTPUTDIR_LHC17p_fast/282030/$fileName $OUTPUTDIR_LHC17p_fast/282031/$fileName

        ls $OUTPUTDIR_LHC17p_fast/GammaConv*Material*.root > fileLHC17p_fast.txt
        fileNumbers=`cat fileLHC17p_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17p_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17p_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi


    OUTPUTDIR_LHC17p_woSDD=$BASEDIR/$TRAINDIR/Material_LHC17p_woSDD_LowInt
    mkdir -p $OUTPUTDIR_LHC17p_woSDD
    if [ $LHC17p_woSDD != "" ]; then
        echo "Downloading " $LHC17p_woSDD " in " $OUTPUTDIR_LHC17p_woSDD
        runs=`cat runlists/runNumbersLHC17p_LowInt.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run
            if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_woSDD/$fileName file:$OUTPUTDIR_LHC17p_woSDD/$run/
            fi
        done;

        hadd -f $OUTPUTDIR_LHC17p_woSDD/$fileName $OUTPUTDIR_LHC17p_woSDD/282008/$fileName $OUTPUTDIR_LHC17p_woSDD/282016/$fileName $OUTPUTDIR_LHC17p_woSDD/282021/$fileName $OUTPUTDIR_LHC17p_woSDD/282025/$fileName $OUTPUTDIR_LHC17p_woSDD/282030/$fileName $OUTPUTDIR_LHC17p_woSDD/282031/$fileName

        ls $OUTPUTDIR_LHC17p_woSDD/GammaConv*Material*.root > fileLHC17p_woSDD.txt
        fileNumbers=`cat fileLHC17p_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17p_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17p_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi


fi

if [ $1 = "mergeTriggers" ]; then

	rm $OUTPUTDIR/GammaConvV1_LHC17pq_fast-pass1_*.root
	ls $OUTPUTDIR/GammaConvV1_LHC17*_fast-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f 8 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC17pq_fast-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root
		fi
	done

	rm $OUTPUTDIR/GammaConvV1_LHC17pq_wSDD-pass1_*.root
	ls $OUTPUTDIR/GammaConvV1_LHC17*_wSDD-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f 8 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC17pq_wSDD-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root
		fi
	done

	rm $OUTPUTDIR/GammaConvV1_LHC17pq_woSDD-pass1_*.root
	ls $OUTPUTDIR/GammaConvV1_LHC17*_woSDD-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f 8 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC17pq_woSDD-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root
		fi
	done

    rm $OUTPUTDIR/GammaConvV1_LHC17p_all-pass1_*.root
	ls $OUTPUTDIR/GammaConvV1_LHC17p_*-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f 8 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root ]  && [ -f $OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC17p_all-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root
		fi
	done

	rm $OUTPUTDIR/GammaConvV1_LHC17q_all-pass1_*.root
	ls $OUTPUTDIR/GammaConvV1_LHC17q_*-pass1_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f 8 | cut -d "_" -f 4 | cut -d "." -f1`
		echo $number
		if [ -f $OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root ]  && [ -f $OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root ] ; then
			hadd -f $OUTPUTDIR/GammaConvV1_LHC17q_all-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root $OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root
		fi
	done

fi

if [ $1 = "mergeAllData" ]; then

    OUTPUTDIR_LHC17p_fast=$BASEDIR/$TRAINDIR/LHC17p_fast
    OUTPUTDIR_LHC17p_wSDD=$BASEDIR/$TRAINDIR/LHC17p_wSDD
    OUTPUTDIR_LHC17p_woSDD=$BASEDIR/$TRAINDIR/LHC17p_woSDD
    mkdir -p $OUTPUTDIR_LHC17p_fast
    mkdir -p $OUTPUTDIR_LHC17p_wSDD
    mkdir -p $OUTPUTDIR_LHC17p_woSDD

    OUTPUTDIR_LHC17q_fast=$BASEDIR/$TRAINDIR/LHC17q_fast
    OUTPUTDIR_LHC17q_wSDD=$BASEDIR/$TRAINDIR/LHC17q_wSDD
    OUTPUTDIR_LHC17q_woSDD=$BASEDIR/$TRAINDIR/LHC17q_woSDD
    mkdir -p $OUTPUTDIR_LHC17q_fast
    mkdir -p $OUTPUTDIR_LHC17q_wSDD
    mkdir -p $OUTPUTDIR_LHC17q_woSDD

    fileName="GammaConvV1_400.root"
    runs=`cat runlists/runNumbersLHC17p_all.txt`
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17p_wSDD/$run/$fileName $OUTPUTDIR_LHC17p_wSDD/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17p_wSDD/$fileName $OUTPUTDIR_LHC17p_wSDD/intermediate.root $OUTPUTDIR_LHC17p_wSDD/$run/$fileName
#                 mv $OUTPUTDIR_LHC17p_wSDD/$fileName $OUTPUTDIR_LHC17p_wSDD/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17p_wSDD/intermediate.root $OUTPUTDIR_LHC17p_wSDD/$fileName
#     ls $OUTPUTDIR_LHC17p_wSDD/GammaConvV1_*.root > fileLHC17p_wSDD.txt
#     fileNumbers=`cat fileLHC17p_wSDD.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_wSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_wSDD_$number.log\"\)
#     done;

    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17p_woSDD/$run/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17p_woSDD/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17p_woSDD/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$fileName
    ls $OUTPUTDIR_LHC17p_woSDD/GammaConvV1_*.root > fileLHC17p_woSDD.txt
    fileNumbers=`cat fileLHC17p_woSDD.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_woSDD_$number.log\"\)
    done;

#     runs=`cat runlists/runNumbersLHC17p_fast.txt`
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17p_fast/$run/$fileName $OUTPUTDIR_LHC17p_fast/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17p_fast/$fileName $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$run/$fileName
#                 mv $OUTPUTDIR_LHC17p_fast/$fileName $OUTPUTDIR_LHC17p_fast/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17p_fast/intermediate.root $OUTPUTDIR_LHC17p_fast/$fileName
#     ls $OUTPUTDIR_LHC17p_fast/GammaConvV1_*.root > fileLHC17p_fast.txt
#     fileNumbers=`cat fileLHC17p_fast.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_fast_$number.log\"\)
#     done;
#
#
    runs=`cat runlists/runNumbersLHC17q_all.txt`
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17q_wSDD/$run/$fileName $OUTPUTDIR_LHC17q_wSDD/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17q_wSDD/$fileName $OUTPUTDIR_LHC17q_wSDD/intermediate.root $OUTPUTDIR_LHC17q_wSDD/$run/$fileName
#                 mv $OUTPUTDIR_LHC17q_wSDD/$fileName $OUTPUTDIR_LHC17q_wSDD/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17q_wSDD/intermediate.root $OUTPUTDIR_LHC17q_wSDD/$fileName
#     ls $OUTPUTDIR_LHC17q_wSDD/GammaConvV1_*.root > fileLHC17q_wSDD.txt
#     fileNumbers=`cat fileLHC17q_wSDD.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_wSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17q_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17q_wSDD_$number.log\"\)
#     done;
#
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17q_woSDD/$run/$fileName $OUTPUTDIR_LHC17q_woSDD/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17q_woSDD/$fileName $OUTPUTDIR_LHC17q_woSDD/intermediate.root $OUTPUTDIR_LHC17q_woSDD/$run/$fileName
#                 mv $OUTPUTDIR_LHC17q_woSDD/$fileName $OUTPUTDIR_LHC17q_woSDD/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17q_woSDD/intermediate.root $OUTPUTDIR_LHC17q_woSDD/$fileName
#     ls $OUTPUTDIR_LHC17q_woSDD/GammaConvV1_*.root > fileLHC17q_woSDD.txt
#     fileNumbers=`cat fileLHC17q_woSDD.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17q_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17q_woSDD_$number.log\"\)
#     done;
#

    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17q_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17q_fast/$run/$fileName $OUTPUTDIR_LHC17q_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17q_fast/$fileName $OUTPUTDIR_LHC17q_fast/intermediate.root $OUTPUTDIR_LHC17q_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17q_fast/$fileName $OUTPUTDIR_LHC17q_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17q_fast/intermediate.root $OUTPUTDIR_LHC17q_fast/$fileName
    ls $OUTPUTDIR_LHC17q_fast/GammaConvV1_*.root > fileLHC17q_fast.txt
    fileNumbers=`cat fileLHC17q_fast.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17q_fast_$number.log\"\)
    done;

fi



if [ $2 = "no" ]; then
   echo "Not dowloading MC";

elif [ $2 = "yes" ]; then

    OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/LHC17l3b_fast
    OUTPUTDIR_LHC17l3b_cent=$BASEDIR/$TRAINDIR/LHC17l3b_wSDD
    OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/LHC17l3b_woSDD
    mkdir -p $OUTPUTDIR_LHC17l3b_fast
    mkdir -p $OUTPUTDIR_LHC17l3b_cent
    mkdir -p $OUTPUTDIR_LHC17l3b_woSDD

    OUTPUTDIR_LHC17l4b_fast=$BASEDIR/$TRAINDIR/LHC17l4b_fast
    OUTPUTDIR_LHC17l4b_cent=$BASEDIR/$TRAINDIR/LHC17l4b_wSDD
    OUTPUTDIR_LHC17l4b_woSDD=$BASEDIR/$TRAINDIR/LHC17l4b_woSDD
    mkdir -p $OUTPUTDIR_LHC17l4b_fast
    mkdir -p $OUTPUTDIR_LHC17l4b_cent
    mkdir -p $OUTPUTDIR_LHC17l4b_woSDD

    if [ $3 = "runwise" ]; then

        ################################ LHC17l3b_fast
        rm runNumbersLHC17l3bfastNotMerged.txt
        echo $OUTPUTDIR_LHC17l3b_fast
        runs=`cat runlists/runNumbersLHC17l3b_fast.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run
            if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/$fileName file:$OUTPUTDIR_LHC17l3b_fast/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17l3bfastNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17l3bfastNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            rm $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17l3b_fast/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName $OUTPUTDIR_LHC17l3b_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17l3b_fast/$run/$fileName $OUTPUTDIR_LHC17l3b_fast/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17l3b_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
            fi
        done;


        ################################ LHC17l3b_cent
        rm runNumbersLHC17l3bwSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17l3b_cent
        runs=`cat runlists/runNumbersLHC17l3b_wSDD.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l3b_cent/$run
            if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/$fileName file:$OUTPUTDIR_LHC17l3b_cent/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17l3bwSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17l3bwSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            rm $OUTPUTDIR_LHC17l3b_cent/$run/$fileName
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17l3b_cent/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17l3b_cent/$run/$fileName $OUTPUTDIR_LHC17l3b_cent/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17l3b_cent/$run/$fileName $OUTPUTDIR_LHC17l3b_cent/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17l3b_cent/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_cent/$run/$fileName
            fi
        done;


        ################################ LHC17l3b_woSDD
        rm runNumbersLHC17l3bwoSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17l3b_woSDD
        runs=`cat runlists/runNumbersLHC17l3b_woSDD.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run
            if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l3b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/$fileName file:$OUTPUTDIR_LHC17l3b_woSDD/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17l3bwoSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17l3bwoSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l3b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            rm $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17l3b_woSDD/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17l3b_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
            fi
        done;


        ################################ LHC17l4b_fast
        rm runNumbersLHC17qfastNotMerged.txt
        echo $OUTPUTDIR_LHC17l4b_fast
        runs=`cat runlists/runNumbersLHC17l4b_fast.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l4b_fast/$run
            if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/$fileName file:$OUTPUTDIR_LHC17l4b_fast/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17qfastNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17qfastNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            rm $OUTPUTDIR_LHC17l4b_fast/$run/$fileName
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17l4b_fast/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17l4b_fast/$run/$fileName $OUTPUTDIR_LHC17l4b_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17l4b_fast/$run/$fileName $OUTPUTDIR_LHC17l4b_fast/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17l4b_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_fast/$run/$fileName
            fi
        done;


        ################################ LHC17l4b_cent
        rm runNumbersLHC17qwSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17l4b_cent
        runs=`cat runlists/runNumbersLHC17l4b_wSDD.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l4b_cent/$run
            if [ -f $OUTPUTDIR_LHC17l4b_cent/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l4b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/$fileName file:$OUTPUTDIR_LHC17l4b_cent/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17l4b_cent/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17qwSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17qwSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/sim/2017/LHC17l4b_cent/$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l4b_cent/$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            rm $OUTPUTDIR_LHC17l4b_cent/$run/$fileName
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17l4b_cent/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17l4b_cent/$run/$fileName $OUTPUTDIR_LHC17l4b_cent/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17l4b_cent/$run/$fileName $OUTPUTDIR_LHC17l4b_cent/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17l4b_cent/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_cent/$run/$fileName
            fi
        done;


        ################################ LHC17l4b_cent_woSDD
        rm runNumbersLHC17qwoSDDNotMerged.txt
        echo $OUTPUTDIR_LHC17l4b_woSDD
        runs=`cat runlists/runNumbersLHC17l4b_woSDD.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l4b_woSDD/$run
            if [ -f $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/$fileName file:$OUTPUTDIR_LHC17l4b_woSDD/$run/
            fi
            if [ -f $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                echo $run >> runNumbersLHC17qwoSDDNotMerged.txt
            fi
        done;

        NotMergedruns=`cat runNumbersLHC17qwoSDDNotMerged.txt`
        for run in $NotMergedruns; do
            echo "copying stage_1 output for " $run
            mkdir -p $run
            stageOutputs=`alien_ls /alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/Stage_1/`
            for stageOutput in $stageOutputs; do
                mkdir -p $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput
                if [ -f $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/Stage_1/$stageOutput/$fileName file:$OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/
                fi
            done;
            counter=0;
            rm $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName
            if [ $fileName != "AnalysisResults.root" ]; then
                for stageOutput in $stageOutputs; do
                    if [ -f $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/$fileName ]; then
                        if [ $counter = 0 ]; then
                            cp $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/$fileName $OUTPUTDIR_LHC17l4b_woSDD/$run/intermediate_$run.root
                            counter=$(($counter+1));
                            echo $counter;
                        else
                            hadd -f $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l4b_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/$fileName
                            mv $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l4b_woSDD/$run/intermediate_$run.root
                        fi
                    fi
                done;
                mv $OUTPUTDIR_LHC17l4b_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName
            fi
        done;


        if [ $fileName != "AnalysisResults.root" ]; then

            runs=`cat runlists/runNumbersLHC17l3b_fast.txt`
            for run in $runs; do
                ls $OUTPUTDIR_LHC17l3b_fast/$run/GammaConvV1_*.root > fileLHC17l3b_fast.txt
                fileNumbers=`cat fileLHC17l3b_fast.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_fast/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_fast/$run/GammaConvV1_LHC17l3b_fast_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l3b_fast/$run/GammaConvV1_LHC17l3b_fast_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_fast/$run/CutSelection_LHC17l3b_fast_$number.log\"\)
                done;
            done;

            runs=`cat runlists/runNumbersLHC17l3b_wSDD.txt`
            for run in $runs; do
                ls $OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_*.root > fileLHC17l3b_cent.txt
                fileNumbers=`cat fileLHC17l3b_cent.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_LHC17l3b_cent_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_LHC17l3b_cent_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_cent/$run/CutSelection_LHC17l3b_cent_$number.log\"\)
                done;
            done;

            runs=`cat runlists/runNumbersLHC17l3b_woSDD.txt`
            for run in $runs; do
                ls $OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_*.root > fileLHC17l3b_woSDD.txt
                fileNumbers=`cat fileLHC17l3b_woSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_LHC17l3b_woSDD_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_LHC17l3b_woSDD_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_woSDD/$run/CutSelection_LHC17l3b_woSDD_$number.log\"\)
                done;
            done;

            runs=`cat runlists/runNumbersLHC17l4b_fast.txt`
            for run in $runs; do
                ls $OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_*.root > fileLHC17l4b_fast.txt
                fileNumbers=`cat fileLHC17l4b_fast.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_LHC17l4b_fast_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_LHC17l4b_fast_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_fast/$run/CutSelection_LHC17l4b_fast_$number.log\"\)
                done;
            done;

            runs=`cat runlists/runNumbersLHC17l4b_wSDD.txt`
            for run in $runs; do
                ls $OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_*.root > fileLHC17l4b_wSDD.txt
                fileNumbers=`cat fileLHC17l4b_wSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_LHC17l4b_wSDD_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_LHC17l4b_wSDD_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_cent/$run/CutSelection_LHC17l4b_wSDD_$number.log\"\)
                done;
            done;

            runs=`cat runlists/runNumbersLHC17l4b_woSDD.txt`
            for run in $runs; do
                ls $OUTPUTDIR_LHC17l4b_woSDD/$run/GammaConvV1_*.root > fileLHC17l4b_woSDD.txt
                fileNumbers=`cat fileLHC17l4b_woSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_woSDD/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_woSDD/$run/GammaConvV1_LHC17l4b_woSDD_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l4b_woSDD/$run/GammaConvV1_LHC17l4b_woSDD_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_woSDD/$run/CutSelection_LHC17l4b_woSDD_$number.log\"\)
                done;
            done;
        fi
    fi


    if [ $LHC17l3b_fast != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/merge/GammaConvV1* file:$OUTPUTDIR_LHC17l3b_fast/
        ls $OUTPUTDIR_LHC17l3b_fast/GammaConvV1_*.root > fileLHC17l3b_fast.txt
        fileNumbers=`cat fileLHC17l3b_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_fast_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_fast_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l3b_fast_$number.log\"\)
        done;
    fi

    if [ $LHC17l3b_cent != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/merge/GammaConvV1* file:$OUTPUTDIR_LHC17l3b_cent/
        ls $OUTPUTDIR_LHC17l3b_cent/GammaConvV1_*.root > fileLHC17l3b_cent.txt
        fileNumbers=`cat fileLHC17l3b_cent.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_cent/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_cent_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_cent_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l3b_cent_$number.log\"\)
        done;
    fi

    if [ $LHC17l3b_woSDD != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/merge/GammaConvV1* file:$OUTPUTDIR_LHC17l3b_woSDD/
        ls $OUTPUTDIR_LHC17l3b_woSDD/GammaConvV1_*.root > fileLHC17l3b_woSDD.txt
        fileNumbers=`cat fileLHC17l3b_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l3b_woSDD_$number.log\"\)
        done;
    fi

    if [ $LHC17l4b_fast != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/merge/GammaConvV1* file:$OUTPUTDIR_LHC17l4b_fast/
        ls $OUTPUTDIR_LHC17l4b_fast/GammaConvV1_*.root > fileLHC17l4b_fast.txt
        fileNumbers=`cat fileLHC17l4b_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l4b_fast_$number.log\"\)
        done;
    fi

    if [ $LHC17l4b_cent != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/merge/GammaConvV1* file:$OUTPUTDIR_LHC17l4b_cent/
        ls $OUTPUTDIR_LHC17l4b_cent/GammaConvV1_*.root > fileLHC17l4b_cent.txt
        fileNumbers=`cat fileLHC17l4b_cent.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_cent/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_cent_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_cent_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l4b_cent_$number.log\"\)
        done;
    fi

    if [ $LHC17l4b_woSDD != "" ]; then
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/merge/GammaConvV1* file:$OUTPUTDIR_LHC17l4b_woSDD/
        ls $OUTPUTDIR_LHC17l4b_woSDD/GammaConvV1_*.root > fileLHC17l4b_woSDD.txt
        fileNumbers=`cat fileLHC17l4b_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_woSDD_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_woSDD_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l4b_woSDD_$number.log\"\)
        done;
    fi

fi

if [ $2 = "mergefills" ]; then

    OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/LHC17l3b_fast
    OUTPUTDIR_LHC17l3b_cent=$BASEDIR/$TRAINDIR/LHC17l3b_wSDD
    OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/LHC17l3b_woSDD
    mkdir -p $OUTPUTDIR_LHC17l3b_fast
    mkdir -p $OUTPUTDIR_LHC17l3b_cent
    mkdir -p $OUTPUTDIR_LHC17l3b_woSDD

    OUTPUTDIR_LHC17l4b_fast=$BASEDIR/$TRAINDIR/LHC17l4b_fast
    OUTPUTDIR_LHC17l4b_cent=$BASEDIR/$TRAINDIR/LHC17l4b_wSDD
    OUTPUTDIR_LHC17l4b_woSDD=$BASEDIR/$TRAINDIR/LHC17l4b_woSDD
    mkdir -p $OUTPUTDIR_LHC17l4b_fast
    mkdir -p $OUTPUTDIR_LHC17l4b_cent
    mkdir -p $OUTPUTDIR_LHC17l4b_woSDD

    fileName="GammaConvV1_400.root"
    fileName2="GammaConvV1_400"
    ######################################################################### fill 1
    runs=`cat runlists/runNumbersLHC17p_fill1.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill1.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill1.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill1.root
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_fast/$run/$fileName $OUTPUTDIR_LHC17l3b_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill1.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill1.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill1.root

    ######################################################################### fill 2
    runs=`cat runlists/runNumbersLHC17p_fill2.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill2.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill2.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill2.root
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_fast/$run/$fileName $OUTPUTDIR_LHC17l3b_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill2.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill2.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill2.root

    ######################################################################### fill 3
    runs=`cat runlists/runNumbersLHC17p_fill3.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill3.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill3.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill3.root
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_fast/$run/$fileName $OUTPUTDIR_LHC17l3b_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill3.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill3.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill3.root

    ######################################################################### fill 4
    runs=`cat runlists/runNumbersLHC17q_fast.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill4.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill4.root $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$fileName2\_fill4.root
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_fast/$run/$fileName $OUTPUTDIR_LHC17l3b_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill4.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill4.root $OUTPUTDIR_LHC17l3b_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$fileName2\_fill4.root
fi

if [ $2 = "material" ]; then

    OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/Material_LHC17l3b_fast_LowInt
    mkdir -p $OUTPUTDIR_LHC17l3b_fast
    if [ $LHC17l3b_fast != "" ]; then
        echo "Downloading " $LHC17l3b_fast " in " $OUTPUTDIR_LHC17l3b_fast
        runs=`cat runlists/runNumbersLHC17l3b_LowInt.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run
            if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/$fileName file:$OUTPUTDIR_LHC17l3b_fast/$run/
            fi
        done;

        hadd -f $OUTPUTDIR_LHC17l3b_fast/$fileName $OUTPUTDIR_LHC17l3b_fast/282008/$fileName $OUTPUTDIR_LHC17l3b_fast/282016/$fileName $OUTPUTDIR_LHC17l3b_fast/282021/$fileName $OUTPUTDIR_LHC17l3b_fast/282025/$fileName $OUTPUTDIR_LHC17l3b_fast/282031/$fileName

        ls $OUTPUTDIR_LHC17l3b_fast/GammaConv*Material*.root > fileLHC17l3b_fast.txt
        fileNumbers=`cat fileLHC17l3b_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17l3b_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17l3b_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi


    OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/Material_LHC17l3b_woSDD_LowInt
    mkdir -p $OUTPUTDIR_LHC17l3b_woSDD
    if [ $LHC17l3b_woSDD != "" ]; then
        echo "Downloading " $LHC17l3b_woSDD " in " $OUTPUTDIR_LHC17l3b_woSDD
        runs=`cat runlists/runNumbersLHC17l3b_LowInt.txt`
        for run in $runs; do
            echo $run
            mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run
            if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
                echo "file has already been copied"
            else
                alien_cp alien:/alice/sim/2017/LHC17l3b_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/$fileName file:$OUTPUTDIR_LHC17l3b_woSDD/$run/
            fi
        done;

        hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$fileName $OUTPUTDIR_LHC17l3b_woSDD/282008/$fileName $OUTPUTDIR_LHC17l3b_woSDD/282016/$fileName $OUTPUTDIR_LHC17l3b_woSDD/282021/$fileName $OUTPUTDIR_LHC17l3b_woSDD/282025/$fileName $OUTPUTDIR_LHC17l3b_woSDD/282031/$fileName

        ls $OUTPUTDIR_LHC17l3b_woSDD/GammaConv*Material*.root > fileLHC17l3b_woSDD.txt
        fileNumbers=`cat fileLHC17l3b_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17l3b_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

fi

if [ $2 = "mergeTriggers" ]; then

	rm $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_all*.root
	ls $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f 8 | cut -d "_" -f 5 | cut -d "." -f1`
    done;
    echo $number
    if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_fast_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_cent_$number.root ]  && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD_$number.root ] ; then
        hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_all_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_fast_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_cent_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD_$number.root
    fi

	rm $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_all_*.root
	ls $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_*.root > filesForMerging.txt
	filesForMerging=`cat filesForMerging.txt`
	for fileName in $filesForMerging; do
		echo $fileName
		number=`echo $fileName  | cut -d "/" -f 8 | cut -d "_" -f 5 | cut -d "." -f1`
    done
    echo $number
    if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_cent_$number.root ]  && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_woSDD_$number.root ] ; then
        hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_all_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_cent_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC17l4b_woSDD_$number.root
    fi

fi

if [ $2 = "jetjet" ]; then

    OUTPUTDIR_LHC18b8_fast=$BASEDIR/$TRAINDIR/LHC18b8_fast
    OUTPUTDIR_LHC18b8_cent=$BASEDIR/$TRAINDIR/LHC18b8_wSDD
    OUTPUTDIR_LHC18b8_woSDD=$BASEDIR/$TRAINDIR/LHC18b8_woSDD
    mkdir -p $OUTPUTDIR_LHC18b8_fast
    mkdir -p $OUTPUTDIR_LHC18b8_cent
    mkdir -p $OUTPUTDIR_LHC18b8_woSDD

            ################################ LHC18b8_fast
            rm runNumbersLHC18b8fastNotMerged.txt
            echo $OUTPUTDIR_LHC18b8_fast
            runs=`cat runlists/runNumbersLHC18b8_fast.txt`
            for run in $runs; do
                echo $run
                binNumbersJJ=`cat runlists/binsJetJetLHC18b8.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    mkdir -p $OUTPUTDIR_LHC18b8_fast/$binNumber/$run
                    if [ -f $OUTPUTDIR_LHC18b8_fast/$binNumber/$run/$fileName ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2018/LHC18b8_fast/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_fast/$fileName file:$OUTPUTDIR_LHC18b8_fast/$binNumber/$run
                    fi
                    if [ -f $OUTPUTDIR_LHC18b8_fast/$binNumber/$run/$fileName ]; then
                        echo "file has already been copied"
                    else
                        echo $run >> runNumbersLHC18b8fastNotMerged.txt
                    fi
#                 done;

#                 mkdir -p $OUTPUTDIR_LHC18b8_fast/$run/$binNumber
#                     if [ -f $OUTPUTDIR_LHC18b8_fast/$run/$binNumber/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         alien_cp alien:/alice/sim/2018/LHC18b8_fast/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_fast/$fileName file:$OUTPUTDIR_LHC18b8_fast/$run/$binNumber\_$fileName
#                     fi
#                     if [ -f $OUTPUTDIR_LHC18b8_fast/$run/$binNumber/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         echo $run >> runNumbersLHC18b8fastNotMerged.txt
#                     fi
                done;
            done;

            rm $OUTPUTDIR_LHC18b8_fast/GammaConv*.root
            firstrunNumber=`head -n1 runlists/runNumbersLHC18b8_fast.txt`
            firstbinNumber=`head -n1 runlists/binsJetJetLHC18b8.txt`
            ls $OUTPUTDIR_LHC18b8_fast/$firstbinNumber/$firstrunNumber/GammaConvV1_*.root > fileJJmergeLHC18b8.txt
            fileNumbers=`cat fileJJmergeLHC18b8.txt`
            MergeAccordingToSpecificRunlist fileJJmergeLHC18b8.txt $OUTPUTDIR_LHC18b8_fast $NSlashes3 GammaConvV1 fast runlists/runNumbersLHC18b8_fast.txt runlists/binsJetJetLHC18b8.txt
#
#             ################################ LHC18b8_cent
#             rm runNumbersLHC18b8wSDDNotMerged.txt
#             echo $OUTPUTDIR_LHC18b8_cent
#             runs=`cat runlists/runNumbersLHC18b8_wSDD.txt`
#             for run in $runs; do
#                 echo $run
#                 binNumbersJJ=`cat runlists/binsJetJetLHC18b8.txt`
#                 for binNumber in $binNumbersJJ; do
#                     echo $binNumber
#                     mkdir -p $OUTPUTDIR_LHC18b8_cent/$binNumber/$run
#                     if [ -f $OUTPUTDIR_LHC18b8_cent/$binNumber/$run/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                     alien_cp alien:/alice/sim/2018/LHC18b8_cent/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_cent/$fileName file:$OUTPUTDIR_LHC18b8_cent/$binNumber/$run
#                     fi
#                     if [ -f $OUTPUTDIR_LHC18b8_cent/$binNumber/$run/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         echo $run >> runNumbersLHC18b8fastNotMerged.txt
#                     fi
# #                 done;
# #             done;
#                     mkdir -p $OUTPUTDIR_LHC18b8_cent/$run/$binNumber
#                     if [ -f $OUTPUTDIR_LHC18b8_cent/$run/$binNumber/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                     alien_cp alien:/alice/sim/2018/LHC18b8_cent/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_cent/$fileName file:$OUTPUTDIR_LHC18b8_cent/$run/$binNumber\_$fileName
#                     fi
#                     if [ -f $OUTPUTDIR_LHC18b8_cent/$run/$binNumber/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         echo $run >> runNumbersLHC18b8fastNotMerged.txt
#                     fi
#                 done;
#             done;
#
#             rm $OUTPUTDIR_LHC18b8_cent/GammaConv*.root
#             firstrunNumber=`head -n1 runlists/runNumbersLHC18b8_wSDD.txt`
#             firstbinNumber=`head -n1 runlists/binsJetJetLHC18b8.txt`
#             ls $OUTPUTDIR_LHC18b8_cent/$firstbinNumber/$firstrunNumber/GammaConvV1_*.root > fileJJmergeLHC18b8.txt
#             fileNumbers=`cat fileJJmergeLHC18b8.txt`
#             MergeAccordingToSpecificRunlist fileJJmergeLHC18b8.txt $OUTPUTDIR_LHC18b8_cent $NSlashes3 GammaConvV1 wSDD runlists/runNumbersLHC18b8_wSDD.txt runlists/binsJetJetLHC18b8.txt
#
#
#             ################################ LHC18b8_woSDD
#             rm runNumbersLHC18b8woSDDNotMerged.txt
#             echo $OUTPUTDIR_LHC18b8_woSDD
#             runs=`cat runlists/runNumbersLHC18b8_woSDD.txt`
#             for run in $runs; do
#                 echo $run
#                 binNumbersJJ=`cat runlists/binsJetJetLHC18b8.txt`
#                 for binNumber in $binNumbersJJ; do
#                     echo $binNumber
#                     mkdir -p $OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run
#                     if [ -f $OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         alien_cp alien:/alice/sim/2018/LHC18b8_woSDD/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_woSDD/$fileName file:$OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run
#                     fi
#                     if [ -f $OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         echo $run >> runNumbersLHC18b8fastNotMerged.txt
#                     fi
# #                 done;
# #             done;
#
#                     mkdir -p $OUTPUTDIR_LHC18b8_woSDD/$run/$binNumber
#                     if [ -f $OUTPUTDIR_LHC18b8_woSDD/$run/$binNumber/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         alien_cp alien:/alice/sim/2018/LHC18b8_cent_woSDD/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_woSDD/$fileName file:$OUTPUTDIR_LHC18b8_woSDD/$run/$binNumber\_$fileName
#                     fi
#                     if [ -f $OUTPUTDIR_LHC18b8_woSDD/$run/$binNumber/$fileName ]; then
#                         echo "file has already been copied"
#                     else
#                         echo $run >> runNumbersLHC18b8fastNotMerged.txt
#                     fi
#                 done;
#             done;
#
#             rm $OUTPUTDIR_LHC18b8_woSDD/GammaConv*.root
#             firstrunNumber=`head -n1 runlists/runNumbersLHC18b8_woSDD.txt`
#             firstbinNumber=`head -n1 runlists/binsJetJetLHC18b8.txt`
#             ls $OUTPUTDIR_LHC18b8_woSDD/$firstbinNumber/$firstrunNumber/GammaConvV1_*.root > fileJJmergeLHC18b8.txt
#             fileNumbers=`cat fileJJmergeLHC18b8.txt`
#             MergeAccordingToSpecificRunlist fileJJmergeLHC18b8.txt $OUTPUTDIR_LHC18b8_woSDD $NSlashes3 GammaConvV1 woSDD runlists/runNumbersLHC18b8_woSDD.txt runlists/binsJetJetLHC18b8.txt

fi

if [ $2 = "mergeAllMC" ]; then

    fileName="GammaConvV1_400.root"
    OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/LHC17l3b_fast
    OUTPUTDIR_LHC17l3b_cent=$BASEDIR/$TRAINDIR/LHC17l3b_wSDD
    OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/LHC17l3b_woSDD
    mkdir -p $OUTPUTDIR_LHC17l3b_fast
    mkdir -p $OUTPUTDIR_LHC17l3b_cent
    mkdir -p $OUTPUTDIR_LHC17l3b_woSDD

    OUTPUTDIR_LHC17l4b_fast=$BASEDIR/$TRAINDIR/LHC17l4b_fast
    OUTPUTDIR_LHC17l4b_cent=$BASEDIR/$TRAINDIR/LHC17l4b_wSDD
    OUTPUTDIR_LHC17l4b_woSDD=$BASEDIR/$TRAINDIR/LHC17l4b_woSDD
    mkdir -p $OUTPUTDIR_LHC17l4b_fast
    mkdir -p $OUTPUTDIR_LHC17l4b_cent
    mkdir -p $OUTPUTDIR_LHC17l4b_woSDD
#
#     runs=`cat runlists/runNumbersLHC17l3b_wSDD.txt`
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17l3b_cent/$run/$fileName $OUTPUTDIR_LHC17l3b_cent/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17l3b_cent/$fileName $OUTPUTDIR_LHC17l3b_cent/intermediate.root $OUTPUTDIR_LHC17l3b_cent/$run/$fileName
#                 mv $OUTPUTDIR_LHC17l3b_cent/$fileName $OUTPUTDIR_LHC17l3b_cent/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17l3b_cent/intermediate.root $OUTPUTDIR_LHC17l3b_cent/$fileName
#     ls $OUTPUTDIR_LHC17l3b_cent/GammaConvV1_*.root > fileLHC17l3b_cent.txt
#     fileNumbers=`cat fileLHC17l3b_cent.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_cent/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l3b_wSDD_$number.log\"\)
#     done;
#
#
#     runs=`cat runlists/runNumbersLHC17l3b_woSDD.txt`
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
#                 mv $OUTPUTDIR_LHC17l3b_woSDD/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$fileName
#     ls $OUTPUTDIR_LHC17l3b_woSDD/GammaConvV1_*.root > fileLHC17l3b_woSDD.txt
#     fileNumbers=`cat fileLHC17l3b_woSDD.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l3b_woSDD_$number.log\"\)
#     done;
#
#     counter=0;
#     runs=`cat runlists/runNumbersLHC17l3b_fast.txt`
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17l3b_fast/$run/$fileName $OUTPUTDIR_LHC17l3b_fast/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17l3b_fast/$fileName $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
#                 mv $OUTPUTDIR_LHC17l3b_fast/$fileName $OUTPUTDIR_LHC17l3b_fast/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17l3b_fast/intermediate.root $OUTPUTDIR_LHC17l3b_fast/$fileName
#     ls $OUTPUTDIR_LHC17l3b_fast/GammaConvV1_*.root > fileLHC17l3b_fast.txt
#     fileNumbers=`cat fileLHC17l3b_fast.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l3b_fast_$number.log\"\)
#     done;
#
#
#     counter=0;
#     runs=`cat runlists/runNumbersLHC17l4b_wSDD.txt`
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17l4b_cent/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17l4b_cent/$run/$fileName $OUTPUTDIR_LHC17l4b_cent/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17l4b_cent/$fileName $OUTPUTDIR_LHC17l4b_cent/intermediate.root $OUTPUTDIR_LHC17l4b_cent/$run/$fileName
#                 mv $OUTPUTDIR_LHC17l4b_cent/$fileName $OUTPUTDIR_LHC17l4b_cent/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17l4b_cent/intermediate.root $OUTPUTDIR_LHC17l4b_cent/$fileName
#     ls $OUTPUTDIR_LHC17l4b_cent/GammaConvV1_*.root > fileLHC17l4b_cent.txt
#     fileNumbers=`cat fileLHC17l4b_cent.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_cent/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_wSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_wSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l4b_wSDD_$number.log\"\)
#     done;
#
#     counter=0;
#     runs=`cat runlists/runNumbersLHC17l4b_woSDD.txt`
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l4b_woSDD/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17l4b_woSDD/$fileName $OUTPUTDIR_LHC17l4b_woSDD/intermediate.root $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName
#                 mv $OUTPUTDIR_LHC17l4b_woSDD/$fileName $OUTPUTDIR_LHC17l4b_woSDD/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17l4b_woSDD/intermediate.root $OUTPUTDIR_LHC17l4b_woSDD/$fileName
#     ls $OUTPUTDIR_LHC17l4b_woSDD/GammaConvV1_*.root > fileLHC17l4b_woSDD.txt
#     fileNumbers=`cat fileLHC17l4b_woSDD.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l4b_woSDD_$number.log\"\)
#     done;

    counter=0;
    runs=`cat runlists/runNumbersLHC17l4b_fast.txt`
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l4b_fast/$run/$fileName $OUTPUTDIR_LHC17l4b_fast/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l4b_fast/$fileName $OUTPUTDIR_LHC17l4b_fast/intermediate.root $OUTPUTDIR_LHC17l4b_fast/$run/$fileName
                mv $OUTPUTDIR_LHC17l4b_fast/$fileName $OUTPUTDIR_LHC17l4b_fast/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l4b_fast/intermediate.root $OUTPUTDIR_LHC17l4b_fast/$fileName
    ls $OUTPUTDIR_LHC17l4b_fast/GammaConvV1_*.root > fileLHC17l4b_fast.txt
    fileNumbers=`cat fileLHC17l4b_fast.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l4b_fast_$number.log\"\)
    done;

fi
