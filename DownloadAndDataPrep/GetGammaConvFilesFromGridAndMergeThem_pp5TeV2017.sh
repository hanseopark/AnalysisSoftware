#! /bin/bash
source basicFunction.sh

# download script for pp 5TeV from 2017
BASEDIR=/home/admin1/leardini/GridOutput/pp
mkdir -p $BASEDIR
TRAINDIR=Legotrain-vAN-20181127-1_increasedMCStatQA;

### Data
LHC17pqMETA=2532_20181127-1634;
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
LHC17lMETA=1161_20181130-1442; #997_20180831-1159;
LHC17l3b_fast=$LHC17lMETA\_child_1
LHC17l3b_woSDD=$LHC17lMETA\_child_2
LHC17l3b_cent=;
LHC17l4b_fast=;
LHC17l4b_woSDD=;
LHC17l4b_cent=;
LHC18j2_fast=$LHC17lMETA\_child_3;
LHC18j2_woSDD=$LHC17lMETA\_child_4;

# LHC17l3b_fast=$LHC17lMETA\_child_1
# LHC17l3b_cent=$LHC17lMETA\_child_2
# LHC17l3b_woSDD=$LHC17lMETA\_child_3
# LHC17l4b_fast=$LHC17lMETA\_child_4
# LHC17l4b_cent=$LHC17lMETA\_child_5
# LHC17l4b_woSDD=$LHC17lMETA\_child_6

########### additional Pythia ########
LHC18d6cMETA=;
LHC18d6c_fast=; #$LHC18d6cMETA\_child_1;
LHC18d6c_woSDD=; #$LHC18d6cMETA\_child_2;
LHC18d6c_woSDD_part2=; #$LHC18d6cMETA\_child_3;
LHC18d6c_wSDD=; #$LHC18d6cMETA\_child_4;

###### all low intensity pythia ######
LowIntPythiaMETA=; #3619_20181127-1812; #3477_20180823-1142;
LHC17b3lLowInt_fast=$LowIntPythiaMETA\_child_1;
LHC17b3lLowInt_woSDD=$LowIntPythiaMETA\_child_2;
LHC17b3lLowInt_wSDD=$LowIntPythiaMETA\_child_3;
LHC18d6cLowInt_fast=$LowIntPythiaMETA\_child_4;
LHC18d6cLowInt_woSDD=$LowIntPythiaMETA\_child_5;
LHC18d6cLowInt_woSDDpart2=$LowIntPythiaMETA\_child_6;
LHC18d6cLowInt_wSDD=$LowIntPythiaMETA\_child_7;


############## Phojet ################
LHC18d6bMETA=; #3622_20181127-1813; #3477_20180823-1142;
LHC18d6b_fast=$LHC18d6bMETA\_child_1;
LHC18d6b_woSDD=$LHC18d6bMETA\_child_2;
LHC18d6b_fastpart2=$LHC18d6bMETA\_child_3;
LHC18d6b_woSDDpart2=$LHC18d6bMETA\_child_4;
LHC18d6b_wSDD=;


### MC Jet Jet
LHC18b8META=; #974_20180817-1718
LHC18b8_fast=$LHC18b8META\_child_1
LHC18b8_cent=$LHC18b8META\_child_2
LHC18b8_woSDD=$LHC18b8META\_child_3

###################### Data 2015 #####################
LHC15n_pass4=;
LHC16k5a_pass3=;
LHC16k5b_pass3=;
LHC17e2_pass4=;
LHC16h3_JJ_pass4=;
######################################################

OUTPUTDIR=$BASEDIR/$TRAINDIR
TRAINPATHData=GA_pp #_AOD
AODflag=;
if [ $2 = "jetjet" ]; then
    TRAINPATHMC=GA_pp_MC_AOD
if [ $2 = "AOD" ]; then
    TRAINPATHMC=GA_pp_MC_AOD
    AODflag=AOD209/
else
    TRAINPATHMC=GA_pp_MC
fi

NSlashes=10 #without mergelist = 9
NSlashes2=10
NSlashes3=11

mergeFolder=merge_runlist_1;

#QA 2017:
#data 2454_20180815-1616
#phojet onfly 3453_20180816-1714
#phojet offline 3455_20180816-1718
#pythia offline 3454_20180816-1717
#pythia onfly 3449_20180815-1258

#QA 2015:
#data 2455_20180815-1305
#MC 16k5a 3444_20180815-1258
#16k5b
#17e2 3446_20180815-1252


fileToDownloadQAPhoton=AnalysisResults.root;
fileToDownloadQAEvent=GammaConvV1_400.root;
# fileName=GammaConvV1_406.root;

# fileToDownload=GammaConv_Material_121.root;

if [ $1 = "no" ]; then
   echo "Not dowloading data";
elif [ $1 = "yes" ]; then

    if [ $3 = "runwise" ]; then

        ################################ LHC17p_fast
        if [ $LHC17p_fast != "" ]; then
            rm runNumbersLHC17pfastNotMergedQAP.txt
            rm runNumbersLHC17pfastNotMergedQAE.txt
            echo $OUTPUTDIR_LHC17p_fast
            runs=`cat runlists/runNumbersLHC17p_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_fast/$run
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17p_fast/$run/
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17p_fast/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileToDownloadQAEvent ]; then
                echo "file " $OUTPUTDIR_LHC17p_fast/$run/$fileToDownloadQAEvent "has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17pfastNotMergedQAE.txt
            fi
            NotMergedruns=`cat runNumbersLHC17pfastNotMergedQAE.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                if [ $fileToDownloadQAEvent != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17p_fast.txt
                        fileNames =`cat fileStagedToMergeLHC17p_fast.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17p_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_fast/$run/$fileName
                fi
            done;


            if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileToDownloadQAPhoton ]; then
                echo "file " $OUTPUTDIR_LHC17p_fast/$run/$fileToDownloadQAPhoton "has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17pfastNotMergedQAP.txt
            fi
            NotMergedruns=`cat runNumbersLHC17pfastNotMergedQAP.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/Stage_1/$stageOutput/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/
                done;
            done;
        fi

        ################################ LHC17p_wSDD
        if [ $LHC17p_wSDD != "" ]; then
            rm runNumbersLHC17pwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17p_wSDD
            runs=`cat runlists/runNumbersLHC17p_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_wSDD/$run
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/$fileToDownload file:$OUTPUTDIR_LHC17p_wSDD/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/$fileToDownload ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17pwSDDNotMerged.txt
            fi

            NotMergedruns=`cat runNumbersLHC17pwSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/Stage_1/$stageOutput/$fileToDownload file:$OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17p_wSDD.txt
                        fileNames =`cat fileStagedToMergeLHC17p_wSDD.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17p_wSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_wSDD/$run/$fileName
                fi
            done;
        fi

        ################################ LHC17p_woSDD
        if [ $LHC17p_woSDD != "" ]; then
            rm runNumbersLHC17pwoSDDNotMergedQAP.txt
            rm runNumbersLHC17pwoSDDNotMergedQAE.txt
            echo $OUTPUTDIR_LHC17p_woSDD
            runs=`cat runlists/runNumbersLHC17p_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17p_woSDD/$run/
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17p_woSDD/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileToDownloadQAEvent ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17pwoSDDNotMergedQAE.txt
            fi
            NotMergedruns=`cat runNumbersLHC17pwoSDDNotMergedQAE.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                if [ $fileToDownloadQAEvent != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17p_woSDD.txt
                        fileNames =`cat fileStagedToMergeLHC17p_woSDD.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17p_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17p_woSDD/$run/$fileName
                fi
            done;

            if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileToDownloadQAPhoton ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17pwoSDDNotMergedQAP.txt
            fi
            NotMergedruns=`cat runNumbersLHC17pwoSDDNotMergedQAP.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/Stage_1/$stageOutput/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/
                done;
            done;
        fi

        ################################ LHC17q_fast
        if [ $LHC17q_fast != "" ]; then
            rm runNumbersLHC17qfastNotMergedQAE.txt
            rm runNumbersLHC17qfastNotMergedQAP.txt
            echo $OUTPUTDIR_LHC17q_fast
            runs=`cat runlists/runNumbersLHC17q_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17q_fast/$run
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17q_fast/$run/
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17q_fast/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17q_fast/$run/$fileToDownloadQAEvent ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17qfastNotMergedQAE.txt
            fi
            NotMergedruns=`cat runNumbersLHC17qfastNotMergedQAE.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                if [ $fileToDownloadQAEvent != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17q_fast.txt
                        fileNames =`cat fileStagedToMergeLHC17q_fast.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17q_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_fast/$run/$fileName
                fi
            done;

            if [ -f $OUTPUTDIR_LHC17q_fast/$run/$fileToDownloadQAPhoton ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17qfastNotMergedQAP.txt
            fi
            NotMergedruns=`cat runNumbersLHC17qfastNotMergedQAP.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/Stage_1/$stageOutput/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/
                done;
            done;

        fi

        ################################ LHC17q_wSDD
        if [ $LHC17q_wSDD != "" ]; then
            rm runNumbersLHC17qwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17q_wSDD
            runs=`cat runlists/runNumbersLHC17q_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17q_wSDD/$run
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/$fileToDownload file:$OUTPUTDIR_LHC17q_wSDD/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/$fileToDownload ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17qwSDDNotMerged.txt
            fi

            NotMergedruns=`cat runNumbersLHC17qwSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/Stage_1/$stageOutput/$fileToDownload file:$OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17q_wSDD.txt
                        fileNames =`cat fileStagedToMergeLHC17q_wSDD.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17q_wSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_wSDD/$run/$fileName
                fi
            done;
        fi

        ################################ LHC17q_woSDD
        if [ $LHC17q_woSDD != "" ]; then
            rm runNumbersLHC17qwoSDDNotMergedQAE.txt
            rm runNumbersLHC17qwoSDDNotMergedQAP.txt
            echo $OUTPUTDIR_LHC17q_woSDD
            runs=`cat runlists/runNumbersLHC17q_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17q_woSDD/$run
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17q_woSDD/$run/
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17q_woSDD/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/$fileToDownloadQAEvent ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17qwoSDDNotMergedQAE.txt
            fi
            NotMergedruns=`cat runNumbersLHC17qwoSDDNotMergedQAE.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                if [ $fileToDownloadQAEvent != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17q_woSDD.txt
                        fileNames =`cat fileStagedToMergeLHC17q_woSDD.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17q_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17q_woSDD/$run/$fileName
                fi
            done;

            if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/$fileToDownloadQAPhoton ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17qwoSDDNotMergedQAP.txt
            fi
            NotMergedruns=`cat runNumbersLHC17qwoSDDNotMergedQAP.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/Stage_1/$stageOutput/$fileToDownloadQAPhoton file:$OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/
                done;
            done;

        fi

        if [ $fileToDownloadQAEvent != "AnalysisResults.root" ]; then
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

    if [ $3 = "QA" ]; then

        ################################ LHC17p_fast
        if [ $LHC17p_fast != "" ]; then
            rm runNumbersLHC17pfastNotMerged.txt
            echo $OUTPUTDIR_LHC17p_fast
            runs=`cat runlists/runNumbersLHC17p_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_fast/$run
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/AnalysisResults.root file:$OUTPUTDIR_LHC17p_fast/$run/
                if [ -f $OUTPUTDIR_LHC17p_fast/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/Stage_1/$stageOutput/AnalysisResults.root file:$OUTPUTDIR_LHC17p_fast/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17p_wSDD
        if [ $LHC17p_wSDD != "" ]; then
            rm runNumbersLHC17pwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17p_wSDD
            runs=`cat runlists/runNumbersLHC17p_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_wSDD/$run
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/AnalysisResults.root file:$OUTPUTDIR_LHC17p_wSDD/$run/
                if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/Stage_1/$stageOutput/AnalysisResults.root file:$OUTPUTDIR_LHC17p_wSDD/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17p_woSDD
        if [ $LHC17p_woSDD != "" ]; then
            rm runNumbersLHC17pwoSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17p_woSDD
            runs=`cat runlists/runNumbersLHC17p_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run
                alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/AnalysisResults.root file:$OUTPUTDIR_LHC17p_woSDD/$run/
                if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17p_woSDD/Stage_1/$stageOutput/AnalysisResults.root file:$OUTPUTDIR_LHC17p_woSDD/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17q_fast
        if [ $LHC17q_fast != "" ]; then
            rm runNumbersLHC17qfastNotMerged.txt
            echo $OUTPUTDIR_LHC17q_fast
            runs=`cat runlists/runNumbersLHC17q_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17q_fast/$run
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/AnalysisResults.root file:$OUTPUTDIR_LHC17q_fast/$run/
                if [ -f $OUTPUTDIR_LHC17q_fast/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17q_fast/Stage_1/$stageOutput/AnalysisResults.root file:$OUTPUTDIR_LHC17q_fast/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17q_wSDD
        if [ $LHC17q_wSDD != "" ]; then
            rm runNumbersLHC17qwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17q_wSDD
            runs=`cat runlists/runNumbersLHC17q_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17q_wSDD/$run
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/AnalysisResults.root file:$OUTPUTDIR_LHC17q_wSDD/$run/
                if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17q_wSDD/Stage_1/$stageOutput/AnalysisResults.root file:$OUTPUTDIR_LHC17q_wSDD/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17q_woSDD
        if [ $LHC17q_woSDD != "" ]; then
            rm runNumbersLHC17qwoSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17q_woSDD
            runs=`cat runlists/runNumbersLHC17q_all.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17q_woSDD/$run
                alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/AnalysisResults.root file:$OUTPUTDIR_LHC17q_woSDD/$run/
                if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/data/2017/LHC17q/000$run/pass1_CENT_woSDD/PWGGA/$TRAINPATHData/$LHC17q_woSDD/Stage_1/$stageOutput/AnalysisResults.root file:$OUTPUTDIR_LHC17q_woSDD/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi
    fi


    if [ $LHC17p_fast != "" ]; then
        OUTPUTDIR_LHC17p_fast=$BASEDIR/$TRAINDIR/LHC17p_fast/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17p_fast
        echo "Downloading " $LHC17p_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_fast/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17p_fast/
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
        OUTPUTDIR_LHC17p_wSDD=$BASEDIR/$TRAINDIR/LHC17p_wSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17p_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_wSDD/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17p_wSDD/
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
        OUTPUTDIR_LHC17p_woSDD=$BASEDIR/$TRAINDIR/LHC17p_woSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17p_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_woSDD/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17p_woSDD/
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
        OUTPUTDIR_LHC17q_fast=$BASEDIR/$TRAINDIR/LHC17q_fast/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17q_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17q_fast/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17q_fast/
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
        OUTPUTDIR_LHC17q_wSDD=$BASEDIR/$TRAINDIR/LHC17q_wSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17q_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17q_wSDD/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17q_wSDD/
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
        OUTPUTDIR_LHC17q_woSDD=$BASEDIR/$TRAINDIR/LHC17q_woSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17q_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17q_woSDD/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17q_woSDD/
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

    if [ $LHC17p_fast != "" ]; then
        OUTPUTDIR_LHC17p_fast=$BASEDIR/$TRAINDIR/Material_LHC17p_fast_LowInt
        mkdir -p $OUTPUTDIR_LHC17p_fast
        echo "Downloading " $LHC17p_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_fast/merge_runlist_3/GammaConv_Material_* file:$OUTPUTDIR_LHC17p_fast/
        ls $OUTPUTDIR_LHC17p_fast/GammaConv*Material*.root > fileLHC17p_fast.txt
        fileNumbers=`cat fileLHC17p_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17p_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17p_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC17p_woSDD != "" ]; then
        OUTPUTDIR_LHC17p_woSDD=$BASEDIR/$TRAINDIR/Material_LHC17p_woSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC17p_woSDD
        echo "Downloading " $LHC17p_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_woSDD/merge_runlist_3/GammaConv_Material_*  file:$OUTPUTDIR_LHC17p_woSDD/
        ls $OUTPUTDIR_LHC17p_woSDD/GammaConv*Material*.root > fileLHC17p_woSDD.txt
        fileNumbers=`cat fileLHC17p_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17p_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17p_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC17p_wSDD != "" ]; then
        OUTPUTDIR_LHC17p_wSDD=$BASEDIR/$TRAINDIR/Material_LHC17p_wSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC17p_wSDD
        echo "Downloading " $LHC17p_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC17p_wSDD/merge_runlist_3/GammaConv_Material_* file:$OUTPUTDIR_LHC17p_wSDD/
        ls $OUTPUTDIR_LHC17p_wSDD/GammaConv*Material*.root > fileLHC17p_wSDD.txt
        fileNumbers=`cat fileLHC17p_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17p_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17p_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $3 = "runwise" ]; then

        if [ $LHC17p_fast != "" ]; then
            echo "Downloading " $LHC17p_fast " in " $OUTPUTDIR_LHC17p_fast
            runs=`cat runlists/runNumbersLHC17p_LowInt.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_fast/$run
                if [ -f $OUTPUTDIR_LHC17p_fast/$run/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_fast/$fileToDownload file:$OUTPUTDIR_LHC17p_fast/$run/
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

        if [ $LHC17p_woSDD != "" ]; then
            echo "Downloading " $LHC17p_woSDD " in " $OUTPUTDIR_LHC17p_woSDD
            runs=`cat runlists/runNumbersLHC17p_LowInt.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_woSDD/$run
                if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_FAST/PWGGA/$TRAINPATHData/$LHC17p_woSDD/$fileToDownload file:$OUTPUTDIR_LHC17p_woSDD/$run/
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

        if [ $LHC17p_wSDD != "" ]; then
            echo "Downloading " $LHC17p_wSDD " in " $OUTPUTDIR_LHC17p_wSDD
            runs=`cat runlists/runNumbersLHC17p_LowInt.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17p_wSDD/$run
                if [ -f $OUTPUTDIR_LHC17p_wSDD/$run/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/data/2017/LHC17p/000$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHData/$LHC17p_wSDD/$fileToDownload file:$OUTPUTDIR_LHC17p_wSDD/$run/
                fi
            done;

            hadd -f $OUTPUTDIR_LHC17p_wSDD/$fileName $OUTPUTDIR_LHC17p_wSDD/282008/$fileName $OUTPUTDIR_LHC17p_wSDD/282016/$fileName $OUTPUTDIR_LHC17p_wSDD/282021/$fileName $OUTPUTDIR_LHC17p_wSDD/282025/$fileName $OUTPUTDIR_LHC17p_wSDD/282030/$fileName $OUTPUTDIR_LHC17p_wSDD/282031/$fileName

            ls $OUTPUTDIR_LHC17p_wSDD/GammaConv*Material*.root > fileLHC17p_wSDD.txt
            fileNumbers=`cat fileLHC17p_wSDD.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
                echo $number
                root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17p_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17p_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
            done;
        fi
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

    fileName="GammaConvV1_400.root"
    OUTPUTDIR_LHC17pq_woSDD=$BASEDIR/$TRAINDIR/LHC17pq_woSDD
    echo $OUTPUTDIR_LHC17pq_woSDD
    runs=`cat runlists/runNumbersLHC17p_all.txt`
    counter=0;
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17pq_woSDD/$run/$fileName ]; then
            echo `ls $OUTPUTDIR_LHC17pq_woSDD/$run/$fileName`
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17pq_woSDD/$run/$fileName $OUTPUTDIR_LHC17pq_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17pq_woSDD/$fileName $OUTPUTDIR_LHC17pq_woSDD/intermediate.root $OUTPUTDIR_LHC17pq_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17pq_woSDD/$fileName $OUTPUTDIR_LHC17pq_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17pq_woSDD/intermediate.root $OUTPUTDIR_LHC17pq_woSDD/$fileName
    ls $OUTPUTDIR_LHC17pq_woSDD/GammaConvV1_*.root > fileLHC17pq_woSDD.txt
    fileNumbers=`cat fileLHC17p_woSDD.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17pq_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_woSDD_$number.log\"\)
    done;


#     OUTPUTDIR_LHC17p_fast=$BASEDIR/$TRAINDIR/LHC17p_fast
#     OUTPUTDIR_LHC17p_wSDD=$BASEDIR/$TRAINDIR/LHC17p_wSDD
#     OUTPUTDIR_LHC17p_woSDD=$BASEDIR/$TRAINDIR/LHC17p_woSDD
#     mkdir -p $OUTPUTDIR_LHC17p_fast
#     mkdir -p $OUTPUTDIR_LHC17p_wSDD
#     mkdir -p $OUTPUTDIR_LHC17p_woSDD

#     OUTPUTDIR_LHC17q_fast=$BASEDIR/$TRAINDIR/LHC17q_fast
#     OUTPUTDIR_LHC17q_wSDD=$BASEDIR/$TRAINDIR/LHC17q_wSDD
#     OUTPUTDIR_LHC17q_woSDD=$BASEDIR/$TRAINDIR/LHC17q_woSDD
#     mkdir -p $OUTPUTDIR_LHC17q_fast
#     mkdir -p $OUTPUTDIR_LHC17q_wSDD
#     mkdir -p $OUTPUTDIR_LHC17q_woSDD


#
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
#     runs=`cat runlists/runNumbersLHC17p_all.txt`
#     counter=0;
#     fileName="GammaConvV1_408.root"
#
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
#
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17p_woSDD/$run/$fileName ]; then
#             echo `ls $OUTPUTDIR_LHC17p_woSDD/$run/$fileName`
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17p_woSDD/$run/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17p_woSDD/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$run/$fileName
#                 mv $OUTPUTDIR_LHC17p_woSDD/$fileName $OUTPUTDIR_LHC17p_woSDD/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17p_woSDD/intermediate.root $OUTPUTDIR_LHC17p_woSDD/$fileName
#     ls $OUTPUTDIR_LHC17p_woSDD/GammaConvV1_*.root > fileLHC17p_woSDD.txt
#     fileNumbers=`cat fileLHC17p_woSDD.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17p_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17p_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17p_woSDD_$number.log\"\)
#     done;

#     runs=`cat runlists/runNumbersLHC17q_all.txt`

#     fileName="GammaConvV1_408.root"
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

#     fileName="GammaConvV1_408.root"
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
#
#     fileName="GammaConvV1_408.root"
#     counter=0;
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17q_fast/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17q_fast/$run/$fileName $OUTPUTDIR_LHC17q_fast/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17q_fast/$fileName $OUTPUTDIR_LHC17q_fast/intermediate.root $OUTPUTDIR_LHC17q_fast/$run/$fileName
#                 mv $OUTPUTDIR_LHC17q_fast/$fileName $OUTPUTDIR_LHC17q_fast/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17q_fast/intermediate.root $OUTPUTDIR_LHC17q_fast/$fileName
#     ls $OUTPUTDIR_LHC17q_fast/GammaConvV1_*.root > fileLHC17q_fast.txt
#     fileNumbers=`cat fileLHC17q_fast.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17q_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC17q_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC17q_fast_$number.log\"\)
#     done;

fi



if [ $2 = "no" ]; then
   echo "Not dowloading MC";

elif [ $2 = "yes" ]; then

    if [ $3 = "runwise" ]; then

        ################################ LHC17l3b_fast
        if [ $LHC17l3b_fast != "" ]; then
            rm runNumbersLHC17l3bfastNotMerged.txt
            echo $OUTPUTDIR_LHC17l3b_fast
            runs=`cat runlists/runNumbersLHC17l3b_fast.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run
                alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_fast/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileToDownloadQAEvent ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17l3bfastNotMerged.txt
            fi

            NotMergedruns=`cat runNumbersLHC17l3bfastNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                rm $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17l3b_fast.txt
                        fileNames =`cat fileStagedToMergeLHC17l3b_fast.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17l3b_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_fast/$run/$fileName
                fi
            done;
        fi

        ################################ LHC17l3b_cent
        if [ $LHC17l3b_cent != "" ]; then
            rm runNumbersLHC17l3bwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l3b_cent
            runs=`cat runlists/runNumbersLHC17l3b_wSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_cent/$run
                alien_cp alien:/alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/GammaConv* file:$OUTPUTDIR_LHC17l3b_cent/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/$fileToDownloadQAEvent ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17l3bwSDDNotMerged.txt
            fi

            NotMergedruns=`cat runNumbersLHC17l3bwSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                rm $OUTPUTDIR_LHC17l3b_cent/$run/$fileName
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17l3b_cent.txt
                        fileNames =`cat fileStagedToMergeLHC17l3b_cent.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17l3b_cent/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_cent/$run/$fileName
                fi
            done;
        fi

        ################################ LHC17l3b_woSDD
        if [ $LHC17l3b_woSDD != "" ]; then
            rm runNumbersLHC17l3bwoSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l3b_woSDD
            runs=`cat runlists/runNumbersLHC17l3b_woSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run
                alien_cp alien:/alice/sim/2017/LHC17l3b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_woSDD/$run/
            done;
            if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileToDownloadQAEvent ]; then
                echo "file has already been copied for sucessfully for run " $runNumber
            else
                echo $run >> runNumbersLHC17l3bwoSDDNotMerged.txt
            fi

            NotMergedruns=`cat runNumbersLHC17l3bwoSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/sim/2017/LHC17l3b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                rm $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17l3b_woSDD.txt
                        fileNames =`cat fileStagedToMergeLHC17l3b_woSDD.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17l3b_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
                fi
            done;
        fi

        ################################ LHC17l4b_fast
        if [ $LHC17l4b_fast != "" ]; then
            rm runNumbersLHC17qfastNotMerged.txt
            echo $OUTPUTDIR_LHC17l4b_fast
            runs=`cat runlists/runNumbersLHC17l4b_fast.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l4b_fast/$run
                alien_cp alien:/alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/GammaConv* file:$OUTPUTDIR_LHC17l4b_fast/$run/
            done;

            NotMergedruns=`cat runNumbersLHC17qfastNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                rm $OUTPUTDIR_LHC17l4b_fast/$run/$fileName
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17l4b_fast.txt
                        fileNames =`cat fileStagedToMergeLHC17l4b_fast.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17l4b_fast/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_fast/$run/$fileName
                fi
            done;
        fi

        ################################ LHC17l4b_cent
        if [ $LHC17l4b_cent != "" ]; then
            rm runNumbersLHC17qwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l4b_cent
            runs=`cat runlists/runNumbersLHC17l4b_wSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l4b_cent/$run
                alien_cp alien:/alice/sim/2017/LHC17l4b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/GammaConv* file:$OUTPUTDIR_LHC17l4b_cent/$run/
            done;

            NotMergedruns=`cat runNumbersLHC17qwSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l4b_cent/$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/sim/2017/LHC17l4b_cent/$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                rm $OUTPUTDIR_LHC17l4b_cent/$run/$fileName
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17l4b_cent.txt
                        fileNames =`cat fileStagedToMergeLHC17l4b_cent.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17l4b_cent/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_cent/$run/$fileName
                fi
            done;
        fi

        ################################ LHC17l4b_cent_woSDD
        if [ $LHC17l4b_cent_woSDD != "" ]; then
            rm runNumbersLHC17qwoSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l4b_woSDD
            runs=`cat runlists/runNumbersLHC17l4b_woSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l4b_woSDD/$run
                alien_cp alien:/alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/GammaConv* file:$OUTPUTDIR_LHC17l4b_woSDD/$run/
            done;

            NotMergedruns=`cat runNumbersLHC17qwoSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput
                    alien_cp alien:/alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/
                done;
                counter=0;
                rm $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName
                if [ $fileName != "AnalysisResults.root" ]; then
                    for stageOutput in $stageOutputs; do
                        ls $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/GammaConvV1_*.root > fileStagedToMergeLHC17l4b_woSDD.txt
                        fileNames =`cat fileStagedToMergeLHC17l4b_woSDD.txt`
                        for fileName in $fileNames; do
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
                    done;
                    mv $OUTPUTDIR_LHC17l4b_woSDD/$run/intermediate_$run.root $OUTPUTDIR_LHC17l4b_woSDD/$run/$fileName
                fi
            done;
        fi

        if [ $fileName != "AnalysisResults.root" ]; then

            runs=`cat runlists/runNumbersLHC17l_all.txt`
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

                ls $OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_*.root > fileLHC17l3b_cent.txt
                fileNumbers=`cat fileLHC17l3b_cent.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_LHC17l3b_cent_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l3b_cent/$run/GammaConvV1_LHC17l3b_cent_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_cent/$run/CutSelection_LHC17l3b_cent_$number.log\"\)
                done;

                ls $OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_*.root > fileLHC17l3b_woSDD.txt
                fileNumbers=`cat fileLHC17l3b_woSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_LHC17l3b_woSDD_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/$run/GammaConvV1_LHC17l3b_woSDD_$number.root\"\,\"$OUTPUTDIR_LHC17l3b_woSDD/$run/CutSelection_LHC17l3b_woSDD_$number.log\"\)
                done;

                ls $OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_*.root > fileLHC17l4b_fast.txt
                fileNumbers=`cat fileLHC17l4b_fast.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_LHC17l4b_fast_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l4b_fast/$run/GammaConvV1_LHC17l4b_fast_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_fast/$run/CutSelection_LHC17l4b_fast_$number.log\"\)
                done;

                ls $OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_*.root > fileLHC17l4b_wSDD.txt
                fileNumbers=`cat fileLHC17l4b_wSDD.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_LHC17l4b_wSDD_$number.root\"\,\"GammaConvV1_$number\"\)
                    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR_LHC17l4b_cent/$run/GammaConvV1_LHC17l4b_wSDD_$number.root\"\,\"$OUTPUTDIR_LHC17l4b_cent/$run/CutSelection_LHC17l4b_wSDD_$number.log\"\)
                done;

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

    if [ $3 = "QA" ]; then

        ################################ LHC17l3b_fast
        if [ $LHC17l3b_fast != "" ]; then
            OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/LHC17l3b_fast
            mkdir -p $OUTPUTDIR_LHC17l3b_fast
            rm runNumbersLHC17l3bfastNotMerged.txt
            echo $OUTPUTDIR_LHC17l3b_fast
            runs=`cat runlists/runNumbersLHC17l3b_fast.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run
                alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC17l3b_fast/AnalysisResults.root file:$OUTPUTDIR_LHC17l3b_fast/$run/
                if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/AnalysisResults.root ]; then
                    echo "file has already been copied"
                else
                    echo $run >> runNumbersLHC17l3bfastNotMerged.txt
                fi
            done;

            NotMergedruns=`cat runNumbersLHC17l3bfastNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_fast/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC17l3b_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput
                    if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC17l3b_fast/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_fast/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17l3b_cent
        if [ $LHC17l3b_cent != "" ]; then
            OUTPUTDIR_LHC17l3b_cent=$BASEDIR/$TRAINDIR/LHC17l3b_wSDD
            mkdir -p $OUTPUTDIR_LHC17l3b_cent
            rm runNumbersLHC17l3bwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l3b_cent
            runs=`cat runlists/runNumbersLHC17l3b_wSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_cent/$run
                alien_cp alien:/alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/AnalysisResults.root file:$OUTPUTDIR_LHC17l3b_cent/$run/
                if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2017/LHC17l3b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_cent/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17l3b_woSDD
        if [ $LHC17l3b_woSDD != "" ]; then
            OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/LHC17l3b_woSDD
            mkdir -p $OUTPUTDIR_LHC17l3b_woSDD
            rm runNumbersLHC17l3bwoSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l3b_woSDD
            runs=`cat runlists/runNumbersLHC17l3b_woSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run
                alien_cp alien:/alice/sim/2017/LHC17l3b_cent_woSDD/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/AnalysisResults.root file:$OUTPUTDIR_LHC17l3b_woSDD/$run/
                if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/AnalysisResults.root ]; then
                    echo "file has already been copied"
                else
                    echo $run >> runNumbersLHC17l3bwoSDDNotMerged.txt
                fi
            done;

            NotMergedruns=`cat runNumbersLHC17l3bwoSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2017/LHC17l3b_cent_woSDD/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput
                    if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2017/LHC17l3b_cent_woSDD/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_woSDD/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17l4b_fast
        if [ $LHC17l4b_fast != "" ]; then
            rm runNumbersLHC17qfastNotMerged.txt
            echo $OUTPUTDIR_LHC17l4b_fast
            runs=`cat runlists/runNumbersLHC17l4b_fast.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l4b_fast/$run
                alien_cp alien:/alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/AnalysisResults.root file:$OUTPUTDIR_LHC17l4b_fast/$run/
                if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2017/LHC17l4b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l4b_fast/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17l4b_cent
        if [ $LHC17l4b_cent != "" ]; then
            rm runNumbersLHC17qwSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l4b_cent
            runs=`cat runlists/runNumbersLHC17l4b_wSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l4b_cent/$run
                alien_cp alien:/alice/sim/2017/LHC17l4b_cent/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/AnalysisResults.root file:$OUTPUTDIR_LHC17l4b_cent/$run/
                if [ -f $OUTPUTDIR_LHC17l4b_cent/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2017/LHC17l4b_cent/$run/pass1_CENT_wSDD/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l4b_cent/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        ################################ LHC17l4b_cent_woSDD
        if [ $LHC17l4b_cent_woSDD != "" ]; then
            rm runNumbersLHC17qwoSDDNotMerged.txt
            echo $OUTPUTDIR_LHC17l4b_woSDD
            runs=`cat runlists/runNumbersLHC17l4b_woSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l4b_woSDD/$run
                alien_cp alien:/alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/AnalysisResults.root file:$OUTPUTDIR_LHC17l4b_woSDD/$run/
                if [ -f $OUTPUTDIR_LHC17l4b_woSDD/$run/AnalysisResults.root ]; then
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
                    if [ -f $OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2017/LHC17l4b_cent_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l4b_woSDD/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi

        if [ $LHC18j2_fast != "" ]; then
            OUTPUTDIR_LHC18j2_fast=$BASEDIR/$TRAINDIR/LHC18j2_fast
            mkdir -p $OUTPUTDIR_LHC18j2_fast
            rm runNumbersLHC18j2fastNotMerged.txt
            echo $OUTPUTDIR_LHC18j2_fast
            runs=`cat runlists/runNumbersLHC17l3b_fast.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC18j2_fast/$run
                alien_cp alien:/alice/sim/2018/LHC18j2_fast/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC18j2_fast/AnalysisResults.root file:$OUTPUTDIR_LHC18j2_fast/$run/
                if [ -f $OUTPUTDIR_LHC18j2_fast/$run/AnalysisResults.root ]; then
                    echo "file has already been copied"
                else
                    echo $run >> runNumbersLHC18j2fastNotMerged.txt
                fi
            done;

            NotMergedruns=`cat runNumbersLHC18j2fastNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2018/LHC18j2_fast/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC18j2_fast/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC18j2_fast/$run/Stage_1/$stageOutput
                    if [ -f $OUTPUTDIR_LHC18j2_fast/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2018/LHC18j2_fast/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC18j2_fast/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC18j2_fast/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi


        if [ $LHC18j2_woSDD != "" ]; then
            OUTPUTDIR_LHC18j2_woSDD=$BASEDIR/$TRAINDIR/LHC18j2_woSDD
            mkdir -p $OUTPUTDIR_LHC18j2_woSDD
            rm runNumbersLHC18j2woSDDNotMerged.txt
            echo $OUTPUTDIR_LHC18j2_woSDD
            runs=`cat runlists/runNumbersLHC17l3b_woSDD.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC18j2_woSDD/$run
                alien_cp alien:/alice/sim/2018/LHC18j2_woSDD/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC18j2_woSDD/AnalysisResults.root file:$OUTPUTDIR_LHC18j2_woSDD/$run/
                if [ -f $OUTPUTDIR_LHC18j2_woSDD/$run/AnalysisResults.root ]; then
                    echo "file has already been copied"
                else
                    echo $run >> runNumbersLHC18j2woSDDNotMerged.txt
                fi
            done;

            NotMergedruns=`cat runNumbersLHC18j2woSDDNotMerged.txt`
            for run in $NotMergedruns; do
                echo "copying stage_1 output for " $run
                mkdir -p $run
                stageOutputs=`alien_ls /alice/sim/2018/LHC18j2_woSDD/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC18j2_woSDD/Stage_1/`
                for stageOutput in $stageOutputs; do
                    mkdir -p $OUTPUTDIR_LHC18j2_woSDD/$run/Stage_1/$stageOutput
                    if [ -f $OUTPUTDIR_LHC18j2_woSDD/$run/Stage_1/$stageOutput/AnalysisResults.root ]; then
                        echo "file has already been copied"
                    else
                        alien_cp alien:/alice/sim/2018/LHC18j2_woSDD/$run/$AODflagPWGGA/$TRAINPATHMC/$LHC18j2_woSDD/Stage_1/$stageOutput/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC18j2_woSDD/$run/Stage_1/$stageOutput/
                    fi
                done;
            done;
        fi


    fi


    if [ $LHC17l3b_fast != "" ]; then
        OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/LHC17l3b_fast/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17l3b_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17l3b_fast/
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
        OUTPUTDIR_LHC17l3b_cent=$BASEDIR/$TRAINDIR/LHC17l3b_wSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17l3b_cent
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_cent/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17l3b_cent/
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
        OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/LHC17l3b_woSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17l3b_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17l3b_woSDD/
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
        OUTPUTDIR_LHC17l4b_fast=$BASEDIR/$TRAINDIR/LHC17l4b_fast/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17l4b_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l4b_fast/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17l4b_fast/
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
        OUTPUTDIR_LHC17l4b_cent=$BASEDIR/$TRAINDIR/LHC17l4b_wSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17l4b_cent
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l4b_cent/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17l4b_cent/
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
        OUTPUTDIR_LHC17l4b_woSDD=$BASEDIR/$TRAINDIR/LHC17l4b_woSDD/$mergeFolder
        mkdir -p $OUTPUTDIR_LHC17l4b_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l4b_woSDD/$mergeFolder/GammaConv* file:$OUTPUTDIR_LHC17l4b_woSDD/
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

    if [ $LHC18d6b_fast != "" ]; then
        OUTPUTDIR_LHC18d6b_fast=$BASEDIR/$TRAINDIR/LHC18d6b_fast_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_fast/merge/GammaConv* file:$OUTPUTDIR_LHC18d6b_fast/
        ls $OUTPUTDIR_LHC18d6b_fast/GammaConv*.root > fileLHC18d6b_fast.txt
        fileNumbers=`cat fileLHC18d6b_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18d6b_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18d6b_fast_LowInt_$number.root\"\,\"GammaConvV1_$number\"\)
        done;
    fi

    if [ $LHC18d6b_woSDD != "" ]; then
        OUTPUTDIR_LHC18d6b_woSDD=$BASEDIR/$TRAINDIR/LHC18d6b_woSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_woSDD/merge/GammaConv* file:$OUTPUTDIR_LHC18d6b_woSDD/
        ls $OUTPUTDIR_LHC18d6b_woSDD/GammaConv*.root > fileLHC18d6b_woSDD.txt
        fileNumbers=`cat fileLHC18d6b_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18d6b_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18d6b_woSDD_LowInt_$number.root\"\,\"GammaConvV1_$number\"\)
        done;
    fi

    if [ $LHC18d6b_fastpart2 != "" ]; then
        OUTPUTDIR_LHC18d6b_fastpart2=$BASEDIR/$TRAINDIR/LHC18d6b_fastpart2_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_fastpart2
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_fastpart2/merge/GammaConv* file:$OUTPUTDIR_LHC18d6b_fastpart2/
        ls $OUTPUTDIR_LHC18d6b_fastpart2/GammaConv*.root > fileLHC18d6b_fastpart2.txt
        fileNumbers=`cat fileLHC18d6b_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18d6b_fastpart2/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18d6b_fastpart2_LowInt_$number.root\"\,\"GammaConvV1_$number\"\)
        done;
    fi

    if [ $LHC18d6b_woSDDpart2 != "" ]; then
        OUTPUTDIR_LHC18d6b_woSDDpart2=$BASEDIR/$TRAINDIR/LHC18d6b_woSDDpart2_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_woSDDpart2
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_woSDDpart2/merge/GammaConv* file:$OUTPUTDIR_LHC18d6b_woSDDpart2/
        ls $OUTPUTDIR_LHC18d6b_woSDDpart2/GammaConv*.root > fileLHC18d6b_woSDDpart2.txt
        fileNumbers=`cat fileLHC18d6b_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18d6b_woSDDpart2/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18d6b_woSDDpart2_LowInt_$number.root\"\,\"GammaConvV1_$number\"\)
        done;
    fi

    if [ $LHC18j2_fast != "" ]; then
        OUTPUTDIR_LHC18j2_fast=$BASEDIR/$TRAINDIR/LHC18j2_fast
        mkdir -p $OUTPUTDIR_LHC18j2_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18j2_fast/merge/GammaConv* file:$OUTPUTDIR_LHC18j2_fast/
        ls $OUTPUTDIR_LHC18j2_fast/GammaConv*.root > fileLHC18j2_fast.txt
        fileNumbers=`cat fileLHC18j2_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18j2_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18j2_fast_$number.root\"\,\"GammaConvV1_$number\"\)
        done;
    fi

    if [ $LHC18j2_woSDD != "" ]; then
        OUTPUTDIR_LHC18j2_woSDD=$BASEDIR/$TRAINDIR/LHC18j2_woSDD
        mkdir -p $OUTPUTDIR_LHC18j2_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18j2_woSDD/merge/GammaConv* file:$OUTPUTDIR_LHC18j2_woSDD/
        ls $OUTPUTDIR_LHC18j2_woSDD/GammaConv*.root > fileLHC18j2_woSDD.txt
        fileNumbers=`cat fileLHC18j2_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18j2_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18j2_woSDD_$number.root\"\,\"GammaConvV1_$number\"\)
        done;
    fi



#     if [ $LHC18d6b_wSDD != "" ]; then
#         OUTPUTDIR_LHC18d6b_wSDD=$BASEDIR/$TRAINDIR/LHC18d6b_wSDD_LowInt
#         mkdir -p $OUTPUTDIR_LHC18d6b_wSDD
#         alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_wSDD/merge/GammaConv* file:$OUTPUTDIR_LHC18d6b_wSDD/
#         ls $OUTPUTDIR_LHC18d6b_wSDD/GammaConv*.root > fileLHC18d6b_wSDD.txt
#         fileNumbers=`cat fileLHC18d6b_wSDD.txt`
#         for fileName in $fileNumbers; do
#             echo $fileName
#             number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
#             echo $number
#             root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18d6b_wSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18d6b_wSDD_LowInt_$number.root\"\,\"GammaConvV1_$number\"\)
#         done;
#     fi
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

    if [ $LHC17l3b_fast != "" ]; then
        OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/Material_LHC17l3b_fast_LowInt
        mkdir -p $OUTPUTDIR_LHC17l3b_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/merge_runlist_2/GammaConv_Material_* file:$OUTPUTDIR_LHC17l3b_fast/
        ls $OUTPUTDIR_LHC17l3b_fast/GammaConv*Material*.root > fileLHC17l3b_fast.txt
        fileNumbers=`cat fileLHC17l3b_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17l3b_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17l3b_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC17l3b_woSDD != "" ]; then
        OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/Material_LHC17l3b_woSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC17l3b_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/merge_runlist_2/GammaConv_Material_* file:$OUTPUTDIR_LHC17l3b_woSDD/
        ls $OUTPUTDIR_LHC17l3b_woSDD/GammaConv*Material*.root > fileLHC17l3b_woSDD.txt
        fileNumbers=`cat fileLHC17l3b_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17l3b_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    LHC17l3b_wSDD=$LHC17l3b_cent
    if [ $LHC17l3b_wSDD != "" ]; then
        OUTPUTDIR_LHC17l3b_wSDD=$BASEDIR/$TRAINDIR/Material_LHC17l3b_wSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC17l3b_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17l3b_wSDD/merge_runlist_2/GammaConv_Material_* file:$OUTPUTDIR_LHC17l3b_wSDD/
        ls $OUTPUTDIR_LHC17l3b_wSDD/GammaConv*Material*.root > fileLHC17l3b_wSDD.txt
        fileNumbers=`cat fileLHC17l3b_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17l3b_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17l3b_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6b_fast != "" ]; then
        OUTPUTDIR_LHC18d6b_fast=$BASEDIR/$TRAINDIR/Material_LHC18d6b_fast_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_fast/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6b_fast/
        ls $OUTPUTDIR_LHC18d6b_fast/GammaConv*Material*.root > fileLHC18d6b_fast.txt
        fileNumbers=`cat fileLHC18d6b_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6b_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6b_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6b_woSDD != "" ]; then
        OUTPUTDIR_LHC18d6b_woSDD=$BASEDIR/$TRAINDIR/Material_LHC18d6b_woSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_woSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6b_woSDD/
        ls $OUTPUTDIR_LHC18d6b_woSDD/GammaConv*Material*.root > fileLHC18d6b_woSDD.txt
        fileNumbers=`cat fileLHC18d6b_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6b_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6b_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6b_fastpart2 != "" ]; then
        OUTPUTDIR_LHC18d6b_fastpart2=$BASEDIR/$TRAINDIR/Material_LHC18d6b_fastpart2_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_fastpart2
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_fastpart2/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6b_fastpart2/
        ls $OUTPUTDIR_LHC18d6b_fastpart2/GammaConv*Material*.root > fileLHC18d6b_fastpart2.txt
        fileNumbers=`cat fileLHC18d6b_fastpart2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6b_fastpart2/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6b_fastpart2_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6b_woSDDpart2 != "" ]; then
        OUTPUTDIR_LHC18d6b_woSDDpart2=$BASEDIR/$TRAINDIR/Material_LHC18d6b_woSDDpart2_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_woSDDpart2
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_woSDDpart2/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6b_woSDDpart2/
        ls $OUTPUTDIR_LHC18d6b_woSDDpart2/GammaConv*Material*.root > fileLHC18d6b_woSDDpart2.txt
        fileNumbers=`cat fileLHC18d6b_woSDDpart2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6b_woSDDpart2/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6b_woSDDpart2_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6b_wSDD != "" ]; then
        OUTPUTDIR_LHC18d6b_wSDD=$BASEDIR/$TRAINDIR/Material_LHC18d6b_wSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6b_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6b_wSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6b_wSDD/
        ls $OUTPUTDIR_LHC18d6b_wSDD/GammaConv*Material*.root > fileLHC18d6b_wSDD.txt
        fileNumbers=`cat fileLHC18d6b_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6b_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6b_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi


    if [ $LHC18d6c_fast != "" ]; then
        OUTPUTDIR_LHC18d6c_fast=$BASEDIR/$TRAINDIR/Material_LHC18d6c_fast_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6c_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6c_fast/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6c_fast/
        ls $OUTPUTDIR_LHC18d6c_fast/GammaConv*Material*.root > fileLHC18d6c_fast.txt
        fileNumbers=`cat fileLHC18d6c_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6c_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6c_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6c_woSDD != "" ]; then
        OUTPUTDIR_LHC18d6c_woSDD=$BASEDIR/$TRAINDIR/Material_LHC18d6c_woSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6c_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6c_woSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6c_woSDD/
        ls $OUTPUTDIR_LHC18d6c_woSDD/GammaConv*Material*.root > fileLHC18d6c_woSDD.txt
        fileNumbers=`cat fileLHC18d6c_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6c_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6c_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6c_woSDD_part2 != "" ]; then
        OUTPUTDIR_LHC18d6c_woSDD=$BASEDIR/$TRAINDIR/Material_LHC18d6c_woSDD_LowInt_part2
        mkdir -p $OUTPUTDIR_LHC18d6c_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6c_woSDD_part2/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6c_woSDD/
        ls $OUTPUTDIR_LHC18d6c_woSDD/GammaConv*Material*.root > fileLHC18d6c_woSDD.txt
        fileNumbers=`cat fileLHC18d6c_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6c_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6c_woSDDpart2_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6c_wSDD != "" ]; then
        OUTPUTDIR_LHC18d6c_wSDD=$BASEDIR/$TRAINDIR/Material_LHC18d6c_wSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6c_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6c_wSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6c_wSDD/
        ls $OUTPUTDIR_LHC18d6c_wSDD/GammaConv*Material*.root > fileLHC18d6c_wSDD.txt
        fileNumbers=`cat fileLHC18d6c_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6c_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6c_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC17b3lLowInt_fast != "" ]; then
        OUTPUTDIR_LHC17b3lLowInt_fast=$BASEDIR/$TRAINDIR/Material_LHC17b3lLowInt_fast_LowInt
        mkdir -p $OUTPUTDIR_LHC17b3lLowInt_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17b3lLowInt_fast/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC17b3lLowInt_fast/
        ls $OUTPUTDIR_LHC17b3lLowInt_fast/GammaConv*Material*.root > fileLHC17b3lLowInt_fast.txt
        fileNumbers=`cat fileLHC17b3lLowInt_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17b3lLowInt_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17b3lLowInt_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC17b3lLowInt_woSDD != "" ]; then
        OUTPUTDIR_LHC17b3lLowInt_woSDD=$BASEDIR/$TRAINDIR/Material_LHC17b3lLowInt_woSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC17b3lLowInt_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17b3lLowInt_woSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC17b3lLowInt_woSDD/
        ls $OUTPUTDIR_LHC17b3lLowInt_woSDD/GammaConv*Material*.root > fileLHC17b3lLowInt_woSDD.txt
        fileNumbers=`cat fileLHC17b3lLowInt_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17b3lLowInt_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17b3lLowInt_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC17b3lLowInt_wSDD != "" ]; then
        OUTPUTDIR_LHC17b3lLowInt_wSDD=$BASEDIR/$TRAINDIR/Material_LHC17b3lLowInt_wSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC17b3lLowInt_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17b3lLowInt_wSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC17b3lLowInt_wSDD/
        ls $OUTPUTDIR_LHC17b3lLowInt_wSDD/GammaConv*Material*.root > fileLHC17b3lLowInt_wSDD.txt
        fileNumbers=`cat fileLHC17b3lLowInt_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17b3lLowInt_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17b3lLowInt_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6cLowInt_fast != "" ]; then
        OUTPUTDIR_LHC18d6cLowInt_fast=$BASEDIR/$TRAINDIR/Material_LHC18d6cLowInt_fast_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6cLowInt_fast
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6cLowInt_fast/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6cLowInt_fast/
        ls $OUTPUTDIR_LHC18d6cLowInt_fast/GammaConv*Material*.root > fileLHC18d6cLowInt_fast.txt
        fileNumbers=`cat fileLHC18d6cLowInt_fast.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6cLowInt_fast/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6cLowInt_fast_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

        if [ $LHC18d6cLowInt_woSDD != "" ]; then
        OUTPUTDIR_LHC18d6cLowInt_woSDD=$BASEDIR/$TRAINDIR/Material_LHC18d6cLowInt_woSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6cLowInt_woSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6cLowInt_woSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6cLowInt_woSDD/
        ls $OUTPUTDIR_LHC18d6cLowInt_woSDD/GammaConv*Material*.root > fileLHC18d6cLowInt_woSDD.txt
        fileNumbers=`cat fileLHC18d6cLowInt_woSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6cLowInt_woSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6cLowInt_woSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6cLowInt_woSDDpart2 != "" ]; then
        OUTPUTDIR_LHC18d6cLowInt_woSDDpart2=$BASEDIR/$TRAINDIR/Material_LHC18d6cLowInt_woSDDpart2_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6cLowInt_woSDDpart2
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6cLowInt_woSDDpart2/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6cLowInt_woSDDpart2/
        ls $OUTPUTDIR_LHC18d6cLowInt_woSDDpart2/GammaConv*Material*.root > fileLHC18d6cLowInt_woSDDpart2.txt
        fileNumbers=`cat fileLHC18d6cLowInt_woSDDpart2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6cLowInt_woSDDpart2/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6cLowInt_woSDDpart2_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC18d6cLowInt_wSDD != "" ]; then
        OUTPUTDIR_LHC18d6cLowInt_wSDD=$BASEDIR/$TRAINDIR/Material_LHC18d6cLowInt_wSDD_LowInt
        mkdir -p $OUTPUTDIR_LHC18d6cLowInt_wSDD
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC18d6cLowInt_wSDD/merge/GammaConv_Material_* file:$OUTPUTDIR_LHC18d6cLowInt_wSDD/
        ls $OUTPUTDIR_LHC18d6cLowInt_wSDD/GammaConv*Material*.root > fileLHC18d6cLowInt_wSDD.txt
        fileNumbers=`cat fileLHC18d6cLowInt_wSDD.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC18d6cLowInt_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC18d6cLowInt_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi


    if [ $3 = "runwise" ]; then

        if [ $LHC17l3b_fast != "" ]; then
            echo "Downloading " $LHC17l3b_fast " in " $OUTPUTDIR_LHC17l3b_fast
            runs=`cat runlists/runNumbersLHC17l3b_LowInt.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_fast/$run
                if [ -f $OUTPUTDIR_LHC17l3b_fast/$run/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l3b_fast/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_fast/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_fast/$run/
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

        if [ $LHC17l3b_woSDD != "" ]; then
            echo "Downloading " $LHC17l3b_woSDD " in " $OUTPUTDIR_LHC17l3b_woSDD
            runs=`cat runlists/runNumbersLHC17l3b_LowInt.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_woSDD/$run
                if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l3b_woSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_woSDD/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_woSDD/$run/
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

        if [ $LHC17l3b_wSDD != "" ]; then
            echo "Downloading " $LHC17l3b_wSDD " in " $OUTPUTDIR_LHC17l3b_wSDD
            runs=`cat runlists/runNumbersLHC17l3b_LowInt.txt`
            for run in $runs; do
                echo $run
                mkdir -p $OUTPUTDIR_LHC17l3b_wSDD/$run
                if [ -f $OUTPUTDIR_LHC17l3b_wSDD/$run/$fileName ]; then
                    echo "file has already been copied"
                else
                    alien_cp alien:/alice/sim/2017/LHC17l3b_wSDD/$run/PWGGA/$TRAINPATHMC/$LHC17l3b_wSDD/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC17l3b_wSDD/$run/
                fi
            done;

            hadd -f $OUTPUTDIR_LHC17l3b_wSDD/$fileName $OUTPUTDIR_LHC17l3b_wSDD/282008/$fileName $OUTPUTDIR_LHC17l3b_wSDD/282016/$fileName $OUTPUTDIR_LHC17l3b_wSDD/282021/$fileName $OUTPUTDIR_LHC17l3b_wSDD/282025/$fileName $OUTPUTDIR_LHC17l3b_wSDD/282031/$fileName

            ls $OUTPUTDIR_LHC17l3b_wSDD/GammaConv*Material*.root > fileLHC17l3b_wSDD.txt
            fileNumbers=`cat fileLHC17l3b_wSDD.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 3 | cut -d "." -f1`
                echo $number
                root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17l3b_wSDD/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17l3b_wSDD_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
            done;
        fi
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

#     fileName=GammaConvV1_400.root;

    OUTPUTDIR_LHC18b8_fast=$BASEDIR/$TRAINDIR/LHC18b8_fast
    mkdir -p $OUTPUTDIR_LHC18b8_fast
    OUTPUTDIR_LHC18b8_woSDD=$BASEDIR/$TRAINDIR/LHC18b8_woSDD
    mkdir -p $OUTPUTDIR_LHC18b8_woSDD
#     OUTPUTDIR_LHC18b8_cent=$BASEDIR/$TRAINDIR/LHC18b8_wSDD
#     mkdir -p $OUTPUTDIR_LHC18b8_cent

#     ############################### LHC18b8_fast
#     rm runNumbersLHC18b8fastNotMerged.txt
#     echo $OUTPUTDIR_LHC18b8_fast
#     runs=`cat runlists/runNumbersLHC18b8_fast.txt`
#     for run in $runs; do
#         echo $run
#         binNumbersJJ=`cat runlists/binsJetJetLHC18b8.txt`
#         for binNumber in $binNumbersJJ; do
#             echo $binNumber
#             mkdir -p $OUTPUTDIR_LHC18b8_fast/$binNumber/$run
#             OutputFiles=`alien_ls /alice/sim/2018/LHC18b8_fast/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_fast/GammaConv*`
#             for file in $OutputFiles; do
#                 if [ -f $OUTPUTDIR_LHC18b8_fast/$binNumber/$run/$file ]; then
#                     echo "file " $OUTPUTDIR_LHC18b8_fast/$binNumber/$run/$file " has already been copied"
#                 else
#                     alien_cp alien:/alice/sim/2018/LHC18b8_fast/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_fast/GammaConv* file:$OUTPUTDIR_LHC18b8_fast/$binNumber/$run
#                 fi
#             done;
#         done;
#     done;
#
#     ################################ LHC18b8_woSDD
#     rm runNumbersLHC18b8woSDDNotMerged.txt
#     echo $OUTPUTDIR_LHC18b8_woSDD
#     runs=`cat runlists/runNumbersLHC18b8_woSDD.txt`
#     for run in $runs; do
#         echo $run
#         binNumbersJJ=`cat runlists/binsJetJetLHC18b8.txt`
#         for binNumber in $binNumbersJJ; do
#             echo $binNumber
#             mkdir -p $OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run
#             OutputFiles=`alien_ls /alice/sim/2018/LHC18b8_cent_woSDD/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_woSDD/GammaConv*`
#             for file in $OutputFiles; do
#                 if [ -f $OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run/$file ]; then
#                     echo "file " $OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run/$file " has already been copied"
#                 else
#                     alien_cp alien:/alice/sim/2018/LHC18b8_cent_woSDD/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_woSDD/GammaConv* file:$OUTPUTDIR_LHC18b8_woSDD/$binNumber/$run
#                 fi
#             done;
#         done;
#     done;
#
#     ################################ LHC18b8_cent
#     rm runNumbersLHC18b8wSDDNotMerged.txt
#     echo $OUTPUTDIR_LHC18b8_cent
#     runs=`cat runlists/runNumbersLHC18b8_wSDD.txt`
#     for run in $runs; do
#         echo $run
#         binNumbersJJ=`cat runlists/binsJetJetLHC18b8.txt`
#         for binNumber in $binNumbersJJ; do
#             echo $binNumber
#             mkdir -p $OUTPUTDIR_LHC18b8_cent/$binNumber/$run
#             if [ -f $OUTPUTDIR_LHC18b8_cent/$binNumber/$run/$fileName ]; then
#                 echo "file has already been copied"
#             else
#             alien_cp alien:/alice/sim/2018/LHC18b8_cent/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC18b8_cent/$fileToDownloadQAEvent file:$OUTPUTDIR_LHC18b8_cent/$binNumber/$run
#             fi
#             if [ -f $OUTPUTDIR_LHC18b8_cent/$binNumber/$run/$fileName ]; then
#                 echo "file has already been copied"
#             else
#                 echo $run >> runNumbersLHC18b8fastNotMerged.txt
#             fi


    rm $OUTPUTDIR_LHC18b8_fast/GammaConv*.root
    firstrunNumber=`head -n1 runlists/runNumbersLHC18b8_fast.txt`
    firstbinNumber=`head -n1 runlists/binsJetJetLHC18b8.txt`
    ls $OUTPUTDIR_LHC18b8_fast/$firstbinNumber/$firstrunNumber/GammaConvV1_*.root > fileJJmergeLHC18b8.txt
    fileNumbers=`cat fileJJmergeLHC18b8.txt`
    MergeAccordingToSpecificRunlist fileJJmergeLHC18b8.txt $OUTPUTDIR_LHC18b8_fast $NSlashes3 GammaConvV1 fast runlists/runNumbersLHC18b8_fast.txt runlists/binsJetJetLHC18b8.txt

    ls $OUTPUTDIR_LHC18b8_fast/GammaConv*.root > fileLHC18b8_fast.txt
    fileNumbers=`cat fileLHC18b8_fast.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18b8_fast/GammaConvV1-fast_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18b8_fast_$number.root\"\,\"GammaConvV1_$number\"\)
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC18b8_fast_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC18b8_fast_$number.log\"\)
    done;


    rm $OUTPUTDIR_LHC18b8_woSDD/GammaConv*.root
    firstrunNumber=`head -n1 runlists/runNumbersLHC18b8_woSDD.txt`
    firstbinNumber=`head -n1 runlists/binsJetJetLHC18b8.txt`
    ls $OUTPUTDIR_LHC18b8_woSDD/$firstbinNumber/$firstrunNumber/GammaConvV1_*.root > fileJJmergeLHC18b8.txt
    fileNumbers=`cat fileJJmergeLHC18b8.txt`
    MergeAccordingToSpecificRunlist fileJJmergeLHC18b8.txt $OUTPUTDIR_LHC18b8_woSDD $NSlashes3 GammaConvV1 woSDD runlists/runNumbersLHC18b8_woSDD.txt runlists/binsJetJetLHC18b8.txt

    ls $OUTPUTDIR_LHC18b8_woSDD/GammaConv*.root > fileLHC18b8_woSDD.txt
    fileNumbers=`cat fileLHC18b8_woSDD.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC18b8_woSDD/GammaConvV1-woSDD_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC18b8_woSDD_$number.root\"\,\"GammaConvV1_$number\"\)
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC18b8_woSDD_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC18b8_woSDD_$number.log\"\)
    done;


#     rm $OUTPUTDIR_LHC18b8_cent/GammaConv*.root
#     firstrunNumber=`head -n1 runlists/runNumbersLHC18b8_wSDD.txt`
#     firstbinNumber=`head -n1 runlists/binsJetJetLHC18b8.txt`
#     ls $OUTPUTDIR_LHC18b8_cent/$firstbinNumber/$firstrunNumber/GammaConvV1_*.root > fileJJmergeLHC18b8.txt
#     fileNumbers=`cat fileJJmergeLHC18b8.txt`
#     MergeAccordingToSpecificRunlist fileJJmergeLHC18b8.txt $OUTPUTDIR_LHC18b8_cent $NSlashes3 GammaConvV1 wSDD runlists/runNumbersLHC18b8_wSDD.txt runlists/binsJetJetLHC18b8.txt
#


fi

if [ $2 = "mergeAllMC" ]; then

    fileName="GammaConvV1_400.root"
#     OUTPUTDIR_LHC17l3b_fast=$BASEDIR/$TRAINDIR/LHC17l3b_fast
#     OUTPUTDIR_LHC17l3b_cent=$BASEDIR/$TRAINDIR/LHC17l3b_wSDD
    OUTPUTDIR_LHC17l3b_woSDD=$BASEDIR/$TRAINDIR/LHC17l3b_woSDD
#     mkdir -p $OUTPUTDIR_LHC17l3b_fast
#     mkdir -p $OUTPUTDIR_LHC17l3b_cent
#     mkdir -p $OUTPUTDIR_LHC17l3b_woSDD

#     OUTPUTDIR_LHC17l4b_fast=$BASEDIR/$TRAINDIR/LHC17l4b_fast
#     OUTPUTDIR_LHC17l4b_cent=$BASEDIR/$TRAINDIR/LHC17l4b_wSDD
#     OUTPUTDIR_LHC17l4b_woSDD=$BASEDIR/$TRAINDIR/LHC17l4b_woSDD
#     mkdir -p $OUTPUTDIR_LHC17l4b_fast
#     mkdir -p $OUTPUTDIR_LHC17l4b_cent
#     mkdir -p $OUTPUTDIR_LHC17l4b_woSDD

    runs=`cat runlists/runNumbersLHC17l3b_woSDD.txt`
    counter=0;
    fileName="GammaConvV1_400.root"
    for run in $runs; do
        echo $run
        if [ -f $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName ]; then
            if [ $counter = 0 ]; then
                cp $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
                counter=$(($counter+1));
                echo $counter;
            else
                hadd -f $OUTPUTDIR_LHC17l3b_woSDD/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$run/$fileName
                mv $OUTPUTDIR_LHC17l3b_woSDD/$fileName $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root
            fi
        fi
    done;
    mv $OUTPUTDIR_LHC17l3b_woSDD/intermediate.root $OUTPUTDIR_LHC17l3b_woSDD/$fileName
    ls $OUTPUTDIR_LHC17l3b_woSDD/GammaConvV1_*.root > fileLHC17l3b_woSDD.txt
    fileNumbers=`cat fileLHC17l3b_woSDD.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l3b_woSDD/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l3b_woSDD-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l3b_woSDD_$number.log\"\)
    done;

#     runs=`cat runlists/runNumbersLHC17l3b_fast.txt`
#     counter=0;
#     fileName="GammaConv_Material_121.root"
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
#     runs=`cat runlists/runNumbersLHC17l3b_wSDD.txt`
#     counter=0;
#     fileName="GammaConvV1_440.root"
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
#
#     runs=`cat runlists/runNumbersLHC17l4b_wSDD.txt`
#     counter=0;
#     fileName="GammaConvV1_440.root"
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
#     fileName="GammaConvV1_440.root"
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
#
#     runs=`cat runlists/runNumbersLHC17l4b_fast.txt`
#
#     counter=0;
#     fileName="GammaConvV1_440.root"
#     for run in $runs; do
#         echo $run
#         if [ -f $OUTPUTDIR_LHC17l4b_fast/$run/$fileName ]; then
#             if [ $counter = 0 ]; then
#                 cp $OUTPUTDIR_LHC17l4b_fast/$run/$fileName $OUTPUTDIR_LHC17l4b_fast/intermediate.root
#                 counter=$(($counter+1));
#                 echo $counter;
#             else
#                 hadd -f $OUTPUTDIR_LHC17l4b_fast/$fileName $OUTPUTDIR_LHC17l4b_fast/intermediate.root $OUTPUTDIR_LHC17l4b_fast/$run/$fileName
#                 mv $OUTPUTDIR_LHC17l4b_fast/$fileName $OUTPUTDIR_LHC17l4b_fast/intermediate.root
#             fi
#         fi
#     done;
#     mv $OUTPUTDIR_LHC17l4b_fast/intermediate.root $OUTPUTDIR_LHC17l4b_fast/$fileName
#     ls $OUTPUTDIR_LHC17l4b_fast/GammaConvV1_*.root > fileLHC17l4b_fast.txt
#     fileNumbers=`cat fileLHC17l4b_fast.txt`
#     for fileName in $fileNumbers; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
#         echo $number
#         root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17l4b_fast/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast-pass1_$number.root\"\,\"GammaConvV1_$number\"\)
#         root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC17l4b_fast-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC17l4b_fast_$number.log\"\)
#     done;

fi

if [ $1 = "2015" ]; then

    if [ $LHC15n_pass4 != "" ]; then
        OUTPUTDIR_LHC15n=$BASEDIR/$TRAINDIR/LHC15n_pass4 #_LowInt (runlist 5)
        mkdir -p $OUTPUTDIR_LHC15n

        echo "Downloading 2015 data, pass4: " $LHC15n_pass4
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC15n_pass4/merge_runlist_1/GammaConv* file:$OUTPUTDIR_LHC15n/
        ls $OUTPUTDIR_LHC15n/GammaConv*.root > fileLHC15n_pass4.txt
        fileNumbers=`cat fileLHC15n_pass4.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
#             number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC15n/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_LHC15n_pass4_$number.root\"\,\"GammaConvV1_$number\"\)
#             root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC15n/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC15n_pass4_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)
        done;
    fi

    if [ $LHC16k5a_pass3 != "" ]; then
        OUTPUTDIR_LHC16k5a=$BASEDIR/$TRAINDIR/LHC16k5a_pass3
        mkdir -p $OUTPUTDIR_LHC16k5a

        echo "Downloading 2015 MC, pass3: " $LHC16k5a_pass3
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC16k5a_pass3/merge/GammaConv* file:$OUTPUTDIR_LHC16k5a/
        ls $OUTPUTDIR_LHC16k5a/GammaConv*.root > fileLHC16k5a_pass3.txt
        fileNumbers=`cat fileLHC16k5a_pass3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
#             number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16k5a/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC16k5a_pass3_$number.root\"\,\"GammaConvV1_$number\"\)
#             root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC16k5a/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC16k5a_pass3_$number.root\"\,\"GammaConvMaterial_$number\"\)

        done;
    fi

    if [ $LHC16k5b_pass3 != "" ]; then
        OUTPUTDIR_LHC16k5b=$BASEDIR/$TRAINDIR/LHC16k5b_pass3
        mkdir -p $OUTPUTDIR_LHC16k5b

        echo "Downloading 2015 MC, pass3: " $LHC16k5b_pass3
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC16k5b_pass3/merge/GammaConv* file:$OUTPUTDIR_LHC16k5b/
        ls $OUTPUTDIR_LHC16k5b/GammaConv*.root > fileLHC16k5b_pass3.txt
        fileNumbers=`cat fileLHC16k5b_pass3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
#             number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16k5b/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC16k5b_pass3_$number.root\"\,\"GammaConvV1_$number\"\)
#             root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC16k5b/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC16k5b_pass3_$number.root\"\,\"GammaConvMaterial_$number\"\)

        done;
    fi

    if [ $LHC17e2_pass4 != "" ]; then
        OUTPUTDIR_LHC17e2=$BASEDIR/$TRAINDIR/LHC17e2_pass4 #_LowInt (runlist 2)
        mkdir -p $OUTPUTDIR_LHC17e2

        echo "Downloading 2015 MC, pass4: " $LHC17e2_pass4
        alien_cp alien:/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC17e2_pass4/merge_runlist_1/GammaConv* file:$OUTPUTDIR_LHC17e2/
        ls $OUTPUTDIR_LHC17e2/GammaConv*.root > fileLHC17e2_pass4.txt
        fileNumbers=`cat fileLHC17e2_pass4.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
#             number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC17e2/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC17e2_pass4_$number.root\"\,\"GammaConvV1_$number\"\)
#             root -l -b -q -x ChangeStructureToStandardMaterial.C\(\"$OUTPUTDIR_LHC17e2/GammaConv_Material_$number.root\"\,\"$OUTPUTDIR/MaterialBudget_LHC17e2_pass4_LowInt_$number.root\"\,\"GammaConvMaterial_$number\"\)

        done;
    fi

    if [ $LHC16h3_JJ_pass4 != "" ]; then

        OUTPUTDIR_LHC16h3_JJ_pass4=$BASEDIR/$TRAINDIR/LHC16h3_JJ_pass4
        mkdir -p $OUTPUTDIR_LHC16h3_JJ_pass4

        rm runNumbersLHC16h3NotMerged.txt
        echo $OUTPUTDIR_LHC16h3_JJ_pass4
        runs=`cat runlists/runNumbersLHC16h3.txt`
        for run in $runs; do
            echo $run
            binNumbersJJ=`cat runlists/binsJetJetLHC18b8.txt`
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                mkdir -p $OUTPUTDIR_LHC16h3_JJ_pass4/$binNumber/$run
                OutputFiles=`alien_ls /alice/sim/2016/LHC16h3/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC16h3_JJ_pass4/GammaConv*`
                for file in $OutputFiles; do
                    if [ -f $OUTPUTDIR_LHC16h3_JJ_pass4/$binNumber/$run/$file ]; then
                        echo "file " $OUTPUTDIR_LHC16h3_JJ_pass4/$binNumber/$run/$file " has already been copied"
                    else
                        alien_cp alien:/alice/sim/2016/LHC16h3/$binNumber/$run/PWGGA/$TRAINPATHMC/$LHC16h3_JJ_pass4/GammaConv* file:$OUTPUTDIR_LHC16h3_JJ_pass4/$binNumber/$run
                    fi
                done;
            done;
        done;

        rm $OUTPUTDIR_LHC16h3_JJ_pass4/GammaConv*.root
        firstrunNumber=`head -n1 runlists/runNumbersLHC16h3.txt`
        firstbinNumber=`head -n1 runlists/binsJetJetLHC18b8.txt`
        ls $OUTPUTDIR_LHC16h3_JJ_pass4/$firstbinNumber/$firstrunNumber/GammaConvV1_*.root > fileJJmergeLHC16h3.txt
        fileNumbers=`cat fileJJmergeLHC16h3.txt`
        MergeAccordingToSpecificRunlist fileJJmergeLHC16h3.txt $OUTPUTDIR_LHC16h3_JJ_pass4 $NSlashes3 GammaConvV1 pass4 runlists/runNumbersLHC16h3.txt runlists/binsJetJetLHC18b8.txt

        ls $OUTPUTDIR_LHC16h3_JJ_pass4/GammaConv*.root > fileLHC16h3_JJ_pass4.txt
        fileNumbers=`cat fileLHC16h3_JJ_pass4.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f 9 | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR_LHC16h3_JJ_pass4/GammaConvV1-pass4_$number.root\"\,\"$OUTPUTDIR/GammaConvV1_MC_LHC16h3_JJ_pass4_$number.root\"\,\"GammaConvV1_$number\"\)
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC16h3_JJ_pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_MC_LHC16h3_JJ_pass4_$number.log\"\)
        done;

    fi

fi