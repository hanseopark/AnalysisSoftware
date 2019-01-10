#! /bin/bash
# $1  is data file
# $2  is Pythia file
# $3  is Phojet file
# $4  is cocktail file
# $5  cut
# $6  energy
# $7 finder type 1:onfly; 2: offline


echo "The data file is $1 "
echo "The Pythia file is $2 "
echo "The Phojet file is $3 "
echo "The cocktial file is $4"
echo "The cut to be used is $5"
echo "The energy is  $6"
echo "The finder selected is $7"


DATAROOTFILE=$1
PYTHIAROOTFILE=$2
PHOJETROOTFILE=$3
COCKTAILROOTFILE=$4
CUT=$5
ENERGY=$6
FINDERTYPE=$7

CUT1=""
CUT2=""
FINDERNAME1=""
FINDERNAME2=""

NAMEPYTHIA=""
NAMEPHOJET=""

 if [ $ENERGY == "13TeV" ]; then
   NAMEPYTHIA=LHC16_Pythia
   NAMEPHOJET=LHC16_Phojet
 fi

if [ $ENERGY == "5TeV" ]; then
   NAMEPYTHIA=LHC17_Pythia
   NAMEPHOJET=LHC17_Phojet
fi

if [ $FINDERTYPE -eq 1 ]; then
   FINDERNAME=MaterialBudgetHistosOnfly
   FINDERNAME1=$FINDERNAME
   CUT1=$CUT
   rm weightStudies.log
fi

if [ $FINDERTYPE -eq 2 ]; then
   FINDERNAME=MaterialBudgetHistosOffline
   FINDERNAME2=$FINDERNAME
   CUT2=$CUT
fi





if [ $FINDERTYPE -eq 1 ]; then
    echo "    "
    echo " Running prepare secondaries for Onfly"
    root -l -b -x -q TaskV1/PrepareSecondaries.C\+\(\"Pi0\",\"$COCKTAILROOTFILE\",\"pdf\",\"$CUT1\_0152103500000000\",\"$ENERGY\",\"No\",\"0.80\",\"\",70,0,kFALSE\)
    echo "  "
    echo "  "
    echo "Onfly V0 finder running Pythia"
    # running with Pythia Onfly
    root -b -q -x -l TaskV1/AnalyseMaterialHistosV2.C\+\(\"$DATAROOTFILE\"\,\"$PYTHIAROOTFILE\"\,\"$CUT1\"\,\"$ENERGY\"\,\"$NAMEPYTHIA\"\,\"pdf\",\"$FINDERNAME\"\)
    # running with Phojet Onfly
    echo "  "
    echo "Onfly V0 finder running Phojet"
    root -b -q -x -l TaskV1/AnalyseMaterialHistosV2.C\+\(\"$DATAROOTFILE\"\,\"$PHOJETROOTFILE\"\,\"$CUT1\"\,\"$ENERGY\"\,\"$NAMEPHOJET\"\,\"pdf\",\"$FINDERNAME\"\)

    ls MCInputFileMaterialBudgetWeights$NAMEPYTHIA\_$CUT1\.root > weightStudies.log
    ls MCInputFileMaterialBudgetWeights$NAMEPHOJET\_$CUT1\.root >> weightStudies.log

fi

if [ $FINDERTYPE -eq 2 ]; then
    echo "    "
    echo " Running prepare secondaries for Offline"
    root -l -b -x -q TaskV1/PrepareSecondaries.C\+\(\"Pi0\",\"$COCKTAILROOTFILE\",\"pdf\",\"$CUT2\_0152103500000000\",\"$ENERGY\",\"No\",\"0.80\",\"\",70,0,kFALSE\)
    echo "  "
    echo "  "
    echo "Offline V0 finder running Pythia"
    # running with Pythia Offline
    root -b -q -x -l TaskV1/AnalyseMaterialHistosV2.C\+\(\"$DATAROOTFILE\"\,\"$PYTHIAROOTFILE\"\,\"$CUT2\"\,\"$ENERGY\"\,\"$NAMEPYTHIA\"\,\"pdf\",\"$FINDERNAME\"\)
    echo "  "
    echo "Offline V0 finder running Phojet"
    # running with Phojet Offline
    root -b -q -x -l TaskV1/AnalyseMaterialHistosV2.C\+\(\"$DATAROOTFILE\"\,\"$PHOJETROOTFILE\"\,\"$CUT2\"\,\"$ENERGY\"\,\"$NAMEPHOJET\"\,\"pdf\",\"$FINDERNAME\"\)
    ls MCInputFileMaterialBudgetWeights$NAMEPYTHIA\_$CUT2\.root >> weightStudies.log
    ls MCInputFileMaterialBudgetWeights$NAMEPHOJET\_$CUT2\.root >> weightStudies.log
fi

 echo "  "
 echo "  "
 echo "  "
 echo "  "

if [ $FINDERTYPE -eq 2 ]; then

    echo "  "
    echo "Running systematic evaluation"

    root -l -b -x -q TaskV1/WeightStudiesOverview.C\+\(\"weightStudies.log\"\,\"pdf\"\,\"$ENERGY\"\,\"PythiaPhojetMultWeightsPtWeights$ENERGY\",4,0,0,\"kFALSE\"\)  
fi


