#! /bin/bash
#
#
#
# This script gests as input a directory where the GammaConvV1 root file is stored,
# it also needs the desired output directory where the produced root files are put.
# If nothing is given it will use ./ for the input directory and ./Output for the output
#
#Input 1: Root file to analyze Default: AnalyisResults
#Input 2: Input directory  Default:$PWD
#Input 3: Output directory Default: $PWD/Results  (directory will be created if it does not exist)
#

PROGNAME=$0
ROOTVERSION=`echo $ROOT_RELEASE`
ISROOT6=0

# running modes
ONLYCORRECTION=0
ONLYCUTS=0

# general settings
SUFFIX=eps
NORMALCUTS=0
DATAFILE=1
DATAFILEOK=0
MCFILE=1
NAMECUTSTUDIES="none"
PERIODNAME="No"
USETHNSPARSE=0
BINSPTGAMMA=0
BINSPTPI0=0
BINSPTETA=0
BINSPTETAPRIME=0
ENERGY=""
MODE=-1
DIRECTPHOTON="No"

# fileNames
DATAROOTFILE=""
MCROOTFILE=""
MCROOTFILEADDSIG=""
MCROOTFILEADDSIGETA=""
MCROOTFILEMERGEA=""
MCROOTFILEMERGEB=""
MCROOTFILEGJ=""

# special settings
MERGINGMC=0
ADVMESONQA=""
NEWGAMMAMACROS=1
ADDEDSIG=0
ADDPILEUPCORR=kFALSE
FILEWEIGHTGAMMA=""

# Boolean for switching on and of different meson processings
DOGAMMA=1
DOPI0=1
DOETA=1
DOPI0INETABINS=1
    # Heavy meson configuration set based on mode>100,
    # search for "Set heavy meson configuration" in this script
DOETAPRIME=0
DOPI0TAG=0



# switch for turning of ToyMC
DISABLETOYMC=0
NEVTSTOY=1e7
MINPTTOY=0
MAXPTTOY=70
MAXPTTOYLAMBDA=90
EXTINPUTFILE=""

# include functions for setting number of bins
source start_FullMesonAnalysis_helperPP.sh
source start_FullMesonAnalysis_helperPPb.sh
source start_FullMesonAnalysis_helperPbPb.sh

function ExtractSignal()
{
    root -x -q -l -b  TaskV1/ExtractSignalV2.C\+\($1\,$MODE\,$USETHNSPARSE\)
}

function ExtractSignalGamma()
{
    if [ $ISROOT6 -eq 0 ]; then
        root -x -l -b -q TaskV1/ExtractGammaSignal.C\+\($1\,$ADDEDSIG\,$MODE\)
    else
        root -x -q -l -b TaskV1/ExtractGammaSignal.C\($1\,$ADDEDSIG\,$MODE\)
    fi
}

function ExtractSignalGammaV2()
{
    if [ $ISROOT6 -eq 0 ]; then
        root -x -l -b -q TaskV1/ExtractGammaSignalV2.C\+\($1\,0\,\"\"\,\"$2\"\)
    else
        root -x -l -b -q TaskV1/ExtractGammaSignalV2.C\($1\,0\,\"\"\,\"$2\"\)
    fi
}


function CorrectSignal()
{
    root -x -l -b -q TaskV1/CorrectSignalV2.C\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$ENERGY\"\,\"\"\,\"$ESTIMATEPILEUP\"\,kFALSE\,$MODE\)

}

function CorrectSignalGamma()
{
    if [ $ISROOT6 -eq 0 ]; then
        root -x -b -q -l CompileCorrectGamma.C\+\+
        root -x -l -b -q TaskV1/CorrectGamma.C\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$ENERGY\"\,\"\"\,\"$ESTIMATEPILEUP\"\,$MODE\)
    else
        root -x -b -q -l CompileCorrectGamma.C
        root -x -l -b -q TaskV1/CorrectGamma.C\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$ENERGY\"\,\"\"\,\"$ESTIMATEPILEUP\"\,$MODE\)
    fi
}

function CorrectSignalGammaV2()
{
    if [ $ISROOT6 -eq 0 ]; then
        root -x -b -q -l CompileCorrectGammaV2.C\+\+
        root -x -l -b -q TaskV1/CorrectGammaV2.C\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$ENERGY\"\,\"\"\,\"$ESTIMATEPILEUP\"\,$MODE\,\"$7\"\)
    else
        root -x -b -q -l CompileCorrectGammaV2.C
        root -x -l -b -q TaskV1/CorrectGammaV2.C\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$ENERGY\"\,\"\"\,\"$ESTIMATEPILEUP\"\,$MODE\,\"$7\"\)
    fi
}

function CreateGammaFinalResultsV3()
{
    if [ $ISROOT6 -eq 0 ]; then
        root -x -l -b -q TaskV1/CalculateGammaToPi0V4.C\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$7\"\,\"$ENERGY\"\,$MODE\,\"$FILEWEIGHTGAMMA\"\)
    else
        root -x -l -b -q TaskV1/CalculateGammaToPi0V4.C\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$7\"\,\"$ENERGY\"\,$MODE\,\"$FILEWEIGHTGAMMA\"\)
    fi
}

function RunPi0Tagging()
{
    if [ $ISROOT6 -eq 0 ]; then
        root -x -l -b -q TaskV1/Pi0Tagging.C\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$7\"\,\"$ENERGY\"\,$MODE\)
    else
        root -x -l -b -q TaskV1/Pi0Tagging.C\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$7\"\,\"$ENERGY\"\,$MODE\)
    fi
}

function AskForTHnSparseOption()
{
    echo "Do you want to use THnSparse for the background? Yes/No?";
    read answer
    if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
        echo -e "--> Will use THnSparse for the background ...\n";
        USETHNSPARSE=1
    elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
        echo -e "--> Will NOT use THnSparse for the background ...\n";
        USETHNSPARSE=0
    else
        echo "--> Command \"$answer\" not found. Please try again."
    fi
}


function AskForMinBiasEffiOnly()
{
    echo "Do you want to use MinBias Efficiencies only? Yes/No?";
    read answer
    if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
        OPTMINBIASEFF="MinBiasEffOnly"
        echo -e "--> Calculating MinBias Efficiecy only ...\n";
        CORRECT=1
    elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
        OPTMINBIASEFF=""
        echo -e "--> Nothing to be done ...\n";
        CORRECT=1
    else
        echo "--> Command \"$answer\" not found. Please try again."
    fi
}

function Usage()
{
    echo -e "
This Script is provided by the Gamma Conversion Group of ALICE
the main developers are
    Friederike Bock \t friederike.bock@cern.ch
    Ana Marin \t\t marin@physi.uni-heidelberg.de

If there are any complications with running this script do not hesitate to connect them.


How to use this script?
$PROGNAME -h  \t\t\t\t\t usage will be displayed
$PROGNAME --help \t\t\t\t usage will be displayed
$PROGNAME -c data.root suffix \t\t\t Will only correct and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -d data.root suffix \t\t\t Will only fullfil cutstudies and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -m DirectoryOfMergedOutputs suffix \t Automatic merging for LHC10bc and LHC10de for efficiencies will be done \n \t\t\t\t\t\t\t\t\t according to the fraction in Data, be careful files have to be arranged in a certain structure \n
    \t\t\t\t\t\t\t\t DATAROOTFILE=DirectoryOfMergedOutputs/mergedALL/GammaConvV1Data.root
    \t\t\t\t\t\t\t\t MCROOTFILE=DirectoryOfMergedOutputs/mergedALL/GammaConvV1MC.root
    \t\t\t\t\t\t\t\t MCROOTFILEMERGEA=DirectoryOfMergedOutputs/mergedBC/GammaConvV1MC.root
    \t\t\t\t\t\t\t\t MCROOTFILEMERGEB=DirectoryOfMergedOutputs/mergedDE/GammaConvV1MC.root\n
$PROGNAME -r data.root suffix \t\t\t Will only execute the production of the final results and \n \t\t\t\t\t\t\t\t\t produce the graphical output in the specified suffix-format.\n
$PROGNAME data.root MC.root suffix \t \t Will execute Gamma Conversion Analysis for data.root file and MC.root file \n \t\t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME  *-*gammaOff* \t\t\t gamma calculation switched off \n
$PROGNAME  *-*gammaOnly* \t\t\t gamma calculation only\n
$PROGNAME  *-*pi0Only* \t\t\t\t pi0 calculation only\n
$PROGNAME  *-*pi0etaOnly* \t\t\t pi0 in eta binnin calculation only\n
$PROGNAME  *-*etaOnly* \t\t\t\t eta calculation only\n
$PROGNAME  *-*etaOff* \t\t\t\t eta calculation switched off\n

This script needs as a basis the output from the GammaConversion-Software which is provided with the Aliroot software, both the 'data.root' and the 'MC.root' have to be output files of this. The script then will check which cutnumbers are in these files and will ask you which exactly you want to analyse, furthermore you have to set a standard cut which has to be always the first in the CutLogFile. Because for this cut the systematic errors will be calculated. Not only one analysis will be fullfiled, you can additionally choose to do a Alpha studie or Chi2 of the meson as well which will give you the oportunity to calculate two indepent error due to cut variation. Additionally the CorrectSignal.C will correct the spectrum for all possible contributions and afterwards calculate the systematic error due to yield extraction. All these error will then enter the final results where you will have plots with only statistical errors as well as systematic + static errors. Several data output files are created, and
stored in the corresponding cut-directory or the working directory for the cutsstudies file.

    "
    exit
}

if [ "$#" == "0" ]; then
    echo -e "
$PROGNAME -h  \t\t\t\t\t usage will be displayed
$PROGNAME --help \t\t\t\t usage will be displayed
$PROGNAME -c data.root suffix \t\t\t Will only correct and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -d data.root suffix \t\t\t Will only fullfil cutstudies and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -m DirectoryOfMergedOutputs suffix \t Automatic merging for LHC10bc and LHC10de for efficiencies will be done \n \t\t\t\t\t\t\t\t\t according to the fraction in Data, be careful files have to be arranged in a certain structure \n
    \t\t\t\t\t\t\t\t DATAROOTFILE=DirectoryOfMergedOutputs/mergedALL/GammaConvV1Data.root
    \t\t\t\t\t\t\t\t MCROOTFILE=DirectoryOfMergedOutputs/mergedALL/GammaConvV1MC.root
    \t\t\t\t\t\t\t\t MCROOTFILEMERGEA=DirectoryOfMergedOutputs/mergedBC/GammaConvV1MC.root
    \t\t\t\t\t\t\t\t MCROOTFILEMERGEB=DirectoryOfMergedOutputs/mergedDE/GammaConvV1MC.root\n
$PROGNAME -r data.root suffix \t\t\t Will only execute the production of the final results and \n \t\t\t\t\t\t\t\t\t produce the graphical output in the specified suffix-format.\n
$PROGNAME data.root MC.root suffix \t \t Will execute Gamma Conversion Analysis for data.root file and MC.root file \n \t\t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME  *-*gammaOff* \t\t\t gamma calculation switched off \n
$PROGNAME  *-*gammaOnly* \t\t\t gamma calculation only\n
$PROGNAME  *-*pi0Only* \t\t\t\t pi0 calculation only\n
$PROGNAME  *-*pi0etaOnly* \t\t\t pi0 in eta binnin calculation only\n
$PROGNAME  *-*etaOnly* \t\t\t\t eta calculation only\n
$PROGNAME  *-*etaOff* \t\t\t\t eta calculation switched off\n
"
    exit
fi

echo $1
echo $2
echo $3
echo $4


if [[ "$ROOTVERSION" == *v6* ]]; then
    ISROOT6=1
    echo "root 6 version detected"
fi

if [[ "$1" == *-*gammaOff* ]]; then
    DOGAMMA=0
    echo "gamma calculation switched off"
fi
etaOnly=0


if [[ "$1" == *-*etaOff* ]]; then
    DOPI0=1
    DOETA=0
    DOPI0INETABINS=0
    echo "Eta calculation switched off"
fi

if [[ "$1" == *-*etaOnly* ]]; then
    DOPI0=0
    DOETA=1
    DOPI0INETABINS=0
    echo "Eta calculation only"
fi

if [[ "$1" == *-*pi0etaOnly* ]]; then
    DOPI0=0
    DOETA=0
    DOPI0INETABINS=1
    echo "Pi0 in Eta binning calculation only"
fi

if [[ "$1" == *-*pi0Only* ]]; then
    DOPI0=1
    DOETA=0
    DOPI0INETABINS=0
    echo "Pi0 calculation only"
fi

if [[ "$1" == *-*gammaOnly* ]]; then
    DOPI0=0
    DOETA=0
    DOPI0INETABINS=0
    DOGAMMA=1
    echo "gamma calculation only"
fi

if [[ "$1" == *-*ToyOff* ]]; then
    DISABLETOYMC=1
    echo "toy MC switched off"
fi


if [[ "$1" == *-*aPUC* ]]; then
    ADDPILEUPCORR=kTRUE
    echo "additionally pileup correction on"
fi

if [[ "$1" == *-h* ]] ; then
    Usage
elif [[ "$1" == *-mAddSig2760GeV* ]] ; then
    MERGINGMC=1
    DATAROOTFILE=$2
    MCROOTFILE=$3
    MCROOTFILEADDSIG=$3
    MCROOTFILEADDSIGETA=$3
    SUFFIX=$4
    ADDEDSIG=1
    echo ""
elif [[ "$1" == *-mMerged11a* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    CONFIG=$3
    SUFFIX=$4
    USETHNSPARSE=0
    minPtMergePi0=7
    minPtMergeEta=20
    DATAROOTFILE=$DIRECTORY/GammaCaloMerged_LHC11a-pass4_$CONFIG.root
    MCROOTFILE=$DIRECTORY/GammaCaloMerged_MC_LHC15g1aFinerPtHardBins_$CONFIG.root
    MCROOTFILEGJ=$DIRECTORY/GammaCaloMerged_MC_LHC15g1b_$CONFIG.root
elif [[ "$1" == *-mMerged13g* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    CONFIG=$3
    SUFFIX=$4
    USETHNSPARSE=0
    minPtMergePi0=7
    minPtMergeEta=20
    DATAROOTFILE=$DIRECTORY/GammaCaloMerged_LHC13g-pass1_$CONFIG.root
    MCROOTFILE=$DIRECTORY/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_$CONFIG.root
    MCROOTFILEGJ=$DIRECTORY/GammaCaloMerged_MC_LHC15a3b_$CONFIG.root
elif [[ "$1" == *-mJJ5TeVA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConv_Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConv_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConv_MC_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConv_MC_A.root
    minPtMergePi0=7
    minPtMergeEta=8
    ADDEDSIG=2
elif [[ "$1" == *-mAddSigPbPbLHC13d2A* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC13d2_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_MC_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_MC_A.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC13d2bA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC13d2b_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2b_MC_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2b_MC_A.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC13d2xA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_A.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC13d2B* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_B.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC13d2_MC_B.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_MC_B.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_MC_B.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC13d2bB* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_B.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC13d2b_MC_B.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2b_MC_B.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2b_MC_B.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC13d2xB* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_B.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC13d2_LHC13d2b_MC_B.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_B.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC13d2_LHC13d2b_MC_B.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC14a1aA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Pi0_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Eta_A.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC14a1bA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC14a1b_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Pi0_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Eta_A.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPbPbLHC14a1abA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC14a1a_LHC14a1b_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1a_LHC14a1b_MC_Pi0_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1a_LHC14a1b_MC_Eta_A.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigPi0ForGammaPbPbLHC14a1aA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1a_MC_Pi0_A.root
    ADDEDSIG=1
    DOPI0=1
    DOETA=0
    DOPI0INETABINS=0
    DOGAMMA=0
elif [[ "$1" == *-mAddSigPi0ForGammaPbPbLHC14a1bB* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_B.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC14a1b_MC_B.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_LHC14a1b_MC_Pi0_B.root
    ADDEDSIG=1
    DOPI0=1
    DOETA=0
    DOPI0INETABINS=0
    DOGAMMA=0
elif [[ "$1" == *-mAddSigForGammaPbPbLHC14a1aA* ]] ; then
    MERGINGMC=0
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC14a1a_MC_A.root
    ADDEDSIG=0
    DOPI0=0
    DOETA=0
    DOPI0INETABINS=0
    DOGAMMA=1
elif [[ "$1" == *-mAddSigForGammaPbPbLHC14a1bB* ]] ; then
    MERGINGMC=0
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_B.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_LHC14a1b_MC_B.root
    ADDEDSIG=0;
    DOPI0=0
    DOETA=0
    DOPI0INETABINS=0
    DOGAMMA=1
elif [[ "$1" == *-mAddSigpPbHIJINGA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_HIJING_MC_A.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_A.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_A.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigpPbHIJINGB* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_B.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_HIJING_MC_B.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_B.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_B.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigpPbHIJINGC* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1Data_C.root
    MCROOTFILE=$DIRECTORY/mergedMinBias/GammaConvV1_HIJING_MC_C.root
    MCROOTFILEADDSIG=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_C.root
    MCROOTFILEADDSIGETA=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_C.root
    ADDEDSIG=1
elif [[ "$1" == *-mAddSigpPb* ]] ; then
    MERGINGMC=1
    SUFFIX=$5
    DATAROOTFILE=$2
    MCROOTFILE=$3
    MCROOTFILEADDSIG=$3
    MCROOTFILEADDSIGETA=$4
    ADDEDSIG=1
elif [[ "$1" == *-mAddSig* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedAll/GammaConvV1Data.root
    MCROOTFILE=$DIRECTORY/mergedAll/GammaConvV1MC.root
    MCROOTFILEMERGEA=$DIRECTORY/WithoutAddedSignals/GammaConvV1MC.root
    MCROOTFILEMERGEB=$DIRECTORY/WithAddedSignals/GammaConvV1MC.root
elif [[ "$1" == *-m* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    SUFFIX=$3
    DATAROOTFILE=$DIRECTORY/mergedALL/GammaConvV1Data.root
    MCROOTFILE=$DIRECTORY/mergedALL/GammaConvV1MC.root
    MCROOTFILEMERGEA=$DIRECTORY/mergedBC/GammaConvV1MC.root
    MCROOTFILEMERGEB=$DIRECTORY/mergedDE/GammaConvV1MC.root
elif [[ "$1" == *-c* ]] ; then
    ONLYCORRECTION=1
    DATAROOTFILE=$2
    SUFFIX=$3;
elif [[ "$1" == -d* ]] ; then
    ONLYCORRECTION=1
    ONLYCUTS=1
    DATAROOTFILE=$2
    SUFFIX=$3;
elif [[ "$1" != -* ]] ; then
    DATAROOTFILE=$1
    MCROOTFILE=$2
    SUFFIX=$3;
else
    DATAROOTFILE=$2
    MCROOTFILE=$3
    SUFFIX=$4;
fi


# test if files have been set correctly
if [ -f $DATAROOTFILE ]; then
    DATAFILEOK=1
    echo "The data file specified is $DATAROOTFILE"
else
    echo "No data file specified, analysis can not be fullfiled."
#    exit
fi
if [ -f $MCROOTFILE ]; then
    echo "The MC file specified is $MCROOTFILE"
    if [ "$MCROOTFILEMERGEA" ]; then
        echo "The MC file for period merge A specified is $MCROOTFILEMERGEA"
    fi
    if [ "$MCROOTFILEMERGEB" ]; then
        echo "The MC file for period merge A specified is $MCROOTFILEMERGEB"
    fi
    if [ "$MCROOTFILEADDSIG" ]; then
           echo "The MC file specified for added signals is $MCROOTFILEADDSIG"
    fi
    if [ "$MCROOTFILEADDSIGETA" ]; then
        echo "The MC file specified for added signals eta only is $MCROOTFILEADDSIGETA"
    fi
    if [ "$MCROOTFILEGJ" ]; then
        echo "The MC file specified for GJ-MC is $MCROOTFILEGJ"
    fi
else
    echo "No MC file specified, analysis will only made paritally, please be careful with the results."
    MCFILE=0
fi


# Special case for only cutstudies overview running for systematics or other comparsions
if [ $ONLYCUTS -eq 1 ]; then
    CORRECT=0
    while [ $CORRECT -eq 0 ]
    do
    echo "Which name should your CutStudies have? None/*name*"
    read answer
    if [ $answer = "None" ]; then
        echo -e "--> No name for CutStudies selected...\n";
        CORRECT=1
    else
        NAMECUTSTUDIES=$answer
        echo -e "--> You have selcted the name: $NAMECUTSTUDIES\n";
        CORRECT=1
    fi
    done

    CORRECT=0
    while [ $CORRECT -eq 0 ]
    do
    echo "Are these CutStudies for a specific data period? No/*name*"
    read answer
    if [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
        CORRECT=1
    else
        PERIODNAME=$answer
        echo "The period name was set to: " $PERIODNAME;
        CORRECT=1
    fi
    done
fi

echo -e "\n____________________________"
echo -e "READING ANSWERS TO QUESTIONS\n"

# Inquire the mode
CORRECT=0
while [ $CORRECT -eq 0 ]
do
  echo "Which mode are you running?
        0 (PCM-PCM)
        1 (PCM-Dalitz)
        2 (PCM-EMCAL)
        3 (PCM-PHOS)
        4 (EMCAL-EMCAL)
        5 (PHOS-PHOS)
        9 (old files)
        10 (EMC-merged)
        11 (PHOS-merged)
        12 (DCal-DCal)
        13 (PCM-DCal)
        and add 100 for heavy mesons (e.g. 102 for PCM-EMCAL)"
    read answer
    # Set heavy meson configuration
    if [ $answer -ge 100 ]; then
        DOPI0=0
        DOETA=1
        DOETAPRIME=1
        DOPI0INETABINS=0
        echo "  Set options for EtaPrime mode (mode >100)"
        echo "    DOPI0=$DOPI0"
        echo "    DOETA=$DOETA"
        echo "    DOETAPRIME=$DOETAPRIME"
        echo "    DOPI0INETABINS=$DOPI0INETABINS"
    fi
    if [ $answer = "0" ] || [ $answer = "100" ]; then
        echo -e "--> You are analysing PCM-PCM output\n";
        NEVTSTOY=1e7
        MINPTTOY=0
        MAXPTTOY=50
        MODE=$answer
        CORRECT=1
    elif [ $answer = "1" ] || [ $answer = "101" ]; then
        echo "--> You are trying to analyse PCM-Dalitz output, this is the wrong script, please use another one.";
        MODE=$answer
        NEVTSTOY=1e7
        MINPTTOY=0
        MAXPTTOY=50
        CORRECT=0
    elif [ $answer = "2" ] || [ $answer = "102" ]; then
        echo -e "--> You are analysing PCM-EMCAL output\n";
        MODE=$answer
#        ADVMESONQA="AdvancedMesonQA"
        CORRECT=1
    elif [ $answer = "3" ] || [ $answer = "103" ]; then
        echo -e "--> You are analysing PCM-PHOS output\n";
        MODE=$answer
        ADVMESONQA="AdvancedMesonQA"
        CORRECT=1
    elif [ $answer = "4" ] || [ $answer = "104" ]; then
        echo -e "--> You are analysing EMCAL-EMCAL output\n";
        MODE=$answer
        ADVMESONQA="AdvancedMesonQA"
        CORRECT=1
    elif [ $answer = "5" ] || [ $answer = "105" ]; then
        echo -e "--> You are analysing PHOS-PHOS output\n";
        MODE=$answer
        ADVMESONQA="AdvancedMesonQA"
        CORRECT=1
    elif [ $answer = "10" ] || [ $answer = "110" ]; then
        echo -e "--> You are analysing EMC-merged output\n";
        MODE=$answer
        CORRECT=1
        DOETA=0;
        DOPI0INETABINS=0;
    elif [ $answer = "11" ] || [ $answer = "111" ]; then
        echo -e "--> You are analysing PHOS-merged output\n";
        MODE=$answer
        CORRECT=1
        DOETA=0;
        DOPI0INETABINS=0;
    elif [ $answer = "9" ] || [ $answer = "109" ]; then
        echo -e "--> You are analysing the old output of PCM-PCM\n";
        MODE=$answer
        CORRECT=1
    elif [ $answer = "12" ] || [ $answer = "112" ]; then
        echo -e "--> You are analysing DCAL-DCAL output\n";
        MODE=$answer
        ADVMESONQA=""#"AdvancedMesonQA"
        CORRECT=1
    elif [ $answer = "13" ] || [ $answer = "113" ]; then
        echo -e "--> You are analysing PCM-DCAL output\n";
        MODE=$answer
        ADVMESONQA=""#"AdvancedMesonQA"
        CORRECT=1
    else
        echo -e "--> Command \"$answer\" not found. Please try again.\n"
    fi
done

# Inquire cut
CORRECT=0
while [ $CORRECT -eq 0 ]
do
    echo "Do you want to take an already exitsting CutSelection.log-file. Yes/No"
    read answer
    if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
        echo -e "--> Chosen already existing logfile ...\n"
        cat CutSelection.log
        CORRECT=1
    elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
        echo -e "--> Generating new CutSelection.log file ...\n"
        if [ $DATAFILEOK -eq 1 ] ; then
            root -b -q -x -l TaskV1/MakeCutLog.C\(\"$DATAROOTFILE\"\,\"CutSelection.log\"\,$MODE\)
            CORRECT=1
        else
            root -b -q -x -l TaskV1/MakeCutLog.C\(\"$MCROOTFILE\"\,\"CutSelection.log\"\,$MODE\)
            CORRECT=1
        fi
    else
        echo "--> Command \"$answer\" not found. Please try again."
    fi
done

# Inquire collision system
CORRECT=0
while [ $CORRECT -eq 0 ]
do
    echo -e "\nWhich collision system do you want to process?
          pp Systems:
          \t 900GeV (pp@900GeV), 2.76TeV (pp@2.76TeV), 5TeV (pp@5.02TeV), 5TeV2017 (2017 pp@5.02TeV), 7TeV (pp@7TeV),  8TeV (pp@8TeV), 13TeV (pp@13TeV), 13TeVLowB (pp@13TeV), 13TeVRBins (pp@13TeV)
          pPb Systems:
          \t pPb_5.023TeV (pPb@5.023TeV), pPb_8TeV (pPb@8TeV)
          A-A Systems:
          \t PbPb_2.76TeV (PbPb@2.76TeV), PbPb_5.02TeV (PbPb@5.02TeV), XeXe_5.44TeV(XeXe@5.44TeV)"
    read answer
    if [ $answer = "900GeV" ] || [ $answer = "900" ] || [ $answer = "9" ] || [ $answer = "0.9" ]; then
        ENERGY="900GeV";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "2.76TeV" ] || [ $answer = "2" ] || [ $answer = "2.76" ]; then
        ENERGY="2.76TeV";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "5TeV" ] || [ $answer = "5.02TeV" ] || [ $answer = "5" ] || [ $answer = "5.02" ]; then
        ENERGY="5TeV";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "5TeV2017" ] || [ $answer = "5.02TeV2017" ] || [ $answer = "52" ] || [ $answer = "5.022017" ]; then
        ENERGY="5TeV2017";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "7TeV" ] || [ $answer = "7" ]; then
        ENERGY="7TeV";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "8TeV" ] || [ $answer = "8" ]; then
        ENERGY="8TeV";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "13TeV" ] || [ $answer = "13" ]; then
        ENERGY="13TeV";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "13TeVRBins" ] ; then
        ENERGY="13TeVRBins";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "13TeVLowB" ]; then
        ENERGY="13TeVLowB";
        EXTINPUTFILE="ExternalInput/IdentifiedCharged/ChargedIdentifiedSpectraPP_2016_08_14.root";
    elif [ $answer = "PbPb_2.76TeV" ] || [ $answer = "PbPb_2.76" ] || [ $answer = "PbPb2" ] || [ $answer = "Pb2" ]; then
        ENERGY="PbPb_2.76TeV";
        EXTINPUTFILE="";
    elif [ $answer = "PbPb_5.02TeV" ] || [ $answer = "PbPb_5.02" ] || [ $answer = "PbPb5" ] || [ $answer = "Pb5" ]; then
        ENERGY="PbPb_5.02TeV";
        EXTINPUTFILE="";
    elif [ $answer = "XeXe_5.44TeV" ] || [ $answer = "XeXe_5.44" ] || [ $answer = "XeXe5" ] || [ $answer = "Xe5" ]; then
        ENERGY="XeXe_5.44TeV";
        EXTINPUTFILE="";
    elif [ $answer = "pPb_5.023TeV" ] || [ $answer = "pPb_5.023" ] || [ $answer = "pPb5" ];  then
        ENERGY="pPb_5.023TeV";
        EXTINPUTFILE="";
    elif [ $answer = "pPb_5.023TeVCent" ] || [ $answer = "pPb_5.023TeVCent" ] || [ $answer = "pPb5MB" ];  then
        ENERGY="pPb_5.023TeVCent";
        EXTINPUTFILE="";
    elif [ $answer = "pPb_5.023TeVRun2" ] || [ $answer = "pPb_5.023R2" ] || [ $answer = "pPb5R2" ];  then
        ENERGY="pPb_5.023TeVRun2";
        EXTINPUTFILE="";
    elif [ $answer = "pPb_8TeV" ] || [ $answer = "pPb_8" ] || [ $answer = "pPb8" ];  then
        ENERGY="pPb_8TeV";
        EXTINPUTFILE="";
    fi
    echo -e "--> The collision system has been selected to be $ENERGY.\n"

    echo "Is a cocktail file available? Yes/No?"
    read answer
    if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
        echo "Please enter the filepath of the cocktail file."
            read COCKROOTFILE
            if [ -f $COCKROOTFILE ]; then
                echo -e "--> The cocktail file specified is $COCKROOTFILE\n"
                USECOCK=1
                    echo "Please enter the rapidity used in the cocktail, e.g. 0.80"
                    read COCKRAP
                    echo -e "--> Rapidity of $COCKRAP has been chosen.\n"
            else
                echo -e "--> No cocktail file specified, it will not be used.\n"
                USECOCK=0
            fi
    elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
        echo -e "--> Will not use cocktail input for secondary correction or double ratio.\n"
        USECOCK=0
    fi

    #######################################################################################################
    # Set Special variable for direct photon running
    #######################################################################################################
    echo "Do you want to run only direct photon studies? Yes/No?";
    read answer
    if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
        DIRECTPHOTON="directPhoton"
        DOPI0=1
        DOETA=0
        DOPI0INETABINS=0
        DOETAPRIME=0
        DOGAMMA=1
        echo -e "--> Switching off eta, pi0 in eta binnings ...\n";

        if [ $MODE = 2 ] || [ $MODE = 3 ]; then
            echo "Do you want to run Pi0-tagging instead of DR? Yes/No?";
            read answer
            if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
                echo -e "--> Running pi0-tagging ...\n";
                DIRECTPHOTON="directPhotonTagging"
                DOPI0TAG=1
                DOPI0=0
                DOGAMMA=1
            elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
                echo -e "--> Running standard DR ...\n";
                DOPI0TAG=0
            else
                echo "--> Command \"$answer\" not found. Please try again."
            fi
        fi
    elif [ $answer = "YesPCMEMC" ] || [ $answer = "YPCMEMC" ] || [ $answer = "yPCMEMC" ] || [ $answer = "yesPCMEMC" ]; then
        echo -e "--> Will produce Direct Photon plots with special PCMEMC binning...\n";
        DIRECTPHOTON="directPhotonA"
        DOPI0=1
        DOETA=0
        DOPI0INETABINS=0
        DOETAPRIME=0
        DOGAMMA=1
        echo -e "--> Switching off eta, pi0 in eta binnings ...\n";
    elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
        echo -e "--> No special modifications have been done to the settings.\n"
    else
        echo "--> Command \"$answer\" not found. Please try again."
    fi

    #######################################################################################################
    # Set TNSPARSE off for Calo modes by default
    #######################################################################################################
    if [ $MODE = 2 ] || [ $MODE = 3 ] || [ $MODE = 4 ] || [ $MODE = 5 ] || [ $MODE = 10 ] || [ $MODE = 11 ] || [ $MODE = 12 ] || [ $MODE = 13 ]; then
        USETHNSPARSE=0
        ADVMESONQA="AdvancedMesonQA"
    fi

    #######################################################################################################
    # Obtain desired number of bins for different energies
    #######################################################################################################
    if [ $ONLYCORRECTION -eq 0 ]; then
        # PP systems
        if [ $ENERGY = "900GeV" ]; then
            GiveBinning900GeV
        elif [ $ENERGY = "2.76TeV" ]; then
            # hack to run dir gamma in parallel if wanted
            if [ $MODE -eq 4 ]  || [ $MODE -eq 12 ] ; then
                DIRECTPHOTON="No"
            else
                DIRECTPHOTON="Gamma"
            fi
            GiveBinning2760GeV
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi
        elif [ $ENERGY = "5TeV" ]; then
            GiveBinning5TeV
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi
        elif [ $ENERGY = "5TeV2017" ]; then
            GiveBinning5TeV2017
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi
        elif [ $ENERGY = "7TeV" ]; then
            GiveBinning7TeV
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi
        elif [ $ENERGY = "8TeV" ]; then
            GiveBinning8TeV
        elif [ $ENERGY = "13TeV" ]; then
            GiveBinning13TeV
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi

        ## PPb systems
        elif [ $ENERGY = "pPb_5.023TeV" ] || [ $ENERGY = "pPb_5.023TeVCent" ] || [ $ENERGY = "pPb_5.023TeVRun2" ]  ; then
            GiveBinningpPb5
            AskForMinBiasEffiOnly
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi
        elif [ $ENERGY = "pPb_8TeV" ]  ; then
            GiveBinningpPb8
            AskForMinBiasEffiOnly
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi

        ## A-A systems
        elif [ $ENERGY = "PbPb_2.76TeV" ]; then
            GiveBinningPbPb2760
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi
        elif [ $ENERGY = "PbPb_5.02TeV" ]; then
            GiveBinningPbPb5TeV
            if [ $MODE -lt 10 ] || [ $MODE -gt 11 ] ; then
                AskForTHnSparseOption
            fi
        elif [ $ENERGY = "XeXe_5.44TeV" ]; then
            GiveBinningXeXe5440GeV
        fi

    fi
    CORRECT=1

done

if [ $MODE -eq 0 ] || [ $MODE -eq 9 ] ; then
    CORRECT=0
    while [ $CORRECT -eq 0 ]
    do
        echo "Please don't forget to run the TaskV1/AnalyseDCADist.C separately to estimate train pileup?";
        echo "Do the output of this macro exist already? Yes/No?";
        read answer
        if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
            ESTIMATEPILEUP="EstimateTrainPileUp"
            echo -e "--> Running with additional histos ...\n";
            CORRECT=1
        elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
            ESTIMATEPILEUP=""
            echo -e "--> Running in Normal mode ...\n";
            CORRECT=1
        else
            echo "--> Command \"$answer\" not found. Please try again."
        fi
    done
fi

echo "Checking if mode $MODE is standard mode...";
if [ $MODE -lt 10 ]  || [ $MODE = 12 ] ||  [ $MODE = 13 ] || [ $MODE -ge 100 ]; then
    echo -e "--> I went into standard modes\n";
    if [ $ONLYCORRECTION -eq 0 ];  then
        # echo "Extraction will be done using modified Gaussian.";
        # crystal=Gaussian

        CORRECT=0
        while [ $CORRECT -eq 0 ]
        do
            echo "Which fit do you want to do? CrystalBall or gaussian convoluted with an exponential function? CrystalBall/Gaussian?";
            read answer
            if [ $answer = "CrystalBall" ] || [ $answer = "C" ] || [ $answer = "c" ]; then
                echo -e "--> CrystalBall chosen ...\n";
                CORRECT=1
                crystal=CrystalBall
            elif [ $answer = "Gaussian" ] || [ $answer = "G" ] || [ $answer = "g" ]; then
                echo -e "--> Gaussian chosen ...\n";
                CORRECT=1
                crystal=Gaussian
            else
                echo "--> Command \"$answer\" not found. Please try again."
            fi
        done
    fi

    # echo "Hauptroutine stimmt"
    CORRECT=0
    while [ $CORRECT -eq 0 ]
    do
        echo "Please check that you really want to process all cuts, otherwise change the CutSelection.log. Remember at first all gamma cutstudies will be carried out. Make sure that the standard cut is the first in the file. Continue? Yes/No?";
        read answer
        if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
            echo -e "--> Continuing ...\n";
            CORRECT=1
        elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
            echo -e "--> Aborting ...\n";
            exit
        else
            echo "--> Command \"$answer\" not found. Please try again."
        fi
    done


# Read the different cuts form the Cut selection log file
    CutSelections=`cat CutSelection.log`
    for CUTSELECTION in $CutSelections; do
        if [ -d $CUTSELECTION ]; then
            echo "CutSelection $CUTSELECTION directory already exists, all files will be overwritten ";
            mkdir $CUTSELECTION/$ENERGY
        else
            mkdir $CUTSELECTION
            mkdir $CUTSELECTION/$ENERGY
        fi

        if [ $ONLYCUTS -eq 0 ]; then
            if [ -d $CUTSELECTION/$ENERGY/$SUFFIX ]; then
                echo "Graphical Output $SUFFIX directory already exists, all files will be overwritten ";
            else
                mkdir $CUTSELECTION/$ENERGY/$SUFFIX
            fi

            if [ $DISABLETOYMC -eq 0 ] && [ $USECOCK -eq 0 ] && [ $ONLYCORRECTION -eq 0 ]; then
                rm ToyMCOutputs.txt
                if [ $ISROOT6 -eq 0 ]; then
                    root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\+\($NEVTSTOY,0,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                    root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\+\($NEVTSTOY,1,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                    root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\+\($NEVTSTOY,2,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                else
                    root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\($NEVTSTOY,0,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                    root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\($NEVTSTOY,1,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                    root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\($NEVTSTOY,2,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                fi
            fi
            if [ $USECOCK -eq 1 ] && [ $ONLYCORRECTION -eq 0 ]; then
                root -b -x -l -q TaskV1/PrepareSecondaries.C\+\(\"Pi0\"\,\"$COCKROOTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,\"$COCKRAP\"\,\"\"\,$BINSPTPI0\,$MODE,kFALSE\)
            fi

            if [ $ONLYCORRECTION -eq 0 ]; then
                echo "CutSelection is $CUTSELECTION";
                if [ $DOPI0 -eq 1 ]; then
                    echo -e "\n\n_________________________"
                    echo -e "EXTRACTING SIGNAL FOR PI0\n"
                    OPTIONSPI0DATA=\"Pi0\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
                    echo $DATAROOTFILE
                    if [ -f $DATAROOTFILE ]; then
                        ExtractSignal $OPTIONSPI0DATA
                    fi
                    PI0DATARAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_data_GammaConvV1WithoutCorrection_*.root`
                    if [ $MCFILE -eq 1 ]; then
                        OPTIONSPI0MC=\"Pi0\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
                        ExtractSignal $OPTIONSPI0MC
                        PI0MCRAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION.root`
                        PI0MCCORRFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION.root`
                        if [ $MERGINGMC -eq 1 ]; then
                            if [ $ADDEDSIG -eq 1 ]; then
                                OPTIONSPI0MC2=\"Pi0\"\,\"$MCROOTFILEADDSIG\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"AddSig\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kTRUE
                                ExtractSignal $OPTIONSPI0MC2
                                PI0MCCORRECTION=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistos_*.root`
                                PI0MCCORRECTIONADDSIG=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosAddSig_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\+\(\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0MCCORRECTION\"\,\"$CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0MCCORRECTIONADDSIG\"\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\(\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0MCCORRECTION\"\,\"$CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0MCCORRECTIONADDSIG\"\)
                                fi
                            elif [ $ADDEDSIG -eq 2 ]; then
                                OPTIONSPI0MC2=\"Pi0\"\,\"$MCROOTFILEADDSIG\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"JetJetMC\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
                                ExtractSignal $OPTIONSPI0MC2
                                PI0MCCORRECTION=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistos_*.root`
                                PI0MCCORRECTIONADDSIG=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosJetJetMC_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\+\(\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0MCCORRECTION\"\,\"$CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0MCCORRECTIONADDSIG\"\,$minPtMergePi0\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\(\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0MCCORRECTION\"\,\"$CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0MCCORRECTIONADDSIG\"\,$minPtMergePi0\)
                                fi
                            else
                                OPTIONSPI0MCA=\"Pi0\"\,\"$MCROOTFILEMERGEA\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"A\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
                                echo $OPTIONSPI0MCA
                                ExtractSignal $OPTIONSPI0MCA
                                OPTIONSPI0MCB=\"Pi0\"\,\"$MCROOTFILEMERGEB\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"B\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,kFALSE
                                echo $OPTIONSPI0MCB
                                ExtractSignal $OPTIONSPI0MCB
                                PI0MCCORRECTIONFILEA=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosA_*.root`
                                PI0MCCORRECTIONFILEB=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistosB_*.root`
                                root -x -q -l -b  TaskV1/MergeEffiWithProperWeighting.C\(\"$DIRECTORY\"\,\"$CUTSELECTION\"\,\"Pi0\",\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0MCCORRFILE\"\,\"$Pi0MCCORRECTIONBCFILE\",\"$Pi0MCCORRECTIONDFILE\"\)
                            fi
                        fi

                        if [ $ISROOT6 -eq 0 ]; then
                            root -x -l -b -q TaskV1/CompareMesonQuantities.C\+\(\"$PI0DATARAWFILE\"\,\"$PI0MCRAWFILE\"\,\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,$BINSPTPI0\,$MODE\)
                        else
                            root -x -l -b -q TaskV1/CompareMesonQuantities.C\(\"$PI0DATARAWFILE\"\,\"$PI0MCRAWFILE\"\,\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,$BINSPTPI0\,$MODE\)
                        fi
                    fi
                fi
                if [ $DOGAMMA -eq 1 ]; then
                    echo -e "\n\n_________________________"
                    echo -e "EXTRACTING SIGNAL FOR GAMMA\n"
                    OPTIONSPI0DATA=\"Pi0\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ESTIMATEPILEUP\"\,$BINSPTPI0\,kFALSE
                    if [ $NEWGAMMAMACROS == 0 ]; then
                        ExtractSignalGamma $OPTIONSPI0DATA
                    else
                        OPTIONSGAMMADATA=\"Pi0\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,\"\"\,$BINSPTGAMMA\,kFALSE\,$MODE
                        ExtractSignalGammaV2 $OPTIONSGAMMADATA $OPTMINBIASEFF
                    fi
                    if [ $MCFILE -eq 1 ]; then
                        OPTIONSPI0MC=\"Pi0\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ESTIMATEPILEUP\"\,$BINSPTPI0
                        if [ $NEWGAMMAMACROS == 0 ]; then
                            ExtractSignalGamma $OPTIONSPI0MC
                        else
                            OPTIONSGAMMAMC=\"Pi0\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,\"\"\,$BINSPTGAMMA\,kFALSE\,$MODE
                            ExtractSignalGammaV2 $OPTIONSGAMMAMC $OPTMINBIASEFF
                        fi
                    fi
                fi
                if [ $DOPI0INETABINS -eq 1 ]; then
                    echo -e "\n\n____________________________"
                    echo -e "EXTRACTING SIGNAL FOR PI0-ETA\n"
                    if [ -f $DATAROOTFILE ]; then
                        OPTIONSPI0ETADATA=\"Pi0EtaBinning\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                        ExtractSignal $OPTIONSPI0ETADATA
                    fi
                    PI0ETADATARAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_data_GammaConvV1WithoutCorrection_$CUTSELECTION.root`
                    if [ $MCFILE -eq 1 ]; then
                        OPTIONSPI0ETAMC=\"Pi0EtaBinning\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                        ExtractSignal $OPTIONSPI0ETAMC
                        PI0ETAMCRAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION.root`
                        PI0ETAMCCORRFILE=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION.root`
                        if [ $MERGINGMC -eq 1 ]; then
                            if [ $ADDEDSIG -eq 1 ]; then
                                OPTIONSPI0ETAMC2=\"Pi0EtaBinning\"\,\"$MCROOTFILEADDSIG\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"AddSig\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kTRUE
                                ExtractSignal $OPTIONSPI0ETAMC2
                                PI0ETAMCCORR=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_*.root`
                                PI0ETAMCCORRADDSIG=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosAddSig_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\+\(\"$CUTSELECTION\"\,\"Pi0EtaBinning\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0ETAMCCORRADDSIG\"\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\(\"$CUTSELECTION\"\,\"Pi0EtaBinning\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0ETAMCCORRADDSIG\"\)
                                fi
                            elif [ $ADDEDSIG -eq 2 ]; then
                                OPTIONSPI0ETAMC2=\"Pi0EtaBinning\"\,\"$MCROOTFILEADDSIG\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"JetJetMC\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                                ExtractSignal $OPTIONSPI0ETAMC2
                                PI0ETAMCCORR=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_*.root`
                                PI0ETAMCCORRADDSIG=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosJetJetMC_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\+\(\"$CUTSELECTION\"\,\"Pi0EtaBinning\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0ETAMCCORRADDSIG\"\,$minPtMergePi0\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\(\"$CUTSELECTION\"\,\"Pi0EtaBinning\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$PI0ETAMCCORRADDSIG\"\,$minPtMergePi0\)
                                fi
                            else
                                OPTIONSPI0ETAMCA=\"Pi0EtaBinning\"\,\"$MCROOTFILEMERGEA\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"BC\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                                echo $OPTIONSPI0ETAMCA
                                ExtractSignal $OPTIONSPI0ETAMCA
                                OPTIONSPI0ETAMCB=\"Pi0EtaBinning\"\,\"$MCROOTFILEMERGEB\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"D\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                                echo $OPTIONSPI0ETAMCB
                                ExtractSignal $OPTIONSPI0ETAMCB

                                PI0ETAMCCORRFILEA=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosBC_*.root`
                                PI0ETAMCCORRFILEB=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosD_*.root`
                                root -x -q -l -b  TaskV1/MergeEffiWithProperWeighting.C\(\"$DIRECTORY\"\,\"$CUTSELECTION\"\,\"Pi0EtaBinning\",\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0ETAMCCORRFILE\"\,\"$PI0ETAMCCORRFILEA\",\"$PI0ETAMCCORRFILEB\"\)
                            fi
                        fi

                        echo -e "\n________________________"
                        echo -e "COMPARE MESON QUANTITIES\n"
                        root -x -l -b -q TaskV1/CompareMesonQuantities.C\+\(\"$PI0ETADATARAWFILE\"\,\"$PI0ETAMCRAWFILE\"\,\"$CUTSELECTION\"\,\"Pi0EtaBinning\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,$BINSPTETA\,$MODE\)
                    fi
                fi
                if [ $DOETA -eq 1 ]; then
                    echo -e "\n\n_________________________"
                    echo -e "EXTRACTING SIGNAL FOR ETA\n"
                    if [ -f $DATAROOTFILE ]; then
                        OPTIONSETADATA=\"Eta\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                        ExtractSignal $OPTIONSETADATA
                    fi
                    ETADATARAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_data_GammaConvV1WithoutCorrection_*$CUTSELECTION.root`
                    if [ $MCFILE -eq 1 ]; then
                        OPTIONSETAMC=\"Eta\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                        ExtractSignal $OPTIONSETAMC
                        ETAMCRAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION.root`
                        ETAMCCORRFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION.root`

                        if [ $MERGINGMC -eq 1 ]; then
                            if [ $ADDEDSIG -eq 1 ]; then
                                OPTIONSETAMC2=\"Eta\"\,\"$MCROOTFILEADDSIGETA\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"AddSig\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kTRUE
                                ExtractSignal $OPTIONSETAMC2
                                ETAMCCORR=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistos_*.root`
                                ETAMCCORRADDSIG=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosAddSig_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\+\(\"$CUTSELECTION\"\,\"Eta\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAMCCORRADDSIG\"\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\(\"$CUTSELECTION\"\,\"Eta\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAMCCORRADDSIG\"\)
                                fi
                            elif [ $ADDEDSIG -eq 2 ]; then
                                OPTIONSETAMC2=\"Eta\"\,\"$MCROOTFILEADDSIG\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"JetJetMC\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                                ExtractSignal $OPTIONSETAMC2
                                ETAMCCORR=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistos_*.root`
                                ETAMCCORRADDSIG=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosJetJetMC_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\+\(\"$CUTSELECTION\"\,\"Eta\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAMCCORRADDSIG\"\,$minPtMergeEta\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\(\"$CUTSELECTION\"\,\"Eta\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAMCCORR\"\,\"$CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAMCCORRADDSIG\"\,$minPtMergeEta\)
                                fi
                            else
                                OPTIONSETAMCA=\"Eta\"\,\"$MCROOTFILEMERGEA\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"A\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                                ExtractSignal $OPTIONSETAMCA
                                OPTIONSETAMCB=\"Eta\"\,\"$MCROOTFILEMERGEB\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"B\"\,\"$ADVMESONQA\"\,$BINSPTETA\,kFALSE
                                ExtractSignal $OPTIONSETAMCB
                                ETAMCCORRFILEA=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosA_*.root`
                                ETAMCCORRFILEB=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistosB_*.root`

                                root -x -q -l -b  TaskV1/MergeEffiWithProperWeighting.C\(\"$DIRECTORY\"\,\"$CUTSELECTION\"\,\"Eta\",\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAMCCORRFILE\"\,\"$ETAMCCORRFILEA\",\"$ETAMCCORRFILEB\"\)
                            fi
                        fi
                        if [ $ISROOT6 -eq 0 ]; then
                            root -x -l -b -q TaskV1/CompareMesonQuantities.C\+\(\"$ETADATARAWFILE\"\,\"$ETAMCRAWFILE\"\,\"$CUTSELECTION\"\,\"Eta\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,$BINSPTETA\,$MODE\)
                        else
                            root -x -l -b -q TaskV1/CompareMesonQuantities.C\(\"$ETADATARAWFILE\"\,\"$ETAMCRAWFILE\"\,\"$CUTSELECTION\"\,\"Eta\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,$BINSPTETA\,$MODE\)
                        fi
                    fi
                fi
                if [ $DOETAPRIME -eq 1 ]; then
                    echo -e "\n\n_______________________________"
                    echo -e "EXTRACTING SIGNAL FOR ETA PRIME\n"
                    if [ -f $DATAROOTFILE ]; then
                        OPTIONSETAPRIMEDATA=\"EtaPrime\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETAPRIME\,kFALSE
                        ExtractSignal $OPTIONSETAPRIMEDATA
                    fi
                    ETAPRIMEDATARAWFILE=`ls $CUTSELECTION/$ENERGY/EtaPrime_data_GammaConvV1WithoutCorrection_*$CUTSELECTION.root`
                    if [ $MCFILE -eq 1 ]; then
                        OPTIONSETAPRIMEMC=\"EtaPrime\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETAPRIME\,kFALSE
                        ExtractSignal $OPTIONSETAPRIMEMC
                        ETAPRIMEMCRAWFILE=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION.root`
                        ETAPRIMEMCCORRFILE=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION.root`

                        if [ $MERGINGMC -eq 1 ]; then
                            if [ $ADDEDSIG -eq 1 ]; then
                                OPTIONSETAPRIMEMC2=\"EtaPrime\"\,\"$MCROOTFILEADDSIGETAPrime\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"AddSig\"\,\"$ADVMESONQA\"\,$BINSPTETAPRIME\,kTRUE
                                ExtractSignal $OPTIONSETAPRIMEMC2
                                ETAPRIMEMCCORR=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistos_*.root`
                                ETAPRIMEMCCORRADDSIG=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosAddSig_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\+\(\"$CUTSELECTION\"\,\"EtaPrime\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAPRIMEMCCORR\"\,\"$CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAPRIMEMCCORRADDSIG\"\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiWithProperWeighting2760GeV.C\(\"$CUTSELECTION\"\,\"EtaPrime\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAPRIMEMCCORR\"\,\"$CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAPRIMEMCCORRADDSIG\"\)
                                fi
                            elif [ $ADDEDSIG -eq 2 ]; then
                                OPTIONSETAPRIMEMC2=\"EtaPrime\"\,\"$MCROOTFILEADDSIG\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"JetJetMC\"\,\"$ADVMESONQA\"\,$BINSPTETAPRIME\,kFALSE
                                ExtractSignal $OPTIONSETAPRIMEMC2
                                ETAPRIMEMCCORR=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistos_*.root`
                                ETAPRIMEMCCORRADDSIG=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosJetJetMC_*.root`
                                if [ $ISROOT6 -eq 0 ]; then
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\+\(\"$CUTSELECTION\"\,\"EtaPrime\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAPRIMEMCCORR\"\,\"$CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAPRIMEMCCORRADDSIG\"\,$minPtMergeEtaPrime\)
                                else
                                    root -b -x -q -l TaskV1/MergeEffiJetJetMC.C\(\"$CUTSELECTION\"\,\"EtaPrime\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAPRIMEMCCORR\"\,\"$CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosMinBias_$CUTSELECTION.root\"\,\"$ETAPRIMEMCCORRADDSIG\"\,$minPtMergeEtaPrime\)
                                fi
                            else
                                OPTIONSETAPRIMEMCA=\"EtaPrime\"\,\"$MCROOTFILEMERGEA\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"A\"\,\"$ADVMESONQA\"\,$BINSPTETAPRIME\,kFALSE
                                ExtractSignal $OPTIONSETAPRIMEMCA
                                OPTIONSETAPRIMEMCB=\"EtaPrime\"\,\"$MCROOTFILEMERGEB\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$crystal\"\,\"$DIRECTPHOTON\"\,\"$OPTMINBIASEFF\"\,\"B\"\,\"$ADVMESONQA\"\,$BINSPTETAPRIME\,kFALSE
                                ExtractSignal $OPTIONSETAPRIMEMCB
                                ETAPRIMEMCCORRFILEA=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosA_*.root`
                                ETAPRIMEMCCORRFILEB=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistosB_*.root`

                                root -x -q -l -b  TaskV1/MergeEffiWithProperWeighting.C\(\"$DIRECTORY\"\,\"$CUTSELECTION\"\,\"EtaPrime\",\"$SUFFIX\"\,\"$ENERGY\"\,\"$ETAPRIMEMCCORRFILE\"\,\"$ETAPRIMEMCCORRFILEA\",\"$ETAPRIMEMCCORRFILEB\"\)
                            fi
                        fi
                        if [ $ISROOT6 -eq 0 ]; then
                            root -x -l -b -q TaskV1/CompareMesonQuantities.C\+\(\"$ETAPRIMEDATARAWFILE\"\,\"$ETAPRIMEMCRAWFILE\"\,\"$CUTSELECTION\"\,\"EtaPrime\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,$BINSPTETAPRIME\,$MODE\)
                        else
                            root -x -l -b -q TaskV1/CompareMesonQuantities.C\(\"$ETAPRIMEDATARAWFILE\"\,\"$ETAPRIMEMCRAWFILE\"\,\"$CUTSELECTION\"\,\"EtaPrime\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,$BINSPTETAPRIME\,$MODE\)
                        fi
                    fi
                fi
                if [ "$SUFFIX" == "pdf" ]; then
                    rm -f $CUTSELECTION/$ENERGY/$SUFFIX/ExtractSignal/ExtractSignal_all.pdf
                    pdfunite $CUTSELECTION/$ENERGY/$SUFFIX/ExtractSignal/*.pdf $CUTSELECTION/$ENERGY/$SUFFIX/ExtractSignal/*/*.pdf $CUTSELECTION/$ENERGY/$SUFFIX/ExtractSignal/ExtractSignal_all.pdf
                fi
            fi

            PI0DATARAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_data_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            PI0MCRAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            PI0MCCORRFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION*.root`
            PI0ETADATARAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_data_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            PI0ETAMCRAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            PI0ETAMCCORRFILE=`ls $CUTSELECTION/$ENERGY/Pi0EtaBinning_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION*.root`
            ETADATARAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_data_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            ETAMCRAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            ETAMCCORRFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION*.root`
            ETAPRIMEDATARAWFILE=`ls $CUTSELECTION/$ENERGY/EtaPrime_data_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            ETAPRIMEMCRAWFILE=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1WithoutCorrection_*$CUTSELECTION*.root`
            ETAPRIMEMCCORRFILE=`ls $CUTSELECTION/$ENERGY/EtaPrime_MC_GammaConvV1CorrectionHistos_*$CUTSELECTION*.root`

            if [ $DOPI0 -eq 1 ]; then
                if [ -f $PI0DATARAWFILE ] && [ -f $PI0MCCORRFILE ]; then
                    CorrectSignal $PI0DATARAWFILE $PI0MCCORRFILE $CUTSELECTION $SUFFIX Pi0 kFALSE $ESTIMATEPILEUP $DIRECTPHOTON
                fi
                if [ -f $PI0MCRAWFILE ] && [ -f $PI0MCCORRFILE ]; then
                    CorrectSignal $PI0MCRAWFILE $PI0MCCORRFILE $CUTSELECTION $SUFFIX Pi0 kTRUE $ESTIMATEPILEUP $DIRECTPHOTON
                fi
            fi
            if [ $DOGAMMA -eq 1 ]; then

                if [ -f $PI0DATARAWFILE ] && [ -f $PI0MCCORRFILE ]; then
                    if [ $NEWGAMMAMACROS == 0 ]; then
                        CorrectSignalGamma $PI0DATARAWFILE $PI0MCCORRFILE $CUTSELECTION $SUFFIX Pi0 kFALSE
                    else
                        CorrectSignalGammaV2 $PI0DATARAWFILE $PI0MCCORRFILE $CUTSELECTION $SUFFIX Pi0 kFALSE $DIRECTPHOTON
                    fi
                fi
                if [ -f $PI0MCRAWFILE ] && [ -f $PI0MCCORRFILE ]; then
                    if [ $NEWGAMMAMACROS == 0 ]; then
                        CorrectSignalGamma $PI0MCRAWFILE $PI0MCCORRFILE $CUTSELECTION $SUFFIX Pi0 kTRUE
                    else
                        CorrectSignalGammaV2 $PI0MCRAWFILE $PI0MCCORRFILE $CUTSELECTION $SUFFIX Pi0 kTRUE $DIRECTPHOTON
                    fi
                fi
                Pi0dataCorrFILE=`ls $CUTSELECTION/$ENERGY/Pi0_data_GammaConvV1Correction_*.root`
                GammaPi0dataCorrFILE=`ls $CUTSELECTION/$ENERGY/Gamma_Pi0_data_GammaConvV1Correction_*.root`
                Pi0MCCorrFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1Correction_*.root`
                GammaPi0MCCorrFILE=`ls $CUTSELECTION/$ENERGY/Gamma_Pi0_MC_GammaConvV1Correction_*.root`
                if [ $USECOCK -eq 1 ] && [ -f $Pi0dataCorrFILE ]; then
                    root -b -x -l -q TaskV1/PrepareCocktail.C\+\(\"$COCKROOTFILE\"\,\"$Pi0dataCorrFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,\"$COCKRAP\"\,\"\"\,$BINSPTPI0\,$MODE\)
                fi
                GammaCocktailFile=`ls $CUTSELECTION/$ENERGY/GammaCocktail_$COCKRAP*.root`
                if [ $USECOCK  ]; then
                    if [ $DOPI0TAG = 1  ]; then
                      Pi0dataUnCorrFILE=`ls $CUTSELECTION/$ENERGY/Pi0_data_GammaConvV1WithoutCorrection_*.root`
                      Pi0MCUnCorrFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaConvV1WithoutCorrection_*.root`

                      RunPi0Tagging $GammaPi0dataCorrFILE $Pi0dataUnCorrFILE $PI0MCCORRFILE $GammaCocktailFile $CUTSELECTION $SUFFIX kFALSE;
                      RunPi0Tagging $GammaPi0MCCorrFILE $Pi0MCUnCorrFILE $PI0MCCORRFILE $GammaCocktailFile $CUTSELECTION $SUFFIX kTRUE;
                    else
                      CreateGammaFinalResultsV3 $GammaPi0dataCorrFILE $Pi0dataCorrFILE $GammaCocktailFile $CUTSELECTION $SUFFIX Pi0 kFALSE;
                      CreateGammaFinalResultsV3 $GammaPi0MCCorrFILE $Pi0MCCorrFILE $GammaCocktailFile $CUTSELECTION $SUFFIX Pi0 kTRUE;
                    fi
                fi
            fi

            if [ $DOPI0INETABINS -eq 1 ]; then
                if [ -f $PI0ETADATARAWFILE ] && [ -f $PI0ETAMCCORRFILE ] ; then
                    CorrectSignal $PI0ETADATARAWFILE $PI0ETAMCCORRFILE $CUTSELECTION $SUFFIX Pi0EtaBinning kFALSE $ESTIMATEPILEUP $DIRECTPHOTON
                fi
                if [ -f $PI0ETAMCRAWFILE ] && [ -f $PI0ETAMCCORRFILE ] ; then
                    CorrectSignal $PI0ETAMCRAWFILE $PI0ETAMCCORRFILE $CUTSELECTION $SUFFIX Pi0EtaBinning kTRUE $ESTIMATEPILEUP $DIRECTPHOTON
                fi
            fi
            if [ $DOETA -eq 1 ]; then
                if [ -f $ETADATARAWFILE ] && [ -f $ETAMCCORRFILE ]; then
                    CorrectSignal $ETADATARAWFILE $ETAMCCORRFILE $CUTSELECTION $SUFFIX Eta kFALSE $ESTIMATEPILEUP $DIRECTPHOTON
                fi
                if [ -f $ETAMCRAWFILE ] && [ -f $ETAMCCORRFILE ]; then
                    CorrectSignal $ETAMCRAWFILE $ETAMCCORRFILE $CUTSELECTION $SUFFIX Eta kTRUE $ESTIMATEPILEUP $DIRECTPHOTON
                fi

            fi
            if [ $DOETAPRIME -eq 1 ]; then
                if [ -f $ETAPRIMEDATARAWFILE ] && [ -f $ETAPRIMEMCCORRFILE ]; then
                    CorrectSignal $ETAPRIMEDATARAWFILE $ETAPRIMEMCCORRFILE $CUTSELECTION $SUFFIX EtaPrime kFALSE $ESTIMATEPILEUP $DIRECTPHOTON
                fi
                if [ -f $ETAPRIMEMCRAWFILE ] && [ -f $ETAPRIMEMCCORRFILE ]; then
                    CorrectSignal $ETAPRIMEMCRAWFILE $ETAPRIMEMCCORRFILE $CUTSELECTION $SUFFIX EtaPrime kTRUE $ESTIMATEPILEUP $DIRECTPHOTON
                fi
            fi
            if [ "$SUFFIX" == "pdf" ]; then
                rm -f $CUTSELECTION/$ENERGY/$SUFFIX/CorrectSignalV2/CorrectSignal_all.pdf
                pdfunite $CUTSELECTION/$ENERGY/$SUFFIX/CorrectSignalV2/*.pdf $CUTSELECTION/$ENERGY/$SUFFIX/CorrectSignalV2/CorrectSignal_all.pdf
            fi

        fi
        NORMALCUTS=`expr $NORMALCUTS + 1`
    done

    if [ $NEWGAMMAMACROS == 1 ]; then
        if [ $DOGAMMA == 1 ]; then
            root -x -q -l -b TaskV1/GammaCutStudiesV3.C\+\(\"CutSelection.log\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,\"$SUFFIX\"\,$MODE\)
        fi
        DOGAMMA=0;
    fi

    if [ $DOPI0 -eq 1 ]; then
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Pi0\"\,\"kFALSE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,$DOGAMMA\,\"\"\,\"$PERIODNAME\"\,$MODE\)
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Pi0\"\,\"kTRUE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,$DOGAMMA\,\"\"\,\"$PERIODNAME\"\,$MODE\)
    fi
    if [ $DOPI0INETABINS -eq 1 ]; then
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Pi0EtaBinning\"\,\"kFALSE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"\"\,\"$PERIODNAME\"\,$MODE\)
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Pi0EtaBinning\"\,\"kTRUE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"\"\,\"$PERIODNAME\"\,$MODE\)
    fi
    if [ $DOETA -eq 1 ]; then
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Eta\"\,\"kFALSE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"\"\,\"$PERIODNAME\"\,$MODE\)
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Eta\"\,\"kTRUE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"\"\,\"$PERIODNAME\"\,$MODE\)
    fi
    if [ $DOETAPRIME -eq 1 ]; then
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"EtaPrime\"\,\"kFALSE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"\"\,\"$PERIODNAME\"\,$MODE\)
        root -x -q -l -b TaskV1/CutStudiesOverview.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"EtaPrime\"\,\"kTRUE\"\,\"$OPTMINBIASEFF\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"\"\,\"$PERIODNAME\"\,$MODE\)
    fi
else
    CORRECT=0
    while [ $CORRECT -eq 0 ]
    do
        echo "Please check that you really want to process all cuts, otherwise change the CutSelection.log. Make sure that the standard cut is the first in the file. Continue? Yes/No?";
        read answer
        if [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "y" ] || [ $answer = "yes" ]; then
            echo -e "--> Continuing ...\n";
            correct=1
        elif [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
            echo -e "--> Aborting ...\n";
            exit
        else
            echo "--> Command \"$answer\" not found. Please try again."
        fi
    done


echo -e "\n_____________________"
echo -e "STARTING ANALYSIS ...\n"
echo "  DOPI0=$DOPI0"
echo "  DOETA=$DOETA"
echo "  DOETAPRIME=$DOETAPRIME"
echo "  DOPI0INETABINS=$DOPI0INETABINS"
echo ""

    #Read the different cuts form the Cut selection log file
    CutSelections=`cat CutSelection.log`
    for CUTSELECTION in $CutSelections; do
        if [ -d $CUTSELECTION ]; then
            echo "CutSelection $CUTSELECTION directory already exists, all files will be overwritten ";
            mkdir $CUTSELECTION/$ENERGY
        else
            mkdir $CUTSELECTION
            mkdir $CUTSELECTION/$ENERGY
        fi

        if [ $ONLYCUTS -eq 0 ]; then
            if [ -d $CUTSELECTION/$ENERGY/$SUFFIX ]; then
                echo "Graphical Output $SUFFIX directory already exists, all files will be overwritten ";
            else
                mkdir $CUTSELECTION/$ENERGY/$SUFFIX
            fi

            if [ $DISABLETOYMC -eq 0 ] && [ $USECOCK -eq 0 ] && [ $ONLYCORRECTION -eq 0 ]; then
                rm ToyMCOutputs.txt
                root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\+\($NEVTSTOY,0,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\+\($NEVTSTOY,1,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
                root -b -x -l -q ToyModels/ModelSecondaryDecaysToPi0.C\+\($NEVTSTOY,2,\"$ENERGY\"\,$MINPTTOY\,$MAXPTTOY\,\"$EXTINPUTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,$MODE\)
            fi
            if [ $USECOCK -eq 1 ] && [ $ONLYCORRECTION -eq 0 ]; then
                root -b -x -l -q TaskV1/PrepareSecondaries.C\+\(\"Pi0\"\,\"$COCKROOTFILE\"\,\"$SUFFIX\"\,\"$CUTSELECTION\"\,\"$ENERGY\"\,\"$DIRECTPHOTON\"\,\"$COCKRAP\"\,\"\"\,$BINSPTPI0\,$MODE,kFALSE\)
            fi

            if [ $ONLYCORRECTION -eq 0 ]; then
                echo "CutSelection is $CUTSELECTION";
                if [ $DOPI0 -eq 1 ]; then
                    if [ -f $DATAROOTFILE ]; then
                        root -b -x -q -l TaskV1/ExtractSignalMergedMesonV2.C\+\(\"Pi0\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,kFALSE\,\"$ENERGY\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,$MODE\)
                    fi
                    PI0DATARAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_data_GammaMergedWithoutCorrection_*.root`
                    if [ $MCFILE -eq 1 ]; then
                        root -b -x -q -l TaskV1/ExtractSignalMergedMesonV2.C\+\(\"Pi0\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,kTRUE\,\"$ENERGY\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,$MODE\)
                        PI0MCRAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaMergedWithoutCorrection_*$CUTSELECTION*.root`
                        PI0MCCORRFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaMergedCorrectionHistos_*$CUTSELECTION*.root`
                        if [ $MERGINGMC -eq 1 ]; then
                            root -b -x -q -l TaskV1/ExtractSignalMergedMesonV2.C\+\(\"Pi0\"\,\"$MCROOTFILEGJ\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,kTRUE\,\"$ENERGY\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTPI0\,$MODE\,1\)
                            PI0MCCORRFILEJJG=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaMergedCorrectionHistosJJGammaTrigg_*.root`
                            root -b -x -q -l TaskV1/MergeCorrFactorsJJandJJGammaTrigMergedCluster.C\+\(\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,\"$PI0MCCORRFILE\"\,\"$CUTSELECTION/$ENERGY/Pi0_MC_GammaMergedCorrectionHistosJJ_$CUTSELECTION.root\"\,\"$PI0MCCORRFILEJJG\"\)

                        fi
                        echo $PI0DATARAWFILE
                        echo $PI0MCRAWFILE
                        echo $CUTSELECTION
                        root -b -x -q -l TaskV1/CompareShapeMergedClusterQuantities.C\+\(\"$PI0DATARAWFILE\"\,\"$PI0MCRAWFILE\"\,\"$CUTSELECTION\"\,\"Pi0\"\,\"$SUFFIX\"\,\"$ENERGY\"\,$BINSPTPI0\,$MODE\)
                    fi
                fi

                if [ $DOETA -eq 1 ]; then
                    if [ -f $DATAROOTFILE ]; then
                        root -b -x -q -l TaskV1/ExtractSignalMergedMesonV2.C\+\(\"Eta\"\,\"$DATAROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,kFALSE\,\"$ENERGY\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,$MODE\)
                    fi
                    ETADATARAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_data_GammaMergedWithoutCorrection_*.root`
                    if [ $MCFILE -eq 1 ]; then
                        root -b -x -q -l TaskV1/ExtractSignalMergedMesonV2.C\+\(\"Eta\"\,\"$MCROOTFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,kTRUE\,\"$ENERGY\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,$MODE\)
                        ETAMCRAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaMergedWithoutCorrection_*$CUTSELECTION*.root`
                        ETAMCCORRFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaMergedCorrectionHistos_*$CUTSELECTION*.root`
                        if [ $MERGINGMC -eq 1 ]; then
                            root -b -x -q -l TaskV1/ExtractSignalMergedMesonV2.C\+\(\"Eta\"\,\"$MCROOTFILEGJ\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,kTRUE\,\"$ENERGY\"\,\"\"\,\"$ADVMESONQA\"\,$BINSPTETA\,$MODE\,1\)
                        fi
                        root -b -x -q -l TaskV1/CompareShapeMergedClusterQuantities.C\+\(\"$ETADATARAWFILE\"\,\"$ETAMCRAWFILE\"\,\"$CUTSELECTION\"\,\"Eta\"\,\"$SUFFIX\"\,\"$ENERGY\"\,$BINSPTETA\,$MODE\)
                    fi
                fi
            fi

            if [ $DOPI0 -eq 1 ]; then
                PI0DATARAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_data_GammaMergedWithoutCorrection_$CUTSELECTION*.root`
                PI0MCRAWFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaMergedWithoutCorrection_$CUTSELECTION*.root`
                PI0MCCORRFILE=`ls $CUTSELECTION/$ENERGY/Pi0_MC_GammaMergedCorrectionHistos_$CUTSELECTION*.root`
            fi
            if [ $DOETA -eq 1 ]; then
                ETADATARAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_data_GammaMergedWithoutCorrection_$CUTSELECTION*.root`
                ETAMCRAWFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaMergedWithoutCorrection_$CUTSELECTION*.root`
                ETAMCCORRFILE=`ls $CUTSELECTION/$ENERGY/Eta_MC_GammaMergedCorrectionHistos_$CUTSELECTION*.root`
            fi

            if [ $DOPI0 -eq 1 ]; then
                if [ -f $PI0DATARAWFILE ] && [ -f $PI0MCCORRFILE ]; then
                      root -b -x -q -l TaskV1/CorrectSignalMergedV2.C\+\(\"$PI0DATARAWFILE\"\,\"$PI0MCCORRFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"Pi0\"\,kFALSE\,\"$ENERGY\"\,\"$PERIODNAME\"\,10\)
                fi
                if [ -f $PI0MCRAWFILE ] && [ -f $PI0MCCORRFILE ]; then
                    root -b -x -q -l TaskV1/CorrectSignalMergedV2.C\+\(\"$PI0MCRAWFILE\"\,\"$PI0MCCORRFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"Pi0\"\,kTRUE\,\"$ENERGY\"\,\"$PERIODNAME\"\,10\)
                fi
            fi
            if [ $DOETA -eq 1 ]; then
                if [ -f $ETADATARAWFILE ] && [ -f $ETAMCCORRFILE ]; then
                    root -b -x -q -l TaskV1/CorrectSignalMergedV2.C\+\(\"$ETADATARAWFILE\"\,\"$ETAMCCORRFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"Eta\"\,kFALSE\,\"$ENERGY\"\,\"$PERIODNAME\"\,10\)
                fi
                if [ -f $ETAMCRAWFILE ] && [ -f $ETAMCCORRFILE ]; then
                    root -b -x -q -l TaskV1/CorrectSignalMergedV2.C\+\(\"$ETAMCRAWFILE\"\,\"$ETAMCCORRFILE\"\,\"$CUTSELECTION\"\,\"$SUFFIX\"\,\"Eta\"\,kTRUE\,\"$ENERGY\"\,\"$PERIODNAME\"\,10\)
                fi
            fi
        fi
        NORMALCUTS=`expr $NORMALCUTS + 1`
    done

    if [ $DOPI0 -eq 1 ]; then
        root -l -b -x -q TaskV1/CutStudiesOverviewMerged.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Pi0\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,\"$PERIODNAME\"\,$MODE\)
        root -l -b -x -q TaskV1/CutStudiesOverviewMerged.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Pi0\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,\"$PERIODNAME\"\,$MODE\)
    fi
    if [ $DOETA -eq 1 ]; then
        root -l -b -x -q TaskV1/CutStudiesOverviewMerged.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Eta\"\,\"kFALSE\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,\"$PERIODNAME\"\,$MODE\)
        root -l -b -x -q TaskV1/CutStudiesOverviewMerged.C\+\(\"CutSelection.log\"\,\"$SUFFIX\"\,\"Eta\"\,\"kTRUE\"\,\"$ENERGY\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,\"$PERIODNAME\"\,$MODE\)
    fi
fi
