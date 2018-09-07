#! /bin/bash
#
#
#
# This script gests as input a directory where the GammaConvV1 root file is stored
# , it also needs the desired output directory where the produced root files are put.
# If nothing is given it will use ./ for the input directory and ./Output for the output
#
#Input 1: Root file to analyze Default: AnalyisResults
#Input 2: Input directory  Default:$PWD 
#Input 3: Output directory Default: $PWD/Results  (directory will be created if it does not exist)
#

PROGNAME=$0

function GiveBinningDalitz7TeV()
{
	echo "How many p_T bins do you want to use for the Pi0? 19(4GeV), 21(5GeV), 24(8GeV), 26(15GeV)";
	read answer
		if [ $answer = 19 ]; then
		    echo "19 Bins --> Max p_T = 4 GeV ...";
		    correctPi0=1
		    BinsPtPi0=19
		elif [ $answer = 21 ]; then
		    echo "21 Bins --> Max p_T = 5 GeV ...";
		    correctPi0=1
		    BinsPtPi0=21
		elif [ $answer = 24 ]; then
		    echo "24 Bins --> Max p_T = 8 GeV ...";
		    correctPi0=1
		    BinsPtPi0=24
		elif [ $answer = 26 ]; then
		    echo "26 Bins --> Max p_T = 15 GeV ...";
		    correctPi0=1
		    BinsPtPi0=26
		else
		    echo "Pi0 Binning was not set correctly. Please try again.";
		    correctPi0=0
		fi
        echo "How many p_t bins do you want to use for the eta meson? 7 (4.4GeV), 8 (6. GeV), 9 (10 GeV), 10 (14 GeV)"
 	read answer
		if [ $answer = 7 ]; then
		    echo "7 Bins --> Max p_T = 4.4 GeV ...";
		    correctEta=1
		    BinsPtEta=7
		elif [ $answer = 8 ]; then
		    echo "8 Bins --> Max p_T = 6 GeV ...";
		    correctEta=1
		    BinsPtEta=8
		elif [ $answer = 9 ]; then
		    echo "9 Bins --> Max p_T = 10 GeV ...";
		    correctEta=1
		    BinsPtEta=9
		elif [ $answer = 10 ]; then
		    echo "10 Bins --> Max p_T = 14 GeV ...";
		    correctEta=1
		    BinsPtEta=10
		else
		    echo "Eta Binning was not set correctly. Please try again.";
		    correctEta=0
		fi

}

function GiveBinningDalitz5TeV2017()
{
	echo "How many p_T bins do you want to use for the Pi0? 19(5GeV), 20(6GeV), 21(8GeV), 22(10GeV)";
	read answer
		if [ $answer = 19 ]; then
		    echo "19 Bins --> Max p_T = 5 GeV ...";
		    correctPi0=1
		    BinsPtPi0=19
		elif [ $answer = 20 ]; then
		    echo "20 Bins --> Max p_T = 6 GeV ...";
		    correctPi0=1
		    BinsPtPi0=20
		elif [ $answer = 21 ]; then
		    echo "21 Bins --> Max p_T = 8 GeV ...";
		    correctPi0=1
		    BinsPtPi0=21
		elif [ $answer = 22 ]; then
		    echo "22 Bins --> Max p_T = 10 GeV ...";
		    correctPi0=1
		    BinsPtPi0=22
		else
		    echo "Pi0 Binning was not set correctly. Please try again.";
		    correctPi0=0
		fi
 	echo "How many p_t bins do you want to use for the eta meson? 7 (4.4GeV), 8 (6. GeV), 9 (10 GeV)"
 	read answer
		if [ $answer = 7 ]; then
		    echo "7 Bins --> Max p_T = 4.4 GeV ...";
		    correctEta=1
		    BinsPtEta=7
		elif [ $answer = 8 ]; then
		    echo "8 Bins --> Max p_T = 6 GeV ...";
		    correctEta=1
		    BinsPtEta=8
		elif [ $answer = 9 ]; then
		    echo "9 Bins --> Max p_T = 10 GeV ...";
		    correctEta=1
		    BinsPtEta=9
		else
		    echo "Eta Binning was not set correctly. Please try again.";
		    correctEta=0
		fi

}

function GiveBinningDirectPhotonHI()
{
    echo "How many p_T bins do you want to use for the Pi0? 17(8GeV), 18(11GeV) 19 (14GeV)";
    read answer
    if [ $answer = 17 ]; then
	echo "17 Bins --> Max p_T = 8 GeV ...";
	correctPi0=1
	BinsPtPi0=17
    elif [ $answer = 18 ]; then
	echo "18 Bins --> Max p_T = 11 GeV ...";
	correctPi0=1
	BinsPtPi0=18
    elif [ $answer = 19 ]; then
	echo "19 Bins --> Max p_T = 20 GeV ...";
	correctPi0=1
	BinsPtPi0=18
    else
   echo "Pi0 Binning was not set correctly. Please try again.";
   correctPi0=0
    fi
    echo "How many p_T bins do you want to use for the Eta? 2(2GeV), 3(4GeV), 4(7GeV)";
    read answer
    if [ $answer = 2 ]; then
   echo "2 Bins --> Max p_T = 2 GeV ...";
   correctEta=1
   BinsPtEta=2
    elif [ $answer = 3 ]; then
   echo "3 Bins --> Max p_T = 4 GeV ...";
   correctEta=1
   BinsPtEta=3
    elif [ $answer = 4 ]; then
   echo "4 Bins --> Max p_T = 7 GeV ...";
   correctEta=1
   BinsPtEta=4
    else
   echo "Eta Binning was not set correctly. Please try again.";
   correctEta=0
    fi
}

function GiveBinningDirectPhoton2760TeV()
{
    echo "How many p_T bins do you want to use for the Pi0? 13(5.5GeV), 14(7GeV) 15 (10GeV)";
    read answer
    if [ $answer = 13 ]; then
	echo "13 Bins --> Max p_T = 5.5 GeV ...";
	correctPi0=1
	BinsPtPi0=13
    elif [ $answer = 14 ]; then
	echo "14 Bins --> Max p_T = 7 GeV ...";
	correctPi0=1
	BinsPtPi0=14
    elif [ $answer = 15 ]; then
	echo "15 Bins --> Max p_T = 10 GeV ...";
	correctPi0=1
	BinsPtPi0=15
    else
	echo "Pi0 Binning was not set correctly. Please try again.";
	correctPi0=0
    fi
    echo "How many p_T bins do you want to use for the Eta? 5(2.5GeV), 6(4GeV), 7(6GeV)";
    read answer
    if [ $answer = 5 ]; then
	echo "5 Bins --> Max p_T = 2.5 GeV ...";
	correctEta=1
	BinsPtEta=5
    elif [ $answer = 6 ]; then
	echo "6 Bins --> Max p_T = 4 GeV ...";
	correctEta=1
	BinsPtEta=6
    elif [ $answer = 7 ]; then
	echo "7 Bins --> Max p_T = 6 GeV ...";
	correctEta=1
	BinsPtEta=7
    else
	echo "Eta Binning was not set correctly. Please try again.";
	correctEta=0
    fi
}

function GiveBinningDirectPhoton7TeV()
{
    echo "How many p_T bins do you want to use for the Pi0? 19(10GeV), 20(12GeV), 21(16GeV), 22(20GeV) 23 (25GeV)";
    read answer
    if [ $answer = 19 ]; then
	echo "19 Bins --> Max p_T = 8 GeV ...";
	correctPi0=1
	BinsPtPi0=19
    elif [ $answer = 20 ]; then
	echo "20 Bins --> Max p_T = 12 GeV ...";
	correctPi0=1
	BinsPtPi0=20
    elif [ $answer = 21 ]; then
	echo "21 Bins --> Max p_T = 16 GeV ...";
	correctPi0=1
	BinsPtPi0=21
    elif [ $answer = 22 ]; then
	echo "22 Bins --> Max p_T = 20 GeV ...";
	correctPi0=1
	BinsPtPi0=22
    elif [ $answer = 23 ]; then
	echo "23 Bins --> Max p_T = 25 GeV ...";
	correctPi0=1
	BinsPtPi0=23
    else
	echo "Pi0 Binning was not set correctly. Please try again.";
	correctPi0=0
    fi
    echo "How many p_T bins do you want to use for the Eta? 10(4GeV), 11(6GeV), 12(8GeV)";
    read answer
    if [ $answer = 10 ]; then
	echo "10 Bins --> Max p_T = 4 GeV ...";
	correctEta=1
	BinsPtEta=10
    elif [ $answer = 11 ]; then
	echo "11 Bins --> Max p_T = 6 GeV ...";
	correctEta=1
	BinsPtEta=11
    elif [ $answer = 12 ]; then
	echo "12 Bins --> Max p_T = 8 GeV ...";
	correctEta=1
	BinsPtEta=12
    else
	echo "Eta Binning was not set correctly. Please try again.";
	correctEta=0
    fi
}

function GiveBinningNormal()
{
    echo "How many p_T bins do you want to use for the Pi0? 22(7GeV), 27(8GeV), 28(10GeV), 29(12GeV) 30 (16GeV) 32 (25GeV)";
    
    read answer
    if [ $answer = 26 ]; then
	echo "26 Bins --> Max p_T = 7 GeV ...";
	correctPi0=1
	BinsPtPi0=26
    elif [ $answer = 27 ]; then
	echo "27 Bins --> Max p_T = 8 GeV ...";
	correctPi0=1
	BinsPtPi0=27
    elif [ $answer = 28 ]; then
	echo "28 Bins --> Max p_T = 10 GeV ...";
	correctPi0=1
	BinsPtPi0=28
    elif [ $answer = 29 ]; then
	echo "29 Bins --> Max p_T = 12 GeV ...";
	correctPi0=1
	BinsPtPi0=29
    elif [ $answer = 30 ]; then
	echo "30 Bins --> Max p_T = 16 GeV ...";
	correctPi0=1
	BinsPtPi0=30
    elif [ $answer = 32 ]; then
	echo "32 Bins --> Max p_T = 25 GeV ...";
	correctPi0=1
	BinsPtPi0=32
    else
	echo "Pi0 Binning was not set correctly. Please try again.";
	correctPi0=0
    fi
    echo "How many p_T bins do you want to use for the Eta? 10(4GeV), 11(6GeV), 12(8GeV)";
    read answer
    if [ $answer = 10 ]; then
	echo "10 Bins --> Max p_T = 4 GeV ...";
	correctEta=1
	BinsPtEta=10
    elif [ $answer = 11 ]; then
	echo "11 Bins --> Max p_T = 6 GeV ...";
	correctEta=1
	BinsPtEta=11
    elif [ $answer = 12 ]; then
	echo "12 Bins --> Max p_T = 8 GeV ...";
	correctEta=1
	BinsPtEta=12
    else
	echo "Eta Binning was not set correctly. Please try again.";
	correctEta=0
    fi

}
function GiveBinning900GeV()
{
    echo "How many p_T bins do you want to use? 10(3GeV), 11(4GeV)";
    read answer
    if [ $answer = 10 ]; then
	echo "10 Bins --> Max p_T = 3 GeV ...";
	correctPi0=1
	BinsPtPi0=10
    elif [ $answer = 11 ]; then
	echo "11 Bins --> Max p_T = 4 GeV ...";
	correctPi0=1
	BinsPtPi0=11
    else
	echo "Pi0 Binning was not set correctly. Please try again.";
    fi
    echo "How many p_t bins do you want to use for the eta meson? 2 (1.8 GeV), 3 (4 GeV)"
    read answer
    if [ $answer = 2 ]; then
	echo "2 Bins --> Max p_T = 1.8 GeV ...";
	correctEta=1
	BinsPtEta=2
    elif [ $answer = 3 ]; then
	echo "3 Bins --> Max p_T = 4 GeV ...";
	correctEta=1
	BinsPtEta=3
    else
	echo "Eta Binning was not set correctly. Please try again.";
	correctEta=0
    fi

}

function GiveBinning2760GeV()
{
    echo "How many p_T bins do you want to use?  6 (5GeV), 7 (10GeV)";
    read answer
    if [ $answer = 6 ]; then
	echo "6 Bins --> Max p_T = 5 GeV ...";
	correctPi0=1
	BinsPtPi0=6
    elif [ $answer = 7 ]; then
	echo "7 Bins --> Max p_T = 10 GeV ...";
	correctPi0=1
	BinsPtPi0=7
    else
	echo "Pi0 Binning was not set correctly. Please try again.";
    fi
    echo "How many p_t bins do you want to use for the eta meson? 6 (4. GeV), 7 (6 GeV)"
    read answer
    if [ $answer = 6 ]; then
	echo "6 Bins --> Max p_T = 4. GeV ...";
	correctEta=1
	BinsPtEta=6
    elif [ $answer = 7 ]; then
	echo "7 Bins --> Max p_T = 6 GeV ...";
	correctEta=1
	BinsPtEta=7
    else
	echo "Eta Binning was not set correctly. Please try again.";
	correctEta=0
    fi
}

function GiveBinningHI()
{
    echo "How many p_T bins do you want to use for Pi0? 7(6GeV), 8(10GeV)";
    read answer
    if [ $answer = 7 ]; then
	echo "7 Bins --> Max p_T = 6 GeV ...";
	correctPi0=1
	BinsPtPi0=7
    elif [ $answer = 8 ]; then
	echo "8 Bins --> Max p_T = 10 GeV ...";
	correctPi0=1
	BinsPtPi0=8
    else
	echo "Pi0 Binning was not set correctly. Please try again.";
	correctPi0=0
    fi
    if [ $DoEtaHI -eq 1 ]; then
	echo "How many p_T bins do you want to use for Eta? 3(4GeV), 4(7GeV)";
	read answer
	if [ $answer = 3 ]; then
	    echo "3 Bins --> Max p_T = 4 GeV ...";
	    correctEta=1
	    BinsPtEta=3
	elif [ $answer = 4 ]; then
	    echo "4 Bins --> Max p_T = 7 GeV ...";
	    correctEta=1
	    BinsPtEta=4
	else
	    echo "Pi0 Binning was not set correctly. Please try again.";
	    correctEta=0
	fi
    fi
}

function GiveBinningDalitz5023GeV()
{
	echo "How many p_T bins do you want to use for the Pi0? 19(6GeV), 21(10GeV), 22(15GeV)";
	read answer
		if [ $answer = 19 ]; then
		    echo "19 Bins --> Max p_T = 6 GeV ...";
		    correctPi0=1
		    BinsPtPi0=19
		elif [ $answer = 21 ]; then
		    echo "21 Bins --> Max p_T = 10 GeV ...";
		    correctPi0=1
		    BinsPtPi0=21
		elif [ $answer = 22 ]; then
		    echo "22 Bins --> Max p_T = 15 GeV ...";
		    correctPi0=1
		    BinsPtPi0=22
		else
		    echo "Pi0 Binning was not set correctly. Please try again.";
		    correctPi0=0
		fi
 	echo "How many p_t bins do you want to use for the eta meson? 7 (4.4GeV), 8 (6. GeV), 9 (10 GeV)"
 	read answer
		if [ $answer = 7 ]; then
		    echo "7 Bins --> Max p_T = 4.4 GeV ...";
		    correctEta=1
		    BinsPtEta=7
		elif [ $answer = 8 ]; then
		    echo "8 Bins --> Max p_T = 6 GeV ...";
		    correctEta=1
		    BinsPtEta=8
		elif [ $answer = 9 ]; then
		    echo "9 Bins --> Max p_T = 10 GeV ...";
		    correctEta=1
		    BinsPtEta=9
		else
		    echo "Eta Binning was not set correctly. Please try again.";
		    correctEta=0
		fi

}

function GiveBinningpPb()
{
    echo "How many p_T bins do you want to use for Pi0? 26(7GeV), 27(8GeV), 28(10GeV), 29(15GeV)";
    read answer
     if [ $answer = 26 ]; then
      echo "26 Bins --> Max p_T = 7 GeV ...";
      correctPi0=1
      BinsPtPi0=26
   elif [ $answer = 27 ]; then
   	echo "27 Bins --> Max p_T = 8 GeV ...";
   	correctPi0=1
  	BinsPtPi0=27
    elif [ $answer = 28 ]; then
  	echo "28 Bins --> Max p_T = 10 GeV ...";
   	correctPi0=1
   	BinsPtPi0=28
    elif [ $answer = 29 ]; then
  	echo "29 Bins --> Max p_T = 15 GeV ...";
   	correctPi0=1
   	BinsPtPi0=29
    else
   	echo "Pi0 Binning was not set correctly. Please try again.";
   	correctPi0=0
    fi
    if [ $DoEtaHI -eq 1 ]; then 
    	echo "How many p_t bins do you want to use for the eta meson? 6 (4. GeV), 7 (6 GeV), 8 (8 GeV)";
    	read answer
    	if [ $answer = 6 ]; then
        	echo "6 Bins --> Max p_T = 4. GeV ...";
       		correctEta=1
       		BinsPtEta=6
    	elif [ $answer = 7 ]; then
         echo "7 Bins --> Max p_T = 6 GeV ...";
         correctEta=1
         BinsPtEta=7
      elif [ $answer = 8 ]; then
         echo "8 Bins --> Max p_T = 6 GeV ...";
         correctEta=1
         BinsPtEta=8
      else
        	echo "Eta Binning was not set correctly. Please try again.";
        	correctEta=0
    	fi
    fi
}


function ExtractSignal()
{
    root -x -q -l -b  TaskV1/ExtractSignalV2.C\+\+\($1\)  
    if [ $2 = "Pi0" ]; then
      if [ $GammaOn -eq 1 ]; then
         root -x -l -b -q TaskV1/ExtractGammaSignal.C\($1\)
      fi
    fi
}

function CorrectSignal()
{
	
	 #root -x -l -b -q TaskV1/CorrectSignalDalitzV2.C\+\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$energy\"\,\"$7\"\,\"$8\"\,\"kTRUE\"\,$9\,$mode\)
	 #root -x -l -b -q TaskV1/CorrectSignalV2.C\+\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$energy\"\,\"\"\,\"$ESTIMATEPILEUP\"\,kFALSE\,$7\,$mode\)
	  root -x -l -b -q TaskV1/CorrectSignalDalitzV2.C\+\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$energy\"\,\"\"\,\"$ESTIMATEPILEUP\"\,kTRUE\,\"$7\"\,$mode\)
	 	  
}


function CreateFinalResults() 
{
	root -x -l -b -q TaskV1/ProduceFinalResults.C\+\+\($1\,\"Dalitz\"\,\"kTRUE\"\,$AddPileUpCorr\,$mode\)
}

function ProduceFinalResultspPb() 
{
    root -x -l -b -q TaskV1/ProduceFinalResultspPb.C\+\+\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$7\"\,kFALSE\,\"kTRUE\"\,$mode\)
    
}

function CreateGammaFinalResults() 
{
    if [ $5 = "Pi0" ]; then
	root -x -l -b -q TaskV1/CalculateGammaToPi0.C++\(\"$1\"\,\"$2\"\,\"$3\"\,\"$4\"\,\"$5\"\,\"$6\"\,\"$energy\"\,\"$7\"\)
    fi
}


function Usage() 
{
    echo -e "
This Script is provided by the Gamma Conversion Group of ALICE 
the main developers are 
	Friederike Bock \t friederike.bock@cern.ch
	Kathrin Koch \t\t kkoch@physi.uni-heidelberg.de
	Ana Marin \t\t marin@physi.uni-heidelberg.de  
	Kenneth Aamoth \t\t work.kenneth.aamodt@gmail.com
			
If there are any complications with running this script do not hesitate to connect them.
			

How to use this script? 
$PROGNAME -h  \t\t\t\t\t usage will be displayed 
$PROGNAME --help \t\t\t\t usage will be displayed 
$PROGNAME -c data.root suffix \t\t\t Will only correct and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -d data.root suffix \t\t\t Will only fullfil cutstudies and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -m DirectoryOfMergedOutputs suffix \t Automatic merging for LHC10bc and LHC10de for efficiencies will be done \n \t\t\t\t\t\t\t\t according to the fraction in Data, be careful files have to be arranged in a certain structure \n
	\t\t\t\t\t\t\t\t DataRootFile=DirectoryOfMergedOutputs/mergedALL/GammaConvV1Data.root
	\t\t\t\t\t\t\t\t MCRootFile=DirectoryOfMergedOutputs/mergedALL/GammaConvV1MC.root
	\t\t\t\t\t\t\t\t MCRootFileBC=DirectoryOfMergedOutputs/mergedBC/GammaConvV1MC.root
	\t\t\t\t\t\t\t\t MCRootFileD=DirectoryOfMergedOutputs/mergedDE/GammaConvV1MC.root\n 
$PROGNAME -r data.root suffix \t\t\t Will only execute the production of the final results and \n \t\t\t\t\t\t\t\t produce the graphical output in the specified suffix-format.\n
$PROGNAME data.root MC.root suffix \t \t Will execute Gamma Conversion Analysis for data.root file and MC.root file \n \t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n

This script needs as a basis the output from the GammaConversion-Software which is provided with the Aliroot software, both the 'data.root' and the 'MC.root' have to be output files of this. The script then will check which cutnumbers are in these files and will ask you which exactly you want to analyse, furthermore you have to set a standard cut which has to be always the first in the CutLogFile. Because for this cut the systematic errors will be calculated. Not only one analysis will be fullfiled, you can additionally choose to do a Alpha studie or Chi2 of the meson as well which will give you the oportunity to calculate two indepent error due to cut variation. Additionally the CorrectSignal.C will correct the spectrum for all possible contributions and afterwards calculate the systematic error due to yield extraction. All these error will then enter the final results where you will have plots with only statistical errors as well as systematic + static errors. Several data output files are created, and 
stored in the corresponding cut-directory or the working directory for the cutsstudies file. 
		
	"
    exit 
}

if [ "$#" == "0" ]; then
    echo -e "
$PROGNAME -h  \t\t\t\t\t usage will be displayed 
$PROGNAME --help \t\t\t\t usage will be displayed 
$PROGNAME -c data.root suffix \t\t\t Will only correct and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -d data.root suffix \t\t\t Will only fullfil cutstudies and produce the final results for already existing RAWdata \n \t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
$PROGNAME -m DirectoryOfMergedOutputs suffix \t Automatic merging for LHC10bc and LHC10de for efficiencies will be done \n \t\t\t\t\t\t\t\t according to the fraction in Data, be careful files have to be arranged in a certain structure \n
	\t\t\t\t\t\t\t\t DataRootFile=DirectoryOfMergedOutputs/mergedALL/GammaConvV1Data.root
	\t\t\t\t\t\t\t\t MCRootFile=DirectoryOfMergedOutputs/mergedALL/GammaConvV1MC.root
	\t\t\t\t\t\t\t\t MCRootFileBC=DirectoryOfMergedOutputs/mergedBC/GammaConvV1MC.root
	\t\t\t\t\t\t\t\t MCRootFileD=DirectoryOfMergedOutputs/mergedDE/GammaConvV1MC.root\n 
$PROGNAME -r data.root suffix \t\t\t Will only execute the production of the final results and \n \t\t\t\t\t\t\t\t produce the graphical output in the specified suffix-format.\n
$PROGNAME data.root MC.root suffix \t \t Will execute Gamma Conversion Analysis for data.root file and MC.root file \n \t\t\t\t\t\t\t\t and produce the graphical output in the specified suffix-format.\n
"
    exit
fi

ONLYCORRECTION=0
ONLYRESULTS=0
HIRUN=0
ONLYCUTS=0
Suffix=eps
MERGINGMC=0
PARTLY=0
DATAFILE=1
MCFILE=1
addedSig=0
pp2760=0
GammaOn=1
AddPileUpCorr=kFALSE
NAMECUTSTUDIES="none"
DoEtaHI=1
#ComputeGGCont=0
#LHC14b2MC=0
MCSample="Sample"
#setBR=0
mode=0

if [[ "$1" == *-*gammaOff* ]]; then
    GammaOn=0
    echo "gamma calculation switched off"
fi

if [[ "$1" == *-*aPUC* ]]; then
    AddPileUpCorr=kTRUE
    echo "additionally pileup correction on"
fi

if [[ "$1" == *-h* ]] ; then
    Usage
elif [[ "$1" == *-mAddSig2760GeVA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedMinBias/GammaConvV1Data.root
    MCRootFile=$DIRECTORY/mergedMinBias/GammaConvV1MC.root
    MCRootFileAddSig=$DIRECTORY/mergedAddSignal/GammaConvV1MCAddSigA.root
	pp2760=1
    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
elif [[ "$1" == *-mAddSig2760GeVB* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedMinBias/GammaConvV1Data.root
    MCRootFile=$DIRECTORY/mergedMinBias/GammaConvV1MC.root
    MCRootFileAddSig=$DIRECTORY/mergedAddSignal/GammaConvV1MCAddSigB.root
	pp2760=1
    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
elif [[ "$1" == *-mAddSigWithSDD2760GeVA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedMinBias/GammaConvV1Data_with_SDD.root
    MCRootFile=$DIRECTORY/mergedMinBias/GammaConvV1MC_with_SDD.root
    MCRootFileAddSig=$DIRECTORY/mergedAddSignal/GammaConvV1MCAddSigA_with_SDD.root
	 pp2760=1
    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
elif [[ "$1" == *-mAddSigWithSDD2760GeVB* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedMinBias/GammaConvV1Data_with_SDD.root
    MCRootFile=$DIRECTORY/mergedMinBias/GammaConvV1MC_with_SDD.root
    MCRootFileAddSig=$DIRECTORY/mergedAddSignal/GammaConvV1MCAddSigB_with_SDD.root
    pp2760=1
    if [ -n $DataRootFile ]; then
   echo "The data file specified is $DataRootFile"
    else 
   echo "No data file specified, analysis can not be fullfiled."
   exit
    fi
    if [ -n $MCRootFile ]; then
   echo "The MC file specified is $MCRootFile"
    else 
   echo "No MC file specified, analysis will only made paritally, please be careful with the results."
   PARTLY=1 
   MCFILE=0
    fi
elif [[ "$1" == *-mAddSigpPbHIJINGA* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedMinBias/GammaConvV1Data_A.root
    MCRootFile=$DIRECTORY/mergedMinBias/GammaConvV1_HIJING_MC_A.root
    MCRootFileAddSig=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_A.root
    addedSig=1
    if [ -f $DataRootFile ]; then
   echo "The data file specified is $DataRootFile"
    else 
   echo "No data file specified, analysis can not be fullfiled."
#  exit
    fi
    if [ -f $MCRootFile ]; then
   echo "The MC file specified is $MCRootFile"
    else 
   echo "No MC file specified, analysis will only made paritally, please be careful with the results."
   PARTLY=1 
   MCFILE=0
    fi                    
elif [[ "$1" == *-mAddSigpPbHIJINGB* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedMinBias/GammaConvV1Data_B.root
    MCRootFile=$DIRECTORY/mergedMinBias/GammaConvV1_HIJING_MC_B.root
    MCRootFileAddSig=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_B.root
   addedSig=1
    if [ -f $DataRootFile ]; then
   echo "The data file specified is $DataRootFile"
    else 
   echo "No data file specified, analysis can not be fullfiled."
#  exit
    fi
    if [ -f $MCRootFile ]; then
   echo "The MC file specified is $MCRootFile"
    else 
   echo "No MC file specified, analysis will only made paritally, please be careful with the results."
   PARTLY=1 
   MCFILE=0
    fi                        
elif [[ "$1" == *-mAddSigpPbHIJINGC* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedMinBias/GammaConvV1Data_C.root
    MCRootFile=$DIRECTORY/mergedMinBias/GammaConvV1_HIJING_MC_C.root
    MCRootFileAddSig=$DIRECTORY/mergedAddSignal/GammaConvV1_HIJING_MC_C.root
   addedSig=1
    if [ -f $DataRootFile ]; then
   echo "The data file specified is $DataRootFile"
    else 
   echo "No data file specified, analysis can not be fullfiled."
#  exit
    fi
    if [ -f $MCRootFile ]; then
   echo "The MC file specified is $MCRootFile"
    else 
   echo "No MC file specified, analysis will only made paritally, please be careful with the results."
   PARTLY=1 
   MCFILE=0
    fi                            
elif [[ "$1" == *-mAddSig* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedAll/GammaConvV1Data.root
    MCRootFile=$DIRECTORY/mergedAll/GammaConvV1MC.root
    MCRootFileBC=$DIRECTORY/WithoutAddedSignals/GammaConvV1MC.root
    MCRootFileD=$DIRECTORY/WithAddedSignals/GammaConvV1MC.root

    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
elif [[ "$1" == *-m* ]] ; then
    MERGINGMC=1
    DIRECTORY=$2
    Suffix=$3
    DataRootFile=$DIRECTORY/mergedALL/GammaConvV1Data.root
    MCRootFile=$DIRECTORY/mergedALL/GammaConvV1MC.root
    MCRootFileBC=$DIRECTORY/mergedBC/GammaConvV1MC.root
    MCRootFileD=$DIRECTORY/mergedDE/GammaConvV1MC.root
    DataRootFileBC=$DIRECTORY/mergedBC/GammaConvV1Data.root
    DataRootFileD=$DIRECTORY/mergedDE/GammaConvV1Data.root

    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
elif [[ "$1" == *-c* ]] ; then
    ONLYCORRECTION=1
    DataRootFile=$2
    Suffix=$3;
elif [[ "$1" == *-d* ]] ; then
    ONLYCORRECTION=1
    ONLYCUTS=1
    DataRootFile=$2
    Suffix=$3;
elif [[ "$1" == *-e* ]] ; then
    ONLYCORRECTION=1
    ONLYCUTS=1
    Suffix=$2;
elif [[ "$1" == *-r* ]] ; then
    ONLYCORRECTION=1
    ONLYRESULTS=1
    DataRootFile=$2
    Suffix=$3;
elif [[ "$1" == *-HI* ]] ; then
    HIRUN=1
    DataRootFile=$2
    MCRootFile=$3
    Suffix=$4;
elif [[ "$1" != -* ]] ; then
    DataRootFile=$1
    MCRootFile=$2
    Suffix=$3;
    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
elif [[ "$1" == *-LHC14b2* ]]; then
     #ComputeGGCont=1
     #setBR=1
     #MCLHC14b2=1
     #MCLHC14b2
     MCSample="LHC14b2"
     
     echo "MC (LHC14b2)The gamma gamma contamination will be scaled"
     DataRootFile=$2
     MCRootFile=$3
     Suffix=$4
    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
else 
    DataRootFile=$2
    MCRootFile=$3
    Suffix=$4;
    if [ -n $DataRootFile ]; then
	echo "The data file specified is $DataRootFile"
    else 
	echo "No data file specified, analysis can not be fullfiled."
	exit
    fi
    if [ -n $MCRootFile ]; then
	echo "The MC file specified is $MCRootFile"
    else 
	echo "No MC file specified, analysis will only made paritally, please be careful with the results."
	PARTLY=1	
	MCFILE=0
    fi
fi

Noeta=0
NoPi0Eta=1;
NORMALCUTS=0
# ALPHACUTS=0

if [ $ONLYCUTS -eq 1 ]; then 
    correct=0
    while [ $correct -eq 0 ]
    do
	echo "Which name should your CutStudies have? None/*name*"
	read answer
	if [ $answer = "None" ]; then
	    correct=1
	else
	    NAMECUTSTUDIES=$answer
	    echo "You have selcted the name: " $NAMECUTSTUDIES;
	    correct=1
	fi
    done
fi


correct=0
while [ $correct -eq 0 ]
do
	echo "Which mode are you running? 9 (Dalitz old files), 1 (New PCM-Dalitz), 6 (New Dalitz-EMCAL), 7 (Dalitz-PHOS)"
	read answer
	if [ $answer = "9" ]; then
		echo "You are analysing Dalitz old output";
		mode=9
		correct=1
	elif [ $answer = "1" ]; then
		echo "You are analysing Dalitz new output.";
		mode=1
		correct=1
	elif [ $answer = "6" ]; then
		echo "You are analysing Dalitz with gamma reconstructed in EMCAL";
		mode=6
		correct=1
	elif [ $answer = "7" ]; then
		echo "You are analysing Dalitz with gamma reconstructed in PHOS";
		mode=7
		correct=1
	else
		echo "Command not found. Please try again.";
	fi
	
done


correct=0
while [ $correct -eq 0 ]
do
    echo "Do you want to take an already exitsting CutSelection.log-file. Yes/No"
    read answer
    if [ $answer = "Yes" ]; then
	echo "Chosen already existing logfile ...";
	cat CutSelection.log
	correct=1
    elif [ $answer = "No" ]; then
	if [ $pp2760 -eq 0 ]; then
	   root -b -q -x -l TaskV1/MakeCutLog.C\(\"$DataRootFile\"\,\"CutSelection.log\"\,$mode\) 
	   correct=1
	else 
 	  root -b -q -x -l TaskV1/MakeCutLog.C\(\"$MCRootFileAddSig\"\,\"CutSelection.log\"\,$mode\)
	  correct=1
	fi
    else
	echo "Command not found. Please try again.";
    fi
done

echo "On the top you see all cuts being available, please tell me the standard cut for Meson analysis.";
read standardCutMeson

correct=0
while [ $correct -eq 0 ]
do
    echo "The standard cut was set to $standardCutMeson. Is this correct? Yes/No?";
    read answer
    if [ $answer = "Yes" ]; then
	echo "Continuing ...";
	correct=1
    elif [ $answer = "No" ]; then
	echo "Aborting ...";
	exit
    else
	echo "Command not found. Please try again.";
    fi
done



correct=0
while [ $correct -eq 0 ]
do
    echo "Which collision system do you want to process? 7TeV (pp@7TeV), 5TeV2017 (pp@5TeV2017), 900GeV (pp@900GeV), 2.76TeV (pp@2.76TeV), PbPb_2.76TeV (PbPb@2.76TeV), pPb_5.023TeV (PbPb@2.76TeV)"
    read energy
    if [ $energy = "7TeV" ]; then
      
         Conference="No"
         Con=0
         if [ $ONLYCORRECTION -eq 0 ]; then
	 GiveBinningDalitz7TeV 
         correctPi0=1
         correctEta=1
         fi
         #if [ $correctPi0 -eq 0 ]; then
         #correct=0
         #elif [ $correctEta -eq 0 ]; then
         #correct=0
         #else 
         correct=1
         #fi
      #else
      #   echo "Command not found. Please try again.";
      #fi
    elif [ $energy = "5TeV2017" ]; then
      
         Conference="No"
         Con=0
         if [ $ONLYCORRECTION -eq 0 ]; then
         GiveBinningDalitz5TeV2017
         correctPi0=1
         correctEta=1
         fi
         #if [ $correctPi0 -eq 0 ]; then
         #correct=0
         #elif [ $correctEta -eq 0 ]; then
         #correct=0
         #else 
         correct=1
         #fi
      #else
      #   echo "Command not found. Please try again.";
      #fi

    elif [ $energy = "900GeV" ]; then
      
         directphoton="No"
         Con=0
         if [ $ONLYCORRECTION -eq 0 ]; then
         GiveBinning900GeV
         correctPi0=1
         correctEta=1
         fi
         #if [ $correctPi0 -eq 0 ]; then
         #correct=0
         #elif [ $correctEta -eq 0 ]; then
         #correct=0
         #else 
         correct=1
         #fi
      #else
      #   echo "Command not found. Please try again.";
      #fi
    elif [ $energy = "2.76TeV" ]; then
      
               
         directphoton="No"
         Con=0
         if [ $ONLYCORRECTION -eq 0 ]; then
         GiveBinning2760GeV
         correctPi0=1
         correctEta=1
         fi
         #if [ $correctPi0 -eq 0 ]; then
         #correct=0
         #elif [ $correctEta -eq 0 ]; then
         #correct=0
         #else 
         correct=1
         #fi
         #else
	 #  echo "Command not found. Please try again.";
	#fi
    elif [ $energy = "PbPb_2.76TeV" ]; then
      
         directphoton="No"
         Con=0
         if [ $ONLYCORRECTION -eq 0 ]; then
            GiveBinningHI
            correctPi0=1
            correctEta=1
         fi
         #if [ $correctPi0 -eq 0 ]; then
         #   correct=0
         #elif [ $correctEta -eq 0 ]; then
         #   correct=0
         #else 
            correct=1
         #fi
     
    elif [ $energy = "pPb_5.023TeV" ]; then
      
         directphoton="No"
         Con=0
         if [ $ONLYCORRECTION -eq 0 ]; then
            GiveBinningDalitz5023GeV
            correctPi0=1
            correctEta=1
         fi
         #if [ $correctPi0 -eq 0 ]; then
         #   correct=0
         #elif [ $correctEta -eq 0 ]; then
         #   correct=0
         #else 
            correct=1
         #fi
     else
	echo "Command not found. Please try again.";
     fi
    
done

correct=0
while [ $correct -eq 0 ]
do
    echo "Do you want to run in thesis mode? Yes/No?";
    read answer
    if [ $answer = "Yes" ]; then
	THESIS="thesis"
	echo "Running in Thesis mode ...";
	correct=1
    elif [ $answer = "No" ]; then
	THESIS="normal"
	echo "Running in Normal mode ...";
	correct=1
    else
	echo "Command not found. Please try again.";
    fi
done


correct=0
while [ $correct -eq 0 ]
do
    echo "Do you want to use MinBias Efficiencies only? Yes/No?";
    read answer
    if [ $answer = "Yes" ]; then
	OPTMINBIASEFF="MinBiasEffOnly"
	echo "Calculating MinBias Efficiecy only ...";
	correct=1
    elif [ $answer = "No" ]; then
	OPTMINBIASEFF="No"
	echo "Nothing to be done ...";
	correct=1
    else
	echo "Command not found. Please try again.";
    fi
done

echo "bla";

if [ $ONLYRESULTS -eq 0 ]; then
	echo "Hauptroutine stimmt -eq On the right path"
	if [ $ONLYCORRECTION -eq 0 ];  then
		correct=0
		while [ $correct -eq 0 ]
		do
			echo "Which fit do you want to do? CrystalBall or gaussian convoluted with an exponential function? CrystalBall/Gaussian?";
			read answer
			if [ $answer = "CrystalBall" ]; then
				echo "CrystalBall chosen ...";
				correct=1
				crystal=CrystalBall
			elif [ $answer = "Gaussian" ]; then
				echo "Gaussian chosen ...";
				correct=1	
				crystal=Gaussian
			else
				echo "Command not found. Please try again.";
			fi
		done
	fi
	
	echo "Hauptroutine stimmt -eq On the right path"
	correct=0
	while [ $correct -eq 0 ]
	do
		echo "Please check that you really want to process all cuts, otherwise change the CutSelection.log. Remember at first all gamma cutstudies will be carried out. Make sure that the standard cut is the first in the file. Continue? Yes/No?";
		read answer
		if [ $answer = "Yes" ]; then
			echo "Continuing ...";
			correct=1
		elif [ $answer = "No" ]; then
			echo "Aborting ...";
			exit
		else
			echo "Command not found. Please try again.";
		fi
	done

#Read the different cuts form the Cut selection log file


	CutSelections=`cat CutSelection.log`
	for cutSelection in $CutSelections; do
		if [ -d $cutSelection ]; then
			echo "CutSelection $cutSelection directory already exists, all files will be overwritten ";
			mkdir $cutSelection/$energy
		else
			mkdir $cutSelection
			mkdir $cutSelection/$energy
		fi

		if [ $ONLYCUTS -eq 0 ]; then 
			if [ -d $cutSelection/$energy/$Suffix ]; then
				echo "Graphical Output $Suffix directory already exists, all files will be overwritten ";
			else
				mkdir $cutSelection/$energy/$Suffix
			fi

			if [ $ONLYCORRECTION -eq 0 ]; then
				echo "CutSelection is $cutSelection";
				
				root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0\"\,\"$DataRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kFALSE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$THESIS\"\,$BinsPtPi0\,kFALSE\,$mode\)  
				Pi0dataRAWFILE=`ls $cutSelection/$energy/Pi0_data_GammaConvV1DalitzWithoutCorrection_*.root`           
	            
				if [ $MCFILE -eq 1 ]; then 
					
					root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0\"\,\"$MCRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$THESIS\"\,$BinsPtPi0\,kFALSE\,$mode\)   
					Pi0MCRAWFILE=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzWithoutCorrection_*.root`
					Pi0MCcorrectionFILE=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzCorrectionHistos_*.root`
					
				  
					if [ $MERGINGMC -eq 1 ]; then
					     if [ $addedSig -eq 1 ]; then
					     
							root -x -q -l -b TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0\"\,\"$MCRootFileAddSig\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"AddSig\"\,\"$THESIS\"\,$BinsPtPi0\,kTRUE\,$mode\)
							Pi0MCcorrection=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzCorrectionHistos_*.root`
							Pi0MCcorrectionAddSig=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzCorrectionHistosAddSig_*.root`
							root -b -x -q -l TaskV1/MergeEffiWithProperWeightingpPbDalitz.C\+\+\(\"$cutSelection\"\,\"Pi0\"\,\"$Suffix\"\,\"$energy\"\,\"$Pi0MCcorrection\"\,\"$cutSelection/$energy/Pi0_MC_GammaConvV1CorrectionHistosMinBias_$cutSelection.root\"\,\"$Pi0MCcorrectionAddSig\"\)
							
					      else
						root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0\"\,\"$MCRootFileBC\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"BC\"\,\"$THESIS\"\,$BinsPtPi0\,kFALSE\,$mode\)
						root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0\"\,\"$MCRootFileD\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"D\"\,\"$THESIS\"\,$BinsPtPi0\,kFALSE\,$mode\)
						Pi0MCcorrectionBCFILE=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzCorrectionHistosBC_*.root`
						Pi0MCcorrectionDFILE=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzCorrectionHistosD_*.root`
 						root -x -q -l -b  TaskV1/MergeEffiWithProperWeightingDalitz.C\(\"$DIRECTORY\"\,\"$cutSelection\"\,\"Pi0\",\"$Suffix\"\,\"$energy\"\,\"$Pi0MCcorrectionFILE\"\,\"$Pi0MCcorrectionBCFILE\",\"$Pi0MCcorrectionDFILE\"\)
					      fi
					fi
					
					root -x -q -l -b TaskV1/CompareMesonQuantitiesDalitz.C\+\+\(\"$Pi0dataRAWFILE\"\,\"$Pi0MCRAWFILE\"\,\"$cutSelection\"\,\"Pi0\"\,\"$Suffix\"\,\"$energy\"\)
				fi	
				
				if [ $Noeta -eq 0  ] && [ "$energy" != "PbPb_2.76TeV" ] && [ "$energy" != "pPb_5.023TeV" ]; then
					if [ $NoPi0Eta -eq 0 ]; then
						
						root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0EtaBinning\"\,\"$DataRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kFALSE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)               
						Pi0EtadataRAWFILE=`ls $cutSelection/$energy/Pi0EtaBinning_data_GammaConvV1DalitzWithoutCorrection_*.root`
					  
						if [ $MCFILE -eq 1 ]; then 

							
							root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0EtaBinning\"\,\"$MCRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)           
							Pi0EtaMCRAWFILE=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzWithoutCorrection_*.root`   
							Pi0EtaMCcorrectionFILE=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzCorrectionHistos_*.root`	
							
							if [ $MERGINGMC -eq 1 ]; then
							
								if [ $addedSig -eq 1 ]; then
									root -x -q -l -b TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0EtaBinning\"\,\"$MCRootFileAddSig\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"AddSig\"\,\"$THESIS\"\,$BinsPtPi0\,kTRUE\,$mode\)
									Pi0EtaMCcorrection=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzCorrectionHistos_*.root`
									Pi0EtaMCcorrectionAddSig=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzCorrectionHistosAddSig_*.root`
									root -b -x -q -l TaskV1/MergeEffiWithProperWeightingpPbDalitz.C\+\+\(\"$cutSelection\"\,\"Pi0EtaBinning\"\,\"$Suffix\"\,\"$energy\"\,\"$Pi0EtaMCcorrection\"\,\"$cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1CorrectionHistosMinBias_$cutSelection.root\"\,\"$Pi0EtaMCcorrectionAddSig\"\)
								else 
								root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0EtaBinning\"\,\"$MCRootFileBC\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"BC\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)
								root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Pi0EtaBinning\"\,\"$MCRootFileD\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"D\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)          
								Pi0EtaMCcorrectionBCFILE=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzCorrectionHistosBC_*.root`
								Pi0EtaMCcorrectionDFILE=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzCorrectionHistosD_*.root`
								root -x -q -l -b  TaskV1/MergeEffiWithProperWeightingDalitz.C\(\"$DIRECTORY\"\,\"$cutSelection\"\,\"Pi0EtaBinning\",\"$Suffix\"\,\"$energy\"\,\"$Pi0EtaMCcorrectionFILE\"\,\"$Pi0EtaMCcorrectionBCFILE\",\"$Pi0EtaMCcorrectionDFILE\"\)
								
								fi

							fi
							 root -x -q -l -b TaskV1/CompareMesonQuantitiesDalitz.C\+\+\(\"$Pi0EtadataRAWFILE\"\,\"$Pi0EtaMCRAWFILE\"\,\"$cutSelection\"\,\"Pi0EtaBinning\"\,\"$Suffix\"\,\"$energy\"\)
						fi
					fi
					
					root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Eta\"\,\"$DataRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kFALSE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)               
					EtadataRAWFILE=`ls $cutSelection/$energy/Eta_data_GammaConvV1DalitzWithoutCorrection_*.root`		

					if [ $MCFILE -eq 1 ]; then 
						root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Eta\"\,\"$MCRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)
						EtaMCRAWFILE=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzWithoutCorrection_*.root`
						EtaMCcorrectionFILE=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzCorrectionHistos_*.root`
									
						if [ $MERGINGMC -eq 1 ]; then
						
							if [ $addedSig -eq 1 ]; then
								root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Eta\"\,\"$MCRootFileAddSig\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"AddSig\"\,$BinsPtEta\,kTRUE\,$mode\)
								EtaMCcorrection=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzCorrectionHistos_*.root`
								EtaMCcorrectionAddSig=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzCorrectionHistosAddSig_*.root`
								root -b -x -q -l TaskV1/MergeEffiWithProperWeightingpPbDalitz.C\+\+\(\"$cutSelection\"\,\"Eta\"\,\"$Suffix\"\,\"$energy\"\,\"$EtaMCcorrection\"\,\"$cutSelection/$energy/Eta_MC_GammaConvV1CorrectionHistosMinBias_$cutSelection.root\"\,\"$EtaMCcorrectionAddSig\"\)
							else
						
							root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Eta\"\,\"$MCRootFileBC\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"BC\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)
							root -x -q -l -b  TaskV1/ExtractSignalDalitz.C\+\+\(\"Eta\"\,\"$MCRootFileD\"\,\"$cutSelection\"\,\"$Suffix\"\,\"kTRUE\"\,\"$energy\"\,\"$crystal\"\,\"$OPTMINBIASEFF\"\,\"D\"\,\"$THESIS\"\,$BinsPtEta\,kFALSE\,$mode\)               
							EtaMCcorrectionBCFILE=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzCorrectionHistosBC_*.root`
							EtaMCcorrectionDFILE=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzCorrectionHistosD_*.root`
							root -x -q -l -b  TaskV1/MergeEffiWithProperWeightingDalitz.C\(\"$DIRECTORY\"\,\"$cutSelection\"\,\"Eta\",\"$Suffix\"\,\"$energy\"\,\"$EtaMCcorrectionFILE\"\,\"$EtaMCcorrectionBCFILE\",\"$EtaMCcorrectionDFILE\"\)
							
							fi
						fi
		 				root -x -q -l -b TaskV1/CompareMesonQuantitiesDalitz.C\+\+\(\"$EtadataRAWFILE\"\,\"$EtaMCRAWFILE\"\,\"$cutSelection\"\,\"Eta\"\,\"$Suffix\"\,\"$energy\"\)
					fi
				
				fi
				if [ $MCFILE -eq 1 ] && [ $energy != "PbPb_2.76TeV" ]; then 
	  
				    root -x -q -l -b TaskQA/ElectronQAv1.C\+\+\(\"$DataRootFile\"\,\"$MCRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"$energy\"\,\"\"\,$mode\)


				    if [ $MERGINGMC -eq 1 ]; then
					#DataRootFileBC
					
					root -x -q -l -b TaskQA/ElectronQAv1.C\+\+\(\"$DataRootFileBC\"\,\"$MCRootFileBC\"\,\"$cutSelection\"\,\"$Suffix\"\,\"$energy\"\,\"BC\"\,$mode\)
					root -x -q -l -b TaskQA/ElectronQAv1.C\+\+\(\"$DataRootFileD\"\,\"$MCRootFileD\"\,\"$cutSelection\"\,\"$Suffix\"\,\"$energy\"\,\"D\"\,$mode\)
				
				
				    fi				
				fi  	
			fi

				
					
# 			
			Pi0dataRAWFILE=`ls $cutSelection/$energy/Pi0_data_GammaConvV1DalitzWithoutCorrection_*.root`             
			Pi0MCRAWFILE=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzWithoutCorrection_*.root`
			Pi0MCcorrectionFILE=`ls $cutSelection/$energy/Pi0_MC_GammaConvV1DalitzCorrectionHistos_*.root`

			if [ $NoPi0Eta -eq 0 ] && [ "$energy" != "PbPb_2.76TeV" ] && [ "$energy" != "pPb_5.023TeV" ]; then

			  Pi0EtadataRAWFILE=`ls $cutSelection/$energy/Pi0EtaBinning_data_GammaConvV1DalitzWithoutCorrection_*.root`
			  Pi0EtaMCRAWFILE=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzWithoutCorrection_*.root`   
			  Pi0EtaMCcorrectionFILE=`ls $cutSelection/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzCorrectionHistos_*.root`

			fi

			if [ $Noeta -eq 0 ] && [ "$energy" != "PbPb_2.76TeV" ] && [ "$energy" != "pPb_5.023TeV" ]; then
			
			EtadataRAWFILE=`ls $cutSelection/$energy/Eta_data_GammaConvV1DalitzWithoutCorrection_*.root`
			EtaMCRAWFILE=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzWithoutCorrection_*.root`
			EtaMCcorrectionFILE=`ls $cutSelection/$energy/Eta_MC_GammaConvV1DalitzCorrectionHistos_*.root`
			
			fi

			
			if [ -n $Pi0dataRAWFILE ] && [ -n $Pi0MCcorrectionFILE ]; then
			
			       	CorrectSignal $Pi0dataRAWFILE $Pi0MCcorrectionFILE $cutSelection $Suffix Pi0 kFALSE $MCSample
			else 
				PARTLY=1
			fi

			if [ -n $Pi0MCRAWFILE ] && [ -n $Pi0MCcorrectionFILE ]; then
				CorrectSignal $Pi0MCRAWFILE $Pi0MCcorrectionFILE $cutSelection $Suffix Pi0 kTRUE $MCSample
			else 
				PARTLY=1
			fi


			if [ $Noeta -eq 0 ] && [ "$energy" != "PbPb_2.76TeV" ] && [ "$energy" != "pPb_5.023TeV" ]; then
				
				if [ -n $Pi0EtadataRAWFILE ] && [ -n $Pi0EtaMCcorrectionFILE ] && [ $NoPi0Eta -eq 0 ]; then
					CorrectSignal $Pi0EtadataRAWFILE $Pi0EtaMCcorrectionFILE $cutSelection $Suffix Pi0EtaBinning kFALSE $MCSample
				else 
						PARTLY=1
				fi

				if [ -n $Pi0EtaMCRAWFILE ] && [ -n $Pi0EtaMCcorrectionFILE ] && [ $NoPi0Eta -eq 0 ]; then
					CorrectSignal $Pi0EtaMCRAWFILE $Pi0EtaMCcorrectionFILE $cutSelection $Suffix Pi0EtaBinning kTRUE $MCSample
				else 
						PARTLY=1
				fi

				if [ -n $EtadataRAWFILE ] && [ -n $EtaMCcorrectionFILE ]; then
					CorrectSignal $EtadataRAWFILE $EtaMCcorrectionFILE $cutSelection $Suffix Eta kFALSE $MCSample
				else 
						PARTLY=1
				fi
				if [ -n $EtaMCRAWFILE ] && [ -n $EtaMCcorrectionFILE ]; then
					CorrectSignal $EtaMCRAWFILE $EtaMCcorrectionFILE $cutSelection $Suffix Eta kTRUE $MCSample
				else 
						PARTLY=1
				fi
						
			fi
                        if [ "$energy" == "5TeV2017" ]; then
                                root -x -q -l -b TaskQA/ElectronQAv1.C\+\+\(\"$DataRootFile\"\,\"$MCRootFile\"\,\"$cutSelection\"\,\"$Suffix\"\,\"$energy\"\,\"\"\,$mode\)
                        fi

		fi

		NORMALCUTS=`expr $NORMALCUTS + 1`
	done

	root -x -q -l -b TaskV1/CutStudiesDalitzOverview.C\(\"CutSelection.log\"\,\"$Suffix\"\,\"Pi0\"\,\"kFALSE\"\,\"$OPTMINBIASEFF\"\,\"$energy\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"Dalitz\"\,$mode\)
	root -x -q -l -b TaskV1/CutStudiesDalitzOverview.C\(\"CutSelection.log\"\,\"$Suffix\"\,\"Pi0\"\,\"kTRUE\"\,\"$OPTMINBIASEFF\"\,\"$energy\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"Dalitz\"\,$mode\)
	

	if [ $Noeta -eq 0 ] && [ "$energy" != "PbPb_2.76TeV" ] && [ "$energy" != "pPb_5.023TeV" ]; then
		if [ $NoPi0Eta -eq 0 ]; then
			root -x -q -l -b TaskV1/CutStudiesDalitzOverview.C\(\"CutSelection.log\"\,\"$Suffix\"\,\"Pi0EtaBinning\"\,\"kFALSE\"\,\"$OPTMINBIASEFF\"\,\"$energy\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"Dalitz\"\,$mode\)
			root -x -q -l -b TaskV1/CutStudiesDalitzOverview.C\(\"CutSelection.log\"\,\"$Suffix\"\,\"Pi0EtaBinning\"\,\"kTRUE\"\,\"$OPTMINBIASEFF\"\,\"$energy\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"Dalitz\"\,$mode\)
		fi
		root -x -q -l -b TaskV1/CutStudiesDalitzOverview.C\(\"CutSelection.log\"\,\"$Suffix\"\,\"Eta\"\,\"kFALSE\"\,\"$OPTMINBIASEFF\"\,\"$energy\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"Dalitz\"\,$mode\)
		root -x -q -l -b TaskV1/CutStudiesDalitzOverview.C\(\"CutSelection.log\"\,\"$Suffix\"\,\"Eta\"\,\"kTRUE\"\,\"$OPTMINBIASEFF\"\,\"$energy\"\,\"$NAMECUTSTUDIES\"\,$NORMALCUTS\,0\,\"Dalitz\"\,$mode\)
	fi
fi		


if [ "$energy" != "PbPb_2.76TeV" ]; then
   if [ "$energy" != "pPb_5.023TeV" ]; then

      

      Pi0dataCorr=`ls $standardCutMeson/$energy/Pi0_data_GammaConvV1DalitzCorrection_*.root`
      Pi0MCCorr=`ls $standardCutMeson/$energy/Pi0_MC_GammaConvV1DalitzCorrection_*.root`


      if [ $Noeta -eq 0 ]; then

	  if [ $NoPi0Eta -eq 0 ]; then
	    Pi0EtadataCorr=`ls $standardCutMeson/$energy/Pi0EtaBinning_data_GammaConvV1DalitzCorrection_*.root`
	    Pi0EtaMCCorr=`ls $standardCutMeson/$energy/Pi0EtaBinning_MC_GammaConvV1DalitzCorrection_*.root`
	  fi

	EtadataCorr=`ls $standardCutMeson/$energy/Eta_data_GammaConvV1DalitzCorrection_*.root`
	EtaMCCorr=`ls $standardCutMeson/$energy/Eta_MC_GammaConvV1DalitzCorrection_*.root`
         
      fi
      
      CallPi0Data=\"$Pi0dataCorr\"\,\"$EtadataCorr\"\,\"$standardCutMeson\"\,\"$Suffix\"\,\"kFALSE\"\,\"Levy\"\,\"\"\,\"$energy\"\,\"$directphoton\"\,\"$OPTMINBIASEFF\"\,\"$THESIS\"
      CallPi0MC=\"$Pi0MCCorr\"\,\"$EtaMCCorr\"\,\"$standardCutMeson\"\,\"$Suffix\"\,\"kTRUE\"\,\"Levy\"\,\"\"\,\"$energy\"\,\"$directphoton\"\,\"$OPTMINBIASEFF\"\,\"$THESIS\"
      CallPi0EtaData=\"$Pi0EtadataCorr\"\,\"$EtadataCorr\"\,\"$standardCutMeson\"\,\"$Suffix\"\,\"kFALSE\"\,\"Levy\"\,\"same\"\,\"$energy\"\,\"$directphoton\"\,\"$OPTMINBIASEFF\"\,\"$THESIS\"
      CallPi0EtaMC=\"$Pi0EtaMCCorr\"\,\"$EtaMCCorr\"\,\"$standardCutMeson\"\,\"$Suffix\"\,\"kTRUE\"\,\"Levy\"\,\"same\"\,\"$energy\"\,\"$directphoton\"\,\"$OPTMINBIASEFF\"\,\"$THESIS\"


      if [ -n $EtadataCorr ] && [ -n $Pi0dataCorr ] && [ -n $SysErrFiledata ]; then

         CreateFinalResults $CallPi0Data ;
   
      else 
         PARTLY=1
      fi

      echo "have done this" 

      if [ -n $EtaMCCorr ] && [ -n $Pi0MCCorr ] && [ -n $SysErrFileMC ]; then
         CreateFinalResults $CallPi0MC
         #CreateGammaFinalResults $GammaPi0MCCorr $Pi0MCCorr $standardCutMeson $Suffix Pi0 kTRUE $directphoton;
      else 
         PARTLY=1
      fi


      if [ $Noeta -eq 0 ]; then
											 
         if [ -n $EtadataCorr ] && [ -n $Pi0EtadataCorr ] && [ -n $SysErrFiledata ] && [ $NoPi0Eta -eq 0 ]; then
            CreateFinalResults $CallPi0EtaData
      # 	if [ $Noeta -eq 0 ]; then CreateGammaFinalResults $GammaPi0EtadataCorr $Pi0EtadataCorr $standardCutMeson $Suffix Pi0EtaBinning kFALSE $directphoton; fi
         else 
            PARTLY=1
         fi
         
         if [ -n $EtaMCCorr ] && [ -n $Pi0EtaMCCorr ] && [ -n $SysErrFileMC ] && [ $NoPi0Eta -eq 0 ]; then
            CreateFinalResults $CallPi0EtaMC
      # if [ $Noeta -eq 0 ]; then CreateGammaFinalResults $GammaPi0EtaMCCorr $Pi0EtaMCCorr $standardCutMeson $Suffix Pi0EtaBinning kTRUE $directphoton; fi
         else 
            PARTLY=1
         fi
      fi
   fi
fi
if [ "$energy" == "pPb_5.023TeV" ]; then

   Pi0dataCorr=`ls $standardCutMeson/$energy/Pi0_data_GammaConvV1DalitzCorrection_*.root`

   #if [ $Noeta -eq 0 ]; then
   #   Pi0EtadataCorr=`ls $standardCutMeson/$energy/Pi0EtaBinning_data_GammaConvV1DalitzCorrection_*.root`
   #   EtadataCorr=`ls $standardCutMeson/$energy/Eta_data_GammaConvV1DalitzCorrection_*.root`
   #fi
   
   ProduceFinalResultspPb $Pi0dataCorr $standardCutMeson $Suffix Pi0 Levy $energy
   # if [ $Noeta -eq 0 ]; then
   #	ProduceFinalResultspPb $Pi0EtadataCorr $standardCutMeson $Suffix Pi0EtaBinning Levy $energy $Pi0EtadataCorr
   #   ProduceFinalResultspPb $EtadataCorr $standardCutMeson $Suffix Eta Levy $energy $Pi0EtadataCorr
   #fi
	
fi 

exit 
