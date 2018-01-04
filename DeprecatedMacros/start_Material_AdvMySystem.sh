#! /bin/bash
#
#
#
# This script gests as input a directory where the AnalysisResults root file is stored
# , it also needs the desired output directory where the produced root files are put.
# If nothing is given it will use ./ for the input directory and ./Output for the output
#
#Input 1: Root file to analyze Default: AnalyisResults
#Input 2: Input directory  Default:$PWD 
#Input 3: Output directory Default: $PWD/Results  (directory will be created if it does not exist)
#


PROGNAME=$0


# function Usage () 
# {
# 	echo -e "
# This Script is provided by the Gamma Conversion Group of ALICE 
# the main developers are 
# 	Friederike Bock \t friederike.bock@cern.ch
# 	Kathrin Koch \t\t kkoch@physi.uni-heidelberg.de
# 	Ana Marin \t\t marin@physi.uni-heidelberg.de  
# 	Kenneth Aamoth \t\t work.kenneth.aamodt@gmail.com
# 			
# If there are any complications with running this script do not hesitate to connect them.
# 			
# 
# How to use this script? 
# $PROGNAME -h  \t\t\t\t\t usage will be displayed 
# $PROGNAME --help \t\t\t\t usage will be displayed 
# 
# This script needs as a basis the output from the GammaConversion-Software which is provided with the Aliroot software, both the 'data.root' and the 'MC.root' have to be output files of this. The script then will check which cutnumbers are in these files and will ask you which exactly you want to analyse, furthermore you have to set a standard cut which has to be always the first in the CutLogFile. Afterwards it will be distinguished between Phojet and Pythia, for which the corresponding directories will be created. One of the main outputs will be the  
# 		
# 	"
# 	exit 
# }

if [ "$#" == "0" ]; then
	echo -e "Please call either" 
	echo -e "$PROGNAME -h  \t\t\t\t\t usage will be displayed"
	echo -e "$PROGNAME --help \t\t\t\t usage will be displayed"
	exit
fi

Suffix=eps

if [ $1 = "-h" ] || [ $1 = "--help" ]; then
	Usage
else 
	DataRootFile=$1
	MCRootFile=$2
	DataRootFileSec=$3
	MCRootFileSec=$4
	Suffix=$5;
fi

PARTLY=0
DATAFILE=1
MCFILE=1
RES=0

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

ALLFILES=1

echo "Testing whether all files are there ..."

if [ -f PlottingGammaConversionAdditional.h ]; then
	echo -e "PlottingGammaConversionAdditional.h exists ... \t\t yes" 
else
	echo -e "PlottingGammaConversionAdditional.h exists ... \t\t no" 	
	ALLFILES=0
fi

if [ -f PlottingGammaConversionHistos.h ]; then
	echo -e "PlottingGammaConversionHistos.h exists ... \t\t yes" 
else
	echo -e "PlottingGammaConversionHistos.h exists ... \t\t no" 	
	ALLFILES=0
fi

if [ -f PhotonCharacteristicsAdv.C ]; then
	echo -e "PhotonCharacteristicsAdv.C exists ... \t\t yes" 
else
	echo -e "PhotonCharacteristicsAdv.C exists ... \t\t no" 	
	ALLFILES=0
fi

if [ -f PhotonResolutionAdv.C ]; then
	echo -e "PhotonResolutionAdv.C exists ... \t\t yes" 
else
	echo -e "PhotonResolutionAdv.C exists ... \t\t no" 	
	ALLFILES=0
fi

if [ -f MappingMaterialAdv.C ]; then
	echo -e "MappingMaterialAdv.C exists ... \t\t yes" 
else
	echo -e "MappingMaterialAdv.C exists ... \t\t no" 	
	ALLFILES=0
fi

if [ -f FittingGammaConversion.h ]; then
	echo -e "FittingGammaConversion.h exists ... \t\t\t\t yes" 
else
	echo -e "FittingGammaConversion.h exists ... \t\t\t\t no" 	
	ALLFILES=0
fi

if [ -f ALICE_logo.eps ]; then
	echo -e "ALICE_logo.eps exists ... \t\t\t\t yes" 
else
	echo -e "ALICE_logo.eps exists ... \t\t\t\t no" 	
	ALLFILES=0
fi

if [ $ALLFILES -eq 0 ]; then
	echo "There are files missing, please put them in the working directory."
	exit
else
	echo "All files exist. Continuing ..."
fi

# . alilogin trunk

correct=0
while [ $correct -eq 0 ]
do
	echo "Do you want to take an already exitsting CutSelection.log-file. Yes/No"
	read answer
	if [ $answer = "Yes" ]; then
		echo Chosen already existing logfile ...;
		cat CutSelection.log
		cat CutSelectionSecHad.log
		correct=1
	elif [ $answer = "No" ]; then
		echo Aborting ...;
		root -b -q MakeCutLogGamConv.C\(\"$DataRootFile\"\,\"CutSelection.log\"\)
		root -b -q MakeCutLogSecHad.C\(\"$DataRootFileSec\"\,\"CutSelectionSecHad.log\"\)
		correct=1
	else
		echo Command not found. Please try again.;
	fi
done

correct=0
while [ $correct -eq 0 ]
do
	echo "Which collision system do you want to process? 7TeV (pp@7TeV), 2.76TeV (pp@2.76TeV), 900GeV (pp@900GeV), HI (PbPb@2.76GeV)"
	read energy
	if [ $energy = "7TeV" ]; then
		echo Do you want to produce conference plots? Yes/No?;
		read answer
		if [ $answer = "Yes" ]; then
			echo Will produce conference plots ...;
			Conference="conference"
			Con=1
			correct=1
		elif [ $answer = "No" ]; then
			echo No conference plots will be produced ...;
			Conference="No"
			Con=0
			correct=1
		else
			echo Command not found. Please try again.;
		fi
	elif [ $energy = "2.76TeV" ]; then
		Noeta=1
		echo Do you want to produce conference plots? Yes/No?;
		read answer
		if [ $answer = "Yes" ]; then
			echo Will produce conference plots ...;
			Conference="conference"
			Con=1
			correct=1
		elif [ $answer = "No" ]; then
			echo No conference plots will be produced ...;
			Conference="No"
			Con=0
			correct=1
		else
			echo Command not found. Please try again.;
		fi
	elif [ $energy = "900GeV" ]; then
		Noeta=1
		echo Do you want to produce conference plots? Yes/No?;
		read answer
		if [ $answer = "Yes" ]; then
			echo Will produce conference plots ...;
			Conference="conference"
			Con=1
			correct=1
		elif [ $answer = "No" ]; then
			echo No conference plots will be produced ...;
			Conference="No"
			Con=0
			correct=1
		else
			echo Command not found. Please try again.;
		fi
	elif [ $energy = "HI" ]; then
		HIRUN=1
		correct=1
	else
	    echo Command not found. Please try again.;
	fi
done

correct=0
while [ $correct -eq 0 ]
do
	echo "Do you want to run resolution? Yes/No";
	read answer
	if [ $answer = "Yes" ]; then
		echo "Resolution will be run."
		RES=1
		correct=1
	elif [ $answer = "No" ]; then
		echo "Resolution will not be run."
		RES=0
		correct=1
	else
		echo Command not found. Please try again.;
		
	fi
done

SECHAD=0
correct=0
while [ $correct -eq 0 ]
do
	echo "Do you want to secondary hadronic interactions? Yes/No";
	read answer
	if [ $answer = "Yes" ]; then
		echo "Secondary hadronic interactions will be run."
		SECHAD=1
		correct=1
	elif [ $answer = "No" ]; then
		echo "Secondary hadronic interactions will not be run."
		SECHAD=0
		correct=1
	else
		echo Command not found. Please try again.;
		
	fi
done


correct=0
while [ $correct -eq 0 ]
do
	echo "What kind of MC file is it? Pythia6/Phojet/Pythia8/Hijing/None?";
	read answer
	if [ $answer = "Phojet" ]; then
		echo "Phojet has been specified."
		Generator=Phojet
		correct=1
	elif [ $answer = "Pythia6" ]; then
		echo "Pythia has been specified."
		Generator=Pythia6
		correct=1
	elif [ $answer = "Pythia8" ]; then
		echo "Pythia has been specified."
		Generator=Pythia8
		correct=1
	elif [ $answer = "Hijing" ]; then
		echo "Hijing has been specified."
		Generator=Hijing
		correct=1
	elif [ $answer = "None" ]; then
		echo "None has been specified."
		Generator=None
		correct=1
	else
		echo Command not found. Please try again.;
	fi
done

correct=0
while [ $correct -eq 0 ]
do
	echo "Specify the Periods your using! All/LHC???";
	read answer
	if [ $answer = "All" ]; then
		echo "Merged Periods"
		Period=All
		correct=1
	elif [ $answer = "" ]; then
		echo "No valid answer."
		correct=0
	else 
		echo "You specified periods $answer."
		Period=$answer
		correct=1
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
	echo "Please check that you really want to process all cuts, otherwise change the CutSelection.log. Remember at first all gamma cutstudies will be carried out. Make sure that the standard cut is the first in the file. Continue? Yes/No?";
	read answer
	if [ $answer = "Yes" ]; then
		echo Continuing ...;
		correct=1
	elif [ $answer = "No" ]; then
		echo Aborting ...;
		exit
	else
		echo Command not found. Please try again.;
	fi
done

if [ -d $Generator ]; then
	echo "Directory for specified generator already exists"
else 
	mkdir $Generator
fi	


cp -d PlottingGammaConversionAdditional.h $Generator/
cp -d PlottingGammaConversionHistos.h $Generator/
cp -d PhotonResolutionAdv.C $Generator/
cp -d MappingMaterialTree.C $Generator/
cp -d MappingMaterialHadronicAdv.C $Generator/
cp -d ALICE_logo.eps $Generator/
cp -d CutSelection.log $Generator/
cp -d CutSelectionSecHad.log $Generator/
cp -d FittingGammaConversion.h $Generator/

cd $Generator
#Read the different cuts form the Cut selection log file
CutSelections=`cat CutSelection.log`
for cutSelection in $CutSelections; do
	if [ -d $cutSelection ]; then
		echo CutSelection $cutSelection directory already exists, all files will be overwritten ;
		mkdir $cutSelection/$energy
	else
		mkdir $cutSelection
		mkdir $cutSelection/$energy
	fi

	if [ -d $cutSelection/$energy/$Suffix ]; then
		echo Graphical Output $Suffix directory already exists, all files will be overwritten ;
	else
		mkdir $cutSelection/$energy/$Suffix
	fi

	if [ $MCFILE -eq 1 ]; then
		root -l -x -q -b  MappingMaterialTree.C\+\+\(\"$DataRootFile\"\,\"$MCRootFile\"\,\"$cutSelection\"\,\"MappingMaterialOverviewAdv_$cutSelection\"\,\"$Suffix\"\,\"$Conference\"\,\"$energy\"\,\"$Generator\"\,\"$Period\"\,\"$THESIS\"\)		
	fi 
done

CutSelections=`cat CutSelectionSecHad.log`
for cutSelection in $CutSelections; do
	if [ -d $cutSelection ]; then
		echo CutSelection $cutSelection directory already exists, all files will be overwritten ;
		mkdir $cutSelection/$energy
	else
		mkdir $cutSelection
		mkdir $cutSelection/$energy
	fi

	if [ -d $cutSelection/$energy/$Suffix ]; then
		echo Graphical Output $Suffix directory already exists, all files will be overwritten ;
	else
		mkdir $cutSelection/$energy/$Suffix
	fi



	if [ $SECHAD -eq 1 ]; then 
		root -x -q -l -b MappingMaterialHadronicAdv.C\+\+\(\"$DataRootFileSec\"\,\"$MCRootFileSec\"\,\"$cutSelection\"\,\"MappingMaterialHadronic3_$cutSelection\"\,\"$Suffix\"\,\"\"\,\"$energy\"\,\"$Generator\"\,\"$Period\"\)
	fi 

done

