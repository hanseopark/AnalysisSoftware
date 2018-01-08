#! /bin/bash

# This macro copies runwise GammaConvV1_<trainConfig> and/or AnalysisResults.root files from grid and puts them in the correct directories
#
# the files will be saved in these places:
#    $DIR/GammaConv/$SYSTEM/$PERIOD/$run/GammaConvV1_<trainconfig>.root
#    $DIR/PhotonQA/$SYSTEM/$PERIOD/$run/AnalysisResults.root
# 
#                or     .../<runNumber>/Stage*/*/<file>
#
# In the basic settings, one has to specify the base directory "DIR"
# and the current information on the datasets
#
# Available functions are:
#   MergeList()
#   ChangeStructure()
#   MainRoutine()
#
# Have to check by hand the merging success; if stage_i output exists it doesn't mean that stage_i merging was completely successful
# Merged file: include pass/runlist in name
#
# TO DO:
# dry-run option
# check if file is corrupted
# Error / warning and algorithm if folder exists already
# automatic train selection like in old macro
# create info.txt 
# meaningful filename for overall merged file via ChangeStructureToStandard
# possibility to download part of the data and merge afterwards the whole runlist
#

MergeList(){  # $1:list   $2:target path   $3:target name    
    
    fileList=`cat $1`
    hadd -f $2/$3 $fileList

}

ChangeStructure(){  # $1: file path + name
#"$BASEDIR/$SYSTEM/$PERIOD/$fileName" 

if [ -f "$1" ]; then
    echo "Change Structure to standard..."
    config=`echo $1 | cut -d "_" -f 2 | cut -d "." -f 1`
    echo "detected trainconfig: " $config
    root -l -b -q -x ChangeStructureToStandard.C\(\"$1\"\,\"$1\"\,\"GammaConvV1_$config\"\)
else
    echo "no file for ChangeStructureToStandard"
fi

}

DownloadMerged(){ # $1: source, $2:target , $3:change structure to standard

    echo "copy $1 to $2 ..."
    alien_cp alien:$1 file:$2
    if [ -f "$2" ]; then
	if [ "$3" == "yes" ];then
	    ChangeStructure $2
	fi
    else
	echo "download failed"
    fi
}



DownloadRunwise(){

    # arguments:
    # $1:  "GammaConv" or "PhotonQA"
    # $2:  "ESD" or "AOD"
    # $3:  "LHC15o"
    # $4:  "337_20171030-1322"
    # $5:  "data" or "sim"
    # $6:  "pass1" / "pass1_pidfix" / "pass3_lowIR_pidfix" / ""

    # merge all runs in the end?
    declare -i doMerge=1      # merge runs: 1 for yes, 0 for no
    declare -i doMergeRun=1   # merge sub-run files: 1 for yes, 0 for no

    DIR=/home/meike/analysis/data/GridOutput
    TASK=$1
    BASEDIR=$DIR/$1
    DATATYPE=$2
    SYSTEM=PbPb/$2
    PERIOD=$3
    TRAIN=$4
    TYPE=$5
    PASS=$6

    # settings according to the pass and runlist
    if [ "$PASS" == "pass1" ];then
	RUN=(246994 246991 246989 246984 246982 246980 246948 246945 246928 246871 246870 246867 246865 246851 246847 246846 246845 246844 246810 246809 246808 246807 246805 246804 246766 246765 246763 246760 246759 246758 246757 246751 246750 246676 246675 246495 246493 246488 246487 246434 246431 246428 246424 246276 246275 246272 246271 246225 246222 246217 246185 246182 246181 246180 246178 246153 246152 246151 246115 246113 246089 246087 246053 246052 246049 246048 246042 246037 246036 246012 246003 246001 245954 245952 245949 245923 245833 245831 245829 245793 245785 245775 245766 245759 245752 245738 245731 245729 245705 245702 245700 245692 245683);
    elif [ "$PASS" == "pass1_pidfix" ];then
	RUN=(245554 245545 245544 245543 245542 245540 245535 245507 245505 245504 245501 245497 245496 245454 245452 245450 245446 245441 245439 245411 245410 245409 245407 245401 245397 245396 245353 245349 245347 245346 245345 245343 245259 245232 245231 245152 245151 245146 245145);	
    elif [ "$PASS" == "pass3_lowIR_pidfix" ];then
	RUN=(246392 246391 246390 245068 245066 245064 244983 244982 244980 244975 244918);
    elif [ "$PASS" == "" ];then
	echo "no PASS definded because using MC"
    else
	echo "PASS $PASS not defined"
    fi

    # settings according to data or MC
    if [ "$TYPE" == "data" ];then
	YEAR="2015"
	PREFIX="000"
	if [ "$DATATYPE" == "ESD" ];then
	    AODNO=""
	    TRAINDIR="GA_PbPb"
	elif [ "$DATATYPE" == "AOD" ];then
	    AODNO="AOD186"
	    TRAINDIR="GA_PbPb_AOD"
	else
	    echo "DATATYPE $DATATYPE not defined"
	fi
    elif [ "$TYPE" == "sim" ];then
	YEAR="2016"
	PREFIX=""
	PASS=""
	if [ "$DATATYPE" == "ESD" ];then
	    AODNO=""
	    TRAINDIR="GA_PbPb_MC"
	elif [ "$DATATYPE" == "AOD" ];then
	    AODNO="AOD188"
	    TRAINDIR="GA_PbPb_MC_AOD"
	else
	    echo "DATATYPE $DATATYPE not defined"
	fi
    else
	echo "TYPE $TYPE not defined"
    fi

    # settings according to which task was run and which trainconfigs were used
    if [ "$TASK" == "GammaConv" ];then
	FILENAMES=("GammaConvV1_246.root" "GammaConvV1_247.root");   #  ("GammaConvV1_246.root") or ("GammaConvV1_245.root" "GammaConvV1_246.root")
	declare -i doChStr=1                  # 1 for yes, 0 for no
    elif [ "$TASK" == "PhotonQA" ];then
	FILENAMES=("AnalysisResults.root");
	declare -i doChStr=0   
    else
	echo "TASK $TASK not defined"
    fi

    # end of settings

    # print settings ...
    echo "YOU HAVE CHOSEN THE FOLLOWING SETTINGS"
    echo "DIR = $DIR"
    echo "BASEDIR = $BASEDIR"
    echo "DATATYPE = $DATATYPE"
    echo "SYSTEM = $SYSTEM"
    echo "PERIOD = $PERIOD"
    echo "TRAIN = $TRAIN"
    echo "TYPE = $TYPE"
    echo "TRAINDIR = $TRAINDIR"
    echo "YEAR = $YEAR"
    echo "PREFIX = $PREFIX"
    echo "PASS = $PASS"
    echo "AODNO = $AODNO"
    echo "FILENAMES = ${FILENAMES[@]}"
    echo "doMerge = $doMerge"
    echo "doMergeRun = $doMergeRun"
    echo "Runs: ${RUN[@]}"
    echo "END OF SETTINGS"

    declare -i r=0  # number of runs
    declare -i k=0  # number of downloaded runs of one filename
    declare -i l=0  # number of unmerged files for one run
    
    for FILENAME in "${FILENAMES[@]}"; do
	
	if [ -f FilesToMerge.txt ]; then rm FilesToMerge.txt; fi 
	if [ -f UnmergedFilesToMerge.txt ]; then rm UnmergedFilesToMerge.txt; fi
	if [ -f MergedRuns.txt ]; then rm MergedRuns.txt; fi
	
	r=0
	k=0
	
	for run in "${RUN[@]}"; do
	    r=r+1
	    l=0     # for one specific run: number of unmerged files
	    if [ -f UnmergedFilesToMerge.txt ]; then rm UnmergedFilesToMerge.txt; fi # clear list before every new run
	    echo "**********************************************************************************************************"
	    echo "Run $r of ${#RUN[@]}: $run";
	    SOURCEDIR=/alice/$TYPE/$YEAR/$PERIOD/$PREFIX$run/$PASS/$AODNO/PWGGA/$TRAINDIR/$TRAIN
	    echo "source directory:    $SOURCEDIR"
	    echo "search for file:     $FILENAME";
	    fileName=`alien_ls $SOURCEDIR/$FILENAME`                  
	    if [ "$fileName" = "$FILENAME" ]; then
		echo "merged file exists:  $fileName"
		OUTPUTDIR=$BASEDIR/$SYSTEM/$PERIOD/$run/
		echo "output directory:    $OUTPUTDIR"
		mkdir -p $OUTPUTDIR
		alien_cp alien:$SOURCEDIR/$fileName file:$OUTPUTDIR/   
		if [ -f $OUTPUTDIR/$fileName ]; then  # only if ile exists
		    echo "copied               $SOURCEDIR/$fileName "
		    echo "to                   $OUTPUTDIR" 
		    echo -n " $OUTPUTDIR/$fileName" >> FilesToMerge.txt   
		    k=k+1
		    echo $run >> MergedRuns.txt
		    if [ "$doChStr" -eq 1 ]; then
			ChangeStructure "$OUTPUTDIR/$fileName" 
		    fi
		else
		    echo "ERROR: download failed"
		fi
	    else
		echo "   no merged file found"
		for((i=5; i>=0; ((i=i-1)) )); do       # search first stage_5, last stage_1
		    echo "   search Stage $i ..."
		    if [ "$i" -eq 0 ]; then SOURCEDIR="/alice/$TYPE/$YEAR/$PERIOD/$PREFIX$run/$PASS/$AODNO/PWGGA/$TRAINDIR/$TRAIN/"
		    else SOURCEDIR="/alice/$TYPE/$YEAR/$PERIOD/$PREFIX$run/$PASS/$AODNO/PWGGA/$TRAINDIR/$TRAIN/Stage_$i"; fi
		    unmergedFolders=`alien_ls $SOURCEDIR`
		    notFoundString="no such file or directory"
		    if [[ $unmergedFolders = *$notFoundString* ]]; then
			echo "   stage $i not found"
		    else   # if stage_i folder exists, it does not mean that merging stage_i has completed !!!
			echo "   folders found in stage $i "
			for unmergedFolder in $unmergedFolders; do
			    if [ "${unmergedFolder//[0-9]}" = "" ]; then #  consists only of digits. Can be also .xml file
				echo "      download from $unmergedFolder ..."
				if [ "$i" -eq 0 ]; then SOURCEDIR="/alice/$TYPE/$YEAR/$PERIOD/$PREFIX$run/$PASS/$AODNO/PWGGA/$TRAINDIR/$TRAIN/$unmergedFolder"
				else SOURCEDIR="/alice/$TYPE/$YEAR/$PERIOD/$PREFIX$run/$PASS/$AODNO/PWGGA/$TRAINDIR/$TRAIN/Stage_$i/$unmergedFolder"; fi
				echo "      source directory:    $SOURCEDIR"
				fileName=`alien_ls $SOURCEDIR/$FILENAME`
				if [ "$fileName" = "$FILENAME" ]; then
				    echo "      unmerged file exists:$fileName"
				    OUTPUTDIR=$BASEDIR/$SYSTEM/$PERIOD/$run/$unmergedFolder
				    echo "      output directory:    $OUTPUTDIR"
				    mkdir -p $OUTPUTDIR
				    alien_cp alien:$SOURCEDIR/$fileName file:$OUTPUTDIR/
				    if [ -f $OUTPUTDIR/$fileName ]; then  # only if file exists:
					echo "      copied               $SOURCEDIR/$fileName "
					echo "      to                   $OUTPUTDIR"
					echo -n " $OUTPUTDIR/$fileName" >> UnmergedFilesToMerge.txt
					l=l+1
				    else
					echo "      ERROR: copying failed"
				    fi
				else
				    echo "      ERROR: found $fileName, not $FILENAME"
				fi
			    else
				echo "      do not download from $unmergedFolder ..."
			    fi
			done
			break   # download only the 'highest stage'
		    fi
		done
		echo "downloaded $l unmerged file(s) for run $run"
		if [ "$l" -gt 1 ]; then
		    if [ "$doMergeRun" -eq 1 ]; then
			echo "merge $l files for run $run ..."
			MergeList "UnmergedFilesToMerge.txt" $BASEDIR/$SYSTEM/$PERIOD/$run $fileName
			if [ -f $BASEDIR/$SYSTEM/$PERIOD/$run/$fileName ]; then  # only if ile exists
			    echo -n " $BASEDIR/$SYSTEM/$PERIOD/$run/$fileName" >> FilesToMerge.txt
			    k=k+1
			    echo $run >> MergedRuns.txt
			    if [ "$doChStr" -eq 1 ]; then
				ChangeStructure "$BASEDIR/$SYSTEM/$PERIOD/$run/$fileName"
			    fi
			else
			    echo "ERROR: merging failed"
			fi
		    fi
		elif [ "$l" -eq 1 ]; then
		    echo "Only one file, nothing to merge"
		    singleFile=`cat UnmergedFilesToMerge.txt`
		    mv $singleFile $BASEDIR/$SYSTEM/$PERIOD/$run/
		    if [ -f $BASEDIR/$SYSTEM/$PERIOD/$run/$fileName ]; then  # only if ile exists
			echo "moved $singleFile "
			echo "to    $BASEDIR/$SYSTEM/$PERIOD/$run/$fileName"
			echo -n " $BASEDIR/$SYSTEM/$PERIOD/$run/$fileName" >> FilesToMerge.txt
			k=k+1
			echo $run >> MergedRuns.txt
			if [ "$doChStr" -eq 1 ]; then
			    ChangeStructure "$BASEDIR/$SYSTEM/$PERIOD/$run/$fileName"
			fi
		    fi
		else
		    echo "nothing to merge"
		fi
	    fi
	done # end runs loop
	
	
	echo "**********************************************************************************************************"
	echo "**********************************************************************************************************"
	echo "have $k of ${#RUN[@]} run files $FILENAME"
	
	if [ "$doMerge" -eq 1 ]; then
	    if [ "$k" -gt 1 ]; then
		echo "merge files from $k runs..."
		MergeList "FilesToMerge.txt" $BASEDIR/$SYSTEM/$PERIOD $fileName	
		textFileName="mergedRuns_$fileName.txt"
		mv MergedRuns.txt $BASEDIR/$SYSTEM/$PERIOD/$textFileName
		echo "moved MergedRuns.txt to $BASEDIR/$SYSTEM/$PERIOD/$textFileName";
	    elif [ "$k" -eq 1 ]; then
		    echo "Only one file, nothing to merge"
		    singleFile=`cat FilesToMerge.txt`
		    mv $singleFile $BASEDIR/$SYSTEM/$PERIOD/$fileName
		    echo "moved $singleFile "
		    echo "to    $BASEDIR/$SYSTEM/$PERIOD/$fileName"
		    textFileName="mergedRuns_$fileName.txt"
		    mv MergedRuns.txt $BASEDIR/$SYSTEM/$PERIOD/$textFileName
		    echo "moved MergedRuns.txt to $BASEDIR/$SYSTEM/$PERIOD/$textFileName";
	    else
		echo "nothing to merge"
	    fi
	    if [ -f $BASEDIR/$SYSTEM/$PERIOD/$fileName ]; then
		echo "OK; merged file exists"
	    else
		echo "ERROR: no merged file"
	    fi
	fi
	
    done # end filenames loop
    
    echo "**********************************************************************************************************"
    echo "**********************************************************************************************************"
    
    # clean-up if necessary
    if [ -f FilesToMerge.txt ]; then rm FilesToMerge.txt; fi
    if [ -f UnmergedFilesToMerge.txt ]; then rm UnmergedFilesToMerge.txt; fi
    if [ -f MergedRuns.txt ]; then rm MergedRuns.txt; fi
}



#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "337_20171030-1322" "data" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "338_20171030-1331" "data" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "339_20171030-1322" "data" "pass3_lowIR_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1" "646_20171030-1323" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1" "647_20171030-1324" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1" "648_20171030-1324" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1a" "649_20171030-1326" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1a" "650_20171030-1326" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1a" "651_20171031-0950" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1b" "652_20171030-1328" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1b" "653_20171030-1329" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1b" "654_20171030-1605" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1c" "655_20171030-1329" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1c" "656_20171030-1329" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1c" "657_20171030-1532" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "666_20171106-1018" "sim" "pass1" # train configs 246 & 247
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "669_20171106-1019" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "668_20171106-1019" "sim" "pass3_lowIR_pidfix"

DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/352_20171130-1612/merge_runlist_4/GammaConvV1_254.root" "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_254_list1_train352.root" "yes"
DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/353_20171130-1607/merge/GammaConvV1_254.root" "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_254_list2_train353.root" "yes"
DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/354_20171130-1608/merge/GammaConvV1_254.root" "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_254_list3_train354.root" "yes"
