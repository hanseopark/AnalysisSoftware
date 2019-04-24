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
    fn=`echo $1 | rev | cut -d "/" -f 1 | rev`
    part2=`echo $fn | cut -d "_" -f 2`
    if [[ "$part2" == *".root" ]];
    then
    	config=`echo $part2 | cut -d "." -f 1`
    else
    	config=$part2
    fi
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
	    echo "Check for missing tracks and V0s..."
	    root -l -b -q ../TaskQA/CheckForMissingTracksAndV0sPeriodwise.C\(\"$2\"\)
	fi
    else
	echo "Warning: download failed. Try once more."
	#alien_cp alien:$1 file:$2
	if [ -f "$2" ]; then
	    if [ "$3" == "yes" ];then
		ChangeStructure $2
		echo "Check for missing tracks and V0s..."
		root -l -b -q ../TaskQA/CheckForMissingTracksAndV0sPeriodwise.C\(\"$2\"\)
	    fi
	else
	    echo "ERROR: download failed. Giving up."
	fi
    fi
}



DownloadRunwise(){

    # arguments:
    # $1:  "GammaConvN" "GammaConvA" or "PhotonQA"
    # $2:  "ESD" or "AOD"
    # $3:  "LHC15o"
    # $4:  "337_20171030-1322"
    # $5:  "data" or "sim"
    # $6:  "pass1" / "pass1_pidfix" / "pass3_lowIR_pidfix" / ""
    
    # merge all runs in the end?
    declare -i doMerge=0      # merge runs: 1 for yes, 0 for no
    declare -i doMergeRun=1   # merge sub-run files: 1 for yes, 0 for no

    DIR=/home/meike/analysis/data/GridOutput
    TASK=$1   
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
	    AODNO="AOD194" 
	    TRAINDIR="GA_PbPb_AOD"
	else
	    echo "DATATYPE $DATATYPE not defined"
	fi
    elif [ "$TYPE" == "sim" ];then
	#YEAR="2016"
	YEAR="2018"
	PREFIX=""
	PASS=""
	if [ "$DATATYPE" == "ESD" ];then
	    AODNO=""
	    TRAINDIR="GA_PbPb_MC"
	elif [ "$DATATYPE" == "AOD" ];then
	    #if [[ "$PERIOD" == *"extra" ]]; then	   
	    #AODNO=""   
	    #else
	    AODNO="AOD198"   
	    #fi
	    TRAINDIR="GA_PbPb_MC_AOD"
	else
	    echo "DATATYPE $DATATYPE not defined"
	fi
    else
	echo "TYPE $TYPE not defined"
    fi

    # settings according to which task was run and which trainconfigs were used
    if [ "$TASK" == "GammaConvC" ];then
	FILENAMES=("GammaConvV1_288.root" "GammaConvV1_289.root" "GammaConvV1_290.root" "GammaConvV1_291.root" "GammaConvV1_292.root" "GammaConvV1_293.root");  
	declare -i doChStr=1                
	BASEDIR=$DIR/GammaConv/
    elif [ "$TASK" == "GammaConvN" ];then
	FILENAMES=("GammaConvV1_568.root" "GammaConvV1_569.root");
	declare -i doChStr=1               
	BASEDIR=$DIR/GammaConv/
    elif [ "$TASK" == "GammaConvA" ];then
	FILENAMES=("GammaConvV1_570.root" "GammaConvV1_571.root");
	declare -i doChStr=1
	BASEDIR=$DIR/GammaConv/
    elif [ "$TASK" == "PhotonQA" ];then
	FILENAMES=("AnalysisResults.root");
	declare -i doChStr=0
	BASEDIR=$DIR/PhotonQA/
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
		    echo "WARNING: download failed. Try once more."
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
			echo "ERROR: download failed. Giving up."
		    fi
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
			# mkdir -p $BASEDIR/$SYSTEM/$PERIOD/$run/Stage_$i and also change below
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
					echo "      WARNING: copying failed. Try once more."
					alien_cp alien:$SOURCEDIR/$fileName file:$OUTPUTDIR/
					if [ -f $OUTPUTDIR/$fileName ]; then  # only if file exists:
					    echo "      copied               $SOURCEDIR/$fileName "
					    echo "      to                   $OUTPUTDIR"
					    echo -n " $OUTPUTDIR/$fileName" >> UnmergedFilesToMerge.txt
					    l=l+1
					else
					    echo "      ERROR: copying failed. Giving up."
					fi
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

#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/352_20171130-1612/merge_runlist_4/GammaConvV1_254.root" "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_254_list1_train352.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/353_20171130-1607/merge/GammaConvV1_254.root" "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_254_list2_train353.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/354_20171130-1608/merge/GammaConvV1_254.root" "/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_254_list3_train354.root" "yes"

#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "365_20180116-1257" "data" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "366_20180116-1257" "data" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "367_20180116-1257" "data" "pass3_lowIR_pidfix"

#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "365_20180116-1257" "data" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "366_20180116-1257" "data" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC15o" "367_20180116-1257" "data" "pass3_lowIR_pidfix"

#DownloadRunwise "PhotonQA" "ESD" "LHC16g1" "735_20180117-1125" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1" "736_20180117-1126" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1" "737_20180117-1127" "sim" "pass3_lowIR_pidfix"

#DownloadRunwise "PhotonQA" "ESD" "LHC16g1a" "738_20180117-1127" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1a" "739_20180117-1127" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1a" "740_20180117-1127" "sim" "pass3_lowIR_pidfix"

#DownloadRunwise "PhotonQA" "ESD" "LHC16g1b" "742_20180117-1126" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1b" "741_20180117-1127" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1b" "743_20180117-1126" "sim" "pass3_lowIR_pidfix"

#DownloadRunwise "PhotonQA" "ESD" "LHC16g1c" "744_20180117-1126" "sim" "pass1"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1c" "745_20180117-1126" "sim" "pass1_pidfix"
#DownloadRunwise "PhotonQA" "ESD" "LHC16g1c" "746_20180117-1126" "sim" "pass3_lowIR_pidfix"

#DownloadRunwise "GammaConv" "ESD" "LHC16g1" "726_20171215-1812" "sim" "pass1"                   
#DownloadRunwise "GammaConv" "ESD" "LHC16g1" "727_20171215-1811" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16g1" "728_20171215-1812" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16g1a" "717_20171215-1810" "sim" "pass1" 
#DownloadRunwise "GammaConv" "ESD" "LHC16g1a" "718_20171215-1810" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16g1a" "719_20171215-1810" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16g1b" "721_20171215-1810" "sim" "pass1" 
#DownloadRunwise "GammaConv" "ESD" "LHC16g1b" "720_20171215-1810" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16g1b" "722_20171215-1811" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16g1c" "723_20171215-1811" "sim" "pass1" 
#DownloadRunwise "GammaConv" "ESD" "LHC16g1c" "724_20171215-1811" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16g1c" "725_20171215-1811" "sim" "pass3_lowIR_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "729_20171215-1812" "sim" "pass1" 
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "730_20171215-1812" "sim" "pass1_pidfix"         
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "731_20171215-1812" "sim" "pass3_lowIR_pidfix"  
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "732_20171215-1813" "sim" "pass1"                 
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "733_20171215-1812" "sim" "pass1_pidfix"         
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "734_20171215-1813" "sim" "pass3_lowIR_pidfix"   

#for trainConfig in 44 45 46 47; do  
    #filename=GammaConvV1_${trainConfig}_list1_train356.root
    #if [ -f MergedFilesToMerge.txt ]; then rm MergedFilesToMerge.txt; fi 
    #for folder in 001 002 003 005 006 007 008 009 010 011 012; do
	#mkdir -p $path/Stage_1_train356/$folder
	#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/356_20171209-1124/merge_runlist_4/Stage_1/$folder/GammaConvV1_$trainConfig.root" "$path/Stage_1_train356/$folder/$filename" "yes"
	#echo -n " $path/Stage_1_train356/$folder/$filename" >> MergedFilesToMerge.txt
    #done
    #MergeList "MergedFilesToMerge.txt" $path $filename
    #ChangeStructure $path/$filename
#    if [ -f MergedFilesToMerge.txt ]; then rm MergedFilesToMerge.txt; fi 
    #filename=GammaConvV1_${trainConfig}_list2_train359.root
    #for folder in 001 002 003 004 005 006 007; do
	#mkdir -p $path/Stage_1_train359/$folder
	#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/359_20171209-1128/merge/Stage_1/$folder/GammaConvV1_$trainConfig.root" "$path/Stage_1_train359/$folder/$filename" "yes"
#	echo -n " $path/Stage_1_train359/$folder/$filename" >> MergedFilesToMerge.txt
    #done
#    MergeList MergedFilesToMerge.txt $path $filename
#    ChangeStructure $path/$filename
#    if [ -f MergedFilesToMerge.txt ]; then rm MergedFilesToMerge.txt; fi 
#    # DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/362_20171209-1131/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train362.root" "yes"
#done

#for trainConfig in 48 49 50; do
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/357_20171209-1125/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list1_train357.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/360_20171209-1129/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list2_train360.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/363_20171209-1131/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train363.root" "yes"
#done

#LHC16g1  
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1/"
#for trainConfig in 266 268 270 272; do
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/726_20171215-1812/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list1_train726.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/727_20171215-1811/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list2_train727.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/728_20171215-1812/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train728.root" "yes"
#done
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/704_20171210-0003/merge_runlist_1/AnalysisResults.root" "$path/AnalysisResults_list1_train704.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/705_20171210-0004/merge/AnalysisResults.root" "$path/AnalysisResults_list2_train705.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/706_20171210-0004/merge/AnalysisResults.root" "$path/AnalysisResults_list3_train706.root" "no"

#LHC16g1a
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1a/"
#for trainConfig in 266 268 270 272; do
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/717_20171215-1810/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list1_train717.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/718_20171215-1810/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list2_train718.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/719_20171215-1810/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train719.root" "yes"
#done
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/695_20171209-1211/merge_runlist_1/AnalysisResults.root" "$path/AnalysisResults_list1_train695.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/696_20171210-0117/merge/AnalysisResults.root" "$path/AnalysisResults_list2_train696.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/697_20171210-1308/merge/AnalysisResults.root" "$path/AnalysisResults_list3_train697.root" "no"

#LHC16g1b
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1b/"
#for trainConfig in 266 268 270 272; do
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/721_20171215-1810/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list1_train721.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/720_20171215-1810/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list2_train720.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/722_20171215-1811/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train722.root" "yes"
#done
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/699_20171209-1149/merge_runlist_1/AnalysisResults.root" "$path/AnalysisResults_list1_train699.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/698_20171209-1148/merge/AnalysisResults.root" "$path/AnalysisResults_list2_train698.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/700_20171209-1150/merge/AnalysisResults.root" "$path/AnalysisResults_list3_train700.root" "no"

#LHC16g1c 
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1c/"
#for trainConfig in 266 268 270 272; do
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/723_20171215-1811/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list1_train723.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/724_20171215-1811/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list2_train724.root" "yes"
#     DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/725_20171215-1811/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train725.root" "yes"
#done
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/701_20171209-1150/merge_runlist_1/AnalysisResults.root" "$path/AnalysisResults_list1_train701.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/702_20171210-0000/merge/AnalysisResults.root" "$path/AnalysisResults_list2_train702.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/703_20171209-1152/merge/AnalysisResults.root" "$path/AnalysisResults_list3_train703.root" "no"

#LHC16h4 
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16h4/"
#for trainConfig in 266 268 270 272; do
     #DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/729_20171215-1812/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list1_train729.root" "yes"
     #merging#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/730_20171215-1812/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list2_train730.root" "yes"
     #DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/731_20171215-1812/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train731.root" "yes"
#done

#for trainConfig in 267 269 271 273; do
    #DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/732_20171215-1813/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list1_train732.root" "yes"
    #merging#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/733_20171215-1812/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list2_train733.root" "yes"
    #DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/734_20171215-1813/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_list3_train734.root" "yes"
#done

#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/707_20171210-1309/merge_runlist_1/AnalysisResults.root" "$path/AnalysisResults_list1_train707.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/708_20171210-1300/merge/AnalysisResults.root" "$path/AnalysisResults_list2_train708.root" "no"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/709_20171210-1259/merge/AnalysisResults.root" "$path/AnalysisResults_list3_train709.root" "no"

#LHC16i1a 
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i1a/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/792_20180125-1532/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train792.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/792_20180125-1532/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train792.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/792_20180125-1532/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train792.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/792_20180125-1532/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train792.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/792_20180125-1532/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train792.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/792_20180125-1532/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train792.root" "yes"

#LHC16i1b
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i1b/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/793_20180125-1532/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train793.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/793_20180125-1532/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train793.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/793_20180125-1532/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train793.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/793_20180125-1532/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train793.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/793_20180125-1532/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train793.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/793_20180125-1532/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train793.root" "yes"

#LHC16i1c: 
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i1c/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/794_20180125-1532/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train794.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/794_20180125-1532/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train794.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/794_20180125-1532/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train794.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/794_20180125-1532/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train794.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/794_20180125-1532/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train794.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/794_20180125-1532/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train794.root" "yes"

#LHC16i2a  
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i2a/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/795_20180125-1532/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train795.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/795_20180125-1532/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train795.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/795_20180125-1532/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train795.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/795_20180125-1532/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train795.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/795_20180125-1532/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train795.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/795_20180125-1532/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train795.root" "yes"

#LHC16i2b 
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i2b/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/796_20180125-1532/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train796.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/796_20180125-1532/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train796.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/796_20180125-1532/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train796.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/796_20180125-1532/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train796.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/796_20180125-1532/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train796.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/796_20180125-1532/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train796.root" "yes"

#LHC16i2c  
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i2c/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/797_20180125-1532/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train797.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/797_20180125-1532/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train797.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/797_20180125-1532/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train797.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/797_20180125-1532/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train797.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/797_20180125-1532/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train797.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/797_20180125-1532/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train797.root" "yes"

#LHC16i3a  
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i3a/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/798_20180125-1532/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train798.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/798_20180125-1532/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train798.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/798_20180125-1532/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train798.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/798_20180125-1532/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train798.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/798_20180125-1532/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train798.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/798_20180125-1532/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train798.root" "yes"

#LHC16i3b 
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i3b/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/799_20180129-1437/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train799.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/799_20180129-1437/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train799.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/799_20180129-1437/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train799.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/799_20180129-1437/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train799.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/799_20180129-1437/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train799.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/799_20180129-1437/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train799.root" "yes"

#LHC16i3c 
#path="/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i3c/"
#mkdir -p $path
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/800_20180125-1533/merge_runlist_1/GammaConvV1_270.root" "$path/GammaConvV1_270_list1_train800.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/800_20180125-1533/merge_runlist_1/GammaConvV1_272.root" "$path/GammaConvV1_272_list1_train800.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/800_20180125-1533/merge_runlist_3/GammaConvV1_270.root" "$path/GammaConvV1_270_list2_train800.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/800_20180125-1533/merge_runlist_3/GammaConvV1_272.root" "$path/GammaConvV1_272_list2_train800.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/800_20180125-1533/merge_runlist_4/GammaConvV1_270.root" "$path/GammaConvV1_270_list3_train800.root" "yes"
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/800_20180125-1533/merge_runlist_4/GammaConvV1_272.root" "$path/GammaConvV1_272_list3_train800.root" "yes"




# =================================== FEB 2018 ===========================================================
# ==================================  Data AOD  ==========================================================

#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#for trainConfig in 266 268 274 276; do
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/442_20180217-2113/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_AOD194_list1_train442.root" "yes"  
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/443_20180217-2113/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_AOD194_list2_train443.root" "yes"  
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/444_20180217-2113/merge/GammaConvV1_$trainConfig.root" "$path/GammaConvV1_${trainConfig}_AOD194_list3_train444.root" "yes"  
#done
#DownloadRunwise "GammaConv" "AOD" "LHC15o" "442_20180217-2113" "data" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC15o" "443_20180217-2113" "data" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC15o" "444_20180217-2113" "data" "pass3_lowIR_pidfix"


#for trainConfig in 51 52 53 54 55; do
#    for folder in 001 002 003 004 005 006; do
#	path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/dca/Stage1/$folder
#	mkdir -p $path
#	filename=GammaConvV1_${trainConfig}_AOD194_list1_train445.root
#	DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/445_20180218-1024/merge/Stage_1/$folder/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#   done
#done

#for trainConfig in 51 52 53 54 55; do
#    for folder in 001 002 003 004 005; do
#	path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/dca/Stage1/$folder
#	mkdir -p $path
#	filename=GammaConvV1_${trainConfig}_AOD194_list2_train446.root
#	DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/446_20180218-1024/merge/Stage_1/$folder/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    done
#done

#for trainConfig in 51 52 53 54 55; do
#	path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/dca/
#	mkdir -p $path
#	filename=GammaConvV1_${trainConfig}_AOD194_list3_train447.root
#	DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/447_20180218-1024/merge/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    done
#done

# errors
#/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/dca/Stage1/005/GammaConvV1_51_AOD194_list1_train445.root



# ==================================  Data ESD  ==========================================================

#for trainConfig in 51 52 53 54 55; do
#    #list1:
#    for folder in 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016; do
#	path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/dca/Stage1/$folder
#	mkdir -p $path
#	filename=GammaConvV1_${trainConfig}_list1_train389.root
#	DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/389_20180219-0050/merge_runlist_4/Stage_1/$folder/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#   done
#    #list2: 
#    for folder in 001 002 003 004 005 006 007 008 009 010; do
#	path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/dca/Stage1/$folder
#	mkdir -p $path
#	filename=GammaConvV1_${trainConfig}_list2_train390.root
#	DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/390_20180218-1819/merge/Stage_1/$folder/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    done
#    #list3:
#    path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/dca/    
#    mkdir -p $path
#    filename=GammaConvV1_${trainConfig}_list3_train391.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb/391_20180218-1819/merge/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done



# ================================== MARCH 2018 ==========================================================
# ==================================   MC AOD   ==========================================================

# LHC16g1 - train 885- download done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train885.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/885_20180309-1206/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train885.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/885_20180309-1206/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train885.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/885_20180309-1206/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16g1" "885_20180309-1206" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1" "885_20180309-1206" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1" "885_20180309-1206" "sim" "pass3_lowIR_pidfix"

# LHC16g1a - train 903
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/    
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train903.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/903_20180326-1023/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train903.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/903_20180326-1023/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train903.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/903_20180326-1023/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16g1a" "903_20180326-1023" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1a" "903_20180326-1023" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1a" "903_20180326-1023" "sim" "pass3_lowIR_pidfix"


# LHC16g1b - train 887 - download done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train887.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/887_20180309-1232/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train887.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/887_20180309-1232/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train887.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/887_20180309-1232/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b" "887_20180309-1232" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b" "887_20180309-1232" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b" "887_20180309-1232" "sim" "pass3_lowIR_pidfix"

##LHC16g1b_extra - train 888
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD_list1_train888.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/888_20180309-1207/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list2_train888.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/888_20180309-1207/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list3_train888.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/888_20180309-1207/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b_extra" "888_20180309-1207" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b_extra" "888_20180309-1207" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b_extra" "888_20180309-1207" "sim" "pass3_lowIR_pidfix"

##LHC16g1c - train 889
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train889.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/889_20180309-1208/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train889.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/889_20180309-1208/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train889.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/889_20180309-1208/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c" "889_20180309-1208" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c" "889_20180309-1208" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c" "889_20180309-1208" "sim" "pass3_lowIR_pidfix"

##LHC16g1c_extra - train 890
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD_list1_train890.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/890_20180309-1208/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list2_train890.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/890_20180309-1208/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list3_train890.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/890_20180309-1208/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c_extra" "890_20180309-1208" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c_extra" "890_20180309-1208" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c_extra" "890_20180309-1208" "sim" "pass3_lowIR_pidfix"

# LHC16h4 - train 891
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train891.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/891_20180309-1233/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train891.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/891_20180309-1233/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train891.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/891_20180309-1233/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16h4" "891_20180309-1233" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16h4" "891_20180309-1233" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16h4" "891_20180309-1233" "sim" "pass3_lowIR_pidfix"


# LHC16i1a - train 892 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train892.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/892_20180309-1233/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train892.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/892_20180309-1233/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train892.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/892_20180309-1233/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i1a" "892_20180309-1233" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1a" "892_20180309-1233" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1a" "892_20180309-1233" "sim" "pass3_lowIR_pidfix"


# LHC16i1b - train 893  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train893.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/893_20180309-1233/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train893.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/893_20180309-1233/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train893.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/893_20180309-1233/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i1b" "893_20180309-1233" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1b" "893_20180309-1233" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1b" "893_20180309-1233" "sim" "pass3_lowIR_pidfix"


# LHC16i1c - train 894  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train894.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/894_20180309-1233/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train894.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/894_20180309-1233/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train894.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/894_20180309-1233/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i1c" "894_20180309-1233" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1c" "894_20180309-1233" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1c" "894_20180309-1233" "sim" "pass3_lowIR_pidfix"


# LHC16i2a - train 895  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train895.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/895_20180309-1233/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train895.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/895_20180309-1233/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train895.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/895_20180309-1233/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i2a" "895_20180309-1233" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2a" "895_20180309-1233" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2a" "895_20180309-1233" "sim" "pass3_lowIR_pidfix"



# LHC16i2b - train 896  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train896.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/896_20180309-1233/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train896.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/896_20180309-1233/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train896.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/896_20180309-1233/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i2b" "896_20180309-1233" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2b" "896_20180309-1233" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2b" "896_20180309-1233" "sim" "pass3_lowIR_pidfix"

# LHC16i2c - train 897  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train897.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/897_20180309-1233/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train897.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/897_20180309-1233/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train897.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/897_20180309-1233/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i2c" "897_20180309-1233" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2c" "897_20180309-1233" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2c" "897_20180309-1233" "sim" "pass3_lowIR_pidfix"


# LHC16i3a - train 898  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train898.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/898_20180309-1234/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train898.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/898_20180309-1234/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train898.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/898_20180309-1234/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i3a" "898_20180309-1234" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3a" "898_20180309-1234" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3a" "898_20180309-1234" "sim" "pass3_lowIR_pidfix"

# LHC16i3b - train 899  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train899.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/899_20180309-1234/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train899.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/899_20180309-1234/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train899.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/899_20180309-1234/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i3b" "899_20180309-1234" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3b" "899_20180309-1234" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3b" "899_20180309-1234" "sim" "pass3_lowIR_pidfix"

# LHC16i3c - train 900 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train900.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/900_20180309-1234/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train900.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/900_20180309-1234/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train900.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/900_20180309-1234/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "AOD" "LHC16i3c" "900_20180309-1234" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3c" "900_20180309-1234" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3c" "900_20180309-1234" "sim" "pass3_lowIR_pidfix"


# LHC16g1 - dca train 902 - not runwise 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/dca/   
#for trainConfig in 51 52 53 54 55; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train902.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/902_20180311-1323/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train902.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/902_20180311-1323/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train902.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/902_20180311-1323/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done


# ==================================   MC ESD   ==========================================================

# LHC16g1 - train 859 859_20180312-0841 was killed -> 899_20180317-1830 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train899.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/899_20180317-1830/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train899.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/899_20180317-1830/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train899.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/899_20180317-1830/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

# LHC16g1a - train 854 854_20180312-0833 was killed -> 898_20180317-1740 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1a/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train898.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/898_20180317-1740/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train898.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/898_20180317-1740/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train898.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/898_20180317-1740/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done


# LHC16g1b & LHC16g1b_extra - train 884 884_20180312-0833 was killed -> 904_20180317-1743 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train904.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/904_20180317-1743_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train904.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/904_20180317-1743_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train904.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/904_20180317-1743_child_1/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    #LHC16g1b_extra
#    filename=GammaConvV1_${trainConfig}_list1_train904_extra.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/904_20180317-1743_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train904_extra.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/904_20180317-1743_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train904_extra.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/904_20180317-1743_child_2/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
# runwise: 883 ! files 266 and 268 already exist and will be overwritten ! need two files, one from LHC16g1b and one from LHC16g1b_extra


# LHC16g1c & LHC16g1c_extra - train 883 883_20180312-0832 was killed -> 905_20180317-1745
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train905.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/905_20180317-1745_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train905.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/905_20180317-1745_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train905.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/905_20180317-1745_child_1/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    #LHC16g1b_extra
#    filename=GammaConvV1_${trainConfig}_list1_train905_extra.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/905_20180317-1745_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train905_extra.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/905_20180317-1745_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train905_extra.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/905_20180317-1745_child_2/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
# runwise: 884   ! files 266 and 268 already exist and will be overwritten !  need two files, one from LHC16g1b and one from LHC16g1b_extra


# LHC16h4 - train 860 860_20180312-0846 was killed  -> 860_20180312-0846 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16h4/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train860.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/860_20180312-0846/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train860.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/860_20180312-0846/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train860.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/860_20180312-0846/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "860_20180312-0846" "sim" "pass1"                  # ! files 266 and 268 already exist and will be overwritten !
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "860_20180312-0846" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16h4" "860_20180312-0846" "sim" "pass3_lowIR_pidfix"

# LHC16i1a - train 861  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i1a/   "pass3_lowIR_pidfix"
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train861.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/861_20180312-0836/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train861.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/861_20180312-0836/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train861.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/861_20180312-0836/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16i1a" "861_20180312-0836" "sim" "pass1"
#DownloadRunwise "GammaConv" "ESD" "LHC16i1a" "861_20180312-0836" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16i1a" "861_20180312-0836" "sim" "pass3_lowIR_pidfix"

# LHC16i1b - train 862 862_20180312-0834 was killed  -> 902_20180317-1742
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i1b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train902.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/902_20180317-1742/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train902.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/902_20180317-1742/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train902.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/902_20180317-1742/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

# LHC16i1c - train 863 863_20180312-0832 was killed -> 903_20180317-1743 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i1c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train903.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/903_20180317-1743/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train903.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/903_20180317-1743/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train903.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/903_20180317-1743/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done


# LHC16i2a - train 864 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i2a/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train864.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/864_20180309-1202/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train864.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/864_20180309-1202/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train864.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/864_20180309-1202/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16i2a" "864_20180309-1202" "sim" "pass1"
#DownloadRunwise "GammaConv" "ESD" "LHC16i2a" "864_20180309-1202" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16i2a" "864_20180309-1202" "sim" "pass3_lowIR_pidfix"

# LHC16i2b - train 865 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i2b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train865.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/865_20180309-1203/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train865.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/865_20180309-1203/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train865.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/865_20180309-1203/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16i2b" "865_20180309-1203" "sim" "pass1"
#DownloadRunwise "GammaConv" "ESD" "LHC16i2b" "865_20180309-1203" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16i2b" "865_20180309-1203" "sim" "pass3_lowIR_pidfix"

# LHC16i2c - train 866  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i2c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train866.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/866_20180309-1204/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train866.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/866_20180309-1204/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train866.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/866_20180309-1204/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16i2c" "866_20180309-1204" "sim" "pass1"
#DownloadRunwise "GammaConv" "ESD" "LHC16i2c" "866_20180309-1204" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16i2c" "866_20180309-1204" "sim" "pass3_lowIR_pidfix"

# LHC16i3a - train 867 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i3a/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train867.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/867_20180309-1204/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train867.root
#   DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/867_20180309-1204/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train867.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/867_20180309-1204/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16i3a" "867_20180309-1204" "sim" "pass1"
#DownloadRunwise "GammaConv" "ESD" "LHC16i3a" "867_20180309-1204" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16i3a" "867_20180309-1204" "sim" "pass3_lowIR_pidfix"

# LHC16i3b - train 868  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i3b/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train868.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/868_20180309-1205/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train868.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/868_20180309-1205/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train868.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/868_20180309-1205/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16i3b" "868_20180309-1205" "sim" "pass1"
#DownloadRunwise "GammaConv" "ESD" "LHC16i3b" "868_20180309-1205" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16i3b" "868_20180309-1205" "sim" "pass3_lowIR_pidfix"

# LHC16i3c - train  869  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16i3c/   
#for trainConfig in 266 268 274 276; do
#    filename=GammaConvV1_${trainConfig}_list1_train869.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/869_20180309-1205/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train869.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/869_20180309-1205/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train869.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/869_20180309-1205/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done
#DownloadRunwise "GammaConv" "ESD" "LHC16i3c" "869_20180309-1205" "sim" "pass1"
#DownloadRunwise "GammaConv" "ESD" "LHC16i3c" "869_20180309-1205" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "ESD" "LHC16i3c" "869_20180309-1205" "sim" "pass3_lowIR_pidfix"



# LHC16g1 - dca train 870  - not runwise  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC16g1/dca/   
#for trainConfig in 51 52 53 54 55; do
#    filename=GammaConvV1_${trainConfig}_list1_train870.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/870_20180309-1908/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list2_train870.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/870_20180309-1908/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_list3_train870.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC/870_20180309-1908/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done


# ================================== JUNE 2018 ==========================================================
# ==================================   MC AOD   ==========================================================

# LHC16g1 - train 921
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train921.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/921_20180529-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train921.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/921_20180529-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train921.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/921_20180529-1513/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

# LHC16g1a - train 922 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/    
#for trainConfig in 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train922.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/922_20180529-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train922.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/922_20180529-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train922.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/922_20180529-1513/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

# LHC16g1b - train 923 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/   
#for trainConfig in 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train923.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/923_20180529-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train923.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/923_20180529-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train923.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/923_20180529-1513/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

##LHC16g1b_extra - train 924 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/   
#for trainConfig in 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD_list1_train924.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/924_20180529-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list2_train924.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/924_20180529-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list3_train924.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/924_20180529-1513/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

##LHC16g1c - train 925
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/   
#for trainConfig in 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD198_list1_train925.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/925_20180529-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list2_train925.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/925_20180529-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD198_list3_train925.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/925_20180529-1513/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done


##LHC16g1c_extra - train 926 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/   
#for trainConfig in 274 276; do
#    filename=GammaConvV1_${trainConfig}_AOD_list1_train926.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/926_20180529-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list2_train926.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/926_20180529-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD_list3_train926.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/926_20180529-1513/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

# ==================================   Data AOD   ==========================================================

## LHC15o - train 452
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/   
#for trainConfig in 266 268 270 272 274 276 278 280; do
#    filename=GammaConvV1_${trainConfig}_AOD194_list1_train452.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/452_20180530-1122_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD194_list2_train452.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/452_20180530-1122_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#    filename=GammaConvV1_${trainConfig}_AOD194_list3_train452.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/452_20180530-1122_child_3/merge/GammaConvV1_$trainConfig.root" "$path/$filename" "yes"
#done

## LHC15o - train 458
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#filename1=GammaConvV1_282_AOD194_list1_train458.root
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/458_20180619-1034_child_1/merge_runlist_1/GammaConvV1_282.root" "$path/$filename1" "yes"
#filename2=GammaConvV1_282_AOD194_list2_train458.root
#DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/458_20180619-1034_child_2/merge/GammaConvV1_282.root" "$path/$filename2" "yes"
#hadd $path/GammaConvV1_282_AOD194_highIR_train458.root $path/$filename1 $path/$filename2


# ================================== JULY 2018 ==========================================================
# ==================================   Data AOD   ==========================================================

## LHC15o - train 461
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#for trainConfig in 283 284 285 286; do
#    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train461.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/461_20180629-1021_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train461.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/461_20180629-1021_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD194_highIR_train461.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD194_list3_train461.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/461_20180629-1021_child_3/merge/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC15o" "461_20180629-1021_child_1" "data" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC15o" "461_20180629-1021_child_2" "data" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC15o" "461_20180629-1021_child_3" "data" "pass3_lowIR_pidfix"


# ==================================   MC AOD   ==========================================================

## LHC16g1 - train 942
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train942.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/942_20180629-1521/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train942.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/942_20180629-1521/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train942.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train942.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/942_20180629-1521/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16g1a - train 943  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train943.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/943_20180629-1521/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train943.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/943_20180629-1521/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train943.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train943.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/943_20180629-1521/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16g1b - train 944   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train944.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/944_20180629-1521/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train944.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/944_20180629-1521/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train944.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train944.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/944_20180629-1521/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16g1b_extra - train 945   ++++ renamed files later from 'AOD198' to 'AOD' ++++++
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train945.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/945_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train945.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/945_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train945.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train945.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/945_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16g1c - train 946   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train946.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/946_20180629-1520/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train946.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/946_20180629-1520/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train946.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train946.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/946_20180629-1520/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16g1c_extra - train 947      ++++ renamed files later from 'AOD198' to 'AOD' ++++++
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train947.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/947_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train947.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/947_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train947.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train947.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/947_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16h4 - train 948  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 283 284 285 286; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train948.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/948_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train948.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/948_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train948.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train948.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/948_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i1a - train 949  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train949.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/949_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train949.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/949_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train949.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train949.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/949_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i1b - train 950  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train950.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/950_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train950.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/950_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train950.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train950.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/950_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i1c - train 951  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train951.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/951_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train951.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/951_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train951.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train951.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/951_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i2a - train 952 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train952.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/952_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train952.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/952_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train952.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train952.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/952_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i2b - train 953  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train953.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/953_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train953.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/953_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train953.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train953.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/953_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i2c - train 954  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train954.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/954_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train954.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/954_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train954.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train954.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/954_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i3a - train 955  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train955.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/955_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train955.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/955_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train955.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train955.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/955_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i3b - train 956  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train956.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/956_20180629-1523/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train956.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/956_20180629-1523/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train956.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train956.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/956_20180629-1523/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

## LHC16i3c - train 957  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train957.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/957_20180629-1522/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train957.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/957_20180629-1522/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train957.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train957.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/957_20180629-1522/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# ==================================   MC AOD   ==========================================================

## LHC16g1 - train 958_20180706-1107
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train958.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/958_20180706-1107/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train958.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/958_20180706-1107/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train958.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train958.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/958_20180706-1107/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16g1" "958_20180706-1107" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1" "958_20180706-1107" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1" "958_20180706-1107" "sim" "pass3_lowIR_pidfix"

# LHC16g1a - train 959_20180706-1105   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train959.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/959_20180706-1105/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train959.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/959_20180706-1105/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train959.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train959.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/959_20180706-1105/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16g1a" "959_20180706-1105" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1a" "959_20180706-1105" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1a" "959_20180706-1105" "sim" "pass3_lowIR_pidfix"

## LHC16g1b - train 960_20180706-1105   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train960.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/960_20180706-1105/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train960.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/960_20180706-1105/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train960.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train960.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/960_20180706-1105/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16g1b" "960_20180706-1105" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b" "960_20180706-1105" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b" "960_20180706-1105" "sim" "pass3_lowIR_pidfix"

## LHC16g1b_extra - train 961_20180706-1105  ++++ renamed files later from 'AOD198' to 'AOD' ++++++    
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train961.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/961_20180706-1105/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train961.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/961_20180706-1105/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train961.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train961.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/961_20180706-1105/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16g1b_extra" "961_20180706-1105" "sim" "pass1"
#ownloadRunwise "GammaConv" "AOD" "LHC16g1b_extra" "961_20180706-1105" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1b_extra" "961_20180706-1105" "sim" "pass3_lowIR_pidfix"

## LHC16g1c - train 962_20180706-1105   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train962.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/962_20180706-1105/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train962.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/962_20180706-1105/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train962.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train962.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/962_20180706-1105/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16g1c" "962_20180706-1105" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c" "962_20180706-1105" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c" "962_20180706-1105" "sim" "pass3_lowIR_pidfix"

## LHC16g1c_extra - train 963_20180706-1105   ++++ renamed files later from 'AOD198' to 'AOD' ++++++
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train963.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/963_20180706-1105/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train963.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/963_20180706-1105/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train963.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train963.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/963_20180706-1105/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16g1c_extra" "963_20180706-1105" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c_extra" "963_20180706-1105" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16g1c_extra" "963_20180706-1105" "sim" "pass3_lowIR_pidfix"

## LHC16h4 - train 964_20180706-1425 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 283 284 285 286; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train964.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/964_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train964.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/964_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train964.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train964.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/964_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16h4" "964_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16h4" "964_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16h4" "964_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i1a - train 965_20180706-1425 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train965.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/965_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train965.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/965_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train965.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train965.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/965_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i1a" "965_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1a" "965_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1a" "965_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i1b - train 966_20180706-1425  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train966.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/966_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train966.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/966_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train966.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train966.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/966_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i1b" "966_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1b" "966_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1b" "966_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i1c - train 967_20180706-1425   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train967.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/967_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train967.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/967_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train967.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train967.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/967_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i1c" "967_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1c" "967_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i1c" "967_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i2a - train 968_20180706-1425  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train968.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/968_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train968.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/968_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train968.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train968.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/968_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i2a" "968_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2a" "968_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2a" "968_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i2b - train 969_20180706-1425 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train969.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/969_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train969.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/969_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train969.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train969.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/969_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i2b" "969_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2b" "969_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2b" "969_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i2c - train 970_20180706-1425   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train970.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/970_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train970.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/970_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train970.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train970.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/970_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i2c" "970_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2c" "970_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i2c" "970_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i3a - train 971_20180706-1425 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train971.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/971_20180706-1425/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train971.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/971_20180706-1425/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train971.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train971.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/971_20180706-1425/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i3a" "971_20180706-1425" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3a" "971_20180706-1425" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3a" "971_20180706-1425" "sim" "pass3_lowIR_pidfix"

## LHC16i3b - train 972_20180706-1426   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train972.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/972_20180706-1426/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train972.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/972_20180706-1426/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train972.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train972.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/972_20180706-1426/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i3b" "972_20180706-1426" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3b" "972_20180706-1426" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3b" "972_20180706-1426" "sim" "pass3_lowIR_pidfix"

## LHC16i3c - train 973_20180706-1426  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train973.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/973_20180706-1426/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train973.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/973_20180706-1426/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD198_highIR_train973.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train973.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/973_20180706-1426/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#DownloadRunwise "GammaConv" "AOD" "LHC16i3c" "973_20180706-1426" "sim" "pass1"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3c" "973_20180706-1426" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConv" "AOD" "LHC16i3c" "973_20180706-1426" "sim" "pass3_lowIR_pidfix"


#for trainConfig in 283 285; do

    ## merge MCs with weighting
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_MB_AOD198_highIR_WMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/GammaConvV1_${trainConfig}_AOD198_list1_train958.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/GammaConvV1_${trainConfig}_AOD198_list2_train958.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/GammaConvV1_${trainConfig}_AOD198_list1_train959.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/GammaConvV1_${trainConfig}_AOD198_list2_train959.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/GammaConvV1_${trainConfig}_AOD198_list1_train960.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/GammaConvV1_${trainConfig}_AOD198_list2_train960.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/GammaConvV1_${trainConfig}_AOD_list1_train961.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/GammaConvV1_${trainConfig}_AOD_list2_train961.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/GammaConvV1_${trainConfig}_AOD198_list1_train962.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/GammaConvV1_${trainConfig}_AOD198_list2_train962.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/GammaConvV1_${trainConfig}_AOD_list1_train963.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/GammaConvV1_${trainConfig}_AOD_list2_train963.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/GammaConvV1_${trainConfig}_AOD198_list1_train964.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/GammaConvV1_${trainConfig}_AOD198_list2_train964.root   
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_LF_AOD198_highIR_WMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/GammaConvV1_${trainConfig}_AOD198_list1_train965.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/GammaConvV1_${trainConfig}_AOD198_list2_train965.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/GammaConvV1_${trainConfig}_AOD198_list1_train966.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/GammaConvV1_${trainConfig}_AOD198_list2_train966.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/GammaConvV1_${trainConfig}_AOD198_list1_train967.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/GammaConvV1_${trainConfig}_AOD198_list2_train967.root    
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_HF_AOD198_highIR_WMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/GammaConvV1_${trainConfig}_AOD198_list1_train968.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/GammaConvV1_${trainConfig}_AOD198_list2_train968.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/GammaConvV1_${trainConfig}_AOD198_list1_train969.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/GammaConvV1_${trainConfig}_AOD198_list2_train969.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/GammaConvV1_${trainConfig}_AOD198_list1_train970.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/GammaConvV1_${trainConfig}_AOD198_list2_train970.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/GammaConvV1_${trainConfig}_AOD198_list1_train971.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/GammaConvV1_${trainConfig}_AOD198_list2_train971.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/GammaConvV1_${trainConfig}_AOD198_list1_train972.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/GammaConvV1_${trainConfig}_AOD198_list2_train972.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/GammaConvV1_${trainConfig}_AOD198_list1_train973.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/GammaConvV1_${trainConfig}_AOD198_list2_train973.root
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_ALL_AOD198_highIR_WMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_MB_AOD198_highIR_WMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_LF_AOD198_highIR_WMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_HF_AOD198_highIR_WMW.root

    ## merge MCs without weighting
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_MB_AOD198_highIR_WoMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/GammaConvV1_${trainConfig}_AOD198_list1_train942.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/GammaConvV1_${trainConfig}_AOD198_list2_train942.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/GammaConvV1_${trainConfig}_AOD198_list1_train943.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/GammaConvV1_${trainConfig}_AOD198_list2_train943.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/GammaConvV1_${trainConfig}_AOD198_list1_train944.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/GammaConvV1_${trainConfig}_AOD198_list2_train944.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/GammaConvV1_${trainConfig}_AOD_list1_train945.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/GammaConvV1_${trainConfig}_AOD_list2_train945.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/GammaConvV1_${trainConfig}_AOD198_list1_train946.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/GammaConvV1_${trainConfig}_AOD198_list2_train946.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/GammaConvV1_${trainConfig}_AOD_list1_train947.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/GammaConvV1_${trainConfig}_AOD_list2_train947.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/GammaConvV1_${trainConfig}_AOD198_list1_train948.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/GammaConvV1_${trainConfig}_AOD198_list2_train948.root
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_LF_AOD198_highIR_WoMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/GammaConvV1_${trainConfig}_AOD198_list1_train949.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/GammaConvV1_${trainConfig}_AOD198_list2_train949.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/GammaConvV1_${trainConfig}_AOD198_list1_train950.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/GammaConvV1_${trainConfig}_AOD198_list2_train950.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/GammaConvV1_${trainConfig}_AOD198_list1_train951.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/GammaConvV1_${trainConfig}_AOD198_list2_train951.root
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_HF_AOD198_highIR_WoMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/GammaConvV1_${trainConfig}_AOD198_list1_train952.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/GammaConvV1_${trainConfig}_AOD198_list2_train952.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/GammaConvV1_${trainConfig}_AOD198_list1_train953.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/GammaConvV1_${trainConfig}_AOD198_list2_train953.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/GammaConvV1_${trainConfig}_AOD198_list1_train954.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/GammaConvV1_${trainConfig}_AOD198_list2_train954.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/GammaConvV1_${trainConfig}_AOD198_list1_train955.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/GammaConvV1_${trainConfig}_AOD198_list2_train955.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/GammaConvV1_${trainConfig}_AOD198_list1_train956.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/GammaConvV1_${trainConfig}_AOD198_list2_train956.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/GammaConvV1_${trainConfig}_AOD198_list1_train957.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/GammaConvV1_${trainConfig}_AOD198_list2_train957.root
    #hadd /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_ALL_AOD198_highIR_WoMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_MB_AOD198_highIR_WoMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_LF_AOD198_highIR_WoMW.root /home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/merged/GammaConvV1_${trainConfig}_MC_HF_AOD198_highIR_WoMW.root

#done

#============== pT weighting ========================
# August 2018

## LHC16h4 added particles
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 283 284 285 286; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train991.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/991_20180817-1409/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train991.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/991_20180817-1409/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train991.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/991_20180817-1409/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done



# normal particles from all MCs

# LHC16g1
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train992.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/992_20180827-1700/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train992.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/992_20180827-1700/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train992.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/992_20180827-1700/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16g1a
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train993.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/993_20180827-1700/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train993.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/993_20180827-1700/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train993.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/993_20180827-1700/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16g1b + LHC16g1b_extra
# was: two different trains
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train994.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/994_20180827-1703_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train994.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/994_20180827-1703_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train994.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/994_20180827-1703_child_1/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train994.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/994_20180827-1703_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train994.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/994_20180827-1703_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD_list3_train994.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/994_20180827-1703_child_2/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16g1c + LHC16g1c_extra  
# was: two different trains
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train995.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/995_20180827-1701_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train995.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/995_20180827-1701_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train995.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/995_20180827-1701_child_1/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train995.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/995_20180827-1701_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train995.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/995_20180827-1701_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD_list3_train995.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/995_20180827-1701_child_2/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


## LHC16h4
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train996.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/996_20180827-1701/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train996.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/996_20180827-1701/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train996.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/996_20180827-1701/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16i1a 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train997.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/997_20180827-1702/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train997.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/997_20180827-1702/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train997.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/997_20180827-1702/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i1b 998_20180827-1702 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train998.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/998_20180827-1702/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train998.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/998_20180827-1702/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train998.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/998_20180827-1702/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i1c  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train999.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/999_20180827-1702/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train999.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/999_20180827-1702/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train999.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/999_20180827-1702/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16i2a  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1000.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1000_20180827-1702/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1000.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1000_20180827-1702/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1000.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1000_20180827-1702/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16i2b 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1001.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1001_20180827-1702/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1001.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1001_20180827-1702/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1001.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1001_20180827-1702/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i2c  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1002.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1002_20180827-1702/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1002.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1002_20180827-1702/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1002.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1002_20180827-1702/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i3a  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1003.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1003_20180827-1703/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1003.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1003_20180827-1703/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1003.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1003_20180827-1703/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i3b  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1004.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1004_20180827-1703/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1004.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1004_20180827-1703/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1004.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1004_20180827-1703/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i3c  
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1005.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1005_20180827-1703/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1005.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1005_20180827-1703/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1005.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1005_20180827-1703/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


#============== pT weighting 2nd iteration ========================
# October 2018


## LHC16h4 added particles  1037_20181022-1120
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 284 286; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1037.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1037_20181022-1120/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1037.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1037_20181022-1120/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1037.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1037_20181022-1120/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# normal particles from all MCs

# LHC16g1   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1024.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1024_20181022-1003/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1024.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1024_20181022-1003/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1024.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1024_20181022-1003/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
# Runwise
#DownloadRunwise "GammaConvN" "AOD" "LHC16g1" "1024_20181022-1003" "sim" "pass1"
#DownloadRunwise "GammaConvN" "AOD" "LHC16g1" "1024_20181022-1003" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConvN" "AOD" "LHC16g1" "1024_20181022-1003" "sim" "pass3_lowIR_pidfix"

# LHC16g1a   1025_20181022-1008 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1025.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1025_20181022-1008/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1025.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1025_20181022-1008/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1025.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1025_20181022-1008/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16g1b + LHC16g1b_extra   1026_20181022-1022
# was: two different trains
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1026.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1026_20181022-1022_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1026.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1026_20181022-1022_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1026.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1026_20181022-1022_child_1/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1026.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1026_20181022-1022_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1026.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1026_20181022-1022_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD_list3_train1026.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1026_20181022-1022_child_2/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16g1c + LHC16g1c_extra    1027_20181022-1039 
# was: two different trains
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1027.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1027_20181022-1039_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1027.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1027_20181022-1039_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1027.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1027_20181022-1039_child_1/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1027.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1027_20181022-1039_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1027.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1027_20181022-1039_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD_list3_train1027.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1027_20181022-1039_child_2/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


## LHC16h4    1037_20181022-1120
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1037.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1037_20181022-1120/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1037.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1037_20181022-1120/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1037.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1037_20181022-1120/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
# Runwise
#DownloadRunwise "GammaConvN" "AOD" "LHC16h4" "1037_20181022-1120" "sim" "pass1"
#DownloadRunwise "GammaConvN" "AOD" "LHC16h4" "1037_20181022-1120" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConvN" "AOD" "LHC16h4" "1037_20181022-1120" "sim" "pass3_lowIR_pidfix"

# LHC16i1a     1028_20181022-1039 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1028.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1028_20181022-1039/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1028.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1028_20181022-1039/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1028.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1028_20181022-1039/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
# Runwise
#DownloadRunwise "GammaConvN" "AOD" "LHC16i1a" "1028_20181022-1039" "sim" "pass1"
#DownloadRunwise "GammaConvN" "AOD" "LHC16i1a" "1028_20181022-1039" "sim" "pass1_pidfix"
#DownloadRunwise "GammaConvN" "AOD" "LHC16i1a" "1028_20181022-1039" "sim" "pass3_lowIR_pidfix"

# LHC16i1b 1029_20181022-1100
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1029.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1029_20181022-1100/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1029.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1029_20181022-1100/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1029.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1029_20181022-1100/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i1c     1030_20181022-1100 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1030.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1030_20181022-1100/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1030.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1030_20181022-1100/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1030.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1030_20181022-1100/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

# LHC16i2a     1031_20181022-1005 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1031.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1031_20181022-1005/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1031.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1031_20181022-1005/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1031.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1031_20181022-1005/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i2b    1032_20181022-1006 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1032.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1032_20181022-1006/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1032.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1032_20181022-1006/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1032.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1032_20181022-1006/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i2c    1033_20181022-1006 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1033.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1033_20181022-1006/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1033.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1033_20181022-1006/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1033.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1033_20181022-1006/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i3a   1034_20181022-1006
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1034.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1034_20181022-1006/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1034.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1034_20181022-1006/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1034.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1034_20181022-1006/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i3b     1035_20181022-1006 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1035.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1035_20181022-1006/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1035.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1035_20181022-1006/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1035.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1035_20181022-1006/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done


# LHC16i3c   1036_20181022-1006   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
#for trainConfig in 283 285; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1036.root
#   DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1036_20181022-1006/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1036.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1036_20181022-1006/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1036.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1036_20181022-1006/merge_runlist_4/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done

#====================== data =============

## LHC15o - train 477_20181029-1012 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#for trainConfig in 288 289 290 291 292 293; do
#    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train477.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/477_20181029-1012_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train477.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/477_20181029-1012_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd $path/GammaConvV1_${trainConfig}_AOD194_highIR_train477.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD194_list3_train477.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/477_20181029-1012_child_3/merge/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done
 
#DownloadRunwise "GammaConvC" "AOD" "LHC15o" "477_20181029-1012_child_1" "data" "pass1"
#DownloadRunwise "GammaConvC" "AOD" "LHC15o" "477_20181029-1012_child_2" "data" "pass1_pidfix"
#DownloadRunwise "GammaConvC" "AOD" "LHC15o" "477_20181029-1012_child_3" "data" "pass3_lowIR_pidfix"


# LHC15o with new event mixing - 479_20181109-1254 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#for trainConfig in 294 295 500 501 504 505; do
#    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train479.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/479_20181109-1254_child_1/merge/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train479.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/479_20181109-1254_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_highIR_train479.root $path/$filename1 $path/$filename2
#    filename3=GammaConvV1_${trainConfig}_AOD194_list3_train479.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/479_20181109-1254_child_3/merge/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#done



# MC without mult weighting, without pT weighting, final cent classes ========================

# LHC16g1 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1060.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1060_20181115-0904/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1060.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1060_20181115-0904/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1060.root $path/$filename1 $path/$filename2
#done
# Runwise
#DownloadRunwise "GammaConvN" "AOD" "LHC16g1" "1060_20181115-0904" "sim" "pass1"
#DownloadRunwise "GammaConvN" "AOD" "LHC16g1" "1060_20181115-0904" "sim" "pass1_pidfix"


## LHC16h4 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 294 295 296 297; do   # normal and added particles
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1050.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1050_20181115-0903/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1050.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1050_20181115-0903/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1050.root $path/$filename1 $path/$filename2
#done

# LHC16g1a   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1047.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1047_20181115-0906/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1047.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1047_20181115-0906/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1047.root $path/$filename1 $path/$filename2
#done

# LHC16g1b + LHC16g1b_extra   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1048.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1048_20181115-0905_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1048.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1048_20181115-0905_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1048.root $path/$filename1 $path/$filename2
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1048.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1048_20181115-0905_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1048.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1048_20181115-0905_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1048.root $path/$filename1 $path/$filename2
#done

# LHC16g1c + LHC16g1c_extra     
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1049.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1049_20181115-0905_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1049.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1049_20181115-0905_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1049.root $path/$filename1 $path/$filename2
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1049.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1049_20181115-0905_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1049.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1049_20181115-0905_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1049.root $path/$filename1 $path/$filename2
#done


# LHC16i1a    
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1051.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1051_20181115-0905/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1051.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1051_20181115-0905/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1051.root $path/$filename1 $path/$filename2
#done

# LHC16i1b    
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1052.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1052_20181115-0905/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1052.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1052_20181115-0905/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1052.root $path/$filename1 $path/$filename2
#done


# LHC16i1c     
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1053.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1053_20181115-0905/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1053.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1053_20181115-0905/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1053.root $path/$filename1 $path/$filename2
#done

# LHC16i2a   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1054.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1054_20181115-0905/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1054.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1054_20181115-0905/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1054.root $path/$filename1 $path/$filename2
#done


# LHC16i2b      
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1055.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1055_20181115-0905/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1055.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1055_20181115-0905/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1055.root $path/$filename1 $path/$filename2
#done


# LHC16i2c   
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1056.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1056_20181115-0904/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1056.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1056_20181115-0904/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1056.root $path/$filename1 $path/$filename2
#done


# LHC16i3a 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1057.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1057_20181115-0904/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1057.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1057_20181115-0904/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1057.root $path/$filename1 $path/$filename2
#done


# LHC16i3b      
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1058.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1058_20181115-0904/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1058.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1058_20181115-0904/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1058.root $path/$filename1 $path/$filename2
#done


# LHC16i3c     
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1059.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1059_20181115-0904/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1059.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1059_20181115-0904/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1059.root $path/$filename1 $path/$filename2
#done


# MC with mult weighting, without pT weighting, final cent classes ========================

# LHC16g1 1066
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1066.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1066_20181121-1841/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1066.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1066_20181121-1841/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1066.root $path/$filename1 $path/$filename2
#done

# LHC16g1 1069 (update on mult weights)
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1069.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1069_20181205-1121/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1069.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1069_20181205-1121/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1069.root $path/$filename1 $path/$filename2
#done
# Runwise ! ATTENTION ! DONT OVERWRITE OLD FILES WITHOUT WEIGHTING BUT WITH QA PLOTS !
#DownloadRunwise "GammaConvN" "AOD" "LHC16g1" "1069_20181205-1121" "sim" "pass1"
#DownloadRunwise "GammaConvN" "AOD" "LHC16g1" "1069_20181205-1121" "sim" "pass1_pidfix"

# LHC16g1a   1070_20181206-1704 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1070.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1070_20181206-1704/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1070.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1070_20181206-1704/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1070.root $path/$filename1 $path/$filename2
#done

# LHC16g1b + LHC16g1b_extra 1071_20181206-1704    
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1071.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1071_20181206-1704_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1071.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1071_20181206-1704_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1071.root $path/$filename1 $path/$filename2
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1071.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1071_20181206-1704_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1071.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1071_20181206-1704_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1071.root $path/$filename1 $path/$filename2
#done

# LHC16g1c + LHC16g1c_extra   1072_20181206-1705      
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1072.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1072_20181206-1705_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1072.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1072_20181206-1705_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1072.root $path/$filename1 $path/$filename2
#done
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1072.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1072_20181206-1705_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1072.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1072_20181206-1705_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1072.root $path/$filename1 $path/$filename2
#done

## LHC16h4   1073_20181206-1705
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
#for trainConfig in 294 295 296 297; do   # normal and added particles
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1073.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1073_20181206-1705/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1073.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1073_20181206-1705/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1073.root $path/$filename1 $path/$filename2
#done

# LHC16i1a     1074_20181206-1705 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1074.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1074_20181206-1705/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1074.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1074_20181206-1705/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1074.root $path/$filename1 $path/$filename2
#done

# LHC16i1b    1075_20181206-1705 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1075.root
#   DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1075_20181206-1705/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1075.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1075_20181206-1705/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1075.root $path/$filename1 $path/$filename2
#done


# LHC16i1c      1076_20181206-1705 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1076.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1076_20181206-1705/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1076.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1076_20181206-1705/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1076.root $path/$filename1 $path/$filename2
#done

# LHC16i2a   1077_20181206-1705 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1077.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1077_20181206-1705/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1077.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1077_20181206-1705/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1077.root $path/$filename1 $path/$filename2
#done


# LHC16i2b      1078_20181206-1706 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1078.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1078_20181206-1706/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1078.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1078_20181206-1706/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1078.root $path/$filename1 $path/$filename2
#done


# LHC16i2c   1079_20181206-1706 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1079.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1079_20181206-1706/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1079.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1079_20181206-1706/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1079.root $path/$filename1 $path/$filename2
#done


# LHC16i3a   1080_20181206-1706 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1080.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1080_20181206-1706/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1080.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1080_20181206-1706/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1080.root $path/$filename1 $path/$filename2
#done


# LHC16i3b   1081_20181206-1706    
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1081.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1081_20181206-1706/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1081.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1081_20181206-1706/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1081.root $path/$filename1 $path/$filename2
#done


# LHC16i3c     1082_20181206-1707
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1082.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1082_20181206-1707/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1082.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1082_20181206-1707/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1082.root $path/$filename1 $path/$filename2
#done

# MC with mult weighting, with pT weighting, final cent classes ========================

# train configs 294,295, 500,501,504,505
: '
# LHC16g1   1084_20181213-1707 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1084.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1084_20181213-1707/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1084.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1084_20181213-1707/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1084.root $path/$filename1 $path/$filename2
done

# LHC16g1a   1085_20181213-1708  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1085.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1085_20181213-1708/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1085.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1085_20181213-1708/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1085.root $path/$filename1 $path/$filename2
done

# LHC16g1b + LHC16g1b_extra 1086_20181213-1708    
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1086.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1086_20181213-1708_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1086.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1086_20181213-1708_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1086.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1086.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1086_20181213-1708_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1086.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1086_20181213-1708_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1086.root $path/$filename1 $path/$filename2
done

# LHC16g1c + LHC16g1c_extra   1087_20181214-1609       
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1087.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1087_20181214-1609_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1087.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1087_20181214-1609_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1087.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1087.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1087_20181214-1609_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1087.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1087_20181214-1609_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1087.root $path/$filename1 $path/$filename2
done

## LHC16h4   1097_20181213-1717
# train configs 294,295,296,297, 500,501,502,503,504,505,506,507
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
for trainConfig in 500 501 502 503 504 505 506 507; do   # normal and added particles
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1097.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1097_20181213-1717/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1097.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1097_20181213-1717/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1097.root $path/$filename1 $path/$filename2
done

# LHC16i1a     1088_20181213-1710  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1088.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1088_20181213-1710/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1088.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1088_20181213-1710/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1088.root $path/$filename1 $path/$filename2
done

# LHC16i1b    1089_20181213-1710  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1089.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1089_20181213-1710/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1089.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1089_20181213-1710/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1089.root $path/$filename1 $path/$filename2
done


# LHC16i1c      1090_20181213-1710  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1090.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1090_20181213-1710/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1090.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1090_20181213-1710/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1090.root $path/$filename1 $path/$filename2
done

# LHC16i2a   1091_20181213-1716 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1091.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1091_20181213-1716/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1091.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1091_20181213-1716/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1091.root $path/$filename1 $path/$filename2
done


# LHC16i2b      1092_20181213-1716 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1092.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1092_20181213-1716/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1092.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1092_20181213-1716/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1092.root $path/$filename1 $path/$filename2
done


# LHC16i2c  1093_20181213-1716  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1093.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1093_20181213-1716/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1093.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1093_20181213-1716/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1093.root $path/$filename1 $path/$filename2
done


# LHC16i3a   1094_20181213-1716  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1094.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1094_20181213-1716/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1094.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1094_20181213-1716/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1094.root $path/$filename1 $path/$filename2
done


# LHC16i3b   1095_20181213-1716   
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1095.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1095_20181213-1716/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1095.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1095_20181213-1716/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1095.root $path/$filename1 $path/$filename2
done


# LHC16i3c     1096_20181213-1717 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 500 501 504 505; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1096.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1096_20181213-1717/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1096.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1096_20181213-1717/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1096.root $path/$filename1 $path/$filename2
done
'

#============= cut variations ==============================================
# standard and PsiPair variation
# normal particles:
# train configs 294 295 508 509 512 513
# added particles:
# train configs 296 297 510 511 514 515

: '
## LHC16h4
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
for trainConfig in 294 295 508 509 512 513; do  # normal particles
    # without added particles 1126_20190110-1729
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1126.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1126_20190110-1729/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1126.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1126_20190110-1729/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1126.root $path/$filename1 $path/$filename2
done
for trainConfig in 296 297 510 511 514 515; do  # added particles
    # with added pi0 and eta  1129_20190110-1730 
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1129.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1129_20190110-1730/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1129.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1129_20190110-1730/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1129.root $path/$filename1 $path/$filename2
    # with added eta          1127_20190110-1729
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1127.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1127_20190110-1729/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1127.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1127_20190110-1729/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1127.root $path/$filename1 $path/$filename2
    # with added pi0          1128_20190110-1729
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1128.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1128_20190110-1729/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1128.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1128_20190110-1729/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1128.root $path/$filename1 $path/$filename2
done

# LHC16g1    1113_20190110-1725
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1113.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1113_20190110-1725/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1113.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1113_20190110-1725/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1113.root $path/$filename1 $path/$filename2
done

# LHC16g1a  1114_20190110-1726
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1114.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1114_20190110-1726/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1114.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1114_20190110-1726/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1114.root $path/$filename1 $path/$filename2
done

# LHC16g1b + LHC16g1b_extra    1115_20190110-1726 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1115.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1115_20190110-1726_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1115.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1115_20190110-1726_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1115.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1115.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1115_20190110-1726_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1115.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1115_20190110-1726_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1115.root $path/$filename1 $path/$filename2
done

# LHC16g1c + LHC16g1c_extra      1116_20190110-1726    
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1116.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1116_20190110-1726_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1116.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1116_20190110-1726_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1116.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1116.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1116_20190110-1726_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1116.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1116_20190110-1726_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1116.root $path/$filename1 $path/$filename2
done

# LHC16i1a      1117_20190110-1727 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1117.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1117_20190110-1727/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1117.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1117_20190110-1727/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1117.root $path/$filename1 $path/$filename2
done

# LHC16i1b   1118_20190110-1727 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1118.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1118_20190110-1727/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1118.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1118_20190110-1727/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1118.root $path/$filename1 $path/$filename2
done


# LHC16i1c      1119_20190110-1727 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1119.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1119_20190110-1727/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1119.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1119_20190110-1727/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1119.root $path/$filename1 $path/$filename2
done

# LHC16i2a  1120_20190110-1727 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1120.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1120_20190110-1727/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1120.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1120_20190110-1727/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1120.root $path/$filename1 $path/$filename2
done


# LHC16i2b     1121_20190110-1727 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1121.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1121_20190110-1727/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1121.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1121_20190110-1727/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1121.root $path/$filename1 $path/$filename2
done


# LHC16i2c    1122_20190110-1727
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1122.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1122_20190110-1727/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1122.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1122_20190110-1727/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1122.root $path/$filename1 $path/$filename2
done


# LHC16i3a    1123_20190110-1728 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1123.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1123_20190110-1728/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1123.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1123_20190110-1728/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1123.root $path/$filename1 $path/$filename2
done


# LHC16i3b   1124_20190110-1728 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1124.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1124_20190110-1728/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1124.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1124_20190110-1728/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1124.root $path/$filename1 $path/$filename2
done


# LHC16i3c      1125_20190110-1728 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 294 295 508 509 512 513; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1125.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1125_20190110-1728/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1125.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1125_20190110-1728/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1125.root $path/$filename1 $path/$filename2
done
'

# LHC15o    492_20190110-1713 
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#for trainConfig in 294 295 508 509 512 513; do
#    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train492.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/492_20190110-1713_child_1/merge/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train492.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/492_20190110-1713_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_highIR_train492.root $path/$filename1 $path/$filename2
#done




#============= cut variations ==============================================
# pT cut variation
# normal particles: 516, 517, 520, 521, 524, 525
# added particles: 518, 519, 522, 523, 526, 527

: '
# LHC15o    
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train493.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/493_20190118-1732_child_1/merge/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train493.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/493_20190118-1732_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_highIR_train493.root $path/$filename1 $path/$filename2
done

## LHC16h4
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
for trainConfig in 516 517 520 521 524 525; do  # normal particles
    # without added particles 1143_20190118-1745 
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1143.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1143_20190118-1745/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1143.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1143_20190118-1745/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1143.root $path/$filename1 $path/$filename2
done
for trainConfig in 518 519 522 523 526 527; do  # added particles
    # with added pi0 and eta   1146_20190118-1743
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1146.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1146_20190118-1743/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1146.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1146_20190118-1743/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1146.root $path/$filename1 $path/$filename2
    # with added eta          1144_20190118-1743 
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1144.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1144_20190118-1743/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1144.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1144_20190118-1743/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1144.root $path/$filename1 $path/$filename2
    # with added pi0          1145_20190118-1743 
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1145.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1145_20190118-1743/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1145.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1145_20190118-1743/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1145.root $path/$filename1 $path/$filename2
done

# LHC16g1    1130_20190118-1738 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1130.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1130_20190118-1738/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1130.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1130_20190118-1738/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1130.root $path/$filename1 $path/$filename2
done

# LHC16g1a   1131_20190118-1739 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1131.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1131_20190118-1739/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1131.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1131_20190118-1739/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1131.root $path/$filename1 $path/$filename2
done

# LHC16g1b + LHC16g1b_extra    1132_20190118-1739
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1132.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1132_20190118-1739_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1132.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1132_20190118-1739_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1132.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1132.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1132_20190118-1739_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1132.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1132_20190118-1739_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1132.root $path/$filename1 $path/$filename2
done

# LHC16g1c + LHC16g1c_extra         1133_20190118-1739 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1133.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1133_20190118-1739_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1133.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1133_20190118-1739_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1133.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1133.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1133_20190118-1739_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1133.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1133_20190118-1739_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1133.root $path/$filename1 $path/$filename2
done

# LHC16i1a    1134_20190118-1741 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1134.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1134_20190118-1741/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1134.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1134_20190118-1741/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1134.root $path/$filename1 $path/$filename2
done

# LHC16i1b   1135_20190118-1742
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1135.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1135_20190118-1742/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1135.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1135_20190118-1742/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1135.root $path/$filename1 $path/$filename2
done


# LHC16i1c     1136_20190118-1742 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1136.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1136_20190118-1742/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1136.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1136_20190118-1742/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1136.root $path/$filename1 $path/$filename2
done

# LHC16i2a    1137_20190118-1743 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1137.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1137_20190118-1743/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1137.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1137_20190118-1743/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1137.root $path/$filename1 $path/$filename2
done


# LHC16i2b     1138_20190118-1743 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1138.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1138_20190118-1743/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1138.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1138_20190118-1743/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1138.root $path/$filename1 $path/$filename2
done


# LHC16i2c   1139_20190118-1743
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1139.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1139_20190118-1743/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1139.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1139_20190118-1743/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1139.root $path/$filename1 $path/$filename2
done


# LHC16i3a     1140_20190118-1744 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1140.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1140_20190118-1744/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1140.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1140_20190118-1744/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1140.root $path/$filename1 $path/$filename2
done


# LHC16i3b    1141_20190118-1744 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1141.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1141_20190118-1744/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1141.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1141_20190118-1744/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1141.root $path/$filename1 $path/$filename2
done


# LHC16i3c      1142_20190118-1744 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 516 517 520 521 524 525; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1142.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1142_20190118-1744/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1142.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1142_20190118-1744/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1142.root $path/$filename1 $path/$filename2
done
'

#============= cut variations ==============================================
# TPC cluster cut variation       + standard
# normal particles: 528 529 532 533 294 295 
# added particles:  530 531 534 535 296 297

: '
# LHC15o    494_20190124-1507
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train494.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/494_20190124-1507_child_1/merge/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train494.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/494_20190124-1507_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_highIR_train494.root $path/$filename1 $path/$filename2
done

## LHC16h4 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
for trainConfig in 528 529 532 533 294 295; do  # normal particles
    # without added particles   1161_20190125-1722
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1161.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1161_20190125-1722/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1161.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1161_20190125-1722/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1161.root $path/$filename1 $path/$filename2
done
for trainConfig in 530 531 534 535 296 297; do  # added particles
    # with added pi0 and eta    1163_20190129-1252 
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1163.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1163_20190129-1252/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1163.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1163_20190129-1252/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1163.root $path/$filename1 $path/$filename2
    # with added eta         1162_20190129-1252   
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1162.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1162_20190129-1252/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1162.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1162_20190129-1252/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1162.root $path/$filename1 $path/$filename2
    # with added pi0       1161_20190125-1722    
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1161.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1161_20190125-1722/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1161.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1161_20190125-1722/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1161.root $path/$filename1 $path/$filename2
done

# LHC16g1   1148_20190124-1511 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1148.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1148_20190124-1511/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1148.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1148_20190124-1511/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1148.root $path/$filename1 $path/$filename2
done

# LHC16g1a   1149_20190124-1512 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1a/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1149.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1149_20190124-1512/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1149.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1149_20190124-1512/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1149.root $path/$filename1 $path/$filename2
done

# LHC16g1b + LHC16g1b_extra   1150_20190124-1512 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1150.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1150_20190124-1512_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1150.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1150_20190124-1512_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1150.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1b_extra/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1150.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1150_20190124-1512_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1150.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1150_20190124-1512_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1150.root $path/$filename1 $path/$filename2
done

# LHC16g1c + LHC16g1c_extra     1151_20190124-1512     
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1151.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1151_20190124-1512_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1151.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1151_20190124-1512_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1151.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16g1c_extra/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD_list1_train1151.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1151_20190124-1512_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD_list2_train1151.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1151_20190124-1512_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD_highIR_train1151.root $path/$filename1 $path/$filename2
done

# LHC16i1a     1152_20190124-1512 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1152.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1152_20190124-1512/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1152.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1152_20190124-1512/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1152.root $path/$filename1 $path/$filename2
done

# LHC16i1b   1153_20190124-1512
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1153.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1153_20190124-1512/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1153.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1153_20190124-1512/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1153.root $path/$filename1 $path/$filename2
done


# LHC16i1c     1154_20190124-1512 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1154.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1154_20190124-1512/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1154.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1154_20190124-1512/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1154.root $path/$filename1 $path/$filename2
done

# LHC16i2a    1155_20190124-1512 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1155.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1155_20190124-1512/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1155.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1155_20190124-1512/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1155.root $path/$filename1 $path/$filename2
done


# LHC16i2b     1156_20190124-1513 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1156.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1156_20190124-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1156.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1156_20190124-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1156.root $path/$filename1 $path/$filename2
done


# LHC16i2c   1157_20190124-1513 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1157.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1157_20190124-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1157.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1157_20190124-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1157.root $path/$filename1 $path/$filename2
done


# LHC16i3a     1158_20190124-1513
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1158.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1158_20190124-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1158.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1158_20190124-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1158.root $path/$filename1 $path/$filename2
done


# LHC16i3b   1159_20190124-1513 

path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1159.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1159_20190124-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1159.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1159_20190124-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1159.root $path/$filename1 $path/$filename2
done


# LHC16i3c     1160_20190124-1513 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 528 529 532 533 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1160.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1160_20190124-1513/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1160.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1160_20190124-1513/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1160.root $path/$filename1 $path/$filename2
done
'

#============= cut variations ==============================================
# Chi2 variations + small variations
# Attention with MCs: cross check without pT weighting
# normal particles: 294 295 500 501 504 505 564 565 
# added particles:  296 297 502 503 506 507 566 567  

: '
# LHC15o    497_20190208-1200 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train497.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/497_20190208-1200_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train497.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/497_20190208-1200_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_highIR_train497.root $path/$filename1 $path/$filename2
done

# LHC16i1a     1177_20190212-1413 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1177.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1177_20190212-1413/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1177.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1177_20190212-1413/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1177.root $path/$filename1 $path/$filename2
done

# LHC16i1b    1178_20190212-1413
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1178.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1178_20190212-1413/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1178.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1178_20190212-1413/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1178.root $path/$filename1 $path/$filename2
done


# LHC16i1c     1179_20190212-1409 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1179.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1179_20190212-1409/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1179.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1179_20190212-1409/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1179.root $path/$filename1 $path/$filename2
done


# LHC16i2a   1191_20190212-1407  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1191.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1191_20190212-1407/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1191.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1191_20190212-1407/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1191.root $path/$filename1 $path/$filename2
done


# LHC16i2b    1181_20190212-1412 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1181.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1181_20190212-1412/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1181.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1181_20190212-1412/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1181.root $path/$filename1 $path/$filename2
done

# LHC16i2c    1182_20190212-1412 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1182.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1182_20190212-1412/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1182.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1182_20190212-1412/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1182.root $path/$filename1 $path/$filename2
done

# LHC16i3a     1183_20190212-1412 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1183.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1183_20190212-1412/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1183.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1183_20190212-1412/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1183.root $path/$filename1 $path/$filename2
done


# LHC16i3b     1184_20190212-1412 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1184.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1184_20190212-1412/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1184.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1184_20190212-1412/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1184.root $path/$filename1 $path/$filename2
done

# LHC16i3c     1185_20190212-1412
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 294 295 500 501 504 505 564 565; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1185.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1185_20190212-1412/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1185.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1185_20190212-1412/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1185.root $path/$filename1 $path/$filename2
done
'

#============= new cent class 0-20% for eta ====================================
# normal particles: 294 295

# LHC15o   501_20190222-1453
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#for trainConfig in 294 295; do
#    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train501.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/501_20190222-1453_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train501.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/501_20190222-1453_child_2/merge/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_highIR_train501.root $path/$filename1 $path/$filename2
#done

: '
#============= new MCs without weights ====================================
# normal particles: 294 295


# LHC18e1 1197_20190221-1136
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1197.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1197_20190221-1136/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1197.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1197_20190221-1136/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1197.root $path/$filename1 $path/$filename2
done

# LHC18e1a   1198_20190221-1338
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1198.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1198_20190221-1338/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1198.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1198_20190221-1338/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1198.root $path/$filename1 $path/$filename2
done


# LHC18e1b + LHC18e1b_extra   1199_20190221-1150 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1199.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1199_20190221-1150_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1199.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1199_20190221-1150_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1199.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b_extra/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1199.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1199_20190221-1150_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1199.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1199_20190221-1150_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1199.root $path/$filename1 $path/$filename2
done

# LHC18e1c + LHC18e1c_extra     1200_20190221-1150    
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1200.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1200_20190221-1150_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1200.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1200_20190221-1150_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1200.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c_extra/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1200.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1200_20190221-1150_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1200.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1200_20190221-1150_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1200.root $path/$filename1 $path/$filename2
done

'

#============= new MCs with mult weights ====================================
# normal particles: 294 295

: '

# LHC18e1 1218_20190306-1441
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1218.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1218_20190306-1441/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1218.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1218_20190306-1441/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1218.root $path/$filename1 $path/$filename2
done

# LHC18e1a  1214_20190306-0935
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1214.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1214_20190306-0935/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1214.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1214_20190306-0935/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1214.root $path/$filename1 $path/$filename2
done


# LHC18e1b + LHC18e1b_extra    1215_20190306-1441 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1215.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1215_20190306-1441_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1215.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1215_20190306-1441_child_1/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1215.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b_extra/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1215.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1215_20190306-1441_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1215.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1215_20190306-1441_child_2/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1215.root $path/$filename1 $path/$filename2
done

# LHC18e1c + LHC18e1c_extra       1216_20190306-0936 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1216.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1216_20190306-0936_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1216.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1216_20190306-0936_child_1/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1216.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c_extra/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1216.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1216_20190306-0936_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1216.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1216_20190306-0936_child_2/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1216.root $path/$filename1 $path/$filename2
done

## LHC16h4 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
for trainConfig in 294 295; do  # normal particles
    # without added particles   1217_20190307-0509
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1217.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1217_20190307-0509/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1217.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1217_20190307-0509/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1217.root $path/$filename1 $path/$filename2
done
for trainConfig in 296 297; do  # added particles
    # with added pi0 and eta    1217_20190307-0509
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1217.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1217_20190307-0509/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1217.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1217_20190307-0509/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1217.root $path/$filename1 $path/$filename2
    # with added eta           1220_20190307-0509
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1220.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1220_20190307-0509/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1220.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1220_20190307-0509/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1220.root $path/$filename1 $path/$filename2
    # with added pi0         1219_20190307-0509
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1219.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1219_20190307-0509/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1219.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1219_20190307-0509/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1219.root $path/$filename1 $path/$filename2
done

# LHC16i1a    1221_20190308-0943 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1221.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1221_20190308-0943/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1221.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1221_20190308-0943/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1221.root $path/$filename1 $path/$filename2
done

# LHC16i1b    1222_20190308-0943 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1222.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1222_20190308-0943/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1222.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1222_20190308-0943/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1222.root $path/$filename1 $path/$filename2
done


# LHC16i1c    1223_20190308-0943 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1223.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1223_20190308-0943/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1223.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1223_20190308-0943/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1223.root $path/$filename1 $path/$filename2
done


# LHC16i2a    1224_20190308-0947 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1224.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1224_20190308-0947/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1224.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1224_20190308-0947/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1224.root $path/$filename1 $path/$filename2
done


# LHC16i2b   1225_20190308-0947 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1225.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1225_20190308-0947/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1225.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1225_20190308-0947/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1225.root $path/$filename1 $path/$filename2
done

# LHC16i2c   1226_20190308-0947 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1226.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1226_20190308-0947/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1226.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1226_20190308-0947/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1226.root $path/$filename1 $path/$filename2
done

# LHC16i3a     1227_20190308-0944 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1227.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1227_20190308-0944/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1227.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1227_20190308-0944/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1227.root $path/$filename1 $path/$filename2
done


# LHC16i3b    1228_20190308-0944 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1228.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1228_20190308-0944/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1228.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1228_20190308-0944/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1228.root $path/$filename1 $path/$filename2
done

# LHC16i3c     1229_20190308-0944 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1229.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1229_20190308-0944/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1229.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1229_20190308-0944/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1229.root $path/$filename1 $path/$filename2
done

'

#============= new MCs with mult weights and pT weights it #1 ====================================
# normal particles: 294 295
# added particles:  296 297

: '

# LHC18e1      1235_20190327-1626
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1235.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1235_20190327-1626/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1235.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1235_20190327-1626/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1235.root $path/$filename1 $path/$filename2
done

# LHC18e1a  1236_20190328-1606 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1236.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1236_20190328-1606/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1236.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1236_20190328-1606/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1236.root $path/$filename1 $path/$filename2
done


# LHC18e1b + LHC18e1b_extra    1237_20190327-1626
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1237.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1237_20190327-1626_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1237.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1237_20190327-1626_child_1/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1237.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b_extra/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1237.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1237_20190327-1626_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1237.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1237_20190327-1626_child_2/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1237.root $path/$filename1 $path/$filename2
done

# LHC18e1c + LHC18e1c_extra       1238_20190327-1626 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1238.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1238_20190327-1626_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1238.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1238_20190327-1626_child_1/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1238.root $path/$filename1 $path/$filename2
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c_extra/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1238.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1238_20190327-1626_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1238.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1238_20190327-1626_child_2/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1238.root $path/$filename1 $path/$filename2
done

## LHC16h4 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
for trainConfig in 294 295; do  # normal particles
    # without added particles   1248_20190327-1625
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1248.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1248_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1248.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1248_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1248.root $path/$filename1 $path/$filename2
done
for trainConfig in 296 297; do  # added particles
    # with added pi0 and eta    1248_20190327-1625
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1248.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1248_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1248.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1248_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1248.root $path/$filename1 $path/$filename2
    # with added eta           1250_20190327-1625 
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1250.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1250_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1250.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1250_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1250.root $path/$filename1 $path/$filename2
    # with added pi0        1249_20190327-1625
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1249.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1249_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1249.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1249_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1249.root $path/$filename1 $path/$filename2
done

# LHC16i1a    1239_20190328-1314  
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1239.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1239_20190328-1314/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1239.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1239_20190328-1314/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1239.root $path/$filename1 $path/$filename2
done

# LHC16i1b     1240_20190327-1626 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1240.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1240_20190327-1626/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1240.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1240_20190327-1626/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1240.root $path/$filename1 $path/$filename2
done


# LHC16i1c   1241_20190327-1626 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1241.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1241_20190327-1626/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1241.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1241_20190327-1626/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1241.root $path/$filename1 $path/$filename2
done


# LHC16i2a     1242_20190328-1315 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1242.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1242_20190328-1315/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1242.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1242_20190328-1315/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1242.root $path/$filename1 $path/$filename2
done


# LHC16i2b   1243_20190327-1625
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1243.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1243_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1243.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1243_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1243.root $path/$filename1 $path/$filename2
done

# LHC16i2c   1244_20190327-1625 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1244.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1244_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1244.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1244_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1244.root $path/$filename1 $path/$filename2
done

# LHC16i3a      1245_20190328-1315
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1245.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1245_20190328-1315/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1245.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1245_20190328-1315/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1245.root $path/$filename1 $path/$filename2
done


# LHC16i3b    1246_20190327-1625 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1246.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1246_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1246.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1246_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1246.root $path/$filename1 $path/$filename2
done

# LHC16i3c     1247_20190327-1625 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 294 295; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1247.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1247_20190327-1625/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1247.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1247_20190327-1625/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1247.root $path/$filename1 $path/$filename2
done

'

#============= new MCs with mult weights and pT weights it #1 ====================================
#                   with TOF cut    without TOF cut     with material budget weights
# normal particles: 294 295         568 569             572 573  
# added particles:  296 297         570 571             574 575


# LHC15o   505_20190411-0852
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC15o/
#for trainConfig in 294 295 568 569 572 573; do
#    filename1=GammaConvV1_${trainConfig}_AOD194_list1_train505.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/505_20190411-0852_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD194_list2_train505.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/505_20190411-0852_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD194_list3_train505.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_AOD/505_20190411-0852_child_3/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_highIR_train505.root $path/$filename1 $path/$filename2
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD194_train505.root $path/$filename1 $path/$filename2 $path/$filename3
#done

#DownloadRunwise "GammaConvN" "AOD" "LHC15o" "505_20190411-0852_child_1" "data" "pass1"
#DownloadRunwise "GammaConvN" "AOD" "LHC15o" "505_20190411-0852_child_2" "data" "pass1_pidfix"
#DownloadRunwise "GammaConvN" "AOD" "LHC15o" "505_20190411-0852_child_3" "data" "pass3_lowIR_pidfix"

DownloadRunwise "GammaConvN" "AOD" "LHC18e1" "1261_20190411-1045" "sim" "pass1"
DownloadRunwise "GammaConvN" "AOD" "LHC18e1" "1261_20190411-1045" "sim" "pass1_pidfix"
DownloadRunwise "GammaConvN" "AOD" "LHC18e1" "1261_20190411-1045" "sim" "pass3_lowIR_pidfix"

: '

# LHC18e1        1261_20190411-1045
#path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1/
#for trainConfig in 294 295 568 569 572 573; do
#    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1261.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1261_20190411-1045/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
#    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1261.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1261_20190411-1045/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
#    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1261.root
#    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1261_20190411-1045/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1261.root $path/$filename1 $path/$filename2
#    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1261.root $path/$filename1 $path/$filename2 $path/$filename3    
#done

# LHC18e1a          1262_20190411-1045 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1a/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1262.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1262_20190411-1045/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1262.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1262_20190411-1045/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1262.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1262_20190411-1045/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1262.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1262.root $path/$filename1 $path/$filename2 $path/$filename3
done


# LHC18e1b + LHC18e1b_extra       1263_20190411-1046 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1263.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1263_20190411-1046_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1263.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1263_20190411-1046_child_1/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1263.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1263_20190411-1046_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1263.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1263.root $path/$filename1 $path/$filename2 $path/$filename3
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1b_extra/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1263.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1263_20190411-1046_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1263.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1263_20190411-1046_child_2/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1263.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1263_20190411-1046_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1263.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1263.root $path/$filename1 $path/$filename2 $path/$filename3
done

# LHC18e1c + LHC18e1c_extra         1264_20190411-1046
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1264.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1264_20190411-1046_child_1/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1264.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1264_20190411-1046_child_1/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1264.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1264_20190411-1046_child_1/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1264.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1264.root $path/$filename1 $path/$filename2 $path/$filename3
done
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC18e1c_extra/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1264.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1264_20190411-1046_child_2/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1264.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1264_20190411-1046_child_2/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1264.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1264_20190411-1046_child_2/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1264.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198__train1264.root $path/$filename1 $path/$filename2 $path/$filename3
done

## LHC16h4      1251_20190411-1044 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16h4/
for trainConfig in 294 295 568 569 572 573; do  # normal particles
    # without added particles   
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1251.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1251_20190411-1044/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1251.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1251_20190411-1044/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1251.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1251_20190411-1044/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1251.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1251.root $path/$filename1 $path/$filename2 $path/$filename3
done
for trainConfig in 296 297 570 571 574 575; do  # added particles
    # with added pi0 and eta    
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1251.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1251_20190411-1044/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1251.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1251_20190411-1044/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1251.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1251_20190411-1044/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1251.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1251.root $path/$filename1 $path/$filename2 $path/$filename3
    # with added eta           
    # with added pi0      
done

# LHC16i1a     1252_20190411-1044 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1a/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1252.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1252_20190411-1044/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1252.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1252_20190411-1044/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1252.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1252_20190411-1044/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1252.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1252.root $path/$filename1 $path/$filename2 $path/$filename3
done

# LHC16i1b     1253_20190411-1044 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1b/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1253.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1253_20190411-1044/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1253.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1253_20190411-1044/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1253.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1253_20190411-1044/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1253.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1253.root $path/$filename1 $path/$filename2 $path/$filename3
done


# LHC16i1c   1254_20190411-1044 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i1c/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1254.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1254_20190411-1044/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1254.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1254_20190411-1044/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1254.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1254_20190411-1044/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1254.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1254.root $path/$filename1 $path/$filename2 $path/$filename3
done


# LHC16i2a      1255_20190411-1044
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2a/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1255.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1255_20190411-1044/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1255.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1255_20190411-1044/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1255.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1255_20190411-1044/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"    
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1255.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1255.root $path/$filename1 $path/$filename2 $path/$filename3
done


# LHC16i2b     1256_20190411-1044 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2b/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1256.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1256_20190411-1044/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1256.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1256_20190411-1044/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1256.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1256_20190411-1044/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1256.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1256.root $path/$filename1 $path/$filename2 $path/$filename3
done

# LHC16i2c     1257_20190411-1045 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i2c/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1257.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1257_20190411-1045/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1257.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1257_20190411-1045/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1257.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1257_20190411-1045/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1257.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1257.root $path/$filename1 $path/$filename2 $path/$filename3
done

# LHC16i3a     1258_20190411-1045 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3a/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1258.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1258_20190411-1045/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1258.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1258_20190411-1045/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1258.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1258_20190411-1045/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1258.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1258.root $path/$filename1 $path/$filename2 $path/$filename3
done


# LHC16i3b    1259_20190411-1045 
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3b/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1259.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1259_20190411-1045/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1259.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1259_20190411-1045/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1259.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1259_20190411-1045/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1259.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1259.root $path/$filename1 $path/$filename2 $path/$filename3
done

# LHC16i3c    1260_20190411-1045
path=/home/meike/analysis/data/GridOutput/GammaConv/PbPb/AOD/LHC16i3c/
for trainConfig in 294 295 568 569 572 573; do
    filename1=GammaConvV1_${trainConfig}_AOD198_list1_train1260.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1260_20190411-1045/merge_runlist_1/GammaConvV1_$trainConfig.root" "$path/$filename1" "yes"
    filename2=GammaConvV1_${trainConfig}_AOD198_list2_train1260.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1260_20190411-1045/merge_runlist_2/GammaConvV1_$trainConfig.root" "$path/$filename2" "yes"
    filename3=GammaConvV1_${trainConfig}_AOD198_list3_train1260.root
    DownloadMerged "/alice/cern.ch/user/a/alitrain/PWGGA/GA_PbPb_MC_AOD/1260_20190411-1045/merge_runlist_3/GammaConvV1_$trainConfig.root" "$path/$filename3" "yes"
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_highIR_train1260.root $path/$filename1 $path/$filename2
    hadd -f $path/GammaConvV1_${trainConfig}_AOD198_train1260.root $path/$filename1 $path/$filename2 $path/$filename3
done

'
