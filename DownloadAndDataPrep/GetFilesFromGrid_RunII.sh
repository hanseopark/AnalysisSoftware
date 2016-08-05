#! /bin/bash

# This macro copies GammaConvV1_<trainConfig> or AnalysisResults.root files from grid and puts them in the correct directories

#****** arguments ******************************************************************************************************************
# $1 = period 
# $2 = "G"       (Download GammaConv files)
#      "A"       (Download AnalysisResults.root files)
# $3 = "runwise" (download runwise. If $2="G": optional merge afterwards) 
#      "merge"   (merge ONLY. Only for $2="G")

#***** variables to be set manually ******************************************************************************************************************
# specify where the files should be saved:
#     $DIR/GammaConv/$SYSTEM/$PERIOD/$run/GammaConvV1_<trainconfig>.root
# or  $DIR/PhotonQA/$SYSTEM/$PERIOD/$run/AnalysisResults.root
DIR=/home/meike/analysis/data/GridOutput

#****** automatic settings according to period given as argument ****************
# update when there is a new production or a new production pass

case "$1" in    # default: echo "no valid period chosen";

    "LHC15f")
	PERIOD=LHC15f
	YEAR=2015
	isMC=FALSE
	SYSTEM=pp
	PASS=pass2
	# listPCMgood_kINT7:
	RUN=(225026 225031 225035 225037 225041 225043 225050 225051 225052 225106 225305 225307 225309 225310 225313 225314 225315 225322 225576 225578 225579 225586 225587 225705 225707 225708 225709 225710 225716 225717 225719 225753 225757 225762 225763 225766 225768 226062 226170 226220 226225 226444 226445 226452 226466 226468 226472 226476 226483 226495 226500);
	AODNO=171
	;;
    "LHC15g3a3") # anchored to 15f	
	PERIOD=LHC15g3a3
	YEAR=2015
	isMC=TRUE
	SYSTEM=pp
	# no pass
	AODNO=176
	# listPCMgood_kINT7: 
	RUN=(225026 225031 225035 225037 225041 225043 225050 225051 225052 225106 225305 225307 225309 225310 225313 225314 225315 225322 225576 225578 225579 225586 225587 225705 225707 225708 225709 225710 225716 225717 225719 225753 225757 225762 225763 225766 225768 226062 226170 226220 226225 226444 226445 226452 226466 226468 226472 226476 226483 226495 226500);
	;;
    "LHC15g3c3") # anchored to 15f
	PERIOD=LHC15g3c3
	YEAR=2015
	isMC=TRUE
	SYSTEM=pp
	# no pass
	AODNO=176
	# listPCMgood_kINT7:
	RUN=(225026 225031 225035 225037 225041 225043 225050 225051 225052 225106 225305 225307 225309 225310 225313 225314 225315 225322 225576 225578 225579 225586 225587 225705 225707 225708 225709 225710 225716 225717 225719 225753 225757 225762 225763 225766 225768 226062 226170 226220 226225 226444 226445 226452 226466 226468 226472 226476 226483 226495 226500);
	;;
    "LHC15g")
	PERIOD=LHC15g
	YEAR=2015
	isMC=FALSE
	SYSTEM=pp
	PASS=pass1_megarun
	;;
    "LHC15h")
	PERIOD=LHC15h
	YEAR=2015
	isMC=FALSE
	SYSTEM=pp
	PASS=pass1
	RUN=(232914 232916 232986 232993 232995 233020 233059 233061 233093 233116 233120 233144 233169 233217 233232 233239 233242 233361 233465 233472 233614 233621 233623 233627 233678 233686 233692 233696 233697 233698 233700 233716 233719 233720 233721 233743 233799 233828 233830 233837 233858 233912 233969 233971 233972 233973 233974 233975 233976 233977 233978 234031 234039 234040 234043 234045 234048 234049 234050);
	;;
    "LHC15i")
	PERIOD=LHC15i
	YEAR=2015
	isMC=FALSE
	SYSTEM=pp
	PASS=pass1
	RUN=(236866 236863 236862 236860 236850 236848 236835 236824 236822 236821 236816 236815 236814 236813 236569 236565 236564 236563 236562 236558 236556 236444 236443 236441  236393 236389 236360 236359 236357 236356 236354 236353 236352 236349 236348 236337 236334 236331 236284 236281 236248 236246 236244 236242 236240 236238 236234 236227 236222 236204 236203 236164 236163 236161 236158 236153 236151 236150 236138 236137 236062 235898 235897 235896 235895 235893 235891 235890 235889 235888 235886 235841 235839 235811 235759 235721 235714 235710 235694 235684 235683 235573 235547 235462 235459 235454 235450 235443 235436 235435 235432 235423 235383 235380 235364 235362 235356 235347 235345 235344 235245 235242 235226 235204 235203 235201 235196);
	;;
    "LHC15oLowIR")
	PERIOD=LHC15o
	YEAR=2015
	isMC=FALSE
	SYSTEM=PbPb
	PASS=pass2_lowIR
	# all 13 runs:
	RUN=(244917 244918 244975 244980 244982 244983 245061 245064 245066 245068 246390 246391 246392);
	;;
    "LHC15k1a1") # anchored to LHC15oLowIR
	PERIOD=LHC15k1a1
	YEAR=2015
	isMC=TRUE
	SYSTEM=PbPb
	# no pass
	# 10 runs:
	RUN=(244918 244975 244982 244983 245064 245066 245068 246390 246391 246392);
	;;
    "LHC15k1a2") # anchored to LHC15oLowIR
	PERIOD=LHC15k1a2
	YEAR=2015
	isMC=TRUE
	SYSTEM=PbPb
	# no pass
	# 10 runs:
	RUN=(244918 244975 244982 244983 245064 245066 245068 246390 246391 246392);
	;;
    "LHC15k1a3") # anchored to LHC15oLowIR
	PERIOD=LHC15k1a3
	YEAR=2015
	isMC=TRUE
	SYSTEM=PbPb
	# no pass
	# 10 runs:
	RUN=(244918 244975 244982 244983 245064 245066 245068 246390 246391 246392);
	;;
    "LHC15oHighIR")
	PERIOD=LHC15o
	YEAR=2015
	isMC=FALSE
	SYSTEM=PbPb
	PASS=pass1
	# all 31 runs:
	RUN=(245683 246807 246808 246809 246810 246844 246845 246846 246847 246851 246855 246858 246859 246864 246865 246867 246870 246871 246928 246930 246937 246942 246945 246948 246949 246980 246982 246984 246989 246991 246994);
	;;
esac

echo "chose period: " $PERIOD", year: "$YEAR", MC: "$isMC", pass: "$PASS", collision system: "$SYSTEM", refiltering number if available: "$AODNO;

if [ $2 = "G" ]; then
    echo "Chose GammaConvV1_<trainconfig>.root";
    FILETYPE=G
    BASEDIR=$DIR/GammaConv/
elif [ $2 = "A" ]; then
    echo "Chose AnalysisResults.root"
    FILETYPE=A
    BASEDIR=$DIR/PhotonQA/
else
    echo "second argument should be G or A";
fi

# ********* common functions **********************************************************************************************************
Structure() {
    if [ "$FILETYPE" == "G" ]; then
	ls $OUTPUTDIR/GammaConvV1_*.root > fileData.txt
	fileNumbers=`cat fileData.txt`
	for fileName in $fileNumbers; do
	    number=`echo $fileName | cut -d "_" -f 2 | cut -d "." -f 1`
	    echo "trainconfig: " $number
	    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1.root\"\,\"GammaConvV1_$number\"\)
	done
    fi
}

Merge() {

	echo "merge ${#RUN[@]} files";               # not the actual number of downloaded runs!    PRODUCE FILE WITH MERGED RUN NUMBERS!!!
	declare -i i=0
	declare -i actualNoRuns=${#RUN[@]}
	for run in "${RUN[@]}"; do
	    i=i+1
	    echo "Run $i of ${#RUN[@]} ...";
	    if [ -f $BASEDIR/$SYSTEM/$PERIOD/$run/GammaConvV1.root ]; then
		if [ $i = 1 ]; then
		    mkdir -p $BASEDIR/$SYSTEM/$PERIOD/mymerge/${#RUN[@]}runs/
		    cp $BASEDIR/$SYSTEM/$PERIOD/$run/GammaConvV1.root $BASEDIR/$SYSTEM/$PERIOD/mymerge/${#RUN[@]}runs/GammaConvV1.root
		else
		    hadd -f ./GammaConvV1.root $BASEDIR/$SYSTEM/$PERIOD/mymerge/${#RUN[@]}runs/GammaConvV1.root $BASEDIR/$SYSTEM/$PERIOD/$run/GammaConvV1.root    # merge every file to mymerge/GammaConv/GammaConvV1.root
		    mv GammaConvV1.root $BASEDIR/$SYSTEM/$PERIOD/mymerge/${#RUN[@]}runs/
		fi
	    else
		    echo "file "$BASEDIR/$SYSTEM/$PERIOD/$run/"GammaConvV1.root does not exist!";
		    actualNoRuns=actualNoRuns-1
	    fi      
	done
	echo $actualNoRuns" runs merged";
	#if [ $actualNoRuns -ne ${#RUN[@]} ]; then
	mv $BASEDIR/$SYSTEM/$PERIOD/mymerge/${#RUN[@]}runs/Gam* $BASEDIR/$SYSTEM/$PERIOD/GammaConvV1_$actualNoRuns"runs.root"
	rm -r $BASEDIR/$SYSTEM/$PERIOD/mymerge/
	#fi
}

WriteInfo() {
    InfoExists=`ls $BASEDIR/$SYSTEM/$PERIOD/info.txt`
    if [ "$InfoExists" == "" ]; then
	echo -e "$TRAIN\n$PERIOD\n$PASS" > $BASEDIR/$SYSTEM/$PERIOD/info.txt
	echo "info.txt created"
	cat $BASEDIR/$SYSTEM/$PERIOD/info.txt
    else
	echo "Do you want to overwrite Info.txt? yes / no";
        cat $BASEDIR/$SYSTEM/$PERIOD/info.txt
        read answerInfo
	if [ "$answerInfo" != "yes" ] && [ "$answerInfo" != "no" ]; then
	    echo "no valid answer has been given. Will assume no.";
	    answerInfo=no;
	fi
	if [ "$answerInfo" == "yes" ]; then
	   echo -e "$TRAIN\n$PERIOD\n$PASS" > $BASEDIR/$SYSTEM/$PERIOD/info.txt  # write also ESD/AOD, AliPhysics version!!
	   echo "info.txt written"
	   cat $BASEDIR/$SYSTEM/$PERIOD/info.txt
	fi
    fi
}

OldData() {
    if [ -f $BASEDIR/$SYSTEM/$PERIOD/info.txt ]; then
	echo "Data for this period already exists: "
	cat $BASEDIR/$SYSTEM/$PERIOD/info.txt
	echo "What do you want to do? move and download (m) download (d) abort (a)";
	read answerMove  
	if [ "$answerMove" == "m" ]; then
	    line=$(head -n 1 $BASEDIR/$SYSTEM/$PERIOD/info.txt)
	    mkdir -p $BASEDIR/$SYSTEM/$PERIOD/old/$PASS/$line                      
	    mv $BASEDIR/$SYSTEM/$PERIOD/* $BASEDIR/$SYSTEM/$PERIOD/old/$PASS/$line
	    echo "moved old data to $BASEDIR/$PERIOD/old/$pass/$line.";
	    echo "Go on with download...";
	elif [ "$answerMove" == "d" ]; then
	    echo "go on with download...";
	elif [ "$answerMove" == "a" ]; then
	    echo "Abort.";
	    exit
	else
	    echo "No valid answer has been given. Abort.";
	    exit
	fi
    fi
} # additional variables: answerMove line

Show() {   # -> split into two functions, one being GetSourceDir
# shows runwise files

if [ "$isMC" == "FALSE" ];then # data
    if [ "$isAOD" == "FALSE" ]; then
	TRAINDIR=GA_$SYSTEM                         # don't need TRAINDIR AND AOD, define SOURCEDIR right away ?      
	AOD=
    else
	TRAINDIR=GA_"$SYSTEM"_AOD
	AOD=/AOD$AODNO
    fi
else # MC
    if [ "$isAOD" == "FALSE" ]; then
	TRAINDIR=GA_"$SYSTEM"_MC
	AOD=
    else
	TRAINDIR=GA_"$SYSTEM"_MC_AOD
	AOD=/AOD$AODNO
    fi
fi

for run in "${RUN[@]}"; do
    if [ "$isMC" == "FALSE" ];then # data
	SOURCEDIR=/alice/data/$YEAR/$PERIOD/000$run/$PASS$AOD/PWGGA/$TRAINDIR
    else
	SOURCEDIR=/alice/sim/$YEAR/$PERIOD/$run$AOD/PWGGA/$TRAINDIR
    fi
    if [ "$FILETYPE" == "G" ]; then
	expr=$SOURCEDIR/[[:digit:]]*_[[:digit:]]*-[[:digit:]]*/GammaConvV1_[[:digit:]]*.root
	echo "search "$SOURCEDIR;
	alien_find $SOURCEDIR/*/ GammaConvV1_*.root | grep $expr
    else
	expr=$SOURCEDIR/[[:digit:]]*_[[:digit:]]*-[[:digit:]]*/AnalysisResults.root
	echo "search "$SOURCEDIR;
	alien_find $SOURCEDIR/*/ AnalysisResults.root | grep $expr
    fi        
done

}

#******* MAIN ROUTINE ***********************************************************************************************************************

# ********************* automatic settings to according to second argument ***************************
if [ $3 = "runwise" ]; then
    echo "chose runwise download";
elif [ $3 = "merge" ]; then
    Merge
    exit
else
    echo "third argument should be runwise or merge";
fi

    #******* SELECT TRAIN **********          if not $3 = merged!
    echo "available data from ESDs: "
    isAOD=FALSE
    Show 
    echo "available data from AODs: "
    isAOD=TRUE
    Show 
    echo "which train do you want to use?";
    read TRAIN
    echo "chose train: "$TRAIN;

    # defined: BASEDIR 
    #          PERIOD YEAR isMC SYSTEM PASS RUN AODNO
    #          TRAIN
    # not usable: TRAINDIR SOURCEDIR AOD expr isAOD

    # need automatic setting: $TRAIN -> isAOD 
    echo "isAOD: TRUE or FALSE?"
    read isAOD
    echo "chose isAOD="$isAOD;


    #**************** make settings according to data/MC and ESD/AOD ***************************************
    if [ "$isMC" == "TRUE" ]; then
	TYPE=sim
	if [ "$isAOD" == "TRUE" ]; then
	    #SOURCEDIR=/alice/$TYPE/$YEAR/$PERIOD/$run/AOD$AODNO/PWGGA/GA_"$SYSTEM"_MC_AOD/$TRAIN    # MC AOD
	    PREFIX=
	    TRAINDIR=GA_"$SYSTEM"_MC_AOD
	    AOD=AOD$AODNO
	    PASS=
	elif [ "$isAOD" == "FALSE" ]; then
	    #SOURCEDIR=/alice/$TYPE/$YEAR/$PERIOD/$run/PWGGA/GA_"$SYSTEM"_MC/$TRAIN                  # MC ESD
	    PREFIX=
	    TRAINDIR=GA_"$SYSTEM"_MC
	    PASS=
	    AOD=
	else
	    echo "isAOD must be TRUE or FALSE. Is "$isAOD;
	fi
    elif [ "$isMC" == "FALSE" ]; then   
	TYPE=data
	if [ "$isAOD" == "TRUE" ]; then
	    #SOURCEDIR=/alice/$TYPE/$YEAR/$PERIOD/000$run/$PASS/AOD$AODNO/PWGGA/GA_"$SYSTEM"_AOD/$TRAIN       # data AOD
	    PREFIX=000
	    TRAINDIR=GA_"$SYSTEM"_AOD
	    AOD=AOD$AODNO
	elif [ "$isAOD" == "FALSE" ]; then
	    #SOURCEDIR=/alice/$TYPE/$YEAR/$PERIOD/000$run/$PASS/PWGGA/GA_"$SYSTEM"/$TRAIN                     # data ESD
	    PREFIX=000
	    TRAINDIR=GA_"$SYSTEM"
	    AOD=
	else
	    echo "isAOD must be TRUE or FALSE. Is "$isAOD;
	fi
    else
	echo "isMC mus be TRUE or FALSE. Is "$isMC;
    fi

    # defined: BASEDIR 
    #          PERIOD YEAR isMC SYSTEM PASS RUN AODNO
    #          TRAIN isAOD TRAINDIR PREFIX AOD
    # not usable: SOURCEDIR expr 

    #check if there is already data for this period        if not $3 = merged!!
    OldData

    #***** RUNWISE ************* download runwise and ask afterwards for merging
    # possible scenarios: no data / part of the data available already

	declare -i i=0
	for run in "${RUN[@]}"; do
	    i=i+1
	    echo "Run $i of ${#RUN[@]} ...";
	    SOURCEDIR=/alice/$TYPE/$YEAR/$PERIOD/$PREFIX$run/$PASS/$AOD/PWGGA/$TRAINDIR/$TRAIN 
	    OUTPUTDIR=$BASEDIR/$SYSTEM/$PERIOD/$run/
	    if [ $2 = "G" ]; then
		alien_ls $SOURCEDIR/GammaConvV1_*.root > filesToCopy.txt
	    else
		alien_ls $SOURCEDIR/AnalysisResults.root > filesToCopy.txt
	    fi
	    files=`cat filesToCopy.txt`
	    if [ "$files" = "" ]; then
		continue
	    fi
	    mkdir -p $OUTPUTDIR  # move below ?
	    EXISTS=`ls $OUTPUTDIR`
	    if [ "$EXISTS" != "" ] && [ "$answer" != "all" ] && [ "$answer" != "none" ]; then
	        echo "files already exist for this run: "$EXISTS" download anyways? yes / all / no / none";  # !!! Diese Abfrage funktioniert noch nicht
		read answer
		if [ "$answer" != "yes" ] && [ "$answer" != "no" ] && [ "$answer" != "all" ] && [ "$answer" != "none" ]; then
		    echo "no valid answer has been given";
		    answer=no                                                                                 
		fi
	    fi
	    for fileToCopy in $files; do
		if [ -f $OUTPUTDIR/$fileToCopy ]; then
		    echo "already copied $OUTPUTDIR_Data/$fileToCopy ";
		elif [ "$answer" == "no" ] || [ "$answer" == "none" ]; then
		    echo "don't copy "$SOURCEDIR/$fileToCopy;
		else
		    alien_cp alien:$SOURCEDIR/$fileToCopy file:$OUTPUTDIR/
		    echo "copied " $SOURCEDIR/$fileToCopy " to " $OUTPUTDIR;
		    #Change structure to standard
		    Structure
		fi
	    done  
	done
	WriteInfo
	if [ $2 = "G" ]; then
	    echo "Do you want to merge the given run list? yes / no";
	    read answerMerge
	    if [ "$answerMerge" != "yes" ] && [ "$answerMerge" != "no" ]; then
		echo "no valid answer has been given. Will assume no.";
		answerMerge=no;
	    fi
	    if [ "$answerMerge" == "yes" ]; then 
		Merge
	    fi
	fi
	
# remove help files
if [ -f fileData.txt ]; then rm fileData.txt; fi
if [ -f filesToCopy.txt ]; then rm filesToCopy.txt; fi


