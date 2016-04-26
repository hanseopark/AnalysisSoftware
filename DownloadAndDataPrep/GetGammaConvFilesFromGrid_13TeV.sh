#! /bin/bash
# creates directories
# copies files from grid
# $1 = data / MC
# $2 = runwise / merged 
# specify period, (pass,) train, run(range) below ! ( see 2x #--> )

BASEDIR=/home/ceres/danisch/GridOutputGammaConv/pp;

Structure() {

        # change structure to standard 
	ls $OUTPUTDIR/GammaConvV1_*.root > fileData.txt
	fileNumbers=`cat fileData.txt`
	for fileName in $fileNumbers; do
	    number=`echo $fileName | cut -d "_" -f 3 | cut -d "." -f 1`
	    echo "trainconfig: " $number
	    root -l -b -q -x ChangeStructureToStandard.C\(\"$OUTPUTDIR/GammaConvV1_$number.root\"\,\"$OUTPUTDIR/GammaConvV1.root\"\,\"GammaConvV1_$number\"\)
	done
}

# -------------- data -------------------------------------------------------------

if [ $1 = "data" ]; then
#-->
    #PERIOD=LHC15f;
    #PASS=pass1;
    #TRAIN=812_20150727-1717;       # output not runwise        

    #PERIOD=LHC15f;
    #PASS=pass2;
    #TRAIN=1008_20150924-1808;

    #PERIOD=LHC15f;                  # output not merged
    #PASS=pass2;
    #TRAIN=1023_20151010-1141;
    #RUN=(225000 225011 225016 225026 225031 225035 225037 225041 225043 225050 225051 225052 225106 225305 225307 225309 225310 225313 225314 225315 225322 225576 225578 225579 225580 225582 225586 225587 225705 225707 225708 225709 225710 225716 225717 225719 225753 225757 225762 225763 225766 225768 226062 226170 226220 226225 226444 226445 226452 226466 226468 226472 226476 226483 226495 226500);    

    #PERIOD=LHC15h;                
    #PASS=pass1;
    #TRAIN=1045_20151026-1652;
    #RUN=(234050 234049 234048 234039 234031 233978 233977 233976 233975 233974 233973 233972 233971 233969 233912);

    #PERIOD=LHC15f;                   # output not merged           
    #PASS=pass2;
    #TRAIN=1044_20151026-1550;  
    #below 56-run-list
    #RUN=(225000 225011 225016 225026 225031 225035 225037 225041 225043 225050 225051 225052 225106 225305 225307 225309 225310 225313 225314 225315 225322 225576 225578 225579 225580 225582 225586 225587 225705 225707 225708 225709 225710 225716 225717 225719 225753 225757 225762 225763 225766 225768 226062 226170 226220 226225 226444 226445 226452 226466 226468 226472 226476 226483 226495 226500);    

    #PERIOD=LHC15f;                 
    #PASS=pass2;
    #TRAIN=1061_20151103-0940; #61 runs below
    #RUN=(225000 225011 225016 225026 225031 225035 225037 225041 225043 225050 225051 225052 225105 225106 225305 225307 225309 225310 225313 225314 225315 225322 225576 225578 225579 225580 225582 225586 225587 225589 225609 225611 225705 225707 225708 225709 225710 225716 225717 225719 225753 225757 225762 225763 225766 225768 226062 226170 226175 226220 226225 226444 226445 226452 226466 226468 226472 226476 226483 226495 226500);    

    #PERIOD=LHC15h;                 
    #PASS=pass1;
    #TRAIN=1062_20151103-0937; #62 runs below
    #RUN=(234050 234049 234048 234045 234043 234040 234039 234031 233978 233977 233976 233975 233974 233973 233972 233971 233969 233912 233858 233837 233830 233828 233799 233743 233721 233720 233719 233718 233716 233700 233698 233697 233696 233692 233686 233678 233627 233623 233621 233614 233613 233472 233465 233361 233242 233239 233232 233217 233172 233169 233144 233120 233116 233093 233061 233059 233020 232995 232993 232986 232916 232914);
    
    PERIOD=LHC15i;                 
    PASS=pass1;
    TRAIN=1063_20151103-0937 ; #107 runs below
    RUN=(236866 236863 236862 236860 236850 236848 236835 236824 236822 236821 236816 236815 236814 236813 236569 236565 236564 236563 236562 236558 236556 236444 236443 236441  236393 236389 236360 236359 236357 236356 236354 236353 236352 236349 236348 236337 236334 236331 236284 236281 236248 236246 236244 236242 236240 236238 236234 236227 236222 236204 236203 236164 236163 236161 236158 236153 236151 236150 236138 236137 236062 235898 235897 235896 235895 235893 235891 235890 235889 235888 235886 235841 235839 235811 235759 235721 235714 235710 235694 235684 235683 235573 235547 235462 235459 235454 235450 235443 235436 235435 235432 235423 235383 235380 235364 235362 235356 235347 235345 235344 235245 235242 235226 235204 235203 235201 235196);
    

 # -- merged -------------------
     if [ $2 = "merged" ]; then
	OUTPUTDIR=$BASEDIR/$PERIOD/$PASS/$TRAIN/merged/
	SOURCEDIR=/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$TRAIN/merge
	mkdir -p $OUTPUTDIR                                                                                         
	alien_ls $SOURCEDIR/Gam* > filesToCopy.txt
	files=`cat filesToCopy.txt`
	for fileToCopy in $files; do
	    if [ -f $OUTPUTDIR/$fileToCopy ]; then
		echo "already copied $OUTPUTDIR_Data/$fileToCopy ";
	    else 
		alien_cp alien:$SOURCEDIR/$fileToCopy file:$OUTPUTDIR/  
		echo "copied " $SOURCEDIR/$fileToCopy " to " $OUTPUTDIR
	    fi
	done
        Structure                                                                      #change structure to standard

# -- runwise ----------------------
     elif [ $2 = "runwise" ]; then
	declare -i i=0
	if [ "${#RUN[@]}" -gt 1 ]; then
	    echo "${#RUN[@]} runs given. Do you want to merge the runs? yes or no?";
	    read answer
	    if [ "$answer" != "yes" ] && [ "$answer" != "no" ]; then 
		echo "no valid answer has been given"
		exit
	    fi
	fi
	for run in "${RUN[@]}"; do
	    i=i+1
	    echo "Run $i of ${#RUN[@]} ...";
	    OUTPUTDIR=$BASEDIR/$PERIOD/$PASS/$TRAIN/$run/
	    SOURCEDIR=/alice/data/2015/$PERIOD/000$run/$PASS/PWGGA/GA_pp/$TRAIN
	    mkdir -p $OUTPUTDIR                                                                                       
	    alien_ls $SOURCEDIR/Gam* > filesToCopy.txt
	    files=`cat filesToCopy.txt`
	    for fileToCopy in $files; do
		if [ -f $OUTPUTDIR/$fileToCopy ]; then
		    echo "already copied $OUTPUTDIR_Data/$fileToCopy ";
		else 
		    alien_cp alien:$SOURCEDIR/$fileToCopy file:$OUTPUTDIR/         # copy files to directory
		    echo "copied " $SOURCEDIR/$fileToCopy " to " $OUTPUTDIR
		fi
	    done
            Structure                                                              # change structure to standard 
	    if [ "$answer" = "yes" ]; then       
		echo "merge files";
		if [ $i = 1 ]; then
		    mkdir -p $BASEDIR/$PERIOD/$PASS/$TRAIN/mymerge/${#RUN[@]}runs/
		    cp $BASEDIR/$PERIOD/$PASS/$TRAIN/$run/GammaConvV1.root $BASEDIR/$PERIOD/$PASS/$TRAIN/mymerge/${#RUN[@]}runs/GammaConvV1.root 
		else
		    hadd -f ./GammaConvV1.root $BASEDIR/$PERIOD/$PASS/$TRAIN/mymerge/${#RUN[@]}runs/GammaConvV1.root $BASEDIR/$PERIOD/$PASS/$TRAIN/$run/GammaConvV1.root    # merge every file to mymerge/GammaConv/GammaConvV1.root
		    mv GammaConvV1.root $BASEDIR/$PERIOD/$PASS/$TRAIN/mymerge/${#RUN[@]}runs/
		fi
	    fi    
	done
     else
	echo "second argument should be runwise or merged";
     fi




# ----- MC -----------------------------------------------------------------------------------

elif [ $1 = "MC" ]; then
#-->

    #PERIOD=LHC15g3;
    #TRAIN=872_20150814-1240;
    #RUN=(226062 225717 225106); 

    #PERIOD=LHC15g3c2;
    #TRAIN=1094_20150924-1802;
    #RUN=(226062);

    #PERIOD=LHC15g3c2;
    #TRAIN=1219_20151010-1140;
    #RUN=(226062);

    PERIOD=LHC15g3;
    TRAIN=1255_20151026-1612;
    RUN=(225310 225309 225307 225305 225106 225052 225051 225050 225043 225041 225037 225035 225031 225026 225016 225011 225000);

# -- merged ------------------
    if [ $2 = "merged" ]; then
	OUTPUTDIR=$BASEDIR/$PERIOD/$TRAIN/merged/
	SOURCEDIR=/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$TRAIN/merge
	mkdir -p $OUTPUTDIR                                                                                         
	alien_ls $SOURCEDIR/Gam* > filesToCopy.txt
	files=`cat filesToCopy.txt`
	for fileToCopy in $files; do
	    if [ -f $OUTPUTDIR/$fileToCopy ]; then
		echo "already copied $OUTPUTDIR_Data/$fileToCopy ";
	    else 
		alien_cp alien:$SOURCEDIR/$fileToCopy file:$OUTPUTDIR/       
		echo "copied " $SOURCEDIR/$fileToCopy " to " $OUTPUTDIR
	    fi
	done

        Structure                                                               # change structure to standard 


# -- runwise ---------------------
    elif [ $2 = "runwise" ]; then
	declare -i i=0
	if [ "${#RUN[@]}" -gt 1 ]; then
	    echo "${#RUN[@]} runs given. Do you want to merge the runs? yes or no?";
	    read answer
	    if [ "$answer" != "yes" ] && [ "$answer" != "no" ]; then 
		echo "no valid answer has been given"
		exit
	    fi
	fi
	for run in "${RUN[@]}"; do
	    i=i+1
	    echo "Run $i of ${#RUN[@]} ...";
	    OUTPUTDIR=$BASEDIR/$PERIOD/$TRAIN/$run/
	    SOURCEDIR=/alice/sim/2015/$PERIOD/$run/PWGGA/GA_pp_MC/$TRAIN
	    mkdir -p $OUTPUTDIR                                                                                         
	    alien_ls $SOURCEDIR/Gam* > filesToCopy.txt
	    files=`cat filesToCopy.txt`
	    for fileToCopy in $files; do
		if [ -f $OUTPUTDIR/$fileToCopy ]; then
		    echo "already copied $OUTPUTDIR_Data/$fileToCopy ";
		else 
		    alien_cp alien:$SOURCEDIR/$fileToCopy file:$OUTPUTDIR/       
		    echo "copied " $SOURCEDIR/$fileToCopy " to " $OUTPUTDIR
		fi
	    done
            Structure  
	    if [ "$answer" = "yes" ]; then       
		echo "merge files";
		if [ $i = 1 ]; then
		    mkdir -p $BASEDIR/$PERIOD/$TRAIN/mymerge/${#RUN[@]}runs/
		    cp $BASEDIR/$PERIOD/$TRAIN/$run/GammaConvV1.root $BASEDIR/$PERIOD/$TRAIN/mymerge/${#RUN[@]}runs/GammaConvV1.root 
		else
		    hadd -f ./GammaConvV1.root $BASEDIR/$PERIOD/$TRAIN/mymerge/${#RUN[@]}runs/GammaConvV1.root $BASEDIR/$PERIOD/$TRAIN/$run/GammaConvV1.root    # merge every file to mymerge/GammaConv/GammaConvV1.root
		    mv GammaConvV1.root $BASEDIR/$PERIOD/$TRAIN/mymerge/${#RUN[@]}runs/
		fi
	    fi                                                                       # change structure to standard 
	done
    else
	echo "secon argument should be runwise or merged";
    fi

else 
    echo "first argument must be data or MC";

fi

if [ -f fileData.txt ]; then rm fileData.txt; fi
if [ -f filesToCopy.txt ]; then rm filesToCopy.txt; fi


