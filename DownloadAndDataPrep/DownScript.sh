#! /bin/bash

# This Script is intended to automize the download of train outputs

# Version: V3
echo  -e "\e[1;36m+++++++++++++++++++++++++++++++++++++\e[21;39m"
echo "DownScript.sh Version: V3"

# Author: Adrian Mechler (mechler@ikf.uni-frankfurt.de)

# to use this you have to download 2 packages html-xml-utils and w3m
# additional you need an alien token 
# and you have to enter the path to your certificat at line ~120

# Use: bash DownScript.sh [TrainPage] [Trainnumber] -[OPTIONS]

# actually all input is also asked interactive in case you forget something

# 	|-> TrainPage (e.g. GA_pp_AOD)
# 	|-> TrainNumber (e.g +768)
#                            YOU CAN SELECT multiple trains at once e.g. "+768 +775"

#####################################

# OPTIONS:
# 	|-> OutName (is optional)
# 	|-> Search (search for string in filename to be downloaded e.g GammaCalo, default is .root)
#							 YOU CAN SELECT multiple at once e.g. "?_GammaCalo ?_GammaConv"
# 	|-> OptRunlistName (select Runlist to be downloaded e.g. listDPGCalo, default is all availible)
#                            YOU CAN SELECT multiple at once e.g. "-RL_listDPGCalo -RL_listDPGEDC"
# 	|-> mergechilds (merges all childs)
# 	|-> runwise (download all runs)
# 	|-> childsareperiods

# Examlple: bash DownScript_V2.sh 776 GA_pp_AOD -GammaConvCalo -RL_listDPGCalo -RL_listDPGEDC -mergechilds -childsareperiods -runwise

#####################################
#####################################

# we need a certificat to get info from alimonitor 
pathtocert=""
UserName=""
alienUserName=""
BASEDIR=${PWD}
thisuser=`echo ${USER}`
if [[ $thisuser = "adrian" || $thisuser = "amechler" ]] 
then
	BASEDIR=${PWD}
	UserName="Adrian Mechler";
	pathtocert="~/.globus"
	alienUserName="amechler"
else
	echo
	echo -e "This feature is not supported for user:$thisuser."
	echo -e "You have to add your certificat path to folder."
	echo
	exit
fi


#####################################
#####################################



# Fuction to check if file is up to date and download it if not
function GetFile_check()
{
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then 
	echo -e "\e[1;33mGetFile(): \e[21;39m alien:/$1  file:/$2" 
	fi
	if [ -f $2 ]
	then
		# Check if older than one day (or not)
		if [ "$(( $(date +"%s") - $(stat -c "%Y" $2) ))" -gt "86400" ]
		then
			echo -e "\e[1;33m|-> \e[21;39mFile $2 exists. \e[1;33m|-> \e[21;39mbut is older then 1 day! download new one from alien:/$1"
			alien_cp alien:$1 file:$2 | tee -a $3
			touch "$2";
		else
			echo -e "\e[1;33m|-> \e[21;39mFile $2 exists. \e[1;33m|-> \e[21;39mand is up-to-date"
		fi
	else
		echo -e "\e[1;33m|-> \e[21;39mDownloading file alien:$1"
		alien_cp alien:$1 file:$2 | tee -a $3
		touch "$2";
	fi
}


function GetWebpage()
{
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then 
		echo -e "\e[1;33mGetWebpage(): \e[21;39m\t$1\t$2\t$3" 
	fi
	certpath="$3"
	Webpage="$2"
	# we download the html to grep needed information from there
	if [ -f $1 ]
	then # Check if older than one day (or not)
if [ "$(( $(date +"%s") - $(stat -c "%Y" $1) ))" -gt "86400" ]
then
	echo -e "\e[1;33m|-> \e[21;39mFile $1 exists. \e[1;33m|-> \e[21;39mbut is older then 1 day! download new one from alien:/$1"
	cmd="curl -s '${Webpage}' --key $certpath/key.pem --cacert $certpath/ca.pem --cert $certpath/client.pem &> $1"
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then 
		echo $cmd 
	fi
	eval $cmd
	touch "$1";
else
	echo -e "\e[1;33m|-> \e[21;39mFile $1 exists. \e[1;33m|-> \e[21;39mand is up-to-date"
fi
else
	echo -e "\e[1;33m|-> \e[21;39mDownloading ${Webpage}"
	cmd="curl -s '${Webpage}' --key $certpath/key.pem --cacert $certpath/ca.pem --cert $certpath/client.pem &> $1"
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then 
		echo $cmd 
	fi
	eval $cmd
	touch "$1";
fi
}

# Fuction to check if file is there and download it if not
function GetFile()
{ 
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then 
		echo -e "\e[1;33mGetFile(): \e[21;39m alien:/$1  file:/$2" 
	fi

	if [[ -f $2 ]] # Check if file there
	then
		echo -e "\e[1;33m|-> \e[21;39mFile $2 exists."
	else
		echo -e "\e[1;33m|-> \e[21;39mDownloading file alien:/$1"
		alien_cp -o  alien:/$1 file:/$2 | tee -a $3
	fi
}


#####################################
#####################################

echo;echo;

echo  -e "\e[1;36m------------------------------------\e[21;39m" 
echo -e "Using settings for: $UserName"
echo  -e "\e[1;36m------------------------------------\e[21;39m"

# chack if valid token is active
if [[ `alien-token-info` =~ "No Token found!" ]]
then 
	echo -e "\e[1;31mWARNING:\e[21;39m No alien token"
	alien-token-init $alienUserName
	echo  -e "\e[1;36m------------------------------------\e[21;39m"
	echo;echo;
fi


################
#There are Settings to be used beforehand
TrainPage=""
TrainPageNum=""
TrainNumber=""
MergeTrains=0
DoDown=1
MergeTrainsOutname=0
TrainNumberFile="TrainNumberFile.txt"
OutName=""
Search=".root"
Searchfile="Searchfile.txt"
OptRunlistName=list
OptRunlistNameSet=0
OptRunlistNamefile=OptRunlistNames.txt
OptAllRunlists=0
UseMerge=0
Userunwise=0
Usechildsareperiods=0
re='^[0-9]+$'

if [[ -f $TrainNumberFile ]]
then
	rm $TrainNumberFile
fi
if [[ -f $Searchfile ]]
then
	rm $Searchfile
fi
if [[ -f $OptRunlistNamefile ]]
then
	rm $OptRunlistNamefile
fi

for setting in "$@"
do
	# echo $setting
	if [[ $setting = "GA_"* ]]
	then
		TrainPage="$setting"
	elif [[ $setting =~ $re ]]
	then
		TrainNumbertmp=${setting}
		echo $TrainNumbertmp >> $TrainNumberFile
		if [[ $TrainNumber = "" ]]
		then
			TrainNumber=$TrainNumbertmp
			MergeTrains=0
		else
			TrainNumber=$TrainNumber'+'$TrainNumbertmp
			MergeTrains=1
		fi
	elif [[ $setting = "-Name_"* ]]
	then
		OutName=${setting#*-Name_}
	elif [[ $setting = "?_"* ]]
	then
		Searchtmp=${setting#*\?_}
		echo $Searchtmp >> $Searchfile
		if [[ $Search = ".root" ]]
		then
			Search=*$Searchtmp*
		else
			Search=$Search'|'*$Searchtmp*
		fi
	elif [[ $setting = "-RL_"* ]]
	then
		OptRunlistNametmp=${setting#*-RL_}
		echo $OptRunlistNametmp >> $OptRunlistNamefile
		if [[ $OptRunlistName = "list" ]]
		then
			OptRunlistName=*$OptRunlistNametmp*
		else
			OptRunlistName=$OptRunlistName'|'*$OptRunlistNametmp*
		fi
	elif [[ $setting = "-mergechilds" ]]
	then
		UseMerge=1
	elif [[ $setting = "-runwise" ]]
	then
		Userunwise=1
	elif [[ $setting = "-childsareperiods" ]]
	then
		Usechildsareperiods=1
	elif [[ $setting = "-noDown" ]]
	then
		DoDown=0
	fi
done

if [[ ! -f $OptRunlistNamefile ]]
then
	OptAllRunlists=1
fi
################
# Required Settings can also be Set interactive
if [[ $TrainPage = "" ]]
then
	echo -e "\e[1;33m|-> \e[21;39m TrainPage not set! What shall be used? GA_pp, GA_pp_MC, GA_pp_AOD, GA_pp_MC_AOD"
	read TrainPage
fi
if [[ $TrainPage = "GA_pp" ]]
then
	TrainPageNum="23"
elif [[ $TrainPage = "GA_pp_MC" ]]
then
	TrainPageNum="43"
elif [[ $TrainPage = "GA_pp_AOD" ]]
then
	TrainPageNum="62"
elif [[ $TrainPage = "GA_pp_MC_AOD" ]]
then  
	TrainPageNum="66"
else
	echo -e "\e[1;31m|-> $TrainPage\e[21;39m not supported"
fi
if [[ $TrainNumber = "" ]]
then
	echo -e "\e[1;33m|-> \e[21;39m TrainNumber not set! What shall be used? "
	read TrainNumber
fi

##################
# Print Settings
if [[ $OutName = "" ]]
then
	OutName=$TrainPage-$TrainNumber
	echo -e "\e[1;33m|-> \e[21;39m OutName not set! using TrainNumber: OutName = $OutName"
else
	echo -e "\e[1;33m|-> \e[21;39m OutName = $OutName"
fi
if [[ $MergeTrains = 0 ]]
then
	echo -e "\e[1;33m|-> \e[21;39m TrainNumber = $TrainNumber"
elif [[ $MergeTrains = 1 ]]
then
	echo -e "\e[1;33m|-> \e[21;39m TrainNumber set: $TrainNumber"
	MergeTrainsOutname=$TrainPage-$TrainNumber
fi
echo -e "\e[1;33m|-> \e[21;39m TrainPage = $TrainPage"
echo -e "\e[1;33m|-> \e[21;39m Search = $Search"
echo -e "\e[1;33m|-> \e[21;39m RunlistName = $OptRunlistName"
echo -e "\e[1;33m|-> \e[21;39m Usechildsareperiods = $Usechildsareperiods"
if [[ $DoDown = 0 ]]
then 
	echo -e "\e[1;33m|-> \e[21;39m Downloading disabled"
else
	if [[ $UseMerge = 1 ]]
	then 
		echo -e "\e[1;33m|-> \e[21;39m Merge files"
	fi
	if [[ $Usechildsareperiods = 1 ]]
	then
		echo -e "\e[1;33m|-> \e[21;39m Childs are periods"		
	fi
fi
echo  -e "\e[1;36m------------------------------------\e[21;39m" 
echo -e "\e[1;33mIs this Setup correct?\e[21;39m (y/n)"
printf "Answer: "
read -e setupcorrect
if [ ! "$setupcorrect" = "y" ]
then
  echo ""
  echo "please try again   :-("
  echo ""
  exit
fi
echo  -e "\e[1;36m------------------------------------\e[21;39m" 
echo ""
echo ""

if [[ $DoDown = 1 ]]
then
	for TrainNumber in `cat $TrainNumberFile`
	do
		echo;echo;echo;echo;
		echo  -e "\e[1;36m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\e[21;39m" 
		if [[ $OutName = "" ]] || [[ $MergeTrains = 1 ]]
		then
			OutName=$TrainPage-$TrainNumber
		fi
		echo -e "\e[1;36m |-> \e[21;39m OutName = $OutName\e[1;36m | \e[21;39mTrainNumber = $TrainNumber\e[1;36m | \e[21;39mTrainPage = $TrainPage\e[1;36m | \e[21;39mSearch = $Search\e[1;36m | \e[21;39mRunlistName = $OptRunlistName\e[1;36m | \e[21;39mUsechildsareperiods = $Usechildsareperiods\e[1;36m | \e[21;39mUseMerge = $UseMerge"
		echo  -e "\e[1;36m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\e[21;39m" 
		echo;echo;

		#####################################
		#####################################


		AlienDir="/alice/cern.ch/user/a/alitrain/PWGGA/$TrainPage/"
		OUTPUTDIR=$BASEDIR/$OutName
		List="listGrid.txt" 
		GlobalVariablesFile="globalvariables.C"
		envFile="env.sh"
		FileList="$OUTPUTDIR/FileList.txt"
		TrainPageHTML="$TrainPage.html"
		ErrorLog="$OUTPUTDIR/Errors.log"
		trainwebpage="https://alimonitor.cern.ch/trains/train.jsp?train_id=${TrainPageNum}"
		PeriodName=""


		if [[ -f $ErrorLog ]]
		then
			rm $ErrorLog
		fi


		# prepare directory
		mkdir -p $OUTPUTDIR


		#####################################
		#####################################

		echo 
		GetWebpage $TrainPageHTML $trainwebpage $pathtocert
		echo 
		echo  -e "\e[1;36m------------------------------------\e[21;39m" 
		echo 

		# get all childs involved in train run
		alien_ls $AlienDir | grep $TrainNumber\_ &> $List


		# get one env.sh to get basic information about the train
		Childeone=`sed '2q;d' $List`
		echo -e "\e[1;33m|->\e[21;39m $Childeone"
		if [[ $debug = 1 ]] || [[ $debug = 2 ]]
		then 
			GetFile $AlienDir$Childeone/$envFile $OUTPUTDIR/.$envFile /dev/null
		else 
			GetFile $AlienDir$Childeone/$envFile $OUTPUTDIR/.$envFile /dev/null &> /dev/null
		fi
		PERIODNAME=`grep "PERIOD_NAME=" $OUTPUTDIR/.$envFile | awk -F "=" '{print $2}' | awk -F ";" '{print $1}'`
		echo "PERIODNAME = $PERIODNAME"
		NCHilds=`grep "CHILD_DATASETS=" $OUTPUTDIR/.$envFile | awk -F "=" '{print $2}' | awk -F ";" '{print $1}'`
		echo "NCHilds = $NCHilds"
		echo 



		#####################################
		#####################################


		# Loop over all childs 
		for child in `sort ${List}`;
		do
			echo 
			echo  -e "\e[1;35m=====================================\e[21;39m" 
			echo 
			echo -e "\e[1;35m $child \e[21;39m" 

			childID=${child##*child_}

			foundperiod=0
			foundchild=0
			lookchildID=0	
			foundRunlists=0
			RunlistID=0
			ChildName=""

			# prepare directory
			mkdir -p "$OUTPUTDIR/.$child"

			# get one globalvariables.C to get information about the child
			GlobalVariablesPath="$OUTPUTDIR/.$child/$GlobalVariablesFile"
			if [[ $debug = 1 ]] || [[ $debug = 2 ]]
			then 
				GetFile "$AlienDir$child/$GlobalVariablesFile" "$GlobalVariablesPath" /dev/null &> /dev/null
			else 
				GetFile "$AlienDir$child/$GlobalVariablesFile" "$GlobalVariablesPath" /dev/null &> /dev/null
			fi

			ChildName=`grep "periodName =" "$GlobalVariablesPath" | awk -F "= " '{print $2}' | awk -F ";" '{print $1}' | awk -F "\"" '{print $2}'`
			echo -e " periodName = $ChildName"


			#####################################
			#####################################


			# parse html to get infos
			for line in `cat $TrainPageHTML`
			do
				if [[ $line = *"editPeriod"* ]]
				then 
					if [[ $found = 1 ]]
					then
						found=0
					fi
					if [[ $line = *"$PERIODNAME"* ]]
					then				
						found=1
						if [[ $debug = 2 ]]
						then 
							echo -e "\e[1;33m|-> \e[21;39m found period" 
						fi
					fi
				fi
				if [[ $found = 1 ]]
				then
					if [[ $lookchildID = 1 ]]
					then	
						if [[ $debug = 2 ]]
						then 
							echo -e "\e[1;33m|-> \e[21;39m lookchildID $line" 
						fi							
						foundchild=0
						if [[ $line = "$childID" ]]
						then 
							if [[ $Usechildsareperiods = 0 ]]
							then
								echo -e "\e[1;33m|->\e[21;39m ChildName = $ChildName"
							else
								echo -e "\e[1;33m|->\e[21;39m will load child/period name later"
							fi
							echo -e "\e[1;33m|->\e[21;39m childID = $childID"
							foundchild=1
						fi
						lookchildID=0	
					fi
					if [[ $line = *"child"* ]]
					then 
						lookchildID=1				
					fi
					if [[ $line = *"</tr>"* ]]
					then 
						foundchild=0			
					fi
				fi 
				if [[ $foundchild = 1 ]]
				then
					if [[ $debug = 2 ]]
					then 
						echo -e "\e[1;33m|-> \e[21;39m found child: $line" 
					fi	
					if [[ $foundRunlists = 1 ]]
					then
						RunlistName=${line%%:}
						RunlistID=$(( $RunlistID +1))

						if [[ $OptAllRunlists = 1 ]]
						then 
							echo "$RunlistName" > $OptRunlistNamefile
						fi

						for OptRunlistName in `cat $OptRunlistNamefile`
						do
							if [[ ! "$RunlistName" = "$OptRunlistName" ]]
							then 
								continue
							fi

							# look for availible Files
							if [ -f $FileList ]; then
								rm $FileList
							fi
							for Search in `cat $Searchfile`
							do
								alien_ls $AlienDir$child/merge_runlist_$RunlistID/ | grep "$Search" >> $FileList
							done

							echo;
							echo  -e "\e[1;36m------------------------------------\e[21;39m" 
							echo -e "\e[1;36m|-> Runlist $RunlistID\tName = $RunlistName \e[21;39m"

							# prepare directory
							mkdir -p $OUTPUTDIR/.$child/$RunlistName

							if [[ $Usechildsareperiods = 1 ]]
							then
								outFile="$OUTPUTDIR/.$child/$RunlistName/Stage_1.xml"
								inFile="$AlienDir$child/merge_runlist_$RunlistID/Stage_1.xml"
								if [[ $debug = 1 ]] || [[ $debug = 2 ]]
								then 
									GetFile "$inFile" "$outFile" /dev/null 
								else 
									GetFile "$inFile" "$outFile" /dev/null &> /dev/null
								fi
								grep "turl" $outFile | awk -F "turl" '{print $2}' | awk -F "\"" '{print $2}' > RunlistLinks.txt
								ChildName=`echo $(head -n 1 RunlistLinks.txt) | awk -F "/" '{print $7}' `
								echo -e "\e[1;33m|->\e[21;39m Changed periodName = $ChildName"
							fi

							mkdir -p $OUTPUTDIR/$ChildName 
							mkdir -p $OUTPUTDIR/$ChildName/$RunlistName

							# mv $OUTPUTDIR/$child/* $OUTPUTDIR/$ChildName/. &> /dev/null
							# rm -r $OUTPUTDIR/$child &> /dev/null

							if [[ $UseMerge = 1 ]]
							then 
								mkdir -p $OUTPUTDIR/$RunlistName
							# rm $OUTPUTDIR/$RunlistName/*
							fi
							# Download all relevant files 
							for filename in `cat $FileList`
							do
								echo
								echo -e "Found files:\e[1;33m|->\e[21;39m Filename = $filename"
								outFile=$OUTPUTDIR/$ChildName/$RunlistName/$filename
								inFile=$AlienDir$child/merge_runlist_$RunlistID/$filename
								downlogFile=$OUTPUTDIR/$ChildName/$RunlistName/${filename%%.root}.downlog
								GetFile $inFile $outFile $downlogFile

								if [[ $UseMerge = 1 ]]
								then 
									alreadyMerged=$OUTPUTDIR/$ChildName/$RunlistName/.${filename%%.root}.merged
									logFile=$OUTPUTDIR/$ChildName/$RunlistName/${filename%%.root}.log
									mergedFile=$OUTPUTDIR/$RunlistName/$filename
									if [[ -f $alreadyMerged ]]
									then
										echo -e "\e[1;33m|-> \e[21;39m$outFile already merged" 
									else
										if [[ -f $logFile ]]
										then
											rm $logFile
										fi
										if [[ -f $outFile ]]
										then
											if [[ -f $mergedFile ]]
											then
												echo -e "\e[1;33m|->\e[21;39m merging $outFile"  #(log: $logFile)" 
												hadd -k $mergedFile.tmp $mergedFile $outFile  &> $logFile
												if [[ $debug = 1 ]] || [[ $debug = 2 ]]
												then 
													echo " Log:"
													cat $logFile 
													echo;echo;
												fi
												if [[ $debug = 1 ]] || [[ $debug = 2 ]]
												then 
													echo "rm $mergedFile" 
												fi
												if [[ `wc -l $logFile` > 9 ]]
												then
													echo -e "\e[1;31mError\e[21;39m $outFile not merged correctly" | tee -a $ErrorLog
													cat $logFile >> $ErrorLog
													rm $mergedFile.tmp
												else
													rm $mergedFile  
													if [[ $debug = 1 ]] || [[ $debug = 2 ]]
													then 
														echo "mv $mergedFile.tmp $mergedFile " 
													fi
													mv $mergedFile.tmp $mergedFile
													touch $alreadyMerged
												fi
											else
												echo -e "\e[1;33m|->\e[21;39m Copy: $outFile is first"
												echo "Copyed to $mergedFile" > $logFile
												cp $outFile $mergedFile
												touch $alreadyMerged
											fi
										else
											echo -e "\e[1;31m|->\e[21;39m missing $outFile  " 
										fi
									fi
								fi
								if [[ $Userunwise = 1 ]]
								then
								# Download all runs
									for runfilename in `cat RunlistLinks.txt`
									do
										runName=`grep "$runfilename" RunlistLinks.txt | awk -F "/" '{print $8}'`
										printf "Download Run:  $runName " 
										runDir=$OUTPUTDIR/$ChildName/$RunlistName/$runName
										mkdir -p $runDir &> /dev/null
										runoutFile=$runDir/$filename
										runinFile=${runfilename%%root_archive.zip}/$filename
										rundownlogFile=$runDir/${filename%%.root}.downlog
										if [[ -f $rundownlogFile ]]
										then
											rm $rundownlogFile
										fi
										GetFile $runinFile $runoutFile $rundownlogFile #>> $rundownlogFile
									done
								fi
							done
							if [[ $debug = 1 ]] || [[ $debug = 2 ]]
							then 
								echo "skipping $RunlistName (searching for *$OptRunlistName*)" 
							fi

						done
						foundRunlists=0
					fi	
					if [[ $line = *"Runlist"* ]]
					then 					
						foundRunlists=1
					fi
					if [[ $line = *"table_row"* ]]
					then 
						RunlistID=0
					fi
				fi
			done		
		done



		rm $List
		rm $FileList
		# rm $OUTPUTDIR/$envFile
		# rm $TrainPageHTML
		echo;echo;echo;echo;

	done
fi

if [[ $MergeTrains = 1 ]]
then
	echo  -e "\e[1;36m------------------------------------\e[21;39m" 
	echo "Start merging trains: $MergeTrainsOutname"
	echo  -e "\e[1;36m------------------------------------\e[21;39m" 
	for TrainNumber in `cat $TrainNumberFile`
	do
		echo;echo;
		echo  -e "\e[1;36m|-> \e[21;39m $TrainNumber" 
		singletrainDir=$TrainPage-$TrainNumber				
		for RunlistName in `cat $OptRunlistNamefile`
		do
			echo  -e "\t\e[1;36m|-> \e[21;39m $RunlistName" 
			mkdir -p $MergeTrainsOutname/$RunlistName

			# look for availible Files
			FileList="$MergeTrainsOutname/FileList.txt"
			if [ -f $FileList ]; then
				rm $FileList
			fi
			ls $singletrainDir/$RunlistName/*.root >> $FileList

			for filenametmp in `cat $FileList`
			do
				filename=${filenametmp#*/*/}
				echo  -e "\t\t\e[1;36m|-> \e[21;39m $filename" 
				outFile=$singletrainDir/$RunlistName/$filename
				alreadyMerged=$singletrainDir/$RunlistName/.${filename%%.root}.merged
				logFile=$singletrainDir/$RunlistName/${filename%%.root}.log
				mergedFile=$MergeTrainsOutname/$filename
				if [[ -f $alreadyMerged ]]
				then
					echo -e "\t\t\e[1;33m|-> \e[21;39m$outFile already merged" 
				else
					if [[ -f $logFile ]]
					then
						rm $logFile
					fi
					if [[ -f $outFile ]]
					then
						if [[ -f $mergedFile ]]
						then
							echo -e "\t\t\e[1;33m|->\e[21;39m merging $outFile"  #(log: $logFile)" 
							hadd -k $mergedFile.tmp $mergedFile $outFile  &> $logFile
							if [[ $debug = 1 ]] || [[ $debug = 2 ]]
							then 
								echo " Log:"
								cat $logFile 
								echo;echo;
							fi
							if [[ $debug = 1 ]] || [[ $debug = 2 ]]
							then 
								echo "rm $mergedFile" 
							fi
							if [[ `wc -l $logFile` > 9 ]]
							then
								echo -e "\t\t\e[1;31mError\e[21;39m $outFile not merged correctly" | tee -a $ErrorLog
								cat $logFile >> $ErrorLog
								rm $mergedFile.tmp
							else
								rm $mergedFile  
								if [[ $debug = 1 ]] || [[ $debug = 2 ]]
								then 
									echo "mv $mergedFile.tmp $mergedFile " 
								fi
								mv $mergedFile.tmp $mergedFile
								touch $alreadyMerged
							fi
						else
							echo -e "\t\t\e[1;33m|->\e[21;39m Copy: $outFile is first"
							echo "Copyed to $mergedFile" > $logFile
							cp $outFile $mergedFile
							touch $alreadyMerged
						fi
					else
						echo -e "\t\t\t\e[1;31m|->\e[21;39m missing $outFile  " 
					fi
				fi
			done
		done
	done
fi


rm $TrainNumberFile
rm $Searchfile
rm $OptRunlistNamefile


echo;echo;echo;echo;