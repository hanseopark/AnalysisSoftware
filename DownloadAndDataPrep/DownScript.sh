#! /bin/bash
debug=0
# This Script is intended to automize the download of train outputs

# Version: V3.3
echo  -e "\e[36m+++++++++++++++++++++++++++++++++++++\e[0m"
echo "DownScript.sh Version: V3.4"

# Author: Adrian Mechler (mechler@ikf.uni-frankfurt.de)

# to use this you have to download 2 packages html-xml-utils and w3m
# additional you need an alien token
# and you have to enter the path to your certificat at line ~120

# Use: bash DownScript.sh [TrainPage] [Trainnumber] -[OPTIONS]

# actually all inputs are also asked interactive in case you forget something

# 	|-> TrainPage (e.g. GA_pp_AOD)
# 	|-> TrainNumber (e.g 768)
#                            YOU CAN SELECT multiple trains at once e.g. "768 775"

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
	BASEDIR="/media/adrian/Elements/grid_data"
	if [[ $I_AM_ALIDOCK = "1" ]]
	then
		BASEDIR="/home/alidock/media/adrian/Elements/grid_data"
	fi
	UserName="Adrian Mechler";
	pathtocert="~/.globus"
	if [[ $I_AM_ALIDOCK = "1" ]]
	then
		pathtocert="~/.globus"
	fi
	alienUserName="amechler"
	key="key.pem"
	cacert="ca.pem"
	clientca="client.pem"
else
	echo
	echo -e "This feature is not supported for user:$thisuser."
	echo -e "You have to add your certificat path to folder."
	echo
	exit
fi

# cd $BASEDIR


#####################################
#####################################



# Fuction to check if file is up to date and download it if not
function GetFile_check()
{
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then
	echo -e "\e[33mGetFile(): \e[0m alien:/$1  file:/$2"
	fi
	if [ -f $2 ]
	then
		# Check if older than one day (or not)
		if [ "$(( $(date +"%s") - $(stat -c "%Y" $2) ))" -gt "86400" ]
		then
			echo -e "\e[33m|-> \e[0mFile $2 exists. \e[33m|-> \e[0mbut is older then 1 day! download new one from alien:/$1"
			alien_cp alien:$1 file:$2 | tee -a $3
			touch "$2";
		else
			echo -e "\e[33m|-> \e[0mFile $2 exists. \e[33m|-> \e[0mand is up-to-date"
		fi
	else
		echo -e "\e[33m|-> \e[0mDownloading file alien:$1"
		alien_cp alien:$1 file:$2 | tee -a $3
		touch "$2";
	fi
}


function GetWebpage()
{
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then
		echo -e "\e[33mGetWebpage(): \e[0m\t$1\t$2\t$3"
	fi
	certpath="$3"
	Webpage="$2"

	# we download the html to grep needed information from there
	if [ -f $1 ]
	then # Check if older than one day (or not)
	if [[ `wc -l $1` < 100 ]]
	then
		rm $1
	fi
	if [ "$(( $(date +"%s") - $(stat -c "%Y" $1) ))" -gt "86400" ]
	then
		echo -e "\e[33m|-> \e[0mFile $1 exists. \e[33m|-> \e[0mbut is older then 1 Day! download new one from alien:/$1"
		# cmd="curl -s '${Webpage}' --key $certpath/$key --cacert $certpath/$cacert --cert $certpath/$clientca &> $1"
		cmd="curl '${Webpage}' --key $certpath/$key -k --cert $certpath/$clientca &> $1"
		if [[ $debug = 1 ]] || [[ $debug = 2 ]]
		then
			echo $cmd
		fi
		eval $cmd
		touch "$1";
	else
		echo -e "\e[33m|-> \e[0mFile $1 exists. \e[33m|-> \e[0mand is up-to-date"
	fi
	else
		echo -e "\e[33m|-> \e[0mDownloading ${Webpage}"
		# cmd="curl '${Webpage}' --key $certpath/$key --cacert $certpath/$cacert --cert $certpath/$clientca &> $1"
		cmd="curl '${Webpage}' --key $certpath/$key -k --cert $certpath/$clientca &> $1"
		if [[ $debug = 1 ]] || [[ $debug = 2 ]]
		then
			echo $cmd
		fi
		eval $cmd
		touch "$1";
	fi
}
# # if this is not working try
# openssl pkcs12 -in myCertificate.p12 -out newcert.pem
# openssl pkcs12 -in myCertificate.p12 -out $cacert -cacerts -nokeys
# openssl pkcs12 -in myCertificate.p12 -out $clientca -clcerts -nokeys
# openssl pkcs12 -in myCertificate.p12 -out $key -nocerts



# Fuction to check if file is there and download it if not
function GetFile()
{
	if [[ $debug = 1 ]] || [[ $debug = 2 ]]
	then
		echo -e "\e[33mGetFile(): \e[0m alien:/$1  file:/$2"
	fi

	if [[ -f $2 ]] # Check if file there
	then
		echo -e "\e[33m|-> \e[0mFile $2 exists."
	else
		echo -e "\e[33m|-> \e[0mDownloading file alien:/$1"
		alien_cp -o  alien:/$1 file:/$2 | tee -a $3
	fi
}


#####################################
#####################################

echo;echo;

echo  -e "\e[36m------------------------------------\e[0m"
echo -e "Using settings for: $UserName"
echo  -e "\e[36m------------------------------------\e[0m"

# chack if valid token is active
if [[ `alien-token-info` =~ "No Token found!" ]]
then
	echo -e "\e[31mWARNING:\e[0m No alien token"
	alien-token-init $alienUserName
	echo  -e "\e[36m------------------------------------\e[0m"
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
SetOutName=""
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
		SetOutName=${setting#*-Name_}
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
	echo -e "\e[33m|-> \e[0m TrainPage not set! What shall be used? GA_pp, GA_pp_MC, GA_pp_AOD, GA_pp_MC_AOD"
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
	echo -e "\e[31m|-> $TrainPage\e[0m not supported"
fi
if [[ $TrainNumber = "" ]]
then
	echo -e "\e[33m|-> \e[0m TrainNumber not set! What shall be used? "
	read TrainNumber
fi

##################
# Print Settings
OutName=.$TrainPage-$TrainNumber
if [[ $SetOutName = "" ]]
then
	# SetOutName=$OutName
	# echo -e "\e[33m|-> \e[0m OutName not set! using TrainNumber: OutName = $OutName"
	echo -e "\e[33m|-> \e[0m OutName = $OutName"
else
	echo -e "\e[33m|-> \e[0m OutName = $SetOutName"
fi
if [[ $MergeTrains = 0 ]]
then
	echo -e "\e[33m|-> \e[0m TrainNumber = $TrainNumber"
elif [[ $MergeTrains = 1 ]]
then
	echo -e "\e[33m|-> \e[0m TrainNumber set: $TrainNumber"
	MergeTrainsOutname=$SetOutName
fi
echo -e "\e[33m|-> \e[0m TrainPage = $TrainPage"
echo -e "\e[33m|-> \e[0m Search = $Search"
echo -e "\e[33m|-> \e[0m RunlistName = $OptRunlistName"
echo -e "\e[33m|-> \e[0m Usechildsareperiods = $Usechildsareperiods"
if [[ $DoDown = 0 ]]
then
	echo -e "\e[33m|-> \e[0m Downloading disabled"
else
	if [[ $UseMerge = 1 ]]
	then
		echo -e "\e[33m|-> \e[0m Merge files"
	fi
	if [[ $Usechildsareperiods = 1 ]]
	then
		echo -e "\e[33m|-> \e[0m Childs are periods"
	fi
fi
echo  -e "\e[36m------------------------------------\e[0m"
echo -e "\e[33mIs this Setup correct?\e[0m (y/n)"
printf "Answer: "
read -e setupcorrect
if [ ! "$setupcorrect" = "y" ]
then
  echo ""
  echo "please try again   :-("
  echo ""
  exit
fi
echo  -e "\e[36m------------------------------------\e[0m"
echo ""
echo ""

if [[ $DoDown = 1 ]]
then
	for TrainNumber in `cat $TrainNumberFile`
	do
		echo;echo;echo;echo;
		echo  -e "\e[36m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\e[0m"
		OutName=.$TrainPage-$TrainNumber
		echo -e "\e[36m |-> \e[0m OutName = $OutName\e[36m | \e[0mTrainNumber = $TrainNumber\e[36m | \e[0mTrainPage = $TrainPage\e[36m | \e[0mSearch = $Search\e[36m | \e[0mRunlistName = $OptRunlistName\e[36m | \e[0mUsechildsareperiods = $Usechildsareperiods\e[36m | \e[0mUseMerge = $UseMerge"
		echo  -e "\e[36m~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\e[0m"
		echo;echo;

		#####################################
		#####################################


		AlienDir="/alice/cern.ch/user/a/alitrain/PWGGA/$TrainPage/"
		OUTPUTDIR=$BASEDIR/$OutName
		List="listGrid.txt"
		GlobalVariablesFile="globalvariables.C"
		envFile="env.sh"
		FileList="$OUTPUTDIR/FileList.txt"
		TrainPageHTML="$BASEDIR/$TrainPage.html"
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
		echo  -e "\e[36m------------------------------------\e[0m"
		echo

		# get all childs involved in train run
		alien_ls $AlienDir | grep $TrainNumber\_2 &> $List


		# get one env.sh to get basic information about the train
		Childeone=`sed '2q;d' $List`
		echo -e "\e[33m|->\e[0m $Childeone"
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
			echo  -e "\e[35m=====================================\e[0m"
			echo
			echo -e "\e[35m $child \e[0m"

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
				GetFile "$AlienDir$child/$GlobalVariablesFile" "$GlobalVariablesPath" /dev/null
			else
				GetFile "$AlienDir$child/$GlobalVariablesFile" "$GlobalVariablesPath" /dev/null &> /dev/null
			fi

			ChildName=`grep "periodName =" "$GlobalVariablesPath" | awk -F "= " '{print $2}' | awk -F ";" '{print $1}' | awk -F "\"" '{print $2}'`
			echo -e " periodName = $ChildName"


			#####################################
			#####################################


			if [[ `wc -l $TrainPageHTML` < 1 ]]
			then
				echo -e "\e[31mError\e[0m $TrainPageHTML not found! aborting" | tee -a $ErrorLog
				continue
			fi

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
							echo -e "\e[33m|-> \e[0m found period"
						fi
					fi
				fi
				if [[ $found = 1 ]]
				then
					if [[ $lookchildID = 1 ]]
					then
						if [[ $debug = 2 ]]
						then
							echo -e "\e[33m|-> \e[0m lookchildID $line"
						fi
						foundchild=0
						if [[ $line = "$childID" ]]
						then
							if [[ $Usechildsareperiods = 0 ]]
							then
								echo -e "\e[33m|->\e[0m ChildName = $ChildName"
							else
								echo -e "\e[33m|->\e[0m will load child/period name later"
							fi
							echo -e "\e[33m|->\e[0m childID = $childID"
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
						echo -e "\e[33m|-> \e[0m found child: $line"
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
								cmd="alien_ls $AlienDir$child/merge_runlist_$RunlistID/ | grep "$Search" | tee -a $FileList"
								echo $cmd
								eval $cmd
							done

							echo;
							echo  -e "\e[36m------------------------------------\e[0m"
							echo -e "\e[36m|-> Runlist $RunlistID\tName = $RunlistName \e[0m"

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
								echo -e "\e[33m|->\e[0m Changed periodName = $ChildName"
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
								echo -e "Found files:\e[33m|->\e[0m Filename = $filename"
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
										echo -e "\e[33m|-> \e[0m$outFile already merged"
									else
										if [[ -f $logFile ]]
										then
											rm $logFile
										fi
										if [[ -f $outFile ]]
										then
											if [[ -f $mergedFile ]]
											then
												echo -e "\e[33m|->\e[0m merging $outFile"  #(log: $logFile)"
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
													echo -e "\e[31mError\e[0m $outFile not merged correctly" | tee -a $ErrorLog
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
												echo -e "\e[33m|->\e[0m Copy: $outFile is first"
												echo "Copyed to $mergedFile" > $logFile
												cp $outFile $mergedFile
												touch $alreadyMerged
											fi
										else
											echo -e "\e[31m|->\e[0m missing $outFile  "
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
										runDir=$OUTPUTDIR/$RunlistName/$ChildName/$runName
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


if [[ $NTrainFile = 1 ]] #&& [[ ! $TrainNumber = *"+"* ]]
then
	cmd="ln -sf $BASEDIR/$OutName $BASEDIR/$SetOutName"
	echo $cmd
	eval $cmd
fi


if [[ $MergeTrains = 1 ]]
then
	echo  -e "\e[36m------------------------------------\e[0m"
	echo "Start merging trains: $MergeTrainsOutname"
	echo  -e "\e[36m------------------------------------\e[0m"
	for TrainNumber in `cat $TrainNumberFile`
	do
		echo;echo;
		echo  -e "\e[36m|-> \e[0m $TrainNumber"
		singletrainDir=.$TrainPage-$TrainNumber
		for RunlistName in `cat $OptRunlistNamefile`
		do
			echo  -e "\t\e[36m|-> \e[0m $RunlistName"
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
				echo  -e "\t\t\e[36m|-> \e[0m $filename"
				outFile=$singletrainDir/$RunlistName/$filename
				alreadyMerged=$MergeTrainsOutname/$RunlistName/.$TrainNumber-${filename%%.root}.merged
				logFile=$MergeTrainsOutname/$RunlistName/$TrainNumber-${filename%%.root}.log
				mergedFile=$MergeTrainsOutname/$RunlistName/$filename
				if [[ -f $alreadyMerged ]]
				then
					echo -e "\t\t\e[33m|-> \e[0m$outFile already merged"
				else
					if [[ -f $logFile ]]
					then
						rm $logFile
					fi
					if [[ -f $outFile ]]
					then
						if [[ -f $mergedFile ]]
						then
							echo -e "\t\t\e[33m|->\e[0m merging $outFile"  #(log: $logFile)"
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
								echo -e "\t\t\e[31mError\e[0m $outFile not merged correctly" | tee -a $ErrorLog
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
							echo -e "\t\t\e[33m|->\e[0m Copy: $outFile is first"
							echo "Copyed to $mergedFile" > $logFile
							cp $outFile $mergedFile
							touch $alreadyMerged
						fi
					else
						echo -e "\t\t\t\e[31m|->\e[0m missing $outFile  "
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
