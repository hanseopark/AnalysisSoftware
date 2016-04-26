#jobnumbers which need to be processed
jobnumbersFailed=`cat $1`;

#path definition
numberOfSlashesToRunnumber=$2
numberOfSlashesToStage=$3
numberOfSlashesToMergingJob=$4
# delim=`echo $path | tr -dc "/" | wc -c`

#find position of run number and subdirectory
# if [ "$3" = "data" ]; then 
# 	nSlashesToRunnumber=`expr $delim + 2`
# 	additionalPath=$4
# 	delim2=`echo $additionalPath | tr -dc "/" | wc -c`
# 	echo $delim2
# 	nSlashesToSubdirectory=`expr $delim + 2 + $delim2`
# 	echo $nSlashesToSubdirectory
# 	
# else 
# 	nSlashesToRunnumber=`expr $delim + 2`
# 	nSlashesToSubdirectory=`expr $delim + 3`
# 	additionalPath="/"
# fi 

# which kind of data are you processing
dataType="/root_archive.zip"
dataType2="root_archive.zip"
#dataType="/AliAODs.root"
maxJobnumbersFailed=`cat $1 | wc -l`;
counterJob=1;

# Loop for all jobumbers
for jobnumberFailed in $jobnumbersFailed; do
	# which job are you currently analysing
	echo "$jobnumberFailed, number $counterJob/ $maxJobnumbersFailed";
	# get the jdl from the job number
	alien_ps -jdl $jobnumberFailed > jdl.txt
	
	# take out all unwanted characters from jdl
	sed 's/{//g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/}//g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/"//g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/&//g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/\[//g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/\]//g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/ //g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/LF://g' jdl.txt > jdlmod.txt	
	mv jdlmod.txt jdl.txt 
	sed 's/,nodownload//g' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed '/LPM/d' jdl.txt > jdlmod.txt
	mv jdlmod.txt jdl.txt 
	sed 's/;//g' jdl.txt > jdlmod.txt	
	mv jdlmod.txt jdl.txt 
	sed 's/,//g' jdl.txt > jdlmod.txt	
	mv jdlmod.txt jdl.txt 
	
	outputName=`cat jdl.txt | grep "OutputDir="`
	echo $outputName
	
	dirNumber=`echo $outputName | cut -d "/" -f $numberOfSlashesToMergingJob`
	echo $dirNumber;
	runNumber=`echo $outputName | cut -d "/" -f $numberOfSlashesToRunnumber`
	echo $runNumber;
	stage=`echo $outputName | cut -d "/" -f $numberOfSlashesToStage`
	echo $stage;
	
	mkdir $runNumber/$stage/$dirNumber
	
	# find processed files 		
	cat jdl.txt | grep $dataType > fileNames.txt
	# remove duplicates
	sort fileNames.txt | uniq -u
	fileNames=`cat fileNames.txt`
# 	echo $fileNames

	fileCounter=1;
	#check each of the files in file name
	for fileName in $fileNames; do
		echo current $fileName
		echo alien:$fileName;
		# create local directory
		mkdir -p $runNumber/$stage/$dirNumber/$fileCounter
# 		rm -r $runNumber/$stage/$dirNumber/$fileCounter/*.zip
		rm -r $runNumber/$stage/$dirNumber/$fileCounter/*.root

		`alien_cp -o alien:$fileName file:$runNumber/$stage/$dirNumber/$fileCounter/$dataType2`
 
# 		cd $runNumber/$stage/$dirNumber/$fileCounter
		unzip -u $runNumber/$stage/$dirNumber/$fileCounter/$dataType2 -d $runNumber/$stage/$dirNumber/$fileCounter
# 		cd -
		fileCounter=$((fileCounter+1));
	done
	FilesToMerge=`ls $runNumber/$stage/$dirNumber/1/*.root`
	for FileToMerge in $FilesToMerge; do
		echo $FileToMerge
		FileToMerge=`echo $FileToMerge | cut -d "/" -f 5`
		echo $FileToMerge;
		hadd -f $runNumber/$stage/$dirNumber/$FileToMerge $runNumber/$stage/$dirNumber/*/$FileToMerge
	done
	
	counterJob=$((counterJob+1));
	
done

directory=`echo $outputName | cut -d "/" -f 2-$numberOfSlashesToStage`
echo $directory;

successfulJobOutputs=`alien_ls /$directory/`

maxJobnumbersSucc=`cat $successfulJobOutputs | wc -l`;
counterJob=1;

# # Loop for all jobumbers
# for successfulJobOutput in $successfulJobOutputs; do
# 	echo "$successfulJobOutput, number $counterJob/ $maxJobnumbersSucc";
# 
# 	mkdir $runNumber/$stage/$successfulJobOutput
# 	`alien_cp -o alien:/$directory/$successfulJobOutput/$dataType2 file:$runNumber/$stage/$successfulJobOutput/$dataType2`
# 	cd $runNumber/$stage/$successfulJobOutput
# 	unzip $dataType2
# 	cd -
# 	
# 	counterJob=$((counterJob+1));
# done	

FilesToMerge=`ls $runNumber/$stage/001/*.root`
for FileToMerge in $FilesToMerge; do
	echo $FileToMerge
	FileToMerge=`echo $FileToMerge | cut -d "/" -f 4`
	echo $FileToMerge;
	hadd -f $runNumber/$stage/$FileToMerge $runNumber/$stage/*/$FileToMerge
done

