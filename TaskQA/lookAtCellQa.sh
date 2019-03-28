#script to make sorting detailed cell output easier
# first line cd into directory you want to sort
# --action1 " "--action2 " " ... --action9 "" map commands to the keys 1to9 (change direcotries as desired)
# last line specifies a directory to copy pictures to when using the "f" key. Even if you don't use this feature, i think
# that feh needs this argument.
#
# for more info see man feh
location=$1
location2=$2
suffix=$3

echo "*******************************************************************************************"
echo "************************************ SETTINGS *********************************************"
echo "*******************************************************************************************"
echo "The folder you are gonna analyse is: " $location
echo "The corresponding bad/ok/maybe classification will be moved to: " $location2
echo "You are looking for $suffix"
echo "*******************************************************************************************"
echo "*******************************************************************************************"
echo "*******************************************************************************************"

mkdir -p $location2/bad
mkdir -p $location2/maybe
mkdir -p $location2/ok
basedir=`pwd`
echo $basedir
NCellsToBeLookedAt=`ls $location/Cell*.$suffix | wc -l`

echo "*******************************************************************************************"
echo "In the analysis folder are: "
if [ $NCellsToBeLookedAt -gt 0 ]; then
	echo "$NCellsToBeLookedAt cells to be investigated"
else
	echo "no $suffix's were found in $location/. Exiting ..."
	exit
fi

echo "*******************************************************************************************"
echo "For each of the cells you will be asked to say whether its ok (3), maybe ok (2) or bad (1)";
echo "Are you ready to proceed? Yes/No"
read answer
if [ $answer = "No" ] || [ $answer = "N" ] || [ $answer = "no" ] || [ $answer = "n" ]; then
        echo "Come back later, if you are up to it!";
	exit;
elif [ $answer = "Yes" ] || [ $answer = "Y" ] || [ $answer = "yes" ] || [ $answer = "y" ]; then
  echo "Then lets start! Happy answering ;)";
else
        echo "your answer was not valid"
  exit;
fi

NCellsGood=`ls $location/Good/Cell*.png | wc -l`
if [ $NCellsToBeLookedAt -gt 0 ]; then
	echo "Let me show you some good examples first. Use Key X to finish viewing examples ............ "
	sleep 5s
	cd $location/Good
        feh -F --cycle-once
fi
echo "*******************************************************************************************"
echo "*******************************************************************************************"


echo "Now lets start classifying (use the 'd' key to toggle display or filename)............................................."
sleep 5s
cd ..
echo $location $basedir/$location2 $suffix
pwd
feh -F --cycle-once --action1 "mv -f '%F' $basedir/$location2/bad/%F " \
--action2  "mv -f '%F' $basedir/$location2/maybe/%F" \
--action3  "mv -f '%F' $basedir/$location2/ok/%F" \
*.$suffix

echo "*******************************************************************************************"
echo "Great you have done it! Congrats!!"
echo "*******************************************************************************************"
echo "*******************************************************************************************"
exit

