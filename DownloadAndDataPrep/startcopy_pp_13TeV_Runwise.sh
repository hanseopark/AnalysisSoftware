#!/bin/bash

# script to download runwise output
# $1 = PERIOD
# $2 MC or Data
# needs GridJobFileListDL

##################################################################
##############             DATA              #####################
##################################################################


FILENAME="root_archive.zip"



YEAR="2016"
PASS="pass1"
DATA="sim"                                                # data or sim
GROUP="GA_pp_MC_AOD"
TRAINNUM="1276"
TRAINNUMLONG="1276_20190128-1835"
LOCALPATH1="/media/gustav/external_drive/Data"            #general location of data
LOCALPATH2_1="pp_13TeV_2016_AOD_sim_1276"                 # specified folder for this data set (specific train run)

if [ $2 = "MC" ]; then
YEAR="2017"
DATA="sim"
GROUP="GA_pp_MC_AOD"
TRAINNUM="1276"
TRAINNUMLONG="1276_20190128-1835"
fi


#################################################################### DATA
if [ $2 = "data" ]; then

if [ $1 = "LHC16d" ]; then
PERIOD=""
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_1"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16d_EMC.txt"

elif [ $1 = "LHC16e" ]; then
PERIOD="$1"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_10"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16e_EMC.txt"

elif [ $1 = "LHC16g" ]; then
PERIOD="$1"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_2"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16g_EMC.txt"

elif [ $1 = "LHC16h" ]; then
PERIOD="$1"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_3"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16h_EMC.txt"

elif [ $1 = "LHC16i" ]; then
PERIOD="$1"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_4"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16i_EMC.txt"

elif [ $1 = "LHC16j" ]; then
PERIOD="$1"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_5"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16j_EMC.txt"

elif [ $1 = "LHC16k" ]; then
PERIOD="$1"
PASS="pass2"
TRAIN="$TRAINNUMLONG""_child_6"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16k_EMC.txt"

elif [ $1 = "LHC16l" ]; then
PERIOD="$1"
PASS="pass2"
TRAIN="$TRAINNUMLONG""_child_7"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16l_EMC.txt"

elif [ $1 = "LHC16o" ]; then
PERIOD="$1"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_8"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16oEMC.txt"

elif [ $1 = "LHC16p" ]; then
PERIOD="$1"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_9"
LOCALPATH2="$LOCALPATH2_1/$1"
RUNLIST="./runlists/runNumbersLHC16p_EMC.txt"


fi
fi




################################################# MC

if [ $2 = "MC" ]; then

if [ $1 = "LHC16d" ]; then
PERIOD="LHC17f6"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_1"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16d_EMC.txt"

elif [ $1 = "LHC16e" ]; then
PERIOD="LHC17f9"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_2"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16e_EMC.txt"

elif [ $1 = "LHC16g" ]; then
PERIOD="LHC17d17"
PASS="pass1"
TRAIN="$TRAINNUMLONG_""_child_3"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16g_EMC.txt"

elif [ $1 = "LHC16h" ]; then
PERIOD="LHC17f5"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_4"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16h_EMC.txt"

elif [ $1 = "LHC16i" ]; then
PERIOD="LHC17d3"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_5"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16i_EMC.txt"

elif [ $1 = "LHC16j" ]; then
PERIOD="LHC17e5"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_6"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16j_EMC.txt"

elif [ $1 = "LHC16k" ]; then
YEAR="2018"
PERIOD="LHC18f1"
PASS="pass2"
TRAIN="$TRAINNUMLONG""_child_9"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16k_EMC.txt"

elif [ $1 = "LHC16l" ]; then
YEAR="2018"
PERIOD="LHC18d8"
PASS="pass2"
TRAIN="$TRAINNUMLONG""_child_10"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16l_EMC.txt"

elif [ $1 = "LHC16o" ]; then
PERIOD="LHC17d16"
PASS="pass1"
TRAIN="$TRAINNUMLONG""_child_7"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16o_EMC.txt"

elif [ $1 = "LHC16p" ]; then
PERIOD="LHC17d18"
PASS="pass1"
child="_child_12"
TRAIN="$TRAINNUMLONG""_child_8"
LOCALPATH2="$LOCALPATH2_1/"$1"_MC"
RUNLIST="./runlists/runNumbersLHC16p_EMC.txt"

fi
fi


RUNNUM=()

#READ RUNNUMBERS
while read -r line || [[ -n "$line" ]]; do
    RUNNUM+="$line "
    echo "Text read from file: $line"
done < "$RUNLIST"



root -q ConnectToAlien.C


checkfiles=0
loopcounter=0

while [ $checkfiles -eq 0 ] && [ $loopcounter -lt 10 ]
do

echo "trying download"
echo "loop: $loopcounter"


checkfiles=0
loopcounter=$((loopcounter + 1))

for i in `echo $RUNNUM`
do
echo "Runnum: $i"


files=($LOCALPATH1/$LOCALPATH2/$i/*)
if [ ${#files[@]} -gt 1 ]; then continue;            ################## MUST BE 1 NORMALLY
if [ -e ""$LOCALPATH1"/"$LOCALPATH2"/"$i"/$FILENAME" ]; then echo ""$LOCALPATH1"/"$LOCALPATH2"/"$i"/$FILENAME"; continue;
else

echo "Run $i is not downloaded"

checkfiles=1



if [ $2 = "MC" ]; then
GRIDPATH="/alice/$DATA/$YEAR/$PERIOD/$i/AOD209/PWGGA/$GROUP/$TRAIN/"$FILENAME""
elif [ $2 = "data" ]; then
GRIDPATH="/alice/$DATA/$YEAR/$PERIOD/000$i/$PASS/PWGGA/$GROUP/$TRAIN/"$FILENAME""
fi

echo $GRIDPATH > GridPathFile.txt

root -q -b GridJobFileListDL.C\(\"GridPathFile.txt\",\"$LOCALPATH2\",\"$LOCALPATH1\",\""$FILENAME"\"\)


if [[ -d "$LOCALPATH1/$LOCALPATH2/$i" ]]; then
    echo "directory exists"
else
    mkdir "$LOCALPATH1/$LOCALPATH2/$i"
fi
cp -R "$LOCALPATH1/$LOCALPATH2/0/." "$LOCALPATH1/$LOCALPATH2/$i/"
rm -r "$LOCALPATH1/$LOCALPATH2/0"

fi

done
done
