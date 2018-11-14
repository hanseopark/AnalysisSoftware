dataOrMC=data #data or MC
ESDOrAOD=ESD #ESD or AOD
#energy=pPb_5TeV
energy=pp_5TeV

# set these if you want to download data
#yearData=2016
#periodData=LHC16q
#pass=pass1_CENT_woSDD
yearData=2017
periodData=LHC17p
pass=pass1_CENT_woSDD

# set these if you want to download MC
yearMC=2018
periodMC=LHC18f3_cent_woSDD_2
periodMCshort=LHC18f3_2
AODFILTER=202

# must be set either way
#RUN=265500
RUN=282078
downloadFile=root_archive.zip

RUNDL=$RUN
passDL=$pass

# give number of maxim files to be downloaded for the data set
NMaxFiles=20
counter=0

if [ $dataOrMC = "MC" ]; then
    if [ $ESDOrAOD = "AOD" ]; then
        RUNDL="$RUN/AOD$AODFILTER"
        downloadFile=root_archive.zip
    fi
    mkdir -p $energy/$periodMCshort/$ESDOrAOD/$RUN
    for PARTS in $( alien_ls /alice/sim/$yearMC/$periodMC/$RUNDL/ )
    do
        if [ -f $energy/$periodMCshort/$ESDOrAOD/$RUN/$PARTS/$downloadFile ]; then
            echo "file $downloadFile has already been copied for run " $RUN "and part" $PARTS
        else
            mkdir $energy/$periodMCshort/$ESDOrAOD/$RUN/$PARTS
            echo /alice/sim/$yearMC/$periodMC/$RUNDL/$PARTS/$downloadFile
            alien_cp alien:/alice/sim/$yearMC/$periodMC/$RUNDL/$PARTS/$downloadFile file:$energy/$periodMCshort/$ESDOrAOD/$RUN/$PARTS/$downloadFile
            unzip -t $energy/$periodMCshort/$ESDOrAOD/$RUN/$PARTS/$downloadFile
            md5sum $energy/$periodMCshort/$ESDOrAOD/$RUN/$PARTS/$downloadFile
        fi
        counter=`expr $counter + 1`
        if [ $counter = $NMaxFiles ]; then
          exit
        fi
    done
fi


if [ $dataOrMC = "data" ]; then
    if [ $ESDOrAOD = "AOD" ]; then
        passDL="$pass/AOD"
        downloadFile=aod_archive.zip
    fi
    mkdir -p $energy/$periodData/$pass/$ESDOrAOD/$RUN/
    for PARTS in $( alien_ls /alice/data/$yearData/$periodData/000$RUN/$passDL )
    do
        if [ -f $energy/$periodData/$pass/$ESDOrAOD/$RUN/$PARTS/$downloadFile ]; then
            echo "file $downloadFile has already been copied for run " $RUN "and part" $PARTS
        else
            mkdir $energy/$periodData/$pass/$ESDOrAOD/$RUN/$PARTS
            echo /alice/data/$yearData/$periodData/000$RUN/$passDL/$PARTS/$downloadFile
            alien_cp alien:/alice/data/$yearData/$periodData/000$RUN/$passDL/$PARTS/$downloadFile file:$energy/$periodData/$pass/$ESDOrAOD/$RUN/$PARTS/$downloadFile
            unzip -t $energy/$periodData/$pass/$ESDOrAOD/$RUN/$PARTS/$downloadFile
            md5sum $energy/$periodData/$pass/$ESDOrAOD/$RUN/$PARTS/$downloadFile
        fi
        counter=`expr $counter + 1`
        if [ $counter = $NMaxFiles ]; then
          exit
        fi
    done
fi
