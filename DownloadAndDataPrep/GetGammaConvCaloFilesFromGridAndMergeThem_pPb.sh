#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash

# This script has to be run with "bash"
source basicFunction.sh

# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1
MERGEONMC=1
SINGLERUN=1
SEPARATEON=0
MERGEONSINGLEData=0
MERGEONSINGLEMC=1
SPECIALMERGE=0
ISADDDOWNLOAD=0
ADDDOWNLOADALREADY=0
CLEANUPADDMERGE=1
CLEANUP=1
CLEANUPMAYOR=$2

function ParseJDLFilesDownloadAndMerge()
{
    echo $1
    maxJobnumbers=`cat $1 | wc -l`;
    if [ $maxJobnumbers > 0 ]; then
        jobnumbers=`cat $1`;
        echo $jobnumbers
        counterJob=1
        # Loop for all jobumbers
        rm fileNamesAddDownload.txt
        for jobnumber in $jobnumbers; do
            # which job are you currently analysing
            echo "$jobnumber, number $counterJob/ $maxJobnumbers";
            # get the jdl from the job number
            alien_ps -jdl $jobnumber > jdl.txt

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
            sed 's/,//g' jdl.txt > jdlmod.txt
            mv jdlmod.txt jdl.txt
            # find processed files
            cat jdl.txt | grep '/root_archive.zip' >> fileNamesAddDownload.txt
            counterJob=$((counterJob+1));
        done
        cat fileNamesAddDownload.txt
        mkdir -p $2/inputFailedStage1Merge
        echo "$2/inputFailedStage1Merge"
        maxNFiles=`cat fileNamesAddDownload.txt | wc -l`;
        sed 's/\/root_archive.zip//g' fileNamesAddDownload.txt > fileNamesAddDownloadMod.txt
        cat fileNamesAddDownloadMod.txt
        fileNames=`cat fileNamesAddDownloadMod.txt`;

        if [ $ADDDOWNLOADALREADY == 0 ]; then
            counterFiles=1
            ISADDDOWNLOAD=1
            # Loop for all files
            for fileName in $fileNames; do
                # which job are you currently analysing
                echo "number $counterFiles/ $maxNFiles";
                echo $fileName;

                CopyFileIfNonExisitent $2/inputFailedStage1Merge/$counterFiles $fileName 10
                counterFiles=$((counterFiles+1));
            done
            ISADDDOWNLOAD=0

            text=`ls $2/inputFailedStage1Merge/1/*.zip`
            ls $2/inputFailedStage1Merge/1/*.root > fileListAdd.txt
            NCurrSlash=`tr -dc '/' <<<"$text" | awk '{ print length; }'`
            ADDFILESTOMERGE=`cat fileListAdd.txt`;

#             echo $ADDFILESTOMERGE
            for TOMERGE in $ADDFILESTOMERGE; do
                echo $TOMERGE
                TOMERGE2=$TOMERGE;
                TOMERGE=`echo $TOMERGE  | cut -d "/" -f $((NCurrSlash+1)) `
                echo $TOMERGE
                TOMERGE2=`echo $TOMERGE2  | cut -d "." -f 1 `
                echo $TOMERGE2
                hadd -n 10 -f $2/inputFailedStage1Merge/$TOMERGE $2/inputFailedStage1Merge/*/$TOMERGE
                if [ $SEPARATEON == 1 ]; then
                    if [[ $TOMERGE == *"GammaConvCalo_"* ]]; then
                        SeparateCutsIfNeeded $TOMERGE2 2
                    elif [[ $TOMERGE == *"GammaConvV1_"* ]]; then
                        SeparateCutsIfNeededConv $TOMERGE2 0
                    elif [[ $TOMERGE == *"GammaCalo_"* ]]; then
                        SeparateCutsIfNeededCalo $TOMERGE2 4
                    fi
                fi
            done
        fi

        if [ $CLEANUPADDMERGE == 1 ]; then
            counterFiles=1
            # Loop for all files
            for fileName in $fileNames; do
                # which job are you currently analysing
                echo "cleaning $counterFiles/ $maxNFiles";
                echo $2/inputFailedStage1Merge/$counterFiles;
                rm -rf $2/inputFailedStage1Merge/$counterFiles;
                counterFiles=$((counterFiles+1));
            done
        fi
   fi
}

# default trainconfigurations
LHC13bcData="";
LHC13bData="";
LHC13cData=""; #ESD
LHC13deData=""; #ESD
LHC13dData=""; #ESD
LHC13eData=""; #ESD
LHC13fData=""; #ESD
LHC13e7MC="";
LHC13b2_efix="";
LHC13b2_efix_p1MC="";
LHC13b2_efix_p2MC="";
LHC13b2_efix_p3MC="" ;
LHC13b2_efix_p4MC="";
LHC13b4fixMC="";
LHC13b4plusMC="";
LHC16c3aMC="";
LHC16c3a2MC="";
LHC16c3bMC="";
LHC16c3b2MC="";

passNr="2";

NSlashes=10;

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pPb
elif [ $1 = "fbockExt" ]; then
    BASEDIR=/media/fbock/Elements/OutputLegoTrains/pPb
elif [ $1 = "fbockGSI" ]; then
    BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pPb
elif [ $1 = "leardini" ]; then
    BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then
    BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
elif [ $1 = "leardiniGSI" ]; then
    BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/pPb
elif [ $1 = "passfeld" ]; then
    BASEDIR=~/work/Gridoutput/pPb
elif [ $1 = "passfeldMAF" ]; then
    BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then
    BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pPb
elif [ $1 = "amarin" ]; then
    BASEDIR=/Users/marin/
elif [ $1 = "amarinGSI" ]; then
    BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/pPb
elif [ $1 = "amarinALICESERV1" ]; then
    BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/
elif [ $1 = "mwilde" ]; then
    BASEDIR=~/work/GridOutput
elif [ $1 = "mwildeGSI" ]; then
    BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pPb
elif [ $1 = "pgonzales" ]; then
    BASEDIR=~/work/GridOutput
elif [ $1 = "pgonzalesGSI" ]; then
    BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pPb
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
fi
# Definitition of number of slashes in your path to different depths
NSlashesBASE=`tr -dc '/' <<<"$BASEDIR" | wc -c`
NSlashes=`expr $NSlashesBASE + 4`
NSlashes2=`expr $NSlashes - 1`
NSlashes3=`expr $NSlashes + 1`
NSlashes4=`expr $NSlashes + 2`
echo "$NSlashesBASE $NSlashes $NSlashes2 $NSlashes3 $NSlashes4"


# TRAINDIR=Legotrain-vAN20170525FF-newDefaultPlusSys
# LHC13bData="637"; #pass 3
# LHC13cData="638"; #pass 2
# LHC13b2_efix_p1MC="978";
# LHC13b2_efix_p2MC="982";
# LHC13b2_efix_p3MC="980";
# LHC13b2_efix_p4MC="981";
# LHC13e7MC="977";
# LHC13bData="640"; #pass 3
# LHC13cData="641"; #pass 2
# LHC13b2_efix_p1MC="973";
# LHC13b2_efix_p2MC="974";
# LHC13b2_efix_p3MC="975";
# LHC13b2_efix_p4MC="976";

# LHC13b2_efix_p1MC="983";
# LHC13b2_efix_p2MC="984";
# LHC13b2_efix_p3MC="985";
# LHC13b2_efix_p4MC="986";

# LHC13b2_efix_p1MC="987";
# LHC13b2_efix_p2MC="988";
# LHC13b2_efix_p3MC="989";
# LHC13b2_efix_p4MC="990";

# LHC13b2_efix_p1MC="999";
# LHC13b2_efix_p2MC="1000";
# LHC13b2_efix_p3MC="993";
# LHC13b2_efix_p4MC="994";
# LHC13b2_efix_p1MC="995";
# LHC13b2_efix_p2MC="996";
# LHC13b2_efix_p3MC="997";
# LHC13b2_efix_p4MC="998";

# LHC13b2_efix_p1MC="1001";
# LHC13b2_efix_p2MC="1002";
# LHC13b2_efix_p3MC="1003";
# LHC13b2_efix_p4MC="1004";
# LHC13bData="647"; #pass 3
# LHC13cData="648"; #pass 2
# LHC13b2_efix_p1MC="1026";
# LHC13b2_efix_p2MC="1046";
# LHC13b2_efix_p3MC="1047";
# LHC13b2_efix_p4MC="1048";
# LHC13bData="660"; #pass 3
# LHC13cData="661"; #pass 2

# TRAINDIR=Legotrain-vAN20180607-trigg2
# LHC13bcData="765"; #pass 3
# LHC13bData="child_1"; #pass 3
# LHC13cData="child_2"; #pass 2
# LHC13deData="766"; #pass 3
# LHC13dData="child_1"; #pass 3
# LHC13eData="child_2"; #pass 2
# LHC13fData="761"; #pass 2

# TRAINDIR=Legotrain-vAN20180618-dirGamma
# LHC13bcData="759"; #pass 3
# LHC13bData="child_1"; #pass 3
# LHC13cData="child_2"; #pass 2
# # LHC13deData="741"; #pass 3
# # LHC13dData="child_1"; #pass 3
# # LHC13eData="child_2"; #pass 2
#
# LHC13b2_efix="1279"
# LHC13b2_efix_p1MC="child_1";
# LHC13b2_efix_p2MC="child_2";
# LHC13b2_efix_p3MC="child_3";
# LHC13b2_efix_p4MC="child_4";

TRAINDIR=Legotrain-vAN20180718-triggQA
# LHC13bcData="765"; #pass 3
# LHC13bData="child_1"; #pass 3
# LHC13cData="child_2"; #pass 2
# LHC13deData="763"; #pass 3
# LHC13dData="child_1"; #pass 4
# LHC13eData="child_2"; #pass 4
# LHC13fData="764"; #pass 4
# LHC13b4fixMC="1286";
# LHC13b4plusMC="1287";
# LHC16c3aMC="1290";
LHC16c3a2MC="1292";
# LHC16c3bMC="1291";
# LHC16c3b2MC="1293";

OUTPUTDIR=$BASEDIR/$TRAINDIR
OUTPUTDIRMC=$BASEDIR/$TRAINDIR/GA_pPb_MC
OUTPUTDIRData=$BASEDIR/$TRAINDIR/GA_pPb
ALIENDIRData="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/"
ALIENDIRMC="/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/"
mkdir -p $OUTPUTDIR/CutSelections
mkdir -p $OUTPUTDIR/SingleFiles
mkdir -p $OUTPUTDIR/JJMCSingleBins

# parse grid directories for correct train output dir for LHC13bc
FindCorrectTrainDirectory $LHC13bData $OUTPUTDIRData $ALIENDIRData $LHC13bcData
HAVELHC13b=$tempBool
LHC13bData=$tempDir
OUTPUTDIR_LHC13b=$tempPath
echo "$HAVELHC13b $LHC13bData $OUTPUTDIR_LHC13b"

FindCorrectTrainDirectory $LHC13cData $OUTPUTDIRData $ALIENDIRData $LHC13bcData
HAVELHC13c=$tempBool
LHC13cData=$tempDir
OUTPUTDIR_LHC13c=$tempPath
echo "$HAVELHC13c $LHC13cData $OUTPUTDIR_LHC13c"

# parse grid directories for correct train output dir for LHC13de
FindCorrectTrainDirectory $LHC13dData $OUTPUTDIRData $ALIENDIRData $LHC13deData
HAVELHC13d=$tempBool
LHC13dData=$tempDir
OUTPUTDIR_LHC13d=$tempPath
echo "$HAVELHC13d $LHC13dData $OUTPUTDIR_LHC13d"

FindCorrectTrainDirectory $LHC13eData $OUTPUTDIRData $ALIENDIRData $LHC13deData
HAVELHC13e=$tempBool
LHC13eData=$tempDir
OUTPUTDIR_LHC13e=$tempPath
echo "$HAVELHC13e $LHC13eData $OUTPUTDIR_LHC13e"

# parse grid directories for correct train output dir for LHC13f
FindCorrectTrainDirectory $LHC13fData $OUTPUTDIRData $ALIENDIRData
HAVELHC13f=$tempBool
LHC13fData=$tempDir
OUTPUTDIR_LHC13f=$tempPath
echo "$HAVELHC13f $LHC13fData $OUTPUTDIR_LHC13f"

# Settting directories for LHC13b2_efix productions
FindCorrectTrainDirectory $LHC13b2_efix_p1MC $OUTPUTDIRMC $ALIENDIRMC $LHC13b2_efix
HAVELHC13b2efixp1=$tempBool
LHC13b2_efix_p1MC=$tempDir
OUTPUTDIR_LHC13b2_efix_p1=$tempPath
echo "$HAVELHC13b2efixp1 $LHC13b2_efix_p1MC $OUTPUTDIR_LHC13b2_efix_p1"

FindCorrectTrainDirectory $LHC13b2_efix_p2MC $OUTPUTDIRMC $ALIENDIRMC $LHC13b2_efix
HAVELHC13b2efixp2=$tempBool
LHC13b2_efix_p2MC=$tempDir
OUTPUTDIR_LHC13b2_efix_p2=$tempPath
echo "$HAVELHC13b2efixp2 $LHC13b2_efix_p2MC $OUTPUTDIR_LHC13b2_efix_p2"

FindCorrectTrainDirectory $LHC13b2_efix_p3MC $OUTPUTDIRMC $ALIENDIRMC $LHC13b2_efix
HAVELHC13b2efixp3=$tempBool
LHC13b2_efix_p3MC=$tempDir
OUTPUTDIR_LHC13b2_efix_p3=$tempPath
echo "$HAVELHC13b2efixp3 $LHC13b2_efix_p3MC $OUTPUTDIR_LHC13b2_efix_p3"

FindCorrectTrainDirectory $LHC13b2_efix_p4MC $OUTPUTDIRMC $ALIENDIRMC $LHC13b2_efix
HAVELHC13b2efixp4=$tempBool
LHC13b2_efix_p4MC=$tempDir
OUTPUTDIR_LHC13b2_efix_p4=$tempPath
echo "$HAVELHC13b2efixp4 $LHC13b2_efix_p4MC $OUTPUTDIR_LHC13b2_efix_p4"

# Settting directories for LHC13e7 production
FindCorrectTrainDirectory $LHC13e7MC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC13e7=$tempBool
LHC13e7MC=$tempDir
OUTPUTDIR_LHC13e7=$tempPath
echo "$HAVELHC13e7 $LHC13e7MC $OUTPUTDIR_LHC13e7"

# Settting directories for JJ production
FindCorrectTrainDirectory $LHC13b4fixMC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC13b4fix=$tempBool
LHC13b4fixMC=$tempDir
OUTPUTDIR_LHC13b4fix=$tempPath
echo "$HAVELHC13b4fix $LHC13b4fixMC $OUTPUTDIR_LHC13b4fix"

FindCorrectTrainDirectory $LHC13b4plusMC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC13b4plus=$tempBool
LHC13b4plusMC=$tempDir
OUTPUTDIR_LHC13b4plus=$tempPath
echo "$HAVELHC13b4plus $LHC13b4plusMC $OUTPUTDIR_LHC13b4plus"

FindCorrectTrainDirectory $LHC16c3aMC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC16c3a=$tempBool
LHC16c3aMC=$tempDir
OUTPUTDIR_LHC16c3a=$tempPath
echo "$HAVELHC16c3a $LHC16c3aMC $OUTPUTDIR_LHC16c3a"

FindCorrectTrainDirectory $LHC16c3a2MC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC16c3a2=$tempBool
LHC16c3a2MC=$tempDir
OUTPUTDIR_LHC16c3a2=$tempPath
echo "$HAVELHC16c3a2 $LHC16c3a2MC $OUTPUTDIR_LHC16c3a2"

FindCorrectTrainDirectory $LHC16c3bMC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC16c3b=$tempBool
LHC16c3bMC=$tempDir
OUTPUTDIR_LHC16c3b=$tempPath
echo "$HAVELHC16c3b $LHC16c3bMC $OUTPUTDIR_LHC16c3b"

FindCorrectTrainDirectory $LHC16c3b2MC $OUTPUTDIRMC $ALIENDIRMC
HAVELHC16c3b2=$tempBool
LHC16c3b2MC=$tempDir
OUTPUTDIR_LHC16c3b2=$tempPath
echo "$HAVELHC16c3b2 $LHC16c3b2MC $OUTPUTDIR_LHC16c3b2"

if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC13b == 1 ]; then
        echo "downloading LHC13b"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b/$runNumber "/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/GA_pPb/$LHC13bData" $NSlashes3
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC13b/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b.txt`
                ls $OUTPUTDIR_LHC13b/$firstrunNumber/GammaConvCalo_*.root > fileLHC13b.txt
                fileNumbers=`cat fileLHC13b.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13b/GammaConvCalo_$number.root $OUTPUTDIR_LHC13b/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b "$ALIENDIRData$LHC13bData/merge_runlist_1" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b "$ALIENDIRData$LHC13bData/merge_runlist_2" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi

    if [ $HAVELHC13c == 1 ]; then
        echo "downloading LHC13c"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13c.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13c/$runNumber "/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/GA_pPb/$LHC13cData" $NSlashes3
            done;
            if [ $MERGEONSINGLEData == 1 ]  && [ ! -f $OUTPUTDIR_LHC13c/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13c.txt`
                ls $OUTPUTDIR_LHC13c/$firstrunNumber/GammaConvCalo_*.root > fileLHC13c.txt
                fileNumbers=`cat fileLHC13c.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13c/GammaConvCalo_$number.root $OUTPUTDIR_LHC13c/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13c/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13c "$ALIENDIRData$LHC13cData/merge_runlist_1" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13c "$ALIENDIRData$LHC13cData/merge_runlist_2" PHOSGood $NSlashes3 "" kFALSE
        fi

    fi
    if [ $HAVELHC13d == 1 ]; then
        echo "downloading LHC13d"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13d_pass4.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13d/$runNumber "/alice/data/2013/LHC13d/000$runNumber/pass4/PWGGA/GA_pPb/$LHC13dData" $NSlashes3
            done;
            if [ $MERGEONSINGLEData == 1 ]  && [ ! -f $OUTPUTDIR_LHC13d/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13d_pass4.txt`
                ls $OUTPUTDIR_LHC13d/$firstrunNumber/GammaConvCalo_*.root > fileLHC13d.txt
                fileNumbers=`cat fileLHC13d.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13d/GammaConvCalo_$number.root $OUTPUTDIR_LHC13d/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13d/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13d "$ALIENDIRData$LHC13dData/merge" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13d "$ALIENDIRData$LHC13dData/merge" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi

    if [ $HAVELHC13e == 1 ]; then
        echo "downloading LHC13e"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13e_pass4_DPGTracks.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13e/$runNumber "/alice/data/2013/LHC13e/000$runNumber/pass4/PWGGA/GA_pPb/$LHC13eData" $NSlashes3
            done;
            if [ $MERGEONSINGLEData == 1 ]  && [ ! -f $OUTPUTDIR_LHC13e/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13e_pass4_DPGTracks.txt`
                ls $OUTPUTDIR_LHC13e/$firstrunNumber/GammaConvCalo_*.root > fileLHC13e.txt
                fileNumbers=`cat fileLHC13e.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13e/GammaConvCalo_$number.root $OUTPUTDIR_LHC13e/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13e/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13e "$ALIENDIRData$LHC13eData/merge_runlist_2" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13e "$ALIENDIRData$LHC13eData/merge_runlist_2" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi
    if [ $HAVELHC13f == 1 ]; then
        echo "downloading LHC13f"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13f_pass4.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13f/$runNumber "/alice/data/2013/LHC13f/000$runNumber/pass4/PWGGA/GA_pPb/$LHC13fData" $NSlashes3
            done;
            if [ $MERGEONSINGLEData == 1 ]  && [ ! -f $OUTPUTDIR_LHC13f/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13f_pass4.txt`
                ls $OUTPUTDIR_LHC13f/$firstrunNumber/GammaConvCalo_*.root > fileLHC13f.txt
                fileNumbers=`cat fileLHC13f.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13f/GammaConvCalo_$number.root $OUTPUTDIR_LHC13f/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13f/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13f "$ALIENDIRData$LHC13fData/merge" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13f "$ALIENDIRData$LHC13fData/merge" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi

    if [ $HAVELHC13b2efixp1 == 1 ]; then
        echo "downloading LHC13b2_efix_p1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1/$runNumber "/alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p1/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix1.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p1/$firstrunNumber/GammaConvCalo_*.root > fileLHC13b2_efix_p1.txt
                fileNumbers=`cat fileLHC13b2_efix_p1.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p1/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p1 "$ALIENDIRMC$LHC13b2_efix_p1MC/merge_runlist_1" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p1 "$ALIENDIRMC$LHC13b2_efix_p1MC/merge_runlist_2" PHOSGood $NSlashes3 "" kFALSE
        fi
        ParseJDLFilesDownloadAndMerge missingMergesLHC13b2_efix_p1.txt  $OUTPUTDIR_LHC13b2_efix_p1

    fi
    if [ $HAVELHC13b2efixp2 == 1 ]; then
        echo "downloading LHC13b2_efix_p2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix2.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2/$runNumber "/alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p2/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix2.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p2/$firstrunNumber/GammaConvCalo_*.root > fileLHC13b2_efix_p2.txt
                fileNumbers=`cat fileLHC13b2_efix_p2.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p2/mergedAllConvCalo.txt
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p2/GammaConvCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvCalo_$number.root
                done;
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p2 "$ALIENDIRMC$LHC13b2_efix_p2MC/merge_runlist_1" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p2 "$ALIENDIRMC$LHC13b2_efix_p2MC/merge_runlist_2" PHOSGood $NSlashes3 "" kFALSE
        fi
        ParseJDLFilesDownloadAndMerge missingMergesLHC13b2_efix_p2.txt  $OUTPUTDIR_LHC13b2_efix_p2
    fi

    if [ $HAVELHC13b2efixp3 == 1 ]; then
        echo "downloading LHC13b2_efix_p3"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix3.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3/$runNumber "/alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p3/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix3.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p3/$firstrunNumber/GammaConvCalo_*.root > fileLHC13b2_efix_p3.txt
                fileNumbers=`cat fileLHC13b2_efix_p3.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p3/GammaConvCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p3/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p3 "$ALIENDIRMC$LHC13b2_efix_p3MC/merge_runlist_1" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p3 "$ALIENDIRMC$LHC13b2_efix_p3MC/merge_runlist_2" PHOSGood $NSlashes3 "" kFALSE
        fi
        ParseJDLFilesDownloadAndMerge missingMergesLHC13b2_efix_p3.txt  $OUTPUTDIR_LHC13b2_efix_p3
    fi
    if [ $HAVELHC13b2efixp4 == 1 ]; then
        echo "downloading LHC13b2_efix_p4"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix4.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4/$runNumber "/alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p4/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix4.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p4/$firstrunNumber/GammaConvCalo_*.root > fileLHC13b2_efix_p4.txt
                fileNumbers=`cat fileLHC13b2_efix_p4.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p4/GammaConvCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p4/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p4 "$ALIENDIRMC$LHC13b2_efix_p4MC/merge_runlist_1" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b2_efix_p4 "$ALIENDIRMC$LHC13b2_efix_p4MC/merge_runlist_2" PHOSGood $NSlashes3 "" kFALSE
        fi
        ParseJDLFilesDownloadAndMerge missingMergesLHC13b2_efix_p4.txt  $OUTPUTDIR_LHC13b2_efix_p4
    fi
    if [ $HAVELHC13e7 == 1 ]; then
        echo "downloading LHC13e7"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13e7.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7/$runNumber "/alice/sim/2013/LHC13e7/$runNumber/PWGGA/GA_pPb_MC/$LHC13e7MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13e7/mergedAllConvCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13e7.txt`
                ls $OUTPUTDIR_LHC13e7/$firstrunNumber/GammaConvCalo_*.root > fileLHC13e7.txt
                fileNumbers=`cat fileLHC13e7.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                        number=$number\_$alpha
                        echo $number
                    fi
                    echo $number
                    hadd -f $OUTPUTDIR_LHC13e7/GammaConvCalo_$number.root $OUTPUTDIR_LHC13e7/*/GammaConvCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13e7/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13e7 "$ALIENDIRMC$LHC13e7MC/merge" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13e7 "$ALIENDIRMC$LHC13e7MC/merge" PHOSGood $NSlashes3 "" kFALSE
        fi
        ParseJDLFilesDownloadAndMerge missingMergesLHC13e7.txt  $OUTPUTDIR_LHC13e7
    fi

    currentDir=$PWD
    if [ $HAVELHC13b4fix == 1 ]; then
        echo "downloading LHC13b4fix"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b4_fix.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binsJetJetLHC13b4.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
#                     CopyFileIfNonExisitent $OUTPUTDIR_LHC13b4fix/$binNumber/$runNumber "/alice/sim/2013/LHC13b4_fix/$runNumber/$binNumber/PWGGA/GA_pPb_MC/$LHC13b4fixMC" $NSlashes3 "/alice/sim/2013/LHC13b4_fix/$runNumber/$binNumber/PWGGA/GA_pPb_MC/$LHC13b4fixMC/" kTRUE
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b4fix/$binNumber/$runNumber "/alice/sim/2013/LHC13b4_fix/$runNumber/$binNumber/PWGGA/GA_pPb_MC/$LHC13b4fixMC" $NSlashes3 "none" kTRUE
                done;
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b4fix/mergedAllConvCalo.txt ]; then
                echo "HERERRE"
                cd $currentDir
                rm $OUTPUTDIR_LHC13b4fix/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b4_fix.txt`
                firstbinNumber=`head -n1 runlists/binsJetJetLHC13b4.txt`
                ls $OUTPUTDIR_LHC13b4fix/$firstbinNumber/$firstrunNumber/GammaConvCalo_*.root > fileLHC13b4.txt
                fileNumbers=`cat fileLHC13b4.txt`
                MergeAccordingToSpecificRunlist fileLHC13b4.txt $OUTPUTDIR_LHC13b4fix $NSlashes4 GammaConvCalo EMCandPCMGood runlists/runNumbersLHC13b4_fix.txt runlists/binsJetJetLHC13b4.txt
                MergeAccordingToSpecificRunlist fileLHC13b4.txt $OUTPUTDIR_LHC13b4fix $NSlashes4 GammaConvCalo PHOSGood runlists/runNumbersLHC13b4_fix.txt runlists/binsJetJetLHC13b4.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b4fix "$ALIENDIRMC$LHC13b4_fix/merge" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b4fix "$ALIENDIRMC$LHC13b4_fix/merge" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi
    if [ $HAVELHC13b4plus == 1 ]; then
        echo "downloading LHC13b4plus"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b4_plus.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binsJetJetLHC13b4.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
#                     CopyFileIfNonExisitent $OUTPUTDIR_LHC13b4plus/$binNumber/$runNumber "/alice/sim/2013/LHC13b4_plus/$runNumber/$binNumber/PWGGA/GA_pPb_MC/$LHC13b4plusMC" $NSlashes3 "/alice/sim/2013/LHC13b4_plus/$runNumber/$binNumber/PWGGA/GA_pPb_MC/$LHC13b4plusMC/" kTRUE
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b4plus/$binNumber/$runNumber "/alice/sim/2013/LHC13b4_plus/$runNumber/$binNumber/PWGGA/GA_pPb_MC/$LHC13b4plusMC" $NSlashes3 "none" kTRUE
                done;
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b4plus/mergedAllConvCalo.txt ]; then
                echo "HERERRE"
                cd $currentDir
                rm $OUTPUTDIR_LHC13b4plus/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b4_plus.txt`
                firstbinNumber=`head -n1 runlists/binsJetJetLHC13b4.txt`
                ls $OUTPUTDIR_LHC13b4plus/$firstbinNumber/$firstrunNumber/GammaConvCalo_*.root > fileLHC13b4.txt
                fileNumbers=`cat fileLHC13b4.txt`
                MergeAccordingToSpecificRunlist fileLHC13b4.txt $OUTPUTDIR_LHC13b4plus $NSlashes4 GammaConvCalo EMCandPCMGood runlists/runNumbersLHC13b4_plus.txt runlists/binsJetJetLHC13b4.txt
                MergeAccordingToSpecificRunlist fileLHC13b4.txt $OUTPUTDIR_LHC13b4plus $NSlashes4 GammaConvCalo PHOSGood runlists/runNumbersLHC13b4_plus.txt runlists/binsJetJetLHC13b4.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b4plus "$ALIENDIRMC$LHC13b4_plus/merge" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC13b4plus "$ALIENDIRMC$LHC13b4_plus/merge" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi

    if [ $HAVELHC16c3a == 1 ]; then
        echo "downloading LHC16c3a"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16c3a.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binsJetJetLHC16c3a.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
#                     CopyFileIfNonExisitent $OUTPUTDIR_LHC16c3a/$binNumber/$runNumber "/alice/sim/2016/LHC16c3a/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC16c3aMC" $NSlashes3 "/alice/sim/2016/LHC16c3a/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC16c3aMC/" kTRUE
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC16c3a/$binNumber/$runNumber "/alice/sim/2016/LHC16c3a/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC16c3aMC" $NSlashes3 "none" kTRUE
                done;
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC16c3a/mergedAllConvCalo.txt ]; then
                echo "HERERRE"
                cd $currentDir
                rm $OUTPUTDIR_LHC16c3a/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16c3a.txt`
                firstbinNumber=`head -n1 runlists/binsJetJetLHC16c3a.txt`
                ls $OUTPUTDIR_LHC16c3a/$firstbinNumber/$firstrunNumber/GammaConvCalo_*.root > fileLHC16c3a.txt
                fileNumbers=`cat fileLHC16c3a.txt`
                MergeAccordingToSpecificRunlist fileLHC16c3a.txt $OUTPUTDIR_LHC16c3a $NSlashes4 GammaConvCalo EMCandPCMGood runlists/runNumbersLHC16c3a.txt runlists/binsJetJetLHC16c3a.txt
                MergeAccordingToSpecificRunlist fileLHC16c3a.txt $OUTPUTDIR_LHC16c3a $NSlashes4 GammaConvCalo PHOSGood runlists/runNumbersLHC16c3a.txt runlists/binsJetJetLHC16c3a.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16c3a "$ALIENDIRMC$LHC16c3a/merge" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16c3a "$ALIENDIRMC$LHC16c3a/merge" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi


    if [ $HAVELHC16c3a2 == 1 ]; then
        echo "downloading LHC16c3a2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC16c3a2.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binsJetJetLHC16c3a.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
#                     CopyFileIfNonExisitent $OUTPUTDIR_LHC16c3a2/$binNumber/$runNumber "/alice/sim/2016/LHC16c3a2/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC16c3a2MC" $NSlashes3 "/alice/sim/2016/LHC16c3a2/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC16c3a2MC/" kTRUE
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC16c3a2/$binNumber/$runNumber "/alice/sim/2016/LHC16c3a2/$binNumber/$runNumber/PWGGA/GA_pPb_MC/$LHC16c3a2MC" $NSlashes3 "none" kTRUE
                done;
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC16c3a2/mergedAllConvCalo.txt ]; then
                echo "HERERRE"
                cd $currentDir
                rm $OUTPUTDIR_LHC16c3a2/GammaConvCalo*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC16c3a2.txt`
                firstbinNumber=`head -n1 runlists/binsJetJetLHC16c3a.txt`
                ls $OUTPUTDIR_LHC16c3a2/$firstbinNumber/$firstrunNumber/GammaConvCalo_*.root > fileLHC16c3a2.txt
                fileNumbers=`cat fileLHC16c3a2.txt`
                MergeAccordingToSpecificRunlist fileLHC16c3a2.txt $OUTPUTDIR_LHC16c3a2 $NSlashes4 GammaConvCalo EMCandPCMGood runlists/runNumbersLHC16c3a2.txt runlists/binsJetJetLHC16c3a.txt
                MergeAccordingToSpecificRunlist fileLHC16c3a2.txt $OUTPUTDIR_LHC16c3a2 $NSlashes4 GammaConvCalo PHOSGood runlists/runNumbersLHC16c3a2.txt runlists/binsJetJetLHC16c3a.txt
            fi
        else
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16c3a2 "$ALIENDIRMC$LHC16c3a2/merge" EMCandPCMGood $NSlashes3 "" kFALSE
            CopyFileIfNonExisitentDiffList $OUTPUTDIR_LHC16c3a2 "$ALIENDIRMC$LHC16c3a2/merge" PHOSGood $NSlashes3 "" kFALSE
        fi
    fi

    echo -e "EMCandPCMGood\nPHOSGood" > runlistsToMerge.txt
    listsToMerge=`cat runlistsToMerge.txt`
    for runListName in $listsToMerge; do
        if [ $HAVELHC13b == 1 ]; then
            ls $OUTPUTDIR_LHC13b/GammaConvCalo-$runListName\_*.root > fileLHC13b.txt
            fileNumbers=`cat fileLHC13b.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass$passNr-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13c == 1 ]; then
            ls $OUTPUTDIR_LHC13c/GammaConvCalo-$runListName\_*.root > fileLHC13c.txt
            fileNumbers=`cat fileLHC13c.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass$passNr-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13d == 1 ]; then
            ls $OUTPUTDIR_LHC13d/GammaConvCalo-$runListName\_*.root > fileLHC13d.txt
            fileNumbers=`cat fileLHC13d.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13d $NSlashes "LHC13d-pass$passNr-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13e == 1 ]; then
            ls $OUTPUTDIR_LHC13e/GammaConvCalo-$runListName\_*.root > fileLHC13e.txt
            fileNumbers=`cat fileLHC13e.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13e $NSlashes "LHC13e-pass$passNr-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13f == 1 ]; then
            ls $OUTPUTDIR_LHC13f/GammaConvCalo-$runListName\_*.root > fileLHC13f.txt
            fileNumbers=`cat fileLHC13f.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13f $NSlashes "LHC13f-pass$passNr-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13b2efixp1 == 1 ]; then
            ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvCalo-$runListName\_*.root > fileLHC13b2efixp1.txt
            fileNumbers=`cat fileLHC13b2efixp1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b2_efix_p1 $NSlashes "MC_LHC13b2_efix_p1-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13b2efixp2 == 1 ]; then
            echo "entered"
            ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvCalo-$runListName\_*.root
            ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvCalo-$runListName\_*.root > fileLHC13b2efixp2.txt
            fileNumbers=`cat fileLHC13b2efixp2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b2_efix_p2 $NSlashes "MC_LHC13b2_efix_p2-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13b2efixp3 == 1 ]; then
            ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvCalo-$runListName\_*.root > fileLHC13b2efixp3.txt
            fileNumbers=`cat fileLHC13b2efixp3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b2_efix_p3 $NSlashes "MC_LHC13b2_efix_p3-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13b2efixp4 == 1 ]; then
            ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvCalo-$runListName\_*.root > fileLHC13b2efixp4.txt
            fileNumbers=`cat fileLHC13b2efixp4.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b2_efix_p4 $NSlashes "MC_LHC13b2_efix_p4-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13e7 == 1 ]; then
            ls $OUTPUTDIR_LHC13e7/GammaConvCalo_*.root > fileLHC13e7.txt
            fileNumbers=`cat fileLHC13e7.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13e7 $NSlashes "MC_LHC13e7-$runListName" "-$runListName"
            done;
        fi

        if [ $HAVELHC13b4fix == 1 ]; then
            ls $OUTPUTDIR_LHC13b4fix/GammaConvCalo-$runListName\_*.root > fileLHC13b4fix.txt
            fileNumbers=`cat fileLHC13b4fix.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b4fix $NSlashes "MC_LHC13b4fix-$runListName" "-$runListName"
            done;
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                ls $OUTPUTDIR_LHC13b4fix/GammaConvCalo-$runListName\_*.root > fileLHC13b4fix.txt
                fileNumbers=`cat fileLHC13b4fix.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes)) 2
                    cp $OUTPUTDIR_LHC13b4fix/$binNumber/GammaConvCalo-$runListName\_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC13b4fix-$binNumber\_$runListName\_$number.root
                done
            done;
        fi
        if [ $HAVELHC13b4plus == 1 ]; then
            ls $OUTPUTDIR_LHC13b4plus/GammaConvCalo-$runListName\_*.root > fileLHC13b4plus.txt
            fileNumbers=`cat fileLHC13b4plus.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC13b4plus $NSlashes "MC_LHC13b4plus-$runListName" "-$runListName"
            done;
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                ls $OUTPUTDIR_LHC13b4plus/GammaConvCalo-$runListName\_*.root > fileLHC13b4plus.txt
                fileNumbers=`cat fileLHC13b4plus.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes)) 2
                    cp $OUTPUTDIR_LHC13b4plus/$binNumber/GammaConvCalo-$runListName\_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC13b4plus-$binNumber\_$runListName\_$number.root
                done
            done;
        fi

        if [ $HAVELHC16c3a == 1 ]; then
            ls $OUTPUTDIR_LHC16c3a/GammaConvCalo-$runListName\_*.root > fileLHC16c3a.txt
            fileNumbers=`cat fileLHC16c3a.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededPCMCalo $fileName $OUTPUTDIR_LHC16c3a $NSlashes "MC_LHC16c3a-$runListName" "-$runListName"
            done;
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                ls $OUTPUTDIR_LHC16c3a/GammaConvCalo-$runListName\_*.root > fileLHC16c3a.txt
                fileNumbers=`cat fileLHC16c3a.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes)) 2
                    cp $OUTPUTDIR_LHC16c3a/$binNumber/GammaConvCalo-$runListName\_$number.root  $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_LHC16c3a-$binNumber\_$runListName\_$number.root
                done
            done;
        fi

    done
    echo "--> rewrite of files done"

    if [ $MERGEON == 1 ]; then
        echo -e "EMCandPCMGood\nPHOSGood" > runlistsToMerge.txt
        listsToMerge=`cat runlistsToMerge.txt`
        for runListName in $listsToMerge; do
            ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr-$runListName\_*.root > filesForMerging.txt
            filesForMerging=`cat filesForMerging.txt`
            periodList=`echo -e "b\nc"`

            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                echo $number
                nameOut=""
                rm listCurrMerge.txt
                echo $fileName
                for periodID in $periodList; do
                    echo $periodID
                    currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-pass$passNr-$runListName\_$number.root
                    if [ -f $currFile ]; then
                        nameOut+=$periodID
                        echo -e "$currFile\n" >> listCurrMerge.txt
                    fi
                done
                MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass$passNr-$runListName\_$number.root
            done

            if [ $HAVELHC13d == 1 ] && [ $HAVELHC13e == 1 ]; then
                ls $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                periodList=`echo -e "d\ne"`
                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass$passNr-$runListName\_$number.root
                done
                ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                periodList=`echo -e "b\n\nc\nd\ne"`
                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass$passNr-$runListName\_$number.root
                done
            fi

            if [ $HAVELHC13d == 1 ] && [ $HAVELHC13e == 1 ] && [ $HAVELHC13f == 1 ]; then
                periodList=`echo -e "d\ne\nf"`

                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass$passNr-$runListName\_$number.root
                done

                periodList=`echo -e "b\n\nc\nd\ne\nf"`
                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                    echo $number
                    nameOut=""
                    rm listCurrMerge.txt
                    echo $fileName
                    for periodID in $periodList; do
                        echo $periodID
                        currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-pass$passNr-$runListName\_$number.root
                        if [ -f $currFile ]; then
                            nameOut+=$periodID
                            echo -e "$currFile\n" >> listCurrMerge.txt
                        fi
                    done
                    MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_LHC13$nameOut-pass$passNr-$runListName\_$number.root
                done
            fi
            ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr-$runListName\_*.root > filesForMerging.txt
            periodList=`echo -e "b\nc\nd\ne\nf"`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 3 "bla" 1
                echo $number
                rm listCurrMerge.txt
                echo $fileName
                for periodID in $periodList; do
                    echo $periodID
                    currFile=$OUTPUTDIR/GammaConvCalo_LHC13$periodID-pass$passNr-$runListName\_$number.root
                    if [ -f $currFile ]; then
                        echo -e "$currFile\n" >> listCurrMerge.txt
                    fi
                done
                fileList=`cat listCurrMerge.txt`
                for currFile in $fileList; do
                    mv $currFile $OUTPUTDIR/SingleFiles/
                done
            done
        done

        if [ $MERGEONMC -eq 1 ]; then
            echo -e "EMCandPCMGood\nPHOSGood" > runlistsToMerge.txt
            listsToMerge=`cat runlistsToMerge.txt`
            for runListName in $listsToMerge; do
                periodList=`echo -e "p1\np2\np3\np4"`
                ls $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1-$runListName\_*.root  > filesForMergingMC.txt
                filesForMerging=`cat filesForMergingMC.txt`
                listsToMerge=`cat runlistsToMerge.txt`
                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 6 "bla" 1
                    echo "number:"$number
                        rm listCurrMerge.txt
                        nameOut=""
                        for periodID in $periodList; do
                            echo $periodID
                            currFile=$OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_$periodID-$runListName\_$number.root
                            if [ -f $currFile ]; then
                                nameOut+="_"$periodID
                                echo -e "$currFile\n" >> listCurrMerge.txt
                            fi
                        done
                        MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix$nameOut-$runListName\_$number.root
                        fileList=`cat listCurrMerge.txt`
                        for currFile in $fileList; do
                            mv $currFile $OUTPUTDIR/SingleFiles/
                        done
                done
                ls $OUTPUTDIR/GammaConvCalo_MC_LHC13e7-$runListName\_*.root > filesForMerging.txt
                filesForMerging=`cat filesForMerging.txt`
                for fileName in $filesForMerging; do
                    echo $fileName
                    GetFileNumberMerging $fileName $((NSlashes-1)) 4 "bla" 1
                    echo $number
                    if [ -f  $OUTPUTDIR/GammaConvCalo_MC_LHC13e7-$runListName\_$number.root ] ; then
                        hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_p4_LHC13e7-$runListName_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_p4-$runListName\_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13e7_$number.root
                    fi
                done

                if [ $HAVELHC13b4fix == 1 ] && [ $HAVELHC13b4plus == 1 ]; then
                    ls $OUTPUTDIR/GammaConvCalo_MC_LHC13b4fix-$runListName\_*.root > filesForMerging.txt
                    filesForMerging=`cat filesForMerging.txt`
                    periodList=`echo -e "LHC13b4fix\nLHC13b4plus"`
                    for fileName in $filesForMerging; do
                        echo $fileName
                        GetFileNumberMerging $fileName $((NSlashes-1)) 4 "bla" 1
                        echo "number:"$number
                        rm listCurrMerge.txt
                        nameOut=""
                        for periodID in $periodList; do
                            echo $periodID
                            currFile=$OUTPUTDIR/GammaConvCalo_MC_$periodID-$runListName\_$number.root
                            if [ -f $currFile ]; then
                                nameOut+="_"$periodID
                                echo -e "$currFile\n" >> listCurrMerge.txt
                            fi
                        done
                        MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/GammaConvCalo_MC$nameOut-$runListName\_$number.root
                        fileList=`cat listCurrMerge.txt`
                        for currFile in $fileList; do
                            mv $currFile $OUTPUTDIR/SingleFiles/
                        done

                        binNumbersJJ=`cat runlists/binsJetJetLHC13b4.txt`
                        for binNumber in $binNumbersJJ; do
                            echo $binNumber
                            nameOut=""
                            rm listCurrMerge.txt
                            for periodID in $periodList; do
#                                 echo $periodID
                                currFile=$OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC_$periodID-$binNumber\_$runListName\_$number.root
#                                 echo $currFile
                                if [ -f $currFile ]; then
                                    nameOut+="_"$periodID
                                    echo -e "$currFile\n" >> listCurrMerge.txt
                                fi
                            done
                            MergeAccordingToList listCurrMerge.txt $OUTPUTDIR/JJMCSingleBins/GammaConvCalo_MC$nameOut-$binNumber\_$runListName\_$number.root
                        done

                    done
                fi
            done
        fi
        echo "--> Merging done"
    fi
else
    if [ $HAVELHC13b == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13b";
        rm $OUTPUTDIR_LHC13b/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13c == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13c";
        rm $OUTPUTDIR_LHC13c/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13d == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13d";
        rm $OUTPUTDIR_LHC13d/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13e == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13e";
        rm $OUTPUTDIR_LHC13e/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13f == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13f";
        rm $OUTPUTDIR_LHC13f/*/GammaConvCalo_*.root
    fi

    if [ $HAVELHC13b2efixp1 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13b2efix1";
        rm $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13b2efixp2 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13b2efix2";
        rm $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13b2efixp3 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13b2efix3";
        rm $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13b2efixp4 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13b2efix4";
        rm $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvCalo_*.root
    fi
    if [ $HAVELHC13e7 == 1 ]; then
        echo "removing all GammaConvCalo files in runFolders for LHC13e7";
        rm $OUTPUTDIR_LHC13e7/*/GammaConvCalo_*.root
    fi
fi
