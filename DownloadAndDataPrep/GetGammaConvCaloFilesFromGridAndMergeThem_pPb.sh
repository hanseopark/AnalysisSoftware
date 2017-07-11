#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash

# This script has to be run with "bash"

# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1
SINGLERUN=1
SEPARATEON=1
MERGEONSINGLEData=1
MERGEONSINGLEMC=1
SPECIALMERGE=0
ISADDDOWNLOAD=0
ADDDOWNLOADALREADY=0
CLEANUPADDMERGE=1
CLEANUP=1
CLEANUPMAYOR=$2

function GetFileNumberList()
{
    echo $1
    ls $1/GammaConvCalo_*.root > filesTemp.txt
    cat filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm fileNumbers.txt
    for fileName in $fileNumbers; do
        number=`echo $fileName  | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $3
    cat $3
}

function SeparateCutsIfNeeded()
{
    if [ -f $1\_A.root ]; then
        echo "separated file $1.root already"
        if [ $CLEANUP == 1 ]; then
            echo "removing $1.root as already separated"
            rm $1.root
        fi
    else
        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$1.root\"\,\"$1\"\,$2\,kTRUE\)
        if [ $CLEANUP == 1 ]; then
            if [ -f $1\_A.root ]; then
                echo "removing $1.root after successful separation"
                rm $1.root
            fi
        fi
    fi
}

function CopyFileIfNonExisitent()
{
    if [ $DOWNLOADON == 1 ]; then
        if [ -f $1/root_archive.zip ] && [ -s $1/root_archive.zip ]; then
            echo "$1/root_archive.zip exists";
        else
            mkdir -p $1
            alien_cp alien:$2/root_archive.zip file:$1/
        fi
        unzip -u $1/root_archive.zip -d $1/
    fi

    if [ $SEPARATEON == 1 ] && [ $ISADDDOWNLOAD == 0 ]; then
        GetFileNumberList $1 $3 fileNumbers2.txt
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaConvCalo_$fileNumber 2
        done;
    fi
}

function ChangeStructureIfNeeded()
{
    if [[ $1 == *"Basic.root"* ]]; then
        echo "Nothing to be done"
    else
        number1=`echo $1  | cut -d "/" -f $3 | cut -d "_" -f 2 | cut -d "." -f1`
        number2=`echo $1  | cut -d "/" -f $3 | cut -d "_" -f 3 | cut -d "." -f1`
        if [ -z "$number2" ]; then
            number=$number1
        else
            echo $number2
            number=$number1\_$number2
        fi
        echo $number
        cp $2/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_$4\_$number.root

        if [ $MERGEONSINGLEMC == 1 ]; then
        if [ -f $2/inputFailedStage1Merge/GammaConvCalo_$number.root ]; then
            echo "adding failed merges"
            hadd -f $OUTPUTDIR/GammaConvCalo_$4\_$number.root $2/GammaConvCalo_$number.root $2/inputFailedStage1Merge/GammaConvCalo_$number.root
        fi
        fi
        if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log ]; then
            echo "nothing to be done";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log\"\,2\)
        fi
    fi
}

function GetFileNumberMerging()
{
    echo $1
    NCurrSub=$3
    number1=`echo $1  | cut -d "/" -f $2 | cut -d "_" -f $NCurrSub | cut -d "." -f1`
    number2=`echo $1  | cut -d "/" -f $2 | cut -d "_" -f $((NCurrSub+1)) | cut -d "." -f1`
    if [ -z "$number2" ]; then
        number=$number1
    else
        echo $number2
        number=$number1\_$number2
    fi
}

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


# check if train configuration has actually been given
HAVELHC13b=1
HAVELHC13c=1
HAVELHC13d=1
HAVELHC13e=1
HAVELHC13f=1
HAVELHC13b2efixp1=1
HAVELHC13b2efixp2=1
HAVELHC13b2efixp3=1
HAVELHC13b2efixp4=1
HAVELHC13e7=1

# default trainconfigurations
LHC13bData="";
LHC13cData=""; #ESD
LHC13dData=""; #ESD
LHC13eData=""; #ESD
LHC13fData=""; #ESD
LHC13e7MC="";
LHC13b2_efix_p1MC="";
LHC13b2_efix_p2MC="";
LHC13b2_efix_p3MC="" ;
LHC13b2_efix_p4MC="";

passNr="2";

NSlashes=10;

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pPb
    NSlashes=8
    NSlashes2=7
    NSlashes3=9
    NSlashes4=10
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
    NSlashes=9;
fi


# TRAINDIR=Legotrain-Sys-ConvCalo
# LHC13bData="527";
# LHC13cData="528";
# LHC13dData="529";
# LHC13eData="530";
# LHC13fData="531";

# LHC13bData="532";
# LHC13cData="533";

# LHC13bData="534";
# LHC13cData="535";

# LHC13e7MC="708";
# LHC13b2_efix_p1MC="704";
# LHC13b2_efix_p2MC="705";
# LHC13b2_efix_p3MC="706" ;
# LHC13b2_efix_p4MC="707";

# LHC13b2_efix_p1MC="709";
# LHC13b2_efix_p2MC="710";
# LHC13b2_efix_p3MC="711" ;
# LHC13b2_efix_p4MC="712";

# LHC13b2_efix_p1MC="713";
# LHC13b2_efix_p2MC="714";
# LHC13b2_efix_p3MC="715" ;
# LHC13b2_efix_p4MC="716";

# LHC13b2_efix_p1MC="717";
# LHC13b2_efix_p2MC="718";
# LHC13b2_efix_p3MC="719" ;
# LHC13b2_efix_p4MC="720";

# TRAINDIR=Legotrain-Sys-ConvCalo-100ns
# LHC13bData="541";
# LHC13cData="542";

# TRAINDIR=Legotrain-Sys-ConvCalo-200ns
# LHC13bData="543";
# LHC13cData="544";

# TRAINDIR=Legotrain-vAN20160901-sys-PCMEMC
# LHC13bData="563";
# LHC13cData="564";

# LHC13b2_efix_p1MC="763";
# LHC13b2_efix_p2MC="767";
# LHC13b2_efix_p3MC="768";
# LHC13b2_efix_p4MC="769";

# LHC13b2_efix_p1MC="770";
# LHC13b2_efix_p2MC="771";
# LHC13b2_efix_p3MC="772";
# LHC13b2_efix_p4MC="773";

# LHC13b2_efix_p1MC="774";
# LHC13b2_efix_p2MC="775";
# LHC13b2_efix_p3MC="776";
# LHC13b2_efix_p4MC="777";

# TRAINDIR=Legotrain-vAN20161211-sys-PCMEMC
# LHC13bData="571";
# LHC13cData="572";
# LHC13dData="573";
# LHC13eData="574";
# LHC13fData="575";
#
# LHC13b2_efix_p1MC="807";
# LHC13b2_efix_p2MC="808";
# LHC13b2_efix_p3MC="809";
# LHC13b2_efix_p4MC="810";
# LHC13e7MC="806";

# TRAINDIR=Legotrain-vAN20170411-dirGamma
# LHC13bData="593";
# LHC13cData="595";
# # LHC13dData="573";
# # LHC13eData="574";
# # LHC13fData="575";
#
# # LHC13b2_efix_p1MC="889";
# # LHC13b2_efix_p2MC="890";
# # LHC13b2_efix_p3MC="891";
# # LHC13b2_efix_p4MC="892";
# # LHC13e7MC="893";
# LHC13b2_efix_p1MC="894";
# LHC13b2_efix_p2MC="895";
# LHC13b2_efix_p3MC="896";
# LHC13b2_efix_p4MC="897";
# LHC13e7MC="898";

# TRAINDIR=Legotrain-vAN20170417-Weighting
# LHC13bData="599"; #pass 3
# LHC13cData="601"; #pass 2
# LHC13b2_efix_p1MC="904";
# LHC13b2_efix_p2MC="905";
# LHC13b2_efix_p3MC="906";
# LHC13b2_efix_p4MC="907";
# LHC13e7MC="908";

# TRAINDIR=Legotrain-vAN20170406-PHOSrelated
# LHC13bData="593"; #pass 3
# LHC13cData="595"; #pass 2
# LHC13b2_efix_p1MC="894";
# LHC13b2_efix_p2MC="895";
# LHC13b2_efix_p3MC="896";
# LHC13b2_efix_p4MC="897";
# LHC13e7MC="898";

# TRAINDIR=Legotrain-pPb-CalibMike
# LHC13bData="448"; #pass 3
# LHC13cData="449"; #pass 2
# LHC13dData="450"; #pass 2
# LHC13eData="451"; #pass 2
# LHC13b2_efix_p1MC="596";
# LHC13b2_efix_p2MC="597";
# LHC13b2_efix_p3MC="598";
# LHC13b2_efix_p4MC="599";

# TRAINDIR=Legotrain-vAN20170501-EMCwoCalib
# LHC13bData="615"; #pass 3
# LHC13cData="616"; #pass 2
# LHC13dData="617"; #pass 2
# LHC13eData="618"; #pass 2
# LHC13b2_efix_p1MC="909";
# LHC13b2_efix_p2MC="911";
# LHC13b2_efix_p3MC="913";
# LHC13b2_efix_p4MC="915";
# LHC13e7MC="917";

# TRAINDIR=Legotrain-vAN20170506-EMCrecalibAndWithMultWeights
# LHC13bData="625"; #pass 3
# LHC13cData="626"; #pass 2
# LHC13dData="617"; #pass 2
# LHC13eData="618"; #pass 2
# LHC13b2_efix_p1MC="920";
# LHC13b2_efix_p2MC="921";
# LHC13b2_efix_p3MC="922";
# LHC13b2_efix_p4MC="923";
# LHC13e7MC="924";
# LHC13b2_efix_p1MC="925";
# LHC13b2_efix_p2MC="926";
# LHC13b2_efix_p3MC="927";
# LHC13b2_efix_p4MC="928";
# LHC13e7MC="929";
#
# TRAINDIR=Legotrain-vAN20170511-EMCrecalibIte1
# LHC13b2_efix_p1MC="930";
# LHC13b2_efix_p2MC="931";
# LHC13b2_efix_p3MC="932";
# LHC13b2_efix_p4MC="933";
# LHC13e7MC="934";

# TRAINDIR=Legotrain-vAN20170511-EMCrecalibIte2
# LHC13b2_efix_p1MC="935";
# LHC13b2_efix_p2MC="936";
# LHC13b2_efix_p3MC="937";
# LHC13b2_efix_p4MC="938";
# LHC13e7MC="945";
# LHC13b2_efix_p1MC="940";
# LHC13b2_efix_p2MC="941";
# LHC13b2_efix_p3MC="942";
# LHC13b2_efix_p4MC="943";
# LHC13e7MC="944";

# TRAINDIR=Legotrain-vAN20170511-EMCrecalibIte3
# LHC13b2_efix_p1MC="946";
# LHC13b2_efix_p2MC="947";
# LHC13b2_efix_p3MC="948";
# LHC13b2_efix_p4MC="949";
# LHC13e7MC="950";

# TRAINDIR=Legotrain-vAN20170519-EMCrecalibIte4
# LHC13bData="629"; #pass 3
# LHC13cData="630"; #pass 2
# LHC13b2_efix_p1MC="951";
# LHC13b2_efix_p2MC="953";
# LHC13b2_efix_p3MC="955";
# LHC13b2_efix_p4MC="957";
# LHC13e7MC="959";
# LHC13b2_efix_p1MC="952";
# LHC13b2_efix_p2MC="954";
# LHC13b2_efix_p3MC="956";
# LHC13b2_efix_p4MC="958";
# LHC13e7MC="960";


TRAINDIR=Legotrain-vAN20170525FF-newDefaultPlusSys
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
LHC13b2_efix_p1MC="1026";
LHC13b2_efix_p2MC="1046";
LHC13b2_efix_p3MC="1047";
LHC13b2_efix_p4MC="1048";


OUTPUTDIR=$BASEDIR/$TRAINDIR
mkdir -p $OUTPUTDIR/CutSelections

if [ "$LHC13bData" == "" ]; then
    HAVELHC13b=0;
fi
if [ "$LHC13cData" = "" ]; then
    HAVELHC13c=0;
fi
if [ "$LHC13dData" = "" ]; then
    HAVELHC13d=0;
fi
if [ "$LHC13eData" = "" ]; then
    HAVELHC13e=0;
fi
if [ "$LHC13fData" = "" ]; then
    HAVELHC13f=0;
fi

if [ "$LHC13b2_efix_p1MC" = "" ]; then
    HAVELHC13b2efixp1=0;
fi
if [ "$LHC13b2_efix_p2MC" = "" ]; then
    HAVELHC13b2efixp2=0;
fi
if [ "$LHC13b2_efix_p3MC" = "" ]; then
    HAVELHC13b2efixp3=0;
fi
if [ "$LHC13b2_efix_p4MC" = "" ]; then
    HAVELHC13b2efixp4=0;
fi
if [ "$LHC13e7MC" = "" ]; then
    HAVELHC13e7=0;
fi

if [ $HAVELHC13b == 1 ]; then
    LHC13bData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13bData\_`
    if [ "$LHC13bData" == "" ]; then
        HAVELHC13b=0;
    else
        OUTPUTDIR_LHC13b=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13bData
    fi
fi
if [ $HAVELHC13c == 1 ]; then
    LHC13cData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13cData\_`
    if [ "$LHC13cData" == "" ]; then
        HAVELHC13c=0;
    else
        OUTPUTDIR_LHC13c=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13cData
    fi
fi
if [ $HAVELHC13d == 1 ]; then
    LHC13dData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13dData\_`
    if [ "$LHC13dData" == "" ]; then
        HAVELHC13d=0;
    else
        OUTPUTDIR_LHC13d=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13dData
    fi
fi
if [ $HAVELHC13e == 1 ]; then
    LHC13eData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13eData\_`
    if [ "$LHC13eData" == "" ]; then
        HAVELHC13e=0;
    else
        OUTPUTDIR_LHC13e=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13eData
    fi
fi
if [ $HAVELHC13f == 1 ]; then
    LHC13fData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13fData\_`
    if [ "$LHC13fData" == "" ]; then
        HAVELHC13f=0;
    else
        OUTPUTDIR_LHC13f=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13fData
    fi
fi

if [ $HAVELHC13b2efixp1 == 1 ]; then
    LHC13b2_efix_p1MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p1MC\_`
    if [ "$LHC13b2_efix_p1MC" == "" ]; then
        HAVELHC13b2efixp1=0;
    else
        OUTPUTDIR_LHC13b2_efix_p1=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p1MC
    fi
fi
if [ $HAVELHC13b2efixp2 == 1 ]; then
    LHC13b2_efix_p2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p2MC\_`
    if [ "$LHC13b2_efix_p2MC" == "" ]; then
        HAVELHC13b2efixp2=0;
    else
        OUTPUTDIR_LHC13b2_efix_p2=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p2MC
    fi
fi
if [ $HAVELHC13b2efixp3 == 1 ]; then
    LHC13b2_efix_p3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p3MC\_`
    if [ "$LHC13b2_efix_p3MC" == "" ]; then
        HAVELHC13b2efixp3=0;
    else
        OUTPUTDIR_LHC13b2_efix_p3=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p3MC
    fi
fi
if [ $HAVELHC13b2efixp4 == 1 ]; then
    LHC13b2_efix_p4MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p4MC\_`
    if [ "$LHC13b2_efix_p4MC" == "" ]; then
        HAVELHC13b2efixp4=0;
    else
        OUTPUTDIR_LHC13b2_efix_p4=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p4MC
    fi
fi
if [ $HAVELHC13e7 == 1 ]; then
    LHC13e7MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13e7MC\_`
    if [ "$LHC13e7MC" == "" ]; then
        HAVELHC13e7=0;
    else
        OUTPUTDIR_LHC13e7=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13e7MC
    fi
fi

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
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13bData/merge_runlist_1" $NSlashes
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
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13c "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge_runlist_1" $NSlashes
        fi

    fi
    if [ $HAVELHC13d == 1 ]; then
        echo "downloading LHC13d"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13d "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13dData/merge" $NSlashes
    fi
    if [ $HAVELHC13e == 1 ]; then
        echo "downloading LHC13e"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13e "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13eData/merge_runlist_2" $NSlashes
    fi
    if [ $HAVELHC13f == 1 ]; then
        echo "downloading LHC13f"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13f "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13fData/merge" $NSlashes
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
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/merge" $NSlashes
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
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p2/mergedAllConvCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/merge" $NSlashes
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
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/merge" $NSlashes
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
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/merge" $NSlashes
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
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge" $NSlashes
        fi
        ParseJDLFilesDownloadAndMerge missingMergesLHC13e7.txt  $OUTPUTDIR_LHC13e7
    fi


    if [ $HAVELHC13b == 1 ]; then
        ls $OUTPUTDIR_LHC13b/GammaConvCalo_*.root > fileLHC13b.txt
        fileNumbers=`cat fileLHC13b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass$passNr"
        done;
    fi

    if [ $HAVELHC13c == 1 ]; then
        ls $OUTPUTDIR_LHC13c/GammaConvCalo_*.root > fileLHC13c.txt
        fileNumbers=`cat fileLHC13c.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass$passNr"
        done;

    fi

    if [ $HAVELHC13d == 1 ]; then
        ls $OUTPUTDIR_LHC13d/GammaConvCalo_*.root > fileLHC13d.txt
        fileNumbers=`cat fileLHC13d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13d $NSlashes "LHC13d-pass$passNr"
        done;
    fi

    if [ $HAVELHC13e == 1 ]; then
        ls $OUTPUTDIR_LHC13e/GammaConvCalo_*.root > fileLHC13e.txt
        fileNumbers=`cat fileLHC13e.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13e $NSlashes "LHC13e-pass$passNr"
        done;
    fi

    if [ $HAVELHC13f == 1 ]; then
        ls $OUTPUTDIR_LHC13f/GammaConvCalo_*.root > fileLHC13f.txt
        fileNumbers=`cat fileLHC13f.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13f $NSlashes "LHC13f-pass$passNr"
        done;
    fi

    if [ $HAVELHC13b2efixp1 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvCalo_*.root > fileLHC13b2efixp1.txt
        fileNumbers=`cat fileLHC13b2efixp1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p1 $NSlashes "MC_LHC13b2_efix_p1"
        done;
    fi

    if [ $HAVELHC13b2efixp2 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvCalo_*.root > fileLHC13b2efixp2.txt
        fileNumbers=`cat fileLHC13b2efixp2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p2 $NSlashes "MC_LHC13b2_efix_p2"
        done;
    fi

    if [ $HAVELHC13b2efixp3 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvCalo_*.root > fileLHC13b2efixp3.txt
        fileNumbers=`cat fileLHC13b2efixp3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p3 $NSlashes "MC_LHC13b2_efix_p3"
        done;
    fi

    if [ $HAVELHC13b2efixp4 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvCalo_*.root > fileLHC13b2efixp4.txt
        fileNumbers=`cat fileLHC13b2efixp4.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p4 $NSlashes "MC_LHC13b2_efix_p4"
        done;
    fi

    if [ $HAVELHC13e7 == 1 ]; then
        ls $OUTPUTDIR_LHC13e7/GammaConvCalo_*.root > fileLHC13e7.txt
        fileNumbers=`cat fileLHC13e7.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13e7 $NSlashes "MC_LHC13e7"
        done;
    fi

    if [ $MERGEON == 1 ]; then
#         ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_*.root > filesForMerging.txt
#         filesForMerging=`cat filesForMerging.txt`
#         for fileName in $filesForMerging; do
#             echo $fileName
#             GetFileNumberMerging $fileName $((NSlashes-1)) 3
#             echo $number
#             ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root
#             ls $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root
#             if [ -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root ] ; then
#                 hadd -f $OUTPUTDIR/GammaConvCalo_LHC13bc-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root
#             fi
#         done
#
#         if [ $HAVELHC13d == 1 ] && [ $HAVELHC13e == 1 ]; then
#             ls $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_*.root > filesForMerging.txt
#             filesForMerging=`cat filesForMerging.txt`
#             for fileName in $filesForMerging; do
#                 echo $fileName
#                 GetFileNumberMerging $fileName $((NSlashes-1)) 3
#                 echo $number
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
#                 if [ -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root ] &&  [ -f $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root ] ; then
#                     hadd -f $OUTPUTDIR/GammaConvCalo_LHC13bcde-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
#                 fi
#
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
#                 if  [ -f $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root ] ; then
#                     hadd -f $OUTPUTDIR/GammaConvCalo_LHC13cde-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
#                 fi
#
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root
#                 ls $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
#                 if  [ -f $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root ] ; then
#                     hadd -f $OUTPUTDIR/GammaConvCalo_LHC13de-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
#                 fi
#             done
#         fi

        ls $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_*.root | sed 's/.*p1_p2.*//' | sed '/^$/d' > filesForMergingMC.txt
        filesForMerging=`cat filesForMergingMC.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 6
            echo "number:"$number
            if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p4_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p2_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p3_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p4_$number.root
            fi
        done
        ls $OUTPUTDIR/GammaConvCalo_MC_LHC13e7_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 4
            echo $number
            if [ -f  $OUTPUTDIR/GammaConvCalo_MC_LHC13e7_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_p4_LHC13e7_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13e7_$number.root
            fi
        done
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
