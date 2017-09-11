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
SEPARATEON=0
MERGEONSINGLEData=1
MERGEONSINGLEMC=1
CLEANUP=1
CLEANUPMAYOR=$2

function GetFileNumberList()
{
    ls $1/GammaCalo_*.root > filesTemp.txt
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

    if [ $SEPARATEON == 1 ] ; then
        GetFileNumberList $1 $3 fileNumbers2.txt
        fileNumbers=`cat fileNumbers2.txt`
        echo fileNumbers
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaCalo_$fileNumber 4
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
        cp $2/GammaCalo_$number.root $OUTPUTDIR/GammaCalo_$4\_$number.root

        if [ -f $2/inputFailedStage1Merge/GammaCalo_$number.root ]; then
            echo "adding failed merges"
            hadd -f $OUTPUTDIR/GammaCalo_$4\_$number.root $2/GammaCalo_$number.root $2/inputFailedStage1Merge/GammaCalo_$number.root
        fi

        if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log ]; then
            echo "nothing to be done";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log\"\,4\)
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
HAVETOBUILDMC=0
HAVETOBUILDData=0
# default trainconfigurations
LHC13bData="";
LHC13cData=""; #ESD
LHC13dData=""; #ESD
LHC13eData=""; #ESD
LHC13fData=""; #ESD
LHC13bcData="";
LHC13e7MC="";
LHC13b2_efix_p1MC="";
LHC13b2_efix_p2MC="";
LHC13b2_efix_p3MC="" ;
LHC13b2_efix_p4MC="";
LHC13b2_efix_MC="";

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

# TRAINDIR=Legotrain-vAN20170411-dirGamma
# LHC13bData="593";
# LHC13cData="595";
# LHC13dData="573";
# LHC13eData="574";
# LHC13fData="575";

# LHC13b2_efix_p1MC="889";
# LHC13b2_efix_p2MC="890";
# LHC13b2_efix_p3MC="891";
# LHC13b2_efix_p4MC="892";
# LHC13e7MC="893";
# LHC13b2_efix_p1MC="894";
# LHC13b2_efix_p2MC="895";
# LHC13b2_efix_p3MC="896";
# LHC13b2_efix_p4MC="897";
# LHC13e7MC="898";

# TRAINDIR=Legotrain-vAN20170416-Weighting
# LHC13bData="599"; #pass 3
# LHC13cData="601"; #pass 2
# LHC13b2_efix_p1MC="899";
# LHC13b2_efix_p2MC="900";
# LHC13b2_efix_p3MC="901";
# LHC13b2_efix_p4MC="902";
# LHC13e7MC="903";

# TRAINDIR=Legotrain-vAN20170417-Weighting
# LHC13bData="599"; #pass 3
# LHC13cData="601"; #pass 2
# LHC13b2_efix_p1MC="904";
# LHC13b2_efix_p2MC="905";
# LHC13b2_efix_p3MC="906";
# LHC13b2_efix_p4MC="907";
# LHC13e7MC="908";

# TRAINDIR=Legotrain-pPb-CalibMike
# LHC13bData="448"; #pass 3
# LHC13cData="449"; #pass 2
# LHC13DData="450"; #pass 2
# LHC13EData="451"; #pass 2
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

# TRAINDIR=Legotrain-vAN20170525FF-newDefaultPlusSys
# LHC13bData="640"; #pass 3
# LHC13cData="641"; #pass 2
# LHC13b2_efix_p1MC="973";
# LHC13b2_efix_p2MC="974";
# LHC13b2_efix_p3MC="975";
# LHC13b2_efix_p4MC="976";
# LHC13e7MC="977";
# LHC13b2_efix_p1MC="978";
# LHC13b2_efix_p2MC="982";
# LHC13b2_efix_p3MC="980";
# LHC13b2_efix_p4MC="981";
# LHC13b2_efix_p1MC="1026";
# LHC13b2_efix_p2MC="1046";
# LHC13b2_efix_p3MC="1047";
# LHC13b2_efix_p4MC="1048";
# LHC13bData="660"; #pass 3
# LHC13cData="661"; #pass 2
# LHC13b2_efix_MC="1065"
# LHC13b2_efix_MC="1066"
# LHC13b2_efix_p1MC="child_1";
# LHC13b2_efix_p2MC="child_2";
# LHC13b2_efix_p3MC="child_3";
# LHC13b2_efix_p4MC="child_4";

TRAINDIR=Legotrain-vAN20170905-dirGamma
LHC13bcData="672"; #pass 2
LHC13bData="child_1"; #pass 3
LHC13cData="child_2"; #pass 2
LHC13b2_efix_MC="1071"
LHC13b2_efix_p1MC="child_1";
LHC13b2_efix_p2MC="child_2";
LHC13b2_efix_p3MC="child_3";
LHC13b2_efix_p4MC="child_4";

OUTPUTDIR=$BASEDIR/$TRAINDIR


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
if [ "$LHC13bcData" != "" ]; then
    HAVETOBUILDData=1;
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
if [ "$LHC13b2_efix_MC" != "" ]; then
    HAVETOBUILDMC=1;
fi
if [ "$LHC13e7MC" = "" ]; then
    HAVELHC13e7=0;
fi

if [ $HAVELHC13b == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC13bData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13bcData\_ | grep $LHC13bData`
    else
        LHC13bData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13bData\_`
    fi
    if [ "$LHC13bData" == "" ]; then
        HAVELHC13b=0;
    else
        OUTPUTDIR_LHC13b=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13bData
    fi
    echo $OUTPUTDIR_LHC13b
fi
if [ $HAVELHC13c == 1 ]; then
    if [ $HAVETOBUILDData == 1 ]; then
        LHC13cData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13bcData\_ | grep $LHC13cData`
    else
        LHC13cData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/ | grep $LHC13cData\_`
    fi
    if [ "$LHC13cData" == "" ]; then
        HAVELHC13c=0;
    else
        OUTPUTDIR_LHC13c=$BASEDIR/$TRAINDIR/GA_pPb-$LHC13cData
    fi
    echo $OUTPUTDIR_LHC13c
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
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC13b2_efix_p1MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_MC\_ | grep $LHC13b2_efix_p1MC `
    else
        LHC13b2_efix_p1MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p1MC\_`
    fi
    if [ "$LHC13b2_efix_p1MC" == "" ]; then
        HAVELHC13b2efixp1=0;
    else
        OUTPUTDIR_LHC13b2_efix_p1=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p1MC
    fi
    echo $OUTPUTDIR_LHC13b2_efix_p1
fi
if [ $HAVELHC13b2efixp2 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC13b2_efix_p2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_MC\_ | grep $LHC13b2_efix_p2MC `
    else
        LHC13b2_efix_p2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p2MC\_`
    fi
    if [ "$LHC13b2_efix_p2MC" == "" ]; then
        HAVELHC13b2efixp2=0;
    else
        OUTPUTDIR_LHC13b2_efix_p2=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p2MC
    fi
    echo $OUTPUTDIR_LHC13b2_efix_p2
fi
if [ $HAVELHC13b2efixp3 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC13b2_efix_p3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_MC\_ | grep $LHC13b2_efix_p3MC `
    else
        LHC13b2_efix_p3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p3MC\_`
    fi
    if [ "$LHC13b2_efix_p3MC" == "" ]; then
        HAVELHC13b2efixp3=0;
    else
        OUTPUTDIR_LHC13b2_efix_p3=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p3MC
    fi
    echo $OUTPUTDIR_LHC13b2_efix_p3
fi
if [ $HAVELHC13b2efixp4 == 1 ]; then
    if [ $HAVETOBUILDMC == 1 ]; then
        LHC13b2_efix_p4MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_MC\_ | grep $LHC13b2_efix_p4MC `
    else
        LHC13b2_efix_p4MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13b2_efix_p4MC\_`
    fi
    if [ "$LHC13b2_efix_p4MC" == "" ]; then
        HAVELHC13b2efixp4=0;
    else
        OUTPUTDIR_LHC13b2_efix_p4=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13b2_efix_p4MC
    fi
    echo $OUTPUTDIR_LHC13b2_efix_p4
fi
if [ $HAVELHC13e7 == 1 ]; then
    LHC13e7MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/ | grep $LHC13e7MC\_`
    if [ "$LHC13e7MC" == "" ]; then
        HAVELHC13e7=0;
    else
        OUTPUTDIR_LHC13e7=$BASEDIR/$TRAINDIR/GA_pPb_MC-$LHC13e7MC
    fi
fi

mkdir -p $OUTPUTDIR/CutSelections


if [ $CLEANUPMAYOR == 0 ]; then
    if [ $HAVELHC13b == 1 ]; then
        echo "downloading LHC13b"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b/$runNumber "/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/GA_pPb/$LHC13bData" $NSlashes3
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC13b/mergedAllCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b.txt`
                ls $OUTPUTDIR_LHC13b/$firstrunNumber/GammaCalo_*.root > fileLHC13b.txt
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
                    hadd -f $OUTPUTDIR_LHC13b/GammaCalo_$number.root $OUTPUTDIR_LHC13b/*/GammaCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b/mergedAllCalo.txt
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
            if [ $MERGEONSINGLEData == 1 ]  && [ ! -f $OUTPUTDIR_LHC13c/mergedAllCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13c.txt`
                ls $OUTPUTDIR_LHC13c/$firstrunNumber/GammaCalo_*.root > fileLHC13c.txt
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
                    hadd -f $OUTPUTDIR_LHC13c/GammaCalo_$number.root $OUTPUTDIR_LHC13c/*/GammaCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13c/mergedAllCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13c "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge_runlist_1" $NSlashes
        fi
    fi
    if [ $HAVELHC13d == 1 ]; then
        echo "downloading LHC13d"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13d "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13dData/merge"
    fi
    if [ $HAVELHC13e == 1 ]; then
        echo "downloading LHC13e"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13e "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13eData/merge_runlist_2"
    fi
    if [ $HAVELHC13f == 1 ]; then
        echo "downloading LHC13f"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13f "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13fData/merge"
    fi


    if [ $HAVELHC13b2efixp1 == 1 ]; then
        echo "downloading LHC13b2_efix_p1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1/$runNumber "/alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p1/mergedAllCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix1.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p1/$firstrunNumber/GammaCalo_*.root > fileLHC13b2_efix_p1.txt
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p1/GammaCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p1/*/GammaCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p1/mergedAllCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/merge" $NSlashes
        fi

    fi
    if [ $HAVELHC13b2efixp2 == 1 ]; then
        echo "downloading LHC13b2_efix_p2"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix2.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2/$runNumber "/alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p2/mergedAllCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix2.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p2/$firstrunNumber/GammaCalo_*.root > fileLHC13b2_efix_p2.txt
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p2/GammaCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p2/*/GammaCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p2/mergedAllCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/merge" $NSlashes
        fi
    fi
    if [ $HAVELHC13b2efixp3 == 1 ]; then
        echo "downloading LHC13b2_efix_p3"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix3.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3/$runNumber "/alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p3/mergedAllCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix3.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p3/$firstrunNumber/GammaCalo_*.root > fileLHC13b2_efix_p3.txt
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p3/GammaCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p3/*/GammaCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p3/mergedAllCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/merge" $NSlashes
        fi
    fi
    if [ $HAVELHC13b2efixp4 == 1 ]; then
        echo "downloading LHC13b2_efix_p4"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix4.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4/$runNumber "/alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p4/mergedAllCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix4.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p4/$firstrunNumber/GammaCalo_*.root > fileLHC13b2_efix_p4.txt
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p4/GammaCalo_$number.root $OUTPUTDIR_LHC13b2_efix_p4/*/GammaCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p4/mergedAllCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/merge" $NSlashes
        fi
    fi
    if [ $HAVELHC13e7 == 1 ]; then
        echo "downloading LHC13e7"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13e7.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7/$runNumber "/alice/sim/2013/LHC13e7/$runNumber/PWGGA/GA_pPb_MC/$LHC13e7MC" $NSlashes3
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13e7/mergedAllCalo.txt ]; then
                firstrunNumber=`head -n1 runlists/runNumbersLHC13e7.txt`
                ls $OUTPUTDIR_LHC13e7/$firstrunNumber/GammaCalo_*.root > fileLHC13e7.txt
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
                    hadd -f $OUTPUTDIR_LHC13e7/GammaCalo_$number.root $OUTPUTDIR_LHC13e7/*/GammaCalo_$number.root
                done;
                echo "done" > $OUTPUTDIR_LHC13e7/mergedAllCalo.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge" $NSlashes
        fi
    fi

    if [ $HAVELHC13b == 1 ]; then
        ls $OUTPUTDIR_LHC13b/GammaCalo_*.root > fileLHC13b.txt
        fileNumbers=`cat fileLHC13b.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass$passNr"
        done;
    fi

    if [ $HAVELHC13c == 1 ]; then
        ls $OUTPUTDIR_LHC13c/GammaCalo_*.root > fileLHC13c.txt
        fileNumbers=`cat fileLHC13c.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass$passNr"
        done;
    fi

    if [ $HAVELHC13d == 1 ]; then
        ls $OUTPUTDIR_LHC13d/GammaCalo_*.root > fileLHC13d.txt
        fileNumbers=`cat fileLHC13d.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13d $NSlashes "LHC13d-pass$passNr"
        done;
    fi

    if [ $HAVELHC13e == 1 ]; then
        ls $OUTPUTDIR_LHC13e/GammaCalo_*.root > fileLHC13e.txt
        fileNumbers=`cat fileLHC13e.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13e $NSlashes "LHC13e-pass$passNr"
        done;
    fi

    if [ $HAVELHC13f == 1 ]; then
        ls $OUTPUTDIR_LHC13f/GammaCalo_*.root > fileLHC13f.txt
        fileNumbers=`cat fileLHC13f.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13f $NSlashes "LHC13f-pass$passNr"
        done;
    fi

    if [ $HAVELHC13b2efixp1 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p1/GammaCalo_*.root > fileLHC13b2efixp1.txt
        fileNumbers=`cat fileLHC13b2efixp1.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p1 $NSlashes "MC_LHC13b2_efix_p1"
        done;
    fi

    if [ $HAVELHC13b2efixp2 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p2/GammaCalo_*.root > fileLHC13b2efixp2.txt
        fileNumbers=`cat fileLHC13b2efixp2.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p2 $NSlashes "MC_LHC13b2_efix_p2"
        done;
    fi

    if [ $HAVELHC13b2efixp3 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p3/GammaCalo_*.root > fileLHC13b2efixp3.txt
        fileNumbers=`cat fileLHC13b2efixp3.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p3 $NSlashes "MC_LHC13b2_efix_p3"
        done;
    fi

    if [ $HAVELHC13b2efixp4 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p4/GammaCalo_*.root > fileLHC13b2efixp4.txt
        fileNumbers=`cat fileLHC13b2efixp4.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p4 $NSlashes "MC_LHC13b2_efix_p4"
        done;
    fi

    if [ $HAVELHC13e7 == 1 ]; then
        ls $OUTPUTDIR_LHC13e7/GammaCalo_*.root > fileLHC13e7.txt
        fileNumbers=`cat fileLHC13e7.txt`
        for fileName in $fileNumbers; do
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13e7 $NSlashes "MC_LHC13e7"
        done;
    fi


    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaCalo_LHC13b-pass$passNr\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            ls $OUTPUTDIR/GammaCalo_LHC13b-pass$passNr\_$number.root
            ls $OUTPUTDIR/GammaCalo_LHC13c-pass$passNr\_$number.root
            if [ -f $OUTPUTDIR/GammaCalo_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC13c-pass$passNr\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCalo_LHC13bc-pass$passNr\_$number.root $OUTPUTDIR/GammaCalo_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaCalo_LHC13c-pass$passNr\_$number.root
            fi
        done
#
#         ls $OUTPUTDIR/GammaCalo_LHC13d-pass$passNr\_*.root > filesForMerging.txt
#         filesForMerging=`cat filesForMerging.txt`
#         for fileName in $filesForMerging; do
#             echo $fileName
#             number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
#             echo $number
#             if [ -f $OUTPUTDIR/GammaCalo_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC13e-pass$passNr\_$number.root ] ; then
#                 hadd -f $OUTPUTDIR/GammaCalo_LHC13de-pass$passNr\_$number.root $OUTPUTDIR/GammaCalo_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaCalo_LHC13e-pass$passNr\_$number.root
#             fi
#         done
#
#         filesForMerging=`cat filesForMerging.txt`
#         for fileName in $filesForMerging; do
#             echo $fileName
#             number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 5 | cut -d "." -f1`
#             echo $number
#             if [ -f $OUTPUTDIR/GammaCalo_LHC13bc-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_LHC13de-pass$passNr\_$number.root ] ; then
#                 hadd -f $OUTPUTDIR/GammaCalo_LHC13bcde-pass$passNr\_$number.root $OUTPUTDIR/GammaCalo_LHC13bc-pass$passNr\_$number.root $OUTPUTDIR/GammaCalo_LHC13de-pass$passNr\_$number.root
#             fi
#         done

        ls $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_*.root | sed 's/.*p1_p2.*//' | sed '/^$/d' > filesForMergingMC.txt
        filesForMerging=`cat filesForMergingMC.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 6 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p2_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p3_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p4_$number.root ]; then
                hadd -f $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p2_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p3_$number.root $OUTPUTDIR/GammaCalo_MC_LHC13b2_efix_p4_$number.root
            fi
        done
    fi
else
    if [ $HAVELHC13b == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13b";
        rm $OUTPUTDIR_LHC13b/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13c == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13c";
        rm $OUTPUTDIR_LHC13c/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13d == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13d";
        rm $OUTPUTDIR_LHC13d/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13e == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13e";
        rm $OUTPUTDIR_LHC13e/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13f == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13f";
        rm $OUTPUTDIR_LHC13f/*/GammaCalo_*.root
    fi

    if [ $HAVELHC13b2efixp1 == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13b2efix1";
        rm $OUTPUTDIR_LHC13b2_efix_p1/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13b2efixp2 == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13b2efix2";
        rm $OUTPUTDIR_LHC13b2_efix_p2/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13b2efixp3 == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13b2efix3";
        rm $OUTPUTDIR_LHC13b2_efix_p3/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13b2efixp4 == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13b2efix4";
        rm $OUTPUTDIR_LHC13b2_efix_p4/*/GammaCalo_*.root
    fi
    if [ $HAVELHC13e7 == 1 ]; then
        echo "removing all GammaCalo files in runFolders for LHC13e7";
        rm $OUTPUTDIR_LHC13e7/*/GammaCalo_*.root
    fi
fi

