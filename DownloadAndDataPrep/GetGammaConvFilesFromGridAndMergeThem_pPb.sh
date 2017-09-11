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
MERGEONSINGLEMC=0
CLEANUP=1
CLEANUPMAYOR=$2
ENABLEFLOW=0

function GetFileNumberList()
{
    ls $1/GammaConvV1_*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm fileNumbers.txt
    for fileName in $fileNumbers; do
        number=`echo $fileName  | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $3
    cat $3
}

function GetFileNumberListFlow()
{
    ls $1/GammaConvFlow_*.root > filesTemp.txt
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
        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$1.root\"\,\"$1\"\,0\,kTRUE\)
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

    if [ $SEPARATEON == 1 ]; then
        GetFileNumberList $1 $3 fileNumbers2.txt
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaConvV1_$fileNumber
        done;
        if [ $ENABLEFLOW == 1 ]; then
            GetFileNumberListFlow $1 $3 fileNumbers2.txt
            fileNumbers=`cat fileNumbers2.txt`
            for fileNumber in $fileNumbers; do
                echo $fileNumber
                SeparateCutsIfNeeded $1/GammaConvFlow_$fileNumber
            done;
        fi
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
        cp $2/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_$4\_$number.root
        if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log ]; then
            echo "nothing to be done";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log\"\,0\)
        fi
   fi
}

function ChangeStructureIfNeededFlow()
{
    number1=`echo $1  | cut -d "/" -f $3 | cut -d "_" -f 2 | cut -d "." -f1`
    number2=`echo $1  | cut -d "/" -f $3 | cut -d "_" -f 3 | cut -d "." -f1`
    if [ -z "$number2" ]; then
        number=$number1
    else
        echo $number2
        number=$number1\_$number2
    fi
    echo $number
    cp $2/GammaConvFlow_$number.root $OUTPUTDIR/GammaConvFlow_$4\_$number.root
    if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaConvFlow_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaConvFlow_$4_$number.log ]; then
        echo "nothing to be done";
    else
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvFlow_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaConvFlow_$4_$number.log\"\,0\)
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

# check if train configuration has actually been given
HAVELHC13b=1
HAVELHC13b=1
HAVELHC13c=1
HAVELHC13d=1
HAVELHC13e=1
HAVELHC13f=1
HAVETOBUILDData=0
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
LHC13bcData="";
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
elif [ $1 = "passfeld" ]; then
    BASEDIR=~/work/Gridoutput/pPb
elif [ $1 = "passfeldMAF" ]; then
    BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then
    BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pPb
elif [ $1 = "pgonzales" ]; then
    BASEDIR=~/work/GridOutput
elif [ $1 = "pgonzalesGSI" ]; then
    BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pPb
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
    NSlashes=9;
fi



# TRAINDIR=Legotrain-vAN20170123-sys-PCM-dirGamma
# LHC13bData="583"; #pass 3
# LHC13cData="584"; #pass 2
# LHC13bData="582"; #pass 4
# LHC13cData="585"; #pass 4
# LHC13dData="573";
# LHC13eData="574";
# LHC13fData="575";

# LHC13b2_efix_p1MC="820";
# LHC13b2_efix_p2MC="821";
# LHC13b2_efix_p3MC="822";
# LHC13b2_efix_p4MC="824";
# LHC13e7MC="825";
# LHC13b2_efix_p1MC="826";
# LHC13b2_efix_p2MC="827";
# LHC13b2_efix_p3MC="828";
# LHC13b2_efix_p4MC="829";
# LHC13e7MC="830";

# TRAINDIR=Legotrain-vAN20170205-sys-PCM-dirGammaWoWeight
# # LHC13b2_efix_p1MC="831";
# # LHC13b2_efix_p2MC="832";
# # LHC13b2_efix_p3MC="833";
# # LHC13b2_efix_p4MC="834";
# LHC13b2_efix_p1MC="835";
# LHC13b2_efix_p2MC="836";
# LHC13b2_efix_p3MC="837";
# LHC13b2_efix_p4MC="838";

# TRAINDIR=Legotrain-vAN20170215-sys-PCM-dirGammaCentAndPurity
# LHC13bData="589"; #pass 3
# LHC13cData="590"; #pass 2
# LHC13b2_efix_p1MC="839";
# LHC13b2_efix_p2MC="840";
# LHC13b2_efix_p3MC="841";
# LHC13b2_efix_p4MC="842";
# LHC13b2_efix_p1MC="843";
# LHC13b2_efix_p2MC="844";
# LHC13b2_efix_p3MC="845";
# LHC13b2_efix_p4MC="846";

# TRAINDIR=Legotrain-vAN20170318-testFixAddSig
# LHC13e7MC="847";

# TRAINDIR=Legotrain-vAN20170329-MCfix
# LHC13b2_efix_p1MC="875";
# LHC13b2_efix_p2MC="876";
# LHC13b2_efix_p3MC="877";
# LHC13b2_efix_p4MC="878";
# LHC13e7MC="879";
# LHC13b2_efix_p1MC="880";
# LHC13b2_efix_p2MC="881";
# LHC13b2_efix_p3MC="882";
# LHC13b2_efix_p4MC="883";
# # LHC13e7MC="884";
# LHC13b2_efix_p1MC="885";
# LHC13b2_efix_p2MC="886";
# LHC13b2_efix_p3MC="887";
# LHC13b2_efix_p4MC="888";
# LHC13e7MC="874";

# TRAINDIR=Legotrain-vAN20170417-Weighting
# LHC13bData="599"; #pass 3
# LHC13cData="601"; #pass 2
# LHC13b2_efix_p1MC="904";
# LHC13b2_efix_p2MC="905";
# LHC13b2_efix_p3MC="906";
# LHC13b2_efix_p4MC="907";
# LHC13e7MC="908";

# TRAINDIR=Legotrain-vAN20170525FF-newDefaultPlusSys
# LHC13bData="631"; #pass 3
# LHC13cData="632"; #pass 2
# LHC13b2_efix_p1MC="961";
# LHC13b2_efix_p2MC="962";
# LHC13b2_efix_p3MC="963";
# LHC13b2_efix_p4MC="964";
# LHC13e7MC="965";
# LHC13e7MC="966";
# LHC13bData="639"; #pass 3
# LHC13cData="636"; #pass 2
# LHC13b2_efix_p1MC="967";
# LHC13b2_efix_p2MC="968";
# LHC13b2_efix_p3MC="969";
# LHC13b2_efix_p4MC="970";
# LHC13e7MC="971";
# LHC13e7MC="972";
# LHC13bData="647"; #pass 3
# LHC13cData="648"; #pass 2

TRAINDIR=Legotrain-vAN20170905-dirGamma
LHC13bcData="672"; #pass 2
LHC13bData="child_1"; #pass 3
LHC13cData="child_2"; #pass 2


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
if [ "$LHC13e7MC" = "" ]; then
    HAVELHC13e7=0;
fi

mkdir -p $OUTPUTDIR/CutSelections


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
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC13b/mergedAllConv.txt.txt ]; then
                rm $OUTPUTDIR_LHC13b/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b.txt`
                ls $OUTPUTDIR_LHC13b/$firstrunNumber/GammaConvV1_*.root > fileLHC13b.txt
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
                    hadd -f $OUTPUTDIR_LHC13b/GammaConvV1_$number.root $OUTPUTDIR_LHC13b/*/GammaConvV1_$number.root
                done;

                if [ $ENABLEFLOW == 1 ]; then
                    ls $OUTPUTDIR_LHC13b/$firstrunNumber/GammaConvFlow_*.root > fileLHC13b.txt
                    fileNumbers=`cat fileLHC13b.txt`
                    for fileName in $fileNumbers; do
                        echo $fileName
                        alpha=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f3 | cut -d "." -f1`
                        if [ -z "$alpha" ]; then
                            number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2`
                            number=$number\_$alpha
                        else
                            number=`echo $fileName  | cut -d "/" -f $NSlashes3 | cut -d "_" -f2 | cut -d "." -f1`
                        fi
                        echo $number
                        hadd -f $OUTPUTDIR_LHC13b/GammaConvFlow_$number.root $OUTPUTDIR_LHC13b/*/GammaConvFlow_$number.root
                    done;
                fi
                echo "done" > $OUTPUTDIR_LHC13b/mergedAllConv.txt.txt
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
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $OUTPUTDIR_LHC13c/mergedAllConv.txt.txt ]; then
                rm $OUTPUTDIR_LHC13c/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13c.txt`
                ls $OUTPUTDIR_LHC13c/$firstrunNumber/GammaConvV1_*.root > fileLHC13c.txt
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
                    hadd -f $OUTPUTDIR_LHC13c/GammaConvV1_$number.root $OUTPUTDIR_LHC13c/*/GammaConvV1_$number.root
                done;
                if [ $ENABLEFLOW == 1 ]; then
                    ls $OUTPUTDIR_LHC13c/$firstrunNumber/GammaConvFlow_*.root > fileLHC13c.txt
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
                        hadd -f $OUTPUTDIR_LHC13c/GammaConvFlow_$number.root $OUTPUTDIR_LHC13c/*/GammaConvFlow_$number.root
                    done;
                fi
                echo "done" > $OUTPUTDIR_LHC13c/mergedAllConv.txt.txt
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
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13b2_efix_p1/mergedAllConv.txt.txt ]; then
                rm $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix1.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p1/$firstrunNumber/GammaConvV1_*.root > fileLHC13b2efixp1.txt
                fileNumbers=`cat fileLHC13b2efixp1.txt`
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvV1_$number.root
                done;
                if [ $ENABLEFLOW == 1 ]; then
                    ls $OUTPUTDIR_LHC13b2_efix_p1/$firstrunNumber/GammaConvFlow_*.root > fileLHC13b2efixp1.txt
                    fileNumbers=`cat fileLHC13b2efixp1.txt`
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
                        hadd -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvFlow_$number.root $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvFlow_$number.root
                    done;
                fi
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p1/mergedAllConv.txt.txt
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
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC13b2_efix_p2/mergedAllConv.txt.txt ]; then
                rm $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix2.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p2/$firstrunNumber/GammaConvV1_*.root > fileLHC13b2efixp2.txt
                fileNumbers=`cat fileLHC13b2efixp2.txt`
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvV1_$number.root
                done;
                if [ $ENABLEFLOW == 1 ]; then
                    ls $OUTPUTDIR_LHC13b2_efix_p2/$firstrunNumber/GammaConvFlow_*.root > fileLHC13b2efixp2.txt
                    fileNumbers=`cat fileLHC13b2efixp2.txt`
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
                        hadd -f $OUTPUTDIR_LHC13b2_efix_p2/GammaConvFlow_$number.root $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvFlow_$number.root
                    done;
                fi
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p2/mergedAllConv.txt.txt
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
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC13b2_efix_p3/mergedAllConv.txt.txt ]; then
                rm $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix3.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p3/$firstrunNumber/GammaConvV1_*.root > fileLHC13b2efixp3.txt
                fileNumbers=`cat fileLHC13b2efixp3.txt`
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvV1_$number.root
                done;
                if [ $ENABLEFLOW == 1 ]; then
                    ls $OUTPUTDIR_LHC13b2_efix_p3/$firstrunNumber/GammaConvFlow_*.root > fileLHC13b2efixp3.txt
                    fileNumbers=`cat fileLHC13b2efixp3.txt`
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
                        hadd -f $OUTPUTDIR_LHC13b2_efix_p3/GammaConvFlow_$number.root $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvFlow_$number.root
                    done;
                fi
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p3/mergedAllConv.txt.txt
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
            if [ $MERGEONSINGLEMC == 1 ]  && [ ! -f $OUTPUTDIR_LHC13b2_efix_p4/mergedAllConv.txt.txt ]; then
                rm $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13b2_efix4.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p4/$firstrunNumber/GammaConvV1_*.root > fileLHC13b2efixp4.txt
                fileNumbers=`cat fileLHC13b2efixp4.txt`
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
                    hadd -f $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvV1_$number.root
                done;
                if [ $ENABLEFLOW == 1 ]; then
                    ls $OUTPUTDIR_LHC13b2_efix_p4/$firstrunNumber/GammaConvFlow_*.root > fileLHC13b2efixp4.txt
                    fileNumbers=`cat fileLHC13b2efixp4.txt`
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
                        hadd -f $OUTPUTDIR_LHC13b2_efix_p4/GammaConvFlow_$number.root $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvFlow_$number.root
                    done;
                fi
                echo "done" > $OUTPUTDIR_LHC13b2_efix_p4/mergedAllConv.txt.txt
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
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $OUTPUTDIR_LHC13e7/mergedAllConv.txt.txt ]; then
                rm $OUTPUTDIR_LHC13e7/GammaConvV1*.root*
                firstrunNumber=`head -n1 runlists/runNumbersLHC13e7.txt`
                ls $OUTPUTDIR_LHC13e7/$firstrunNumber/GammaConvV1_*.root > fileLHC13e7.txt
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
                    hadd -f $OUTPUTDIR_LHC13e7/GammaConvV1_$number.root $OUTPUTDIR_LHC13e7/*/GammaConvV1_$number.root
                done;
                if [ $ENABLEFLOW == 1 ]; then
                    ls $OUTPUTDIR_LHC13e7/$firstrunNumber/GammaConvFlow_*.root > fileLHC13e7.txt
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
                        hadd -f $OUTPUTDIR_LHC13e7/GammaConvFlow_$number.root $OUTPUTDIR_LHC13e7/*/GammaConvFlow_$number.root
                    done;
                fi
                echo "done" > $OUTPUTDIR_LHC13e7/mergedAllConv.txt.txt
            fi
        else
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge" $NSlashes
        fi
    fi

    if [ $HAVELHC13b == 1 ]; then
        ls $OUTPUTDIR_LHC13b/GammaConvV1_*.root > fileLHC13b.txt
        fileNumbers=`cat fileLHC13b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass$passNr"
        done;

        if [ $ENABLEFLOW == 1 ]; then
            ls $OUTPUTDIR_LHC13b/GammaConvFlow_*.root > fileLHC13b.txt
            fileNumbers=`cat fileLHC13b.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass$passNr"
            done;
        fi
    fi

    if [ $HAVELHC13c == 1 ]; then
        ls $OUTPUTDIR_LHC13c/GammaConvV1_*.root > fileLHC13c.txt
        fileNumbers=`cat fileLHC13c.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass$passNr"
        done;

        if [ $ENABLEFLOW == 1 ]; then
            ls $OUTPUTDIR_LHC13c/GammaConvFlow_*.root > fileLHC13c.txt
            fileNumbers=`cat fileLHC13c.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass$passNr"
            done;
        fi
    fi

    if [ $HAVELHC13d == 1 ]; then
        ls $OUTPUTDIR_LHC13d/GammaConvV1_*.root > fileLHC13d.txt
        fileNumbers=`cat fileLHC13d.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13d $NSlashes "LHC13d-pass$passNr"
        done;
    fi

    if [ $HAVELHC13e == 1 ]; then
        ls $OUTPUTDIR_LHC13e/GammaConvV1_*.root > fileLHC13e.txt
        fileNumbers=`cat fileLHC13e.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13e $NSlashes "LHC13e-pass$passNr"
        done;
    fi

    if [ $HAVELHC13f == 1 ]; then
        ls $OUTPUTDIR_LHC13f/GammaConvV1_*.root > fileLHC13f.txt
        fileNumbers=`cat fileLHC13f.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13f $NSlashes "LHC13f-pass$passNr"
        done;
    fi

    if [ $HAVELHC13b2efixp1 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_*.root > fileLHC13b2efixp1.txt
        fileNumbers=`cat fileLHC13b2efixp1.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p1 $NSlashes "MC_LHC13b2_efix_p1"
        done;
        if [ $ENABLEFLOW == 1 ]; then
            ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvFlow_*.root > fileLHC13b2efixp1.txt
            fileNumbers=`cat fileLHC13b2efixp1.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p1 $NSlashes "MC_LHC13b2_efix_p1"
            done;
        fi
    fi

    if [ $HAVELHC13b2efixp2 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_*.root > fileLHC13b2efixp2.txt
        fileNumbers=`cat fileLHC13b2efixp2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p2 $NSlashes "MC_LHC13b2_efix_p2"
        done;
        if [ $ENABLEFLOW == 1 ]; then
            ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvFlow_*.root > fileLHC13b2efixp2.txt
            fileNumbers=`cat fileLHC13b2efixp2.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p2 $NSlashes "MC_LHC13b2_efix_p2"
            done;
        fi
    fi

    if [ $HAVELHC13b2efixp3 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_*.root > fileLHC13b2efixp3.txt
        fileNumbers=`cat fileLHC13b2efixp3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p3 $NSlashes "MC_LHC13b2_efix_p3"
        done;
        if [ $ENABLEFLOW == 1 ]; then
            ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvFlow_*.root > fileLHC13b2efixp3.txt
            fileNumbers=`cat fileLHC13b2efixp3.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p3 $NSlashes "MC_LHC13b2_efix_p3"
            done;
        fi
    fi

    if [ $HAVELHC13b2efixp4 == 1 ]; then
        ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_*.root > fileLHC13b2efixp4.txt
        fileNumbers=`cat fileLHC13b2efixp4.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p4 $NSlashes "MC_LHC13b2_efix_p4"
        done;
        if [ $ENABLEFLOW == 1 ]; then
            ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvFlow_*.root > fileLHC13b2efixp4.txt
            fileNumbers=`cat fileLHC13b2efixp4.txt`
            for fileName in $fileNumbers; do
                ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p4 $NSlashes "MC_LHC13b2_efix_p4"
            done;
        fi
    fi

    if [ $HAVELHC13e7 == 1 ]; then
        ls $OUTPUTDIR_LHC13e7/GammaConvV1_*.root > fileLHC13e7.txt
        fileNumbers=`cat fileLHC13e7.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13e7 $NSlashes "MC_LHC13e7"
        done;
        if [ $ENABLEFLOW == 1 ]; then
            ls $OUTPUTDIR_LHC13e7/GammaConvFlow_*.root > fileLHC13e7.txt
            fileNumbers=`cat fileLHC13e7.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13e7 $NSlashes "MC_LHC13e7"
            done;
        fi
    fi

    if [ $MERGEON == 1 ]; then
        ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 3
            echo $number
            ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root
            ls $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root
            if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_LHC13bc-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root
            fi
        done

        if [ $ENABLEFLOW == 1 ]; then
            rm $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_LHC13c-pass$passNr\_*.root
            ls $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_*.root > filesForMerging.txt
            filesForMerging=`cat filesForMerging.txt`
            for fileName in $filesForMerging; do
                echo $fileName
                GetFileNumberMerging $fileName $((NSlashes-1)) 3
                echo $number
                ls $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_$number.root
                ls $OUTPUTDIR/GammaConvFlow_LHC13c-pass$passNr\_$number.root
                if [ -f $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_LHC13c-pass$passNr\_$number.root ] ; then
                    hadd -f $OUTPUTDIR/GammaConvFlow_LHC13bc-pass$passNr\_$number.root $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvFlow_LHC13c-pass$passNr\_$number.root
                fi
            done
        fi

#         ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_*.root | sed 's/.*p1_p2.*//' | sed '/^$/d' > filesForMergingMC.txt
#         filesForMerging=`cat filesForMergingMC.txt`
#         for fileName in $filesForMerging; do
#             echo $fileName
#             GetFileNumberMerging $fileName $((NSlashes-1)) 6
#             echo $number
#             if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root ]; then
#                     hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root
#             fi
#             if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13e7_$number.root ] ; then
#                 hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_LHC13e7_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13e7_$number.root
#             fi
#         done
#
#
#         if [ $ENABLEFLOW == 1 ]; then
#             ls $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_*.root | sed 's/.*p1_p2.*//' | sed '/^$/d' > filesForMergingMC.txt
#             filesForMerging=`cat filesForMergingMCFlow.txt`
#             for fileName in $filesForMerging; do
#                 echo $fileName
#                 GetFileNumberMerging $fileName $((NSlashes-1)) 6
#                 echo $number
#                 if [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p4_$number.root ]; then
#                     hadd -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p2_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p3_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p4_$number.root
#                 fi
#                 if [ -f $OUTPUTDIR_LHC13e7/GammaConvFlow_$number.root ] ; then
#                     hadd -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_p2_p3_p4_LHC13e7_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13e7_$number.root
#                 fi
#             done
#         fi
    fi
else
    if [ $HAVELHC13b == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13b";
        rm $OUTPUTDIR_LHC13b/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13b/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13c == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13c";
        rm $OUTPUTDIR_LHC13c/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13c/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13d == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13d";
        rm $OUTPUTDIR_LHC13d/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13d/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13e == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13e";
        rm $OUTPUTDIR_LHC13e/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13e/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13f == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13f";
        rm $OUTPUTDIR_LHC13f/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13f/*/GammaConvFlow_*.root
    fi

    if [ $HAVELHC13b2efixp1 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13b2efix1";
        rm $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13b2efixp2 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13b2efix2";
        rm $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13b2efixp3 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13b2efix3";
        rm $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13b2efixp4 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13b2efix4";
        rm $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvFlow_*.root
    fi
    if [ $HAVELHC13e7 == 1 ]; then
        echo "removing all GammaConv files in runFolders for LHC13e7";
        rm $OUTPUTDIR_LHC13e7/*/GammaConvV1_*.root
        rm $OUTPUTDIR_LHC13e7/*/GammaConvFlow_*.root
    fi
fi
