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
SINGLEJOB=1
SEPARATEON=1
MERGEONSINGLEData=0
MERGEONSINGLEMC=0
SPECIALMERGE=0

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
    else 
        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$1.root\"\,\"$1\"\,0\)
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
        GetFileNumberListFlow $1 $3 fileNumbers2.txt                    
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaConvFlow_$fileNumber
        done;
    fi
}



function ChangeStructureIfNeeded()
{
    if [ $SPECIALMERGE == 1 ]; then
        number=$1
        echo $number
        cp $2/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR/GammaConvV1_$4\_$number.root 
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
    fi
    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_$4_$number.log\"\,0\)
}

function ChangeStructureIfNeededFlow()
{
    if [ $SPECIALMERGE == 1 ]; then
        number=$1
        echo $number
        cp $2/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR/GammaConvFlow_$4\_$number.root
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
        cp $2/GammaConvFlow_$number.root $OUTPUTDIR/GammaConvFlow_$4\_$number.root 
    fi
    root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvFlow_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvFlow_$4_$number.log\"\,0\)
}


function GetFileNumberMerging()
{
    if [ $SPECIALMERGE == 1 ]; then
        number=$1
    else
        NCurrSub=$3
        number1=`echo $1  | cut -d "/" -f $2 | cut -d "_" -f $NCurrSub | cut -d "." -f1`
        number2=`echo $1  | cut -d "/" -f $2 | cut -d "_" -f $((NCurrSub+1)) | cut -d "." -f1`
        if [ -z "$number2" ]; then
            number=$number1
        else 
            echo $number2
            number=$number1\_$number2
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

TRAINDIR=Legotrain-vAN20170329-MCfix
# LHC13b2_efix_p1MC="857"; 
# LHC13b2_efix_p2MC="858"; 
# LHC13b2_efix_p3MC="859"; 
LHC13b2_efix_p4MC="860";
# LHC13e7MC="870";
# LHC13b2_efix_p1MC="861"; 
# LHC13b2_efix_p2MC="862"; 
# LHC13b2_efix_p3MC="863"; 
# LHC13b2_efix_p4MC="864";
# LHC13e7MC="871";


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
    
if [ $HAVELHC13b == 1 ]; then
    echo "downloading LHC13b"
    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13bData/merge_runlist_1" $NSlashes
    if [ $SINGLERUN == 1 ]; then
        runNumbers=`cat runlists/runNumbersLHC13b.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b/$runNumber "/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/GA_pPb/$LHC13bData" $NSlashes3
        done;    
        if [ $MERGEONSINGLEData == 1 ]; then
            firstrunNumber=`head -n1 runlists/runNumbersLHC13b.txt`
            ls $OUTPUTDIR_LHC13b/$firstrunNumber/GammaConvV1_*.root > fileLHC13b.txt
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
                hadd -f $OUTPUTDIR_LHC13b/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b/*/GammaConvV1_$number.root
            done;
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
                hadd -f $OUTPUTDIR_LHC13b/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b/*/GammaConvFlow_$number.root
            done;
        fi    
    fi
fi    
if [ $HAVELHC13c == 1 ]; then
    echo "downloading LHC13c"
    CopyFileIfNonExisitent $OUTPUTDIR_LHC13c "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge_runlist_1" $NSlashes
    if [ $SINGLERUN == 1 ]; then
        runNumbers=`cat runlists/runNumbersLHC13c.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13c/$runNumber "/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/GA_pPb/$LHC13cData" $NSlashes3
        done;    
        if [ $MERGEONSINGLEData == 1 ]; then
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
                hadd -f $OUTPUTDIR_LHC13c/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13c/*/GammaConvV1_$number.root
            done;
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
                hadd -f $OUTPUTDIR_LHC13c/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13c/*/GammaConvFlow_$number.root
            done;
        fi    
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
    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/merge" $NSlashes
    if [ $SINGLERUN == 1 ]; then
        runNumbers=`cat runlists/runNumbersLHC13b2_efix1.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            if [ $SINGLEJOB == 1 ]; then 
                if [ $SEPARATEON == 1 ]; then 
                    SEPARATEONTemp=1
                    SEPARATEON=0
                fi
                alien_ls /alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/ > singleJobOutput.txt 
                cat singleJobOutput.txt
                filelocations=`cat singleJobOutput.txt`
                for t in $filelocations; do
                    dummy=`echo $t | cut -d "." -f 2`
                    filenumber=`echo $t | cut -d "." -f 1`
                    if [ $dummy == "$filenumber" ]; then 
                    echo "dir found" $filenumber;
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1/$runNumber/$filenumber "alien:/alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/$filenumber" $NSlashes4
                    else 
                    echo "that a file not a dir" $dummy;
                    fi
                done
                firstFileNumber=`head -n1 singleJobOutput.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p1/$runNumber/$firstFileNumber/GammaConvV1_*.root > fileLHC13b2efixp1.txt
                fileNumbers=`cat fileLHC13b2efixp1.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        echo $number
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2`
                        number=$number\_$alpha
                    fi
                    echo $number
                    hadd -n 10 -f $OUTPUTDIR_LHC13b2_efix_p1/$runNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p1/$runNumber/*/GammaConvV1_$number.root
                done;
                if [ $SEPARATEONTemp == 1 ]; then 
                    SEPARATEON=1
                fi
                rm -rf $OUTPUTDIR_LHC13b2_efix_p1/$runNumber/*/GammaConvV1_*.root
            else 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1/$runNumber "/alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC" $NSlashes3
            fi  
        done;
        if [ $MERGEONSINGLEMC == 1 ]; then
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvV1_$number.root
            done;
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p1/*/GammaConvFlow_$number.root
            done;
        fi    
    fi
fi    
if [ $HAVELHC13b2efixp2 == 1 ]; then
    echo "downloading LHC13b2_efix_p2"
    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/merge" $NSlashes
    if [ $SINGLERUN == 1 ]; then
        runNumbers=`cat runlists/runNumbersLHC13b2_efix2.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            if [ $SINGLEJOB == 1 ]; then 
                if [ $SEPARATEON == 1 ]; then 
                    SEPARATEONTemp=1
                    SEPARATEON=0
                fi
                alien_ls /alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/ > singleJobOutput.txt 
                cat singleJobOutput.txt
                filelocations=`cat singleJobOutput.txt`
                for t in $filelocations; do
                    dummy=`echo $t | cut -d "." -f 2`
                    filenumber=`echo $t | cut -d "." -f 1`
                    if [ $dummy == "$filenumber" ]; then 
                    echo "dir found" $filenumber;
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/$filenumber "alien:/alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/$filenumber" $NSlashes4
                    else 
                    echo "that a file not a dir" $dummy;
                    fi
                done
                firstFileNumber=`head -n1 singleJobOutput.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/$firstFileNumber/GammaConvV1_*.root > fileLHC13b2efixp2.txt
                fileNumbers=`cat fileLHC13b2efixp2.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        echo $number
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2`
                        number=$number\_$alpha
                    fi
                    echo $number
                    hadd -n 10 -f $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/*/GammaConvV1_$number.root
                    rm -rf $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/*/GammaConvV1_$number.root
                done;
                if [ $SEPARATEONTemp == 1 ]; then 
                    SEPARATEON=1
                fi
                rm -rf $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/*/GammaConvV1_*.root
            else 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2/$runNumber "/alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC" $NSlashes3
            fi
        done;
        if [ $MERGEONSINGLEMC == 1 ]; then
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvV1_$number.root
            done;
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p2/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p2/*/GammaConvFlow_$number.root
            done;
        fi    
    fi
fi    
if [ $HAVELHC13b2efixp3 == 1 ]; then
    echo "downloading LHC13b2_efix_p3"
    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/merge" $NSlashes
    if [ $SINGLERUN == 1 ]; then
        runNumbers=`cat runlists/runNumbersLHC13b2_efix3.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            if [ $SINGLEJOB == 1 ]; then 
                if [ $SEPARATEON == 1 ]; then 
                    SEPARATEONTemp=1
                    SEPARATEON=0
                fi
                alien_ls /alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/ > singleJobOutput.txt 
                cat singleJobOutput.txt
                filelocations=`cat singleJobOutput.txt`
                for t in $filelocations; do
                    dummy=`echo $t | cut -d "." -f 2`
                    filenumber=`echo $t | cut -d "." -f 1`
                    if [ $dummy == "$filenumber" ]; then 
                    echo "dir found" $filenumber;
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3/$runNumber/$filenumber "alien:/alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/$filenumber" $NSlashes4
                    else 
                    echo "that a file not a dir" $dummy;
                    fi
                done
                firstFileNumber=`head -n1 singleJobOutput.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p3/$runNumber/$firstFileNumber/GammaConvV1_*.root > fileLHC13b2efixp3.txt
                fileNumbers=`cat fileLHC13b2efixp3.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        echo $number
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2`
                        number=$number\_$alpha
                    fi
                    echo $number
                    hadd -n 10 -f $OUTPUTDIR_LHC13b2_efix_p3/$runNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p3/$runNumber/*/GammaConvV1_$number.root
                    rm -rf $OUTPUTDIR_LHC13b2_efix_p3/$runNumber/*/GammaConvV1_*.root
                done;
                if [ $SEPARATEONTemp == 1 ]; then 
                    SEPARATEON=1
                fi
            else 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3/$runNumber "/alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC" $NSlashes3
            fi
        done;
        if [ $MERGEONSINGLEMC == 1 ]; then
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvV1_$number.root
            done;
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p3/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p3/*/GammaConvFlow_$number.root
            done;
        fi   
    fi
fi    
if [ $HAVELHC13b2efixp4 == 1 ]; then
    echo "downloading LHC13b2_efix_p4"
    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/merge" $NSlashes
    if [ $SINGLERUN == 1 ]; then
        runNumbers=`cat runlists/runNumbersLHC13b2_efix4.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            if [ $SINGLEJOB == 1 ]; then 
                if [ $SEPARATEON == 1 ]; then 
                    SEPARATEONTemp=1
                    SEPARATEON=0
                fi
                alien_ls /alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/ > singleJobOutput.txt 
                cat singleJobOutput.txt
                filelocations=`cat singleJobOutput.txt`
                for t in $filelocations; do
                    dummy=`echo $t | cut -d "." -f 2`
                    filenumber=`echo $t | cut -d "." -f 1`
                    if [ $dummy == "$filenumber" ]; then 
                    echo "dir found" $filenumber;
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4/$runNumber/$filenumber "alien:/alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/$filenumber" $NSlashes4
                    else 
                    echo "that a file not a dir" $dummy;
                    fi
                done
                firstFileNumber=`head -n1 singleJobOutput.txt`
                ls $OUTPUTDIR_LHC13b2_efix_p4/$runNumber/$firstFileNumber/GammaConvV1_*.root > fileLHC13b2efixp4.txt
                fileNumbers=`cat fileLHC13b2efixp4.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        echo $number
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2`
                        number=$number\_$alpha
                    fi
                    echo $number
                    hadd -n 10 -f $OUTPUTDIR_LHC13b2_efix_p4/$runNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC13b2_efix_p4/$runNumber/*/GammaConvV1_$number.root
                    rm -rf $OUTPUTDIR_LHC13b2_efix_p4/$runNumber/*/GammaConvV1_*.root
                done;
                if [ $SEPARATEONTemp == 1 ]; then 
                    SEPARATEON=1
                fi
            else 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4/$runNumber "/alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC" $NSlashes3
            fi
        done;
        if [ $MERGEONSINGLEMC == 1 ]; then
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvV1_$number.root
            done;
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
                hadd -f $OUTPUTDIR_LHC13b2_efix_p4/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p4/*/GammaConvFlow_$number.root
            done;
        fi
    fi
fi    
if [ $HAVELHC13e7 == 1 ]; then
    echo "downloading LHC13e7"
    CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge" $NSlashes
    if [ $SINGLERUN == 1 ]; then
        runNumbers=`cat runlists/runNumbersLHC13e7.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            if [ $SINGLEJOB == 1 ]; then 
                if [ $SEPARATEON == 1 ]; then 
                    SEPARATEONTemp=1
                    SEPARATEON=0
                fi
                alien_ls /alice/sim/2013/LHC13e7/$runNumber/PWGGA/GA_pPb_MC/$LHC13e7MC/ > singleJobOutput.txt 
                cat singleJobOutput.txt
                filelocations=`cat singleJobOutput.txt`
                for t in $filelocations; do
                    dummy=`echo $t | cut -d "." -f 2`
                    filenumber=`echo $t | cut -d "." -f 1`
                    if [ $dummy == "$filenumber" ]; then 
                    echo "dir found" $filenumber;
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7/$runNumber/$filenumber "alien:/alice/sim/2013/LHC13e7/$runNumber/PWGGA/GA_pPb_MC/$LHC13e7MC/$filenumber" $NSlashes4
                    else 
                    echo "that a file not a dir" $dummy;
                    fi
                done
                firstFileNumber=`head -n1 singleJobOutput.txt`
                ls $OUTPUTDIR_LHC13e7/$runNumber/$firstFileNumber/GammaConvV1_*.root > fileLHC13b2efixp1.txt
                fileNumbers=`cat fileLHC13b2efixp1.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    alpha=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f3 | cut -d "." -f1`
                    if [ -z "$alpha" ]; then
                        echo $alpha
                        echo $number
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2 | cut -d "." -f1`
                        echo $number
                    else
                        echo $alpha
                        number=`echo $fileName  | cut -d "/" -f $NSlashes4 | cut -d "_" -f2`
                        number=$number\_$alpha
                    fi
                    echo $number
                    hadd -n 10 -f $OUTPUTDIR_LHC13e7/$runNumber/GammaConvV1_$number.root $OUTPUTDIR_LHC13e7/$runNumber/*/GammaConvV1_$number.root
                done;
                if [ $SEPARATEONTemp == 1 ]; then 
                    SEPARATEON=1
                fi
                rm -rf $OUTPUTDIR_LHC13e7/$runNumber/*/GammaConvV1_*.root
            else 
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7/$runNumber "/alice/sim/2013/LHC13e7/$runNumber/PWGGA/GA_pPb_MC/$LHC13e7MC" $NSlashes3
            fi
        done;
        if [ $MERGEONSINGLEMC == 1 ]; then
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
                hadd -f $OUTPUTDIR_LHC13e7/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13e7/*/GammaConvV1_$number.root
            done;
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
                hadd -f $OUTPUTDIR_LHC13e7/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13e7/*/GammaConvFlow_$number.root
            done;
        fi
    fi
fi        
    
if [ $HAVELHC13b == 1 ]; then
    ls $OUTPUTDIR_LHC13b/GammaConvV1_*.root > fileLHC13b.txt
    if [ $SPECIALMERGE == 1 ]; then
        fileNumbers=`cat file_mergeSpecial.txt`
    else 
        fileNumbers=`cat fileLHC13b.txt`
    fi
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass$passNr"
    done;
    
    ls $OUTPUTDIR_LHC13b/GammaConvFlow_*.root > fileLHC13b.txt
    if [ $SPECIALMERGE == 1 ]; then
        fileNumbers=`cat file_mergeSpecial2.txt`
    else 
        fileNumbers=`cat fileLHC13b.txt`
    fi
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b $NSlashes "LHC13b-pass$passNr"
    done;    
fi

if [ $HAVELHC13c == 1 ]; then
    ls $OUTPUTDIR_LHC13c/GammaConvV1_*.root > fileLHC13c.txt
    if [ $SPECIALMERGE == 1 ]; then
        fileNumbers=`cat file_mergeSpecial.txt`
    else 
        fileNumbers=`cat fileLHC13c.txt`
    fi
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass$passNr"
    done;
    
    ls $OUTPUTDIR_LHC13c/GammaConvFlow_*.root > fileLHC13c.txt
    if [ $SPECIALMERGE == 1 ]; then
        fileNumbers=`cat file_mergeSpecial2.txt`
    else 
        fileNumbers=`cat fileLHC13c.txt`
    fi
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13c $NSlashes "LHC13c-pass$passNr"
    done;
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
    ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvFlow_*.root > fileLHC13b2efixp1.txt
    fileNumbers=`cat fileLHC13b2efixp1.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p1 $NSlashes "MC_LHC13b2_efix_p1"
    done;
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
    ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvFlow_*.root > fileLHC13b2efixp2.txt
    fileNumbers=`cat fileLHC13b2efixp2.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p2 $NSlashes "MC_LHC13b2_efix_p2"
    done;
fi

if [ $HAVELHC13b2efixp3 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_*.root > fileLHC13b2efixp3.txt
    fileNumbers=`cat fileLHC13b2efixp3.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p3 $NSlashes "MC_LHC13b2_efix_p3"
    done;
    ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvFlow_*.root > fileLHC13b2efixp3.txt
    fileNumbers=`cat fileLHC13b2efixp3.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p3 $NSlashes "MC_LHC13b2_efix_p3"    
    done;
fi

if [ $HAVELHC13b2efixp4 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_*.root > fileLHC13b2efixp4.txt
    fileNumbers=`cat fileLHC13b2efixp4.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13b2_efix_p4 $NSlashes "MC_LHC13b2_efix_p4"
    done;
    ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvFlow_*.root > fileLHC13b2efixp4.txt
    fileNumbers=`cat fileLHC13b2efixp4.txt`
    for fileName in $fileNumbers; do
        ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13b2_efix_p4 $NSlashes "MC_LHC13b2_efix_p4"
    done;
fi

if [ $HAVELHC13e7 == 1 ]; then
    ls $OUTPUTDIR_LHC13e7/GammaConvV1_*.root > fileLHC13e7.txt
    fileNumbers=`cat fileLHC13e7.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeeded $fileName $OUTPUTDIR_LHC13e7 $NSlashes "MC_LHC13e7"
    done;
    ls $OUTPUTDIR_LHC13e7/GammaConvFlow_*.root > fileLHC13e7.txt
    fileNumbers=`cat fileLHC13e7.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        ChangeStructureIfNeededFlow $fileName $OUTPUTDIR_LHC13e7 $NSlashes "MC_LHC13e7"
    done;
fi

mkdir -p $OUTPUTDIR/CutSelections
mv $OUTPUTDIR/CutSelection_*.log $OUTPUTDIR/CutSelections/

if [ $MERGEON == 1 ]; then
    if [ $SPECIALMERGE == 0 ]; then 
        rm $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_*.root
        ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_*.root > filesForMerging.txt
    fi
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        GetFileNumberMerging $fileName $((NSlashes-1)) 3
        echo $number
        ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root
        ls $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root
        if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root
        fi
    done

    if [ $SPECIALMERGE == 0 ]; then 
        rm $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_LHC13c-pass$passNr\_*.root
        ls $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_*.root > filesForMerging.txt
    fi
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        GetFileNumberMerging $fileName $((NSlashes-1)) 3
        echo $number
        ls $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_$number.root
        ls $OUTPUTDIR/GammaConvFlow_LHC13c-pass$passNr\_$number.root
        if [ -f $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_LHC13c-pass$passNr\_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvFlow_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvFlow_LHC13c-pass$passNr\_$number.root
        fi
    done
    
    if [ $SPECIALMERGE == 1 ]; then 
        filesForMerging=`cat filesForMergingMC.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 7
            echo $number
            if [ -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number\_mergedByHand.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_$number\_mergedByHand.root
            fi    
        done

        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 7
            echo $number
            if [ -f $OUTPUTDIR_LHC13e7/GammaConvV1_$number\_mergedByHand.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_LHC13e7_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR_LHC13e7/GammaConvV1_$number\_mergedByHand.root 
            fi
        done
    else 
        ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 6
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root
            fi
        done    
    fi

    if [ $SPECIALMERGE == 1 ]; then 
        filesForMerging=`cat filesForMergingMCFlow.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 7
            echo $number
            if [ -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvFlow_$number\_mergedByHand.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR_LHC13b2_efix_p1/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p2/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p3/GammaConvFlow_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p4/GammaConvFlow_$number\_mergedByHand.root
            fi    
        done

        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 7
            echo $number
            if [ -f $OUTPUTDIR_LHC13e7/GammaConvFlow_$number\_mergedByHand.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_p2_p3_p4_LHC13e7_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR_LHC13e7/GammaConvFlow_$number\_mergedByHand.root 
            fi
        done
    else 
        ls $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            GetFileNumberMerging $fileName $((NSlashes-1)) 6
            echo $number
            if [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p4_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p2_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p3_$number.root $OUTPUTDIR/GammaConvFlow_MC_LHC13b2_efix_p4_$number.root
            fi
        done    
    fi
fi
