#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

#! /bin/bash

# This script has to be run with "bash"

function CopyFileIfNonExisitent()
{
    if [ -f $1/root_archive.zip ] && [ -s $1/root_archive.zip ]; then 
        echo "$1/root_archive.zip exists";
    else     
        mkdir -p $1
        alien_cp alien:$2/root_archive.zip file:$1/
    fi    
    unzip -u $1/root_archive.zip -d $1/
}

function ChangeStructureIfNeeded()
{
    cp $1 $2
#     if [ -f $2 ]; then 
#         echo "already changed"
#     else
#         root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$1\"\,\"$2\"\,\"GammaConvV1_$3\"\)
#     fi    
}

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



# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1
SINGLERUN=1
SEPARATEON=1
MERGEONSINGLEData=0
MERGEONSINGLEMC=1

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

TRAINDIR=Legotrain-vAN20170205-sys-PCM-dirGammaWoWeight
# LHC13b2_efix_p1MC="831"; 
# LHC13b2_efix_p2MC="832"; 
# LHC13b2_efix_p3MC="833"; 
# LHC13b2_efix_p4MC="834";
LHC13b2_efix_p1MC="835"; 
LHC13b2_efix_p2MC="836"; 
LHC13b2_efix_p3MC="837"; 
LHC13b2_efix_p4MC="838";

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
    
if [ $DOWNLOADON == 1 ]; then
    if [ $HAVELHC13b == 1 ]; then
        echo "downloading LHC13b"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13bData/merge_runlist_1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b/$runNumber "/alice/data/2013/LHC13b/000$runNumber/ESDs/pass3/PWGGA/GA_pPb/$LHC13bData"
                if [ $SEPARATEON == 1 ]; then
                    GetFileNumberList $OUTPUTDIR_LHC13b/$runNumber $NSlashes3 fileNumbers2.txt                    
                    fileNumbers=`cat fileNumbers2.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$OUTPUTDIR_LHC13b/$runNumber/GammaConvV1_$fileNumber.root\"\,\"$OUTPUTDIR_LHC13b/$runNumber/GammaConvV1_$fileNumber\"\,0\)
                    done;
                fi    
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
            fi    
        fi
    fi    
    if [ $HAVELHC13c == 1 ]; then
        echo "downloading LHC13c"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13c "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge_runlist_1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13c.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13c/$runNumber "/alice/data/2013/LHC13c/000$runNumber/ESDs/pass2/PWGGA/GA_pPb/$LHC13cData"
                if [ $SEPARATEON == 1 ]; then
                    GetFileNumberList $OUTPUTDIR_LHC13c/$runNumber $NSlashes3 fileNumbers2.txt                    
                    fileNumbers=`cat fileNumbers2.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$OUTPUTDIR_LHC13c/$runNumber/GammaConvV1_$fileNumber.root\"\,\"$OUTPUTDIR_LHC13c/$runNumber/GammaConvV1_$fileNumber\"\,0\)
                    done;
                fi    
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
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC/merge"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p1/$runNumber "/alice/sim/2013/LHC13b2_efix_p1/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p1MC"
                if [ $SEPARATEON == 1 ]; then
                    GetFileNumberList $OUTPUTDIR_LHC13b2_efix_p1/$runNumber $NSlashes3 fileNumbers2.txt                    
                    fileNumbers=`cat fileNumbers2.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        if [ -f $OUTPUTDIR_LHC13b2_efix_p1/$runNumber/GammaConvV1_$fileNumber\_A.root ]; then 
                            echo "separated $OUTPUTDIR_LHC13b2_efix_p1/$runNumber/GammaConvV1_$fileNumber.root already"
                        else 
                            root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$OUTPUTDIR_LHC13b2_efix_p1/$runNumber/GammaConvV1_$fileNumber.root\"\,\"$OUTPUTDIR_LHC13b2_efix_p1/$runNumber/GammaConvV1_$fileNumber\"\,0\)
                        fi
                    done;
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
            fi    
        fi
    fi    
    if [ $HAVELHC13b2efixp2 == 1 ]; then
        echo "downloading LHC13b2_efix_p2"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/merge"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix2.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2/$runNumber "/alice/sim/2013/LHC13b2_efix_p2/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC"
                if [ $SEPARATEON == 1 ]; then
                    GetFileNumberList $OUTPUTDIR_LHC13b2_efix_p2/$runNumber $NSlashes3 fileNumbers2.txt                    
                    fileNumbers=`cat fileNumbers2.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        if [ -f $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/GammaConvV1_$fileNumber\_A.root ]; then 
                            echo "separated $OUTPUTDIR_LHC13b2_efix_p2/$runNumber/GammaConvV1_$fileNumber.root already"
                        else 
                            root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$OUTPUTDIR_LHC13b2_efix_p2/$runNumber/GammaConvV1_$fileNumber.root\"\,\"$OUTPUTDIR_LHC13b2_efix_p2/$runNumber/GammaConvV1_$fileNumber\"\,0\)
                        fi   
                    done;
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
            fi    
        fi
    fi    
    if [ $HAVELHC13b2efixp3 == 1 ]; then
        echo "downloading LHC13b2_efix_p3"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/merge"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix3.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3/$runNumber "/alice/sim/2013/LHC13b2_efix_p3/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC"
                if [ $SEPARATEON == 1 ]; then
                    GetFileNumberList $OUTPUTDIR_LHC13b2_efix_p3/$runNumber $NSlashes3 fileNumbers2.txt                    
                    fileNumbers=`cat fileNumbers2.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        if [ -f $OUTPUTDIR_LHC13b2_efix_p3/$runNumber/GammaConvV1_$fileNumber\_A.root ]; then 
                            echo "separated $OUTPUTDIR_LHC13b2_efix_p3/$runNumber/GammaConvV1_$fileNumber.root already"
                        else 
                            root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$OUTPUTDIR_LHC13b2_efix_p3/$runNumber/GammaConvV1_$fileNumber.root\"\,\"$OUTPUTDIR_LHC13b2_efix_p3/$runNumber/GammaConvV1_$fileNumber\"\,0\)
                        fi
                    done;
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
            fi   
        fi
    fi    
    if [ $HAVELHC13b2efixp4 == 1 ]; then
        echo "downloading LHC13b2_efix_p4"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/merge"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13b2_efix4.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4/$runNumber "/alice/sim/2013/LHC13b2_efix_p4/$runNumber/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC"
                if [ $SEPARATEON == 1 ]; then
                    GetFileNumberList $OUTPUTDIR_LHC13b2_efix_p4/$runNumber $NSlashes3 fileNumbers2.txt                    
                    fileNumbers=`cat fileNumbers2.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        if [ -f $OUTPUTDIR_LHC13b2_efix_p4/$runNumber/GammaConvV1_$fileNumber\_A.root ]; then 
                            echo "separated $OUTPUTDIR_LHC13b2_efix_p4/$runNumber/GammaConvV1_$fileNumber.root already"
                        else 
                            root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$OUTPUTDIR_LHC13b2_efix_p4/$runNumber/GammaConvV1_$fileNumber.root\"\,\"$OUTPUTDIR_LHC13b2_efix_p4/$runNumber/GammaConvV1_$fileNumber\"\,0\)
                        fi
                    done;
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
            fi
        fi
    fi    
    if [ $HAVELHC13e7 == 1 ]; then
        echo "downloading LHC13e7"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbersLHC13e7.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7/$runNumber "/alice/sim/2013/LHC13e7/$runNumber/PWGGA/GA_pPb_MC/$LHC13e7MC"
                if [ $SEPARATEON == 1 ]; then
                    GetFileNumberList $OUTPUTDIR_LHC13e7/$runNumber $NSlashes3 fileNumbers2.txt                    
                    fileNumbers=`cat fileNumbers2.txt`
                    for fileNumber in $fileNumbers; do
                        echo $fileNumber
                        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$OUTPUTDIR_LHC13e7/$runNumber/GammaConvV1_$fileNumber.root\"\,\"$OUTPUTDIR_LHC13e7/$runNumber/GammaConvV1_$fileNumber\"\,0\)
                    done;
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
            fi
        fi
    fi        
fi
    
if [ $HAVELHC13b == 1 ]; then
    ls $OUTPUTDIR_LHC13b/GammaConvV1_*.root > fileLHC13b.txt
    fileNumbers=`cat file_mergeSpecial.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=$fileName
#         number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root $number
#         ChangeStructureIfNeeded $OUTPUTDIR_LHC13b/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_LHC13b_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13c == 1 ]; then
    ls $OUTPUTDIR_LHC13c/GammaConvV1_*.root > fileLHC13c.txt
    fileNumbers=`cat file_mergeSpecial.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=$fileName
#       number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13c/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_LHC13c_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13d == 1 ]; then
    ls $OUTPUTDIR_LHC13d/GammaConvV1_*.root > fileLHC13d.txt
    fileNumbers=`cat fileLHC13d.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13d/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_LHC13d_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13e == 1 ]; then
    ls $OUTPUTDIR_LHC13e/GammaConvV1_*.root > fileLHC13e.txt
    fileNumbers=`cat fileLHC13e.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13e/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_LHC13e-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13e-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_LHC13e_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13f == 1 ]; then
    ls $OUTPUTDIR_LHC13f/GammaConvV1_*.root > fileLHC13f.txt
    fileNumbers=`cat fileLHC13f.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13f/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_LHC13f-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_LHC13f-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_LHC13f_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13b2efixp1 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_*.root > fileLHC13b2efixp1.txt
    fileNumbers=`cat fileLHC13b2efixp1.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_MC_LHC13b2_efix_p1_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13b2efixp2 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_*.root > fileLHC13b2efixp2.txt
    fileNumbers=`cat fileLHC13b2efixp2.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_MC_LHC13b2_efix_p2_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13b2efixp3 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_*.root > fileLHC13b2efixp3.txt
    fileNumbers=`cat fileLHC13b2efixp3.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_MC_LHC13b2_efix_p3_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13b2efixp4 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_*.root > fileLHC13b2efixp4.txt
    fileNumbers=`cat fileLHC13b2efixp4.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_MC_LHC13b2_efix_p4_$number.log\"\,0\)
    done;
fi

if [ $HAVELHC13e7 == 1 ]; then
    ls $OUTPUTDIR_LHC13e7/GammaConvV1_*.root > fileLHC13e7.txt
    fileNumbers=`cat fileLHC13e7.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13e7/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13e7_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_MC_LHC13e7_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvV1_MC_LHC13e7_$number.log\"\,0\)
    done;
fi

mkdir -p $OUTPUTDIR/CutSelections
mv $OUTPUTDIR/CutSelection_*.log $OUTPUTDIR/CutSelections/

if [ $MERGEON == 1 ]; then
#     rm $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_*.root
#     ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_*.root > filesForMerging.txt
#     filesForMerging=`cat filesForMerging.txt`
#     for fileName in $filesForMerging; do
#         echo $fileName
#         number=$fileName
# #         number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
#         echo $number
#         ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root
#         ls $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root
#         if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root ] ; then
#             hadd -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13c-pass$passNr\_$number.root
#         fi
#     done
#     
#     ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_*.root > filesForMerging.txt
#     filesForMerging=`cat filesForMerging.txt`
#     for fileName in $filesForMerging; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 4 | cut -d "." -f1`
#         echo $number
#         if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_$number.root ] ; then
#             hadd -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_$number.root
#         fi
#     done
#     
#     ls $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_*.root > filesForMerging.txt
#     filesForMerging=`cat filesForMerging.txt`
#     for fileName in $filesForMerging; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 5 | cut -d "." -f1`
#         echo $number
#         if [ -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13e-pass$passNr\_$number.root ] ; then
#             hadd -f $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_LHC13e-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13e-pass$passNr\_$number.root
#         fi
#     done
# 
#     rm $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_LHC13e-pass$passNr\_*.root
#     ls $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_*.root > filesForMerging.txt
#     filesForMerging=`cat filesForMerging.txt`
#     for fileName in $filesForMerging; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
#         echo $number
#         if [ -f $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_LHC13e-pass$passNr\_$number.root ] ; then
#             hadd -f $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_LHC13e-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvV1_LHC13e-pass$passNr\_$number.root
#         fi
#     done

#     rm $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_*.root
#     ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_*.root > filesForMerging.txt
#     filesForMerging=`cat filesForMerging.txt`
#     for fileName in $filesForMerging; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 6 | cut -d "." -f1`
#         echo $number
#         if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root ] ; then
#             hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p2_$number.root
#         fi
#     done
#     #    

    
#     ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMergingMC.txt`
#     filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
          number=$fileName
#         number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 7 | cut -d "." -f1`
        echo $number
#         if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root ] ; then
#             hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p3_$number.root
#         fi
        if [ -f $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number\_mergedByHand.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR_LHC13b2_efix_p1/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p2/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p3/GammaConvV1_$number\_mergedByHand.root $OUTPUTDIR_LHC13b2_efix_p4/GammaConvV1_$number\_mergedByHand.root
        fi    
    done

    for fileName in $filesForMerging; do
        echo $fileName
        number=$fileName
        echo $number
        if [ -f $OUTPUTDIR_LHC13e7/GammaConvV1_$number\_mergedByHand.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_LHC13e7_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR_LHC13e7/GammaConvV1_$number\_mergedByHand.root 
        fi
    done
    
    #     
#     ls $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_*.root > filesForMerging.txt
#     filesForMerging=`cat filesForMerging.txt`
#     for fileName in $filesForMerging; do
#         echo $fileName
#         number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 8 | cut -d "." -f1`
#         echo $number
#         if [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root ] ; then
#             hadd -f $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaConvV1_MC_LHC13b2_efix_p4_$number.root
#         fi
#     done;
fi
