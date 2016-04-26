#! /bin/bash

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
    if [ -f $2 ]; then 
        echo "already changed"
    else
        root -l -b -q -x ChangeStructureToStandard.C\(\"$1\"\,\"$2\"\,\"GammaConvV1_$3\"\)
    fi    
}



# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs

# switches to enable/disable certain procedures
DOWNLOADON=0
MERGEON=0

# check if train configuration has actually been given
HAVELHC10h=1
HAVELHC13d2=1
HAVELHC13d2b=1

LHC10hData=""
LHC13d2MC="";
LHC13d2bMC=""; 


NSlashes=10
NSlashes2=9

if [ $1 = "fbock" ]; then 
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/PbPb
    NSlashes=8
    NSlashes2=7
elif [ $1 = "fbockGSI" ]; then 
    BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/PbPb
elif [ $1 = "leardini" ]; then 
    BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then 
    BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
elif [ $1 = "leardiniGSI" ]; then 
    BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/PbPb
elif [ $1 = "passfeld" ]; then 
    BASEDIR=~/work/Gridoutput/PbPb
elif [ $1 = "passfeldMAF" ]; then 
    BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then  
    BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/PbPb
elif [ $1 = "amarin" ]; then     
    BASEDIR=/Users/marin/cd 
elif [ $1 = "amarinGSI" ]; then     
    BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/PbPb 
elif [ $1 = "amarinALICESERV1" ]; then     
    BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/   
elif [ $1 = "mwilde" ]; then        
    BASEDIR=~/work/GridOutput 
elif [ $1 = "mwildeGSI" ]; then           
    BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/PbPb 
elif [ $1 = "pgonzales" ]; then     
    BASEDIR=~/work/GridOutput 
elif [ $1 = "pgonzalesGSI" ]; then        
    BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/PbPb
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
    NSlashes=9
    NSlashes2=8
elif [ $1 = "bsahlmul" ]; then
    BASEDIR=/Users/sahlmul/Documents/ALICE/mycode/AnalysisSoftware/grid
    NSlashes=11
    NSlashes2=10
fi
mkdir -p $BASEDIR
  
TRAINDIR=DirectPhotonTest_Kappa
LHC10hData=203; #ESD
LHC13d2MC=250;

  
OUTPUTDIR=$BASEDIR/$TRAINDIR
if [ $2 = "AOD" ]; then
   TRAINPATHData=GA_PbPb_AOD
else
   TRAINPATHData=GA_PbPb
fi   

if [ $2 = "AOD" ]; then
   TRAINPATHMC=GA_PbPb_MC_AOD
else
   TRAINPATHMC=GA_PbPb_MC
fi   

if [ "$LHC10hData" == "" ]; then 
    HAVELHC10h=0;
fi
if [ "$LHC13d2" = "" ]; then 
    HAVELHC13d2MC=0; 
fi
if [ "$LHC13d2bMC" = "" ]; then 
    HAVELHC13d2b=0; 
fi


if [ $HAVELHC10h == 1 ]; then
    LHC10hData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/ | grep $LHC10hData\_`
    if [ "$LHC10hData" == "" ]; then
        HAVELHC10h=0;
    else 
        OUTPUTDIR_LHC10h=$BASEDIR/$TRAINDIR/$TRAINPATHData-$LHC10hData
    fi
fi

if [ $HAVELHC13d2 == 1 ]; then
    LHC13d2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/ | grep $LHC13d2MC\_`
    if [ "$LHC13d2MC" == "" ]; then
        HAVELHC13d2=0;
    else 
        OUTPUTDIR_LHC13d2=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13d2MC
    fi    
fi  

if [ $HAVELHC13d2b == 1 ]; then
    LHC13d2bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/ | grep $LHC13d2bMC\_`
    if [ "$LHC13d2bMC" == "" ]; then
        HAVELHC13d2b=0;
    else 
        OUTPUTDIR_LHC13d2b=$BASEDIR/$TRAINDIR/$TRAINPATHMC-$LHC13d2bMC
    fi    
fi

if [ $DOWNLOADON == 1 ]; then
    if [ $HAVELHC10h == 1 ]; then
        echo "downloading LHC10h"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC10h "/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHData/$LHC10hData/merge"
    fi    
    if [ $HAVELHC13d2 == 1 ]; then
        echo "downloading LHC13d2"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13d2 "/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC13d2MC/merge"
    fi    
    if [ $HAVELHC13d2b == 1 ]; then
        echo "downloading LHC13d2b"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13d2b "/alice/cern.ch/user/a/alitrain/PWGGA/$TRAINPATHMC/$LHC13d2b/merge"
    fi    
fi
    
if [ $2 = "AOD" ]; then
 
    if [ $HAVELHC13d2 == 1 ]; then
        ls $OUTPUTDIR_LHC13d2/GammaConvV1_*.root > fileLHC13d2.txt
        fileNumbers=`cat fileLHC13d2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC13d2/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13d2_AOD_$number.log\"\)

        done;
    fi
    
    if [ $HAVELHC13d2b == 1 ]; then
        ls $OUTPUTDIR_LHC13d2b/GammaConvV1_*.root > fileLHC13d2b.txt
        fileNumbersb=`cat fileLHC13d2b.txt`
        for fileName in $fileNumbersb; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC13d2b/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13d2b_AOD_$number.log\"\)
        done;
    fi
    
    if [ $HAVELHC10h == 1 ]; then
        ls $OUTPUTDIR_LHC10h/GammaConvV1_*.root > fileLHC10h.txt
        fileNumbersData=`cat fileLHC10h.txt`
        for fileName in $fileNumbersData; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC10h/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_LHC10h-pass2_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_LHC10h-pass2_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC10h_$number.log\"\)
        done;
    fi
    
    if [ $MERGEON == 1 ]; then
        rm $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_*.root
        ls $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 7 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_LHC13d2b_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2b_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_AOD_LHC13d2_$number.root
            fi
        done
    fi
else 

    if [ $HAVELHC13d2 == 1 ]; then
        ls $OUTPUTDIR_LHC13d2/GammaConvV1_*.root > fileLHC13d2.txt
        fileNumbers=`cat fileLHC13d2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC13d2/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13d2_$number.log\"\)

        done;
    fi

    if [ $HAVELHC13d2b == 1 ]; then
        ls $OUTPUTDIR_LHC13d2b/GammaConvV1_*.root > fileLHC13d2b.txt
        fileNumbersb=`cat fileLHC13d2b.txt`
        for fileName in $fileNumbersb; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC13d2b/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2b_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC13d2b_$number.log\"\)
        done;
    fi
   
    if [ $HAVELHC13d2 == 1 ]; then
        ls $OUTPUTDIR_LHC10h/GammaConvV1_*.root > fileLHC10h.txt
        fileNumbersData=`cat fileLHC10h.txt`
        for fileName in $fileNumbersData; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC10h/GammaConvV1_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_LHC10h-pass2_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_GA_PbPb_LHC10h-pass2_$number.root\"\,\"$OUTPUTDIR/CutSelection_LHC10h_$number.log\"\)
        done;
    fi
   
    if [ $MERGEON == 1 ]; then
        rm $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_*.root
        ls $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 6 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2b_$number.root ] && [ -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2_LHC13d2b_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2b_$number.root $OUTPUTDIR/GammaConvV1_GA_PbPb_MC_LHC13d2_$number.root
            fi
        done
    fi
fi