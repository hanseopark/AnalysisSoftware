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
#         root -l -b -q -x ChangeStructureToStandardConvCalo.C\(\"$1\"\,\"$2\"\,\"GammaConvCalo_$3\"\)
#     fi    
}


# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1

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

TRAINDIR=Legotrain-vAN20170411-dirGamma
LHC13bData="593"; 
LHC13cData="595"; 
# LHC13dData="573"; 
# LHC13eData="574"; 
# LHC13fData="575"; 

# LHC13b2_efix_p1MC="889"; 
# LHC13b2_efix_p2MC="890"; 
# LHC13b2_efix_p3MC="891"; 
# LHC13b2_efix_p4MC="892";
# LHC13e7MC="893";
LHC13b2_efix_p1MC="894"; 
LHC13b2_efix_p2MC="895"; 
LHC13b2_efix_p3MC="896"; 
LHC13b2_efix_p4MC="897";
LHC13e7MC="898";

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
    HAVELHC13e7MC=0; 
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
    fi    
    if [ $HAVELHC13c == 1 ]; then
        echo "downloading LHC13c"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13c "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb/$LHC13cData/merge_runlist_1"
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
    fi    
    if [ $HAVELHC13b2efixp2 == 1 ]; then
        echo "downloading LHC13b2_efix_p2"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p2MC/merge"
    fi    
    if [ $HAVELHC13b2efixp3 == 1 ]; then
        echo "downloading LHC13b2_efix_p3"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p3MC/merge"
    fi    
    if [ $HAVELHC13b2efixp4 == 1 ]; then
        echo "downloading LHC13b2_efix_p4"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13b2_efix_p4 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13b2_efix_p4MC/merge"
    fi    
    if [ $HAVELHC13e7 == 1 ]; then
        echo "downloading LHC13e7"
        CopyFileIfNonExisitent $OUTPUTDIR_LHC13e7 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pPb_MC/$LHC13e7MC/merge"
    fi        
fi
    
if [ $HAVELHC13b == 1 ]; then
    ls $OUTPUTDIR_LHC13b/GammaConvCalo_*.root > fileLHC13b.txt
    fileNumbers=`cat fileLHC13b.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_LHC13b_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13c == 1 ]; then
    ls $OUTPUTDIR_LHC13c/GammaConvCalo_*.root > fileLHC13c.txt
    fileNumbers=`cat fileLHC13c.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13c/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_LHC13c_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13d == 1 ]; then
    ls $OUTPUTDIR_LHC13d/GammaConvCalo_*.root > fileLHC13d.txt
    fileNumbers=`cat fileLHC13d.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13d/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_LHC13d_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13e == 1 ]; then
    ls $OUTPUTDIR_LHC13e/GammaConvCalo_*.root > fileLHC13e.txt
    fileNumbers=`cat fileLHC13e.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13e/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_LHC13e_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13f == 1 ]; then
    ls $OUTPUTDIR_LHC13f/GammaConvCalo_*.root > fileLHC13f.txt
    fileNumbers=`cat fileLHC13f.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13f/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_LHC13f-pass$passNr\_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_LHC13f-pass$passNr\_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_LHC13f_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13b2efixp1 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p1/GammaConvCalo_*.root > fileLHC13b2efixp1.txt
    fileNumbers=`cat fileLHC13b2efixp1.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p1/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC13b2_efix_p1_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13b2efixp2 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p2/GammaConvCalo_*.root > fileLHC13b2efixp2.txt
    fileNumbers=`cat fileLHC13b2efixp2.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p2/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p2_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p2_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC13b2_efix_p2_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13b2efixp3 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p3/GammaConvCalo_*.root > fileLHC13b2efixp3.txt
    fileNumbers=`cat fileLHC13b2efixp3.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p3/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p3_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p3_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC13b2_efix_p3_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13b2efixp4 == 1 ]; then
    ls $OUTPUTDIR_LHC13b2_efix_p4/GammaConvCalo_*.root > fileLHC13b2efixp4.txt
    fileNumbers=`cat fileLHC13b2efixp4.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13b2_efix_p4/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p4_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p4_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC13b2_efix_p4_$number.log\"\,2\)
    done;
fi

if [ $HAVELHC13e7 == 1 ]; then
    ls $OUTPUTDIR_LHC13e7/GammaConvCalo_*.root > fileLHC13e7.txt
    fileNumbers=`cat fileLHC13e7.txt`
    for fileName in $fileNumbers; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number
        ChangeStructureIfNeeded $OUTPUTDIR_LHC13e7/GammaConvCalo_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13e7_$number.root $number
        root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_MC_LHC13e7_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaConvCalo_MC_LHC13e7_$number.log\"\,2\)
    done;
fi

mkdir -p CutSelections
mv CutSelection_*.log CutSelections/

if [ $MERGEON == 1 ]; then
    rm $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_*.root
    ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
        echo $number
        ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root
        ls $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root
        if [ -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13c-pass$passNr\_$number.root
        fi
    done
    
    ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 4 | cut -d "." -f1`
        echo $number
        if [ -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root
        fi
    done
    
    ls $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 5 | cut -d "." -f1`
        echo $number
        if [ -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_LHC13e-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13b-pass$passNr\_LHC13c-pass$passNr\_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
        fi
    done

    rm $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_LHC13e-pass$passNr\_*.root
    ls $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 3 | cut -d "." -f1`
        echo $number
        if [ -f $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_LHC13e-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13d-pass$passNr\_$number.root $OUTPUTDIR/GammaConvCalo_LHC13e-pass$passNr\_$number.root
        fi
    done

    rm $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_*.root
    ls $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 6 | cut -d "." -f1`
        echo $number
        if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p2_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p2_$number.root
        fi
    done
    #    
    ls $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 7 | cut -d "." -f1`
        echo $number
        if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p3_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p3_$number.root
        fi
    done
    
    ls $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    for fileName in $filesForMerging; do
        echo $fileName
        number=`echo $fileName  | cut -d "/" -f $((NSlashes-1)) | cut -d "_" -f 8 | cut -d "." -f1`
        echo $number
        if [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_$number.root ] && [ -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p4_$number.root ] ; then
            hadd -f $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_p4_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p1_p2_p3_$number.root $OUTPUTDIR/GammaConvCalo_MC_LHC13b2_efix_p4_$number.root
        fi
    done

    
fi

   
