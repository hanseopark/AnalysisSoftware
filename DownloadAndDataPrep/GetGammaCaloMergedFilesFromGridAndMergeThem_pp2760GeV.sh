#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure
# merges files according to the pPb needs
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
}

NSlashes=10
NSlashes2=9

# switches to enable/disable certain procedures
DOWNLOADON=0
MERGEON=1
MERGEONBINSSingle=0
MERGEONBINS=0

# check if train configuration has actually been given
HAVELHC11a=1
HAVELHC15g1a=1
HAVELHC15g1b=1
HAVELHC13g=1
HAVELHC15a3a=1
HAVELHC15a3aplus=1
HAVELHC15a3b=1
HAVELHC12f1a=1
HAVELHC12f1b=1
HAVELHC15g2=1

# default trainconfigurations
LHC11aData="";
LHC12f1aMC="";
LHC12f1bMC="";
LHC15g1aMC="";
LHC15g1bMC="";
LHC13gData="";
LHC15a3aMC="";
LHC15a3aplusMC="";
LHC15a3bMC="";
LHC15g2MC="";

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp
    NSlashes=8
    NSlashes2=7
elif [ $1 = "dmuhlheim" ]; then
   BASEDIR=/home/daniel/Desktop/Grid
   NSlashes=9
   NSlashes=8
fi

# TRAINDIR=Legotrain-vAN20161029_TMEffi
# LHC11aData="1905";
# LHC15g1aMC="2604";
#
# LHC13gData="1907";
# LHC15a3aMC="2608";
# LHC15a3aplusMC="2609";

# TRAINDIR=Legotrain-vAN20161111_TMEffi
# LHC11aData="1905";
# LHC11aData="1945";
# LHC11aData="1949";
# LHC11aData="1965";
# LHC11aData="1894";
# LHC15g1aMC="2653";
# LHC15g1aMC="2669";
# LHC15g1aMC="2670";
# LHC15g1aMC="2689";
# LHC15g1aMC="2697";
# LHC15g1aMC="2722";

# LHC13gData="1907";
# LHC13gData="1946";
# LHC13gData="1951";
# LHC13gData="1961";
# LHC13gData="1895";
# LHC15a3aMC="2650";
# LHC15a3aplusMC="2651";
# LHC15a3aMC="2673";
# LHC15a3aplusMC="2675";
# LHC15a3aMC="2674";
# LHC15a3aplusMC="2676";
# LHC15a3aMC="2687";
# LHC15a3aplusMC="2690";
# LHC15a3aMC="2698";
# LHC15a3aplusMC="2699";
# LHC15a3aMC="2710";
# LHC15a3aplusMC="2711";

TRAINDIR=Legotrain-vAN20170329_TMEffiMCfix
LHC15g1aMC="2866";
LHC12f1aMC="2870";
LHC12f1bMC="2869";

# LHC15g2MC="2841";
LHC15g2MC="2842";
LHC15a3aMC="2867";
LHC15a3aplusMC="2868";


OUTPUTDIR=$BASEDIR/$TRAINDIR

if [ $2 = "LHC11a" ]; then

    if [ "$LHC11aData" == "" ]; then
        HAVELHC11a=0;
    fi
    if [ "$LHC15g1aMC" = "" ]; then
        HAVELHC15g1a=0;
    fi
    if [ "$LHC15g1bMC" = "" ]; then
        HAVELHC15g1b=0;
    fi
    if [ "$LHC12f1aMC" = "" ]; then
        HAVELHC12f1a=0;
    fi
    if [ "$LHC12f1bMC" = "" ]; then
        HAVELHC12f1b=0;
    fi


    if [ $HAVELHC11a == 1 ]; then
        LHC11aData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC11aData\_`
        if [ "$LHC11aData" == "" ]; then
            HAVELHC11a=0;
        else
            OUTPUTDIR_LHC11a=$BASEDIR/$TRAINDIR/GA_pp-$LHC11aData
        fi
    fi
    if [ $HAVELHC15g1a == 1 ]; then
        LHC15g1aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g1aMC\_`
        if [ "$LHC15g1aMC" == "" ]; then
            HAVELHC15g1a=0;
        else
            OUTPUTDIR_LHC15g1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g1aMC
        fi
    fi
    if [ $HAVELHC15g1b == 1 ]; then
       LHC15g1bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g1bMC\_`
        if [ "$LHC15g1bMC" == "" ]; then
            HAVELHC15g1b=0;
        else
            OUTPUTDIR_LHC15g1b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g1bMC
        fi
    fi
    if [ $HAVELHC12f1a == 1 ]; then
        LHC12f1aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12f1aMC\_`
        if [ "$LHC12f1aMC" == "" ]; then
            HAVELHC12f1a=0;
        else
            OUTPUTDIR_LHC12f1a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1aMC
        fi
    fi

    if [ $HAVELHC12f1b == 1 ]; then
        LHC12f1bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12f1bMC\_`
        if [ "$LHC12f1bMC" == "" ]; then
            HAVELHC12f1b=0;
        else
            OUTPUTDIR_LHC12f1b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12f1bMC
        fi
    fi


    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC11a == 1 ]; then
            echo "downloading LHC11a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC11a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC11aData/merge_runlist_1"
#             runNumbers=`cat runlists/runNumbersLHC11a_pass4_wSDD.txt`
#             echo $runNumbers
#             for runNumber in $runNumbers; do
#                 echo $runNumber
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC11a/$runNumber "/alice/data/2011/LHC11a/000$runNumber/ESDs/pass4_with_SDD/PWGGA/GA_pp/$LHC11aData"
#             done;
        fi
        if [ $HAVELHC12f1a == 1 ]; then
            echo "downloading LHC12f1a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1aMC/merge_runlist_1"
        fi
        if [ $HAVELHC12f1b == 1 ]; then
            echo "downloading LHC12f1b"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1bMC/merge_runlist_1"
        fi

        if [ $HAVELHC15g1a == 1 ]; then
            echo "downloading LHC15g1a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g1aMC/merge"

            echo "copying LHC11aJetJet"
            runNumbers=`cat runlists/runNumbersLHC11aJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1a/$binNumber/$runNumber "/alice/sim/2015/LHC15g1a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15g1aMC"
                done;
            done;
        fi
        if [ $HAVELHC15g1b == 1 ]; then
            echo "downloading LHC15g1b"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g1bMC/merge"
            echo "copying LHC11aJetJet - gamma"
            runNumbers=`cat runlists/runNumbersLHC11aJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJGamma=`cat runlists/binNumbersJJGamma.txt`
                for binNumber in $binNumbersJJGamma; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1b/$binNumber/$runNumber "/alice/sim/2015/LHC15g1b/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15g1bMC"
                done;
            done;
        fi

    fi

    if [ $HAVELHC11a == 1 ]; then
        ls $OUTPUTDIR_LHC11a/GammaCaloMerged_*.root > fileLHC11a.txt
        fileNumbers=`cat fileLHC11a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC11a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_LHC11a-pass4_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_LHC11a-pass4_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaMergedCalo_LHC11a_$number.log\"\,10\)
        done;
    fi

    if [ $HAVELHC12f1a == 1 ]; then
        ls $OUTPUTDIR_LHC12f1a/GammaCaloMerged_*.root > fileLHC12f1a.txt
        fileNumbers=`cat fileLHC12f1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC12f1a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC12f1a_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC12f1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCaloMerged_MC_LHC12f1a_$number.log\"\,10\)
        done;
    fi

    if [ $HAVELHC12f1b == 1 ]; then
        ls $OUTPUTDIR_LHC12f1b/GammaCaloMerged_*.root > fileLHC12f1b.txt
        fileNumbers=`cat fileLHC12f1b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC12f1b/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC12f1b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC12f1b_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCaloMerged_MC_LHC12f1b_$number.log\"\,10\)
        done;
    fi


    if [ $HAVELHC15g1a == 1 ]; then
        ls $OUTPUTDIR_LHC15g1a/GammaCaloMerged_*.root > fileLHC15g1a.txt
        fileNumbers=`cat fileLHC15g1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1a_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15g1a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaMergedCalo_MC_LHC15g1a_$number.log\"\,10\)
        done;

        binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        if [ $MERGEONBINSSingle = 1 ]; then
            ls $OUTPUTDIR_LHC15g1a/GammaCaloMerged_*.root > filetemp.txt
            mkdir $OUTPUTDIR/LHC15g1aFineBins
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    hadd -f $OUTPUTDIR_LHC15g1a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15g1a/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$binNumber\_$number.root $number
                done;
            done;
        fi
    fi

    if [ $HAVELHC15g1b == 1 ]; then
        ls $OUTPUTDIR_LHC15g1b/GammaCaloMerged_*.root > fileLHC15g1b.txt
        fileNumbers=`cat fileLHC15g1b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1b/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15g1b_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaMergedCalo_MC_LHC15g1b_$number.log\"\,10\)
        done;
        if [ $MERGEONBINSSingle = 1 ]; then
            binNumbersJJGamma=`cat runlists/binNumbersJJGamma.txt`
            ls $OUTPUTDIR_LHC15g1b/GammaCaloMerged_*.root > filetemp.txt
            mkdir $OUTPUTDIR/LHC15g1aFineBins
            for binNumber in $binNumbersJJGamma; do
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    hadd -f $OUTPUTDIR_LHC15g1b/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15g1b/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15g1b/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1b$binNumber\_$number.root $number
                done;
            done;
        fi
    fi

    if [ $MERGEONBINS == 1 ]; then
        rm $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1aFinerPtHardBins_*.root
        ls $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15g1aFineBins/GammaCaloMerged_MC_LHC15g1a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15g1aFinerPtHardBins_$number.root $TOMERGE
        done;
    fi

elif [ $2 = "LHC13g" ]; then

    if [ "$LHC13gData" == "" ]; then
        HAVELHC13g=0;
    fi
    if [ "$LHC15a3aMC" = "" ]; then
        HAVELHC15a3a=0;
    fi
    if [ "$LHC15a3aplusMC" = "" ]; then
        HAVELHC15a3aplus=0;
    fi
    if [ "$LHC15a3bMC" = "" ]; then
        HAVELHC15a3b=0;
    fi
    if [ "$LHC15g2MC" == "" ]; then
        HAVELHC15g2=0;
    fi


    if [ $HAVELHC13g == 1 ]; then
        LHC13gData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC13gData\_`
        if [ "$LHC13gData" == "" ]; then
            HAVELHC13g=0;
        else
            OUTPUTDIR_LHC13g=$BASEDIR/$TRAINDIR/GA_pp-$LHC13gData
        fi
    fi

    if [ $HAVELHC15g2 == 1 ]; then
        LHC15g2MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15g2MC\_`
        if [ "$LHC15g2MC" == "" ]; then
            HAVELHC15g2=0;
        else
            OUTPUTDIR_LHC15g2=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15g2MC
        fi
    fi

    if [ $HAVELHC15a3a == 1 ]; then
        LHC15a3aMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3aMC\_`
        if [ "$LHC15a3aMC" == "" ]; then
            HAVELHC15a3a=0;
        else
            OUTPUTDIR_LHC15a3a=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aMC
        fi
    fi

    if [ $HAVELHC15a3aplus == 1 ]; then
        LHC15a3aplusMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3aplusMC\_`
        if [ "$LHC15a3aplusMC" == "" ]; then
            HAVELHC15a3aplus=0;
        else
            OUTPUTDIR_LHC15a3aplus=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3aplusMC
        fi
    fi

    if [ $HAVELHC15a3b == 1 ]; then
        LHC15a3bMC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC15a3bMC\_`
        if [ "$LHC15a3bMC" == "" ]; then
            HAVELHC15a3b=0;
        else
            OUTPUTDIR_LHC15a3b=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC15a3bMC
        fi
    fi

    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC13g == 1 ]; then
            echo "copying LHC13g"
#             CopyFileIfNonExisitent $OUTPUTDIR_LHC13g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC13gData/merge"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC13gData/merge_runlist_1"
            runNumbers=`cat runlists/runNumbersLHC13g_pass1.txt`
            echo $runNumbers
#             for runNumber in $runNumbers; do
#                 echo $runNumber
#                 CopyFileIfNonExisitent $OUTPUTDIR_LHC13g/$runNumber "/alice/data/2013/LHC13g/000$runNumber/pass1/PWGGA/GA_pp/$LHC13gData"
#             done;
        fi

        if [ $HAVELHC15g2 == 1 ]; then
            echo "downloading LHC15g2"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g2MC/merge"
        fi

        if [ $HAVELHC15a3a == 1 ]; then
            echo "copying LHC15a3a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aMC/merge"

            echo "copying LHC15a3a single bins"
            runNumbers=`cat runlists/runNumbersLHC13gJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aMC"
                done;
            done;

        fi

        if [ $HAVELHC15a3aplus == 1 ]; then
            echo "copying LHC15a3a_plus"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aplusMC/merge"
            echo "copying LHC15a3a_plus single bins"
            runNumbers=`cat runlists/runNumbersLHC13gJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a_plus/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aplusMC"
                done;
            done;
        fi

        if [ $HAVELHC15a3b == 1 ]; then
            echo "copying LHC15a3b"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3bMC/merge"
            echo "copying LHC15a3b single bins"
            runNumbers=`cat runlists/runNumbersLHC13gJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJGamma=`cat runlists/binNumbersJJGamma.txt`
                for binNumber in $binNumbersJJGamma; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3b/$binNumber/$runNumber "/alice/sim/2015/LHC15a3b/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3bMC"
                done;
            done;
        fi
    fi

    if [ $HAVELHC13g == 1 ]; then
        ls $OUTPUTDIR_LHC13g/GammaCaloMerged_*.root > fileLHC13g.txt
        fileNumbers=`cat fileLHC13g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC13g/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_LHC13g-pass1_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_LHC13g-pass1_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaMergedCalo_LHC13g_$number.log\"\,10\)
        done;
    fi

    if [ $HAVELHC15g2 == 1 ]; then
        ls $OUTPUTDIR_LHC15g2/GammaCaloMerged_*.root > fileLHC15g2.txt
        fileNumbers=`cat fileLHC15g2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15g2/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15g2_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15g2_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaCaloMerged_MC_LHC15g2_$number.log\"\,10\)
        done;
    fi


    if [ $HAVELHC15a3a == 1 ]; then
        binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        ls $OUTPUTDIR_LHC15a3a/GammaCaloMerged_*.root > filetemp.txt
        mkdir -p $OUTPUTDIR/LHC15a3aXFineBins
        if [ $MERGEONBINSSingle = 1 ]; then
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number

                    hadd -f $OUTPUTDIR_LHC15a3a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15a3a/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3a/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$binNumber\_$number.root $number
                done;
            done;
        fi

        ls $OUTPUTDIR_LHC15a3a/GammaCaloMerged_*.root > fileLHC15a3a.txt
        fileNumbers=`cat fileLHC15a3a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3a/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root  $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaMergedCalo_MC_LHC15a3a_$number.log\"\,10\)
        done;
    fi

    if [ $HAVELHC15a3aplus == 1 ]; then
        binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        mkdir -p $OUTPUTDIR/LHC15a3aXFineBins
        ls $OUTPUTDIR_LHC15a3aplus/GammaCaloMerged_*.root > filetemp.txt
        if [ $MERGEONBINSSingle = 1 ]; then
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    hadd -f $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR_LHC15a3aplus/$binNumber/*/GammaCaloMerged_$number.root
                    ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaCaloMerged_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$binNumber\_$number.root $number
                done;
            done;
        fi

        ls $OUTPUTDIR_LHC15a3aplus/GammaCaloMerged_*.root > fileLHC15a3aplus.txt
        fileNumbers=`cat fileLHC15a3aplus.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3aplus/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaMergedCalo_MC_LHC15a3aplus_$number.log\"\,10\)
        done;
    fi

    if [ $HAVELHC15a3b == 1 ]; then
        ls $OUTPUTDIR_LHC15a3b/GammaCaloMerged_*.root > fileLHC15a3b.txt
        fileNumbers=`cat fileLHC15a3b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeeded $OUTPUTDIR_LHC15a3b/GammaCaloMerged_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3b_$number.root $number
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_MC_LHC15a3b_$number.root\"\,\"$OUTPUTDIR/CutSelection_GammaMergedCalo_MC_LHC15a3b_$number.log\"\,10\)
        done;
    fi

    if [ $MERGEON = 1 ]; then
        rm $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_LHC15a3aplus_*.root
        ls $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root ] && [ -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_LHC15a3aplus_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplus_$number.root
            fi
        done
    fi

    if [ $MERGEONBINS = 1 ]; then
        rm $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_*.root
        ls $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_$number.root $TOMERGE

            TOMERGE="";
            for bin in $binsForMerging; do
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15a3aXFineBins/GammaCaloMerged_MC_LHC15a3aplus$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplusFinerPtHardBins_$number.root $TOMERGE

            hadd -f $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaCaloMerged_MC_LHC15a3aFinerPtHardBins_$number.root
        done;
    fi
fi