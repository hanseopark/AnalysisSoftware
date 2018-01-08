#! /bin/bash

# copies files from grid
# creates directory
# changes internal structure

#! /bin/bash
source basicFunction.sh

# switches to enable/disable certain procedures
DOWNLOADON=1
MERGEON=1
MERGEONBINSSingle=1
MERGEONBINS=1
SEPARATEON=0

# check if train configuration has actually been given
HAVELHC11a=1
HAVELHC12f1a=1
HAVELHC12f1b=1
HAVELHC12i3=1
HAVELHC15g1a=1
HAVELHC13g=1
HAVELHC15g2=1
HAVELHC15a3a=1
HAVELHC15a3aplus=1

# default trainconfigurations
LHC11aData="";
LHC12f1aMC="";
LHC12f1bMC="";
LHC12i3MC="";
LHC15g1aMC="";
LHC13gData="";
LHC15a3aMC="";
LHC15a3aplusMC="";
LHC15g2MC="";


NSlashes=10
NSlashes2=9
NSlashes3=11

if [ $1 = "fbock" ]; then
    BASEDIR=/mnt/additionalStorage/OutputLegoTrains/pp
    NSlashes=8
    NSlashes2=7
elif [ $1 = "fbockGSI" ]; then
    BASEDIR=/hera/alice/fbock/Grid/OutputLegoTrains/pp
elif [ $1 = "leardini" ]; then
    BASEDIR=/Users/lucy/
elif [ $1 = "leardiniALICESERV1" ]; then
    BASEDIR=/alidata50/alice_u/leardini/GridOutput/PbPb/
elif [ $1 = "leardiniGSI" ]; then
    BASEDIR=/hera/alice/leardini/Grid/OutputLegoTrains/pp
elif [ $1 = "passfeld" ]; then
    BASEDIR=~/work/Gridoutput/pp
elif [ $1 = "passfeldMAF" ]; then
    BASEDIR=/data9/a_pass02/gamma_test/AnalysisSoftware/LegoTrain/
elif [ $1 = "passfeldGSI" ]; then
    BASEDIR=/hera/alice/passfeld/Grid/OutputLegoTrains/pp
elif [ $1 = "amarin" ]; then
    BASEDIR=/Users/marin/
elif [ $1 = "amarinGSI" ]; then
    BASEDIR=/hera/alice/marin/Grid/OutputLegoTrains/pp
elif [ $1 = "amarinALICESERV1" ]; then
    BASEDIR=/alidata50/alice_u/amarin/GridOutput/PbPb/
elif [ $1 = "mwilde" ]; then
    BASEDIR=~/work/GridOutput
elif [ $1 = "mwildeGSI" ]; then
    BASEDIR=/hera/alice/mwilde/Grid/OutputLegoTrains/pp
elif [ $1 = "pgonzales" ]; then
    BASEDIR=~/work/GridOutput
elif [ $1 = "pgonzalesGSI" ]; then
    BASEDIR=/hera/alice/pgonzales/Grid/OutputLegoTrains/pp
elif [ $1 = "dmuhlhei" ]; then
    BASEDIR=~/data/work/Grid
    NSlashes=9
    NSlashes2=8
elif [ $1 = "bsahlmul" ]; then
    BASEDIR=/Users/sahlmul/Documents/ALICE/mycode/AnalysisSoftware/grid
    NSlashes=11
    NSlashes2=10
fi


# TRAINDIR=Legotrain-vAN20170417_DirGamma
# LHC11aData="2046";
# LHC12f1aMC="2883";
# LHC12f1bMC="2884";

# TRAINDIR=Legotrain-vAN20170830_DirGamma
# LHC11aData="2239";
# LHC12f1aMC="3107";
# LHC12f1bMC="3109";
# LHC12f1aMC="3110";
# LHC12f1bMC="3111";
# LHC12f1aMC="3120";
# LHC12f1bMC="3122";
# LHC12f1aMC="3121";
# LHC12f1bMC="3123";

# TRAINDIR=Legotrain-vAN20170905_DirGamma
# LHC11aData="2241";
# LHC12f1aMC="3128";
# LHC12f1bMC="3129";

TRAINDIR=Legotrain-vAN20171204_DirGamma
# LHC11aData="2266";
# LHC12f1aMC="3189";
# LHC12f1bMC="3190";
# LHC11aData="2268";
# LHC12f1aMC="3194";
# LHC12f1bMC="3195";
# LHC15g1aMC="3196";
# LHC15g1aMC="3200";
LHC11aData="2272";
LHC12f1aMC="3201";
LHC12f1bMC="3202";


OUTPUTDIR=$BASEDIR/$TRAINDIR
mkdir -p $OUTPUTDIR/CutSelections

echo "************************************************************************************************";
echo "Directory to be writen into " $OUTPUTDIR
echo "************************************************************************************************";
echo ""
echo ""


if [ $2 = "LHC11a" ]; then
    echo "************************************************************************************************";
    echo "********************************* Copying LHC11a ***********************************************";
    echo "************************************************************************************************";

    if [ "$LHC11aData" == "" ]; then
        HAVELHC11a=0;
    fi
    if [ "$LHC12f1aMC" = "" ]; then
        HAVELHC12f1a=0;
    fi
    if [ "$LHC12f1bMC" = "" ]; then
        HAVELHC12f1b=0;
    fi
    if [ "$LHC12i3MC" = "" ]; then
        HAVELHC12i3=0;
    fi
    if [ "$LHC15g1aMC" = "" ]; then
        HAVELHC15g1a=0;
    fi

    if [ $HAVELHC11a == 1 ]; then
        LHC11aData=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/ | grep $LHC11aData\_`
        if [ "$LHC11aData" == "" ]; then
            HAVELHC11a=0;
        else
            OUTPUTDIR_LHC11a=$BASEDIR/$TRAINDIR/GA_pp-$LHC11aData
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

    if [ $HAVELHC12i3 == 1 ]; then
        LHC12i3MC=`alien_ls /alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/ | grep $LHC12i3MC\_`
        if [ "$LHC12i3MC" == "" ]; then
            HAVELHC12i3=0;
        else
            OUTPUTDIR_LHC12i3=$BASEDIR/$TRAINDIR/GA_pp_MC-$LHC12i3MC
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

    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC11a == 1 ]; then
            echo "downloading LHC11a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC11a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC11aData/merge_runlist_1" $NSlashes3 "none" kFALSE
        fi
        if [ $HAVELHC12f1a == 1 ]; then
            echo "downloading LHC12f1a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1aMC/merge_runlist_1" $NSlashes3 "none" kFALSE
        fi
        if [ $HAVELHC12f1b == 1 ]; then
            echo "downloading LHC12f1b"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12f1b "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12f1bMC/merge_runlist_1" $NSlashes3 "none" kFALSE
        fi
        if [ $HAVELHC12i3 == 1 ]; then
            echo "LHC12i3"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC12i3 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC12i3MC/merge_runlist_1" $NSlashes3 "none" kFALSE
        fi
        if [ $HAVELHC15g1a == 1 ]; then
            echo "LHC15g1a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g1aMC/merge" $NSlashes3 "none" kFALSE

            echo "copying LHC11aJetJet"
            runNumbers=`cat runlists/runNumbersLHC11aJetJet.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                echo $runNumber
                binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    CopyFileIfNonExisitent $OUTPUTDIR_LHC15g1a/$binNumber/$runNumber "/alice/sim/2015/LHC15g1a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15g1aMC" $NSlashes3 "none" kFALSE
                done;
            done;
        fi
    fi

    if [ $HAVELHC11a == 1 ]; then
        ls $OUTPUTDIR_LHC11a/GammaCalo_*.root > fileLHC11a.txt
        fileNumbers=`cat fileLHC11a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC11a $NSlashes "LHC11a-pass4" ""
        done;
    fi

    if [ $HAVELHC12f1a == 1 ]; then
        ls $OUTPUTDIR_LHC12f1a/GammaCalo_*.root > fileLHC12f1a.txt
        fileNumbers=`cat fileLHC12f1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC12f1a $NSlashes "MC_LHC12f1a" ""
        done;
    fi

    if [ $HAVELHC12f1b == 1 ]; then
        ls $OUTPUTDIR_LHC12f1b/GammaCalo_*.root > fileLHC12f1b.txt
        fileNumbers=`cat fileLHC12f1b.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC12f1b $NSlashes "MC_LHC12f1b" ""
        done;
    fi

    if [ $HAVELHC12i3 == 1 ]; then
        ls $OUTPUTDIR_LHC12i3/GammaCalo_*.root > fileLHC12i3.txt
        fileNumbers=`cat fileLHC12i3.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC12i3 $NSlashes "MC_LHC12i3" ""
        done;
    fi

    if [ $HAVELHC15g1a == 1 ]; then
        ls $OUTPUTDIR_LHC15g1a/GammaCalo_*.root > fileLHC15g1a.txt
        fileNumbers=`cat fileLHC15g1a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC15g1a $NSlashes "MC_LHC15g1a" ""
        done;

        if [ $MERGEONBINSSingle = 1 ]; then
            binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
            echo $binNumbersJJ
            ls $OUTPUTDIR_LHC15g1a/GammaCalo_*.root > filetemp.txt
            mkdir -p $OUTPUTDIR/LHC15g1aFineBins
            for binNumber in $binNumbersJJ; do
                echo $binNumber
                fileNumbers=`cat filetemp.txt`
                for fileName in $fileNumbers; do
                    echo $fileName
                    number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                    echo $number
                    hadd -f $OUTPUTDIR_LHC15g1a/$binNumber/GammaCalo_$number.root $OUTPUTDIR_LHC15g1a/$binNumber/*/GammaCalo_$number.root
                    cp $OUTPUTDIR_LHC15g1a/$binNumber/GammaCalo_$number.root $OUTPUTDIR/LHC15g1aFineBins/GammaCalo_MC_LHC15g1a$binNumber\_$number.root $number
                done;
            done;
        fi
    fi

    if [ $MERGEON == 1 ]; then
        rm $OUTPUTDIR/GammaCalo_MC_LHC12f1a_LHC12f1b_*.root
        ls $OUTPUTDIR/GammaCalo_MC_LHC12f1a_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
        filesForMerging=`cat filesForMerging.txt`
            echo $number
            if [ -f $OUTPUTDIR/GammaCalo_MC_LHC12f1a_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC12f1b_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCalo_MC_LHC12f1a_LHC12f1b_$number.root $OUTPUTDIR/GammaCalo_MC_LHC12f1a_$number.root $OUTPUTDIR/GammaCalo_MC_LHC12f1b_$number.root
            fi
        done

        ls $OUTPUTDIR/GammaCalo_MC_LHC12f1a_*.root > filesForMerging.txt
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaCalo_MC_LHC12f1a_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC12i3_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCalo_MC_LHC12f1a_LHC12i3_$number.root $OUTPUTDIR/GammaCalo_MC_LHC12i3_$number.root $OUTPUTDIR/GammaCalo_MC_LHC12f1a_$number.root
            fi
        done


        ls $OUTPUTDIR/GammaCalo_MC_LHC12f1a_LHC12f1b_*.root > filesForMerging.txt
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 5 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaCalo_MC_LHC12f1a_LHC12f1b_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC12i3_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCalo_MC_LHC12f1a_LHC12f1b_LHC12i3_$number.root $OUTPUTDIR/GammaCalo_MC_LHC12f1a_LHC12f1b_$number.root $OUTPUTDIR/GammaCalo_MC_LHC12i3_$number.root
            fi
        done
    fi

    ls $OUTPUTDIR/GammaCalo_MC_LHC15g1a_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    if [ $MERGEONBINS == 1 ]; then
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15g1aFineBins/GammaCalo_MC_LHC15g1a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15g1aFineBins/GammaCalo_MC_LHC15g1a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15g1aFineBins/GammaCalo_MC_LHC15g1a$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15g1aFineBins/GammaCalo_MC_LHC15g1a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCalo_MC_LHC15g1aFinerPtHardBins_$number.root $TOMERGE
        done;
    fi

elif [ $2 = "LHC13g" ]; then

    if [ "$LHC13gData" == "" ]; then
        HAVELHC13g=0;
    fi

    if [ "$LHC15a3aMC" == "" ]; then
        HAVELHC15a3a=0;
    fi

    if [ "$LHC15a3aplusMC" == "" ]; then
        HAVELHC15a3aplus=0;
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
        echo $LHC15g2MC "   -> full dir"
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

    if [ $DOWNLOADON == 1 ]; then
        if [ $HAVELHC13g == 1 ]; then
            echo "downloading LHC13g"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC13g "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp/$LHC13gData/merge_runlist_1" $NSlashes3 "none" kFALSE
            runNumbers=`cat runNumbersLHC13g_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $OUTPUTDIR_LHC13g/$runNumber "/alice/data/2013/LHC13g/000$runNumber/pass1/PWGGA/GA_pp/$LHC13gData" $NSlashes3 "none" kFALSE
            done;
        fi

        if [ $HAVELHC15g2 == 1 ]; then
            echo "downloading LHC15g2"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15g2 "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15g2MC/merge" $NSlashes3 "none" kFALSE
        fi

        if [ $HAVELHC15a3a == 1 ]; then
            echo "downloading LHC15a3a"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aMC/merge" $NSlashes3 "none" kFALSE
        fi

        if [ $HAVELHC15a3aplus == 1 ]; then
            echo "downloading LHC15a3a_plus"
            CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus "/alice/cern.ch/user/a/alitrain/PWGGA/GA_pp_MC/$LHC15a3aplusMC/merge" $NSlashes3 "none" kFALSE
        fi

        echo "copying LHC13gJetJet in bins"
        runNumbers=`cat runlists/runNumbersLHC13gJetJet.txt`
        echo $runNumbers
        for runNumber in $runNumbers; do
            echo $runNumber
            binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
            if [ $DOWNLOADON = 1 ]; then
                for binNumber in $binNumbersJJ; do
                    echo $binNumber
                    if [ $HAVELHC15a3a == 1 ]; then
                        CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3a/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aMC" $NSlashes3 "none" kFALSE
                    fi
                    if [ $HAVELHC15a3aplus == 1 ]; then
                        mkdir -p $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber
                        CopyFileIfNonExisitent $OUTPUTDIR_LHC15a3aplus/$binNumber/$runNumber "/alice/sim/2015/LHC15a3a_plus/$binNumber/$runNumber/PWGGA/GA_pp_MC/$LHC15a3aplusMC" $NSlashes3 "none" kFALSE
                    fi
                done;
            fi
        done;
    fi

    if [ $HAVELHC13g == 1 ]; then
        ls $OUTPUTDIR_LHC13g/GammaCalo_*.root > fileLHC13g.txt
        fileNumbers=`cat fileLHC13g.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
            echo $number
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC13g $NSlashes "LHC13g-pass1" ""
            mkdir -p $OUTPUTDIR/LHC13gRunWise
            runNumbers=`cat runlists/runNumbersLHC13g_pass1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                cp $OUTPUTDIR_LHC13g/$runNumber/GammaCalo_$number.root $OUTPUTDIR/LHC13gRunWise/GammaCalo_LHC13g-pass1_$runNumber\_$number.root $number
            done;

        done;
    fi

    if [ $HAVELHC15g2 == 1 ]; then
        ls $OUTPUTDIR_LHC15g2/GammaCalo_*.root > fileLHC15g2.txt
        fileNumbers=`cat fileLHC15g2.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC15g2 $NSlashes "MC_LHC15g2" ""
        done;
    fi

    if [ $HAVELHC15a3a == 1 ]; then
        ls $OUTPUTDIR_LHC15a3a/GammaCalo_*.root > fileLHC15a3a.txt
        fileNumbers=`cat fileLHC15a3a.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC15a3a $NSlashes "MC_LHC15a3a" ""
        done;
    fi

    if [ $HAVELHC15a3aplus == 1 ]; then
        ls $OUTPUTDIR_LHC15a3aplus/GammaCalo_*.root > fileLHC15a3aplus.txt
        fileNumbers=`cat fileLHC15a3aplus.txt`
        for fileName in $fileNumbers; do
            echo $fileName
            ChangeStructureIfNeededCalo $fileName $OUTPUTDIR_LHC15a3aplus $NSlashes "MC_LHC15a3aplus" ""
        done;
    fi

    if [ $MERGEONBINSSingle = 1 ]; then
        binNumbersJJ=`cat runlists/binNumbersJJToMerge.txt`
        echo $binNumbersJJ
        ls $OUTPUTDIR_LHC15a3a/GammaCalo_*.root > filetemp.txt
        mkdir -p $OUTPUTDIR/LHC15a3aXFineBins
        for binNumber in $binNumbersJJ; do
            echo $binNumber
            fileNumbers=`cat filetemp.txt`
            for fileName in $fileNumbers; do
                echo $fileName
                number=`echo $fileName  | cut -d "/" -f $NSlashes | cut -d "_" -f 2 | cut -d "." -f1`
                echo $number
                if [ $HAVELHC15a3a == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC15a3a/$binNumber/GammaCalo_$number.root $OUTPUTDIR_LHC15a3a/$binNumber/*/GammaCalo_$number.root
                    cp $OUTPUTDIR_LHC15a3a/$binNumber/GammaCalo_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3a$binNumber\_$number.root $number
                fi
                if [ $HAVELHC15a3aplus == 1 ]; then
                    hadd -f $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaCalo_$number.root $OUTPUTDIR_LHC15a3aplus/$binNumber/*/GammaCalo_$number.root
                    cp $OUTPUTDIR_LHC15a3aplus/$binNumber/GammaCalo_$number.root $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3aplus$binNumber\_$number.root $number
                fi
            done;
        done;
    fi


    rm $OUTPUTDIR/GammaCalo_LHC13gRed_*.root
    ls $OUTPUTDIR/GammaCalo_LHC13g-pass1_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    if [ $MERGEON == 1 ]; then
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 3 | cut -d "." -f1`
            echo $number
            runsForMerging=`cat runlists/runNumbersLHC13g_pass1_reduced.txt`
            TOMERGE="";
            for run in $runsForMerging; do
                echo $OUTPUTDIR/LHC13gRunWise/GammaCalo_LHC13g-pass1_$run\_$number.root
                if [ -f $OUTPUTDIR/LHC13gRunWise/GammaCalo_LHC13g-pass1_$run\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC13gRunWise/GammaCalo_LHC13g-pass1_$run"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for run $run, number $number, $OUTPUTDIR/LHC13gRunWise/GammaCalo_LHC13g-pass1_$run\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCalo_LHC13g-pass1-reducedRunList_$number.root $TOMERGE
        done

    fi

    rm $OUTPUTDIR/GammaCalo_MC_LHC15a3a_LHC15a3aplus_*.root
    ls $OUTPUTDIR/GammaCalo_MC_LHC15a3a_*.root > filesForMerging.txt
    filesForMerging=`cat filesForMerging.txt`
    if [ $MERGEON == 1 ]; then
        for fileName in $filesForMerging; do
            echo $fileName
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            echo $number
            if [ -f $OUTPUTDIR/GammaCalo_MC_LHC15a3a_$number.root ] && [ -f $OUTPUTDIR/GammaCalo_MC_LHC15a3aplus_$number.root ] ; then
                hadd -f $OUTPUTDIR/GammaCalo_MC_LHC15a3a_LHC15a3aplus_$number.root $OUTPUTDIR/GammaCalo_MC_LHC15a3a_$number.root $OUTPUTDIR/GammaCalo_MC_LHC15a3aplus_$number.root
            fi
        done
    fi

    if [ $MERGEONBINS == 1 ]; then
        rm $OUTPUTDIR/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_*.root
        filesForMerging=`cat filesForMerging.txt`
        for fileName in $filesForMerging; do
            binsForMerging=`cat runlists/binNumbersJJToMerge.txt`
            number=`echo $fileName  | cut -d "/" -f $NSlashes2 | cut -d "_" -f 4 | cut -d "." -f1`
            TOMERGE="";
            for bin in $binsForMerging; do
                echo $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3a$bin\_$number.root
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3a$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3a$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number, $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3a$bin\_$number.root";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCalo_MC_LHC15a3aFinerPtHardBins_$number.root $TOMERGE

            TOMERGE="";
            for bin in $binsForMerging; do
                if [ -f $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3aplus$bin\_$number.root ]; then
                    TOMERGE="$TOMERGE $OUTPUTDIR/LHC15a3aXFineBins/GammaCalo_MC_LHC15a3aplus$bin"
                    TOMERGE+="_$number.root"
                else
                    echo "I couldn't find the file for bin $bin, number $number";
                fi
            done;
            hadd -f $OUTPUTDIR/GammaCalo_MC_LHC15a3aplusFinerPtHardBins_$number.root $TOMERGE

            hadd -f $OUTPUTDIR/GammaCalo_MC_LHC15a3aFinerPtHardBins_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaCalo_MC_LHC15a3aplusFinerPtHardBins_$number.root $OUTPUTDIR/GammaCalo_MC_LHC15a3aFinerPtHardBins_$number.root
        done;
    fi
fi