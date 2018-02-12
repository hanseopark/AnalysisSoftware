#! /bin/bash

# This script has to be run with "bash"

# switches to enable/disable certain procedures
DOWNLOADON=0
MERGEON=0
SINGLERUN=0
SEPARATEON=0
MERGEONSINGLEData=0
MERGEONSINGLEMC=0
CLEANUP=0
CLEANUPMAYOR=0
number=""

function GetFileNumberListPCM()
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


function GetFileNumberListPCMCalo()
{
    ls $1/GammaConvCalo_*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm fileNumbers.txt
    for fileName in $fileNumbers; do
        number=`echo $fileName  | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $3
    cat $3
}



function GetFileNumberListCalo()
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

function GetFileNumberListCaloMerged()
{
    ls $1/GammaCaloMerged_*.root > filesTemp.txt
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
    if [ -f $1\_A.root ] || [ -f $1\_AWOTree.root ] ; then
        echo "separated file $1.root already"
        if [ $CLEANUP == 1 ]; then
            echo "removing $1.root as already separated"
            rm $1.root
        fi
    else
        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$1.root\"\,\"$1\"\,$2\,kTRUE\,$3\)
        if [ $CLEANUP == 1 ]; then
            if [ -f $1\_A.root ] || [ -f $1\_AWOTree.root ] ; then
                echo "removing $1.root after successful separation"
                rm $1.root
            fi
        fi
    fi
}


function CopyFileIfNonExisitent()
{
    echo "file path: "$1
    echo "alien path: " $2
    echo "number of subdirectories: " $3
    echo "detailed path subdirs: " $4
    echo "separate also the tree: " $5
    if [ $DOWNLOADON == 1 ]; then
        if [ -f $1/root_archive.zip ] && [ -s $1/root_archive.zip ]; then
            echo "$1/root_archive.zip exists";
        else
            mkdir -p $1
            alien_cp alien:$2/root_archive.zip file:$1/
            if [ -f $1/root_archive.zip ] && [ -s $1/root_archive.zip ]; then
                echo "copied correctly"
            else
                rm locallog.txt
                if [ $4 != "none" ]; then
                    if [ $4 == "*Stage*" ]; then
                        alien_ls -F $4/ | grep "/" > locallog.txt
                    else
                        alien_ls -F $4/ | grep "/" | grep -v "Stage" > locallog.txt
                    fi
                    echo "********** log output ***********"
                    cat locallog.txt
                    dirNumbers=`cat locallog.txt`
                    firstDir=`head -n1 locallog.txt`
                    for dirName in $dirNumbers; do
                        mkdir -p $1/Stage1/$dirName
                        if [ -f $1/Stage1/$dirName/root_archive.zip ] && [ -s $1/Stage1/$dirName/root_archive.zip ]; then
                            echo "$1/Stage1/$dirName/root_archive.zip exists";
                        else
                            alien_cp alien:$4/$dirName/root_archive.zip file:$1/Stage1/$dirName/
                            unzip -u $1/Stage1/$dirName/root_archive.zip -d $1/Stage1/$dirName/
                        fi
                    done
                    rm locallog.txt
                    BASEDIR=$PWD
                    cd $1/Stage1/$firstDir/
                    ls G*.root > $BASEDIR/locallog.txt
                    cd -
                    echo "********** log output 2 for merging***********"
                    cat $BASEDIR/locallog.txt
                    fileNumbers=`cat $BASEDIR/locallog.txt`
                    for fileName in $fileNumbers; do
                        hadd -n 10 -f $1/$fileName $1/Stage1/*/$fileName
                    done
                    cd $1/
                    zip root_archive.zip *.root
                    cd -
                fi
            fi
      fi
      unzip -u $1/root_archive.zip -d $1/
    fi

    if [ $SEPARATEON == 1 ]; then
        rm fileNumbers2.txt
        GetFileNumberListPCM $1 $3 fileNumbers2.txt
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaConvV1_$fileNumber 0 $5
        done;
        rm fileNumbers2.txt
        GetFileNumberListPCMCalo $1 $3 fileNumbers2.txt
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaConvCalo_$fileNumber 2 $5
        done;
        rm fileNumbers2.txt
        GetFileNumberListCalo $1 $3 fileNumbers2.txt
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaCalo_$fileNumber 4 $5
        done;
        rm fileNumbers2.txt
    fi
}



function ChangeStructureIfNeededPCM()
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
        cp $2/GammaConvV1$5_$number.root $OUTPUTDIR/GammaConvV1_$4\_$number.root
        if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log ]; then
            echo "nothing to be done";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log\"\,0\)
        fi
   fi
}

function ChangeStructureIfNeededPCMCalo()
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
        cp $2/GammaConvCalo$5_$number.root $OUTPUTDIR/GammaConvCalo_$4\_$number.root
        if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log ]; then
            echo "nothing to be done";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log\"\,2\)
        fi
   fi
}

function ChangeStructureIfNeededCalo()
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
        cp $2/GammaCalo$5_$number.root $OUTPUTDIR/GammaCalo_$4\_$number.root
        if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log ]; then
            echo "nothing to be done";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log\"\,4\)
        fi
   fi
}

function ChangeStructureIfNeededCaloMerged()
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
        cp $2/GammaCaloMerged$5_$number.root $OUTPUTDIR/GammaCaloMerged_$4\_$number.root
        if [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaCaloMerged_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaCaloMerged_$4_$number.log ]; then
            echo "nothing to be done";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaCaloMerged_$4_$number.log\"\,10\)
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

function GetFileNumberMerging2()
{
    echo $1
    NCurrSub=$2
    alpha=`echo $1  | cut -d "/" -f $NCurrSub | cut -d "_" -f3 | cut -d "." -f1`
    number=`echo $1  | cut -d "/" -f $NCurrSub | cut -d "_" -f2 | cut -d "." -f1`
    if [ -z "$alpha" ]; then
        echo $number
    else
        number1=`echo $1  | cut -d "/" -f $NCurrSub | cut -d "_" -f2`
        number=$number1\_$alpha
        echo $number
    fi
    echo $number

}

function MergeAccordingToSpecificRunlist()
{
    echo "filenumbers $1"
    echo "path to subdir $2"
    echo "number of slashes $3"
    echo "general name root file $4"
    echo "name list for output $5"
    echo "runlist $6"

    fileNumbers=`cat $1`
    for fileName in $fileNumbers; do
        echo $fileName
        GetFileNumberMerging2 $fileName $3
        echo $number
        if [[ $number == *"WTree"* ]]; then
            echo "won't merge these";
        else
            TOMERGE="";
            runs=`cat $6`
            for run in $runs; do
                nameCurrFile=`echo $2/$run/$4\_$number.root`
                if [ -f $nameCurrFile ]; then
                    TOMERGE="$TOMERGE $nameCurrFile"
                else
                    echo "I couldn't find the file for bin $run, $nameCurrFile";
                fi
            done
            hadd -n 10 -f $2/$4-$5_$number.root $TOMERGE
        fi

    done;
    echo "done" > $2/mergedAll$4.txt
}

function MergeAccordingToList()
{
    echo "fileList $1"
    echo "final-Output $2"
    fileNamesMerging=`cat $1`
    TOMERGE="";
    counter=0
    for fileName in $fileNamesMerging; do
        if [ -f $fileName ]; then
            TOMERGE="$TOMERGE $fileName"
            counter=$((counter+1))
        else
            echo "I couldn't find the file for bin $fileName";
        fi
    done;
    if [ $counter -gt 1 ]; then
        hadd -n 10 -f $2 $TOMERGE
    else
        cp $TOMERGE $2
    fi
}

