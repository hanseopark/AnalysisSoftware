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
minFileSize=1000
SEPARATEONLYConv=1

tempDir=""
tempPath=""
tempBool=1

function FindCorrectTrainDirectory()
{
  tempDir=""
  tempPath=""
  tempBool=1
  echo $1 $2 $3 $4
  if [ "$1" == "" ]; then
    tempBool=0;
  fi

  testSub=0
  if [ $# -gt 3 ]; then
    if [ "$4" != "" ]; then
      testSub=1;
    fi
  fi

  if [ $tempBool == 1 ]; then
      if [ $testSub == 1 ]; then
#           tempDir=`alien_ls $3 | grep $4\_ | grep $1`
          alien_ls $3 | grep $4\_ > listGrid.txt
          folderNames=`cat listGrid.txt`
          rm listGrid.txt
          for folderName in $folderNames; do
            testing=`echo $folderName  | cut -d "_" -f1`
#             echo $folderName " " $4 " " $testing;
            if [[ $testing == $4 ]]; then
                echo "$folderName" >> listGrid.txt
            fi
          done
          InterMediate=`head -n1 listGrid.txt`\_$1
          InterMediate=${InterMediate/$1\_$1/$1}
          InterMediate=${InterMediate/child_1_$1/$1}
          InterMediateExists="$( cat listGrid.txt | grep -w "$InterMediate" )"
          tempDir="$InterMediate"
          rm listGrid.txt
      else
          tempDir=`alien_ls $3 | grep $1\_`
      fi
      if [ "$tempDir" == "" ]; then
          tempBool=0;
      else
          tempPath="$2-$tempDir"
      fi
  fi
}

function GetFileNumberMerging()
{
    NCurrSub=$3
    number1=`echo $1  | cut -d "/" -f $2 | cut -d "_" -f $NCurrSub | cut -d "." -f1`
    number2=`echo $1  | cut -d "/" -f $2 | cut -d "_" -f $((NCurrSub+1)) | cut -d "." -f1`
    if [ -z "$number2" ]; then
        number=$number1
    else
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

function GetFileNumberListSpecial()
{
    ls $1/$2_*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm -f fileNumbers.txt
    for fileName in $fileNumbers; do
        if [ $4 -eq 0 ]; then
            number=`echo $fileName  | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        else
            GetFileNumberMerging $fileName $2 2
        fi
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $5
    cat $3
}


function GetFileNumberListPCM()
{
    ls $1/GammaConvV1_*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm -f fileNumbers.txt
    for fileName in $fileNumbers; do
        if [ $4 -eq 0 ]; then
            number=`echo $fileName | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        else
            GetFileNumberMerging $fileName $2 2
        fi
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $3
    cat $3
}


function GetFileNumberListPCMCalo()
{
    ls $1/GammaConvCalo*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm -f fileNumbers.txt
    for fileName in $fileNumbers; do
        if [ $4 -eq 0 ]; then
            number=`echo $fileName | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        else
            GetFileNumberMerging $fileName $2 2
        fi
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $3
    cat $3
}

function GetFileNumberListCalo()
{
    ls $1/GammaCalo*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm -f fileNumbers.txt
    for fileName in $fileNumbers; do
        if [ $4 -eq 0 ]; then
            number=`echo $fileName | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        else
            GetFileNumberMerging $fileName $2 2
        fi
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $3
    cat $3
}

function GetFileNumberListHeavy()
{
    ls $1/HeavyNeutralMesonToGG*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm -f fileNumbers.txt
    for fileName in $fileNumbers; do
        number=`echo $fileName  | cut -d "/" -f $2 | cut -d "_" -f 2 | cut -d "." -f1`
        number+=\_`echo $fileName  | cut -d "/" -f $2 | cut -d "_" -f 3 | cut -d "." -f1`
        echo $number >> fileNumbers.txt
    done
    sort -u fileNumbers.txt > $3
    cat $3
}

function GetFileNumberListCaloMerged()
{
    ls $1/GammaCaloMerged*.root > filesTemp.txt
    fileNumbers=`cat filesTemp.txt`
    rm -f fileNumbers.txt
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
            rm -f $1.root
        fi
    else
        root -b -x -q -l SeparateDifferentCutnumbers.C\+\(\"$1.root\"\,\"$1\"\,$2\,kTRUE\,$3\)
        if [ $CLEANUP == 1 ]; then
            if [ -f $1\_A.root ] || [ -f $1\_AWOTree.root ] ; then
                echo "removing $1.root after successful separation"
                rm -f $1.root
            fi
        fi
    fi
}

function CopyRunwiseAndMergeAccordingToRunlistMC()
{
    if [ $2 == 1 ]; then
        echo "downloading $1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbers$1${11}.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
#                 CopyFileIfNonExisitent $3/$runNumber "$7/$1/$runNumber/$5/$4" $8 "$7/$1/$runNumber/$5/$4/Stage_1/" kTRUE
                CopyFileIfNonExisitent $3/$runNumber "$7/$1/$runNumber/$5/$4" $8 "$7/$1/$runNumber/$5/$4/" kTRUE
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $3/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $3/${10}*.root*
                echo runlists/runNumbers$1${11}.txt
                firstrunNumber=`head -n1 runlists/runNumbers$1${11}.txt`
                echo $firstrunNumber
                ls $3/$firstrunNumber/${10}\_*.root > file$1.txt

                listsToMerge=`cat $9`
                for runListName in $listsToMerge; do
                    MergeAccordingToSpecificRunlist file$1.txt $3 $8 ${10} $runListName runlists/runNumbers$1${11}_$runListName.txt "no"
                done
            fi
        else
            CopyFileIfNonExisitent $3 "/alice/cern.ch/user/a/alitrain/PWGGA/$6/$4/merge" $NSlashes "" kTRUE
        fi
    fi
}

function CopyRunwiseAndMergeAccordingToRunlistJJMC()
{
    if [ $2 == 1 ]; then
        echo "downloading $1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbers$1.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                for pThard in {1..20}; do
                    CopyFileIfNonExisitent $3/$pThard/$runNumber "$7/$1/$pThard/$runNumber/$5/$4" $8 "$7/$1/$runNumber/$5/$4/Stage_1/" kTRUE
                done
            done;
            if [ $MERGEONSINGLEMC == 1 ] && [ ! -f $3/mergedAllCalo.txt ]; then
                cd $currentDir
                rm $3/${10}*.root*
                echo runlists/runNumbers$1.txt
                firstrunNumber=`head -n1 runlists/runNumbers$1.txt`
                echo $firstrunNumber
                ls $3/$firstrunNumber/${10}\_*.root > file$1.txt

                listsToMerge=`cat $9`
                for runListName in $listsToMerge; do
                    MergeAccordingToSpecificRunlist file$1.txt $3 $8 ${10} $runListName runlists/runNumbers$1_$runListName.txt "no"
                done
            fi
        else
            CopyFileIfNonExisitent $3 "/alice/cern.ch/user/a/alitrain/PWGGA/$6/$4/merge" $NSlashes "" kTRUE
        fi
    fi
}

function CopyRunwiseAndMergeAccordingToRunlistData()
{
    if [ $2 == 1 ]; then
        echo "downloading $1"
        if [ $SINGLERUN == 1 ]; then
            runNumbers=`cat runlists/runNumbers$1\_${10}.txt`
            echo $runNumbers
            for runNumber in $runNumbers; do
                CopyFileIfNonExisitent $3/$runNumber "$7/$1/000$runNumber/$5/$4" $8 "$7/$1/000$runNumber/$5/$4" kTRUE
            done;
            if [ $MERGEONSINGLEData == 1 ] && [ ! -f $3/mergedAllCalo.txt ]; then
                rm $3/${11}*.root*
                firstrunNumber=`head -n1 runlists/runNumbers$1\_${10}.txt`
                ls $3/$firstrunNumber/${11}\_*.root > file$1.txt
                listsToMerge=`cat $9`
                for runListName in $listsToMerge; do
                    MergeAccordingToSpecificRunlist file$1.txt $3 $8 ${11} $runListName runlists/runNumbers$1\_${10}\_$runListName.txt "no"
                done
            fi
        else
            CopyFileIfNonExisitent $3 "/alice/cern.ch/user/a/alitrain/PWGGA/$6/$4/merge_runlist_1" $NSlashes "" kTRUE
        fi
    fi

}

function CopyFileIfNonExisitent()
{
    echo "Arguments for CopyFileIfNonExisitent:"
    echo "  file path: "$1
    echo "  alien path: " $2
    echo "  number of subdirectories: " $3
    echo "  detailed path subdirs: " $4
    echo "  separate also the tree: " $5
    if [ $DOWNLOADON == 1 ]; then
        if [ -f $1/root_archive.zip ] && [ -s $1/root_archive.zip ]; then
            echo "$1/root_archive.zip exists";
        else
            mkdir -p $1
            alien_cp alien:$2/root_archive.zip file:$1/
            if [ -f $1/root_archive.zip ] && [ -s $1/root_archive.zip ]; then
                echo "copied correctly"
            else
                echo "prolems in this function: CopyFileIfNonExisitent"
                rm -f locallog.txt
                if [ "$4" != "none" ]; then
                    if [ "$4" == "*Stage*" ]; then
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
                    rm -f locallog.txt
                    BASEDIR=$PWD
                    cd $1/Stage1/$firstDir/
                    ls G*.root > $BASEDIR/locallog.txt
                    ls H*.root >> $BASEDIR/locallog.txt
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
        rm -f fileNumbers2.txt
        GetFileNumberListPCM $1 $3 fileNumbers2.txt 0
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/GammaConvV1_$fileNumber 0 $5
        done;
        if [ $SEPARATEONLYConv == 0 ]; then
            rm -f fileNumbers2.txt
            GetFileNumberListPCMCalo $1 $3 fileNumbers2.txt 0
            fileNumbers=`cat fileNumbers2.txt`
            for fileNumber in $fileNumbers; do
                echo $fileNumber
                SeparateCutsIfNeeded $1/GammaConvCalo_$fileNumber 2 $5
            done;
            rm -f fileNumbers2.txt
            # GetFileNumberListCalo $1 $3 fileNumbers2.txt 0
            # fileNumbers=`cat fileNumbers2.txt`
            # for fileNumber in $fileNumbers; do
            #     echo $fileNumber
            #     SeparateCutsIfNeeded $1/GammaCalo_$fileNumber 4 $5
            # done;
            # rm -f fileNumbers2.txt
        fi
    fi
}


function CopyFileIfNonExisitentDiffList()
{
    # Print function arguments
    echo "Arguments for CopyFileIfNonExisitentDiffList:"
    echo "  file path:                \"$1\""
    echo "  alien path:               \"$2\""
    echo "  list:                     \"$3\""
    echo "  number of subdirectories: \"$4\""
    echo "  detailed path subdirs:    \"$5\""
    echo "  separate also the tree:   \"$6\""

    # Download and copy with proper file name
    if [ $DOWNLOADON == 1 ]; then
        if [ -f $1/$3/root_archive.zip ] && [ -s $1/$3/root_archive.zip ]; then
            echo "$1/$3/root_archive.zip exists";
        else
            echo "Couldn't find \"$1/$3/root_archive.zip\". Copying from the GRID..."
            mkdir -p $1/$3
            alien_ls $2/root_archive.zip > templog.txt
            if  grep -q "no such file" templog.txt ;
            then
                echo "no root_archive.zip has been found, copying single files"
                alien_cp alien:$2/Gamma* file:$1/$3/
                cd $1/$3
                zip root_archive.zip *.root
                cd -
            else
                alien_cp alien:$2/root_archive.zip file:$1/$3/
            fi
            if [ -f $1/$3/root_archive.zip ] && [ -s $1/root_archive.zip ]; then
                echo "copied correctly"
            else
                rm -f locallog.txt
                if [ "$5" != "none" ]; then
                    if [ "$5" == "*Stage*" ]; then
                        alien_ls -F $5/ | grep "/" > locallog.txt
                    else
                        alien_ls -F $5/ | grep "/" | grep -v "Stage" > locallog.txt
                    fi
                    echo "********** log output ***********"
                    cat locallog.txt
                    dirNumbers=`cat locallog.txt`
                    firstDir=`head -n1 locallog.txt`
                    for dirName in $dirNumbers; do
                        mkdir -p $1/$3/Stage1/$dirName
                        if [ -f $1/$3/Stage1/$dirName/root_archive.zip ] && [ -s $1/$3/Stage1/$dirName/root_archive.zip ]; then
                            echo "$1/$3/Stage1/$dirName/root_archive.zip exists";
                        else
                            alien_cp alien:$5/$dirName/root_archive.zip file:$1/$3/Stage1/$dirName/
                            unzip -u $1/$3/Stage1/$dirName/root_archive.zip -d $1/$3/Stage1/$dirName/
                        fi
                    done
                    rm -f locallog.txt
                    BASEDIR=$PWD
                    cd $1/$3/Stage1/$firstDir/
                    ls G*.root > $BASEDIR/locallog.txt
                    ls H*.root >> $BASEDIR/locallog.txt
                    cd -
                    echo "********** log output 2 for merging***********"
                    cat $BASEDIR/locallog.txt
                    fileNumbers=`cat $BASEDIR/locallog.txt`
                    for fileName in $fileNumbers; do
                        hadd -n 10 -f $1/$3/$fileName $1/$3/Stage1/*/$fileName
                    done
                    cd $1/$3/
                    zip root_archive.zip *.root
                    cd -
                fi
            fi
        fi
        unzip -u $1/$3/root_archive.zip -d $1/$3/
    fi

    # Seperate if required
    if [ $SEPARATEON == 1 ]; then
        rm -f fileNumbers2.txt
        GetFileNumberListPCM $1/$3 $4 fileNumbers2.txt 0
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/$3/GammaConvV1_$fileNumber 0 $6
        done;
        rm -f fileNumbers2.txt
        GetFileNumberListPCMCalo $1/$3 $4 fileNumbers2.txt 0
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            echo $fileNumber
            SeparateCutsIfNeeded $1/$3/GammaConvCalo_$fileNumber 2 $6
        done;
        # rm -f fileNumbers2.txt
        # GetFileNumberListCalo $1/$3 $4 fileNumbers2.txt 0
        # fileNumbers=`cat fileNumbers2.txt`
        # for fileNumber in $fileNumbers; do
        #     echo $fileNumber
        #     SeparateCutsIfNeeded $1/$3/GammaCalo_$fileNumber 4 $6
        # done;
        rm -f fileNumbers2.txt
        GetFileNumberListHeavy $1/$3 $4 fileNumbers2.txt 0
        fileNumbers=`cat fileNumbers2.txt`
        for fileNumber in $fileNumbers; do
            SeparateCutsIfNeeded $1/$3/HeavyNeutralMesonToGG_$fileNumber 0 $6
        done;
        rm -f fileNumbers2.txt
    fi

    # Copy to main directory
    rm -f fileNumbers2.txt
    GetFileNumberListPCM $1/$3 $4 fileNumbers2.txt 1
    for fileNumber in `cat fileNumbers2.txt`; do
        cp $1/$3/GammaConvV1_$fileNumber.root $1/GammaConvV1-$3\_$fileNumber.root
    done;
    rm -f fileNumbers2.txt
    GetFileNumberListPCMCalo $1/$3 $4 fileNumbers2.txt 1
    for fileNumber in `cat fileNumbers2.txt`; do
        cp $1/$3/GammaConvCalo_$fileNumber.root $1/GammaConvCalo-$3\_$fileNumber.root
    done;
    rm -f fileNumbers2.txt
    GetFileNumberListCalo $1/$3 $4 fileNumbers2.txt 1
    for fileNumber in `cat fileNumbers2.txt`; do
        cp $1/$3/GammaCalo_$fileNumber.root $1/GammaCalo-$3\_$fileNumber.root
    done;
    rm -f fileNumbers2.txt
    GetFileNumberListHeavy $1/$3 $4 fileNumbers2.txt 1
    for fileNumber in `cat fileNumbers2.txt`; do
        cp $1/$3/HeavyNeutralMesonToGG_$fileNumber.root $1/HeavyNeutralMesonToGG-$3\_$fileNumber.root
    done
    rm -f fileNumbers2.txt
}


function ChangeStructureIfNeededPCM()
{
    if [[ $1 == *"Basic.root"* ]]; then
        echo "Nothing to be done: Basis.root"
    else
        GetFileNumberMerging $1 $3 2
        cp $2/GammaConvV1$5_$number.root $OUTPUTDIR/GammaConvV1_$4\_$number.root
        if  [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log ] &&
            [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log ]; then
            echo "Nothing to be done: \"CutSelection_GammaConvV1_$4_$number.log\" already exists";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvV1_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaConvV1_$4_$number.log\"\,0\)
        fi
    fi
}

function ChangeStructureIfNeededPCMCalo()
{
    if [[ $1 == *"Basic.root"* ]]; then
        echo "Nothing to be done: Basis.root"
    else
        GetFileNumberMerging $1 $3 2
        cp $2/GammaConvCalo$5_$number.root $OUTPUTDIR/GammaConvCalo_$4\_$number.root
        if  [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log ] &&
            [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log ]; then
            echo "Nothing to be done: \"CutSelection_GammaConvCalo_$4_$number.log\" already exists";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaConvCalo_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaConvCalo_$4_$number.log\"\,2\)
        fi
    fi
}

function ChangeStructureIfNeededCalo()
{
    if [[ $1 == *"Basic.root"* ]]; then
        echo "Nothing to be done: Basis.root"
    else
        GetFileNumberMerging $1 $3 2
        cp $2/GammaCalo$5_$number.root $OUTPUTDIR/GammaCalo_$4\_$number.root
        if  [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log ] &&
            [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log ]; then
            echo "Nothing to be done: \"CutSelection_GammaCalo_$4_$number.log\" already exists";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCalo_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaCalo_$4_$number.log\"\,4\)
        fi
    fi
}

function ChangeStructureIfNeededCaloMerged()
{
    if [[ $1 == *"Basic.root"* ]]; then
        echo "Nothing to be done: Basis.root"
    else
        GetFileNumberMerging $1 $3 2
        cp $2/GammaCaloMerged$5_$number.root $OUTPUTDIR/GammaCaloMerged_$4\_$number.root
        if  [ -f $OUTPUTDIR/CutSelections/CutSelection_GammaCaloMerged_$4_$number.log ] &&
            [ -s $OUTPUTDIR/CutSelections/CutSelection_GammaCalomerged_$4_$number.log ]; then
            echo "Nothing to be done: \"CutSelection_GammaCaloMerged_$4_$number.log\" already exists";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/GammaCaloMerged_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_GammaCaloMerged_$4_$number.log\"\,10\)
        fi
    fi
}

function ChangeStructureIfNeededHeavy()
{
    if [[ $1 == *"Basic.root"* ]]; then
        echo "Nothing to be done: Basis.root"
    else
        number1=`echo $1  | cut -d "/" -f $3 | cut -d "_" -f 2 | cut -d "." -f1`
        number2=`echo $1  | cut -d "/" -f $3 | cut -d "_" -f 3 | cut -d "." -f1`
        if [ -z "$number2" ]; then
            number=$number1
        else
            # echo $number2
            number=$number1\_$number2
        fi
        # echo $number
        cp $2/HeavyNeutralMesonToGG$5_$number.root $OUTPUTDIR/HeavyNeutralMesonToGG_$4\_$number.root
        if [ -f $OUTPUTDIR/CutSelections/CutSelection_HeavyNeutralMesonToGG_$4_$number.log ] &&  [ -s $OUTPUTDIR/CutSelections/CutSelection_HeavyNeutralMesonToGG_$4_$number.log ]; then
            echo "Nothing to be done (CutSelection_HeavyNeutralMesonToGG_$4_$number)";
        else
            root -b -l -q -x ../TaskV1/MakeCutLog.C\(\"$OUTPUTDIR/HeavyNeutralMesonToGG_$4\_$number.root\"\,\"$OUTPUTDIR/CutSelections/CutSelection_HeavyNeutralMesonToGG_$4_$number.log\"\,100\)
        fi
    fi
}


function MergeAccordingToSpecificRunlist()
{
    echo "filenumbers $1"
    echo "path to subdir $2"
    echo "number of slashes $3"
    echo "general name root file $4"
    echo "name list for output $5"
    echo "runlist $6"
    echo "binlist $7"

    if [ -e $6 ] && [ -f $6 ]; then

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
                if [ "$7" == "no" ]; then
                    for run in $runs; do
                        nameCurrFile=`echo $2/$run/$4\_$number.root`
                        if [ -f $nameCurrFile ]; then
                            TOMERGE="$TOMERGE $nameCurrFile"
                        else
                            echo "I couldn't find the file for bin $run, $nameCurrFile";
                        fi
                    done
                    hadd -n 10 -f $2/$4-$5_$number.root $TOMERGE
                else
                    bins=`cat $7`
                    TOMERGEJJ="";
                    for bin in $bins; do
                        TOMERGE="";
                        nameCurrJJFile=`echo "$2/$bin/$4-$5""_$number.root"`
                        if [ -f $nameCurrJJFile ] && [ $(stat -c %s $nameCurrJJFile) -gt $minFileSize ]; then
                            echo "file $nameCurrJJFile has been created already nothing to be done, will be added to list merge"
                            TOMERGEJJ="$TOMERGEJJ $nameCurrJJFile"
                        else
                            for run in $runs; do
                                nameCurrFile=`echo $2/$bin/$run/$4\_$number.root`
                                if [ -f $nameCurrFile ]; then
                                    TOMERGE="$TOMERGE $nameCurrFile"
                                else
                                    echo "I couldn't find the file for bin $bin run $run, $nameCurrFile";
                                fi
                            done
                            hadd -n 10 -f $nameCurrJJFile $TOMERGE
                        fi
                    done
                    hadd -n 10 -f "$2/$4-$5""_$number.root" $TOMERGEJJ

                    for run in $runs; do
                        nameCurrJJFileRun=`echo "$2/$run/$4""_$number.root"`
                        if [ -f $nameCurrJJFileRun ] && [ $(stat -c %s $nameCurrJJFileRun) -gt $minFileSize ]; then
                            echo "file $nameCurrJJFileRun has been created already nothing to be done"
                        else
                            TOMERGERun="";
                            for bin in $bins; do
                                nameCurrFile=`  echo $2/$bin/$run/$4\_$number.root`
                                if [ -f $nameCurrFile ]; then
                                    TOMERGERun="$TOMERGERun $nameCurrFile"
                                else
                                    echo "I couldn't find the file for bin $bin run $run, $nameCurrFile";
                                fi
                            done
                            mkdir -p $2/$run/
                            hadd -n 10 -f $nameCurrJJFileRun $TOMERGERun
                        fi
                    done

                fi
            fi
        done;
        echo "done" > $2/mergedAll$4.txt
    else
        echo "runlist: $6 not found"
    fi
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
    elif [ $counter -eq 1 ]; then
        cp $TOMERGE $2
    fi
}

# MergeFilesAccrosYears merges all files in a folder that have the same period and the same train config number
# NOTE $1: input folder
# NOTE $2: file name in that folder that will be used to find all runlists and all train config numbers
function MergeFilesAccrossYears()
{

    # * Import function arguments *
        cd "$1"
        local filenameToUse="$2" # e.g."HeavyNeutralMesonToGG_LHC17e-pass1-DPGTracks_1_400.root

    # * Variable names *
        local currentDir="$(pwd)"
        local outputFolder="merged"
        if [ ! -f "$filenameToUse" ]; then
            echo "ERROR: File \"$filenameToUse\" does not exist!"
            echo "--> Check this script or check folder \"$1\""
            exit
        fi
        local analysisType=${filenameToUse/_LHC1*-pass?-*_*.root}
        local periodToCheck=${filenameToUse/${analysisType}_}
        periodToCheck=${periodToCheck/-pass?-*_*.root}
        local passNrToCheck=${filenameToUse/${analysisType}_${periodToCheck}-pass}
        passNrToCheck=${passNrToCheck/-*_*.root}
        local runlistToCheck=${filenameToUse/${analysisType}_${periodToCheck}-pass${passNrToCheck}-}
        runlistToCheck=${runlistToCheck/_*.root}
        local trainconfigToCheck=${filenameToUse/${analysisType}_${periodToCheck}-pass${passNrToCheck}-${runlistToCheck}_}
        trainconfigToCheck=${trainconfigToCheck/.root}
        echo -e "Merged files will be stored in:"
        echo -e "  \"`pwd`/$outputFolder\"\n"
        echo -e "The following test file name is used for merging:"
        echo -e "  \"$filenameToUse\""
        echo -e "  (modify this one in the script if needed)\n"
        echo -e "Derived parameters:"
        echo -e " - Analysis type:  \"$analysisType\""
        echo -e " - Period name:    \"$periodToCheck\""
        echo -e " - Pass number:    \"$passNrToCheck\""
        echo -e " - Train config:   \"$trainconfigToCheck\""
        echo -e " - Runlist:        \"$runlistToCheck\""
    #

    # * Create list of runlists to be merged *
        rm -f temp_RunLists.txt
        local beginString="${analysisType}_${periodToCheck}-pass?-"
        local endString="_${trainconfigToCheck}.root"
        local newString=""
        for i in `ls $beginString*$endString`; do
            newString=$i
            newString=${newString/$endString}
            newString=${newString/$beginString}
            echo "$newString" >> temp_RunLists.txt
        done
        runLists=`cat temp_RunLists.txt`
        nRunlists=`cat temp_RunLists.txt | wc -l`
        echo -e "\nFound the following $nRunlists runlists for \"$beginString*$endString\":"
        echo -e "$runLists"
    #

    # * Create list of train configs to be merged *
        rm -f temp_TrainConfigs.txt
        beginString="${analysisType}_${periodToCheck}-pass?-${runlistToCheck}_"
        endString=".root"
        for i in `ls $beginString*$endString`; do
            newString=$i
            newString=${newString/$endString}
            newString=${newString/$beginString}
            echo "$newString" >> temp_TrainConfigs.txt
        done
        trainconfigList=`cat temp_TrainConfigs.txt`
        nTrainConfigs=`cat temp_TrainConfigs.txt | wc -l`
        echo -e "\nFound the following $nTrainConfigs train configs for \"$beginString*$endString\":"
        echo -e "$trainconfigList"
        echo ""
        read -p "Press enter if this is ok..."
        echo ""
    #

    # * Make new output directory for merged files *
        if [ ! -d $outputFolder ]; then
            mkdir $outputFolder
        fi

    # * Get period names from files *
        beginString="${analysisType}_LHC"
        local fileNames
        local mergePeriods
        local outputName
        for runList in $runLists; do
            for trainConfig in $trainconfigList; do
                echo -e "\n\n---=== Adding runlist \"${runList}\", train config \"${trainConfig}\" ===---"
                endString="-pass?-${runList}_${trainConfig}.root"
                fileNames=`ls ${beginString}*${endString}`
                mergePeriods=""
                echo "" > temp_YearList.txt # clean file
                for fileName in $fileNames; do
                    # * Get period name from fileName
                    newString=${fileName/$beginString}
                    newString=${newString/$endString}
                    # * Remove year
                    for year in `cat temp_YearList.txt`; do
                        newString=${newString/$year}
                    done
                    # * Add year to temp_YearList.txt if not yet already there
                    if [[ ${newString:0:2} =~ ^[0-9]+$ ]]; then
                        echo "${newString:0:2}" >> temp_YearList.txt
                    fi
                    # * Add period name to mergePeriods
                    if [ -n newString ]; then
                        mergePeriods+=$newString
                    fi
                done
                outputName=$beginString$mergePeriods$endString
                outputName=${outputName/-pass?} # remove pass number
                # echo "Output merged file name: \"$outputName\""
                # * Merge ROOT files
                hadd -f $outputFolder/$outputName $fileNames
            done
        done

    # * Remove temporary list files *
        rm -f temp_RunLists.txt
        rm -f temp_TrainConfigs.txt
        rm -f temp_YearList.txt

    # * Go back to former working directory
        cd "${currentDir}"
}
