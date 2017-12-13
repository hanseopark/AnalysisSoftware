#! /bin/bash
echo "Getting List of runNumber Files"
varListOfRunlists=`ls runNumbersLHC*.txt`
echo "Sorting Files"
for runNumberFile in $varListOfRunlists; do
    sort $runNumberFile -o $runNumberFile
done
echo "Done"
