#local running with this code requires testfiles downloaded and a testSampleESD.txt or testSampleAOD.txt text file with the input files stored in, for example pPb_5TeV/LHC16q/testSampleESD.txt

runMode="MS"
generator="MCGenPythia8TuneX"
numEvents=10000
# not yet used
ecms=5020
tune=14
mkdir -p $generator/$runMode
cd $generator/$runMode
chunk="MCGenPythia8TuneX"
#valgrind \
#  --tool=memcheck \
#  --error-limit=no \
#  --trace-children=yes\
#  --leak-check=full \
#  --max-stackframe=3060888 \
#  --suppressions=$ROOTSYS/etc/valgrind-root.supp \
#  --num-callers=40 \
#  --log-file=/tmp/valgrind_memory.log \
aliroot -x -l -b -q '../../runLocalMCGen.C("'$runmode'",'$numEvents',"'$chunk'",'$ecms','$tune')'


