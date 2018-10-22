#local running with this code requires testfiles downloaded and a testSampleESD.txt or testSampleAOD.txt text file with the input files stored in, for example pPb_5TeV/LHC16q/testSampleESD.txt

#energy="pPb_5TeV"
energy="pp_5TeV"
intMCrunning=0 #0: data, 1: MC, 2: JJ MC
#collsys=2 #0: pp, 1: PbPb, 2: pPb
collsys=0 #0: pp, 1: PbPb, 2: pPb
runPeriod="LHC17p"
#runPeriod="LHC18f3_2"
#runPeriodData="LHC16q"
runPeriodData="LHC17p"
dataType=ESD #ESD or AOD
runMode="P" #switch for which tasks to run: QA (photon and cluster QA), P (PCM), C (Calo [EMC, DMC, PHOS]), H (hybrid PCM-Calo), M (merged EMC), S (skimming ESD or AOD)
recoPassData=1
tenderPassData="pass1"
useCorrTask="kTRUE"
aodConversionCutnumber="80000003_06000008400000001000000000";
numLocalFiles=5

mkdir -p $energy/$runPeriod/$runMode$dataType
cd $energy/$runPeriod/$runMode$dataType
aliroot -x -l -b -q '../../../runLocalAnalysisROOT6.C('$intMCrunning','$collsys', "'$runPeriod'", "'$runPeriodData'", "'$dataType'", "'$runMode'", '$recoPassData', "'$tenderPassData'", '$useCorrTask', "'$aodConversionCutnumber'", '$isRun2', '$numLocalFiles')'
