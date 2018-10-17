Local traintests can be run with the code provided in this folder.
Follow these steps to use it:
1. Download testfiles AOD or ESD with the download script (bash downloadTestFiles.sh).
    Set in the script via dataOrMC if you want data or MC
    Set in the script via ESDOrAOD if you wand ESDs or AODs
    Set energy, yearData/MC, periodData/MC, the AODFILTER and RUN to your needs
2. Once you have downloaded files, create a text file containing the full path to the input files inside the folder energy/periodname, e.g: pPb_5TeV/LHC16q/testSampleESD.txt or pPb_5TeV/LHC16q/testSampleAOD.txt
3. Adjust in runLocalAnalysisROOT6.C the wagons you want to run.
4. Run the local traintest via startANAROOT6.sh.
    Set according to your downloaded files energy, runPeriod and runPeriodData, dataType (ESD or AOD), set recoPassData and if you want to use the correction framework via useCorrTask
