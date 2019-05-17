# NOTE: Macro parameters
    GPPcompile=0
    fileNamePCM="pdf/13TeV/2018_09_16/MC_PCMResultsFullCorrection_pPb.root"
    fileNamePCMEMCAL="pdf/13TeV/2018_09_16/MC_PCM-EMCALResultsFullCorrection_pPb.root"
    fileNamePCMPHOS="pdf/13TeV/2018_09_16/MC_PCM-PHOSResultsFullCorrection_pPb.root"
    fileNameEMCAL="pdf/13TeV/2018_09_16/MC_EMCAL-EMCALResultsFullCorrection_pPb.root"
    fileNamePHOS="pdf/13TeV/2018_09_16/MC_PHOS-PHOSResultsFullCorrection_pPb.root"
    suffix="pdf"
    isMC=kTRUE
    bWCorrection=""
    fileNameCorrFactors=""
#

# * Compile with ROOT :(( *
    if [ $GPPcompile -eq 0 ]; then
        root -l -x -q -b "CombineEtaPrimeMeasurements5TeV_pPb.C+(\"$fileNamePCM\",\"$fileNamePCMEMCAL\",\"$fileNamePCMPHOS\",\"$fileNameEMCAL\",\"$fileNamePHOS\",\"$suffix\",$isMC,\"$bWCorrection\",\"$fileNameCorrFactors\")"
# * Compile with G++ \('.')/ *
    else
        echo -n "Compiling \"CombineEtaPrimeMeasurements5TeV_pPb\" ..."
        compile CombineEtaPrimeMeasurements5TeV_pPb.C
        if [ $? -eq 0 ]; then
            echo -e "\rRunning \"CombineEtaPrimeMeasurements5TeV_pPb\" ...\n"
            ./CombineEtaPrimeMeasurements5TeV_pPb $fileNamePCM $fileNamePCMEMCAL $fileNamePCMPHOS $fileNameEMCAL $fileNamePHOS $suffix $isMC $bWCorrection $fileNameCorrFactors
        fi
    fi
#