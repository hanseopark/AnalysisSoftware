if [ $1 = "msas" ]; then 
   PHOTONCONVDIR=/home/mike/git_afterburner/AnalysisSoftware
fi

echo $PHOTONFLOWDIR

mkdir -p CommonHeaders
mkdir -p TaskFlow
mkdir -p TaskFlow/Results
mkdir -p TaskFlow/Software
mkdir -p TaskFlow/Theory
mkdir -p TaskFlow/V2dir_calculation_PCM_PHOS_combined
mkdir -p TaskFlow/V2dir_calculation_PCM_PHOS_combined/output

ln -sf $PHOTONCONVDIR/Produce_FinalPlots_GammaFlowDir.C .
ln -sf $PHOTONCONVDIR/CommonHeaders/*.h CommonHeaders/
ln -sf $PHOTONCONVDIR/TaskFlow/*.txt TaskFlow/
ln -sf $PHOTONCONVDIR/TaskFlow/*.sh TaskFlow/
ln -sf $PHOTONCONVDIR/TaskFlow/Results/* TaskFlow/Results/
ln -sf $PHOTONCONVDIR/TaskFlow/Software/*.C TaskFlow/Software/
ln -sf $PHOTONCONVDIR/TaskFlow/Software/*.h TaskFlow/Software/
ln -sf $PHOTONCONVDIR/TaskFlow/Theory/* TaskFlow/Theory/
ln -sf $PHOTONCONVDIR/TaskFlow/V2dir_calculation_PCM_PHOS_combined/* TaskFlow/V2dir_calculation_PCM_PHOS_combined/
