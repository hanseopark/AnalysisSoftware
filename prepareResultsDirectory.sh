if [ $1 = "fbock" ]; then 
   PHOTONCONVDIR=/home/fbock/Photon/Software/PCGGIT/AnalysisSoftware
elif [ $1 = "leardini" ]; then 
   PHOTONCONVDIR=/home/admin1/leardini/photonconv/AnalysisSoftware
elif [ $1 = "dmuhlhei" ]; then
   PHOTONCONVDIR=/home/daniel/data/work/photonconv/AnalysisSoftware
fi

echo $PHOTONCONVDIR

mkdir TaskV1
mkdir TaskQA
mkdir CommonHeaders
mkdir SystematicErrorsNew
mkdir ExternalInput
mkdir ExternalInputPbPb
mkdir ExternalInputpPb
mkdir FinalResults
mkdir CocktailInput
mkdir RooUnfold
mkdir DownloadAndDataPrep
ln -sf $PHOTONCONVDIR/SystematicErrorsNew/* SystematicErrorsNew/
ln -sf $PHOTONCONVDIR/CommonHeaders/*.h CommonHeaders/
ln -sf $PHOTONCONVDIR/TaskV1/*.h TaskV1/
ln -sf $PHOTONCONVDIR/TaskV1/*.C TaskV1/
ln -sf $PHOTONCONVDIR/TaskQA/*.h TaskQA/
ln -sf $PHOTONCONVDIR/TaskQA/*.C TaskQA/
ln -sf $PHOTONCONVDIR/ExternalInput/* ExternalInput/
ln -sf $PHOTONCONVDIR/ExternalInputPbPb/* ExternalInputPbPb/
ln -sf $PHOTONCONVDIR/ExternalInputpPb/* ExternalInputpPb/
ln -sf $PHOTONCONVDIR/FinalResults/* FinalResults/
ln -sf $PHOTONCONVDIR/CocktailInput/* CocktailInput/
ln -sf $PHOTONCONVDIR/RooUnfold/* RooUnfold/
ln -sf $PHOTONCONVDIR/DownloadAndDataPrep/* DownloadAndDataPrep/

ln -sf $PHOTONCONVDIR/*.eps .
ln -sf $PHOTONCONVDIR/*.C .
ln -sf $PHOTONCONVDIR/*.h .
ln -sf $PHOTONCONVDIR/*.sh .

if [ $2 = "pp2760GeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *7TeV*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10*
fi