if [ $1 = "fbock" ]; then 
   PHOTONCONVDIR=/home/fbock/Photon/Software/photonconv/AnalysisSoftware
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

ln -sf $PHOTONCONVDIR/*.eps .
ln -sf $PHOTONCONVDIR/*.root .
ln -sf $PHOTONCONVDIR/*.C .
ln -sf $PHOTONCONVDIR/*.dat .
ln -sf $PHOTONCONVDIR/*.h .
ln -sf $PHOTONCONVDIR/*.sh .
ln -sf $PHOTONCONVDIR/*.pl .
ln -sf $PHOTONCONVDIR/binNumbers*.txt .
ln -sf $PHOTONCONVDIR/runNumbers*.txt .

if [ $2 = "pp2760GeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *7TeV*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10*
    rm *LHC13b*
    rm *LHC13e7*
    rm *LHC14*
fi