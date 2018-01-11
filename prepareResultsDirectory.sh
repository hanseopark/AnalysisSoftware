if [ $1 = "fbock" ]; then
   PHOTONCONVDIR=/home/fbock/Photon/Software/PCGGIT
elif [ $1 = "leardini" ]; then
   PHOTONCONVDIR=/home/admin1/leardini/newSoftware/AnalysisSoftware
elif [ $1 = "dmuhlhei" ]; then
   PHOTONCONVDIR=/home/daniel/data/work/pcgGit/AnalysisSoftware
elif [ $1 = "msas" ]; then
   PHOTONCONVDIR=/home/mike/git_afterburner/AnalysisSoftware
elif [ $1 = "jlueh" ]; then
   PHOTONCONVDIR=/home/jens/Cloud/Sciebo/Linux_Arbeitsbereich/Programme_und_Skripte/alice-pcg/AnalysisSoftware
elif [ $1 = "jokonig" ]; then
   PHOTONCONVDIR=PHOTONCONVDIR=/home/joshua/AnalysisSoftware
elif [ $1 = "mdanisch" ]; then
   PHOTONCONVDIR=/home/meike/analysis/software/photonconv/AnalysisSoftware
fi

echo $PHOTONCONVDIR

mkdir -p TaskV1
mkdir -p TaskQA
mkdir -p CommonHeaders
mkdir -p SystematicErrorsNew
mkdir -p ExternalInput
mkdir -p ExternalInputPbPb
mkdir -p ExternalInputpPb
mkdir -p FinalResults
mkdir -p CocktailInput
mkdir -p RooUnfold
mkdir -p DownloadAndDataPrep
mkdir -p ToyModels
mkdir -p SimulationStudies

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
ln -sf $PHOTONCONVDIR/ToyModels/* ToyModels/
ln -sf $PHOTONCONVDIR/SimulationStudies/*.C SimulationStudies/
ln -sf $PHOTONCONVDIR/SimulationStudies/*.h SimulationStudies/

if [ $1 = "dmuhlhei" ]; then
	ln -sf $PHOTONCONVDIR/DataQA DataQA
fi

if [ $1 = "leardini" ]; then
    mkdir -p LHC11hExternalInputs
	ln -sf $PHOTONCONVDIR/LHC11hExternalInputs/* LHC11hExternalInputs/
fi

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
elif [ $2 = "pPb5TeV" ]; then
    rm *PbPb*
    rm *pp*
    rm *2760*
    rm *7TeV*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ $2 = "pp7TeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *2760*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10h*
    rm *LHC12*
elif [ $2 = "pp8TeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *2760*
    rm *7TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ $2 = "pp13TeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *2760*
    rm *7TeV*
    rm *8TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ $2 = "PbPb5TeV" ]; then
    rm *PP*
    rm *pp*
    rm *pPb*
    rm *900G*
    rm *2760*
    rm *7TeV*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
fi
