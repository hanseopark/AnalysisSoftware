if [ $1 = "fbock" ]; then
    PHOTONCONVDIR=/home/fbock/Photon/Software/PCGGIT
elif [ $1 = "fbockNew" ]; then
    PHOTONCONVDIR=/mnt/dae8586b-f761-47b0-b8ef-93a00caf0e11/Photon/Software/PCGGIT
elif [ $1 = "leardini" ]; then
    PHOTONCONVDIR=/home/admin1/leardini/newSoftware/AnalysisSoftware
elif [ $1 = "dmuhlhei" ]; then
    PHOTONCONVDIR=/home/daniel/data/work/pcgGit/AnalysisSoftware
elif [ $1 = "msas" ]; then
    PHOTONCONVDIR=/home/mike/git_afterburner/AnalysisSoftware
elif [ $1 = "jlueh" ]; then
    PHOTONCONVDIR=/home/jens/Programme_und_Skripte/alice-pcg/AnalysisSoftware
elif [ $1 = "jokonig" ]; then
    PHOTONCONVDIR=/home/gustav/AnalysisSoftware/AnalysisSoftware
elif [ $1 = "mdanisch" ]; then
    PHOTONCONVDIR=/home/meike/analysis/software/photonconv/AnalysisSoftware
elif [ $1 = "loizides" ]; then
    PHOTONCONVDIR=/home/loizides/sw/pcm_git
elif [ $1 = "amechler" ]; then
    PHOTONCONVDIR=/home/adrian/git-Framework/AnalysisSoftware
elif [ $1 = "hannahbossi" ]; then
    PHOTONCONVDIR=/Users/hannahbossi/analysis/AnalysisSoftware
elif [ $1 = "ahornung" ]; then
    PHOTONCONVDIR=/Users/andreahornung/pcg-software/AnalysisSoftware
elif [ $1 = "nschmidt" ]; then
    PHOTONCONVDIR=/home/nschmidt/AnalysisSoftware
elif [ $1 = "fjonas" ]; then
    PHOTONCONVDIR=/home/florianjonas/tools/alice/AnalysisSoftware
elif [ $1 = "amarin" ]; then
    PHOTONCONVDIR=/Users/marin/analysis/GIT/AnalysisSoftware
elif [ $1 = "amarinDock" ]; then
    PHOTONCONVDIR=/mnt/myhome/analysis/GIT/AnalysisSoftware
elif [ $1 = "amarinPI" ]; then
    PHOTONCONVDIR=/misc/alidata110/alice_u/marin/GIT/AnalysisSoftware
elif [ $1 = "amechler" ]; then
    PHOTONCONVDIR=/home/adrian/git-Framework/AnalysisSoftware
elif [ $1 = "amechlerAlidock" ]; then
    PHOTONCONVDIR=~/adrian/git-Framework/AnalysisSoftware
elif [ $1 = "redeboer" ]; then
    PHOTONCONVDIR=~/alice/AnalysisSoftware
elif [ $1 = "jlietave" ]; then
    PHOTONCONVDIR=/home/jakub/AnalysisSoftware
elif [ $1 = "llamers" ]; then
    PHOTONCONVDIR=/home/liz/git_afterburner/AnalysisSoftware
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
mkdir -p SupportingMacros
mkdir -p LocalTrainTests

ln -sf $PHOTONCONVDIR/SystematicErrorsNew/* SystematicErrorsNew/
ln -sf $PHOTONCONVDIR/CommonHeaders/*.h CommonHeaders/
ln -sf $PHOTONCONVDIR/TaskV1/*.h TaskV1/
ln -sf $PHOTONCONVDIR/TaskV1/*.C TaskV1/
ln -sf $PHOTONCONVDIR/TaskQA/*.h TaskQA/
ln -sf $PHOTONCONVDIR/TaskQA/*.C TaskQA/
ln -sf $PHOTONCONVDIR/TaskQA/*.sh TaskQA/
ln -sf $PHOTONCONVDIR/ExternalInput/* ExternalInput/
ln -sf $PHOTONCONVDIR/ExternalInputPbPb/* ExternalInputPbPb/
ln -sf $PHOTONCONVDIR/ExternalInputpPb/* ExternalInputpPb/
ln -sf $PHOTONCONVDIR/CocktailInput/* CocktailInput/
ln -sf $PHOTONCONVDIR/RooUnfold/* RooUnfold/
ln -sf $PHOTONCONVDIR/DownloadAndDataPrep/* DownloadAndDataPrep/
ln -sf $PHOTONCONVDIR/ToyModels/* ToyModels/
ln -sf $PHOTONCONVDIR/SimulationStudies/*.C SimulationStudies/
ln -sf $PHOTONCONVDIR/SimulationStudies/*.h SimulationStudies/
ln -sf $PHOTONCONVDIR/SupportingMacros/*.C SupportingMacros/
ln -sf $PHOTONCONVDIR/SupportingMacros/*.h SupportingMacros/
ln -sf $PHOTONCONVDIR/LocalTrainTests/*.C LocalTrainTests/
ln -sf $PHOTONCONVDIR/LocalTrainTests/*.h LocalTrainTests/
ln -sf $PHOTONCONVDIR/LocalTrainTests/*.sh LocalTrainTests/

if [ $1 = "dmuhlhei" ]; then
	ln -sf $PHOTONCONVDIR/DataQA DataQA
fi

if [ $1 = "leardini" ]; then
    mkdir -p LHC11hExternalInputs
	ln -sf $PHOTONCONVDIR/LHC11hExternalInputs/* LHC11hExternalInputs/
fi

if [ $1 = "mdanisch" ]; then
    mkdir -p TaskQA/ExampleConfigurations
    ln -sf $PHOTONCONVDIR/TaskQA/ExampleConfigurations/*LHC15o*.txt TaskQA/ExampleConfigurations/
    ln -sf $PHOTONCONVDIR/0_configurationFile5TeVPbPbPCM*.txt ./
    ln -sf $PHOTONCONVDIR/CompareCentralities*.txt ./
fi

if [ $1 = "fjonas" ]; then
    mkdir -p TaskQA/ExampleConfigurations
    ln -sf $PHOTONCONVDIR/TaskQA/ExampleConfigurations/*.txt TaskQA/ExampleConfigurations/
    ln -sf $PHOTONCONVDIR/DataQA DataQA
fi

if [ $1 = "fjonas" ] || [ $1 = "nschmidt" ] || [ $1 = "jlietave" ] || [ $1 = "redeboer" ]; then
    ln -sf $PHOTONCONVDIR/TaskV1/* TaskV1/
    ln -sf $PHOTONCONVDIR/TaskQA/* TaskQA/
    ln -sf $PHOTONCONVDIR/PublishedResultsAndConfigurations PublishedResultsAndConfigurations
    ln -sf $PHOTONCONVDIR/PreliminaryResultsAndConfigurations PreliminaryResultsAndConfigurations
fi
ln -sf $PHOTONCONVDIR/*.eps .
ln -sf $PHOTONCONVDIR/*.C .
ln -sf $PHOTONCONVDIR/*.h .
ln -sf $PHOTONCONVDIR/*.sh .

if [ "$2" = "pp2760GeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *7TeV*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ "$2" = "pPb5TeV" ]; then
    rm *PbPb*
    rm *pp*
    rm *2760*
    rm *7TeV*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ "$2" = "pp7TeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *2760*
    rm *8TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10h*
    rm *LHC12*
elif [ "$2" = "pp8TeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *2760*
    rm *7TeV*
    rm *13TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ "$2" = "pp13TeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *2760*
    rm *7TeV*
    rm *8TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ "$2" = "pp5TeV" ]; then
    rm *PbPb*
    rm *pPb*
    rm *2760*
    rm *7TeV*
    rm *8TeV*
    rm *LHC11h*
    rm *LHC10*
elif [ "$2" = "PbPb5TeV" ]; then
    rm *PP*
#   rm *pp* need function in PrepareChargedPionDataALICE_pp.C
    rm -r *pPb*
    rm *900G*
#    rm *2760*
    rm *7TeV*
    rm *8TeV*
    rm *13TeV*
#    rm *LHC11h*
elif [ "$2" = "pPb8TeV" ]; then
    rm *PP*
    rm *PbPb*
    rm *pp*
    rm *900G*
    rm *2760*
    rm *7TeV*
    rm *5TeV*
    rm *13TeV*
    rm *LHC11h*
fi
