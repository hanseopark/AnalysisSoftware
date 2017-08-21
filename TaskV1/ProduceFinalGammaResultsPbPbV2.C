#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
//#include "AliHEPDataParser.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "CalculateGammaToPi0V3.h"

void ProduceFinalGammaResultsPbPbV2(TString cutSel        = "",
                                    TString optionEnergy  = "PbPb_2.76TeV",
                                    TString suffix        = "eps",
                                    Int_t mode            = 0)
{

    StyleSettingsThesis();
    SetPlotStyle();

    TString dateForOutput           = ReturnDateStringForOutput();
    TString collisionSystem         = ReturnFullCollisionsSystem(optionEnergy);
    TString collisionSystemOutput   = ReturnCollisionEnergyOutputString(optionEnergy);
    TString centrality              = GetCentralityString(cutSel);
    TString centralityW0Per         = GetCentralityStringWoPer(cutSel);
    TString detectionProcess        = ReturnFullTextReconstructionProcess(mode);

    Color_t colorNLOcalc          = kBlue-2;
    Color_t colorCent             = GetColorDefaultColor(optionEnergy, "", GetCentralityString(cutSel));
    Color_t colorCentNotPi0Fitted = GetColorDefaultColor(optionEnergy, "", GetCentralityString(cutSel)) -7;
    Color_t colorErrA             = kGreen-8;
    Color_t colorErrB             = kRed-8;
    Color_t colorErrC             = kBlue-8;

    //******************************************************
    Int_t nLinesNLOLegends  = 2;
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0)
        nLinesNLOLegends    = 3;
    Double_t minPt   = 0.6;
    Double_t maxPt   = 40.;
    Double_t minYDR  = 0.75;
    Double_t maxYDR  = 2.;
    Double_t minYIR  = 0.;
    Double_t maxYIR  = 2.;

    Int_t binOffset = 4;

    //******************************************************* defining output directory *******************************************************
    TString outputDir            = Form("%s/%s/FinalGammaResults_%s_%s", suffix.Data(), dateForOutput.Data(), centralityW0Per.Data(), collisionSystemOutput.Data() );
    gSystem->Exec("mkdir -p "+outputDir);

    TString inputFileName        = Form("%s/%s/Gamma_Pi0_data_GammaConvV1_InclusiveRatio.root", cutSel.Data(), optionEnergy.Data());
    cout << "trying to read: " << inputFileName.Data() << endl;
    TFile *fileInput             = new TFile(inputFileName.Data());
    if (fileInput->IsZombie()) {
        cout << "file couldn't be read, aborting....";
        return;
    }

    // read stat error hists
    TH1D* histoIncGamma                     = (TH1D*) fileInput->Get("histoGammaSpecCorrPurity");
    TH1D* histoPi0Spectrum                  = (TH1D*) fileInput->Get("CorrectedYieldTrueEff");
    TH1D* histoIncRatio                     = (TH1D*) fileInput->Get("IncRatioPurity_trueEff");
    TH1D* histoIncRatioPi0Fit               = (TH1D*) fileInput->Get("histoIncRatioFitPurity");
    TH1D* histoDR                           = (TH1D*) fileInput->Get("DoubleRatioTrueEffPurity");
    TH1D* histoDRFit                        = (TH1D*) fileInput->Get("DoubleRatioFitPurity");
    TH1D* histoDirGamma                     = (TH1D*) fileInput->Get("histoDirectPhotonSpectrum");
    // read sys error graphs
//     TGraphAsymmErrors* graphIncGammaSysErr      = (TGraphAsymmErrors*) fileInput->Get("histoGammaSpecCorrPurity_SystErr");
//     TGraphAsymmErrors* graphIncRatioSysErr      = (TGraphAsymmErrors*) fileInput->Get("IncRatioPurity_trueEff_SystErr");
//     TGraphAsymmErrors* graphIncRatioPi0FitSysErr= (TGraphAsymmErrors*) fileInput->Get("histoIncRatioFitPurity_SystErr");
//     TGraphAsymmErrors* graphDRSysErr            = (TGraphAsymmErrors*) fileInput->Get("DoubleRatioTrueEffPurity_SystErr");
//     TGraphAsymmErrors* graphDRPi0FitSysErr      = (TGraphAsymmErrors*) fileInput->Get("DoubleRatioFitPurity_SystErr");
    // read theory graphs
    TGraphAsymmErrors* graphNLODR               = (TGraphAsymmErrors*) fileInput->Get("graphNLODoubleRatio");
    TGraphAsymmErrors* graphNLOGammaDir         = (TGraphAsymmErrors*) fileInput->Get("graphNLODirGamma");
    TGraphAsymmErrors* graphNLOGammaPrompt      = (TGraphAsymmErrors*) fileInput->Get("graphPromptPhotonNLO");
    TGraphAsymmErrors* graphNLOGammaFrag        = (TGraphAsymmErrors*) fileInput->Get("graphFragmentationPhotonNLO");

    TH1D *histoPi0SpectrumFit                   = (TH1D*) fileInput->Get("fitPi0YieldC");
    TH1D *histoPi0SpectrumFitA                  = (TH1D*) fileInput->Get("fitPi0YieldA");
    TH1D *histoPi0SpectrumFitB                  = (TH1D*) fileInput->Get("fitPi0YieldB");
    TH1D *histoPi0SpectrumFitC                  = (TH1D*) fileInput->Get("fitPi0YieldC");

    TH1D *histoCocktailAllGamma                 = (TH1D*) fileInput->Get("Gamma_Pt");
    Double_t paramsCocktail[6];
    GetFitParameter("qcd",GetCentralityString(cutSel),paramsCocktail);
    cout << paramsCocktail[0] << "\t" << paramsCocktail[1] << endl;
    TF1 *cocktailFitAllGammaForNLO                                 = (TF1*)  FitObject("qcd","cocktailFit","Pi0",histoCocktailAllGamma,2.0,11,paramsCocktail,"QNRME+");
    cocktailFitAllGammaForNLO->SetRange(0,20);
    cout << WriteParameterToFile(cocktailFitAllGammaForNLO)<< endl;

    Double_t fNcoll    = 0;
    Double_t fNcollErr = 0;
        fNcoll                             = GetNCollFromCutNumber(cutSel,optionEnergy);
        fNcollErr                          = GetNCollErrFromCutNumber(cutSel,optionEnergy);
    cout << "Ncoll = " << fNcoll << " +/- " << fNcollErr << endl;

    TString PubCombGammaFileName = "ExternalInputPbPb/CombDirGamma/Gamma_CombResults_PbPb_2.76TeV_20150729_Pub2015.root";
    TString PubPCMGammaFileName  = "ExternalInputPbPb/PCM/Gamma_PCMResults_PbPb_2.76TeV_20150729_Pub2015.root";
    TFile *filePubGammaPCM       = new TFile(PubPCMGammaFileName.Data());
    TFile *filePubGammaComb      = new TFile(PubCombGammaFileName.Data());
    if (filePubGammaPCM->IsZombie() || filePubGammaComb->IsZombie()) {
        cout << "published gamma file couldn't be read, aborting....";
        return;
    }

    TDirectoryFile* directoryCombGamma_0020                         = (TDirectoryFile*)filePubGammaComb->Get("Gamma_PbPb_2.76TeV_0-20%");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumStat_0020    = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DirGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumSyst_0020    = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DirGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumStat_0020   = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("IncGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumSyst_0020   = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("IncGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioStat_0020         = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DR_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioSyst_0020         = (TGraphAsymmErrors*)directoryCombGamma_0020->Get("DR_comb_SysErr");

    TDirectoryFile* directoryCombGamma_2040                         = (TDirectoryFile*)filePubGammaComb->Get("Gamma_PbPb_2.76TeV_20-40%");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumStat_2040    = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DirGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDirGammaSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DirGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumStat_2040   = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("IncGammaSpec_comb_StatErr");
        TGraphAsymmErrors* graphPubCombInclGammaSpectrumSyst_2040   = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("IncGammaSpec_comb_SysErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioStat_2040         = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DR_comb_StatErr");
        TGraphAsymmErrors* graphPubCombDoubleRatioSyst_2040         = (TGraphAsymmErrors*)directoryCombGamma_2040->Get("DR_comb_SysErr");

    TDirectoryFile* directoryPCMGamma_0020                          = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_0-20%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_0020     = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_0020     = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_0020                 = (TH1D*)directoryPCMGamma_0020->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_0020    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_0020);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_0020    = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_0020                 = (TH1D*)directoryPCMGamma_0020->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_0020    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_0020);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_0020    = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("IncRatioSystError");
        TH1D *histoPubPCMDoubleRatioStat_0020                       = (TH1D*)directoryPCMGamma_0020->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_0020          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_0020);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_0020          = (TGraphAsymmErrors*)directoryPCMGamma_0020->Get("DoubleRatioSystError");

    TDirectoryFile* directoryPCMGamma_0010                         = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_0-10%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_0010   = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_0010                 = (TH1D*)directoryPCMGamma_0010->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_0010    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_0010);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_0010                 = (TH1D*)directoryPCMGamma_0010->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_0010    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_0010);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_0010    = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("IncRatioSystError");
        TH1D *histoPubPCMDoubleRatioStat_0010                       = (TH1D*)directoryPCMGamma_0010->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_0010          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_0010);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_0010          = (TGraphAsymmErrors*)directoryPCMGamma_0010->Get("DoubleRatioSystError");

    TDirectoryFile* directoryPCMGamma_2040                         = (TDirectoryFile*)filePubGammaPCM->Get("Gamma_PbPb_2.76TeV_20-40%");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumStat_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("graphDirGammaSpectrumStat");
        TGraphAsymmErrors* graphPubPCMDirGammaSpectrumSyst_2040   = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("graphDirGammaSpectrumSyst");
        TH1D *histoPubPCMInclGammaSpectrumStat_2040                 = (TH1D*)directoryPCMGamma_2040->Get("IncGammaStatError");
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumStat_2040    = new TGraphAsymmErrors(histoPubPCMInclGammaSpectrumStat_2040);
        TGraphAsymmErrors* graphPubPCMInclGammaSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("IncGammaSystError");
        TH1D *histoPubPCMInclRatioSpectrumStat_2040                 = (TH1D*)directoryPCMGamma_2040->Get("IncRatioStatError");
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumStat_2040    = new TGraphAsymmErrors(histoPubPCMInclRatioSpectrumStat_2040);
        TGraphAsymmErrors* graphPubPCMInclRatioSpectrumSyst_2040    = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("IncRatioSystError");
        TH1D *histoPubPCMDoubleRatioStat_2040                       = (TH1D*)directoryPCMGamma_2040->Get("DoubleRatioStatError");
        TGraphAsymmErrors* graphPubPCMDoubleRatioStat_2040          = new TGraphAsymmErrors(histoPubPCMDoubleRatioStat_2040);
        TGraphAsymmErrors* graphPubPCMDoubleRatioSyst_2040          = (TGraphAsymmErrors*)directoryPCMGamma_2040->Get("DoubleRatioSystError");

    TGraphAsymmErrors* graphPublishedDirGammaSpectrumStat            = NULL;
    TGraphAsymmErrors* graphPublishedDirGammaSpectrumSyst           = NULL;
    TGraphAsymmErrors* graphPublishedInclGammaSpectrumStat            = NULL;
    TGraphAsymmErrors* graphPublishedInclGammaSpectrumSyst           = NULL;
    TGraphAsymmErrors* graphPublishedInclRatioSpectrumStat           = NULL;
    TGraphAsymmErrors* graphPublishedInclRatioSpectrumSyst          = NULL;
    TGraphAsymmErrors* graphPublishedDoubleRatioStat                 = NULL;
    TGraphAsymmErrors* graphPublishedDoubleRatioSyst                 = NULL;
    if(centrality.CompareTo("0-10%")==0){
//         graphPublishedDirGammaSpectrumStat                           = (TGraphAsymmErrors*)graphPubPCMDirGammaSpectrumStat_0020->Clone("graphPublishedDirGammaSpectrumStat");
//         graphPublishedDirGammaSpectrumSyst                          = (TGraphAsymmErrors*)graphPubPCMDirGammaSpectrumSyst_0020->Clone("graphPublishedDirGammaSpectrumSyst");
//         graphPublishedInclGammaSpectrumStat                          = (TGraphAsymmErrors*)graphPubPCMInclGammaSpectrumStat_0020->Clone("graphPublishedInclGammaSpectrumStat");
//         graphPublishedInclGammaSpectrumSyst                         = (TGraphAsymmErrors*)graphPubPCMInclGammaSpectrumSyst_0020->Clone("graphPublishedInclGammaSpectrumSyst");
//         graphPublishedDoubleRatioStat                                = (TGraphAsymmErrors*)graphPubPCMDoubleRatioStat_0020->Clone("graphPublishedDoubleRatioStat");
//         graphPublishedDoubleRatioSyst                                = (TGraphAsymmErrors*)graphPubPCMDoubleRatioSyst_0020->Clone("graphPublishedDoubleRatioSyst");
        graphPublishedDirGammaSpectrumStat                           = (TGraphAsymmErrors*)graphPubPCMDirGammaSpectrumStat_0010->Clone("graphPublishedDirGammaSpectrumStat");
        graphPublishedDirGammaSpectrumSyst                          = (TGraphAsymmErrors*)graphPubPCMDirGammaSpectrumSyst_0010->Clone("graphPublishedDirGammaSpectrumSyst");
        graphPublishedInclGammaSpectrumStat                           = (TGraphAsymmErrors*)graphPubPCMInclGammaSpectrumStat_0010->Clone("graphPublishedInclGammaSpectrumStat");
        graphPublishedInclGammaSpectrumSyst                          = (TGraphAsymmErrors*)graphPubPCMInclGammaSpectrumSyst_0010->Clone("graphPublishedInclGammaSpectrumSyst");
        graphPublishedInclRatioSpectrumStat                          = (TGraphAsymmErrors*)graphPubPCMInclRatioSpectrumStat_0010->Clone("graphPublishedInclRatioSpectrumStat");
        graphPublishedInclRatioSpectrumSyst                         = (TGraphAsymmErrors*)graphPubPCMInclRatioSpectrumSyst_0010->Clone("graphPublishedInclRatioSpectrumSyst");
        graphPublishedDoubleRatioStat                                = (TGraphAsymmErrors*)graphPubPCMDoubleRatioStat_0010->Clone("graphPublishedDoubleRatioStat");
        graphPublishedDoubleRatioSyst                                = (TGraphAsymmErrors*)graphPubPCMDoubleRatioSyst_0010->Clone("graphPublishedDoubleRatioSyst");
    } else {
        graphPublishedDirGammaSpectrumStat                           = (TGraphAsymmErrors*)graphPubPCMDirGammaSpectrumStat_2040->Clone("graphPublishedDirGammaSpectrumStat");
        graphPublishedDirGammaSpectrumSyst                          = (TGraphAsymmErrors*)graphPubPCMDirGammaSpectrumSyst_2040->Clone("graphPublishedDirGammaSpectrumSyst");
        graphPublishedInclGammaSpectrumSyst                          = (TGraphAsymmErrors*)graphPubPCMInclGammaSpectrumSyst_2040->Clone("graphPublishedInclGammaSpectrumStat");
        graphPublishedInclGammaSpectrumStat                          = (TGraphAsymmErrors*)graphPubPCMInclGammaSpectrumStat_2040->Clone("graphPublishedInclGammaSpectrumSyst");
        graphPublishedInclRatioSpectrumStat                          = (TGraphAsymmErrors*)graphPubPCMInclRatioSpectrumStat_2040->Clone("graphPublishedInclRatioSpectrumStat");
        graphPublishedInclRatioSpectrumSyst                         = (TGraphAsymmErrors*)graphPubPCMInclRatioSpectrumSyst_2040->Clone("graphPublishedInclRatioSpectrumSyst");
        graphPublishedDoubleRatioStat                                = (TGraphAsymmErrors*)graphPubPCMDoubleRatioStat_2040->Clone("graphPublishedDoubleRatioStat");
        graphPublishedDoubleRatioSyst                                = (TGraphAsymmErrors*)graphPubPCMDoubleRatioSyst_2040->Clone("graphPublishedDoubleRatioSyst");
    }
    DrawGammaSetMarkerTGraphAsym(graphPublishedDirGammaSpectrumStat , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPublishedDirGammaSpectrumSyst , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPublishedInclGammaSpectrumStat , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPublishedInclGammaSpectrumSyst , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPublishedInclRatioSpectrumStat , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPublishedInclRatioSpectrumSyst , 20, 2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPublishedDoubleRatioStat , 20,2, kGray+1, kGray+1, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPublishedDoubleRatioSyst , 20, 2, kGray+1, kGray+1, 1, kTRUE);


    TString fileNameSysErrDoubleRatio       = Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveragedSepErrType_DoubleRatio_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    TString fileNameSysErrDoubleRatioPi0Fit = Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveragedSepErrType_DoubleRatio_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    //= Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveraged_DoubleRatioPi0Fit_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    TString fileNameSysErrIncGamma          = Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveragedSepErrType_Gamma_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    TString fileNameSysErrIncRatio          = Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveragedSepErrType_IncRatio_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    TString fileNameSysErrIncRatioPi0Fit    = Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveragedSepErrType_IncRatio_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    //= Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveraged_IncRatioPi0Fit_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    TString fileNameSysErrPi0               = Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveragedSepErrType_Pi0_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());
    TString fileNameSysErrPi0Fit            = Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveragedSepErrType_Pi0_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());//= Form("GammaSystematicErrorsCalculated_2017_08_11/SystematicErrorAveraged_Pi0Fit_PbPb2760GeV%s_2017_08_11.dat",centralityW0Per.Data());

    // ******************************************************************
    // *********** reading systematic errors for double ratio ***********
    // ******************************************************************
    ifstream         fileSysErrDoubleRatio;
    fileSysErrDoubleRatio.open(fileNameSysErrDoubleRatio,ios_base::in);
    Double_t         relSystErrorDoubleRatioUp[50];
    Double_t         relSystErrorDoubleRatioDown[50];
    Double_t         relSystErrorWOMaterialDoubleRatioUp[50];
    Double_t         relSystErrorWOMaterialDoubleRatioDown[50];
    Double_t         relSystErrorADoubleRatioUp[50];
    Double_t         relSystErrorADoubleRatioDown[50];
    Double_t         relSystErrorBDoubleRatioUp[50];
    Double_t         relSystErrorBDoubleRatioDown[50];
    Double_t         relSystErrorCDoubleRatioUp[50];
    Double_t         relSystErrorCDoubleRatioDown[50];

    cout << fileNameSysErrDoubleRatio << endl;
    Int_t nPointsErrors=0;
    while(!fileSysErrDoubleRatio.eof() && nPointsErrors < 100){
        fileSysErrDoubleRatio   >> ptSysDoubleRatio[nPointsErrors]
                                >> relSystErrorDoubleRatioDown[nPointsErrors]                 >> relSystErrorDoubleRatioUp[nPointsErrors]
                                >> relSystErrorWOMaterialDoubleRatioDown[nPointsErrors]       >> relSystErrorWOMaterialDoubleRatioUp[nPointsErrors]
                                >> relSystErrorADoubleRatioDown[nPointsErrors]                >> relSystErrorADoubleRatioUp[nPointsErrors]
                                >> relSystErrorBDoubleRatioDown[nPointsErrors]                >> relSystErrorBDoubleRatioUp[nPointsErrors]
                                >> relSystErrorCDoubleRatioDown[nPointsErrors]                >> relSystErrorCDoubleRatioUp[nPointsErrors];
        cout << nPointsErrors <<  "\t"  << ptSysDoubleRatio[nPointsErrors] << "\t"  << relSystErrorDoubleRatioDown[nPointsErrors] << "\t" << relSystErrorDoubleRatioUp[nPointsErrors] << "\t"
             << relSystErrorWOMaterialDoubleRatioDown[nPointsErrors] << "\t"  <<relSystErrorWOMaterialDoubleRatioUp[nPointsErrors] << endl;;
        nPointsErrors++;
    }
    fileSysErrDoubleRatio.close();

    // ******************************************************************
    // ****** reading systematic errors for double ratio pi0 fit ********
    // ******************************************************************
    ifstream         fileSysErrDoubleRatioPi0Fit;
    fileSysErrDoubleRatioPi0Fit.open(fileNameSysErrDoubleRatioPi0Fit, ios_base::in);
    Double_t         relSystErrorDoubleRatioPi0FitUp[50];
    Double_t         relSystErrorDoubleRatioPi0FitDown[50];
    Double_t         relSystErrorWOMaterialDoubleRatioPi0FitUp[50];
    Double_t         relSystErrorWOMaterialDoubleRatioPi0FitDown[50];
    Double_t         relSystErrorADoubleRatioPi0FitUp[50];
    Double_t         relSystErrorADoubleRatioPi0FitDown[50];
    Double_t         relSystErrorBDoubleRatioPi0FitUp[50];
    Double_t         relSystErrorBDoubleRatioPi0FitDown[50];
    Double_t         relSystErrorCDoubleRatioPi0FitUp[50];
    Double_t         relSystErrorCDoubleRatioPi0FitDown[50];

    cout << fileNameSysErrDoubleRatioPi0Fit << endl;
    nPointsErrors=0;
    while(!fileSysErrDoubleRatioPi0Fit.eof() && nPointsErrors < 100){
        fileSysErrDoubleRatioPi0Fit     >> ptSysDoubleRatioFit[nPointsErrors]
                                        >> relSystErrorDoubleRatioPi0FitDown[nPointsErrors]                 >> relSystErrorDoubleRatioPi0FitUp[nPointsErrors]
                                        >> relSystErrorWOMaterialDoubleRatioPi0FitDown[nPointsErrors]         >> relSystErrorWOMaterialDoubleRatioPi0FitUp[nPointsErrors]
                                        >> relSystErrorADoubleRatioPi0FitDown[nPointsErrors]                 >> relSystErrorADoubleRatioPi0FitUp[nPointsErrors]
                                        >> relSystErrorBDoubleRatioPi0FitDown[nPointsErrors]                 >> relSystErrorBDoubleRatioPi0FitUp[nPointsErrors]
                                        >> relSystErrorCDoubleRatioPi0FitDown[nPointsErrors]                 >> relSystErrorCDoubleRatioPi0FitUp[nPointsErrors];
        cout << nPointsErrors << "\t"  << relSystErrorDoubleRatioPi0FitDown[nPointsErrors] << "\t" << relSystErrorDoubleRatioPi0FitUp[nPointsErrors] << "\t"
             << relSystErrorWOMaterialDoubleRatioPi0FitDown[nPointsErrors] << "\t"  <<relSystErrorWOMaterialDoubleRatioPi0FitUp[nPointsErrors] << endl;;
        nPointsErrors++;
    }
    fileSysErrDoubleRatioPi0Fit.close();

    // ******************************************************************
    // ** reading systematic errors for inclusive ratio with fitted pi0 *
    // ******************************************************************
    ifstream         fileSysErrIncRatioPi0Fit;
    fileSysErrIncRatioPi0Fit.open(fileNameSysErrIncRatioPi0Fit,ios_base::in);
    Double_t         relSystErrorIncRatioPi0FitUp[50];
    Double_t         relSystErrorIncRatioPi0FitDown[50];
    Double_t         relSystErrorWOMaterialIncRatioPi0FitUp[50];
    Double_t         relSystErrorWOMaterialIncRatioPi0FitDown[50];
    Double_t         relSystErrorAIncRatioPi0FitUp[50];
    Double_t         relSystErrorAIncRatioPi0FitDown[50];
    Double_t         relSystErrorBIncRatioPi0FitUp[50];
    Double_t         relSystErrorBIncRatioPi0FitDown[50];
    Double_t         relSystErrorCIncRatioPi0FitUp[50];
    Double_t         relSystErrorCIncRatioPi0FitDown[50];

    cout << fileNameSysErrIncRatioPi0Fit << endl;
    nPointsErrors=0;
    while(!fileSysErrIncRatioPi0Fit.eof() && nPointsErrors < 100){
        fileSysErrIncRatioPi0Fit    >> ptSysInclRatioFit[nPointsErrors]
                                    >> relSystErrorIncRatioPi0FitDown[nPointsErrors]                 >> relSystErrorIncRatioPi0FitUp[nPointsErrors]
                                    >> relSystErrorWOMaterialIncRatioPi0FitDown[nPointsErrors]         >> relSystErrorWOMaterialIncRatioPi0FitUp[nPointsErrors]
                                    >> relSystErrorAIncRatioPi0FitDown[nPointsErrors]                 >> relSystErrorAIncRatioPi0FitUp[nPointsErrors]
                                    >> relSystErrorBIncRatioPi0FitDown[nPointsErrors]                 >> relSystErrorBIncRatioPi0FitUp[nPointsErrors]
                                    >> relSystErrorCIncRatioPi0FitDown[nPointsErrors]                 >> relSystErrorCIncRatioPi0FitUp[nPointsErrors];
        cout << nPointsErrors << "\t"  << relSystErrorIncRatioPi0FitDown[nPointsErrors] << "\t"  << relSystErrorIncRatioPi0FitUp[nPointsErrors] << "\t"
             << relSystErrorWOMaterialIncRatioPi0FitDown[nPointsErrors] << "\t" << relSystErrorWOMaterialIncRatioPi0FitUp[nPointsErrors] << endl;;
        nPointsErrors++;
    }
    fileSysErrIncRatioPi0Fit.close();

    // ******************************************************************
    // ********* reading systematic errors for inclusive ratio **********
    // ******************************************************************
    ifstream         fileSysErrIncRatio;
    fileSysErrIncRatio.open(fileNameSysErrIncRatio,ios_base::in);
    Double_t         relSystErrorIncRatioUp[50];
    Double_t         relSystErrorIncRatioDown[50];
    Double_t         relSystErrorWOMaterialIncRatioUp[50];
    Double_t         relSystErrorWOMaterialIncRatioDown[50];
    Double_t         relSystErrorAIncRatioUp[50];
    Double_t         relSystErrorAIncRatioDown[50];
    Double_t         relSystErrorBIncRatioUp[50];
    Double_t         relSystErrorBIncRatioDown[50];
    Double_t         relSystErrorCIncRatioUp[50];
    Double_t         relSystErrorCIncRatioDown[50];

    cout << fileNameSysErrIncRatio << endl;
    nPointsErrors=0;
    while(!fileSysErrIncRatio.eof() && nPointsErrors < 100){
        fileSysErrIncRatio         >> ptSysInclRatio[nPointsErrors]
                                >> relSystErrorIncRatioDown[nPointsErrors]                     >> relSystErrorIncRatioUp[nPointsErrors]
                                >> relSystErrorWOMaterialIncRatioDown[nPointsErrors]         >> relSystErrorWOMaterialIncRatioUp[nPointsErrors]
                                >> relSystErrorAIncRatioDown[nPointsErrors]                 >> relSystErrorAIncRatioUp[nPointsErrors]
                                >> relSystErrorBIncRatioDown[nPointsErrors]                 >> relSystErrorBIncRatioUp[nPointsErrors]
                                >> relSystErrorCIncRatioDown[nPointsErrors]                 >> relSystErrorCIncRatioUp[nPointsErrors];
        cout << nPointsErrors << "\t"  << relSystErrorIncRatioDown[nPointsErrors] << "\t"  << relSystErrorIncRatioUp[nPointsErrors] << "\t"
             << relSystErrorWOMaterialIncRatioDown[nPointsErrors] << "\t" << relSystErrorWOMaterialIncRatioUp[nPointsErrors] << endl;;
        nPointsErrors++;
    }
    fileSysErrIncRatio.close();

    // ******************************************************************
    // ********* reading systematic errors for inclusive gamma **********
    // ******************************************************************
    ifstream         fileSysErrIncGamma;
    fileSysErrIncGamma.open(fileNameSysErrIncGamma,ios_base::in);
    Double_t         relSystErrorIncGammaUp[50];
    Double_t         relSystErrorIncGammaDown[50];
    Double_t         relSystErrorWOMaterialIncGammaUp[50];
    Double_t         relSystErrorWOMaterialIncGammaDown[50];
    Double_t         relSystErrorAIncGammaUp[50];
    Double_t         relSystErrorAIncGammaDown[50];
    Double_t         relSystErrorBIncGammaUp[50];
    Double_t         relSystErrorBIncGammaDown[50];
    Double_t         relSystErrorCIncGammaUp[50];
    Double_t         relSystErrorCIncGammaDown[50];

    cout << fileNameSysErrIncGamma << endl;
    nPointsErrors=0;
    while(!fileSysErrIncGamma.eof() && nPointsErrors < 100){
        fileSysErrIncGamma      >> ptSysGamma[nPointsErrors]
                                >> relSystErrorIncGammaDown[nPointsErrors]                  >> relSystErrorIncGammaUp[nPointsErrors]
                                >> relSystErrorWOMaterialIncGammaDown[nPointsErrors]        >> relSystErrorWOMaterialIncGammaUp[nPointsErrors]
                                >> relSystErrorAIncGammaDown[nPointsErrors]                 >> relSystErrorAIncGammaUp[nPointsErrors]
                                >> relSystErrorBIncGammaDown[nPointsErrors]                 >> relSystErrorBIncGammaUp[nPointsErrors]
                                >> relSystErrorCIncGammaDown[nPointsErrors]                 >> relSystErrorCIncGammaUp[nPointsErrors];
        cout << nPointsErrors << "\t"  << relSystErrorIncGammaDown[nPointsErrors] << "\t"  <<relSystErrorIncGammaUp[nPointsErrors] << "\t"
             << relSystErrorWOMaterialIncGammaDown[nPointsErrors] << "\t" << relSystErrorWOMaterialIncGammaUp[nPointsErrors] << endl;;
        nPointsErrors++;
    }
    fileSysErrIncGamma.close();
    nPointsErrors-2;

    // ******************************************************************
    // *************** reading systematic errors for pi0 ****************
    // ******************************************************************
    ifstream         fileSysErrPi0;
    fileSysErrPi0.open(fileNameSysErrPi0,ios_base::in);
    Double_t         relSystErrorPi0Up[50];
    Double_t         relSystErrorPi0Down[50];
    Double_t         relSystErrorWOMaterialPi0Up[50];
    Double_t         relSystErrorWOMaterialPi0Down[50];

    cout << fileNameSysErrPi0 << endl;
    nPointsErrors=0;
    while(!fileSysErrPi0.eof() && nPointsErrors < 100){
        fileSysErrPi0           >> ptSysPi0[nPointsErrors]
                                >> relSystErrorPi0Down[nPointsErrors]               >> relSystErrorPi0Up[nPointsErrors]
                                >> relSystErrorWOMaterialPi0Down[nPointsErrors]     >> relSystErrorWOMaterialPi0Up[nPointsErrors] ;
        cout << nPointsErrors << "\t"  << relSystErrorPi0Down[nPointsErrors] << "\t"  <<relSystErrorPi0Up[nPointsErrors] << "\t"
             << relSystErrorWOMaterialPi0Down[nPointsErrors] << "\t" << relSystErrorWOMaterialPi0Up[nPointsErrors] << endl;;
        nPointsErrors++;
    }
    fileSysErrIncGamma.close();

    // ******************************************************************
    // *********** reading systematic errors for pi0s fitted ************
    // ******************************************************************
    ifstream         fileSysErrPi0Fit;
    fileSysErrPi0Fit.open(fileNameSysErrPi0Fit,ios_base::in);
    Double_t         relSystErrorPi0FitUp[50];
    Double_t         relSystErrorPi0FitDown[50];
    Double_t         relSystErrorWOMaterialPi0FitUp[50];
    Double_t         relSystErrorWOMaterialPi0FitDown[50];

    cout << fileNameSysErrPi0Fit << endl;
    nPointsErrors=0;
    while(!fileSysErrPi0Fit.eof()){
        fileSysErrPi0Fit         >> ptSysPi0Fit[nPointsErrors]
                                 >> relSystErrorPi0FitDown[nPointsErrors]               >> relSystErrorPi0FitDown[nPointsErrors]
                                 >> relSystErrorWOMaterialPi0FitDown[nPointsErrors]     >> relSystErrorWOMaterialPi0FitUp[nPointsErrors] ;
        cout << nPointsErrors << "\t"  << relSystErrorPi0FitDown[nPointsErrors] << "\t"  <<relSystErrorPi0FitUp[nPointsErrors] << "\t"
             << relSystErrorWOMaterialPi0FitDown[nPointsErrors] << "\t" << relSystErrorWOMaterialPi0FitUp[nPointsErrors] << endl;;
        nPointsErrors++;
    }
    nPointsErrors--;
    fileSysErrPi0Fit.close();

    // ******************************************************************
    // ****************** reading theory graphs *************************
    // ******************************************************************
    TString fileNameTheoryPbPb                                  = "ExternalInputPbPb/Theory/TheoryCompilationPbPb.root";
    TFile* fileTheoryPbPb                           = new TFile( fileNameTheoryPbPb.Data());
    TDirectory* directoryTheoryGamma                = (TDirectory*)fileTheoryPbPb->Get("DirectPhoton");

    Double_t xSection2760GeV        = 55.416*1e-3;      // V0OR
    Double_t recalcBarn             = 1e12; //NLO in pbarn!!!!
    Double_t xSection2760GeVppINEL = 62.8*1e9;
    TGraphAsymmErrors *graphInvYieldPbPbTheoryCTEQ61EPS09       = NULL;
    TGraphAsymmErrors *graphInvYieldPPTheoryCT10BFG2_pdfErr     = NULL;
    TGraphAsymmErrors *graphInvYieldPbPbTheoryEPS09             = NULL;
    TGraphAsymmErrors *graphInvYieldPPTheoryCT10BFG2_scale      = NULL;
    TGraphAsymmErrors* graphDRPbPbCT10BFG2_pdfErr               = NULL;
    TGraphAsymmErrors* graphDRPbPbCTEQ61EPS09                   = NULL;
    graphInvYieldPbPbTheoryCTEQ61EPS09                          = (TGraphAsymmErrors*)directoryTheoryGamma->Get("PbPb276CTEQ61EPS09BFG2_sum_pdferr_InvYield");
    graphInvYieldPbPbTheoryCTEQ61EPS09 = (TGraphAsymmErrors*)ScaleGraph(graphInvYieldPbPbTheoryCTEQ61EPS09, fNcoll/(recalcBarn*xSection2760GeV));
    graphInvYieldPPTheoryCT10BFG2_pdfErr                        = (TGraphAsymmErrors*)directoryTheoryGamma->Get("pp276CT10BFG2_sum_pdferr_InvYield");
    graphInvYieldPPTheoryCT10BFG2_pdfErr = (TGraphAsymmErrors*)ScaleGraph(graphInvYieldPPTheoryCT10BFG2_pdfErr, fNcoll/(recalcBarn*xSection2760GeV));

    graphInvYieldPbPbTheoryEPS09                                = (TGraphAsymmErrors*)directoryTheoryGamma->Get("PbPb276EPS09BFG2_sum_scale_InvYield");
    graphInvYieldPbPbTheoryEPS09 = (TGraphAsymmErrors*)ScaleGraph(graphInvYieldPbPbTheoryEPS09, fNcoll/(recalcBarn*xSection2760GeV));
    graphInvYieldPPTheoryCT10BFG2_scale                         = (TGraphAsymmErrors*)directoryTheoryGamma->Get("pp276CT10BFG2_sum_scale_InvYield");
    graphInvYieldPPTheoryCT10BFG2_scale = (TGraphAsymmErrors*)ScaleGraph(graphInvYieldPPTheoryCT10BFG2_scale, fNcoll/(recalcBarn*xSection2760GeV));

    Bool_t activateTheoryPbPb                                   = kFALSE;
    if (graphInvYieldPbPbTheoryCTEQ61EPS09 && graphInvYieldPPTheoryCT10BFG2_pdfErr && graphInvYieldPbPbTheoryEPS09 &&   graphInvYieldPPTheoryCT10BFG2_scale)
        activateTheoryPbPb                                      = kTRUE;
//     graphInvYieldPbPbTheoryCTEQ61EPS09->RemovePoint(graphInvYieldPbPbTheoryCTEQ61EPS09->GetN()-1);
//     graphInvYieldPPTheoryCT10BFG2_pdfErr->RemovePoint(graphInvYieldPPTheoryCT10BFG2_pdfErr->GetN()-1);
//
//     graphInvYieldPbPbTheoryEPS09->RemovePoint(graphInvYieldPbPbTheoryEPS09->GetN()-1);
//     graphInvYieldPPTheoryCT10BFG2_scale->RemovePoint(graphInvYieldPPTheoryCT10BFG2_scale->GetN()-1);

    if (activateTheoryPbPb){
        for(Int_t i = 0;i<graphInvYieldPbPbTheoryCTEQ61EPS09->GetN();i++){
            Double_t yerrlow1         = graphInvYieldPbPbTheoryCTEQ61EPS09->GetErrorYlow(i);
            Double_t yerrlow2         = 1*graphInvYieldPbPbTheoryEPS09->GetErrorYlow(i);
            Double_t yerrhigh1         = graphInvYieldPbPbTheoryCTEQ61EPS09->GetErrorYhigh(i);
            Double_t yerrhigh2         = 1*graphInvYieldPbPbTheoryEPS09->GetErrorYhigh(i);
            Double_t xerrlow         = graphInvYieldPbPbTheoryEPS09->GetErrorXhigh(i);

            graphInvYieldPbPbTheoryCTEQ61EPS09->SetPointError( i,xerrlow,xerrlow,
                                                        sqrt(yerrlow1*yerrlow1+ yerrlow2*yerrlow2),
                                                        sqrt(yerrhigh1*yerrhigh1+ yerrhigh2*yerrhigh2));
        }

        cout << "======================================================================= DR EPS09 =================================================================" << endl;
        graphDRPbPbCTEQ61EPS09                                     = (TGraphAsymmErrors*)graphInvYieldPbPbTheoryCTEQ61EPS09->Clone("graphDRPbPbCTEQ61EPS09");
        for (Int_t bin = 0; bin < graphDRPbPbCTEQ61EPS09->GetN(); bin++){
            Double_t cocktailIntegral     = cocktailFitAllGammaForNLO->Integral(graphDRPbPbCTEQ61EPS09->GetX()[bin]-graphDRPbPbCTEQ61EPS09->GetErrorXlow(bin),
                                                                            graphDRPbPbCTEQ61EPS09->GetX()[bin]+graphDRPbPbCTEQ61EPS09->GetErrorXhigh(bin))
                                        /(graphDRPbPbCTEQ61EPS09->GetErrorXlow(bin)+graphDRPbPbCTEQ61EPS09->GetErrorXhigh(bin));
            Double_t DRtheo                 = 1 + (graphDRPbPbCTEQ61EPS09->GetY()[bin]/ cocktailIntegral);
            Double_t DRtheoErrDown        = graphDRPbPbCTEQ61EPS09->GetErrorYlow(bin)/ cocktailIntegral;
            Double_t DRtheoErrUp        = graphDRPbPbCTEQ61EPS09->GetErrorYhigh(bin)/ cocktailIntegral;
            cout << graphDRPbPbCTEQ61EPS09->GetY()[bin] << endl;
            cout << graphDRPbPbCTEQ61EPS09->GetX()[bin]<< "\t" << cocktailIntegral << "\t" << DRtheo << "\t +" << DRtheoErrUp << "\t -" <<    DRtheoErrDown << endl;
            graphDRPbPbCTEQ61EPS09->SetPoint(bin, graphDRPbPbCTEQ61EPS09->GetX()[bin], DRtheo );
            graphDRPbPbCTEQ61EPS09->SetPointError(bin, graphDRPbPbCTEQ61EPS09->GetErrorXlow(bin), graphDRPbPbCTEQ61EPS09->GetErrorXhigh(bin), DRtheoErrDown, DRtheoErrUp );
        }
        graphDRPbPbCTEQ61EPS09->Print();

        cout << "======================================================================= DR CT10BFG2 =================================================================" << endl;
        for(Int_t i = 0;i<graphInvYieldPPTheoryCT10BFG2_scale->GetN();i++){
            Double_t yerrlow1         = graphInvYieldPPTheoryCT10BFG2_pdfErr->GetErrorYlow(i);
            Double_t yerrlow2         = graphInvYieldPPTheoryCT10BFG2_scale->GetErrorYlow(i);
            Double_t yerrhigh1         = graphInvYieldPPTheoryCT10BFG2_pdfErr->GetErrorYhigh(i);
            Double_t yerrhigh2         = graphInvYieldPPTheoryCT10BFG2_scale->GetErrorYhigh(i);
            Double_t xerrlow         = graphInvYieldPPTheoryCT10BFG2_scale->GetErrorXhigh(i);

            graphInvYieldPPTheoryCT10BFG2_pdfErr->SetPointError(i,xerrlow,xerrlow,
                                                                sqrt(yerrlow1*yerrlow1+ yerrlow2*yerrlow2),
                                                                sqrt(yerrhigh1*yerrhigh1+ yerrhigh2*yerrhigh2));
        }
        graphDRPbPbCT10BFG2_pdfErr                                 = (TGraphAsymmErrors*)graphInvYieldPPTheoryCT10BFG2_pdfErr->Clone("graphDRPbPbCT10BFG2_pdfErr");
        for (Int_t bin = 0; bin < graphDRPbPbCT10BFG2_pdfErr->GetN(); bin++){
            Double_t cocktailIntegral     = cocktailFitAllGammaForNLO->Integral(graphDRPbPbCT10BFG2_pdfErr->GetX()[bin]-graphDRPbPbCT10BFG2_pdfErr->GetErrorXlow(bin),
                                                                            graphDRPbPbCT10BFG2_pdfErr->GetX()[bin]+graphDRPbPbCT10BFG2_pdfErr->GetErrorXhigh(bin))
                                        /(graphDRPbPbCT10BFG2_pdfErr->GetErrorXlow(bin)+graphDRPbPbCT10BFG2_pdfErr->GetErrorXhigh(bin));
            Double_t DRtheo                 = 1 + (graphDRPbPbCT10BFG2_pdfErr->GetY()[bin]/ cocktailIntegral);
            Double_t DRtheoErrDown        = graphDRPbPbCT10BFG2_pdfErr->GetErrorYlow(bin)/ cocktailIntegral;
            Double_t DRtheoErrUp        = graphDRPbPbCT10BFG2_pdfErr->GetErrorYhigh(bin)/ cocktailIntegral;
    //         cout << graphDRPbPbCT10BFG2_pdfErr->GetX()[bin]<< "\t" << cocktailIntegral << "\t" << DRtheo << "\t +" << DRtheoErrUp << "\t -" <<    DRtheoErrDown << endl;
            graphDRPbPbCT10BFG2_pdfErr->SetPoint(bin, graphDRPbPbCT10BFG2_pdfErr->GetX()[bin], DRtheo );
            graphDRPbPbCT10BFG2_pdfErr->SetPointError(bin, graphDRPbPbCT10BFG2_pdfErr->GetErrorXlow(bin), graphDRPbPbCT10BFG2_pdfErr->GetErrorXhigh(bin), DRtheoErrDown, DRtheoErrUp );
        }
        graphDRPbPbCT10BFG2_pdfErr->Print();
    }

    // ******************************************************************
    // ***** read NL0 calculations and put them in proper format ********
    // ******************************************************************
    TGraphErrors *NLODoubleRatio     = (TGraphErrors*) fileInput->Get("graphNLODoubleRatio");
    TGraphErrors *NLO                = (TGraphErrors*) fileInput->Get("graphNLODirGamma");
    SetStyleGammaNLOTGraphWithBand( NLODoubleRatio, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);
    SetStyleGammaNLOTGraphWithBand( NLO, 3.0, 1, colorNLOcalc, 1001, colorNLOcalc, 0);

    // **************************************************************************
    // ******************** Draw Final Gamma Pictures ***************************
    // **************************************************************************
    TF1 *lineOne = new TF1("lineOne","1",0,16);
    lineOne->SetLineWidth(1.2);
    lineOne->SetLineColor(kGray+2);

    // **************************************************************************
    // ***************** Calculate systematic error graphs **********************
    // **************************************************************************
    TGraphAsymmErrors* graphDoubleRatioFitPi0SysErr     = CalculateSysErrFromRelSysHisto( histoDRFit, "DoubleRatioPi0FitSystError", relSystErrorDoubleRatioPi0FitDown, relSystErrorDoubleRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphDoubleRatioFitPi0SysErrA    = CalculateSysErrFromRelSysHisto( histoDRFit, "DoubleRatioPi0FitSystErrorA", relSystErrorADoubleRatioPi0FitDown, relSystErrorADoubleRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphDoubleRatioFitPi0SysErrB    = CalculateSysErrFromRelSysHisto( histoDRFit, "DoubleRatioPi0FitSystErrorB", relSystErrorBDoubleRatioPi0FitDown, relSystErrorBDoubleRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphDoubleRatioFitPi0SysErrC    = CalculateSysErrFromRelSysHisto( histoDRFit, "DoubleRatioPi0FitSystErrorC", relSystErrorCDoubleRatioPi0FitDown, relSystErrorCDoubleRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphDoubleRatioSysErr  = CalculateSysErrFromRelSysHisto( histoDR, "DoubleRatioSystError", relSystErrorDoubleRatioDown, relSystErrorDoubleRatioUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphDoubleRatioSysErrA = CalculateSysErrFromRelSysHisto( histoDR, "DoubleRatioSystErrorA", relSystErrorADoubleRatioDown, relSystErrorADoubleRatioUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphDoubleRatioSysErrB = CalculateSysErrFromRelSysHisto( histoDR, "DoubleRatioSystErrorB", relSystErrorBDoubleRatioDown, relSystErrorBDoubleRatioUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphDoubleRatioSysErrC = CalculateSysErrFromRelSysHisto( histoDR, "DoubleRatioSystErrorC", relSystErrorCDoubleRatioDown, relSystErrorCDoubleRatioUp, binOffset, nPointsErrors);

    //************************************************************************
    //******************* Calculate error graph for inclusive ratio **********
    //************************************************************************
    TGraphAsymmErrors* graphIncRatioFitPi0SysErr      = CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystError", relSystErrorIncRatioPi0FitDown, relSystErrorIncRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncRatioFitPi0SysErrA     = CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystErrorA", relSystErrorAIncRatioPi0FitDown, relSystErrorAIncRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncRatioFitPi0SysErrB     = CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystErrorB", relSystErrorBIncRatioPi0FitDown, relSystErrorBIncRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncRatioFitPi0SysErrC     = CalculateSysErrFromRelSysHisto( histoIncRatioPi0Fit, "IncRatioPi0FitSystErrorC", relSystErrorCIncRatioPi0FitDown, relSystErrorCIncRatioPi0FitUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncRatioSysErr  = CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystError", relSystErrorIncRatioDown, relSystErrorIncRatioUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncRatioSysErrA = CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystErrorA", relSystErrorAIncRatioDown, relSystErrorAIncRatioUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncRatioSysErrB = CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystErrorB", relSystErrorBIncRatioDown, relSystErrorBIncRatioUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncRatioSysErrC = CalculateSysErrFromRelSysHisto( histoIncRatio, "IncRatioSystErrorC", relSystErrorCIncRatioDown, relSystErrorCIncRatioUp, binOffset, nPointsErrors);

    // **************************************************************************
    // ***************** Calculate systematic error graphs **********************
    // **************************************************************************
    TGraphAsymmErrors* graphIncGammaSysErr      = CalculateSysErrFromRelSysHisto( histoIncGamma, "IncGammaSystError", relSystErrorIncGammaDown, relSystErrorIncGammaUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncGammaSysErrA     = CalculateSysErrFromRelSysHisto( histoIncGamma, "IncGammaSystErrorA", relSystErrorAIncGammaDown, relSystErrorAIncGammaUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncGammaSysErrB     = CalculateSysErrFromRelSysHisto( histoIncGamma, "IncGammaSystErrorB", relSystErrorBIncGammaDown, relSystErrorBIncGammaUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncGammaSysErrC     = CalculateSysErrFromRelSysHisto( histoIncGamma, "IncGammaSystErrorC", relSystErrorCIncGammaDown, relSystErrorCIncGammaUp, binOffset, nPointsErrors);
    TGraphAsymmErrors* graphIncGammaSysErrW0Mat = CalculateSysErrFromRelSysHisto( histoIncGamma, "IncGammaSystErrorWOMat", relSystErrorWOMaterialIncGammaDown, relSystErrorWOMaterialIncGammaUp, binOffset, nPointsErrors);

    //************************************************************************
    //******************* Calculate error graph for pi0 **********************
    //************************************************************************
    TGraphAsymmErrors* graphPi0SysErr    = CalculateSysErrFromRelSysHisto( histoPi0Spectrum, "Pi0SystError", relSystErrorPi0Down, relSystErrorPi0Up, binOffset, nPointsErrors);
//     TGraphAsymmErrors* graphPi0FitSysErr = CalculateSysErrFromRelSysHisto( histoPi0SpectrumFit, "Pi0FitSystError", relSystErrorPi0FitDown, relSystErrorPi0FitUp, binOffset, nPointsErrors);

    cout << "*************************************************************************************************************" << endl;

    //****************************************************************************
    //******************* draw double ratio with pi0 fitted **********************
    //****************************************************************************
    TCanvas *canvasDoubleRatio = GetAndSetCanvas("canvasDoubleRatioFinal");
    TH2D *dummyDR = new TH2D("dummyDR", "dummyDR", 120, 0., 16, 1000., minYDR, maxYDR);
    SetStyleHistoTH2ForGraphs( dummyDR, "#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})",
                               0.045, 0.05, 0.045, 0.05, 0.85, 0.85);
    dummyDR->GetXaxis()->SetRangeUser(minPt, maxPt);
    dummyDR->DrawCopy();

//     for (Int_t i = 1; i< histoDRFit->GetNbinsX()+1; i++){
//         cout << i << "\t" << histoDRFit->GetBinCenter(i) << "\t" << histoDRFit->GetBinContent(i) << "\t" << histoDRFit->GetBinError(i) << endl;
//     }
//     graphDoubleRatioFitPi0SysErr->Print();

    DrawGammaSetMarker(histoDRFit, 20, 2., colorCent, colorCent);
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErr , 20, 2, colorCent, colorCent, 1, kTRUE);
    lineOne->Draw("same");
    NLODoubleRatio->Draw("p3lsame");

    graphDoubleRatioFitPi0SysErr->Draw("E2same");
    histoDRFit->DrawCopy("same");

    TLegend* legendDoubleRatio = GetAndSetLegend(0.15,0.75,4);
    legendDoubleRatio->AddEntry(graphDoubleRatioFitPi0SysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
    legendDoubleRatio->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
    legendDoubleRatio->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
    legendDoubleRatio->Draw();

    canvasDoubleRatio->Print(Form("%s/DoubleRatioPi0Fitted_%s.eps",outputDir.Data(),centrality.Data()));

    //****************************************************************************
    // draw double ratio with pi0 fitted and compare to double ratio with pure pi0
    //****************************************************************************
    canvasDoubleRatio->cd();
    dummyDR->DrawCopy();

    lineOne->Draw("same");
    NLODoubleRatio->Draw("p3lsame");

    DrawGammaSetMarker(histoDR, 20, 2., colorCentNotPi0Fitted, colorCentNotPi0Fitted);
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioSysErr , 20, 2, colorCentNotPi0Fitted, colorCentNotPi0Fitted, 1, kTRUE);
    graphDoubleRatioSysErr->Draw("E2same");
    histoDR->DrawCopy("same");

    graphDoubleRatioFitPi0SysErr->Draw("E2same");
    histoDRFit->DrawCopy("same");

    TLegend* legendDoubleRatio2 = GetAndSetLegend(0.15,0.75,4);
    legendDoubleRatio2->AddEntry(graphDoubleRatioFitPi0SysErr,Form("PCM, %s, #pi^{0} fitted",collisionSystem.Data()),"pf");
    legendDoubleRatio2->AddEntry(graphDoubleRatioSysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
    legendDoubleRatio2->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
    legendDoubleRatio2->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
    legendDoubleRatio2->Draw();

    canvasDoubleRatio->Print(Form("%s/DoubleRatio_ComparedPi0FittedAndNot_%s.eps",outputDir.Data(),centrality.Data()));

    canvasDoubleRatio->cd();
    dummyDR->DrawCopy();

    graphDoubleRatioSysErr->Draw("E2same");
    histoDR->DrawCopy("same");

    lineOne->Draw("same");
    NLODoubleRatio->Draw("p3lsame");

    graphPublishedDoubleRatioSyst->Draw("E2same");
    graphPublishedDoubleRatioStat->Draw("p,e1,same");

    TLegend* legendDoubleRatioWP = GetAndSetLegend(0.15,0.75,4);
    legendDoubleRatioWP->AddEntry(graphDoubleRatioSysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
    legendDoubleRatioWP->AddEntry(graphPublishedDoubleRatioSyst,Form("PCM published, %s",collisionSystem.Data()),"pf");
    legendDoubleRatioWP->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
    legendDoubleRatioWP->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
    legendDoubleRatioWP->Draw();

    canvasDoubleRatio->Print(Form("%s/DoubleRatioPi0_%s_withPublished.eps",outputDir.Data(),centrality.Data()));


    //****************************************************************************
    //*************** Double Ratio with pi0 fitted with ABC errors ***************
    //****************************************************************************
    TCanvas *canvasDoubleRatioABC = GetAndSetCanvas("canvasDoubleRatioFinal");
    dummyDR->GetXaxis()->SetRangeUser(minPt, maxPt);
    dummyDR->DrawCopy();

    lineOne->Draw("same");
    NLODoubleRatio->Draw("p3lsame");
    histoDRFit->DrawCopy("same");
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErrA , 20, 1, colorErrA, colorErrA, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErrB , 20, 1, colorErrB, colorErrB, 1, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphDoubleRatioFitPi0SysErrC , 20, 1, colorErrC, colorErrC, 1, kTRUE);

    graphDoubleRatioFitPi0SysErrA->Draw("E2same");
    graphDoubleRatioFitPi0SysErrB->Draw("E2same");
    graphDoubleRatioFitPi0SysErrC->Draw("E2same");
    graphDoubleRatioFitPi0SysErr->Draw("E2same");

    TLegend* legendDoubleRatioABC = GetAndSetLegend(0.15,0.65,6);
    legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
    legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErrA,"Error Type A","f");
    legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErrB,"Error Type B","f");
    legendDoubleRatioABC->AddEntry(graphDoubleRatioFitPi0SysErrC,"Error Type C","f");
    legendDoubleRatioABC->AddEntry(NLODoubleRatio,"NLO prediction: 1 + (#it{N}_{coll}#it{N}_{#gamma_{direct,pp,NLO}}/#it{N}_{#gamma_{decay}})","l");
    legendDoubleRatioABC->AddEntry((TObject*)0, "for #mu = 0.5 to 2.0 #it{p}_{T}", "");
    legendDoubleRatioABC->Draw();

    canvasDoubleRatioABC->Print(Form("%s/DoubleRatioPi0Fitted_ABCErrors_%s.eps",outputDir.Data(),centrality.Data()));


    // **************************************************************************
    // ******************** plotting inclusive ratio ****************************
    // **************************************************************************
    TCanvas *canvasIncRatio = GetAndSetCanvas("canvasIncRatioFinal");
    TH2D *dummyIncR = new TH2D("dummyIncR", "dummyIncR", 120, 0., 16, 1000., minYIR, maxYIR);
    SetStyleHistoTH2ForGraphs( dummyIncR, "#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}",0.045, 0.05, 0.045, 0.05, 0.85, 0.85);
    dummyIncR->GetXaxis()->SetRangeUser(minPt, maxPt);
    dummyIncR->DrawCopy();

    DrawGammaSetMarker(histoIncRatio, 20, 2., colorCentNotPi0Fitted, colorCentNotPi0Fitted);
    DrawGammaSetMarkerTGraphAsym(graphIncRatioSysErr , 20, 2, colorCentNotPi0Fitted, colorCentNotPi0Fitted, 1, kTRUE);
    graphIncRatioSysErr->Draw("E2same");
    histoIncRatio->DrawCopy("same");

    TLegend* legendIncRatio = GetAndSetLegend(0.12,0.2,1.5);
    legendIncRatio->AddEntry(graphIncRatioSysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
    legendIncRatio->Draw();

    canvasIncRatio->Print(Form("%s/IncRatio_%s.eps",outputDir.Data(),centrality.Data()));

    canvasIncRatio->cd();
    dummyIncR->DrawCopy();

    graphIncRatioSysErr->Draw("E2same");
    histoIncRatio->DrawCopy("same");

    graphPublishedInclRatioSpectrumSyst->Draw("E2same");
    graphPublishedInclRatioSpectrumStat->Draw("p,e1,same");

    TLegend* legendIncRatioWP = GetAndSetLegend(0.12,0.2,1.5);
    legendIncRatioWP->AddEntry(graphIncRatioSysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
    legendIncRatioWP->AddEntry(graphPublishedInclRatioSpectrumSyst,Form("PCM published, %s",collisionSystem.Data()),"pf");
    legendIncRatioWP->Draw();

    canvasIncRatio->Print(Form("%s/IncRatio_%s_withPublished.eps",outputDir.Data(),centrality.Data()));


    // **************************************************************************
    // ****************** plotting inclusive ratio with pi0 fitted **************
    // **************************************************************************
    TCanvas *canvasIncRatioPi0Fit = GetAndSetCanvas("canvasIncRatioPi0FitFinal");
    dummyIncR->DrawCopy();

    DrawGammaSetMarker(histoIncRatioPi0Fit, 20, 2., colorCent, colorCent);
    DrawGammaSetMarkerTGraphAsym(graphIncRatioFitPi0SysErr , 20, 2, colorCent, colorCent, 1, kTRUE);
    graphIncRatioFitPi0SysErr->Draw("E2same");
    histoIncRatioPi0Fit->DrawCopy("same");

    TLegend* legendIncRatioFit = GetAndSetLegend(0.12,0.2,1.5);
    legendIncRatioFit->AddEntry(graphIncRatioFitPi0SysErr,Form("PCM, %s",collisionSystem.Data()),"pf");
    legendIncRatioFit->Draw();

    canvasIncRatioPi0Fit->Print(Form("%s/IncRatioPi0Fitted_%s.eps",outputDir.Data(),centrality.Data()));

    // **************************************************************************
    // ******************** plotting pi0 spectrum *******************************
    // **************************************************************************
    TCanvas *canvasPi0 = GetAndSetCanvas("canvasPi0Final");
    canvasPi0->SetLeftMargin(0.14);
    canvasPi0->SetLogy();
    canvasPi0->SetLogx();

    TH2D *dummyPi0 = new TH2D("dummyPi0", "dummyPi0", 120, 0., 16, 1000.,1e-10,1e4);
    SetStyleHistoTH2ForGraphs( dummyPi0, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.045, 0.05, 0.045, 0.05, 0.85, 1.2);
    dummyPi0->GetXaxis()->SetLabelOffset(-0.015);
    dummyPi0->GetXaxis()->SetRangeUser(minPt+0.3, maxPt+2);
    dummyPi0->DrawCopy();

    DrawGammaSetMarker(histoPi0Spectrum, 20, 2., colorCentNotPi0Fitted, colorCentNotPi0Fitted);
    DrawGammaSetMarkerTGraphAsym(graphPi0SysErr , 20, 2, colorCentNotPi0Fitted, colorCentNotPi0Fitted, 1, kTRUE);
    graphPi0SysErr->Draw("E2same");
    histoPi0Spectrum->DrawCopy("same");

    TLegend* legendPi0 = GetAndSetLegend(0.18,0.2,1.5);
    legendPi0->AddEntry(graphPi0SysErr,Form("#pi^{0} PCM, %s",collisionSystem.Data()),"pf");
    legendPi0->Draw();

    canvasPi0->Print(Form("%s/Pi0Spectrum_%s.eps",outputDir.Data(),centrality.Data()));

    // **************************************************************************
    // ********************** plotting fiited pi0 spectrum **********************
    // **************************************************************************
//     dummyPi0->DrawCopy();

//     DrawGammaSetMarker(histoPi0SpectrumFit, 20, 2., colorCent, colorCent);
//     DrawGammaSetMarkerTGraphAsym(graphPi0FitSysErr , 20, 2, colorCent, colorCent, 1, kTRUE);
//     graphPi0FitSysErr->Draw("E2same");
//     histoPi0SpectrumFit->DrawCopy("same");

//     TLegend* legendPi0Fit = GetAndSetLegend(0.18,0.2,1.5);
//     legendPi0Fit->AddEntry(graphPi0FitSysErr,Form("#pi^{0} PCM, %s",collisionSystem.Data()),"pf");

//     legendPi0Fit->Draw();

//     canvasPi0->Print(Form("%s/Pi0SprectrumFitted_%s.eps",outputDir.Data(),centrality.Data()));


    // **************************************************************************
    // ******************* Plotting inclusive Gamma Spectrum ********************
    // **************************************************************************
    TCanvas *canvasIncGamma = GetAndSetCanvas("canvasIncGamma");
    canvasIncGamma->SetLeftMargin(0.14);
    canvasIncGamma->SetLogy();
    canvasIncGamma->SetLogx();

    TH2D *dummyGamma ;
    dummyGamma = new TH2D("dummyGamma", "dummyGamma", 120, 0., 16, 1000.,1e-10,1e4);
    SetStyleHistoTH2ForGraphs( dummyGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.045, 0.05, 0.045, 0.05, 0.85, 1.2);
    dummyGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyGamma->GetXaxis()->SetRangeUser(minPt, maxPt);
    dummyGamma->GetYaxis()->SetRangeUser(1e-6,1e2);
    dummyGamma->DrawCopy();


    DrawGammaSetMarker(histoIncGamma, 20, 2., colorCent, colorCent);
    DrawGammaSetMarkerTGraphAsym(graphIncGammaSysErr , 20, 2, colorCent, colorCent, 1, kTRUE);
    graphIncGammaSysErr->Draw("E2same");
    histoIncGamma->DrawCopy("same");

    TLegend* legendIncGamma = GetAndSetLegend(0.18,0.2,1.5);
    legendIncGamma->AddEntry(graphIncGammaSysErr,Form("#gamma_{inc} PCM, %s",collisionSystem.Data()),"pf");
    legendIncGamma->Draw();

    canvasIncGamma->Print(Form("%s/IncGammaSpectrum_%s.eps",outputDir.Data(),centrality.Data()));

   canvasIncGamma->cd();
    dummyGamma->DrawCopy();

    graphIncGammaSysErr->Draw("E2same");
    histoIncGamma->DrawCopy("same");

    graphPublishedInclGammaSpectrumSyst->Draw("E2same");
    graphPublishedInclGammaSpectrumStat->Draw("p,same,e1");

    TLegend* legendIncGammaWP = GetAndSetLegend(0.18,0.2,1.5);
    legendIncGammaWP->AddEntry(graphIncGammaSysErr,Form("#gamma_{inc} PCM, %s",collisionSystem.Data()),"pf");
    legendIncGammaWP->AddEntry(graphPublishedInclGammaSpectrumSyst,Form("#gamma_{inc} PCM published, %s",collisionSystem.Data()),"pf");
    legendIncGammaWP->Draw();

    canvasIncGamma->Print(Form("%s/IncGammaSpectrum_%s_withPublished.eps",outputDir.Data(),centrality.Data()));


    // **************************************************************************
    // ***************** Calculate direct photon spectrum ***********************
    // **************************************************************************

    //_______________________ copy inclusive photon spectra _____________________
    TH1D *histoDirGammaSpectrumError               = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrum");
    TH1D *histoDirGammaSpectrumErrorSummed         = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrumSummed");
    TH1D *histoDirGammaSpectrumErrorSyst           = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrumSyst");
    TH1D *histoDirGammaSpectrumErrorSystA          = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrumSystA");
    TH1D *histoDirGammaSpectrumErrorSystB          = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrumSystB");
    TH1D *histoDirGammaSpectrumErrorSystC          = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrumSystC");
    TH1D *histoDirGammaSpectrumErrorSystWithoutMat = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrumSyst");
    TH1D *histoDirGammaSpectrumErrorStat           = (TH1D*)histoIncGamma->Clone("DirectPhotonSpectrumStat");

    //_______________________ get arrays of double ratio errors __________________
    Double_t *systErrorsDoubleRatio                = new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
    Double_t *systErrorsDoubleRatioA               = new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
    Double_t *systErrorsDoubleRatioB               = new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
    Double_t *systErrorsDoubleRatioC               = new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
    Double_t *systErrorsDoubleRatioX               = new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];
    Double_t *systErrorsWithoutMatDoubleRatio      = new Double_t[graphDoubleRatioFitPi0SysErr->GetN()];

    for (Int_t i = 0; i< graphDoubleRatioFitPi0SysErr->GetN(); i++){
        systErrorsDoubleRatio[i]                   = relSystErrorDoubleRatioPi0FitUp[i];
        systErrorsDoubleRatioA[i]                  = relSystErrorADoubleRatioPi0FitUp[i];
        systErrorsDoubleRatioB[i]                  = relSystErrorBDoubleRatioPi0FitUp[i];
        systErrorsDoubleRatioC[i]                  = relSystErrorCDoubleRatioPi0FitUp[i];
        systErrorsWithoutMatDoubleRatio[i]         = relSystErrorWOMaterialDoubleRatioPi0FitUp[i];
    }

    systErrorsDoubleRatioX =  graphDoubleRatioFitPi0SysErr->GetX();

    //_______________________ copy inclusive photon spectra _____________________
    TH1D* histoDRWithSummedErrors           = (TH1D*) histoDRFit->Clone("DoubleRatioWithSummedErrors");
    TH1D* histoDRWithStatErrors             = (TH1D*) histoDRFit->Clone("DoubleRatioWithStatErrors");
    TH1D* histoDRWithSystErrors             = (TH1D*) histoDRFit->Clone("DoubleRatioWithSystErrors");
    TH1D* histoDRWithSystErrorsA            = (TH1D*) histoDRFit->Clone("DoubleRatioWithSystErrorsA");
    TH1D* histoDRWithSystErrorsB            = (TH1D*) histoDRFit->Clone("DoubleRatioWithSystErrorsB");
    TH1D* histoDRWithSystErrorsC            = (TH1D*) histoDRFit->Clone("DoubleRatioWithSystErrorsC");
    TH1D* histoDRWithSystErrorsWithoutMat   = (TH1D*) histoDRFit->Clone("DoubleRatioWithSystErrors");

    for(Int_t i = 0; i<histoDRWithSummedErrors->GetNbinsX();i++){
        //cout<<systErrorsDoubleRatioX[i]<<"  "<<histoDRWithSummedErrors->GetBinCenter(i+binOffset)<<endl;
        Double_t binErrorSummed = sqrt( pow( (histoDRWithSummedErrors->GetBinError(i+binOffset)/histoDRWithSummedErrors->GetBinContent(i+binOffset))*100,2) + pow(systErrorsDoubleRatio[i],2) );
        Double_t binErrorSyst = systErrorsDoubleRatio[i];
        Double_t binErrorSystA = systErrorsDoubleRatioA[i];
        Double_t binErrorSystB = systErrorsDoubleRatioB[i];
        Double_t binErrorSystC = systErrorsDoubleRatioC[i];
        Double_t binErrorSystWitoutMat = systErrorsWithoutMatDoubleRatio[i];
        Double_t binErrorStat = (histoDRWithSummedErrors->GetBinError(i+binOffset)/histoDRWithSummedErrors->GetBinContent(i+binOffset))*100;

        histoDRWithSummedErrors->SetBinError(i+binOffset,(binErrorSummed/100)*histoDRWithSummedErrors->GetBinContent(i+binOffset));
        histoDRWithSystErrors->SetBinError(i+binOffset,(binErrorSyst/100)*histoDRWithSystErrors->GetBinContent(i+binOffset));
        histoDRWithSystErrorsA->SetBinError(i+binOffset,(binErrorSystA/100)*histoDRWithSystErrorsA->GetBinContent(i+binOffset));
        histoDRWithSystErrorsB->SetBinError(i+binOffset,(binErrorSystB/100)*histoDRWithSystErrorsB->GetBinContent(i+binOffset));
        histoDRWithSystErrorsC->SetBinError(i+binOffset,(binErrorSystC/100)*histoDRWithSystErrorsC->GetBinContent(i+binOffset));
        histoDRWithSystErrorsWithoutMat->SetBinError(i+binOffset,(binErrorSystWitoutMat/100)*histoDRWithSystErrorsWithoutMat->GetBinContent(i+binOffset));
//         histoDRWithStatErrors->SetBinError(i+binOffset,(binErrorStat/100)*histoDRWithStatErrors->GetBinContent(i+binOffset));

    }

    for(Int_t i = 0; i<histoIncGamma->GetNbinsX();i++){
        histoDirGammaSpectrumError->SetBinContent(i+binOffset,-1);
        histoDirGammaSpectrumErrorSummed->SetBinContent(i+binOffset,-1);
        histoDirGammaSpectrumErrorSyst->SetBinContent(i+binOffset,-1);
        histoDirGammaSpectrumErrorSystA->SetBinContent(i+binOffset,-1);
        histoDirGammaSpectrumErrorSystB->SetBinContent(i+binOffset,-1);
        histoDirGammaSpectrumErrorSystC->SetBinContent(i+binOffset,-1);
        histoDirGammaSpectrumErrorStat->SetBinContent(i+binOffset,-1);

        histoDirGammaSpectrumError->SetBinError(i+binOffset,0);
        histoDirGammaSpectrumErrorSummed->SetBinError(i+binOffset,0);
        histoDirGammaSpectrumErrorSyst->SetBinError(i+binOffset,0);
        histoDirGammaSpectrumErrorSystA->SetBinError(i+binOffset,0);
        histoDirGammaSpectrumErrorSystB->SetBinError(i+binOffset,0);
        histoDirGammaSpectrumErrorSystC->SetBinError(i+binOffset,0);
        histoDirGammaSpectrumErrorStat->SetBinError(i+binOffset,0);
    }

    // get the binning of the direct photons from the DR
    TH1D*histoDirGammaSpectrumSyst                 = (TH1D*)histoDRWithSystErrors->Clone();
         histoDirGammaSpectrumSyst->Reset();
    TH1D*histoDirGammaSpectrumStat                 = (TH1D*)histoDRWithStatErrors->Clone();
         histoDirGammaSpectrumStat->Reset();
    TH1D*histoDirGammaSpectrumSummed               = (TH1D*)histoDRWithSummedErrors->Clone();
         histoDirGammaSpectrumSummed->Reset();
    TH1D*histoDirGammaSpectrumSystWithoutMat       = (TH1D*)histoDRWithSystErrorsWithoutMat->Clone();
         histoDirGammaSpectrumSystWithoutMat->Reset();

    for(Int_t i = 0; i<histoDRWithSystErrors->GetNbinsX(); i++){
        // obtain common quantities
        Double_t Rgamma     = histoDRWithSystErrors->GetBinContent(i+binOffset);
        Double_t nIncGamma  = graphIncGammaSysErr->GetY()[i];
        //cout << "bin " << histoDRWithSystErrors->GetBinCenter(i+binOffset) << " " << graphIncGammaSysErr->GetX()[i] << endl;
        //cout << "Rgamma " << Rgamma << " nIncGamma " << nIncGamma << endl;

        // calculating systematics graph
        Double_t errRgamma  = histoDRWithSystErrors->GetBinError(i+binOffset);
        Double_t errNIncGam = graphIncGammaSysErr->GetEYhigh()[i];
        Double_t q1         = 1 - 1/ Rgamma;

        Double_t q1Error    = errRgamma/(Rgamma*Rgamma);
        Double_t content    = nIncGamma * ( 1 - 1/ Rgamma);
        Double_t error      = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        Double_t errDR      = content - error;
        histoDirGammaSpectrumSyst->SetBinError(i+binOffset, error);
        histoDirGammaSpectrumSyst->SetBinContent(i+binOffset, content);
        histoDirGammaSpectrumErrorSyst->SetBinContent(i+binOffset, errDR);

        // calculating stat graphs
        errRgamma           = histoDRWithStatErrors->GetBinError(i+binOffset);
        errNIncGam          = histoIncGamma->GetBinError(i+binOffset);
        q1                  = 1 - 1/ Rgamma;
        q1Error             = errRgamma/(Rgamma*Rgamma);
        content             = nIncGamma * ( 1 - 1/ Rgamma);
        error               = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR               = content - error;
        histoDirGammaSpectrumStat->SetBinError(i+binOffset, error);
        histoDirGammaSpectrumStat->SetBinContent(i+binOffset, content);
        histoDirGammaSpectrumErrorStat->SetBinContent(i+binOffset, errDR);

        // calculating summed error graphs
        errRgamma           = histoDRWithSummedErrors->GetBinError(i+binOffset);
        errNIncGam          = sqrt( pow( histoIncGamma->GetBinError(i+binOffset),2) + pow( graphIncGammaSysErr->GetEYhigh()[i], 2) );
        q1                  = 1 - 1/ Rgamma;
        q1Error             = errRgamma/(Rgamma*Rgamma);
        content             = nIncGamma * ( 1 - 1/ Rgamma);
        error               = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR               = content - error;
        histoDirGammaSpectrumSummed->SetBinError(i+binOffset, error);
        histoDirGammaSpectrumSummed->SetBinContent(i+binOffset, content);
        histoDirGammaSpectrumErrorSummed->SetBinContent(i+binOffset, errDR);

        // calculating sys err without error graph
        errRgamma           = histoDRWithSystErrorsWithoutMat->GetBinError(i+binOffset);
        errNIncGam          = graphIncGammaSysErrW0Mat->GetEYhigh()[i];
        q1                  = 1 - 1/ Rgamma;
        q1Error             = errRgamma/(Rgamma*Rgamma);
        content             = nIncGamma * ( 1 - 1/ Rgamma);
        error               = sqrt( pow( q1 * errNIncGam ,2) + pow( q1Error * nIncGamma ,2));
        errDR               = content - error;
        histoDirGammaSpectrumSystWithoutMat->SetBinError(i+binOffset, error);
        histoDirGammaSpectrumSystWithoutMat->SetBinContent(i+binOffset, content);
        histoDirGammaSpectrumErrorSystWithoutMat->SetBinContent(i+binOffset, errDR);

    }

    // purely calculating points based on all systematic errors
    TGraphAsymmErrors *graphDirGammaSpectrumSyst = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSyst,histoDirGammaSpectrumStat,0,0.5,3);
    if(graphDirGammaSpectrumSyst){
        graphDirGammaSpectrumSyst->SetName("graphDirGammaSpectrumSyst");
        graphDirGammaSpectrumSyst->Print();
    }
    // purely calculating points based on statistical errors
    TGraphAsymmErrors *graphDirGammaSpectrumStat = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorStat,histoDirGammaSpectrumStat,0,0.5,3);
    if(graphDirGammaSpectrumStat){
        graphDirGammaSpectrumStat->SetName("graphDirGammaSpectrumStat");
        graphDirGammaSpectrumStat->Print();
    }
    // purely calculating points based on all systematic + statistical errors
    TGraphAsymmErrors *graphDirGammaSpectrumSummed = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,0,0.5,3);
    if(graphDirGammaSpectrumSummed)graphDirGammaSpectrumSummed->SetName("graphDirGammaSpectrumSummed");
    // calculate points above confidence level summed errors
    TGraphAsymmErrors *graphDirGammaSpectrumSummedConfi = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,2,0.5,3);
    if(graphDirGammaSpectrumSummedConfi)graphDirGammaSpectrumSummedConfi->SetName("graphDirGammaSpectrumSummedConfi");
    // calculate upperlimits summed errors
    TGraphAsymmErrors *graphDirGammaSpectrumSummedUL = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,4,0.5,3);
    if(graphDirGammaSpectrumSummedUL){
        graphDirGammaSpectrumSummedUL->SetName("graphDirGammaSpectrumSummedUL");
        DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumSummedUL , 20, 0, colorCent, colorCent, 1, kTRUE);
        graphDirGammaSpectrumSummedUL->Draw("||");
        PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphDirGammaSpectrumSummedUL);
        graphDirGammaSpectrumSummedUL->Print();
    }
    // calculate arrows for points with 0, error summed
    TGraphAsymmErrors *graphDirGammaSpectrumSummedAr = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,5,0.5,3);
    if(graphDirGammaSpectrumSummedAr){
        graphDirGammaSpectrumSummedAr->SetName("graphDirGammaSpectrumSummedAr");
        DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumSummedAr , 1, 2, colorCent, colorCent, 1, kTRUE);
        graphDirGammaSpectrumSummedAr->Draw("|>");
        PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphDirGammaSpectrumSummedAr);
        graphDirGammaSpectrumSummedAr->Print();
    }
    // calculate points below confidence level summed errors
    TGraphAsymmErrors *graphDirGammaSpectrumSummedULConfi = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,6,0.5,3);
    if(graphDirGammaSpectrumSummedULConfi) graphDirGammaSpectrumSummedULConfi->SetName("graphDirGammaSpectrumSummedULConfi");
    // calculate points below confidence level summed errors with arrows
    TGraphAsymmErrors *graphDirGammaSpectrumSummedArConfi = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSummed,histoDirGammaSpectrumStat,7,0.5,3);
    if(graphDirGammaSpectrumSummedArConfi){
        graphDirGammaSpectrumSummedArConfi->SetName("graphDirGammaSpectrumSummedArConfi");
        DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumSummedArConfi , 1, 2, colorCent, colorCent, 1, kTRUE);
        graphDirGammaSpectrumSummedArConfi->Draw("|>");
        PlotErrorBarAtUpperEdgeOfTGraphAsymErr(graphDirGammaSpectrumSummedArConfi);
        graphDirGammaSpectrumSummedArConfi->Print();
    }
    // purely calculating points based on systematic errors without material budget
    TGraphAsymmErrors *graphDirGammaSpectrumSystwoMat = CalculateDirectPhotonPointsAndUpperLimits(histoDirGammaSpectrumErrorSystWithoutMat,histoDirGammaSpectrumStat,0,0.5,3);
    if(graphDirGammaSpectrumSystwoMat) graphDirGammaSpectrumSystwoMat->SetName("graphDirGammaSpectrumSystWOMat");

    // style setting of direct photon graphs
    DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumSyst , 20,2, colorCent, colorCent, 1, kTRUE);
    graphDirGammaSpectrumSyst->Draw("Z2");
    DrawGammaSetMarkerTGraphAsym(graphDirGammaSpectrumStat , 20, 2, colorCent, colorCent, 1, kTRUE);
    graphDirGammaSpectrumStat->Draw("pE1");

    // **************************************************************************
    // ******************* Plotting direct Gamma Spectrum ***********************
    // **************************************************************************
    TCanvas *canvasDirGamma = GetAndSetCanvas("canvasDirGamma");
    canvasDirGamma->SetLeftMargin(0.14);
    canvasDirGamma->SetLogy();
    canvasDirGamma->SetLogx();

    TH2D *dummyDirGamma = new TH2D("dummyDirGamma", "dummyDirGamma", 1200, 0., 16, 1000., 1e-9,1e4);
    SetStyleHistoTH2ForGraphs( dummyDirGamma, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})",
                               0.045, 0.05, 0.045, 0.05, 0.85, 1.2);
    dummyDirGamma->GetXaxis()->SetLabelOffset(-0.015);
    dummyDirGamma->GetXaxis()->SetRangeUser(minPt, maxPt);
    dummyDirGamma->GetYaxis()->SetRangeUser(1e-8,1e1);
    dummyDirGamma->GetYaxis()->SetRangeUser(1e-6,1e1);
    dummyDirGamma->DrawCopy();

    graphDirGammaSpectrumSyst->Draw("Z2,same");
    graphDirGammaSpectrumStat->Draw("zp,same,e1");

    TLegend* legendDirGamma = GetAndSetLegend(0.18,0.2,1.5);
    legendDirGamma->AddEntry(graphDirGammaSpectrumSyst,Form("#gamma_{dir} PCM, %s",collisionSystem.Data()),"pf");
    legendDirGamma->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_%s.eps",outputDir.Data(),centrality.Data()));

    canvasDirGamma->cd();
    dummyDirGamma->DrawCopy();

        graphDirGammaSpectrumSyst->Draw("Z2,same");
        graphDirGammaSpectrumStat->Draw("zp,same,e1");

        graphPublishedDirGammaSpectrumSyst->Draw("Z2,same");
        graphPublishedDirGammaSpectrumStat->Draw("zp,same,e1");

        TLegend* legendDirGammaWP = GetAndSetLegend(0.18,0.2,1.5);
        legendDirGammaWP->AddEntry(graphDirGammaSpectrumSyst,Form("#gamma_{dir} PCM, %s",collisionSystem.Data()),"pf");
        legendDirGammaWP->AddEntry(graphPublishedDirGammaSpectrumSyst,Form("#gamma_{dir} PCM published, %s",collisionSystem.Data()),"pf");
        legendDirGammaWP->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_%s_withPublished.eps",outputDir.Data(),centrality.Data()));

    canvasDirGamma->cd();
    dummyDirGamma->DrawCopy();

        graphDirGammaSpectrumSyst->Draw("Z2,same");
        graphDirGammaSpectrumStat->Draw("zpE1,same");
        if(graphDirGammaSpectrumSummedUL)graphDirGammaSpectrumSummedUL->Draw("||,same");
        if(graphDirGammaSpectrumSummedAr)graphDirGammaSpectrumSummedAr->Draw("|>,same");

        TLegend* legendDirGammaUL = GetAndSetLegend(0.18,0.2,1.5);
        legendDirGammaUL->AddEntry(graphDirGammaSpectrumSyst,Form("#gamma_{dir} PCM, %s",collisionSystem.Data()),"pf");
        legendDirGammaUL->Draw();

    canvasDirGamma->Print(Form("%s/DirGammaSpectrum_%s_UpperLimits.eps",outputDir.Data(),centrality.Data()));

    // **************************************************************************
    // ****************** Write data to output file *****************************
    // **************************************************************************
    TString optionOutput = optionEnergy;
    const char* fileNameOutputComp = Form("Gamma_PCMResults_%s.root",optionOutput.Data());
    cout << "Creating output file " << fileNameOutputComp << endl;
    TFile* fileGammaSpectrum = new TFile(fileNameOutputComp,"UPDATE");
        fileGammaSpectrum->mkdir(Form("Gamma_%s_%s",optionEnergy.Data(),GetCentralityString(cutSel).Data()));
        fileGammaSpectrum->cd(Form("Gamma_%s_%s", optionEnergy.Data(),GetCentralityString(cutSel).Data()));

            histoDR->SetName("DoubleRatioStatError");
            SetHistogramm(histoDR,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
            histoDR->Write("DoubleRatioStatError",TObject::kOverwrite);
            graphDoubleRatioSysErr->Write("DoubleRatioSystError",TObject::kOverwrite);
            graphDoubleRatioSysErrA->Write("DoubleRatioSystErrorA",TObject::kOverwrite);
            graphDoubleRatioSysErrB->Write("DoubleRatioSystErrorB",TObject::kOverwrite);
            graphDoubleRatioSysErrC->Write("DoubleRatioSystErrorC",TObject::kOverwrite);
            histoDRFit->SetName("DoubleRatioPi0FitStatError");
            SetHistogramm(histoDRFit,"#it{p}_{T} (GeV/#it{c})", "(#it{N}_{#gamma_{inc}}/#it{N}_{#pi^{0}})/(#it{N}_{#gamma_{decay}}/#it{N}_{#pi^{0}})");
            histoDRFit->Write("DoubleRatioPi0FitStatError",TObject::kOverwrite);
            graphDoubleRatioFitPi0SysErr->Write("DoubleRatioPi0FitSystError",TObject::kOverwrite);
            graphDoubleRatioFitPi0SysErrA->Write("DoubleRatioPi0FitSystErrorA",TObject::kOverwrite);
            graphDoubleRatioFitPi0SysErrB->Write("DoubleRatioPi0FitSystErrorB",TObject::kOverwrite);
            graphDoubleRatioFitPi0SysErrC->Write("DoubleRatioPi0FitSystErrorC",TObject::kOverwrite);

            histoIncRatio->SetName("IncRatioStatError");
            SetHistogramm(histoIncRatio,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
            histoIncRatio->Write("IncRatioStatError",TObject::kOverwrite);
            graphIncRatioSysErr->Write("IncRatioSystError",TObject::kOverwrite);
            graphIncRatioSysErrA->Write("IncRatioSystErrorA",TObject::kOverwrite);
            graphIncRatioSysErrB->Write("IncRatioSystErrorB",TObject::kOverwrite);
            graphIncRatioSysErrC->Write("IncRatioSystErrorC",TObject::kOverwrite);
            histoIncRatioPi0Fit->SetName("IncRatioPi0FitStatError");
            SetHistogramm(histoIncRatioPi0Fit,"#it{p}_{T} (GeV/#it{c})", "#gamma_{inc}/#pi^{0}");
            histoIncRatioPi0Fit->Write("IncRatioPi0FitStatError",TObject::kOverwrite);
            graphIncRatioFitPi0SysErr->Write("IncRatioPi0FitSystError",TObject::kOverwrite);
            graphIncRatioFitPi0SysErrA->Write("IncRatioPi0FitSystErrorA",TObject::kOverwrite);
            graphIncRatioFitPi0SysErrB->Write("IncRatioPi0FitSystErrorB",TObject::kOverwrite);
            graphIncRatioFitPi0SysErrC->Write("IncRatioPi0FitSystErrorC",TObject::kOverwrite);

            histoIncGamma->SetName("IncGammaStatError");
            SetHistogramm(histoIncGamma,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#gamma}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
            histoIncGamma->Write("IncGammaStatError",TObject::kOverwrite);
            graphIncGammaSysErr->Write("IncGammaSystError",TObject::kOverwrite);
            graphIncGammaSysErrA->Write("IncGammaSystErrorA",TObject::kOverwrite);
            graphIncGammaSysErrB->Write("IncGammaSystErrorB",TObject::kOverwrite);
            graphIncGammaSysErrC->Write("IncGammaSystErrorC",TObject::kOverwrite);

            histoPi0Spectrum->SetName("Pi0StatError");
            SetHistogramm(histoPi0Spectrum,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
            histoPi0Spectrum->Write("Pi0StatError",TObject::kOverwrite);
            graphPi0SysErr->Write("Pi0SystError",TObject::kOverwrite);
//             histoPi0SpectrumFit->SetName("Pi0FitStatError");
//             SetHistogramm(histoPi0SpectrumFit,"#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}N_{#pi^{0}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV^{-2}#it{c}^{2})");
//             histoPi0SpectrumFit->Write("Pi0FitStatError",TObject::kOverwrite);
//             graphPi0FitSysErr->Write("Pi0FitSystError",TObject::kOverwrite);

            if(graphDirGammaSpectrumSyst)            graphDirGammaSpectrumSyst->Write("graphDirGammaSpectrumSyst",TObject::kOverwrite);
            if(graphDirGammaSpectrumStat)            graphDirGammaSpectrumStat->Write("graphDirGammaSpectrumStat",TObject::kOverwrite);
            if(graphDirGammaSpectrumSummed)            graphDirGammaSpectrumSummed->Write("graphDirGammaSpectrumSummed",TObject::kOverwrite);
            if(graphDirGammaSpectrumSummedConfi)    graphDirGammaSpectrumSummedConfi->Write("graphDirGammaSpectrumSummedConfi",TObject::kOverwrite);
            if(graphDirGammaSpectrumSummedUL)        graphDirGammaSpectrumSummedUL->Write("graphDirGammaSpectrumSummedUL",TObject::kOverwrite);
            if(graphDirGammaSpectrumSummedAr)        graphDirGammaSpectrumSummedAr->Write("graphDirGammaSpectrumSummedAr",TObject::kOverwrite);
            if(graphDirGammaSpectrumSummedULConfi)    graphDirGammaSpectrumSummedULConfi->Write("graphDirGammaSpectrumSummedULConfi",TObject::kOverwrite);
            if(graphDirGammaSpectrumSummedArConfi)    graphDirGammaSpectrumSummedArConfi->Write("graphDirGammaSpectrumSummedAr",TObject::kOverwrite);
            if(graphDirGammaSpectrumSystwoMat)         graphDirGammaSpectrumSystwoMat->Write("graphDirGammaSpectrumSystWOMat",TObject::kOverwrite);

            fileGammaSpectrum->mkdir("Theory");
            fileGammaSpectrum->cd("Theory");
                if (optionOutput.CompareTo("pp") == 0 ) NLODoubleRatio->Write(Form("NLODoubleRatio_%s_%s",optionEnergy.Data(),GetCentralityString(cutSel).Data()),TObject::kOverwrite);
                    else NLODoubleRatio->Write(Form("NLODoubleRatio_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
                if (optionOutput.CompareTo("pp") == 0 ) NLO->Write(Form("NLO_%s_%s",optionEnergy.Data(),GetCentralityString(cutSel).Data()),TObject::kOverwrite);
                    else NLO->Write(Form("NLO_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
                if (graphInvYieldPbPbTheoryCTEQ61EPS09) graphInvYieldPbPbTheoryCTEQ61EPS09->Write(Form("EPS09_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
                if (graphInvYieldPPTheoryCT10BFG2_pdfErr) graphInvYieldPPTheoryCT10BFG2_pdfErr->Write(Form("CT10BF_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
                if (graphDRPbPbCTEQ61EPS09) graphDRPbPbCTEQ61EPS09->Write(Form("EPS09DoubleRatio_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
                if (graphDRPbPbCT10BFG2_pdfErr) graphDRPbPbCT10BFG2_pdfErr->Write(Form("CT10BFDoubleRatio_%s",GetCentralityString(cutSel).Data()),TObject::kOverwrite);
    fileGammaSpectrum->Write();
    fileGammaSpectrum->Close();


}

