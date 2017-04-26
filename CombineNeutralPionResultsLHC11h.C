// ***************************************************************/
// ******      Friederike Bock, friederike.bock@cern.ch       ****
// ******      Lucia Leardini, lucia.leardini@cern.ch 		  ****
// ***************************************************************/
#include <Riostream.h>
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
#include "TGraphAsymmErrors.h" 
#include "TGaxis.h"
#include "TMarker.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CombineMesonMeasurementsPbPbLHC11hV1.h"


extern TRandom*   gRandom;
extern TBenchmark*   gBenchmark;
extern TSystem*   gSystem;
extern TMinuit*   gMinuit;

void ScaleMCYield(TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson, TString optionDalitz ){
   histoCorrectedToBeScaled->Sumw2();
   histoCorrectedToBeScaled->Scale(1./deltaRapid);
   histoCorrectedToBeScaled->Scale(scaling);
   histoCorrectedToBeScaled->Scale(1./nEvtMC);
   for (Int_t i = 1; i < histoCorrectedToBeScaled->GetNbinsX()+1 ; i++){
      Double_t newBinContent = histoCorrectedToBeScaled->GetBinContent(i)/histoCorrectedToBeScaled->GetBinCenter(i);
      Double_t newBinError = histoCorrectedToBeScaled->GetBinError(i)/histoCorrectedToBeScaled->GetBinCenter(i);
      histoCorrectedToBeScaled->SetBinContent(i,newBinContent);
      histoCorrectedToBeScaled->SetBinError(i,newBinError);
   }
   if (nameMeson.CompareTo("Pi0") == 0 ||nameMeson.CompareTo("Pi0EtaBinning") == 0 ){
      if (optionDalitz.CompareTo("kFALSE")==0){
         histoCorrectedToBeScaled->Scale(1./0.98798);
      } else {
         histoCorrectedToBeScaled->Scale(1./0.01198);
      }
   }else{
      if (optionDalitz.CompareTo("kFALSE")==0){
         histoCorrectedToBeScaled->Scale(1./0.3931);
      } else {
         histoCorrectedToBeScaled->Scale(1./6.8e-5);
      }
      
   }
}


void CombineNeutralPionResultsLHC11h(TString suffix = "pdf",
                                     TString nameFilePbPbLHC11h = "data_PCMResults_PbPb_2.76TeV",
                                     Bool_t runDrawReweighted = kTRUE,
                                     Bool_t runPPplotting = kFALSE,
                                     TString thisthesis="ALICE work in progress" //thisthesis.Data()
//                                      Bool_t thesisPlotting = kFALSE

){

	gROOT->Reset();   
	gROOT->SetStyle("Plain");
	
	TString dateForOutput = ReturnDateStringForOutput();

    //file from soon to be published paper, Fredi email from 4th Oct 2016
    //saved in ExternalInputPbPb/NeutralMesonpp2760GeVReference folder
    TString nameFilePP = "ExternalInputPbPb/NeutralMesonspp276GeVReference/CombinedResultsPaperPP2760GeV_2016_10_03_FrediV2Clusterizer.root";
	TString nameFilePi0PbPbLHC10h = "ExternalInputPbPb/data_PCMResults_PbPb_2.76TeV_LHC10h.root";
    TString nameFilePi0PbPbPaper = "FinalResults/CombinedResultsPbPb_ShiftedX_PaperRAA_13_Aug_2014.root";
    TString nameFileChargedPionKaonSpectra = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Kaon_08052014.root";
    TString nameFileChargedPionKaonRAA = "ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Pion_08052014.root";
    TString nameFileOtherExp = "ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesPbPbWithEta.root";
    TString nameFileCocktailInput = "/home/admin1/leardini/newSoftware/cocktail_input/CocktailInputPbPb_Param_PbPb_2.76TeV_20April2017.root";

	TString outputDir = Form("%s/%s/NeutralMesonsPCMResults",suffix.Data(),dateForOutput.Data());
	gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec(Form("cp %s %s/InputFilePionPPPaper.root ",nameFilePP.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFilePCMPionPbPbLHC10h.root ",nameFilePi0PbPbLHC10h.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFilePionPbPbPaper.root ",nameFilePi0PbPbPaper.Data(),outputDir.Data() ));
	gSystem->Exec(Form("cp %s %s/InputFilePCMPbPbLHC11h.root ",nameFilePbPbLHC11h.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileChargedSpectra.root ",nameFileChargedPionKaonSpectra.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileChargedRaa.root ",nameFileChargedPionKaonRAA.Data(),outputDir.Data() ));
    gSystem->Exec(Form("cp %s %s/InputFileOtherExp.root ",nameFileOtherExp.Data(),outputDir.Data() ));
   
	StyleSettingsThesis();  
	SetPlotStyle();

    TString cent0010 = "0#font[122]{-}10%";
    TString cent2050 = "20#font[122]{-}50%";
    TString cent2040 = "20#font[122]{-}40%";

	Color_t  colorCombPP             = kBlack;
	Color_t  colorCombPbPb0005          = kRed+1;
	Color_t  colorCombPbPb0010          = kRed+1;
	Color_t  colorCombPbPb0510          = 807;
	Color_t  colorCombPbPb1020          = 800;
	Color_t  colorCombPbPb2040          = kAzure+1;//-5;//kGreen+2;
	Color_t  colorCombPbPb4060          = kCyan+2;
	Color_t  colorCombPbPb4050          = kTeal+2;
	Color_t  colorCombPbPb5060          = kAzure+5;
	Color_t  colorCombPbPb6070          = kAzure-5;
	Color_t  colorCombPbPb6080          = kBlue+1;
	Color_t  colorCombPbPb7080          = kBlue-1;
	Color_t  colorCombPbPb8090          = kViolet+3;
	Color_t  colorCombPbPb7590          = kViolet-3;
	
	Color_t  colorCombMCPbPb0005           = kRed+3;
	Color_t  colorCombMCPbPb0010           = kRed+3;
	Color_t  colorCombMCPbPb0510           = 807+2;
	Color_t  colorCombMCPbPb1020           = 800+2;
	Color_t  colorCombMCPbPb2040           = kGreen+4;
	Color_t  colorCombMCPbPb4060           = kCyan+4;
	Color_t  colorCombMCPbPb6080           = kBlue+3;

	Style_t  markerStylePP     = 33 ;
	Style_t  markerStylePbPb0005  = 20 ;
	Style_t  markerStylePbPb0010  = 20 ;
	Style_t  markerStylePbPb0510  = 21 ;
	Style_t  markerStylePbPb1020  = 29 ;
	Style_t  markerStylePbPb2040  = 21 ;
	Style_t  markerStylePbPb4050  = 24 ;
	Style_t  markerStylePbPb4060  = 20 ;
	Style_t  markerStylePbPb5060  = 24 ;
	Style_t  markerStylePbPb6070  = 25 ;
	Style_t  markerStylePbPb6080  = 21 ;
	Style_t  markerStylePbPb7080  = 25 ;
	Style_t  markerStylePbPb8090  = 27 ;
	Style_t  markerStylePbPb7590  = 30 ;
	
	Style_t  markerStylePbPb0005MC   = 24 ;
	Style_t  markerStylePbPb0010MC   = 24 ;
	Style_t  markerStylePbPb0510MC   = 25 ;
	Style_t  markerStylePbPb1020MC   = 30 ;
	Style_t  markerStylePbPb2040MC   = 27 ;
	Style_t  markerStylePbPb4060MC   = 24 ;
	Style_t  markerStylePbPb6080MC   = 25 ;
	
	Size_t   markerSizePP      = 2.5;
	Size_t   markerSizePbPb0005   = 2.;
	Size_t   markerSizePbPb0010   = 2.;
	Size_t   markerSizePbPb0510   = 2.;
	Size_t   markerSizePbPb1020   = 2.5;
	Size_t   markerSizePbPb2040   = 2.;
	Size_t   markerSizePbPb4060   = 2.;
	Size_t   markerSizePbPb6080   = 2.;
	
	Color_t  colorPi0900GeV          = kRed +2;
	Color_t  colorPi02760GeV         = kMagenta+2;
	Color_t  colorPi07TeV            = kBlue+2;
	Color_t  colorPi0900GeVBox = colorPi0900GeV-10;
	Color_t  colorPi02760GeVBox = colorPi02760GeV-10;
	Color_t  colorPi07TeVBox = colorPi07TeV-10;

	Color_t  colorMCPythiaPP900GeV   = colorPi0900GeV-4;
	Color_t  colorMCPythiaPP2760GeV = colorPi02760GeV+2;
	Color_t  colorMCPythiaPP7TeV  = colorPi07TeV+3;
	Color_t  colorMCPhojetPP900GeV   = colorPi0900GeV+2;
	Color_t  colorMCPhojetPP2760GeV = colorPi02760GeV-4;
	Color_t  colorMCPhojetPP7TeV  = colorPi07TeV-3;

	Style_t  markerStyleSpectrum7TeVMC  = 24 ;
	Style_t  markerStyleSpectrum900GeVMC = 25 ;
	Style_t  markerStyleSpectrum2760GeVMC = 30 ;
	Style_t  markerStyleSpectrum7TeV    = 20 ;
	Style_t  markerStyleSpectrum900GeV = 21 ;
	Style_t  markerStyleSpectrum2760GeV = 29 ;
	
	Double_t xSection7TeVppINEL = 73.2*1e9;
	Double_t xSection2760GeVppINEL = 62.8*1e9;
	Double_t xSection900GeVppINEL = 52.5*1e9;
	
	Style_t  markerStyleMCPP7TeV  = 24 ;
	Style_t  markerStyleMCPP900GeV   = 25 ;
	Style_t  markerStyleMCPP2760GeV  = 30 ;
	
	Size_t   markerSizePi0PP7TeV  = 1.8;
	Size_t   markerSizePi0PP900GeV = 1.8;
	Size_t   markerSizePi0PP2760GeV  = 2.2;
	
    Color_t colorSec[4]                             = {kRed+2, 807, kCyan+2, kBlue};
    Style_t markerStyleSec[4]                       = {21, 33, 29, 34};
    Size_t markerSizeSec[4]                         = {1.5, 1.75, 2., 1.5};
    Color_t colorSecFromToy[3]                      = {kRed-2, kOrange+1, kCyan-2};
    Style_t markerStyleSecFromToy[3]                = {20, 33, 29};
    Size_t markerSizeSecFromToy[3]                  = {1.7, 2, 2.2};

	TString collisionSystemPP2760GeV = "pp #sqrt{#it{s}} = 2.76 TeV";      
	TString collisionSystemPP7TeV = "pp #sqrt{#it{s}} = 7 TeV";      
	TString collisionSystemPP900GeV = "pp #sqrt{#it{s}} = 0.9 TeV";     
    TString collisionSystemPbPb = "Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb0005 = "0#font[122]{-}5% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb0010 = "0#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb0510 = "5#font[122]{-}10% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb1020 = "10#font[122]{-}20% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb2030 = "20#font[122]{-}30% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb3040 = "30#font[122]{-}40% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb2040 = "20#font[122]{-}40% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb2050 = "20#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb3050 = "30#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb4060 = "40#font[122]{-}60% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb6080 = "60#font[122]{-}80% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb6070 = "60#font[122]{-}70% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb7080 = "70#font[122]{-}80% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb7590 = "75#font[122]{-}90% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb8090 = "80#font[122]{-}90% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb4050 = "40#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb5060 = "50#font[122]{-}60% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb0020 = "0#font[122]{-}20% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb0080 = "0#font[122]{-}80% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb0040 = "0#font[122]{-}40% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	TString collisionSystemPbPb4080 = "40#font[122]{-}80% Pb#font[122]{-}Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";
	
	Size_t markerSizeComparison = 0.5;
	Double_t maxPtMesonEffFit = 12.;
	Double_t minPtMesonEffFit = 1.2;
	Int_t offsetCorrectionHighPt= 1;
	TF1* fitTrueEffi = new TF1("EffiFitDummy","1 - [0]*exp([1]*x)+[1]");
	fitTrueEffi->SetRange(minPtMesonEffFit,maxPtMesonEffFit);
   
    Bool_t thesisPlotting;
    if(thisthesis.CompareTo("")==0) thesisPlotting = kFALSE;
    else thesisPlotting = kTRUE;

    Double_t textSize = 0.035;
    TLatex *thesisLabelTopLeft = new TLatex(0.13,0.9,thisthesis.Data());
    SetStyleTLatex( thesisLabelTopLeft, textSize,4);

    TLatex *thesisLabel = new TLatex(0.6,0.93,thisthesis.Data());
    SetStyleTLatex( thesisLabel, textSize,4);

    TLatex *thesisLabel2 = new TLatex(0.8,0.15,thisthesis.Data());
    SetStyleTLatex( thesisLabel2, textSize,4);

    TLatex *thesisLabel3 = new TLatex(0.12,0.9,thisthesis.Data());
    SetStyleTLatex( thesisLabel3, textSize,4);

	TLatex *labelRawPi0PbPb = new TLatex(0.6,0.9,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelRawPi0PbPb, textSize,4);

	TLatex *labelRawEtaPbPb = new TLatex(0.6,0.9,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelRawEtaPbPb,textSize,4);

	TLatex *labelMassEtaPP = new TLatex(0.2,0.9,collisionSystemPP7TeV.Data());
	SetStyleTLatex( labelMassEtaPP, 0.062,4);

	TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
	SetStyleTLatex( labelLegendAMass, 0.08,4);

    TH1D*    histoSecAcceptance[5][4] = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL} }; 
    TH1D*    histoSecTrueEffi[5][4] = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL} }; 
    TH1D*    histoYieldSecMeson[5][4] = { { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL}, { NULL, NULL, NULL, NULL} }; 
    TString nameSecMeson[4]                         = {"K0S", "Lambda", "K0L", "Rest"};
    TString nameSecMesonPlot[4]                     = {"K_{s}^{0}", "#Lambda", "K_{l}^{0}", "Rest"};


    //PbPb PCM only - 2011 data
    TFile*   filePCMPbPbLHC11h = new TFile(nameFilePbPbLHC11h.Data());
    cout << "Pi0 0-10%" << endl;
    TDirectory* directoryPi0PbPbLHC11h0010 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_0-10%");
        TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h0010                 = (TH1D*)directoryPi0PbPbLHC11h0010->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h0010 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0010->Get("Pi0SystError");

        TH1D* histoPCMPi0MassDataPbPbLHC11h0010         = (TH1D*)directoryPi0PbPbLHC11h0010->Get("MassPi0");
        histoPCMPi0MassDataPbPbLHC11h0010->Scale(1000.);
        TH1D* histoPCMPi0MassMCPbPbLHC11h0010           = (TH1D*)directoryPi0PbPbLHC11h0010->Get("TrueMassPi0");
        histoPCMPi0MassMCPbPbLHC11h0010->Scale(1000.);

        for (Int_t j = 0; j < 4; j++){        
            histoSecAcceptance[2][j] 	= (TH1D*)directoryPi0PbPbLHC11h0010->Get(Form("Pi0From%sAcceptance",nameSecMeson[j].Data()));
            histoSecTrueEffi[2][j] 			= (TH1D*)directoryPi0PbPbLHC11h0010->Get(Form("Pi0From%sEfficiency",nameSecMeson[j].Data()));
            histoYieldSecMeson[2][j] 	    = (TH1D*)directoryPi0PbPbLHC11h0010->Get(Form("Pi0From%sRawYieldPerEvent",nameSecMeson[j].Data()));
        }
        
        TH1D* histoPCMPi0WidthDataPbPbLHC11h0010        = (TH1D*)directoryPi0PbPbLHC11h0010->Get("FWHMPi0MeV");
        TH1D* histoPCMPi0WidthMCPbPbLHC11h0010          = (TH1D*)directoryPi0PbPbLHC11h0010->Get("TrueFWHMPi0MeV");
        TH1D* histoPi0TrueEffiPtPbPbLHC11h0010          = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_Efficiency");
        TH1D* histoPi0AcceptPtPbPbLHC11h0010            = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_Acceptance");
        TH1D* histoMCYieldPi0PtPbPbLHC11h0010           = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input_Reweighted");
        TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights0010  = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input");
        //    TH1D* histoPi0WeightsPbPbLHC11h0010       = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Weights");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig0010   = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0010 = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input_AddedSig");
        TH1D* histoPi0RawYieldPbPbLHC11h0010            = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_RawYieldPerEvent");

        TGraphAsymmErrors* graphPi0RAA0010              = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0010->Get("Pi0RAA");
        TGraphAsymmErrors* graphPi0RAASys0010           = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0010->Get("Pi0RAASys");


    cout << "0-5%" << endl;
    TDirectory* directoryPi0PbPbLHC11h0005 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_0-5%");

        TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h0005                 = (TH1D*)directoryPi0PbPbLHC11h0005->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h0005 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0005->Get("Pi0SystError");

        TH1D* histoPCMPi0MassDataPbPbLHC11h0005         = (TH1D*)directoryPi0PbPbLHC11h0005->Get("MassPi0");
        histoPCMPi0MassDataPbPbLHC11h0005->Scale(1000.);
        TH1D* histoPCMPi0MassMCPbPbLHC11h0005           = (TH1D*)directoryPi0PbPbLHC11h0005->Get("TrueMassPi0");
        histoPCMPi0MassMCPbPbLHC11h0005->Scale(1000.);
        TH1D* histoPCMPi0WidthDataPbPbLHC11h0005        = (TH1D*)directoryPi0PbPbLHC11h0005->Get("FWHMPi0MeV");
        TH1D* histoPCMPi0WidthMCPbPbLHC11h0005          = (TH1D*)directoryPi0PbPbLHC11h0005->Get("TrueFWHMPi0MeV");
        TH1D* histoPi0TrueEffiPtPbPbLHC11h0005          = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_Efficiency");
        TH1D* histoPi0AcceptPtPbPbLHC11h0005            = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_Acceptance");
        TH1D* histoMCYieldPi0PtPbPbLHC11h0005           = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input_Reweighted");
        TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights0005  = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input");
        //    TH1D* histoPi0WeightsPbPbLHC11h0005       = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Weights");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig0005   = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0005 = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input_AddedSig");
        TH1D* histoPi0RawYieldPbPbLHC11h0005            = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_RawYieldPerEvent");

        TGraphAsymmErrors* graphPi0RAA0005              = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0005->Get("Pi0RAA");
        TGraphAsymmErrors* graphPi0RAASys0005           = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0005->Get("Pi0RAASys");

    cout << "5-10%" << endl;
    TDirectory* directoryPi0PbPbLHC11h0510 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_5-10%");

        TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h0510                 = (TH1D*)directoryPi0PbPbLHC11h0510->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h0510 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0510->Get("Pi0SystError");

        TH1D* histoPCMPi0MassDataPbPbLHC11h0510         = (TH1D*)directoryPi0PbPbLHC11h0510->Get("MassPi0");
        histoPCMPi0MassDataPbPbLHC11h0510->Scale(1000.);
        TH1D* histoPCMPi0MassMCPbPbLHC11h0510           = (TH1D*)directoryPi0PbPbLHC11h0510->Get("TrueMassPi0");
        histoPCMPi0MassMCPbPbLHC11h0510->Scale(1000.);

        TH1D* histoPCMPi0WidthDataPbPbLHC11h0510        = (TH1D*)directoryPi0PbPbLHC11h0510->Get("FWHMPi0MeV");
        TH1D* histoPCMPi0WidthMCPbPbLHC11h0510          = (TH1D*)directoryPi0PbPbLHC11h0510->Get("TrueFWHMPi0MeV");
        TH1D* histoPi0TrueEffiPtPbPbLHC11h0510          = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_Efficiency");
        TH1D* histoPi0AcceptPtPbPbLHC11h0510            = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_Acceptance");
        TH1D* histoMCYieldPi0PtPbPbLHC11h0510           = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input_Reweighted");
        TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights0510  = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input");
        //    TH1D* histoPi0WeightsPbPbLHC11h0510       = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Weights");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig0510   = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0510 = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input_AddedSig");
        TH1D* histoPi0RawYieldPbPbLHC11h0510            = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_RawYieldPerEvent");

    cout << "20-40%" << endl;
    TDirectory* directoryPi0PbPbLHC11h2040 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_20-40%");

        TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h2040                 = (TH1D*)directoryPi0PbPbLHC11h2040->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h2040 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2040->Get("Pi0SystError");

        TH1D* histoPCMPi0MassDataPbPbLHC11h2040         = (TH1D*)directoryPi0PbPbLHC11h2040->Get("MassPi0");
        histoPCMPi0MassDataPbPbLHC11h2040->Scale(1000.);
        TH1D* histoPCMPi0MassMCPbPbLHC11h2040           = (TH1D*)directoryPi0PbPbLHC11h2040->Get("TrueMassPi0");
        histoPCMPi0MassDataPbPbLHC11h2040->Scale(1000.);

        TH1D* histoPCMPi0WidthDataPbPbLHC11h2040        = (TH1D*)directoryPi0PbPbLHC11h2040->Get("FWHMPi0MeV");
        TH1D* histoPCMPi0WidthMCPbPbLHC11h2040          = (TH1D*)directoryPi0PbPbLHC11h2040->Get("TrueFWHMPi0MeV");
        TH1D* histoPi0TrueEffiPtPbPbLHC11h2040          = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_Efficiency");
        TH1D* histoPi0AcceptPtPbPbLHC11h2040            = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_Acceptance");
        TH1D* histoMCYieldPi0PtPbPbLHC11h2040           = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input_Reweighted");
        TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights2040  = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input");
        //    TH1D* histoPi0WeightsPbPbLHC11h2040       = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Weights");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig2040   = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2040 = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input_AddedSig");
        TH1D* histoPi0RawYieldPbPbLHC11h2040            = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_RawYieldPerEvent");

        TGraphAsymmErrors* graphPi0RAA2040              = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2040->Get("Pi0RAA");
        TGraphAsymmErrors* graphPi0RAASys2040           = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2040->Get("Pi0RAASys");

    cout << "20-50%" << endl;
    TDirectory* directoryPi0PbPbLHC11h2050 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_20-50%");

        TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h2050                 = (TH1D*)directoryPi0PbPbLHC11h2050->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h2050 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2050->Get("Pi0SystError");

        TH1D* histoPCMPi0MassDataPbPbLHC11h2050         = (TH1D*)directoryPi0PbPbLHC11h2050->Get("MassPi0");
        histoPCMPi0MassDataPbPbLHC11h2050->Scale(1000.);
        TH1D* histoPCMPi0MassMCPbPbLHC11h2050           = (TH1D*)directoryPi0PbPbLHC11h2050->Get("TrueMassPi0");
        histoPCMPi0MassMCPbPbLHC11h2050->Scale(1000.);

        for (Int_t j = 0; j < 4; j++){        
            histoSecAcceptance[4][j] 	    	= (TH1D*)directoryPi0PbPbLHC11h0010->Get(Form("Pi0From%sAcceptance",nameSecMeson[j].Data()));
            histoSecTrueEffi[4][j] 			= (TH1D*)directoryPi0PbPbLHC11h0010->Get(Form("Pi0From%sEfficiency",nameSecMeson[j].Data()));
            histoYieldSecMeson[4][j] 	    = (TH1D*)directoryPi0PbPbLHC11h0010->Get(Form("Pi0From%sRawYieldPerEvent",nameSecMeson[j].Data()));
        }

        TH1D* histoPCMPi0WidthDataPbPbLHC11h2050        = (TH1D*)directoryPi0PbPbLHC11h2050->Get("FWHMPi0MeV");
        TH1D* histoPCMPi0WidthMCPbPbLHC11h2050          = (TH1D*)directoryPi0PbPbLHC11h2050->Get("TrueFWHMPi0MeV");
        TH1D* histoPi0TrueEffiPtPbPbLHC11h2050          = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_Efficiency");
        TH1D* histoPi0AcceptPtPbPbLHC11h2050            = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_Acceptance");
        TH1D* histoMCYieldPi0PtPbPbLHC11h2050           = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input_Reweighted");
        TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights2050  = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input");
        //    TH1D* histoPi0WeightsPbPbLHC11h2050       = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Weights");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig2050   = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2050 = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input_AddedSig");
        TH1D* histoPi0RawYieldPbPbLHC11h2050            = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_RawYieldPerEvent");

        TGraphAsymmErrors* graphPi0RAA2050              = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2050->Get("Pi0RAA");
        TGraphAsymmErrors* graphPi0RAASys2050           = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2050->Get("Pi0RAASys");


    cout << "Eta 0-10%" << endl;
    TDirectory* directoryEtaPbPbLHC11h0010 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_0-10%");

        TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h0010                 = (TH1D*)directoryEtaPbPbLHC11h0010->Get("CorrectedYieldEta");
        TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h0010 =        (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0010->Get("EtaSystError");

        TH1D* histoPCMEtaMassDataPbPbLHC11h0010         = (TH1D*)directoryEtaPbPbLHC11h0010->Get("MassEta");
        histoPCMEtaMassDataPbPbLHC11h0010->Scale(1000.);
        TH1D* histoPCMEtaMassMCPbPbLHC11h0010           = (TH1D*)directoryEtaPbPbLHC11h0010->Get("TrueMassEta");
        histoPCMEtaMassMCPbPbLHC11h0010->Scale(1000.);

        TH1D* histoPCMEtaWidthDataPbPbLHC11h0010        = (TH1D*)directoryEtaPbPbLHC11h0010->Get("FWHMEtaMeV");
        TH1D* histoPCMEtaWidthMCPbPbLHC11h0010          = (TH1D*)directoryEtaPbPbLHC11h0010->Get("TrueFWHMEtaMeV");
        TH1D* histoEtaTrueEffiPtPbPbLHC11h0010          = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_Efficiency");
        TH1D* histoEtaAcceptPtPbPbLHC11h0010            = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_Acceptance");
        TH1D* histoMCYieldEtaPtPbPbLHC11h0010           = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input_Reweighted");
        TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights0010  = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input");
        //    TH1D* histoEtaWeightsPbPbLHC11h0010       = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Weights");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig0010   = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0010 = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input_AddedSig");
        TH1D* histoEtaRawYieldPbPbLHC11h0010            = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_RawYieldPerEvent");

        TGraphAsymmErrors* graphEtaRAA0010              = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0010->Get("EtaRAA");
        TGraphAsymmErrors* graphEtaRAASys0010           = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0010->Get("EtaRAASys");

        TH1D* histoPCMEtaToPi0RatioPbPb0010 =            (TH1D*)directoryEtaPbPbLHC11h0010->Get("EtatoPi0Ratio");
        TGraphAsymmErrors *graphPCMEtaToPi0RatioPbPb0010 = new TGraphAsymmErrors(histoPCMEtaToPi0RatioPbPb0010);
        TGraphAsymmErrors* graphPCMEtaToPi0RatioSysErrPbPb0010=    (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0010->Get("EtatoPi0RatioSys");

    cout << "0-5%" << endl;
    TDirectory* directoryEtaPbPbLHC11h0005 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_0-5%");

        TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h0005                 = (TH1D*)directoryEtaPbPbLHC11h0005->Get("CorrectedYieldEta");
        TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h0005 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0005->Get("EtaSystError");

        TH1D* histoPCMEtaMassDataPbPbLHC11h0005         = (TH1D*)directoryEtaPbPbLHC11h0005->Get("MassEta");
        histoPCMEtaMassDataPbPbLHC11h0005->Scale(1000.);
        TH1D* histoPCMEtaMassMCPbPbLHC11h0005           = (TH1D*)directoryEtaPbPbLHC11h0005->Get("TrueMassEta");
        histoPCMEtaMassMCPbPbLHC11h0005->Scale(1000.);

        TH1D* histoPCMEtaWidthDataPbPbLHC11h0005        = (TH1D*)directoryEtaPbPbLHC11h0005->Get("FWHMEtaMeV");
        TH1D* histoPCMEtaWidthMCPbPbLHC11h0005          = (TH1D*)directoryEtaPbPbLHC11h0005->Get("TrueFWHMEtaMeV");
        TH1D* histoEtaTrueEffiPtPbPbLHC11h0005          = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_Efficiency");
        TH1D* histoEtaAcceptPtPbPbLHC11h0005            = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_Acceptance");
        TH1D* histoMCYieldEtaPtPbPbLHC11h0005           = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input_Reweighted");
        TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights0005  = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input");
        //    TH1D* histoEtaWeightsPbPbLHC11h0005 = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Weights");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig0005   = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0005 = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input_AddedSig");
        TH1D* histoEtaRawYieldPbPbLHC11h0005            = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_RawYieldPerEvent");

        TGraphAsymmErrors* graphEtaRAA0005              = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0005->Get("EtaRAA");
        TGraphAsymmErrors* graphEtaRAASys0005           = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0005->Get("EtaRAASys");

    cout << "5-10%" << endl;
    TDirectory* directoryEtaPbPbLHC11h0510 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_5-10%");

        TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h0510                 = (TH1D*)directoryEtaPbPbLHC11h0510->Get("CorrectedYieldEta");
        TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h0510 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0510->Get("EtaSystError");

        TH1D* histoPCMEtaMassDataPbPbLHC11h0510         = (TH1D*)directoryEtaPbPbLHC11h0510->Get("MassEta");
        histoPCMEtaMassDataPbPbLHC11h0510->Scale(1000.);
        TH1D* histoPCMEtaMassMCPbPbLHC11h0510           = (TH1D*)directoryEtaPbPbLHC11h0510->Get("TrueMassEta");
        histoPCMEtaMassMCPbPbLHC11h0510->Scale(1000.);

        TH1D* histoPCMEtaWidthDataPbPbLHC11h0510        = (TH1D*)directoryEtaPbPbLHC11h0510->Get("FWHMEtaMeV");
        TH1D* histoPCMEtaWidthMCPbPbLHC11h0510          = (TH1D*)directoryEtaPbPbLHC11h0510->Get("TrueFWHMEtaMeV");
        TH1D* histoEtaTrueEffiPtPbPbLHC11h0510          = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_Efficiency");
        TH1D* histoEtaAcceptPtPbPbLHC11h0510            = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_Acceptance");
        TH1D* histoMCYieldEtaPtPbPbLHC11h0510           = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input_Reweighted");
        TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights0510  = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input");
        //    TH1D* histoEtaWeightsPbPbLHC11h0510       = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Weights");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig0510   = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0510 = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input_AddedSig");
        TH1D* histoEtaRawYieldPbPbLHC11h0510            = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_RawYieldPerEvent");

    cout << "20-40%" << endl;
    TDirectory* directoryEtaPbPbLHC11h2040 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_20-40%");

        TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h2040                 = (TH1D*)directoryEtaPbPbLHC11h2040->Get("CorrectedYieldEta");
        TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h2040 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2040->Get("EtaSystError");

        TH1D* histoPCMEtaMassDataPbPbLHC11h2040         = (TH1D*)directoryEtaPbPbLHC11h2040->Get("MassEta");
        histoPCMEtaMassDataPbPbLHC11h2040->Scale(1000.);
        TH1D* histoPCMEtaMassMCPbPbLHC11h2040           = (TH1D*)directoryEtaPbPbLHC11h2040->Get("TrueMassEta");
        histoPCMEtaMassMCPbPbLHC11h2040->Scale(1000.);

        TH1D* histoPCMEtaWidthDataPbPbLHC11h2040        = (TH1D*)directoryEtaPbPbLHC11h2040->Get("FWHMEtaMeV");
        TH1D* histoPCMEtaWidthMCPbPbLHC11h2040          = (TH1D*)directoryEtaPbPbLHC11h2040->Get("TrueFWHMEtaMeV");
        TH1D* histoEtaTrueEffiPtPbPbLHC11h2040          = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_Efficiency");
        TH1D* histoEtaAcceptPtPbPbLHC11h2040            = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_Acceptance");
        TH1D* histoMCYieldEtaPtPbPbLHC11h2040           = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input_Reweighted");
        TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights2040  = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input");
        //    TH1D* histoEtaWeightsPbPbLHC11h0010       = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Weights");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig2040   = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2040 = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input_AddedSig");
        TH1D* histoEtaRawYieldPbPbLHC11h2040            = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_RawYieldPerEvent");

        TGraphAsymmErrors* graphEtaRAA2040              = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2040->Get("EtaRAA");
        TGraphAsymmErrors* graphEtaRAASys2040           = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2040->Get("EtaRAASys");

        TH1D* histoPCMEtaToPi0RatioPbPb2040 =            (TH1D*)directoryEtaPbPbLHC11h2040->Get("EtatoPi0Ratio");
        TGraphAsymmErrors *graphPCMEtaToPi0RatioPbPb2040 = new TGraphAsymmErrors(histoPCMEtaToPi0RatioPbPb2040);
        TGraphAsymmErrors* graphPCMEtaToPi0RatioSysErrPbPb2040=    (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2040->Get("EtatoPi0RatioSys");

    cout << "20-50%" << endl;
    TDirectory* directoryEtaPbPbLHC11h2050 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_20-50%");

        TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h2050                 = (TH1D*)directoryEtaPbPbLHC11h2050->Get("CorrectedYieldEta");
        TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h2050 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2050->Get("EtaSystError");

        TH1D* histoPCMEtaMassDataPbPbLHC11h2050         = (TH1D*)directoryEtaPbPbLHC11h2050->Get("MassEta");
        histoPCMEtaMassDataPbPbLHC11h2050->Scale(1000.);
        TH1D* histoPCMEtaMassMCPbPbLHC11h2050           = (TH1D*)directoryEtaPbPbLHC11h2050->Get("TrueMassEta");
        histoPCMEtaMassMCPbPbLHC11h2050->Scale(1000.);

        TH1D* histoPCMEtaWidthDataPbPbLHC11h2050        = (TH1D*)directoryEtaPbPbLHC11h2050->Get("FWHMEtaMeV");
        TH1D* histoPCMEtaWidthMCPbPbLHC11h2050          = (TH1D*)directoryEtaPbPbLHC11h2050->Get("TrueFWHMEtaMeV");
        TH1D* histoEtaTrueEffiPtPbPbLHC11h2050          = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_Efficiency");
        TH1D* histoEtaAcceptPtPbPbLHC11h2050            = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_Acceptance");
        TH1D* histoMCYieldEtaPtPbPbLHC11h2050           = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input_Reweighted");
        TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights2050  = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig2050   = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input_Reweighted_AddedSig");
        TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2050 = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input_AddedSig");
        TH1D* histoEtaRawYieldPbPbLHC11h2050            = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_RawYieldPerEvent");

        TGraphAsymmErrors* graphEtaRAA2050              = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2050->Get("EtaRAA");
        TGraphAsymmErrors* graphEtaRAASys2050           = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2050->Get("EtaRAASys");

        TH1D* histoPCMEtaToPi0RatioPbPb2050 =            (TH1D*)directoryEtaPbPbLHC11h2050->Get("EtatoPi0Ratio");
        TGraphAsymmErrors *graphPCMEtaToPi0RatioPbPb2050 = new TGraphAsymmErrors(histoPCMEtaToPi0RatioPbPb2050);
        TGraphAsymmErrors* graphPCMEtaToPi0RatioSysErrPbPb2050=    (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2050->Get("EtatoPi0RatioSys");


      TString fileNamePP2760GeVpublished      = "ExternalInputPbPb/NeutralMesonspp276GeVReference/CombinedResultsPaperPP2760GeV_2017_01_26_FrediV2Clusterizer.root";
       TFile *fileFinalResultsPP2760GeV  =        new TFile(fileNamePP2760GeVpublished.Data());
        TDirectoryFile* directoryEtaPP2760GeV         = (TDirectoryFile*)fileFinalResultsPP2760GeV->Get("Eta2.76TeV");

        TGraphAsymmErrors* graphRatioEtaToPi0PCM2760GeVStatErr = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphRatioEtaToPi0Comb2760GeVStatErr");
        TGraphAsymmErrors* graphRatioEtaToPi0PCM2760GeVSysErr = (TGraphAsymmErrors*)directoryEtaPP2760GeV->Get("graphRatioEtaToPi0Comb2760GeVStatErr");



    //Pi0 PbPb paper results (combined PCM + PHOS)
    TFile *filePublished    = new TFile(nameFilePi0PbPbPaper.Data());
    TGraphAsymmErrors* graphYieldsPublished0010     = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMStatErr_0010");
    TGraphAsymmErrors* graphYieldsSysPublished0010  = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMSysErr_0010");
    TGraphAsymmErrors* graphYieldsPublished2040     = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMStatErr_2040");
    TGraphAsymmErrors* graphYieldsSysPublished2040  = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMSysErr_2040");

    TGraphAsymmErrors* graphRAAPublished0005        = (TGraphAsymmErrors*)filePublished->Get("graphRAAStatErr_0510");
	TGraphAsymmErrors* graphRAASysPublished0005     = (TGraphAsymmErrors*)filePublished->Get("graphRAASysErr_0510");
    TGraphAsymmErrors* graphRAAPublished0510        = (TGraphAsymmErrors*)filePublished->Get("graphRAAStatErr_0510");
    TGraphAsymmErrors* graphRAASysPublished0510     = (TGraphAsymmErrors*)filePublished->Get("graphRAASysErr_0510");
    TGraphAsymmErrors* graphRAAPublished0010        = (TGraphAsymmErrors*)filePublished->Get("graphRAAStatErr_0010");
	TGraphAsymmErrors* graphRAASysPublished0010     = (TGraphAsymmErrors*)filePublished->Get("graphRAASysErr_0010");
	TGraphAsymmErrors* graphRAAPublished2040        = (TGraphAsymmErrors*)filePublished->Get("graphRAAStatErr_2040");
	TGraphAsymmErrors* graphRAASysPublished2040     = (TGraphAsymmErrors*)filePublished->Get("graphRAASysErr_2040");

	
    //Pi0 PbPb PCM only results
    TFile *filePCMPbPbLHC10h    = new TFile(nameFilePi0PbPbLHC10h.Data());
    cout << "Pi0 0-5%" << endl;
    TDirectory* directoryPi0PbPbLHC10h0005 = (TDirectory*)filePCMPbPbLHC10h->Get("Pi0_PbPb_2.76TeV_0-5%");
        TH1D* histoPCMPi0CorrectedSpecPbPbLHC10h0005                 = (TH1D*)directoryPi0PbPbLHC10h0005->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC10h0005 = (TGraphAsymmErrors*)directoryPi0PbPbLHC10h0005->Get("Pi0SystErrorA");
        TH1D* histoPi0TrueEffiPtPbPbLHC10h0005                       = (TH1D*)directoryPi0PbPbLHC10h0005->Get("Pi0_Efficiency");
        TH1D* histoPi0RawYieldPbPbLHC10h0005                         = (TH1D*)directoryPi0PbPbLHC10h0005->Get("Pi0_RawYieldPerEvent");

    cout << "Pi0 5-10%" << endl;
    TDirectory* directoryPi0PbPbLHC10h0510 = (TDirectory*)filePCMPbPbLHC10h->Get("Pi0_PbPb_2.76TeV_5-10%");
        TH1D* histoPCMPi0CorrectedSpecPbPbLHC10h0510                 = (TH1D*)directoryPi0PbPbLHC10h0510->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC10h0510 = (TGraphAsymmErrors*)directoryPi0PbPbLHC10h0510->Get("Pi0SystErrorA");
        TH1D* histoPi0TrueEffiPtPbPbLHC10h0510                       = (TH1D*)directoryPi0PbPbLHC10h0510->Get("Pi0_Efficiency");
        TH1D* histoPi0RawYieldPbPbLHC10h0510                         = (TH1D*)directoryPi0PbPbLHC10h0510->Get("Pi0_RawYieldPerEvent");

    cout << "Pi0 0-10%" << endl;
    TDirectory* directoryPi0PbPbLHC10h0010 = (TDirectory*)filePCMPbPbLHC10h->Get("Pi0_PbPb_2.76TeV_0-10%");
        TH1D* histoPCMPi0CorrectedSpecPbPbLHC10h0010                 = (TH1D*)directoryPi0PbPbLHC10h0010->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC10h0010 = (TGraphAsymmErrors*)directoryPi0PbPbLHC10h0010->Get("Pi0SystErrorA");
        TH1D* histoPi0TrueEffiPtPbPbLHC10h0010                       = (TH1D*)directoryPi0PbPbLHC10h0010->Get("Pi0_Efficiency");
        TH1D* histoPi0RawYieldPbPbLHC10h0010                         = (TH1D*)directoryPi0PbPbLHC10h0010->Get("Pi0_RawYieldPerEvent");

    cout << "20-40%" << endl;
    TDirectory* directoryPi0PbPbLHC10h2040 = (TDirectory*)filePCMPbPbLHC10h->Get("Pi0_PbPb_2.76TeV_20-40%");
        TH1D* histoPCMPi0CorrectedSpecPbPbLHC10h2040                 = (TH1D*)directoryPi0PbPbLHC10h2040->Get("CorrectedYieldPi0");   
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC10h2040 = (TGraphAsymmErrors*)directoryPi0PbPbLHC10h2040->Get("Pi0SystErrorA"); 
        TH1D* histoPi0TrueEffiPtPbPbLHC10h2040                       = (TH1D*)directoryPi0PbPbLHC10h2040->Get("Pi0_Efficiency");
        TH1D* histoPi0RawYieldPbPbLHC10h2040                         = (TH1D*)directoryPi0PbPbLHC10h2040->Get("Pi0_RawYieldPerEvent");


    //Charged pion and kaon spectra
	TFile *fileChargedKaonRAA   = new TFile(nameFileChargedPionKaonSpectra.Data());
	TH1D *histoStatChargedKaon0005                  = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Stat_0_5");
	TH1D *histoSystChargedKaon0005                  = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Syst_0_5");
	TGraphAsymmErrors* graphChargedKaonRAA0005      = new TGraphAsymmErrors(histoStatChargedKaon0005);
	TGraphAsymmErrors* graphChargedKaonRAASys0005   = new TGraphAsymmErrors(histoSystChargedKaon0005);
	
	TH1D *histoStatChargedKaon0510                  = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Stat_5_10");
	TH1D *histoSystChargedKaon0510                  = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Syst_5_10");
	TGraphAsymmErrors* graphChargedKaonRAA0510      = new TGraphAsymmErrors(histoStatChargedKaon0510);
	TGraphAsymmErrors* graphChargedKaonRAASys0510   = new TGraphAsymmErrors(histoSystChargedKaon0510);

	TH1D* histoStatChargedKaon0010                  = (TH1D*)histoStatChargedKaon0510->Clone("histoStatChargedKaon0010");
	TH1D* histoSystChargedKaon0010                  = (TH1D*)histoSystChargedKaon0510->Clone("histoSystChargedKaon0010");
	histoStatChargedKaon0010->Add(histoStatChargedKaon0005);
	histoSystChargedKaon0010->Add(histoSystChargedKaon0005);
	histoStatChargedKaon0010->Scale(0.5);
	histoSystChargedKaon0010->Scale(0.5);
	for (Int_t i = 1; i < histoSystChargedKaon0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoSystChargedKaon0005->GetBinContent(i) != 0){
			relErrLowerCent= histoSystChargedKaon0005->GetBinError(i)/histoSystChargedKaon0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoSystChargedKaon0510->GetBinContent(i) != 0){
			relErrHigherCent = histoSystChargedKaon0510->GetBinError(i)/histoSystChargedKaon0510->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoSystChargedKaon0010->SetBinError(i, histoSystChargedKaon0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoSystChargedKaon0010->SetBinError(i, histoSystChargedKaon0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   
	TGraphAsymmErrors* graphChargedKaonRAA0010      = new TGraphAsymmErrors(histoStatChargedKaon0010);
	TGraphAsymmErrors* graphChargedKaonRAASys0010   = new TGraphAsymmErrors(histoSystChargedKaon0010);

	TH1D *histoStatChargedKaon2040                  = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Stat_20_40");
	TH1D *histoSystChargedKaon2040                  = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Syst_20_40");
	TGraphAsymmErrors* graphChargedKaonRAA2040      = new TGraphAsymmErrors(histoStatChargedKaon2040);
	TGraphAsymmErrors* graphChargedKaonRAASys2040   = new TGraphAsymmErrors(histoSystChargedKaon2040);


    //Charged pion and kaon Raa
	TFile *fileChargedPionRAA                       = new TFile(nameFileChargedPionKaonRAA.Data());
	TH1D *histoStatChargedPion0005                  = (TH1D*)fileChargedPionRAA->Get("RAAPion_Stat_0_5");
	TH1D *histoSystChargedPion0005                  = (TH1D*)fileChargedPionRAA->Get("RAAPion_Syst_0_5");
	TGraphAsymmErrors* graphChargedPionRAA0005      = new TGraphAsymmErrors(histoStatChargedPion0005);
	TGraphAsymmErrors* graphChargedPionRAASys0005   = new TGraphAsymmErrors(histoSystChargedPion0005);
	
	TH1D *histoStatChargedPion0510                  = (TH1D*)fileChargedPionRAA->Get("RAAPion_Stat_5_10");
	TH1D *histoSystChargedPion0510                  = (TH1D*)fileChargedPionRAA->Get("RAAPion_Syst_5_10");
	TGraphAsymmErrors* graphChargedPionRAA0510      = new TGraphAsymmErrors(histoStatChargedPion0510);
	TGraphAsymmErrors* graphChargedPionRAASys0510   = new TGraphAsymmErrors(histoSystChargedPion0510);

	TH1D* histoStatChargedPion0010                  = (TH1D*)histoStatChargedPion0510->Clone("histoStatChargedPion0010");
	TH1D* histoSystChargedPion0010                  = (TH1D*)histoSystChargedPion0510->Clone("histoSystChargedPion0010");
	histoStatChargedPion0010->Add(histoStatChargedPion0005);
	histoSystChargedPion0010->Add(histoSystChargedPion0005);
	histoStatChargedPion0010->Scale(0.5);
	histoSystChargedPion0010->Scale(0.5);
	for (Int_t i = 1; i < histoSystChargedPion0010->GetNbinsX()+1; i++){
		Double_t relErrLowerCent = 0;
		if (histoSystChargedPion0005->GetBinContent(i) != 0){
			relErrLowerCent= histoSystChargedPion0005->GetBinError(i)/histoSystChargedPion0005->GetBinContent(i)*100 ;
		}
		Double_t relErrHigherCent = 0;
		if (histoSystChargedPion0510->GetBinContent(i) != 0){
			relErrHigherCent = histoSystChargedPion0510->GetBinError(i)/histoSystChargedPion0510->GetBinContent(i)*100 ;
		}
		
		if (relErrHigherCent > relErrLowerCent){
			histoSystChargedPion0010->SetBinError(i, histoSystChargedPion0010->GetBinContent(i)*relErrHigherCent/100);
		} else {
			histoSystChargedPion0010->SetBinError(i, histoSystChargedPion0010->GetBinContent(i)*relErrLowerCent/100);
		}         
	}   
	TGraphAsymmErrors* graphChargedPionRAA0010      = new TGraphAsymmErrors(histoStatChargedPion0010);
	TGraphAsymmErrors* graphChargedPionRAASys0010   = new TGraphAsymmErrors(histoSystChargedPion0010);
	
	TH1D *histoStatChargedPion2040                  = (TH1D*)fileChargedPionRAA->Get("RAAPion_Stat_20_40");
	TH1D *histoSystChargedPion2040                  = (TH1D*)fileChargedPionRAA->Get("RAAPion_Syst_20_40");
	TGraphAsymmErrors* graphChargedPionRAA2040      = new TGraphAsymmErrors(histoStatChargedPion2040);
	TGraphAsymmErrors* graphChargedPionRAASys2040   = new TGraphAsymmErrors(histoSystChargedPion2040);
	

	//RHIC and SPS results
	TFile *filePHENIX   = new TFile(nameFileOtherExp.Data());
	TGraphErrors* 	graphWA98_17_3GeVRAA_0013           = (TGraphErrors*)filePHENIX->Get("graphWA98RAA_0013");
	TGraphErrors* graphPHENIX200GeVRAA_0010             = (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVRAA_0010");
	TGraphErrors* graphPHENIX200GeVRAA_2040             = (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVRAA_2040");
	TGraphErrors* graphPHENIX39GeVRAA_0010              = (TGraphErrors*)filePHENIX->Get("graphPHENIX39GeVRAA_0010");
	TGraphErrors* graphPHENIX39GeVRAA_2040              = (TGraphErrors*)filePHENIX->Get("graphPHENIX39GeVRAA_2040");
	TGraphErrors* graphPHENIX62GeVRAA_0010              = (TGraphErrors*)filePHENIX->Get("graphPHENIX62GeVRAA_0010");
	TGraphErrors* graphPHENIX62GeVRAA_2040              = (TGraphErrors*)filePHENIX->Get("graphPHENIX62GeVRAA_2040");
	TGraphErrors* graphPHENIX200GeVRAAvsNPartAt7GeVc    = (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVRAAvsNPartAt7GeVc");
	TGraphErrors* graphPHENIX62GeVRAAvsNPartAt7GeVc     = (TGraphErrors*)filePHENIX->Get("graphPHENIX62GeVRAAvsNPartAt7GeVc");
	TGraphErrors* graphPHENIX39GeVRAAvsNPartAt7GeVc     = (TGraphErrors*)filePHENIX->Get("graphPHENIX39GeVRAAvsNPartAt7GeVc");

	TGraphErrors* graphPHENIX200GeVEtaRAA_0010          = (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVEtaRAA_0010");
	TGraphErrors* graphPHENIX200GeVEtaRAA_2040          = (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVEtaRAA_2040");
	TGraphErrors* graphPHENIX200GeVEtaRAA_2060          = (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVEtaRAA_2060");

	Double_t normErr0010 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0010/tAA0010),2)+pow((commonCentralityErr0010/100),2),0.5);
	Double_t normErr2040 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
	Double_t normErr2050 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
	
	
    //PP PCM only
//     TFile* fileNeutralPionDataPP = new TFile(nameFilePP.Data());
// 	TDirectory* directoryPi07TeV =      (TDirectory*)fileNeutralPionDataPP->Get("Pi07TeV");
// 	TH1D* histoAccPi07TeV =             (TH1D*)directoryPi07TeV->Get("AcceptancePi0");
// 	TH1D* histoTrueEffPtPi07TeV =       (TH1D*)directoryPi07TeV->Get("EfficiencyPi0");
// 	TH1D* histoRawYieldPi07TeV =        (TH1D*)directoryPi07TeV->Get("RAWYieldPerEventsPi0");
// 	TDirectory* directoryEta7TeV =      (TDirectory*)fileNeutralPionDataPP->Get("Eta7TeV");
// 	TH1D* histoCorrectedYieldEta7TeV = (TH1D*)directoryEta7TeV->Get("CorrectedYieldEta");
// 	TGraphAsymmErrors* graphCorrectedYieldSysEta7TeV =   (TGraphAsymmErrors*)directoryEta7TeV->Get("EtaSystError");
// 	TH1D* histoAccEta7TeV =             (TH1D*)directoryEta7TeV->Get("AcceptanceEta");
// 	TH1D* histoTrueEffPtEta7TeV =       (TH1D*)directoryEta7TeV->Get("EfficiencyEta");
// 	TH1D* histoRawYieldEta7TeV =        (TH1D*)directoryEta7TeV->Get("RAWYieldPerEventsEta");
// 	TH1D* histoEtaMassData7TeV =        (TH1D*)directoryEta7TeV->Get("MassEta");
// 	TH1D* histoEtaWidthData7TeV =       (TH1D*)directoryEta7TeV->Get("FWHMEtaMeV");
// 	TH1D* histoEtaMassMC7TeV =          (TH1D*)directoryEta7TeV->Get("TrueMassEta");
// 	TH1D* histoEtaWidthMC7TeV =         (TH1D*)directoryEta7TeV->Get("TrueFWHMEtaMeV");
// 	TH1D* histoRatioEtaPi07TeV=        (TH1D*)directoryEta7TeV->Get("EtatoPi0RatioConversion");
// 	TGraphAsymmErrors* graphRatioEtaPi0SystErr7TeV=             (TGraphAsymmErrors*)directoryEta7TeV->Get("EtatoPi0RatioConversionSys");
//
// 	histoEtaMassData7TeV->Scale(1000.);
// 	histoEtaMassMC7TeV->Scale(1000.);
// 	histoEtaMassData7TeV->SetBinContent(histoEtaMassData7TeV->GetNbinsX(),0);
// 	histoEtaMassMC7TeV->SetBinContent(histoEtaMassMC7TeV->GetNbinsX(),0);
// 	histoEtaWidthData7TeV->SetBinContent(histoEtaWidthData7TeV->GetNbinsX(),10000.);
// 	histoEtaWidthMC7TeV->SetBinContent(histoEtaWidthMC7TeV->GetNbinsX(),10000.);
//
// 	TDirectory* directoryPi0900GeV =    (TDirectory*)fileNeutralPionDataPP->Get("Pi0900GeV");
// 	TH1D* histoAccPi0900GeV =           (TH1D*)directoryPi0900GeV->Get("AcceptancePi0");
// 	TH1D* histoTrueEffPtPi0900GeV =  (TH1D*)directoryPi0900GeV->Get("EfficiencyPi0");
// 	TH1D* histoRawYieldPi0900GeV =      (TH1D*)directoryPi0900GeV->Get("RAWYieldPerEventsPi0");
// 	TDirectory* directoryEta900GeV =    (TDirectory*)fileNeutralPionDataPP->Get("Eta900GeV");
// 	TH1D* histoCorrectedYieldEta900GeV =     (TH1D*)directoryEta900GeV->Get("CorrectedYieldEta");
// 	TGraphAsymmErrors* graphCorrectedYieldSysEta900GeV =  (TGraphAsymmErrors*)directoryEta900GeV->Get("EtaSystError");
// 	TH1D* histoAccEta900GeV =           (TH1D*)directoryEta900GeV->Get("AcceptanceEta");
// 	TH1D* histoTrueEffPtEta900GeV =  (TH1D*)directoryEta900GeV->Get("EfficiencyEta");
// 	TH1D* histoRawYieldEta900GeV =      (TH1D*)directoryEta900GeV->Get("RAWYieldPerEventsEta");
// 	TH1D* histoRatioEtaPi0900GeV=      (TH1D*)directoryEta900GeV->Get("EtatoPi0RatioConversion");
// 	TGraphAsymmErrors* graphRatioEtaPi0SystErr900GeV=           (TGraphAsymmErrors*)directoryEta900GeV->Get("EtatoPi0RatioConversionSys");
	
// 	TDirectory* directoryPi02760GeV =   (TDirectory*)fileNeutralPionDataPP->Get("Pi02.76TeV");
// 	TH1D* histoAccPi02760GeV =          (TH1D*)directoryPi02760GeV->Get("AcceptancePi0");
// 	TH1D* histoTrueEffPtPi02760GeV =    (TH1D*)directoryPi02760GeV->Get("EfficiencyPi0");
// 	TH1D* histoRawYieldPi02760GeV =  (TH1D*)directoryPi02760GeV->Get("RAWYieldPerEventsPi0");
// 	TH1D* graphInvSectionPCMPi02760GeV =  (TH1D*)directoryPi02760GeV->Get("CorrectedYieldPi0");
// 	TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeV = (TGraphAsymmErrors*)directoryPi02760GeV->Get("Pi0SystError");
// 	TDirectory* directoryEta2760GeV =   (TDirectory*)fileNeutralPionDataPP->Get("Eta2.76TeV");
// 	TH1D* histoAccEta2760GeV =          (TH1D*)directoryEta2760GeV->Get("AcceptanceEta");
// 	TH1D* histoTrueEffPtEta2760GeV =    (TH1D*)directoryEta2760GeV->Get("EfficiencyEta");
// 	TH1D* graphInvSectionPCMEta2760GeV =  (TH1D*)directoryEta2760GeV->Get("CorrectedYieldEta");
// 	TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeV = (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtaSystError");
// 	TH1D* histoRawYieldEta2760GeV =  (TH1D*)directoryEta2760GeV->Get("RAWYieldPerEventsEta");
// 	TH1D* histoRatioEtaPi02760GeV=     (TH1D*)directoryEta2760GeV->Get("EtatoPi0RatioConversion");
// 	TGraphAsymmErrors* graphRatioEtaPi0SystErr2760GeV=          (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtatoPi0RatioConversionSys");
//
// 	TH1D* histoEtaPCMMassData2760GeV =     (TH1D*)directoryEta2760GeV->Get("MassEta");
// 	TH1D* histoEtaPCMWidthData2760GeV =    (TH1D*)directoryEta2760GeV->Get("FWHMEtaMeV");
// 	TH1D* histoEtaPCMMassMC2760GeV =       (TH1D*)directoryEta2760GeV->Get("TrueMassEta");
// 	TH1D* histoEtaPCMWidthMC2760GeV =      (TH1D*)directoryEta2760GeV->Get("TrueFWHMEtaMeV");
// 	histoEtaPCMMassData2760GeV->Scale(1000.);
// 	histoEtaPCMMassMC2760GeV->Scale(1000.);
// 	histoEtaPCMMassData2760GeV->SetBinContent(histoEtaPCMMassData2760GeV->GetNbinsX(),0);
// 	histoEtaPCMMassMC2760GeV->SetBinContent(histoEtaPCMMassMC2760GeV->GetNbinsX(),0);
// 	histoEtaPCMWidthData2760GeV->SetBinContent(histoEtaPCMWidthData2760GeV->GetNbinsX(),10000.);
// 	histoEtaPCMWidthMC2760GeV->SetBinContent(histoEtaPCMWidthMC2760GeV->GetNbinsX(),10000.);
//
// 	TH1D* histoPi0PCMMassData2760GeV =      (TH1D*)directoryPi02760GeV->Get("MassPi0");
// 	TH1D* histoPi0PCMWidthData2760GeV =     (TH1D*)directoryPi02760GeV->Get("FWHMPi0MeV");
// 	TH1D* histoPi0PCMMassMC2760GeV =        (TH1D*)directoryPi02760GeV->Get("TrueMassPi0");
// 	TH1D* histoPi0PCMWidthMC2760GeV =          (TH1D*)directoryPi02760GeV->Get("TrueFWHMPi0MeV");
// 	histoPi0PCMMassData2760GeV->Scale(1000.);
// 	histoPi0PCMMassMC2760GeV->Scale(1000.);
// 	histoPi0PCMMassData2760GeV->SetBinContent(histoPi0PCMMassData2760GeV->GetNbinsX(),0);
// 	histoPi0PCMMassMC2760GeV->SetBinContent(histoPi0PCMMassMC2760GeV->GetNbinsX(),0);
// 	histoPi0PCMWidthData2760GeV->SetBinContent(histoPi0PCMWidthData2760GeV->GetNbinsX(),10000.);
// 	histoPi0PCMWidthMC2760GeV->SetBinContent(histoPi0PCMWidthMC2760GeV->GetNbinsX(),10000.);

// 	TFile*   fileCocktail =                new TFile("CocktailInput/cocktail_allCentpluspp.root");
// 	TDirectory* directoryCocktailpp2760GeV =           (TDirectory*)fileCocktail->Get("cocktail_pp_2760GeV_qcd");
// 	TH1D* histoEtaFromCocktailpp2760GeV = (TH1D*)directoryCocktailpp2760GeV->Get("ptEta");
// 	TDirectory* directoryCocktail0010 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_0010_qcd");
// 	TH1D* histoEtaFromCocktail0010 = (TH1D*)directoryCocktail0010->Get("ptEta");
// 	TDirectory* directoryCocktail0510 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_0510_qcd");
// 	TH1D* histoEtaFromCocktail0510 = (TH1D*)directoryCocktail0510->Get("ptEta");
// 	TDirectory* directoryCocktail1020 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_1020_qcd");
// 	TH1D* histoEtaFromCocktail1020 = (TH1D*)directoryCocktail1020->Get("ptEta");
// 	TDirectory* directoryCocktail2040 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_2040_qcd");
// 	TH1D* histoEtaFromCocktail2040 = (TH1D*)directoryCocktail2040->Get("ptEta");
// 	TDirectory* directoryCocktail4060 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_4060_qcd");
// 	TH1D* histoEtaFromCocktail4060 = (TH1D*)directoryCocktail4060->Get("ptEta");
// 	TDirectory* directoryCocktail6080 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_6080_qcd");
// 	TH1D* histoEtaFromCocktail6080 = (TH1D*)directoryCocktail6080->Get("ptEta");
// 	cout << "here 7TeV" << endl;
	
	TFile* fileConversionsPP = new TFile(nameFilePP.Data());
	cout << "Pi0 in pp 2.76TeV" << endl;
    TDirectory *folderConversionsPi0PP = (TDirectory*)fileConversionsPP->Get("Pi02.76TeV");
    TF1* fitInvCrossSectionPi0Comb2760GeV = (TF1*)folderConversionsPi0PP->Get("TsallisFitPi0");
    fitInvCrossSectionPi0Comb2760GeV->SetParameter(0, fitInvCrossSectionPi0Comb2760GeV->GetParameter(0)*1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMPi02760GeV=       (TGraphAsymmErrors*)folderConversionsPi0PP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
    graphInvSectionPCMPi02760GeV = ScaleGraph(graphInvSectionPCMPi02760GeV,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMSysPi02760GeV=        (TGraphAsymmErrors*)folderConversionsPi0PP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
    graphInvSectionPCMSysPi02760GeV = ScaleGraph(graphInvSectionPCMSysPi02760GeV,1./(xSection2760GeVpp*recalcBarn)*factorToInel);

    TDirectory *folderConversionsPi0PPSupporting = (TDirectory*)folderConversionsPi0PP->Get("Supporting");
    TH1D* histoTrueEffPtPi02760GeV =    (TH1D*)folderConversionsPi0PPSupporting->Get("Pi0CorrectionFactorPCM");
    TH1D* histoPi0PCMMassData2760GeV =     (TH1D*)folderConversionsPi0PPSupporting->Get("Pi0MassDataPCM");
    TH1D* histoPi0PCMWidthData2760GeV =    (TH1D*)folderConversionsPi0PPSupporting->Get("Pi0WidthDataPCM");
    TH1D* histoPi0PCMMassMC2760GeV =       (TH1D*)folderConversionsPi0PPSupporting->Get("Pi0MassMCPCM");
    TH1D* histoPi0PCMWidthMC2760GeV =      (TH1D*)folderConversionsPi0PPSupporting->Get("Pi0WidthMCPCM");
    histoPi0PCMMassData2760GeV->Scale(1000.);
    histoPi0PCMMassMC2760GeV->Scale(1000.);

    cout << "Eta in pp 2.76TeV" << endl;
    TDirectory *folderConversionsEtaPP = (TDirectory*)fileConversionsPP->Get("Eta2.76TeV");
    TF1* fitInvCrossSectionEtaComb2760GeV = (TF1*)folderConversionsEtaPP->Get("TsallisFitEta");
    fitInvCrossSectionEtaComb2760GeV->SetParameter(0, fitInvCrossSectionEtaComb2760GeV->GetParameter(0)*1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMEta2760GeV=       (TGraphAsymmErrors*)folderConversionsEtaPP->Get("graphInvCrossSectionEtaPCM2760GeVStatErr");
    graphInvSectionPCMEta2760GeV = ScaleGraph(graphInvSectionPCMEta2760GeV,1./(xSection2760GeVpp*recalcBarn)*factorToInel);
    TGraphAsymmErrors* graphInvSectionPCMSysEta2760GeV=        (TGraphAsymmErrors*)folderConversionsEtaPP->Get("graphInvCrossSectionEtaPCM2760GeVSysErr");
    graphInvSectionPCMSysEta2760GeV = ScaleGraph(graphInvSectionPCMSysEta2760GeV,1./(xSection2760GeVpp*recalcBarn)*factorToInel);

    TDirectory *folderConversionsEtaPPSupporting = (TDirectory*)folderConversionsEtaPP->Get("Supporting");
    TH1D* histoTrueEffPtEta2760GeV =    (TH1D*)folderConversionsEtaPPSupporting->Get("EtaCorrectionFactorPCM");
    TH1D* histoEtaPCMMassData2760GeV =     (TH1D*)folderConversionsEtaPPSupporting->Get("EtaMassDataPCM");
    TH1D* histoEtaPCMWidthData2760GeV =    (TH1D*)folderConversionsEtaPPSupporting->Get("EtaWidthDataPCM");
    TH1D* histoEtaPCMMassMC2760GeV =       (TH1D*)folderConversionsEtaPPSupporting->Get("EtaMassMCPCM");
    TH1D* histoEtaPCMWidthMC2760GeV =      (TH1D*)folderConversionsEtaPPSupporting->Get("EtaWidthMCPCM");
    histoEtaPCMMassData2760GeV->Scale(1000.);
    histoEtaPCMMassMC2760GeV->Scale(1000.);

		
	Width_t  widthLinesBoxes            = 1.4;
	Width_t  widthCommonFit             = 2.;
	Width_t  widthStatErrBars           = 1.5;
	Width_t  widthCommonErrors          = 1.1;
	Width_t  widthCommonSpectrumBoxes         = 0.99;
	if (suffix.CompareTo("eps")==0){
		widthLinesBoxes            = 1.4;
		widthCommonFit             = 2.;
		widthStatErrBars           = 1.5;
		widthCommonErrors          = 1.1;
		widthCommonSpectrumBoxes         = 0.99;
	} else {
		widthLinesBoxes            = 2.3;
		widthCommonFit             = 2.6;
		widthStatErrBars           = 2.6;
		widthCommonErrors          = 2.;
		widthCommonSpectrumBoxes         = 2.3;
	}

	
	// *******************************************************************************************************
	// ************************** 		Acceptance x Efficiency			**************************************
	// *******************************************************************************************************   
    TH1D* histoPi0AccEffPtPbPbLHC11h0010            = (TH1D*)histoPi0AcceptPtPbPbLHC11h0010->Clone("histoPi0AccEffPtPbPbLHC11h0010");
    histoPi0AccEffPtPbPbLHC11h0010->Multiply(histoPi0TrueEffiPtPbPbLHC11h0010);
	TH1D* histoPi0AccEffPtPbPbLHC11h2050            = (TH1D*)histoPi0AcceptPtPbPbLHC11h2050->Clone("histoPi0AccEffPtPbPbLHC11h2050");
    histoPi0AccEffPtPbPbLHC11h2050->Multiply(histoPi0TrueEffiPtPbPbLHC11h2050);
	TH1D* histoEtaAccEffPtPbPbLHC11h0010            = (TH1D*)histoEtaAcceptPtPbPbLHC11h0010->Clone("histoEtaAccEffPtPbPbLHC11h0010");
    histoEtaAccEffPtPbPbLHC11h0010->Multiply(histoEtaTrueEffiPtPbPbLHC11h0010);
	TH1D* histoEtaAccEffPtPbPbLHC11h2050            = (TH1D*)histoEtaAcceptPtPbPbLHC11h2050->Clone("histoEtaAccEffPtPbPbLHC11h2050");
    histoEtaAccEffPtPbPbLHC11h2050->Multiply(histoEtaTrueEffiPtPbPbLHC11h2050);

    TCanvas* canvasAccEff = new TCanvas("canvasAccEff","",200,10,1350,1350);  // gives the page size
    DrawGammaCanvasSettings( canvasAccEff, 0.1, 0.02, 0.035, 0.09);
    canvasAccEff->SetLogy();
    
        TH2F * histo2DAccEff = new TH2F("histo2DAccEff","histo2DAccEff",1000,0.,15,2000,2e-6,1e-2 );
        SetStyleHistoTH2ForGraphs(histo2DAccEff, "#it{p}_{T} (GeV/#it{c})","Acceptance x efficiency",0.035,0.04, 0.035,0.04, .9,1.2);
        histo2DAccEff->Draw("copy");
        
        DrawGammaSetMarker(histoPi0AccEffPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);   
        histoPi0AccEffPtPbPbLHC11h0010->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoPi0AccEffPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
        histoPi0AccEffPtPbPbLHC11h2050->DrawCopy("e1,same"); 

        DrawGammaSetMarker(histoEtaAccEffPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, kPink-7, kPink-7);   
        histoEtaAccEffPtPbPbLHC11h0010->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoEtaAccEffPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, kAzure+4,  kAzure+4);   
        histoEtaAccEffPtPbPbLHC11h2050->DrawCopy("e1,same"); 

        TLegend* legendAccepEff = new TLegend(0.5,0.12,0.9,0.3); //0.25,0.13,0.93,0.43);
        legendAccepEff->SetFillColor(0);
        legendAccepEff->SetLineColor(0);
        legendAccepEff->SetTextSize(0.035);
        legendAccepEff->SetTextFont(42);
        legendAccepEff->SetMargin(0.17);
        legendAccepEff->SetNColumns(2);
        legendAccepEff->SetHeader(collisionSystemPbPb.Data());
        legendAccepEff->AddEntry(histoPi0AccEffPtPbPbLHC11h0010,"#pi^{0}, 0#font[122]{-}10%","p");
        legendAccepEff->AddEntry(histoEtaAccEffPtPbPbLHC11h0010,"#eta, 0#font[122]{-}10%","p");
        legendAccepEff->AddEntry(histoPi0AccEffPtPbPbLHC11h2050,"#pi^{0}, 20#font[122]{-}50%","p");
        legendAccepEff->AddEntry(histoEtaAccEffPtPbPbLHC11h2050,"#eta, 20#font[122]{-}50%","p");
        legendAccepEff->Draw();
        if(thesisPlotting) thesisLabelTopLeft->Draw();

	canvasAccEff->SaveAs(Form("%s/AcceEffPi0PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));

    

	// *******************************************************************************************************
	// ************************** 				Acceptance 				**************************************
	// *******************************************************************************************************   
	TCanvas* canvasAcceptance = new TCanvas("canvasAcceptance","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasAcceptance, 0.1, 0.02, 0.035, 0.09);
    
        TH2F * histo2DAcceptanceEta = new TH2F("histo2DAcceptanceEta","histo2DAcceptanceEta",1000,0.001,11,2000,0.5,1.02 );
        SetStyleHistoTH2ForGraphs(histo2DAcceptanceEta, "#it{p}_{T} (GeV/#it{c})","A_{#eta} in |y| < 0.85",0.035,0.04, 0.035,0.04, .9,1.2);
        histo2DAcceptanceEta->Draw("copy");

    // 	DrawGammaSetMarker(histoAccEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
    // 	histoAccEta2760GeV->DrawCopy("e1,same");  

        DrawGammaSetMarker(histoEtaAcceptPtPbPbLHC11h0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);   
    // 	histoEtaAcceptPtPbPbLHC11h0010->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoEtaAcceptPtPbPbLHC11h0510, markerStylePbPb0510, markerSizePbPb0510, colorCombPbPb0510, colorCombPbPb0510);   
    // 	histoEtaAcceptPtPbPbLHC11h0010->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoEtaAcceptPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);   
        histoEtaAcceptPtPbPbLHC11h0010->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoEtaAcceptPtPbPbLHC11h2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
    // 	histoEtaAcceptPtPbPbLHC11h2040->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoEtaAcceptPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
        histoEtaAcceptPtPbPbLHC11h2050->DrawCopy("e1,same"); 

        TLegend* legendAcceptance = new TLegend(0.6,0.12,0.9,0.3); //0.25,0.13,0.93,0.43);
        legendAcceptance->SetFillColor(0);
        legendAcceptance->SetLineColor(0);
        legendAcceptance->SetTextSize(0.035);
        legendAcceptance->SetTextFont(42);
        legendAcceptance->SetMargin(0.17);
        legendAcceptance->SetHeader(collisionSystemPbPb.Data());
    // 	legendAcceptance->AddEntry(histoEtaAcceptPtPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
    // 	legendAcceptance->AddEntry(histoEtaAcceptPtPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
        legendAcceptance->AddEntry(histoEtaAcceptPtPbPbLHC11h0010,"0#font[122]{-}10%","p");
        // 	legendAcceptance->AddEntry(histoEtaAcceptPtPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
    // 	legendAcceptance->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
        legendAcceptance->AddEntry(histoEtaAcceptPtPbPbLHC11h2050,"20#font[122]{-}50%","p");
    // 	legendAcceptance->AddEntry(histoAccEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
    // 	legendAcceptance->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.8","");
        legendAcceptance->Draw();
        if(thesisPlotting) thesisLabelTopLeft->Draw();
        
	canvasAcceptance->SaveAs(Form("%s/AcceptanceEtaPbPb2760GeV.%s",outputDir.Data(),suffix.Data()));
	
    canvasAcceptance->cd();

        TH2F *histo2DAcceptancePi0 = new TH2F("histo2DAcceptancePi0","histo2DAcceptancePi0",1000,0.001,15.,2000,0.7,1.05);
        SetStyleHistoTH2ForGraphs(histo2DAcceptancePi0, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}} in |y| < 0.85",0.035,0.04, 0.035,0.04, 0.9,1.2);
        histo2DAcceptancePi0->Draw("copy");

    // 	DrawGammaSetMarker(histoAccPi02760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
    // 	histoAccPi02760GeV->DrawCopy("e1,same");  

        DrawGammaSetMarker(histoPi0AcceptPtPbPbLHC11h0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);   
    // 	histoPi0AcceptPtPbPbLHC11h0005->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoPi0AcceptPtPbPbLHC11h0510, markerStylePbPb0510, markerSizePbPb0510, colorCombPbPb0510, colorCombPbPb0510);   
    // 	histoPi0AcceptPtPbPbLHC11h0510->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoPi0AcceptPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);   
        histoPi0AcceptPtPbPbLHC11h0010->DrawCopy("e1,same");
        DrawGammaSetMarker(histoPi0AcceptPtPbPbLHC11h2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
    // 	histoPi0AcceptPtPbPbLHC11h2040->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoPi0AcceptPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
        histoPi0AcceptPtPbPbLHC11h2050->DrawCopy("e1,same");

        legendAcceptance->Draw();
        if(thesisPlotting) thesisLabelTopLeft->Draw();
        
	canvasAcceptance->SaveAs(Form("%s/AcceptancePi0PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));

    canvasAcceptance->cd();

        histo2DAcceptancePi0->Draw("copy");
        
        DrawGammaSetMarker(histoPi0AcceptPtPbPbLHC11h0010, 20, 1.7, kBlack, kBlack);   
        histoPi0AcceptPtPbPbLHC11h0010->DrawCopy("e1,same");
        TLegend* legendAccSec = new TLegend(0.6, 0.12, 0.9, 0.32);
        legendAccSec->SetFillColor(0);
        legendAccSec->SetLineColor(0);
        legendAccSec->SetTextSize(0.035);
        legendAccSec->SetTextFont(42);
        legendAccSec->SetMargin(0.17);;
        legendAccSec->SetHeader(collisionSystemPbPb.Data());
        legendAccSec->AddEntry(histoPi0AcceptPtPbPbLHC11h0010,"primary #pi^{0}");
        for (Int_t j = 0; j < 4; j++){
            if (histoSecAcceptance[2][j]){
                DrawGammaSetMarker(histoSecAcceptance[2][j], markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoSecAcceptance[2][j]->Draw("e1,same");
                legendAccSec->AddEntry(histoSecAcceptance[2][j],Form("sec. from %s", nameSecMesonPlot[j].Data()));
            }    
        
        }
        legendAccSec->Draw();

    canvasAcceptance->Update();
    canvasAcceptance->SaveAs(Form("%s/Pi0WithSecAcceptance_0010.%s",outputDir.Data(),suffix.Data()));

    canvasAcceptance->cd();

        histo2DAcceptancePi0->Draw("copy");
        
        DrawGammaSetMarker(histoPi0AcceptPtPbPbLHC11h2050, 20, 1.7, kBlack, kBlack);   
        histoPi0AcceptPtPbPbLHC11h2050->DrawCopy("e1,same");
        for (Int_t j = 0; j < 4; j++){
            if (histoSecAcceptance[4][j]){
                DrawGammaSetMarker(histoSecAcceptance[4][j], markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoSecAcceptance[4][j]->Draw("e1,same");
            }    
        }
        legendAccSec->Draw();

    canvasAcceptance->Update();
    canvasAcceptance->SaveAs(Form("%s/Pi0WithSecAcceptance_2050.%s",outputDir.Data(),suffix.Data()));
	
	// *******************************************************************************************************
	// ************************** 				Efficiency 				**************************************
	// *******************************************************************************************************   
   
	TCanvas* canvasEfficiency = new TCanvas("canvasEfficiency","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasEfficiency, 0.1, 0.02, 0.035, 0.09);
	canvasEfficiency->SetLogy();
        
        TH2F * histo2DEffEta = new TH2F("histo2DEffEta","histo2DEffEta",1000,0,11,2000,2e-4,3e-3 );
        SetStyleHistoTH2ForGraphs(histo2DEffEta, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #eta}  ",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DEffEta->Draw("copy");
        
        DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005); 
    // 	histoEtaTrueEffiPtPbPbLHC11h0005->DrawCopy("e1,same");    
        DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h0010, markerStylePbPb0510, markerSizePbPb0510, colorCombPbPb0510, colorCombPbPb0510); 
    // 	histoEtaTrueEffiPtPbPbLHC11h0010->DrawCopy("e1,same");    
        DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010); 
        histoEtaTrueEffiPtPbPbLHC11h0010->DrawCopy("e1,same");     
        DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040); 
    // 	histoEtaTrueEffiPtPbPbLHC11h2040->DrawCopy("e1,same"); 
        DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040); 
        histoEtaTrueEffiPtPbPbLHC11h2050->DrawCopy("e1,same");  

        TLegend* legendEfficiency = new TLegend(0.6,0.12,0.9,0.3); //0.25,0.13,0.93,0.43);
        legendEfficiency->SetFillColor(0);
        legendEfficiency->SetLineColor(0);
        legendEfficiency->SetTextSize(0.035);
        legendEfficiency->SetTextFont(42);
        legendEfficiency->SetMargin(0.17);
        legendEfficiency->SetHeader(collisionSystemPbPb.Data());
    // 	legendEfficiency->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
    // 	legendEfficiency->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
        legendEfficiency->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0010,"0#font[122]{-}10%","p");
    // 	legendEfficiency->AddEntry(histoEtaTrueEffiPtPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
        legendEfficiency->AddEntry(histoEtaTrueEffiPtPbPbLHC11h2050,"20#font[122]{-}50%","p");
        if(runPPplotting){
        DrawGammaSetMarker(histoTrueEffPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
        histoTrueEffPtEta2760GeV->DrawCopy("pe1,same");
        legendEfficiency->AddEntry(histoTrueEffPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
        }
        legendEfficiency->Draw();

        if(thesisPlotting) thesisLabelTopLeft->Draw();
        
	canvasEfficiency->SaveAs(Form("%s/EfficiencyEtaPbPb2760GeV.%s",outputDir.Data(),suffix.Data()));
	
    canvasEfficiency->cd();
    
        TH2F * histo2DEffPi0 = new TH2F("histo2DEffPi0","histo2DEffPi0",1000,0.,15.,2000,2e-6,1.e-2   );
        SetStyleHistoTH2ForGraphs(histo2DEffPi0, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}  ",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DEffPi0->Draw("copy");

        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010); 
        histoPi0TrueEffiPtPbPbLHC11h0010->DrawCopy("e1,same");
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005); 
    // 	histoPi0TrueEffiPtPbPbLHC11h0005->DrawCopy("e1,same");    
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510); 
    // 	histoPi0TrueEffiPtPbPbLHC11h0510->DrawCopy("e1,same");    
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040); 
        histoPi0TrueEffiPtPbPbLHC11h2050->DrawCopy("e1,same");  
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040); 
    // 	histoPi0TrueEffiPtPbPbLHC11h2040->DrawCopy("e1,same");    

        if(runPPplotting){
        DrawGammaSetMarker(histoTrueEffPtPi02760GeV, markerStyleSpectrum2760GeVMC-1, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
        histoTrueEffPtPi02760GeV->DrawCopy("pe1,same");
        legendEfficiency->AddEntry(histoTrueEffPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
        }
        legendEfficiency->Draw();

        if(thesisPlotting) thesisLabelTopLeft->Draw();

	canvasEfficiency->SaveAs(Form("%s/EfficiencyPi0PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));

    
	canvasEfficiency->cd();
        histo2DEffPi0->GetYaxis()->SetTitle("#epsilon_{reco}  ");
        histo2DEffPi0->Draw("copy");
	
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010); 
        histoPi0TrueEffiPtPbPbLHC11h0010->DrawCopy("e1,same");
        DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h0010, markerStylePbPb0010+1, markerSizePbPb0010, colorCombPbPb1020, colorCombPbPb1020); 
        histoEtaTrueEffiPtPbPbLHC11h0010->DrawCopy("e1,same");  

        TLegend* legendEffi0010 = new TLegend(0.45, 0.12, 0.8, 0.25);
        legendEffi0010->SetFillColor(0);
        legendEffi0010->SetLineColor(0);
        legendEffi0010->SetTextSize(0.035);
        legendEffi0010->SetTextFont(42);
        legendEffi0010->SetMargin(0.17);
        legendEffi0010->SetNColumns(2);
        legendEffi0010->SetHeader("0#font[122]{-}10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV");
        legendEffi0010->AddEntry(histoPi0TrueEffiPtPbPbLHC11h0010,"#pi^{0}","p");
        legendEffi0010->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0010,"#eta","p");
        legendEffi0010->Draw();
	
	canvasEfficiency->SaveAs(Form("%s/MesonsEfficiency0010PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	canvasEfficiency->cd();
        histo2DEffPi0->Draw("copy");
	
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040); 
        histoPi0TrueEffiPtPbPbLHC11h2050->DrawCopy("e1,same");  
        DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h2050, markerStylePbPb2040+1, markerSizePbPb2040, colorCombPbPb4050, colorCombPbPb4050); 
        histoEtaTrueEffiPtPbPbLHC11h2050->DrawCopy("e1,same");  
	
        TLegend* legendEffi2050 = new TLegend(0.45, 0.12, 0.8, 0.25);
        legendEffi2050->SetFillColor(0);
        legendEffi2050->SetLineColor(0);
        legendEffi2050->SetTextSize(0.035);
        legendEffi2050->SetTextFont(42);
        legendEffi2050->SetMargin(0.17);
        legendEffi2050->SetNColumns(2);
        legendEffi2050->SetHeader("20#font[122]{-}50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV");
        legendEffi2050->AddEntry(histoPi0TrueEffiPtPbPbLHC11h2050,"#pi^{0}","p");
        legendEffi2050->AddEntry(histoEtaTrueEffiPtPbPbLHC11h2050,"#eta","p");
        legendEffi2050->Draw();
        
	canvasEfficiency->SaveAs(Form("%s/MesonsEfficiency2050PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));

    canvasEfficiency->cd();
        TH2F * histo2DEffPi0WithSec = new TH2F("histo2DEffPi0WithSec","histo2DEffPi0WithSec",1000,0.,15.,2000,1e-7,1.e-2);
        SetStyleHistoTH2ForGraphs(histo2DEffPi0WithSec, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco}  ",0.035,0.04, 0.035,0.04, 1.,1.);
        histo2DEffPi0WithSec->DrawCopy();
    
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h0010, 20, 1.7, kBlack, kBlack);   
        histoPi0TrueEffiPtPbPbLHC11h0010->DrawCopy("same,e1,p");
        
        TLegend* legendEffWithSec = new TLegend(0.6, 0.12, 0.9, 0.32);
        legendEffWithSec->SetFillColor(0);
        legendEffWithSec->SetLineColor(0);
        legendEffWithSec->SetTextSize(0.035);
        legendEffWithSec->SetTextFont(42);
        legendEffWithSec->SetMargin(0.17);
        legendEffWithSec->SetHeader(collisionSystemPbPb.Data());
        legendEffWithSec->AddEntry(histoPi0TrueEffiPtPbPbLHC11h0010,"val. prim","p");
        for (Int_t j = 0; j < 4; j++){
            if (histoSecTrueEffi[2][j]){
                DrawGammaSetMarker(histoSecTrueEffi[2][j],  markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoSecTrueEffi[2][j]->DrawCopy("same,e1");  
                legendEffWithSec->AddEntry(histoSecTrueEffi[2][j],Form("val. #pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
            }
        }
        legendEffWithSec->Draw();
    
    canvasEfficiency->Update();
    canvasEfficiency->SaveAs(Form("%s/Pi0WithSecTrueEff_0010.%s",outputDir.Data(),suffix.Data()));     
    
    canvasEfficiency->cd();
        histo2DEffPi0WithSec->Draw("copy");
    
        DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h2050, 20, 1.7, kBlack, kBlack);   
        histoPi0TrueEffiPtPbPbLHC11h2050->DrawCopy("same,e1,p");
        
        for (Int_t j = 0; j < 4; j++){
            if (histoSecTrueEffi[4][j]){
                DrawGammaSetMarker(histoSecTrueEffi[4][j],  markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoSecTrueEffi[4][j]->DrawCopy("same,e1");  
            }
        }
        legendEffWithSec->Draw();
    
    canvasEfficiency->Update();
    canvasEfficiency->SaveAs(Form("%s/Pi0WithSecTrueEff_2050.%s",outputDir.Data(),suffix.Data()));     

	Double_t arrayBoundariesX1_4[3];
	Double_t arrayBoundariesY1_4[2];
	Double_t relativeMarginsX[3];
	Double_t relativeMarginsY[3];
	ReturnCorrectValuesForCanvasScaling(1000,500, 2, 1,0.06, 0.005, 0.005,0.09,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);
	Double_t margin = relativeMarginsX[0]*0.8*1000;
	Double_t textsizeLabels1 = 0;
	Double_t textsizeFac1 = 0;
	Double_t textsizeLabels2 = 0;
	Double_t textsizeFac2 = 0;

	TCanvas* canvasEffRatio = new TCanvas("canvasEffRatio","",200,10,1000,500);  // gives the page size
	DrawGammaCanvasSettings( canvasEffRatio,  0.15, 0.05, 0.08, 0.09);

        TPad* pad2PartEffiRatio1 = new TPad("pad2PartEffiRatio1", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
        DrawGammaPadSettings( pad2PartEffiRatio1, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2]);
        pad2PartEffiRatio1->Draw();
        TPad* pad2PartEffiRatio3 = new TPad("pad2PartEffiRatio3", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
        DrawGammaPadSettings( pad2PartEffiRatio3, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[2]);
        pad2PartEffiRatio3->Draw();

        ReturnCorrectValuesTextSize(pad2PartEffiRatio1,textsizeLabels1, textsizeFac1, 22, margin);
        ReturnCorrectValuesTextSize(pad2PartEffiRatio3,textsizeLabels2, textsizeFac2, 22, margin);
        
        TH2F * histo2DEffRatio2 = new TH2F("histo2DEffRatio2","histo2DEffRatio2",1000,0.001,15.,1000,0.5,2.);
        SetStyleHistoTH2ForGraphs(histo2DEffRatio2, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco} #frac{20#font[122]{-}50% }{0#font[122]{-}10%}  ",0.85*textsizeLabels1, textsizeLabels1,
                                    0.85*textsizeLabels1, textsizeLabels1, 0.8,0.21/(textsizeFac1*margin));
        TH2F* histo2DEffRatio = new TH2F("histo2DEffRatio","histo2DEffRatio",1000,0.001,11,1000,0.5,2.);
        SetStyleHistoTH2ForGraphs(histo2DEffRatio, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco} #frac{20#font[122]{-}50% }{0#font[122]{-}10%}", 0.85*textsizeLabels2, textsizeLabels2,
                                    0.85*textsizeLabels2, textsizeLabels2, 0.8,0.25/(textsizeFac2*margin));
        histo2DEffRatio2->GetYaxis()->SetRangeUser(0.8,1.8);	
        histo2DEffRatio->GetYaxis()->SetRangeUser(0.8,1.8);
        
        pad2PartEffiRatio1->cd();
        histo2DEffRatio2->DrawCopy();

        TH1D* ratioTrueEffiPtPi0 = (TH1D*)histoPi0TrueEffiPtPbPbLHC11h2050->Clone("ratioTrueEffiPtPi0");
        ratioTrueEffiPtPi0->Divide(ratioTrueEffiPtPi0,histoPi0TrueEffiPtPbPbLHC11h0010,1.,1.,"B");
        DrawGammaSetMarker(ratioTrueEffiPtPi0, markerStylePbPb0010, 1.5, kBlack, kBlack); 
        ratioTrueEffiPtPi0->DrawCopy("e1,same");

        TLatex *labelSystem = new TLatex(0.6,0.9,collisionSystemPbPb.Data());
        SetStyleTLatex( labelSystem, 0.035,4);
        labelSystem->Draw();
        TLegend* legendEffiRatioPi0 = new TLegend(0.6,0.8,0.8,0.87);
        legendEffiRatioPi0->SetFillColor(0);
        legendEffiRatioPi0->SetLineColor(0);
        legendEffiRatioPi0->SetTextSize(0.035);
        legendEffiRatioPi0->AddEntry(ratioTrueEffiPtPi0,"#pi^{0}","p");
        legendEffiRatioPi0->Draw();
        DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);
	
        histo2DEffRatio2->Draw("axis,same");
        pad2PartEffiRatio1->Update();
        pad2PartEffiRatio3->cd();
        histo2DEffRatio->DrawCopy();
                
        TH1D* ratioTrueEffiPtEta = (TH1D*)histoEtaTrueEffiPtPbPbLHC11h2050->Clone("ratioTrueEffiPtEta");
        ratioTrueEffiPtEta->Divide(ratioTrueEffiPtEta,histoEtaTrueEffiPtPbPbLHC11h0010,1.,1.,"B");
        DrawGammaSetMarker(ratioTrueEffiPtEta, markerStylePbPb2040, 1.5, kBlack, kBlack); 
        ratioTrueEffiPtEta->DrawCopy("e1,same");
        labelSystem->Draw();

        TLegend* legendEffiRatioEta = new TLegend(0.6,0.8,0.8,0.87);
        legendEffiRatioEta->SetFillColor(0);
        legendEffiRatioEta->SetLineColor(0);
        legendEffiRatioEta->SetTextSize(0.035);
        legendEffiRatioEta->AddEntry(ratioTrueEffiPtEta,"#eta","p");
        legendEffiRatioEta->Draw();
        DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);	
            
        histo2DEffRatio->Draw("axis,same");
        pad2PartEffiRatio3->Update();

    canvasEffRatio->Update();	
    canvasEffRatio->SaveAs(Form("%s/EfficienciesRatio.%s",outputDir.Data(),suffix.Data()));

	
	// *******************************************************************************************************
	// ************************** 				Raw yields				**************************************
	// *******************************************************************************************************
	
    TCanvas* canvasYields = new TCanvas("canvasYields","",200,10,1350,1350*1.15);  // gives the page size
    DrawGammaCanvasSettings( canvasYields, 0.16, 0.02, 0.02, 0.09);
    canvasYields->SetLogx();
    canvasYields->SetLogy();

        TH2F *histo2DRawPi0PbPb = new TH2F("histo2DRawPi0PbPb","histo2DRawPi0PbPb",11000,0.2,20,1000,1.e-7,1);
            SetStyleHistoTH2ForGraphs(histo2DRawPi0PbPb, "#it{p}_{T} (GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N^{raw}_{#pi^{0}}}}{d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.6);
        histo2DRawPi0PbPb->Draw("copy");

        DrawGammaSetMarker(histoPi0RawYieldPbPbLHC11h0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);   
    // 	histoPi0RawYieldPbPbLHC11h0005->DrawCopy("e1,same"); 
        DrawGammaSetMarker(histoPi0RawYieldPbPbLHC11h0510, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0510, colorCombPbPb0510);   
    // 	histoPi0RawYieldPbPbLHC11h0510->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoPi0RawYieldPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);   
        histoPi0RawYieldPbPbLHC11h0010->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoPi0RawYieldPbPbLHC11h2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
    // 	histoPi0RawYieldPbPbLHC11h2040->DrawCopy("e1,same");   
        DrawGammaSetMarker(histoPi0RawYieldPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
        histoPi0RawYieldPbPbLHC11h2050->DrawCopy("e1,same");   

        if(thesisPlotting) thesisLabel->Draw();
        labelRawPi0PbPb->Draw();
        TLegend* legendRawYields = new TLegend(0.2,0.12,0.45,0.25);
        legendRawYields->SetFillColor(0);
        legendRawYields->SetLineColor(0);
        legendRawYields->SetTextSize(0.035);
        legendRawYields->SetTextFont(42);
        legendRawYields->SetHeader(collisionSystemPbPb.Data());
    // 	legendRawYields->AddEntry(histoPi0RawYieldPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
    // 	legendRawYields->AddEntry(histoPi0RawYieldPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
        legendRawYields->AddEntry(histoPi0RawYieldPbPbLHC11h0010,"0#font[122]{-}10%","p");
    // 	legendRawYields->AddEntry(histoPi0RawYieldPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
        legendRawYields->AddEntry(histoPi0RawYieldPbPbLHC11h2050,"20#font[122]{-}50%","p");
        legendRawYields->Draw();
	
	canvasYields->SaveAs(Form("%s/RawYieldsPi0PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));
	
    canvasYields->cd();

        TH2F *histo2DRawEtaPbPb = new TH2F("histo2DRawEtaPbPb","histo2DRawEtaPbPb",11000,0.7,20,1000,1.e-7,2.e-2);
            SetStyleHistoTH2ForGraphs(histo2DRawEtaPbPb, "#it{p}_{T} (GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N^{raw}_{#eta}}}{d#it{p}_{T}} ",0.035,0.04, 0.035,0.04, 1.,1.6);
        histo2DRawEtaPbPb->Draw("copy");

	DrawGammaSetMarker(histoEtaRawYieldPbPbLHC11h0510, markerStylePbPb0510, markerSizePbPb0510, colorCombPbPb0510, colorCombPbPb0510);   
// 	histoEtaRawYieldPbPbLHC11h0510->DrawCopy("e1,same");   
	DrawGammaSetMarker(histoEtaRawYieldPbPbLHC11h0005, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);   
// 	histoEtaRawYieldPbPbLHC11h0005->DrawCopy("e1,same");   
	DrawGammaSetMarker(histoEtaRawYieldPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);   
	histoEtaRawYieldPbPbLHC11h0010->DrawCopy("e1,same");   
	DrawGammaSetMarker(histoEtaRawYieldPbPbLHC11h2040, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
// 	histoEtaRawYieldPbPbLHC11h2040->DrawCopy("e1,same");
	DrawGammaSetMarker(histoEtaRawYieldPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040);   
	histoEtaRawYieldPbPbLHC11h2050->DrawCopy("e1,same");   

        if(thesisPlotting) thesisLabel->Draw();
        labelRawEtaPbPb->Draw();
        legendRawYields->Draw();
	
	canvasYields->SaveAs(Form("%s/RawYieldsEtaPbPb2760GeV.%s",outputDir.Data(),suffix.Data()));

    canvasYields->cd();
        TH2F *histo2DRawPi0WithSecPbPb = new TH2F("histo2DRawPi0WithSecPbPb","histo2DRawPi0WithSecPbPb",11000,0.2,20,1000,1.e-9,4e-1);
            SetStyleHistoTH2ForGraphs(histo2DRawPi0WithSecPbPb, "#it{p}_{T} (GeV/#it{c})","#frac{1}{#it{N}_{ev}} #frac{d#it{N^{raw}_{#pi^{0}}}}{d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.6);
        histo2DRawPi0WithSecPbPb->Draw("copy");

        TLegend* legendSecRAWYield = new TLegend(0.2,0.12,0.45,0.25);
        legendSecRAWYield->SetFillColor(0);
        legendSecRAWYield->SetLineColor(0);
        legendSecRAWYield->SetTextSize(0.035);
        legendSecRAWYield->SetTextFont(42);
        legendSecRAWYield->SetHeader(collisionSystemPbPb.Data());
        
        DrawGammaSetMarker(histoPi0RawYieldPbPbLHC11h0010, 20,1.7, kBlack, kBlack);   
        histoPi0RawYieldPbPbLHC11h0010->Draw("same,e1");
        legendSecRAWYield->AddEntry(histoPi0RawYieldPbPbLHC11h0010,"#pi^{0} raw yield","p");
        for (Int_t j = 0; j < 4; j++){
            if (histoYieldSecMeson[2][j]){
                DrawGammaSetMarker(histoYieldSecMeson[2][j],  markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoYieldSecMeson[2][j]->DrawCopy("same,e1");  
                legendSecRAWYield->AddEntry(histoYieldSecMeson[2][j],Form("#pi^{0} from %s",nameSecMesonPlot[j].Data()),"p");
            }    
        }
        legendSecRAWYield->Draw();
        
    canvasYields->Update();
    canvasYields->SaveAs(Form("%s/Pi0WithSecRAWYield_0010.%s",outputDir.Data(),suffix.Data()));
        
    canvasYields->cd();
        histo2DRawPi0WithSecPbPb->Draw("copy");

        DrawGammaSetMarker(histoPi0RawYieldPbPbLHC11h2050, 20, 1.7, kBlack, kBlack);   
        histoPi0RawYieldPbPbLHC11h2050->Draw("same,e1");
        for (Int_t j = 0; j < 4; j++){
            if (histoYieldSecMeson[4][j]){
                DrawGammaSetMarker(histoYieldSecMeson[4][j],  markerStyleSec[j] , markerSizeSec[j], colorSec[j], colorSec[j]);  
                histoYieldSecMeson[4][j]->DrawCopy("same,e1");  
            }    
        }
        legendSecRAWYield->Draw();
        
    canvasYields->Update();
    canvasYields->SaveAs(Form("%s/Pi0WithSecRAWYield_2050.%s",outputDir.Data(),suffix.Data()));
	
	// *******************************************************************************************************
	// ************************** 			Mass and Width together  	**************************************
	// *******************************************************************************************************   
	Double_t mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	Double_t mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();
	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegend[4]  = {0.,0.18,0.47,0.75};
	Double_t rowsLegend[6]     = {0.88,0.75,0.57,0.4,0.22,0.05}; //with EMCAL {0.88,0.75,0.57,0.4,0.22,0.05};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumn  = 0.13;
	Size_t textSizeTopRow      = 0.13; 
	Size_t textSizeSecondRow   = 0.11;
	//******************* Offsets ***********************
	Double_t offsetSystColumn  = 0.15;
	Double_t offsetMarkerX     = 0.1;
	Double_t offsetMarkerY     = 0.05;
	Double_t offsetBoxSizeY    = 0.05;
	Double_t offsetFit         = 0.04;
	//****************** Scale factors ******************
	Double_t scaleWidthLine       = 0.8;
		
	TCanvas * canvas6PartMassWidth = new TCanvas("canvas6PartMassWidth","",10,10,2400,1300);  // gives the page size     
	canvas6PartMassWidth->cd();
	DrawGammaCanvasSettings( canvas6PartMassWidth, 0., 0., 0., 0.);

	TPad* pad6PartMassWidth1 = new TPad("pad6PartMassWidth1", "", 0., 0.54, 0.525, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth1, 0.13, 0.0, 0.04, 0.);
	pad6PartMassWidth1->Draw();
	TPad* pad6PartMassWidth2 = new TPad("pad6PartMassWidth2", "", 0., 0., 0.525, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth2, 0.13, 0.0, 0., 0.16);
	pad6PartMassWidth2->Draw();

	TPad* pad6PartMassWidth3 = new TPad("pad6PartMassWidth3", "", 0.525, 0.54, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth3, 0.0, 0.04, 0.04, 0.);
	pad6PartMassWidth3->Draw();
	TPad* pad6PartMassWidth4 = new TPad("pad6PartMassWidth4", "", 0.525, 0., 1., 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth4, 0.0, 0.04, 0., 0.16);
	pad6PartMassWidth4->Draw();
    
	TH2D *histo2DPi0FWHM = new TH2D("histo2DPi0FWHM", "histo2DPi0FWHM", 20,0.3,20. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DPi0FWHM, "#it{p}_{T} (GeV/#it{c})","FWHM/2.36 (MeV/#it{c}^{2})", 0.035,0.05, 0.065,0.08, .95,.8, 515, 504); 
    histo2DPi0FWHM->GetYaxis()->SetRangeUser(0.2,15);
    histo2DPi0FWHM->GetYaxis()->SetLabelOffset(0.01);

	TH2D *histo2DPi0Mass = new TH2D("histo2DPi0Mass", "histo2DPi0Mass", 20,0.3,20. ,1000.,125.,150);
	SetStyleHistoTH2ForGraphs(histo2DPi0Mass, "#it{p}_{T} (GeV/#it{c})","m_{#pi^{0}} (MeV/#it{c}^{2})", 0.062,0.076, 0.062,0.078, .95,.8, 515, 510); 
	histo2DPi0Mass->GetYaxis()->SetRangeUser(130.,140.5);
	histo2DPi0Mass->GetXaxis()->SetLabelOffset(-0.02);
	histo2DPi0Mass->GetYaxis()->SetLabelOffset(0.01);
	pad6PartMassWidth1->cd();
	pad6PartMassWidth1->SetLogx();
        histo2DPi0FWHM->DrawCopy();
                
        DrawGammaSetMarker(histoPCMPi0WidthDataPbPbLHC11h0010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
        histoPCMPi0WidthDataPbPbLHC11h0010->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMPi0WidthMCPbPbLHC11h0010, markerStylePbPb0005MC, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
        histoPCMPi0WidthMCPbPbLHC11h0010->DrawCopy("same,p");

        TLatex *labelMassPbPb0010 = new TLatex(0.16,0.87,collisionSystemPbPb0010.Data());
        SetStyleTLatex( labelMassPbPb0010, 0.062,4);
        labelMassPbPb0010->Draw();

        TLatex *labelLegendTypeData = new TLatex(0.16,0.8,"Data: full points");
        SetStyleTLatex( labelLegendTypeData, 0.062,4);
        labelLegendTypeData->Draw();
        TLatex *labelLegendTypeMC = new TLatex(0.16,0.74,"MC: empty points");
        SetStyleTLatex( labelLegendTypeMC, 0.062,4);
        labelLegendTypeMC->Draw();

	pad6PartMassWidth1->Update();
	pad6PartMassWidth2->cd();
	pad6PartMassWidth2->SetLogx();
        histo2DPi0Mass->DrawCopy();

        DrawGammaSetMarker(histoPCMPi0MassDataPbPbLHC11h0010 , markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);               
        histoPCMPi0MassDataPbPbLHC11h0010->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMPi0MassMCPbPbLHC11h0010 , markerStylePbPb0005MC , markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);                
        histoPCMPi0MassMCPbPbLHC11h0010->DrawCopy("same,p");
        DrawGammaLines(0.3, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,kGray+2,3);
 
    pad6PartMassWidth2->Update();
	pad6PartMassWidth3->cd();
	pad6PartMassWidth3->SetLogx();
        histo2DPi0FWHM->DrawCopy();
 
        DrawGammaSetMarker(histoPCMPi0WidthDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
        histoPCMPi0WidthDataPbPbLHC11h2050->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMPi0WidthMCPbPbLHC11h2050, markerStylePbPb6080MC, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
        histoPCMPi0WidthMCPbPbLHC11h2050->DrawCopy("same,p");

        TLatex *labelMassPbPb2050 = new TLatex(0.05,0.87,collisionSystemPbPb2050.Data());
        SetStyleTLatex( labelMassPbPb2050, 0.062,4);
        labelMassPbPb2050->Draw();
        
    pad6PartMassWidth3->Update();
    pad6PartMassWidth4->cd();
    pad6PartMassWidth4->SetLogx();
        histo2DPi0Mass->DrawCopy();

        DrawGammaSetMarker(histoPCMPi0MassDataPbPbLHC11h2050 , markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);               
        histoPCMPi0MassDataPbPbLHC11h2050->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMPi0MassMCPbPbLHC11h2050, markerStylePbPb0005MC , markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);                
        histoPCMPi0MassMCPbPbLHC11h2050->DrawCopy("same,p");
        DrawGammaLines(0.3, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,kGray+2,3);

        TLatex *thesisLabelMassWidth = new TLatex(0.75,0.2,thisthesis.Data());
        SetStyleTLatex( thesisLabelMassWidth, 0.062,4);
        if(thesisPlotting) thesisLabelMassWidth->Draw();

    pad6PartMassWidth4->Update();
	canvas6PartMassWidth->Update();
	canvas6PartMassWidth->SaveAs(Form("%s/MassWidthPi0PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));

    canvas6PartMassWidth->cd();
    
	TH2D *histo2DEtaFWHM = new TH2D("histo2DEtaFWHM", "histo2DEtaFWHM", 20,0.5,20. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DEtaFWHM, "#it{p}_{T} (GeV/#it{c})","FWHM/2.36 (MeV/#it{c}^{2})", 0.035,0.05, 0.065,0.08, .95,.8, 515, 504); 
	histo2DEtaFWHM->GetYaxis()->SetRangeUser(0.5,15);
	histo2DEtaFWHM->GetYaxis()->SetLabelOffset(0.01);

	TH2D *histo2DEtaMass = new TH2D("histo2DEtaMass", "histo2DEtaMass", 20,0.5,20. ,1000.,520.,570);
	SetStyleHistoTH2ForGraphs(histo2DEtaMass, "#it{p}_{T} (GeV/#it{c})","m_{#eta} (MeV/#it{c}^{2})", 0.062,0.076, 0.062,0.078, .95,.8, 515, 504); 
	histo2DEtaMass->GetYaxis()->SetRangeUser(540.,562);
	histo2DEtaMass->GetXaxis()->SetLabelOffset(-0.02);
	histo2DEtaMass->GetYaxis()->SetLabelOffset(0.01);
	pad6PartMassWidth1->cd();
	pad6PartMassWidth1->SetLogx();
        histo2DEtaFWHM->DrawCopy();
                
        DrawGammaSetMarker(histoPCMEtaWidthDataPbPbLHC11h0010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
        histoPCMEtaWidthDataPbPbLHC11h0010->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMEtaWidthMCPbPbLHC11h0010, markerStylePbPb0005MC, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
        histoPCMEtaWidthMCPbPbLHC11h0010->DrawCopy("same,p");

        labelMassPbPb0010->Draw();
        labelLegendTypeData->Draw();
        labelLegendTypeMC->Draw();
    
	pad6PartMassWidth1->Update();
	pad6PartMassWidth2->cd();
	pad6PartMassWidth2->SetLogx();
        histo2DEtaMass->DrawCopy();

        DrawGammaSetMarker(histoPCMEtaMassDataPbPbLHC11h0010 , markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);               
        histoPCMEtaMassDataPbPbLHC11h0010->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMEtaMassMCPbPbLHC11h0010 , markerStylePbPb0005MC , markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);                
        histoPCMEtaMassMCPbPbLHC11h0010->DrawCopy("same,p");
        DrawGammaLines(0.5, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,kGray+2,3);
 
    pad6PartMassWidth2->Update();
	pad6PartMassWidth3->cd();
	pad6PartMassWidth3->SetLogx();
        histo2DEtaFWHM->DrawCopy();
 
        DrawGammaSetMarker(histoPCMEtaWidthDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
        histoPCMEtaWidthDataPbPbLHC11h2050->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMEtaWidthMCPbPbLHC11h2050, markerStylePbPb6080MC, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
        histoPCMEtaWidthMCPbPbLHC11h2050->DrawCopy("same,p");

        labelMassPbPb2050->Draw();
        
    pad6PartMassWidth3->Update();
    pad6PartMassWidth4->cd();
    pad6PartMassWidth4->SetLogx();
        histo2DEtaMass->DrawCopy();

        DrawGammaSetMarker(histoPCMEtaMassDataPbPbLHC11h2050 , markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);               
        histoPCMEtaMassDataPbPbLHC11h2050->DrawCopy("same,p");
        DrawGammaSetMarker(histoPCMEtaMassMCPbPbLHC11h2050, markerStylePbPb0005MC , markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);                
        histoPCMEtaMassMCPbPbLHC11h2050->DrawCopy("same,p");
        DrawGammaLines(0.5, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,kGray+2,3);

        if(thesisPlotting) thesisLabelMassWidth->Draw();

    pad6PartMassWidth4->Update();
	canvas6PartMassWidth->Update();
	canvas6PartMassWidth->SaveAs(Form("%s/MassWidthEtaPbPb2760GeV.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartMassWidth1; 
	delete pad6PartMassWidth2; 
	delete pad6PartMassWidth3; 
	delete pad6PartMassWidth4; 
	delete canvas6PartMassWidth;  
          
		
	// *******************************************************************************************************
	// ************************** 			Mass resolution syst error	**************************************
	// *******************************************************************************************************   
	TCanvas * canvas6PartMassResol = new TCanvas("canvas6PartMassResol","",10,10,2400,1300);  // gives the page size     
	canvas6PartMassResol->cd();
	DrawGammaCanvasSettings( canvas6PartMassResol, 0., 0., 0., 0.);

	TPad* pad6PartMassResol1 = new TPad("pad6PartMassResol1", "", 0., 0.54, 0.525, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassResol1, 0.13, 0.0, 0.04, 0.);
	pad6PartMassResol1->Draw();
	TPad* pad6PartMassResol2 = new TPad("pad6PartMassResol2", "", 0., 0., 0.525, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassResol2, 0.13, 0.0, 0., 0.16);
	pad6PartMassResol2->Draw();

	TPad* pad6PartMassResol3 = new TPad("pad6PartMassResol3", "", 0.525, 0.54, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassResol3, 0.0, 0.04, 0.04, 0.);
	pad6PartMassResol3->Draw();
	TPad* pad6PartMassResol4 = new TPad("pad6PartMassResol4", "", 0.525, 0., 1., 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassResol4, 0.0, 0.04, 0., 0.16);
	pad6PartMassResol4->Draw();
    
	TH2D *histo2DPi0 = new TH2D("histo2DPi0", "histo2DPi0", 20,0.3,20. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DPi0, "#it{p}_{T} (GeV/#it{c})","#pi^{0} mass resolution",0.035,0.05, 0.065,0.08, .95,.8, 515, 504); 
	histo2DPi0->GetYaxis()->SetRangeUser(-2,3);
	histo2DPi0->GetYaxis()->SetLabelOffset(0.01);

	TH2D *histo2DEta = new TH2D("histo2DEta", "histo2DEta", 20,0.3,20. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DEta, "#it{p}_{T} (GeV/#it{c})","#eta mass resolution",  0.062,0.076, 0.062,0.078, .95,.8, 515, 504); 
	histo2DEta->GetYaxis()->SetRangeUser(-1.5, 4);
	histo2DEta->GetYaxis()->SetLabelOffset(0.01);
    
    TH1D* resolutionMassPi00010 = (TH1D*)histoPCMPi0MassMCPbPbLHC11h0010->Clone();
    resolutionMassPi00010->Add(histoPCMPi0MassDataPbPbLHC11h0010,-1);
    resolutionMassPi00010->Scale(1./(mesonMassExpectPi0*1000.));
    resolutionMassPi00010->Scale(400.);
    TH1D* resolutionMassPi02050 = (TH1D*)histoPCMPi0MassMCPbPbLHC11h2050->Clone();
    resolutionMassPi02050->Add(histoPCMPi0MassDataPbPbLHC11h2050,-1);
    resolutionMassPi02050->Scale(1./(mesonMassExpectPi0*1000.));
    resolutionMassPi02050->Scale(400.);
    TH1D* resolutionMassEta0010 = (TH1D*)histoPCMEtaMassMCPbPbLHC11h0010->Clone();
    resolutionMassEta0010->Add(histoPCMEtaMassDataPbPbLHC11h0010,-1);
    resolutionMassEta0010->Scale(1./(mesonMassExpectEta*1000));
    resolutionMassEta0010->Scale(400.);
    TH1D* resolutionMassEta2050 = (TH1D*)histoPCMEtaMassMCPbPbLHC11h2050->Clone();
    resolutionMassEta2050->Add(histoPCMEtaMassDataPbPbLHC11h2050,-1);
    resolutionMassEta2050->Scale(1./(mesonMassExpectEta*1000));
    resolutionMassEta2050->Scale(400.);    
    
	pad6PartMassResol1->cd();
	pad6PartMassResol1->SetLogx();
        histo2DPi0->DrawCopy();
                
        DrawGammaLines(0.3, 20. , 1.5, 1.5, 1., kGray+2, 3);
        DrawGammaLines(0.3, 20. , 1., 1., 1., kGray+2, 1);
        DrawGammaLines(0.3, 20. , .5, .5, 1., kGray+2, 3);
        labelMassPbPb0010->Draw();

        DrawGammaSetMarker(resolutionMassPi00010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
        resolutionMassPi00010->DrawCopy("same,p");

        TF1* polPi00010 = new TF1("polPi00010","[0]",0.4,14);
        polPi00010->SetLineColor(kRed);
        polPi00010->SetLineWidth(1);    
        resolutionMassPi00010->Fit(polPi00010,"NRMEX0+");
        polPi00010->Draw("same");
        
    pad6PartMassResol1->Update();
	pad6PartMassResol2->cd();
	pad6PartMassResol2->SetLogx();
        histo2DEta->DrawCopy();

        labelMassPbPb0010->Draw();
        DrawGammaLines(0.3, 20. , 1.5, 1.5, 1., kGray+2, 3);
        DrawGammaLines(0.3, 20. , 1., 1., 1., kGray+2, 1);
        DrawGammaLines(0.3, 20. , .5, .5, 1., kGray+2, 3);
 
        DrawGammaSetMarker(resolutionMassEta0010 , markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);               
        resolutionMassEta0010->DrawCopy("same,p");

        TF1* polEta0010 = new TF1("polEta0010","[0]",1.,10);
        polEta0010->SetLineColor(kRed);
        polEta0010->SetLineWidth(1);    
        resolutionMassEta0010->Fit(polEta0010,"NRMEX0+");
        polEta0010->Draw("same");
        
    pad6PartMassResol2->Update();
	pad6PartMassResol3->cd();
	pad6PartMassResol3->SetLogx();
        histo2DPi0->DrawCopy();
 
        DrawGammaLines(0.3, 20. , 1.5, 1.5, 1., kGray+2, 3);
        DrawGammaLines(0.3, 20. , 1., 1., 1., kGray+2, 1);
        DrawGammaLines(0.3, 20. , .5, .5, 1., kGray+2, 3);
        labelMassPbPb2050->Draw();
        
        DrawGammaSetMarker(resolutionMassPi02050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
        resolutionMassPi02050->DrawCopy("same,p");
        
        TF1* polPi02050 = new TF1("polPi02050","[0]",0.4,14);
        polPi02050->SetLineColor(kRed);
        polPi02050->SetLineWidth(1);    
        resolutionMassPi02050->Fit(polPi02050,"NRMEX0+");
        polPi02050->Draw("same");
        
    pad6PartMassResol3->Update();
    pad6PartMassResol4->cd();
    pad6PartMassResol4->SetLogx();
        histo2DEta->DrawCopy();

        labelMassPbPb2050->Draw();
        DrawGammaLines(0.3, 20. , 1.5, 1.5, 1., kGray+2, 3);
        DrawGammaLines(0.3, 20. , 1., 1., 1., kGray+2, 1);
        DrawGammaLines(0.3, 20. , .5, .5, 1., kGray+2, 3);
        //         if(thesisPlotting) thesisLabelMassResol->Draw();
        
        DrawGammaSetMarker(resolutionMassEta2050 , markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);               
        resolutionMassEta2050->DrawCopy("same,p");

        TF1* polEta2050 = new TF1("polEta2050","[0]",1.,10);
        polEta2050->SetLineColor(kRed);
        polEta2050->SetLineWidth(1);    
        resolutionMassEta2050->Fit(polEta2050,"NRMEX0+");
        polEta2050->Draw("same");

    pad6PartMassResol4->Update();
	canvas6PartMassResol->Update();
	canvas6PartMassResol->SaveAs(Form("%s/MassResolutionPbPb2760GeV.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartMassResol1; 
	delete pad6PartMassResol2; 
	delete pad6PartMassResol3; 
	delete pad6PartMassResol4; 
	delete canvas6PartMassResol;  
    
    
/*	
	// *******************************************************************************************************
	// ************************** 			Mass and Width together  	**************************************
	// *******************************************************************************************************   
	
	Double_t mesonMassExpectPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
	Double_t mesonMassExpectEta = TDatabasePDG::Instance()->GetParticle(221)->Mass();

	//********************************** Defintion of the Legend **************************************************   
	Double_t columnsLegend[4]  = {0.,0.18,0.47,0.75};
	Double_t rowsLegend[6]     = {0.88,0.75,0.57,0.4,0.22,0.05}; //with EMCAL {0.88,0.75,0.57,0.4,0.22,0.05};
	//******************* Text sizes *******************
	Size_t textSizeLeftColumn  = 0.13;
	Size_t textSizeTopRow      = 0.13; 
	Size_t textSizeSecondRow   = 0.11;
	//******************* Offsets ***********************
	Double_t offsetSystColumn  = 0.15;
	Double_t offsetMarkerX     = 0.1;
	Double_t offsetMarkerY     = 0.05;
	Double_t offsetBoxSizeY    = 0.05;
	Double_t offsetFit         = 0.04;
	//****************** Scale factors ******************
	Double_t scaleWidthLine       = 0.8;
		
	TCanvas * canvas6PartMassWidth = new TCanvas("canvas6PartMassWidth","",10,10,2400,1300);  // gives the page size     
	canvas6PartMassWidth->cd();
	DrawGammaCanvasSettings( canvas6PartMassWidth, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartMassWidth1 = new TPad("pad6PartMassWidth1", "", 0., 0.54, 0.37, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth1, 0.16, 0.0, 0.02, 0.);
	pad6PartMassWidth1->Draw();
	TPad* pad6PartMassWidth2 = new TPad("pad6PartMassWidth2", "", 0., 0., 0.37, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth2, 0.16, 0.0, 0., 0.14);
	pad6PartMassWidth2->Draw();

	TPad* pad6PartMassWidth3 = new TPad("pad6PartMassWidth3", "", 0.37, 0.54, 0.685, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth3, 0.0, 0.0, 0.02, 0.);
	pad6PartMassWidth3->Draw();
	TPad* pad6PartMassWidth4 = new TPad("pad6PartMassWidth4", "", 0.37, 0., 0.685, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth4, 0.0, 0.0, 0., 0.14);
	pad6PartMassWidth4->Draw();

	TPad* pad6PartMassWidth5 = new TPad("pad6PartMassWidth5", "", 0.685, 0.54, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth5, 0.0, 0.02, 0.02, 0.);
	pad6PartMassWidth5->Draw();
	TPad* pad6PartMassWidth6 = new TPad("pad6PartMassWidth6", "", 0.685, 0., 1., 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartMassWidth6, 0.0, 0.02, 0., 0.14);
	pad6PartMassWidth6->Draw();

	TPad* padMassLegend1 = new TPad("padMassLegend1", "", 0.07, 0.09, 0.20, 0.185,-1, -1, -2);
	DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
	padMassLegend1->Draw();

	TPad* padFWHMLegend1 = new TPad("padFWHMLegend1", "", 0.07, 0.85, 0.28, 0.94,-1, -1, -2);
	DrawGammaPadSettings( padFWHMLegend1, 0., 0., 0., 0.);
	padFWHMLegend1->Draw();


	TH2D *histo2DPi0FWHM;
	histo2DPi0FWHM = new TH2D("histo2DPi0FWHM", "histo2DPi0FWHM", 20,0.35,20. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DPi0FWHM, "#it{p}_{T} (GeV/#it{c})","FWHM/2.36 (MeV/#it{c}^{2})", 0.035,0.05, 0.065,0.08, 1,1, 515, 504); 
	histo2DPi0FWHM->GetYaxis()->SetRangeUser(-0.5,10);
	histo2DPi0FWHM->GetYaxis()->SetLabelOffset(0.01);

	TH2D *histo2DPi0Mass;
	histo2DPi0Mass = new TH2D("histo2DPi0Mass", "histo2DPi0Mass", 20,0.35,20. ,1000.,125.,150);
	SetStyleHistoTH2ForGraphs(histo2DPi0Mass, "#it{p}_{T} (GeV/#it{c})","m_{#pi^{0}} (MeV/#it{c}^{2})", 0.062,0.076, 0.062,0.078, 0.8,1, 515, 510); 
	histo2DPi0Mass->GetYaxis()->SetRangeUser(130.,140.5);
	histo2DPi0Mass->GetXaxis()->SetLabelOffset(-0.02);
	histo2DPi0Mass->GetYaxis()->SetLabelOffset(0.01);
	pad6PartMassWidth1->cd();
	pad6PartMassWidth1->SetLogx();
	histo2DPi0FWHM->DrawCopy();
			
	DrawGammaSetMarker(histoPi0PCMWidthData2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
	histoPi0PCMWidthData2760GeV->DrawCopy("same,p");
	DrawGammaSetMarker(histoPi0PCMWidthMC2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoPi0PCMWidthMC2760GeV->DrawCopy("same,p");

	TLatex *labelMassPi0PP = new TLatex(0.2,0.9,collisionSystemPP2760GeV.Data());
	SetStyleTLatex( labelMassPi0PP, 0.062,4);
	labelMassPi0PP->Draw();
// 	TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
// 	SetStyleTLatex( labelLegendAMass, 0.08,4);
	labelLegendAMass->Draw();

	pad6PartMassWidth1->Update();
	pad6PartMassWidth2->cd();
	pad6PartMassWidth2->SetLogx();
	histo2DPi0Mass->DrawCopy();

	DrawGammaSetMarker(histoPi0PCMMassData2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
	histoPi0PCMMassData2760GeV->DrawCopy("same,p");
	DrawGammaSetMarker(histoPi0PCMMassMC2760GeV , markerStyleSpectrum2760GeVMC , markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoPi0PCMMassMC2760GeV->DrawCopy("same,p");

	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1.);
	TLatex *labelLegendDMass = new TLatex(0.92,0.9,"d)");
	SetStyleTLatex( labelLegendDMass, 0.075,4);
	labelLegendDMass->Draw();

	pad6PartMassWidth2->Update();

	pad6PartMassWidth5->cd();
	pad6PartMassWidth5->SetLogx();
	histo2DPi0FWHM->DrawCopy();
			
	DrawGammaSetMarker(histoPCMPi0WidthDataPbPbLHC11h0010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
	histoPCMPi0WidthDataPbPbLHC11h0010->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMPi0WidthMCPbPbLHC11h0010, markerStylePbPb0005MC, markerSizePbPb0005, colorCombMCPbPb0005, colorCombMCPbPb0005);
	histoPCMPi0WidthMCPbPbLHC11h0010->DrawCopy("same,p");

	TLatex *labelMassPi0PbPb0005 = new TLatex(0.05,0.9,collisionSystemPbPb0010.Data());
	SetStyleTLatex( labelMassPi0PbPb0005, 0.062,4);
	labelMassPi0PbPb0005->Draw();
	TLatex *labelLegendBMass = new TLatex(0.89,0.88,"c)");
	SetStyleTLatex( labelLegendBMass, 0.08,4);
	labelLegendBMass->Draw();

	TLatex *labelLegendTypeData = new TLatex(0.05,0.83,"Data: full points");
	SetStyleTLatex( labelLegendTypeData, 0.062,4);
	labelLegendTypeData->Draw();
	TLatex *labelLegendTypeMC = new TLatex(0.05,0.77,"MC: empty points");
	SetStyleTLatex( labelLegendTypeMC, 0.062,4);
	labelLegendTypeMC->Draw();

	
	pad6PartMassWidth5->Update();
	pad6PartMassWidth6->cd();
	pad6PartMassWidth6->SetLogx();
	histo2DPi0Mass->DrawCopy();

	DrawGammaSetMarker(histoPCMPi0MassDataPbPbLHC11h0010 , markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);               
	histoPCMPi0MassDataPbPbLHC11h0010->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMPi0MassMCPbPbLHC11h0010 , markerStylePbPb0005MC , markerSizePbPb0005, colorCombMCPbPb0005, colorCombMCPbPb0005);                
	histoPCMPi0MassMCPbPbLHC11h0010->DrawCopy("same,p");
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1);
	TLatex *labelLegendFMass = new TLatex(0.89,0.9,"f)");
	SetStyleTLatex( labelLegendFMass, 0.075,4);
	labelLegendFMass->Draw();

	pad6PartMassWidth6->Update();

	pad6PartMassWidth3->cd();
	pad6PartMassWidth3->SetLogx();
	histo2DPi0FWHM->DrawCopy();
		
	DrawGammaSetMarker(histoPCMPi0WidthDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
	histoPCMPi0WidthDataPbPbLHC11h2050->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMPi0WidthMCPbPbLHC11h2050, markerStylePbPb6080MC, markerSizePbPb6080, colorCombMCPbPb6080, colorCombMCPbPb6080);
	histoPCMPi0WidthMCPbPbLHC11h2050->DrawCopy("same,p");

	TLatex *labelMassPi0PbPb6080 = new TLatex(0.05,0.9,collisionSystemPbPb2050.Data());
	SetStyleTLatex( labelMassPi0PbPb6080, 0.062,4);
	labelMassPi0PbPb6080->Draw();
	TLatex *labelLegendCMass = new TLatex(0.91,0.88,"b)");
	SetStyleTLatex( labelLegendCMass, 0.08,4);
	labelLegendCMass->Draw();

	pad6PartMassWidth3->Update();
	pad6PartMassWidth4->cd();
	pad6PartMassWidth4->SetLogx();
	histo2DPi0Mass->DrawCopy();
		
	DrawGammaSetMarker(histoPCMPi0MassDataPbPbLHC11h2050 , markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);               
	histoPCMPi0MassDataPbPbLHC11h2050->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMPi0MassMCPbPbLHC11h2050, markerStylePbPb0005MC , markerSizePbPb6080, colorCombMCPbPb6080, colorCombMCPbPb6080);                
	histoPCMPi0MassMCPbPbLHC11h2050->DrawCopy("same,p");
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1);
	TLatex *labelLegendEMass = new TLatex(0.91,0.9,"e)");
	SetStyleTLatex( labelLegendEMass, 0.075,4);
	labelLegendEMass->Draw();

	pad6PartMassWidth4->Update();

    if(thesisPlotting) thesisLabel->Draw();
	canvas6PartMassWidth->Update();
	canvas6PartMassWidth->SaveAs(Form("%s/MassWidth_Pi0.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartMassWidth1; 
	delete pad6PartMassWidth2; 
	delete pad6PartMassWidth3; 
	delete pad6PartMassWidth4; 
	delete canvas6PartMassWidth;  
   

    TCanvas * canvas6PartEtaMassWidthLHC11h = new TCanvas("canvas6PartEtaMassWidthLHC11h","",10,10,2400,1300);  // gives the page size     
	canvas6PartEtaMassWidthLHC11h->cd();
	DrawGammaCanvasSettings(canvas6PartEtaMassWidthLHC11h, 0.13, 0.0, 0.02, 0.09);

	TPad* pad6PartEtaMassWidth1LHC11h = new TPad("pad6PartEtaMassWidth1LHC11h", "", 0., 0.54, 0.37, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth1LHC11h, 0.16, 0.0, 0.02, 0.);
	pad6PartEtaMassWidth1LHC11h->Draw();
	TPad* pad6PartEtaMassWidth2LHC11h = new TPad("pad6PartEtaMassWidth2LHC11h", "", 0., 0., 0.37, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth2LHC11h, 0.16, 0.0, 0., 0.14);
	pad6PartEtaMassWidth2LHC11h->Draw();

	TPad* pad6PartEtaMassWidth3LHC11h = new TPad("pad6PartEtaMassWidth3LHC11h", "", 0.37, 0.54, 0.685, 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth3LHC11h, 0.0, 0.0, 0.02, 0.);
	pad6PartEtaMassWidth3LHC11h->Draw();
	TPad* pad6PartEtaMassWidth4LHC11h = new TPad("pad6PartEtaMassWidth4LHC11h", "", 0.37, 0., 0.685, 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth4LHC11h, 0.0, 0.0, 0., 0.14);
	pad6PartEtaMassWidth4LHC11h->Draw();

	TPad* pad6PartEtaMassWidth5LHC11h = new TPad("pad6PartEtaMassWidth5LHC11h", "", 0.685, 0.54, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth5LHC11h, 0.0, 0.03, 0.02, 0.);
	pad6PartEtaMassWidth5LHC11h->Draw();
	TPad* pad6PartEtaMassWidth6LHC11h = new TPad("pad6PartEtaMassWidth6LHC11h", "", 0.685, 0., 1., 0.54,-1, -1, -2);
	DrawGammaPadSettings( pad6PartEtaMassWidth6LHC11h, 0.0, 0.03, 0., 0.14);
	pad6PartEtaMassWidth6LHC11h->Draw();

	TH2D *histo2DEtaFWHMLHC11h;
	histo2DEtaFWHMLHC11h = new TH2D("histo2DEtaFWHMLHC11h", "histo2DEtaFWHMLHC11h", 20,0.35,12. ,1000.,-30,40);
	SetStyleHistoTH2ForGraphs(histo2DEtaFWHMLHC11h, "#it{p}_{T} (GeV/#it{c})","FWHM/2.36 (MeV/#it{c}^{2})", 0.035,0.05, 0.065,0.08, 1,1, 515, 504); 
	histo2DEtaFWHMLHC11h->GetYaxis()->SetRangeUser(-0.5,20);
	histo2DEtaFWHMLHC11h->GetYaxis()->SetLabelOffset(0.01);

	TH2D *histo2DEtaMassLHC11h;
	histo2DEtaMassLHC11h = new TH2D("histo2DEtaMassLHC11h", "histo2DEtaMassLHC11h", 20,0.35,12. ,1000.,540.,559.5);
	SetStyleHistoTH2ForGraphs(histo2DEtaMassLHC11h, "#it{p}_{T} (GeV/#it{c})","m_{#eta} (MeV/#it{c}^{2})", 0.062,0.076, 0.062,0.078, 0.8,1, 515, 504); 
	// histo2DEtaMassLHC11h->GetYaxis()->SetRangeUser(130.,140.5);
	histo2DEtaMassLHC11h->GetXaxis()->SetLabelOffset(-0.02);
	histo2DEtaMassLHC11h->GetYaxis()->SetLabelOffset(0.01);
	pad6PartEtaMassWidth1LHC11h->cd();
	pad6PartEtaMassWidth1LHC11h->SetLogx();
	histo2DEtaFWHMLHC11h->DrawCopy();
			
	DrawGammaSetMarker(histoEtaPCMWidthData2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
	histoEtaPCMWidthData2760GeV->DrawCopy("same,p");
	DrawGammaSetMarker(histoEtaPCMWidthMC2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoEtaPCMWidthMC2760GeV->DrawCopy("same,p");

// 	labelMassPi0PP->Draw();
// 	TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
// 	SetStyleTLatex( labelLegendAMass, 0.08,4);

// 	labelLegendAMass->Draw();

	pad6PartEtaMassWidth1LHC11h->Update();
	pad6PartEtaMassWidth2LHC11h->cd();
	pad6PartEtaMassWidth2LHC11h->SetLogx();
	histo2DEtaMassLHC11h->DrawCopy();

	DrawGammaSetMarker(histoEtaPCMMassData2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
	histoEtaPCMMassData2760GeV->DrawCopy("same,p");
	DrawGammaSetMarker(histoEtaPCMMassMC2760GeV , markerStyleSpectrum2760GeVMC , markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoEtaPCMMassMC2760GeV->DrawCopy("same,p");

	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1.);
// 	labelLegendDMass->Draw();

	pad6PartEtaMassWidth2LHC11h->Update();

	pad6PartEtaMassWidth5LHC11h->cd();
	pad6PartEtaMassWidth5LHC11h->SetLogx();
	histo2DEtaFWHMLHC11h->DrawCopy();
			
	DrawGammaSetMarker(histoPCMEtaWidthDataPbPbLHC11h0010,markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
	histoPCMEtaWidthDataPbPbLHC11h0010->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMEtaWidthMCPbPbLHC11h0010, markerStylePbPb0005MC, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
	histoPCMEtaWidthMCPbPbLHC11h0010->DrawCopy("same,p");

	TLatex *labelMassEtaPbPbLHC11h0010 = new TLatex(0.05,0.9,collisionSystemPbPb0010.Data());
	SetStyleTLatex( labelMassEtaPbPbLHC11h0010, 0.062,4);
	labelMassEtaPbPbLHC11h0010->Draw();
// 	labelLegendBMass->Draw();
	
// 	labelLegendTypeData->Draw();
// 	labelLegendTypeMC->Draw();

	pad6PartEtaMassWidth5LHC11h->Update();
	pad6PartEtaMassWidth6LHC11h->cd();
	pad6PartEtaMassWidth6LHC11h->SetLogx();
	histo2DEtaMassLHC11h->DrawCopy();

	DrawGammaSetMarker(histoPCMEtaMassDataPbPbLHC11h0010,markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);
	histoPCMEtaMassDataPbPbLHC11h0010->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMEtaMassMCPbPbLHC11h0010, markerStylePbPb0010MC, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);
	histoPCMEtaMassMCPbPbLHC11h0010->DrawCopy("same,p");
	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1);
// 	labelLegendFMass->Draw();

	pad6PartEtaMassWidth6LHC11h->Update();

	pad6PartEtaMassWidth3LHC11h->cd();
	pad6PartEtaMassWidth3LHC11h->SetLogx();
	histo2DEtaFWHMLHC11h->DrawCopy();
	DrawGammaSetMarker(histoPCMEtaWidthDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
	histoPCMEtaWidthDataPbPbLHC11h2050->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMEtaWidthMCPbPbLHC11h2050, markerStylePbPb6080MC, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
	histoPCMEtaWidthMCPbPbLHC11h2050->DrawCopy("same,p");

	TLatex *labelMassEtaPbPb2050 = new TLatex(0.05,0.9,collisionSystemPbPb2050.Data());
	SetStyleTLatex( labelMassEtaPbPb2050, 0.062,4);
	labelMassEtaPbPb2050->Draw();
// 	labelLegendCMass->Draw();

	pad6PartEtaMassWidth3LHC11h->Update();
	pad6PartEtaMassWidth4LHC11h->cd();
	pad6PartEtaMassWidth4LHC11h->SetLogx();
	histo2DEtaMassLHC11h->DrawCopy();
	
	DrawGammaSetMarker(histoPCMEtaMassDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);              
	histoPCMEtaMassDataPbPbLHC11h2050->DrawCopy("same,p");
	DrawGammaSetMarker(histoPCMEtaMassMCPbPbLHC11h2050, markerStylePbPb6080MC , markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);            
	histoPCMEtaMassMCPbPbLHC11h2050->DrawCopy("same,p");
	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1);
// 	labelLegendEMass->Draw();

	pad6PartEtaMassWidth4LHC11h->Update();
    if(thesisPlotting) thesisLabel->Draw();

	canvas6PartEtaMassWidthLHC11h->Update();  
	canvas6PartEtaMassWidthLHC11h->SaveAs(Form("%s/MassWidth_Eta.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartEtaMassWidth1LHC11h; 
	delete pad6PartEtaMassWidth2LHC11h; 
	delete pad6PartEtaMassWidth3LHC11h; 
	delete pad6PartEtaMassWidth4LHC11h; 
	delete canvas6PartEtaMassWidthLHC11h;*/


    gStyle->SetErrorX(0.);
    gStyle->SetEndErrorSize(0); //or Z option
	
	// *******************************************************************************************************
	// ************************** 			Invariant yields			**************************************
	// *******************************************************************************************************   
	TCanvas* canvasEtaSpectra = new TCanvas("canvasEtaSpectra","",200,10,1350,1350*1.15); // gives the page size
	DrawGammaCanvasSettings( canvasEtaSpectra,  0.16, 0.02, 0.02, 0.09);
	canvasEtaSpectra->SetLogy();
	canvasEtaSpectra->SetLogx();
    TH2F *histo2DInvYieldSectionEtaLHC11h = new TH2F("histo2DInvYieldSectionEtaLHC11h","histo2DInvYieldSectionEtaLHC11h",11000,0.7,15,1000,1e-7,1e2);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionEtaLHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#eta}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",textSize,0.04, textSize,0.04, 1.,1.6);
//     histo2DInvYieldSectionEtaLHC11h->GetXaxis()->SetMoreLogLabels();
    histo2DInvYieldSectionEtaLHC11h->GetXaxis()->SetLabelOffset(-0.01);
    histo2DInvYieldSectionEtaLHC11h->Draw("copy");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
//  histoPCMEtaCorrectedSpecPbPbLHC11h0005->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005, widthLinesBoxes, kTRUE);
// 	graphPCMEtaCorrectedSpecSysPbPbLHC11h0005->Draw("E2same");
	
	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510);
//     histoPCMEtaCorrectedSpecPbPbLHC11h0510->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510, widthLinesBoxes, kTRUE);
// 	graphPCMEtaCorrectedSpecSysPbPbLHC11h0510->Draw("E2same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
	histoPCMEtaCorrectedSpecPbPbLHC11h0010->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010, widthLinesBoxes, kTRUE);
	graphPCMEtaCorrectedSpecSysPbPbLHC11h0010->Draw("E2same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
//     histoPCMEtaCorrectedSpecPbPbLHC11h2040->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
// 	graphPCMEtaCorrectedSpecSysPbPbLHC11h2040->Draw("E2same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
	histoPCMEtaCorrectedSpecPbPbLHC11h2050->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
	graphPCMEtaCorrectedSpecSysPbPbLHC11h2050->Draw("E2same");

	TLegend* legendEtaSpectra = new TLegend(0.2,0.14,0.6,0.28);// with PP points 0.2,0.12,0.8,0.3);
	legendEtaSpectra->SetFillColor(0);
	legendEtaSpectra->SetLineColor(0);
    legendEtaSpectra->SetTextFont(42);
	legendEtaSpectra->SetTextSize(textSize);
	legendEtaSpectra->SetMargin(0.17);
    legendEtaSpectra->SetHeader(collisionSystemPbPb.Data());
// 	legendEtaSpectra->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h0005,"PCM measured","pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h0010
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb0005.Data(),"");
// 	legendEtaSpectra->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h0510,"PCM measured","pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h0010
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb0510.Data(),"");
	legendEtaSpectra->AddEntry(graphPCMEtaCorrectedSpecSysPbPbLHC11h0010," 0#font[122]{-}10%","pf");//collisionSystemPbPb0010.Data(),"pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h0010
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb0010.Data(),"");
	legendEtaSpectra->AddEntry(graphPCMEtaCorrectedSpecSysPbPbLHC11h2050,"20#font[122]{-}50%","pf");//collisionSystemPbPb2050.Data(),"pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h2050
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb2050.Data(),"");
//  legendEtaSpectra->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h2040,"measured","pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h2040
//  legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb2040.Data(),"");

    if(runPPplotting){
      graphInvSectionPCMSysEta2760GeV = ScaleGraph(graphInvSectionPCMSysEta2760GeV,1e1);
      graphInvSectionPCMEta2760GeV = ScaleGraph(graphInvSectionPCMEta2760GeV,1e1);
      DrawGammaSetMarkerTGraphAsym(graphInvSectionPCMSysEta2760GeV, markerStylePP,markerSizePP, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
      graphInvSectionPCMSysEta2760GeV->Draw("E2same");
      DrawGammaSetMarkerTGraphAsym(graphInvSectionPCMEta2760GeV, markerStylePP,markerSizePP, kBlack , kBlack);
      graphInvSectionPCMEta2760GeV->Draw("p,same");

//      histoEtaFromCocktailpp2760GeV->SetLineWidth(2.);
//      histoEtaFromCocktailpp2760GeV->SetLineColor(kBlack);
//      histoEtaFromCocktailpp2760GeV->Draw("hist, l,same");

      legendEtaSpectra->AddEntry(graphInvSectionPCMSysEta2760GeV,"measured","pf");
//	legendEtaSpectra->AddEntry(histoEtaFromCocktailpp2760GeV,"m_{T} scaled","l");
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPP2760GeV.Data(),"");
//  legendEtaSpectra->AddEntry(graphInvSectionCombSysPi02760GeVPlot,collisionSystemPP.Data(),"pf");
//  legendEtaSpectra->AddEntry(fitInvCrossSectionPi0Comb2760GeV,"Tsallis Fit","l");
//  legendEtaSpectra->AddEntry(fitInvCrossSectionPi0Comb2760GeVPow,"Powerlaw Fit","l");
    }
    labelRawEtaPbPb->Draw();
    if(thesisPlotting) thesisLabel->Draw();
	legendEtaSpectra->Draw();
    histo2DInvYieldSectionEtaLHC11h->Draw("axis,same");
	canvasEtaSpectra->Update();
	canvasEtaSpectra->Print(Form("%s/SpectraEtaPbPb2760GeV.%s",outputDir.Data(),suffix.Data()));
	

	TCanvas* canvasPi0Spectra = new TCanvas("canvasPi0Spectra","",200,10,1350,1350*1.15);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0Spectra, 0.16, 0.02, 0.02, 0.09);
	canvasPi0Spectra->SetLogy();
	canvasPi0Spectra->SetLogx();
    TH2F * histo2DInvYieldSectionPi0LHC11h = new TH2F("histo2DInvYieldSectionPi0LHC11h","histo2DInvYieldSectionPi0LHC11h",11000,0.2,20,1000,1e-7,1e3);
    SetStyleHistoTH2ForGraphs(histo2DInvYieldSectionPi0LHC11h, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N_{#pi^{0}}}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}",textSize,0.04, textSize,0.04, 1.,1.5);
//     histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetMoreLogLabels();
    histo2DInvYieldSectionPi0LHC11h->GetXaxis()->SetLabelOffset(-0.01);
    histo2DInvYieldSectionPi0LHC11h->Draw("copy");

	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
//  histoPCMPi0CorrectedSpecPbPbLHC11h0005->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005, widthLinesBoxes, kTRUE);
// 	graphPCMPi0CorrectedSpecSysPbPbLHC11h0005->Draw("E2same");
	
	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510);
//     histoPCMPi0CorrectedSpecPbPbLHC11h0510->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510, widthLinesBoxes, kTRUE);
// 	graphPCMPi0CorrectedSpecSysPbPbLHC11h0510->Draw("E2same");

	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
	histoPCMPi0CorrectedSpecPbPbLHC11h0010->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010, widthLinesBoxes, kTRUE);
	graphPCMPi0CorrectedSpecSysPbPbLHC11h0010->Draw("E2same");

	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
//     histoPCMPi0CorrectedSpecPbPbLHC11h2040->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
// 	graphPCMPi0CorrectedSpecSysPbPbLHC11h2040->Draw("E2same");

	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
	histoPCMPi0CorrectedSpecPbPbLHC11h2050->Draw("p,same");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
	graphPCMPi0CorrectedSpecSysPbPbLHC11h2050->Draw("E2same");

	TLatex *labelSpectraPi0LabelPbPb = new TLatex(0.6,0.89,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraPi0LabelPbPb, 0.035  ,4);
	labelSpectraPi0LabelPbPb->Draw();

	TLegend* legendPi0Spectra = new TLegend(0.2,0.14,0.6,0.28);// with PP points 0.2,0.12,0.8,0.3);
	legendPi0Spectra->SetFillColor(0);
	legendPi0Spectra->SetLineColor(0);
    legendPi0Spectra->SetTextFont(42);
	legendPi0Spectra->SetTextSize(textSize);
	legendPi0Spectra->SetMargin(0.17);
    legendPi0Spectra->SetHeader(collisionSystemPbPb.Data());
// 	legendPi0Spectra->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h0005,"PCM measured","pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h0005
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb0005.Data(),"");
// 	legendPi0Spectra->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h0510,"PCM measured","pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h0510
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb0510.Data(),"");
	legendPi0Spectra->AddEntry(graphPCMPi0CorrectedSpecSysPbPbLHC11h0010," 0#font[122]{-}10%","pf");//collisionSystemPbPb0010.Data(),"pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h0010
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb0010.Data(),"");
// 	legendPi0Spectra->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h2040,"PCM measured","pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h2040
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb2040.Data(),"");
	legendPi0Spectra->AddEntry(graphPCMPi0CorrectedSpecSysPbPbLHC11h2050,"20#font[122]{-}50%","pf");//collisionSystemPbPb2050.Data(),"pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h2050
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb2050.Data(),"");

    if(runPPplotting){

        graphInvSectionPCMSysPi02760GeV->RemovePoint(17);
        DrawGammaSetMarkerTGraphAsym(graphInvSectionPCMSysPi02760GeV, markerStylePP,markerSizePP, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
        graphInvSectionPCMSysPi02760GeV->Draw("E2same");

        graphInvSectionPCMPi02760GeV->RemovePoint(19);
        DrawGammaSetMarkerTGraphAsym(graphInvSectionPCMPi02760GeV, markerStylePP,markerSizePP, kBlack , kBlack);
        graphInvSectionPCMPi02760GeV->Draw("p,same");

        legendPi0Spectra->AddEntry(graphInvSectionPCMSysPi02760GeV,"PCM","pf");
        legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPP2760GeV.Data(),"");
    }

    if(thesisPlotting) thesisLabel->Draw();
	legendPi0Spectra->Draw();
    histo2DInvYieldSectionPi0LHC11h->Draw("axis,same");
	canvasPi0Spectra->Update();
	canvasPi0Spectra->Print(Form("%s/SpectraPi0PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));


	// *******************************************************************************************************
	// ************************** 			Eta to Pi0 ratio			**************************************
	// *******************************************************************************************************   
    TFile *fileCocktail = new TFile(nameFileCocktailInput.Data());
    TList *folder0010 = (TList*)fileCocktail->Get("PbPb_2.76TeV_0010");
        TF1* paramPi0PCM0010 = (TF1*)folder0010->FindObject("NPionPCMStat_Fit");
        TF1* paramEtaPCM0010 = (TF1*)folder0010->FindObject("EtaPCMStat_Fit");
        TF1* paramEtaToPi0RatioPCM0010 = (TF1*)folder0010->FindObject("EtaToNPionPCMStat_Fit");
        TF1* paramCombEtaToPi0RatioPCM0010 = (TF1*)folder0010->FindObject("EtaToNPionCombStat_Fit");

        TF1 *mTScaledEtaFromPi0PbPb2760GeV_0010 = (TF1*)MtScaledParam(paramPi0PCM0010, 221, 0.476);
        TF1* etapi0RatioFromMtScalingPbPb2760GeV_0010  = DivideTF1(mTScaledEtaFromPi0PbPb2760GeV_0010, paramPi0PCM0010, "etapi0RatioFromMtScalingPbPb2760GeV_0010");
        TF1* etapi0RatioFromParamPbPb2760GeV_0010  = DivideTF1(paramEtaPCM0010, paramPi0PCM0010, "etapi0RatioFromParamPbPb2760GeV_0010");
        
        TH1D* histoetatopi0ParamToRatio_0010 = (TH1D*)paramEtaToPi0RatioPCM0010->GetHistogram();
        TH1D* histoetatopi0FromMtScaling_0010 = (TH1D*)etapi0RatioFromMtScalingPbPb2760GeV_0010->GetHistogram();
        TH1D* histoetatopi0FromRatioParams_0010 = (TH1D*)etapi0RatioFromParamPbPb2760GeV_0010->GetHistogram();
        
    TList *folder2050 = (TList*)fileCocktail->Get("PbPb_2.76TeV_2050");
        TF1* paramPi0PCM2050 = (TF1*)folder2050->FindObject("NPionPCMStat_Fit");
        TF1* paramEtaPCM2050 = (TF1*)folder2050->FindObject("EtaPCMStat_Fit");
        TF1* paramEtaToPi0RatioPCM2050 = (TF1*)folder2050->FindObject("EtaToNPionPCMStat_Fit");
        TF1* paramCombEtaToPi0RatioPCM2050 = (TF1*)folder2050->FindObject("EtaToNPionCombStat_Fit");

        TF1 *mTScaledEtaFromPi0PbPb2760GeV_2050 = (TF1*)MtScaledParam(paramPi0PCM2050, 221, 0.476);
        TF1* etapi0RatioFromMtScalingPbPb2760GeV_2050  = DivideTF1(mTScaledEtaFromPi0PbPb2760GeV_2050, paramPi0PCM2050, "etapi0RatioFromMtScalingPbPb2760GeV_2050");
        TF1* etapi0RatioFromParamPbPb2760GeV_2050  = DivideTF1(paramEtaPCM2050, paramPi0PCM2050, "etapi0RatioFromParamPbPb2760GeV_2050");
        
        TH1D* histoetatopi0ParamToRatio_2050 = (TH1D*)paramEtaToPi0RatioPCM2050->GetHistogram();
        TH1D* histoetatopi0FromMtScaling_2050 = (TH1D*)etapi0RatioFromMtScalingPbPb2760GeV_2050->GetHistogram();
        TH1D* histoetatopi0FromRatioParams_2050 = (TH1D*)etapi0RatioFromParamPbPb2760GeV_2050->GetHistogram();

    TCanvas* canvasEtatoPi0combo = new TCanvas("canvasEtatoPi0combo","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasEtatoPi0combo, 0.09, 0.03, 0.035, 0.1);
//     canvasEtatoPi0combo->SetLogx();

    TH2F * histo2DEtatoPi0combo = new TH2F("histo2DEtatoPi0combo","histo2DEtatoPi0combo",11000,0./*7*/,10.8/*15*/,1000,0.01,1.2);
    SetStyleHistoTH2ForGraphs(histo2DEtatoPi0combo, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}",0.035,0.04, 0.035,0.04, 1.2,1.);
//     histo2DEtatoPi0combo->GetXaxis()->SetMoreLogLabels();
//     histo2DEtatoPi0combo->GetXaxis()->SetLabelOffset(-0.01);
    histo2DEtatoPi0combo->GetYaxis()->SetRangeUser(0.01,1.2/*05*/);
    histo2DEtatoPi0combo->Draw("copy");

//     TCanvas* canvasRatioEtaToPi0ALICEPbPb = new TCanvas("canvasRatioEtaToPi0ALICEPbPb","",200,10,1350,900);  // gives the page size
//     DrawGammaCanvasSettings( canvasRatioEtaToPi0ALICEPbPb, 0.09, 0.01, 0.015, 0.115);
//
//     TH2D *histo2DRatioEtaToPi0ALICEPbPb = new TH2D("histo2DRatioEtaToPi0ALICEPbPb", "histo2DRatioEtaToPi0ALICEPbPb", 20,0.1,15.01,1000.,-0.4,1.3);
//     SetStyleHistoTH2ForGraphs(histo2DRatioEtaToPi0ALICEPbPb, "#it{p}_{T} (GeV/#it{c})","#eta/#pi^{0}", 0.046,0.058, 0.046,0.058, 0.85,0.74, 510, 510);
//     histo2DRatioEtaToPi0ALICEPbPb->GetYaxis()->SetRangeUser(0.,1.32);
//     histo2DRatioEtaToPi0ALICEPbPb->Draw("copy");

      if(thesisPlotting) thesisLabel3->Draw();


      DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioSysErrPbPb0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010, 2, kTRUE);
      DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioPbPb0010,markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
      ProduceGraphAsymmWithoutXErrors(graphPCMEtaToPi0RatioPbPb0010);
      DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioSysErrPbPb2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, 2, kTRUE);
      DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioPbPb2050,  markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
      ProduceGraphAsymmWithoutXErrors(graphPCMEtaToPi0RatioPbPb2050);
      DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioSysErrPbPb2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, 2, kTRUE);
      DrawGammaSetMarkerTGraphAsym(graphPCMEtaToPi0RatioPbPb2040,  markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
      ProduceGraphAsymmWithoutXErrors(graphPCMEtaToPi0RatioPbPb2040);

      graphPCMEtaToPi0RatioSysErrPbPb0010->Draw("E2same");
      graphPCMEtaToPi0RatioPbPb0010->Draw("p,same");
      graphPCMEtaToPi0RatioSysErrPbPb2050->Draw("E2same");
      graphPCMEtaToPi0RatioPbPb2050->Draw("p,same");

      TLegend* legendEtatoPi0combo_onlyPbPb = new TLegend(0.12,0.73,0.53,0.89);
      legendEtatoPi0combo_onlyPbPb->SetFillColor(0);
      legendEtatoPi0combo_onlyPbPb->SetLineColor(0);
      legendEtatoPi0combo_onlyPbPb->SetTextFont(42);
      legendEtatoPi0combo_onlyPbPb->SetTextSize(textSize);
      legendEtatoPi0combo_onlyPbPb->SetMargin(0.17);
      legendEtatoPi0combo_onlyPbPb->SetHeader(collisionSystemPbPb.Data());
  //     legendRatioALICEdata2->AddEntry(graphCombEtaToPi0RatioSysErrpp7TeV,Form("ALICE, %s",collisionSystemPP7TeV.Data()),"pf");
  //     legendRatioALICEdata2->AddEntry((TObject*)0,"Phys.Lett. B717 (2012) 162-172","");
      legendEtatoPi0combo_onlyPbPb->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb0010," 0#font[122]{-}10%","pf");//collisionSystemPbPb0010.Data(),"pf");
      legendEtatoPi0combo_onlyPbPb->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb2050,"20#font[122]{-}50%","pf");//collisionSystemPbPb2050.Data(),"pf");
      legendEtatoPi0combo_onlyPbPb->Draw();

    canvasEtatoPi0combo->Update();  //qui
    canvasEtatoPi0combo->SaveAs(Form("%s/EtaToPi0RatioPbPb2760GeV.%s",outputDir.Data(), suffix.Data()));


    canvasEtatoPi0combo->cd();
        histo2DEtatoPi0combo->Draw("copy");

        graphPCMEtaToPi0RatioSysErrPbPb0010->Draw("E2same");
        graphPCMEtaToPi0RatioPbPb0010->Draw("p,same");
        graphPCMEtaToPi0RatioSysErrPbPb2050->Draw("E2same");
        graphPCMEtaToPi0RatioPbPb2050->Draw("p,same");

        DrawGammaSetMarkerTGraphAsym(graphRatioEtaToPi0PCM2760GeVSysErr, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);
        DrawGammaSetMarkerTGraphAsym(graphRatioEtaToPi0PCM2760GeVStatErr, markerStylepp, markerSizepp, kBlack, kBlack, 1, kTRUE);
        graphRatioEtaToPi0PCM2760GeVStatErr->Draw("p,same");
        graphRatioEtaToPi0PCM2760GeVSysErr->Draw("E2same");

        legendEtatoPi0combo_onlyPbPb->Draw();
        if(thesisPlotting) thesisLabel3->Draw();

        TLegend* legendEtatoPi0combo_withPP276GeV = new TLegend(0.55,0.15,0.95,0.26);
        legendEtatoPi0combo_withPP276GeV->SetFillColor(0);
        legendEtatoPi0combo_withPP276GeV->SetLineColor(0);
        legendEtatoPi0combo_withPP276GeV->SetTextFont(42);
        legendEtatoPi0combo_withPP276GeV->SetTextSize(textSize);
        legendEtatoPi0combo_withPP276GeV->SetMargin(0.17);
        legendEtatoPi0combo_withPP276GeV->SetHeader(collisionSystemPP2760GeV.Data());
        legendEtatoPi0combo_withPP276GeV->AddEntry(graphRatioEtaToPi0PCM2760GeVSysErr,"arXiv: 1702.00917","fp");
        legendEtatoPi0combo_withPP276GeV->Draw();

    canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioPbPb2760GeV_WithPP.%s",outputDir.Data(),suffix.Data()));

    canvasEtatoPi0combo->cd();
        histo2DEtatoPi0combo->Draw("copy");

        graphPCMEtaToPi0RatioSysErrPbPb0010->Draw("E2same");
        graphPCMEtaToPi0RatioPbPb0010->Draw("p,same");

        paramEtaToPi0RatioPCM0010->SetLineColor(kRed+2);
        paramEtaToPi0RatioPCM0010->SetLineStyle(2);
        paramEtaToPi0RatioPCM0010->SetLineWidth(2);
        paramEtaToPi0RatioPCM0010->Draw("c,histo,same");

        paramCombEtaToPi0RatioPCM0010->SetLineColor(kMagenta+2);
        paramCombEtaToPi0RatioPCM0010->SetLineStyle(4);
        paramCombEtaToPi0RatioPCM0010->SetLineWidth(2);
        paramCombEtaToPi0RatioPCM0010->Draw("c,histo,same");
        
        etapi0RatioFromMtScalingPbPb2760GeV_0010->SetLineColor(kBlue+1);
        etapi0RatioFromMtScalingPbPb2760GeV_0010->SetLineStyle(3);
        etapi0RatioFromMtScalingPbPb2760GeV_0010->SetLineWidth(2);
        //         etapi0RatioFromMtScalingPbPb2760GeV_0010->GetXaxis()->SetRangeUser(1,11);
        etapi0RatioFromMtScalingPbPb2760GeV_0010->Draw("c,histo,same");

        etapi0RatioFromParamPbPb2760GeV_0010->SetLineColor(kGreen+2);
        etapi0RatioFromParamPbPb2760GeV_0010->SetLineStyle(5);
        etapi0RatioFromParamPbPb2760GeV_0010->SetLineWidth(2);
        //         etapi0RatioFromParamPbPb2760GeV_0010->GetXaxis()->SetRangeUser(0.8,11);
        etapi0RatioFromParamPbPb2760GeV_0010->Draw("c,histo,same");

        if(thesisPlotting) thesisLabel3->Draw();

        TLegend* legendEtatoPi0_withCurves = new TLegend(0.12,0.67,0.53,0.89);
        legendEtatoPi0_withCurves->SetFillColor(0);
        legendEtatoPi0_withCurves->SetLineColor(0);
        legendEtatoPi0_withCurves->SetTextFont(42);
        legendEtatoPi0_withCurves->SetTextSize(textSize);
        legendEtatoPi0_withCurves->SetMargin(0.17);
        legendEtatoPi0_withCurves->SetHeader(collisionSystemPbPb.Data());
        legendEtatoPi0_withCurves->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb0010," 0#font[122]{-}10%","pf");//collisionSystemPbPb0010.Data(),"pf");
        legendEtatoPi0_withCurves->AddEntry(paramEtaToPi0RatioPCM0010,"param to PCM ratio","l");
        legendEtatoPi0_withCurves->AddEntry(paramCombEtaToPi0RatioPCM0010,"param to combined PCM-EMCal ratio","l");
        legendEtatoPi0_withCurves->AddEntry(etapi0RatioFromParamPbPb2760GeV_0010,"PCM spectra param.","l");
        legendEtatoPi0_withCurves->AddEntry(etapi0RatioFromMtScalingPbPb2760GeV_0010,"m_{T}-scaled #eta from PCM #pi^{0} param.","l");
        legendEtatoPi0_withCurves->Draw();

        TLegend* legendEtatoPi0combo_withCurves = new TLegend(0.55,0.15,0.95,0.26);
        legendEtatoPi0combo_withCurves->SetFillColor(0);
        legendEtatoPi0combo_withCurves->SetLineColor(0);
        legendEtatoPi0combo_withCurves->SetTextFont(42);
        legendEtatoPi0combo_withCurves->SetTextSize(textSize);
        legendEtatoPi0combo_withCurves->SetMargin(0.17);
        legendEtatoPi0combo_withCurves->SetHeader(collisionSystemPP2760GeV.Data());
        legendEtatoPi0combo_withCurves->AddEntry(graphRatioEtaToPi0PCM2760GeVSysErr,"arXiv: 1702.00917","fp");
        //         legendEtatoPi0combo_withCurves->Draw();

    canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioPbPb2760GeV_WithCurves_0010.%s",outputDir.Data(),suffix.Data()));

   canvasEtatoPi0combo->cd();
        histo2DEtatoPi0combo->Draw("copy");

        graphPCMEtaToPi0RatioSysErrPbPb2050->Draw("E2same");
        graphPCMEtaToPi0RatioPbPb2050->Draw("p,same");

        paramEtaToPi0RatioPCM2050->SetLineColor(kRed+2);
        paramEtaToPi0RatioPCM2050->SetLineStyle(2);
        paramEtaToPi0RatioPCM2050->SetLineWidth(2);
        paramEtaToPi0RatioPCM2050->Draw("c,histo,same");

        paramCombEtaToPi0RatioPCM2050->SetLineColor(kMagenta+2);
        paramCombEtaToPi0RatioPCM2050->SetLineStyle(4);
        paramCombEtaToPi0RatioPCM2050->SetLineWidth(2);
        paramCombEtaToPi0RatioPCM2050->Draw("c,histo,same");
        
        etapi0RatioFromMtScalingPbPb2760GeV_2050->SetLineColor(kBlue+1);
        etapi0RatioFromMtScalingPbPb2760GeV_2050->SetLineStyle(3);
        etapi0RatioFromMtScalingPbPb2760GeV_2050->SetLineWidth(2);
        //         etapi0RatioFromMtScalingPbPb2760GeV_2050->GetXaxis()->SetRangeUser(1,11);
        etapi0RatioFromMtScalingPbPb2760GeV_2050->Draw("c,histo,same");

        etapi0RatioFromParamPbPb2760GeV_2050->SetLineColor(kGreen+2);
        etapi0RatioFromParamPbPb2760GeV_2050->SetLineStyle(5);
        etapi0RatioFromParamPbPb2760GeV_2050->SetLineWidth(2);
        //         etapi0RatioFromParamPbPb2760GeV_2050->GetXaxis()->SetRangeUser(0.8,11);
        etapi0RatioFromParamPbPb2760GeV_2050->Draw("c,histo,same");

        if(thesisPlotting) thesisLabel3->Draw();

        TLegend* legendEtatoPi0_withCurves2 = new TLegend(0.12,0.67,0.53,0.89);
        legendEtatoPi0_withCurves2->SetFillColor(0);
        legendEtatoPi0_withCurves2->SetLineColor(0);
        legendEtatoPi0_withCurves2->SetTextFont(42);
        legendEtatoPi0_withCurves2->SetTextSize(textSize);
        legendEtatoPi0_withCurves2->SetMargin(0.17);
        legendEtatoPi0_withCurves2->SetHeader(collisionSystemPbPb.Data());
        legendEtatoPi0_withCurves2->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb2050,"20#font[122]{-}50%","pf");//collisionSystemPbPb2050.Data(),"pf");
        legendEtatoPi0_withCurves2->AddEntry(paramEtaToPi0RatioPCM2050,"param to PCM ratio","l");
        legendEtatoPi0_withCurves2->AddEntry(paramCombEtaToPi0RatioPCM2050,"param to combined PCM-EMCal ratio","l");
        legendEtatoPi0_withCurves2->AddEntry(etapi0RatioFromParamPbPb2760GeV_2050,"PCM spectra param.","l");
        legendEtatoPi0_withCurves2->AddEntry(etapi0RatioFromMtScalingPbPb2760GeV_2050,"m_{T}-scaled #eta from PCM #pi^{0} param.","l");
        legendEtatoPi0_withCurves2->Draw();

//         TLegend* legendEtatoPi0combo_withCurves = new TLegend(0.55,0.15,0.95,0.26);
//         legendEtatoPi0combo_withCurves->SetFillColor(0);
//         legendEtatoPi0combo_withCurves->SetLineColor(0);
//         legendEtatoPi0combo_withCurves->SetTextFont(42);
//         legendEtatoPi0combo_withCurves->SetTextSize(textSize);
//         legendEtatoPi0combo_withCurves->SetMargin(0.17);
//         legendEtatoPi0combo_withCurves->SetHeader(collisionSystemPP2760GeV.Data());
//         legendEtatoPi0combo_withCurves->AddEntry(graphRatioEtaToPi0PCM2760GeVSysErr,"arXiv: 1702.00917","fp");
        //         legendEtatoPi0combo_withCurves->Draw();

    canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioPbPb2760GeV_WithCurves_2050.%s",outputDir.Data(),suffix.Data()));

    
    TFile *filePhotons0010 = new TFile("/home/admin1/leardini/Results/DirectPhotonsAnalysisApril2017/50100013_00200009847005008750404000_0152501500000000/PbPb_2.76TeV/Gamma_Pi0_MC_GammaConvV1Correction_50100013_00200009847005008750404000_0152501500000000.root");
    TFile *filePhotons2040 = new TFile("/home/admin1/leardini/Results/DirectPhotonsAnalysisApril2017/52400013_00200009847005008750404000_0152501500000000/PbPb_2.76TeV/Gamma_Pi0_MC_GammaConvV1Correction_52400013_00200009847005008750404000_0152501500000000.root");
    
    TH1D* purity_0010 = (TH1D*)filePhotons0010->Get("GammaPurityWOSec_Pt");
    TH1D* purity_2040 = (TH1D*)filePhotons2040->Get("GammaPurityWOSec_Pt");
    
    TCanvas *canvasPurity                       = GetAndSetCanvas("canvasPurity");
    TH2F * histo2DPurity = new TH2F("histo2DPurity","histo2DPurity",11000,0.,16,1000,0.,1.2);
    SetStyleHistoTH2ForGraphs(histo2DPurity, "#it{p}_{T} (GeV/#it{c})","#gamma purity",0.035,0.04, 0.035,0.04, 1.1,1.);
    histo2DPurity->GetYaxis()->SetRangeUser(0.805,1.05);
    histo2DPurity->GetXaxis()->SetRangeUser(0.5,13.5);
    histo2DPurity->Draw("copy");

        DrawGammaSetMarker(purity_0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
        purity_0010->Draw("same");
        DrawGammaSetMarker(purity_2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
        purity_2040->Draw("same");
        DrawGammaLines(0.5, 13.5 , 1, 1 ,1,kGray,2);

        cout << "fitting purity 0010" << endl;
        TF1* polPurity0010 = new TF1("polPurity0010","[0]",2.,14);
        polPurity0010->SetLineColor(kRed+2);
        polPurity0010->SetLineWidth(1);    
        purity_0010->Fit(polPurity0010,"NRMEX0+");
//         polPurity0010->Draw("same");

        cout << "fitting purity 2040" << endl;
        TF1* polPurity2040 = new TF1("polPurity2040","[0]",2.,14);
        polPurity2040->SetLineColor(kBlue+1);
        polPurity2040->SetLineWidth(1);    
        purity_2040->Fit(polPurity2040,"NRMEX0+");
//         polPurity2040->Draw("same");
        
        TLegend* legendPurity = GetAndSetLegend(0.15,0.15,3);
        legendPurity->SetHeader(collisionSystemPbPb.Data());
        legendPurity->AddEntry(purity_0010,"0#font[122]{-}10%","p");
        legendPurity->AddEntry(purity_2040,"20#font[122]{-}40%","p");
        legendPurity->Draw();
        TLatex *thesisLabelPurity = new TLatex(0.16,0.3,thisthesis.Data());
        SetStyleTLatex( thesisLabelPurity, textSize,4);
        if(thesisPlotting) thesisLabelPurity->Draw();

    canvasPurity->SaveAs(Form("%s/PurityGamma.%s",outputDir.Data(),suffix.Data()));
    delete canvasPurity;
    
    
//       TFile *fileChargedRatios = new TFile("pdf/2017_01_29/CombineMesonMeasurementsPbPb2760GeVX/InputALICEResultsPbPb2760GeV_2017_01_29.root");
//       if(fileChargedRatios){
//         TDirectoryFile* directoryChargedPbPb = (TDirectoryFile*)fileChargedRatios->Get("ChargedParticles_PbPb_2.76TeV");
//         TGraphAsymmErrors* graphChargedRatioKaonToPion0010      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion0010");
//         TGraphAsymmErrors* graphChargedRatioKaonToPionSys0010   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys0010");
//         TGraphAsymmErrors* graphChargedRatioKaonToPion2040      = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPion2040");
//         TGraphAsymmErrors* graphChargedRatioKaonToPionSys2040   = (TGraphAsymmErrors*)directoryChargedPbPb->Get("graphChargedRatioKaonToPionSys2040");
//
//
//         canvasEtatoPi0combo->cd();
//             histo2DEtatoPi0combo->Draw("copy");
//             histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
//             histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");
//
//             DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPion0010, 24,markerSizeComb, colorCharged, colorCharged);
//             DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPionSys0010,24,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
//             graphChargedRatioKaonToPionSys0010->Draw("2same");
//             graphChargedRatioKaonToPion0010->Draw("p,same");
//
//             TLegend* legendChargedRatio0010 = new TLegend(0.12,0.7,0.5,0.88);
//             legendChargedRatio0010->SetFillColor(0);
//             legendChargedRatio0010->SetLineColor(0);
//             legendChargedRatio0010->SetTextFont(42);
//             legendChargedRatio0010->SetTextSize(textSize);
//             legendChargedRatio0010->SetMargin(0.17);
//             legendChargedRatio0010->SetHeader(collisionSystemPbPb.Data());
//   //           legendChargedRatio0010->SetNColumns(2);
//             legendChargedRatio0010->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb0010,Form("#eta/#pi^{0},   %s",cent0010.Data()),"fp");
//             legendChargedRatio0010->AddEntry(graphChargedRatioKaonToPionSys0010,Form("K^{#pm}/#pi^{#pm}, %s",cent0010.Data()),"fp");
//             legendChargedRatio0010->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
//             legendChargedRatio0010->Draw();
//
//             if(thesisPlotting) thesisLabel3->Draw();
//             graphPCMEtaToPi0RatioSysErrPbPb0010->Draw("E2same");
//             graphPCMEtaToPi0RatioPbPb0010->Draw("p,same");
//
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions_0010.%s",outputDir.Data(),suffix.Data()));
//
//         canvasEtatoPi0combo->cd();
//             histo2DEtatoPi0combo->Draw("copy");
//             histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
//             histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");
//
//             DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPion2040, 25,markerSizeComb, colorCharged, colorCharged);
//             DrawGammaSetMarkerTGraphAsym(graphChargedRatioKaonToPionSys2040,25,markerSizeComb, colorCharged, colorCharged, 1, kTRUE);
//             graphChargedRatioKaonToPionSys2040->Draw("2same");
//             graphChargedRatioKaonToPion2040->Draw("p,same");
//
//             TLegend* legendChargedRatio2040 = new TLegend(0.12,0.7,0.5,0.88);
//             legendChargedRatio2040->SetFillColor(0);
//             legendChargedRatio2040->SetLineColor(0);
//             legendChargedRatio2040->SetTextFont(42);
//             legendChargedRatio2040->SetTextSize(textSize);
//             legendChargedRatio2040->SetMargin(0.17);
//             legendChargedRatio2040->SetHeader(collisionSystemPbPb.Data());
//   //           legendChargedRatio2040->SetNColumns(2);
// //             legendChargedRatio2040->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb2040,Form("#eta/#pi^{0},   %s",cent2040.Data()),"fp");
//             legendChargedRatio2040->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb2040,Form("#eta/#pi^{0},   %s",cent2040.Data()),"fp");
//             legendChargedRatio2040->AddEntry(graphChargedRatioKaonToPionSys2040,Form("K^{#pm}/#pi^{#pm}, %s",cent2040.Data()),"fp");
//             legendChargedRatio2040->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
//             legendChargedRatio2040->Draw();
//
//             if(thesisPlotting) thesisLabel3->Draw();
//             graphPCMEtaToPi0RatioSysErrPbPb2040->Draw("E2same");
//             graphPCMEtaToPi0RatioPbPb2040->Draw("p,same");
//
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions_2040.%s",outputDir.Data(),suffix.Data()));
//
//         canvasEtatoPi0combo->cd();
//             histo2DEtatoPi0combo->Draw("copy");
//             histo2DEtatoPi0combo->GetYaxis()->SetTitleOffset(1.1);
//             histo2DEtatoPi0combo->GetYaxis()->SetTitle("Particle ratio");
//
//             graphChargedRatioKaonToPionSys2040->Draw("2same");
//             graphChargedRatioKaonToPion2040->Draw("p,same");
//
//             TLegend* legendChargedRatio2050 = new TLegend(0.12,0.7,0.5,0.88);
//             legendChargedRatio2050->SetFillColor(0);
//             legendChargedRatio2050->SetLineColor(0);
//             legendChargedRatio2050->SetTextFont(42);
//             legendChargedRatio2050->SetTextSize(textSize);
//             legendChargedRatio2050->SetMargin(0.17);
//             legendChargedRatio2050->SetHeader(collisionSystemPbPb.Data());
//   //           legendChargedRatio2050->SetNColumns(2);
// //             legendChargedRatio2050->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb2050,Form("#eta/#pi^{0},   %s",cent2050.Data()),"fp");
//             legendChargedRatio2050->AddEntry(graphPCMEtaToPi0RatioSysErrPbPb2050,Form("#eta/#pi^{0},   %s",cent2050.Data()),"fp");
//             legendChargedRatio2050->AddEntry(graphChargedRatioKaonToPionSys2040,Form("K^{#pm}/#pi^{#pm}, %s",cent2040.Data()),"fp");
//             legendChargedRatio2050->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
//             legendChargedRatio2050->Draw();
//
//             if(thesisPlotting) thesisLabel3->Draw();
//             graphPCMEtaToPi0RatioSysErrPbPb2050->Draw("E2same");
//             graphPCMEtaToPi0RatioPbPb2050->Draw("p,same");
//
//         canvasEtatoPi0combo->SaveAs(Form("%s/EtatoPi0RatioCombined_withKaonsToPions_2050.%s",outputDir.Data(),suffix.Data()));
//
//
//       }
//
	// *******************************************************************************************************
	// ************************** 			RCP			**************************************
	// *******************************************************************************************************   
/*	
	TGraphAsymmErrors *graphPi0CentralInvYield = new TGraphAsymmErrors(histoPCMPi0CorrectedSpecPbPbLHC11h0010);
	TGraphAsymmErrors *graphPi0SemicentralInvYield = new TGraphAsymmErrors(histoPCMPi0CorrectedSpecPbPbLHC11h2050);

	TGraphAsymmErrors *graphEtaCentralInvYield = new TGraphAsymmErrors(histoPCMEtaCorrectedSpecPbPbLHC11h0010);
	TGraphAsymmErrors *graphEtaSemicentralInvYield = new TGraphAsymmErrors(histoPCMEtaCorrectedSpecPbPbLHC11h2050);
	
	TGraphAsymmErrors* graphPi0RCP;
	TGraphAsymmErrors* graphPi0RCPSys;
	TGraphAsymmErrors* graphEtaRCP;
	TGraphAsymmErrors* graphEtaRCPSys;
	
	Double_t NcollC = 1500; 
	Double_t NcollP = 349.3;
	
	CalcRcp( graphPi0CentralInvYield, graphPi0SemicentralInvYield, graphPCMPi0CorrectedSpecSysPbPbLHC11h0010, graphPCMPi0CorrectedSpecSysPbPbLHC11h2050, 
			 &graphPi0RCP, &graphPi0RCPSys, NcollC, NcollP, "Pi0");

	CalcRcp( graphEtaCentralInvYield, graphEtaSemicentralInvYield, graphPCMEtaCorrectedSpecSysPbPbLHC11h0010,graphPCMEtaCorrectedSpecSysPbPbLHC11h2050, 
			 &graphEtaRCP, &graphEtaRCPSys, NcollC, NcollP, "Eta");
	
    TCanvas* canvasRCP = new TCanvas("canvasRCP","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRCP,  0.13, 0.01, 0.015, 0.08);
		
	TH2F * histo2DRCP;
	histo2DRCP = new TH2F("histo2DRCP","histo2DRCP",1000,0.,14.5,1000,0,2.);
	SetStyleHistoTH2ForGraphs(histo2DRCP, "#it{p}_{T} (GeV/#it{c})","#it{R}_{CP}", 0.03,0.04, 0.03,0.04, 0.83,1.4);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{d#it{N}^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{N_{evt}^{pp}}#frac{d#it{N}^{pp}}{d#it{p}_{T} d#it{y}}}
// 	histo2DRCP->GetXaxis()->SetLabelOffset(0.01);
// 	histo2DRCP->GetYaxis()->SetLabelOffset(0.01);
	histo2DRCP->DrawCopy(); 

	graphPi0RCP->Print();
	graphPi0RCPSys->Print();
	graphEtaRCP->Print();
	graphEtaRCPSys->Print();

	DrawGammaSetMarkerTGraphAsym(graphPi0RCP, 20,1, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RCPSys,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
	graphPi0RCPSys->Draw("2same");
	graphPi0RCP->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphEtaRCP, 20,1, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRCPSys,20,1, kRed, kRed, 2, kTRUE, kRed-9);
	graphEtaRCPSys->Draw("2same");
	graphEtaRCP->Draw("p,same");
	
	TLegend* legendRCP = new TLegend(0.5,0.75,0.9,0.9);  //0.16,0.05,0.73,0.2);
	legendRCP->SetFillColor(0);
	legendRCP->SetFillStyle(0);
	legendRCP->SetLineColor(0);
	legendRCP->SetTextSize(0.035);
// 	legendRCP->SetNColumns(2);
	legendRCP->SetMargin(0.2);
	legendRCP->AddEntry(graphPi0RCPSys,"PCM measured #pi^{0} R_{CP} scaled by N_{coll}"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h0010
	legendRCP->AddEntry(graphEtaRCPSys,"PCM measured #eta R_{CP} scaled by N_{coll}"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h2050
	
	legendRCP->Draw();
	canvasRCP->Update();
	canvasRCP->Print(Form("%s/RCP_LHC11h.%s",outputDir.Data(),suffix.Data()));*/

	
	// *******************************************************************************************************
	// *****************************			RAA 				******************************************
	// *******************************************************************************************************   
    TCanvas* canvasRAA = new TCanvas("canvasRAA","",200,10,1200,1100);  //200,10,1350,900);  // gives the page size
        DrawGammaCanvasSettings( canvasRAA, 0.08, 0.02, 0.035, 0.09);

    TH2F * histo2DRAA = new TH2F("histo2DRAA","histo2DRAA",11000,0.,70.,1000,-0.5,2.);
    SetStyleHistoTH2ForGraphs(histo2DRAA, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,.9);
    histo2DRAA->GetYaxis()->SetRangeUser(0.,1.4);
    histo2DRAA->GetXaxis()->SetRangeUser(0.,15);
    histo2DRAA->Draw("copy");

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0005, 20 ,1, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0005,20 ,1, kRed, kRed, 2, kTRUE, kRed-9);
// 	graphPi0RAASys0005->Draw("2same");
// 	graphPi0RAA0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 20,1, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
// 	graphPi0RAASys2040->Draw("2same");
// 	graphPi0RAA2040->Draw("p,same");

    DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPi0RAA0010,markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
    ProduceGraphAsymmWithoutXErrors(graphPi0RAA0010);

    DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPi0RAA2050,  markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
    ProduceGraphAsymmWithoutXErrors(graphPi0RAA2050);

    DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040,  markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
    ProduceGraphAsymmWithoutXErrors(graphPi0RAA2040);

    graphPi0RAASys2050->Draw("E2same");
    graphPi0RAA2050->Draw("p,same");
    graphPi0RAASys0010->Draw("E2same");
    graphPi0RAA0010->Draw("p,same");

    TLegend* legendPi0RAA = new TLegend(0.12,0.7,0.5,0.88);
    legendPi0RAA->SetFillColor(0);
    legendPi0RAA->SetLineColor(0);
    legendPi0RAA->SetTextFont(42);
    legendPi0RAA->SetTextSize(textSize);
    legendPi0RAA->SetMargin(0.17);
    legendPi0RAA->SetHeader(Form("#pi^{0} %s",collisionSystemPbPb.Data()));
    legendPi0RAA->AddEntry(graphPi0RAASys0010," 0#font[122]{-}10%","pf");
    legendPi0RAA->AddEntry(graphPi0RAASys2050,"20#font[122]{-}50%","pf");
    legendPi0RAA->Draw();

    if(thesisPlotting) thesisLabel3->Draw();
// 	if(thesisPlotting) thesisLabel->Draw();
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0PbPb2760GeV.%s",outputDir.Data(),suffix.Data()));


    canvasRAA->cd();
    histo2DRAA->Draw("copy");
    histo2DRAA->GetXaxis()->SetRangeUser(0.,10.3);

    DrawGammaSetMarkerTGraphAsym(graphEtaRAASys0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphEtaRAA0010,markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
    ProduceGraphAsymmWithoutXErrors(graphEtaRAA0010);

    DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphEtaRAA2050,  markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
    ProduceGraphAsymmWithoutXErrors(graphEtaRAA2050);

    DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphEtaRAA2040,  markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
    ProduceGraphAsymmWithoutXErrors(graphEtaRAA2040);

    graphEtaRAASys2050->Draw("E2same");
    graphEtaRAA2050->Draw("p,same");
    graphEtaRAASys0010->Draw("E2same");
    graphEtaRAA0010->Draw("p,same");


    TLegend* legendEtaRAA = new TLegend(0.12,0.7,0.5,0.88);
    legendEtaRAA->SetFillColor(0);
    legendEtaRAA->SetLineColor(0);
    legendEtaRAA->SetTextFont(42);
    legendEtaRAA->SetTextSize(textSize);
    legendEtaRAA->SetMargin(0.17);
    legendEtaRAA->SetHeader(Form("#eta %s",collisionSystemPbPb.Data()));
    legendEtaRAA->AddEntry(graphEtaRAASys0010," 0#font[122]{-}10%","pf");
    legendEtaRAA->AddEntry(graphEtaRAASys2050,"20#font[122]{-}50%","pf");
    legendEtaRAA->Draw();

    if(thesisPlotting) thesisLabel3->Draw();

    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAEtaPbPb2760GeV.%s",outputDir.Data(),suffix.Data()));


// 	canvasRAA->cd();
// 	histo2DRAA->Draw("copy");
//
// 	DrawGammaSetMarkerTGraphAsym(graphRAAPublished0005, 20,2, kRed, kRed, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished0005,20,2, kRed-9, kRed-9, 2, kTRUE, kRed-9);
// 	graphRAASysPublished0005->SetFillStyle(0);
// // 	graphRAASysPublished0005->Draw("2same");
// // 	graphRAAPublished0005->Draw("p,same");
//
// 	graphRAASysPublished0010->Draw("2same");
// 	graphRAAPublished0010->Draw("p,same");
//
// 	DrawGammaSetMarkerTGraphAsym(graphRAAPublished2040, 20,2, kGreen+2, kGreen+2, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished2040,20,2, kGreen-6, kGreen-6, 2, kTRUE, kGreen-6);
// 	graphRAASysPublished2040->SetFillStyle(0);
// // 	graphRAASysPublished2040->Draw("2same");
// // 	graphRAAPublished2040->Draw("p,same");
//
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0005, 21,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0005,21,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
// 	graphPi0RAASys0005->SetFillStyle(0);
// // 	graphPi0RAASys0005->Draw("2same");
// // 	graphPi0RAA0005->Draw("p,same");
//
// 	graphPi0RAASys0010->Draw("2same");
// 	graphPi0RAA0010->Draw("p,same");
//
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2050, 20,1, kBlue, kBlue, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2050,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
// 	graphPi0RAASys2050->SetFillStyle(3008);
// // 	graphPi0RAASys2050->Draw("2same");
// // 	graphPi0RAA2050->Draw("p,same");
//
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 21,2, kBlue,kBlue, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,21,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
// 	graphPi0RAASys2040->SetFillStyle(0);
// // 	graphPi0RAASys2040->Draw("2same");
// // 	graphPi0RAA2040->Draw("p,same");
//
// 	TLegend* legendRAAPi0AllPCM0010 = new TLegend(0.2,0.75,0.7,0.9);  //0.16,0.05,0.73,0.2);
// 	legendRAAPi0AllPCM0010->SetFillColor(0);
// 	legendRAAPi0AllPCM0010->SetFillStyle(0);
// 	legendRAAPi0AllPCM0010->SetLineColor(0);
// 	legendRAAPi0AllPCM0010->SetTextSize(0.035);
// 	legendRAAPi0AllPCM0010->SetMargin(0.2);
// 	legendRAAPi0AllPCM0010->SetHeader("0#font[122]{-}10% PbPb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
// // 	legendRAAPi0AllPCM0010->AddEntry(graphRAAPublished0005,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108) - 0-5%");
// 	legendRAAPi0AllPCM0010->AddEntry(graphRAAPublished0010,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108)");
// // 	legendRAAPi0AllPCM0010->AddEntry(graphRAAPublished2040,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108) - 20-40%");
// // 	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAA0005,"PCM #pi^{0} R_{AA} (2010) - 0-5%", "pl");
// 	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAA0010,"PCM #pi^{0} R_{AA} (2010)","pl");
// // 	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAA2040,"PCM #pi^{0} R_{AA} (2010) - 20-40%", "pl");
// // 	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAASys2050,"PCM #pi^{0} R_{AA} (2010) - 20#font[122]{-}50%");
// 	legendRAAPi0AllPCM0010->Draw();
//
// 	canvasRAA->Update();
// 	canvasRAA->SaveAs(Form("%s/RAAPi0WithPublished0010.%s",outputDir.Data(),suffix.Data()));
	
	
// 	canvasRAA->cd();
// 	histo2DRAA->Draw("copy");
//
// 	DrawGammaSetMarkerTGraphAsym(graphRAAPublished2040, 20,2, kGreen+2, kGreen+2, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished2040,20,2, kGreen-6, kGreen-6, 2, kTRUE, kGreen-6);
// 	graphRAASysPublished2040->SetFillStyle(0);
// 	graphRAASysPublished2040->Draw("2same");
// 	graphRAAPublished2040->Draw("p,same");
//
//
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 21,2, kBlue,kBlue, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,21,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
// 	graphPi0RAASys2040->SetFillStyle(0);
// 	graphPi0RAASys2040->Draw("2same");
// 	graphPi0RAA2040->Draw("p,same");
//
// 	TLegend* legendRAAPi0AllPCM2040 = new TLegend(0.2,0.75,0.7,0.9);  //0.16,0.05,0.73,0.2);
// 	legendRAAPi0AllPCM2040->SetFillColor(0);
// 	legendRAAPi0AllPCM2040->SetFillStyle(0);
// 	legendRAAPi0AllPCM2040->SetLineColor(0);
// 	legendRAAPi0AllPCM2040->SetTextSize(0.035);
// 	legendRAAPi0AllPCM2040->SetMargin(0.2);
// 	legendRAAPi0AllPCM2040->SetHeader("20-40% PbPb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
// 	legendRAAPi0AllPCM2040->AddEntry(graphRAAPublished2040,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108)");
// 	legendRAAPi0AllPCM2040->AddEntry(graphPi0RAA2040,"PCM #pi^{0} R_{AA} (2010)", "pl");
// 	legendRAAPi0AllPCM2040->Draw();
//
// 	canvasRAA->Update();
// 	canvasRAA->SaveAs(Form("%s/RAAPi0WithPublished2040.%s",outputDir.Data(),suffix.Data()));
	


//==========================================================================================
//     TString folderComparison = Form("%s/ComparisonPublished",outputDir.Data());
//     gSystem->Exec("mkdir -p "+folderComparison);
//
//     TCanvas* canvasRatiowithPub = new TCanvas("canvasRatiowithPub","",200,10,1350,900);  // gives the page size
//     DrawGammaCanvasSettings( canvasRatiowithPub, 0.09, 0.01, 0.015, 0.115);
//
//     TH2D *histo2DCompCombinedRatioPub = new TH2D("histo2DCompCombinedRatioPub", "histo2DCompCombinedRatioPub", 20,0.1,12.01,1000.,-0.4,2.);
//     SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatioPub, "#it{p}_{T} (GeV/#it{c})","#pi^{0} #frac{2010}{2011}", 0.035,0.04, 0.035,0.04, 1.,1.);
//     histo2DCompCombinedRatioPub->GetYaxis()->SetRangeUser(0.,2.5);
//     histo2DCompCombinedRatioPub->Draw("copy");
//
//     TGraphErrors* grapha    = NULL;
//     TGraphErrors* graphb    = NULL;
//     TGraphErrors* graphc    = NULL;
//     TGraphErrors* graphd     = NULL;
//     TGraphErrors* graph1    = NULL;
//     TGraphErrors* graph2     = NULL;
//     TGraphErrors* graph3    = NULL;
//     TGraphErrors* graph4     = NULL;
//
//
//     TGraphErrors* RatioPi0YieldsPublishedTo2011_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMPi0CorrectedSpecPbPbLHC10h0010, graphPCMPi0CorrectedSpecSysPbPbLHC10h0010,
//                                                                                                      histoPCMPi0CorrectedSpecPbPbLHC11h0010, graphPCMPi0CorrectedSpecSysPbPbLHC11h0010,
//                                                                                                      kTRUE,  kTRUE,
//                                                                                                      &graphc, &graphd,
//                                                                                                      &grapha, &graphb ) ;
// //     RatioPi0YieldsPublishedTo2011_0010->Print();
//
//     TGraphErrors* RatioPi0YieldsPublishedTo2011_2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMPi0CorrectedSpecPbPbLHC10h2040, graphPCMPi0CorrectedSpecSysPbPbLHC10h2040,
//                                                                                                      histoPCMPi0CorrectedSpecPbPbLHC11h2040, graphPCMPi0CorrectedSpecSysPbPbLHC11h2040,
//                                                                                                      kTRUE,  kTRUE,
//                                                                                                      &graph3, &graph4,
//                                                                                                      &graph1, &graph2 ) ;
// //     RatioPi0YieldsPublishedTo2011_2040->Print();
//
//     DrawGammaSetMarkerTGraphErr(RatioPi0YieldsPublishedTo2011_0010, 20., 2., kRed, kRed);
//     DrawGammaSetMarkerTGraphErr(RatioPi0YieldsPublishedTo2011_2040, 20., 2., kAzure+1, kAzure+1);
//
//     RatioPi0YieldsPublishedTo2011_0010->Draw("same,p");
//     RatioPi0YieldsPublishedTo2011_2040->Draw("same,p");
//
//     TLegend* legendRatio1 = new TLegend(0.15,0.13,0.32,0.32);
//     legendRatio1->SetTextSize(0.035);
//     legendRatio1->SetFillColor(0);
//     legendRatio1->SetFillStyle(0);
//     legendRatio1->SetBorderSize(0);
//     legendRatio1->AddEntry(RatioPi0YieldsPublishedTo2011_0010," 0#font[122]{-}10%","pl");
//     legendRatio1->AddEntry(RatioPi0YieldsPublishedTo2011_2040," 20-40%","pl");
//     legendRatio1->Draw();
//
//     canvasRatiowithPub->SaveAs(Form("%s/Pi0yields_RatioWithPublished.%s",folderComparison.Data(),suffix.Data()));
//
//     canvasRatiowithPub->cd();
//     histo2DCompCombinedRatioPub->Draw("copy");
//
//     grapha    = NULL;
//     graphb    = NULL;
//     graphc    = NULL;
//     graphd     = NULL;
//     graph1    = NULL;
//     graph2     = NULL;
//     graph3    = NULL;
//     graph4     = NULL;
//
//     TGraphErrors* RatioPi0RawYieldsPublishedTo2011_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPi0RawYieldPbPbLHC10h0010, histoPi0RawYieldPbPbLHC10h0010,
//                                                                                                      histoPi0RawYieldPbPbLHC11h0010, histoPi0RawYieldPbPbLHC11h0010,
//                                                                                                      kTRUE,  kTRUE,
//                                                                                                      &graphc, &graphd,
//                                                                                                      &grapha, &graphb ) ;
// //     RatioPi0YieldsPublishedTo2011_0010->Print();
//
//     TGraphErrors* RatioPi0RawYieldsPublishedTo2011_2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPi0RawYieldPbPbLHC10h2040, histoPi0RawYieldPbPbLHC10h2040,
//                                                                                                      histoPi0RawYieldPbPbLHC11h2040, histoPi0RawYieldPbPbLHC11h2040,
//                                                                                                      kTRUE,  kTRUE,
//                                                                                                      &graph3, &graph4,
//                                                                                                      &graph1, &graph2 ) ;
// //     RatioPi0YieldsPublishedTo2011_2040->Print();
//
//     DrawGammaSetMarkerTGraphErr(RatioPi0RawYieldsPublishedTo2011_0010, 20., 2., kRed, kRed);
//     DrawGammaSetMarkerTGraphErr(RatioPi0RawYieldsPublishedTo2011_2040, 20., 2., kAzure+1, kAzure+1);
//
//     RatioPi0RawYieldsPublishedTo2011_0010->Draw("same,p");
//     RatioPi0RawYieldsPublishedTo2011_2040->Draw("same,p");
//
//     legendRatio1->Draw();
//
//     canvasRatiowithPub->SaveAs(Form("%s/Pi0Raw_RatioWithPublished.%s",folderComparison.Data(),suffix.Data()));
//
//     canvasRatiowithPub->cd();
// 	histo2DCompCombinedRatioPub->Draw("copy");
//
//     grapha    = NULL;
//     graphb    = NULL;
//     graphc    = NULL;
//     graphd     = NULL;
//     graph1    = NULL;
//     graph2     = NULL;
//     graph3    = NULL;
//     graph4     = NULL;
//
//     TGraphErrors* RatioPi0EffiPublishedTo2011_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPi0TrueEffiPtPbPbLHC10h0010, histoPi0TrueEffiPtPbPbLHC10h0010,
//                                                                                                      histoPi0TrueEffiPtPbPbLHC11h0010, histoPi0TrueEffiPtPbPbLHC11h0010,
//                                                                                                      kTRUE,  kTRUE,
//                                                                                                      &graphc, &graphd,
//                                                                                                      &grapha, &graphb ) ;
// //     RatioPi0YieldsPublishedTo2011_0010->Print();
//
//     TGraphErrors* RatioPi0EffiPublishedTo2011_2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPi0TrueEffiPtPbPbLHC10h2040, histoPi0TrueEffiPtPbPbLHC10h2040,
//                                                                                                      histoPi0TrueEffiPtPbPbLHC11h2040, histoPi0TrueEffiPtPbPbLHC11h2040,
//                                                                                                      kTRUE,  kTRUE,
//                                                                                                      &graph3, &graph4,
//                                                                                                      &graph1, &graph2 ) ;
// //     RatioPi0YieldsPublishedTo2011_2040->Print();
//
//     DrawGammaSetMarkerTGraphErr(RatioPi0EffiPublishedTo2011_0010, 20., 2., kRed, kRed);
//     DrawGammaSetMarkerTGraphErr(RatioPi0EffiPublishedTo2011_2040, 20., 2., kAzure+1, kAzure+1);
//
//     RatioPi0EffiPublishedTo2011_0010->Draw("same,p");
//     RatioPi0EffiPublishedTo2011_2040->Draw("same,p");
//
//     legendRatio1->Draw();
//
//     canvasRatiowithPub->SaveAs(Form("%s/Pi0Eff_RatioWithPublished.%s",folderComparison.Data(),suffix.Data()));
//
// 	TCanvas* canvasRatioRaa = new TCanvas("canvasRatioRaa","",200,10,1350,900);  // gives the page size
// 	DrawGammaCanvasSettings( canvasRatioRaa, 0.09, 0.01, 0.015, 0.115);
//
// 	TH2F * histo2DCompCombinedRatio2;
// 	histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.,12.,1000,0.2,2);
// 	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0} #it{R}_{AA} #frac{2010}{2011}", 0.03,0.04, 0.03,0.04, 0.83,1);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{d#it{N}^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{N_{evt}^{pp}}#frac{d#it{N}^{pp}}{d#it{p}_{T} d#it{y}}}
// 	histo2DCompCombinedRatio2->DrawCopy();
//
//
// 	TGraphErrors* graphRaaStatPi0Pt2010_0010 	= NULL;
// 	TGraphErrors* graphRaaSysPi0Pt2010_0010 	= NULL;
// 	TGraphErrors* graphRaaStatPi0Pt2011_0010	= NULL;
// 	TGraphErrors* graphRaaSysPi0Pt2011_0010 	= NULL;
// 	TGraphErrors* graphRaaStatPi0Pt2010_2040 	= NULL;
// 	TGraphErrors* graphRaaSysPi0Pt2010_2040 	= NULL;
// 	TGraphErrors* graphRaaStatPi0Pt2011_2040 	= NULL;
// 	TGraphErrors* graphRaaSysPi0Pt2011_2040 	= NULL;
//
// 	TGraphErrors* RatioPi0RaaPublishedTo2011_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphRAAPublished0010, graphRAASysPublished0010,
// 																									 graphPi0RAA0010, graphPi0RAASys0010,
// 																									 kTRUE,  kTRUE,
// 																									 &graphRaaStatPi0Pt2011_0010, &graphRaaSysPi0Pt2011_0010,
// 																									 &graphRaaStatPi0Pt2010_0010, &graphRaaSysPi0Pt2010_0010 ) ;
// // 	RatioPi0RaaPublishedTo2011_0010->Print();
//
// 	TGraphErrors* RatioPi0RaaPublishedTo2011_2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphRAAPublished2040, graphRAASysPublished2040,
// 																									 graphPi0RAA2040, graphPi0RAASys2040,
// 																									 kTRUE,  kTRUE,
// 																									 &graphRaaStatPi0Pt2011_2040, &graphRaaSysPi0Pt2011_2040,
// 																									 &graphRaaStatPi0Pt2010_2040, &graphRaaSysPi0Pt2010_2040 ) ;
// // 	RatioPi0RaaPublishedTo2011_2040->Print();
//
// 	DrawGammaSetMarkerTGraphErr(RatioPi0RaaPublishedTo2011_0010, 20., 1.5, kRed, kRed);
// 	DrawGammaSetMarkerTGraphErr(RatioPi0RaaPublishedTo2011_2040, 20., 1.5, kAzure+1, kAzure+1);
//
// 	RatioPi0RaaPublishedTo2011_0010->Draw("same,p");
// 	RatioPi0RaaPublishedTo2011_2040->Draw("same,p");
//
//
// 	TLegend* legendRatio = new TLegend(0.25,0.13,0.32,0.25);
// 	legendRatio->SetTextSize(0.035);
// 	legendRatio->SetFillColor(0);
// 	legendRatio->SetFillStyle(0);
// 	legendRatio->SetBorderSize(0);
// 	legendRatio->AddEntry(RatioPi0RaaPublishedTo2011_0010," 0#font[122]{-}10%","pl");
// 	legendRatio->AddEntry(RatioPi0RaaPublishedTo2011_2040," 20-40%","pl");
// 	legendRatio->Draw();
//
// 	canvasRatioRaa->SaveAs(Form("%s/Pi0RAA_RatioWithPublished.%s",folderComparison.Data(),suffix.Data()));
//==========================================================================================
	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");

	//Comparison with charged kaons and pions:
	//PbPb276.fullpT.RATIOS.20140329.root RAA_Pion_08052014.root RAA_Kaon_08052014.root

    while(graphChargedPionRAA0010->GetX()[0]<=0.2){
		graphChargedPionRAA0010->RemovePoint(0);
		graphChargedPionRAASys0010->RemovePoint(0);
        graphChargedPionRAA2040->RemovePoint(0);
        graphChargedPionRAASys2040->RemovePoint(0);
	}

    DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0010, markerStylePbPb0010,markerSizePbPb0010, kGray+1 , kGray+1, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0010,  markerStylePbPb0010,markerSizePbPb0010, kGray+1 , kGray+1);
	graphChargedPionRAASys0010->Draw("2same");
	graphChargedPionRAA0010->Draw("p,same");

// 	DrawGammaSetMarkerTGraphAsym(graphRAAPublished0010, 20,2, kRed, kRed, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished0010,20,2, kRed-9, kRed-9, 2, kTRUE, kRed-9);
// 	graphRAASysPublished0010->SetFillStyle(0);
// 	graphRAASysPublished0010->Draw("2same");
// 	graphRAAPublished0010->Draw("p,same");

    graphPi0RAASys0010->Draw("2same");
	graphPi0RAA0010->Draw("p,same"); 

		
// 	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040, 20,1, kGreen+2, kGreen+2, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040,20,1, kGreen+2, kGreen+2, 2, kTRUE, kGreen-6);
// 	graphChargedKaonRAASys2040->Draw("2same");
// 	graphChargedKaonRAA2040->Draw("p,same"); 
// 		
// 	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA2040, 20,1, kBlue, kBlue, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys2040,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
// 	graphPi0RAASys2050->SetFillStyle(3008);
// 	graphChargedPionRAASys2040->Draw("2same");
// 	graphChargedPionRAA2040->Draw("p,same");
// 
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 20,1, kTeal-6, kTeal-6, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,20,1, kTeal-6, kTeal-6, 2, kTRUE, kTeal-9);
// 	graphPi0RAASys2040->SetFillStyle(3001);
// 	graphPi0RAASys2040->Draw("2same");
// 	graphPi0RAA2040->Draw("p,same");
	
    TLegend* legendPi0AndChargedRAA_0010 = new TLegend(0.12,0.7,0.5,0.88);
    legendPi0AndChargedRAA_0010->SetFillColor(0);
    legendPi0AndChargedRAA_0010->SetLineColor(0);
    legendPi0AndChargedRAA_0010->SetTextFont(42);
    legendPi0AndChargedRAA_0010->SetTextSize(textSize);
    legendPi0AndChargedRAA_0010->SetMargin(0.17);
    legendPi0AndChargedRAA_0010->SetHeader(collisionSystemPbPb0010.Data());
    legendPi0AndChargedRAA_0010->AddEntry(graphPi0RAASys0010,"#pi^{0}","fp");
    legendPi0AndChargedRAA_0010->AddEntry(graphChargedPionRAA0010,"#pi^{#pm}","fp");
    legendPi0AndChargedRAA_0010->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
    legendPi0AndChargedRAA_0010->Draw();
    legendPi0AndChargedRAA_0010->Draw();

	if(thesisPlotting) thesisLabel3->Draw();
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0PbPb2760GeV_WithChargedPions0010.%s",outputDir.Data(),suffix.Data()));

	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	
    DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys2040, markerStylePbPb2040,markerSizePbPb2040, kGray+1 , kGray+1, 2, kTRUE);
    graphChargedPionRAASys2040->Draw("2same");
    DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA2040,  markerStylePbPb2040,markerSizePbPb2040, kGray+1 , kGray+1);
    graphChargedPionRAA2040->Draw("p,same");

    graphPi0RAASys2040->Draw("2same");
    graphPi0RAA2040->Draw("p,same");

    TLegend* legendPi0AndChargedRAA_2040 = new TLegend(0.12,0.7,0.5,0.88);
    legendPi0AndChargedRAA_2040->SetFillColor(0);
    legendPi0AndChargedRAA_2040->SetLineColor(0);
    legendPi0AndChargedRAA_2040->SetTextFont(42);
    legendPi0AndChargedRAA_2040->SetTextSize(textSize);
    legendPi0AndChargedRAA_2040->SetMargin(0.17);
    legendPi0AndChargedRAA_2040->SetHeader(collisionSystemPbPb2040.Data());
    legendPi0AndChargedRAA_2040->AddEntry(graphPi0RAASys2040,"#pi^{0}","fp");
    legendPi0AndChargedRAA_2040->AddEntry(graphChargedPionRAA2040,"#pi^{#pm}","fp");
    legendPi0AndChargedRAA_2040->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
    legendPi0AndChargedRAA_2040->Draw();
    legendPi0AndChargedRAA_2040->Draw();

	if(thesisPlotting) thesisLabel3->Draw();
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0PbPb2760GeV_WithChargedPions2040.%s",outputDir.Data(),suffix.Data()));

    canvasRAA->cd();
    histo2DRAA->Draw("copy");

    graphChargedPionRAASys2040->Draw("2same");
    graphChargedPionRAA2040->Draw("p,same");

    graphPi0RAASys2050->Draw("2same");
    graphPi0RAA2050->Draw("p,same");

    TLegend* legendPi0AndChargedRAA_2050 = new TLegend(0.12,0.7,0.5,0.88);
    legendPi0AndChargedRAA_2050->SetFillColor(0);
    legendPi0AndChargedRAA_2050->SetLineColor(0);
    legendPi0AndChargedRAA_2050->SetTextFont(42);
    legendPi0AndChargedRAA_2050->SetTextSize(textSize);
    legendPi0AndChargedRAA_2050->SetMargin(0.17);
    legendPi0AndChargedRAA_2050->SetHeader(collisionSystemPbPb.Data());
    legendPi0AndChargedRAA_2050->AddEntry(graphPi0RAASys2050,Form("#pi^{0} %s",cent2050.Data()),"fp");
    legendPi0AndChargedRAA_2050->AddEntry(graphChargedPionRAA2040,Form("#pi^{#pm} %s",cent2040.Data()),"fp");
    legendPi0AndChargedRAA_2050->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
    legendPi0AndChargedRAA_2050->Draw();
    legendPi0AndChargedRAA_2050->Draw();

    if(thesisPlotting) thesisLabel3->Draw();
    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAPi0PbPb2760GeV_WithChargedPions2050.%s",outputDir.Data(),suffix.Data()));




    canvasRAA->cd();
    histo2DRAA->Draw("copy");

    while(graphChargedKaonRAA0010->GetX()[0]<0.2){
      graphChargedKaonRAA0010->RemovePoint(0);
      graphChargedKaonRAASys0010->RemovePoint(0);
      graphChargedKaonRAA2040->RemovePoint(0);
      graphChargedKaonRAASys2040->RemovePoint(0);
    }


    DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0010, markerStylePbPb0010,markerSizePbPb0010, kGray+1 , kGray+1, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0010,  markerStylePbPb0010,markerSizePbPb0010, kGray+1 , kGray+1);
    graphChargedKaonRAASys0010->Draw("2same");
    graphChargedKaonRAA0010->Draw("p,same");

    graphEtaRAASys0010->Draw("2same");
    graphEtaRAA0010->Draw("p,same");

    TLegend* legendEtaAndChargedRAA_0010 = new TLegend(0.12,0.7,0.5,0.88);
    legendEtaAndChargedRAA_0010->SetFillColor(0);
    legendEtaAndChargedRAA_0010->SetLineColor(0);
    legendEtaAndChargedRAA_0010->SetTextFont(42);
    legendEtaAndChargedRAA_0010->SetTextSize(textSize);
    legendEtaAndChargedRAA_0010->SetMargin(0.17);
    legendEtaAndChargedRAA_0010->SetHeader(collisionSystemPbPb0010.Data());
    legendEtaAndChargedRAA_0010->AddEntry(graphEtaRAASys0010,"#eta","fp");
    legendEtaAndChargedRAA_0010->AddEntry(graphChargedKaonRAASys0010,"K^{#pm}","fp");
    legendEtaAndChargedRAA_0010->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
    legendEtaAndChargedRAA_0010->Draw();
    legendEtaAndChargedRAA_0010->Draw();

    if(thesisPlotting) thesisLabel3->Draw();

    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAEtaPbPb2760GeV_WithChargedKaons0010.%s",outputDir.Data(),suffix.Data()));

    canvasRAA->cd();
    histo2DRAA->Draw("copy");

    DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040, markerStylePbPb2040,markerSizePbPb2040, kGray+1 , kGray+1, 2, kTRUE);
    DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040,  markerStylePbPb2040,markerSizePbPb2040, kGray+1 , kGray+1);
    graphChargedKaonRAASys2040->Draw("2same");
    graphChargedKaonRAA2040->Draw("p,same");

    graphEtaRAASys2040->Draw("2same");
    graphEtaRAA2040->Draw("p,same");

    TLegend* legendEtaAndChargedRAA_2040 = new TLegend(0.12,0.7,0.5,0.88);
    legendEtaAndChargedRAA_2040->SetFillColor(0);
    legendEtaAndChargedRAA_2040->SetLineColor(0);
    legendEtaAndChargedRAA_2040->SetTextFont(42);
    legendEtaAndChargedRAA_2040->SetTextSize(textSize);
    legendEtaAndChargedRAA_2040->SetMargin(0.17);
    legendEtaAndChargedRAA_2040->SetHeader(collisionSystemPbPb2040.Data());
    legendEtaAndChargedRAA_2040->AddEntry(graphEtaRAASys2040,"#eta","fp");
    legendEtaAndChargedRAA_2040->AddEntry(graphChargedKaonRAASys2040,"K^{#pm}","fp");
    legendEtaAndChargedRAA_2040->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
    legendEtaAndChargedRAA_2040->Draw();

    if(thesisPlotting) thesisLabel3->Draw();

    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAEtaPbPb2760GeV_WithChargedKaons2040.%s",outputDir.Data(),suffix.Data()));

    canvasRAA->cd();
    histo2DRAA->Draw("copy");

    graphChargedKaonRAASys2040->Draw("2same");
    graphChargedKaonRAA2040->Draw("p,same");

    graphEtaRAASys2050->Draw("2same");
    graphEtaRAA2050->Draw("p,same");

    TLegend* legendEtaAndChargedRAA_2050 = new TLegend(0.12,0.7,0.5,0.88);
    legendEtaAndChargedRAA_2050->SetFillColor(0);
    legendEtaAndChargedRAA_2050->SetLineColor(0);
    legendEtaAndChargedRAA_2050->SetTextFont(42);
    legendEtaAndChargedRAA_2050->SetTextSize(textSize);
    legendEtaAndChargedRAA_2050->SetMargin(0.17);
    legendEtaAndChargedRAA_2050->SetHeader(collisionSystemPbPb.Data());
    legendEtaAndChargedRAA_2050->AddEntry(graphEtaRAASys2050,Form("#eta %s",cent2040.Data()),"fp");
    legendEtaAndChargedRAA_2050->AddEntry(graphChargedKaonRAASys2040,Form("K^{#pm} %s",cent2040.Data()),"fp");
    legendEtaAndChargedRAA_2050->AddEntry((TObject*)0,"PLB 736 (2014) 196","");
    legendEtaAndChargedRAA_2050->Draw();

    if(thesisPlotting) thesisLabel3->Draw();

    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAEtaPbPb2760GeV_WithChargedKaons2050.%s",outputDir.Data(),suffix.Data()));


    TBox* boxErrorNorm0010          = CreateBoxConv(kRed-7, 0.25, 1.-normErr0010 , 0.5, 1.+normErr0010);
    TBox* boxErrorNorm2050          = CreateBoxConv(kAzure-4, 0.5, 1.-normErr2050 , 0.75, 1.+normErr2050);
    TBox* boxErrorNorm0010_Single   = CreateBoxConv(colorComb0005Box, 0.2, 1.-normErr0010 , 0.5, 1.+normErr0010);
    TBox* boxErrorNorm2040_Single   = CreateBoxConv(colorCombo2050, 0.2, 1.-normErr2040 , 0.5, 1.+normErr2040);
    TBox* boxErrorNorm2050_Single   = CreateBoxConv(colorCombo2050, 0.2, 1.-normErr2050 , 0.5, 1.+normErr2050);
    TBox* boxErrorNorm0010Only      = CreateBoxConv(kRed-7, 0.25, 1.-normErr0010 , 0.5, 1.+normErr0010);
    TBox* boxErrorNorm2050Only      = CreateBoxConv(kAzure-4, 0.25, 1.-normErr2050 , 0.5, 1.+normErr2050);
    TBox* boxErrorNorm2050EPOnly    = CreateBoxConv(kBlue-7, 0.25, 1.-normErr2050 , 0.5, 1.+normErr2050);

    Int_t textSizeLabelsPixelRAA = 50;
    Double_t marginRAA = 0.14*1200;
    Double_t textsizeLabelsRAA = 0;
    Double_t textsizeFacRAA = 0;
    if (canvasRAA->XtoPixel(canvasRAA->GetX2()) < canvasRAA->YtoPixel(canvasRAA->GetY1())){
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA->XtoPixel(canvasRAA->GetX2()) ;
        textsizeFacRAA = (Double_t)1./canvasRAA->XtoPixel(canvasRAA->GetX2()) ;
    } else {
        textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAA->YtoPixel(canvasRAA->GetY1());
        textsizeFacRAA = (Double_t)1./canvasRAA->YtoPixel(canvasRAA->GetY1());
    }

    canvasRAA->cd();
      TH2F * histo2DRAAcomboPHENIX = new TH2F("histo2DRAAcomboPHENIX","histo2DRAAcomboPHENIX",11000,0.,70.,1000,-0.5,2.);
      SetStyleHistoTH2ForGraphs(histo2DRAAcomboPHENIX, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}",0.035,0.04, 0.035,0.04, 1.,.92);
      histo2DRAAcomboPHENIX->GetYaxis()->SetRangeUser(0.,2.);
      histo2DRAAcomboPHENIX->GetXaxis()->SetRangeUser(0.,20.01);
      histo2DRAAcomboPHENIX->Draw("copy");
      DrawGammaLines(0., 20 , 1, 1 ,1,kGray,2);
      DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
      DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_0010, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
      DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_0010, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);
      DrawGammaSetMarkerTGraphErr(graphWA98_17_3GeVRAA_0013, markerStyleWA98,markerSizeWA98, kGray+2 , kGray+2);

        graphPi0RAASys0010->Draw("E2same");
        graphPi0RAA0010->Draw("p,same");

        TLatex *labelRAAALICEPbPb0010 = new TLatex(0.41,0.9,"#pi^{0} ALICE 0#font[122]{-}10% Pb#font[122]{-}Pb");
        graphPHENIX200GeVRAA_0010->Draw("p,same");
        graphPHENIX39GeVRAA_0010->Draw("p,same");
        graphPHENIX62GeVRAA_0010->Draw("p,same");
        graphWA98_17_3GeVRAA_0013->Draw("p,same");

        SetStyleTLatex( labelRAAALICEPbPb0010, textSize,4);
        labelRAAALICEPbPb0010->Draw();
        TLegend* legendRAASinglePbPb0010 = new TLegend(0.41,0.84,0.65,0.885);
        legendRAASinglePbPb0010->SetFillColor(0);
        legendRAASinglePbPb0010->SetLineColor(0);
        legendRAASinglePbPb0010->SetNColumns(1);
        legendRAASinglePbPb0010->SetTextFont(42);
        legendRAASinglePbPb0010->SetTextSize(textSize);
        legendRAASinglePbPb0010->SetMargin(0.17);
        if(thesisPlotting)
        legendRAASinglePbPb0010->AddEntry(graphPi0RAASys0010,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV (this thesis)","pf");
        else legendRAASinglePbPb0010->AddEntry(graphPi0RAASys0010,"0#font[122]{-}10% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
        legendRAASinglePbPb0010->Draw();

        TLatex *labelRAAPHENIXPbPb0010 = new TLatex(0.41,0.79,"#pi^{0} PHENIX 0#font[122]{-}10% Au#font[122]{-}Au");
        SetStyleTLatex( labelRAAPHENIXPbPb0010, textSize,4);
        labelRAAPHENIXPbPb0010->Draw();

        TLegend* legendRAARHICPbPb0010 = new TLegend(0.41,0.66,0.95,0.78);
        legendRAARHICPbPb0010->SetFillColor(0);
        legendRAARHICPbPb0010->SetLineColor(0);
        legendRAARHICPbPb0010->SetNColumns(2);
        legendRAARHICPbPb0010->SetTextFont(42);
        legendRAARHICPbPb0010->SetMargin(0.17);
        legendRAARHICPbPb0010->SetTextSize(textSize);
        legendRAARHICPbPb0010->AddEntry(graphPHENIX200GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
        legendRAARHICPbPb0010->AddEntry(graphPHENIX62GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
        legendRAARHICPbPb0010->AddEntry(graphPHENIX39GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
        legendRAARHICPbPb0010->Draw();

        TLatex *labelRAAWA98PbPb0010 = new TLatex(0.41,0.61,"#pi^{0} WA98 0#font[122]{-}13% Pb#font[122]{-}Pb");
        SetStyleTLatex( labelRAAWA98PbPb0010, textSize,4);
        labelRAAWA98PbPb0010->Draw();

        TLegend* legendRAASPSPbPb0010 = new TLegend(0.41,0.55,0.95,0.59);
        legendRAASPSPbPb0010->SetFillColor(0);
        legendRAASPSPbPb0010->SetLineColor(0);
        legendRAASPSPbPb0010->SetNColumns(2);
        legendRAASPSPbPb0010->SetTextFont(42);
        legendRAASPSPbPb0010->SetTextSize(textSize);
        legendRAASPSPbPb0010->SetMargin(0.17);
        legendRAASPSPbPb0010->AddEntry(graphWA98_17_3GeVRAA_0013,"#sqrt{#it{s}_{_{NN}}} = 17.3 GeV","p");
        legendRAASPSPbPb0010->Draw();

        boxErrorNorm0010_Single->Draw();
//         if(thesisPlotting)thesisLabel3->Draw();

    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAPi0PbPb2760GeV_WithWorldData0010.%s",outputDir.Data(),suffix.Data()));

    canvasRAA->cd();
    histo2DRAAcomboPHENIX->DrawCopy("");
    DrawGammaLines(0., 20 , 1, 1 ,1,kGray,2);

        DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_2040, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
        DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_2040, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
        DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_2040, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);

        graphPHENIX200GeVRAA_2040->Draw("p,same");
        graphPHENIX39GeVRAA_2040->Draw("p,same");
        graphPHENIX62GeVRAA_2040->Draw("p,same");

//         graphPi0RAASys2050->Draw("E2same");
//         graphPi0RAA2050->Draw("p,same");
        graphPi0RAASys2040->Draw("E2same");
        graphPi0RAA2040->Draw("p,same");

        TLatex *labelRAAALICEPbPb2040 = new TLatex(0.57,0.87,"#pi^{0} ALICE 20#font[122]{-}50% Pb#font[122]{-}Pb ");
        SetStyleTLatex( labelRAAALICEPbPb2040, textSize,4);
        labelRAAALICEPbPb2040->Draw();
        TLegend* legendRAASinglePbPb2040 = new TLegend(0.57,0.81,0.83,0.85);
        legendRAASinglePbPb2040->SetFillColor(0);
        legendRAASinglePbPb2040->SetLineColor(0);
        legendRAASinglePbPb2040->SetNColumns(1);
        legendRAASinglePbPb2040->SetTextFont(42);
        legendRAASinglePbPb2040->SetTextSize(textSize);
        legendRAASinglePbPb2040->SetMargin(0.17);
        legendRAASinglePbPb2040->AddEntry(graphPi0RAASys2040,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
        legendRAASinglePbPb2040->Draw();

        TLatex *labelRAAPHENIXPbPb2040 = new TLatex(0.57,0.75,"#pi^{0} PHENIX 20#font[122]{-}40% Au#font[122]{-}Au");
        SetStyleTLatex( labelRAAPHENIXPbPb2040,textSize,4);
        labelRAAPHENIXPbPb2040->Draw();

        TLegend* legendRAARHICPbPb2040 = new TLegend(0.57,0.55,0.84,0.73);
        legendRAARHICPbPb2040->SetFillColor(0);
        legendRAARHICPbPb2040->SetLineColor(0);
    //  legendRAARHICPbPb2040->SetNColumns(2);
        legendRAARHICPbPb2040->SetTextFont(42);
        legendRAARHICPbPb2040->SetTextSize(textSize);
        legendRAARHICPbPb2040->SetMargin(0.17);
        legendRAARHICPbPb2040->AddEntry(graphPHENIX200GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
        legendRAARHICPbPb2040->AddEntry(graphPHENIX62GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
        legendRAARHICPbPb2040->AddEntry(graphPHENIX39GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
        legendRAARHICPbPb2040->Draw();

        if(thesisPlotting)thesisLabel3->Draw();

     boxErrorNorm2040_Single->Draw();
//         boxErrorNorm2050_Single->Draw();
    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAPi0PbPb2760GeV_WithWorldData2040.%s",outputDir.Data(),suffix.Data()));



    canvasRAA->cd();
    histo2DRAAcomboPHENIX->GetYaxis()->SetRangeUser(0.,1.4);
    histo2DRAAcomboPHENIX->GetXaxis()->SetRangeUser(0.,18);
    histo2DRAAcomboPHENIX->DrawCopy("");

      DrawGammaLines(0., 18 , 1, 1 ,1,kGray,2);

      graphEtaRAASys0010->Draw("E2same");
      graphEtaRAA0010->Draw("p,same");

      DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
      graphPHENIX200GeVEtaRAA_0010->Draw("p,same");

      TLegend* legendEtaRAAcomp =new TLegend(0.55,0.83,0.8,0.93);
      legendEtaRAAcomp->SetFillColor(0);
      legendEtaRAAcomp->SetLineColor(0);
      legendEtaRAAcomp->SetTextFont(42);
      legendEtaRAAcomp->SetMargin(0.17);
      legendEtaRAAcomp->SetTextSize(textSize);
      legendEtaRAAcomp->SetHeader("#eta ALICE 0#font[122]{-}10% Pb#font[122]{-}Pb");
      legendEtaRAAcomp->AddEntry(graphEtaRAASys0010,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV", "pf");
      legendEtaRAAcomp->Draw();
      TLegend* legendEtaRAAcompPH =new TLegend(0.55,0.72,0.8,0.82);
      legendEtaRAAcompPH->SetFillColor(0);
      legendEtaRAAcompPH->SetLineColor(0);
      legendEtaRAAcompPH->SetTextFont(42);
      legendEtaRAAcompPH->SetMargin(0.17);
      legendEtaRAAcompPH->SetTextSize(textSize);
      legendEtaRAAcompPH->SetHeader("#eta PHENIX 0#font[122]{-}10% Au#font[122]{-}Au");
      legendEtaRAAcompPH->AddEntry(graphPHENIX200GeVEtaRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
      legendEtaRAAcompPH->Draw();

      boxErrorNorm0010_Single->Draw();
        if(thesisPlotting)thesisLabel3->Draw();

    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAEtaPbPb2760GeV_WithWorldData0010.%s",outputDir.Data(),suffix.Data()));


    canvasRAA->cd();
    histo2DRAAcomboPHENIX->DrawCopy("");
      DrawGammaLines(0., 18 , 1, 1 ,1,kGray,2);

      DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2040, 25,2, colorPhenix,colorPhenix, 0.1, kFALSE);
      DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2040,24,2, colorPhenix,colorPhenix, 1, kTRUE,colorPhenix);
      graphPHENIX200GeVEtaRAA_2040->SetFillStyle(3003);
//       graphPHENIX200GeVEtaRAA_2040->Draw("2same");
      graphPHENIX200GeVEtaRAA_2040->Draw("p,same");

      DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2060, 25,2, colorPhenix,colorPhenix, 0.1, kFALSE);
      DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVEtaRAA_2060,25,2, colorPhenix,colorPhenix, 1, kTRUE, colorPhenix);
      graphPHENIX200GeVEtaRAA_2060->SetFillStyle(3002);
//       graphPHENIX200GeVEtaRAA_2060->Draw("2same");
      graphPHENIX200GeVEtaRAA_2060->Draw("p,same");

      graphEtaRAASys2040->Draw("E2same");
      graphEtaRAA2040->Draw("p,same");

      TLegend* legendEtaRAAcomp1 = new TLegend(0.52,0.83,0.8,0.93);
      legendEtaRAAcomp1->SetFillColor(0);
      legendEtaRAAcomp1->SetLineColor(0);
      legendEtaRAAcomp1->SetTextFont(42);
      legendEtaRAAcomp1->SetMargin(0.17);
      legendEtaRAAcomp1->SetTextSize(textSize);
      legendEtaRAAcomp1->SetHeader("#eta ALICE Pb-Pb #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
      legendEtaRAAcomp1->AddEntry(graphEtaRAASys2040,cent2050.Data(),"pf");
      legendEtaRAAcomp1->Draw();

      TLegend* legendEtaRAAcomp2 = new TLegend(0.52,0.72,0.9,0.82);
      legendEtaRAAcomp2->SetFillColor(0);
      legendEtaRAAcomp2->SetLineColor(0);
      legendEtaRAAcomp2->SetTextFont(42);
      legendEtaRAAcomp2->SetMargin(0.17);
      legendEtaRAAcomp2->SetTextSize(textSize);
      legendEtaRAAcomp2->SetNColumns(2);
      legendEtaRAAcomp2->SetHeader("#eta PHENIX Au-Au #sqrt{#it{s}_{_{NN}}} = 200 GeV");
      legendEtaRAAcomp2->AddEntry(graphPHENIX200GeVEtaRAA_2040,"20#font[122]{-}40%","p");
      legendEtaRAAcomp2->AddEntry(graphPHENIX200GeVEtaRAA_2060,"20#font[122]{-}60%","p");
      legendEtaRAAcomp2->Draw();

      //       boxErrorNorm2050Only->Draw();
     boxErrorNorm2040_Single->Draw();
        if(thesisPlotting)thesisLabel3->Draw();

    canvasRAA->Update();
    canvasRAA->SaveAs(Form("%s/RAAEtaPbPb2760GeV_WithWorldData2040.%s",outputDir.Data(),suffix.Data()));





	TString  nameFinalResDat = Form("%s/FitResultsMC.dat",outputDir.Data());
	TString forOutput;
	fstream  fileFinalResults;
	fileFinalResults.open(nameFinalResDat.Data(), ios::out);
	  
	TCanvas* canvasSpectraMC = new TCanvas("canvasSpectraMC","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraMC,  0.13, 0.01, 0.015, 0.08);
	canvasSpectraMC->SetLogy();
	canvasSpectraMC->SetLogx();
	TH2F * histo2DSpectraMC = new TH2F("histo2DSpectraMC","histo2DSpectraMC",1000,0.03,40.,1000,1e-8,2e4 );
	SetStyleHistoTH2ForGraphs(histo2DSpectraMC, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
	
	TCanvas* canvasFraction2 = new TCanvas("canvasFraction2","",1550,1200);  // gives the page size
	canvasFraction2->SetTickx();
	canvasFraction2->SetTicky();
	canvasFraction2->SetGridx(0);
	canvasFraction2->SetGridy(0);
	canvasFraction2->SetLogy(0);
	canvasFraction2->SetLeftMargin(0.13);
	canvasFraction2->SetRightMargin(0.02);
	canvasFraction2->SetTopMargin(0.02);
	canvasFraction2->SetFillColor(0);
	
	//**********************************************************************************
	//**************************** Pi0 reweighting evalulation 0#font[122]{-}10% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDPi0PbPbLHC11h0010 = NULL;
	if (directoryPi0PbPbLHC11h0010){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldPi0PtPbPbLHC11h0010->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldPi0PtPbPbLHC11h0010->Draw("hist,pe1,same");
		
		histoPCMPi0CorrectedSpecPbPbLHC11h0010->Draw("hist,pe1,same");
		
		fitYieldDataQCDPi0PbPbLHC11h0010 = FitObject("qcd","fitYieldDataQCDPi0PbPbLHC11h0010","Pi0",histoPCMPi0CorrectedSpecPbPbLHC11h0010,0.4,14,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPbLHC11h0010, 1, 1.5, colorCombPbPb6080);
        fitYieldDataQCDPi0PbPbLHC11h0010->SetLineWidth(2);
		fitYieldDataQCDPi0PbPbLHC11h0010->Draw("same");
		
		forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPbLHC11h0010);
		fileFinalResults<< forOutput.Data()<< endl;  
        
        TLegend* legendYieldandFit0010 = new TLegend(0.66,0.81,0.9,0.9);
        legendYieldandFit0010->SetFillColor(0);
        legendYieldandFit0010->SetLineColor(0);
        legendYieldandFit0010->SetTextSize(0.035);
        legendYieldandFit0010->SetMargin(0.2);
        legendYieldandFit0010->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h0010,"Data","p");
        legendYieldandFit0010->AddEntry(histoMCYieldPi0PtPbPbLHC11h0010,"MC yield","p");
        legendYieldandFit0010->AddEntry(fitYieldDataQCDPi0PbPbLHC11h0010,"QCD fit to data","l");
        legendYieldandFit0010->Draw();
        
		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Pi0_MCInputSpectraFittedPbPbLHC11h0010.%s",outputDir.Data(),suffix.Data()));

		
// 		histoMCYieldPi0PtPbPbLHC11hAddedSig0010->SetMarkerStyle(markerStylePbPb4060MC);
		TH1D* histoRatioDatatoFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpecPbPbLHC11h0010, fitYieldDataQCDPi0PbPbLHC11h0010);
		TH1D* histoRatioMCtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11h0010, fitYieldDataQCDPi0PbPbLHC11h0010);
// 		TH1D* histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hAddedSig0010, fitYieldDataQCDPi0PbPbLHC11h0010);
		TH1D* histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010 = NULL;
		if (histoMCYieldPi0PtPbPbLHC11hWOWeights0010) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hWOWeights0010, fitYieldDataQCDPi0PbPbLHC11h0010);
		
		canvasFraction2->cd();
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoRatioDatatoFitQCDPbPbLHC11h0010, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPbLHC11h0010,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, -0.5, 8.,
				kTRUE, 0., 13.5);
		histoRatioDatatoFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPbLHC11h0010, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
		histoRatioMCtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		
// 		DrawGammaSetMarker(histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010, markerStylePbPb4060MC,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060);
// 		histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p");
		
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		
		TLegend* legendFit0010 = new TLegend(0.16,0.73,0.4,0.9);
		legendFit0010->SetFillColor(0);
		legendFit0010->SetLineColor(0);
		legendFit0010->SetTextSize(0.035);
		legendFit0010->SetMargin(0.2);
		legendFit0010->AddEntry(histoRatioDatatoFitQCDPbPbLHC11h0010,"Data/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit0010->AddEntry(histoRatioMCtoDataFitQCDPbPbLHC11h0010,"MC weighted/QCD fit to Data (0.4 < pT < 14)","p");
// 		if (runDrawReweighted) legendFit0010->AddEntry(histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010,"MC added sig. weighted/QCD fit to Data (0.4 < pT < 14)","p");
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010) legendFit0010->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010,"MC/QCD fit to Data (0.4 < pT < 14)","p");
		legendFit0010->Draw();
		
		TLatex *labelRatioMCData0010 = new TLatex(0.2,0.92,collisionSystemPbPb0010.Data());
		SetStyleTLatex( labelRatioMCData0010, 0.04,4);
		labelRatioMCData0010->Draw();
		DrawGammaLines(0., 30.,1., 1.,0.1);
		
		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Pi0_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb0010.%s",outputDir.Data(),suffix.Data()));
	}

	//**********************************************************************************
	//**************************** Pi0 reweighting evalulation 0-5% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDPi0PbPbLHC11h0005 = NULL;
	if (directoryPi0PbPbLHC11h0005){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldPi0PtPbPbLHC11h0005->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldPi0PtPbPbLHC11h0005->Draw("hist,pe1,same");

		histoPCMPi0CorrectedSpecPbPbLHC11h0005->Draw("hist,pe1,same");
		
		fitYieldDataQCDPi0PbPbLHC11h0005 = FitObject("qcd","fitYieldDataQCDPi0PbPbLHC11h0005","Pi0",histoPCMPi0CorrectedSpecPbPbLHC11h0005,0.4,14,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPbLHC11h0005, 1, 1.5, colorCombPbPb6080);
		fitYieldDataQCDPi0PbPbLHC11h0005->Draw("same");
		
		forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPbLHC11h0005);
		fileFinalResults<< forOutput.Data()<< endl;  
		
		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Pi0_MCInputSpectraFittedPbPbLHC11h0005.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoRatioDatatoFitQCDPbPbLHC11h0005 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpecPbPbLHC11h0005, fitYieldDataQCDPi0PbPbLHC11h0005);
		TH1D* histoRatioMCtoDataFitQCDPbPbLHC11h0005 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11h0005, fitYieldDataQCDPi0PbPbLHC11h0005);
		TH1D* histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005 = NULL;
		if (histoMCYieldPi0PtPbPbLHC11hWOWeights0005) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hWOWeights0005, fitYieldDataQCDPi0PbPbLHC11h0005);
		
		canvasFraction2->cd();
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoRatioDatatoFitQCDPbPbLHC11h0005, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPbLHC11h0005,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, -0.5, 8.,
				kTRUE, 0., 13.5);
		histoRatioDatatoFitQCDPbPbLHC11h0005->Draw("same,e,p");  

		DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
		histoRatioMCtoDataFitQCDPbPbLHC11h0005->Draw("same,e,p"); 
		
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005->Draw("same,e,p");  
		
		TLegend* legendFit0005 = new TLegend(0.16,0.8,0.4,0.9);
		legendFit0005->SetFillColor(0);
		legendFit0005->SetLineColor(0);
		legendFit0005->SetTextSize(0.035);
		legendFit0005->SetMargin(0.2);
		legendFit0005->AddEntry(histoRatioDatatoFitQCDPbPbLHC11h0005,"Data/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit0005->AddEntry(histoRatioMCtoDataFitQCDPbPbLHC11h0005,"MC weighted/QCD fit to Data (0.4 < pT < 14)","p");
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005) legendFit0005->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0005,"MC/QCD fit to Data (0.4 < pT < 14)","p");
		legendFit0005->Draw();
		
		TLatex *labelRatioMCData0005 = new TLatex(0.2,0.92,collisionSystemPbPb0005.Data());
		SetStyleTLatex( labelRatioMCData0005, 0.04,4);
		labelRatioMCData0005->Draw();
		DrawGammaLines(0., 30.,1., 1.,0.1);
		
		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Pi0_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb0005.%s",outputDir.Data(),suffix.Data()));
	}

	
	//**********************************************************************************
	//**************************** Pi0 reweighting evalulation 5-10% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDPi0PbPbLHC11h0510 = NULL;
	if (directoryPi0PbPbLHC11h0510){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldPi0PtPbPbLHC11h0510->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldPi0PtPbPbLHC11h0510->Draw("hist,pe1,same");

		histoPCMPi0CorrectedSpecPbPbLHC11h0510->Draw("hist,pe1,same");
		
		fitYieldDataQCDPi0PbPbLHC11h0510 = FitObject("qcd","fitYieldDataQCDPi0PbPbLHC11h0510","Pi0",histoPCMPi0CorrectedSpecPbPbLHC11h0510,0.4,14,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPbLHC11h0510, 1, 1.5, colorCombPbPb6080);
		fitYieldDataQCDPi0PbPbLHC11h0510->Draw("same");

		forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPbLHC11h0510);
		fileFinalResults<< forOutput.Data()<< endl;  

		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Pi0_MCInputSpectraFittedPbPbLHC11h0510.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoRatioDatatoFitQCDPbPbLHC11h0510 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpecPbPbLHC11h0510, fitYieldDataQCDPi0PbPbLHC11h0510);
		TH1D* histoRatioMCtoDataFitQCDPbPbLHC11h0510 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11h0510, fitYieldDataQCDPi0PbPbLHC11h0510);
		TH1D* histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510 = NULL;
		if (histoMCYieldPi0PtPbPbLHC11hWOWeights0510) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hWOWeights0510, fitYieldDataQCDPi0PbPbLHC11h0510);

		canvasFraction2->cd();
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoRatioDatatoFitQCDPbPbLHC11h0510, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPbLHC11h0510,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, -0.5, 8.,
				kTRUE, 0., 13.5);
		histoRatioDatatoFitQCDPbPbLHC11h0510->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPbLHC11h0510, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
		histoRatioMCtoDataFitQCDPbPbLHC11h0510->Draw("same,e,p");  
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510->Draw("same,e,p");  
		
		TLegend* legendFit0510 = new TLegend(0.16,0.8,0.4,0.9);
		legendFit0510->SetFillColor(0);
		legendFit0510->SetLineColor(0);
		legendFit0510->SetTextSize(0.035);
		legendFit0510->SetMargin(0.2);
		legendFit0510->AddEntry(histoRatioDatatoFitQCDPbPbLHC11h0510,"Data/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit0510->AddEntry(histoRatioMCtoDataFitQCDPbPbLHC11h0510,"MC weighted/QCD fit to Data (0.4 < pT < 14)","p");
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510) legendFit0510->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0510,"MC/QCD fit to Data (0.4 < pT < 14)","p");
		legendFit0510->Draw();
		
		TLatex *labelRatioMCData0510 = new TLatex(0.2,0.92,collisionSystemPbPb0510.Data());
		SetStyleTLatex( labelRatioMCData0510, 0.04,4);
		labelRatioMCData0510->Draw();
		DrawGammaLines(0., 30.,1., 1.,0.1);
		
		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Pi0_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb0510.%s",outputDir.Data(),suffix.Data()));
	}

	
	//**********************************************************************************
	//**************************** Pi0 reweighting evalulation 20-40% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDPi0PbPbLHC11h2040 = NULL;
	if (directoryPi0PbPbLHC11h2040){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldPi0PtPbPbLHC11h2040->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldPi0PtPbPbLHC11h2040->Draw("hist,pe1,same");

		histoPCMPi0CorrectedSpecPbPbLHC11h2040->Draw("hist,pe1,same");

		fitYieldDataQCDPi0PbPbLHC11h2040 = FitObject("qcd","fitYieldDataQCDPi0PbPbLHC11h2040","Pi0",histoPCMPi0CorrectedSpecPbPbLHC11h2040,0.4,14,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPbLHC11h2040, 1, 1.5, colorCombPbPb6080);
        fitYieldDataQCDPi0PbPbLHC11h2040->SetLineWidth(2);
		fitYieldDataQCDPi0PbPbLHC11h2040->Draw("same");

		forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPbLHC11h2040);
		fileFinalResults<< forOutput.Data()<< endl;  
        
        TLegend* legendYieldandFit2040 = new TLegend(0.66,0.8,0.9,0.9);
        legendYieldandFit2040->SetFillColor(0);
        legendYieldandFit2040->SetLineColor(0);
        legendYieldandFit2040->SetTextSize(0.035);
        legendYieldandFit2040->SetMargin(0.2);
        legendYieldandFit2040->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h2040,"Data","p");
        legendYieldandFit2040->AddEntry(histoMCYieldPi0PtPbPbLHC11h2040,"MC yield","p");
        legendYieldandFit2040->AddEntry(fitYieldDataQCDPi0PbPbLHC11h2040,"QCD fit to data","l");
        legendYieldandFit2040->Draw();

		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Pi0_MCInputSpectraFittedPbPbLHC11h2040.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoRatioDatatoFitQCDPbPbLHC11h2040 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpecPbPbLHC11h2040, fitYieldDataQCDPi0PbPbLHC11h2040);
		TH1D* histoRatioMCtoDataFitQCDPbPbLHC11h2040 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11h2040, fitYieldDataQCDPi0PbPbLHC11h2040);
		TH1D* histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040 = NULL;
		if (histoMCYieldPi0PtPbPbLHC11hWOWeights2040) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hWOWeights2040, fitYieldDataQCDPi0PbPbLHC11h2040);
		
		canvasFraction2->cd();
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoRatioDatatoFitQCDPbPbLHC11h2040, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPbLHC11h2040,
				"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
				kFALSE, 1.5, 0, kTRUE,
				kTRUE, -0.5, 8.,
				kTRUE, 0., 13.5);
		histoRatioDatatoFitQCDPbPbLHC11h2040->Draw("same,e,p");  

		DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb0005 , colorCombPbPb0005);
		histoRatioMCtoDataFitQCDPbPbLHC11h2040->Draw("same,e,p");  
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040->Draw("same,e,p");  
		
		TLegend* legendFit2040 = new TLegend(0.16,0.8,0.4,0.9);
		legendFit2040->SetFillColor(0);
		legendFit2040->SetLineColor(0);
		legendFit2040->SetTextSize(0.035);
		legendFit2040->SetMargin(0.2);
		legendFit2040->AddEntry(histoRatioDatatoFitQCDPbPbLHC11h2040,"Data/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit2040->AddEntry(histoRatioMCtoDataFitQCDPbPbLHC11h2040,"MC weighted/QCD fit to Data (0.4 < pT < 14)","p");
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040) legendFit2040->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2040,"MC/QCD fit to Data (0.4 < pT < 14)","p");
		legendFit2040->Draw();

		TLatex *labelRatioMCData2040 = new TLatex(0.2,0.92,collisionSystemPbPb2040.Data());
		SetStyleTLatex( labelRatioMCData2040, 0.04,4);
		labelRatioMCData2040->Draw();
		DrawGammaLines(0., 30.,1., 1.,0.1);

		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Pi0_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb2040.%s",outputDir.Data(),suffix.Data()));
	}

	//**********************************************************************************
	//**************************** Pi0 reweighting evalulation 20-50 % LHC11h***********
	//**********************************************************************************
	TF1* fitYieldDataQCDPi0PbPbLHC11h2050 = NULL;
	if (directoryPi0PbPbLHC11h2050){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldPi0PtPbPbLHC11h2050->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldPi0PtPbPbLHC11h2050->Draw("hist,pe1,same");

		histoPCMPi0CorrectedSpecPbPbLHC11h2050->Draw("hist,pe1,same");

		fitYieldDataQCDPi0PbPbLHC11h2050 = FitObject("qcd","fitYieldDataQCDPi0PbPbLHC11h2050","Pi0",histoPCMPi0CorrectedSpecPbPbLHC11h2050,0.4,14,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDPi0PbPbLHC11h2050, 1, 1.5, colorCombPbPb6080);
        fitYieldDataQCDPi0PbPbLHC11h2050->SetLineWidth(2);
		fitYieldDataQCDPi0PbPbLHC11h2050->Draw("same");

		forOutput= WriteParameterToFile(fitYieldDataQCDPi0PbPbLHC11h2050);
		fileFinalResults<< forOutput.Data()<< endl;          
        
        TLegend* legendYieldandFit2050 = new TLegend(0.66,0.8,0.9,0.9);
        legendYieldandFit2050->SetFillColor(0);
        legendYieldandFit2050->SetLineColor(0);
        legendYieldandFit2050->SetTextSize(0.035);
        legendYieldandFit2050->SetMargin(0.2);
        legendYieldandFit2050->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h2050,"Data","p");
        legendYieldandFit2050->AddEntry(histoMCYieldPi0PtPbPbLHC11h2050,"MC yield","p");
        legendYieldandFit2050->AddEntry(fitYieldDataQCDPi0PbPbLHC11h2050,"QCD fit to data","l");
        legendYieldandFit2050->Draw();

		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Pi0_MCInputSpectraFittedPbPbLHC11h2050.%s",outputDir.Data(),suffix.Data()));

        histoMCYieldPi0PtPbPbLHC11hAddedSig2050->SetMarkerStyle(markerStylePbPb4060MC);
		TH1D* histoRatioDatatoFitQCDPbPbLHC11h2050 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpecPbPbLHC11h2050, fitYieldDataQCDPi0PbPbLHC11h2050);
		TH1D* histoRatioMCtoDataFitQCDPbPbLHC11h2050 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11h2050, fitYieldDataQCDPi0PbPbLHC11h2050);
		TH1D* histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h2050 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hAddedSig2050, fitYieldDataQCDPi0PbPbLHC11h2050);
		TH1D* histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050 = NULL;
		if (histoMCYieldPi0PtPbPbLHC11hWOWeights2050) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hWOWeights2050, fitYieldDataQCDPi0PbPbLHC11h2050);
		
		canvasFraction2->cd();
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050) DrawGammaSetMarker(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
		DrawGammaSetMarker(histoRatioDatatoFitQCDPbPbLHC11h2050, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoRatioDatatoFitQCDPbPbLHC11h2050,
					"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
					kFALSE, 1.5, 0, kTRUE,
					kTRUE, -0.5, 8.,
					kTRUE, 0., 13.5);
		histoRatioDatatoFitQCDPbPbLHC11h2050->Draw("same,e,p");  

		DrawGammaSetMarker(histoRatioMCtoDataFitQCDPbPbLHC11h2050, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
		histoRatioMCtoDataFitQCDPbPbLHC11h2050->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h2050, markerStylePbPb4060MC,markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);
		histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h2050->Draw("same,e,p");  
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050->Draw("same,e,p");  
		
		TLegend* legendFit2050 = new TLegend(0.16,0.73,0.4,0.9);
		legendFit2050->SetFillColor(0);
		legendFit2050->SetLineColor(0);
		legendFit2050->SetTextSize(0.035);
		legendFit2050->SetMargin(0.2);
		legendFit2050->AddEntry(histoRatioDatatoFitQCDPbPbLHC11h2050,"Data/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit2050->AddEntry(histoRatioMCtoDataFitQCDPbPbLHC11h2050,"MC weighted/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit2050->AddEntry(histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h2050,"MC added sig. weighted/QCD fit to Data (0.4 < pT < 14)","p");
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050) legendFit2050->AddEntry(histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h2050,"MC/QCD fit to Data (0.4 < pT < 14)","p");
		legendFit2050->Draw();
		
		TLatex *labelRatioMCData2050 = new TLatex(0.2,0.92,collisionSystemPbPb2050.Data());
		SetStyleTLatex( labelRatioMCData2050, 0.04,4);
		labelRatioMCData2050->Draw();
		DrawGammaLines(0., 30.,1., 1.,0.1);
		
		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Pi0_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb2050.%s",outputDir.Data(),suffix.Data()));
	}

  
	//**********************************************************************************
	//**************************** Eta reweighting evalulation 0#font[122]{-}10% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDEtaPbPbLHC11h0010 = NULL;
	if (directoryEtaPbPbLHC11h0010){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldEtaPtPbPbLHC11h0010->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldEtaPtPbPbLHC11h0010->Draw("hist,pe1,same");

		histoMCYieldEtaPtPbPbLHC11hAddedSig0010->SetMarkerStyle(markerStylePbPb4060MC);
		histoMCYieldEtaPtPbPbLHC11hAddedSig0010->Draw("hist,pe1,same");

		histoPCMEtaCorrectedSpecPbPbLHC11h0010->Draw("hist,pe1,same");
		
		fitYieldDataQCDEtaPbPbLHC11h0010 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h0010","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h0010,1,10,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDEtaPbPbLHC11h0010, 1, 1.5, colorCombPbPb0010);
		fitYieldDataQCDEtaPbPbLHC11h0010->SetLineWidth(2);
		fitYieldDataQCDEtaPbPbLHC11h0010->Draw("same");
		
		forOutput= WriteParameterToFile(fitYieldDataQCDEtaPbPbLHC11h0010);
		fileFinalResults<< forOutput.Data()<< endl;  
		
		TLegend* legendEtaYieldandFit0010 = new TLegend(0.66,0.81,0.9,0.9);
        legendEtaYieldandFit0010->SetFillColor(0);
        legendEtaYieldandFit0010->SetLineColor(0);
        legendEtaYieldandFit0010->SetTextSize(0.035);
        legendEtaYieldandFit0010->SetMargin(0.2);
        legendEtaYieldandFit0010->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h0010,"Data","p");
        legendEtaYieldandFit0010->AddEntry(histoMCYieldEtaPtPbPbLHC11hAddedSig0010,"MC yield","p");
        legendEtaYieldandFit0010->AddEntry(fitYieldDataQCDEtaPbPbLHC11h0010,"QCD fit to data","l");
        legendEtaYieldandFit0010->Draw();

		
		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Eta_MCInputSpectraFittedPbPbLHC11h0010.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoEtaRatioDatatoFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoPCMEtaCorrectedSpecPbPbLHC11h0010, fitYieldDataQCDEtaPbPbLHC11h0010);
		TH1D* histoEtaRatioMCtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11h0010, fitYieldDataQCDEtaPbPbLHC11h0010);
		TH1D* histoEtaRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11hAddedSig0010, fitYieldDataQCDEtaPbPbLHC11h0010);
		TH1D* histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010 = NULL;
		if (histoMCYieldEtaPtPbPbLHC11hWOWeights0010) histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11hWOWeights0010, fitYieldDataQCDEtaPbPbLHC11h0010);
		
		canvasFraction2->cd();
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoEtaRatioDatatoFitQCDPbPbLHC11h0010, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitQCDPbPbLHC11h0010,
					"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
					kFALSE, 1.5, 0, kTRUE,
					kTRUE, -0.5, 8.,
					kTRUE, 0., 10);
		histoEtaRatioDatatoFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoEtaRatioMCtoDataFitQCDPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
		if (runDrawReweighted) histoEtaRatioMCtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoEtaRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010, markerStylePbPb4060MC,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060);
		if (runDrawReweighted) histoEtaRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010) histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p"); 
		
		TLegend* legendEtaFit0010 = new TLegend(0.16,0.73,0.4,0.9);
		legendEtaFit0010->SetFillColor(0);
		legendEtaFit0010->SetLineColor(0);
		legendEtaFit0010->SetTextSize(0.035);
		legendEtaFit0010->SetMargin(0.2);
		legendEtaFit0010->AddEntry(histoEtaRatioDatatoFitQCDPbPbLHC11h0010,"Data/Tsallis fit to Data (1 < pT < 10)","p");
		if (runDrawReweighted) legendEtaFit0010->AddEntry(histoEtaRatioMCtoDataFitQCDPbPbLHC11h0010,"MC weighted/Tsallis fit to Data (1 < pT < 10)","p");
		if (runDrawReweighted) legendEtaFit0010->AddEntry(histoEtaRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010,"MC added sig. weighted/Tsallis fit to Data (1 < pT < 10)","p");
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010) legendEtaFit0010->AddEntry(histoEtaRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010,"MC/Tsallis fit to Data (1 < pT < 10)","p");
		legendEtaFit0010->Draw();

		TLatex *labelRatioMCData0010 = new TLatex(0.2,0.92,collisionSystemPbPb0010.Data());
		SetStyleTLatex( labelRatioMCData0010, 0.04,4);
		labelRatioMCData0010->Draw();
		DrawGammaLines(0., 30.,1., 1.,0.1);
		
		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Eta_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb0010.%s",outputDir.Data(),suffix.Data()));
	}

	//**********************************************************************************
	//**************************** Eta reweighting evalulation 0-5% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDEtaPbPbLHC11h0005 = NULL;
	if (directoryEtaPbPbLHC11h0005){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy();

		histoMCYieldEtaPtPbPbLHC11h0005->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldEtaPtPbPbLHC11h0005->Draw("hist,pe1,same");

		histoPCMEtaCorrectedSpecPbPbLHC11h0005->Draw("hist,pe1,same");

		fitYieldDataQCDEtaPbPbLHC11h0005 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h0005","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h0005,1.,10,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDEtaPbPbLHC11h0005, 1, 1.5, colorCombPbPb6080);
		fitYieldDataQCDEtaPbPbLHC11h0005->SetLineWidth(2);
		fitYieldDataQCDEtaPbPbLHC11h0005->Draw("same");

		forOutput= WriteParameterToFile(fitYieldDataQCDEtaPbPbLHC11h0005);
		fileFinalResults<< forOutput.Data()<< endl;

		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Eta_MCInputSpectraFittedPbPbLHC11h0005.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoEtaRatioDatatoFitQCDPbPb0005 = CalculateHistoRatioToFit (histoPCMEtaCorrectedSpecPbPbLHC11h0005, fitYieldDataQCDEtaPbPbLHC11h0005);
		TH1D* histoEtaRatioMCtoDataFitQCDPbPb0005 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11h0005, fitYieldDataQCDEtaPbPbLHC11h0005);
		TH1D* histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005 = NULL;
		if (histoMCYieldEtaPtPbPbLHC11hWOWeights0005) histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11hWOWeights0005, fitYieldDataQCDEtaPbPbLHC11h0005);

		canvasFraction2->cd();
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );

		DrawGammaSetMarker(histoEtaRatioDatatoFitQCDPbPb0005, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitQCDPbPb0005,
					"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
					kFALSE, 1.5, 0, kTRUE,
					kTRUE, -0.5, 8.,
					kTRUE, 0., 10);
		histoEtaRatioDatatoFitQCDPbPb0005->Draw("same,e,p");

		DrawGammaSetMarker(histoEtaRatioMCtoDataFitQCDPbPb0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
		if (runDrawReweighted) histoEtaRatioMCtoDataFitQCDPbPb0005->Draw("same,e,p");
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005) histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005->Draw("same,e,p");

		TLegend* legendEtaFit0005 = new TLegend(0.16,0.8,0.4,0.9);
		legendEtaFit0005->SetFillColor(0);
		legendEtaFit0005->SetLineColor(0);
		legendEtaFit0005->SetTextSize(0.035);
		legendEtaFit0005->SetMargin(0.2);
		legendEtaFit0005->AddEntry(histoEtaRatioDatatoFitQCDPbPb0005,"Data/Tsallis fit to Data (1 < pT < 10)","p");
		if (runDrawReweighted) legendEtaFit0005->AddEntry(histoEtaRatioMCtoDataFitQCDPbPb0005,"MC weighted/Tsallis fit to Data (1 < pT < 10)","p");
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005) legendEtaFit0005->AddEntry(histoEtaRatioMCUnweightedtoDataFitQCDPbPb0005,"MC/Tsallis fit to Data (1 < pT < 10)","p");
		legendEtaFit0005->Draw();
		TLatex *labelRatioMCData0005 = new TLatex(0.2,0.92,collisionSystemPbPb0005.Data());
		SetStyleTLatex( labelRatioMCData0005, 0.04,4);
		labelRatioMCData0005->Draw();
		DrawGammaLines(0., 30.,1., 1.,0.1);

		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Eta_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb0005.%s",outputDir.Data(),suffix.Data()));
	}

	
	//**********************************************************************************
	//**************************** Eta reweighting evalulation 5-10% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDEtaPbPbLHC11h0510 = NULL;
	if (directoryEtaPbPbLHC11h0510){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldEtaPtPbPbLHC11h0510->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldEtaPtPbPbLHC11h0510->Draw("hist,pe1,same");

		histoPCMEtaCorrectedSpecPbPbLHC11h0510->Draw("hist,pe1,same");
		
		fitYieldDataQCDEtaPbPbLHC11h0510 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h0510","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h0510,1.,10,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDEtaPbPbLHC11h0510, 1, 1.5, colorCombPbPb6080);
		fitYieldDataQCDEtaPbPbLHC11h0510->Draw("same");
		
		forOutput= WriteParameterToFile(fitYieldDataQCDEtaPbPbLHC11h0510);
		fileFinalResults<< forOutput.Data()<< endl;  
		
		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Eta_MCInputSpectraFittedPbPbLHC11h0510.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoEtaRatioDatatoFitQCDPbPb0510 = CalculateHistoRatioToFit (histoPCMEtaCorrectedSpecPbPbLHC11h0510, fitYieldDataQCDEtaPbPbLHC11h0510);
		TH1D* histoEtaRatioMCtoDataFitQCDPbPb0510 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11h0510, fitYieldDataQCDEtaPbPbLHC11h0510);
		TH1D* histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510 = NULL;
		if (histoMCYieldEtaPtPbPbLHC11hWOWeights0510) histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11hWOWeights0510, fitYieldDataQCDEtaPbPbLHC11h0510);
		
		canvasFraction2->cd();
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoEtaRatioDatatoFitQCDPbPb0510, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitQCDPbPb0510,
					"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
					kFALSE, 1.5, 0, kTRUE,
					kTRUE, -0.5, 8.,
					kTRUE, 0., 10);
		histoEtaRatioDatatoFitQCDPbPb0510->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoEtaRatioMCtoDataFitQCDPbPb0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0005 , colorCombPbPb0005);
		if (runDrawReweighted) histoEtaRatioMCtoDataFitQCDPbPb0510->Draw("same,e,p");  
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510) histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510->Draw("same,e,p");  
		
		TLegend* legendEtaFit0510 = new TLegend(0.16,0.8,0.4,0.9);
		legendEtaFit0510->SetFillColor(0);
		legendEtaFit0510->SetLineColor(0);
		legendEtaFit0510->SetTextSize(0.035);
		legendEtaFit0510->SetMargin(0.2);
		legendEtaFit0510->AddEntry(histoEtaRatioDatatoFitQCDPbPb0510,"Data/Tsallis fit to Data (1 < pT < 10)","p");
		if (runDrawReweighted) legendEtaFit0510->AddEntry(histoEtaRatioMCtoDataFitQCDPbPb0510,"MC weighted/Tsallis fit to Data (1 < pT < 10)","p");
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510) legendEtaFit0510->AddEntry(histoEtaRatioMCUnweightedtoDataFitQCDPbPb0510,"MC/Tsallis fit to Data (1 < pT < 10)","p");
		legendEtaFit0510->Draw();
		TLatex *labelRatioMCData0510 = new TLatex(0.2,0.92,collisionSystemPbPb0510.Data());
		SetStyleTLatex( labelRatioMCData0510, 0.04,4);
		labelRatioMCData0510->Draw();
		
		DrawGammaLines(0., 30.,1., 1.,0.1);
		
		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Eta_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb0510.%s",outputDir.Data(),suffix.Data()));
	}

	
	//**********************************************************************************
	//**************************** Eta reweighting evalulation 20-40% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDEtaPbPbLHC11h2040 = NULL;
	if (directoryEtaPbPbLHC11h2040){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldEtaPtPbPbLHC11h2040->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldEtaPtPbPbLHC11h2040->Draw("hist,pe1,same");

		histoPCMEtaCorrectedSpecPbPbLHC11h2040->Draw("hist,pe1,same");
		
		fitYieldDataQCDEtaPbPbLHC11h2040 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h2040","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h2040,1.,10,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDEtaPbPbLHC11h2040, 1, 1.5, colorCombPbPb6080);
		fitYieldDataQCDEtaPbPbLHC11h2040->Draw("same");
		
		forOutput= WriteParameterToFile(fitYieldDataQCDEtaPbPbLHC11h2040);
		fileFinalResults<< forOutput.Data()<< endl;  
		
		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Eta_MCInputSpectraFittedPbPbLHC11h2040.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoEtaRatioDatatoFitQCDPbPb2040 = CalculateHistoRatioToFit (histoPCMEtaCorrectedSpecPbPbLHC11h2040, fitYieldDataQCDEtaPbPbLHC11h2040);
		TH1D* histoEtaRatioMCtoDataFitQCDPbPb2040 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11h2040, fitYieldDataQCDEtaPbPbLHC11h2040);
		TH1D* histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040 = NULL;
		if (histoMCYieldEtaPtPbPbLHC11hWOWeights2040) histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11hWOWeights2040, fitYieldDataQCDEtaPbPbLHC11h2040);
		
		canvasFraction2->cd();
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoEtaRatioDatatoFitQCDPbPb2040, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitQCDPbPb2040,
					"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
					kFALSE, 1.5, 0, kTRUE,
					kTRUE, -0.5, 8.,
					kTRUE, 0., 10);
		histoEtaRatioDatatoFitQCDPbPb2040->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoEtaRatioMCtoDataFitQCDPbPb2040, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
		if (runDrawReweighted) histoEtaRatioMCtoDataFitQCDPbPb2040->Draw("same,e,p");  
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040) histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040->Draw("same,e,p");  
		
		TLegend* legendEtaFit2040 = new TLegend(0.16,0.8,0.4,0.9);
		legendEtaFit2040->SetFillColor(0);
		legendEtaFit2040->SetLineColor(0);
		legendEtaFit2040->SetTextSize(0.035);
		legendEtaFit2040->SetMargin(0.2);
		legendEtaFit2040->AddEntry(histoEtaRatioDatatoFitQCDPbPb2040,"Data/Tsallis fit to Data (1 < pT < 10)","p");
		if (runDrawReweighted) legendEtaFit2040->AddEntry(histoEtaRatioMCtoDataFitQCDPbPb2040,"MC weighted/Tsallis fit to Data (1 < pT < 10)","p");
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040) legendEtaFit2040->AddEntry(histoEtaRatioMCUnweightedtoDataFitQCDPbPb2040,"MC/Tsallis fit to Data (1 < pT < 10)","p");
		legendEtaFit2040->Draw();
		TLatex *labelRatioMCData2040 = new TLatex(0.2,0.92,collisionSystemPbPb2040.Data());
		SetStyleTLatex( labelRatioMCData2040, 0.04,4);
		labelRatioMCData2040->Draw();

		DrawGammaLines(0., 30.,1., 1.,0.1);
		
		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Eta_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb2040.%s",outputDir.Data(),suffix.Data()));
	}

	
	//**********************************************************************************
	//**************************** Eta reweighting evalulation 20#font[122]{-}50% LHC11h*************
	//**********************************************************************************
	TF1* fitYieldDataQCDEtaPbPbLHC11h2050 = NULL;
	if (directoryEtaPbPbLHC11h2050){
		canvasSpectraMC->cd();
		histo2DSpectraMC->DrawCopy(); 

		histoMCYieldEtaPtPbPbLHC11h2050->SetMarkerStyle(markerStylePbPb6080MC);
		histoMCYieldEtaPtPbPbLHC11h2050->Draw("hist,pe1,same");
		
		histoMCYieldEtaPtPbPbLHC11hAddedSig2050->SetMarkerStyle(markerStylePbPb4060MC);
		histoMCYieldEtaPtPbPbLHC11hAddedSig2050->Draw("hist,pe1,same");

		histoPCMEtaCorrectedSpecPbPbLHC11h2050->Draw("hist,pe1,same");
		
		fitYieldDataQCDEtaPbPbLHC11h2050 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h2050","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h2050,1.,10,NULL,"QNRME+");
		DrawGammaSetMarkerTF1(fitYieldDataQCDEtaPbPbLHC11h2050, 1, 1.5, colorCombPbPb6080);
		fitYieldDataQCDEtaPbPbLHC11h2050->Draw("same");
		
		forOutput= WriteParameterToFile(fitYieldDataQCDEtaPbPbLHC11h2050);
		fileFinalResults<< forOutput.Data()<< endl;  

		TLegend* legendEtaYieldandFit2050 = new TLegend(0.66,0.81,0.9,0.9);
        legendEtaYieldandFit2050->SetFillColor(0);
        legendEtaYieldandFit2050->SetLineColor(0);
        legendEtaYieldandFit2050->SetTextSize(0.035);
        legendEtaYieldandFit2050->SetMargin(0.2);
        legendEtaYieldandFit2050->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h2050,"Data","p");
        legendEtaYieldandFit2050->AddEntry(histoMCYieldEtaPtPbPbLHC11hAddedSig2050,"MC yield","p");
        legendEtaYieldandFit2050->AddEntry(fitYieldDataQCDEtaPbPbLHC11h2050,"QCD fit to data","l");
        legendEtaYieldandFit2050->Draw();
		
		canvasSpectraMC->Update();
		canvasSpectraMC->Print(Form("%s/Eta_MCInputSpectraFittedPbPbLHC11h2050.%s",outputDir.Data(),suffix.Data()));

		TH1D* histoEtaRatioDatatoFitQCDPbPb2050 = CalculateHistoRatioToFit (histoPCMEtaCorrectedSpecPbPbLHC11h2050, fitYieldDataQCDEtaPbPbLHC11h2050);
		TH1D* histoEtaRatioMCtoDataFitQCDPbPb2050 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11h2050, fitYieldDataQCDEtaPbPbLHC11h2050);
		TH1D* histoEtaRatioMCAddedSigtoDataFitQCDPbPb2050 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11hAddedSig2050, fitYieldDataQCDEtaPbPbLHC11h2050);
		TH1D* histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050 = NULL;
		if (histoMCYieldEtaPtPbPbLHC11hWOWeights2050) histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050 = CalculateHistoRatioToFit (histoMCYieldEtaPtPbPbLHC11hWOWeights2050, fitYieldDataQCDEtaPbPbLHC11h2050);
		
		canvasFraction2->cd();
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050) DrawGammaSetMarker(histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510 );
		
		DrawGammaSetMarker(histoEtaRatioDatatoFitQCDPbPb2050, markerStylePbPb1020,markerSizePbPb1020, kBlack , kBlack);
		DrawAutoGammaMesonHistos( histoEtaRatioDatatoFitQCDPbPb2050,
					"", "#it{p}_{T} (GeV/#it{c})", "Spectrum/ fit to Spectrum",
					kFALSE, 1.5, 0, kTRUE,
					kTRUE, -0.5, 8.,
					kTRUE, 0., 10);
		histoEtaRatioDatatoFitQCDPbPb2050->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoEtaRatioMCtoDataFitQCDPbPb2050, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
		if (runDrawReweighted) histoEtaRatioMCtoDataFitQCDPbPb2050->Draw("same,e,p");  
		
		DrawGammaSetMarker(histoEtaRatioMCAddedSigtoDataFitQCDPbPb2050, markerStylePbPb4060MC,markerSizePbPb4060, colorCombPbPb4060, colorCombPbPb4060);
		if (runDrawReweighted) histoEtaRatioMCAddedSigtoDataFitQCDPbPb2050->Draw("same,e,p");  
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050) histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050->Draw("same,e,p");  

		TLegend* legendEtaFit2050 = new TLegend(0.16,0.73,0.4,0.9);
		legendEtaFit2050->SetFillColor(0);
		legendEtaFit2050->SetLineColor(0);
		legendEtaFit2050->SetTextSize(0.035);
		legendEtaFit2050->SetMargin(0.2);
		legendEtaFit2050->AddEntry(histoEtaRatioDatatoFitQCDPbPb2050,"Data/Tsallis fit to Data (1 < pT < 10)","p");
		if (runDrawReweighted) legendEtaFit2050->AddEntry(histoEtaRatioMCtoDataFitQCDPbPb2050,"MC weighted/Tsallis fit to Data (1 < pT < 10)","p");
		if (runDrawReweighted) legendEtaFit2050->AddEntry(histoEtaRatioMCAddedSigtoDataFitQCDPbPb2050,"MC added sig. weighted/Tsallis fit to Data (1 < pT < 10)","p");
		if (histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050) legendEtaFit2050->AddEntry(histoEtaRatioMCUnweightedtoDataFitQCDPbPb2050,"MC/Tsallis fit to Data (1 < pT < 10)","p");
		legendEtaFit2050->Draw();
		TLatex *labelRatioMCData2050 = new TLatex(0.2,0.92,collisionSystemPbPb2050.Data());
		SetStyleTLatex( labelRatioMCData2050, 0.04,4);
		labelRatioMCData2050->Draw();

		DrawGammaLines(0., 30.,1., 1.,0.1);

		canvasFraction2->Update();
		canvasFraction2->SaveAs(Form("%s/Eta_Ratio_MCLHC11h_To_DataFitLHC11h_PbPb2050.%s",outputDir.Data(),suffix.Data()));
	}
	
	
//    TFile fdatainput(nameFilePbPbLHC11h.Data(),"UPDATE");
//    	graphPi0RCP->Write("Pi0RCP",TObject::kOverwrite);
// 	graphPi0RCPSys->Write("Pi0RCPsys",TObject::kOverwrite);
// 	graphEtaRCP->Write("EtaRCP",TObject::kOverwrite);
// 	graphEtaRCPSys->Write("EtaRCPsys",TObject::kOverwrite);
//    fdatainput.Close();
   
//    TFile fEtatoPi0input("EtaToPi0InputsForCombination.root","UPDATE");
// 	
// 	graphPi0RAA0005->Write();
// 	graphPi0RAASys0005->Write();
// 	graphPi0RAA0010->Write();
// 	graphPi0RAASys0010->Write();
// 	graphPi0RAA2040->Write();
// 	graphPi0RAASys2040->Write();
// 	graphPi0RAA2050->Write();
// 	graphPi0RAASys2050->Write();
// 
// 	graphEtaRAA0005->Write();
// 	graphEtaRAASys0005->Write();
// 	graphEtaRAA0010->Write();
// 	graphEtaRAASys0010->Write();
// 	graphEtaRAA2040->Write();
// 	graphEtaRAASys2040->Write();
// 	graphEtaRAA2050->Write();
// 	graphEtaRAASys2050->Write();
// 
// 	graphRAAPublished0005->Write();
// 	graphRAASysPublished0005->Write();
// 	graphRAAPublished0010->Write();
// 	graphRAASysPublished0010->Write();
// 	graphRAAPublished2040->Write();
// 	graphRAASysPublished2040->Write();
// 
// 
// 	graphChargedKaonRAA0005->Write();
// 	graphChargedKaonRAASys0005->Write();
// 	graphChargedKaonRAA2040->Write();
// 	graphChargedKaonRAASys2040->Write();
// 
// 	graphChargedPionRAA0005->Write();
// 	graphChargedPionRAASys0005->Write();
// 	graphChargedPionRAA2040->Write();
// 	graphChargedPionRAASys2040->Write();
// 
//    fEtatoPi0input.Close();
	
	
   TFile fMCSpectraInput("MCSpectraInputLHC11h.root","UPDATE");
      if (fitYieldDataQCDPi0PbPbLHC11h0010){
         fitYieldDataQCDPi0PbPbLHC11h0010->SetRange(0,30);
		 fitYieldDataQCDPi0PbPbLHC11h0010->Write("Pi0_Fit_Data_PbPb_2760GeV_0010V0M",TObject::kOverwrite);
// 		 fitYieldDataQCDPi0PbPbLHC11h0010->Write("Pi0_Fit_Data_PbPbLHC11h_2760GeV_0010V0M",TObject::kOverwrite);
      }
      if (fitYieldDataQCDPi0PbPbLHC11h0005){
         fitYieldDataQCDPi0PbPbLHC11h0005->SetRange(0,30);
		 fitYieldDataQCDPi0PbPbLHC11h0005->Write("Pi0_Fit_Data_PbPb_2760GeV_0005V0M",TObject::kOverwrite);
//          fitYieldDataQCDPi0PbPbLHC11h0005->Write("Pi0_Fit_Data_PbPbLHC11h_2760GeV_0005V0M",TObject::kOverwrite);
      }
      if (fitYieldDataQCDPi0PbPbLHC11h0510){
         fitYieldDataQCDPi0PbPbLHC11h0510->SetRange(0,30);
		 fitYieldDataQCDPi0PbPbLHC11h0510->Write("Pi0_Fit_Data_PbPb_2760GeV_0510V0M",TObject::kOverwrite);
//          fitYieldDataQCDPi0PbPbLHC11h0510->Write("Pi0_Fit_Data_PbPbLHC11h_2760GeV_0510V0M",TObject::kOverwrite);
      }
      if (fitYieldDataQCDPi0PbPbLHC11h2040){
         fitYieldDataQCDPi0PbPbLHC11h2040->SetRange(0,30);
		 fitYieldDataQCDPi0PbPbLHC11h2040->Write("Pi0_Fit_Data_PbPb_2760GeV_2040V0M",TObject::kOverwrite);
//       fitYieldDataQCDPi0PbPbLHC11h2040->Write("Pi0_Fit_Data_PbPbLHC11h_2760GeV_2040V0M",TObject::kOverwrite);
      }
      if (fitYieldDataQCDPi0PbPbLHC11h2050){
         fitYieldDataQCDPi0PbPbLHC11h2050->SetRange(0,30);
		 fitYieldDataQCDPi0PbPbLHC11h2050->Write("Pi0_Fit_Data_PbPb_2760GeV_2050V0M",TObject::kOverwrite);
//       fitYieldDataQCDPi0PbPbLHC11h2050->Write("Pi0_Fit_Data_PbPbLHC11h_2760GeV_2050V0M",TObject::kOverwrite);
      }

      
      if (fitYieldDataQCDEtaPbPbLHC11h0010){
         fitYieldDataQCDEtaPbPbLHC11h0010->SetRange(0,30);
         fitYieldDataQCDEtaPbPbLHC11h0010->Write("Eta_Fit_Data_PbPb_2760GeV_0010V0M",TObject::kOverwrite);
      }
      if (fitYieldDataQCDEtaPbPbLHC11h0005){
         fitYieldDataQCDEtaPbPbLHC11h0005->SetRange(0,30);
		 fitYieldDataQCDEtaPbPbLHC11h0005->Write("Eta_Fit_Data_PbPb_2760GeV_0005V0M",TObject::kOverwrite);
//          fitYieldDataQCDEtaPbPbLHC11h0005->Write("Eta_Fit_Data_PbPbLHC11h_2760GeV_0005V0M",TObject::kOverwrite);
      }
      if (fitYieldDataQCDEtaPbPbLHC11h0510){
         fitYieldDataQCDEtaPbPbLHC11h0510->SetRange(0,30);
		 fitYieldDataQCDEtaPbPbLHC11h0510->Write("Eta_Fit_Data_PbPb_2760GeV_0510V0M",TObject::kOverwrite);
//          fitYieldDataQCDEtaPbPbLHC11h0510->Write("Eta_Fit_Data_PbPbLHC11h_2760GeV_0510V0M",TObject::kOverwrite);
      }
      if (fitYieldDataQCDEtaPbPbLHC11h2040){
         fitYieldDataQCDEtaPbPbLHC11h2040->SetRange(0,30);
		 fitYieldDataQCDEtaPbPbLHC11h2040->Write("Eta_Fit_Data_PbPb_2760GeV_2040V0M",TObject::kOverwrite);
//          fitYieldDataQCDEtaPbPbLHC11h2040->Write("Eta_Fit_Data_PbPbLHC11h_2760GeV_2040V0M",TObject::kOverwrite);
      }
	  if (fitYieldDataQCDEtaPbPbLHC11h2050){
         fitYieldDataQCDEtaPbPbLHC11h2050->SetRange(0,30);
         fitYieldDataQCDEtaPbPbLHC11h2050->Write("Eta_Fit_Data_PbPb_2760GeV_2050V0M",TObject::kOverwrite);
      }

	  cout << "here" << endl;
	  // only true if none of the MCs is mixed for final result
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0010->Write("Eta_Hijing_LHC14a1a_PbPb_2760GeV_0010TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0010->Write("Eta_Hijing_LHC14a1a_addSig_PbPb_2760GeV_0010TPC",TObject::kOverwrite);  
      histoMCYieldPi0PtPbPbLHC11hWOWeights0010->Write("Pi0_Hijing_LHC14a1a_PbPb_2760GeV_0010TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0010->Write("Pi0_Hijing_LHC14a1a_addSig_PbPb_2760GeV_0010TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0010->Write("Eta_Hijing_LHC14a1b_PbPb_2760GeV_0010TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0010->Write("Eta_Hijing_LHC14a1b_addSig_PbPb_2760GeV_0010TPC",TObject::kOverwrite);  
      histoMCYieldPi0PtPbPbLHC11hWOWeights0010->Write("Pi0_Hijing_LHC14a1b_PbPb_2760GeV_0010TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0010->Write("Pi0_Hijing_LHC14a1b_addSig_PbPb_2760GeV_0010TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0010->Write("Eta_Hijing_LHC14a1c_PbPb_2760GeV_0010TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0010->Write("Eta_Hijing_LHC14a1c_addSig_PbPb_2760GeV_0010TPC",TObject::kOverwrite);  
      histoMCYieldPi0PtPbPbLHC11hWOWeights0010->Write("Pi0_Hijing_LHC14a1c_PbPb_2760GeV_0010TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0010->Write("Pi0_Hijing_LHC14a1c_addSig_PbPb_2760GeV_0010TPC",TObject::kOverwrite);
	  cout << "here" << endl;
	  // only true if none of the MCs is mixed for final result
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0005->Write("Eta_Hijing_LHC14a1a_PbPb_2760GeV_0005TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0005->Write("Eta_Hijing_LHC14a1a_addSig_PbPb_2760GeV_0005TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights0005->Write("Pi0_Hijing_LHC14a1a_PbPb_2760GeV_0005TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0005->Write("Pi0_Hijing_LHC14a1a_addSig_PbPb_2760GeV_0005TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0005->Write("Eta_Hijing_LHC14a1b_PbPb_2760GeV_0005TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0005->Write("Eta_Hijing_LHC14a1b_addSig_PbPb_2760GeV_0005TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights0005->Write("Pi0_Hijing_LHC14a1b_PbPb_2760GeV_0005TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0005->Write("Pi0_Hijing_LHC14a1b_addSig_PbPb_2760GeV_0005TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0005->Write("Eta_Hijing_LHC14a1c_PbPb_2760GeV_0005TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0005->Write("Eta_Hijing_LHC14a1c_addSig_PbPb_2760GeV_0005TPC",TObject::kOverwrite);  
      histoMCYieldEtaPtPbPbLHC11hWOWeights0005->Write("Pi0_Hijing_LHC14a1c_PbPb_2760GeV_0005TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0005->Write("Pi0_Hijing_LHC14a1c_addSig_PbPb_2760GeV_0005TPC",TObject::kOverwrite);
	  cout << "here" << endl;
	  // only true if none of the MCs is mixed for final result
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0510->Write("Eta_Hijing_LHC14a1a_PbPb_2760GeV_0510TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0510->Write("Eta_Hijing_LHC14a1a_addSig_PbPb_2760GeV_0510TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights0510->Write("Pi0_Hijing_LHC14a1a_PbPb_2760GeV_0510TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0510->Write("Pi0_Hijing_LHC14a1a_addSig_PbPb_2760GeV_0510TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0510->Write("Eta_Hijing_LHC14a1b_PbPb_2760GeV_0510TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0510->Write("Eta_Hijing_LHC14a1b_addSig_PbPb_2760GeV_0510TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights0510->Write("Pi0_Hijing_LHC14a1b_PbPb_2760GeV_0510TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0510->Write("Pi0_Hijing_LHC14a1b_addSig_PbPb_2760GeV_0510TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights0510->Write("Eta_Hijing_LHC14a1c_PbPb_2760GeV_0510TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0510->Write("Eta_Hijing_LHC14a1c_addSig_PbPb_2760GeV_0510TPC",TObject::kOverwrite);  
      histoMCYieldEtaPtPbPbLHC11hWOWeights0510->Write("Pi0_Hijing_LHC14a1c_PbPb_2760GeV_0510TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0510->Write("Pi0_Hijing_LHC14a1c_addSig_PbPb_2760GeV_0510TPC",TObject::kOverwrite);
	  cout << "here" << endl;
	  // only true if none of the MCs is mixed for final result
	  histoMCYieldEtaPtPbPbLHC11hWOWeights2040->Write("Eta_Hijing_LHC14a1a_PbPb_2760GeV_2040TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2040->Write("Eta_Hijing_LHC14a1a_addSig_PbPb_2760GeV_2040TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights2040->Write("Pi0_Hijing_LHC14a1a_PbPb_2760GeV_2040TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2040->Write("Pi0_Hijing_LHC14a1a_addSig_PbPb_2760GeV_2040TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights2040->Write("Eta_Hijing_LHC14a1b_PbPb_2760GeV_2040TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2040->Write("Eta_Hijing_LHC14a1b_addSig_PbPb_2760GeV_2040TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights2040->Write("Pi0_Hijing_LHC14a1b_PbPb_2760GeV_2040TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2040->Write("Pi0_Hijing_LHC14a1b_addSig_PbPb_2760GeV_2040TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights2040->Write("Eta_Hijing_LHC14a1c_PbPb_2760GeV_2040TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2040->Write("Eta_Hijing_LHC14a1c_addSig_PbPb_2760GeV_2040TPC",TObject::kOverwrite);  
      histoMCYieldEtaPtPbPbLHC11hWOWeights2040->Write("Pi0_Hijing_LHC14a1c_PbPb_2760GeV_2040TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2040->Write("Pi0_Hijing_LHC14a1c_addSig_PbPb_2760GeV_2040TPC",TObject::kOverwrite);
	  cout << "here" << endl;
	  // only true if none of the MCs is mixed for final result
	  histoMCYieldEtaPtPbPbLHC11hWOWeights2050->Write("Eta_Hijing_LHC14a1a_PbPb_2760GeV_2050TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2050->Write("Eta_Hijing_LHC14a1a_addSig_PbPb_2760GeV_2050TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights2050->Write("Pi0_Hijing_LHC14a1a_PbPb_2760GeV_2050TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2050->Write("Pi0_Hijing_LHC14a1a_addSig_PbPb_2760GeV_2050TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights2050->Write("Eta_Hijing_LHC14a1b_PbPb_2760GeV_2050TPC",TObject::kOverwrite);  
	  histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2050->Write("Eta_Hijing_LHC14a1b_addSig_PbPb_2760GeV_2050TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hWOWeights2050->Write("Pi0_Hijing_LHC14a1b_PbPb_2760GeV_2050TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2050->Write("Pi0_Hijing_LHC14a1b_addSig_PbPb_2760GeV_2050TPC",TObject::kOverwrite);
	  histoMCYieldEtaPtPbPbLHC11hWOWeights2050->Write("Eta_Hijing_LHC14a1c_PbPb_2760GeV_2050TPC",TObject::kOverwrite);  
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2050->Write("Eta_Hijing_LHC14a1c_addSig_PbPb_2760GeV_2050TPC",TObject::kOverwrite);  
      histoMCYieldEtaPtPbPbLHC11hWOWeights2050->Write("Pi0_Hijing_LHC14a1c_PbPb_2760GeV_2050TPC",TObject::kOverwrite);
	  histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2050->Write("Pi0_Hijing_LHC14a1c_addSig_PbPb_2760GeV_2050TPC",TObject::kOverwrite);
	  
   fMCSpectraInput.Close();
}

