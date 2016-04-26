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
#include "CombineMesonMeasurementsPbPbLHC11hV2.h"


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


void CombineNeutralPionResultsLHC11h(TString suffix = "pdf", TString nameFilePP = "data_PCMResultsFullCorrection_PP_NoBinShifting.root", TString nameFilePbPbLHC11h = "data_PCMResults_PbPb_2.76TeV", Bool_t runDrawReweighted = kTRUE){

	gROOT->Reset();   
	gROOT->SetStyle("Plain");
	
	TString nameFileEtaPbPbLHC10h = "data_PCMResults_PbPb_2.76TeV_LHC10hEta.root";
	TString nameFilePi0PbPbLHC10h = "data_PCMResults_PbPb_2.76TeV_LHC10h.root";
	
	TString dateForOutput = ReturnDateStringForOutput();
	TString outputDir = Form("%s/%s/CombineNeutralPionResults",suffix.Data(),dateForOutput.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	gSystem->Exec(Form("cp %s %s/InputFilePCMPionPP.root ",nameFilePP.Data(),outputDir.Data() ));
	gSystem->Exec(Form("cp %s %s/InputFilePCMPionPbPbLHC10h.root ",nameFilePi0PbPbLHC10h.Data(),outputDir.Data() ));
	gSystem->Exec(Form("cp %s %s/InputFilePCMEtaPbPbLHC10h.root ",nameFileEtaPbPbLHC10h.Data(),outputDir.Data() ));
	gSystem->Exec(Form("cp %s %s/InputFilePCMPbPbLHC11h.root ",nameFilePbPbLHC11h.Data(),outputDir.Data() ));
	
	gSystem->Exec("mkdir -p "+outputDir);
   
	StyleSettingsThesis();  
	SetPlotStyle();

	Color_t  colorCombPP             = kBlack;
	Color_t  colorCombPbPb0005          = kRed+1;
	Color_t  colorCombPbPb0010          = kRed+1;
	Color_t  colorCombPbPb0510          = 807;
	Color_t  colorCombPbPb1020          = 800;
	Color_t  colorCombPbPb2040          = kAzure-5;//kGreen+2;
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
	Style_t  markerStylePbPb2040  = 33 ;
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
	Size_t   markerSizePbPb2040   = 2.5;
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
	
	TString collisionSystemPP2760GeV = "pp #sqrt{#it{s}} = 2.76 TeV";      
	TString collisionSystemPP7TeV = "pp #sqrt{#it{s}} = 7 TeV";      
	TString collisionSystemPP900GeV = "pp #sqrt{#it{s}} = 0.9 TeV";     
	TString collisionSystemPbPb0005 = "0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";      
	TString collisionSystemPbPb0010 = "0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";      
	TString collisionSystemPbPb0510 = "5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";     
	TString collisionSystemPbPb1020 = "10-20% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb2030 = "20-30% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb3040 = "30-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb2040 = "20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb2050 = "20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb3050 = "30-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb4060 = "40-60% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb6080 = "60-80% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb6070 = "60-70% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb7080 = "70-80% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb7590 = "75-90% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb8090 = "80-90% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb4050 = "40-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb5060 = "50-60% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";    
	TString collisionSystemPbPb0020 = "0-20% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";     
	TString collisionSystemPbPb0080 = "0-80% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";     
	TString collisionSystemPbPb0040 = "0-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";     
	TString collisionSystemPbPb4080 = "40-80% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV";     
	
	Size_t markerSizeComparison = 0.5;
	Double_t maxPtMesonEffFit = 12.;
	Double_t minPtMesonEffFit = 1.2;
	Int_t offsetCorrectionHighPt= 1;
	TF1* fitTrueEffi = new TF1("EffiFitDummy","1 - [0]*exp([1]*x)+[1]");
	fitTrueEffi->SetRange(minPtMesonEffFit,maxPtMesonEffFit);
   
	TLatex *labelRawPi0PbPb = new TLatex(0.6,0.9,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelRawPi0PbPb, 0.038,4);

	TLatex *labelRawEtaPbPb = new TLatex(0.6,0.9,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelRawEtaPbPb, 0.038,4);

	TLatex *labelMassEtaPP = new TLatex(0.2,0.9,collisionSystemPP7TeV.Data());
	SetStyleTLatex( labelMassEtaPP, 0.062,4);

	TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
	SetStyleTLatex( labelLegendAMass, 0.08,4);


/////////////////////////////////////////////////////////////// External file ////////////////////////////////////////////////////////////////////////////


	TFile *filePublished = new TFile("FinalResults/CombinedResultsPbPb_ShiftedX_PaperRAA_13_Aug_2014.root");
    TGraphAsymmErrors* graphYieldsPublished0010 = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMStatErr_0010"); 
    TGraphAsymmErrors* graphYieldsSysPublished0010 = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMSysErr_0010"); 
    
    TGraphAsymmErrors* graphYieldsPublished2040 = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMStatErr_2040"); 
    TGraphAsymmErrors* graphYieldsSysPublished2040 = (TGraphAsymmErrors*)filePublished->Get("InvYieldPbPbPCMSysErr_2040"); 
	
    TGraphAsymmErrors* graphRAAPublished0005 = (TGraphAsymmErrors*)filePublished->Get("graphRAAStatErr_0005"); 
	TGraphAsymmErrors* graphRAASysPublished0005 = (TGraphAsymmErrors*)filePublished->Get("graphRAASysErr_0005"); 
	

    TGraphAsymmErrors* graphRAAPublished0010 = (TGraphAsymmErrors*)filePublished->Get("graphRAAStatErr_0010"); 
	TGraphAsymmErrors* graphRAASysPublished0010 = (TGraphAsymmErrors*)filePublished->Get("graphRAASysErr_0010"); 
	
	TGraphAsymmErrors* graphRAAPublished2040 = (TGraphAsymmErrors*)filePublished->Get("graphRAAStatErr_2040"); 
	TGraphAsymmErrors* graphRAASysPublished2040 = (TGraphAsymmErrors*)filePublished->Get("graphRAASysErr_2040"); 

	
////////////////////////////////////////////////////////////////// LHC10h file //////////////////////////////////////////////////////////////////////////////   
    TFile*   filePCMPbPbLHC10h = new TFile("data_PCMResults_PbPb_2.76TeV_LHC10h.root"); //data_PCMResults_PbPb_2.76TeV_LHC10h.root
     
//////////////////////// Pi0 ////////////////////////////   
    cout << "Pi0 0-10%" << endl;
    TDirectory* directoryPi0PbPbLHC10h0010 = (TDirectory*)filePCMPbPbLHC10h->Get("Pi0_PbPb_2.76TeV_0-10%"); 

        TH1D* histoPCMPi0CorrectedSpecPbPbLHC10h0010                 = (TH1D*)directoryPi0PbPbLHC10h0010->Get("CorrectedYieldPi0");   
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC10h0010 = (TGraphAsymmErrors*)directoryPi0PbPbLHC10h0010->Get("Pi0SystErrorA"); 

    cout << "20-40%" << endl;
    TDirectory* directoryPi0PbPbLHC10h2040 = (TDirectory*)filePCMPbPbLHC10h->Get("Pi0_PbPb_2.76TeV_20-40%"); 
    
        TH1D* histoPCMPi0CorrectedSpecPbPbLHC10h2040                 = (TH1D*)directoryPi0PbPbLHC10h2040->Get("CorrectedYieldPi0");   
        TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC10h2040 = (TGraphAsymmErrors*)directoryPi0PbPbLHC10h2040->Get("Pi0SystErrorA"); 
        
	

		
/////////////////////////////////////////////// Charged pions/kaons file /////////////////////////////////////////////////////////
		
	TFile *fileChargedKaonRAA = new TFile("ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Kaon_08052014.root");
	TH1D *histoStatChargedKaon0005 = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Stat_0_5");
	TH1D *histoSystChargedKaon0005 = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Syst_0_5");
// 	histoStatChargedKaon0005->Scale(0.5);
// 	histoSystChargedKaon0005->Scale(0.5);

	TGraphAsymmErrors* graphChargedKaonRAA0005 = new TGraphAsymmErrors(histoStatChargedKaon0005); 
	TGraphAsymmErrors* graphChargedKaonRAASys0005 = new TGraphAsymmErrors(histoSystChargedKaon0005); 
	
	TH1D *histoStatChargedKaon0510 = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Stat_5_10");
	TH1D *histoSystChargedKaon0510 = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Syst_5_10");
// 	histoStatChargedKaon0510->Scale(0.5);
// 	histoSystChargedKaon0510->Scale(0.5);
	
	TGraphAsymmErrors* graphChargedKaonRAA0510 = new TGraphAsymmErrors(histoStatChargedKaon0510); 
	TGraphAsymmErrors* graphChargedKaonRAASys0510 = new TGraphAsymmErrors(histoSystChargedKaon0510); 


	TH1D* histoStatChargedKaon0010 = (TH1D*)histoStatChargedKaon0510->Clone("histoStatChargedKaon0010");
	TH1D* histoSystChargedKaon0010 = (TH1D*)histoSystChargedKaon0510->Clone("histoSystChargedKaon0010");
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

	TGraphAsymmErrors* graphChargedKaonRAA0010 = new TGraphAsymmErrors(histoStatChargedKaon0010); 
	TGraphAsymmErrors* graphChargedKaonRAASys0010 = new TGraphAsymmErrors(histoSystChargedKaon0010); 

	TH1D *histoStatChargedKaon2040 = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Stat_20_40");
	TH1D *histoSystChargedKaon2040 = (TH1D*)fileChargedKaonRAA->Get("RAAKaon_Syst_20_40");
	TGraphAsymmErrors* graphChargedKaonRAA2040 = new TGraphAsymmErrors(histoStatChargedKaon2040); 
	TGraphAsymmErrors* graphChargedKaonRAASys2040 = new TGraphAsymmErrors(histoSystChargedKaon2040); 


	
	TFile *fileChargedPionRAA = new TFile("ExternalInputPbPb/IdentifiedCharged/JIRA_PWGLF-258/RAA_Pion_08052014.root");
	TH1D *histoStatChargedPion0005 = (TH1D*)fileChargedPionRAA->Get("RAAPion_Stat_0_5");
	TH1D *histoSystChargedPion0005 = (TH1D*)fileChargedPionRAA->Get("RAAPion_Syst_0_5");
// 	histoStatChargedPion0005->Scale(0.5);
// 	histoSystChargedPion0005->Scale(0.5);

	TGraphAsymmErrors* graphChargedPionRAA0005 = new TGraphAsymmErrors(histoStatChargedPion0005); 
	TGraphAsymmErrors* graphChargedPionRAASys0005 = new TGraphAsymmErrors(histoSystChargedPion0005); 
	
	TH1D *histoStatChargedPion0510 = (TH1D*)fileChargedPionRAA->Get("RAAPion_Stat_5_10");
	TH1D *histoSystChargedPion0510 = (TH1D*)fileChargedPionRAA->Get("RAAPion_Syst_5_10");
// 	histoStatChargedPion0510->Scale(0.5);
// 	histoSystChargedPion0510->Scale(0.5);
	
	TGraphAsymmErrors* graphChargedPionRAA0510 = new TGraphAsymmErrors(histoStatChargedPion0510); 
	TGraphAsymmErrors* graphChargedPionRAASys0510 = new TGraphAsymmErrors(histoSystChargedPion0510); 


	TH1D* histoStatChargedPion0010 = (TH1D*)histoStatChargedPion0510->Clone("histoStatChargedPion0010");
	TH1D* histoSystChargedPion0010 = (TH1D*)histoSystChargedPion0510->Clone("histoSystChargedPion0010");
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

	TGraphAsymmErrors* graphChargedPionRAA0010 = new TGraphAsymmErrors(histoStatChargedPion0010); 
	TGraphAsymmErrors* graphChargedPionRAASys0010 = new TGraphAsymmErrors(histoSystChargedPion0010); 
	
	TH1D *histoStatChargedPion2040 = (TH1D*)fileChargedPionRAA->Get("RAAPion_Stat_20_40");
	TH1D *histoSystChargedPion2040 = (TH1D*)fileChargedPionRAA->Get("RAAPion_Syst_20_40");
	TGraphAsymmErrors* graphChargedPionRAA2040 = new TGraphAsymmErrors(histoStatChargedPion2040); 
	TGraphAsymmErrors* graphChargedPionRAASys2040 = new TGraphAsymmErrors(histoSystChargedPion2040); 
	

////////////////////////////////////////////////// PHENIX adn SPS file //////////////////////////////////////////////////////////////////
	
	TFile *filePHENIX = new TFile("ExternalInputPbPb/OtherExperiments/DataCompilationFromOtherEnergiesPbPbWithEta.root");
	TGraphErrors* 	graphWA98_17_3GeVRAA_0013= (TGraphErrors*)filePHENIX->Get("graphWA98RAA_0013");
	TGraphErrors* graphPHENIX200GeVRAA_0010= (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVRAA_0010");
	TGraphErrors* graphPHENIX200GeVRAA_2040= (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVRAA_2040");
	TGraphErrors* graphPHENIX39GeVRAA_0010= (TGraphErrors*)filePHENIX->Get("graphPHENIX39GeVRAA_0010");
	TGraphErrors* graphPHENIX39GeVRAA_2040= (TGraphErrors*)filePHENIX->Get("graphPHENIX39GeVRAA_2040");
	TGraphErrors* graphPHENIX62GeVRAA_0010= (TGraphErrors*)filePHENIX->Get("graphPHENIX62GeVRAA_0010");
	TGraphErrors* graphPHENIX62GeVRAA_2040= (TGraphErrors*)filePHENIX->Get("graphPHENIX62GeVRAA_2040");
	TGraphErrors* graphPHENIX200GeVRAAvsNPartAt7GeVc = (TGraphErrors*)filePHENIX->Get("graphPHENIX200GeVRAAvsNPartAt7GeVc");
	TGraphErrors* graphPHENIX62GeVRAAvsNPartAt7GeVc = (TGraphErrors*)filePHENIX->Get("graphPHENIX62GeVRAAvsNPartAt7GeVc");
	TGraphErrors* graphPHENIX39GeVRAAvsNPartAt7GeVc = (TGraphErrors*)filePHENIX->Get("graphPHENIX39GeVRAAvsNPartAt7GeVc");

	TGraphAsymmErrors* graphEtaRAAPHENIX0010 = (TGraphAsymmErrors*)filePHENIX->Get("graphPHENIX200GeVEtaRAA_0010"); 
// 	TGraphAsymmErrors* graphPi0RAAPHENIX0010 = (TGraphAsymmErrors*)filePHENIX->Get("graphPHENIX200GeVRAA_0010"); 
	TGraphAsymmErrors* graphEtaRAAPHENIX2040 = (TGraphAsymmErrors*)filePHENIX->Get("graphPHENIX200GeVEtaRAA_2040"); 
// 	TGraphAsymmErrors* graphPi0RAAPHENIX2040 = (TGraphAsymmErrors*)filePHENIX->Get("graphPHENIX200GeVRAA_2040"); 
	TGraphAsymmErrors* graphEtaRAAPHENIX2060 = (TGraphAsymmErrors*)filePHENIX->Get("graphPHENIX200GeVEtaRAA_2060"); 

	Double_t normErr0010 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr0010/tAA0010),2)+pow((commonCentralityErr0010/100),2),0.5);
	Double_t normErr2040 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
	Double_t normErr2050 = pow(pow(xSection2760GeVErrpp/(xSection2760GeVpp*1e3),2)+pow((tAAErr2040/tAA2040),2)+pow((commonCentralityErr2040/100),2),0.5);
	
	
////////////////////////////////////////////////////////////////// PP file //////////////////////////////////////////////////////////////////////////////
	TFile* fileNeutralPionDataPP = new TFile(nameFilePP.Data());
	
	TDirectory* directoryPi07TeV =      (TDirectory*)fileNeutralPionDataPP->Get("Pi07TeV"); 
	TH1D* histoAccPi07TeV =             (TH1D*)directoryPi07TeV->Get("AcceptancePi0");
	TH1D* histoTrueEffPtPi07TeV =       (TH1D*)directoryPi07TeV->Get("EfficiencyPi0");
	TH1D* histoRawYieldPi07TeV =        (TH1D*)directoryPi07TeV->Get("RAWYieldPerEventsPi0");
	TDirectory* directoryEta7TeV =      (TDirectory*)fileNeutralPionDataPP->Get("Eta7TeV");
	TH1D* histoCorrectedYieldEta7TeV = (TH1D*)directoryEta7TeV->Get("CorrectedYieldEta");
	TGraphAsymmErrors* graphCorrectedYieldSysEta7TeV =   (TGraphAsymmErrors*)directoryEta7TeV->Get("EtaSystError");
	TH1D* histoAccEta7TeV =             (TH1D*)directoryEta7TeV->Get("AcceptanceEta");
	TH1D* histoTrueEffPtEta7TeV =       (TH1D*)directoryEta7TeV->Get("EfficiencyEta");
	TH1D* histoRawYieldEta7TeV =        (TH1D*)directoryEta7TeV->Get("RAWYieldPerEventsEta");
	TH1D* histoEtaMassData7TeV =        (TH1D*)directoryEta7TeV->Get("MassEta");
	TH1D* histoEtaWidthData7TeV =       (TH1D*)directoryEta7TeV->Get("FWHMEtaMeV");
	TH1D* histoEtaMassMC7TeV =          (TH1D*)directoryEta7TeV->Get("TrueMassEta");
	TH1D* histoEtaWidthMC7TeV =         (TH1D*)directoryEta7TeV->Get("TrueFWHMEtaMeV");
	TH1D* histoRatioEtaPi07TeV=        (TH1D*)directoryEta7TeV->Get("EtatoPi0RatioConversion");   
	TGraphAsymmErrors* graphRatioEtaPi0SystErr7TeV=             (TGraphAsymmErrors*)directoryEta7TeV->Get("EtatoPi0RatioConversionSys");

	histoEtaMassData7TeV->Scale(1000.);
	histoEtaMassMC7TeV->Scale(1000.);
	histoEtaMassData7TeV->SetBinContent(histoEtaMassData7TeV->GetNbinsX(),0);
	histoEtaMassMC7TeV->SetBinContent(histoEtaMassMC7TeV->GetNbinsX(),0);
	histoEtaWidthData7TeV->SetBinContent(histoEtaWidthData7TeV->GetNbinsX(),10000.);
	histoEtaWidthMC7TeV->SetBinContent(histoEtaWidthMC7TeV->GetNbinsX(),10000.);

	TDirectory* directoryPi0900GeV =    (TDirectory*)fileNeutralPionDataPP->Get("Pi0900GeV"); 
	TH1D* histoAccPi0900GeV =           (TH1D*)directoryPi0900GeV->Get("AcceptancePi0");
	TH1D* histoTrueEffPtPi0900GeV =  (TH1D*)directoryPi0900GeV->Get("EfficiencyPi0");
	TH1D* histoRawYieldPi0900GeV =      (TH1D*)directoryPi0900GeV->Get("RAWYieldPerEventsPi0");
	TDirectory* directoryEta900GeV =    (TDirectory*)fileNeutralPionDataPP->Get("Eta900GeV"); 
	TH1D* histoCorrectedYieldEta900GeV =     (TH1D*)directoryEta900GeV->Get("CorrectedYieldEta");
	TGraphAsymmErrors* graphCorrectedYieldSysEta900GeV =  (TGraphAsymmErrors*)directoryEta900GeV->Get("EtaSystError");
	TH1D* histoAccEta900GeV =           (TH1D*)directoryEta900GeV->Get("AcceptanceEta");
	TH1D* histoTrueEffPtEta900GeV =  (TH1D*)directoryEta900GeV->Get("EfficiencyEta");
	TH1D* histoRawYieldEta900GeV =      (TH1D*)directoryEta900GeV->Get("RAWYieldPerEventsEta");
	TH1D* histoRatioEtaPi0900GeV=      (TH1D*)directoryEta900GeV->Get("EtatoPi0RatioConversion");   
	TGraphAsymmErrors* graphRatioEtaPi0SystErr900GeV=           (TGraphAsymmErrors*)directoryEta900GeV->Get("EtatoPi0RatioConversionSys");
	
	TDirectory* directoryPi02760GeV =   (TDirectory*)fileNeutralPionDataPP->Get("Pi02.76TeV"); 
	TH1D* histoAccPi02760GeV =          (TH1D*)directoryPi02760GeV->Get("AcceptancePi0");
	TH1D* histoTrueEffPtPi02760GeV =    (TH1D*)directoryPi02760GeV->Get("EfficiencyPi0");
	TH1D* histoRawYieldPi02760GeV =  (TH1D*)directoryPi02760GeV->Get("RAWYieldPerEventsPi0");
	TH1D* histoCorrectedYieldPi02760GeV =  (TH1D*)directoryPi02760GeV->Get("CorrectedYieldPi0");
	TGraphAsymmErrors* graphCorrectedYieldSysPi02760GeV = (TGraphAsymmErrors*)directoryPi02760GeV->Get("Pi0SystError"); 
	TDirectory* directoryEta2760GeV =   (TDirectory*)fileNeutralPionDataPP->Get("Eta2.76TeV"); 
	TH1D* histoAccEta2760GeV =          (TH1D*)directoryEta2760GeV->Get("AcceptanceEta");
	TH1D* histoTrueEffPtEta2760GeV =    (TH1D*)directoryEta2760GeV->Get("EfficiencyEta");
	TH1D* histoCorrectedYieldEta2760GeV =  (TH1D*)directoryEta2760GeV->Get("CorrectedYieldEta");
	TGraphAsymmErrors* graphCorrectedYieldSysEta2760GeV = (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtaSystError");
	TH1D* histoRawYieldEta2760GeV =  (TH1D*)directoryEta2760GeV->Get("RAWYieldPerEventsEta");
	TH1D* histoEtaMassData2760GeV =     (TH1D*)directoryEta2760GeV->Get("MassEta");
	TH1D* histoEtaWidthData2760GeV =    (TH1D*)directoryEta2760GeV->Get("FWHMEtaMeV");
	TH1D* histoEtaMassMC2760GeV =       (TH1D*)directoryEta2760GeV->Get("TrueMassEta");
	TH1D* histoEtaWidthMC2760GeV =      (TH1D*)directoryEta2760GeV->Get("TrueFWHMEtaMeV");
	TH1D* histoRatioEtaPi02760GeV=     (TH1D*)directoryEta2760GeV->Get("EtatoPi0RatioConversion");   
	TGraphAsymmErrors* graphRatioEtaPi0SystErr2760GeV=          (TGraphAsymmErrors*)directoryEta2760GeV->Get("EtatoPi0RatioConversionSys");

	histoEtaMassData2760GeV->Scale(1000.);
	histoEtaMassMC2760GeV->Scale(1000.);
	histoEtaMassData2760GeV->SetBinContent(histoEtaMassData2760GeV->GetNbinsX(),0);
	histoEtaMassMC2760GeV->SetBinContent(histoEtaMassMC2760GeV->GetNbinsX(),0);
	histoEtaWidthData2760GeV->SetBinContent(histoEtaWidthData2760GeV->GetNbinsX(),10000.);
	histoEtaWidthMC2760GeV->SetBinContent(histoEtaWidthMC2760GeV->GetNbinsX(),10000.);
	
	TH1D* histoPCMMassDataPP2760GeV =      (TH1D*)directoryPi02760GeV->Get("MassPi0");
	TH1D* histoPCMWidthDataPP2760GeV =     (TH1D*)directoryPi02760GeV->Get("FWHMPi0MeV");
	TH1D* histoPCMMassMCPP2760GeV =        (TH1D*)directoryPi02760GeV->Get("TrueMassPi0");
	TH1D* histoPCMWidthMCPP2760GeV =          (TH1D*)directoryPi02760GeV->Get("TrueFWHMPi0MeV");
	histoPCMMassDataPP2760GeV->Scale(1000.);
	histoPCMMassMCPP2760GeV->Scale(1000.);
	histoPCMMassDataPP2760GeV->SetBinContent(histoPCMMassDataPP2760GeV->GetNbinsX(),0);
	histoPCMMassMCPP2760GeV->SetBinContent(histoPCMMassMCPP2760GeV->GetNbinsX(),0);
	histoPCMWidthDataPP2760GeV->SetBinContent(histoPCMWidthDataPP2760GeV->GetNbinsX(),10000.);
	histoPCMWidthMCPP2760GeV->SetBinContent(histoPCMWidthMCPP2760GeV->GetNbinsX(),10000.);

	
	
	
	TFile*   fileCocktail =                new TFile("CocktailInput/cocktail_allCentpluspp.root");
	TDirectory* directoryCocktailpp2760GeV =           (TDirectory*)fileCocktail->Get("cocktail_pp_2760GeV_qcd"); 
	TH1D* histoEtaFromCocktailpp2760GeV = (TH1D*)directoryCocktailpp2760GeV->Get("ptEta"); 
	TDirectory* directoryCocktail0010 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_0010_qcd"); 
	TH1D* histoEtaFromCocktail0010 = (TH1D*)directoryCocktail0010->Get("ptEta");  
	TDirectory* directoryCocktail0510 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_0510_qcd"); 
	TH1D* histoEtaFromCocktail0510 = (TH1D*)directoryCocktail0510->Get("ptEta");  
	TDirectory* directoryCocktail1020 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_1020_qcd"); 
	TH1D* histoEtaFromCocktail1020 = (TH1D*)directoryCocktail1020->Get("ptEta");  
	TDirectory* directoryCocktail2040 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_2040_qcd"); 
	TH1D* histoEtaFromCocktail2040 = (TH1D*)directoryCocktail2040->Get("ptEta");  
	TDirectory* directoryCocktail4060 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_4060_qcd"); 
	TH1D* histoEtaFromCocktail4060 = (TH1D*)directoryCocktail4060->Get("ptEta");  
	TDirectory* directoryCocktail6080 =             (TDirectory*)fileCocktail->Get("cocktail_PbPb_6080_qcd"); 
	TH1D* histoEtaFromCocktail6080 = (TH1D*)directoryCocktail6080->Get("ptEta");  
	cout << "here 7TeV" << endl;
	
	TFile* fileNeutralPionCombDataPP = new TFile("CombinedResultsPaperX_18_Feb_2014.root");
	TGraphAsymmErrors* graphInvYieldPi0Comb7TeV= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeV");
	graphInvYieldPi0Comb7TeV = ScaleGraph(graphInvYieldPi0Comb7TeV,1./xSection7TeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0Comb7TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVStatErr");
	graphInvYieldPi0Comb7TeVStatErr = ScaleGraph(graphInvYieldPi0Comb7TeVStatErr,1./xSection7TeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0Comb7TeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb7TeVSysErr");
	graphInvYieldPi0Comb7TeVSysErr = ScaleGraph(graphInvYieldPi0Comb7TeVSysErr,1./xSection7TeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PCM7TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMStat7TeV");
	graphInvYieldPi0PCM7TeVStatErr = ScaleGraph(graphInvYieldPi0PCM7TeVStatErr,1./xSection7TeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PCM7TeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMSys7TeV");
	graphInvYieldPi0PCM7TeVSysErr = ScaleGraph(graphInvYieldPi0PCM7TeVSysErr,1./xSection7TeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PHOS7TeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSStat7TeV");
	graphInvYieldPi0PHOS7TeVStatErr = ScaleGraph(graphInvYieldPi0PHOS7TeVStatErr,1./xSection7TeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PHOS7TeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSSys7TeV");
	graphInvYieldPi0PHOS7TeVSysErr = ScaleGraph(graphInvYieldPi0PHOS7TeVSysErr,1./xSection7TeVppINEL);
	cout << "here 2.76TeV" << endl;
	TGraphAsymmErrors* graphInvYieldPi0Comb2760GeV= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeV");
	graphInvYieldPi0Comb2760GeV = ScaleGraph(graphInvYieldPi0Comb2760GeV,1./xSection2760GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVStatErr");
	graphInvYieldPi0Comb2760GeVStatErr = ScaleGraph(graphInvYieldPi0Comb2760GeVStatErr,1./xSection2760GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0Comb2760GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb2760GeVSysErr");
	graphInvYieldPi0Comb2760GeVSysErr = ScaleGraph(graphInvYieldPi0Comb2760GeVSysErr,1./xSection2760GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PCM2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCM2760GeVStatErr");
	graphInvYieldPi0PCM2760GeVStatErr = ScaleGraph(graphInvYieldPi0PCM2760GeVStatErr,1./xSection2760GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PCM2760GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCM2760GeVSysErr");
	graphInvYieldPi0PCM2760GeVSysErr = ScaleGraph(graphInvYieldPi0PCM2760GeVSysErr,1./xSection2760GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PHOS2760GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOS2760GeVStatErr");
	graphInvYieldPi0PHOS2760GeVStatErr = ScaleGraph(graphInvYieldPi0PHOS2760GeVStatErr,1./xSection2760GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PHOS2760GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOS2760GeVSysErr");
	graphInvYieldPi0PHOS2760GeVSysErr = ScaleGraph(graphInvYieldPi0PHOS2760GeVSysErr,1./xSection2760GeVppINEL);
	
	cout << "here 900 GeV" << endl;
	TGraphAsymmErrors* graphInvYieldPi0Comb900GeV= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeV");
	graphInvYieldPi0Comb900GeV = ScaleGraph(graphInvYieldPi0Comb900GeV,1./xSection900GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0Comb900GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVStatErr");
	graphInvYieldPi0Comb900GeVStatErr = ScaleGraph(graphInvYieldPi0Comb900GeVStatErr,1./xSection900GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0Comb900GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0Comb900GeVSysErr");
	graphInvYieldPi0Comb900GeVSysErr = ScaleGraph(graphInvYieldPi0Comb900GeVSysErr,1./xSection900GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PCM900GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMStat900GeV");
	graphInvYieldPi0PCM900GeVStatErr = ScaleGraph(graphInvYieldPi0PCM900GeVStatErr,1./xSection900GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PCM900GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PCMSys900GeV");
	graphInvYieldPi0PCM900GeVSysErr = ScaleGraph(graphInvYieldPi0PCM900GeVSysErr,1./xSection900GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PHOS900GeVStatErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSStat900GeV");
	graphInvYieldPi0PHOS900GeVStatErr = ScaleGraph(graphInvYieldPi0PHOS900GeVStatErr,1./xSection900GeVppINEL);
	TGraphAsymmErrors* graphInvYieldPi0PHOS900GeVSysErr= (TGraphAsymmErrors*)fileNeutralPionCombDataPP->Get("graphInvCrossSectionPi0PHOSSys900GeV");
	graphInvYieldPi0PHOS900GeVSysErr = ScaleGraph(graphInvYieldPi0PHOS900GeVSysErr,1./xSection900GeVppINEL);

	cout << "efficiencies Added Sig 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPPEffiAddedSig = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosAddSig_000001100209366300380000000_01631031009000.root");
	TH1D* histoEffiAddedSigPP = (TH1D*)fileNeutralPion2760GeVDataPPEffiAddedSig->Get("TrueMesonEffiPt");
	
	TFile* fileNeutralPion7TeVDataPPEffi2760GeVCut = new TFile("ppAdditionalInput/0000011002093663003800000_01631031009/7TeV/Pi0_MC_GammaConvV1CorrectionHistosD_0000011002093663003800000_01631031009.root");
	TH1D* histoEffi7TeV2760GeVCut = (TH1D*)fileNeutralPion7TeVDataPPEffi2760GeVCut->Get("TrueMesonEffiPt");
	
	TFile* fileNeutralPion2760GeVDataPPEffiWithSDD = new TFile("ppAdditionalInput/0002011002093663003800000_01631031009/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistos_0002011002093663003800000_01631031009.root");
	TH1D* histoEffi2760GeVWithSDD = (TH1D*)fileNeutralPion2760GeVDataPPEffiWithSDD->Get("TrueMesonEffiPt");
	
	cout << "efficiencies MinBias 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPPEffiMinBias = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosMinBias_000001100209366300380000000_01631031009000.root");
	TH1D* histoEffiMinBiasPP = (TH1D*)fileNeutralPion2760GeVDataPPEffiMinBias->Get("TrueMesonEffiPt");

//    cout << "old spectra 2.76TeV GeV" << endl;
//    TFile* fileNeutralPionDataPPOld = new TFile("0000011002093663003800000_01631031009_old/data_GammaConversionResultsFullCorrectionNoBinShifting.root");
//    TDirectory* directoryPi02760GeVOld =         (TDirectory*)fileNeutralPionDataPPOld->Get("Pi02.76TeV"); 
//    TH1D* histoRawYieldPi02760GeVOld =           (TH1D*)directoryPi02760GeVOld->Get("RAWYieldPerEventsPi0");
//    TH1D* histoCorrectedYieldPi02760GeVOld =           (TH1D*)directoryPi02760GeVOld->Get("CorrectedYieldPi0");
   
	cout << "MC spectra 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPP = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_data_GammaConvV1Correction_000001100209366300380000000_01631031009000.root");
	TH1D* histoMCYieldPi02760GeVFinal = (TH1D*)fileNeutralPion2760GeVDataPP->Get("MCYield_Meson_oldBin");
	TFile* fileEta2760GeVDataPP = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Eta_data_GammaConvV1Correction_000001100209366300380000000_01631031009000.root");
	TH1D* histoMCYieldEta2760GeVFinal = (TH1D*)fileEta2760GeVDataPP->Get("MCYield_Meson_oldBin");
	
	cout << "MC spectra added signal 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPPAddSig = new TFile("000001100209366300380000000_01631031009000/2.76TeV/Pi0_MC_GammaConvV1CorrectionHistosAddSig_000001100209366300380000000_01631031009000.root");
	TH1D* histoMCYieldPi02760GeVAddSig = (TH1D*)fileNeutralPion2760GeVDataPPAddSig->Get("MC_Meson_genPt_oldBin");
	TH1F *histoEventQualityMC =         (TH1F*)fileNeutralPion2760GeVDataPPAddSig->Get("NEvents");
	Float_t nEvtMC = GetNEvents(histoEventQualityMC);
	TString rapidityRange = "";
	Double_t deltaRapid =  ReturnRapidityStringAndDouble("01631031009", rapidityRange);
	Double_t scaling = 1./(2.*TMath::Pi());
	ScaleMCYield(histoMCYieldPi02760GeVAddSig,  deltaRapid,  scaling,  nEvtMC,  "Pi0" ,"kFALSE");
	cout << "MC spectra withough weighting added signal 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPPAddSigWOWeighting = new TFile("ppAdditionalInput/0000012002093663003800000_01631031009/2.76TeV/Pi0_data_GammaConvV1Correction_0000012002093663003800000_01631031009.root");
	TH1D* histoMCYieldPi02760GeVAddSigWOWeighting = (TH1D*)fileNeutralPion2760GeVDataPPAddSigWOWeighting->Get("MCYield_Meson_oldBin");
	
	cout << "MC spectra without weighting 2.76TeV GeV" << endl;
	TFile* fileNeutralPion2760GeVDataPPWOWeigthing = new TFile("ppAdditionalInput/0000011002093663003800000_01631031009/2.76TeV/Pi0_data_GammaConvV1Correction_0000011002093663003800000_01631031009.root");
	TH1D* histoMCYieldPi02760GeVWOWeighting = (TH1D*)fileNeutralPion2760GeVDataPPWOWeigthing->Get("MCYield_Meson_oldBin");
	cout << "MC spectra 7TeV GeV" << endl;
	TFile* fileNeutralPion7TeVDataPP = new TFile("900366208010033211360000000900/7TeV/Pi0_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
	TH1D* histoMCYieldPi07TeV = (TH1D*)fileNeutralPion7TeVDataPP->Get("MCYield_Meson_oldBin");
	TFile* fileEta7TeVDataPP = new TFile("900366208010033211360000000900/7TeV/Eta_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
	TH1D* histoMCYieldEta7TeV = (TH1D*)fileEta7TeVDataPP->Get("MCYield_Meson_oldBin");
	cout << "MC spectra 0.9TeV GeV" << endl;
	TFile* fileNeutralPion900GeVDataPP = new TFile("900366208010033211360000000900/900GeV/Pi0_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
	TH1D* histoMCYieldPi0900GeV = (TH1D*)fileNeutralPion900GeVDataPP->Get("MCYield_Meson_oldBin");
	TFile* fileEta900GeVDataPP = new TFile("900366208010033211360000000900/900GeV/Eta_data_AnalysisResultsCorrection_900366208010033211360000000900.root");
	TH1D* histoMCYieldEta900GeV = (TH1D*)fileEta900GeVDataPP->Get("MCYield_Meson_oldBin");

   
   
////////////////////////////////////////////////////////////////// LHC11h file //////////////////////////////////////////////////////////////////////////////   
	TFile*   filePCMPbPbLHC11h = new TFile(nameFilePbPbLHC11h.Data());
     
//////////////////////// Pi0 ////////////////////////////   
	cout << "Pi0 0-10%" << endl;
	TDirectory* directoryPi0PbPbLHC11h0010 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_0-10%"); 

		TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h0010 				 = (TH1D*)directoryPi0PbPbLHC11h0010->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h0010 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0010->Get("Pi0SystError"); 
		
		TH1D* histoPCMPi0MassDataPbPbLHC11h0010 		= (TH1D*)directoryPi0PbPbLHC11h0010->Get("MassPi0");
		histoPCMPi0MassDataPbPbLHC11h0010->Scale(1000.);
		TH1D* histoPCMPi0MassMCPbPbLHC11h0010 			= (TH1D*)directoryPi0PbPbLHC11h0010->Get("TrueMassPi0");
		histoPCMPi0MassMCPbPbLHC11h0010->Scale(1000.);
	
		TH1D* histoPCMPi0WidthDataPbPbLHC11h0010 		= (TH1D*)directoryPi0PbPbLHC11h0010->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0WidthMCPbPbLHC11h0010 			= (TH1D*)directoryPi0PbPbLHC11h0010->Get("TrueFWHMPi0MeV");
		TH1D* histoPi0TrueEffiPtPbPbLHC11h0010 			= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_Efficiency");
		TH1D* histoPi0AcceptPtPbPbLHC11h0010 			= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_Acceptance");
		TH1D* histoMCYieldPi0PtPbPbLHC11h0010		 	= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input_Reweighted");
		TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights0010 	= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input");
		//    TH1D* histoPi0WeightsPbPbLHC11h0010	 	= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Weights");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig0010 	= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0010 = (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Input_AddedSig");
		TH1D* histoPi0RawYieldPbPbLHC11h0010 			= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_RawYieldPerEvent");
	
		TGraphAsymmErrors* graphPi0RAA0010 				= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0010->Get("Pi0RAA"); 
		TGraphAsymmErrors* graphPi0RAASys0010 			= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0010->Get("Pi0RAASys"); 
   
	cout << "0-5%" << endl;
	TDirectory* directoryPi0PbPbLHC11h0005 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_0-5%"); 
	
		TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h0005 				 = (TH1D*)directoryPi0PbPbLHC11h0005->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h0005 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0005->Get("Pi0SystError"); 
		
		TH1D* histoPCMPi0MassDataPbPbLHC11h0005 		= (TH1D*)directoryPi0PbPbLHC11h0005->Get("MassPi0");
		histoPCMPi0MassDataPbPbLHC11h0005->Scale(1000.);
		TH1D* histoPCMPi0MassMCPbPbLHC11h0005 			= (TH1D*)directoryPi0PbPbLHC11h0005->Get("TrueMassPi0");
		histoPCMPi0MassMCPbPbLHC11h0005->Scale(1000.);
		TH1D* histoPCMPi0WidthDataPbPbLHC11h0005 		= (TH1D*)directoryPi0PbPbLHC11h0005->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0WidthMCPbPbLHC11h0005 			= (TH1D*)directoryPi0PbPbLHC11h0005->Get("TrueFWHMPi0MeV");
		TH1D* histoPi0TrueEffiPtPbPbLHC11h0005 			= (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_Efficiency");
		TH1D* histoPi0AcceptPtPbPbLHC11h0005 			= (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_Acceptance");
		TH1D* histoMCYieldPi0PtPbPbLHC11h0005 			= (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input_Reweighted");
		TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights0005 	= (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input");
		//    TH1D* histoPi0WeightsPbPbLHC11h0005 		= (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Weights");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig0005 	= (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0005 = (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_HIJING_Input_AddedSig");
		TH1D* histoPi0RawYieldPbPbLHC11h0005 			= (TH1D*)directoryPi0PbPbLHC11h0005->Get("Pi0_RawYieldPerEvent");

		TGraphAsymmErrors* graphPi0RAA0005 				= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0005->Get("Pi0RAA"); 
		TGraphAsymmErrors* graphPi0RAASys0005 			= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0005->Get("Pi0RAASys"); 
 
	cout << "5-10%" << endl;
	TDirectory* directoryPi0PbPbLHC11h0510 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_5-10%"); 
   
		TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h0510 				 = (TH1D*)directoryPi0PbPbLHC11h0510->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h0510 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h0510->Get("Pi0SystError"); 
		
		TH1D* histoPCMPi0MassDataPbPbLHC11h0510 		= (TH1D*)directoryPi0PbPbLHC11h0510->Get("MassPi0");
		histoPCMPi0MassDataPbPbLHC11h0510->Scale(1000.);
		TH1D* histoPCMPi0MassMCPbPbLHC11h0510 			= (TH1D*)directoryPi0PbPbLHC11h0510->Get("TrueMassPi0");
		histoPCMPi0MassMCPbPbLHC11h0510->Scale(1000.);
		
		TH1D* histoPCMPi0WidthDataPbPbLHC11h0510 		= (TH1D*)directoryPi0PbPbLHC11h0510->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0WidthMCPbPbLHC11h0510 			= (TH1D*)directoryPi0PbPbLHC11h0510->Get("TrueFWHMPi0MeV");
		TH1D* histoPi0TrueEffiPtPbPbLHC11h0510 			= (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_Efficiency");
		TH1D* histoPi0AcceptPtPbPbLHC11h0510 			= (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_Acceptance");
		TH1D* histoMCYieldPi0PtPbPbLHC11h0510 			= (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input_Reweighted");
		TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights0510 	= (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input");
		//    TH1D* histoPi0WeightsPbPbLHC11h0510 		= (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Weights");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig0510 	= (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights0510 = (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_HIJING_Input_AddedSig");
		TH1D* histoPi0RawYieldPbPbLHC11h0510 			= (TH1D*)directoryPi0PbPbLHC11h0510->Get("Pi0_RawYieldPerEvent");
   
	cout << "20-40%" << endl;
	TDirectory* directoryPi0PbPbLHC11h2040 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_20-40%"); 
	
		TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h2040 				 = (TH1D*)directoryPi0PbPbLHC11h2040->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h2040 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2040->Get("Pi0SystError"); 
		
		TH1D* histoPCMPi0MassDataPbPbLHC11h2040 		= (TH1D*)directoryPi0PbPbLHC11h2040->Get("MassPi0");
		histoPCMPi0MassDataPbPbLHC11h2040->Scale(1000.);
		TH1D* histoPCMPi0MassMCPbPbLHC11h2040 			= (TH1D*)directoryPi0PbPbLHC11h2040->Get("TrueMassPi0");
		histoPCMPi0MassDataPbPbLHC11h2040->Scale(1000.);
		
		TH1D* histoPCMPi0WidthDataPbPbLHC11h2040 		= (TH1D*)directoryPi0PbPbLHC11h2040->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0WidthMCPbPbLHC11h2040 			= (TH1D*)directoryPi0PbPbLHC11h2040->Get("TrueFWHMPi0MeV");
		TH1D* histoPi0TrueEffiPtPbPbLHC11h2040 			= (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_Efficiency");
		TH1D* histoPi0AcceptPtPbPbLHC11h2040 			= (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_Acceptance");
		TH1D* histoMCYieldPi0PtPbPbLHC11h2040 			= (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input_Reweighted");
		TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights2040 	= (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input");
		//    TH1D* histoPi0WeightsPbPbLHC11h2040 		= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Weights");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig2040 	= (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2040 = (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_HIJING_Input_AddedSig");
		TH1D* histoPi0RawYieldPbPbLHC11h2040 			= (TH1D*)directoryPi0PbPbLHC11h2040->Get("Pi0_RawYieldPerEvent");

		TGraphAsymmErrors* graphPi0RAA2040 				= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2040->Get("Pi0RAA"); 
		TGraphAsymmErrors* graphPi0RAASys2040 			= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2040->Get("Pi0RAASys"); 
	
	cout << "20-50%" << endl;
	TDirectory* directoryPi0PbPbLHC11h2050 = (TDirectory*)filePCMPbPbLHC11h->Get("Pi0_PbPb_2.76TeV_20-50%"); 
	
		TH1D* histoPCMPi0CorrectedSpecPbPbLHC11h2050 				 = (TH1D*)directoryPi0PbPbLHC11h2050->Get("CorrectedYieldPi0");   
		TGraphAsymmErrors* graphPCMPi0CorrectedSpecSysPbPbLHC11h2050 = (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2050->Get("Pi0SystError"); 
		
		TH1D* histoPCMPi0MassDataPbPbLHC11h2050 		= (TH1D*)directoryPi0PbPbLHC11h2050->Get("MassPi0");
		histoPCMPi0MassDataPbPbLHC11h2050->Scale(1000.);
		TH1D* histoPCMPi0MassMCPbPbLHC11h2050 			= (TH1D*)directoryPi0PbPbLHC11h2050->Get("TrueMassPi0");
		histoPCMPi0MassMCPbPbLHC11h2050->Scale(1000.);
		
		TH1D* histoPCMPi0WidthDataPbPbLHC11h2050 		= (TH1D*)directoryPi0PbPbLHC11h2050->Get("FWHMPi0MeV");
		TH1D* histoPCMPi0WidthMCPbPbLHC11h2050 			= (TH1D*)directoryPi0PbPbLHC11h2050->Get("TrueFWHMPi0MeV");
		TH1D* histoPi0TrueEffiPtPbPbLHC11h2050 			= (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_Efficiency");
		TH1D* histoPi0AcceptPtPbPbLHC11h2050 			= (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_Acceptance");
		TH1D* histoMCYieldPi0PtPbPbLHC11h2050 			= (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input_Reweighted");
		TH1D* histoMCYieldPi0PtPbPbLHC11hWOWeights2050 	= (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input");
		//    TH1D* histoPi0WeightsPbPbLHC11h2050 		= (TH1D*)directoryPi0PbPbLHC11h0010->Get("Pi0_HIJING_Weights");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSig2050 	= (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldPi0PtPbPbLHC11hAddedSigWOWeights2050 = (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_HIJING_Input_AddedSig");
		TH1D* histoPi0RawYieldPbPbLHC11h2050 			= (TH1D*)directoryPi0PbPbLHC11h2050->Get("Pi0_RawYieldPerEvent");
		
		TGraphAsymmErrors* graphPi0RAA2050 				= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2050->Get("Pi0RAA"); 
		TGraphAsymmErrors* graphPi0RAASys2050 			= (TGraphAsymmErrors*)directoryPi0PbPbLHC11h2050->Get("Pi0RAASys"); 
	
   
   
//////////////////////// Eta ////////////////////////////
	cout << "Eta 0-10%" << endl;
	TDirectory* directoryEtaPbPbLHC11h0010 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_0-10%"); 
	
		TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h0010 				 = (TH1D*)directoryEtaPbPbLHC11h0010->Get("CorrectedYieldEta");   
		TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h0010 =        (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0010->Get("EtaSystError"); 
		
		TH1D* histoPCMEtaMassDataPbPbLHC11h0010 		= (TH1D*)directoryEtaPbPbLHC11h0010->Get("MassEta");
		histoPCMEtaMassDataPbPbLHC11h0010->Scale(1000.);
		TH1D* histoPCMEtaMassMCPbPbLHC11h0010	 		= (TH1D*)directoryEtaPbPbLHC11h0010->Get("TrueMassEta");
		histoPCMEtaMassMCPbPbLHC11h0010->Scale(1000.);
		
		TH1D* histoPCMEtaWidthDataPbPbLHC11h0010 		= (TH1D*)directoryEtaPbPbLHC11h0010->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaWidthMCPbPbLHC11h0010 			= (TH1D*)directoryEtaPbPbLHC11h0010->Get("TrueFWHMEtaMeV");
		TH1D* histoEtaTrueEffiPtPbPbLHC11h0010 			= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_Efficiency");
		TH1D* histoEtaAcceptPtPbPbLHC11h0010 			= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_Acceptance");
		TH1D* histoMCYieldEtaPtPbPbLHC11h0010 			= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input_Reweighted");
		TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights0010 	= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input");
		//    TH1D* histoEtaWeightsPbPbLHC11h0010 		= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Weights");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig0010 	= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0010 = (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Input_AddedSig");
		TH1D* histoEtaRawYieldPbPbLHC11h0010 			= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_RawYieldPerEvent");
		
		TGraphAsymmErrors* graphEtaRAA0010 				= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0010->Get("EtaRAA"); 
		TGraphAsymmErrors* graphEtaRAASys0010 			= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0010->Get("EtaRAASys"); 

	cout << "0-5%" << endl;
	TDirectory* directoryEtaPbPbLHC11h0005 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_0-5%"); 

		TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h0005 				 = (TH1D*)directoryEtaPbPbLHC11h0005->Get("CorrectedYieldEta");   
		TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h0005 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0005->Get("EtaSystError"); 

		TH1D* histoPCMEtaMassDataPbPbLHC11h0005 		= (TH1D*)directoryEtaPbPbLHC11h0005->Get("MassEta");
		histoPCMEtaMassDataPbPbLHC11h0005->Scale(1000.);
		TH1D* histoPCMEtaMassMCPbPbLHC11h0005			= (TH1D*)directoryEtaPbPbLHC11h0005->Get("TrueMassEta");
		histoPCMEtaMassMCPbPbLHC11h0005->Scale(1000.);

		TH1D* histoPCMEtaWidthDataPbPbLHC11h0005 		= (TH1D*)directoryEtaPbPbLHC11h0005->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaWidthMCPbPbLHC11h0005 			= (TH1D*)directoryEtaPbPbLHC11h0005->Get("TrueFWHMEtaMeV");
		TH1D* histoEtaTrueEffiPtPbPbLHC11h0005 			= (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_Efficiency");
		TH1D* histoEtaAcceptPtPbPbLHC11h0005 			= (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_Acceptance");
		TH1D* histoMCYieldEtaPtPbPbLHC11h0005 			= (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input_Reweighted");
		TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights0005 	= (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input");
		//    TH1D* histoEtaWeightsPbPbLHC11h0005 = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Weights");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig0005 	= (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0005 = (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_HIJING_Input_AddedSig");
		TH1D* histoEtaRawYieldPbPbLHC11h0005 			= (TH1D*)directoryEtaPbPbLHC11h0005->Get("Eta_RawYieldPerEvent");

		TGraphAsymmErrors* graphEtaRAA0005 				= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0005->Get("EtaRAA"); 
		TGraphAsymmErrors* graphEtaRAASys0005 			= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0005->Get("EtaRAASys"); 

	cout << "5-10%" << endl;
	TDirectory* directoryEtaPbPbLHC11h0510 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_5-10%"); 
	
		TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h0510 				 = (TH1D*)directoryEtaPbPbLHC11h0510->Get("CorrectedYieldEta");   
		TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h0510 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h0510->Get("EtaSystError"); 
		
		TH1D* histoPCMEtaMassDataPbPbLHC11h0510 		= (TH1D*)directoryEtaPbPbLHC11h0510->Get("MassEta");
		histoPCMEtaMassDataPbPbLHC11h0510->Scale(1000.);
		TH1D* histoPCMEtaMassMCPbPbLHC11h0510 			= (TH1D*)directoryEtaPbPbLHC11h0510->Get("TrueMassEta");
		histoPCMEtaMassMCPbPbLHC11h0510->Scale(1000.);
		
		TH1D* histoPCMEtaWidthDataPbPbLHC11h0510 		= (TH1D*)directoryEtaPbPbLHC11h0510->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaWidthMCPbPbLHC11h0510 			= (TH1D*)directoryEtaPbPbLHC11h0510->Get("TrueFWHMEtaMeV");
		TH1D* histoEtaTrueEffiPtPbPbLHC11h0510 			= (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_Efficiency");
		TH1D* histoEtaAcceptPtPbPbLHC11h0510 			= (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_Acceptance");
		TH1D* histoMCYieldEtaPtPbPbLHC11h0510 			= (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input_Reweighted");
		TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights0510 	= (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input");
		//    TH1D* histoEtaWeightsPbPbLHC11h0510 		= (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Weights");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig0510 	= (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights0510 = (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_HIJING_Input_AddedSig");
		TH1D* histoEtaRawYieldPbPbLHC11h0510 			= (TH1D*)directoryEtaPbPbLHC11h0510->Get("Eta_RawYieldPerEvent");
	
	cout << "20-40%" << endl;
	TDirectory* directoryEtaPbPbLHC11h2040 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_20-40%");
	
		TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h2040 				 = (TH1D*)directoryEtaPbPbLHC11h2040->Get("CorrectedYieldEta");   
		TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h2040 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2040->Get("EtaSystError"); 
		
		TH1D* histoPCMEtaMassDataPbPbLHC11h2040 		= (TH1D*)directoryEtaPbPbLHC11h2040->Get("MassEta");
		histoPCMEtaMassDataPbPbLHC11h2040->Scale(1000.);
		TH1D* histoPCMEtaMassMCPbPbLHC11h2040 			= (TH1D*)directoryEtaPbPbLHC11h2040->Get("TrueMassEta");
		histoPCMEtaMassMCPbPbLHC11h2040->Scale(1000.);
		
		TH1D* histoPCMEtaWidthDataPbPbLHC11h2040 		= (TH1D*)directoryEtaPbPbLHC11h2040->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaWidthMCPbPbLHC11h2040 			= (TH1D*)directoryEtaPbPbLHC11h2040->Get("TrueFWHMEtaMeV");
		TH1D* histoEtaTrueEffiPtPbPbLHC11h2040 			= (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_Efficiency");
		TH1D* histoEtaAcceptPtPbPbLHC11h2040 			= (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_Acceptance");
		TH1D* histoMCYieldEtaPtPbPbLHC11h2040 			= (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input_Reweighted");
		TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights2040 	= (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input");
		//    TH1D* histoEtaWeightsPbPbLHC11h0010 		= (TH1D*)directoryEtaPbPbLHC11h0010->Get("Eta_HIJING_Weights");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig2040 	= (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2040 = (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_HIJING_Input_AddedSig");
		TH1D* histoEtaRawYieldPbPbLHC11h2040 			= (TH1D*)directoryEtaPbPbLHC11h2040->Get("Eta_RawYieldPerEvent");
		
		TGraphAsymmErrors* graphEtaRAA2040 				= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2040->Get("EtaRAA"); 
		TGraphAsymmErrors* graphEtaRAASys2040 			= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2040->Get("EtaRAASys"); 
	
	cout << "20-50%" << endl;
	TDirectory* directoryEtaPbPbLHC11h2050 = (TDirectory*)filePCMPbPbLHC11h->Get("Eta_PbPb_2.76TeV_20-50%"); 
	
		TH1D* histoPCMEtaCorrectedSpecPbPbLHC11h2050 				 = (TH1D*)directoryEtaPbPbLHC11h2050->Get("CorrectedYieldEta");   
		TGraphAsymmErrors* graphPCMEtaCorrectedSpecSysPbPbLHC11h2050 = (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2050->Get("EtaSystError"); 

		TH1D* histoPCMEtaMassDataPbPbLHC11h2050 		= (TH1D*)directoryEtaPbPbLHC11h2050->Get("MassEta");
		histoPCMEtaMassDataPbPbLHC11h2050->Scale(1000.);
		TH1D* histoPCMEtaMassMCPbPbLHC11h2050 			= (TH1D*)directoryEtaPbPbLHC11h2050->Get("TrueMassEta");
		histoPCMEtaMassMCPbPbLHC11h2050->Scale(1000.);

		TH1D* histoPCMEtaWidthDataPbPbLHC11h2050 		= (TH1D*)directoryEtaPbPbLHC11h2050->Get("FWHMEtaMeV");
		TH1D* histoPCMEtaWidthMCPbPbLHC11h2050 			= (TH1D*)directoryEtaPbPbLHC11h2050->Get("TrueFWHMEtaMeV");
		TH1D* histoEtaTrueEffiPtPbPbLHC11h2050 			= (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_Efficiency");
		TH1D* histoEtaAcceptPtPbPbLHC11h2050 			= (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_Acceptance");
		TH1D* histoMCYieldEtaPtPbPbLHC11h2050 			= (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input_Reweighted");
		TH1D* histoMCYieldEtaPtPbPbLHC11hWOWeights2050 	= (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSig2050 	= (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input_Reweighted_AddedSig");
		TH1D* histoMCYieldEtaPtPbPbLHC11hAddedSigWOWeights2050 = (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_HIJING_Input_AddedSig");
		TH1D* histoEtaRawYieldPbPbLHC11h2050 			= (TH1D*)directoryEtaPbPbLHC11h2050->Get("Eta_RawYieldPerEvent");

		TGraphAsymmErrors* graphEtaRAA2050 				= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2050->Get("EtaRAA"); 
		TGraphAsymmErrors* graphEtaRAASys2050 			= (TGraphAsymmErrors*)directoryEtaPbPbLHC11h2050->Get("EtaRAASys");    
   
		
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
	// ************************** 				Efficiency 				**************************************
	// *******************************************************************************************************   
   
	TCanvas* canvasEffEta = new TCanvas("canvasEffEta","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasEffEta, 0.1, 0.02, 0.035, 0.09);
	canvasEffEta->SetLogy();
	TH2F * histo2DEffEta = new TH2F("histo2DEffEta","histo2DEffEta",1000,0,10.5,2000,1e-4,1e-2 );
	SetStyleHistoTH2ForGraphs(histo2DEffEta, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #eta}  ",0.03,0.04, 0.03,0.04, 1.,1.);
	histo2DEffEta->Draw("copy");
	
// 	DrawGammaSetMarker(histoTrueEffPtEta2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);  
// 	histoTrueEffPtEta2760GeV->DrawCopy("pe1,same");    
	
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

	TLegend* legendEffiEtaLHC11h = new TLegend(0.35,0.15,0.85,0.3);
	legendEffiEtaLHC11h->SetFillColor(0);
	legendEffiEtaLHC11h->SetLineColor(0);
	legendEffiEtaLHC11h->SetTextSize(0.035);
	legendEffiEtaLHC11h->SetNColumns(1);
// 	legendEffiEtaLHC11h->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendEffiEtaLHC11h->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffiEtaLHC11h->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0010,"0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendEffiEtaLHC11h->AddEntry(histoEtaTrueEffiPtPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffiEtaLHC11h->AddEntry(histoEtaTrueEffiPtPbPbLHC11h2050,"20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendEffiEtaLHC11h->AddEntry(histoTrueEffPtEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
	legendEffiEtaLHC11h->Draw();	
	
	canvasEffEta->SaveAs(Form("%s/EfficiencyEtaLHC11h.%s",outputDir.Data(),suffix.Data()));
	
	
	TCanvas* canvasEffPi0 = new TCanvas("canvasEffPi0","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasEffPi0, 0.1, 0.02, 0.035, 0.09);
	canvasEffPi0->SetLogy();
	TH2F * histo2DEffPi0 = new TH2F("histo2DEffPi0","histo2DEffPi0",1000,0.,15.,2000,2e-6,1.e-2   );
	SetStyleHistoTH2ForGraphs(histo2DEffPi0, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco, #pi^{0}}  ",0.03,0.04, 0.03,0.04, 1.,1.);
	histo2DEffPi0->Draw("copy");

// 	DrawGammaSetMarker(histoTrueEffPtPi02760GeV, markerStyleSpectrum2760GeVMC-1, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);   
// 	histoTrueEffPtPi02760GeV->DrawCopy("pe1,same");
	
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

	TLegend* legendEffiPi0LHC11h = new TLegend(0.35,0.15,0.85,0.3);
	legendEffiPi0LHC11h->SetFillColor(0);
	legendEffiPi0LHC11h->SetLineColor(0);
	legendEffiPi0LHC11h->SetTextSize(0.035);
// 	legendEffiPi0LHC11h->AddEntry(histoPi0TrueEffiPtPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendEffiPi0LHC11h->AddEntry(histoPi0TrueEffiPtPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffiPi0LHC11h->AddEntry(histoPi0TrueEffiPtPbPbLHC11h0010,"0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendEffiPi0LHC11h->AddEntry(histoPi0TrueEffiPtPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffiPi0LHC11h->AddEntry(histoPi0TrueEffiPtPbPbLHC11h2050,"20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendEffiPi0LHC11h->AddEntry(histoTrueEffPtPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
	legendEffiPi0LHC11h->Draw();
// 	legendEffiPi0LHC11h->Draw();

	canvasEffPi0->SaveAs(Form("%s/EfficiencyPi0LHC11h.%s",outputDir.Data(),suffix.Data()));

	
	canvasEffPi0->cd();
	histo2DEffPi0->GetYaxis()->SetTitle("#epsilon_{reco}  ");
	histo2DEffPi0->Draw("copy");
	
	DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h0010, markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010); 
	histoPi0TrueEffiPtPbPbLHC11h0010->DrawCopy("e1,same");
	DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h0010, markerStylePbPb0010+1, markerSizePbPb0010, colorCombPbPb1020, colorCombPbPb1020); 
	histoEtaTrueEffiPtPbPbLHC11h0010->DrawCopy("e1,same");  

	TLegend* legendEffi0010 = new TLegend(0.3,0.15,0.85,0.3);
	legendEffi0010->SetFillColor(0);
	legendEffi0010->SetLineColor(0);
	legendEffi0010->SetTextSize(0.035);
	legendEffi0010->AddEntry(histoPi0TrueEffiPtPbPbLHC11h0010,"#pi^{0}  0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffi0010->AddEntry(histoEtaTrueEffiPtPbPbLHC11h0010,"#eta   0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffi0010->Draw();
	
	canvasEffPi0->SaveAs(Form("%s/EfficiencyPi0andEta0010_LHC11h.%s",outputDir.Data(),suffix.Data()));
	
	canvasEffPi0->cd();
	histo2DEffPi0->Draw("copy");
	
	DrawGammaSetMarker(histoPi0TrueEffiPtPbPbLHC11h2050, markerStylePbPb2040, markerSizePbPb2040, colorCombPbPb2040, colorCombPbPb2040); 
	histoPi0TrueEffiPtPbPbLHC11h2050->DrawCopy("e1,same");  
	DrawGammaSetMarker(histoEtaTrueEffiPtPbPbLHC11h2050, markerStylePbPb2040+1, markerSizePbPb2040, colorCombPbPb4050, colorCombPbPb4050); 
	histoEtaTrueEffiPtPbPbLHC11h2050->DrawCopy("e1,same");  
	
	TLegend* legendEffi2050 = new TLegend(0.3,0.15,0.85,0.3);
	legendEffi2050->SetFillColor(0);
	legendEffi2050->SetLineColor(0);
	legendEffi2050->SetTextSize(0.035);
	legendEffi2050->AddEntry(histoPi0TrueEffiPtPbPbLHC11h2050,"#pi^{0}  20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffi2050->AddEntry(histoEtaTrueEffiPtPbPbLHC11h2050,"#eta   20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffi2050->Draw();
	canvasEffPi0->SaveAs(Form("%s/EfficiencyPi0andEta2050_LHC11h.%s",outputDir.Data(),suffix.Data()));


	Double_t arrayBoundariesX1_4[3];
	Double_t arrayBoundariesY1_4[2];
	Double_t relativeMarginsX[3];
	Double_t relativeMarginsY[3];
	ReturnCorrectValuesForCanvasScaling(1000,500, 2, 1,0.06, 0.005, 0.005,0.09,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

	TCanvas* canvasEffRatio = new TCanvas("canvasEffRatio","",0,0,1000,500);  // gives the page size
	DrawGammaCanvasSettings( canvasEffRatio,  0.13, 0.02, 0.03, 0.06);

	TPad* pad2PartEffiRatio1 = new TPad("pad2PartEffiRatio1", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( pad2PartEffiRatio1, relativeMarginsX[0], relativeMarginsX[1], relativeMarginsY[0], relativeMarginsY[2]);
	pad2PartEffiRatio1->Draw();

	TPad* pad2PartEffiRatio3 = new TPad("pad2PartEffiRatio3", "", arrayBoundariesX1_4[1], arrayBoundariesY1_4[1], arrayBoundariesX1_4[2], arrayBoundariesY1_4[0],-1, -1, -2);
	DrawGammaPadSettings( pad2PartEffiRatio3, relativeMarginsX[1], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[2]);
	pad2PartEffiRatio3->Draw();

	Double_t margin = relativeMarginsX[0]*0.8*1000;
	Double_t textsizeLabels1 = 0;
	Double_t textsizeFac1 = 0;
	Double_t textsizeLabels2 = 0;
	Double_t textsizeFac2 = 0;

	ReturnCorrectValuesTextSize(pad2PartEffiRatio1,textsizeLabels1, textsizeFac1, 22, margin);
	ReturnCorrectValuesTextSize(pad2PartEffiRatio3,textsizeLabels2, textsizeFac2, 22, margin);

	TH2F * histo2DEffRatio2 = new TH2F("histo2DEffRatio2","histo2DEffRatio2",1000,0.,15.,1000,0.5,2.);
	SetStyleHistoTH2ForGraphs(histo2DEffRatio2, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco} #frac{20-50% }{0-10%}  ",0.85*textsizeLabels1, textsizeLabels1,
								  0.85*textsizeLabels1, textsizeLabels1, 0.8,0.2/(textsizeFac1*margin), 512, 505);
	histo2DEffRatio2->GetYaxis()->SetRangeUser(0.5,1.8);
	
	TH2F* histo2DEffRatio = new TH2F("histo2DEffRatio","histo2DEffRatio",1000,0.,10.5,1000,0.5,2.);
	SetStyleHistoTH2ForGraphs(histo2DEffRatio, "#it{p}_{T} (GeV/#it{c})","#epsilon_{reco} #frac{20-50% }{0-10%}", 0.85*textsizeLabels2, textsizeLabels2,
								  0.85*textsizeLabels2, textsizeLabels2, 0.8,0.25/(textsizeFac2*margin), 512, 505);
	histo2DEffRatio->GetYaxis()->SetRangeUser(0.5,1.8);

	pad2PartEffiRatio1->cd();
	histo2DEffRatio2->DrawCopy();

	TH1D* ratioTrueEffiPtPi0 = (TH1D*)histoPi0TrueEffiPtPbPbLHC11h2050->Clone("ratioTrueEffiPtPi0");
	ratioTrueEffiPtPi0->Divide(ratioTrueEffiPtPi0,histoPi0TrueEffiPtPbPbLHC11h0010,1.,1.,"B");

	DrawGammaSetMarker(ratioTrueEffiPtPi0, markerStylePbPb0010, 1.5, colorCombPbPb2040, colorCombPbPb2040); 
	ratioTrueEffiPtPi0->DrawCopy("e1,same");

	TLegend* legendEffiRatioPi0 = new TLegend(0.3,0.13,0.8,0.27);
	legendEffiRatioPi0->SetFillColor(0);
	legendEffiRatioPi0->SetLineColor(0);
	legendEffiRatioPi0->SetTextSize(0.035);
	legendEffiRatioPi0->AddEntry(ratioTrueEffiPtPi0,"#pi^{0}, #epsilon_{reco} #frac{20-50% }{0-10%} Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffiRatioPi0->Draw();

	DrawGammaLines(0., 20.5 , 1, 1 ,1,kGray, 2);
	
	histo2DEffRatio2->Draw("axis,same");
	pad2PartEffiRatio1->Update();
	pad2PartEffiRatio3->cd();
	histo2DEffRatio->DrawCopy();
			
	TH1D* ratioTrueEffiPtEta = (TH1D*)histoEtaTrueEffiPtPbPbLHC11h2050->Clone("ratioTrueEffiPtEta");
	ratioTrueEffiPtEta->Divide(ratioTrueEffiPtEta,histoEtaTrueEffiPtPbPbLHC11h0010,1.,1.,"B");

	DrawGammaSetMarker(ratioTrueEffiPtEta, markerStylePbPb0010, 1.5, colorCombPbPb2040, colorCombPbPb2040); 
	ratioTrueEffiPtEta->DrawCopy("e1,same");

	TLegend* legendEffiRatioEta = new TLegend(0.23,0.13,0.8,0.27);
	legendEffiRatioEta->SetFillColor(0);
	legendEffiRatioEta->SetLineColor(0);
	legendEffiRatioEta->SetTextSize(0.035);
	legendEffiRatioEta->AddEntry(ratioTrueEffiPtEta,"#eta, #epsilon_{reco} #frac{20-50% }{0-10%} Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendEffiRatioEta->Draw();

	DrawGammaLines(0., 19.5 , 1, 1 ,1, kGray, 2);	
		
	histo2DEffRatio->Draw("axis,same");
	pad2PartEffiRatio3->Update();

	canvasEffRatio->Update();	
	canvasEffRatio->SaveAs(Form("%s/EfficiencyRatio_LHC11h.%s",outputDir.Data(),suffix.Data()));
	

	
	// *******************************************************************************************************
	// ************************** 				Acceptance 				**************************************
	// *******************************************************************************************************   
	
	TCanvas* canvasAcceptanceEta = new TCanvas("canvasAcceptanceEta","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasAcceptanceEta, 0.1, 0.02, 0.035, 0.09);
	TH2F * histo2DAcceptanceEta = new TH2F("histo2DAcceptanceEta","histo2DAcceptanceEta",1000,0,10.5,2000,0.4,1.02 );
	SetStyleHistoTH2ForGraphs(histo2DAcceptanceEta, "#it{p}_{T} (GeV/#it{c})","A_{#eta}",0.03,0.04, 0.03,0.04, 1.,1.);
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

	TLegend* legendAcceptanceEta = new TLegend(0.4,0.15,0.9,0.32); //0.25,0.13,0.93,0.43);
	legendAcceptanceEta->SetFillColor(0);
	legendAcceptanceEta->SetLineColor(0);
	legendAcceptanceEta->SetTextSize(0.035);
	legendAcceptanceEta->SetMargin(0.2);
// 	legendAcceptanceEta->AddEntry(histoEtaAcceptPtPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendAcceptanceEta->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
// 	legendAcceptanceEta->AddEntry(histoEtaAcceptPtPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendAcceptanceEta->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
	legendAcceptanceEta->AddEntry(histoEtaAcceptPtPbPbLHC11h0010,"0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendAcceptanceEta->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
	legendAcceptanceEta->AddEntry((TObject*)0, "","");
	// 	legendAcceptanceEta->AddEntry(histoEtaAcceptPtPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendAcceptanceEta->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
	legendAcceptanceEta->AddEntry(histoEtaAcceptPtPbPbLHC11h2050,"20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendAcceptanceEta->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
// 	legendAcceptanceEta->AddEntry(histoAccEta2760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
// 	legendAcceptanceEta->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.8","");
	legendAcceptanceEta->Draw();
	
	canvasAcceptanceEta->SaveAs(Form("%s/AcceptanceEtaLHC11h.%s",outputDir.Data(),suffix.Data()));
	

	TCanvas* canvasAcceptancePi0 = new TCanvas("canvasAcceptancePi0","",200,10,1350,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasAcceptancePi0, 0.1, 0.02, 0.035, 0.09);
	TH2F * histo2DAcceptancePi0;
	histo2DAcceptancePi0 = new TH2F("histo2DAcceptancePi0","histo2DAcceptancePi0",1000,0,15.,2000,0.7,1.02   );
	SetStyleHistoTH2ForGraphs(histo2DAcceptancePi0, "#it{p}_{T} (GeV/#it{c})","A_{#pi^{0}}",0.03,0.04, 0.03,0.04, 1.,1.);
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

	TLegend* legendAcceptancePi0 = new TLegend(0.4,0.15,0.9,0.32); //0.25,0.13,0.93,0.43);
	legendAcceptancePi0->SetFillColor(0);
	legendAcceptancePi0->SetLineColor(0);
	legendAcceptancePi0->SetTextSize(0.035);
	legendAcceptancePi0->SetMargin(0.2);
// 	legendAcceptancePi0->AddEntry(histoPi0AcceptPtPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendAcceptancePi0->AddEntry((TObject*)0, " |y| < 0.6","");
// 	legendAcceptancePi0->AddEntry(histoPi0AcceptPtPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendAcceptancePi0->AddEntry((TObject*)0, " |y| < 0.6","");
	legendAcceptancePi0->AddEntry(histoPi0AcceptPtPbPbLHC11h0010,"0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendAcceptancePi0->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
	legendAcceptancePi0->AddEntry((TObject*)0, "","");
// 	legendAcceptancePi0->AddEntry(histoPi0AcceptPtPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendAcceptancePi0->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
	legendAcceptancePi0->AddEntry(histoPi0AcceptPtPbPbLHC11h2050,"20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendAcceptancePi0->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.85","");
// 	legendAcceptancePi0->AddEntry(histoAccPi02760GeV,"pp #sqrt{#it{s}} = 2.76 TeV","p");
// 	legendAcceptancePi0->AddEntry((TObject*)0, " |#eta| < 0.9, |y| < 0.8","");
	legendAcceptancePi0->Draw();
	
	canvasAcceptancePi0->SaveAs(Form("%s/AcceptancePi0LHC11h.%s",outputDir.Data(),suffix.Data()));

	
   
   
	// *******************************************************************************************************
	// ************************** 				Raw yields				**************************************
	// *******************************************************************************************************
	
	TCanvas* canvasRawYieldPi0PbPb = new TCanvas("canvasRawYieldPi0PbPb","",200,10,1350*1.4,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasRawYieldPi0PbPb, 0.12, 0.02, 0.035, 0.09);
	canvasRawYieldPi0PbPb->SetLogy();
	//canvasRawYieldPi0PbPb->SetLogx();
	TH2F * histo2DRawPi0PbPb = new TH2F("histo2DRawPi0PbPb","histo2DRawPi0PbPb",1000,0.,15.,2000,1.e-7,1e0 );
	SetStyleHistoTH2ForGraphs(histo2DRawPi0PbPb, "#it{p}_{T} (GeV/#it{c})","#frac{N_{raw}^{#pi^{0}}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
	histo2DRawPi0PbPb->GetYaxis()->SetRangeUser(1.e-7,2e-1 );
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

	labelRawPi0PbPb->Draw();
	TLegend* legendRawYieldPi0PbPb = new TLegend(0.55,0.75,0.75,0.88);
	legendRawYieldPi0PbPb->SetFillColor(0);
	legendRawYieldPi0PbPb->SetLineColor(0);
	legendRawYieldPi0PbPb->SetTextSize(0.035);
// 	legendRawYieldPi0PbPb->AddEntry(histoPi0RawYieldPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendRawYieldPi0PbPb->AddEntry(histoPi0RawYieldPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendRawYieldPi0PbPb->AddEntry(histoPi0RawYieldPbPbLHC11h0010,"0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendRawYieldPi0PbPb->AddEntry(histoPi0RawYieldPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendRawYieldPi0PbPb->AddEntry(histoPi0RawYieldPbPbLHC11h2050,"20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendRawYieldPi0PbPb->Draw();
	
	canvasRawYieldPi0PbPb->SaveAs(Form("%s/RawYieldPi0PbPbLHC11h.%s",outputDir.Data(),suffix.Data()));
	

	TCanvas* canvasRawYieldEtaPbPb = new TCanvas("canvasRawYieldEtaPbPb","",200,10,1350*1.4,1350);  // gives the page size
	DrawGammaCanvasSettings( canvasRawYieldEtaPbPb, 0.12, 0.02, 0.035, 0.09);
	canvasRawYieldEtaPbPb->SetLogy();
	//       canvasRawYieldEtaPbPb->SetLogx();
	TH2F * histo2DRawEtaPbPb = new TH2F("histo2DRawEtaPbPb","histo2DRawEtaPbPb",1000,0.,10.5,2000,1.e-7,2e-2 );
	SetStyleHistoTH2ForGraphs(histo2DRawEtaPbPb, "#it{p}_{T} (GeV/#it{c})","#frac{N_{raw}^{#eta}}{N_{evt} d#it{p}_{T}}",0.035,0.04, 0.035,0.04, 1.,1.3);
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

	labelRawEtaPbPb->Draw();
	TLegend* legendRawYieldEtaPbPb = new TLegend(0.55,0.75,0.75,0.88);
	legendRawYieldEtaPbPb->SetFillColor(0);
	legendRawYieldEtaPbPb->SetLineColor(0);
	legendRawYieldEtaPbPb->SetTextSize(0.035);
// 	legendRawYieldEtaPbPb->AddEntry(histoEtaRawYieldPbPbLHC11h0005,"0-5% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendRawYieldEtaPbPb->AddEntry(histoEtaRawYieldPbPbLHC11h0510,"5-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendRawYieldEtaPbPb->AddEntry(histoEtaRawYieldPbPbLHC11h0010,"0-10% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
// 	legendRawYieldEtaPbPb->AddEntry(histoEtaRawYieldPbPbLHC11h2040,"20-40% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendRawYieldEtaPbPb->AddEntry(histoEtaRawYieldPbPbLHC11h2050,"20-50% Pb-Pb, #sqrt{s_{_{NN}}} = 2.76 TeV","p");
	legendRawYieldEtaPbPb->Draw();
	
	canvasRawYieldEtaPbPb->SaveAs(Form("%s/RawYieldEtaPbPbLHC11h.%s",outputDir.Data(),suffix.Data()));
	
	
	
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
			
	DrawGammaSetMarker(histoPCMWidthDataPP2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
	histoPCMWidthDataPP2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMWidthMCPP2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoPCMWidthMCPP2760GeV->DrawCopy("same,e1,p"); 

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

	DrawGammaSetMarker(histoPCMMassDataPP2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);               
	histoPCMMassDataPP2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMMassMCPP2760GeV , markerStyleSpectrum2760GeVMC , markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);               
	histoPCMMassMCPP2760GeV->DrawCopy("same,e1,p"); 

	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1.);
	TLatex *labelLegendDMass = new TLatex(0.92,0.9,"d)");
	SetStyleTLatex( labelLegendDMass, 0.075,4);
	labelLegendDMass->Draw();

	pad6PartMassWidth2->Update();

	pad6PartMassWidth5->cd();
	pad6PartMassWidth5->SetLogx();
	histo2DPi0FWHM->DrawCopy();
			
	DrawGammaSetMarker(histoPCMPi0WidthDataPbPbLHC11h0010, markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
	histoPCMPi0WidthDataPbPbLHC11h0010->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMPi0WidthMCPbPbLHC11h0010, markerStylePbPb0005MC, markerSizePbPb0005, colorCombMCPbPb0005, colorCombMCPbPb0005);
	histoPCMPi0WidthMCPbPbLHC11h0010->DrawCopy("same,e1,p"); 

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
	histoPCMPi0MassDataPbPbLHC11h0010->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMPi0MassMCPbPbLHC11h0010 , markerStylePbPb0005MC , markerSizePbPb0005, colorCombMCPbPb0005, colorCombMCPbPb0005);                
	histoPCMPi0MassMCPbPbLHC11h0010->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1);
	TLatex *labelLegendFMass = new TLatex(0.89,0.9,"f)");
	SetStyleTLatex( labelLegendFMass, 0.075,4);
	labelLegendFMass->Draw();

	pad6PartMassWidth6->Update();

	pad6PartMassWidth3->cd();
	pad6PartMassWidth3->SetLogx();
	histo2DPi0FWHM->DrawCopy();
		
	DrawGammaSetMarker(histoPCMPi0WidthDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
	histoPCMPi0WidthDataPbPbLHC11h2050->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMPi0WidthMCPbPbLHC11h2050, markerStylePbPb6080MC, markerSizePbPb6080, colorCombMCPbPb6080, colorCombMCPbPb6080);
	histoPCMPi0WidthMCPbPbLHC11h2050->DrawCopy("same,e1,p"); 

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
	histoPCMPi0MassDataPbPbLHC11h2050->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMPi0MassMCPbPbLHC11h2050, markerStylePbPb0005MC , markerSizePbPb6080, colorCombMCPbPb6080, colorCombMCPbPb6080);                
	histoPCMPi0MassMCPbPbLHC11h2050->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,1.,1);
	TLatex *labelLegendEMass = new TLatex(0.91,0.9,"e)");
	SetStyleTLatex( labelLegendEMass, 0.075,4);
	labelLegendEMass->Draw();

	pad6PartMassWidth4->Update();

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
			
	DrawGammaSetMarker(histoEtaWidthData2760GeV, markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);
	histoEtaWidthData2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaWidthMC2760GeV, markerStyleSpectrum2760GeVMC, markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);
	histoEtaWidthMC2760GeV->DrawCopy("same,e1,p"); 

	labelMassPi0PP->Draw();
// 	TLatex *labelLegendAMass = new TLatex(0.92,0.88,"a)");
// 	SetStyleTLatex( labelLegendAMass, 0.08,4);

	labelLegendAMass->Draw();

	pad6PartEtaMassWidth1LHC11h->Update();
	pad6PartEtaMassWidth2LHC11h->cd();
	pad6PartEtaMassWidth2LHC11h->SetLogx();
	histo2DEtaMassLHC11h->DrawCopy();

	DrawGammaSetMarker(histoEtaMassData2760GeV , markerStyleSpectrum2760GeV, markerSizePi0PP2760GeV, colorPi02760GeV, colorPi02760GeV);               
	histoEtaMassData2760GeV->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoEtaMassMC2760GeV , markerStyleSpectrum2760GeVMC , markerSizePi0PP2760GeV, colorMCPythiaPP2760GeV, colorMCPythiaPP2760GeV);               
	histoEtaMassMC2760GeV->DrawCopy("same,e1,p"); 

	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1.);
	labelLegendDMass->Draw();

	pad6PartEtaMassWidth2LHC11h->Update();

	pad6PartEtaMassWidth5LHC11h->cd();
	pad6PartEtaMassWidth5LHC11h->SetLogx();
	histo2DEtaFWHMLHC11h->DrawCopy();
			
	DrawGammaSetMarker(histoPCMEtaWidthDataPbPbLHC11h0010,markerStylePbPb0005, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
	histoPCMEtaWidthDataPbPbLHC11h0010->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMEtaWidthMCPbPbLHC11h0010, markerStylePbPb0005MC, markerSizePbPb0005, colorCombPbPb0005, colorCombPbPb0005);
	histoPCMEtaWidthMCPbPbLHC11h0010->DrawCopy("same,e1,p"); 

	TLatex *labelMassEtaPbPbLHC11h0010 = new TLatex(0.05,0.9,collisionSystemPbPb0010.Data());
	SetStyleTLatex( labelMassEtaPbPbLHC11h0010, 0.062,4);
	labelMassEtaPbPbLHC11h0010->Draw();
	labelLegendBMass->Draw();
	
	labelLegendTypeData->Draw();
	labelLegendTypeMC->Draw();

	pad6PartEtaMassWidth5LHC11h->Update();
	pad6PartEtaMassWidth6LHC11h->cd();
	pad6PartEtaMassWidth6LHC11h->SetLogx();
	histo2DEtaMassLHC11h->DrawCopy();

	DrawGammaSetMarker(histoPCMEtaMassDataPbPbLHC11h0010,markerStylePbPb0010, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);
	histoPCMEtaMassDataPbPbLHC11h0010->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMEtaMassMCPbPbLHC11h0010, markerStylePbPb0010MC, markerSizePbPb0010, colorCombPbPb0010, colorCombPbPb0010);
	histoPCMEtaMassMCPbPbLHC11h0010->DrawCopy("same,e1,p");  
	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1);
	labelLegendFMass->Draw();

	pad6PartEtaMassWidth6LHC11h->Update();

	pad6PartEtaMassWidth3LHC11h->cd();
	pad6PartEtaMassWidth3LHC11h->SetLogx();
	histo2DEtaFWHMLHC11h->DrawCopy();
	DrawGammaSetMarker(histoPCMEtaWidthDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
	histoPCMEtaWidthDataPbPbLHC11h2050->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMEtaWidthMCPbPbLHC11h2050, markerStylePbPb6080MC, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);
	histoPCMEtaWidthMCPbPbLHC11h2050->DrawCopy("same,e1,p"); 

	TLatex *labelMassEtaPbPb2050 = new TLatex(0.05,0.9,collisionSystemPbPb2050.Data());
	SetStyleTLatex( labelMassEtaPbPb2050, 0.062,4);
	labelMassEtaPbPb2050->Draw();
	labelLegendCMass->Draw();

	pad6PartEtaMassWidth3LHC11h->Update();
	pad6PartEtaMassWidth4LHC11h->cd();
	pad6PartEtaMassWidth4LHC11h->SetLogx();
	histo2DEtaMassLHC11h->DrawCopy();
	
	DrawGammaSetMarker(histoPCMEtaMassDataPbPbLHC11h2050, markerStylePbPb6080, markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);              
	histoPCMEtaMassDataPbPbLHC11h2050->DrawCopy("same,e1,p"); 
	DrawGammaSetMarker(histoPCMEtaMassMCPbPbLHC11h2050, markerStylePbPb6080MC , markerSizePbPb6080, colorCombPbPb6080, colorCombPbPb6080);            
	histoPCMEtaMassMCPbPbLHC11h2050->DrawCopy("same,e1,p"); 
	DrawGammaLines(0.35, 20. , mesonMassExpectEta*1000., mesonMassExpectEta*1000.,1.,1);
	labelLegendEMass->Draw();

	pad6PartEtaMassWidth4LHC11h->Update();

	canvas6PartEtaMassWidthLHC11h->Update();  
	canvas6PartEtaMassWidthLHC11h->SaveAs(Form("%s/MassWidth_Eta.%s",outputDir.Data(),suffix.Data()));
	delete pad6PartEtaMassWidth1LHC11h; 
	delete pad6PartEtaMassWidth2LHC11h; 
	delete pad6PartEtaMassWidth3LHC11h; 
	delete pad6PartEtaMassWidth4LHC11h; 
	delete canvas6PartEtaMassWidthLHC11h;


	
	// *******************************************************************************************************
	// ************************** 			Invariant yields			**************************************
	// *******************************************************************************************************   
	
	TCanvas* canvasEtaSpectra = new TCanvas("canvasEtaSpectra","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasEtaSpectra,  0.13, 0.01, 0.015, 0.08);
	canvasEtaSpectra->SetLogy();
	canvasEtaSpectra->SetLogx();
	TH2F * histo2DEtaSpectra = new TH2F("histo2DEtaSpectra","histo2DEtaSpectra",1000,0.1,11.,1000,1e-8,2e2 );
	SetStyleHistoTH2ForGraphs(histo2DEtaSpectra, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
	histo2DEtaSpectra->DrawCopy(); 
	
// 	graphCorrectedYieldSysEta2760GeV = ScaleGraph(graphCorrectedYieldSysEta2760GeV,1e1);
// 	histoCorrectedYieldEta2760GeV->Scale(1e1);
// 	DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldSysEta2760GeV, markerStylePP,markerSizePP, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
// 	graphCorrectedYieldSysEta2760GeV->Draw("E2same");
// 	DrawGammaSetMarker(histoCorrectedYieldEta2760GeV, markerStylePP,markerSizePP, kBlack , kBlack);
// 	histoCorrectedYieldEta2760GeV->Draw("p,same,e1");
// 	
// 	histoEtaFromCocktailpp2760GeV->SetLineWidth(2.);
// 	histoEtaFromCocktailpp2760GeV->SetLineColor(kBlack);
// 	histoEtaFromCocktailpp2760GeV->Draw("hist, l,same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
//  histoPCMEtaCorrectedSpecPbPbLHC11h0005->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005, widthLinesBoxes, kTRUE);
// 	graphPCMEtaCorrectedSpecSysPbPbLHC11h0005->Draw("E2same");
	
	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510);
//     histoPCMEtaCorrectedSpecPbPbLHC11h0510->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510, widthLinesBoxes, kTRUE);
// 	graphPCMEtaCorrectedSpecSysPbPbLHC11h0510->Draw("E2same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
	histoPCMEtaCorrectedSpecPbPbLHC11h0010->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010, widthLinesBoxes, kTRUE);
	graphPCMEtaCorrectedSpecSysPbPbLHC11h0010->Draw("E2same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
//     histoPCMEtaCorrectedSpecPbPbLHC11h2040->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
// 	graphPCMEtaCorrectedSpecSysPbPbLHC11h2040->Draw("E2same");

	DrawGammaSetMarker(histoPCMEtaCorrectedSpecPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
	histoPCMEtaCorrectedSpecPbPbLHC11h2050->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMEtaCorrectedSpecSysPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
	graphPCMEtaCorrectedSpecSysPbPbLHC11h2050->Draw("E2same");

	
	TLatex *labelSpectraEtaLabelPbPb = new TLatex(0.6,0.89,"#eta #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraEtaLabelPbPb, 0.05  ,4);
	labelSpectraEtaLabelPbPb->Draw();

	TLegend* legendEtaSpectra = new TLegend(0.2,0.12,0.8,0.3);
	legendEtaSpectra->SetFillColor(0);
	legendEtaSpectra->SetLineColor(0);
	legendEtaSpectra->SetTextSize(0.03);
	legendEtaSpectra->SetMargin(0.2);
// 	legendEtaSpectra->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h0005,"PCM measured","pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h0010
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb0005.Data(),"");
// 	legendEtaSpectra->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h0510,"PCM measured","pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h0010
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb0510.Data(),"");
	legendEtaSpectra->AddEntry(graphPCMEtaCorrectedSpecSysPbPbLHC11h0010,Form("PCM, %s ",collisionSystemPbPb0010.Data()),"pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h0010
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb0010.Data(),"");
	legendEtaSpectra->AddEntry(graphPCMEtaCorrectedSpecSysPbPbLHC11h2050,Form("PCM, %s ",collisionSystemPbPb2050.Data()),"pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h2050
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb2050.Data(),"");
//  legendEtaSpectra->AddEntry(histoPCMEtaCorrectedSpecPbPbLHC11h2040,"measured","pf"); //graphPCMEtaCorrectedSpecSysPbPbLHC11h2040
//  legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPbPb2040.Data(),"");
	
//	legendEtaSpectra->AddEntry(histoEtaFromCocktailpp2760GeV,"m_{T} scaled","l");
// 	legendEtaSpectra->AddEntry(graphCorrectedYieldSysEta2760GeV,"measured","pf");
// 	legendEtaSpectra->AddEntry((TObject*)0, collisionSystemPP2760GeV.Data(),"");
//  legendEtaSpectra->AddEntry(graphInvSectionCombSysPi02760GeVPlot,collisionSystemPP.Data(),"pf");
//  legendEtaSpectra->AddEntry(fitInvCrossSectionPi0Comb2760GeV,"Tsallis Fit","l");
//  legendEtaSpectra->AddEntry(fitInvCrossSectionPi0Comb2760GeVPow,"Powerlaw Fit","l");
	legendEtaSpectra->Draw();

	canvasEtaSpectra->Update();
	canvasEtaSpectra->Print(Form("%s/SpectraEtaPbPbLHC11h.%s",outputDir.Data(),suffix.Data()));
	

	TCanvas* canvasPi0Spectra = new TCanvas("canvasPi0Spectra","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasPi0Spectra,  0.13, 0.01, 0.015, 0.08);
	canvasPi0Spectra->SetLogy();
	canvasPi0Spectra->SetLogx();
	TH2F * histo2DPi0Spectra = new TH2F("histo2DPi0Spectra","histo2DPi0Spectra",1000,0.1,15.,1000,1e-7,2e3 );
	SetStyleHistoTH2ForGraphs(histo2DPi0Spectra, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-2}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
	histo2DPi0Spectra->DrawCopy(); 

	
// 	graphCorrectedYieldSysPi02760GeV->RemovePoint(17);
// 	DrawGammaSetMarkerTGraphAsym(graphCorrectedYieldSysPi02760GeV, markerStylePP,markerSizePP, kBlack , kBlack, widthLinesBoxes, kTRUE);//, colorCombPbPb1020-5);
// 	graphCorrectedYieldSysPi02760GeV->Draw("E2same");   

// 	histoCorrectedYieldPi02760GeV->SetBinContent(19,0);
// 	DrawGammaSetMarker(histoCorrectedYieldPi02760GeV, markerStylePP,markerSizePP, kBlack , kBlack);
// 	histoCorrectedYieldPi02760GeV->Draw("p,same,e1");
	
	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005);
//  histoPCMPi0CorrectedSpecPbPbLHC11h0005->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h0005, markerStylePbPb0005,markerSizePbPb0005, colorCombPbPb0005 , colorCombPbPb0005, widthLinesBoxes, kTRUE);
// 	graphPCMPi0CorrectedSpecSysPbPbLHC11h0005->Draw("E2same");
	
	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510);
//     histoPCMPi0CorrectedSpecPbPbLHC11h0510->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h0510, markerStylePbPb0510,markerSizePbPb0510, colorCombPbPb0510 , colorCombPbPb0510, widthLinesBoxes, kTRUE);
// 	graphPCMPi0CorrectedSpecSysPbPbLHC11h0510->Draw("E2same");

	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010);
	histoPCMPi0CorrectedSpecPbPbLHC11h0010->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h0010, markerStylePbPb0010,markerSizePbPb0010, colorCombPbPb0010 , colorCombPbPb0010, widthLinesBoxes, kTRUE);
	graphPCMPi0CorrectedSpecSysPbPbLHC11h0010->Draw("E2same");

	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
//     histoPCMPi0CorrectedSpecPbPbLHC11h2040->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h2040, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
// 	graphPCMPi0CorrectedSpecSysPbPbLHC11h2040->Draw("E2same");

	DrawGammaSetMarker(histoPCMPi0CorrectedSpecPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040);
	histoPCMPi0CorrectedSpecPbPbLHC11h2050->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphAsym(graphPCMPi0CorrectedSpecSysPbPbLHC11h2050, markerStylePbPb2040,markerSizePbPb2040, colorCombPbPb2040 , colorCombPbPb2040, widthLinesBoxes, kTRUE);
	graphPCMPi0CorrectedSpecSysPbPbLHC11h2050->Draw("E2same");

	TLatex *labelSpectraPi0LabelPbPb = new TLatex(0.6,0.89,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-}");
	SetStyleTLatex( labelSpectraPi0LabelPbPb, 0.04  ,4);
	labelSpectraPi0LabelPbPb->Draw();

	TLegend* legendPi0Spectra = new TLegend(0.2,0.12,0.8,0.3);
	legendPi0Spectra->SetFillColor(0);
	legendPi0Spectra->SetLineColor(0);
	legendPi0Spectra->SetTextSize(0.03);
	legendPi0Spectra->SetMargin(0.2);
// 	legendPi0Spectra->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h0005,"PCM measured","pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h0005
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb0005.Data(),"");
// 	legendPi0Spectra->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h0510,"PCM measured","pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h0510
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb0510.Data(),"");
	legendPi0Spectra->AddEntry(graphPCMPi0CorrectedSpecSysPbPbLHC11h0010,Form("PCM, %s ",collisionSystemPbPb0010.Data()),"pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h0010
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb0010.Data(),"");
// 	legendPi0Spectra->AddEntry(histoPCMPi0CorrectedSpecPbPbLHC11h2040,"PCM measured","pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h2040
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb2040.Data(),"");
	legendPi0Spectra->AddEntry(graphPCMPi0CorrectedSpecSysPbPbLHC11h2050,Form("PCM, %s ",collisionSystemPbPb2050.Data()),"pf"); //graphPCMPi0CorrectedSpecSysPbPbLHC11h2050
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPbPb2050.Data(),"");
// 	legendPi0Spectra->AddEntry(graphCorrectedYieldSysPi02760GeV,"PCM","pf");
// 	legendPi0Spectra->AddEntry((TObject*)0, collisionSystemPP2760GeV.Data(),"");
	legendPi0Spectra->Draw();
	
	canvasPi0Spectra->Update();
	canvasPi0Spectra->Print(Form("%s/SpectraPi0PbPbLHC11h.%s",outputDir.Data(),suffix.Data()));

	
	
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
	
    TCanvas* canvasRAA = new TCanvas("canvasRAA","",200,10,1200,1100);  // gives the page size
    DrawGammaCanvasSettings( canvasRAA,  0.13, 0.05, 0.02, 0.08);
	TH2F * histo2DRAA = new TH2F("histo2DRAA","histo2DRAA",1000,0.,15.,1000,0.,1.);
	SetStyleHistoTH2ForGraphs(histo2DRAA, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.03,0.04, 0.03,0.04, 0.83,1.4);
	histo2DRAA->DrawCopy(); 

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0005, 20 ,1, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0005,20 ,1, kRed, kRed, 2, kTRUE, kRed-9);
// 	graphPi0RAASys0005->Draw("2same");
// 	graphPi0RAA0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0010, 20,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0010,20,2, kRed-9, kRed-9, 2, kTRUE, kRed-9);
	graphPi0RAASys0010->SetFillStyle(0);
	graphPi0RAASys0010->Draw("2same");
	graphPi0RAA0010->Draw("p,same"); 
	
	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 20,1, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
// 	graphPi0RAASys2040->Draw("2same");
// 	graphPi0RAA2040->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2050, 20,2, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2050,20,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
	graphPi0RAASys2050->SetFillStyle(0);
	graphPi0RAASys2050->Draw("2,same");
	graphPi0RAA2050->Draw("p,same");

	
	TLegend* legendPi0RAA = new TLegend(0.45,0.77,0.9,0.95);  //0.16,0.05,0.73,0.2);
	legendPi0RAA->SetFillColor(0);
	legendPi0RAA->SetFillStyle(0);
	legendPi0RAA->SetLineColor(0);
	legendPi0RAA->SetTextSize(0.035);
	legendPi0RAA->SetMargin(0.2);
// 	legendPi0RAA->AddEntry(graphPi0RAASys0005,"PCM #pi^{0} R_{AA} (2010) - 0-5%");
	legendPi0RAA->AddEntry(graphPi0RAA0010,"PCM #pi^{0} R_{AA} (2010) - 0-10%", "pl");
// 	legendPi0RAA->AddEntry(graphPi0RAASys2040,"PCM #pi^{0} R_{AA} (2010) - 20-40%");
	legendPi0RAA->AddEntry(graphPi0RAA2050,"PCM #pi^{0} R_{AA} (2010) - 20-50%", "pl");
	legendPi0RAA->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0PbPbLHC11h.%s",outputDir.Data(),suffix.Data()));

	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");

	DrawGammaSetMarkerTGraphAsym(graphRAAPublished0005, 20,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished0005,20,2, kRed-9, kRed-9, 2, kTRUE, kRed-9);
	graphRAASysPublished0005->SetFillStyle(0);
// 	graphRAASysPublished0005->Draw("2same");
// 	graphRAAPublished0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphRAAPublished0010, 20,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished0010,20,2, kRed-9, kRed-9, 2, kTRUE, kRed-9);
	graphRAASysPublished0010->SetFillStyle(0);
	graphRAASysPublished0010->Draw("2same");
	graphRAAPublished0010->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphRAAPublished2040, 20,2, kGreen+2, kGreen+2, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished2040,20,2, kGreen-6, kGreen-6, 2, kTRUE, kGreen-6);
	graphRAASysPublished2040->SetFillStyle(0);
// 	graphRAASysPublished2040->Draw("2same");
// 	graphRAAPublished2040->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0005, 21,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0005,21,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
	graphPi0RAASys0005->SetFillStyle(0);
// 	graphPi0RAASys0005->Draw("2same");
// 	graphPi0RAA0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0010, 21,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0010,21,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
	graphPi0RAASys0010->SetFillStyle(0);
	graphPi0RAASys0010->Draw("2same");
	graphPi0RAA0010->Draw("p,same"); 
	
	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2050, 20,1, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2050,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
	graphPi0RAASys2050->SetFillStyle(3008);
// 	graphPi0RAASys2050->Draw("2same");
// 	graphPi0RAA2050->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 21,2, kBlue,kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,21,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
	graphPi0RAASys2040->SetFillStyle(0);
// 	graphPi0RAASys2040->Draw("2same");
// 	graphPi0RAA2040->Draw("p,same");
	
	TLegend* legendRAAPi0AllPCM0010 = new TLegend(0.2,0.75,0.7,0.9);  //0.16,0.05,0.73,0.2);
	legendRAAPi0AllPCM0010->SetFillColor(0);
	legendRAAPi0AllPCM0010->SetFillStyle(0);
	legendRAAPi0AllPCM0010->SetLineColor(0);
	legendRAAPi0AllPCM0010->SetTextSize(0.035);
	legendRAAPi0AllPCM0010->SetMargin(0.2);
	legendRAAPi0AllPCM0010->SetHeader("0-10% PbPb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
// 	legendRAAPi0AllPCM0010->AddEntry(graphRAAPublished0005,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108) - 0-5%");
	legendRAAPi0AllPCM0010->AddEntry(graphRAAPublished0010,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108)");
// 	legendRAAPi0AllPCM0010->AddEntry(graphRAAPublished2040,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108) - 20-40%");
// 	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAA0005,"PCM #pi^{0} R_{AA} (2010) - 0-5%", "pl");
	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAA0010,"PCM #pi^{0} R_{AA} (2010)","pl");
// 	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAA2040,"PCM #pi^{0} R_{AA} (2010) - 20-40%", "pl");
// 	legendRAAPi0AllPCM0010->AddEntry(graphPi0RAASys2050,"PCM #pi^{0} R_{AA} (2010) - 20-50%");
	legendRAAPi0AllPCM0010->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0WithPublished0010.%s",outputDir.Data(),suffix.Data()));
	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	
	DrawGammaSetMarkerTGraphAsym(graphRAAPublished2040, 20,2, kGreen+2, kGreen+2, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished2040,20,2, kGreen-6, kGreen-6, 2, kTRUE, kGreen-6);
	graphRAASysPublished2040->SetFillStyle(0);
	graphRAASysPublished2040->Draw("2same");
	graphRAAPublished2040->Draw("p,same"); 


	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 21,2, kBlue,kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,21,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
	graphPi0RAASys2040->SetFillStyle(0);
	graphPi0RAASys2040->Draw("2same");
	graphPi0RAA2040->Draw("p,same");
	
	TLegend* legendRAAPi0AllPCM2040 = new TLegend(0.2,0.75,0.7,0.9);  //0.16,0.05,0.73,0.2);
	legendRAAPi0AllPCM2040->SetFillColor(0);
	legendRAAPi0AllPCM2040->SetFillStyle(0);
	legendRAAPi0AllPCM2040->SetLineColor(0);
	legendRAAPi0AllPCM2040->SetTextSize(0.035);
	legendRAAPi0AllPCM2040->SetMargin(0.2);
	legendRAAPi0AllPCM2040->SetHeader("20-40% PbPb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
	legendRAAPi0AllPCM2040->AddEntry(graphRAAPublished2040,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108)");
	legendRAAPi0AllPCM2040->AddEntry(graphPi0RAA2040,"PCM #pi^{0} R_{AA} (2010)", "pl");
	legendRAAPi0AllPCM2040->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0WithPublished2040.%s",outputDir.Data(),suffix.Data()));
	
	

//==========================================================================================
	
    TCanvas* canvasRatiowithPub = new TCanvas("canvasRatiowithPub","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasRatiowithPub, 0.09, 0.01, 0.015, 0.115);
    
    TH2D *histo2DCompCombinedRatioPub = new TH2D("histo2DCompCombinedRatioPub", "histo2DCompCombinedRatioPub", 20,0.1,12.01,1000.,-0.4,2.);
    SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatioPub, "#it{p}_{T} (GeV/#it{c})","#pi^{0} corrected yield #frac{2010}{2011}", 0.035,0.04, 0.035,0.04, 1.,1.);
    histo2DCompCombinedRatioPub->GetYaxis()->SetRangeUser(0.,2.5);
    histo2DCompCombinedRatioPub->Draw("copy");

    TGraphErrors* grapha    = NULL;
    TGraphErrors* graphb    = NULL;
    TGraphErrors* graphc    = NULL;
    TGraphErrors* graphd     = NULL;
    TGraphErrors* graph1    = NULL;
    TGraphErrors* graph2     = NULL;
    TGraphErrors* graph3    = NULL;
    TGraphErrors* graph4     = NULL;

        
    TGraphErrors* RatioPi0YieldsPublishedTo2011_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMPi0CorrectedSpecPbPbLHC10h0010, graphPCMPi0CorrectedSpecSysPbPbLHC10h0010, 
                                                                                                     histoPCMPi0CorrectedSpecPbPbLHC11h0010, graphPCMPi0CorrectedSpecSysPbPbLHC11h0010, 
                                                                                                     kTRUE,  kTRUE, 
                                                                                                     &graphc, &graphd, 
                                                                                                     &grapha, &graphb ) ;
//     RatioPi0YieldsPublishedTo2011_0010->Print();

    TGraphErrors* RatioPi0YieldsPublishedTo2011_2040 = CalculateRatioBetweenSpectraWithDifferentBinning(histoPCMPi0CorrectedSpecPbPbLHC10h2040, graphPCMPi0CorrectedSpecSysPbPbLHC10h2040,
                                                                                                     histoPCMPi0CorrectedSpecPbPbLHC11h2040, graphPCMPi0CorrectedSpecSysPbPbLHC11h2040,
                                                                                                     kTRUE,  kTRUE, 
                                                                                                     &graph3, &graph4, 
                                                                                                     &graph1, &graph2 ) ;
//     RatioPi0YieldsPublishedTo2011_2040->Print();
    
    DrawGammaSetMarkerTGraphErr(RatioPi0YieldsPublishedTo2011_0010, 20., 2., kRed, kRed);
    DrawGammaSetMarkerTGraphErr(RatioPi0YieldsPublishedTo2011_2040, 20., 2., kAzure+1, kAzure+1);

    RatioPi0YieldsPublishedTo2011_0010->Draw("same,p");
    RatioPi0YieldsPublishedTo2011_2040->Draw("same,p");
    
    TLegend* legendRatio1 = new TLegend(0.15,0.13,0.32,0.32);
    legendRatio1->SetTextSize(0.035);         
    legendRatio1->SetFillColor(0);
    legendRatio1->SetFillStyle(0);
    legendRatio1->SetBorderSize(0);
    legendRatio1->AddEntry(RatioPi0YieldsPublishedTo2011_0010," 0-10%","pl");
    legendRatio1->AddEntry(RatioPi0YieldsPublishedTo2011_2040," 20-40%","pl");
    legendRatio1->Draw();

    canvasRatiowithPub->SaveAs(Form("%s/Pi0yields_RatioWithPublished.%s",outputDir.Data(),suffix.Data()));

	
	TCanvas* canvasRatioRaa = new TCanvas("canvasRatioRaa","",200,10,1350,900);  // gives the page size
	DrawGammaCanvasSettings( canvasRatioRaa, 0.09, 0.01, 0.015, 0.115);
	
	TH2F * histo2DCompCombinedRatio2;
	histo2DCompCombinedRatio2 = new TH2F("histo2DCompCombinedRatio2","histo2DCompCombinedRatio2",1000,0.,12.,1000,0.2,2);
	SetStyleHistoTH2ForGraphs(histo2DCompCombinedRatio2, "#it{p}_{T} (GeV/#it{c})","#pi^{0} #it{R}_{AA} #frac{2010}{2011}", 0.03,0.04, 0.03,0.04, 0.83,1);// 512, 505); //#frac{#frac{1}{N_{evt}^{AA}}#frac{d#it{N}^{AA}}{d#it{p}_{T} d#it{y}}}{ #frac{1}{N_{evt}^{pp}}#frac{d#it{N}^{pp}}{d#it{p}_{T} d#it{y}}}
	histo2DCompCombinedRatio2->DrawCopy(); 

	
	TGraphErrors* graphRaaStatPi0Pt2010_0010 	= NULL;
	TGraphErrors* graphRaaSysPi0Pt2010_0010 	= NULL;
	TGraphErrors* graphRaaStatPi0Pt2011_0010	= NULL;
	TGraphErrors* graphRaaSysPi0Pt2011_0010 	= NULL;
	TGraphErrors* graphRaaStatPi0Pt2010_2040 	= NULL;
	TGraphErrors* graphRaaSysPi0Pt2010_2040 	= NULL;
	TGraphErrors* graphRaaStatPi0Pt2011_2040 	= NULL;
	TGraphErrors* graphRaaSysPi0Pt2011_2040 	= NULL;

	TGraphErrors* RatioPi0RaaPublishedTo2011_0010 = CalculateRatioBetweenSpectraWithDifferentBinning(graphRAAPublished0010, graphRAASysPublished0010, 
																									 graphPi0RAA0010, graphPi0RAASys0010, 
																									 kTRUE,  kTRUE, 
																									 &graphRaaStatPi0Pt2011_0010, &graphRaaSysPi0Pt2011_0010, 
																									 &graphRaaStatPi0Pt2010_0010, &graphRaaSysPi0Pt2010_0010 ) ;
// 	RatioPi0RaaPublishedTo2011_0010->Print();

	TGraphErrors* RatioPi0RaaPublishedTo2011_2040 = CalculateRatioBetweenSpectraWithDifferentBinning(graphRAAPublished2040, graphRAASysPublished2040,
																									 graphPi0RAA2040, graphPi0RAASys2040,
																									 kTRUE,  kTRUE, 
																									 &graphRaaStatPi0Pt2011_2040, &graphRaaSysPi0Pt2011_2040, 
																									 &graphRaaStatPi0Pt2010_2040, &graphRaaSysPi0Pt2010_2040 ) ;
// 	RatioPi0RaaPublishedTo2011_2040->Print();
	
	DrawGammaSetMarkerTGraphErr(RatioPi0RaaPublishedTo2011_0010, 20., 1.5, kRed, kRed);
	DrawGammaSetMarkerTGraphErr(RatioPi0RaaPublishedTo2011_2040, 20., 1.5, kAzure+1, kAzure+1);

	RatioPi0RaaPublishedTo2011_0010->Draw("same,p");
	RatioPi0RaaPublishedTo2011_2040->Draw("same,p");
	

	TLegend* legendRatio = new TLegend(0.25,0.13,0.32,0.25);
	legendRatio->SetTextSize(0.035);			
	legendRatio->SetFillColor(0);
	legendRatio->SetFillStyle(0);
	legendRatio->SetBorderSize(0);
	legendRatio->AddEntry(RatioPi0RaaPublishedTo2011_0010," 0-10%","pl");
	legendRatio->AddEntry(RatioPi0RaaPublishedTo2011_2040," 20-40%","pl");
	legendRatio->Draw();

	canvasRatioRaa->SaveAs(Form("%s/Pi0RAA_RatioWithPublished.%s",outputDir.Data(),suffix.Data()));

//==========================================================================================

	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");

	//Comparison with charged kaons and pions:
	//PbPb276.fullpT.RATIOS.20140329.root RAA_Pion_08052014.root RAA_Kaon_08052014.root

	for(Int_t j = 0; j<6; j++){
		graphChargedPionRAA0010->RemovePoint(0);
		graphChargedPionRAASys0010->RemovePoint(0);
	}

	for(Int_t i = 0; i<3; i++){
		graphChargedPionRAA0010->RemovePoint(graphChargedPionRAA0010->GetN()-1);
		graphChargedPionRAASys0010->RemovePoint(graphChargedPionRAASys0010->GetN()-1);
	}

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0005, 20,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0005,20,2, kRed, kRed, 2, kTRUE, kRed-9);
	graphPi0RAASys0005->SetFillStyle(0);
// 	graphPi0RAASys0005->Draw("2same");
// 	graphPi0RAA0005->Draw("p,same"); 

	graphChargedPionRAA0005->RemovePoint(0);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0005, 21,2, kOrange+1, kOrange+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0005,21,2, kOrange-2, kOrange-2, 2, kTRUE, kOrange-2);
// 	graphChargedPionRAASys0005->Draw("2same");
// 	graphChargedPionRAA0005->Draw("p,same"); 

	
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0005, 20,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0005,20,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
// 	graphChargedKaonRAASys0005->Draw("2same");
// 	graphChargedKaonRAA0005->Draw("p,same"); 

	graphChargedPionRAA0010->RemovePoint(0);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0010, 21,2, kOrange+1, kOrange+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0010,21,2, kOrange-2, kOrange-2, 2, kTRUE, kOrange-2);
	graphChargedPionRAASys0010->Draw("2same");
	graphChargedPionRAA0010->Draw("p,same"); 

	
	DrawGammaSetMarkerTGraphAsym(graphRAAPublished0010, 20,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished0010,20,2, kRed-9, kRed-9, 2, kTRUE, kRed-9);
	graphRAASysPublished0010->SetFillStyle(0);
	graphRAASysPublished0010->Draw("2same");
	graphRAAPublished0010->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0010, 20,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0010,20,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
	graphPi0RAASys0010->SetFillStyle(0);
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
	
	TLegend* legendChargedRAAcomp1 = new TLegend(0.2,0.75,0.7,0.9);  //0.16,0.05,0.73,0.2);
	legendChargedRAAcomp1->SetFillColor(0);
	legendChargedRAAcomp1->SetFillStyle(0);
	legendChargedRAAcomp1->SetLineColor(0);
	legendChargedRAAcomp1->SetTextSize(0.035);
	legendChargedRAAcomp1->SetMargin(0.2);
	legendRAAPi0AllPCM2040->SetHeader("0-10% PbPb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
// 	legendChargedRAAcomp1->AddEntry(graphChargedKaonRAASys0005,"K^{#pm} R_{AA} - 0-5%");
// 	legendChargedRAAcomp1->AddEntry(graphChargedPionRAA0005,"#pi^{#pm} R_{AA} - 0-5%");
// 	legendChargedRAAcomp1->AddEntry(graphPi0RAA0005,"PCM measured #pi^{0} R_{AA} - 0-5%", "pl");
	legendChargedRAAcomp1->AddEntry(graphChargedPionRAA0010,"#pi^{#pm} R_{AA} (JIRA PWGLF-258)");
	legendChargedRAAcomp1->AddEntry(graphRAAPublished0010,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108)", "pl");
	legendChargedRAAcomp1->AddEntry(graphPi0RAA0010,"PCM #pi^{0} R_{AA} (2010)", "pl");
// 	legendChargedRAAcomp->AddEntry(graphChargedKaonRAASys2040,"K^{#pm} R_{AA} - 20-40%");
// 	legendChargedRAAcomp->AddEntry(graphChargedPionRAASys2040,"#pi^{#pm} R_{AA} - 20-40%");
// 	legendChargedRAAcomp->AddEntry(graphPi0RAASys2040,"PCM measured #pi^{0} R_{AA} - 20-40%");
	legendChargedRAAcomp1->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0WithChargedPions0010.%s",outputDir.Data(),suffix.Data()));

	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	
// 	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0005, 20,1, kRed, kRed, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0005,20,1, kRed, kRed, 2, kTRUE, kRed-9);
// 	graphChargedKaonRAASys0005->Draw("2same");
// 	graphChargedKaonRAA0005->Draw("p,same"); 
// 
// 	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0005, 20,1, kOrange, kOrange, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0005,20,1, kOrange, kOrange, 2, kTRUE, kOrange-2);
// 	graphChargedPionRAASys0005->Draw("2same");
// 	graphChargedPionRAA0005->Draw("p,same"); 
// 
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0005, 20,1, kMagenta+1, kMagenta+1, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0005,20,1, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
// 	graphPi0RAASys0005->SetFillStyle(3001);
// 	graphPi0RAASys0005->Draw("2same");
// 	graphPi0RAA0005->Draw("p,same"); 

	for(Int_t j = 0; j<6; j++){
		graphChargedPionRAA2040->RemovePoint(0);
		graphChargedPionRAASys2040->RemovePoint(0);
	}

	for(Int_t i = 0; i<3; i++){
		graphChargedPionRAA2040->RemovePoint(graphChargedPionRAA2040->GetN()-1);
		graphChargedPionRAASys2040->RemovePoint(graphChargedPionRAASys2040->GetN()-1);
	}

	
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040, 20,1, kGreen+2, kGreen+2, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040,20,1, kGreen+2, kGreen+2, 2, kTRUE, kGreen-6);
// 	graphChargedKaonRAASys2040->Draw("2same");
// 	graphChargedKaonRAA2040->Draw("p,same"); 
		
	graphChargedPionRAA2040->RemovePoint(0);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA2040, 21,2, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys2040,2,2, kBlue-10, kBlue-10, 2, kTRUE, kBlue-10);
// 	graphChargedPionRAASys2040->SetFillStyle(0);
	graphChargedPionRAASys2040->Draw("2same");
	graphChargedPionRAA2040->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphRAAPublished2040, 20,2, kGreen+3, kGreen+3, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphRAASysPublished2040,20,2, kGreen+2, kGreen+2, 2, kTRUE, kGreen+2);
	graphRAASysPublished2040->SetFillStyle(0);
	graphRAASysPublished2040->Draw("2same");
	graphRAAPublished2040->Draw("p,same"); 

	
	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, 20,2, kAzure+1, kAzure+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040,20,2, kAzure+1, kAzure+1, 2, kTRUE, kAzure-9);
	graphPi0RAASys2040->SetFillStyle(0);
	graphPi0RAASys2040->Draw("2same");
	graphPi0RAA2040->Draw("p,same");
	

	TLegend* legendChargedRAAcomp = new TLegend(0.2,0.75,0.7,0.9);  //0.16,0.05,0.73,0.2);
	legendChargedRAAcomp->SetFillColor(0);
	legendChargedRAAcomp->SetFillStyle(0);
	legendChargedRAAcomp->SetLineColor(0);
	legendChargedRAAcomp->SetTextSize(0.035);
	legendChargedRAAcomp->SetMargin(0.2);
	legendRAAPi0AllPCM2040->SetHeader("20-40% PbPb, #sqrt{#it{s}_{_{NN}}} = 2.76 TeV");
// 	legendChargedRAAcomp->AddEntry(graphChargedKaonRAASys0005,"K^{#pm} R_{AA} - 0-5%");
// 	legendChargedRAAcomp->AddEntry(graphChargedPionRAASys0005,"#pi^{#pm} R_{AA} - 0-5%");
// 	legendChargedRAAcomp->AddEntry(graphPi0RAASys0005,"PCM measured #pi^{0} R_{AA} - 0-5%");
// 	legendChargedRAAcomp->AddEntry(graphChargedKaonRAASys2040,"K^{#pm} R_{AA} - 20-40%");
	legendChargedRAAcomp->AddEntry(graphChargedPionRAA2040,"#pi^{#pm} R_{AA} (JIRA PWGLF-258)");
	legendChargedRAAcomp->AddEntry(graphRAAPublished2040,"PCM-PHOS #pi^{0} R_{AA} (EPJC 74 (2014) 3108)", "pl");
	legendChargedRAAcomp->AddEntry(graphPi0RAA2040,"PCM #pi^{0} R_{AA} (2010)", "pl");
	legendChargedRAAcomp->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAPi0WithChargedPions2040.%s",outputDir.Data(),suffix.Data()));

	
	
	// comparison Pi0 RAA from 2011 as in published paper plot
	TCanvas* canvasRAAasPaper = new TCanvas("canvasRAAasPaper","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasRAAasPaper,  0.1, 0.01, 0.015, 0.1);

	Int_t textSizeLabelsPixelRAA = 50;
	Double_t marginRAA = 0.14*1200;
	Double_t textsizeLabelsRAA = 0;
	Double_t textsizeFacRAA = 0;

	if (canvasRAAasPaper->XtoPixel(canvasRAAasPaper->GetX2()) < canvasRAAasPaper->YtoPixel(canvasRAAasPaper->GetY1())){
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAasPaper->XtoPixel(canvasRAAasPaper->GetX2()) ;
		textsizeFacRAA = (Double_t)1./canvasRAAasPaper->XtoPixel(canvasRAAasPaper->GetX2()) ;
	} else {
		textsizeLabelsRAA = (Double_t)textSizeLabelsPixelRAA/canvasRAAasPaper->YtoPixel(canvasRAAasPaper->GetY1());
		textsizeFacRAA = (Double_t)1./canvasRAAasPaper->YtoPixel(canvasRAAasPaper->GetY1());
	}

	TH2F * histo2DRAAAll3Up = new TH2F("histo2DRAAAll3Up","histo2DRAAAll3Up",1000,0.,20.5,1000,0.,2.1);
	SetStyleHistoTH2ForGraphs(histo2DRAAAll3Up, "#it{p}_{T} (GeV/#it{c})","#it{R}_{AA}", 0.85*textsizeLabelsRAA,textsizeLabelsRAA,0.85*textsizeLabelsRAA,textsizeLabelsRAA, 0.95,1., 512, 505); 
	histo2DRAAAll3Up->GetXaxis()->SetLabelFont(42);
	histo2DRAAAll3Up->GetYaxis()->SetLabelFont(42);
	histo2DRAAAll3Up->GetYaxis()->SetLabelOffset(0.005);
	histo2DRAAAll3Up->DrawCopy("");

	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0010, markerStyleCommmonSpectrum0010,markerSizeCommonSpectrum0010, colorComb0010 , colorComb0010, widthLinesBoxes, kTRUE, kGray+1);
	graphPi0RAASys0010->Draw("E2same");
	graphPi0RAASys0010->SetFillStyle(0);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0010, markerStyleCommmonSpectrum0010,markerSizeCommonSpectrum0010, colorComb0010 , colorComb0010);
	graphPi0RAA0010->Draw("p,same,e1");
	DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_0010, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_0010, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_0010, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphWA98_17_3GeVRAA_0013, markerStyleWA98,markerSizeWA98, kGray+2 , kGray+2);

	graphPHENIX200GeVRAA_0010->Draw("p,same,e1");	
	graphPHENIX39GeVRAA_0010->Draw("p,same,e1");	
	graphPHENIX62GeVRAA_0010->Draw("p,same,e1");	
	graphWA98_17_3GeVRAA_0013->Draw("p,same,e1");	

	histo2DRAAAll3Up->Draw("axis,same");

	TLatex *labelRAAALICEPbPb0010 = new TLatex(0.35,0.93,"#pi^{0} ALICE    0-10% Pb-Pb (LHC11h)");
	SetStyleTLatex( labelRAAALICEPbPb0010, 0.85*textsizeLabelsRAA,4);
	labelRAAALICEPbPb0010->Draw();
	TLegend* legendRAASinglePbPb0010 = new TLegend(0.35,0.88,0.65,0.92);
	legendRAASinglePbPb0010->SetFillColor(0);
	legendRAASinglePbPb0010->SetLineColor(0);
	legendRAASinglePbPb0010->SetNColumns(1);
	legendRAASinglePbPb0010->SetTextFont(42);
	legendRAASinglePbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
// 	legendRAASinglePbPb0010->AddEntry(graphPi0RAASys0010,"#sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
	legendRAASinglePbPb0010->AddEntry(graphPi0RAASys0010,"0-10% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
	legendRAASinglePbPb0010->Draw();

	TLatex *labelRAAPHENIXPbPb0010 = new TLatex(0.35,0.83,"#pi^{0} PHENIX 0-10% Au-Au");
	SetStyleTLatex( labelRAAPHENIXPbPb0010, 0.85*textsizeLabelsRAA,4);
	labelRAAPHENIXPbPb0010->Draw();

	TLegend* legendRAARHICPbPb0010 = new TLegend(0.35,0.73,0.95,0.82);
	legendRAARHICPbPb0010->SetFillColor(0);
	legendRAARHICPbPb0010->SetLineColor(0);
	legendRAARHICPbPb0010->SetNColumns(2);
	legendRAARHICPbPb0010->SetTextFont(42);
	legendRAARHICPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAARHICPbPb0010->AddEntry(graphPHENIX200GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
	legendRAARHICPbPb0010->AddEntry(graphPHENIX62GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
	legendRAARHICPbPb0010->AddEntry(graphPHENIX39GeVRAA_0010,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
	legendRAARHICPbPb0010->Draw();

	TLatex *labelRAAWA98PbPb0010 = new TLatex(0.35,0.68,"#pi^{0} WA98     0-13% Pb-Pb");
	SetStyleTLatex( labelRAAWA98PbPb0010, 0.85*textsizeLabelsRAA,4);
	labelRAAWA98PbPb0010->Draw();

	TLegend* legendRAASPSPbPb0010 = new TLegend(0.35,0.63,0.95,0.67);
	legendRAASPSPbPb0010->SetFillColor(0);
	legendRAASPSPbPb0010->SetLineColor(0);
	legendRAASPSPbPb0010->SetNColumns(2);
	legendRAASPSPbPb0010->SetTextFont(42);
	legendRAASPSPbPb0010->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAASPSPbPb0010->AddEntry(graphWA98_17_3GeVRAA_0013,"#sqrt{#it{s}_{_{NN}}} = 17.3 GeV","p");
	legendRAASPSPbPb0010->Draw();

	TBox* boxErrorNorm0010_Single = CreateBoxConv(colorComb0005Box, 0.2, 1.-normErr0010 , 0.5, 1.+normErr0010);
	boxErrorNorm0010_Single->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	canvasRAAasPaper->Update();
	canvasRAAasPaper->Print(Form("%s/RAAPi0asPaper0010.%s",outputDir.Data(),suffix.Data()));
	
	
	canvasRAAasPaper->cd();
	histo2DRAAAll3Up->DrawCopy("");

	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040, widthLinesBoxes, kTRUE, kGray+1);
	graphPi0RAASys2040->Draw("E2same");
	graphPi0RAASys2040->SetFillStyle(0);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2040, markerStyleCommmonSpectrum2040,markerSizeCommonSpectrum2040, colorComb2040 , colorComb2040);
	graphPi0RAA2040->Draw("p,same,e1");
	
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2050, markerStyleCommmonSpectrum2040+1,markerSizeCommonSpectrum2040, kBlue , kBlue, widthLinesBoxes, kTRUE, kGray+1);
// 	graphPi0RAASys2050->Draw("E2same");
	graphPi0RAASys2050->SetFillStyle(0);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2050, markerStyleCommmonSpectrum2040+1,markerSizeCommonSpectrum2040, kBlue , kBlue);
// 	graphPi0RAA2050->Draw("p,same,e1");

	DrawGammaSetMarkerTGraphErr(graphPHENIX200GeVRAA_2040, markerStylePHENIX200GeV,markerSizePHENIX200GeV, kBlack , kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX62GeVRAA_2040, markerStylePHENIX62GeV,markerSizePHENIX62GeV, kBlack, kBlack);
	DrawGammaSetMarkerTGraphErr(graphPHENIX39GeVRAA_2040, markerStylePHENIX39GeV,markerSizePHENIX39GeV, kBlack , kBlack);

	graphPHENIX200GeVRAA_2040->Draw("p,same,e1");	
	graphPHENIX39GeVRAA_2040->Draw("p,same,e1");	
	graphPHENIX62GeVRAA_2040->Draw("p,same,e1");	

	histo2DRAAAll3Up->Draw("axis,same");

// 	TLatex *labelRAAALICEPbPb2040 = new TLatex(0.15,0.93,"#pi^{0} ALICE    20-40% Pb-Pb (LHC11h)");
	TLatex *labelRAAALICEPbPb2040 = new TLatex(0.5,0.93,"#pi^{0} ALICE    Pb-Pb (LHC11h)");
	SetStyleTLatex( labelRAAALICEPbPb2040, 0.85*textsizeLabelsRAA,4);
	labelRAAALICEPbPb2040->Draw();
	TLegend* legendRAASinglePbPb2040 = new TLegend(0.5,0.88,0.8,0.92);
	legendRAASinglePbPb2040->SetFillColor(0);
	legendRAASinglePbPb2040->SetLineColor(0);
	legendRAASinglePbPb2040->SetNColumns(1);
	legendRAASinglePbPb2040->SetTextFont(42);
	legendRAASinglePbPb2040->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAASinglePbPb2040->AddEntry(graphPi0RAASys2040,"20-40% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
// 	legendRAASinglePbPb2040->AddEntry(graphPi0RAASys2050,"20-50% #sqrt{#it{s}_{_{NN}}} = 2.76 TeV","pf");
	legendRAASinglePbPb2040->Draw();

	TLatex *labelRAAPHENIXPbPb2040 = new TLatex(0.5,0.83,"#pi^{0} PHENIX 20-40% Au-Au");
	SetStyleTLatex( labelRAAPHENIXPbPb2040, 0.85*textsizeLabelsRAA,4);
	labelRAAPHENIXPbPb2040->Draw();

	TLegend* legendRAARHICPbPb2040 = new TLegend(0.5,0.68,0.95,0.82);
	legendRAARHICPbPb2040->SetFillColor(0);
	legendRAARHICPbPb2040->SetLineColor(0);
// 	legendRAARHICPbPb2040->SetNColumns(2);
	legendRAARHICPbPb2040->SetTextFont(42);
	legendRAARHICPbPb2040->SetTextSize(0.85*textsizeLabelsRAA);
	legendRAARHICPbPb2040->AddEntry(graphPHENIX200GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 200 GeV","p");
	legendRAARHICPbPb2040->AddEntry(graphPHENIX62GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 62.4 GeV","p");
	legendRAARHICPbPb2040->AddEntry(graphPHENIX39GeVRAA_2040,"#sqrt{#it{s}_{_{NN}}} = 39 GeV","p");
	legendRAARHICPbPb2040->Draw();

	TBox* boxErrorNorm2040_Single = CreateBoxConv(colorComb2040Box, 0.2, 1.-normErr2040 , 0.5, 1.+normErr2040);
	boxErrorNorm2040_Single->Draw();
	TBox* boxErrorNorm2050_Single = CreateBoxConv(colorComb2050Box, 0.2, 1.-normErr2050 , 0.5, 1.+normErr2050);
// 	boxErrorNorm2050_Single->Draw();
	DrawGammaLines(0., 19.5 , 1, 1 ,1,kGray);
	
	canvasRAAasPaper->Update();
	canvasRAAasPaper->Print(Form("%s/RAAPi0asPaper2040.%s",outputDir.Data(),suffix.Data()));

	
	
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	histo2DRAA->GetXaxis()->SetRangeUser(0.,10.5);
	
	DrawGammaSetMarkerTGraphAsym(graphEtaRAA0005, 20,1, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys0005,20,1, kRed, kRed, 2, kTRUE, kRed-9);
// 	graphEtaRAASys0005->Draw("2same");
// 	graphEtaRAA0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA0010, 20,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys0010,20,2, kRed-9, kRed-9, 2, kTRUE, kRed-9);
	graphEtaRAASys0010->SetFillStyle(0);
	graphEtaRAASys0010->Draw("2same");
	graphEtaRAA0010->Draw("p,same"); 
	
// 	graphEtaRAASys2050->RemovePoint(0);
// 	graphEtaRAA2050->RemovePoint(0);

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA2050, 20,2, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2050,20,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
	graphEtaRAASys2050->SetFillStyle(0);
	graphEtaRAASys2050->Draw("2,same");
	graphEtaRAA2050->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA2040, 20,1, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2040,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
// 	graphEtaRAASys2040->Draw("2same");
// 	graphEtaRAA2040->Draw("p,same");

	
	TLegend* legendEtaRAA = new TLegend(0.45,0.77,0.9,0.95);  //0.16,0.05,0.73,0.2);
	legendEtaRAA->SetFillColor(0);
	legendEtaRAA->SetFillStyle(0);
	legendEtaRAA->SetLineColor(0);
	legendEtaRAA->SetTextSize(0.035);
	legendEtaRAA->SetMargin(0.2);
// 	legendEtaRAA->AddEntry(graphEtaRAASys0005,"PCM #eta R_{AA} (2011) - 0-5%");
	legendEtaRAA->AddEntry(graphEtaRAA0010,"PCM #eta R_{AA} (2011) - 0-10%","pl");
// 	legendEtaRAA->AddEntry(graphEtaRAASys2040,"PCM #eta R_{AA} (2011) - 20-40%");
	legendEtaRAA->AddEntry(graphEtaRAA2050,"PCM #eta R_{AA} (2011) - 20-50%","pl");
	legendEtaRAA->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAEtaPbPbLHC11h.%s",outputDir.Data(),suffix.Data()));
	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	
	DrawGammaSetMarkerTGraphAsym(graphEtaRAAPHENIX0010, 20,1, kOrange, kOrange, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAAPHENIX0010,20,1, kOrange, kOrange, 2, kTRUE, kOrange-2);
	graphEtaRAAPHENIX0010->Draw("2same");
	graphEtaRAAPHENIX0010->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphEtaRAAPHENIX2040, 20,2, kGreen+2, kGreen+2, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAAPHENIX2040,20,2, kGreen-6, kGreen-6, 2, kTRUE, kGreen-6);
	graphEtaRAAPHENIX2040->SetFillStyle(0);
	graphEtaRAAPHENIX2040->Draw("2same");
	graphEtaRAAPHENIX2040->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphEtaRAAPHENIX2060, 20,2, kGreen+2, kGreen+2, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAAPHENIX2060,20,2, kGreen-6, kGreen-6, 2, kTRUE, kGreen-6);
	graphEtaRAAPHENIX2060->SetFillStyle(0);
// 	graphEtaRAAPHENIX2060->Draw("2same");
// 	graphEtaRAAPHENIX2060->Draw("p,same"); 
	
	graphEtaRAAPHENIX0010->Draw("2same");
	graphEtaRAAPHENIX0010->Draw("p,same"); 
	graphEtaRAAPHENIX2040->Draw("2same");
	graphEtaRAAPHENIX2040->Draw("p,same"); 
	
	graphEtaRAASys0010->Draw("2same");
	graphEtaRAA0010->Draw("p,same"); 
	
	DrawGammaSetMarkerTGraphAsym(graphEtaRAA2040, 20,2, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2040,20,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
	graphEtaRAASys2040->SetFillStyle(0);
	graphEtaRAASys2040->Draw("2same");
	graphEtaRAA2040->Draw("p,same");

	TLegend* legendEtaRAAcomp1 = new TLegend(0.45,0.72,0.9,0.9);  //0.16,0.05,0.73,0.2);
	legendEtaRAAcomp1->SetFillColor(0);
	legendEtaRAAcomp1->SetFillStyle(0);
	legendEtaRAAcomp1->SetLineColor(0);
	legendEtaRAAcomp1->SetTextSize(0.035);
	legendEtaRAAcomp1->SetMargin(0.2);
	legendEtaRAAcomp1->AddEntry(graphEtaRAAPHENIX0010,"PHENIX #eta R_{AA} - 0-10%","pl");
	legendEtaRAAcomp1->AddEntry(graphEtaRAAPHENIX2040,"PHENIX #eta R_{AA} - 20-40%","pl");
	legendEtaRAAcomp1->AddEntry(graphEtaRAASys0010,"PCM #eta R_{AA} (2011) - 0-10%", "pl");
	legendEtaRAAcomp1->AddEntry(graphEtaRAASys2040,"PCM #eta R_{AA} (2011) - 20-40%", "pl");
	legendEtaRAAcomp1->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAEtaWithPHENIX.%s",outputDir.Data(),suffix.Data()));
	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	
	graphEtaRAAPHENIX0010->Draw("2same");
	graphEtaRAAPHENIX0010->Draw("p,same"); 
	graphEtaRAAPHENIX2060->Draw("2same");
	graphEtaRAAPHENIX2060->Draw("p,same"); 
	
	graphEtaRAASys0010->Draw("2same");
	graphEtaRAA0010->Draw("p,same"); 
	
	graphEtaRAASys2050->Draw("2same");
	graphEtaRAA2050->Draw("p,same");

	TLegend* legendEtaRAAcomp2 = new TLegend(0.45,0.72,0.9,0.9);  //0.16,0.05,0.73,0.2);
	legendEtaRAAcomp2->SetFillColor(0);
	legendEtaRAAcomp2->SetFillStyle(0);
	legendEtaRAAcomp2->SetLineColor(0);
	legendEtaRAAcomp2->SetTextSize(0.035);
	legendEtaRAAcomp2->SetMargin(0.2);
	legendEtaRAAcomp2->AddEntry(graphEtaRAAPHENIX0010,"PHENIX #eta R_{AA} - 0-10%","pl");
 	legendEtaRAAcomp2->AddEntry(graphEtaRAAPHENIX2060,"PHENIX #eta R_{AA} - 20-60%","pl");
	legendEtaRAAcomp2->AddEntry(graphEtaRAA0010,"PCM #eta R_{AA} (2011) - 0-10%", "pl");
 	legendEtaRAAcomp2->AddEntry(graphEtaRAA2050,"PCM #eta R_{AA} (2011) - 20-50%","pl");
	legendEtaRAAcomp2->Draw();

	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAEtaWithPHENIX2.%s",outputDir.Data(),suffix.Data()));

	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	
	for(Int_t j = 0; j<6; j++){
		graphChargedKaonRAA0010->RemovePoint(0);
	}

	for(Int_t i = 0; i<8; i++){
		graphChargedKaonRAA0010->RemovePoint(graphChargedKaonRAA0010->GetN()-1);
		graphChargedKaonRAASys0010->RemovePoint(graphChargedKaonRAASys0010->GetN()-1);
	}
	
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0005, 21,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0005,21,2, kRed, kRed, 2, kTRUE, kRed-9);
// 	graphChargedKaonRAASys0005->Draw("2same");
// 	graphChargedKaonRAA0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA0010, 21,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys0010,21,2, kRed, kRed, 2, kTRUE, kRed-9);
	graphChargedKaonRAASys0010->Draw("2same");
	graphChargedKaonRAA0010->Draw("p,same"); 
	
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA0005, 20,1, kOrange, kOrange, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys0005,20,1, kOrange, kOrange, 2, kTRUE, kOrange-2);
// 	graphChargedPionRAASys0005->Draw("2same");
// 	graphChargedPionRAA0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA0005, 20,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys0005,20,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
	graphEtaRAASys0005->SetFillStyle(0);
// 	graphEtaRAASys0005->Draw("2same");
// 	graphEtaRAA0005->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA0010, 20,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys0010,20,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
	graphEtaRAASys0010->SetFillStyle(0);
	graphEtaRAASys0010->Draw("2same");
	graphEtaRAA0010->Draw("p,same"); 
		
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
	
	TLegend* legendChargedRAAcomp3 = new TLegend(0.4,0.75,0.9,0.9);  //0.16,0.05,0.73,0.2);
	legendChargedRAAcomp3->SetFillColor(0);
	legendChargedRAAcomp3->SetFillStyle(0);
	legendChargedRAAcomp3->SetLineColor(0);
	legendChargedRAAcomp3->SetTextSize(0.035);
	legendChargedRAAcomp3->SetMargin(0.2);
// 	legendChargedRAAcomp3->AddEntry(graphChargedKaonRAA0005,"K^{#pm} R_{AA} - 0-5%");
	legendChargedRAAcomp3->AddEntry(graphChargedKaonRAA0010,"K^{#pm} R_{AA} (PWGLF-258) - 0-10%");
// 	legendChargedRAAcomp3->AddEntry(graphChargedPionRAASys0005,"#pi^{#pm} R_{AA} - 0-5%");
// 	legendChargedRAAcomp3->AddEntry(graphEtaRAA0005,"PCM measured #eta R_{AA} - 0-5%", "pl");
	legendChargedRAAcomp3->AddEntry(graphEtaRAA0010,"PCM #eta R_{AA} (2011) - 0-10%", "pl");
// 	legendChargedRAAcomp3->AddEntry(graphChargedKaonRAASys2040,"K^{#pm} R_{AA} - 20-40%");
// 	legendChargedRAAcomp3->AddEntry(graphChargedPionRAASys2040,"#pi^{#pm} R_{AA} - 20-40%");
// 	legendChargedRAAcomp3->AddEntry(graphPi0RAASys2040,"PCM measured #pi^{0} R_{AA} - 20-40%");
	legendChargedRAAcomp3->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAEtaWithChargedKaons0010.%s",outputDir.Data(),suffix.Data()));

	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");


	for(Int_t j = 0; j<6; j++){
		graphChargedKaonRAA2040->RemovePoint(0);
	}

	for(Int_t i = 0; i<8; i++){
		graphChargedKaonRAA2040->RemovePoint(graphChargedKaonRAA2040->GetN()-1);
		graphChargedKaonRAASys2040->RemovePoint(graphChargedKaonRAASys2040->GetN()-1);
	}
	
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAA2040, 20,2, kCyan+1, kCyan+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRAASys2040,20,2, kCyan+2, kCyan+2, 2, kTRUE, kCyan-8);
	graphChargedKaonRAASys2040->Draw("2same");
	graphChargedKaonRAA2040->Draw("p,same"); 
		
// 	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAA2040, 20,1, kBlue, kBlue, 0.1, kFALSE);
// 	DrawGammaSetMarkerTGraphAsym(graphChargedPionRAASys2040,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
// 	graphPi0RAASys2050->SetFillStyle(3008);
// 	graphChargedPionRAASys2040->Draw("2same");
// 	graphChargedPionRAA2040->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA2040, 20,2, kGreen+2, kGreen+2, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2040,20,2, kGreen+2, kGreen+2, 2, kTRUE, kGreen-6);
	graphEtaRAASys2040->SetFillStyle(0);
	graphEtaRAASys2040->Draw("2same");
	graphEtaRAA2040->Draw("p,same");
	
	TLegend* legendChargedRAAcomp4 = new TLegend(0.4,0.75,0.9,0.9);  //0.16,0.05,0.73,0.2);
	legendChargedRAAcomp4->SetFillColor(0);
	legendChargedRAAcomp4->SetFillStyle(0);
	legendChargedRAAcomp4->SetLineColor(0);
	legendChargedRAAcomp4->SetTextSize(0.035);
	legendChargedRAAcomp4->SetMargin(0.2);
	legendChargedRAAcomp4->AddEntry(graphChargedKaonRAA2040,"K^{#pm} R_{AA} (PWGLF-258) - 20-40%");
// 	legendChargedRAAcomp->AddEntry(graphChargedPionRAASys2040,"#pi^{#pm} R_{AA} - 20-40%");
	legendChargedRAAcomp4->AddEntry(graphEtaRAA2040,"PCM #eta R_{AA} (2011) - 20-40%", "lp");
	legendChargedRAAcomp4->Draw();
	
	canvasRAA->Update();
	canvasRAA->SaveAs(Form("%s/RAAEtaWithChargedKaons2040.%s",outputDir.Data(),suffix.Data()));

	
	
	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	
	DrawGammaSetMarkerTGraphAsym(graphPi0RAA0010, 20,2, kOrange, kOrange, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys0010,20,2, kOrange, kOrange, 2, kTRUE, kOrange-2);
	graphPi0RAASys0010->SetFillStyle(0);
	graphPi0RAASys0010->Draw("2same");
	graphPi0RAA0010->Draw("p,same"); 

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA0010, 21,2, kRed, kRed, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys0010,21,2, kRed, kRed, 2, kTRUE, kRed-9);
	graphEtaRAASys0010->SetFillStyle(0);
	graphEtaRAASys0010->Draw("2same");
	graphEtaRAA0010->Draw("p,same"); 
	
	DrawGammaSetMarkerTGraphAsym(graphEtaRAA2050, 20,2, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2050,20,2, kAzure, kAzure, 2, kTRUE, kAzure-8);
	graphEtaRAASys2050->SetFillStyle(0);
	graphEtaRAASys2050->Draw("2,same");
	graphEtaRAA2050->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphPi0RAA2050, 21,2, kTeal-6, kTeal-6, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphPi0RAASys2050,21,2, kTeal-6, kTeal-6, 2, kTRUE, kTeal-9);
	graphPi0RAASys2050->SetFillStyle(0);
	graphPi0RAASys2050->Draw("2same");
	graphPi0RAA2050->Draw("p,same");

	
	TLegend* legendPi0EtaRAA = new TLegend(0.4,0.72,0.9,0.9);  //0.16,0.05,0.73,0.2);
	legendPi0EtaRAA->SetFillColor(0);
	legendPi0EtaRAA->SetFillStyle(0);
	legendPi0EtaRAA->SetLineColor(0);
	legendPi0EtaRAA->SetTextSize(0.035);
	legendPi0EtaRAA->SetMargin(0.2);
	legendPi0EtaRAA->AddEntry(graphPi0RAA0010,"PCM #pi^{0} R_{AA} (2010) - 0-10%", "lp");
	legendPi0EtaRAA->AddEntry(graphEtaRAA0010,"PCM #eta R_{AA} (2011) - 0-10%", "lp");
	legendPi0EtaRAA->AddEntry(graphPi0RAA2050,"PCM #pi^{0} R_{AA} (2010) - 20-50%", "lp");
	legendPi0EtaRAA->AddEntry(graphEtaRAA2050,"PCM #eta R_{AA} (2011) - 20-50%", "lp");
	legendPi0EtaRAA->Draw();
	
	canvasRAA->Update();
// 	canvasRAA->SaveAs(Form("%s/RAAPi0andEtaPbPbLHC11h.%s",outputDir.Data(),suffix.Data()));
	
	

	canvasRAA->cd();
	histo2DRAA->Draw("copy");
	

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA0010, 21,2, kMagenta+1, kMagenta+1, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys0010,21,2, kMagenta+1, kMagenta+1, 2, kTRUE, kMagenta-7);
	graphEtaRAASys0010->SetFillStyle(0);
	graphEtaRAASys0010->Draw("2same");
	graphEtaRAA0010->Draw("p,same"); 


	DrawGammaSetMarkerTGraphAsym(graphEtaRAA2050, 20,1, kBlue, kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2050,20,1, kAzure, kAzure, 2, kTRUE, kAzure-8);
	graphEtaRAASys2050->SetFillStyle(0);
// 	graphEtaRAASys2050->Draw("2same");
// 	graphEtaRAA2050->Draw("p,same");

	DrawGammaSetMarkerTGraphAsym(graphEtaRAA2040, 21,2, kBlue,kBlue, 0.1, kFALSE);
	DrawGammaSetMarkerTGraphAsym(graphEtaRAASys2040,21,2, kAzure-9, kAzure-9, 2, kTRUE, kAzure-9);
	graphEtaRAASys2040->SetFillStyle(0);
	graphEtaRAASys2040->Draw("2same");
	graphEtaRAA2040->Draw("p,same");
	
	graphPi0RAASys0010->Draw("E2same");
	graphPi0RAA0010->Draw("p,same,e1");
	
	graphPi0RAASys2040->Draw("E2same");
	graphPi0RAA2040->Draw("p,same,e1");
	
	graphEtaRAAPHENIX0010->Draw("2same");
	graphEtaRAAPHENIX0010->Draw("p,same"); 
	graphEtaRAAPHENIX2040->Draw("2same");
	graphEtaRAAPHENIX2040->Draw("p,same"); 

	
	TLegend* legendEtaRAAcomp = new TLegend(0.4,0.65,0.9,0.9);  //0.16,0.05,0.73,0.2);
	legendEtaRAAcomp->SetFillColor(0);
	legendEtaRAAcomp->SetFillStyle(0);
	legendEtaRAAcomp->SetLineColor(0);
	legendEtaRAAcomp->SetTextSize(0.035);
	legendEtaRAAcomp->SetMargin(0.2);
	legendEtaRAAcomp->AddEntry(graphPi0RAA0010,"PCM #pi^{0} R_{AA} (2010) - 0-10%", "lp");
	legendEtaRAAcomp->AddEntry(graphPi0RAA2040,"PCM #pi^{0} R_{AA} (2010) - 20-40%","pl");
	legendEtaRAAcomp->AddEntry(graphEtaRAAPHENIX0010,"PHENIX #eta R_{AA} - 0-10%","pl");
	legendEtaRAAcomp->AddEntry(graphEtaRAAPHENIX2040,"PHENIX #eta R_{AA} - 20-40%","pl");
// 	legendEtaRAAcomp->AddEntry(graphEtaRAAPHENIX2060,"PHENIX #eta R_{AA} - 20-60%");
	legendEtaRAAcomp->AddEntry(graphEtaRAA0010,"PCM #eta R_{AA} (2011) - 0-10%", "pl");
	legendEtaRAAcomp->AddEntry(graphEtaRAA2040,"PCM #eta R_{AA} (2011) - 20-40%", "pl");
// 	legendEtaRAAcomp->AddEntry(graphEtaRAA2050,"PCM measured #eta R_{AA} - 20-50%","pl");
	legendEtaRAAcomp->Draw();
	
	canvasRAA->Update();
// 	canvasRAA->SaveAs(Form("%s/RAAPi0andEtaWithPHENIX.%s",outputDir.Data(),suffix.Data()));
	
	

		
	TString  nameFinalResDat = Form("%s/FitResultsMC.dat",outputDir.Data());
	TString forOutput;
	fstream  fileFinalResults;
	fileFinalResults.open(nameFinalResDat.Data(), ios::out);
	  
	TCanvas* canvasSpectraMC = new TCanvas("canvasSpectraMC","",200,10,1200,1100);  // gives the page size
	DrawGammaCanvasSettings( canvasSpectraMC,  0.13, 0.01, 0.015, 0.08);
	canvasSpectraMC->SetLogy();
	canvasSpectraMC->SetLogx();
	TH2F * histo2DSpectraMC = new TH2F("histo2DSpectraMC","histo2DSpectraMC",1000,0.01,30.,1000,1e-8,2e4 );
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
	//**************************** Pi0 reweighting evalulation 0-10% LHC11h*************
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

		
		histoMCYieldPi0PtPbPbLHC11hAddedSig0010->SetMarkerStyle(markerStylePbPb4060MC);
		TH1D* histoRatioDatatoFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoPCMPi0CorrectedSpecPbPbLHC11h0010, fitYieldDataQCDPi0PbPbLHC11h0010);
		TH1D* histoRatioMCtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11h0010, fitYieldDataQCDPi0PbPbLHC11h0010);
		TH1D* histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010 = CalculateHistoRatioToFit (histoMCYieldPi0PtPbPbLHC11hAddedSig0010, fitYieldDataQCDPi0PbPbLHC11h0010);
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
		
		DrawGammaSetMarker(histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010, markerStylePbPb4060MC,markerSizePbPb4060, colorCombPbPb4060 , colorCombPbPb4060);
		histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		
		if (histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010) histoRatioMCUnweightedtoDataFitQCDPbPbLHC11h0010->Draw("same,e,p");  
		
		TLegend* legendFit0010 = new TLegend(0.16,0.73,0.4,0.9);
		legendFit0010->SetFillColor(0);
		legendFit0010->SetLineColor(0);
		legendFit0010->SetTextSize(0.035);
		legendFit0010->SetMargin(0.2);
		legendFit0010->AddEntry(histoRatioDatatoFitQCDPbPbLHC11h0010,"Data/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit0010->AddEntry(histoRatioMCtoDataFitQCDPbPbLHC11h0010,"MC weighted/QCD fit to Data (0.4 < pT < 14)","p");
		if (runDrawReweighted) legendFit0010->AddEntry(histoRatioMCAddedSigtoDataFitQCDPbPbLHC11h0010,"MC added sig. weighted/QCD fit to Data (0.4 < pT < 14)","p");
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
	//**************************** Eta reweighting evalulation 0-10% LHC11h*************
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
		
		fitYieldDataQCDEtaPbPbLHC11h0005 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h0005","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h0005,1.4,10,NULL,"QNRME+");
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
		
		fitYieldDataQCDEtaPbPbLHC11h0510 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h0510","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h0510,1.4,10,NULL,"QNRME+");
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
		
		fitYieldDataQCDEtaPbPbLHC11h2040 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h2040","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h2040,1.4,10,NULL,"QNRME+");
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
	//**************************** Eta reweighting evalulation 20-50% LHC11h*************
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
		
		fitYieldDataQCDEtaPbPbLHC11h2050 = FitObject("l","fitYieldDataQCDEtaPbPbLHC11h2050","Eta",histoPCMEtaCorrectedSpecPbPbLHC11h2050,1.4,10,NULL,"QNRME+");
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

