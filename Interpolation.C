
/****************************************************************************************************************************
******          provided by Gamma Conversion Group, PWG-GA*****
******          Pedro Gonzalez, pedro.gonzalez.zamora@cern.ch
******          Ana Marin, marin@physi.uni-heidelberg.de    
*****           Annika Passfeld, annikapassfeld@uni-muenster,de                                                                                                  *****
******          Friederike Bock, friederike.bock@cern.ch                                                                                                        *****
******          
*****************************************************************************************************************************/

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
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"
#include "CommonHeaders/Interpolation.h"
#include "CommonHeaders/PlottingInterPolationRpPb.h"
#include "TFitResultPtr.h"



/*extern TRandom*         gRandom;	
extern TBenchmark*      gBenchmark;
extern TSystem*         gSystem;
extern TMinuit*         gMinuit;

void CalcRpPb(TGraphAsymmErrors* PPSpectrumSystErr, TGraphAsymmErrors*  PPSpectrumStatErr, TGraphAsymmErrors* pPbSpectrumSystErr, TGraphAsymmErrors* pPbSpectrumStatErr,
	      TGraphAsymmErrors** graphRpPbSystErr, TGraphAsymmErrors** graphRpPbStatErr);
//TGraphAsymmErrors* CalculateErrorsStat(TGraphErrors* spectrum);
TGraphAsymmErrors* CalculateSystErrors(TGraphErrors* spectrum,TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);
TF1* RebinWithFitToTGraph(TGraphAsymmErrors *spectrum, TGraphAsymmErrors** newSpectrum, TGraphAsymmErrors *newBins, TString FitType,Double_t minPt, Double_t maxPt, Double_t* parameters,Double_t probability);
TGraphErrors *GetInterpolSpectrum2D(TGraphErrors *g1, TGraphErrors *g2,Double_t d1, Double_t d2,Double_t dSqrts);
TGraphAsymmErrors* GetChargeParticlesRpPb();
TGraphErrors *ConvertTGraphAsymmErrorstoTGraphErrors(TGraphAsymmErrors* g1);*/



	      




void Interpolation(TString resultsType ="PCMPileUpCorrection", TString suffix="pdf", TString outputDir="OutputRpPb",TString FitFuncName="Bylinkin", TString thesisPlots="" )
{

	
	gROOT->SetStyle("Plain");
        TH1::AddDirectory(kFALSE);
	
	if(thesisPlots.EqualTo("thesis") ){
	  
	StyleSettingsThesis(suffix);
	cout<<"Entre a thesis plots"<<endl;
	}else {
	
	StyleSettings();
	}
	
	
	//Defining colors
	
	colorsArray[kDalitz] = kCyan+2;
	colorsArray[kPCM]    = kBlack;
	colorsArray[kPHOS]   = kRed+2;
	colorsArray[kPion]   = kGreen-3;
	colorsArray[kKaon]   = kOrange-3;
	colorsArray[kProton] = kMagenta-3;
	
	
	
	TString dateForOutput = ReturnDateStringForOutput();
	outputDir = Form("%s/%s/%s/%s/%s",outputDir.Data(),dateForOutput.Data(),suffix.Data(),resultsType.Data(),FitFuncName.Data());
	gSystem->Exec("mkdir -p "+outputDir);
	
	Bool_t PileUpCorrection = kFALSE;
	Bool_t thesis = kFALSE;

	
	TString fileNameNeutralPionCombResultsPP	="FinalResults/CombinedResultsPP_ShiftedX_PaperRAA_16_May_2014.root";
	//TString fileNameNeutralPionPCMResultspPb	="ExternalInputpPb/PCM/data_PCMResults_pPb_20141023.root";
	//TString fileNameNeutralPionPCMResultspPb        ="ExternalInputpPb/PCM/data_PCMResults_pPb_20150804_standard_CatErrors.root";
	TString fileNameNeutralPionPCMResultspPb        ="ExternalInputpPb/PCM/data_PCMResults_pPb_20151111_standard_CatErrors.root";
	TString fileNameNeutralPionDalitzResultspPb	="ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20150806.root";
	//TString fileNameNeutralPionDalitzResultspPb	="ExternalInputpPb/PCM/data_PCMResults_Dalitz_pPb_20150901.root";
	
	
	
	if( resultsType.CompareTo("PCMPileUpCorrection") == 0 ){
	  
	       PileUpCorrection = kTRUE;
	       
	       cout<<"PileUp correction will be applied to PCM alone spetrum. "<<endl;
	       
	     	  
	}  
	
	if( thesisPlots.CompareTo("thesis") == 0 ) {
	  
	      thesis = kTRUE;
	      thesisPlotLabel = "This thesis";
	      
	} else {
	  
	      thesisPlotLabel = "ALICE work in progress";
	}
	  
	
	
	
	
	
        const char *fileNameEPS09sPi0AKK = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_akk_mb.dat";
	const char *fileNameEPS09sPi0DSS = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_fdss_mb.dat";
	const char *fileNameEPS09sPi0KKP = "ExternalInputpPb/Theory/R_pi0_pPb_y0_eps09s/R_pPb5000_y0_eps09s_kkp_mb.dat";
	const char *fileNameEPS09sPi0CGC = "ExternalInputpPb/Theory/ColorGlassCondensate.dat";
	
	
	
	
	
	
	TString fileNameRpPbChargePionsNew  	     	= "ExternalInputpPb/RpPb_502_PiKp_w_Kinks_12062015.root";
	TString histoNameChargedPionsRpPbSystErr     	= "hsys_RpPb_pion";
	TString histoNameChargedPionsRpPbStatErr     	= "hstat_RpPb_pion";
	TString histoNameChargedParticlesRpPbSystErr 	= "hsys_RpPb_charged";
	TString histoNameChargedParticlesRpPbStatErr 	= "hstat_RpPb_charged";
	TString histoNameChargedKaonsRpPbStatErr     	= "hstat_RpPb_kaon";
	TString histoNameChargedKaonsRpPbSystErr 	= "hsys_RpPb_kaon";
	TString histoNameChargedProtonsRpPbStatErr      = "hstat_RpPb_proton";
	TString histoNameChargedProtonsRpPbSystErr      = "hsys_RpPb_proton";
	
	
	 TFile*  fileRpPbChargedPionsNew = new TFile( fileNameRpPbChargePionsNew.Data() );
	
	
	
	TH1D*   histo_pi_RpPb_SystErr     = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedPionsRpPbSystErr.Data());
	TH1D*   histo_pi_RpPb_StatErr     = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedPionsRpPbStatErr.Data());
	TH1D*   histo_kaon_RpPb_StatErr   = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedKaonsRpPbStatErr.Data());
	TH1D*   histo_kaon_RpPb_SystErr   = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedKaonsRpPbSystErr.Data());
	TH1D*   histo_proton_RpPb_StatErr = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedProtonsRpPbStatErr.Data());
	TH1D*   histo_proton_RpPb_SystErr = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedProtonsRpPbSystErr.Data());
	
	TGraphAsymmErrors*  graphChargedPionRpPbSystErr   = new TGraphAsymmErrors(histo_pi_RpPb_SystErr);
	TGraphAsymmErrors*  graphChargedPionRpPbStatErr   = new TGraphAsymmErrors(histo_pi_RpPb_StatErr);
	TGraphAsymmErrors*  graphChargedPionRpPbSystErrClone  = (TGraphAsymmErrors*) graphChargedPionRpPbSystErr->Clone();
	TGraphAsymmErrors*  graphChargedPionRpPbStatErrClone  = (TGraphAsymmErrors*) graphChargedPionRpPbStatErr->Clone();
	
	
	TGraphAsymmErrors*  graphChargedKaonRpPbStatErr   = new TGraphAsymmErrors(histo_kaon_RpPb_StatErr);
	TGraphAsymmErrors*  graphChargedKaonRpPbSystErr   = new TGraphAsymmErrors(histo_kaon_RpPb_SystErr);
	TGraphAsymmErrors*  graphChargedProtonRpPbStatErr = new TGraphAsymmErrors(histo_proton_RpPb_StatErr);
	TGraphAsymmErrors*  graphChargedProtonRpPbSystErr = new TGraphAsymmErrors(histo_proton_RpPb_SystErr);
	
	
	
	
	TGraphAsymmErrors* graphAsymmChargedParticlesRpPbSystErr = GetChargeParticlesRpPb("Syst");
	TGraphAsymmErrors* graphAsymmChargedParticlesRpPbStatErr = GetChargeParticlesRpPb("Stat");
	
	//TH1D*  histo_charged_RpPb_SystErr = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedParticlesRpPbSystErr.Data());
	//TH1D*  histo_charged_RpPb_StatErr = (TH1D*)fileRpPbChargedPionsNew->Get(histoNameChargedParticlesRpPbStatErr.Data());
	
	
	//TGraphAsymmErrors* graphAsymmChargedParticlesRpPbSystErr = new TGraphAsymmErrors(histo_charged_RpPb_SystErr);
	//TGraphAsymmErrors* graphAsymmChargedParticlesRpPbStatErr = new TGraphAsymmErrors(histo_charged_RpPb_StatErr);
	
	//PHOS Input

	TString fileNameRpPbPHOS  	     = "ExternalInputpPb/PHOS/data_PHOSResults_RpPb_20160208.root";
	TString NamePHOSRpPbSystErr     = "graphRpPbPHOSSystErr";
	TString NamePHOSRpPbStatErr     = "graphRpPbPHOSStatErr";

	 TFile*  fileRpPbPHOS = new TFile( fileNameRpPbPHOS.Data() );
       
	TGraphAsymmErrors*  graphRpPbPHOSSystErr = (TGraphAsymmErrors*)fileRpPbPHOS->Get(NamePHOSRpPbSystErr);
	TGraphAsymmErrors*  graphRpPbPHOSStatErr = (TGraphAsymmErrors*)fileRpPbPHOS->Get(NamePHOSRpPbStatErr);
	
	graphRpPbPHOSSystErr->RemovePoint(0);	
	
	
	TString fileNamePythia8 = "ExternalInputpPb/pythia8-5.02TeV-part-1.root";
	
	TString fitType  = "l";
	TString fitNameLabel = "Tsallis";
	
	Double_t *Parameters2760GeV = NULL;
	Double_t *Parameters7TeV    = NULL;
	
	
	if( FitFuncName.CompareTo("Tsallis") == 0 ){
	  
	  fitType = "l";
	  Parameters2760GeV = new Double_t[3];// = {2.4e+10,6.88,0.139}; //2.,5.,0.18
	  Parameters2760GeV[0] = 2.4e+10;
	  Parameters2760GeV[1] = 6.88;
	  Parameters2760GeV[2] = 0.139;
	  
	  Parameters7TeV    = new Double_t[3];// = {13.7e+09,7.,0.13};
	  Parameters7TeV[0] = 13.7e+09;
	  Parameters7TeV[1] = 7.0;
	  Parameters7TeV[2] = 0.13;
	  
	  fitNameLabel = "Tsallis";
	  
	  
	} else if ( FitFuncName.CompareTo("Bylinkin") == 0){
	    
	  fitType = "tcm";
	  
	  
	  Parameters2760GeV = new Double_t[5];// = { 2.4e+10, 0.3, 1e+10,0.3,8};
	  Parameters2760GeV[0] = 2.4e+10;
	  Parameters2760GeV[1] = 0.3;
	  Parameters2760GeV[2] = 1e+10;
	  Parameters2760GeV[3] = 0.3;
	  Parameters2760GeV[4] = 3.8;
	  
	  Parameters7TeV    = new Double_t[5];// = { 7.4e+10, 0.3, 1e+09,0.3,8};
	  Parameters7TeV[0] = 7.4e+10;
	  Parameters7TeV[1] = 0.3;
	  Parameters7TeV[2] = 1e+09;
	  Parameters7TeV[3] = 0.3;
	  Parameters7TeV[4] = 3.8;
	  
	  fitNameLabel = "Bylinkin-Rostovtsev";
	  
         }
	
	
	
        TH1F* ppPythia = GetPPReferenceFromPythia(fileNamePythia8.Data());
	TGraphAsymmErrors*  graphppPythia = new TGraphAsymmErrors(ppPythia);
	
	
	graphppPythia->Print();
	

	
        TFile*  fileNeutralPionCombResultsPP = new TFile( fileNameNeutralPionCombResultsPP.Data() );
	
	
	
	
	TString graphNameInvCrossSectionPi07TeVStatSystErr;
	TString graphNameInvCrossSectionPi07TeVStatErr;
	TString graphNameInvCrossSectionPi07TeVSysErr;
	
	
	TString graphNameInvCrossSectionPi02760GeVStatSystErr;
	TString graphNameInvCrossSectionPi02760GeVStatErr;
	TString graphNameInvCrossSectionPi02760GeVSysErr;
	
	
	TGraphAsymmErrors*  graphInvCrossSectionPi07TeVStatSystErr;
	TGraphAsymmErrors*  graphInvCrossSectionPi07TeVStatErr;
	TGraphAsymmErrors*  graphInvCrossSectionPi07TeVSystErr;
	
	TGraphAsymmErrors*  graphInvCrossSectionPi02760GeVStatSystErr;   	
        TGraphAsymmErrors*  graphInvCrossSectionPi02760GeVStatErr;   	
	TGraphAsymmErrors*  graphInvCrossSectionPi02760GeVSystErr;   	
	
	
	
	
	if ( resultsType.CompareTo("Comb") == 0   ){
	  
	
	graphNameInvCrossSectionPi07TeVStatSystErr 	= Form("graphInvCrossSectionPi0%s7TeV",resultsType.Data()); 
	graphNameInvCrossSectionPi07TeVStatErr     	= Form("graphInvCrossSectionPi0%s7TeVStatErr",resultsType.Data());
	graphNameInvCrossSectionPi07TeVSysErr      	= Form("graphInvCrossSectionPi0%s7TeVSysErr",resultsType.Data());
	
	graphNameInvCrossSectionPi02760GeVStatSystErr 	= Form("graphInvCrossSectionPi0%s2760GeV",resultsType.Data());
	graphNameInvCrossSectionPi02760GeVStatErr 	= Form("graphInvCrossSectionPi0%s2760GeVStatErr",resultsType.Data());
	graphNameInvCrossSectionPi02760GeVSysErr  	= Form("graphInvCrossSectionPi0%s2760GeVSysErr",resultsType.Data());
	
	
	graphInvCrossSectionPi07TeVStatSystErr      	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVStatSystErr.Data());
	graphInvCrossSectionPi07TeVStatErr          	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVStatErr.Data());
	graphInvCrossSectionPi07TeVSystErr          	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVSysErr.Data());
	
	
	cout<<"************************************************Checking if they are shifted Combined**********************************************************"<<endl;
	graphInvCrossSectionPi07TeVStatErr->Print();
	
	cout<<"**************************************************"<<endl;
	
	graphInvCrossSectionPi07TeVSystErr->Print();
	
	
        graphInvCrossSectionPi02760GeVStatSystErr  	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVStatSystErr.Data());	
        graphInvCrossSectionPi02760GeVStatErr      	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVStatErr.Data());	
	graphInvCrossSectionPi02760GeVSystErr      	= (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVSysErr.Data());	
	
	cout<<"Entro a Comb"<<endl;
	
	
	} else if(  resultsType.CompareTo("PCM") == 0 || resultsType.CompareTo("PCMPileUpCorrection") == 0 ){

	
	  
	
	graphNameInvCrossSectionPi07TeVStatErr    	= "graphInvCrossSectionPi0PCMStat7TeV";
	graphNameInvCrossSectionPi07TeVSysErr     	= "graphInvCrossSectionPi0PCMSys7TeV";
	graphNameInvCrossSectionPi02760GeVStatErr 	= "graphInvCrossSectionPi0PCM2760GeVStatErr";
	graphNameInvCrossSectionPi02760GeVSysErr  	= "graphInvCrossSectionPi0PCM2760GeVSysErr";
	
	
	
	
	graphInvCrossSectionPi07TeVStatErr        = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVStatErr.Data());
	graphInvCrossSectionPi07TeVSystErr        = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi07TeVSysErr.Data());
	
	cout<<"************************************************Checking if they are shifted PCM alone**********************************************************"<<endl;
	graphInvCrossSectionPi07TeVStatErr->Print();
	
	cout<<"**************************************************"<<endl;
	
	graphInvCrossSectionPi07TeVSystErr->Print();
	
	
	
	TGraphAsymmErrors* graphInvCrossSectionPi07TeVStatErrNoPileUpCorr = (TGraphAsymmErrors*)graphInvCrossSectionPi07TeVStatErr->Clone();
	
	if( PileUpCorrection == kTRUE ){
	  
	cout<<"entro a corregir"<<endl;  
	  
	TF1* fitCorrectionFactorsHistvsPtRatio 								= new TF1("fitCorrectionFactorsHistvsPtRatio","[0]/pow(x,[1])+[2]");
	fitCorrectionFactorsHistvsPtRatio->SetParameter(0,2.9737546081);
	fitCorrectionFactorsHistvsPtRatio->SetParameter(1,1.4795520406);
	fitCorrectionFactorsHistvsPtRatio->SetParameter(2,2.2652589579);
	
	Double_t* yStatNoPileUpCorr                             = graphInvCrossSectionPi07TeVStatErrNoPileUpCorr->GetY();

	Double_t* yStat 					= graphInvCrossSectionPi07TeVStatErr->GetY();
	Double_t* xStat 					= graphInvCrossSectionPi07TeVStatErr->GetX();
	
	
	Double_t* ySyst 					= graphInvCrossSectionPi07TeVSystErr->GetY();
	Double_t* xSyst 					= graphInvCrossSectionPi07TeVSystErr->GetX();
	
	const Int_t nBins = graphInvCrossSectionPi07TeVStatErr->GetN();
	
	Double_t ratioX[nBins];
	Double_t ratioY[nBins];
	Double_t errorX[nBins];
	Double_t errorY[nBins];
	
	
	cout<<"Applyin pileUp correction"<<endl;
	
	for (Int_t i = 0; i < nBins; i++){
	  
		cout << xStat[i] << "\t" << 100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]) << endl;
		cout << yStat[i] << "\t" << yStat[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]))/100 << "\t" << yStat[i] << endl;
		yStat[i] 	= yStat[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xStat[i]))/100;
		
		cout << xSyst[i] << "\t" << 100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]) << endl;
		cout << ySyst[i] << "\t" << ySyst[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]))/100 << "\t" << ySyst[i] << endl;
		ySyst[i] 	= ySyst[i]*(100-fitCorrectionFactorsHistvsPtRatio->Eval(xSyst[i]))/100;
		
		
		ratioX[i] = xSyst[i];
		errorX[i] = graphInvCrossSectionPi07TeVStatErr->GetErrorXhigh(i);
		ratioY[i] = yStat[i]/yStatNoPileUpCorr[i];
		errorY[i] = TMath::Sqrt(TMath::Power(graphInvCrossSectionPi07TeVStatErr->GetErrorYhigh(i)/ySyst[i],2) +TMath::Power(graphInvCrossSectionPi07TeVStatErrNoPileUpCorr->GetErrorYhigh(i)/yStatNoPileUpCorr[i],2))*ratioY[i];
		cout << "Ratio: " << ratioX[i] << "\t" <<  errorX[i] << "\t" << ratioY[i] << "\t" << errorY[i] << "\t" << errorY[i]/ratioY[i]*100 << "%"<<endl;

		
		
	}   
	cout << "*************************" << endl;
	
	TGraphErrors* RatioNoPileUpAndPileUpCorr =  new TGraphErrors(nBins-1,ratioX,ratioY,errorX,errorY); 
	
	
	
	TCanvas* canvasRatioNoPileUpAndPileUpCorr = new TCanvas("canvasRatioNoPileUpAndPileUpCorr","Ratio between PileUp Correction And No PileUp correction spectra",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioNoPileUpAndPileUpCorr,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioNoPileUpAndPileUpCorr = new TH2F("histo2DRatioNoPileUpAndPileUpCorr","histo2DRatioNoPileUpAndPileUpCorr",1000,0.,15.,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioNoPileUpAndPileUpCorr, "#it{p}_{T} (GeV/#it{c})","PileUp Corr/No Corr", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioNoPileUpAndPileUpCorr->GetYaxis()->SetRangeUser(0.5,1.5);
	histo2DRatioNoPileUpAndPileUpCorr->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(RatioNoPileUpAndPileUpCorr,21,1.5, kBlack , kBlack);
	RatioNoPileUpAndPileUpCorr->Draw("E1psame");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioNoPileUpAndPileUpCorr = new TLatex(0.16,0.9,"pp minimum bias #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( textRatioNoPileUpAndPileUpCorr,0.06,4); 
	textRatioNoPileUpAndPileUpCorr->Draw();
	
	
	TLegend* legendRatioNoPileUpAndPileUpCorr = new TLegend(0.18,0.15,0.9,0.21);
	legendRatioNoPileUpAndPileUpCorr->SetFillColor(0);
	legendRatioNoPileUpAndPileUpCorr->SetLineColor(0);
	legendRatioNoPileUpAndPileUpCorr->SetNColumns(2);
	legendRatioNoPileUpAndPileUpCorr->SetTextSize(0.045);
	legendRatioNoPileUpAndPileUpCorr->AddEntry(histo2DRatioNoPileUpAndPileUpCorr,"PileUp Correction / No PileUp Correction ","p");
	legendRatioNoPileUpAndPileUpCorr->Draw();
	
	
	canvasRatioNoPileUpAndPileUpCorr->SaveAs(Form("%s/Ratio_NoPileUpAndPileUpCorr.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	  
	  
	}
	
	
	graphInvCrossSectionPi07TeVStatSystErr    =  CalculateCombinedSysAndStatError( graphInvCrossSectionPi07TeVStatErr , graphInvCrossSectionPi07TeVSystErr);
	
	graphInvCrossSectionPi02760GeVStatErr     = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVStatErr.Data());	
	graphInvCrossSectionPi02760GeVSystErr     = (TGraphAsymmErrors*)fileNeutralPionCombResultsPP->Get(graphNameInvCrossSectionPi02760GeVSysErr.Data());	
	
	graphInvCrossSectionPi02760GeVStatSystErr =  CalculateCombinedSysAndStatError( graphInvCrossSectionPi02760GeVStatErr , graphInvCrossSectionPi02760GeVSystErr);
	
	cout<<"Entro a PCM"<<endl;
	  
	  
	} else {
	  
	  cout<<"Invalid results type:  "<<resultsType.Data() <<endl;
	  return;
	  
	}
	  
     
	
	cout<<"Llego 0 "<<endl;

        //return;
	
	TFile*  fileNeutralPionPCMResultspPb = new TFile(fileNameNeutralPionPCMResultspPb.Data());
	TDirectory* fNeutralPionPCMResultspPbContainer = (TDirectory*) fileNeutralPionPCMResultspPb->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	if( ! fNeutralPionPCMResultspPbContainer ) {cout<<"TList fNeutralPionPCMResultspPbContainer does not exist: "<<endl; return;}
	
	
	
	TH1F* histoInvYieldPi0PCMpPb5023GeV              	   = (TH1F*)fNeutralPionPCMResultspPbContainer->Get("CorrectedYieldPi0");
        TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVSystErr    = (TGraphAsymmErrors*)fNeutralPionPCMResultspPbContainer->Get("Pi0SystError"); 
	TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVComplErr   = (TGraphAsymmErrors*)fNeutralPionPCMResultspPbContainer->Get("Pi0ComplError");
 	TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVStatErr    = new TGraphAsymmErrors(histoInvYieldPi0PCMpPb5023GeV);
	
	
	//********************************************************************
	//NOTE Temporal this should be fixed in the ProducedFinalResultspPb
	
	graphInvYieldPi0PCMpPb5023GeVStatErr->RemovePoint(0);
	graphInvYieldPi0PCMpPb5023GeVStatErr->RemovePoint(0);
	graphInvYieldPi0PCMpPb5023GeVStatErr->RemovePoint(0);
	
	//cout<<"Annika Annika"<<endl;
	//graphInvYieldPi0PCMpPb5023GeVStatErr->Print();
	//graphInvYieldPi0PCMpPb5023GeVSystErr->Print();
	
	
	graphInvYieldPi0PCMpPb5023GeVSystErr->RemovePoint(0);
	graphInvYieldPi0PCMpPb5023GeVSystErr->RemovePoint(0);

	
	
	graphInvYieldPi0PCMpPb5023GeVComplErr->RemovePoint(0);
	graphInvYieldPi0PCMpPb5023GeVComplErr->RemovePoint(0);
	
	
	cout<<"Annika Annika"<<endl;
	graphInvYieldPi0PCMpPb5023GeVStatErr->Print();
	graphInvYieldPi0PCMpPb5023GeVSystErr->Print();
	graphInvYieldPi0PCMpPb5023GeVComplErr->Print();
	
	
	//**********************************************************************
	
	
	//*********************************************************************
	
	
	cout<<"LLego 1"<<endl;
	
	
	TFile*  fileNeutralPionDalitzResultspPb = new TFile(fileNameNeutralPionDalitzResultspPb.Data());
	
	TDirectory* fNeutralPionDalitzResultspPbContainer = (TDirectory*) fileNeutralPionDalitzResultspPb->GetDirectory("Pi0_pPb_5.023TeV_0-100%");
	
	if( ! fNeutralPionDalitzResultspPbContainer ) {cout<<"TList fNeutralPionDalitzResultspPbContainer does not exist: "<<endl; return;}
	 

	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVComplErr  = (TGraphAsymmErrors*)fNeutralPionDalitzResultspPbContainer->Get("Pi0ComplError");
	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVSystErr   = (TGraphAsymmErrors*)fNeutralPionDalitzResultspPbContainer->Get("Pi0SystError");
	//TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVStatErr   = (TGraphAsymmErrors*)fNeutralPionDalitzResultspPbContainer->Get("Pi0StatError");
	
	TH1F* histoInvYieldPi0DalitzpPb5023GeV  = (TH1F*)fNeutralPionDalitzResultspPbContainer->Get("CorrectedYieldPi0");
	     
	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVStatErr = new TGraphAsymmErrors( histoInvYieldPi0DalitzpPb5023GeV );
	graphInvYieldPi0DalitzpPb5023GeVStatErr->RemovePoint(0);
	//graphInvYieldPi0DalitzpPb5023GeVStatErr->RemovePoint(0);
	
	//graphInvYieldPi0DalitzpPb5023GeVStatErr->Print();
	
	//graphInvYieldPi0DalitzpPb5023GeVSystErr->Print();
	
	
	if (! graphInvYieldPi0DalitzpPb5023GeVComplErr ){
	    
	   cout<<"graphInvYieldPi0DalitzpPb5023GeVComplErr"<<endl;
	   return ;
	} 
	
	if ( !graphInvYieldPi0DalitzpPb5023GeVSystErr ){
	  
	  cout<<"graphInvYieldPi0DalitzpPb5023GeVSystErr"<<endl;
	  
	  return;
	  
	}
	
	if( !graphInvYieldPi0DalitzpPb5023GeVStatErr ) {
	  
	  cout<<"graphInvYieldPi0DalitzpPb5023GeVStatErr"<<endl;
	  return;
	}
	
	cout<<"Llego 2"<<endl;
	
	
	
	
	
	
	
	TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr    = (TGraphAsymmErrors*)graphInvYieldPi0PCMpPb5023GeVSystErr->Clone();
 	TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVXShiftedStatErr    = (TGraphAsymmErrors*)graphInvYieldPi0PCMpPb5023GeVStatErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr   = (TGraphAsymmErrors*)graphInvYieldPi0PCMpPb5023GeVComplErr->Clone();
	
	
	cout<<"Llego 3"<<endl;
	
	
	
	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVXShiftedStatErr  = (TGraphAsymmErrors*) graphInvYieldPi0DalitzpPb5023GeVStatErr->Clone();	
	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr  = (TGraphAsymmErrors*) graphInvYieldPi0DalitzpPb5023GeVSystErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr = (TGraphAsymmErrors*) graphInvYieldPi0DalitzpPb5023GeVComplErr->Clone();
	
	
	
	
	cout<<"Llego 4"<<endl;

	

	
	
	
	cout<<"*******************************BinXshift Dalitz**************************************"<<endl;
	
	
	/*Double_t paramPi0pPb5023GeV[10];
	
	ReturnParameterSetFittingPbPb("8000",paramPi0pPb5023GeV);
	
	
	TF1* fitTsallisPi0pPb5023GeVPtMult = FitObject("rad","tmpt","Pi0");
	
	SetParametersLimitsForFit (fitTsallisPi0pPb5023GeVPtMult, 5, paramPi0pPb5023GeV);
	
	
	 
	Double_t Param[5] = { 0.97524, 2.64925, 0.155715, 3.02966, 5.9877};
   
        fitTsallisPi0pPb5023GeVPtMult->SetParameters(Param);*/
	
	
	Double_t parametersBylinkinpPb5023GeV[5]   =	{0.26,0.5,1.6,0.6,1.8};
	
	TF1 *fitBylinkinPi0pPb5023GeV = FitObject("tcm","Bylinkin","Pi0");
	fitBylinkinPi0pPb5023GeV->SetParameters(parametersBylinkinpPb5023GeV[0],parametersBylinkinpPb5023GeV[1],parametersBylinkinpPb5023GeV[2],parametersBylinkinpPb5023GeV[3],parametersBylinkinpPb5023GeV[4]);
	
	
	
	Double_t paramPi0pPb5023GeV[3]  = {5.0,8.0,0.15};
	TF1* fitTsallisPi0pPb5023GeVPtMult = FitObject("tmpt","tmpt","Pi0");
	fitTsallisPi0pPb5023GeVPtMult->SetParameters(paramPi0pPb5023GeV[0],paramPi0pPb5023GeV[1], paramPi0pPb5023GeV[2]) ; // standard parameter optimize if necessary
	
		
	
	
	
	//graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr  = ApplyXshift(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr ,fitBylinkinPi0pPb5023GeV,"Pi0");
	//graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr, graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr, fitBylinkinPi0pPb5023GeV, 0, graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->GetN());
	//graphInvYieldPi0DalitzpPb5023GeVXShiftedStatErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr, graphInvYieldPi0DalitzpPb5023GeVXShiftedStatErr, fitBylinkinPi0pPb5023GeV, 0, graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->GetN());
	
		
		
	graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr  = ApplyXshift(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr ,fitTsallisPi0pPb5023GeVPtMult,"Pi0");
	graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr, graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr, fitTsallisPi0pPb5023GeVPtMult, 0, graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->GetN());
	graphInvYieldPi0DalitzpPb5023GeVXShiftedStatErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr, graphInvYieldPi0DalitzpPb5023GeVXShiftedStatErr, fitTsallisPi0pPb5023GeVPtMult, 0, graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->GetN());

	
	cout<<"**********************************BinXshift PCM********************************************"<<endl;
	
	
	
	graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr  = ApplyXshift(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr ,fitTsallisPi0pPb5023GeVPtMult,"Pi0");
	graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr, graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr, fitTsallisPi0pPb5023GeVPtMult, 0, graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->GetN());
	graphInvYieldPi0PCMpPb5023GeVXShiftedStatErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr, graphInvYieldPi0PCMpPb5023GeVXShiftedStatErr, fitTsallisPi0pPb5023GeVPtMult, 0, graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->GetN());
	
	
	
	
	//graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr  = ApplyXshift(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr ,fitBylinkinPi0pPb5023GeV,"Pi0");
	//graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr, graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr, fitBylinkinPi0pPb5023GeV, 0, graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->GetN());
	//graphInvYieldPi0PCMpPb5023GeVXShiftedStatErr   = ApplyXshiftIndividualSpectra(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr, graphInvYieldPi0PCMpPb5023GeVXShiftedStatErr, fitBylinkinPi0pPb5023GeV, 0, graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->GetN());
	
	
	
	


	//////////////////////////////////////////Convert to InvYield  7 TeV////////////////////////////////////////////
	
	
	TGraphAsymmErrors* graphInvYieldPi07TeVStatSystErr         = (TGraphAsymmErrors*) graphInvCrossSectionPi07TeVStatSystErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi07TeVStatErr             = (TGraphAsymmErrors*) graphInvCrossSectionPi07TeVStatErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi07TeVSystErr             = (TGraphAsymmErrors*) graphInvCrossSectionPi07TeVSystErr->Clone();

    

	///////////////////////////////////////////Convert to InvYield 2.76 TeV//////////////////////////////////////////
	
	TGraphAsymmErrors* graphInvYieldPi02760GeVStatSystErr      =  (TGraphAsymmErrors*)graphInvCrossSectionPi02760GeVStatSystErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi02760GeVStatErr          =  (TGraphAsymmErrors*)graphInvCrossSectionPi02760GeVStatErr->Clone();
	TGraphAsymmErrors* graphInvYieldPi02760GeVSystErr          =  (TGraphAsymmErrors*)graphInvCrossSectionPi02760GeVSystErr->Clone();


	
	
	
	/////////////////////////////////////////////////Rebinning ////////////////////////////////////////////////////
	
	
	Float_t minPt = 0.4;
	Float_t maxPt = graphInvYieldPi02760GeVStatSystErr->GetXaxis()->GetBinUpEdge(graphInvYieldPi02760GeVStatSystErr->GetXaxis()->GetNbins());
	
	
	
	
	
	TGraphErrors* graphInvYieldPi02760GeVBinDalitzStatSystErr;
	TGraphErrors* graphInvYieldPi02760GeVBinPCMStatSystErr;
	
	
	TGraphAsymmErrors* graphAInvYieldPi02760GeVBinDalitzStatSystErr;
	
	
	TGraphAsymmErrors* graphAInvYieldPi02760GeVBinPCMStatSystErr;

	TGraphAsymmErrors* graphInvYieldPi02760GeVBinDalitzSystErr;

	TGraphAsymmErrors* graphInvYieldPi02760GeVBinPCMSystErr;
	
	cout<<"Rebinning graphInvYieldPi02760GeVStatSystErr to Dalitz"<<endl;
	
	TF1* fitPi02760GeVBinDalitzStatSystErr = RebinWithFitToTGraph(graphInvYieldPi02760GeVStatSystErr,&graphAInvYieldPi02760GeVBinDalitzStatSystErr,graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr,fitType.Data(),minPt,maxPt,Parameters2760GeV);
	
	cout<<"Rebinning graphInvYieldPi02760GeVStatSystErr to PCM"<<endl;
	
	TF1* fitPi02760GeVBinPCMStatSystErr    = RebinWithFitToTGraph(graphInvYieldPi02760GeVStatSystErr,&graphAInvYieldPi02760GeVBinPCMStatSystErr,graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr,fitType.Data(),minPt,maxPt,Parameters2760GeV);
		
	cout<<"Rebinning graphInvYieldPi02760GeVSystErr to Dalitz"<<endl;
	
	TF1* fitPi02760GeVBinDalitzSystErr     = RebinWithFitToTGraph(graphInvYieldPi02760GeVSystErr,&graphInvYieldPi02760GeVBinDalitzSystErr,graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr,fitType.Data(),minPt,maxPt,Parameters2760GeV);
	
	cout<<"Rebinning graphInvYieldPi02760GeVSystErr to PCM"<<endl;
		
	TF1* fitPi02760GeVBinPCMSystErr        = RebinWithFitToTGraph(graphInvYieldPi02760GeVSystErr,&graphInvYieldPi02760GeVBinPCMSystErr,graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr,fitType.Data(),minPt,maxPt,Parameters2760GeV);
	
	
	
	
		
	
	TGraphAsymmErrors* TGraphAymmmRatioToFitPi0PP2760GeV = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphInvYieldPi02760GeVStatSystErr,fitPi02760GeVBinDalitzStatSystErr);
	
	
	
	TCanvas* canvasRatioToFitPi0PP2760GeV = new TCanvas("canvasRatioToFitPi0PP2760GeV","Ratio between PileUp Correction And No PileUp correction spectra",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioToFitPi0PP2760GeV,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioToFitPi0PP2760GeV = new TH2F("histo2DRatioToFitPi0PP2760GeV","histo2DRatioToFitPi0PP2760GeV",1000,0.,15.,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioToFitPi0PP2760GeV, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioToFitPi0PP2760GeV->GetYaxis()->SetRangeUser(0.8,1.2);
	histo2DRatioToFitPi0PP2760GeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(TGraphAymmmRatioToFitPi0PP2760GeV,20,1, kRed , kRed);
	//HistoRatioToFitPi0PP2760GeV->Draw("E1psame");
	//HistoRatioToFitPi0PP2760GeV->Draw("same");
	TGraphAymmmRatioToFitPi0PP2760GeV->Draw("p,same,e");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioToFitPi0PP2760GeV = new TLatex(0.16,0.9,"pp minimum bias #sqrt{#it{s}} = 2.76 TeV");
	SetStyleTLatex( textRatioToFitPi0PP2760GeV,0.06,4); 
	textRatioToFitPi0PP2760GeV->Draw();
	
	TLatex *textRatioToFitPi0PP2760GeVFitFuncName = new TLatex(0.16,0.8,Form("%s Fit",fitNameLabel.Data()));
	SetStyleTLatex( textRatioToFitPi0PP2760GeVFitFuncName,0.06,4); 
	textRatioToFitPi0PP2760GeVFitFuncName->Draw();
	
	
	
	
	
	canvasRatioToFitPi0PP2760GeV->SaveAs(Form("%s/RatiotoFit_Pi0PP2760.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	
	
	minPt = 0.30;
	maxPt = graphInvYieldPi07TeVStatSystErr->GetXaxis()->GetBinUpEdge(graphInvYieldPi07TeVStatSystErr->GetXaxis()->GetNbins());
	

	TGraphErrors* graphInvYieldPi07TeVBinDalitzStatSystErr;
	TGraphErrors* graphInvYieldPi07TeVBinPCMStatSystErr;
	
	TGraphAsymmErrors* graphAInvYieldPi07TeVBinDalitzStatSystErr;
	TGraphAsymmErrors* graphAInvYieldPi07TeVBinPCMStatSystErr;
	
	
	
	TGraphAsymmErrors* graphInvYieldPi07TeVBinDalitzSystErr;
	TGraphAsymmErrors* graphInvYieldPi07TeVBinPCMSystErr;
	
	
	
	cout<<"Rebinning graphInvYieldPi07TeVStatSystErr to Dalitz"<<endl;
	
	
	TF1* fitPi07TeVBinDalitzStatSystErr = RebinWithFitToTGraph(graphInvYieldPi07TeVStatSystErr,&graphAInvYieldPi07TeVBinDalitzStatSystErr,graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr, fitType.Data(),minPt,maxPt,Parameters7TeV);
	
	
	cout<<"Rebinning graphInvYieldPi07TeVStatSystErr to PCM"<<endl;
	
	
	TF1* fitPi07TeVBinPCMStatSystErr    = RebinWithFitToTGraph(graphInvYieldPi07TeVStatSystErr,&graphAInvYieldPi07TeVBinPCMStatSystErr,   graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr,fitType.Data(),minPt,maxPt,Parameters7TeV);
	
	
	cout<<"Rebinning graphInvYieldPi07TeVStatErr to Dalitz"<<endl;
	
		
	cout<<"Rebinning graphInvYieldPi07TeVSystErr to Dalitz"<<endl;
	
	
	TF1* fitPi07TeVBinDalitzSystErr     = RebinWithFitToTGraph(graphInvYieldPi07TeVSystErr,&graphInvYieldPi07TeVBinDalitzSystErr,graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr,fitType.Data(),minPt,maxPt,Parameters7TeV);
	
	cout<<"Rebinning graphInvYieldPi07TeVSystErr to PCM"<<endl;
	
	
	TF1* fitPi07TeVBinPCMSystErr        = RebinWithFitToTGraph(graphInvYieldPi07TeVSystErr,&graphInvYieldPi07TeVBinPCMSystErr,   graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr,fitType.Data(),minPt,maxPt,Parameters7TeV);
	
		
	
	
	cout<<"Start to computing"<<endl;
	
	
	TGraphAsymmErrors* TGraphAymmmRatioToFitPi0PP7TeV = (TGraphAsymmErrors*) CalculateGraphErrRatioToFit(graphInvYieldPi07TeVStatSystErr,fitPi07TeVBinDalitzStatSystErr);
	
	
	
	TCanvas* canvasRatioToFitPi0PP7TeV = new TCanvas("canvasRatioToFitPi0PP7TeV","Ratio between PileUp Correction And No PileUp correction spectra",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioToFitPi0PP7TeV,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioToFitPi0PP7TeV = new TH2F("histo2DRatioToFitPi0PP7TeV","histo2DRatioToFitPi0PP7TeV",1000,0.,15.,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioToFitPi0PP7TeV, "#it{p}_{T} (GeV/#it{c})","Data/Fit", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioToFitPi0PP7TeV->GetYaxis()->SetRangeUser(0.8,1.2);
	histo2DRatioToFitPi0PP7TeV->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(TGraphAymmmRatioToFitPi0PP7TeV,20,1, kRed , kRed);
	//HistoRatioToFitPi0PP7TeV->Draw("E1psame");
	//HistoRatioToFitPi0PP7TeV->Draw("same");
	TGraphAymmmRatioToFitPi0PP7TeV->Draw("p,same,e");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioToFitPi0PP7TeV = new TLatex(0.16,0.9,"pp minimum bias #sqrt{#it{s}} = 7 TeV");
	SetStyleTLatex( textRatioToFitPi0PP7TeV,0.06,4); 
	textRatioToFitPi0PP7TeV->Draw();
	
	TLatex *textRatioToFitPi0PP7TeVFitFuncName = new TLatex(0.16,0.8,Form("%s Fit",fitNameLabel.Data()));
	SetStyleTLatex( textRatioToFitPi0PP7TeVFitFuncName,0.06,4); 
	textRatioToFitPi0PP7TeVFitFuncName->Draw();
	
	
	
	
	canvasRatioToFitPi0PP7TeV->SaveAs(Form("%s/RatiotoFit_Pi0PP7TeV.%s",outputDir.Data(),suffix.Data()));
	
	
		
	
	
	
	graphInvYieldPi02760GeVBinDalitzStatSystErr 		= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi02760GeVBinDalitzStatSystErr);
	graphInvYieldPi02760GeVBinPCMStatSystErr    		= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi02760GeVBinPCMStatSystErr);

	
	
	graphInvYieldPi07TeVBinDalitzStatSystErr 		= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi07TeVBinDalitzStatSystErr);
	graphInvYieldPi07TeVBinPCMStatSystErr    		= ConvertTGraphAsymmErrorstoTGraphErrors(graphAInvYieldPi07TeVBinPCMStatSystErr);
	
	
	
	TGraphErrors* graphErrosInterPolation5023GeVBinDalitz = GetInterpolSpectrum2D(graphInvYieldPi02760GeVBinDalitzStatSystErr, graphInvYieldPi07TeVBinDalitzStatSystErr, 2760,7000,5023, "Dalitz");
	TGraphErrors* graphErrosInterPolation5023GeVBinPCM    = GetInterpolSpectrum2D(graphInvYieldPi02760GeVBinPCMStatSystErr,    graphInvYieldPi07TeVBinPCMStatSystErr,    2760,7000,5023, "PCM");
	
	TString namePtvsSqrtsPlotPCM    = 	Form("%s/Pt_vs_Sqrts_PCM.%s",outputDir.Data(),suffix.Data());	
	TString namePtvsSqrtsPlotDalitz = 	Form("%s/Pt_vs_Sqrts_Dalitz.%s",outputDir.Data(),suffix.Data());	
	
	
	PlotInterpolationPtBins(graphPtvsSqrtsPCM,gPtvsEnergiesPCM,fPowerlawPCM,graphErrosInterPolation5023GeVBinPCM,6,5,namePtvsSqrtsPlotPCM);
       	PlotInterpolationPtBins(graphPtvsSqrtsDalitz,gPtvsEnergiesDalitz,fPowerlawDalitz,graphErrosInterPolation5023GeVBinDalitz,5,4,namePtvsSqrtsPlotDalitz);
	//void PlotAlphavsPt(TGraphErrors* gAlpha,TString method, TString thesisPlotLabel, TString namePlot){
	  
	  cout<<"Aqui debe antes "<<thesisPlotLabel.Data()<<endl;

	PlotAlphavsPt(graphAlphaPCM,    "PCM",    thesisPlotLabel.Data(),  Form("%s/Alpha_vs_Pt_PCM.%s",   outputDir.Data(),suffix.Data()));
	PlotAlphavsPt(graphAlphaDalitz, "Dalitz", thesisPlotLabel.Data(),  Form("%s/Alpha_vs_Pt_Dalitz.%s",outputDir.Data(),suffix.Data()));

	
	
	cout<<"Statistics for Dalitz"<<endl;  //The statistical errors are taken from the interpolated pp reference 
	
	Int_t     DalitznPoints    =   graphErrosInterPolation5023GeVBinDalitz->GetN();
	Double_t *DalitzxValue     =   graphErrosInterPolation5023GeVBinDalitz->GetX();
	Double_t *DalitzxStatErr   =   graphInvYieldPi07TeVBinDalitzStatSystErr->GetEX(); //X errors are taken from the rebbined 7 spectrum
	Double_t *DalitzyValue     =   graphErrosInterPolation5023GeVBinDalitz->GetY();
	Double_t *DalitzyStatErr   =   graphErrosInterPolation5023GeVBinDalitz->GetEY();
	
	
	graphErrosInterPolation5023GeVBinDalitzStatErr = new TGraphAsymmErrors(DalitznPoints,DalitzxValue,DalitzyValue,DalitzxStatErr,DalitzxStatErr,DalitzyStatErr,DalitzyStatErr);
	
	
	cout<<"Systematics for Dalitz"<<endl;
	graphErrosInterPolation5023GeVBinDalitzSystErr = CalculateSystErrors(graphErrosInterPolation5023GeVBinDalitz,graphInvYieldPi02760GeVBinDalitzSystErr,graphInvYieldPi07TeVBinDalitzSystErr);
	
	
	
	
	
	/////////////////////////////////////////////////////////////////////
	
	
		
	cout<<"Statistics for PCM"<<endl;  //The statistical errors are taken from the interpolated pp reference 
	
	Int_t     PCMnPoints    =   graphErrosInterPolation5023GeVBinPCM->GetN();
	
	Double_t *PCMxValue     =   graphErrosInterPolation5023GeVBinPCM->GetX();
	Double_t *PCMxStatErr   =   graphInvYieldPi07TeVBinPCMStatSystErr->GetEX(); //X errors are taken from the rebbined 7 spectrum
	
	Double_t *PCMyValue     =   graphErrosInterPolation5023GeVBinPCM->GetY();
	Double_t *PCMyStatErr   =   graphErrosInterPolation5023GeVBinPCM->GetEY();
	
	
	graphErrosInterPolation5023GeVBinPCMStatErr = new TGraphAsymmErrors(PCMnPoints,PCMxValue,PCMyValue,PCMxStatErr,PCMxStatErr,PCMyStatErr,PCMyStatErr);
	
	
	
	
	cout<<"Systematics for PCM"<<endl;
	graphErrosInterPolation5023GeVBinPCMSystErr    = CalculateSystErrors(graphErrosInterPolation5023GeVBinPCM,graphInvYieldPi02760GeVBinPCMSystErr,graphInvYieldPi07TeVBinPCMSystErr);
	
	cout<<"Binning Interpolation Dalitz"<<endl;
	
	graphErrosInterPolation5023GeVBinPCMStatErr->Print();
		
	cout<<"Binning Interpolation PCM   "<<endl;
	
	graphErrosInterPolation5023GeVBinPCMStatErr->Print();
	

	
	
	
	if( resultsType.CompareTo("PCM") == 0 || resultsType.CompareTo("PCMPileUpCorrection") == 0 ) {
	  
	graphErrosInterPolation5023GeVBinDalitzSystWOMatErr  =  CancelOutMaterialError(graphErrosInterPolation5023GeVBinDalitzSystErr,"PCMDalitz");
	graphInvYieldPi0DalitzpPb5023GeVXShiftedSystWOMatErr =  CancelOutMaterialError(graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr,"Dalitz");
	
	
	graphErrosInterPolation5023GeVBinPCMSystWOMatErr     = CancelOutMaterialError(graphErrosInterPolation5023GeVBinPCMSystErr, "PCMPCM");
	graphInvYieldPi0PCMpPb5023GeVXShiftedSystWOMatErr    = CancelOutMaterialError(graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr,"PCMPCM");
	
	
	CalcRpPb(graphErrosInterPolation5023GeVBinDalitzSystWOMatErr,graphErrosInterPolation5023GeVBinDalitzStatErr,graphInvYieldPi0DalitzpPb5023GeVXShiftedSystWOMatErr, graphInvYieldPi0DalitzpPb5023GeVXShiftedStatErr,&graphRpPbDalitzSystErr,&graphRpPbDalitzStatErr);
	CalcRpPb(graphErrosInterPolation5023GeVBinPCMSystWOMatErr,graphErrosInterPolation5023GeVBinPCMStatErr,graphInvYieldPi0PCMpPb5023GeVXShiftedSystWOMatErr, graphInvYieldPi0PCMpPb5023GeVXShiftedStatErr,&graphRpPbPCMSystErr,&graphRpPbPCMStatErr);
	
	} else {
		
	CalcRpPb(graphErrosInterPolation5023GeVBinDalitzSystErr,graphErrosInterPolation5023GeVBinDalitzStatErr,graphInvYieldPi0DalitzpPb5023GeVXShiftedSystErr, graphInvYieldPi0DalitzpPb5023GeVXShiftedStatErr,&graphRpPbDalitzSystErr,&graphRpPbDalitzStatErr);
	CalcRpPb(graphErrosInterPolation5023GeVBinPCMSystErr,graphErrosInterPolation5023GeVBinPCMStatErr,graphInvYieldPi0PCMpPb5023GeVXShiftedSystErr, graphInvYieldPi0PCMpPb5023GeVXShiftedStatErr,&graphRpPbPCMSystErr,&graphRpPbPCMStatErr);
	
	}
	
	
	
	/////////////Drawing InvYield pp@2.76 GeV//////////////
	
	
	TCanvas* cInvYield2760GeVBinDalitz = new TCanvas("cInvYield2760GeVBinDalitz","Invariant Yield pp@2.76 TeV",200,10,700,550);
    	DrawGammaCanvasSettings( cInvYield2760GeVBinDalitz,  0.15, 0.02, 0.03, 0.1);
	cInvYield2760GeVBinDalitz->SetLogx();
	cInvYield2760GeVBinDalitz->SetLogy();
	
	TPad* padComparisonInvYield2760GeVBinDalitz = new TPad("padComparisonInvYield2760GeVBinDalitz", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYield2760GeVBinDalitz, 0.15, 0.02, 0.03, 0.1);
	padComparisonInvYield2760GeVBinDalitz->Draw();
	
	
	TH2F * histoInvYield2760GeVBinDalitz;
	histoInvYield2760GeVBinDalitz = new TH2F("histoInvYield2760GeVBinDalitz","histoInvYield2760GeVBinDalitz",1000,0.23,30.,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histoInvYield2760GeVBinDalitz, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYield2760GeVBinDalitz->DrawCopy(); 
	
	graphInvYieldPi02760GeVStatSystErr->SetFillColor(0);
	graphAInvYieldPi02760GeVBinDalitzStatSystErr->SetFillColor(0);

	
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi02760GeVStatSystErr, 20,1,kBlue,kBlue);
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi02760GeVBinDalitzStatSystErr, 20,1,kRed,kRed);
	
	
	TLegend* legendInvYields2760GeVBinDalitz = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields2760GeVBinDalitz->SetFillColor(0);
	legendInvYields2760GeVBinDalitz->SetLineColor(0);
	legendInvYields2760GeVBinDalitz->SetNColumns(1);
	legendInvYields2760GeVBinDalitz->SetTextSize(0.03);
	legendInvYields2760GeVBinDalitz->AddEntry(graphInvYieldPi02760GeVStatSystErr,"pp at #sqrt{#it{s}} = 2.76 TeV","pef");
	legendInvYields2760GeVBinDalitz->AddEntry(graphAInvYieldPi02760GeVBinDalitzStatSystErr,"pp at #sqrt{#it{s}} = 2.76 TeV rebinned","pef");
	legendInvYields2760GeVBinDalitz->AddEntry(fitPi02760GeVBinDalitzStatSystErr,Form("%s fit",fitNameLabel.Data()),"pef");
	
	
	
	graphAInvYieldPi02760GeVBinDalitzStatSystErr->Draw("p,same,e1");
	graphInvYieldPi02760GeVStatSystErr->Draw("p,same,e1");
	fitPi02760GeVBinDalitzStatSystErr->Draw("same");
	legendInvYields2760GeVBinDalitz->Draw("sames");
	
	
	cInvYield2760GeVBinDalitz->Update();
	cInvYield2760GeVBinDalitz->SaveAs(Form("%s/InvYield_PP_BinDalitz_2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	TCanvas* cInvYield7TeVBinDalitz = new TCanvas("cInvYield7TeVBinDalitz","Invariant Yield pp@2.76 TeV",200,10,700,550);
    	DrawGammaCanvasSettings( cInvYield7TeVBinDalitz,  0.15, 0.02, 0.03, 0.1);
	cInvYield7TeVBinDalitz->SetLogx();
	cInvYield7TeVBinDalitz->SetLogy();
	
	TPad* padComparisonInvYield7TeVBinDalitz = new TPad("padComparisonInvYield7TeVBinDalitz", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYield7TeVBinDalitz, 0.15, 0.02, 0.03, 0.1);
	padComparisonInvYield7TeVBinDalitz->Draw();
	
	
	TH2F * histoInvYield7TeVBinDalitz;
	histoInvYield7TeVBinDalitz = new TH2F("histoInvYield7TeVBinDalitz","histoInvYield7TeVBinDalitz",1000,0.23,30.,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histoInvYield7TeVBinDalitz, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYield7TeVBinDalitz->DrawCopy(); 
	

	graphInvYieldPi07TeVStatSystErr->SetFillColor(0);
	graphAInvYieldPi07TeVBinDalitzStatSystErr->SetFillColor(0);
	
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi07TeVBinDalitzStatSystErr, 20,1,kBlue,kBlue);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi07TeVStatSystErr, 20,1,kRed,kRed);
	
	
	TLegend* legendInvYields7TeVBinDalitz = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields7TeVBinDalitz->SetFillColor(0);
	legendInvYields7TeVBinDalitz->SetLineColor(0);
	legendInvYields7TeVBinDalitz->SetNColumns(1);
	legendInvYields7TeVBinDalitz->SetTextSize(0.03);
	legendInvYields7TeVBinDalitz->AddEntry(graphInvYieldPi07TeVStatSystErr,"pp at #sqrt{#it{s}} = 7 TeV","pef");
	legendInvYields7TeVBinDalitz->AddEntry(graphAInvYieldPi07TeVBinDalitzStatSystErr,"pp at #sqrt{#it{s}} = 7 TeV rebinned","pef");
	legendInvYields7TeVBinDalitz->AddEntry(fitPi07TeVBinDalitzStatSystErr,Form("%s fit",fitNameLabel.Data()),"pef");
	
	
	
	
	
	
	
	graphAInvYieldPi07TeVBinDalitzStatSystErr->Draw("p,same,e1");
	graphInvYieldPi07TeVStatSystErr->Draw("p,same,e1");
	fitPi07TeVBinDalitzStatSystErr->Draw("same");
	legendInvYields7TeVBinDalitz->Draw();
	
	
	cInvYield7TeVBinDalitz->Update();
	cInvYield7TeVBinDalitz->SaveAs(Form("%s/InvYield_PP_BinDalitz_7TeV.%s",outputDir.Data(),suffix.Data()));
	
	
	
		
	
	
	TCanvas* cInvYield2760GeVBinPCM = new TCanvas("cInvYield2760GeVBinPCM","Invariant Yield pp@2.76 TeV",200,10,700,550);
    	DrawGammaCanvasSettings( cInvYield2760GeVBinPCM,  0.15, 0.02, 0.03, 0.1);
	cInvYield2760GeVBinPCM->SetLogx();
	cInvYield2760GeVBinPCM->SetLogy();
	
	TPad* padComparisonInvYield2760GeVBinPCM = new TPad("padComparisonInvYield2760GeVBinPCM", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYield2760GeVBinPCM, 0.15, 0.02, 0.03, 0.1);
	padComparisonInvYield2760GeVBinPCM->Draw();
	
	
	TH2F * histoInvYield2760GeVBinPCM;
	histoInvYield2760GeVBinPCM = new TH2F("histoInvYield2760GeVBinPCM","histoInvYield2760GeVBinPCM",1000,0.23,30.,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histoInvYield2760GeVBinPCM, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYield2760GeVBinPCM->DrawCopy(); 
	
        graphInvYieldPi02760GeVStatSystErr->SetFillColor(0);
	graphAInvYieldPi02760GeVBinPCMStatSystErr->SetFillColor(0);
	
	
	
	
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi02760GeVBinPCMStatSystErr, 20,1,kBlue,kBlue);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi02760GeVStatSystErr, 20,1,kRed,kRed);
	
	TLegend* legendInvYields2760GeVBinPCM = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields2760GeVBinPCM->SetFillColor(0);
	legendInvYields2760GeVBinPCM->SetLineColor(0);
	legendInvYields2760GeVBinPCM->SetNColumns(1);
	legendInvYields2760GeVBinPCM->SetTextSize(0.03);
	legendInvYields2760GeVBinPCM->AddEntry(graphInvYieldPi02760GeVStatSystErr,"pp at #sqrt{#it{s}} = 2.76 TeV","pef");
	legendInvYields2760GeVBinPCM->AddEntry(graphAInvYieldPi02760GeVBinPCMStatSystErr,"pp at #sqrt{#it{s}} = 2.76 TeV rebinned","pef");
	legendInvYields2760GeVBinPCM->AddEntry(fitPi02760GeVBinPCMStatSystErr,Form("%s fit",fitNameLabel.Data()),"pef");
	
	
	
	
	
	
	
	graphAInvYieldPi02760GeVBinPCMStatSystErr->Draw("p,same,e1");
	graphInvYieldPi02760GeVStatSystErr->Draw("p,same,e1");
	fitPi02760GeVBinPCMStatSystErr->Draw("p,same");
	legendInvYields2760GeVBinPCM->Draw("same");
	
	
	cInvYield2760GeVBinPCM->Update();
	cInvYield2760GeVBinPCM->SaveAs(Form("%s/InvYield_PP_BinPCM_2760GeV.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	TCanvas* cInvYield7TeVBinPCM = new TCanvas("cInvYield7TeVBinPCM","Invariant Yield pp@2.76 TeV",200,10,700,550);
    	DrawGammaCanvasSettings( cInvYield7TeVBinPCM,  0.15, 0.02, 0.03, 0.1);
	cInvYield7TeVBinPCM->SetLogx();
	cInvYield7TeVBinPCM->SetLogy();
	
	TPad* padComparisonInvYield7TeVBinPCM = new TPad("padComparisonInvYield7TeVBinPCM", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padComparisonInvYield7TeVBinPCM, 0.15, 0.02, 0.03, 0.1);
	padComparisonInvYield7TeVBinPCM->Draw();
	
	
	TH2F * histoInvYield7TeVBinPCM;
	histoInvYield7TeVBinPCM = new TH2F("histoInvYield7TeVBinPCM","histoInvYield7TeVBinPCM",1000,0.23,30.,1000,2e0,10e12);
	SetStyleHistoTH2ForGraphs(histoInvYield7TeVBinPCM, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYield7TeVBinPCM->DrawCopy(); 
	
        graphInvYieldPi07TeVStatSystErr->SetFillColor(0);
	graphAInvYieldPi07TeVBinPCMStatSystErr->SetFillColor(0);
	
	
	
	DrawGammaSetMarkerTGraphAsym(graphAInvYieldPi07TeVBinPCMStatSystErr, 20,1,kBlue,kBlue);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi07TeVStatSystErr, 20,1,kRed,kRed);
	
	
	TLegend* legendInvYields7TeVBinPCM = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields7TeVBinPCM->SetFillColor(0);
	legendInvYields7TeVBinPCM->SetLineColor(0);
	legendInvYields7TeVBinPCM->SetNColumns(1);
	legendInvYields7TeVBinPCM->SetTextSize(0.03);
	legendInvYields7TeVBinPCM->AddEntry(graphInvYieldPi07TeVStatSystErr,"pp at #sqrt{#it{s}} = 7 TeV","pef");
	legendInvYields7TeVBinPCM->AddEntry(graphAInvYieldPi07TeVBinPCMStatSystErr,"pp at #sqrt{#it{s}} = 7 TeV rebinned","pef");
	legendInvYields7TeVBinPCM->AddEntry(fitPi07TeVBinPCMStatSystErr,Form("%s fit",fitNameLabel.Data()),"pef");
	
	
	
	
	graphAInvYieldPi07TeVBinPCMStatSystErr->Draw("p,same,e1");
	graphInvYieldPi07TeVStatSystErr->Draw("p,same,e1");
	fitPi07TeVBinPCMStatSystErr->Draw("same");
	legendInvYields7TeVBinPCM->Draw("same");
	
	
	cInvYield7TeVBinPCM->Update();
	cInvYield7TeVBinPCM->SaveAs(Form("%s/InvYield_PP_BinPCM_7TeV.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	////////////////////////////Scaling graphsErrors to be painted////////////////////
	
	
	
	TGraphErrors* graphInvYieldPi07TeVBinDalitzStatSystErrScaled   		   = new TGraphErrors( *graphInvYieldPi07TeVBinDalitzStatSystErr );
	TGraphErrors* graphInvYieldPi02760GeVBinDalitzStatSystErrScaled 	   = new TGraphErrors( *graphInvYieldPi02760GeVBinDalitzStatSystErr);
	TGraphErrors* graphErrosInterPolation5023GeVBinDalitzScaled     	   = new TGraphErrors( *graphErrosInterPolation5023GeVBinDalitz );
	TGraphErrors* graphErrosInterPolation5023GeVBinDalitzTpPb                  = new TGraphErrors( *graphErrosInterPolation5023GeVBinDalitz );
	//TGraphAsymmErrors* graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErrC        = new TGraphAsymmErrors( *graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr );
	
	
	
	
	graphInvYieldPi07TeVBinDalitzStatSystErrScaled        = ScaleGraph(graphInvYieldPi07TeVBinDalitzStatSystErrScaled,4);
	graphErrosInterPolation5023GeVBinDalitzScaled         = ScaleGraph(graphErrosInterPolation5023GeVBinDalitzScaled,2);
	graphInvYieldPi02760GeVBinDalitzStatSystErrScaled     = ScaleGraph(graphInvYieldPi02760GeVBinDalitzStatSystErrScaled,1);
	graphErrosInterPolation5023GeVBinDalitzTpPb           = ScaleGraph(graphErrosInterPolation5023GeVBinDalitzTpPb,fTpPb);
	//graphInvCrossSectionPi0DalitzpPb5023GeVXShiftedComplErr   = ScaleGraph(graphInvCrossSectionPi0DalitzpPb5023GeVXShiftedComplErr,recalcBarn*xSectionpPb5023GeVINEL);
	
	
	
	
	TGraphErrors* graphInvYieldPi07TeVBinPCMStatSystErrScaled       = new TGraphErrors( *graphInvYieldPi07TeVBinPCMStatSystErr );
	TGraphErrors* graphInvYieldPi02760GeVBinPCMScaled  	        = new TGraphErrors( *graphInvYieldPi02760GeVBinPCMStatSystErr);
	TGraphErrors* graphErrosInterPolation5023GeVBinPCMScaled        = new TGraphErrors( *graphErrosInterPolation5023GeVBinPCM );
	
	TGraphErrors* graphErrosInterPolation5023GeVBinPCMTpPb                  = new TGraphErrors( *graphErrosInterPolation5023GeVBinPCM );
	
	
	graphInvYieldPi07TeVBinPCMStatSystErrScaled        = ScaleGraph(graphInvYieldPi07TeVBinPCMStatSystErrScaled,4);
	graphErrosInterPolation5023GeVBinPCMScaled         = ScaleGraph(graphErrosInterPolation5023GeVBinPCMScaled,2);
	graphInvYieldPi02760GeVBinPCMScaled                = ScaleGraph(graphInvYieldPi02760GeVBinPCMScaled,1);
	graphErrosInterPolation5023GeVBinPCMTpPb           = ScaleGraph(graphErrosInterPolation5023GeVBinPCMTpPb,fTpPb);
	
	
	
	
	
	
        		
													 
	//TCanvas* cInvYieldComparison = new TCanvas("cInvYieldComparison","Comparison Invariant Yields",200,10,700,500);
	TCanvas* cInvYieldComparison = new TCanvas("cInvYieldComparison","Comparison Invariant Yields",550,700);
		
	DrawGammaCanvasSettings( cInvYieldComparison,  0.15, 0.02, 0.03, 0.1);
	cInvYieldComparison->SetLogx();
	cInvYieldComparison->SetLogy();
	
	TPad* padInvYieldComparison = new TPad("padInvYieldComparison", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldComparison, 0.15, 0.02, 0.03, 0.1);
	padInvYieldComparison->Draw();
	
	
	TH2F * histoInvYieldComparison;
	histoInvYieldComparison = new TH2F("histoInvYieldComparison","histoInvYieldComparison",1000,0.23,30.,1000,5e3,10e11);
	SetStyleHistoTH2ForGraphs(histoInvYieldComparison, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldComparison->GetXaxis()->SetRangeUser(0.6,11.0);
	histoInvYieldComparison->GetYaxis()->SetRangeUser(1e4,2e11);
	histoInvYieldComparison->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphInvYieldPi07TeVBinDalitzStatSystErrScaled,     	    20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphErr(graphInvYieldPi02760GeVBinDalitzStatSystErrScaled, 	    20, 1, kGreen+2, kGreen+2);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeVBinDalitzScaled,   20, 1, kBlue+1,  kBlue+1);
	
	
	
	graphInvYieldPi07TeVBinDalitzStatSystErrScaled->Draw("p,same,e1");
	graphInvYieldPi02760GeVBinDalitzStatSystErrScaled->Draw("p,same,e1");
	graphErrosInterPolation5023GeVBinDalitzScaled->Draw("p,same,e1");
	
	
	
	
	//TLatex *textRpPbChargedPionsKaonsProtonsDalitz = new TLatex(0.3,0.35,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	//SetStyleTLatex( textRpPbChargedPionsKaonsProtonsDalitz,0.035,4); 
	//textRpPbChargedPionsKaonsProtonsDalitz->Draw();
	
	
	//TLatex *LabelpPbChargedPionsKaonsProtonsDalitz = new TLatex(0.3,0.4,"ALICE work in progress");
	//SetStyleTLatex( LabelpPbChargedPionsKaonsProtonsDalitz, 0.035,4);
	//LabelpPbChargedPionsKaonsProtonsDalitz->Draw();
	
	
	TLatex *LabelpPb = new TLatex(0.45,0.9,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPb, 0.03,4);
	LabelpPb->Draw();
	TLatex *labelSpectraPi0LabelPCMpPb = new TLatex(0.45,0.85,"#pi^{0} meson");
	SetStyleTLatex( labelSpectraPi0LabelPCMpPb, 0.03,4);
	labelSpectraPi0LabelPCMpPb->Draw();

	//TLatex *labelSpectraPi0LabelpPb = new TLatex(0.65,0.8,"|y_{#pi^{0},lab}| < 0.8");
	//SetStyleTLatex( labelSpectraPi0LabelpPb, 0.03,4);
	//labelSpectraPi0LabelpPb->Draw();  

	graphInvYieldPi07TeVBinDalitzStatSystErrScaled->SetFillColor(0);
	graphInvYieldPi02760GeVBinDalitzStatSystErrScaled->SetFillColor(0);
	graphErrosInterPolation5023GeVBinDalitzScaled->SetFillColor(0);
	TLegend* legendInvYields = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYields->SetFillColor(0);
	legendInvYields->SetLineColor(0);
	legendInvYields->SetNColumns(1);
	legendInvYields->SetTextSize(0.03);
	legendInvYields->AddEntry(graphInvYieldPi07TeVBinDalitzStatSystErrScaled,"calculated pp #sqrt{#it{s}} = 7 TeV (x4)");
	legendInvYields->AddEntry(graphErrosInterPolation5023GeVBinDalitzScaled,"pp reference #sqrt{#it{s}} = 5.02 TeV (x2)");
	legendInvYields->AddEntry(graphInvYieldPi02760GeVBinDalitzStatSystErrScaled,"calculated pp #sqrt{#it{s}} = 2.76 TeV (x1)");
	
	
	
	legendInvYields->Draw();

		
	
	cout<<"TGraphErrors pp@5.023 TeV"<<endl;
	graphErrosInterPolation5023GeVBinDalitzScaled->Print();
	
	cInvYieldComparison->Update();
	cInvYieldComparison->SaveAs(Form("%s/InvYield_Comparison_PP_Dalitz_2760GeV_5023GeV_7TeV.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	TCanvas* cppAndpPbComarisonDalitz = new TCanvas("cppAndpPbComarisonDalitz","",1150,1000);
		
	//DrawGammaCanvasSettings( cppAndpPbComarisonDalitz,  0.18, 0.02, 0.03, 0.11);
	DrawGammaCanvasSettings( cppAndpPbComarisonDalitz,  0.13, 0.02, 0.03, 0.1);
	
	cppAndpPbComarisonDalitz->SetLogx();
	cppAndpPbComarisonDalitz->SetLogy();
	
	//TPad* padppAndpPbComarisonDalitz = new TPad("padppAndpPbComarisonDalitz", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	//DrawGammaPadSettings( padppAndpPbComarisonDalitz, 0.24, 0.02, 0.03, 0.1);
	//padppAndpPbComarisonDalitz->Draw();
	
	
	TH2F * histoppAndpPbComparisonDalitz;
	histoppAndpPbComparisonDalitz = new TH2F("histoppAndpPbComparisonDalitz","histoppAndpPbComparisonDalitz",1000,0.2,20.,1000,1e-6,1e2);
	SetStyleHistoTH2ForGraphs(histoppAndpPbComparisonDalitz,pTLabel.Data(), invYieldLabel.Data(), 0.03,0.03, 0.03,0.03, 1.3,1.8);
	histoppAndpPbComparisonDalitz->GetXaxis()->SetRangeUser(0.5,12);
	histoppAndpPbComparisonDalitz->GetYaxis()->SetRangeUser(1e-6,5);
	histoppAndpPbComparisonDalitz->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr,20,1.5,  kBlack,   kBlack);
	
	
	graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->Draw("p,e1");
	graphErrosInterPolation5023GeVBinDalitzTpPb->Draw("XL,same");
	
	TLatex *LabelppAndpPbPCM = new TLatex(0.6,0.9,thesisPlotLabel.Data());
	SetStyleTLatex( LabelppAndpPbPCM, 0.03,4);
	LabelppAndpPbPCM->Draw();
	TLatex *labelSpectraPi0DalitzpPb = new TLatex(0.6,0.85,DalitzLabel.Data());
	SetStyleTLatex( labelSpectraPi0DalitzpPb, 0.03,4);
	labelSpectraPi0DalitzpPb->Draw();


	graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->SetFillColor(0);
	graphErrosInterPolation5023GeVBinDalitzTpPb->SetFillColor(0);
	TLegend* legendppAndpPbComarisonDalitz = new TLegend(0.2,0.25,0.4,0.35);
	legendppAndpPbComarisonDalitz->SetFillColor(0);
	legendppAndpPbComarisonDalitz->SetTextFont(42);
	legendppAndpPbComarisonDalitz->SetLineColor(0);
	legendppAndpPbComarisonDalitz->SetNColumns(1);
	legendppAndpPbComarisonDalitz->SetTextSize(0.03);
	legendppAndpPbComarisonDalitz->AddEntry(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr,Form("%s, |#it{y}_{cms}| < 0.8",energyLabel.Data()));
	legendppAndpPbComarisonDalitz->AddEntry(graphErrosInterPolation5023GeVBinDalitzTpPb,"pp reference x #LT #it{T}_{pPb} #GT");
	
	legendppAndpPbComarisonDalitz->Draw();

		
	
	cout<<"TGraphErrors pp@5.023 TeV"<<endl;
	graphErrosInterPolation5023GeVBinDalitzScaled->Print();
	
	cppAndpPbComarisonDalitz->Update();
	cppAndpPbComarisonDalitz->SaveAs(Form("%s/InvYield_Comparison_PP_pPb_Dalitz.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	TCanvas* cppAndpPbComarisonPCM = new TCanvas("cppAndpPbComarisonPCM","Comparison Invariant Yields",550,700);
		
	DrawGammaCanvasSettings( cppAndpPbComarisonPCM,  0.18, 0.02, 0.03, 0.1);
	cppAndpPbComarisonPCM->SetLogx();
	cppAndpPbComarisonPCM->SetLogy();
	
	TPad* padppAndpPbComarisonPCM = new TPad("padppAndpPbComarisonPCM", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padppAndpPbComarisonPCM, 0.24, 0.02, 0.03, 0.1);
	padppAndpPbComarisonPCM->Draw();
	
	
	TH2F * histoppAndpPbComparisonPCM;
	histoppAndpPbComparisonPCM = new TH2F("histoppAndpPbComparisonPCM","histoppAndpPbComparisonPCM",1000,0.001,10.,1000,1e-6,1e2);
	SetStyleHistoTH2ForGraphs(histoppAndpPbComparisonPCM, "#it{p}_{T} (GeV/#it{c})","#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.9);
	histoppAndpPbComparisonPCM->GetXaxis()->SetRangeUser(0.6,10);
	histoppAndpPbComparisonPCM->GetYaxis()->SetRangeUser(1e-6,5);
	histoppAndpPbComparisonPCM->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr,20,1,  kBlack,   kBlack);
	
	
	graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->Draw("p,e1");
	graphErrosInterPolation5023GeVBinPCMTpPb->Draw("XL,same");
	
	TLatex *LabelppAndpPb = new TLatex(0.45,0.9,thesisPlotLabel.Data());
	SetStyleTLatex( LabelppAndpPb, 0.04,4);
	LabelppAndpPb->Draw();
	TLatex *labelSpectraPi0PCMpPb = new TLatex(0.45,0.85,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-}");
	SetStyleTLatex( labelSpectraPi0PCMpPb, 0.04,4);
	labelSpectraPi0PCMpPb->Draw();


	graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->SetFillColor(0);
	graphErrosInterPolation5023GeVBinPCMTpPb->SetFillColor(0);
	TLegend* legendppAndpPbComarisonPCM = new TLegend(0.2,0.15,0.6,0.25);
	legendppAndpPbComarisonPCM->SetFillColor(0);
	legendppAndpPbComarisonPCM->SetLineColor(0);
	legendppAndpPbComarisonPCM->SetNColumns(1);
	legendppAndpPbComarisonPCM->SetTextSize(0.04);
	legendppAndpPbComarisonPCM->AddEntry(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr,"p-Pb, MB, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	legendppAndpPbComarisonPCM->AddEntry(graphErrosInterPolation5023GeVBinPCMTpPb,"pp reference x #LT #it{T}_{pPb} #GT");
	
	legendppAndpPbComarisonPCM->Draw();

		
	
	cout<<"TGraphErrors pp@5.023 TeV"<<endl;
	graphErrosInterPolation5023GeVBinPCMScaled->Print();
	
	cppAndpPbComarisonPCM->Update();
	cppAndpPbComarisonPCM->SaveAs(Form("%s/InvYield_Comparison_PP_pPb_PCM.%s",outputDir.Data(),suffix.Data()));
	
	
		
	
	
	
	
	TCanvas* cInvYieldComparisonPCM = new TCanvas("cInvYieldComparisonPCM","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldComparisonPCM,  0.15, 0.02, 0.03, 0.1);
	cInvYieldComparisonPCM->SetLogx();
	cInvYieldComparisonPCM->SetLogy();
	
	TPad* padInvYieldComparisonPCM = new TPad("padInvYieldComparisonPCM", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldComparisonPCM, 0.15, 0.02, 0.03, 0.1);
	padInvYieldComparisonPCM->Draw();
	
	
	TH2F * histoInvYieldComparisonPCM;
	histoInvYieldComparisonPCM = new TH2F("histoInvYieldComparisonPCM","histoInvYieldComparisonPCM",1000,0.23,30.,1000,2e2,10e11);
	SetStyleHistoTH2ForGraphs(histoInvYieldComparisonPCM, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldComparisonPCM->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphInvYieldPi07TeVBinPCMStatSystErrScaled,     	20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphErr(graphInvYieldPi02760GeVBinPCMScaled, 	20, 1, kGreen+2, kGreen+2);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeVBinPCMScaled,  20, 1, kBlue+1,  kBlue+1);
	
	
	
	
	
	graphInvYieldPi07TeVBinPCMStatSystErrScaled->Draw("p,same,e1");
	graphInvYieldPi02760GeVBinPCMScaled->Draw("p,same,e1");
	graphErrosInterPolation5023GeVBinPCMScaled->Draw("p,same,e1");
	
	
	TLatex *LabelpPbPCM = new TLatex(0.65,0.9,"ALICE work in progress");
	SetStyleTLatex( LabelpPbPCM, 0.03,4);
	LabelpPbPCM->Draw();
	TLatex *labelSpectraPi0LabelPCMpPb2 = new TLatex(0.65,0.85,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
	SetStyleTLatex( labelSpectraPi0LabelPCMpPb2, 0.03,4);
	labelSpectraPi0LabelPCMpPb2->Draw();

	TLatex *labelSpectraPi0LabelpPbPCM = new TLatex(0.65,0.8,"|#it{y}_{#pi^{0},lab}| < 0.8");
	SetStyleTLatex( labelSpectraPi0LabelpPbPCM, 0.03,4);
	labelSpectraPi0LabelpPbPCM->Draw();  

	graphInvYieldPi07TeVBinPCMStatSystErrScaled->SetFillColor(0);
	graphInvYieldPi02760GeVBinPCMScaled->SetFillColor(0);
	graphErrosInterPolation5023GeVBinPCMScaled->SetFillColor(0);
	TLegend* legendInvYieldsPCM = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYieldsPCM->SetFillColor(0);
	legendInvYieldsPCM->SetLineColor(0);
	legendInvYieldsPCM->SetNColumns(1);
	legendInvYieldsPCM->SetTextSize(0.03);
	legendInvYieldsPCM->AddEntry(graphInvYieldPi07TeVBinPCMStatSystErrScaled,"pp@7 TeV (x4)","pef");
	legendInvYieldsPCM->AddEntry(graphErrosInterPolation5023GeVBinPCMScaled,"pp@5.02 TeV - Interpolation (x2)","pef");
	legendInvYieldsPCM->AddEntry(graphInvYieldPi02760GeVBinPCMScaled,"pp@2.76 TeV (x1)","pef");
	
	
	legendInvYieldsPCM->Draw();

		
	
	cInvYieldComparisonPCM->Update();
	cInvYieldComparisonPCM->SaveAs(Form("%s/InvYield_Comparison_PP_PCM_2760GeV_5023GeV_7TeV.%s",outputDir.Data(),suffix.Data()));
	
	//cout<<"Check where the something strange is"<<endl; 
	
	
	TGraphErrors* graphRebinnedAStat = NULL; 
	TGraphErrors* graphRebinnedASyst = NULL; 
	TGraphErrors* graphRebinnedBStat = NULL; 
	TGraphErrors* graphRebinnedBSyst = NULL;
	
	TGraphErrors* graphRatioPPReferencePCM_Pythia    = CalculateRatioBetweenSpectraWithDifferentBinning(graphppPythia,graphppPythia, graphErrosInterPolation5023GeVBinPCMStatErr, graphErrosInterPolation5023GeVBinPCMSystErr,       kTRUE,  kTRUE,&graphRebinnedAStat,&graphRebinnedASyst,&graphRebinnedBStat,&graphRebinnedBSyst);
	TGraphErrors* graphRatioPPReferenceDalitz_Pythia = CalculateRatioBetweenSpectraWithDifferentBinning(graphppPythia,graphppPythia, graphErrosInterPolation5023GeVBinDalitzStatErr, graphErrosInterPolation5023GeVBinDalitzSystErr, kTRUE,  kTRUE,&graphRebinnedAStat,&graphRebinnedASyst,&graphRebinnedBStat,&graphRebinnedBSyst);
	
	//cout<<"Check where the something strange is"<<endl; 
	
  
	
	TCanvas* cInvYieldComparisonPCM_Pyhtia = new TCanvas("cInvYieldComparisonPCM_Pyhtia","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldComparisonPCM_Pyhtia,  0.15, 0.02, 0.03, 0.1);
	cInvYieldComparisonPCM_Pyhtia->SetLogx();
	cInvYieldComparisonPCM_Pyhtia->SetLogy();
	
	TPad* padInvYieldComparisonPCM_Pyhtia = new TPad("padInvYieldComparisonPCM_Pyhtia", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldComparisonPCM_Pyhtia, 0.15, 0.02, 0.03, 0.1);
	padInvYieldComparisonPCM_Pyhtia->Draw();
	
	
	
	TH2F * histoInvYieldComparisonPCM_Pyhtia;
	histoInvYieldComparisonPCM_Pyhtia = new TH2F("histoInvYieldComparisonPCM_Pyhtia","histoInvYieldComparisonPCM_Pyhtia",1000,0.23,30.,1000,2e0,10e11);
	SetStyleHistoTH2ForGraphs(histoInvYieldComparisonPCM_Pyhtia, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldComparisonPCM_Pyhtia->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphppPythia,20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeVBinPCM,  20, 1, kBlue+1,  kBlue+1);

	graphErrosInterPolation5023GeVBinPCM->Draw("p,same,e1");
	graphppPythia->Draw("p,same,e1");
	
	
	TLatex *LabelpPbPCM_Pyhtia = new TLatex(0.65,0.9,"ALICE work in progress");
	SetStyleTLatex( LabelpPbPCM_Pyhtia, 0.03,4);
	LabelpPbPCM_Pyhtia->Draw();
	TLatex *labelSpectraPi0LabelPCMpPb2_Pythia = new TLatex(0.65,0.85,"#pi^{0} #rightarrow #gamma #gamma #rightarrow e^{+}e^{-} e^{+}e^{-} (PCM)");
	SetStyleTLatex( labelSpectraPi0LabelPCMpPb2_Pythia, 0.03,4);
	labelSpectraPi0LabelPCMpPb2_Pythia->Draw();

	TLatex *labelSpectraPi0LabelpPbPCM_Pythia = new TLatex(0.65,0.8,"|#it{y}_{#pi^{0},lab}| < 0.8");
	SetStyleTLatex( labelSpectraPi0LabelpPbPCM_Pythia, 0.03,4);
	labelSpectraPi0LabelpPbPCM_Pythia->Draw();  

	graphErrosInterPolation5023GeVBinPCM->SetFillColor(0);
	TLegend* legendInvYieldsPCM_Pythia = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYieldsPCM_Pythia->SetFillColor(0);
	legendInvYieldsPCM_Pythia->SetLineColor(0);
	legendInvYieldsPCM_Pythia->SetNColumns(1);
	legendInvYieldsPCM_Pythia->SetTextSize(0.03);
	legendInvYieldsPCM_Pythia->AddEntry(graphErrosInterPolation5023GeVBinPCM,"pp@5.02 TeV - Interpolation","pef");
	legendInvYieldsPCM_Pythia->AddEntry(graphppPythia,"pp@5.02 TeV  Pythia","pef");	
	legendInvYieldsPCM_Pythia->Draw();
	cInvYieldComparisonPCM_Pyhtia->Update();
	cInvYieldComparisonPCM_Pyhtia->SaveAs(Form("%s/InvYield_Comparison_PP_PCM_Pythia.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	TCanvas* canvasRatioPCM_Pythia = new TCanvas("canvasRatioPCM_Pythia","PP reference PCM compared to Pythia",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioPCM_Pythia,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioPCM_Pythia = new TH2F("histo2DRatioPCM_Pythia","histo2DRatioPCM_Pythia",1000,0.,15.,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioPCM_Pythia, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioPCM_Pythia->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphRatioPPReferencePCM_Pythia,21,1.5, kBlack , kBlack);
	graphRatioPPReferencePCM_Pythia->Draw("E1psame");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioPCM_Pythia = new TLatex(0.16,0.9,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRatioPCM_Pythia,0.06,4); 
	textRatioPCM_Pythia->Draw();
	
	
	TLegend* legendRatioPCM_Pythia = new TLegend(0.18,0.15,0.9,0.21);
	legendRatioPCM_Pythia->SetFillColor(0);
	legendRatioPCM_Pythia->SetLineColor(0);
	legendRatioPCM_Pythia->SetNColumns(2);
	legendRatioPCM_Pythia->SetTextSize(0.045);
	legendRatioPCM_Pythia->AddEntry(graphRatioPPReferencePCM_Pythia,"pp reference Pythia / pp reference data (PCM)","p");
	legendRatioPCM_Pythia->Draw();
	
	
	canvasRatioPCM_Pythia->SaveAs(Form("%s/Ratio_PCM_Pythia.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	TCanvas* cInvYieldComparisonDalitz_Pyhtia = new TCanvas("cInvYieldComparisonDalitz_Pyhtia","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldComparisonDalitz_Pyhtia,  0.15, 0.02, 0.03, 0.1);
	cInvYieldComparisonDalitz_Pyhtia->SetLogx();
	cInvYieldComparisonDalitz_Pyhtia->SetLogy();
	
	TPad* padInvYieldComparisonDalitz_Pyhtia = new TPad("padInvYieldComparisonDalitz_Pyhtia", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldComparisonDalitz_Pyhtia, 0.15, 0.02, 0.03, 0.1);
	padInvYieldComparisonDalitz_Pyhtia->Draw();
	
	
	TH2F * histoInvYieldComparisonDalitz_Pyhtia;
	histoInvYieldComparisonDalitz_Pyhtia = new TH2F("histoInvYieldComparisonDalitz_Pyhtia","histoInvYieldComparisonDalitz_Pyhtia",1000,0.23,30.,1000,2e0,10e11);
	SetStyleHistoTH2ForGraphs(histoInvYieldComparisonDalitz_Pyhtia, "#it{p}_{T} (GeV/#it{c})", "#it{E}#frac{d^{2}#sigma}{d#it{p}^{3}}(pb GeV^{-2} #it{c}^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldComparisonDalitz_Pyhtia->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphppPythia,20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphErr(graphErrosInterPolation5023GeVBinDalitz,  20, 1, kBlue+1,  kBlue+1);
	
	
	
	
	
	graphErrosInterPolation5023GeVBinDalitz->Draw("p,same,e1");
	graphppPythia->Draw("p,same,e1");
	
	
	TLatex *LabelpPbDalitz_Pyhtia = new TLatex(0.65,0.9,"ALICE work in progress");
	SetStyleTLatex( LabelpPbDalitz_Pyhtia, 0.03,4);
	LabelpPbDalitz_Pyhtia->Draw();
	TLatex *labelSpectraPi0LabelDalitzpPb2_Pythia = new TLatex(0.65,0.85,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz)");
	SetStyleTLatex( labelSpectraPi0LabelDalitzpPb2_Pythia, 0.03,4);
	labelSpectraPi0LabelDalitzpPb2_Pythia->Draw();

	TLatex *labelSpectraPi0LabelpPbDalitz_Pythia = new TLatex(0.65,0.8,"|#it{y}_{#pi^{0},lab}| < 0.8");
	SetStyleTLatex( labelSpectraPi0LabelpPbDalitz_Pythia, 0.03,4);
	labelSpectraPi0LabelpPbDalitz_Pythia->Draw();  

	graphErrosInterPolation5023GeVBinDalitz->SetFillColor(0);
	TLegend* legendInvYieldsDalitz_Pythia = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYieldsDalitz_Pythia->SetFillColor(0);
	legendInvYieldsDalitz_Pythia->SetLineColor(0);
	legendInvYieldsDalitz_Pythia->SetNColumns(1);
	legendInvYieldsDalitz_Pythia->SetTextSize(0.03);
	legendInvYieldsDalitz_Pythia->AddEntry(graphErrosInterPolation5023GeVBinDalitz,"pp@5.02 TeV - Interpolation","pef");
	legendInvYieldsDalitz_Pythia->AddEntry(graphppPythia,"pp@5.02 TeV  Pythia","pef");
	
	
	legendInvYieldsDalitz_Pythia->Draw();

		
	
	cInvYieldComparisonDalitz_Pyhtia->Update();
	cInvYieldComparisonDalitz_Pyhtia->SaveAs(Form("%s/InvYield_Comparison_PP_Dalitz_Pythia.%s",outputDir.Data(),suffix.Data()));
	
	
	TCanvas* canvasRatioDalitz_Pythia = new TCanvas("canvasRatioDalitz_Pythia","PP reference Dalitz compared to Pythia",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioDalitz_Pythia,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioDalitz_Pythia = new TH2F("histo2DRatioDalitz_Pythia","histo2DRatioDalitz_Pythia",1000,0.,15.,1000,0.3,2.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioDalitz_Pythia, "#it{p}_{T} (GeV/#it{c})","#it{R}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioDalitz_Pythia->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphRatioPPReferenceDalitz_Pythia,21,1.5, kBlack , kBlack);
	graphRatioPPReferenceDalitz_Pythia->Draw("E1psame");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioDalitz_Pythia = new TLatex(0.16,0.9,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRatioDalitz_Pythia,0.06,4); 
	textRatioDalitz_Pythia->Draw();
	
	
	TLegend* legendRatioDalitz_Pythia = new TLegend(0.18,0.15,0.9,0.21);
	legendRatioDalitz_Pythia->SetFillColor(0);
	legendRatioDalitz_Pythia->SetLineColor(0);
	legendRatioDalitz_Pythia->SetNColumns(2);
	legendRatioDalitz_Pythia->SetTextSize(0.045);
	legendRatioDalitz_Pythia->AddEntry(graphRatioPPReferenceDalitz_Pythia,"pp reference Pythia / pp reference data(Dalitz)","p");
	legendRatioDalitz_Pythia->Draw();
	
	
	canvasRatioDalitz_Pythia->SaveAs(Form("%s/Ratio_Dalitz_Pythia.%s",outputDir.Data(),suffix.Data()));
	
		
	
	
	
	//////////////////////////////////////////
	
	
	
	
	
	TCanvas* cInvYieldPi0DalitzpPb5023GeVXShifted = new TCanvas("cInvYieldPi0DalitzpPb5023GeVXShifted","Pi0 Dalitz pPb at 5.23 TeV X shifted",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldPi0DalitzpPb5023GeVXShifted,  0.15, 0.02, 0.03, 0.1);
	cInvYieldPi0DalitzpPb5023GeVXShifted->SetLogx();
	cInvYieldPi0DalitzpPb5023GeVXShifted->SetLogy();
	
	TPad* padInvYieldPi0DalitzpPb5023GeVXShifted = new TPad("padInvYieldPi0DalitzpPb5023GeVXShifted", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPi0DalitzpPb5023GeVXShifted, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPi0DalitzpPb5023GeVXShifted->Draw();
	
	
	TH2F * histoInvYieldPi0DalitzpPb5023GeVXShifted;
	histoInvYieldPi0DalitzpPb5023GeVXShifted = new TH2F("histoInvYieldPi0DalitzpPb5023GeVXShifted","histoInvYieldPi0DalitzpPb5023GeVXShifted",1000,0.23,30.,1000,2e-8,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPi0DalitzpPb5023GeVXShifted, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPi0DalitzpPb5023GeVXShifted->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr, 	20, 1, kRed,   kRed);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0DalitzpPb5023GeVComplErr,  	        20, 1, kBlue,  kBlue);
	
	graphInvYieldPi0DalitzpPb5023GeVComplErr->SetFillColor(0);
	graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->SetFillColor(0);
	
	graphInvYieldPi0DalitzpPb5023GeVComplErr->Draw("p,e1");
	graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr->Draw("p,same,e1");
	
	TLegend* legendInvYieldPi0DalitzpPb5023GeVXShifted = new TLegend(0.20,0.25,0.70,0.40);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->SetFillColor(0);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->SetLineColor(0);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->SetTextSize(0.03);
	legendInvYieldPi0DalitzpPb5023GeVXShifted->AddEntry(graphInvYieldPi0DalitzpPb5023GeVComplErr,"#pi^{0} p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV (Dalitz)");
	legendInvYieldPi0DalitzpPb5023GeVXShifted->AddEntry(graphInvYieldPi0DalitzpPb5023GeVXShiftedComplErr,"#pi^{0} p-Pb at #sqrt{#it{s}_{NN}} = 5.02  XBin shifted (Dalitz)");
	legendInvYieldPi0DalitzpPb5023GeVXShifted->Draw("same");
	
		
	
	cInvYieldPi0DalitzpPb5023GeVXShifted->Update();
	cInvYieldPi0DalitzpPb5023GeVXShifted->SaveAs(Form("%s/InvYield_XShifted_Dalitz_pPb_5023GeV.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	TCanvas* cInvYieldPi0PCMpPb5023GeVXShifted = new TCanvas("cInvYieldPi0PCMpPb5023GeVXShifted","Pi0 PCM pPb at 5.23 TeV X shifted",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldPi0PCMpPb5023GeVXShifted,  0.15, 0.02, 0.03, 0.1);
	cInvYieldPi0PCMpPb5023GeVXShifted->SetLogx();
	cInvYieldPi0PCMpPb5023GeVXShifted->SetLogy();
	
	TPad* padInvYieldPi0PCMpPb5023GeVXShifted = new TPad("padInvYieldPi0PCMpPb5023GeVXShifted", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPi0PCMpPb5023GeVXShifted, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPi0PCMpPb5023GeVXShifted->Draw();
	
	
	TH2F * histoInvYieldPi0PCMpPb5023GeVXShifted;
	histoInvYieldPi0PCMpPb5023GeVXShifted = new TH2F("histoInvYieldPi0PCMpPb5023GeVXShifted","histoInvYieldPi0PCMpPb5023GeVXShifted",1000,0.23,30.,1000,2e-8,5e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPi0PCMpPb5023GeVXShifted, "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev.}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (c/GeV)^{2}", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPi0PCMpPb5023GeVXShifted->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr, 	20, 1, kRed,   kRed);
	DrawGammaSetMarkerTGraphAsym(graphInvYieldPi0PCMpPb5023GeVComplErr,  	        20, 1, kBlue,  kBlue);
	
	
	graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->SetFillColor(0);
	graphInvYieldPi0PCMpPb5023GeVComplErr->SetFillColor(0);
	
	
	graphInvYieldPi0PCMpPb5023GeVComplErr->Draw("p,e1");
	graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr->Draw("p,same,e1");
	
	TLegend* legendInvYieldPi0PCMpPb5023GeVXShifted = new TLegend(0.20,0.25,0.70,0.40);
	legendInvYieldPi0PCMpPb5023GeVXShifted->SetFillColor(0);
	legendInvYieldPi0PCMpPb5023GeVXShifted->SetLineColor(0);
	legendInvYieldPi0PCMpPb5023GeVXShifted->SetTextSize(0.03);
	legendInvYieldPi0PCMpPb5023GeVXShifted->AddEntry(graphInvYieldPi0PCMpPb5023GeVComplErr,"#pi^{0} p-Pb at #sqrt{#it{s}_{NN}} = 5.02 TeV (PCM)");
	legendInvYieldPi0PCMpPb5023GeVXShifted->AddEntry(graphInvYieldPi0PCMpPb5023GeVXShiftedComplErr,"#pi^{0} p-Pb at #sqrt{#it{s}_{NN}} = 5.02  XBin shifted (PCM)");
	legendInvYieldPi0PCMpPb5023GeVXShifted->Draw("same");
	
		
	
	cInvYieldPi0PCMpPb5023GeVXShifted->Update();
	cInvYieldPi0PCMpPb5023GeVXShifted->SaveAs(Form("%s/InvYield_XShifted_PCM_pPb_5023GeV.%s",outputDir.Data(),suffix.Data()));
	
	
		
	
	
	
	TCanvas* canvasRpPbPCM = new TCanvas("canvasRpPbPCM","R_{pPb} PCM",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbPCM,  0.1, 0.01, 0.015, 0.11);
	TH2F * histo2DRpPbPCM = new TH2F("histo2DRpPbPCM","histo2DRpPbPCM",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbPCM,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbPCM->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbPCM->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbPCM->GetYaxis()->CenterTitle(kTRUE);
	
	histo2DRpPbPCM->Draw(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRpPbPCMSystErr,20,1.5, colorsArray[kPCM], colorsArray[kPCM], 1, kTRUE);
	graphRpPbPCMSystErr->Draw("2same");
	
	
	DrawGammaSetMarkerTGraphAsym(graphRpPbPCMStatErr,20,1.5, colorsArray[kPCM],   colorsArray[kPCM]);
 	graphRpPbPCMStatErr->Draw("p,same,e");
	graphRpPbPCMStatErr->SetFillColor(0);
	
	
	
	
	DrawGammaLines(0., 14.5,1., 1.,1.,kBlack,2);
	
	TLatex *labelpPbPCM = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( labelpPbPCM, 0.04,4);
	labelpPbPCM->Draw();
	
	TLatex *labelpPbPCMEnergy =  new TLatex(0.15,0.87, energyLabel.Data());
	SetStyleTLatex( labelpPbPCMEnergy, 0.04,4);
	labelpPbPCMEnergy->Draw();
	
	
	TLegend* legendRpPbPCM = new TLegend(0.15,0.79,0.70,0.84);
	legendRpPbPCM->SetFillColor(0);
	legendRpPbPCM->SetLineColor(0);
	legendRpPbPCM->SetTextFont(42);
	legendRpPbPCM->SetNColumns(3);
	legendRpPbPCM->SetTextSize(0.04);
	legendRpPbPCM->AddEntry(graphRpPbPCMStatErr,Form("%s, |#it{y}| < 0.8",PCMLabel.Data()),"p");
	
	legendRpPbPCM->Draw();
	canvasRpPbPCM->Update();
	
	

	canvasRpPbPCM->Update(); 
	
	canvasRpPbPCM->SaveAs(Form("%s/RpPb_PCM.%s",outputDir.Data(),suffix.Data()));
	

	
	TCanvas* canvasRpPbDalitz = new TCanvas("canvasRpPbDalitz","R_{pPb} Dalitz",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbDalitz,  0.1, 0.01, 0.015, 0.11);
	TH2F * histo2DRpPbDalitz = new TH2F("histo2DRpPbDalitz","histo2DRpPbDalitz",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbDalitz,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbDalitz->GetXaxis()->SetRangeUser(0.0,10.5);
	histo2DRpPbDalitz->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbDalitz->GetYaxis()->CenterTitle(kTRUE);
	
	histo2DRpPbDalitz->Draw(); 
	
	DrawGammaSetMarkerTGraphAsym(graphRpPbDalitzSystErr,20,1.5, colorsArray[kDalitz], colorsArray[kDalitz], 1, kTRUE);
	graphRpPbDalitzSystErr->Draw("2same");
	
	DrawGammaSetMarkerTGraphAsym(graphRpPbDalitzStatErr,20,1.5, colorsArray[kDalitz], colorsArray[kDalitz]);
 	graphRpPbDalitzStatErr->Draw("p,same,e");
	graphRpPbDalitzStatErr->SetFillColor(0);
	
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	TLatex *labelpPb3 = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( labelpPb3, 0.04,4);
	labelpPb3->Draw();
	
	TLatex *labelpPb4 =  new TLatex(0.15,0.87, energyLabel.Data());
	SetStyleTLatex( labelpPb4, 0.04,4);
	labelpPb4->Draw();
	
	
	TLegend* legendRpPbDalitz = new TLegend(0.15,0.79,0.75,0.84);
	legendRpPbDalitz->SetFillColor(0);
	legendRpPbDalitz->SetLineColor(0);
	legendRpPbDalitz->SetTextFont(42);
	legendRpPbDalitz->SetNColumns(3);
	legendRpPbDalitz->SetTextSize(0.04);
	legendRpPbDalitz->AddEntry(graphRpPbDalitzStatErr,Form("%s, |#it{y}| < 0.8",DalitzLabel.Data()),"p");
	
	legendRpPbDalitz->Draw();
	canvasRpPbDalitz->Update();
	
	canvasRpPbDalitz->SaveAs(Form("%s/RpPb_Dalitz.%s",outputDir.Data(),suffix.Data()));
	
        
	
	
	TCanvas* canvasRpPbPCMDalitz = new TCanvas("canvasRpPbPCMDalitz","Combined R_{pPb}",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbPCMDalitz,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * PCMhisto2DRpPbPCMDalitz = new TH2F("PCMhisto2DRpPbPCMDalitz","PCMhisto2DRpPbPCMDalitz",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(PCMhisto2DRpPbPCMDalitz,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	PCMhisto2DRpPbPCMDalitz->GetXaxis()->SetRangeUser(0.0,10.5);
	PCMhisto2DRpPbPCMDalitz->GetYaxis()->SetRangeUser(0.19,1.81);
	PCMhisto2DRpPbPCMDalitz->GetYaxis()->CenterTitle();
	PCMhisto2DRpPbPCMDalitz->Draw();
	
	
	graphRpPbDalitzSystErr->Draw("2same");
	graphRpPbDalitzStatErr->Draw("p,same,e");
	graphRpPbPCMSystErr->Draw("2same");
	graphRpPbPCMStatErr->Draw("p,same,e");
	
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	TLatex *textRpPbPCMDalitz = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbPCMDalitz,0.04,4); 
	textRpPbPCMDalitz->Draw();
	
	TLatex *LabelpPbPCMDalitz = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbPCMDalitz, 0.04,4);
	LabelpPbPCMDalitz->Draw();
	
	
	
	
	TLegend* legendRpPbPCMDalitz = new TLegend(0.3,0.20,0.75,0.35);
	legendRpPbPCMDalitz->SetFillColor(0);
	legendRpPbPCMDalitz->SetLineColor(0);
	legendRpPbPCMDalitz->SetTextSize(0.03);
	legendRpPbPCMDalitz->SetTextFont(42);
	legendRpPbPCMDalitz->SetTextSize(0.04);
	
	legendRpPbPCMDalitz->AddEntry(graphRpPbPCMStatErr,   Form("%s, |#it{y}| < 0.8",PCMLabel.Data()),"p");
	legendRpPbPCMDalitz->AddEntry(graphRpPbDalitzStatErr,Form("%s, |#it{y}| < 0.8",DalitzLabel.Data()),"p");
	
	
	legendRpPbPCMDalitz->Draw();
	
	
	canvasRpPbPCMDalitz->SaveAs(Form("%s/RpPb_Comparison_PCM_Dalitz.%s",outputDir.Data(),suffix.Data()));
	
	
	//////////////////////////////Comparing all GA spectra ///////////////////////////////////////////////	


	TCanvas* canvasRpPbCombineAllGA = new TCanvas("canvasRpPbCombineAllGA","Combined #it{R}_{pPb}",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRpPbCombineAllGA,  0.1, 0.02, 0.02, 0.13);
	
	TH2F * histo2DRpPbCombineAllGA = new TH2F("histo2DRpPbCombineAllGA","histo2DRpPbCombineAllGA",1000,0.,20.,1000,0.3,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRpPbCombineAllGA, "#it{p}_{T} (GeV/#it{c})","#it{R}^{#pi^{0}}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRpPbCombineAllGA->DrawCopy(); 
	
	 DrawGammaSetMarkerTGraphAsym(graphRpPbPHOSSystErr,20,1.5, kRed+1,kRed+1,1.0,kTRUE);
	 DrawGammaSetMarkerTGraphAsym(graphRpPbPHOSStatErr,20,1.5, kRed+1 ,kRed+1);
	
	 DrawGammaSetMarkerTGraphAsym(graphRpPbPCMSystErr,20,1.5, kBlack,kBlack,1.0,kTRUE);
	 DrawGammaSetMarkerTGraphAsym(graphRpPbPCMStatErr,20,1.5, kBlack ,kBlack);
	
	 DrawGammaSetMarkerTGraphAsym(graphRpPbDalitzSystErr,20,1.5, kCyan+2,kCyan+2,1.0,kTRUE);
	 DrawGammaSetMarkerTGraphAsym(graphRpPbDalitzStatErr,20,1.5, kCyan+2 ,kCyan+2);

	graphRpPbDalitzSystErr->Draw("2same");
	graphRpPbDalitzStatErr->Draw("p,same,e");
	graphRpPbPHOSSystErr->Draw("2same");
	graphRpPbPHOSStatErr->Draw("p,same,e");
	graphRpPbPCMSystErr->Draw("2same");
	graphRpPbPCMStatErr->Draw("p,same,e");
	
	DrawGammaLines(0., 20.,1., 1.,2.,kBlack);
	
	TLatex *textRpPbCombineAllGA = new TLatex(0.2,0.35,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRpPbCombineAllGA,0.035,4); 
	textRpPbCombineAllGA->Draw();
	
	
	
	TLegend* legendRpPbCombineAllGA = new TLegend(0.18,0.15,0.5,0.3);
	legendRpPbCombineAllGA->SetFillColor(0);
	legendRpPbCombineAllGA->SetLineColor(0);
	legendRpPbCombineAllGA->SetTextSize(0.03);
	legendRpPbCombineAllGA->AddEntry(graphRpPbPCMStatErr,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PCM)","pef");
	legendRpPbCombineAllGA->AddEntry(graphRpPbPHOSStatErr,"#pi^{0} #rightarrow #gamma#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (PHOS)","pef");
	legendRpPbCombineAllGA->AddEntry(graphRpPbDalitzStatErr,"#pi^{0} #rightarrow e^{+}e^{-}#gamma #rightarrow e^{+}e^{-}e^{+}e^{-} (Dalitz)","pef");	
	legendRpPbCombineAllGA->Draw();
	
	
	canvasRpPbCombineAllGA->SaveAs(Form("%s/RpPb_Comparison_AllGA.%s",outputDir.Data(),suffix.Data()));	
		
	//////////////////////////////Comparing to Charged Particles ///////////////////////////////////////////////
	
	
	
		
	
	DrawGammaSetMarkerTGraphAsym(graphAsymmChargedParticlesRpPbSystErr,20,1.5, kOrange+2, kOrange+2, 1, kTRUE,kOrange-9);

	
	
	DrawGammaSetMarkerTGraphAsym(graphAsymmChargedParticlesRpPbStatErr,20,1.5, kOrange+2, kOrange+2);
	graphAsymmChargedParticlesRpPbStatErr->SetFillColor(0);
	
	
	
	
	
	TCanvas* canvasRpPbChargedParPCM = new TCanvas("canvasRpPbChargedParCombined","Combined R_{pPb} with charged particles",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbChargedParPCM,  0.1, 0.01, 0.015, 0.11);
	
	
	
	
	TH2F * histo2DRpPbChargedParPCM = new TH2F("histo2DRpPbChargedParPCM","histo2DRpPbChargedParPCM",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedParPCM,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbChargedParPCM->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbChargedParPCM->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbChargedParPCM->GetYaxis()->CenterTitle();
	
	histo2DRpPbChargedParPCM->Draw();
	
	
	
	
	graphAsymmChargedParticlesRpPbSystErr->Draw("2p,same,e");
	graphAsymmChargedParticlesRpPbStatErr->Draw("p,same,e");
	
	
	graphRpPbPCMSystErr->Draw("2p,same,e");
	graphRpPbPCMStatErr->Draw("p,same,e");
	
	
	DrawGammaLines(0., 14.5,1., 1.,1.,kBlack,2);
	
	TLatex *textRpPbPCM2 = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbPCM2,0.04,4); 
	textRpPbPCM2->Draw();
	
	TLatex *LabelpPbPCMChargedPar = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbPCMChargedPar, 0.04,4);
	LabelpPbPCMChargedPar->Draw();
	
	
	 TLegend* legendRpPbPCM2 = new TLegend(0.3,0.20,0.75,0.35);
	legendRpPbPCM2->SetFillColor(0);
	legendRpPbPCM2->SetLineColor(0);
	legendRpPbPCM2->SetTextSize(0.03);
	legendRpPbPCM2->SetTextFont(42);
	legendRpPbPCM2->SetTextSize(0.04);
	
	legendRpPbPCM2->AddEntry(graphRpPbPCMStatErr,Form("%s, |#it{y}| < 0.8",PCMLabel.Data()),"p");
	legendRpPbPCM2->AddEntry(graphAsymmChargedParticlesRpPbStatErr,"ALICE, NSD, charged particles, |#eta_{cms}| < 0.3","p");
	
	legendRpPbPCM2->Draw();
	
	
	
	
	
	
	
	
	canvasRpPbChargedParPCM->SaveAs(Form("%s/RpPb_Comparison_PCM_ChargedParticles.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	TCanvas* canvasRpPbChargedParDalitz = new TCanvas("canvasRpPbChargedParDalitz","Dalitz R_{pPb} with charged particles",200,10,1000,700);
	
	
	DrawGammaCanvasSettings( canvasRpPbChargedParDalitz,  0.1, 0.01, 0.015, 0.11);
	
	
	
	
	TH2F * histo2DRpPbChargedParDalitz = new TH2F("histo2DRpPbChargedParDalitz","histo2DRpPbChargedParDalitz",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedParDalitz,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbChargedParDalitz->GetXaxis()->SetRangeUser(0.0,10.5);
	histo2DRpPbChargedParDalitz->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbChargedParDalitz->GetYaxis()->CenterTitle();
	
	histo2DRpPbChargedParDalitz->Draw();
	
	
	
	graphAsymmChargedParticlesRpPbSystErr->Draw("2p,same,e");
	graphAsymmChargedParticlesRpPbStatErr->Draw("p,same,e");
	
	graphRpPbDalitzStatErr->Draw("p,same,e");
	graphRpPbDalitzSystErr->Draw("2p,same,e");
	
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	TLatex *textRpPbDalitz2 = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbDalitz2,0.04,4); 
	textRpPbDalitz2->Draw();
	
	
	

	
	
	TLatex *LabelpPb4 = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPb4, 0.04,4);
	LabelpPb4->Draw();

	
	TLegend* legendRpPbDalitz2 = new TLegend(0.3,0.20,0.75,0.35);
	legendRpPbDalitz2->SetFillColor(0);
	legendRpPbDalitz2->SetLineColor(0);
	legendRpPbDalitz2->SetTextSize(0.03);
	legendRpPbDalitz2->SetTextFont(42);
	legendRpPbDalitz2->SetTextSize(0.04);
	
	legendRpPbDalitz2->AddEntry(graphRpPbDalitzStatErr,Form("%s, |#it{y}| < 0.8",DalitzLabel.Data()),"p");
	legendRpPbDalitz2->AddEntry(graphAsymmChargedParticlesRpPbStatErr,"ALICE, NSD, charged particles, |#eta_{cms}| < 0.3","p");
	
	legendRpPbDalitz2->Draw();
	
	
	canvasRpPbChargedParDalitz->SaveAs(Form("%s/RpPb_Comparison_Dalitz_ChargedParticles.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	cout<<"Llego aqui antes theory"<<endl;
	
	//return;
	
	//****************************** extracting EPS09s predictions**************************************
	
	ifstream 		inDSS;
	ifstream                inAKK;
	ifstream                inKKP;
	ifstream 		inCGC;	
	
	
	Int_t nlinesEPSsPi0fDSS = 0;
	
	inDSS.open(fileNameEPS09sPi0DSS,ios_base::in);
	
	Double_t xEPSsPi0fDSS[100],yEPSsPi0fDSS[100];
	
	Double_t xUpErrorEPSsPi0DSS[100],xDownErrorEPSsPi0DSS[100];
	Double_t yUpErrorEPSsPi0DSS[100],yDownErrorEPSsPi0DSS[100];
	
	
	while(!inDSS.eof()){
			nlinesEPSsPi0fDSS++;
			inDSS >> xEPSsPi0fDSS[nlinesEPSsPi0fDSS]  >> yEPSsPi0fDSS[nlinesEPSsPi0fDSS] >> yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]>>yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS];
			yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = ( yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
			yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = -1*( yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] - yEPSsPi0fDSS[nlinesEPSsPi0fDSS] ) /yEPSsPi0fDSS[nlinesEPSsPi0fDSS];
			xUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]   = 0;
			xDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS] = 0;
			
			    
			cout << nlinesEPSsPi0fDSS << "         "  << xEPSsPi0fDSS[nlinesEPSsPi0fDSS] << "         "  <<yEPSsPi0fDSS[nlinesEPSsPi0fDSS]<<"         "<<yUpErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<<"          "<<yDownErrorEPSsPi0DSS[nlinesEPSsPi0fDSS]<< endl;
	
	}
	inDSS.close();
	
	cout<<"Llego aqui primer theory file"<<endl;
	
	TGraph* graphPi0DSS5000 = new TGraph(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS);
	TGraphAsymmErrors* graphAsymmErrorsPi0DSS5000 = new TGraphAsymmErrors(nlinesEPSsPi0fDSS,xEPSsPi0fDSS,yEPSsPi0fDSS,xDownErrorEPSsPi0DSS,xUpErrorEPSsPi0DSS,yDownErrorEPSsPi0DSS,yUpErrorEPSsPi0DSS);
	
	
	
	
	Int_t nlinesPi0AKK = 0;
		
	inAKK.open(fileNameEPS09sPi0AKK,ios_base::in);
	
	Double_t xESPsPi0AKK[100], yESPsPi0AKK[100];
	
		
	while(!inAKK.eof()){
			nlinesPi0AKK++;
			inAKK >> xESPsPi0AKK[nlinesPi0AKK]  >> yESPsPi0AKK[nlinesPi0AKK];
			cout << nlinesPi0AKK << "         "  << xESPsPi0AKK[nlinesPi0AKK] << "         "  <<xESPsPi0AKK[nlinesPi0AKK]<<endl;
	
	}
	inAKK.close();
	
	TGraph* graphPi0ESP09sPi0AKK = new TGraph(nlinesPi0AKK,xESPsPi0AKK,yESPsPi0AKK);	
	
	
	Int_t nlinesPi0KKP = 0;
		
	inKKP.open(fileNameEPS09sPi0KKP,ios_base::in);
	
	Double_t xESPsPi0KKP[100], yESPsPi0KKP[100];
	
		
	while(!inKKP.eof()){
			nlinesPi0KKP++;
			inKKP >> xESPsPi0KKP[nlinesPi0KKP]  >> yESPsPi0KKP[nlinesPi0KKP];
			cout << nlinesPi0KKP << "         "  << xESPsPi0KKP[nlinesPi0KKP] << "         "  <<xESPsPi0KKP[nlinesPi0KKP]<<endl;
	
	}
	inKKP.close();
	
	
	cout<<"Final 2 "<<endl;
	//return;
	
	TGraph* graphPi0ESP09sPi0KKP = new TGraph(nlinesPi0KKP,xESPsPi0KKP,yESPsPi0KKP);	
	
		
	
	graphPi0ESP09sPi0KKP->RemovePoint(0);
	graphPi0ESP09sPi0AKK->RemovePoint(0);
	graphPi0DSS5000->RemovePoint(0);
	graphAsymmErrorsPi0DSS5000->RemovePoint(0);
	
	graphPi0ESP09sPi0KKP->SetLineColor(kRed);
	graphPi0ESP09sPi0KKP->SetLineWidth(3);
	graphPi0ESP09sPi0KKP->SetLineStyle(10);
	
	
	//graphPi0ESP09sPi0KKP->RemovePoint(0);
	graphPi0ESP09sPi0AKK->SetLineColor(kBlue);
	graphPi0ESP09sPi0AKK->SetLineWidth(3);
	graphPi0ESP09sPi0AKK->SetLineStyle(7);
	
	
	
	graphPi0DSS5000->SetLineColor(8);
	graphPi0DSS5000->SetLineWidth(3);
	
	

	graphAsymmErrorsPi0DSS5000->SetFillColor(kYellow-7);
	graphAsymmErrorsPi0DSS5000->SetFillStyle(1001);
		
	
        


	Int_t nlinesPi0CGC = 0;	
	inCGC.open(fileNameEPS09sPi0CGC,ios_base::in);
	
	Double_t xESPsPi0CGC[100], yESPsPi0CGC[100];
	
		
	while(!inCGC.eof()){
			nlinesPi0CGC++;
			inCGC >> xESPsPi0CGC[nlinesPi0CGC]  >> yESPsPi0CGC[nlinesPi0CGC];
			cout << nlinesPi0CGC << "         "  << xESPsPi0CGC[nlinesPi0CGC] << "         "  <<yESPsPi0CGC[nlinesPi0CGC]<<endl;
	
	}
	inCGC.close();
	

	
	
	
	
	cout<<"Final 3"<<endl;
	
	//return;
	
	
	
	
	
	
	
	TGraph* graphPi0CGC = new TGraph(nlinesPi0CGC,xESPsPi0CGC,yESPsPi0CGC);	
	
	TCanvas* canvasRpPbESP09sCGCPCM = new TCanvas("canvasRpPbESP09sCGCPCM","PCM Combined predictions for R^{#pi0}_{pPb}",200,10,1000,700);
	
	
	DrawGammaCanvasSettings( canvasRpPbESP09sCGCPCM,  0.1, 0.01, 0.015, 0.13);

	
	
	TH2F * histo2DRpPbESP09sPCM = new TH2F("histo2DRpPbESP09sPCM","histo2DRpPbESP09sPCM",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbESP09sPCM,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbESP09sPCM->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbESP09sPCM->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbESP09sPCM->GetYaxis()->CenterTitle();
	
	
	histo2DRpPbESP09sPCM->DrawCopy(); 
	
	
	//TH2F * histo2DRpPbESP09sCGCPCM = new TH2F("histo2DRpPbESP09sCGCPCM","histo2DRpPbESP09sCGCPCM",1500,0.5,15.,1000,0.0,1.5);
	//SetStyleHistoTH2ForGraphs(histo2DRpPbESP09sCGCPCM, "p_{T} (GeV/c)","R^{#pi^{0}}_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.7, 512, 505); 
	//histo2DRpPbESP09sCGCPCM->DrawCopy(); 
	
	
	DrawGammaLines(0., 14.5,1., 1.,2.,kBlack);
		
	
        Float_t offsetYMaxLegendPCM = -0.15;
	Float_t xMinTheoryPCM = 0.30;
	
	TLatex *LabelRpPbESP09sCGCPCM = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelRpPbESP09sCGCPCM, 0.04,4);
	
	
        TLatex *textRpPbComparedPCM = new TLatex(0.15,0.87,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRpPbComparedPCM,0.04,4); 
	
	TLatex *textRpPbEPS09sPCM = new TLatex(xMinTheoryPCM,0.55 + offsetYMaxLegendPCM,"p-Pb, #sqrt{#it{s}_{NN}} = 5.0 TeV");
	SetStyleTLatex( textRpPbEPS09sPCM,0.04,4); 
	
	TLatex *textRpPbEPS09sPCM2 = new TLatex(xMinTheoryPCM,0.5 + offsetYMaxLegendPCM,"pQCD calculations");
	SetStyleTLatex( textRpPbEPS09sPCM2,0.04,4); 
	
	
	TLatex *textRpPbESP09sCGCPCM2= new TLatex(xMinTheoryPCM,0.23,"CGC calculations");
	SetStyleTLatex( textRpPbESP09sCGCPCM2,0.04,4); 	
	
	TLegend* legendRpPbESP09sCGCPCM = new TLegend(0.15,0.79,0.45,0.84);
	legendRpPbESP09sCGCPCM->SetFillColor(0);
	legendRpPbESP09sCGCPCM->SetLineColor(0);
	legendRpPbESP09sCGCPCM->SetTextSize(0.04);
	legendRpPbESP09sCGCPCM->SetTextFont(42);
	legendRpPbESP09sCGCPCM->AddEntry(graphRpPbPCMStatErr,Form("%s, |#it{y}| < 0.8",PCMLabel.Data()),"p");

	
	TLegend* legendRpPbESP09sPCM = new TLegend(xMinTheoryPCM,0.3+offsetYMaxLegendPCM,0.55,0.48+offsetYMaxLegendPCM);
	legendRpPbESP09sPCM->SetFillColor(0);
	legendRpPbESP09sPCM->SetLineColor(0);
	legendRpPbESP09sPCM->SetTextSize(0.04);
	legendRpPbESP09sPCM->SetTextFont(42);
	legendRpPbESP09sPCM->AddEntry(graphPi0ESP09sPi0KKP,"EPS09s KKP NLO");
	legendRpPbESP09sPCM->AddEntry(graphPi0ESP09sPi0AKK,"EPS09s AKK NLO");
	legendRpPbESP09sPCM->AddEntry(graphPi0DSS5000,"EPS09s fDSS NLO");
	legendRpPbESP09sPCM->AddEntry(graphAsymmErrorsPi0DSS5000,"EPS09s fDSS errors","pef");
	
	
	TLegend* legendRpPbCGCPCM = new TLegend(xMinTheoryPCM+.3,0.33+(offsetYMaxLegendPCM/4),0.85,0.48+offsetYMaxLegendPCM);
	legendRpPbCGCPCM->SetFillColor(0);
	legendRpPbCGCPCM->SetLineColor(0);
	legendRpPbCGCPCM->SetTextSize(0.04);
	legendRpPbCGCPCM->SetTextFont(42);
	legendRpPbCGCPCM->AddEntry(graphPi0CGC,"Color Glas Condensate","p");

	graphPi0ESP09sPi0KKP->SetFillColor(0);
	graphPi0ESP09sPi0AKK->SetFillColor(0);
	graphPi0DSS5000->SetFillColor(0);

	graphPi0CGC->SetFillColor(0);
	graphPi0CGC->SetMarkerColor(kMagenta-2);
	graphPi0CGC->SetMarkerStyle(21);
	graphPi0CGC->SetMarkerSize(1.5);
	graphPi0CGC->SetLineColor(kGreen+3);
	graphPi0CGC->SetLineWidth(3);
	graphPi0CGC->SetLineStyle(1);

	graphPi0ESP09sPi0KKP->SetFillColor(0);
	graphPi0ESP09sPi0AKK->SetFillColor(0);
	graphPi0DSS5000->SetFillColor(0);

	
	
	histo2DRpPbESP09sPCM->DrawCopy();
	legendRpPbESP09sCGCPCM->Draw("same");		
	textRpPbComparedPCM->Draw("same");
	textRpPbEPS09sPCM->Draw("same");
 	textRpPbEPS09sPCM2->Draw("same");
	LabelRpPbESP09sCGCPCM->Draw("same");
	
	legendRpPbESP09sPCM->Draw("same");
	legendRpPbCGCPCM->Draw("same");
	

	graphAsymmErrorsPi0DSS5000->Draw("same,E3");
	graphPi0CGC->Draw("same,p");
	graphPi0DSS5000->Draw("same,p,l");
	graphPi0ESP09sPi0AKK->Draw("same,p,l");
	graphPi0ESP09sPi0KKP->Draw("same,p,l");
	graphRpPbPCMSystErr->Draw("2p,same,e");
	graphRpPbPCMStatErr->Draw("p,same,e");
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	

	canvasRpPbESP09sCGCPCM->Update(); 
	
	canvasRpPbESP09sCGCPCM->SaveAs(Form("%s/RpPb_Comparison_PCM_Models.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	
	
	
	TCanvas* canvasRpPbESP09sDalitz = new TCanvas("canvasRpPbESP09sDalitz","Dalitz Combined predictions for R^{#pi0}_{pPb}",200,10,1000,700);
	
	
	
	DrawGammaCanvasSettings( canvasRpPbESP09sDalitz,  0.1, 0.01, 0.015, 0.13);

	
	
	
	TH2F * histo2DRpPbESP09sDalitz = new TH2F("histo2DRpPbESP09sDalitz","histo2DRpPbESP09sDalitz",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbESP09sDalitz,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbESP09sDalitz->GetXaxis()->SetRangeUser(0.0,10.5);
	histo2DRpPbESP09sDalitz->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbESP09sDalitz->GetYaxis()->CenterTitle();
	
	
	histo2DRpPbESP09sDalitz->DrawCopy(); 
	
	
	Float_t offsetYMaxLegend = -0.15;
	Float_t xMinTheory = 0.30;
	
	TLatex *LabelRpPbESP09sCGCDalitz = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelRpPbESP09sCGCDalitz, 0.04,4);
	
	
        TLatex *textRpPbComparedDalitz = new TLatex(0.15,0.87,"p-Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRpPbComparedDalitz,0.04,4); 
	
	TLatex *textRpPbEPS09sDalitz = new TLatex(xMinTheory,0.55 + offsetYMaxLegend,"p-Pb, #sqrt{#it{s}_{NN}} = 5.0 TeV");
	SetStyleTLatex( textRpPbEPS09sDalitz,0.04,4); 
	
	TLatex *textRpPbEPS09sDalitz2 = new TLatex(xMinTheory,0.5 + offsetYMaxLegend,"pQCD calculations");
	SetStyleTLatex( textRpPbEPS09sDalitz2,0.04,4); 
	
	
	TLatex *textRpPbESP09sCGCDalitz2= new TLatex(xMinTheory,0.23,"CGC calculations");
	SetStyleTLatex( textRpPbESP09sCGCDalitz2,0.04,4); 	
	
	TLegend* legendRpPbESP09sCGCDalitz = new TLegend(0.15,0.79,0.45,0.84);
	legendRpPbESP09sCGCDalitz->SetFillColor(0);
	legendRpPbESP09sCGCDalitz->SetLineColor(0);
	legendRpPbESP09sCGCDalitz->SetTextSize(0.04);
	legendRpPbESP09sCGCDalitz->SetTextFont(42);
	legendRpPbESP09sCGCDalitz->AddEntry(graphRpPbDalitzStatErr,Form("%s, |#it{y}| < 0.8",DalitzLabel.Data()),"p");

	
	TLegend* legendRpPbESP09sDalitz = new TLegend(xMinTheory,0.3+offsetYMaxLegend,0.55,0.48+offsetYMaxLegend);
	legendRpPbESP09sDalitz->SetFillColor(0);
	legendRpPbESP09sDalitz->SetLineColor(0);
	legendRpPbESP09sDalitz->SetTextSize(0.04);
	legendRpPbESP09sDalitz->SetTextFont(42);
	legendRpPbESP09sDalitz->AddEntry(graphPi0ESP09sPi0KKP,"EPS09s KKP NLO");
	legendRpPbESP09sDalitz->AddEntry(graphPi0ESP09sPi0AKK,"EPS09s AKK NLO");
	legendRpPbESP09sDalitz->AddEntry(graphPi0DSS5000,"EPS09s fDSS NLO");
	legendRpPbESP09sDalitz->AddEntry(graphAsymmErrorsPi0DSS5000,"EPS09s fDSS errors","pef");
	
	
	TLegend* legendRpPbCGCDalitz = new TLegend(xMinTheory+.3,0.33+(offsetYMaxLegend/4),0.85,0.48+offsetYMaxLegend);
	legendRpPbCGCDalitz->SetFillColor(0);
	legendRpPbCGCDalitz->SetLineColor(0);
	legendRpPbCGCDalitz->SetTextSize(0.04);
	legendRpPbCGCDalitz->SetTextFont(42);
	legendRpPbCGCDalitz->AddEntry(graphPi0CGC,"Color Glas Condensate","p");

	graphPi0ESP09sPi0KKP->SetFillColor(0);
	graphPi0ESP09sPi0AKK->SetFillColor(0);
	graphPi0DSS5000->SetFillColor(0);

	graphPi0CGC->SetFillColor(0);
	graphPi0CGC->SetMarkerColor(kMagenta-2);
	graphPi0CGC->SetMarkerStyle(21);
	graphPi0CGC->SetMarkerSize(1.5);
	graphPi0CGC->SetLineColor(kGreen+3);
	graphPi0CGC->SetLineWidth(3);
	graphPi0CGC->SetLineStyle(1);
	

	histo2DRpPbESP09sDalitz->DrawCopy();
	legendRpPbESP09sCGCDalitz->Draw("same");		
	textRpPbComparedDalitz->Draw("same");
	textRpPbEPS09sDalitz->Draw("same");
 	textRpPbEPS09sDalitz2->Draw("same");
	LabelRpPbESP09sCGCDalitz->Draw("same");
	
	legendRpPbESP09sDalitz->Draw("same");
	legendRpPbCGCDalitz->Draw("same");
	
	
	graphAsymmErrorsPi0DSS5000->Draw("same,E3");
	graphPi0CGC->Draw("same,p");
	graphPi0DSS5000->Draw("same,p,l");
	graphPi0ESP09sPi0AKK->Draw("same,p,l");
	graphPi0ESP09sPi0KKP->Draw("same,p,l");
	graphRpPbDalitzSystErr->Draw("2p,same,e");
	graphRpPbDalitzStatErr->Draw("p,same,e");
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	

	canvasRpPbESP09sDalitz->Update(); 
	
	canvasRpPbESP09sDalitz->SaveAs(Form("%s/RpPb_Comparison_Dalitz_Models.%s",outputDir.Data(),suffix.Data())); //NOTE* /


	
	
	
	
	
	
	
	//////////////////////////Comparison to Charged Pions/////////////////////////////////
	cout<<"Comparing charged pions"<<endl;
	 
	// TFile*  fileRpPbChargedPions 	 = new TFile( fileNameRpPbChargePions.Data() );
	
	//TGraphAsymmErrors*  graphChargedPionPPSystErr   = new TGraphAsymmErrors(histo_pi_pp_SystErr);
	//TGraphAsymmErrors*  graphChargedPionPPStatErr   = new TGraphAsymmErrors(histo_pi_pp_StatErr);
	
	
	
	
	
	
	
	
	
	/*
	
	TGraphAsymmErrors* graphErrosInterPolation5023GeVBinPCMStatErrScaled  =   (TGraphAsymmErrors*)graphErrosInterPolation5023GeVBinPCMStatErr->Clone();
	TGraphAsymmErrors* graphErrosInterPolation5023GeVBinPCMSystErrScaled  =   (TGraphAsymmErrors*)graphErrosInterPolation5023GeVBinPCMSystErr->Clone();
	
	TGraphAsymmErrors* graphErrosInterPolation5023GeVBinDalitzStatErrScaled  =   (TGraphAsymmErrors*)graphErrosInterPolation5023GeVBinDalitzStatErr->Clone();
	TGraphAsymmErrors* graphErrosInterPolation5023GeVBinDalitzSystErrScaled  =   (TGraphAsymmErrors*)graphErrosInterPolation5023GeVBinDalitzSystErr->Clone();
	
	
		
	graphErrosInterPolation5023GeVBinPCMStatErrScaled  = ScaleGraph(graphErrosInterPolation5023GeVBinPCMStatErrScaled, (1./(recalcBarn*1e-3)));
	graphErrosInterPolation5023GeVBinPCMSystErrScaled  = ScaleGraph(graphErrosInterPolation5023GeVBinPCMSystErrScaled, (1./(recalcBarn*1e-3)));
	
	
	graphErrosInterPolation5023GeVBinDalitzStatErrScaled  = ScaleGraph(graphErrosInterPolation5023GeVBinDalitzStatErrScaled, (1./(recalcBarn*1e-3)));
	graphErrosInterPolation5023GeVBinDalitzSystErrScaled  = ScaleGraph(graphErrosInterPolation5023GeVBinDalitzSystErrScaled, (1./(recalcBarn*1e-3)));
	
	*/
	
	
	//graphChargedPionPPStatErr = ScaleGraph(graphChargedPionPPStatErr,0.5);
	//graphChargedPionPPSystErr = ScaleGraph(graphChargedPionPPSystErr,0.5);
	
	
	
	cout<<"//////////////////////////Ratio to pp reference of charged pions////////////////////////////////"<<endl;
	
	
	
	//TGraphErrors* graphRatioPPReferencePCM_ChargedPions    = CalculateRatioBetweenSpectraWithDifferentBinning(graphErrosInterPolation5023GeVBinPCMStatErrScaled, graphErrosInterPolation5023GeVBinPCMSystErrScaled,graphChargedPionPPStatErr,graphChargedPionPPSystErr,  kTRUE,  kTRUE, &graphRebinnedAStat,&graphRebinnedASyst,&graphRebinnedBStat,&graphRebinnedBSyst);
        //TGraphErrors* graphRatioPPReferenceDalitz_ChargedPions = CalculateRatioBetweenSpectraWithDifferentBinning(graphErrosInterPolation5023GeVBinDalitzStatErr, graphErrosInterPolation5023GeVBinDalitzSystErr, graphChargedPionPPStatErr,graphChargedPionPPSystErr,       kTRUE,  kTRUE, &graphRebinnedAStat,&graphRebinnedASyst,&graphRebinnedBStat,&graphRebinnedBSyst);	
	
	//graphErrosInterPolation5023GeVBinPCMStatErrScaled->Print();
	
	//graphChargedPionPPStatErr->Print();
	
	
	/*
	
	TCanvas* cInvYieldPPComparisonPCM_ChargedPions = new TCanvas("cInvYieldPPComparisonPCM_ChargedPions","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldPPComparisonPCM_ChargedPions,  0.15, 0.02, 0.03, 0.1);
	cInvYieldPPComparisonPCM_ChargedPions->SetLogx();
	cInvYieldPPComparisonPCM_ChargedPions->SetLogy();
	
	TPad* padInvYieldPPComparisonPCM_ChargedPions = new TPad("padInvYieldPPComparisonPCM_ChargedPions", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPPComparisonPCM_ChargedPions, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPPComparisonPCM_ChargedPions->Draw();
	
	
	TH2F * histoInvYieldPPComparisonPCM_ChargedPions;
	histoInvYieldPPComparisonPCM_ChargedPions = new TH2F("histoInvYieldPPComparisonPCM_ChargedPions","histoInvYieldPPComparisonPCM_ChargedPions",1000,0.23,30.,1000,2e-12,10e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPPComparisonPCM_ChargedPions, "p_{T} (GeV/c)", "E#frac{d^{2}#sigma}{dp^{3}}(pb GeV^{-2} c^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPPComparisonPCM_ChargedPions->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphErrosInterPolation5023GeVBinPCMStatErrScaled,     	20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionPPStatErr, 	20, 1, kGreen+2, kGreen+2);
	
	
	
	
	
	graphErrosInterPolation5023GeVBinPCMStatErrScaled->Draw("p,same,e1");
	graphChargedPionPPStatErr->Draw("p,same,e1");
	
	
	
	LabelpPbPCM->Draw();
	labelSpectraPi0LabelPCMpPb2->Draw();

	labelSpectraPi0LabelpPbPCM->Draw();  

	graphErrosInterPolation5023GeVBinPCMStatErrScaled->SetFillColor(0);
	graphChargedPionPPStatErr->SetFillColor(0);
	graphErrosInterPolation5023GeVBinPCMScaled->SetFillColor(0);
	TLegend* legendInvYieldPPComparisonPCM_ChargedPions = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYieldPPComparisonPCM_ChargedPions->SetFillColor(0);
	legendInvYieldPPComparisonPCM_ChargedPions->SetLineColor(0);
	legendInvYieldPPComparisonPCM_ChargedPions->SetNColumns(1);
	legendInvYieldPPComparisonPCM_ChargedPions->SetTextSize(0.03);
	legendInvYieldPPComparisonPCM_ChargedPions->AddEntry(graphErrosInterPolation5023GeVBinPCMStatErrScaled,"#pi^{0} pp reference","pef");
	legendInvYieldPPComparisonPCM_ChargedPions->AddEntry(graphChargedPionPPStatErr,"#pi^{#pm} pp reference ","pef");
	
	
	legendInvYieldPPComparisonPCM_ChargedPions->Draw();

		
	
	cInvYieldPPComparisonPCM_ChargedPions->Update();
	cInvYieldPPComparisonPCM_ChargedPions->SaveAs(Form("%s/InvYield_Comparison_PP_PCM_ChargedPions.%s",outputDir.Data(),suffix.Data()));
	
	
	
	TCanvas* cInvYieldPPComparisonDalitz_ChargedPions = new TCanvas("cInvYieldPPComparisonDalitz_ChargedPions","Comparison Invariant Yields",200,10,700,500);
		
	DrawGammaCanvasSettings( cInvYieldPPComparisonDalitz_ChargedPions,  0.15, 0.02, 0.03, 0.1);
	cInvYieldPPComparisonDalitz_ChargedPions->SetLogx();
	cInvYieldPPComparisonDalitz_ChargedPions->SetLogy();
	
	TPad* padInvYieldPPComparisonDalitz_ChargedPions = new TPad("padInvYieldPPComparisonDalitz_ChargedPions", "", 0., 0.42, 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padInvYieldPPComparisonDalitz_ChargedPions, 0.15, 0.02, 0.03, 0.1);
	padInvYieldPPComparisonDalitz_ChargedPions->Draw();
	
	
	TH2F * histoInvYieldPPComparisonDalitz_ChargedPions;
	histoInvYieldPPComparisonDalitz_ChargedPions = new TH2F("histoInvYieldPPComparisonDalitz_ChargedPions","histoInvYieldPPComparisonDalitz_ChargedPions",1000,0.23,30.,1000,2e-12,10e1);
	SetStyleHistoTH2ForGraphs(histoInvYieldPPComparisonDalitz_ChargedPions, "p_{T} (GeV/c)", "E#frac{d^{2}#sigma}{dp^{3}}(pb GeV^{-2} c^{3})", 0.032,0.04, 0.04,0.04, 1,1.55);
	histoInvYieldPPComparisonDalitz_ChargedPions->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphAsym(graphErrosInterPolation5023GeVBinDalitzStatErrScaled,     	20, 1, kRed+1,   kRed+1);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionPPStatErr, 	20, 1, kGreen+2, kGreen+2);
	
	
	
	
	
	graphErrosInterPolation5023GeVBinDalitzStatErrScaled->Draw("p,same,e1");
	graphChargedPionPPStatErr->Draw("p,same,e1");
	
	
	
	
	
	labelSpectraPi0LabelPCMpPb->Draw();
	
	labelSpectraPi0LabelpPb->Draw();
	
	
	

	graphErrosInterPolation5023GeVBinDalitzStatErrScaled->SetFillColor(0);
	graphChargedPionPPStatErr->SetFillColor(0);
	graphErrosInterPolation5023GeVBinDalitzScaled->SetFillColor(0);
	TLegend* legendInvYieldPPComparisonDalitz_ChargedPions = new TLegend(0.2,0.15,0.6,0.35);
	legendInvYieldPPComparisonDalitz_ChargedPions->SetFillColor(0);
	legendInvYieldPPComparisonDalitz_ChargedPions->SetLineColor(0);
	legendInvYieldPPComparisonDalitz_ChargedPions->SetNColumns(1);
	legendInvYieldPPComparisonDalitz_ChargedPions->SetTextSize(0.03);
	legendInvYieldPPComparisonDalitz_ChargedPions->AddEntry(graphErrosInterPolation5023GeVBinDalitzStatErrScaled,"#pi^{0} pp reference","pef");
	legendInvYieldPPComparisonDalitz_ChargedPions->AddEntry(graphChargedPionPPStatErr,"#pi^{#pm} pp reference ","pef");
	
	
	legendInvYieldPPComparisonDalitz_ChargedPions->Draw();

		
	
	cInvYieldPPComparisonDalitz_ChargedPions->Update();
	cInvYieldPPComparisonDalitz_ChargedPions->SaveAs(Form("%s/InvYield_Comparison_PP_Dalitz_ChargedPions.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	
	
  
	
	TCanvas* canvasRatioPCM_ChargedPions = new TCanvas("canvasRatioPCM_ChargedPions","PP reference PCM compared to ChargedPions",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioPCM_ChargedPions,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioPCM_ChargedPions = new TH2F("histo2DRatioPCM_ChargedPions","histo2DRatioPCM_ChargedPions",1000,0.,15.,1000,-0.05,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioPCM_ChargedPions, "p_{T} (GeV/c)","R_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioPCM_ChargedPions->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphRatioPPReferencePCM_ChargedPions,21,1.5, kBlack , kBlack);
	graphRatioPPReferencePCM_ChargedPions->Draw("E1psame");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioPCM_ChargedPions = new TLatex(0.16,0.9,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRatioPCM_ChargedPions,0.06,4); 
	textRatioPCM_ChargedPions->Draw();
	
	
	TLegend* legendRatioPCM_ChargedPions = new TLegend(0.18,0.15,0.9,0.21);
	legendRatioPCM_ChargedPions->SetFillColor(0);
	legendRatioPCM_ChargedPions->SetLineColor(0);
	legendRatioPCM_ChargedPions->SetNColumns(2);
	legendRatioPCM_ChargedPions->SetTextSize(0.045);
	legendRatioPCM_ChargedPions->AddEntry(graphRatioPPReferencePCM_ChargedPions,"#pi^{0} pp reference (PCM) / #pi^{#pm} pp reference","p");
	legendRatioPCM_ChargedPions->Draw();
	
	
	canvasRatioPCM_ChargedPions->SaveAs(Form("%s/Ratio_PCM_ChargedPions.%s",outputDir.Data(),suffix.Data()));
	
	
	
		
	TCanvas* canvasRatioDalitz_ChargedPions = new TCanvas("canvasRatioDalitz_ChargedPions","PP reference Dalitz compared to ChargedPions",200,10,1200,700);
	
	DrawGammaCanvasSettings( canvasRatioDalitz_ChargedPions,  0.1, 0.01, 0.015, 0.13);
	
	TH2F * histo2DRatioDalitz_ChargedPions = new TH2F("histo2DRatioDalitz_ChargedPions","histo2DRatioDalitz_ChargedPions",1000,0.,15.,1000,-0.05,1.5);
	SetStyleHistoTH2ForGraphs(histo2DRatioDalitz_ChargedPions, "p_{T} (GeV/c)","R_{pPb}", 0.05,0.064, 0.05,0.06, 0.8,0.6, 512, 505); 
	histo2DRatioDalitz_ChargedPions->DrawCopy(); 
	
	DrawGammaSetMarkerTGraphErr(graphRatioPPReferenceDalitz_ChargedPions,21,1.5, kBlack , kBlack);
	graphRatioPPReferenceDalitz_ChargedPions->Draw("E1psame");
	
	
		
	DrawGammaLines(0., 15.,1., 1.,2.,kBlack);
	
	TLatex *textRatioDalitz_ChargedPions = new TLatex(0.16,0.9,"p-Pb minimum bias #sqrt{#it{s}_{NN}} = 5.02 TeV");
	SetStyleTLatex( textRatioDalitz_ChargedPions,0.06,4); 
	textRatioDalitz_ChargedPions->Draw();
	
	
	TLegend* legendRatioDalitz_ChargedPions = new TLegend(0.18,0.15,0.9,0.21);
	legendRatioDalitz_ChargedPions->SetFillColor(0);
	legendRatioDalitz_ChargedPions->SetLineColor(0);
	legendRatioDalitz_ChargedPions->SetNColumns(2);
	legendRatioDalitz_ChargedPions->SetTextSize(0.045);
	legendRatioDalitz_ChargedPions->AddEntry(graphRatioPPReferenceDalitz_ChargedPions,"#pi^{0} pp reference (Dalitz) / #pi^{#pm} pp reference","p");
	legendRatioDalitz_ChargedPions->Draw();
	
	
	canvasRatioDalitz_ChargedPions->SaveAs(Form("%s/Ratio_Dalitz_ChargedPions.%s",outputDir.Data(),suffix.Data()));
	
	*/
	
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbSystErr,20,1.5, kRed, kRed, 1, kTRUE,kRed-10);

	
	
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbStatErr,20,1.5, kRed, kRed);
	graphAsymmChargedParticlesRpPbStatErr->SetFillColor(0);
	
	
	
	
	 
	
	
	TCanvas* canvasRpPbChargedPionsPCM = new TCanvas("canvasRpPbChargedPionsPCM","Dalitz R_{pPb} with charged particles",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbChargedPionsPCM,  0.1, 0.01, 0.015, 0.11);
	
	
	TH2F * histo2DRpPbChargedPionsPCM = new TH2F("histo2DRpPbChargedPionsPCM","histo2DRpPbChargedPionsPCM",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedPionsPCM,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbChargedPionsPCM->GetXaxis()->SetRangeUser(0.0,14.5);
	histo2DRpPbChargedPionsPCM->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbChargedPionsPCM->GetYaxis()->CenterTitle();
	
	histo2DRpPbChargedPionsPCM->Draw(); 
	
	
	graphChargedPionRpPbSystErr->Draw("2p,same,e");
	graphChargedPionRpPbStatErr->Draw("p,same,e");
	
	
	graphRpPbPCMSystErr->Draw("2p,same,e");
	graphRpPbPCMStatErr->Draw("p,same,e");
	
	
	
	
	DrawGammaLines(0., 14.5,1., 1.,1.,kBlack,2);
	
	
	
	
	
	TLatex *textRpPbChargedPionsPCM = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbChargedPionsPCM,0.04,4); 
	textRpPbChargedPionsPCM->Draw();
	
	
	

	  
	
	TLatex *LabelpPbChargedPionsPCM = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbChargedPionsPCM, 0.04,4);
	LabelpPbChargedPionsPCM->Draw();
	
	
	

	
	TLegend* legendRpPbChargedPionsPCM = new TLegend(0.3,0.20,0.75,0.35);
	legendRpPbChargedPionsPCM->SetFillColor(0);
	legendRpPbChargedPionsPCM->SetLineColor(0);
	legendRpPbChargedPionsPCM->SetTextSize(0.03);
	legendRpPbChargedPionsPCM->SetTextFont(42);
	legendRpPbChargedPionsPCM->SetTextSize(0.04);
	
	legendRpPbChargedPionsPCM->AddEntry(graphRpPbPCMStatErr,Form("%s, -0.8 < #it{y} < 0.8",PCMLabel.Data()),"p");
	legendRpPbChargedPionsPCM->AddEntry(graphChargedPionRpPbStatErr,Form("%s, -0.5 < #it{y}_{cms} < 0   for #it{p}_{T} < 2.0 GeV/#it{c}",ChargedPionsLabel.Data()),"p");
	legendRpPbChargedPionsPCM->AddEntry((TObject*)0,"            -0.3 < #it{y}_{cms} < 0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	legendRpPbChargedPionsPCM->Draw();
	
	
	
	canvasRpPbChargedPionsPCM->SaveAs(Form("%s/RpPb_Comparison_PCM_ChargedPions.%s",outputDir.Data(),suffix.Data()));
	
	
	
	
	
	
	
	TCanvas* canvasRpPbChargedPionsDalitz = new TCanvas("canvasRpPbChargedPionsDalitz","R_{pPb} Dalitz",200,10,1000,700);
	
	DrawGammaCanvasSettings( canvasRpPbChargedPionsDalitz,  0.1, 0.01, 0.015, 0.11);
	
	
	TH2F * histo2DRpPbChargedPionsDalitz = new TH2F("histo2DRpPbChargedPionsDalitz","histo2DRpPbChargedPionsDalitz",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedPionsDalitz,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbChargedPionsDalitz->GetXaxis()->SetRangeUser(0.0,10.5);
	histo2DRpPbChargedPionsDalitz->GetYaxis()->SetRangeUser(0.19,1.81);
	histo2DRpPbChargedPionsDalitz->GetYaxis()->CenterTitle();
	
	histo2DRpPbChargedPionsDalitz->Draw(); 
	
	
	
	
	graphChargedPionRpPbSystErr->Draw("2p,same,e");
	graphChargedPionRpPbStatErr->Draw("p,same,e");
	
	graphRpPbDalitzSystErr->Draw("2p,same,e");
	graphRpPbDalitzStatErr->Draw("p,same,e");
	
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	
	
	
	
	TLatex *textRpPbChargedPionsDalitz = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbChargedPionsDalitz,0.04,4); 
	textRpPbChargedPionsDalitz->Draw();
	
	
	

	
	
	TLatex *LabelpPbChargedPionsDalitz = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbChargedPionsDalitz, 0.04,4);
	LabelpPbChargedPionsDalitz->Draw();
	
	

	
	TLegend* legendRpPbChargedPionsDalitz = new TLegend(0.3,0.20,0.75,0.35);
	legendRpPbChargedPionsDalitz->SetFillColor(0);
	legendRpPbChargedPionsDalitz->SetLineColor(0);
	legendRpPbChargedPionsDalitz->SetTextSize(0.03);
	legendRpPbChargedPionsDalitz->SetTextFont(42);
	legendRpPbChargedPionsDalitz->SetTextSize(0.04);
	
	legendRpPbChargedPionsDalitz->AddEntry(graphRpPbDalitzStatErr,Form("%s, -0.8 < #it{y}_{cms} < 0.8",DalitzLabel.Data()),"p");
	legendRpPbChargedPionsDalitz->AddEntry(graphChargedPionRpPbStatErr,Form("%s, -0.5 < #it{y}_{cms} < 0   for #it{p}_{T} < 2.0 GeV/#it{c}",ChargedPionsLabel.Data()),"p");
	legendRpPbChargedPionsDalitz->AddEntry((TObject*)0,"            -0.3 < #it{y}_{cms} < 0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	legendRpPbChargedPionsDalitz->Draw();
	
	
	
	
	
	canvasRpPbChargedPionsDalitz->SaveAs(Form("%s/RpPb_Comparison_Dalitz_ChargedPions.%s",outputDir.Data(),suffix.Data()));
	
	
	
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRpPbSystErr,20,1.5, colorsArray[kKaon],colorsArray[kKaon],1.0,kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphChargedKaonRpPbStatErr,20,1.5, colorsArray[kKaon],colorsArray[kKaon]);
	graphChargedKaonRpPbStatErr->SetFillColor(0);
	
	
	DrawGammaSetMarkerTGraphAsym(graphChargedProtonRpPbSystErr,20,1.5, colorsArray[kProton], colorsArray[kProton], 1,kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphChargedProtonRpPbStatErr,20,1.5, colorsArray[kProton], colorsArray[kProton]);
	graphChargedProtonRpPbStatErr->SetFillColor(0);
	
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbSystErrClone,20,1.5, colorsArray[kPion], colorsArray[kPion], 1,kTRUE);
	DrawGammaSetMarkerTGraphAsym(graphChargedPionRpPbStatErrClone,20,1.5, colorsArray[kPion], colorsArray[kPion]);
	graphChargedPionRpPbStatErrClone->SetFillColor(0);
	
	
	
	
	
	TCanvas* canvasRpPbChargedPionsKaonsProtonsDalitz = new TCanvas("canvasRpPbChargedPionsKaonsProtonsDalitz","R_{pPb} Dalitz",200,10,1000,800);
	DrawGammaCanvasSettings( canvasRpPbChargedPionsKaonsProtonsDalitz,  0.1, 0.01, 0.015, 0.11);

	TH2F * histo2DRpPbChargedPionsKaonsProtonsDalitz = new TH2F("histo2DRpPbChargedPionsKaonsProtonsDalitz","histo2DRpPbChargedPionsKaonsProtonsDalitz",1000,0.,20.0,1000,0.0,5.0);
	SetStyleHistoTH2ForGraphs(histo2DRpPbChargedPionsKaonsProtonsDalitz,pTLabel.Data(),RpPbLabel.Data(), 0.04,0.04, 0.04,0.04, 1.2,1.0, 512, 512); 
	
	
	histo2DRpPbChargedPionsKaonsProtonsDalitz->GetXaxis()->SetRangeUser(0.0,10.5);
	histo2DRpPbChargedPionsKaonsProtonsDalitz->GetYaxis()->SetRangeUser(0.19,1.81);
	
	histo2DRpPbChargedPionsKaonsProtonsDalitz->Draw(); 
	
	
	
	
	graphChargedKaonRpPbSystErr->Draw("2p,same,e");
	graphChargedKaonRpPbStatErr->Draw("p,same,e");
	graphChargedProtonRpPbSystErr->Draw("2p,same,e");
	graphChargedProtonRpPbStatErr->Draw("p,same,e");
	graphChargedPionRpPbSystErrClone->Draw("2p,same,e");
	graphChargedPionRpPbStatErrClone->Draw("p,same,e");
	
	
	graphRpPbDalitzSystErr->Draw("2p,same,e");
	graphRpPbDalitzStatErr->Draw("p,same,e");
	
	
	DrawGammaLines(0., 10.5,1., 1.,1.,kBlack,2);
	
	
	
	
	
	TLatex *textRpPbChargedPionsKaonsProtonsDalitz = new TLatex(0.15,0.87,energyLabel.Data());
	SetStyleTLatex( textRpPbChargedPionsKaonsProtonsDalitz,0.04,4); 
	textRpPbChargedPionsKaonsProtonsDalitz->Draw();
	
	
	
	TLatex *LabelpPbChargedPionsKaonsProtonsDalitz = new TLatex(0.15,0.92,thesisPlotLabel.Data());
	SetStyleTLatex( LabelpPbChargedPionsKaonsProtonsDalitz, 0.04,4);
	LabelpPbChargedPionsKaonsProtonsDalitz->Draw();

	
	TLegend* legendRpPbChargedPionsKaonsProtonsDalitz = new TLegend(0.3,0.15,0.75,0.45);
	legendRpPbChargedPionsKaonsProtonsDalitz->SetFillColor(0);
	legendRpPbChargedPionsKaonsProtonsDalitz->SetLineColor(0);
	legendRpPbChargedPionsKaonsProtonsDalitz->SetTextSize(0.03);
	legendRpPbChargedPionsKaonsProtonsDalitz->SetTextFont(42);
	legendRpPbChargedPionsKaonsProtonsDalitz->SetTextSize(0.03);
	
	legendRpPbChargedPionsKaonsProtonsDalitz->AddEntry(graphRpPbDalitzStatErr,DalitzLabel.Data(),"p");
	legendRpPbChargedPionsKaonsProtonsDalitz->AddEntry(graphChargedPionRpPbStatErrClone,"#pi^{+}+#pi^{-} -0.5 < #it{y}_{CMS} < 0 for #it{p}_{T} < 2.0 GeV/#it{c}","p");
	legendRpPbChargedPionsKaonsProtonsDalitz->AddEntry((TObject*)0,"            -0.3 < #it{y}_{CMS}< 0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	legendRpPbChargedPionsKaonsProtonsDalitz->AddEntry(graphChargedKaonRpPbSystErr,"K^{+} + K^{-} -0.5 < #it{y}_{CMS} < 0 for #it{p}_{T} < 2.0 GeV/#it{c}","p");
	legendRpPbChargedPionsKaonsProtonsDalitz->AddEntry((TObject*)0,"            -0.3 < #it{y}_{CMS}< 0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	legendRpPbChargedPionsKaonsProtonsDalitz->AddEntry(graphChargedProtonRpPbSystErr,"p + #bar{p} -0.5 <#it{y}_{CMS}<0 for #it{p}_{T} < 2.0 GeV/#it{c}","p");
	legendRpPbChargedPionsKaonsProtonsDalitz->AddEntry((TObject*)0,"            -0.3<#it{y}_{CMS}<0.3 for #it{p}_{T} > 2.0 GeV/#it{c}","");
	
	legendRpPbChargedPionsKaonsProtonsDalitz->Draw();
	
	
	canvasRpPbChargedPionsKaonsProtonsDalitz->SaveAs(Form("%s/RpPb_Comparison_Dalitz_ChargedPionsKaonsProtons.%s",outputDir.Data(),suffix.Data()));

	
	
  
  
	
	
	
	
	SavingFiles(outputDir);
	
	
	
	cout<<"Llego al final"<<endl;
	
	


}

TH1F* GetPPReferenceFromPythia(TString fileName){
  
   TFile f1(fileName.Data());
   //TFile f1("pythia8-5.02TeV-part-1.root");
  // TFile f1("job_nd.12/pythia_pi0_non_diffractive.root");
   

   Double_t nEvt;
   Double_t xVal,yVal,yErr;
   
   TH1F* mult = (TH1F*)f1.Get("mult");
   TH1F* pi0PtVar = (TH1F*)f1.Get("pi0PtVar");
   TH1F* pi0Pt    = (TH1F*)f1.Get("pi0Pt");

   nEvt=mult->GetEntries()/1.e7;
   cout<<"Nevt="<< nEvt<<endl;
   
   
   pi0PtVar->Scale(1./(2*TMath::Pi()*nEvt));
   TH1F * OverPtpi0Pt=(TH1F*)pi0PtVar->Clone();
  //pi0Pt.Scale(1./(2*TMath::Pi()*nEvt));
  //   TH1F * OverPtpi0Pt=(TH1F*)pi0Pt->Clone();
     //Double_t sigma=6.789e+01;//48.03e-3;// cross section from Pythia8-tune 4C in mbar
     //Double_t sigma= 6.789e-04;
     //Double_t sigma=1.;


   for(Int_t ii=0;ii<pi0PtVar->GetNbinsX()+1;ii++){
   //for(Int_t ii=0;ii<pi0Pt->GetNbinsX()+1;ii++){

    xVal=pi0PtVar->GetBinCenter(ii+1);
    yVal=pi0PtVar->GetBinContent(ii+1);
    yErr=pi0Pt->GetBinError(ii+1);
    //OverPtpi0Pt->SetBinContent(ii+1,(yVal/(xVal*sigma))*recalcBarn);
    //OverPtpi0Pt->SetBinContent(ii+1,(yVal/xVal) * (sigma * recalcBarn) ); 
    OverPtpi0Pt->SetBinContent(ii+1,(yVal/xVal) * 1e09 ); 
    //    cout<< xVal <<" "<< OverPtpi0Pt.GetBinContent(ii+1)<<endl;

    // OverPtpi0Pt.SetBinError(ii+1,yErr/xVal);

    /*
    xVal=pi0Pt.GetBinCenter(ii+1);
    yVal=pi0Pt.GetBinContent(ii+1);
    yErr=pi0Pt.GetBinError(ii+1);
    OverPtpi0Pt.SetBinContent(ii+1,yVal/(xVal*sigma));
    // OverPtpi0Pt.SetBinError(ii+1,yErr/xVal);
    cout<< xVal <<" "<< OverPtpi0Pt.GetBinContent(ii+1)<<endl;
    */

  }
//OverPtpi0Pt->SetLineColor(2);
//OverPtpi0Pt->Draw("same");


 //if( Parameters2760GeV ) delete Parameters2760GeV;
 //if( Parameters7TeV )    delete Parameters7TeV;
	
  //canvasRpPbCombine->SaveAs(Form("%s/RpPb_Comparison_PCM_Dalitz.%s",outputDir.Data(),suffix.Data()));
  
	
 return OverPtpi0Pt;
}



void SavingFiles(TString outputDir){
  
      TString dateForOutput = ReturnDateStringForOutput();
  
      TFile fCombResults(Form("%s/Pi0RpPb_%s.root",outputDir.Data(),dateForOutput.Data()),"RECREATE");
      
      graphErrosInterPolation5023GeVBinDalitzStatErr->Write("Pi0_pp_reference_DalitzBinning_StatErr");
      graphErrosInterPolation5023GeVBinDalitzSystErr->Write("Pi0_pp_reference_DalitzBinning_SystErr");
      graphErrosInterPolation5023GeVBinPCMStatErr->Write("Pi0_pp_reference_PCMBinning_StatErr");
      graphErrosInterPolation5023GeVBinPCMSystErr->Write("Pi0_pp_reference_PCMBinning_SystErr");
      graphRpPbDalitzSystErr->Write("Pi0_RpPb_Dalitz_SystErr");
      graphRpPbDalitzStatErr->Write("Pi0_RpPb_Dalitz_StatErr");
      graphRpPbPCMStatErr->Write("Pi0_RpPb_PCM_StatErr");
      graphRpPbPCMSystErr->Write("Pi0_RpPb_PCM_SystErr");
      
      if( graphAlphaDalitz )  graphAlphaDalitz->Write("Pi0_RpPb_Dalitz_Alpha");
      if( graphAlphaPCM )     graphAlphaPCM->Write("Pi0_RpPb_PCM_Alpha");
       
      if ( graphPtvsSqrtsDalitz ){
      
	for ( Int_t i = 0;  i < graphAlphaDalitz->GetN(); i++)
	  if( graphPtvsSqrtsDalitz[i] ) graphPtvsSqrtsDalitz[i]->Write(Form("Pt_vs_Sqrts_Dalitz_Bin_%d",i));
      }
      
      if ( graphPtvsSqrtsPCM ){
      
	for ( Int_t i = 0;  i < graphAlphaPCM->GetN(); i++)
	  if( graphPtvsSqrtsPCM[i] ) graphPtvsSqrtsPCM[i]->Write(Form("Pt_vs_Sqrts_PCM_Bin_%d",i));
      }
      
  
}



