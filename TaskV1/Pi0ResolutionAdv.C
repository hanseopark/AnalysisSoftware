/***********************************************************************************************
*** provided by Gamma Conversion Group, PWG4, 									******
***	    Friederike Bock, fbock@physi.uni-heidelberg.de ***							******
************************************************************************************************
************************************************************************************************
*** With this Macro Resolution Studies for the Conversion Method can be carried out 		******
*** For this a run with the flag "resolution" of the Gamma Conversion Software needs 	******
*** to be performed, this study can only be carried out on Montecarlo Productions 		******
************************************************************************************************/

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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

TString textGenerator;
TString collisionSystem;
TString textPeriod;
TString textDate;
TString mesonLatex;

void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Int_t logX, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500){
	TCanvas * canvasStandard = new TCanvas("canvasStandard","",10,10,canvasSizeX,canvasSizeY);  // gives the page size		
	canvasStandard->SetLogx(logX);
	canvasStandard->SetLogz(logZ);
	canvasStandard->SetRightMargin(0.12); 		
	canvasStandard->SetLeftMargin(0.12); 		
	canvasStandard->SetBottomMargin(0.1); 		
	canvasStandard->SetTopMargin(0.04); 		
	canvasStandard->cd();
	DrawAutoGammaHisto2D(	histo2D,
									title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX);
 	DrawAliceLogoPerformance(floatLogo[0],floatLogo[1],floatLogo[2],floatLogo[3],0.00, textDate,collisionSystem, textGenerator,textPeriod,canvasSizeX,canvasSizeY);	
	canvasStandard->Update();
	canvasStandard->SaveAs(nameOutput.Data());
	delete canvasStandard;
}




void ProduceProjectionsPlotWithFits( TH1D** histos, TF1** fitData, Double_t* floatArray, Int_t nBins, Int_t nColumns, Int_t nRows, TString xTitle, TString yTitle, TString legendEntry1, TString legendEntry2, TString title, TString nameOutput, Int_t canvasSizeX= 2800, Int_t canvasSizeY = 1800){
										
	TCanvas * canvasProjection = new TCanvas("canvasProjection","",canvasSizeX,canvasSizeY);  // gives the page size
	canvasProjection->SetTopMargin(0.02);
	canvasProjection->SetBottomMargin(0.02);
	canvasProjection->SetRightMargin(0.02);
	canvasProjection->SetLeftMargin(0.02);

	TPad * padProjection = new TPad("padProjection","",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
	padProjection->SetFillColor(0);
	padProjection->GetFrame()->SetFillColor(0);
	padProjection->SetBorderMode(0);
	padProjection->SetLogy(0);
	padProjection->Divide(nColumns,nRows);
	padProjection->Draw();

	TGaxis::SetMaxDigits(3);
	Double_t relWidthLogo=0.4;
	Double_t padXWidth = 2800/nColumns;
	Double_t padYWidth = 1800/nRows;
	
	Int_t place = 0;
	for(Int_t iPt=0;iPt<nBins;iPt++){
		Double_t startPt = floatArray[iPt];
		Double_t endPt = floatArray[iPt+1];

		place = place + 1;						//give the right place in the page
		if (place == nColumns){
			iPt--;
			padProjection->cd(place);
			padProjection->cd(place)->SetTopMargin(0.12);
			padProjection->cd(place)->SetBottomMargin(0.15);
			padProjection->cd(place)->SetRightMargin(0.05);
			padProjection->cd(place)->SetLeftMargin(0.15);

			string textAlice = "ALICE performance";
			Double_t textHeight = 0.07;
			Double_t startTextX = 0.0;
			Double_t startTextY = 0.65;
			Double_t differenceText = textHeight*1.15;
			Double_t coordinatesStartPadX = 0.3;
			Double_t coordinatesStartPadY = 0.0;
			Double_t coordinatesEndPadX = (coordinatesStartPadX*padXWidth+relWidthLogo*padXWidth)/padXWidth;
			Double_t coordinatesEndPadY = (coordinatesStartPadY*padYWidth+relWidthLogo*padXWidth)/padYWidth;

			TLatex *alice = new TLatex(startTextX,(startTextY+(2*differenceText)),Form("%s",textAlice.c_str()));
			TLatex *process = new TLatex(startTextX, (startTextY+differenceText), collisionSystem.Data());
			TLatex *latexDate = new TLatex(startTextX,startTextY,textDate.Data());
			TPad *myPadLogo = new TPad("myPadLogo", "Pad for ALICE Logo", coordinatesStartPadX, coordinatesStartPadY, coordinatesEndPadX ,coordinatesEndPadY, -1, -1, -2);

			if (textGenerator.CompareTo("")!=0 && textPeriod.CompareTo("")!=0){
				TLatex *generator = new TLatex(startTextX,(startTextY-differenceText),Form("%s   %s",textGenerator.Data(),textPeriod.Data())); // Bo: this was modified
				generator->SetNDC();
				generator->SetTextColor(1);
				generator->SetTextFont(62);
				generator->SetTextSize(textHeight);
				generator->SetLineWidth(2);
				generator->Draw("same");	
			} else if (textGenerator.CompareTo("")!=0) {
				TLatex *generator = new TLatex(startTextX,(startTextY-differenceText),Form("%s",textGenerator.Data())); // Bo: this was modified
				generator->SetNDC();
				generator->SetTextColor(1);
				generator->SetTextFont(62);
				generator->SetTextSize(textHeight);
				generator->SetLineWidth(2);
				generator->Draw("same");	
			}

			alice->SetNDC();
			alice->SetTextColor(2);
			alice->SetTextSize(textHeight);
			alice->Draw();

			process->SetNDC(kTRUE);
			process->SetTextSize(textHeight);
			process->Draw();

			latexDate->SetNDC();
			latexDate->SetTextColor(1);
			latexDate->SetTextSize(textHeight);
			latexDate->Draw();
			
			TLegend* legendProjections = new TLegend(0.55,0.5,1.,0.8);
			legendProjections->SetTextSize(0.08);			
			legendProjections->SetFillColor(0);
			legendProjections->SetBorderSize(0);
			legendProjections->AddEntry(histos[3],legendEntry1.Data(),"pe");
			legendProjections->AddEntry(fitData[3],legendEntry2.Data(),"l");
			legendProjections->Draw();

			myPadLogo->SetFillColor(0);
			myPadLogo->SetFrameFillColor(0);
			myPadLogo->SetBorderMode(0);
			myPadLogo->SetBorderSize(0);
			myPadLogo->SetFrameBorderMode(0);
			myPadLogo->SetLeftMargin(0.0);
			myPadLogo->SetTopMargin(0.0);
			myPadLogo->SetBottomMargin(0.0);
			myPadLogo->SetRightMargin(0.0);
			TASImage *myAliceLogo = new TASImage("ALICE_logo.eps");
			myPadLogo->Draw();
			myPadLogo->cd();
			myAliceLogo->Draw();
			
		} else {
			padProjection->cd(place);
			padProjection->cd(place)->SetTopMargin(0.12);
			padProjection->cd(place)->SetBottomMargin(0.15);
			padProjection->cd(place)->SetRightMargin(0.05);
			padProjection->cd(place)->SetLeftMargin(0.15);
			padProjection->cd(place)->SetLogy(0);
         Float_t startRange = -0.5;
 /*        Int_t i = 1;
         while(histos[iPt]->GetBinContent(i) == 0 && i < histos[iPt]->GetNbinsX()){
//             cout << i << "\t" << histos[iPt]->GetBinContent(i) << endl;
            i++;
         }
         startRange= histos[iPt]->GetBinCenter(i-2);
 */
         Float_t endRange = 0.5;
//          Int_t l = histos[iPt]->GetNbinsX();
//          while(histos[iPt]->GetBinContent(l) == 0 && l > 0){
// //             cout << l << "\t" << histos[iPt]->GetBinContent(l) << endl;
//             l--;
//          }
//         endRange= histos[iPt]->GetBinCenter(l+2);
         
         
			Double_t yMin = 1.;
			if ( histos[iPt]->GetEntries() > 0){
				if (title.CompareTo("Pt")==0){
					DrawGammaHistoWithTitle2( histos[iPt],
							Form("%3.2f GeV/c < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt), 
							xTitle.Data(),yTitle.Data(),
							startRange,endRange,yMin);
				}
				if (fitData[iPt]!=0x0){
					fitData[iPt]->SetLineColor(kBlue);
					fitData[iPt]->SetLineWidth(0.8);
					fitData[iPt]->Draw("same");
				}
			} else {
				place--;
			}
		}
	}
	canvasProjection->Print(nameOutput.Data());
	delete padProjection;
	delete canvasProjection;

}


void Pi0ResolutionAdv( TString mesonName                    = "Pi0",
                       const char *fileMonteCarloInput      = "", 
                       TString optionCutSelection           = "", 
                       TString suffix                       = "eps", 
                       TString optEnergy                    = "", 
                       TString optMCGenerator               = "",
                       TString optPeriod                    = "",
                       Int_t mode                           = 0){	
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettings();	
	SetPlotStyle();

	TString collisionSystem 		= ReturnFullCollisionsSystem(optEnergy);
	TString detectionSystem 		= ReturnFullTextReconstructionProcess(mode);
	if (collisionSystem.CompareTo("") == 0){
		cout << "No correct collision system specification, has been given" << endl;
		return;     
	}
	
	TString centrality = "";
	TString firstCutnumber = optionCutSelection(GetEventSystemCutPosition(),1);
	if (firstCutnumber.CompareTo("0") != 0){
		centrality = GetCentralityString(optionCutSelection);
	}   
	cout << centrality.Data() << endl;
   
   
	if(optMCGenerator.CompareTo("") ==0){
		textGenerator = "";
	} else {
		textGenerator = optMCGenerator;
	}

	textDate = ReturnDateString();
	
    if (mesonName.CompareTo("Pi0") == 0)
        mesonLatex ="#pi^{0}"; 
    else if (mesonName.CompareTo("Eta") == 0)
        mesonLatex ="#eta"; // 
	
	TString outputDirectory;
	if(optPeriod.CompareTo("") ==0 || optPeriod.CompareTo("All") ==0){
		textPeriod = "";
		outputDirectory =	 		Form("%s/%s/%s/Pi0Resolution",optionCutSelection.Data(),optEnergy.Data(),suffix.Data());
		gSystem->Exec("mkdir -p "+outputDirectory);
	} else {
		textPeriod = optPeriod;
		outputDirectory =	 		Form("%s/%s/%s/%s/Pi0Resolution",optionCutSelection.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
		gSystem->Exec("mkdir -p "+outputDirectory);
	}
	
	
//**********************************************************************************************************************
//******************************** Definition of some plotting variables ***********************************************
//**********************************************************************************************************************

	//Array defintion for printing Logo 
	
	
	Float_t floatLocationRightDown2D[4]=	{0.65,0.2,0.11, 0.02};
	
//**********************************************************************************************************************
//******************************** Defintion of arrays for rebinning ***************************************************
//**********************************************************************************************************************

	//pt rebinning array
	Double_t ptbinning[19] = {0.4  ,0.6  ,0.8  ,1.   ,1.2,
                             1.4  ,1.6 , 1.8, 2.   ,2.4  ,
                             2.8, 3.2   ,3.6  ,4.,  4.5  ,
                             5. , 6.    ,10., 15.};
	Int_t ptbinningReb[18] = {4  ,4  ,4  ,4 ,4,
                             4  ,4  ,4 ,4  ,4  ,
                             4, 4  ,4  ,4  ,8  ,
                             8,10  ,10 };	
	
//**********************************************************************************************************************
//**************** Definition of maximum and minimum values for histogramms depending on the eta cutnumber *************
//**********************************************************************************************************************

	Double_t arrayPtMeanPt[2] = {-2.2,2.2};

	TLatex *latexEtaRange;
	latexEtaRange = 	new TLatex(0.15,0.92,"|#eta| < 0.9 "); // Bo: this was modified
	latexEtaRange->SetNDC();
	latexEtaRange->SetTextFont(62);
	latexEtaRange->SetTextSize(0.04);
	latexEtaRange->SetLineWidth(6);      
	

//**********************************************************************************************************************
//****************************************** Loading of Histograms *****************************************************
//**********************************************************************************************************************
   
	TFile fileMC(fileMonteCarloInput);  
	
	//************************** Container Loading ********************************************************************
	TString nameOutputContainer = "GammaConvV1";
	if (mode == 2 || mode == 3) nameOutputContainer = "GammaConvCalo";
		else if (mode == 4 || mode == 5) nameOutputContainer = "GammaCalo";
		
	TList *TopDir =(TList*)fileMC.Get(nameOutputContainer.Data());
	if(TopDir == NULL){
		cout<<"ERROR: TopDir not Found"<<endl;
		return;
	}
	TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",optionCutSelection.Data()));
	//    TList *MCContainer = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",optionCutSelection.Data()));
	TList *TrueConversionContainer = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",optionCutSelection.Data()));

		//************************** Histo loading *************************************************************************

		//****************************** Resolution dPt vs Pt **********************************************
		TH2F * Resolution_Pi0_MCPt_ResolPt = (TH2F*)TrueConversionContainer->FindObject(Form("ESD_TruePrimary%s_MCPt_ResolPt",mesonName.Data()) ); //
	cout << "hier noch " << Resolution_Pi0_MCPt_ResolPt<<endl;
	//**********************************************************************************************************************************************
//*********************************** Rebinning of 2D histograms and fitslices *****************************************************************
//**********************************************************************************************************************************************
	Int_t maxYbin = 0;
// 	Int_t maxXbin = 0;
	Double_t precision = 10E-4;
		
		TH2F *Resolution_Pi0_Ptrebin = new TH2F("Resolution_Pi0_Ptrebin", "Photon Resolution dPt vs Pt ", 18,ptbinning, Resolution_Pi0_MCPt_ResolPt->GetNbinsY() ,-1,1) ;	
		cout << "n bins Y  " << Resolution_Pi0_MCPt_ResolPt->GetNbinsY() << endl;
		Resolution_Pi0_MCPt_ResolPt->Sumw2();
		maxYbin = Resolution_Pi0_MCPt_ResolPt->GetNbinsY();
		for (Int_t i = 0; i <Resolution_Pi0_MCPt_ResolPt->GetNbinsX(); i++) {
				TAxis *xaxis_old = 	Resolution_Pi0_MCPt_ResolPt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 0 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_Pi0_MCPt_ResolPt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_Pi0_MCPt_ResolPt->GetBinContent(i, j);
					Resolution_Pi0_Ptrebin->Fill(binoldpoint, biny, value);
				}
		}
		cout << "hier noch" << endl;
// 		for (Int_t i = 0; i < Resolution_Pi0_Ptrebin->GetNbinsX(); i++){
// 			for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
// 				Double_t value =  Resolution_Pi0_Ptrebin->GetBinContent(i, j);
// 				Resolution_Pi0_Ptrebin->SetBinContent(i, j, value);	/// Resolution_Pi0_MCPt_ResolPt->GetYaxis()->GetBinWidth(j)
// 			}	
// 		}
		cout << "hier noch" << endl;
		TH1D* Resolution_Pi0_PtRes[18];
		TF1* 	fitResolution_Pi0_PtRes[18];
		TH1F *Resolution_Pi0_Ptrebin_1 = new TH1F("Resolution_Pi0_Ptrebin_1", "mean Gamma Resolution dPt vs Pt", 18, ptbinning) ;	
		TH1F *Resolution_Pi0_Ptrebin_2 = new TH1F("Resolution_Pi0_Ptrebin_2", "sigma Gamma Resolution dPt vs Pt", 18, ptbinning) ;	

		ResolutionFittingRebined( Resolution_Pi0_Ptrebin, Resolution_Pi0_PtRes , fitResolution_Pi0_PtRes , 18,ptbinningReb, Resolution_Pi0_Ptrebin_1 ,Resolution_Pi0_Ptrebin_2, "gaus",-0.3,0.3, precision , "Resolution_Pi0_PtRes");

		
//  		ResolutionFittingNormalized( Resolution_Pi0_Ptrebin, Resolution_Pi0_PtRes , fitResolution_Pi0_PtRes , deltaPt_Pi0Res, 29, Resolution_Pi0_Ptrebin_1 ,Resolution_Pi0_Ptrebin_2, "gaus",-1.,1., precision , "Resolution_Pi0_PtRes");

		for (Int_t i = 0; i < Resolution_Pi0_Ptrebin_2->GetNbinsX(); i++){
			Double_t newValue = Resolution_Pi0_Ptrebin_2->GetBinContent(i)*100;
			Double_t newError = Resolution_Pi0_Ptrebin_2->GetBinError(i)*100;
// 			Resolution_Pi0_Ptrebin_1->SetBinContent(i, (Resolution_Pi0_Ptrebin_1->GetBinContent(i)-Resolution_Pi0_Ptrebin_1->GetBinCenter(i))/Resolution_Pi0_Ptrebin_2->GetBinCenter(i)*100);
// 			Resolution_Pi0_Ptrebin_1->SetBinError(i, Resolution_Pi0_Ptrebin_1->GetBinError(i)/Resolution_Pi0_Ptrebin_2->GetBinCenter(i)*100);
// 			Double_t newValue1 = Resolution_Pi0_Ptrebin_1->GetBinContent(i)*100;
// 			Double_t newError1 = Resolution_Pi0_Ptrebin_1->GetBinError(i)*100;
			Resolution_Pi0_Ptrebin_2->SetBinContent(i, newValue);
			Resolution_Pi0_Ptrebin_2->SetBinError(i, newError);
		}
   Resolution_Pi0_Ptrebin_1->Scale(100);
		
//********************************************************************************************************
//*************************************** creating ps-file ***********************************************
//********************************************************************************************************  
	
//****************************************************************************************************************
//************************************ single plotting for better resolution *************************************
//****************************************************************************************************************
	TCanvas * canvasEPGdPtVsPt = new TCanvas("canvasEPGdPtVsPt","",10,10,700,1000);  // gives the page size		
	canvasEPGdPtVsPt->cd();
	DrawGammaCanvasSettings( canvasEPGdPtVsPt, 0.13, 0.02, 0.02, 0.09);
	
	TPad* padEPGdPtVsPt1 = new TPad("padEPGdPtVsPt1", "", 0., 0.5, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padEPGdPtVsPt1, 0.12, 0.02, 0.02, 0.);
	padEPGdPtVsPt1->Draw();
	
	TPad* padEPGdPtVsPt2 = new TPad("padEPGdPtVsPt2", "", 0., 0., 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padEPGdPtVsPt2, 0.12, 0.02, 0., 0.13);
	padEPGdPtVsPt2->Draw();
	
	padEPGdPtVsPt1->cd();
	TH2F * histo2DMean = new TH2F("histo2DMean","histo2DMean",1000,0.,15.,1000,-8.5,7.5);
	SetStyleHistoTH2ForGraphs(histo2DMean, "#it{p}_{T,MC} (GeV/#it{c})", Form (" #mu ((#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}) (%)",mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()),  0.04,0.05, 0.04,0.05, 1.1,1.1);
	if (mode == 0){
		histo2DMean->GetYaxis()->SetRangeUser(-1.9,1.9);
		histo2DMean->GetXaxis()->SetRangeUser(0,6);
	} else if (mode == 2 || mode == 3){
		histo2DMean->GetYaxis()->SetRangeUser(-5.4,2.4);
		histo2DMean->GetXaxis()->SetRangeUser(0,10);
	} else if (mode == 4 || mode == 5){
        histo2DMean->GetYaxis()->SetRangeUser(-8.5,6.5);
        histo2DMean->GetXaxis()->SetRangeUser(0,10);        
    }    
	histo2DMean->DrawCopy(); 

		StylingSliceHistos(Resolution_Pi0_Ptrebin_1,0.8);
		Resolution_Pi0_Ptrebin_1->SetMarkerColor(kBlue+2);
		Resolution_Pi0_Ptrebin_1->SetLineColor(kBlue-8);
		Resolution_Pi0_Ptrebin_1->Draw("same,pe");
// 		DrawResolutionGammaHisto( Resolution_Pi0_Ptrebin_1, 
// 						"", "p_{T,MC} (GeV/c)", " #mu ((p^{#pi^{0}}_{T,rec} -p^{#pi^{0}}_{T,MC}))/p^{#pi^{0}}_{T,MC} (%)", 
// 							kFALSE, 10., 140.,
// 							kFALSE, arrayPtMeanPt[0], arrayPtMeanPt[1],
// 							kTRUE, 0., 5.9);
      
// 		TF1 *fCorr  = new TF1("line","[0]+[1]*x",0.3,16.); 
// 		Resolution_Pi0_Ptrebin_1->Fit(fCorr,"RME0");
// 		Double_t parameterProb[2];
// 		fCorr->GetParameters(parameterProb);
// 		TLine *line = new TLine(0.3,parameterProb[0]+0.3*parameterProb[1],16.,parameterProb[0]+16*parameterProb[1]);
// 		line->SetLineWidth(1);
// 		line->Draw();
		TLatex *labelCentrality1 = new TLatex(0.5,0.92,Form("%s %s",centrality.Data(), collisionSystem.Data()  ));
		SetStyleTLatex( labelCentrality1, 0.038,4);
		labelCentrality1->Draw();
		TLatex *labelProcess = new TLatex(0.51,0.88,Form("%s #rightarrow #gamma#gamma",mesonLatex.Data() ));
		SetStyleTLatex( labelProcess, 0.038,4);
		labelProcess->Draw();
		TLatex *labelDetSys = new TLatex(0.51,0.84,detectionSystem.Data());
		SetStyleTLatex( labelDetSys, 0.038,4);
		labelDetSys->Draw();
		
	padEPGdPtVsPt1->Update();
	padEPGdPtVsPt2->cd();
   
	TH2F * histo2DSigma = new TH2F("histo2DSigma","histo2DSigma",1000,0.,15.,1000,-0.1,100);
	SetStyleHistoTH2ForGraphs(histo2DSigma, Form("#it{p}^{%s}_{T,MC} (GeV/#it{c})",mesonLatex.Data()), Form(" #sigma #it{p}^{%s}_{T} (%)",mesonLatex.Data()), 0.04,0.05, 0.04,0.05, 1.1,1.1);
	if (mode == 0){
		histo2DSigma->GetYaxis()->SetRangeUser(0,7.2);
		histo2DSigma->GetXaxis()->SetRangeUser(0,6);
	} else if (mode == 2 || mode == 3){
		histo2DSigma->GetYaxis()->SetRangeUser(0,11.5);
		histo2DSigma->GetXaxis()->SetRangeUser(0,10);		
	} else if (mode == 4 || mode == 5){
        histo2DSigma->GetYaxis()->SetRangeUser(0,17);
        histo2DSigma->GetXaxis()->SetRangeUser(0,10);        
    } 	
	histo2DSigma->DrawCopy(); 

		StylingSliceHistos(Resolution_Pi0_Ptrebin_2,0.8);
		Resolution_Pi0_Ptrebin_2->SetMarkerColor(kBlue+2);
		Resolution_Pi0_Ptrebin_2->SetLineColor(kBlue-8);
		Resolution_Pi0_Ptrebin_2->Draw("same,pe");
// 		DrawResolutionGammaHisto( Resolution_Pi0_Ptrebin_2, 
// 						"", "p^{#pi^{0}}_{T,MC} (GeV/c)", " #sigma p^{#pi^{0}}_{T} (%)", 
// 							kFALSE, 10., 140.,
// 							kTRUE, 0, 10.,
// 							kTRUE, 0., 5.9);
		
		
	padEPGdPtVsPt2->Update();
	canvasEPGdPtVsPt->Update();		
	canvasEPGdPtVsPt->SaveAs(Form("%s/Resolution_%s_Pt.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()));
	delete padEPGdPtVsPt1;	
	delete padEPGdPtVsPt2;	
	delete canvasEPGdPtVsPt;

	ProduceProjectionsPlotWithFits( Resolution_Pi0_PtRes ,fitResolution_Pi0_PtRes, ptbinning, 18 ,5, 4, Form(" (p^{%s}_{T,rec} -p^{%s}_{T,MC})/p^{%s}_{T,MC}", mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()), Form ("N_{%s}", mesonLatex.Data()),"" , Form("fit %s",mesonLatex.Data()), "Pt",Form("%s/ProjectionsFitted%s.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()));


}


