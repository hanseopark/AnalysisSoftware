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

TString textGenerator;
TString collisionSystem;
TString textPeriod;
TString textDate;


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




void ProduceProjectionsPlotWithFits( TH1D** histos, TF1** fitData, Double_t* floatArray, Int_t nBins, Int_t nColumns, Int_t nRows, TString xTitle, Double_t startX, Double_t endX, TString yTitle, TString legendEntry1, TString legendEntry2, TString title, TString nameOutput, Int_t canvasSizeX= 2800, Int_t canvasSizeY = 1800){
										
	TCanvas * canvasProjection = new TCanvas("canvasProjection","",2800,1800);  // gives the page size
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
			TASImage *myAliceLogo = new TASImage("../ALICE_logo.eps");
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

			Double_t yMin = 1.;
			if ( histos[iPt]->GetEntries() > 0){
				if (title.CompareTo("R")==0){
					DrawGammaHistoWithTitle2( histos[iPt],
							Form("%3.2f cm < R < %3.2f cm",startPt,endPt), 
							xTitle.Data(),yTitle.Data(),
							startX,endX,yMin);
				} else if (title.CompareTo("Z")==0){
					DrawGammaHistoWithTitle2( histos[iPt],
							Form("%3.2f cm < Z < %3.2f cm",startPt,endPt), 
							xTitle.Data(),yTitle.Data(),
							startX,endX,yMin);
				} else if (title.CompareTo("Pt")==0){
					DrawGammaHistoWithTitle2( histos[iPt],
							Form("%3.2f GeV/c < p_{T} < %3.2f GeV/c",startPt,endPt), 
							xTitle.Data(),yTitle.Data(),
							startX,endX,yMin);
				} else {
					DrawGammaHistoWithTitle2( histos[iPt],
							Form("%3.2f GeV/c < p < %3.2f GeV/c",startPt,endPt), 
							xTitle.Data(),yTitle.Data(),
							startX,endX,yMin);
				}
				fitData[iPt]->SetLineColor(kBlue);
				fitData[iPt]->SetLineWidth(0.8);
				fitData[iPt]->Draw("same");
			} else {
				place--;
			}
		}
	}
	canvasProjection->Print(nameOutput.Data());
	delete padProjection;
	delete canvasProjection;

}


void PhotonResolutionAdv(const char *fileMonteCarloInput = "", TString optionCutSelection = "", TString suffix = "eps", TString optEnergy="", TString optMCGenerator="", TString optPeriod=""){	
	
	gROOT->Reset();	
	gROOT->SetStyle("Plain");
	
	StyleSettings();	
	SetPlotStyle();

	if(optEnergy.CompareTo("7TeV") == 0){
		collisionSystem = "pp @ #sqrt{#it{s}} = 7 TeV";		
	} else if( optEnergy.CompareTo("2.76TeV") == 0) {
		collisionSystem = "pp @ #sqrt{#it{s}} = 2.76 TeV";
	} else if( optEnergy.CompareTo("900GeV") == 0) {
		collisionSystem = "pp @ #sqrt{#it{s}} = 900 GeV";
	} else if( optEnergy.CompareTo("HI") == 0) {
		collisionSystem ="PbPb @ #sqrt{#it{s}_{_{NN}}} =2.76 TeV";
	} else {
		cout << "No correct collision system specification, has been given" << endl;
		return;		
	}
	
	if(optMCGenerator.CompareTo("") ==0){
		textGenerator = "";
	} else {
		textGenerator = optMCGenerator;
	}

	TDatime 	now;
	int 		iDate = now.GetDate();
	int 		iYear=iDate/10000;
	int 		iMonth=(iDate%10000)/100;
	int 		iDay=iDate%100;
	char* 	cMonth[12]={"Jan","Feb","Mar","Apr","May","Jun",
		    "Jul","Aug","Sep","Oct","Nov","Dec"};
	TString 	textDayth;
	if (iDay%10 == 1){
		textDayth = "st";
	} else if (iDay%10 == 2){
		textDayth = "nd";
	} else if (iDay%10 == 3){
		textDayth = "rd";
	} else {
		textDayth = "th";
	}
	textDate = Form("%i^{%s} %s %i",iDay, textDayth.Data(),cMonth[iMonth-1], iYear);

	
	TString outputDirectory;
	TString outputDirectory1;
	if(optPeriod.CompareTo("") ==0 || optPeriod.CompareTo("All") ==0){
		if (textGenerator.CompareTo("")!=0){
			outputDirectory =	 		Form("%s/%s/%s/%s/PhotonResolution",textGenerator.Data(),optionCutSelection.Data(),optEnergy.Data(),suffix.Data());
			outputDirectory1 = Form("%s/%s",textGenerator.Data(),optionCutSelection.Data());
		} else {
			textPeriod = "";
			outputDirectory =	 		Form("%s/%s/%s/PhotonResolution",optionCutSelection.Data(),optEnergy.Data(),suffix.Data());
			outputDirectory1 = Form("%s",optionCutSelection.Data());
		}
		gSystem->Exec("mkdir -p "+outputDirectory);
	} else if (textGenerator.CompareTo("")!=0){
		textPeriod = optPeriod;
		outputDirectory =	 		Form("%s/%s/%s/%s/%s/PhotonResolution",textGenerator.Data(),optionCutSelection.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
		outputDirectory1 = Form("%s/%s",textGenerator.Data(),optionCutSelection.Data());
		gSystem->Exec("mkdir -p "+outputDirectory);
	} else {
		textPeriod = optPeriod;
		outputDirectory =	 		Form("%s/%s/%s/%s/PhotonResolution",optionCutSelection.Data(),optEnergy.Data(),suffix.Data(),optPeriod.Data());
		outputDirectory1 = Form("%s/%s",textGenerator.Data(),optionCutSelection.Data());
		gSystem->Exec("mkdir -p "+outputDirectory);
	}
	
//**********************************************************************************************************************
//******************************** Definition of some plotting variables ***********************************************
//**********************************************************************************************************************

	//Array defintion for printing Logo 
	Float_t right_up[4]={0.7,0.63,0.15, 0.02};
	Float_t floatLocationRightUp2D[4]=		{0.65,0.84,0.11, 0.02};
	Float_t floatLocationLeftUp2D[4]=		{0.15,0.84,0.11, 0.02};
	Float_t right_down[4]={0.7,0.23,0.15, 0.02};
	Float_t left_up[4]={0.17,0.73, 0.15, 0.02};

	//Defintion of Rebin
	Float_t rebin = 4;	

//**********************************************************************************************************************
//******************************** Defintion of arrays for rebinning ***************************************************
//**********************************************************************************************************************

	//pt rebinning array
	Float_t ptcalcbinning1[5] =  {0.1,0.5 ,1 ,2 ,5};
	Float_t ptcalcbinning2[4] = {14,40,100,200};
	Float_t ptcalcbinning3[6] = {1,11,17,19,21,23};
	Double_t ptbinning[23];
	ptbinning[0] = 0;
	for ( int j = 0; j < 5 ; j ++ ){
		for ( int i = ptcalcbinning3[j]; i<ptcalcbinning3[j+1]; i++){ 
			ptbinning[i] = ptbinning[i-1] +	 ptcalcbinning1[j];
		}
	}
	
	//R rebinning array
	Float_t Rcalcbinning1[18] =  {5, 0.5, 2, 0.5, 4,0.5,1,0.5,3.5,0.5,4,0.5,5,0.5,4,1,11,15};
	Float_t Rcalcbinning3[19] = {1,2,16,17,27,28,38,39,41,43,65,68,71,74,79,80,81,82,87};
	Double_t Rbinning[87];
	Rbinning[0] = 0;
	for ( int j = 0; j < 18 ; j ++ ){
		for ( int i = Rcalcbinning3[j]; i<Rcalcbinning3[j+1]; i++){ 
			Rbinning[i] = Rbinning[i-1] +	 Rcalcbinning1[j];
		}
	}

	//Z rebinning array
	Float_t Zcalcbinning1[9] =  {50, 10, 5, 4, 2,4,5,10,50};
	Float_t Zcalcbinning3[10] = {1,3,5,11,21,31,41,47,49,52};
	Double_t Zbinning[51];
	Zbinning[0] = -200;
	for ( int j = 0; j < 9 ; j ++ ){
		for ( int i = Zcalcbinning3[j]; i<Zcalcbinning3[j+1]; i++){ 
			Zbinning[i] = Zbinning[i-1] +	 Zcalcbinning1[j];
		}
	}


//**********************************************************************************************************************
//**************** Definition of maximum and minimum values for histogramms depending on the eta cutnumber *************
//**********************************************************************************************************************

	Double_t maximumfracPtITSCluster = 32.5;
	Double_t maxdZAbsSigmaR = 10.5;
	Double_t maxdZAbsSigmaZ = 2.5;
	Double_t arraydZAbsMeanR[2] = {-3.1,3.1};
	Double_t arraydZAbsMeanZ[2] = {-2.5,2.5};
	Double_t maxdRAbsSigmaR = 8.25;
	Double_t arraydRAbsMeanR[2] = {-4.1,4.1};
	Double_t maxdPtSigmaR = 10.25;
	Double_t arraydPtMeanR[2] = {-8.2,6.2};
	Double_t maxdPhiAbsSigmaR = 10.5e-3;
	Double_t arraydPhiAbsMeanR[2] = {-0.9e-3,0.9e-3};
	Double_t maxPtSigmaPt = 23.5;
	Double_t arrayPtMeanPt[2] = {-4.2,4.2};
	Double_t maxdPtSigmaPhi = 10.5;
	Double_t arraydPtMeanPhi[2] = {-2.2,3.2};


	TString etaCutNumber = optionCutSelection(1,1);
	TLatex *latexEtaRange;
	if(etaCutNumber.CompareTo("2") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"0.9 < |#eta| < 1.4 "); // Bo: this was modified
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);
			maxdZAbsSigmaZ = 10.;
			arraydZAbsMeanZ[0] = -10. ;
			arraydZAbsMeanZ[1] = 10.;
	} else if(etaCutNumber.CompareTo("0") == 0){
			latexEtaRange = 	new TLatex(0.15,0.92,"|#eta| < 0.9 "); // Bo: this was modified
			latexEtaRange->SetNDC();
			latexEtaRange->SetTextFont(62);
			latexEtaRange->SetTextSize(0.04);
			latexEtaRange->SetLineWidth(6);      
			maxdZAbsSigmaR = 3.1;
			maxdRAbsSigmaR = 2.35;
			arraydRAbsMeanR[0] = -0.95;
			arraydRAbsMeanR[1] = 0.95;
			arraydPhiAbsMeanR[0] = -0.75e-3;
			arraydPhiAbsMeanR[1] = 0.75e-3;
			maxdPhiAbsSigmaR = 7.5e-3;
			arraydZAbsMeanR[0] = -0.26;
			arraydZAbsMeanR[1] = 0.26;
			maxdPtSigmaR = 2.75;
			arraydPtMeanR[0] = -1;
			arraydPtMeanR[1] = 1;
			maxdPtSigmaPhi = 2.1;
			arraydPtMeanPhi[0] = -0.25;
			arraydPtMeanPhi[1] = 1.1;
	} 
	
//**********************************************************************************************************************
//******************************************* Defintion of file name ***************************************************
//**********************************************************************************************************************
	TString nameGammaDirectory;
	TString nameGammaList;
	if(optionCutSelection.CompareTo("") != 0){
		nameGammaDirectory = Form("GammaConv_%s",  optionCutSelection.Data());
		cout << nameGammaDirectory.Data() << endl;
	}	

//**********************************************************************************************************************
//****************************************** Loading of Histograms *****************************************************
//**********************************************************************************************************************

	TFile fileMC(fileMonteCarloInput);  
	
	//************************** Container Loading ********************************************************************
	TDirectory *directoryResolution = new TDirectory(); 
	
	if(!(directoryResolution = (TDirectory*)fileMC.Get(nameGammaDirectory.Data()))) cout <<"PWG4GammConversion TList NOT loaded correctly"<<endl; 
	
	//******************************** Resolution in Absolute ********************************************
	TH2F * Resolution_dZAbs_VS_R = (TH2F*)directoryResolution->Get("Resolution_dZAbs_VS_R" );
	TH2F * Resolution_dZAbs_VS_Z = (TH2F*)directoryResolution->Get("Resolution_dZAbs_VS_Z" );
	TH2F * Resolution_dZAbs_VS_Phi = (TH2F*)directoryResolution->Get("Resolution_dZAbs_VS_Phi" );
	TH2F * Resolution_dZAbs_VS_Eta = (TH2F*)directoryResolution->Get("Resolution_dZAbs_VS_Eta" );
	TH2F * Resolution_dZAbs_VS_Pt = (TH2F*)directoryResolution->Get("Resolution_dZAbs_VS_Pt" );
	
	TH2F * Resolution_dPhiAbs_VS_R = (TH2F*)directoryResolution->Get("Resolution_dPhiAbs_VS_R" );
	TH2F * Resolution_dPhiAbs_VS_Phi = (TH2F*)directoryResolution->Get("Resolution_dPhiAbs_VS_Phi" );
	TH2F * Resolution_dPhiAbs_VS_Z = (TH2F*)directoryResolution->Get("Resolution_dPhiAbs_VS_Z" );
	TH2F * Resolution_dPhiAbs_VS_Eta = (TH2F*)directoryResolution->Get("Resolution_dPhiAbs_VS_Eta" );
	TH2F * Resolution_dPhiAbs_VS_Pt = (TH2F*)directoryResolution->Get("Resolution_dPhiAbs_VS_Pt" );

	TH2F * Resolution_dRAbs_VS_R = (TH2F*)directoryResolution->Get("Resolution_dRAbs_VS_R" );
	TH2F * Resolution_dRAbs_VS_Phi = (TH2F*)directoryResolution->Get("Resolution_dRAbs_VS_Phi" );
	TH2F * Resolution_dRAbs_VS_Z = (TH2F*)directoryResolution->Get("Resolution_dRAbs_VS_Z" );
	TH2F * Resolution_dRAbs_VS_Eta = (TH2F*)directoryResolution->Get("Resolution_dRAbs_VS_Eta" );
	TH2F * Resolution_dRAbs_VS_Pt = (TH2F*)directoryResolution->Get("Resolution_dRAbs_VS_Pt" );

	TH2F * Resolution_dEtaAbs_VS_R = (TH2F*)directoryResolution->Get("Resolution_dEtaAbs_VS_R" );
	TH2F * Resolution_dEtaAbs_VS_Phi = (TH2F*)directoryResolution->Get("Resolution_dEtaAbs_VS_Phi" );
	TH2F * Resolution_dEtaAbs_VS_Z = (TH2F*)directoryResolution->Get("Resolution_dEtaAbs_VS_Z" );
	TH2F * Resolution_dEtaAbs_VS_Eta = (TH2F*)directoryResolution->Get("Resolution_dEtaAbs_VS_Eta" );
	TH2F * Resolution_dEtaAbs_VS_Pt = (TH2F*)directoryResolution->Get("Resolution_dEtaAbs_VS_Pt" );

	TH2F * Resolution_dPt_VS_R = (TH2F*)directoryResolution->Get("Resolution_dPt_VS_R" );
	TH2F * Resolution_dPt_VS_Phi = (TH2F*)directoryResolution->Get("Resolution_dPt_VS_Phi" );
	TH2F * Resolution_dPt_VS_Z = (TH2F*)directoryResolution->Get("Resolution_dPt_VS_Z" );
	TH2F * Resolution_dPt_VS_Eta = (TH2F*)directoryResolution->Get("Resolution_dPt_VS_Eta" );
	TH2F * Resolution_dPt_VS_Pt = (TH2F*)directoryResolution->Get("Resolution_dPt_VS_Pt" );

//**********************************************************************************************************************************************
//*********************************** Rebinning of 2D histograms and fitslices *****************************************************************
//**********************************************************************************************************************************************
	Int_t maxYbin = 0;
	Int_t maxXbin = 0;

		//************************************** Rebinning of dZAbs R **********************************************
		TH2F *Resolution_dZAbs_VS_Rrebin = new TH2F("Resolution_dZAbs_VS_Rrebin", "Photon Resolution dZabs vs R ", 86, Rbinning , 50, -25., 25.) ;	
		Resolution_dZAbs_VS_R->Sumw2();
		maxXbin = Resolution_dZAbs_VS_R->GetNbinsX();		
		maxYbin = Resolution_dZAbs_VS_R->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dZAbs_VS_R->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dZAbs_VS_R->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dZAbs_VS_R->GetBinContent(i, j);
					Resolution_dZAbs_VS_Rrebin->Fill(binoldpoint, biny, value);
				}
		}

	//************************************** Rebinning of dZAbs Z **********************************************
		TH2F *Resolution_dZAbs_VS_Zrebin = new TH2F("Resolution_dZAbs_VS_Zrebin", "Photon Resolution dZabs vs Z ", 49, Zbinning , 50, -25., 25.) ;	
		Resolution_dZAbs_VS_Z->Sumw2();
		maxXbin = Resolution_dZAbs_VS_Z->GetNbinsX();		
		maxYbin = Resolution_dZAbs_VS_Z->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dZAbs_VS_Z->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dZAbs_VS_Z->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dZAbs_VS_Z->GetBinContent(i, j);
					Resolution_dZAbs_VS_Zrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dZAbs vs Pt ********************************************
		TH2F *Resolution_dZAbs_VS_Ptrebin = new TH2F("Resolution_dZAbs_VS_Ptrebin", "Photon Resolution dRabs vs Pt ", 21, ptbinning , 50, -25., 25.) ;	
		Resolution_dZAbs_VS_Pt->Sumw2();
		maxXbin = Resolution_dZAbs_VS_Pt->GetNbinsX();		
		maxYbin = Resolution_dZAbs_VS_Pt->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dZAbs_VS_Pt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dZAbs_VS_Pt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dZAbs_VS_Pt->GetBinContent(i, j);
					Resolution_dZAbs_VS_Ptrebin->Fill(binoldpoint, biny, value);
				}
		}

		
		//************************************** Rebinning of dRAbs Z **********************************************
		TH2F *Resolution_dRAbs_VS_Zrebin = new TH2F("Resolution_dRAbs_VS_Zrebin", "Photon Resolution dRabs vs Z ", 49, Zbinning , 50, -25., 25.) ;	
		Resolution_dRAbs_VS_Z->Sumw2();
		maxXbin = Resolution_dRAbs_VS_Z->GetNbinsX();		
		maxYbin = Resolution_dRAbs_VS_Z->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dRAbs_VS_Z->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dRAbs_VS_Z->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dRAbs_VS_Z->GetBinContent(i, j);
					Resolution_dRAbs_VS_Zrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dRAbs vs Pt ********************************************
		TH2F *Resolution_dRAbs_VS_Ptrebin = new TH2F("Resolution_dRAbs_VS_Ptrebin", "Photon Resolution dRabs vs Pt ", 21, ptbinning , 50, -25., 25.) ;	
		Resolution_dRAbs_VS_Pt->Sumw2();
		maxXbin = Resolution_dRAbs_VS_Pt->GetNbinsX();		
		maxYbin = Resolution_dRAbs_VS_Pt->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dRAbs_VS_Pt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dRAbs_VS_Pt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dRAbs_VS_Pt->GetBinContent(i, j);
					Resolution_dRAbs_VS_Ptrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dRAbs R **********************************************
		TH2F *Resolution_dRAbs_VS_Rrebin = new TH2F("Resolution_dRAbs_VS_Rrebin", "Photon Resolution dRabs vs R ", 86, Rbinning , 50, -25,25) ;	
		Resolution_dRAbs_VS_R->Sumw2();
		maxXbin = Resolution_dRAbs_VS_R->GetNbinsX();		
		maxYbin = Resolution_dRAbs_VS_R->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dRAbs_VS_R->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dRAbs_VS_R->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dRAbs_VS_R->GetBinContent(i, j);
					Resolution_dRAbs_VS_Rrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dPhiAbs Z **********************************************
		TH2F *Resolution_dPhiAbs_VS_Zrebin = new TH2F("Resolution_dPhiAbs_VS_Zrebin", "Photon Resolution dRabs vs Z ", 49, Zbinning , 100, -TMath::Pi()/30, TMath::Pi()/30);
		Resolution_dPhiAbs_VS_Z->Sumw2();
		maxXbin = Resolution_dPhiAbs_VS_Z->GetNbinsX();		
		maxYbin = Resolution_dPhiAbs_VS_Z->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dPhiAbs_VS_Z->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dPhiAbs_VS_Z->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dPhiAbs_VS_Z->GetBinContent(i, j);
					Resolution_dPhiAbs_VS_Zrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dPhiAbs R **********************************************
		TH2F *Resolution_dPhiAbs_VS_Rrebin = new TH2F("Resolution_dPhiAbs_VS_Rrebin", "Photon Resolution dPhiabs vs R ", 86, Rbinning , 100, -TMath::Pi()/30, TMath::Pi()/30) ;	
		Resolution_dPhiAbs_VS_R->Sumw2();
		maxXbin = Resolution_dPhiAbs_VS_R->GetNbinsX();		
		maxYbin = Resolution_dPhiAbs_VS_R->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dPhiAbs_VS_R->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dPhiAbs_VS_R->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dPhiAbs_VS_R->GetBinContent(i, j);
					Resolution_dPhiAbs_VS_Rrebin->Fill(binoldpoint, biny, value);
				}
		}

		//************************************** Rebinning of dRAbs vs Pt ********************************************
		TH2F *Resolution_dPhiAbs_VS_Ptrebin = new TH2F("Resolution_dPhiAbs_VS_Ptrebin", "Photon Resolution dPhiabs vs Pt ", 21, ptbinning , 100, -TMath::Pi()/30, TMath::Pi()/30);	
		Resolution_dPhiAbs_VS_Pt->Sumw2();
		maxXbin = Resolution_dPhiAbs_VS_Pt->GetNbinsX();		
		maxYbin = Resolution_dPhiAbs_VS_Pt->GetNbinsY();
		for (Int_t i = 1; i <(maxXbin+1); i++) {
				TAxis *xaxis_old = 	Resolution_dPhiAbs_VS_Pt->GetXaxis()	;		
				Double_t binoldpoint = xaxis_old->GetBinCenter(i);							
				for (Int_t j = 1 ; j < (maxYbin+1) ; j++){ 		
					Double_t biny = Resolution_dPhiAbs_VS_Pt->GetYaxis()->GetBinCenter(j);								
					Double_t value =  Resolution_dPhiAbs_VS_Pt->GetBinContent(i, j);
					Resolution_dPhiAbs_VS_Ptrebin->Fill(binoldpoint, biny, value);
				}
		}

		
		
		Double_t precision = 10E-5;

//*********************************** Fitting for dRabs vs R **********************************************
		TH1D* Resolution_R_dRRes[86];
		TF1*  fitResolution_R_dRRes[86];
		TH1F *Resolution_dRAbs_VS_Rrebin_1 = new TH1F("Resolution_dRAbs_VS_Rrebin_1", "Mean Photon Resolution dRabs vs R ", 86, Rbinning ) ;	
		TH1F *Resolution_dRAbs_VS_Rrebin_2 = new TH1F("Resolution_dRAbs_VS_Rrebin_2", "Sigma Photon Resolution dRabs vs R ", 86, Rbinning ) ;
		ResolutionFitting( Resolution_dRAbs_VS_Rrebin , Resolution_R_dRRes , fitResolution_R_dRRes ,86, Resolution_dRAbs_VS_Rrebin_1 ,Resolution_dRAbs_VS_Rrebin_2, "gaus",-10.,10.,precision, "fitResolution_R_dRRes");

//********************************************* fitting for dZAbs vs R *************************************************************
		TH1D* Resolution_R_dZRes[86];
		TF1*  fitResolution_R_dZRes[86];
		TH1F *Resolution_dZAbs_VS_Rrebin_1 = new TH1F("Resolution_dZAbs_VS_Rrebin_1", "Mean Photon Resolution dZabs vs R ", 86, Rbinning ) ;	
		TH1F *Resolution_dZAbs_VS_Rrebin_2 = new TH1F("Resolution_dZAbs_VS_Rrebin_2", "Sigma Photon Resolution dZabs vs R ", 86, Rbinning ) ;
		ResolutionFitting( Resolution_dZAbs_VS_Rrebin , Resolution_R_dZRes , fitResolution_R_dZRes ,86, Resolution_dZAbs_VS_Rrebin_1 ,Resolution_dZAbs_VS_Rrebin_2, "gaus",-10.,10.,precision, "fitResolution_R_dZRes");

//********************************************* fitting for dZAbs vs Z *************************************************************
		TH1D* Resolution_Z_dZRes[50];
		TF1*  fitResolution_Z_dZRes[50];
		TH1F *Resolution_dZAbs_VS_Zrebin_1 = new TH1F("Resolution_dZAbs_VS_Zrebin_1", "Mean Photon Resolution dZabs vs Z ", 50, Zbinning ) ;	
		TH1F *Resolution_dZAbs_VS_Zrebin_2 = new TH1F("Resolution_dZAbs_VS_Zrebin_2", "Sigma Photon Resolution dZabs vs Z ", 50, Zbinning ) ;
		ResolutionFitting( Resolution_dZAbs_VS_Zrebin , Resolution_Z_dZRes , fitResolution_Z_dZRes ,50, Resolution_dZAbs_VS_Zrebin_1 ,Resolution_dZAbs_VS_Zrebin_2, "gaus",-10.,10.,precision, "Resolution_Z_dZRes");

		TH1D* Resolution_Z_dRRes[50];
		TF1*  fitResolution_Z_dRRes[50];
		TH1F *Resolution_dRAbs_VS_Zrebin_1 = new TH1F("Resolution_dRAbs_VS_Zrebin_1", "Mean Photon Resolution dZabs vs Z ", 50, Zbinning ) ;	
		TH1F *Resolution_dRAbs_VS_Zrebin_2 = new TH1F("Resolution_dRAbs_VS_Zrebin_2", "Sigma Photon Resolution dZabs vs Z ", 50, Zbinning ) ;
		ResolutionFitting( Resolution_dRAbs_VS_Zrebin , Resolution_Z_dRRes , fitResolution_Z_dRRes ,50, Resolution_dRAbs_VS_Zrebin_1 ,Resolution_dRAbs_VS_Zrebin_2, "gaus",-10.,10.,precision, "Resolution_Z_dRRes");

		TH1D* Resolution_Pt_dRRes[22];
		TF1*  fitResolution_Pt_dRRes[22];
		TH1F *Resolution_dRAbs_VS_Ptrebin_1 = new TH1F("Resolution_dRAbs_VS_Ptrebin_1", "Mean Photon Resolution dRabs vs Pt ", 21, ptbinning ) ;	
		TH1F *Resolution_dRAbs_VS_Ptrebin_2 = new TH1F("Resolution_dRAbs_VS_Ptrebin_2", "Sigma Photon Resolution dRabs vs Pt ", 21, ptbinning ) ;
		ResolutionFitting( Resolution_dRAbs_VS_Ptrebin , Resolution_Pt_dRRes , fitResolution_Pt_dRRes ,21, Resolution_dRAbs_VS_Ptrebin_1 ,Resolution_dRAbs_VS_Ptrebin_2, "gaus",-10.,10.,precision, "Resolution_Pt_dRRes");
		
		TH1D* Resolution_Pt_dZRes[22];
		TF1*  fitResolution_Pt_dZRes[22];
		TH1F *Resolution_dZAbs_VS_Ptrebin_1 = new TH1F("Resolution_dZAbs_VS_Ptrebin_1", "Mean Photon Resolution dRabs vs Pt ", 21, ptbinning ) ;	
		TH1F *Resolution_dZAbs_VS_Ptrebin_2 = new TH1F("Resolution_dZAbs_VS_Ptrebin_2", "Sigma Photon Resolution dRabs vs Pt ", 21, ptbinning ) ;
		ResolutionFitting( Resolution_dZAbs_VS_Ptrebin , Resolution_Pt_dZRes , fitResolution_Pt_dZRes ,21, Resolution_dZAbs_VS_Ptrebin_1 ,Resolution_dZAbs_VS_Ptrebin_2, "gaus",-10.,10.,precision, "Resolution_Pt_dZRes");
		
		Double_t precisionPhi = 10E-10;
		TH1D* Resolution_Z_dPhiRes[50];
		TF1*  fitResolution_Z_dPhiRes[50];
		TH1F *Resolution_dPhiAbs_VS_Zrebin_1 = new TH1F("Resolution_dPhiAbs_VS_Zrebin_1", "Mean Photon Resolution dZabs vs Z ", 50, Zbinning);	
		TH1F *Resolution_dPhiAbs_VS_Zrebin_2 = new TH1F("Resolution_dPhiAbs_VS_Zrebin_2", "Sigma Photon Resolution dZabs vs Z ", 50, Zbinning);
		ResolutionFitting( Resolution_dPhiAbs_VS_Zrebin , Resolution_Z_dPhiRes , fitResolution_Z_dPhiRes ,50, Resolution_dPhiAbs_VS_Zrebin_1 ,Resolution_dPhiAbs_VS_Zrebin_2, "gaus",-10.,10.,precisionPhi, "Resolution_Z_dPhiRes");
		
		
//********************************************* fitting for dPhiAbs vs R *************************************************************
		TH1D* Resolution_R_dPhiAbsRes[86];
		TF1* fitResolution_R_dPhiAbsRes[86];
		TH1F *Resolution_dPhiAbs_VS_Rrebin_1 = new TH1F("Resolution_dPhiAbs_VS_Rrebin_1", "Mean Photon Resolution dPhiabs vs R ", 86, Rbinning);	
		TH1F *Resolution_dPhiAbs_VS_Rrebin_2 = new TH1F("Resolution_dPhiAbs_VS_Rrebin_2", "Sigma Photon Resolution dPhiabs vs R ", 86, Rbinning);
		ResolutionFitting( Resolution_dPhiAbs_VS_Rrebin , Resolution_R_dPhiAbsRes , fitResolution_R_dPhiAbsRes ,86, Resolution_dPhiAbs_VS_Rrebin_1 ,Resolution_dPhiAbs_VS_Rrebin_2, "gaus",-10.,10.,precisionPhi , "Resolution_R_dPhiAbsRes");

		TH1D* Resolution_Phi_dPhiRes[100];
		TF1*  fitResolution_Phi_dPhiRes[100];
		TH1F *Resolution_dPhiAbs_VS_Phi_1 = new TH1F("Resolution_dPhiAbs_VS_Phirebin_1", "Mean Photon Resolution dPhiabs vs Phi ", 100, 0, 2*TMath::Pi());	
		TH1F *Resolution_dPhiAbs_VS_Phi_2 = new TH1F("Resolution_dPhiAbs_VS_Phirebin_2", "Sigma Photon Resolution dPhiabs vs Phi ", 100,0, 2*TMath::Pi());
		ResolutionFitting( Resolution_dPhiAbs_VS_Phi , Resolution_Phi_dPhiRes , fitResolution_Phi_dPhiRes ,100, Resolution_dPhiAbs_VS_Phi_1 ,Resolution_dPhiAbs_VS_Phi_2, "gaus",-10.,10.,precisionPhi, "Resolution_Phi_dPhiRes");

		TH1D* Resolution_Pt_dPhiRes[22];
		TF1*  fitResolution_Pt_dPhiRes[22];
		TH1F *Resolution_dPhiAbs_VS_Ptrebin_1 = new TH1F("Resolution_dPhiAbs_VS_Ptrebin_1", "Mean Photon Resolution dPhiabs vs Pt ", 21, ptbinning ) ;	
		TH1F *Resolution_dPhiAbs_VS_Ptrebin_2 = new TH1F("Resolution_dPhiAbs_VS_Ptrebin_2", "Sigma Photon Resolution dPhiabs vs Pt ", 21, ptbinning ) ;
		ResolutionFitting( Resolution_dPhiAbs_VS_Ptrebin , Resolution_Pt_dPhiRes , fitResolution_Pt_dPhiRes ,21, Resolution_dPhiAbs_VS_Ptrebin_1 ,Resolution_dPhiAbs_VS_Ptrebin_2, "gaus",-10.,10.,precisionPhi, "Resolution_Pt_dPhiRes");
		
		TH1D* Resolution_Phi_dZRes[100];
		TF1*  fitResolution_Phi_dZRes[100];
		TH1F *Resolution_dZAbs_VS_Phi_1 = new TH1F("Resolution_dZAbs_VS_Phirebin_1", "Mean Photon Resolution dZabs vs Phi ", 100, 0, 2*TMath::Pi());	
		TH1F *Resolution_dZAbs_VS_Phi_2 = new TH1F("Resolution_dZAbs_VS_Phirebin_2", "Sigma Photon Resolution dZabs vs Phi ", 100, 0, 2*TMath::Pi());	
		ResolutionFitting( Resolution_dZAbs_VS_Phi , Resolution_Phi_dZRes , fitResolution_Phi_dZRes ,100, Resolution_dZAbs_VS_Phi_1 ,Resolution_dZAbs_VS_Phi_2, "gaus",-10.,10.,precisionPhi, "Resolution_Phi_dZRes");

		TH1D* Resolution_Phi_dRRes[100];
		TF1*  fitResolution_Phi_dRRes[100];
		TH1F *Resolution_dRAbs_VS_Phi_1 = new TH1F("Resolution_dRAbs_VS_Phirebin_1", "Mean Photon Resolution dRabs vs Phi ", 100, 0, 2*TMath::Pi());	
		TH1F *Resolution_dRAbs_VS_Phi_2 = new TH1F("Resolution_dRAbs_VS_Phirebin_2", "Sigma Photon Resolution dRabs vs Phi ", 100,0, 2*TMath::Pi());	
		ResolutionFitting( Resolution_dRAbs_VS_Phi , Resolution_Phi_dRRes , fitResolution_Phi_dRRes ,100, Resolution_dRAbs_VS_Phi_1 ,Resolution_dRAbs_VS_Phi_2, "gaus",-10.,10.,precisionPhi, "Resolution_Phi_dRRes");

// ----------------------------- Resolution dPhiAbs vs R ---------------------------------------
	TCanvas * canvasPhiAbsVsR = new TCanvas("canvasPhiAbsVsR","",10,10,700,1000);  // gives the page size		
	canvasPhiAbsVsR->cd();
	DrawGammaCanvasSettings( canvasPhiAbsVsR, 0.13, 0.02, 0.02, 0.09);
	
	TPad* padPhiAbsVsR1 = new TPad("padPhiAbsVsR1", "", 0., 0.5, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padPhiAbsVsR1, 0.12, 0.02, 0.02, 0.);
	padPhiAbsVsR1->Draw();
	
	TPad* padPhiAbsVsR2 = new TPad("padPhiAbsVsR2", "", 0., 0., 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( padPhiAbsVsR2, 0.12, 0.02, 0., 0.1);
	padPhiAbsVsR2->Draw();
	
	padPhiAbsVsR1->cd();
	padPhiAbsVsR1->SetTopMargin(0.04);
		StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_1,0.8);
		Resolution_dPhiAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
		Resolution_dPhiAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_1, 
							"", "R (cm)", "Peak position d#phi (rad)", 
							kFALSE, 10., 140.,
							kTRUE, arraydPhiAbsMeanR[0] , arraydPhiAbsMeanR[1],
							kTRUE, 0., 120.);
		DrawGammaLines(0., 120, 0,0,0.005);
	padPhiAbsVsR1->Update();
	padPhiAbsVsR2->cd();
		StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_2,0.8);
		Resolution_dPhiAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
		Resolution_dPhiAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_2, 
						"", "R (cm)", "#sigma d#phi (rad)", 
							kFALSE, 10., 140.,
							kTRUE, 0.0, maxdPhiAbsSigmaR,
							kTRUE, 0., 120.);
	padPhiAbsVsR2->Update();
	canvasPhiAbsVsR->Update();	
	canvasPhiAbsVsR->SaveAs(Form("%s/Resolution_dphiAbs_R.%s",outputDirectory.Data(),suffix.Data()));
	delete padPhiAbsVsR1;	
	delete padPhiAbsVsR2;	
	delete canvasPhiAbsVsR;

// ---------------------------------- Resolution dRAbs vs R -------------------------------
	TCanvas * canvasdAbsRVsR = new TCanvas("canvasdAbsRVsR","",10,10,700,1000);  // gives the page size		
	canvasdAbsRVsR->cd();
	DrawGammaCanvasSettings( canvasdAbsRVsR, 0.13, 0.02, 0.02, 0.09);
	
	TPad* paddAbsRVsR1 = new TPad("paddAbsRVsR1", "", 0., 0.5, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( paddAbsRVsR1, 0.12, 0.02, 0.02, 0.);
	paddAbsRVsR1->Draw();
	
	TPad* paddAbsRVsR2 = new TPad("paddAbsRVsR2", "", 0., 0., 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( paddAbsRVsR2, 0.12, 0.02, 0., 0.1);
	paddAbsRVsR2->Draw();

	paddAbsRVsR1->cd();
		StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_1,0.8);
		Resolution_dRAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
		Resolution_dRAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_1, 
							"", "R (cm)", "Peak pos. dR (cm)", 
							kFALSE, 10., 140.,
							kTRUE, arraydRAbsMeanR[0] , arraydRAbsMeanR[1],
							kTRUE, 0., 120.);
		DrawGammaLines(0., 120, 0,0,0.005);
	paddAbsRVsR1->Update();
	paddAbsRVsR2->cd();
		StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_2,0.8);
		Resolution_dRAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
		Resolution_dRAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_2, 
						"", "R (cm)", "#sigma dR (cm)", 
							kFALSE, 10., 140.,
							kTRUE, 0., maxdRAbsSigmaR,
							kTRUE, 0.8, 120.);
	paddAbsRVsR2->Update();
	canvasdAbsRVsR->Update();	
	canvasdAbsRVsR->SaveAs(Form("%s/Resolution_dRAbs_R.%s",outputDirectory.Data(),suffix.Data()));
	delete paddAbsRVsR1;	
	delete paddAbsRVsR2;	
	delete canvasdAbsRVsR;
// -------------------------------- Resolution dZAbs vs R -----------------------------	
	TCanvas * canvasdAbsZVsR = new TCanvas("canvasdAbsZVsR","",10,10,700,1000);  // gives the page size		
	canvasdAbsZVsR->cd();
	DrawGammaCanvasSettings( canvasdAbsZVsR, 0.13, 0.02, 0.02, 0.09);
	
	TPad* paddAbsZVsR1 = new TPad("paddAbsZVsR1", "", 0., 0.5, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( paddAbsZVsR1, 0.12, 0.02, 0.02, 0.);
	paddAbsZVsR1->Draw();
	
	TPad* paddAbsZVsR2 = new TPad("paddAbsZVsR2", "", 0., 0., 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( paddAbsZVsR2, 0.12, 0.02, 0., 0.1);
	paddAbsZVsR2->Draw();

	paddAbsZVsR1->cd();
		StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_1,0.8);
		Resolution_dZAbs_VS_Rrebin_1->SetMarkerColor(kRed+2);
		Resolution_dZAbs_VS_Rrebin_1->SetLineColor(kBlue-8)	;	
		DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_1, 
							"", "R (cm)", "Peak pos. dZ (cm)", 
							kFALSE, 10., 140.,
							kTRUE, arraydZAbsMeanR[0], arraydZAbsMeanR[1],
							kTRUE, 0., 120.);
		DrawGammaLines(0., 120, 0,0,0.005);
	paddAbsZVsR1->Update();
	paddAbsZVsR2->cd();
		StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_2,0.8);
		Resolution_dZAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
		Resolution_dZAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_2, 
						"", "R (cm)", "#sigma dZ (cm)", 
							kFALSE, 10., 140.,
							kTRUE, 0., maxdZAbsSigmaR,
							kTRUE, 0., 120.);
	paddAbsZVsR2->Update();
	canvasdAbsZVsR->Update();	
	canvasdAbsZVsR->SaveAs(Form("%s/Resolution_dZAbs_R.%s",outputDirectory.Data(),suffix.Data()));
	delete paddAbsZVsR1;	
	delete paddAbsZVsR2;	
	delete canvasdAbsZVsR;

// -------------------------------- Resolution dZAbs vs Z -----------------------------	
	TCanvas * canvasdAbsZVsZ = new TCanvas("canvasdAbsZVsZ","",10,10,700,1000);  // gives the page size		
	canvasdAbsZVsZ->cd();
	DrawGammaCanvasSettings( canvasdAbsZVsZ, 0.13, 0.02, 0.02, 0.09);
	
	TPad* paddAbsZVsZ1 = new TPad("paddAbsZVsZ1", "", 0., 0.5, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( paddAbsZVsZ1, 0.12, 0.02, 0.02, 0.);
	paddAbsZVsZ1->Draw();
	
	TPad* paddAbsZVsZ2 = new TPad("paddAbsZVsZ2", "", 0., 0., 1., 0.5,-1, -1, -2);
	DrawGammaPadSettings( paddAbsZVsZ2, 0.12, 0.02, 0., 0.1);
	paddAbsZVsZ2->Draw();

	paddAbsZVsZ1->cd();
		StylingSliceHistos(Resolution_dZAbs_VS_Zrebin_1,0.8);
		Resolution_dZAbs_VS_Zrebin_1->SetMarkerColor(kRed+2);
		Resolution_dZAbs_VS_Zrebin_1->SetLineColor(kBlue-8)	;	
		DrawResolutionGammaHisto( Resolution_dZAbs_VS_Zrebin_1, 
							"", "Z (cm)", "Peak pos. dZ (cm)", 
							kFALSE, 10., 140.,
							kTRUE, arraydZAbsMeanZ[0], arraydZAbsMeanZ[1],
							kTRUE, -200., 200.);
		DrawGammaLines(-200., 200, 0,0,0.005);
	paddAbsZVsZ1->Update();
	paddAbsZVsZ2->cd();
		StylingSliceHistos(Resolution_dZAbs_VS_Zrebin_2,0.8);
		Resolution_dZAbs_VS_Zrebin_2->SetMarkerColor(kRed+2);
		Resolution_dZAbs_VS_Zrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dZAbs_VS_Zrebin_2, 
						"", "Z (cm)", "#sigma dZ (cm)", 
							kFALSE, 10., 140.,
							kTRUE, 0., maxdZAbsSigmaZ,
							kTRUE, -200., 200.);
	paddAbsZVsZ2->Update();
	canvasdAbsZVsZ->Update();	
	canvasdAbsZVsZ->SaveAs(Form("%s/Resolution_dZAbs_Z.%s",outputDirectory.Data(),suffix.Data()));
	delete paddAbsZVsZ1;	
	delete paddAbsZVsZ2;	
	delete canvasdAbsZVsZ;

	TCanvas * canvasdXvsPt = new TCanvas("canvasdXvsPt","",10,10,700,2000);  // gives the page size		
	canvasdXvsPt->cd();
	
	TPad* padResFullvsPt1 = new TPad("padResFullvsPt1", "",  0., 0.76, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padResFullvsPt1, 0.1, 0.02, 0.02, 0.);
	padResFullvsPt1->Draw();
	TPad* padResFullvsPt2 = new TPad("padResFullvsPt2", "",0., 0.515, 1., 0.755,-1, -1, -2);
	DrawGammaPadSettings( padResFullvsPt2, 0.1, 0.02, 0., 0.0);
	padResFullvsPt2->Draw();
	TPad* padResFullvsPt3 = new TPad("padResFullvsPt3", "",0., 0.275, 1., 0.515,-1, -1, -2);
	DrawGammaPadSettings( padResFullvsPt3, 0.1, 0.02, 0., 0.0);
	padResFullvsPt3->Draw();
	TPad* padResFullvsPt4 = new TPad("padResFullvsPt4", "",0., 0., 1., 0.275,-1, -1, -2);
	DrawGammaPadSettings( padResFullvsPt4, 0.1, 0.02, 0., 0.1);
	padResFullvsPt4->Draw();
	
	padResFullvsPt3->cd();
		Resolution_dPhiAbs_VS_Ptrebin_2->Scale(1000.);
		DrawGammaSetMarker(Resolution_dPhiAbs_VS_Ptrebin_2, 21,1., kRed+2 , kRed-9);
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Ptrebin_2, 
							"", "", "#sigma (mrad)", 
							kFALSE, 10., 140.,
							kTRUE, -1 , maxdPhiAbsSigmaR*1000,
							kTRUE, 0., 10.);
// 		DrawGammaLines(0., 120, 0,0,0.005);
	
	padResFullvsPt4->cd();
		DrawGammaSetMarker(Resolution_dRAbs_VS_Ptrebin_2, 20,1., kBlue+2 , kBlue-9);
		DrawGammaSetMarker(Resolution_dZAbs_VS_Ptrebin_2, 22,1., kGreen+2 , kGreen-9);
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Ptrebin_2, 
						"", "p_{t} (GeV/c)", "#sigma (cm)", 
							kFALSE, 10., 140.,
							kTRUE, 0., maxdRAbsSigmaR,
							kTRUE, 0., 10.);
		Resolution_dZAbs_VS_Ptrebin_2->Draw("same,e1");
		Resolution_dRAbs_VS_Ptrebin_2->Draw("same,e1");
		
	canvasdXvsPt->Update();	
	canvasdXvsPt->SaveAs(Form("%s/Resolution_dZAnddRvsPt.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasdXvsPt;

	

	TCanvas * canvasZdist = new TCanvas("canvasZdist","",10,10,700,500);  // gives the page size		
	canvasZdist->cd();
		canvasZdist->SetTopMargin(0+0.02);
		canvasZdist->SetRightMargin(0+0.01);
		StylingSliceHistos(Resolution_dZAbs_VS_Rrebin_2,0.8);
		Resolution_dZAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
		Resolution_dZAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dZAbs_VS_Rrebin_2, 
						"", "R (cm)", "#sigma dZ (cm)", 
							kFALSE, 10., 140.,
							kFALSE, 0., 2.1,
							kFALSE, 0., 120.);
	canvasZdist->Update();	
	canvasZdist->SaveAs(Form("%s/Resolution_dZvsRdistributionsingle.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasZdist;

	TCanvas * canvasdZZdist = new TCanvas("canvasdZZdist","",10,10,700,500);  // gives the page size		
	canvasdZZdist->cd();
		canvasdZZdist->SetTopMargin(0+0.02);
		canvasdZZdist->SetRightMargin(0+0.01);
		StylingSliceHistos(Resolution_dZAbs_VS_Zrebin_2,0.8);
		Resolution_dZAbs_VS_Zrebin_2->SetMarkerColor(kRed+2);
		Resolution_dZAbs_VS_Zrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dZAbs_VS_Zrebin_2, 
						"", "Z (cm)", "#sigma dZ (cm)", 
							kFALSE, 10., 140.,
							kFALSE, 0., 2.1,
							kFALSE, 0., 120.);
	canvasdZZdist->Update();	
	canvasdZZdist->SaveAs(Form("%s/Resolution_dZvsZdistributionsingle.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasdZZdist;

	TCanvas * canvasPhiDist = new TCanvas("canvasPhiDist","",10,10,700,500);  // gives the page size		
	canvasPhiDist->cd();
		canvasPhiDist->SetTopMargin(0+0.055);
		canvasPhiDist->SetRightMargin(0+0.01);
		StylingSliceHistos(Resolution_dPhiAbs_VS_Rrebin_2,0.8);
		Resolution_dPhiAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
		Resolution_dPhiAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_2, 
						"", "R (cm)", "#sigma d#phi (rad)", 
							kFALSE, 10., 140.,
							kFALSE, 0.0, 0.01,
							kFALSE, 0., 120.);
	canvasPhiDist->Update();	
	canvasPhiDist->SaveAs(Form("%s/Resolution_Phidistributionsingle.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasPhiDist;

	TCanvas * canvasRDist = new TCanvas("canvasRDist","",10,10,700,500);  // gives the page size		
	canvasRDist->cd();
		canvasRDist->SetTopMargin(0+0.02);
		canvasRDist->SetRightMargin(0+0.01);
		StylingSliceHistos(Resolution_dRAbs_VS_Rrebin_2,0.8);
		Resolution_dRAbs_VS_Rrebin_2->SetMarkerColor(kRed+2);
		Resolution_dRAbs_VS_Rrebin_2->SetLineColor(kBlue-8);
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_2, 
						"", "R (cm)", "#sigma dR (cm)", 
							kFALSE, 10., 140.,
							kFALSE, 0., 3.1,
							kFALSE, 0., 120.);
	canvasRDist->Update();	
	canvasRDist->SaveAs(Form("%s/Resolution_Rdistributionsingle.%s",outputDirectory.Data(),suffix.Data()));
	delete canvasRDist;

	
	// ----------------------------- Resolution dPhiAbs vs R ---------------------------------------
	TCanvas * canvasResFullSigma = new TCanvas("canvasResFullSigma","",10,10,1800,1000);  // gives the page size		
	canvasResFullSigma->cd();
	DrawGammaCanvasSettings( canvasResFullSigma, 0.13, 0.0, 0.02, 0.09);
	
	TPad* padResFull1 = new TPad("padResFull1", "",  0., 0.52, 0.35, 1.,-1, -1, -2);
	DrawGammaPadSettings( padResFull1, 0.1, 0.0, 0.0, 0.);
	padResFull1->Draw();
	TPad* padResFull2 = new TPad("padResFull2", "",0., 0., 0.35, 0.52,-1, -1, -2);
	DrawGammaPadSettings( padResFull2, 0.1, 0.00, 0., 0.1);
	padResFull2->Draw();
	
	TPad* padResFull3 = new TPad("padResFull3", "", 0.35, 0.52, 0.675, 1.,-1, -1, -2);
	DrawGammaPadSettings( padResFull3, 0.0, 0.0, 0.0, 0.);
	padResFull3->Draw();
	TPad* padResFull4 = new TPad("padResFull4", "", 0.35, 0., 0.675, 0.52,-1, -1, -2);
	DrawGammaPadSettings( padResFull4, 0.0, 0.0, 0., 0.1);
	padResFull4->Draw();
	
	TPad* padResFull5 = new TPad("padResFull5", "", 0.675, 0.52, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padResFull5, 0.0, 0.02, 0.0, 0.);
	padResFull5->Draw();
	TPad* padResFull6 = new TPad("padResFull6", "", 0.675, 0., 1., 0.52,-1, -1, -2);
	DrawGammaPadSettings( padResFull6, 0.0, 0.02, 0., 0.1);
	padResFull6->Draw();
	
	
	padResFull1->cd();
// 	padResFull1->SetTopMargin(0.05);
		DrawGammaSetMarker(Resolution_dPhiAbs_VS_Rrebin_2, 21,1., kRed+2 , kRed-9);
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_2, 
							"", "R (cm)", "#sigma (#mu rad)", 
							kFALSE, 10., 140.,
							kTRUE, -0.001 , maxdPhiAbsSigmaR,
							kTRUE, 0., 120.);
// 		DrawGammaLines(0., 120, 0,0,0.005);
	padResFull1->Update();
	padResFull2->cd();
		DrawGammaSetMarker(Resolution_dRAbs_VS_Rrebin_2, 20,1., kBlue+2 , kBlue-9);
		DrawGammaSetMarker(Resolution_dZAbs_VS_Rrebin_2, 22,1., kGreen+2 , kGreen-9);
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_2, 
						"", "R (cm)", "#sigma (cm)", 
							kFALSE, 10., 140.,
							kTRUE, 0.0, maxdRAbsSigmaR,
							kTRUE, 0., 120.);
		Resolution_dZAbs_VS_Rrebin_2->Draw("same,e1");
		Resolution_dRAbs_VS_Rrebin_2->Draw("same,e1");
	padResFull2->Update();

	padResFull3->cd();
// 	padResFull3->SetTopMargin(0.05);
	
	TH1D * histoDummyZResPhi = new TH1D("histoDummyZResPhi","histoDummyZResPhi",500,-250., 250.);
	TH2D * histoDummyZResPhi2 = new TH2D("histoDummyZResPhi","histoDummyZResPhi",500,-220., 220., 500, -0.001, maxdPhiAbsSigmaR);
	
	histoDummyZResPhi2->GetYaxis()->SetTitleSize(0.055);	
	histoDummyZResPhi2->GetYaxis()->SetLabelSize(0.045);
	histoDummyZResPhi2->GetYaxis()->SetDecimals();
	histoDummyZResPhi2->GetYaxis()->SetTitleOffset(0.9);
	histoDummyZResPhi2->GetXaxis()->SetTitleOffset(0.9);
	histoDummyZResPhi2->GetXaxis()->SetTitleSize(0.055);
	histoDummyZResPhi2->GetXaxis()->SetLabelSize(0.045);	
	histoDummyZResPhi2->SetTitle("");
	histoDummyZResPhi2->Draw();
	
	DrawGammaSetMarker(Resolution_dPhiAbs_VS_Zrebin_2, 21,1., kRed+2 , kRed-9);
// 	DrawResolutionGammaHisto( histoDummyZResPhi, 
// 							"", "Z (cm)", "#sigma (rad)", 
// 							kFALSE, 10., 140.,
// 							kTRUE, -0.002 , maxdPhiAbsSigmaR,
// 							kTRUE, -220., 220.);
// 		DrawGammaLines(0., 120, 0,0,0.005);
	Resolution_dPhiAbs_VS_Zrebin_2->Draw("same,e1");
	padResFull3->Update();
	padResFull4->cd();
		DrawGammaSetMarker(Resolution_dRAbs_VS_Zrebin_2, 20,1., kBlue+2 , kBlue-9);
		DrawGammaSetMarker(Resolution_dZAbs_VS_Zrebin_2, 22,1., kGreen+2 , kGreen-9);
		DrawResolutionGammaHisto( histoDummyZResPhi, 
						"", "Z (cm)", "#sigma (cm)", 
							kFALSE, 10., 140.,
							kTRUE, 0., maxdRAbsSigmaR,
							kTRUE, -220., 220.);
		Resolution_dZAbs_VS_Zrebin_2->Draw("same,e1");
		Resolution_dRAbs_VS_Zrebin_2->Draw("same,e1");
	padResFull4->Update();

	padResFull5->cd();
// 	padResFull5->SetTopMargin(0.05);
		TMarker* markerdZRes = CreateMarkerFromHisto(Resolution_dZAbs_VS_Zrebin_2,1,1 ,1.8);
		TMarker* markerdRRes = CreateMarkerFromHisto(Resolution_dRAbs_VS_Zrebin_2,1,1 ,1.8);
		TMarker* markerdPhiRes = CreateMarkerFromHisto(Resolution_dPhiAbs_VS_Zrebin_2,1,1 ,1.8);
		
		TLegend* legendAllSigma = new TLegend( 0.6,0.65,0.92,0.9);
		legendAllSigma->SetTextSize(0.06);			
		legendAllSigma->SetFillColor(0);
		legendAllSigma->SetLineColor(0);
		legendAllSigma->AddEntry(markerdRRes,"#sigma_{dR} (cm)","p");
		legendAllSigma->AddEntry(markerdZRes,"#sigma_{dZ} (cm)","p");
		legendAllSigma->AddEntry(markerdPhiRes,"#sigma_{d#phi} (#mu rad)", "p");		

		DrawGammaSetMarker(Resolution_dPhiAbs_VS_Phi_2, 21,1., kRed+2 , kRed-9);
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Phi_2, 
							"", "Z (cm)", "#sigma (rad)", 
							kFALSE, 10., 140.,
							kTRUE, -0.001 , maxdPhiAbsSigmaR,
							kTRUE, 0, 2*TMath::Pi());	
// 		DrawGammaLines(0., 120, 0,0,0.005);
		legendAllSigma->Draw("same");
	padResFull5->Update();
	padResFull6->cd();
		DrawGammaSetMarker(Resolution_dRAbs_VS_Phi_2, 20,1., kBlue+2 , kBlue-9);
		DrawGammaSetMarker(Resolution_dZAbs_VS_Phi_2, 22,1., kGreen+2 , kGreen-9);
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Phi_2, 
						"", "#phi (rad)", "#sigma (cm)", 
							kFALSE, 10., 140.,
							kTRUE, 0., maxdRAbsSigmaR,
							kTRUE, 0, 2*TMath::Pi());	
		Resolution_dZAbs_VS_Phi_2->Draw("same,e1");
		Resolution_dRAbs_VS_Phi_2->Draw("same,e1");
	padResFull6->Update();

	canvasResFullSigma->Update();	
	canvasResFullSigma->SaveAs(Form("%s/Resolution_Sigma_AllVsAll.%s",outputDirectory.Data(),suffix.Data()));
	delete padResFull1;	
	delete padResFull2;	
	delete padResFull3;	
	delete padResFull4;	
	delete padResFull5;	
	delete padResFull6;	
	delete canvasResFullSigma;

	TCanvas * canvasResFullMean = new TCanvas("canvasResFullMean","",10,10,1800,1000);  // gives the page size		
	canvasResFullMean->cd();
	DrawGammaCanvasSettings( canvasResFullMean, 0.13, 0.0, 0.0, 0.09);
	
	TPad* padResFullMean1 = new TPad("padResFullMean1", "", 0., 0.52, 0.35, 1.,-1, -1, -2);
	DrawGammaPadSettings( padResFullMean1, 0.1, 0.0, 0.0, 0.);
	padResFullMean1->Draw();
	TPad* padResFullMean2 = new TPad("padResFullMean2", "", 0., 0., 0.35, 0.52,-1, -1, -2);
	DrawGammaPadSettings( padResFullMean2, 0.1, 0.00, 0., 0.1);
	padResFullMean2->Draw();
	
	TPad* padResFullMean3 = new TPad("padResFullMean3", "", 0.35, 0.52, 0.675, 1.,-1, -1, -2);
	DrawGammaPadSettings( padResFullMean3, 0.0, 0.0, 0.0, 0.);
	padResFullMean3->Draw();
	TPad* padResFullMean4 = new TPad("padResFullMean4", "", 0.35, 0., 0.675, 0.52,-1, -1, -2);
	DrawGammaPadSettings( padResFullMean4, 0.0, 0.0, 0., 0.1);
	padResFullMean4->Draw();
	
	TPad* padResFullMean5 = new TPad("padResFullMean5", "", 0.675, 0.52, 1., 1.,-1, -1, -2);
	DrawGammaPadSettings( padResFullMean5, 0.0, 0.00, 0.0, 0.);
	padResFullMean5->Draw();
	TPad* padResFullMean6 = new TPad("padResFullMean6", "", 0.675, 0., 1., 0.52,-1, -1, -2);
	DrawGammaPadSettings( padResFullMean6, 0.0, 0.0, 0., 0.1);
	padResFullMean6->Draw();
	
	
	padResFullMean1->cd();
// 	padResFullMean1->SetTopMargin(0.05);
		DrawGammaSetMarker(Resolution_dPhiAbs_VS_Rrebin_1, 21,1., kRed+2 , kRed-9);
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Rrebin_1, 
							"", "R (cm)", "#mu (#mu rad)", 
							kFALSE,  10., 140.,
							kTRUE, arraydPhiAbsMeanR[0] , arraydPhiAbsMeanR[1],
							kTRUE, 0., 120.);
// 		DrawGammaLines(0., 120, 0,0,0.005);
	padResFullMean1->Update();
	padResFullMean2->cd();
		DrawGammaSetMarker(Resolution_dRAbs_VS_Rrebin_1, 20,1., kBlue+2 , kBlue-9);
		DrawGammaSetMarker(Resolution_dZAbs_VS_Rrebin_1, 22,1., kGreen+2 , kGreen-9);
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Rrebin_1, 
						"", "R (cm)", "#mu (cm)", 
							kFALSE, 10., 140.,
							kTRUE, arraydRAbsMeanR[0] , arraydRAbsMeanR[1],
							kTRUE, 0., 120.);
		Resolution_dZAbs_VS_Rrebin_1->Draw("same,e1");
		Resolution_dRAbs_VS_Rrebin_1->Draw("same,e1");
	padResFullMean2->Update();

	padResFullMean3->cd();
// 	padResFullMean3->SetTopMargin(0.05);
	
	TH2D * histoDummyZResPhi4 = new TH2D("histoDummyZResPhi4","histoDummyZResPhi4",500,-220., 220.,500, arraydRAbsMeanR[0] , arraydRAbsMeanR[1]);
	TH2D * histoDummyZResPhi3 = new TH2D("histoDummyZResPhi3","histoDummyZResPhi3",500,-220., 220., 500, arraydPhiAbsMeanR[0] , arraydPhiAbsMeanR[1]);
	
	histoDummyZResPhi3->GetYaxis()->SetTitleSize(0.055);	
	histoDummyZResPhi3->GetYaxis()->SetLabelSize(0.045);
	histoDummyZResPhi3->GetYaxis()->SetDecimals();
	histoDummyZResPhi3->GetYaxis()->SetTitleOffset(0.9);
	histoDummyZResPhi3->GetXaxis()->SetTitleOffset(0.9);
	histoDummyZResPhi3->GetXaxis()->SetTitleSize(0.055);
	histoDummyZResPhi3->GetXaxis()->SetLabelSize(0.045);	
	histoDummyZResPhi3->SetTitle("");
	histoDummyZResPhi3->Draw();
	
	TLegend* legendAllMean = new TLegend( 0.6,0.65,0.92,0.9);
	legendAllMean->SetTextSize(0.06);			
	legendAllMean->SetFillColor(0);
	legendAllMean->SetLineColor(0);
	legendAllMean->AddEntry(markerdRRes,"#mu_{dR} (cm)","p");
	legendAllMean->AddEntry(markerdZRes,"#mu_{dZ} (cm)","p");
	legendAllMean->AddEntry(markerdPhiRes,"#mu_{d#phi} (#mu rad)", "p");		

	DrawGammaSetMarker(Resolution_dPhiAbs_VS_Zrebin_1, 21,1., kRed+2 , kRed-9);
// 	DrawResolutionGammaHisto( histoDummyZResPhi, 
// 							"", "Z (cm)", "#sigma (rad)", 
// 							kFALSE, 10., 140.,
// 							kTRUE, -0.002 , maxdPhiAbsSigmaR,
// 							kTRUE, -220., 220.);
// 		DrawGammaLines(0., 120, 0,0,0.005);
	Resolution_dPhiAbs_VS_Zrebin_1->Draw("same,e1");
	legendAllMean->Draw("same");
	padResFullMean3->Update();
	padResFullMean4->cd();
	
	histoDummyZResPhi4->GetYaxis()->SetTitleSize(0.055);	
	histoDummyZResPhi4->GetYaxis()->SetLabelSize(0.045);
	histoDummyZResPhi4->GetYaxis()->SetDecimals();
	histoDummyZResPhi4->GetYaxis()->SetTitleOffset(0.9);
	histoDummyZResPhi4->GetXaxis()->SetTitleOffset(0.9);
	histoDummyZResPhi4->GetXaxis()->SetTitleSize(0.055);
	histoDummyZResPhi4->GetXaxis()->SetLabelSize(0.045);
	histoDummyZResPhi4->GetXaxis()->SetTitle("Z (cm)");
	histoDummyZResPhi4->GetYaxis()->SetTitle("#mu (cm)");
	histoDummyZResPhi4->SetTitle("");
	histoDummyZResPhi4->Draw();

		DrawGammaSetMarker(Resolution_dRAbs_VS_Zrebin_1, 20,1., kBlue+2 , kBlue-9);
		DrawGammaSetMarker(Resolution_dZAbs_VS_Zrebin_1, 22,1., kGreen+2 , kGreen-9);
	/*	DrawResolutionGammaHisto( histoDummyZResPhi4, 
						"", "Z (cm)", "#mu (cm)", 
							kFALSE, 10., 140.,
							kTRUE, arraydRAbsMeanR[0] , arraydRAbsMeanR[1],
							kTRUE, -220., 220.);
	*/	Resolution_dZAbs_VS_Zrebin_1->Draw("same,e1");
		Resolution_dRAbs_VS_Zrebin_1->Draw("same,e1");
	padResFullMean4->Update();

	padResFullMean5->cd();
// 	padResFullMean5->SetTopMargin(0.05);
		
		DrawGammaSetMarker(Resolution_dPhiAbs_VS_Phi_1, 21,1., kRed+2 , kRed-9);
		DrawResolutionGammaHisto( Resolution_dPhiAbs_VS_Phi_1, 
							"", "Z (cm)", "#mu (#mu rad)", 
							kFALSE, 10., 140.,
							kTRUE,arraydPhiAbsMeanR[0] , arraydPhiAbsMeanR[1],	
							kTRUE, 0, 2*TMath::Pi());	
// 		DrawGammaLines(0., 120, 0,0,0.005);

	padResFullMean5->Update();
	padResFullMean6->cd();
		DrawGammaSetMarker(Resolution_dRAbs_VS_Phi_1, 20,1., kBlue+2 , kBlue-9);
		DrawGammaSetMarker(Resolution_dZAbs_VS_Phi_1, 22,1., kGreen+2 , kGreen-9);
		DrawResolutionGammaHisto( Resolution_dRAbs_VS_Phi_1, 
						"", "#phi (rad)", "#mu (cm)", 
							kFALSE, 10., 140.,
							kTRUE, arraydRAbsMeanR[0] , arraydRAbsMeanR[1],
							kTRUE, 0, 2*TMath::Pi());	
		Resolution_dZAbs_VS_Phi_1->Draw("same,e1");
		Resolution_dRAbs_VS_Phi_1->Draw("same,e1");
	padResFullMean6->Update();

	canvasResFullMean->Update();	
	canvasResFullMean->SaveAs(Form("%s/Resolution_Mean_AllVsAll.%s",outputDirectory.Data(),suffix.Data()));
	delete padResFullMean1;	
	delete padResFullMean2;	
	delete padResFullMean3;	
	delete padResFullMean4;	
	delete padResFullMean5;	
	delete padResFullMean6;	
	delete canvasResFullMean;

	
	ProduceProjectionsPlotWithFits( Resolution_Z_dPhiRes ,fitResolution_Z_dPhiRes, Zbinning, 50 ,9,6, "dPhi_{abs}", -10., 10., "N_{#gamma}","#gamma res dZ" , "fit dZ", "Z",Form("%s/ProjectionsFitteddPhiAbsZ.%s",outputDirectory.Data(),suffix.Data()),4*2800, 4*1800);

	ProduceProjectionsPlotWithFits( Resolution_R_dRRes ,fitResolution_R_dRRes, Rbinning, 85 ,11,8 , "dR_{abs}", -10., 10., "N_{#gamma}","#gamma res dR" , "fit dR", "R",Form("%s/ProjectionsFitteddRAbsR.%s",outputDirectory.Data(),suffix.Data()),4*2800, 4*1800);
	ProduceProjectionsPlotWithFits( Resolution_R_dPhiAbsRes ,fitResolution_R_dPhiAbsRes, Rbinning, 85 ,11,8 , "d#Phi_{abs}", -10., 10., "N_{#gamma}","#gamma res d#Phi" , "fit d#Phi", "R",Form("%s/ProjectionsFitteddPhiAbsR.%s",outputDirectory.Data(),suffix.Data()),4*2800, 4*1800);
	ProduceProjectionsPlotWithFits( Resolution_R_dZRes ,fitResolution_R_dZRes, Rbinning, 85 ,11,8 , "dZ_{abs}", -10., 10., "N_{#gamma}","#gamma res dZ" , "fit dZ", "R",Form("%s/ProjectionsFitteddZAbsR.%s",outputDirectory.Data(),suffix.Data()),4*2800, 4*1800);
		
	ProduceProjectionsPlotWithFits( Resolution_Z_dZRes ,fitResolution_Z_dZRes, Zbinning, 50 ,9,6, "dZ_{abs}", -10., 10., "N_{#gamma}","#gamma res dZ" , "fit dZ", "Z",Form("%s/ProjectionsFitteddZAbsZ.%s",outputDirectory.Data(),suffix.Data()),4*2800, 4*1800);

	ProduceProjectionsPlotWithFits( Resolution_Pt_dRRes ,fitResolution_Pt_dRRes, ptbinning, 21 ,6,4 , "dR_{abs}", -10., 10., "N_{#gamma}","#gamma res dR" , "fit dR", "Pt",Form("%s/ProjectionsFitteddRAbsPt.%s",outputDirectory.Data(),suffix.Data()),4*2800, 4*1800);
	ProduceProjectionsPlotWithFits( Resolution_Pt_dZRes ,fitResolution_Pt_dZRes, ptbinning, 21 ,6,4 , "dZ_{abs}", -10., 10., "N_{#gamma}","#gamma res dZ" , "fit dZ", "Pt",Form("%s/ProjectionsFitteddZAbsPt.%s",outputDirectory.Data(),suffix.Data()),4*2800, 4*1800);

	
}


