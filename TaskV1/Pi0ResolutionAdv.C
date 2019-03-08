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
TString detectionProcess;

//**********************************************************************************
//******************** Helper function for 2D histo plotting ***********************
//**********************************************************************************
void DrawAutoGammaHistoPaper2D( TH2* histo1,
                    TString Title, TString XTitle, TString YTitle,
                    Bool_t YRangeMax, Double_t YMaxFactor, Double_t YMinimum,
                    Bool_t YRange, Double_t YMin ,Double_t YMax,
                    Bool_t XRange, Double_t XMin, Double_t XMax, Double_t xOffset, Double_t yOffset) {
    if (YRangeMax && !XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
    }
    if (YRangeMax && XRange){
        YRange = kFALSE;
        Double_t maxRangeR = histo1->GetMaximum();
        Double_t minRangeR = histo1->GetMinimum();
        if(YMinimum > minRangeR){minRangeR = YMinimum;}
        histo1->GetYaxis()->SetRangeUser(minRangeR, maxRangeR*YMaxFactor);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (YRange && XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }
    if (!YRangeMax && !YRange && XRange){
        histo1->GetXaxis()->SetRangeUser(XMin, XMax);
    }

    if (YRange && !XRange){
        histo1->GetYaxis()->SetRangeUser(YMin, YMax);
    }

    histo1->SetTitle(Title.Data());

    if(XTitle.CompareTo("") != 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.CompareTo("") != 0){
        histo1->SetYTitle(YTitle.Data());
    }

    histo1->GetYaxis()->SetLabelFont(42);
    histo1->GetXaxis()->SetLabelFont(42);
    histo1->GetYaxis()->SetTitleFont(62);
    histo1->GetXaxis()->SetTitleFont(62);
    histo1->GetYaxis()->SetLabelSize(0.045);
    histo1->GetYaxis()->SetTitleSize(0.05);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetTitleOffset(xOffset);
    histo1->GetYaxis()->SetTitleOffset(yOffset);
    histo1->GetXaxis()->SetTitleSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.045);
}

//**********************************************************************************
//******************** Plotting 2D resolution matrix *******************************
//**********************************************************************************
void PlotStandard2D( TH2* histo2D, 
                     TString nameOutput,
                     TString title, 
                     TString xTitle, 
                     TString yTitle, 
                     Bool_t kRangeYAutoMax,
                     Double_t maxFacY,
                     Double_t startYAuto,
                     Bool_t kRangeY, 
                     Double_t startY, 
                     Double_t endY, 
                     Bool_t kRangeX, 
                     Double_t startX, 
                     Double_t endX,                
                     Double_t startZ,
                     Int_t logX, 
                     Int_t logZ, 
                     Float_t* floatLogo, 
                     Double_t offsetX,
                     Double_t offsetY,
                     TString additionalLabel  = "",
                     Int_t canvasSizeX        = 500, 
                     Int_t canvasSizeY        = 500
                   ){
  
    TCanvas *canvasStandard = new TCanvas("canvasStandard","canvasStandard",canvasSizeX,canvasSizeY);
    DrawGammaCanvasSettings( canvasStandard, 0.115, 0.11, 0.02, 0.11); 
    canvasStandard->SetLogx(logX);
    canvasStandard->SetLogz(logZ);
    canvasStandard->cd();
    
    DrawAutoGammaHistoPaper2D(	histo2D,
                                title.Data(), xTitle.Data(), yTitle.Data(),
                                kRangeYAutoMax, maxFacY, startYAuto, 
                                kRangeY, startY, endY, 
                                kRangeX, startX, endX,
                                offsetX, offsetY                               
                             );
    histo2D->GetZaxis()->SetRangeUser(startZ, histo2D->GetMaximum());

    histo2D->SetTitle(title.Data());
    histo2D->Draw("colz");
    
    PutProcessLabelAndEnergyOnPlot(floatLogo[0], floatLogo[1], 28, collisionSystem.Data(), detectionProcess.Data(), additionalLabel.Data(), 63, 0.03);

    canvasStandard->Update();
    canvasStandard->SaveAs(nameOutput.Data());
    delete canvasStandard;
}

//**********************************************************************************
//******************** Plotting resolution slices **********************************
//**********************************************************************************
void ProduceProjectionsPlotWithFits(  TH1D** histos, 
                                      TF1** fitData, 
                                      Double_t* floatArray, 
                                      Int_t nBins, 
                                      Int_t nColumns, 
                                      Int_t nRows, 
                                      TString xTitle, 
                                      TString yTitle, 
                                      TString legendEntry1, 
                                      TString legendEntry2, 
                                      TString title, 
                                      TString nameOutput, 
                                      Int_t mode              = 0,
                                      Int_t canvasSizeX       = 2800,
                                      Int_t canvasSizeY       = 1800
                                      
                                   ){

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
            alice->SetTextColor(1);
            alice->SetTextSize(textHeight);
            alice->Draw();

            process->SetNDC(kTRUE);
            process->SetTextSize(textHeight);
            process->Draw();

            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textHeight);
            latexDate->Draw();
            
            TLegend* legendProjections = new TLegend(0.0,0.1,1.,0.3);
            legendProjections->SetTextSize(0.08);			
            legendProjections->SetFillColor(0);
            legendProjections->SetBorderSize(0);
            legendProjections->AddEntry(histos[3],legendEntry1.Data(),"pe");
            legendProjections->AddEntry(fitData[3],legendEntry2.Data(),"l");
            legendProjections->Draw();            
        } else {
            padProjection->cd(place);
            padProjection->cd(place)->SetTopMargin(0.12);
            padProjection->cd(place)->SetBottomMargin(0.15);
            padProjection->cd(place)->SetRightMargin(0.05);
            padProjection->cd(place)->SetLeftMargin(0.15);
            padProjection->cd(place)->SetLogy(0);
            Float_t startRange = -0.5;
            Float_t endRange = 0.5;
            
            Double_t yMin = 1.;
            if (mode == 10) yMin = 0;
            if ( histos[iPt]->GetEntries() > 0){
                if (title.CompareTo("Pt")==0){
                    DrawGammaHistoWithTitle2( histos[iPt],
                        Form("%3.2f GeV/c < #it{p}_{T} < %3.2f GeV/#it{c}",startPt,endPt), 
                        xTitle.Data(),yTitle.Data(),
                        startRange,endRange,yMin);
                }
                if (fitData[iPt]!=0x0){
                    fitData[iPt]->SetLineColor(kBlue);
                    fitData[iPt]->SetLineWidth(1);
                    fitData[iPt]->Draw("same");
                } else {
                   cout << "fit failed" << endl;
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


//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn2D(TH2* histo){
    Double_t minimum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        for (Int_t j = 1; j<histo->GetNbinsY(); j++){
            if (histo->GetBinContent(i,j) < minimum && histo->GetBinContent(i,j) > 0){
                minimum = histo->GetBinContent(i,j);
            }
        }
    }
    return minimum;
}

//**********************************************************************************
//********************** Main function for resolution ******************************
//**********************************************************************************
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

    if (mesonName.CompareTo("Pi0") == 0)
        mesonLatex ="#pi^{0}"; 
    else if (mesonName.CompareTo("Eta") == 0)
        mesonLatex ="#eta"; // 
    
    collisionSystem                 = ReturnFullCollisionsSystem(optEnergy);
    detectionProcess                = ReturnFullTextReconstructionProcess(mode);
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;     
    }
    TString fTextMeasurement    = Form("%s #rightarrow #gamma#gamma", mesonLatex.Data());
    
    TString centrality = "";
    TString firstCutnumber = optionCutSelection(GetEventSystemCutPosition(),1);
    if (firstCutnumber.CompareTo("0") != 0){
        centrality = GetCentralityString(optionCutSelection);
        collisionSystem = Form("%s %s", centrality.Data(), collisionSystem.Data());
    }   
    cout << centrality.Data() << endl;
    
    
    if(optMCGenerator.CompareTo("") ==0){
        textGenerator = "";
    } else {
        textGenerator = optMCGenerator;
    }

    textDate = ReturnDateString();
    
    
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
    Float_t floatLocationRightDown2D[4]=	{0.15,0.25,0.11, 0.02};

    //**********************************************************************************************************************
    //******************************** Defintion of arrays for rebinning ***************************************************
    //**********************************************************************************************************************

    //pt rebinning array
    
    Int_t maxNBinsPt        = 50;
    Int_t nColumns          = 9;
    Int_t nRows             = 6;
    Double_t ptbinning[50];
    Int_t ptbinningReb[49];
    
    if (mode != 10 ){
        maxNBinsPt          = 18;   
        nColumns            = 5;
        nRows               = 4;
        Double_t ptbinningStandard[19]  = { 0.4, 0.6, 0.8, 1.0, 1.2,
                                            1.4, 1.6, 1.8, 2.0, 2.4,
                                            2.8, 3.2, 3.6, 4.0, 4.5,
                                            5.0, 6.0, 10., 15.};
        Int_t ptbinningRebStandard[18]  = { 4, 4, 4, 4, 4,
                                            4, 4, 4, 4, 4,
                                            4, 4, 4, 4, 8,
                                            8, 10,10 };  
        for (Int_t i = 0; i < maxNBinsPt+1; i++){
            if (i < maxNBinsPt) 
                ptbinningReb[i]      = ptbinningRebStandard[i];
            ptbinning[i]             = ptbinningStandard[i];
        }            
    } else {
        maxNBinsPt          = 23;
        nColumns            = 5;
        nRows               = 5;
        Double_t ptbinningStandard[24]  = { 10., 11., 12., 13., 14.,
                                            15., 16., 18., 20., 22., 
                                            24., 26., 28., 30., 35., 
                                            40., 45., 50 , 60., 70.,
                                            80., 100.,150.,200};
        Int_t ptbinningRebStandard[23]  = { 10, 10, 10, 10, 10,
                                            10, 10, 10, 10, 10,
                                            10, 10, 10, 10, 10,
                                            10, 10, 10, 10, 10,
                                            10, 10, 10};
        for (Int_t i = 0; i < maxNBinsPt+1; i++){
            if (i < maxNBinsPt) 
                ptbinningReb[i]      = ptbinningRebStandard[i];
            ptbinning[i]             = ptbinningStandard[i];
        }                            
    }    
                              
                              
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
    TFile* fileMC                   = new TFile(fileMonteCarloInput);  
    TString autoDetectedMainDir     = AutoDetectMainTList(mode , fileMC);
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }

    //************************** Container Loading ********************************************************************
    TList *TopDir                   = (TList*)fileMC->Get(autoDetectedMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    TList *HistosGammaConversion    = (TList*)TopDir->FindObject(Form("Cut Number %s",optionCutSelection.Data()));
    //    TList *MCContainer = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",optionCutSelection.Data()));
    TList *TrueConversionContainer  = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",optionCutSelection.Data()));

    //************************** Histo loading *************************************************************************
//             nameResolHist[0]            = Form("ESD_TruePrimary%sPureMerged_MCPt_ResolPt",outputlabel.Data());
//             nameResolHist[1]            = Form("ESD_TruePrimary%sMergedPartConv_MCPt_ResolPt",outputlabel.Data());
//             nameResolHist[2]            = Form("ESD_TruePrimary%s1Gamma_MCPt_ResolPt",outputlabel.Data());
//             nameResolHist[3]            = Form("ESD_TruePrimary%s1Electron_MCPt_ResolPt",outputlabel.Data());

    //****************************** Resolution dPt vs Pt **********************************************
    TH2F * Resolution_Pi0_MCPt_ResolPt              = (TH2F*)TrueConversionContainer->FindObject(Form("ESD_TruePrimary%s_MCPt_ResolPt",mesonName.Data()) ); //
    if (Resolution_Pi0_MCPt_ResolPt) Resolution_Pi0_MCPt_ResolPt->Sumw2();
    
    TH2F* Resolution_Pi0_MCPt_ResolPt_Merged        = NULL;
    TH2F* Resolution_Pi0_MCPt_ResolPt_PartMerged    = NULL;
    TH2F* Resolution_Pi0_MCPt_ResolPt_OneGamma         = NULL;
    TH2F* Resolution_Pi0_MCPt_ResolPt_OneElectron      = NULL;
    if ( mode == 10 || mode == 11 ){
        Resolution_Pi0_MCPt_ResolPt_Merged          = (TH2F*)TrueConversionContainer->FindObject(Form("ESD_TruePrimary%sPureMerged_MCPt_ResolPt",mesonName.Data()) ); //
        Resolution_Pi0_MCPt_ResolPt_Merged->Sumw2();
        Resolution_Pi0_MCPt_ResolPt_PartMerged      = (TH2F*)TrueConversionContainer->FindObject(Form("ESD_TruePrimary%sMergedPartConv_MCPt_ResolPt",mesonName.Data()));
        Resolution_Pi0_MCPt_ResolPt_PartMerged->Sumw2();
        Resolution_Pi0_MCPt_ResolPt_OneGamma        = (TH2F*)TrueConversionContainer->FindObject(Form("ESD_TruePrimary%s1Gamma_MCPt_ResolPt",mesonName.Data()));
        Resolution_Pi0_MCPt_ResolPt_OneGamma->Sumw2();
        Resolution_Pi0_MCPt_ResolPt_OneElectron     = (TH2F*)TrueConversionContainer->FindObject(Form("ESD_TruePrimary%s1Electron_MCPt_ResolPt",mesonName.Data()));
        Resolution_Pi0_MCPt_ResolPt_OneElectron->Sumw2();
        Resolution_Pi0_MCPt_ResolPt                 = (TH2F*)Resolution_Pi0_MCPt_ResolPt_Merged->Clone(Form("ESD_TruePrimary%s_MCPt_ResolPt",mesonName.Data()) ); //
        Resolution_Pi0_MCPt_ResolPt->Sumw2();
        Resolution_Pi0_MCPt_ResolPt->Add(Resolution_Pi0_MCPt_ResolPt_PartMerged);
        Resolution_Pi0_MCPt_ResolPt->Add(Resolution_Pi0_MCPt_ResolPt_OneGamma);
        Resolution_Pi0_MCPt_ResolPt->Add(Resolution_Pi0_MCPt_ResolPt_OneElectron);
    }    
    cout << "hier noch " << Resolution_Pi0_MCPt_ResolPt<<endl;
     
    //**********************************************************************************************************************************************
    //*********************************** Rebinning of 2D histograms and fitslices *****************************************************************
    //**********************************************************************************************************************************************
    Int_t maxYbin = 0;
    // 	Int_t maxXbin = 0;
    Double_t precision = 10E-4;
      
    TH1D* Resolution_Pi0_PtRes[maxNBinsPt];
    TF1*  fitResolution_Pi0_PtRes[maxNBinsPt];
    TH1F *Resolution_Meson_Mean     = new TH1F("Resolution_Meson_Mean", "mean meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ;	
    TH1F *Resolution_Meson_Sigma    = new TH1F("Resolution_Meson_Sigma", "sigma meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ;	
    TH1F *Resolution_Meson_Mean2    = new TH1F("Resolution_Meson_Mean2", "mean meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ;  
    TH1F *Resolution_Meson_Sigma2   = new TH1F("Resolution_Meson_Sigma2", "RMS meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ; 

    Double_t minFitRange    = -0.5;
    Double_t maxFitRange    = 0.5;
    if (mode == 10){
        minFitRange         = -0.25;
        maxFitRange         = 0.25;       
    }   
    
    ResolutionFittingRebined( Resolution_Pi0_MCPt_ResolPt, Resolution_Pi0_PtRes, fitResolution_Pi0_PtRes, maxNBinsPt, ptbinning, ptbinningReb, Resolution_Meson_Mean, Resolution_Meson_Sigma, "gaus", 
                              minFitRange, maxFitRange, precision, "Resolution_Pi0_PtRes");

    Resolution_Meson_Sigma->Scale(100);
    Resolution_Meson_Mean->Scale(100);
      
    for (Int_t i = 0; i < maxNBinsPt; i++){
        Resolution_Meson_Mean2->SetBinContent(i+1, Resolution_Pi0_PtRes[i]->GetMean());
        Resolution_Meson_Mean2->SetBinError(i+1, Resolution_Pi0_PtRes[i]->GetMeanError());
        Resolution_Meson_Sigma2->SetBinContent(i+1, Resolution_Pi0_PtRes[i]->GetRMS());
        Resolution_Meson_Sigma2->SetBinError(i+1, Resolution_Pi0_PtRes[i]->GetRMSError());
    }
    Resolution_Meson_Sigma2->Scale(100.);
    Resolution_Meson_Mean2->Scale(100.);
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
    TH2F * histo2DMean = new TH2F("histo2DMean","histo2DMean",5000,0.,200.,1000,-50,50);
    SetStyleHistoTH2ForGraphs( histo2DMean, "#it{p}_{T,MC} (GeV/#it{c})", Form (" #mu ((#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}) (%s)",mesonLatex.Data(), mesonLatex.Data(),
                               mesonLatex.Data(),"%"),  0.04,0.05, 0.04,0.05, 1.1,1.1);
    if (mode == 0){
        histo2DMean->GetYaxis()->SetRangeUser(-1.9,1.9);
        histo2DMean->GetXaxis()->SetRangeUser(0,6);
    } else if (mode == 2 || mode == 3){
        histo2DMean->GetYaxis()->SetRangeUser(-5.4,2.4);
        if (optEnergy.CompareTo("pPb_5.023TeV") == 0)
            histo2DMean->GetYaxis()->SetRangeUser(-10.4,5.4);
        histo2DMean->GetXaxis()->SetRangeUser(0,10);
    } else if (mode == 4 || mode == 5){
        histo2DMean->GetYaxis()->SetRangeUser(-8.5,6.5);
        histo2DMean->GetXaxis()->SetRangeUser(0,10);        
    } else if (mode == 10 || mode == 11){
        histo2DMean->GetYaxis()->SetRangeUser(-33.5,29.5);
        histo2DMean->GetXaxis()->SetRangeUser(9.5,200);        
        padEPGdPtVsPt1->SetLogx();
        padEPGdPtVsPt2->SetLogx();
        histo2DMean->GetXaxis()->SetMoreLogLabels();
    }    
    histo2DMean->DrawCopy(); 

        StylingSliceHistos(Resolution_Meson_Mean,0.8);
        Resolution_Meson_Mean->SetMarkerColor(kBlue+2);
        Resolution_Meson_Mean->SetLineColor(kBlue-8);
        Resolution_Meson_Mean->Draw("same,pe");
        StylingSliceHistos(Resolution_Meson_Mean2,0.8);
        Resolution_Meson_Mean2->SetMarkerStyle(24);
        Resolution_Meson_Mean2->SetMarkerColor(kGreen+2);
        Resolution_Meson_Mean2->SetLineColor(kGreen-8);
        Resolution_Meson_Mean2->Draw("same,pe");
        
        PutProcessLabelAndEnergyOnPlot(0.95, 0.95, 0.04, collisionSystem.Data(), fTextMeasurement.Data(), detectionProcess.Data(),62, 0.03, "", 1, 1.25, 31);
        TLegend* legendMean = GetAndSetLegend2(0.45, 0.92-(0.04*2), 0.55, 0.92, 0.04, 1, "", 42, 0.35);
        legendMean->AddEntry(Resolution_Meson_Mean,"#mu Gauss");
        legendMean->AddEntry(Resolution_Meson_Mean2,"#mu Histo");
        legendMean->Draw();
        
    padEPGdPtVsPt1->Update();
    padEPGdPtVsPt2->cd();
    
    TH2F * histo2DSigma = new TH2F("histo2DSigma","histo2DSigma",5000,0.,200.,1000,-0.1,100);
    SetStyleHistoTH2ForGraphs(histo2DSigma, Form("#it{p}^{%s}_{T,MC} (GeV/#it{c})",mesonLatex.Data()), Form(" #sigma #it{p}^{%s}_{T} (%s)",mesonLatex.Data(),"%"), 0.04,0.05, 0.04,0.05, 1.1,1.1);
    if (mode == 0){
      histo2DSigma->GetYaxis()->SetRangeUser(0,7.2);
      histo2DSigma->GetXaxis()->SetRangeUser(0,6);
    } else if (mode == 2 || mode == 3){
      histo2DSigma->GetYaxis()->SetRangeUser(0,11.5);
        if (optEnergy.CompareTo("pPb_5.023TeV") == 0)
            histo2DSigma->GetYaxis()->SetRangeUser(0,20.4);
      histo2DSigma->GetXaxis()->SetRangeUser(0,10);		
    } else if (mode == 4 || mode == 5){
          histo2DSigma->GetYaxis()->SetRangeUser(0,17);
          histo2DSigma->GetXaxis()->SetRangeUser(0,10);        
    } else if (mode == 10 || mode == 11){
        histo2DSigma->GetYaxis()->SetRangeUser(0,58.5);
        histo2DSigma->GetXaxis()->SetRangeUser(9.5,200);                
        histo2DSigma->GetXaxis()->SetMoreLogLabels();
    } 
    histo2DSigma->DrawCopy(); 

        StylingSliceHistos(Resolution_Meson_Sigma,0.8);
        Resolution_Meson_Sigma->SetMarkerColor(kBlue+2);
        Resolution_Meson_Sigma->SetLineColor(kBlue-8);
        Resolution_Meson_Sigma->Draw("same,pe");
        StylingSliceHistos(Resolution_Meson_Sigma2,0.8);
        Resolution_Meson_Sigma2->SetMarkerStyle(24);
        Resolution_Meson_Sigma2->SetMarkerColor(kGreen+2);
        Resolution_Meson_Sigma2->SetLineColor(kGreen-8);
        Resolution_Meson_Sigma2->Draw("same,pe");

        TLegend* legendSigma = GetAndSetLegend2(0.78, 0.96-(0.04*2), 0.95, 0.96, 0.04, 1, "", 42, 0.15);
        legendSigma->AddEntry(Resolution_Meson_Sigma,"#sigma Gauss");
        legendSigma->AddEntry(Resolution_Meson_Sigma2,"RMS Histo");
        legendSigma->Draw();

    padEPGdPtVsPt2->Update();
    canvasEPGdPtVsPt->Update();		
    canvasEPGdPtVsPt->SaveAs(Form("%s/Resolution_%s_Pt.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()));


    if ( mode == 10 || mode == 11){
        TH1D* Resolution_Pi0Merged_PtRes[maxNBinsPt];
        TF1*  fitResolution_Pi0Merged_PtRes[maxNBinsPt];
        TH1F *Resolution_Meson_MeanMerged = new TH1F("Resolution_Meson_MeanMerged", "mean meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ;    
        TH1F *Resolution_Meson_SigmaMerged = new TH1F("Resolution_Meson_SigmaMerged", "sigma meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ; 
        minFitRange = -0.2;
        maxFitRange = 0.2;       
        
        ResolutionFittingRebined( Resolution_Pi0_MCPt_ResolPt_Merged, Resolution_Pi0Merged_PtRes, fitResolution_Pi0Merged_PtRes, maxNBinsPt, ptbinning, ptbinningReb, Resolution_Meson_MeanMerged,
                                  Resolution_Meson_SigmaMerged, "gaus", minFitRange, maxFitRange, precision, "Resolution_Pi0Merged_PtRes");

         for (Int_t i = 0; i < maxNBinsPt; i++){
            Resolution_Meson_MeanMerged->SetBinContent(i+1, Resolution_Pi0Merged_PtRes[i]->GetMean());
            Resolution_Meson_MeanMerged->SetBinError(i+1, Resolution_Pi0Merged_PtRes[i]->GetMeanError());
            Resolution_Meson_SigmaMerged->SetBinContent(i+1, Resolution_Pi0Merged_PtRes[i]->GetRMS());
            Resolution_Meson_SigmaMerged->SetBinError(i+1, Resolution_Pi0Merged_PtRes[i]->GetRMSError());
        }
 
        Resolution_Meson_SigmaMerged->Scale(100);
        Resolution_Meson_MeanMerged->Scale(100);
        
        TH1D* Resolution_Pi0PartMerged_PtRes[maxNBinsPt];
        TF1*  fitResolution_Pi0PartMerged_PtRes[maxNBinsPt];
        TH1F *Resolution_Meson_MeanPartMerged = new TH1F("Resolution_Meson_MeanPartMerged", "mean meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ;    
        TH1F *Resolution_Meson_SigmaPartMerged = new TH1F("Resolution_Meson_SigmaPartMerged", "sigma meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ; 
        
        ResolutionFittingRebined( Resolution_Pi0_MCPt_ResolPt_PartMerged, Resolution_Pi0PartMerged_PtRes, fitResolution_Pi0PartMerged_PtRes, maxNBinsPt, ptbinning, ptbinningReb, Resolution_Meson_MeanPartMerged,
                                  Resolution_Meson_SigmaPartMerged, "gaus", minFitRange, maxFitRange, precision, "Resolution_Pi0PartMerged_PtRes");

        for (Int_t i = 0; i < maxNBinsPt; i++){
            Resolution_Meson_MeanPartMerged->SetBinContent(i+1, Resolution_Pi0PartMerged_PtRes[i]->GetMean());
            Resolution_Meson_MeanPartMerged->SetBinError(i+1, Resolution_Pi0PartMerged_PtRes[i]->GetMeanError());
            Resolution_Meson_SigmaPartMerged->SetBinContent(i+1, Resolution_Pi0PartMerged_PtRes[i]->GetRMS());
            Resolution_Meson_SigmaPartMerged->SetBinError(i+1, Resolution_Pi0PartMerged_PtRes[i]->GetRMSError());
        }
 
        Resolution_Meson_SigmaPartMerged->Scale(100);
        Resolution_Meson_MeanPartMerged->Scale(100);

        TH1D* Resolution_Pi0OneGamma_PtRes[maxNBinsPt];
        TF1*  fitResolution_Pi0OneGamma_PtRes[maxNBinsPt];
        TH1F *Resolution_Meson_MeanOneGamma = new TH1F("Resolution_Meson_MeanOneGamma", "mean meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ;    
        TH1F *Resolution_Meson_SigmaOneGamma = new TH1F("Resolution_Meson_SigmaOneGamma", "sigma meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ; 
        
        ResolutionFittingRebined( Resolution_Pi0_MCPt_ResolPt_OneGamma, Resolution_Pi0OneGamma_PtRes, fitResolution_Pi0OneGamma_PtRes, maxNBinsPt, ptbinning, ptbinningReb, Resolution_Meson_MeanOneGamma,
                                  Resolution_Meson_SigmaOneGamma, "gaus", minFitRange, maxFitRange, precision, "Resolution_Pi0OneGamma_PtRes");

        for (Int_t i = 0; i < maxNBinsPt; i++){
            Resolution_Meson_MeanOneGamma->SetBinContent(i+1, Resolution_Pi0OneGamma_PtRes[i]->GetMean());
            Resolution_Meson_MeanOneGamma->SetBinError(i+1, Resolution_Pi0OneGamma_PtRes[i]->GetMeanError());
            Resolution_Meson_SigmaOneGamma->SetBinContent(i+1, Resolution_Pi0OneGamma_PtRes[i]->GetRMS());
            Resolution_Meson_SigmaOneGamma->SetBinError(i+1, Resolution_Pi0OneGamma_PtRes[i]->GetRMSError());
        }
 
        Resolution_Meson_SigmaOneGamma->Scale(100);
        Resolution_Meson_MeanOneGamma->Scale(100);

        TH1D* Resolution_Pi0OneElectron_PtRes[maxNBinsPt];
        TF1*  fitResolution_Pi0OneElectron_PtRes[maxNBinsPt];
        TH1F *Resolution_Meson_MeanOneElectron = new TH1F("Resolution_Meson_MeanOneElectron", "mean meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ;    
        TH1F *Resolution_Meson_SigmaOneElectron = new TH1F("Resolution_Meson_SigmaOneElectron", "sigma meson Resolution dPt vs Pt", maxNBinsPt,  ptbinning) ; 
        
        ResolutionFittingRebined( Resolution_Pi0_MCPt_ResolPt_OneElectron, Resolution_Pi0OneElectron_PtRes, fitResolution_Pi0OneElectron_PtRes, maxNBinsPt, ptbinning, ptbinningReb, Resolution_Meson_MeanOneElectron,
                                  Resolution_Meson_SigmaOneElectron, "gaus", minFitRange, maxFitRange, precision, "Resolution_Pi0OneElectron_PtRes");

        for (Int_t i = 0; i < maxNBinsPt; i++){
            Resolution_Meson_MeanOneElectron->SetBinContent(i+1, Resolution_Pi0OneElectron_PtRes[i]->GetMean());
            Resolution_Meson_MeanOneElectron->SetBinError(i+1, Resolution_Pi0OneElectron_PtRes[i]->GetMeanError());
            Resolution_Meson_SigmaOneElectron->SetBinContent(i+1, Resolution_Pi0OneElectron_PtRes[i]->GetRMS());
            Resolution_Meson_SigmaOneElectron->SetBinError(i+1, Resolution_Pi0OneElectron_PtRes[i]->GetRMSError());
        }
 
        Resolution_Meson_SigmaOneElectron->Scale(100);
        Resolution_Meson_MeanOneElectron->Scale(100);

        
        padEPGdPtVsPt1->cd();
        histo2DMean->GetYaxis()->SetRangeUser(-38.5,17.5);
        histo2DMean->GetXaxis()->SetRangeUser(9.5,200);        
        histo2DMean->DrawCopy(); 

            DrawGammaSetMarker(Resolution_Meson_MeanMerged, 21, 1, kAzure, kAzure);  
            Resolution_Meson_MeanMerged->Draw("same,pe");
            DrawGammaSetMarker(Resolution_Meson_MeanPartMerged, 20, 1, kRed+2, kRed+2);  
            Resolution_Meson_MeanPartMerged->Draw("same,pe");
            DrawGammaSetMarker(Resolution_Meson_MeanOneGamma, 25, 1, kGreen+2, kGreen+2);  
            Resolution_Meson_MeanOneGamma->Draw("same,pe");
            DrawGammaSetMarker(Resolution_Meson_MeanOneElectron, 24, 1, kViolet+2, kViolet+2);          
            Resolution_Meson_MeanOneElectron->Draw("same,pe");
            
            PutProcessLabelAndEnergyOnPlot(0.17, 0.95, 0.04, collisionSystem.Data(), fTextMeasurement.Data(), detectionProcess.Data());
            TLegend* legendMean2 = GetAndSetLegend2(0.17, 0.78-(0.04*4), 0.35, 0.78, 0.04, 1, "", 42, 0.15);
            legendMean2->AddEntry(Resolution_Meson_MeanMerged,"#mu Histo, pure merged");
            legendMean2->AddEntry(Resolution_Meson_MeanPartMerged,"#mu Histo, part. merged");
            legendMean2->AddEntry(Resolution_Meson_MeanOneGamma,"#mu Histo, only 1 #gamma");
            legendMean2->AddEntry(Resolution_Meson_MeanOneElectron,"#mu Histo, only 1 e^{#pm}");
            legendMean2->Draw();
            
//         padEPGdPtVsPt1->Update();
        padEPGdPtVsPt2->cd();
        histo2DSigma->GetYaxis()->SetRangeUser(0,65.5);
        histo2DSigma->GetXaxis()->SetRangeUser(9.5,200);                
        histo2DSigma->DrawCopy(); 

            DrawGammaSetMarker(Resolution_Meson_SigmaMerged, 21, 1, kAzure, kAzure);  
            Resolution_Meson_SigmaMerged->Draw("same,pe");
            DrawGammaSetMarker(Resolution_Meson_SigmaPartMerged, 20, 1, kRed+2, kRed+2);  
            Resolution_Meson_SigmaPartMerged->Draw("same,pe");
            DrawGammaSetMarker(Resolution_Meson_SigmaOneGamma, 25, 1, kGreen+2, kGreen+2);  
            Resolution_Meson_SigmaOneGamma->Draw("same,pe");
            DrawGammaSetMarker(Resolution_Meson_SigmaOneElectron, 24, 1, kViolet+2, kViolet+2);          
            Resolution_Meson_SigmaOneElectron->Draw("same,pe");

            TLegend* legendSigma2 = GetAndSetLegend2(0.17, 0.96-(0.04*4), 0.35, 0.96, 0.04, 1, "", 42, 0.15);
            legendSigma2->AddEntry(Resolution_Meson_SigmaMerged,"RMS, pure merged");
            legendSigma2->AddEntry(Resolution_Meson_SigmaPartMerged,"RMS, part. merged");
            legendSigma2->AddEntry(Resolution_Meson_SigmaOneGamma,"RMS, only 1 #gamma");
            legendSigma2->AddEntry(Resolution_Meson_SigmaOneElectron,"RMS, only 1 e^{#pm}");
            
            legendSigma2->Draw();

        padEPGdPtVsPt2->Update();
        canvasEPGdPtVsPt->Update();     
        canvasEPGdPtVsPt->SaveAs(Form("%s/Resolution_Split_%s_Pt.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()));
        
        ProduceProjectionsPlotWithFits( Resolution_Pi0Merged_PtRes ,fitResolution_Pi0Merged_PtRes, ptbinning, maxNBinsPt ,nColumns, nRows, Form(" (p^{%s}_{T,rec} -p^{%s}_{T,MC})/p^{%s}_{T,MC}", mesonLatex.Data(),
                                        mesonLatex.Data(), mesonLatex.Data()), Form ("N_{%s}", mesonLatex.Data()), "MC" , Form("fit %s",mesonLatex.Data()), "Pt",Form("%s/ProjectionsFitted_PureMerged%s.%s", 
                                        outputDirectory.Data(), mesonName.Data(), suffix.Data()), mode);

        PlotStandard2D( Resolution_Pi0_MCPt_ResolPt_Merged, Form("%s/Resolution2D_PureMerged%s.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()), "", 
                        "#it{p}_{T,MC} (GeV/#it{c})", Form ("(#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}",mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()), 
                        0, 0, 0,
                        0, -10, 10, 
                        0, 0, 50, FindSmallestEntryIn2D(Resolution_Pi0_MCPt_ResolPt_Merged),
                        0, 1, floatLocationRightDown2D, 1, 1., Form("%s fully merged", mesonLatex.Data()), 1000, 800);
        
        ProduceProjectionsPlotWithFits( Resolution_Pi0PartMerged_PtRes ,fitResolution_Pi0PartMerged_PtRes, ptbinning, maxNBinsPt ,nColumns, nRows, Form(" (p^{%s}_{T,rec} -p^{%s}_{T,MC})/p^{%s}_{T,MC}",
                                        mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()), Form ("N_{%s}", mesonLatex.Data()), "MC" , Form("fit %s",mesonLatex.Data()),
                                        "Pt",Form("%s/ProjectionsFitted_PartMerged%s.%s", outputDirectory.Data(), mesonName.Data(), suffix.Data()), mode);

        PlotStandard2D( Resolution_Pi0_MCPt_ResolPt_PartMerged, Form("%s/Resolution2D_PartMerged%s.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()), "", 
                        "#it{p}_{T,MC} (GeV/#it{c})", Form ("(#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}",mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()), 
                        0, 0, 0,
                        0, -10, 10, 
                        0, 0, 50, FindSmallestEntryIn2D(Resolution_Pi0_MCPt_ResolPt_Merged),
                        0, 1, floatLocationRightDown2D, 1, 1., Form("%s merged part conv.", mesonLatex.Data()), 1000, 800);
        PlotStandard2D( Resolution_Pi0_MCPt_ResolPt_OneGamma, Form("%s/Resolution2D_OneGamma%s.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()), "", 
                        "#it{p}_{T,MC} (GeV/#it{c})", Form ("(#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}",mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()), 
                        0, 0, 0,
                        0, -10, 10, 
                        0, 0, 50, FindSmallestEntryIn2D(Resolution_Pi0_MCPt_ResolPt_OneGamma),
                        0, 1, floatLocationRightDown2D, 1, 1., Form("%s only 1 #gamma", mesonLatex.Data()), 1000, 800);
        PlotStandard2D( Resolution_Pi0_MCPt_ResolPt_OneElectron, Form("%s/Resolution2D_OneElectron%s.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()), "", 
                        "#it{p}_{T,MC} (GeV/#it{c})", Form ("(#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}",mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()), 
                        0, 0, 0,
                        0, -10, 10, 
                        0, 0, 50, FindSmallestEntryIn2D(Resolution_Pi0_MCPt_ResolPt_OneElectron),
                        0, 1, floatLocationRightDown2D, 1, 1., Form("%s only 1 e^{#pm}", mesonLatex.Data()), 1000, 800);    
    }

    ProduceProjectionsPlotWithFits( Resolution_Pi0_PtRes ,fitResolution_Pi0_PtRes, ptbinning, maxNBinsPt ,nColumns, nRows, Form(" (p^{%s}_{T,rec} -p^{%s}_{T,MC})/p^{%s}_{T,MC}", mesonLatex.Data(),
                                    mesonLatex.Data(), mesonLatex.Data()), Form ("N_{%s}", mesonLatex.Data()), "MC" , Form("fit %s",mesonLatex.Data()), "Pt",Form("%s/ProjectionsFitted%s.%s", 
                                    outputDirectory.Data(), mesonName.Data(), suffix.Data()), mode);

    PlotStandard2D( Resolution_Pi0_MCPt_ResolPt, Form("%s/Resolution2D%s.%s",outputDirectory.Data(),mesonName.Data(),suffix.Data()), "", 
                     "#it{p}_{T,MC} (GeV/#it{c})", Form ("(#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}",mesonLatex.Data(), mesonLatex.Data(), mesonLatex.Data()), 
                    0, 0, 0,
                    0, -10, 10, 
                    0, 0, 50, FindSmallestEntryIn2D(Resolution_Pi0_MCPt_ResolPt),
                    0, 1, floatLocationRightDown2D, 1, 1., "", 1000, 800);

    delete padEPGdPtVsPt1;  
    delete padEPGdPtVsPt2;  
    delete canvasEPGdPtVsPt;

}


