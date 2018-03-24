// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

#include <stdlib.h>
#include <iostream>
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
#include "../CommonHeaders/PlottingMeson.h"
#include "TLatex.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
#include "AnalyseDCADist.h"
#include "THnSparse.h"
#include "TTree.h"

TH1F* fHistDCAZUnderMeson_MesonPt_AllCat_Data[30];
TH1F* fHistDCAZUnderMeson_MesonPt_AllCat_MC[30];

void PlotDCADistPtBinDataAndMC(TString namePlot, TString nameCanvas, TString namePad, Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins, Double_t* fRangeBinsPt, TString dateDummy);

void DrawGammaDCAHistoDataAndMC( TH1* , TH1*,
             TString , TString , TString ,
	     Float_t , Float_t , Double_t );

// Main Function
void PlotDCADataAndMC(  TString meson           ="",
                        TString cutSelection    ="",
                        TString suffix          = "",
                        TString optionEnergy    = "",
                        TString optionPeriod    = "",
                        Int_t numberOfBins      = 10,
                        Int_t mode              = 9
                    ) {

    gROOT->Reset();

    // Set default plotting styles
    StyleSettingsThesis();
    SetPlotStyle();

    // Reading out cutSelection
    TString fEventCutSelection              ="";
    TString fGammaCutSelection              ="";
    TString fClusterCutSelection            ="";
    TString fElectronCutSelection           ="";
    TString fMesonCutSelection              ="";
    fCutSelection                           = cutSelection;
    if (mode == 9){
        ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection, fMesonCutSelection);
        fEventCutSelection                  = fGammaCutSelection(0, 7);
        fGammaCutSelection                  = fGammaCutSelection(7, fGammaCutSelection.Length()-7);
        cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
    } else {
        ReturnSeparatedCutNumberAdvanced(cutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    }

    // Read out centrality string
    TString intermediate                    = GetCentralityString(fEventCutSelection);
    fTextCent                                 = "";
    if (intermediate.CompareTo("pp")==0){
        fTextCent                           = "MinBias";
        intermediate                        = "";
    } else {
        fTextCent                           = Form("%s central", intermediate.Data());
        if (intermediate.CompareTo("0-100%") == 0)
            intermediate                    = "";
    }

    // Set energy
    fEnergyText                             = ReturnFullCollisionsSystem(optionEnergy);
    if (intermediate.CompareTo("") != 0 )
        fEnergyText                         = Form("%s %s", intermediate.Data(), fEnergyText.Data());
    if (optionPeriod.CompareTo("") != 0){
        fEnergyText                         = Form("%s, %s", fEnergyText.Data(), optionPeriod.Data());
    }
    fEnergyFlag                             = optionEnergy;

    // Set DCAz cut string
    TString fDCAZCut                        = fGammaCutSelection(GetPhotonDcaZPrimVtxCutPosition(fGammaCutSelection), 1);
    fMaxDcaZPhoton                          = AnalyseDCAZPhotonCutValue(fDCAZCut.Atoi());
    if (fMaxDcaZPhoton> 10) fMaxDcaZPhoton  = 10.;

    // Set meson type
    fMesonType                              = "#pi^{0}";
    if (meson.CompareTo("Eta") == 0)
        fMesonType                          = "#eta";

    // Set date
    fdate                                   = ReturnDateString();

    // Initialize histogram arrays & different ranges for fitting, extraction ...
    InitializeBinning(meson, numberOfBins, fEnergyFlag, "", mode, fEventCutSelection, fClusterCutSelection, -1, kTRUE, intermediate, optionPeriod);

    // Set Output directory
    TString outputDir                       = Form("%s/%s/%s/%s/AnalyseDCADist", cutSelection.Data(), optionEnergy.Data(), optionPeriod.Data(), suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // Input: output file of AnalyseDCADist.C
    const char* nameOutputData = Form("%s/%s/%s_Data_GammaConvV1DCATestAnalysed%s.root",fCutSelection.Data(),optionEnergy.Data(),meson.Data(),optionPeriod.Data());
    const char* nameOutputMC = Form("%s/%s/%s_MC_GammaConvV1DCATestAnalysed%s.root",fCutSelection.Data(),optionEnergy.Data(),meson.Data(),optionPeriod.Data());
    TFile*   fOutputData = new TFile(nameOutputData);
    TFile*   fOutputMC = new TFile(nameOutputMC);


    TH1D* fEventQualityData = (TH1D*)fOutputData->Get("NEvents");
    if (optionEnergy.Contains("PbPb") || optionEnergy.Contains("pPb")){
      fNEvents = fEventQualityData->GetBinContent(1);
    } else {
      fNEvents = GetNEvents(fEventQualityData);
    }

    // Set expected mass
    fMesonMassExpect                        = TDatabasePDG::Instance()->GetParticle(fMesonId)->Mass();

 
    //#############################
    //### GET HISTOS FROM FILES ###
    //#############################


    for (Int_t j = fStartPtBin; j < fNBinsPt ; j++){

      fHistDCAZUnderMeson_MesonPt_AllCat_Data[j] = (TH1F*) fOutputData->Get(Form("HistDCAZUnderMesonAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]));
      fHistDCAZUnderMeson_MesonPt_AllCat_MC[j] = (TH1F*) fOutputMC->Get(Form("HistDCAZUnderMesonAllCat_MesonPt_%3.2f-%3.2f", fBinsPt[j], fBinsPt[j+1]));
     
    }


    //############################
    //######### PLOTS ############
    //############################
    
    PlotDCADistPtBinDataAndMC(Form("%s/%s_DCAProjectionsDataAndMC_AllCat.%s", outputDir.Data(), meson.Data(), suffix.Data()),
    "canvas_Allcat", "pad_Allcat",   fRow, fColumn, fStartPtBin, fNBinsPt, fBinsPt, fdate);

    //PlotDCADistPtBinDataAndMCCat


}

//############################
//###### FUNCTIONS ###########
//############################

//______________________________ DCAz photon under Meson Peak together with Fit and Histo estimate of BG __________________________________
void PlotDCADistPtBinDataAndMC(TString namePlot, TString nameCanvas, TString namePad,
                                        Int_t fRowPlot, Int_t fColumnPlot, Int_t fStartBinPtRange, Int_t fNumberPtBins,
                                        Double_t* fRangeBinsPt, TString dateDummy){
    TCanvas *canvas          = new TCanvas(nameCanvas.Data(),"",2800,1800);  // gives the page size
    canvas->SetTopMargin(0.0);
    canvas->SetBottomMargin(0.0);
    canvas->SetRightMargin(0.0);
    canvas->SetLeftMargin(0.0);

    TPad * padDataFit                 = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas
    padDataFit->SetFillColor(0);
    padDataFit->GetFrame()->SetFillColor(0);
    padDataFit->SetBorderMode(0);
    padDataFit->Divide(fColumnPlot,fRowPlot,0.0,0.0);
// 	padDataFit->SetLeftMargin(0.2);
    padDataFit->Draw();

    Int_t place = 0;
    for(Int_t iPt = fStartBinPtRange; iPt < fNumberPtBins; iPt++){
        Double_t startPt 			= fRangeBinsPt[iPt];
        Double_t endPt                 = fRangeBinsPt[iPt+1];

        place = place + 1;                  //give the right place in the page
        if (place == fColumnPlot){
            iPt--;
            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.001);


            string textAlice 			= "ALICE performance";
            string textEvents			= "Data";
            Double_t textHeight 	= 0.055;
            Double_t startTextX 	= 0.05;
            Double_t startTextY 	= 0.75;
            Double_t differenceText = textHeight*1.25;

            TLatex *alice 			= new TLatex(startTextX, (startTextY+(2*differenceText)), Form("%s", textAlice.c_str()));
            TLatex *energy 			= new TLatex(startTextX, (startTextY-differenceText), Form("%s", fEnergyText.Data()));
            TLatex *events 			= new TLatex(startTextX, startTextY, Form("%s: %2.1e events", textEvents.c_str(), fNEvents));
            TLatex *latexDate 		= new TLatex(startTextX, startTextY+differenceText, dateDummy.Data());
            TLatex *latexCategory 	= new TLatex(startTextX, startTextY-(3*differenceText), Form("All Meson Categories"));

            alice->SetNDC();
            alice->SetTextColor(1);
            alice->SetTextSize(textHeight);
            alice->Draw();

            energy->SetNDC();
            energy->SetTextColor(1);
            energy->SetTextSize(textHeight);
            energy->Draw();

            events->SetNDC();
            events->SetTextColor(1);
            events->SetTextSize(textHeight);
            events->Draw();

            latexDate->SetNDC();
            latexDate->SetTextColor(1);
            latexDate->SetTextSize(textHeight);
            latexDate->Draw();

            latexCategory->SetNDC();
            latexCategory->SetTextColor(1);
            latexCategory->SetTextSize(textHeight*1.5);
            latexCategory->Draw();

            TLegend* legend 		= new TLegend(0.05, 0.1, 1., 0.5);
            legend->SetTextSize(textHeight);
            legend->SetFillColor(0);
            legend->SetLineColor(0);
            legend->SetNColumns(1);
            legend->SetMargin(0.1);
            if (fHistDCAZUnderMeson_MesonPt_AllCat_Data[fStartBinPtRange])
                legend->AddEntry(fHistDCAZUnderMeson_MesonPt_AllCat_Data[fStartBinPtRange], "Data", "p");
            if (fHistDCAZUnderMeson_MesonPt_AllCat_MC[fStartBinPtRange])
	        fHistDCAZUnderMeson_MesonPt_AllCat_MC[fStartBinPtRange]->SetLineColor(kRed);
                legend->AddEntry(fHistDCAZUnderMeson_MesonPt_AllCat_MC[fStartBinPtRange], "MC scaled to data maximum", "l");
            legend->Draw();

        } else {
            padDataFit->cd(place);
            padDataFit->cd(place)->SetTopMargin(0.12);
            padDataFit->cd(place)->SetBottomMargin(0.15);
            padDataFit->cd(place)->SetRightMargin(0.001);

            padDataFit->cd(place)->SetLogy();

            int remaining = (place-1)%fColumnPlot;
            if (remaining > 0) padDataFit->cd(place)->SetLeftMargin(0.15);
            else padDataFit->cd(place)->SetLeftMargin(0.25);

            if (fEnergyFlag.CompareTo("PbPb_2.76TeV")==0 || fEnergyFlag.CompareTo("pPb_5.023TeV")==0){
                if (fHistDCAZUnderMeson_MesonPt_AllCat_Data[iPt] && fHistDCAZUnderMeson_MesonPt_AllCat_MC[iPt]){
		  DrawGammaDCAHistoDataAndMC( fHistDCAZUnderMeson_MesonPt_AllCat_Data[iPt], fHistDCAZUnderMeson_MesonPt_AllCat_MC[iPt], 
                                    Form("%3.2f GeV/c < #it{p}_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
                                    "dca_{z} #gamma (cm)", Form("d(dca_{z}), all cat"),
                                    -10, 10, 0.5e5);
                }
            } else {
                if (fHistDCAZUnderMeson_MesonPt_AllCat_Data[iPt] && fHistDCAZUnderMeson_MesonPt_AllCat_MC[iPt]){
		  DrawGammaDCAHistoDataAndMC( fHistDCAZUnderMeson_MesonPt_AllCat_Data[iPt], fHistDCAZUnderMeson_MesonPt_AllCat_MC[iPt], 
                                    Form("%3.2f GeV/c < #it{p}_{T,%s} < %3.2f GeV/c", startPt, fMesonType.Data(), endPt),
                                    "dca_{z} #gamma (cm)", Form("d(dca_{z}), all cat"),
                                    -10, 10, 0.5e4);
                }
            }
        }
    }
    canvas->Print(namePlot.Data());
    delete padDataFit;
    delete canvas;
}


void DrawGammaDCAHistoDataAndMC( TH1* histo1, TH1* histo2,
                        TString Title, TString XTitle, TString YTitle,
				 Float_t xMin, Float_t xMax, Double_t numberOfOrders ) {

    histo1->GetXaxis()->SetRangeUser(xMin, xMax);

    Double_t yMin = 0;
    Double_t yMax = 0;
    for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
      if (histo1->GetBinContent(i) > yMax){
	yMax = histo1->GetBinContent(i);
      }
    }
    yMin = yMax/numberOfOrders;
    if (yMin < 1e-1) yMin = 1e-1;
    histo1->GetYaxis()->SetRangeUser(yMin, 2*yMax);
    
    
    if(XTitle.Length() > 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.Length() > 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetYaxis()->SetLabelSize(0.02);
    histo1->GetYaxis()->SetTitleSize(0.025);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetYaxis()->SetTitleOffset(0.5);
    histo1->GetXaxis()->SetTitleSize(0.025);
    histo1->GetXaxis()->SetLabelSize(0.02);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerColor(1);
    histo1->SetLineColor(1);
    histo1->SetLineWidth(0.5);
    histo1->SetMarkerSize(0.2);
    histo1->SetTitleOffset(1.2,"xy");
    histo1->SetTitleSize(0.05,"xy");
    histo1->GetYaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetLabelSize(0.05);
    histo1->GetXaxis()->SetNdivisions(507,kTRUE);


    TH1* histo2copy = (TH1*)histo2->Clone("histo2copy");
    histo2copy->Scale( histo1->Integral()  / histo2->Integral() );


    histo1->SetTitle("");
    histo1->DrawCopy("e1,p");
    histo2copy->SetLineColor(kRed);
    histo2copy->DrawCopy("histo, same");
    if(Title.Length() > 0){
      histo1->SetTitle("");
      TLatex *alice = new TLatex(0.1,0.95,Form("%s",Title.Data())); 
      alice->SetNDC();
      alice->SetTextColor(1);
      alice->SetTextSize(0.062);
      alice->Draw();
    }
    
}
