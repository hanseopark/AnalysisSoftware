/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                                  *****
 ******     Friederike Bock, friederike.bock@cern.ch                                                                    *****
 *****************************************************************************************************************************/

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


//*************************************************************************************************
//*********************** Main function ***********************************************************
//*************************************************************************************************
void MergeCorrFactorsJJandJJGammaTrigMergedCluster( TString fCutSelection, 
                                                    TString mesonType, 
                                                    TString fSuffix, 
                                                    TString fEnergyFlag, 
                                                    TString nameFileCorrectionFileFull, 
                                                    TString nameFileJJ, 
                                                    TString nameFileJJGammaTrigg
                                                  ){
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    gSystem->Exec("cp "+nameFileCorrectionFileFull+" "+nameFileJJ);
    TString outputDir                       = Form("%s/%s/%s/ComparisonMCs", fCutSelection.Data(), fEnergyFlag.Data(), fSuffix.Data());
    gSystem->Exec("mkdir "+outputDir);
   
    //*************************************************************************************************
    //****************************** read ***********************************************
    //*************************************************************************************************
    // read new output file
    TFile* file                             = new TFile (nameFileCorrectionFileFull);
    TH1D *histoTrueEffiPt                   = (TH1D*)file->Get("TrueEfficiencyMergedPt"); 
    TH1D *histoTrueEffiPrimPt               = (TH1D*)file->Get("TrueEfficiencyPrimMesonPt"); 
    TH1D *histoAcceptPt                     = (TH1D*)file->Get("AcceptancePt"); 
    TH1D *histoInputRebined                 = (TH1D*)file->Get("MC_Meson_WithinAccept_Rebin"); 
    Double_t maxPt                          = histoTrueEffiPt->GetXaxis()->GetBinUpEdge(histoTrueEffiPt->GetNbinsX());
    
    // read pure JJ MC file
    TFile* fileJJ                           = new TFile (nameFileJJ);
    TH1F *histoEventQualityJJ               = (TH1F*)fileJJ->Get("NEvents");
    TH1D *histoTrueEffiPtJJ                 = (TH1D*)fileJJ->Get("TrueEfficiencyMergedPt"); 
    TH1D *histoTrueEffiPrimPtJJ             = (TH1D*)fileJJ->Get("TrueEfficiencyPrimMesonPt"); 
    TH1D *histoAcceptPtJJ                   = (TH1D*)fileJJ->Get("AcceptancePt"); 
    TH1D *histoInputRebinedJJ               = (TH1D*)fileJJ->Get("MC_Meson_WithinAccept_Rebin"); 
    Double_t nEvtMCJJ                       = GetNEvents(histoEventQualityJJ);
    histoInputRebinedJJ->Sumw2();
    histoInputRebinedJJ->Scale(1./nEvtMCJJ);
    
    // read gamma triggered JJ MC file
    TFile* fileJJGammaTrigg                 = new TFile (nameFileJJGammaTrigg);
    TH1F *histoEventQualityJJGammaTrigg     = (TH1F*)fileJJGammaTrigg->Get("NEvents");
    TH1D *histoTrueEffiPtJJGammaTrigg       = (TH1D*)fileJJGammaTrigg->Get("TrueEfficiencyMergedPt"); 
    TH1D *histoTrueEffiPrimPtJJGammaTrigg   = (TH1D*)fileJJGammaTrigg->Get("TrueEfficiencyPrimMesonPt"); 
    TH1D *histoAcceptPtJJGammaTrigg         = (TH1D*)fileJJGammaTrigg->Get("AcceptancePt"); 
    TH1D *histoInputRebinedJJGammaTrigg     = (TH1D*)fileJJGammaTrigg->Get("MC_Meson_WithinAccept_Rebin"); 
    Double_t nEvtMCJJGammaTrigg             = GetNEvents(histoEventQualityJJGammaTrigg);
    histoInputRebinedJJGammaTrigg->Sumw2();
    histoInputRebinedJJGammaTrigg->Scale(1./nEvtMCJJGammaTrigg);
    //*************************************************************************************************
    //************************* calculate weights and merged histos ***********************************
    //*************************************************************************************************    
    for (Int_t i = 1; i < histoTrueEffiPtJJ->GetNbinsX()+1 ; i++){
        Double_t relErrJJ      = histoTrueEffiPtJJ->GetBinError(i)/histoTrueEffiPtJJ->GetBinContent(i)*100;
        Double_t relErrJJGammaTrigg  = histoTrueEffiPtJJGammaTrigg->GetBinError(i)/histoTrueEffiPtJJGammaTrigg->GetBinContent(i)*100;
                
        Double_t weightJJ           = 1/TMath::Power(histoTrueEffiPtJJ->GetBinError(i),2);
        Double_t weightJJGammaTrigg = 1/TMath::Power(histoTrueEffiPtJJGammaTrigg->GetBinError(i),2);
        Double_t weightSum          = weightJJ + weightJJGammaTrigg;
        Double_t weightedEffi       = (weightJJ*histoTrueEffiPtJJ->GetBinContent(i) + weightJJGammaTrigg * histoTrueEffiPtJJGammaTrigg->GetBinContent(i))/weightSum;
        Double_t weightedEffiErr    = pow((weightJJ +  weightJJGammaTrigg),-0.5);
        
        if (isfinite(weightedEffi) && isfinite(weightedEffiErr)){
            histoTrueEffiPt->SetBinContent(i, weightedEffi);
            histoTrueEffiPt->SetBinError(i, weightedEffiErr);
        }

        Double_t relErrPrimJJ       = histoTrueEffiPrimPtJJ->GetBinError(i)/histoTrueEffiPrimPtJJ->GetBinContent(i)*100;
        Double_t relErrPrimJJGammaTrigg  = histoTrueEffiPrimPtJJGammaTrigg->GetBinError(i)/histoTrueEffiPrimPtJJGammaTrigg->GetBinContent(i)*100;
                
        Double_t weightPrimJJ           = 1/TMath::Power(histoTrueEffiPrimPtJJ->GetBinError(i),2);
        Double_t weightPrimJJGammaTrigg = 1/TMath::Power(histoTrueEffiPrimPtJJGammaTrigg->GetBinError(i),2);
        Double_t weightPrimSum          = weightPrimJJ + weightPrimJJGammaTrigg;
        Double_t weightedEffiPrim       = (weightPrimJJ*histoTrueEffiPrimPtJJ->GetBinContent(i) + weightPrimJJGammaTrigg * histoTrueEffiPrimPtJJGammaTrigg->GetBinContent(i))/weightPrimSum;
        Double_t weightedEffiPrimErr    = pow((weightPrimJJ +  weightPrimJJGammaTrigg),-0.5);
        
        if (isfinite(weightedEffiPrim) && isfinite(weightedEffiPrimErr)){
            histoTrueEffiPrimPt->SetBinContent(i, weightedEffiPrim);
            histoTrueEffiPrimPt->SetBinError(i, weightedEffiPrimErr);
        }
        
        Double_t relErrAJJ          = histoAcceptPtJJ->GetBinError(i)/histoAcceptPtJJ->GetBinContent(i)*100;
        Double_t relErrAJJGammaTrigg= histoAcceptPtJJGammaTrigg->GetBinError(i)/histoAcceptPtJJGammaTrigg->GetBinContent(i)*100;
                
        Double_t weightAJJ          = 1/TMath::Power(histoAcceptPtJJ->GetBinError(i),2);
        Double_t weightAJJGammaTrigg= 1/TMath::Power(histoAcceptPtJJGammaTrigg->GetBinError(i),2);
        Double_t weightASum         = weightJJ + weightJJGammaTrigg;
        Double_t weightedAcc        = (weightJJ*histoAcceptPtJJ->GetBinContent(i) + weightJJGammaTrigg * histoAcceptPtJJGammaTrigg->GetBinContent(i))/weightSum;
        Double_t weightedAccErr     = pow((weightJJ +  weightJJGammaTrigg),-0.5);
        
        if (isfinite(weightedAcc) && isfinite(weightedAccErr)){
            histoAcceptPt->SetBinContent(i, weightedAcc);
            histoAcceptPt->SetBinError(i, weightedAccErr);
        }
    }
    
    //*************************************************************************************************
    //****************************** calculate ratios *************************************************
    //*************************************************************************************************
    TH1D* ratioJJFinal              = (TH1D*) histoTrueEffiPtJJ->Clone("ratioJJFinal");
    ratioJJFinal->Divide(histoTrueEffiPtJJ,histoTrueEffiPt , 1.,1.,"B");
    TH1D* ratioJJGammaTriggFinal    = (TH1D*) histoTrueEffiPtJJGammaTrigg->Clone("ratioJJGammaTriggFinal");
    ratioJJGammaTriggFinal->Divide(histoTrueEffiPtJJGammaTrigg,histoTrueEffiPt , 1.,1.,"B");
    TH1D* ratioJJJJGammaTrigg       = (TH1D*) histoTrueEffiPtJJ->Clone("ratioJJGammaTriggFinal");
    ratioJJJJGammaTrigg->Divide(histoTrueEffiPtJJ,histoTrueEffiPtJJGammaTrigg , 1.,1.,"B");

    TH1D* ratioAccJJFinal           = (TH1D*) histoAcceptPtJJ->Clone("ratioAccJJFinal");
    ratioAccJJFinal->Divide(histoAcceptPtJJ,histoAcceptPt , 1.,1.,"B");
    TH1D* ratioAccJJGammaTriggFinal = (TH1D*) histoAcceptPtJJGammaTrigg->Clone("ratioAccJJGammaTriggFinal");
    ratioAccJJGammaTriggFinal->Divide(histoAcceptPtJJGammaTrigg,histoAcceptPt , 1.,1.,"B");
    TH1D* ratioAccJJJJGammaTrigg    = (TH1D*) histoAcceptPtJJ->Clone("ratioAccJJJJGammaTrigg");
    ratioAccJJJJGammaTrigg->Divide(histoAcceptPtJJ,histoAcceptPtJJGammaTrigg , 1.,1.,"B");
    
    TH1D* ratioInputs           = (TH1D*) histoInputRebinedJJ->Clone("ratioInputs");
    ratioInputs->Divide(histoInputRebinedJJ,histoInputRebinedJJGammaTrigg , 1.,1.,"");

    //***************************************************************************************************************
    //************************************Plotting Meson input yield *****************************************
    //***************************************************************************************************************
    TCanvas* canvasInputSpectra = new TCanvas("canvasInputSpectra","",0,0,1000,1350);// gives the page size
    DrawGammaCanvasSettings( canvasInputSpectra, 0.15, 0.015, 0.015, 0.07);
    canvasInputSpectra->SetLogy();
//     canvasInputSpectra->SetGridx();
    Double_t minCorrYield       = 2e-12;
    Double_t maxCorrYield       = 1e-3;
    
    DrawAutoGammaMesonHistos( histoInputRebinedJJ,
                    "", "#it{p}_{T} (GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}",
                    kTRUE, 5., 10e-10, kTRUE,
                    kFALSE, 0.5, 1.5,
                    kFALSE, 0., 7.9);
    histoInputRebinedJJ->GetYaxis()->SetTitleOffset(1.5);
    DrawGammaSetMarker(histoInputRebinedJJ, 25, 1., kBlue, kBlue);                                        
    histoInputRebinedJJ->Draw("p,E1");
    DrawGammaSetMarker(histoInputRebinedJJGammaTrigg, 26, 1., kRed, kRed);                                      
    histoInputRebinedJJGammaTrigg->Draw("p,E1,same");

    TLegend* legendInput = GetAndSetLegend2(0.6,0.85,0.8,0.95, 28);
    legendInput->SetMargin(0.2);
    legendInput->AddEntry(histoInputRebinedJJ,"Jet-Jet");
    legendInput->AddEntry(histoInputRebinedJJGammaTrigg,"Jet-Jet #gamma triggered");
    legendInput->Draw();
    
    canvasInputSpectra->Update();
    canvasInputSpectra->SaveAs(Form("%s/%s_InputYield.%s",outputDir.Data(),mesonType.Data(),fSuffix.Data()));

    //*************************************************************************************************
    //*************************** ratio inputs**********************************************************
    //*************************************************************************************************
    
    TCanvas* canvasInputRatio = new TCanvas("canvasInputRatio","",200,10,1350,900);
    DrawGammaCanvasSettings( canvasInputRatio, 0.06, 0.015, 0.015, 0.08);
    
    DrawGammaSetMarker(ratioInputs, 24, 1., kBlue+2, kBlue+2);
    DrawAutoGammaMesonHistos( ratioInputs,
                    "", "#it{p}_{T} (GeV/#it{c})", "JJ/ JJ gamma trigg",
                    kFALSE, 5., 10e-10, kTRUE,
                    kTRUE, 0.2, 1.,
                    kFALSE, 0., 7.9);
    ratioInputs->GetYaxis()->SetTitleOffset(0.7);
    ratioInputs->Draw("e1");
    canvasInputRatio->Update();
    DrawGammaLines(0., 50.,0.8, 0.8,0.1);
    DrawGammaLines(0., 50.,0.6, 0.6,0.1);
    
    canvasInputRatio->SaveAs(Form("%s/%s_MC_RatioInput_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    
    //*************************************************************************************************
    //*************************** Ratio effi weighted *************************************************
    //*************************************************************************************************
    TCanvas* canvasFraction2 = new TCanvas("canvasFraction2","",200,10,1350,900);
    DrawGammaCanvasSettings( canvasFraction2, 0.06, 0.015, 0.015, 0.08);
    
    DrawGammaSetMarker(ratioJJFinal, 24, 1., kBlue+2, kBlue+2);
    DrawAutoGammaMesonHistos( ratioJJFinal,
                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff,A}/ #epsilon_{eff,B}",
                    kFALSE, 5., 10e-10, kTRUE,
                    kTRUE, 0.5, 1.5,
                    kFALSE, 0., 7.9);
    ratioJJFinal->GetYaxis()->SetTitleOffset(0.7);
    ratioJJFinal->Draw("e1");
    DrawGammaSetMarker(ratioJJGammaTriggFinal, 25, 1., kRed, kRed);
    ratioJJGammaTriggFinal->Draw("same");
//     DrawGammaSetMarker(ratioJJJJGammaTrigg, 26, 1., kGreen+2, kGreen+2);
//     ratioJJJJGammaTrigg->Draw("same");
    canvasFraction2->Update();
    DrawGammaLines(0., maxPt,1., 1.,0.1);

    TLegend* legendRatioWeightPP = GetAndSetLegend2(0.15,0.85,0.4,0.95, 28);
    legendRatioWeightPP->AddEntry(ratioJJGammaTriggFinal,"Jet-Jet gamma trigg / final","p");
//     legendRatioWeightPP->AddEntry(ratioJJJJGammaTrigg,"Min Bias/ added Signal","p");
    legendRatioWeightPP->AddEntry(ratioJJFinal, "Jet-Jet / final","p");
    legendRatioWeightPP->Draw();
    
    canvasFraction2->SaveAs(Form("%s/%s_MC_RatioComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    //*************************************************************************************************
    //*************************** Drawing different Efficiencies Log **********************************
    //*************************************************************************************************
    TCanvas* canvasEffSimple = new TCanvas("canvasEffSimple","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasEffSimple, 0.065, 0.01, 0.01, 0.08);
    canvasEffSimple->SetLogy(1);

    DrawAutoGammaMesonHistos( histoTrueEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kTRUE, 2., 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., maxPt);
    DrawGammaSetMarker(histoTrueEffiPt, 24, 1., kBlack, kBlack);                                         
    histoTrueEffiPt->GetYaxis()->SetTitleOffset(0.7);
    histoTrueEffiPt->DrawCopy("e1");    
    DrawGammaSetMarker(histoTrueEffiPtJJ, 25, 1., kBlue, kBlue);                                        
    histoTrueEffiPtJJ->DrawCopy("e1,same");    
    DrawGammaSetMarker(histoTrueEffiPtJJGammaTrigg, 26, 1., kRed, kRed);                                      
    histoTrueEffiPtJJGammaTrigg->DrawCopy("e1,same");        

    TLegend* legendEff = GetAndSetLegend2(0.4,0.125,0.6,0.205, 28);
    legendEff->SetMargin(0.2);
    legendEff->AddEntry(histoTrueEffiPt,"combined ");
    legendEff->AddEntry(histoTrueEffiPtJJ,"Jet-Jet");
    legendEff->AddEntry(histoTrueEffiPtJJGammaTrigg,"Jet-Jet #gamma triggered");
    legendEff->Draw();
        
    canvasEffSimple->Update();
    canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonEffiMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    //*************************************************************************************************
    //*************************** Drawing different Efficiencies Lin **********************************
    //*************************************************************************************************
    canvasEffSimple->SetLogy(0);    
    DrawAutoGammaMesonHistos( histoTrueEffiPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "#epsilon_{eff}", 
                                    kTRUE, 0.75, 3e-6, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., maxPt);
            
    histoTrueEffiPt->DrawCopy("e1");    
    histoTrueEffiPtJJ->DrawCopy("e1,same");    
    histoTrueEffiPtJJGammaTrigg->DrawCopy("e1,same");            
    legendEff->Draw();
        
    canvasEffSimple->Update();
    canvasEffSimple->SaveAs(Form("%s/%s_MC_ComparisonEffiMerged_linY_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    delete canvasEffSimple;

    //*************************************************************************************************
    //*************************** Ratio effi weighted *************************************************
    //*************************************************************************************************
    canvasFraction2->cd();
    DrawGammaSetMarker(ratioAccJJFinal, 24, 1., kBlue+2, kBlue+2);
    DrawAutoGammaMesonHistos( ratioAccJJFinal,
                    "", "#it{p}_{T} (GeV/#it{c})", "A_{A}/ A_{B}",
                    kFALSE, 5., 10e-10, kTRUE,
                    kTRUE, 0.5, 1.5,
                    kFALSE, 0., maxPt);
    ratioAccJJFinal->GetYaxis()->SetTitleOffset(0.7);
    ratioAccJJFinal->Draw("e1");
    DrawGammaSetMarker(ratioAccJJGammaTriggFinal, 25, 1., kRed, kRed);
    ratioAccJJGammaTriggFinal->Draw("same");
//     DrawGammaSetMarker(ratioAccJJJJGammaTrigg, 26, 1., kGreen+2, kGreen+2);
//     ratioAccJJJJGammaTrigg->Draw("same");
    canvasFraction2->Update();
    DrawGammaLines(0., maxPt,1., 1.,0.1);

    legendRatioWeightPP->Draw();
    
    canvasFraction2->SaveAs(Form("%s/%s_MC_RatioComparisonAcceptanceMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));

    //*************************************************************************************************
    //*************************** Drawing different Acceptance Lin **********************************
    //*************************************************************************************************
    TCanvas* canvasAccSimple = new TCanvas("canvasAccSimple","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings( canvasAccSimple, 0.065, 0.01, 0.01, 0.08);
    canvasAccSimple->SetLogy(0);

    DrawAutoGammaMesonHistos( histoAcceptPt, 
                                    "", "#it{p}_{T} (GeV/#it{c})", "A", 
                                    kTRUE, 1.25, 0, kFALSE,
                                    kFALSE, 0., 0.7, 
                                    kFALSE, 0., maxPt);
    DrawGammaSetMarker(histoAcceptPt, 24, 1., kBlack, kBlack);                                         
    histoAcceptPt->GetYaxis()->SetTitleOffset(0.7);
    histoAcceptPt->DrawCopy("e1");    
    DrawGammaSetMarker(histoAcceptPtJJ, 25, 1., kBlue, kBlue);                                        
    histoAcceptPtJJ->DrawCopy("e1,same  ");    
    DrawGammaSetMarker(histoAcceptPtJJGammaTrigg, 26, 1., kRed, kRed);                                      
    histoAcceptPtJJGammaTrigg->DrawCopy("e1,same");        

    TLegend* legendAcc = GetAndSetLegend2(0.4,0.125,0.6,0.205, 28);
    legendAcc->SetMargin(0.2);
    legendAcc->AddEntry(histoAcceptPt,"combined ");
    legendAcc->AddEntry(histoAcceptPtJJ,"Jet-Jet");
    legendAcc->AddEntry(histoAcceptPtJJGammaTrigg,"Jet-Jet #gamma triggered");
    legendAcc->Draw();
        
    canvasAccSimple->Update();
    canvasAccSimple->SaveAs(Form("%s/%s_MC_ComparisonAcceptanceMerged_%s.%s",outputDir.Data(),mesonType.Data(),fCutSelection.Data(),fSuffix.Data()));
    
    //*************************************************************************************************
    //**************************** write final merged file ********************************************
    //*************************************************************************************************
    TFile* fileSum = new TFile (nameFileCorrectionFileFull,"UPDATE");
        histoTrueEffiPt->Write("TrueMesonEffiPt",TObject::kOverwrite);
        histoTrueEffiPrimPt->Write("TrueEfficiencyPrimMesonPt",TObject::kOverwrite);
//         histoAcceptPt->Write("AcceptancePt",TObject::kOverwrite);
        histoEventQualityJJGammaTrigg->Write("NEvents_JJGammaTrigg",TObject::kOverwrite);
        fileSum->Write();
    fileSum->Close();

}
