//  **********************************************************************************
//  ******     provided by Gamma Conversion Group, PWGGA,                        *****
//  ******     Lucas Altenkämper, lucas.altenkamper@cern.ch                      *****
//  **********************************************************************************

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
#include "TTree.h"
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ExtractSignalBinning.h"
#include "CommonHeaders/ConversionFunctions.h"

// some definitions ******************************************************************
TString     outputDir                                                       = "";
TString     fDetectionProcess                                               = "";
TString     energyString                                                    = "";
TString     collisionSystem                                                 = "";
TString     centralityString                                                = "";
Double_t    eventNormScalingFactor                                          = 1.;
const Int_t nParticles                                                      = 11;
TString     particleName[nParticles]                                        = {"Pi0", "Eta", "EtaPrime", "Omega", "Pi+-", "K+-", "Phi", "Rho0", "Rho+-", "K0s", "JPsi"};
TString     particleNameCocktailInput[nParticles]                           = {"NPion", "Eta", "EtaPrime", "Omega", "CPion", "CKaon", "Phi", "NRho", "CRho", "NKaonSubS", "JPsi"};
TString     particleNameLatex[nParticles]                                   = {"#pi^{0}", "#eta", "#eta'", "#omega", "(#pi^{+}+#pi^{-})/2", "(K^{+}+K^{-})/2", "#phi", "#rho^{0}", "(#rho^{+}+#rho^{-})/2", "K^{0}_{s}", "J/#psi"};
Int_t       particlePDGCode[nParticles]                                     = {111, 221, 331, 223, 211, 321, 333, 113, 213, 310, 443};
Double_t    particleMtScalingFactor[nParticles]                             = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};     // for scaling from pi0
Double_t    particleMtScalingFactorErr[nParticles]                          = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};     // for scaling from pi0
Double_t    particleMtScalingFactorResFeedDownCorr[nParticles]              = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};     // for scaling from pi0 (res. feed down corrected)
Double_t    particleMtScalingFactorErrResFeedDownCorr[nParticles]           = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};     // for scaling from pi0 (res. feed down corrected)
Double_t    particleMtScalingFactorPhi[nParticles]                          = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};     // for scaling from phi
Double_t    particleMtScalingFactorErrPhi[nParticles]                       = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};     // for scaling from phi

// declarations **********************************************************************
TH1D*   pi0InvYield                                                         = NULL;
TH1D*   pi0Yield                                                            = NULL;
TH1D*   pi0InvYieldWOResFeedDown                                            = NULL;
TH1D*   pi0YieldWOResFeedDown                                               = NULL;
TH1D*   pi0EtaBinningInvYield                                               = NULL;
TH1D*   pi0EtaBinningYield                                                  = NULL;
TH1D*   pi0EtaBinningInvYieldWOResFeedDown                                  = NULL;
TH1D*   pi0EtaBinningYieldWOResFeedDown                                     = NULL;
TH1D*   etaInvYield                                                         = NULL;
TH1D*   etaYield                                                            = NULL;
TH1D*   etaToPi0Ratio                                                       = NULL;
TH1D*   etaToPi0RatioWOResFeedDown                                          = NULL;

TF1*    constantFitPhiToPi0Ratio                                            = NULL;
TF1*    constantFitEtaToPi0WOResFeedDownRatio                               = NULL;
TF1*    constantFitParticleRatio[nParticles*nParticles]                     = {NULL};

TH1D*   particleYield[nParticles]                                           = {NULL};
TH1D*   particleRatio[nParticles*nParticles]                                = {NULL};

TH1D*   mtScalingFactorCocktail                                             = NULL;
TH1D*   mtScalingFactorRecalc                                               = NULL;
TH1D*   mtScalingFactorResFeedDownCorr                                      = NULL;
TH1D*   mtScalingFactorPhi                                                  = NULL;

TF1*    yieldParametrizations[nParticles]                                   = {NULL};
TF1*    yieldParametrizationPi0WOResFeedDown                                = NULL;
TF1*    yieldParametrizationsMtScaled[nParticles]                           = {NULL};
TF1*    yieldParametrizationsMtScaledWOResFeedDown[nParticles]              = {NULL};
TF1*    yieldParametrizationsMtScaledPhi[nParticles]                        = {NULL};

TH1D*   yieldParametrizationsRatioToData[nParticles]                        = {NULL};
TH1D*   yieldParametrizationPi0WOResFeedDownRatioToData                     = NULL;
TH1D*   yieldParametrizationsMtScaledRatioToData[nParticles]                = {NULL};
TH1D*   yieldParametrizationsMtScaledWOResFeedDownRatioToData[nParticles]   = {NULL};
TH1D*   yieldParametrizationsMtScaledPhiRatioToData[nParticles]             = {NULL};

// get iterator for particle ratio ***************************************************
Int_t GetParticleRatioIter(Int_t numPDG, Int_t denomPDG) {
    
    Int_t   ratioIter                                               = 0;
    Bool_t  stop                                                    = kFALSE;
    
    for (Int_t i = 0; i < nParticles; i++) {
        if (stop) break;
        for (Int_t j = 0; j < nParticles; j++) {
            if (particlePDGCode[j] == numPDG && particlePDGCode[i] == denomPDG)
                stop                                                = kTRUE;
            if (stop) break;
            ratioIter++;
        }
    }
    
    return ratioIter;
}

// fct. to transform inv. yield to yield *********************************************
TH1D* TransformInvYield(TH1D* inputHisto, TString outputHistoName) {
    
    if (!inputHisto) return NULL;
    
    TH1D* outputHisto                                               = (TH1D*)inputHisto->Clone(outputHistoName);
    outputHisto->Sumw2();
    
    // multiply with 2*pi
    outputHisto->Scale(2*TMath::Pi());
    
    // multiply with pt
    for (Int_t i = 1; i <= outputHisto->GetNbinsX(); i++) {
        outputHisto->SetBinContent(i, outputHisto->GetBinContent(i)   * outputHisto->GetBinCenter(i));
        outputHisto->SetBinError(i,   outputHisto->GetBinError(i)     * outputHisto->GetBinCenter(i));
    }

    return outputHisto;
}

// fct. to transform graph to hist ***************************************************
TH1D* TransformToTH1D(TObject* inputObject, TString outputHistoName) {
    
    if (!inputObject) return NULL;
    
    TString inputObjectClassName                                    = inputObject->ClassName();
    
    TH1D*               outputHisto                                 = NULL;
    TGraph*             tempGraph                                   = NULL;
    TGraphErrors*       tempGraphErrors                             = NULL;
    TGraphAsymmErrors*  tempGraphAsymmErrors                        = NULL;
    
    if (inputObjectClassName.CompareTo("TGraphErrors") == 0) {
        
        cout << inputObject->GetName() << " is a TGraphErrors" << endl;
        
        tempGraphErrors                                             = (TGraphErrors*)inputObject->Clone("tempGraphErrors");
        Double_t*   xValue                                          = tempGraphErrors->GetX();
        Double_t*   yValue                                          = tempGraphErrors->GetY();
        Double_t*   xError                                          = tempGraphErrors->GetEX();
        Double_t*   yError                                          = tempGraphErrors->GetEY();
        Int_t       nPoints                                         = tempGraphErrors->GetN();

        Double_t*   newBinningX                                     = new Double_t[nPoints+1];
        for(Int_t i = 0; i < nPoints; i++) {
            newBinningX[i]                                          = xValue[i] - xError[i];
        }
        newBinningX[nPoints]                                        = xValue[nPoints-1] + xError[nPoints-1];
        
        outputHisto                                                 = new TH1D(outputHistoName,"",nPoints,newBinningX);
        for(Int_t i = 1; i <= nPoints; i++){
            outputHisto->SetBinContent(i,   yValue[i-1]);
            outputHisto->SetBinError(i,     yError[i-1]);
        }

        delete[] newBinningX;
    } else if (inputObjectClassName.CompareTo("TGraphAsymmErrors") == 0) {
        tempGraphAsymmErrors                                        = (TGraphAsymmErrors*)inputObject->Clone("tempGraphAsymmErrors");
        Double_t*   xValue                                          = tempGraphAsymmErrors->GetX();
        Double_t*   yValue                                          = tempGraphAsymmErrors->GetY();
        Double_t*   xErrorLow                                       = tempGraphAsymmErrors->GetEXlow();
        Double_t*   xErrorHigh                                      = tempGraphAsymmErrors->GetEXhigh();
        Double_t*   yErrorLow                                       = tempGraphAsymmErrors->GetEYlow();
        Double_t*   yErrorHigh                                      = tempGraphAsymmErrors->GetEYhigh();
        Int_t       nPoints                                         = tempGraphAsymmErrors->GetN();

        Double_t*   newBinningX                                     = new Double_t[nPoints+1];
        for(Int_t i = 0; i < nPoints; i++) {
            newBinningX[i]                                          = xValue[i] - xErrorLow[i];
        }
        newBinningX[nPoints]                                        = xValue[nPoints-1] + xErrorHigh[nPoints-1];

        outputHisto                                                 = new TH1D(outputHistoName,"",nPoints,newBinningX);
        for(Int_t i = 1; i <= nPoints; i++){
            outputHisto->SetBinContent(i,   yValue[i-1]);
            outputHisto->SetBinError(i,     (yErrorLow[i-1] + yErrorHigh[i-1]) / 2);
        }

        delete[] newBinningX;
    } else if (inputObjectClassName.Contains("TH1")) {
        outputHisto                                                 = (TH1D*)inputObject->Clone(outputHistoName.Data());
    } else {
        cout << inputObject->GetName() << " type not recognized or not supported" << endl;
        return NULL;
    }
    
    return outputHisto;
}

// fct. to determine the plateua onset of a ratio ************************************
TF1* FitPlateau(TH1D* ratio, Double_t nSigmas = 1, Double_t statisticsThresh = 0.1) {
    
    TF1*        tempConstFit                                        = new TF1("tempConstFit", "[0]",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
    Double_t    plateauOnset                                        = 0.;
    Double_t    tempConst[2]                                        = {0., 0.};
    Double_t    tempConstErr[2]                                     = {0., 0.};
    Int_t       onsetBin                                            = ratio->GetNbinsX();
    for (Int_t i = ratio->GetNbinsX(); i > 0; i--) {
        tempConst[1]                                                = tempConst[0];
        tempConstErr[1]                                             = tempConstErr[0];
        
        ratio->Fit(tempConstFit,"QMNRIE+","",ratio->GetXaxis()->GetBinLowEdge(i),ratio->GetXaxis()->GetXmax());
        tempConst[0]                                                = tempConstFit->GetParameter(0);
        tempConstErr[0]                                             = tempConstFit->GetParError(0);
        
        if (i == ratio->GetNbinsX()) {
            tempConst[1]                                            = tempConst[0];
            tempConstErr[1]                                         = tempConstErr[0];
        }
        
        if (ratio->GetBinError(i)/ratio->GetBinContent(i) > statisticsThresh) continue;
        if ( !( ((tempConst[1] + nSigmas*tempConstErr[1]) > tempConst[0]) && ((tempConst[1] - nSigmas*tempConstErr[1]) < tempConst[0]) ) ) break;
        
        //cout << "bin " << i << " (" << ratio->GetXaxis()->GetBinLowEdge(i) << "):\tcurrent = " << tempConst[0] << " +/- " << tempConstErr[0] << "\tvs. prev = " << tempConst[1] << " +/- " << tempConstErr[1]  << endl;
        
        onsetBin                                                    = i;
    }
    plateauOnset                                                    = ratio->GetXaxis()->GetBinLowEdge(onsetBin);
    
    ratio->Fit(tempConstFit,"QMNRIE+","",plateauOnset,ratio->GetXaxis()->GetXmax());
    return tempConstFit;
}

// fct. to fill mt scaling factors histo *********************************************
TH1D* FillMtScalingFactorsHisto(TString histoName, Double_t factors[nParticles], Double_t errors[nParticles]) {
    
    TH1D* histo                                                     = new TH1D(histoName.Data(), "", nParticles, 0, nParticles);
    for (Int_t i = 1; i < histo->GetNbinsX()+1; i++) {
        histo->GetXaxis()->SetBinLabel( i,  Form("%d",particlePDGCode[i-1]));
        histo->SetBinContent(           i,  factors[i-1]);
        histo->SetBinError(             i,  errors[i-1]);
    }
    return histo;
}

// fct. to calculate ratio TH1D/TF1 **************************************************
TH1D* CalculateRatioToFit(TH1D* histo, TF1* fit, TString name) {
    
    if (!histo || !fit) return NULL;

    // fit range
    Double_t ptMin                      = 0;
    Double_t ptMax                      = 0;
    fit->GetRange(ptMin, ptMax);
    
    // ratio to fit
    TH1D* histoRatioToFit               = (TH1D*)histo->Clone(name);
    histoRatioToFit->Sumw2();
    histoRatioToFit->Divide(fit);
    
    // set bins outside of function range to zero
    for (Int_t i=0; i<histoRatioToFit->GetNbinsX()+1; i++) {
        if (histoRatioToFit->GetXaxis()->GetBinUpEdge(i) <= ptMin || histoRatioToFit->GetXaxis()->GetBinLowEdge(i) >= ptMax) {
            histoRatioToFit->SetBinContent(i,   0);
            histoRatioToFit->SetBinError(i,     0);
        }
    }
    
    return histoRatioToFit;
}

// fct. to plot paricle yield ********************************************************
void PlotYield( TH1D*   yield,
                TF1*    paramYield,
                TH1D*   ratioToParamYield,
                TF1*    paramYieldMt,
                TH1D*   ratioToParamYieldMt,
                TF1*    paramYieldMtFeedDownCorr,
                TH1D*   ratioToParamYieldMtFeedDownCorr,
                TF1*    paramYieldMtPhi,
                TH1D*   ratioToParamYieldMtPhi,
                Int_t   nParticle,
                TString namePlot,
                TString suffix
               ) {
    
    if (!yield) return;
    if (!paramYield && !paramYieldMt && !paramYieldMtFeedDownCorr && !paramYieldMtPhi) return;
    if (!ratioToParamYield && !ratioToParamYieldMt && !ratioToParamYieldMtFeedDownCorr && !ratioToParamYieldMtPhi) return;

    Int_t       nParams                     = 0;
    
    Double_t    xMin                        = yield->GetXaxis()->GetXmin();
    Double_t    xMax                        = yield->GetXaxis()->GetXmax();

    Double_t    yMin                        = yield->GetMinimum(0) * 0.8;
    Double_t    yMax                        = yield->GetMaximum() * 2.0;
    
    TString     xTitle                      = "#it{p}_{T}(GeV/#it{c})";
    TString     yTitle                      = "#frac{1}{#it{N}_{ev}} #frac{d#it{N}^{2}}{d#it{p}_{T}d#it{y}} ((GeV/#it{c})^{-1})";

    // create canvas and pads
    TCanvas* canvas                         = GetAndSetCanvas("canvas");
    DrawGammaCanvasSettings(canvas, 0.11, 0.015, 0.02, 0.09);
    TPad* pad                               = new TPad("pad", "", 0., 0.25, 1., 1.,-1, -1, -2);
    DrawGammaPadSettings(pad, 0.11, 0.015, 0.02, 0.);
    pad->Draw();
    TPad* padRatioToFit                     = new TPad("padRatioToFit", "", 0., 0., 1., 0.25,-1, -1, -2);
    DrawGammaPadSettings(padRatioToFit, 0.11, 0.015, 0.0, 0.3);
    padRatioToFit->Draw();

    // plot yield + params
    pad->cd();
    pad->SetLogy();
    pad->SetLogx();

    SetHistogramm(  yield, xTitle, yTitle, yMin, yMax);
    DrawGammaSetMarker(yield, 24, 2.0, kBlack, kBlack);
    yield->Draw("e1");
    
    if (paramYield) {
        nParams++;
        paramYield->SetLineColor(kRed);
        paramYield->SetLineWidth(4.0);
        paramYield->Draw("e1,same");
    }
    
    if (paramYieldMt) {
        nParams++;
        paramYieldMt->SetLineColor(kBlue);
        paramYieldMt->SetLineWidth(4.0);
        paramYieldMt->Draw("e1,same");
    }
    
    if (paramYieldMtFeedDownCorr) {
        nParams++;
        paramYieldMtFeedDownCorr->SetLineColor(kOrange-3);
        paramYieldMtFeedDownCorr->SetLineWidth(4.0);
        paramYieldMtFeedDownCorr->Draw("e1,same");
    }
    
    if (paramYieldMtPhi) {
        nParams++;
        paramYieldMtPhi->SetLineColor(kAzure+1);
        paramYieldMtPhi->SetLineWidth(4.0);
        paramYieldMtPhi->Draw("e1,same");
    }

    // legend
    TLegend* legend                         = GetAndSetLegend(0.6, 0.65, nParams+1, 1, Form("%s %s", centralityString.Data(), collisionSystem.Data()));
    legend->AddEntry(yield, Form("%s yield", particleNameLatex[nParticle].Data()), "p");
    if (paramYield)                 legend->AddEntry(paramYield,                "yield param.", "l");
    if (paramYieldMt)               legend->AddEntry(paramYieldMt,              "m_{T} scaled #pi^{0}", "l");
    if (paramYieldMtFeedDownCorr)   legend->AddEntry(paramYieldMtFeedDownCorr,  "m_{T} scaled #pi^{0} res. feed down corr.", "l");
    if (paramYieldMtPhi)            legend->AddEntry(paramYieldMtPhi,           "m_{T} scaled #phi", "l");
    legend->Draw();

    // plot ratio
    padRatioToFit->cd();
    padRatioToFit->SetLogx();

    if (ratioToParamYield) {
        SetStyleHistoTH1ForGraphs(ratioToParamYield, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}", 0.14, 0.15, 0.12, 0.10,  0.85, 0.4, 510, 505);
        DrawGammaSetMarker(ratioToParamYield, 20, 2.0, kRed, kRed);
        ratioToParamYield->GetYaxis()->SetRangeUser(0.4,1.6);
        ratioToParamYield->Draw("e1,same");
    }
    
    if (ratioToParamYieldMt) {
        SetStyleHistoTH1ForGraphs(ratioToParamYieldMt, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}", 0.14, 0.15, 0.12, 0.10,  0.85, 0.4, 510, 505);
        DrawGammaSetMarker(ratioToParamYieldMt, 20, 2.0, kBlue, kBlue);
        ratioToParamYieldMt->GetYaxis()->SetRangeUser(0.4,1.6);
        ratioToParamYieldMt->Draw("e1,same");
    }
    
    if (ratioToParamYieldMtFeedDownCorr) {
        SetStyleHistoTH1ForGraphs(ratioToParamYieldMtFeedDownCorr, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}", 0.14, 0.15, 0.12, 0.10,  0.85, 0.4, 510, 505);
        DrawGammaSetMarker(ratioToParamYieldMtFeedDownCorr, 20, 2.0, kOrange-3, kOrange-3);
        ratioToParamYieldMtFeedDownCorr->GetYaxis()->SetRangeUser(0.4,1.6);
        ratioToParamYieldMtFeedDownCorr->Draw("e1,same");
    }
    
    if (ratioToParamYieldMtPhi) {
        SetStyleHistoTH1ForGraphs(ratioToParamYieldMtPhi, "#it{p}_{T} (GeV/#it{c})","#frac{data}{fit}", 0.14, 0.15, 0.12, 0.10,  0.85, 0.4, 510, 505);
        DrawGammaSetMarker(ratioToParamYieldMtPhi, 20, 2.0, kAzure+1, kAzure+1);
        ratioToParamYieldMtPhi->GetYaxis()->SetRangeUser(0.4,1.6);
        ratioToParamYieldMtPhi->Draw("e1,same");
    }

    DrawGammaLines(xMin, xMax, 1.0, 1.0, 2.0, kGray+2, 2);
    DrawGammaLines(xMin, xMax, 1.2, 1.2, 2.0, kGray+2, 8);
    DrawGammaLines(xMin, xMax, 0.8, 0.8, 2.0, kGray+2, 8);

    canvas->SaveAs(Form("%s/%s.%s", outputDir.Data(), namePlot.Data(), suffix.Data()));
}






