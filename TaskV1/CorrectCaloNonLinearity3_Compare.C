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
#include "TProfile2D.h"
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
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

Float_t FunctionNL_kTestBeamv3(Float_t e){
    return ( 0.9615 / ( 0.976941 *( 1. / ( 1. + 0.162310 * exp( -e / 1.08689 ) ) * 1. / ( 1. + 0.0819592 * exp( ( e - 152.338 ) / 30.9594 ) ) ) ) );
}
Float_t FunctionNL_kPi0MCv3(Float_t e){
    return ( 1.0 / ( 0.981039 * ( 1. / ( 1. + 0.113508 * exp( -e / 1.00173 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}
Float_t FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6){
    return ( p6 / ( p0 * ( 1. / ( 1. + p1 * exp( -e / p2 ) ) * 1. / ( 1. + p3 * exp( ( e - p4 ) / p5 ) ) ) ) );
}
Float_t FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2){
    return ( p0 + exp( p1 + ( p2 * e ) ) );
}
Float_t FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5){
  Float_t ret = 1;
  if ((p3 +  p4 * TMath::Power(e,p5 ) ) != 0)
    ret = ( (p0 +  p1 * TMath::Power(e,p2 ) )/(p3 +  p4 * TMath::Power(e,p5 ) ) );
  if (ret != 0.)
    return ret;
  else
    return 1.;
}

//****************************************************************************
//************** Function to compare different CaloNonLinearities ************
//****************************************************************************
void CorrectCaloNonLinearity3_Compare(TString select = "LHC11a-Pythia")
{
    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();

    TString inputDir                = "CorrectCaloNonLinearity3";
    TString outputDir               = Form("%s/Compare",inputDir.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    TString suffix                  = "eps";

    Int_t nNL                       = 0;
    TString *inputFileNames         = 0x0;
    TString *inputFilePaths         = 0x0;
    Bool_t *plotMassData            = 0x0;
    TString xTitle                  = "#it{E}_{Cluster} (GeV)";
    TFile **inputFiles              = 0x0;

    const Int_t nColor              = 11;
    const Int_t nMarkerStyle        = 10;
    Color_t color[nColor]           = { kBlack, 633, 807, /*800,*/ 418, 
                                        kYellow+1, 601, 879, 806, 852,
                                        kCyan+3, 426};
    Int_t markerStyle[nMarkerStyle] = {24, 25, 27, 28, 29, 
                                       30, 31, 20, 33, 34};

    Bool_t isDataVariation          = kFALSE;
    
//*******************************************************************************
// Choosing data set
    if(select.CompareTo("LHC11a-Pythia")==0){
        nNL                         = 2;
        inputFileNames              = new TString[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-ConvCalo";
        inputFileNames[1]           = "LHC11a-Pythia-Calo";
    } else if(select.CompareTo("LHC11a-Phojet")==0){
        nNL                         = 2;
        inputFileNames              = new TString[nNL];
        inputFileNames[0]           = "LHC11a-Phojet-ConvCalo";
        inputFileNames[1]           = "LHC11a-Phojet-Calo";
    } else if(select.CompareTo("ConvCalo")==0){
        nNL                         = 2;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-ConvCalo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-ConvCalo";
        plotMassData[1]             = kFALSE;
    } else if(select.CompareTo("Calo")==0){
        nNL                         = 2;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-Calo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-Calo";
        plotMassData[1]             = kFALSE;
    } else if(select.CompareTo("2.76TeV-OpenTime")==0){
        nNL                         = 6;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-ConvCalo-OpenTime";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Pythia-Calo-OpenTime";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC11a-Phojet-ConvCalo-OpenTime";
        plotMassData[2]             = kTRUE;
        inputFileNames[3]           = "LHC11a-Phojet-Calo-OpenTime";
        plotMassData[3]             = kTRUE;
        inputFileNames[4]           = "LHC13g-Pythia-ConvCalo-OpenTime";
        plotMassData[4]             = kTRUE;
        inputFileNames[5]           = "LHC13g-Pythia-Calo-OpenTime";
        plotMassData[5]             = kTRUE;
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("2.76TeV-JetJet-OpenTime")==0){
        nNL                         = 4;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-JetJet-ConvCalo-OpenTime";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-JetJet-Calo-OpenTime";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC13g-JetJet-ConvCalo-OpenTime";
        plotMassData[2]             = kTRUE;
        inputFileNames[3]           = "LHC13g-JetJet-Calo-OpenTime";
        plotMassData[3]             = kTRUE;
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("LHC11a-ConvCalo")==0){
        nNL                         = 2;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-ConvCalo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-ConvCalo";
        plotMassData[1]             = kFALSE;
    } else if(select.CompareTo("LHC11a-ConvCalo-Pythia-DiffTiming")==0){
        nNL                         = 3;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-ConvCalo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Pythia-ConvCalo-OpenTimeDeltaT";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC11a-Pythia-ConvCalo-OpenTime";
        plotMassData[2]             = kTRUE;
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("LHC11a-ConvCalo-Phojet-DiffTiming")==0){
        nNL                         = 3;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Phojet-ConvCalo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-ConvCalo-OpenTimeDeltaT";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC11a-Phojet-ConvCalo-OpenTime";
        plotMassData[2]             = kTRUE;        
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("LHC11a-Calo-Pythia-DiffTiming")==0){
        nNL                         = 3;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-Calo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Pythia-Calo-OpenTimeDeltaT";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC11a-Pythia-Calo-OpenTime";
        plotMassData[2]             = kTRUE;
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("LHC11a-Calo-Phojet-DiffTiming")==0){
        nNL                         = 3;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Phojet-Calo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-Calo-OpenTimeDeltaT";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC11a-Phojet-Calo-OpenTime";
        plotMassData[2]             = kTRUE;        
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("LHC11a-Pythia-DiffTiming")==0){
        nNL                         = 6;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-Calo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Pythia-ConvCalo";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC11a-Pythia-Calo-OpenTimeDeltaT";
        plotMassData[2]             = kTRUE;
        inputFileNames[3]           = "LHC11a-Pythia-ConvCalo-OpenTimeDeltaT";
        plotMassData[3]             = kTRUE;
        inputFileNames[4]           = "LHC11a-Pythia-Calo-OpenTime";
        plotMassData[4]             = kTRUE;
        inputFileNames[5]           = "LHC11a-Pythia-ConvCalo-OpenTime";
        plotMassData[5]             = kTRUE;
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("LHC11a-Phojet-DiffTiming")==0){
        nNL                         = 6;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Phojet-Calo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-ConvCalo";
        plotMassData[1]             = kTRUE;
        inputFileNames[2]           = "LHC11a-Phojet-Calo-OpenTimeDeltaT";
        plotMassData[2]             = kTRUE;
        inputFileNames[3]           = "LHC11a-Phojet-ConvCalo-OpenTimeDeltaT";
        plotMassData[3]             = kTRUE;
        inputFileNames[4]           = "LHC11a-Phojet-Calo-OpenTime";
        plotMassData[4]             = kTRUE;
        inputFileNames[5]           = "LHC11a-Phojet-ConvCalo-OpenTime";
        plotMassData[5]             = kTRUE;
        isDataVariation             = kTRUE;
    } else if(select.CompareTo("LHC11a-Calo")==0){
        nNL                         = 2;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-Calo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-Calo";
        plotMassData[1]             = kFALSE;
    } else if(select.CompareTo("LHC11a-ConvCalo-Calo")==0){
        nNL                         = 4;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC11a-Pythia-ConvCalo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC11a-Phojet-ConvCalo";
        plotMassData[1]             = kFALSE;
        inputFileNames[2]           = "LHC11a-Pythia-Calo";
        plotMassData[2]             = kTRUE;
        inputFileNames[3]           = "LHC11a-Phojet-Calo";
        plotMassData[3]             = kFALSE;
    } else if(select.CompareTo("LHC10-ConvCalo-Calo")==0){
        nNL                         = 2;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC10-ConvCalo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC10-Calo";
        plotMassData[1]             = kFALSE;
    } else if(select.CompareTo("LHC13bc-ConvCalo-Calo")==0){
        nNL                         = 3;
        inputFileNames              = new TString[nNL];
        plotMassData                = new Bool_t[nNL];
        inputFileNames[0]           = "LHC13bc-ConvCalo";
        plotMassData[0]             = kTRUE;
        inputFileNames[1]           = "LHC13bc-Calo";
        plotMassData[1]             = kFALSE;
        inputFileNames[2]           = "LHC13bc-Calo2";
        plotMassData[2]             = kFALSE;
    } else if(select.CompareTo("LHC12-Pythia")==0){
      nNL                         = 2;
      inputFileNames              = new TString[nNL];
      inputFileNames[0]           = "LHC12-Pythia-ConvCalo";
      inputFileNames[1]           = "LHC12-Pythia-Calo";
    } else if(select.CompareTo("LHC12-Phojet")==0){
      nNL                         = 2;
      inputFileNames              = new TString[nNL];
      inputFileNames[0]           = "LHC12-Phojet-ConvCalo";
      inputFileNames[1]           = "LHC12-Phojet-Calo";
    } else if(select.CompareTo("LHC12-ConvCalo")==0){
      nNL                         = 2;
      inputFileNames              = new TString[nNL];
      plotMassData                = new Bool_t[nNL];
      inputFileNames[0]           = "LHC12-Pythia-ConvCalo";
      plotMassData[0]             = kTRUE;
      inputFileNames[1]           = "LHC12-Phojet-ConvCalo";
      plotMassData[1]             = kFALSE;
    } else if(select.CompareTo("LHC12-Calo")==0){
      nNL                         = 2;
      inputFileNames              = new TString[nNL];
      plotMassData                = new Bool_t[nNL];
      inputFileNames[0]           = "LHC12-Pythia-Calo";
      plotMassData[0]             = kTRUE;
      inputFileNames[1]           = "LHC12-Phojet-Calo";
      plotMassData[1]             = kFALSE;
    } else if(select.CompareTo("LHC12-ConvCalo-Calo")==0){
      nNL                         = 4;
      inputFileNames              = new TString[nNL];
      plotMassData                = new Bool_t[nNL];
      inputFileNames[0]           = "LHC12-Pythia-ConvCalo";
      plotMassData[0]             = kTRUE;
      inputFileNames[1]           = "LHC12-Phojet-ConvCalo";
      plotMassData[1]             = kFALSE;
      inputFileNames[2]           = "LHC12-Pythia-Calo";
      plotMassData[2]             = kTRUE;
      inputFileNames[3]           = "LHC12-Phojet-Calo";
      plotMassData[3]             = kFALSE;
    } else{
        cout << "No valid selection '" << select.Data() << "'' given, returning..." << endl;
        return;
    }

    if(nNL>nColor||nNL>nMarkerStyle){
        cout << "nNL: " << nNL << ", but defined nColor " << nColor << " and defined nMarkerStyle " << nMarkerStyle << endl;
        return;
    }

//*******************************************************************************
// Input

    inputFilePaths = new TString[nNL];
    inputFiles = new TFile*[nNL];
    for(Int_t i=0; i<nNL; i++){
        inputFilePaths[i] = Form("%s/CorrectCaloNonLinearity3_%s.root",inputDir.Data(),inputFileNames[i].Data());
        inputFiles[i] = new TFile(inputFilePaths[i].Data(),"READ");
            if(inputFiles[i]->IsZombie()) {cout << "Info: ROOT file '" << inputFilePaths[i].Data() << "' could not be openend, return!" << endl; return;}
    }

//*******************************************************************************
// Output

    TCanvas *canvas     = new TCanvas("canvas","",200,0,1350,900);  // gives the page size
    Double_t leftMargin = 0.1; Double_t rightMargin = 0.02; Double_t topMargin = 0.06; Double_t bottomMargin = 0.1;
    DrawGammaCanvasSettings(canvas, 0.1, 0.02, 0.02, 0.1);
    canvas->SetLogx(1); canvas->SetLogy(0);
    
    Int_t nColumns      = 1;
    if (nNL > 3) 
        nColumns        = 2;
    TLegend *legend     = GetAndSetLegend2(0.13, 0.95-(nNL/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42);
    legend->SetMargin(0.1);
//*******************************************************************************
// plotting masses MC+Data
    TH1D* histoMassMC[nNL];
    for(Int_t i=0; i<nNL; i++){
        TString drawOption = (i==0)?"p":"p, same";
        histoMassMC[i] = (TH1D*)inputFiles[i]->Get("Mean mass MC");
        if(!histoMassMC[i]){
        cout << "ERROR: Could not find histogram 'Mean mass MC' in " << inputFilePaths[i].Data() << endl;
        continue;
        }
        DrawGammaSetMarker(histoMassMC[i], markerStyle[i], 1, color[i], color[i]);
        histoMassMC[i]->GetYaxis()->SetTitleOffset(1.1);
        histoMassMC[i]->GetXaxis()->SetLabelOffset(-0.01);
        histoMassMC[i]->Draw(drawOption.Data());
        legend->AddEntry(histoMassMC[i],inputFileNames[i].Data());
    }

    legend->Draw("same");
    canvas->Update();
    canvas->SaveAs(Form("%s/Mass_MC_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvas->Clear();

    legend->Clear();
    TH1D* histoMassData[nNL];
    for(Int_t i=0; i<nNL; i++){
        TString plotNameData = inputFileNames[i];
        if(plotMassData && !plotMassData[i]) continue;
        else{
            TObjArray *d = plotNameData.Tokenize("-");
            TObjString* rString = (TObjString*)d->At(0);
            plotNameData = rString->GetString();
            if (!isDataVariation){
                for(Int_t iObjArr=1; iObjArr<d->GetEntries(); iObjArr++){
                    TObjString* tmpObjString = (TObjString*)d->At(iObjArr);
                    TString tmpString = tmpObjString->GetString();
                    if(tmpString.Contains("Calo")){
                        plotNameData+="-";
                        plotNameData+=tmpString;
                    }
                }
            } else {
                for(Int_t iObjArr=1; iObjArr<d->GetEntries(); iObjArr++){
                    TObjString* tmpObjString = (TObjString*)d->At(iObjArr);
                    TString tmpString = tmpObjString->GetString();
                    if( tmpString.CompareTo("Pythia") != 0 && tmpString.CompareTo("Phojet") != 0 &&
                        tmpString.CompareTo("JetJet") != 0 ){
                        plotNameData+="-";
                        plotNameData+=tmpString;
                    }
                }
            }    
        }
        TString drawOption = (i==0)?"p":"p, same";
        histoMassData[i] = (TH1D*)inputFiles[i]->Get("Mean mass Data");
        if(!histoMassData[i]){
        cout << "ERROR: Could not find histogram 'Mean mass Data' in " << inputFilePaths[i].Data() << endl;
        continue;
        }
        DrawGammaSetMarker(histoMassData[i], markerStyle[i], 1, color[i], color[i]);
        histoMassData[i]->GetXaxis()->SetLabelOffset(-0.01);
        histoMassData[i]->GetYaxis()->SetTitleOffset(1.1);
        histoMassData[i]->Draw(drawOption.Data());
        legend->AddEntry(histoMassData[i],plotNameData.Data());
    }

    legend->Draw("same");
    canvas->Update();
    canvas->SaveAs(Form("%s/Mass_Data_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvas->Clear();
    legend->Clear();

//*******************************************************************************
// plotting mass ratios
    TCanvas *canvasRatio     = new TCanvas("canvasRatio","",0,0,750,500);  // gives the page size
    DrawGammaCanvasSettings(canvasRatio, 0.08, 0.02, 0.02, 0.08);
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0);

    TH1D* histoMeanMassRatio[nNL];
    for(Int_t i=0; i<nNL; i++){
        TString drawOption = (i==0)?"p":"p, same";
        histoMeanMassRatio[i] = (TH1D*)inputFiles[i]->Get("MeanMassRatioMCData-noFit");
        if(!histoMeanMassRatio[i]){
        cout << "ERROR: Could not find histogram 'MeanMassRatioMCData-noFit' in " << inputFilePaths[i].Data() << endl;
        continue;
        }
        histoMeanMassRatio[i]->SetTitle(" ");
        DrawGammaSetMarker(histoMeanMassRatio[i], markerStyle[i], 1, color[i], color[i]);
        if(select.Contains("LHC10-")) histoMeanMassRatio[i]->GetYaxis()->SetRangeUser(0.9,1.05);
//         histoMeanMassRatio[i]->GetYaxis()->SetTitleOffset(1.1);
        histoMeanMassRatio[i]->GetXaxis()->SetLabelOffset(-0.01);
        histoMeanMassRatio[i]->GetXaxis()->SetTitleOffset(1.1);
        histoMeanMassRatio[i]->Draw(drawOption.Data());
        legend->AddEntry(histoMeanMassRatio[i],inputFileNames[i].Data());
    }

    legend->Draw("same");
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
    canvasRatio->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasRatio->Clear();
    legend->Clear();

    TF1* histoMeanMassRatioFits[nNL];
    histoMeanMassRatio[0]->Draw("axis");
    for(Int_t i=0; i<nNL; i++){
        histoMeanMassRatioFits[i] = (TF1*)inputFiles[i]->Get("MeanMassRatioMCData-Fit");
        if(!histoMeanMassRatioFits[i]){
        cout << "ERROR: Could not find histogram 'MeanMassRatioMCData-Fit' in " << inputFilePaths[i].Data() << endl;
        continue;
        }
        histoMeanMassRatioFits[i]->SetTitle(" ");
        histoMeanMassRatioFits[i]->SetLineColor(color[i]);
        histoMeanMassRatioFits[i]->GetXaxis()->SetTitle(xTitle.Data());
        if(select.Contains("LHC10-")) histoMeanMassRatioFits[i]->GetYaxis()->SetRangeUser(0.9,1.05);
        histoMeanMassRatioFits[i]->Draw("same");
        legend->AddEntry(histoMeanMassRatioFits[i],inputFileNames[i].Data());
    }

    legend->Draw("same");
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
    canvasRatio->SaveAs(Form("%s/MeanMassRatio_Fits_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasRatio->Clear();
    legend->Clear();


    for(Int_t i=0; i<nNL; i++){
        TString drawOption = (i==0)?"p":"p, same";
        histoMeanMassRatio[i]->Draw(drawOption.Data());
        histoMeanMassRatioFits[i]->Draw("same");
        legend->AddEntry(histoMeanMassRatio[i],inputFileNames[i].Data());
    }

    legend->Draw("same");
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
    canvasRatio->SaveAs(Form("%s/MeanMassRatioAndFits_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasRatio->Clear();
    legend->Clear();

//*******************************************************************************
// plotting total correction
    TH1D* histoTotalCorrection[nNL];
    for(Int_t i=0; i<nNL; i++){
        TString drawOption = (i==0)?"p":"p, same";
        histoTotalCorrection[i] = (TH1D*)inputFiles[i]->Get("Total Correction");
        if(!histoTotalCorrection[i]){
        cout << "ERROR: Could not find histogram 'Total Correction' in " << inputFilePaths[i].Data() << endl;
        continue;
        }
        DrawGammaSetMarker(histoTotalCorrection[i], markerStyle[i], 0.2, color[i], color[i]);
        histoTotalCorrection[i]->Draw(drawOption.Data());
        legend->AddEntry(histoTotalCorrection[i],inputFileNames[i].Data(),"l");
    }

    TH1D* testBeam = (TH1D*) histoTotalCorrection[0]->Clone("testbeam");
    testBeam->Reset("ICE");
    DrawGammaSetMarker(testBeam, markerStyle[nNL], 0.2, color[nNL], color[nNL]);
    for(Int_t iBin = 1; iBin <= testBeam->GetNbinsX()+1; iBin++) {
      Float_t e = testBeam->GetXaxis()->GetBinCenter(iBin);
      Float_t factor = 1;
      factor *= FunctionNL_kPi0MCv3(e);
      factor /= FunctionNL_kTestBeamv3(e);
      testBeam->SetBinContent(iBin,factor);
    }
    legend->AddEntry(testBeam,"kPi0MCv3 / kTestBeamv3","l");
    testBeam->Draw("p,same");

    legend->Draw("same");
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
    canvasRatio->SaveAs(Form("%s/TotalCorrection_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasRatio->Clear();
    legend->Clear();

    if(select.Contains("LHC11a-ConvCalo-Calo")){
        legend->SetNColumns(2);
        legend->SetX2(0.9);

        Int_t iPlot = 0;
        testBeam->Reset("ICE");
        DrawGammaSetMarker(testBeam,  markerStyle[nNL], 0.2, color[nNL], color[nNL]);
        for(Int_t iBin = 1; iBin <= testBeam->GetNbinsX()+1; iBin++) {
            Float_t e = testBeam->GetXaxis()->GetBinCenter(iBin);
            Float_t factor = 1;
            factor *= FunctionNL_kPi0MCv3(e);
            factor /= FunctionNL_kTestBeamv3(e);
            testBeam->SetBinContent(iBin,factor);
        }
        legend->AddEntry(testBeam,"kPi0MCv3 / kTestBeamv3","l");
        testBeam->DrawCopy("p");

        for(Int_t iNL=0; iNL<4; iNL++){
            TH1D* test = (TH1D*) histoTotalCorrection[0]->Clone(Form("NL%i",iNL));
            test->Reset("ICE");
            DrawGammaSetMarker(test, markerStyle[iNL], 0.2, color[iNL], color[iNL]);
            for(Int_t iBin = 1; iBin <= test->GetNbinsX()+1; iBin++) {
                Float_t e = test->GetXaxis()->GetBinCenter(iBin);
                Float_t factor = 1;
                if(iNL==0) factor /= FunctionNL_kSDM(e, 0.984889*0.995*0.9970, -3.65456, -1.12744);
                else if(iNL==1) factor /= FunctionNL_kSDM(e, 0.984384*0.995*0.9970, -3.30287, -1.48516);
                else if(iNL==2) factor /= FunctionNL_kSDM(2.0*e, 0.966151*0.995*0.9981, -2.97974, -0.29463);
                else if(iNL==3) factor /= FunctionNL_kSDM(2.0*e, 0.988814*0.995*0.9981, 0.335011, -4.30322);
                test->SetBinContent(iBin,factor);
            }
            if(iNL==0) legend->AddEntry(test,"ConvCalo-Pythia","l");
            else if(iNL==1) legend->AddEntry(test,"ConvCalo-Phojet","l");
            else if(iNL==2) legend->AddEntry(test,"Calo-Pythia","l");
            else if(iNL==3) legend->AddEntry(test,"Calo-Phojet","l");
            test->DrawCopy("p,same");
        }

        legend->Draw("same");
        canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
        canvasRatio->SaveAs(Form("%s/TotalCorrection_Full_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
        canvasRatio->Clear();
        legend->Clear();
    }

//**************************************************************************************************************************************
//**************************************************************************************************************************************
//**************************************************************************************************************************************

    // Case 1, strict timing - 2.76TeV : Mass Ratio Fitting
    if(select.Contains("2.76TeV")){
      Int_t n=4;
        legend->SetNColumns(2);
        legend->SetY1(1.05);
        legend->SetX2(1.5);

        TH1D* testB = new TH1D(*testBeam);
        testB->Reset("ICE");
        DrawGammaSetMarker(testB,  20, 0.8, kBlack, kBlack);
        for(Int_t iBin = 1; iBin <= testBeam->GetNbinsX()+1; iBin++) {
            if(iBin%10>=1) continue;
            Float_t e = testB->GetXaxis()->GetBinCenter(iBin);
            Float_t factor = 1;
            factor *= FunctionNL_kPi0MCv3(e);
            factor /= FunctionNL_kTestBeamv3(e);
            testB->SetBinContent(iBin,factor);
        }

        testB->GetXaxis()->SetTitleOffset(0.8);
        testB->GetYaxis()->SetRangeUser(0.98,1.1);
        testB->DrawCopy("p");

        Color_t clr[4] = { kAzure, kRed, kAzure+4, kRed+1};

        TH1D* testBarr[n];
        for(Int_t iNL=0; iNL<n; iNL++){
            testBarr[iNL] = (TH1D*) histoTotalCorrection[0]->Clone(Form("NL%i",iNL));
            testBarr[iNL]->Reset("ICE");
            DrawGammaSetMarker(testBarr[iNL], 21, 0.5, clr[iNL], clr[iNL]);
            testBarr[iNL]->SetLineWidth(3);
            for(Int_t iBin = 1; iBin <= testBarr[iNL]->GetNbinsX()+1; iBin++) {
                Float_t energy = testBarr[iNL]->GetXaxis()->GetBinCenter(iBin);
                Float_t factor = 1;
                if(iNL==0) factor /= FunctionNL_kSDM(energy, 0.984889*0.995*0.9970, -3.65456, -1.12744);
                else if(iNL==1) factor /= FunctionNL_kSDM(2.0*energy, 0.966151*0.995*0.9981, -2.97974, -0.29463);
                else if(iNL==2) factor /= FunctionNL_kSDM(energy, 0.983176*0.9945, -3.91107, -0.697613);
                else if(iNL==3) factor /= FunctionNL_kSDM(energy, 0.974358*0.9987, -2.18037, -1.91622);
//                else if(iNL==4) factor /= FunctionNL_kSDM(energy, 0.981892*0.995*0.9970, -5.43438, -1.05468);
//                else if(iNL==5) factor /= FunctionNL_kSDM(2.0*energy, 0.979994*0.995*0.9981, -3.24431, -0.760205);
                if( iNL < 2 && iBin%10>=1) continue;
                testBarr[iNL]->SetBinContent(iBin,factor);
            }
            if(iNL==0) legend->AddEntry(testBarr[iNL],"Pythia8 (LHC11a) - ConvCalo (50ns)","p");
            else if(iNL==1) legend->AddEntry(testBarr[iNL],"Pythia8 (LHC11a) - Calo (50ns)","p");
            else if(iNL==2) legend->AddEntry(testBarr[iNL],"Pythia8 (LHC11a) - ConvCalo (500ns)","l");
            else if(iNL==3) legend->AddEntry(testBarr[iNL],"Pythia8 (LHC11a) - Calo (500ns)","l");
//            else if(iNL==4) legend->AddEntry(testBarr[iNL],"JetJet - ConvCalo","p");
//            else if(iNL==5) legend->AddEntry(testBarr[iNL],"JetJet - Calo","p");
            testBarr[iNL]->DrawCopy("p,same");
        }
        testB->DrawCopy("p,same");
        legend->AddEntry(testB,"kPi0MCv3 / kTestBeamv3","p");
        legend->Draw("same");
        PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 0.03, "pp, #sqrt{#it{s}} = 2.76 TeV", "CCRF","");
        canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
        canvasRatio->SaveAs(Form("%s/TotalCorrection_Full_strict_Case1_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
        canvasRatio->Clear();
        legend->Clear();
    }

    //case 1, open timing - 2.76TeV : Mass Ratio Fitting
    if(select.Contains("2.76TeV")){
      Int_t n=10;
        legend->SetNColumns(2);
        legend->SetY1(1.05);
        legend->SetX2(1.5);

        TH1D* testB = new TH1D(*testBeam);
        testB->Reset("ICE");
        DrawGammaSetMarker(testB,  20, 0.8, kBlack, kBlack);
        for(Int_t iBin = 1; iBin <= testBeam->GetNbinsX()+1; iBin++) {
            if(iBin%10>=1) continue;
            Float_t e = testB->GetXaxis()->GetBinCenter(iBin);
            Float_t factor = 1;
            factor *= FunctionNL_kPi0MCv3(e);
            factor /= FunctionNL_kTestBeamv3(e);
            testB->SetBinContent(iBin,factor);
        }

        testB->GetXaxis()->SetTitleOffset(0.8);
        testB->GetYaxis()->SetRangeUser(0.98,1.1);
        testB->DrawCopy("p");

        Color_t clr[10] = { kAzure, kRed, kAzure+4, kRed-4, kAzure+8, kRed+2, kAzure-3, kRed-2, kAzure-6, kRed-6};

        TH1D* testBarr[n];
        for(Int_t iNL=0; iNL<n; iNL++){
            testBarr[iNL] = (TH1D*) histoTotalCorrection[0]->Clone(Form("NL%i",iNL));
            testBarr[iNL]->Reset("ICE");
            DrawGammaSetMarker(testBarr[iNL], 21, 0.5, clr[iNL], clr[iNL]);
            testBarr[iNL]->SetLineWidth(3);
            for(Int_t iBin = 1; iBin <= testBarr[iNL]->GetNbinsX()+1; iBin++) {
                Float_t energy = testBarr[iNL]->GetXaxis()->GetBinCenter(iBin);
                Float_t factor = 1;
                if(iNL==0) factor /= FunctionNL_kSDM(energy, 0.983176*0.9945, -3.91107, -0.697613);
                else if(iNL==1) factor /= FunctionNL_kSDM(energy, 0.974358*0.9987, -2.18037, -1.91622);
                else if(iNL==2) factor /= FunctionNL_kSDM(energy, 0.972574*0.9942, -3.19191, -0.946239);
                else if(iNL==3) factor /= FunctionNL_kSDM(energy, 0.963307*0.9962, -3.27998, -0.589806);
                else if(iNL==4) factor /= FunctionNL_kSDM(energy, 0.981893*0.9930, -4.05476, -0.710661);
                else if(iNL==5) factor /= FunctionNL_kSDM(energy, 0.97499*0.9995, -0.180148, -4.78066);
                else if(iNL==6) factor /= FunctionNL_kSDM(energy, 0.983176*0.993*0.99, -1.85546, -3.37696);
                else if(iNL==7) factor /= FunctionNL_kSDM(energy, 0.974424*0.998*0.992, -0.533785, -4.06374);
                else if(iNL==8) factor /= FunctionNL_kSDM(energy, 0.977035*0.9835, -3.82187, -1.04332);
                else if(iNL==9) factor /= FunctionNL_kSDM(energy, 0.963307*0.995, -4.01949, -0.38667);
                if(iNL>5 && iBin%10>=1) continue;
                testBarr[iNL]->SetBinContent(iBin,factor);
            }
            if(iNL==0) legend->AddEntry(testBarr[iNL],"Pythia8(LHC11a) - ConvCalo","l");
            else if(iNL==1) legend->AddEntry(testBarr[iNL],"Pythia8(LHC11a) - Calo","l");
            else if(iNL==2) legend->AddEntry(testBarr[iNL],"Pythia8(LHC13g) - ConvCalo","l");
            else if(iNL==3) legend->AddEntry(testBarr[iNL],"Pythia8(LHC13g) - Calo","l");
            else if(iNL==4) legend->AddEntry(testBarr[iNL],"Phojet - ConvCalo","l");
            else if(iNL==5) legend->AddEntry(testBarr[iNL],"Phojet - Calo","l");
            else if(iNL==6) legend->AddEntry(testBarr[iNL],"JetJet(LHC11a) - ConvCalo","p");
            else if(iNL==7) legend->AddEntry(testBarr[iNL],"JetJet(LHC11a) - Calo","p");
            else if(iNL==8) legend->AddEntry(testBarr[iNL],"JetJet(LHC13g) - ConvCalo","p");
            else if(iNL==9) legend->AddEntry(testBarr[iNL],"JetJet(LHC13g) - Calo","p");
            testBarr[iNL]->DrawCopy("p,same");
        }
        testB->DrawCopy("p,same");
        legend->AddEntry(testB,"kPi0MCv3 / kTestBeamv3","p");
        legend->Draw("same");
        PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 0.03, "pp, #sqrt{#it{s}} = 2.76 TeV", "CCRF","");
        canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
        canvasRatio->SaveAs(Form("%s/TotalCorrection_Full_open_Case1_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
        canvasRatio->Clear();
        legend->Clear();
    }

    // Case 2 - 2.76TeV : Mass Fitting, then Ratio
    if(select.Contains("2.76TeV")){
       Int_t n = 10;
        legend->SetNColumns(2);
        legend->SetY1(1.05);
        legend->SetX2(1.5);

        TH1D* testB = new TH1D(*testBeam);
        testB->Reset("ICE");
        DrawGammaSetMarker(testB,  20, 0.8, kBlack, kBlack);
        for(Int_t iBin = 1; iBin <= testB->GetNbinsX()+1; iBin++) {
            if(iBin%10>=1) continue;
            Float_t e = testB->GetXaxis()->GetBinCenter(iBin);
            Float_t factor = 1;
            factor *= FunctionNL_kPi0MCv3(e);
            factor /= FunctionNL_kTestBeamv3(e);
            testB->SetBinContent(iBin,factor);
        }
        testB->GetXaxis()->SetTitleOffset(0.8);
        testB->GetYaxis()->SetRangeUser(0.98,1.1);
        testB->DrawCopy("p");

        Color_t clr[10] = { kAzure, kRed, kAzure+4, kRed-4, kAzure+8, kRed+2, kAzure-3, kRed-2, kAzure-6, kRed-6};

        TH1D* testBarr[n];
        for(Int_t iNL=0; iNL<n; iNL++){
            testBarr[iNL] = (TH1D*) histoTotalCorrection[0]->Clone(Form("NL%i",iNL));
            testBarr[iNL]->Reset("ICE");
            DrawGammaSetMarker(testBarr[iNL], 21, 0.5, clr[iNL], clr[iNL]);
            testBarr[iNL]->SetLineWidth(3);
            for(Int_t iBin = 1; iBin <= testBarr[iNL]->GetNbinsX()+1; iBin++) {
                Float_t energy = testBarr[iNL]->GetXaxis()->GetBinCenter(iBin);
                Float_t factor = 1;
                if(iNL==0) factor /= (FunctionNL_DPOW(energy, 1.0443938253, -0.0691830812, -0.1247555443, 1.1673716264, -0.1853095466, -0.0848801702) - 0.0055);
                else if(iNL==1) factor /= (FunctionNL_DPOW(energy, 0.9980625418, -0.0564782662, -0.5, 1.0383412435, -0.0851830429, -0.4999999996) - 0.00175);
                else if(iNL==2) factor /= (FunctionNL_DPOW(energy, 1.1716155406, -0.1962930603, -0.0193959829, 1.0336659741, -0.0467778485, -0.4407662248) - 0.0055);
                else if(iNL==3) factor /= (FunctionNL_DPOW(energy, 1.0795372569, -0.1347324732, -0.1630736190, 1.1614181498, -0.199995361, -0.1711378093) - 0.0035);
                else if(iNL==4) factor /= (FunctionNL_DPOW(energy, 1.0166321784, -0.0440799552, -0.2611899222, 1.0636538464, -0.0816662488, -0.2173961316) - 0.007);
                else if(iNL==5) factor /= (FunctionNL_DPOW(energy, 1.0232969083, -0.090409434, -0.3592406513, 1.0383412435, -0.0851830429, -0.4999999996) + 0.0007);
                else if(iNL==6) factor /= (FunctionNL_DPOW(energy, 1.1100193881, -0.1389194936, -0.0800000242, 1.1673716264, -0.1853095466, -0.0848801702) - 0.017);
                else if(iNL==7) factor /= (FunctionNL_DPOW(energy, 1.0106037132, -0.0748250591, -0.4999999996, 1.0383412435, -0.0851830429, -0.4999999996) - 0.014);
                else if(iNL==8) factor /= (FunctionNL_DPOW(energy, 1.0520183153, -0.0806102847, -0.1450415920, 1.0336724056, -0.0467844121, -0.4406992764) - 0.016);
                else if(iNL==9) factor /= (FunctionNL_DPOW(energy, 1.0119417393, -0.0755250741, -0.4999999996, 1.1614181498, -0.1999995361, -0.1711378093) - 0.006);
                if(iNL>5 && iBin%10>=1) continue;
                testBarr[iNL]->SetBinContent(iBin,factor);
            }
            if(iNL==0) legend->AddEntry(testBarr[iNL],"Pythia8(LHC11a) - ConvCalo","l");
            else if(iNL==1) legend->AddEntry(testBarr[iNL],"Pythia8(LHC11a) - Calo","l");
            else if(iNL==2) legend->AddEntry(testBarr[iNL],"Pythia8(LHC13g) - ConvCalo","l");
            else if(iNL==3) legend->AddEntry(testBarr[iNL],"Pythia8(LHC13g) - Calo","l");
            else if(iNL==4) legend->AddEntry(testBarr[iNL],"Phojet - ConvCalo","l");
            else if(iNL==5) legend->AddEntry(testBarr[iNL],"Phojet - Calo","l");
            else if(iNL==6) legend->AddEntry(testBarr[iNL],"JetJet(LHC11a) - ConvCalo","p");
            else if(iNL==7) legend->AddEntry(testBarr[iNL],"JetJet(LHC11a) - Calo","p");
            else if(iNL==8) legend->AddEntry(testBarr[iNL],"JetJet(LHC13g) - ConvCalo","p");
            else if(iNL==9) legend->AddEntry(testBarr[iNL],"JetJet(LHC13g) - Calo","p");

            testBarr[iNL]->DrawCopy("p,same");
        }
        testB->DrawCopy("p,same");
        legend->AddEntry(testB,"kPi0MCv3 / kTestBeamv3","p");
        legend->Draw("same");
        PutProcessLabelAndEnergyOnPlot(0.7, 0.25, 0.03, "pp, #sqrt{#it{s}} = 2.76 TeV", "CCMF","");
        canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
        canvasRatio->SaveAs(Form("%s/TotalCorrection_Full_Case2_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
        canvasRatio->Clear();
        legend->Clear();
    }

    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;

    delete canvas;
    delete canvasRatio;

    return;
}
