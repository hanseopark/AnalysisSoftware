/****************************************************************************************************************************
******         provided by Gamma Conversion Group, PWG4,                                                     *****
******        Friederike Bock, friederike.bock@cern.ch                                                    *****
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ExtractSignalBinning.h"
// #include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"
#include "CommonHeaders/CombinationFunctions.h"

extern TRandom*    gRandom;
extern TBenchmark* gBenchmark;
extern TSystem*    gSystem;
extern TMinuit*    gMinuit;


// *** CONTENT
// *** settings
// *** read files: PCM
// *** plot mass & width 
//     - Pi0 PCM
//     - Eta PCM 
// *** plot example invariant mass bins 
//     - PCM
// *** save results   

//____________________________________________________________________________________________________________________________________________
void CombineMesonMeasurementsPbPb5TeV(  TString fileNamePCM         = "/home/meike/analysis/results/photonconvResults/PbPb/LHC15oLowIRp2-LHC15k1a1-LHC15k1a2-LHC15k1a3/11210013_00200009247602008250404000_0652501500000000/PbPb_5.02TeV/Pi0_data_GammaConvV1WithoutCorrection_11210013_00200009247602008250404000_0652501500000000.root",
                                        TString fileNamePCMEMCAL    = "", 
                                        TString fileNameEMCAL       = "",                                          
                                        TString fileNamePHOS        = "",  
                                        TString suffix              = "eps", 
                                        TString isMC                = "", 
                                        TString thesisPlots         = "", 
                                        TString bWCorrection        = "X",
                                        TString fileInputCorrFactors= ""
                                    ){

    TString date = ReturnDateString();
    cout << "date: " << date << endl;

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis();
    SetPlotStyle();
    
    gStyle->SetEndErrorSize(0);
    
    Bool_t plotMassAndWidth = kTRUE;
    Bool_t plotInvMassBins = kTRUE;
    Bool_t saveResults = kFALSE;

    Bool_t havePCM = kFALSE;
    Bool_t havePCMEMCAL = kFALSE;
    Bool_t havePHOS = kFALSE;
    Bool_t haveEMCAL = kFALSE;

    if(fileNamePCM.CompareTo("") != 0)         havePCM = kTRUE;
    if(fileNamePCMEMCAL.CompareTo("") != 0)    havePCMEMCAL = kTRUE;
    if(fileNameEMCAL.CompareTo("") != 0)       haveEMCAL = kTRUE;
    if(fileNamePHOS.CompareTo("") != 0)        havePHOS = kTRUE;

    TString dateForOutput                       = ReturnDateStringForOutput();
    cout << "date for output: " << dateForOutput << endl;

    TString centralityString;
    Int_t exampleBinPi0PCM;
    Int_t yMaxInvMassPi0PCM;
    if(fileNamePCM.Contains("10110013")) {
      centralityString =  "0-10%";
      exampleBinPi0PCM = 3;
      yMaxInvMassPi0PCM = 25000;
    }
    if(fileNamePCM.Contains("11210013")){
      centralityString =  "10-20%";
      exampleBinPi0PCM = 3;
      yMaxInvMassPi0PCM = 15000;
    }
    if(fileNamePCM.Contains("12510013")){
      centralityString =  "20-50%";
      exampleBinPi0PCM = 3;
      yMaxInvMassPi0PCM = 10000;
    }
    if(fileNamePCM.Contains("15910013")){
      centralityString =  "50-90%";
      exampleBinPi0PCM = 3;
      yMaxInvMassPi0PCM = 1000;
    }
    if(fileNamePCM.Contains("10010013")){
      centralityString =  "0-100%";
      exampleBinPi0PCM = 1;
      yMaxInvMassPi0PCM = 600000;
    }
    TString collisionSystemPbPb5TeV;
    collisionSystemPbPb5TeV.Form("pp, #sqrt{#it{s}_{NN}} = 5.02 TeV, %s", centralityString.Data());

    Double_t fBinsPi0HI5020GeVPt[14]  = { 0.0, 1.0, 1.4,
					  1.6, 1.8, 2.0, 2.2,
					  2.4, 2.6, 3.0, 3.5,
					  4.0, 5.0, 7.0};

    TString outputDir                           = Form("%s/%s/CombineMesonMeasurements5TeV%s",suffix.Data(),dateForOutput.Data(),bWCorrection.Data());
    cout << "output directory: "<< outputDir.Data() << endl;

    cout << "PCM file: "<< fileNamePCM.Data() << endl;

    gSystem->Exec("mkdir -p "+outputDir);
    if(havePCM)               gSystem->Exec(Form("cp %s %s/InputPCM.root", fileNamePCM.Data(), outputDir.Data()));
    if(havePCMEMCAL)          gSystem->Exec(Form("cp %s %s/InputPCMEMCAL.root", fileNamePCMEMCAL.Data(), outputDir.Data()));
    if(havePHOS)              gSystem->Exec(Form("cp %s %s/InputPHOS.root", fileNamePHOS.Data(), outputDir.Data()));
    if(haveEMCAL)             gSystem->Exec(Form("cp %s %s/InputEMCAL.root", fileNameEMCAL.Data(), outputDir.Data()));
    
    Bool_t thesis                               = kFALSE;
    if(thesisPlots.CompareTo("thesis") == 0){
        thesis                                  = kTRUE;
    }
    
    TString prefix2                             = "";
    if (isMC.CompareTo("kTRUE")==0){ 
        prefix2                                 = "MC";
    } else {    
        prefix2                                 = "Data";
    }
    
    Double_t mesonMassExpectPi0                 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    Double_t mesonMassExpectEta                 = TDatabasePDG::Instance()->GetParticle(221)->Mass();
    
    Int_t textSizeLabelsPixel                   = 0;
    Double_t textsizeLabelsPP                   = 0.07;
 
    Width_t  widthLinesBoxes                    = 1.4;
    Width_t  widthCommonFit                     = 2;
    
    // Definition of colors, styles and markers sizes
    Color_t  colorComb                          = kMagenta+2;
    Style_t  markerStyleComb                    = 20;
    Size_t   markerSizeComb                     = 2;
    
    Color_t  colorCombLowPt                     = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kFALSE);
    Color_t  colorCombHighPt                    = GetDefaultColorDiffDetectors("Comb", kFALSE, kFALSE, kTRUE);
    Style_t  markerStyleCombLowPt               = 20;
    Style_t  markerStyleCombHighPt              = 20;
    Size_t   markerSizeComparison               = 2;
    
    Color_t colorTrigg      [10]                = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};
    Color_t colorTriggShade [10]                = {kGray+1, kGray, kRed-6, kBlue-6, kGreen-8, kCyan-6, kViolet-8, kMagenta-8,  kRed-8, kBlue-8};
    Marker_t markerTrigg    [10]                = {20, 20, 21, 34, 29, 33, 21, 27, 28, 30 };
    Marker_t markerTriggMC  [10]                = {24, 24, 25, 28, 30, 27, 25, 27, 28, 30 };

    Size_t sizeTrigg        [10]                = {1.5, 1.5, 1.5, 2, 2.2, 2., 1.5, 2., 2.5, 1.5 };

    
    TString  nameMeasGlobal[11]                 = { "PCM", "PHOS", "EMCal", "PCM-PHOS", "PCM-EMCal",
                                                    "PCM-Dalitz", "PHOS-Dalitz", "EMCal-Dalitz", "EMCal high pT", "EMCal merged",
                                                    "PCMOtherDataset"};
    TString  nameTrigger[3]                     = {"MB", "EMC7", "EGA"};
    
    Color_t  colorDet[11];
    Color_t  colorDetMC[11];
    Style_t  markerStyleDet[11];
    Style_t  markerStyleDetMC[11];
    Size_t   markerSizeDet[11];
    Size_t   markerSizeDetMC[11];

    Color_t  colorNLO                           = kAzure-4;
    Style_t  styleMarkerNLOMuHalf               = 24;
    Style_t  styleMarkerNLOMuOne                = 27;
    Style_t  styleMarkerNLOMuTwo                = 30;
    Style_t  styleLineNLOMuHalf                 = 8;
    Style_t  styleLineNLOMuOne                  = 7;
    Style_t  styleLineNLOMuTwo                  = 4;
    Style_t  styleLineNLOMuTwoBKK               = 3;
    Style_t  styleLineNLOMuTwoDSS               = 6;
    Size_t   sizeMarkerNLO                      = 1;
    Width_t  widthLineNLO                       = 2.;

    
    for (Int_t i = 0; i < 11; i++){
        colorDet[i]                             = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kFALSE, kFALSE, kTRUE);
        colorDetMC[i]                           = GetDefaultColorDiffDetectors(nameMeasGlobal[i].Data(), kTRUE, kFALSE, kTRUE);
        markerStyleDet[i]                       = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kFALSE);
        markerStyleDetMC[i]                     = GetDefaultMarkerStyleDiffDetectors(nameMeasGlobal[i].Data(), kTRUE);
        markerSizeDet[i]                        = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kFALSE)*2;
        markerSizeDetMC[i]                      = GetDefaultMarkerSizeDiffDetectors(nameMeasGlobal[i].Data(), kTRUE)*2;
    }    
    

    // *******************************************************************************************************
    // ***************************** read data for PCM *******************************************************
    // *******************************************************************************************************
    Bool_t haveAllEtaInvMassPCM[3]                 = {kFALSE, kFALSE, kFALSE};
    Bool_t haveAllPi0InvMassPCM[3]                 = {kFALSE, kFALSE, kFALSE};
    if(havePCM){
      cout << "read PCM file" << endl;
      TFile* filePCM                                          = new TFile(fileNamePCM.Data());

      // *********
      // ** Pi0 **
      // *********

        TH1D* histoPCMPi0Mass                               = (TH1D*)filePCM->Get("histoMassMeson"); 
	histoPCMPi0Mass->Scale(1000); // [MeV]
        TH1D* histoPCMPi0FWHMMeV                            = (TH1D*)filePCM->Get("histoFWHMMeson");
	histoPCMPi0FWHMMeV->Scale(1000);
	histoPCMPi0FWHMMeV->Scale(1/2.35);  // sigma [MeV]

        TH1D* histoPi0InvMassSigPlusBGPCM[3];
        TH1D* histoPi0InvMassSigPCM[3];
        TH1D* histoPi0InvMassSigRemBGSubPCM[3];
        TH1D* histoPi0InvMassBGPCM[3];
        TH1D* histoPi0InvMassRemBGPCM[3];
        TH1D* histoPi0InvMassBGTotPCM[3];
        TF1* fitPi0InvMassSigPCM[3];
        TF1* fitPi0InvMassBGPCM[3];
        if (plotInvMassBins){
            for (Int_t i = 0; i < 1; i++){
	        histoPi0InvMassSigPCM[i]               = (TH1D*)filePCM->Get(Form("fHistoMappingSignalInvMass_in_Pt_Bin%02d",exampleBinPi0PCM));
                histoPi0InvMassSigPlusBGPCM[i]         = (TH1D*)filePCM->Get(Form("Mapping_GG_InvMass_in_Pt_Bin%02d",exampleBinPi0PCM));
                histoPi0InvMassBGPCM[i]                = (TH1D*)filePCM->Get(Form("Mapping_BckNorm_InvMass_in_Pt_Bin%02d",exampleBinPi0PCM));
                fitPi0InvMassSigPCM[i]                 = (TF1*)filePCM->Get(Form("Signal_InvMassFit_in_Pt_Bin%02d",exampleBinPi0PCM));

                if (histoPi0InvMassSigPCM[i] && histoPi0InvMassSigPlusBGPCM[i] && histoPi0InvMassBGPCM[i] && fitPi0InvMassSigPCM[i]){
                    haveAllPi0InvMassPCM[i]            = kTRUE;
                }


                if (haveAllPi0InvMassPCM[i]){
                    histoPi0InvMassSigPCM[i]->Fit(fitPi0InvMassSigPCM[i],"QRME0");
                    for (Int_t l=0; l < 6; l++){
                        cout << fitPi0InvMassSigPCM[i]->GetParameter(l) << "\t +- " << fitPi0InvMassSigPCM[i]->GetParError(l) << endl;
                    }
                    fitPi0InvMassBGPCM[i]                                  = new TF1("Linearpp","[0]+[1]*x",0.02,0.25);
                    fitPi0InvMassBGPCM[i]->SetParameter(0, fitPi0InvMassSigPCM[i]->GetParameter(4));
                    fitPi0InvMassBGPCM[i]->SetParameter(1, fitPi0InvMassSigPCM[i]->GetParameter(5));
                    TVirtualFitter * fitterPCM                             = TVirtualFitter::GetFitter();
                    Int_t nFreeParPCM                                      = fitPi0InvMassSigPCM[i]->GetNumberFreeParameters();
                    double * covMatrixPCM                                  = fitterPCM->GetCovarianceMatrix();

                    histoPi0InvMassRemBGPCM[i]                             = (TH1D*)histoPi0InvMassBGPCM[i]->Clone(Form("Pi0_InvMassRemBG_Example_%s",nameTrigger[i].Data()));
                    for (Int_t j = 1; j < histoPi0InvMassRemBGPCM[i]->GetNbinsX()+1; j++){
                        histoPi0InvMassRemBGPCM[i]->SetBinContent(j,0);
                        histoPi0InvMassRemBGPCM[i]->SetBinError(j,0);
                    }
                    for (Int_t j = histoPi0InvMassSigPCM[i]->GetXaxis()->FindBin(0.01); j < histoPi0InvMassSigPCM[i]->GetXaxis()->FindBin(0.30)+1; j++){
                        Double_t startBinEdge                                   = histoPi0InvMassSigPCM[i]->GetXaxis()->GetBinLowEdge(j);
                        Double_t endBinEdge                                     = histoPi0InvMassSigPCM[i]->GetXaxis()->GetBinUpEdge(j);
                        Double_t intLinearBack                                  = fitPi0InvMassBGPCM[i]->Integral(startBinEdge, endBinEdge)/(endBinEdge-startBinEdge) ;
                        Double_t errorLinearBck                                 = pow(( pow( (endBinEdge-startBinEdge)*fitPi0InvMassSigPCM[i]->GetParError(4),2) +
                                                                                        pow(0.5*(endBinEdge*endBinEdge-startBinEdge*startBinEdge)*fitPi0InvMassSigPCM[i]->GetParError(5),2)
                                                                                        +2*covMatrixPCM[nFreeParPCM*nFreeParPCM-2]*(endBinEdge-startBinEdge)*0.5*
                                                                                        (endBinEdge*endBinEdge-startBinEdge*startBinEdge)),0.5)/(endBinEdge-startBinEdge);
                        histoPi0InvMassRemBGPCM[i]->SetBinContent(j,intLinearBack);
                        histoPi0InvMassRemBGPCM[i]->SetBinError(j,errorLinearBck);
                    }
                    histoPi0InvMassBGTotPCM[i]         = (TH1D*)histoPi0InvMassBGPCM[i]->Clone(Form("Pi0_InvMassTotBG_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassBGTotPCM[i]->Sumw2();
                    histoPi0InvMassBGTotPCM[i]->Add(histoPi0InvMassRemBGPCM[i]);
                    histoPi0InvMassSigRemBGSubPCM[i]   = (TH1D*)histoPi0InvMassSigPCM[i]->Clone(Form("Pi0_InvMassSigRemBGSub_Example_%s",nameTrigger[i].Data()));
                    histoPi0InvMassSigRemBGSubPCM[i]->Sumw2();
                    histoPi0InvMassSigRemBGSubPCM[i]->Add(histoPi0InvMassRemBGPCM[i],-1);
                    fitPi0InvMassSigPCM[i]->SetParameter(4, 0);
                    fitPi0InvMassSigPCM[i]->SetParameter(5, 0);
                }
                cout << nameTrigger[i].Data() << "\t" << histoPi0InvMassSigPCM[i] << "\t" << histoPi0InvMassSigPlusBGPCM[i] << "\t" << histoPi0InvMassBGPCM[i] << "\t" << fitPi0InvMassSigPCM[i]
                      << "\t" << haveAllPi0InvMassPCM[i] << endl;
            }
        }
    
    }



    // **********************************************************************************************************************
    // **************************************** Mass and width for pi0  *****************************************************
    // **********************************************************************************************************************
    if(plotMassAndWidth){

      Double_t pTMaxPi0PCM = 12.;
      Double_t pTMinPi0PCM = 0.63;

      Double_t arrayBoundariesX1_4[2];
      Double_t arrayBoundariesY1_4[3];
      Double_t relativeMarginsX[3];
      Double_t relativeMarginsY[3];
      textSizeLabelsPixel             = 50;
      ReturnCorrectValuesForCanvasScaling(1350,1250, 1, 2,0.09, 0.005, 0.005,0.085,arrayBoundariesX1_4,arrayBoundariesY1_4,relativeMarginsX,relativeMarginsY);

      TCanvas* canvasMassWidthPi0     = new TCanvas("canvasMassWidthPi0","",0,0,1350,1250);  // gives the page size
      DrawGammaCanvasSettings( canvasMassWidthPi0,  0.13, 0.02, 0.03, 0.06);

      TPad* padWidthPi0               = new TPad("padWidthPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[1], arrayBoundariesX1_4[1], arrayBoundariesY1_4[0],-1, -1, -2);
      DrawGammaPadSettings( padWidthPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[0], relativeMarginsY[1]);
      padWidthPi0->Draw();
      
      TPad* padMassPi0                = new TPad("padMassPi0", "", arrayBoundariesX1_4[0], arrayBoundariesY1_4[2], arrayBoundariesX1_4[1], arrayBoundariesY1_4[1],-1, -1, -2);
      DrawGammaPadSettings( padMassPi0, relativeMarginsX[0], relativeMarginsX[2], relativeMarginsY[1], relativeMarginsY[2]);
      padMassPi0->Draw();
    
      TPad* padMassLegend1            = new TPad("padMassLegend1", "", 0.13, 0.36, 0.52, 0.52,-1, -1, -2);
      DrawGammaPadSettings( padMassLegend1, 0., 0., 0., 0.);
      padMassLegend1->SetFillStyle(0);
      padMassLegend1->Draw();
    
      padWidthPi0->cd();
      padWidthPi0->SetLogx(); 

      Double_t margin                 = relativeMarginsX[0]*2.7*1350;
      Double_t textsizeLabelsWidth    = 0;
      Double_t textsizeFacWidth       = 0;
      if (padWidthPi0->XtoPixel(padWidthPi0->GetX2()) < padWidthPi0->YtoPixel(padWidthPi0->GetY1())){
	textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
	textsizeFacWidth            = (Double_t)1./padWidthPi0->XtoPixel(padWidthPi0->GetX2()) ;
      } else {
	textsizeLabelsWidth         = (Double_t)textSizeLabelsPixel/padWidthPi0->YtoPixel(padWidthPi0->GetY1());
	textsizeFacWidth            = (Double_t)1./padWidthPi0->YtoPixel(padWidthPi0->GetY1());
      }
      cout << textsizeLabelsWidth << endl;
        
      TH2F * histo2DAllPi0FWHM    = new TH2F("histo2DAllPi0FWHM","histo2DAllPi0FWHM", 20, pTMinPi0PCM, pTMaxPi0PCM ,1000., -30, 40);
      SetStyleHistoTH2ForGraphs(histo2DAllPi0FWHM, "#it{p}_{T} (GeV/#it{c})", "Peak width (MeV/#it{c}^{2})", 0.85*textsizeLabelsWidth, textsizeLabelsWidth,
				0.85*textsizeLabelsWidth, textsizeLabelsWidth, 0.8,0.28/(textsizeFacWidth*margin), 512, 505);
      histo2DAllPi0FWHM->GetYaxis()->SetRangeUser(-0.5.,7.5); 
      histo2DAllPi0FWHM->GetYaxis()->SetMoreLogLabels(kTRUE);
      histo2DAllPi0FWHM->GetYaxis()->SetNdivisions(505);
      histo2DAllPi0FWHM->GetYaxis()->SetNoExponent(kTRUE);
      histo2DAllPi0FWHM->GetXaxis()->SetTickLength(0.05);
      histo2DAllPi0FWHM->GetYaxis()->SetTickLength(0.026);
      histo2DAllPi0FWHM->DrawCopy(); 

      DrawGammaSetMarker(histoPCMPi0FWHMMeV, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
      histoPCMPi0FWHMMeV->Draw("p,same,e");
      //DrawGammaSetMarker(histoPCMPi0TrueFWHMMeV, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
      //histoPCMPi0TrueFWHMMeV->Draw("p,same,e");

      TLatex *labelLegendAMass    = new TLatex(0.13,0.06,"a)");
      SetStyleTLatex( labelLegendAMass, textSizeLabelsPixel,4);
      labelLegendAMass->SetTextFont(43);
      labelLegendAMass->Draw();

      TLatex *labelMassPerf       = new TLatex(0.13,0.87,Form("ALICE performance %s",date.Data()));
      SetStyleTLatex( labelMassPerf, textSizeLabelsPixel,4);
      labelMassPerf->SetTextFont(43);
      labelMassPerf->Draw();        
      TLatex *labelMassEnergy     = new TLatex(0.13,0.78,collisionSystemPbPb5TeV.Data());
      SetStyleTLatex( labelMassEnergy, textSizeLabelsPixel,4);
      labelMassEnergy->SetTextFont(43);
      labelMassEnergy->Draw();
      TLatex *labelMassPi0        = new TLatex(0.13,0.69,"#pi^{0} #rightarrow #gamma#gamma");
      SetStyleTLatex( labelMassPi0, textSizeLabelsPixel,4);
      labelMassPi0->SetTextFont(43);
      labelMassPi0->Draw();  
      
      padMassPi0->cd();
      padMassPi0->SetLogx();

      Double_t textsizeLabelsMass         = 0;
      Double_t textsizeFacMass            = 0;
      if (padMassPi0->XtoPixel(padMassPi0->GetX2()) <padMassPi0->YtoPixel(padMassPi0->GetY1()) ){
	textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
	textsizeFacMass                 = (Double_t)1./padMassPi0->XtoPixel(padMassPi0->GetX2()) ;
      } else {
	textsizeLabelsMass              = (Double_t)textSizeLabelsPixel/padMassPi0->YtoPixel(padMassPi0->GetY1());
	textsizeFacMass                 = (Double_t)1./padMassPi0->YtoPixel(padMassPi0->GetY1());
      }
        
      TH2F * histo2DAllPi0Mass            = new TH2F("histo2DAllPi0Mass","histo2DAllPi0Mass",20, pTMinPi0PCM, pTMaxPi0PCM , 1000., 120., 175);
      SetStyleHistoTH2ForGraphs(histo2DAllPi0Mass, "#it{p}_{T} (GeV/#it{c})", "Peak position (MeV/#it{c}^{2})", 0.85*textsizeLabelsMass, textsizeLabelsMass, 0.85*textsizeLabelsMass, 
				textsizeLabelsMass, 0.9, 0.28/(textsizeFacMass*margin), 512, 505);
      histo2DAllPi0Mass->GetYaxis()->SetRangeUser(133.,138.); 
      histo2DAllPi0Mass->GetXaxis()->SetMoreLogLabels(kTRUE);
      histo2DAllPi0Mass->GetYaxis()->SetNdivisions(505);
      histo2DAllPi0Mass->GetYaxis()->SetNoExponent(kTRUE);
      histo2DAllPi0Mass->GetXaxis()->SetTickLength(0.05);
      histo2DAllPi0Mass->GetXaxis()->SetLabelOffset(-0.015);
      histo2DAllPi0Mass->DrawCopy(); 

      DrawGammaSetMarker(histoPCMPi0Mass, markerStyleDet[0], markerSizeDet[0]*0.55, colorDet[0] , colorDet[0]);
      histoPCMPi0Mass->Draw("p,same,e");
      //DrawGammaSetMarker(histoPCMPi0TrueMass, markerStyleDetMC[0], markerSizeDetMC[0]*0.55, colorDetMC[0] , colorDetMC[0]);
      //histoPCMPi0TrueMass->Draw("p,same,e");
      DrawGammaLines(pTMinPi0PCM, pTMaxPi0PCM , mesonMassExpectPi0*1000., mesonMassExpectPi0*1000.,0.1, kGray);

      TLatex *labelLegendBMass            = new TLatex(0.13,0.22,"b)");
      SetStyleTLatex( labelLegendBMass, textSizeLabelsPixel,4);
      labelLegendBMass->SetTextFont(43);
      labelLegendBMass->Draw();
        
      //********************************** Defintion of the Legend **************************************************    
      Double_t columnsLegendMass2[3]      = {0.,0.57,0.84};
      Double_t rowsLegendMass2[4]         = {0.75,0.5,0.25,0.01};
      //******************* Offsets ***********************
      Double_t offsetMarkerXMass2         = 0.1;
      Double_t offsetMarkerYMass2         = 0.1;
      //****************** Scale factors ******************
      Double_t scaleMarkerMass2           = 1.2;
        
      padMassLegend1->cd();
      //****************** first Column **************************************************
      
      TLatex *textMassPCM                 = new TLatex(columnsLegendMass2[0],rowsLegendMass2[1],"PCM");
      SetStyleTLatex( textMassPCM, textSizeLabelsPixel,4);
      textMassPCM->SetTextFont(43);
      textMassPCM->Draw();
    
      //****************** second Column *************************************************
      TLatex *textMassData                = new TLatex(columnsLegendMass2[1],rowsLegendMass2[0] ,"Data");
      SetStyleTLatex( textMassData, textSizeLabelsPixel,4);
      textMassData->SetTextFont(43);
      textMassData->Draw();
      //TLatex *textMassMC                  = new TLatex(columnsLegendMass2[2] ,rowsLegendMass2[0],"MC");
      //SetStyleTLatex( textMassMC, textSizeLabelsPixel,4);
      //textMassMC->SetTextFont(43);
      //textMassMC->Draw();
        

      TMarker* markerPCMPi0Mass        = CreateMarkerFromHisto(histoPCMPi0Mass,columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
      markerPCMPi0Mass->DrawMarker(columnsLegendMass2[1]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2);
      //TMarker* markerPCMPi0MassMC      = CreateMarkerFromHisto(histoPCMPi0TrueMass,columnsLegendMass2[2]+ offsetMarkerXMass2 ,rowsLegendMass2[1]+ offsetMarkerYMass2 ,scaleMarkerMass2);
      //markerPCMPi0MassMC->DrawMarker(columnsLegendMass2[2]+ offsetMarkerXMass2-0.04 ,rowsLegendMass2[1]+ offsetMarkerYMass2);

      canvasMassWidthPi0->Update();
      canvasMassWidthPi0->Print(Form("%s/Pi0_MassAndWidth_%s.%s",outputDir.Data(), centralityString.Data(), suffix.Data()));


    }
    // **********************************************************************************************************************
    // **************************Plot example invariant mass bins ***********************************************************
    // **********************************************************************************************************************
    textSizeLabelsPixel                 = 100*3/5;
    TCanvas* canvasInvMassSamplePlot    = new TCanvas("canvasInvMassSamplePlot","",0,0,1500,1500);  // gives the page size
    DrawGammaCanvasSettings( canvasInvMassSamplePlot,  0.09, 0.01, 0.035, 0.08);

    Style_t markerStyleInvMassSGBG      = 0;
    Size_t markerSizeInvMassSGBG        = 0;
    Color_t markerColorInvMassSGBG      = kBlack;
    Style_t markerStyleInvMassMBG       = 24;
    Size_t markerSizeInvMassMBG         = 1.5;
    Color_t markerColorInvMassMBG       = kGray+2;
    Style_t markerStyleInvMassBG        = 20;
    Size_t markerSizeInvMassBG          = 2;
    Color_t markerColorInvMassBG        = kBlack;
    Style_t markerStyleInvMassSG        = 20;
    Size_t markerSizeInvMassSG          = 3;
    Color_t markerColorInvMassSG        = kRed+2;
    Color_t fitColorInvMassSG           = kAzure+2;

    Double_t marginInvMass          = 0.1*1500;
    Double_t textsizeLabelsInvMass  = 0;
    Double_t textsizeFacInvMass     = 0;

    Double_t InvMassMinPi0 = 0.05;
    Double_t InvMassMaxPi0 = 0.23;
    Double_t InvMassMinEta = 0.41;  
    Double_t InvMassMaxEta = 0.71;

    if (canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) < canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1())){
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->XtoPixel(canvasInvMassSamplePlot->GetX2()) ;
    } else {
        textsizeLabelsInvMass       = (Double_t)textSizeLabelsPixel/canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
        textsizeFacInvMass          = (Double_t)1./canvasInvMassSamplePlot->YtoPixel(canvasInvMassSamplePlot->GetY1());
    }
    cout << textsizeLabelsInvMass << endl;

    TH2F * histo2DPi0InvMassDummy;
    histo2DPi0InvMassDummy             = new TH2F("histo2DPi0InvMassDummy","histo2DPi0InvMassDummy",11000,0.02,0.255,21000,-1000,600000);
    SetStyleHistoTH2ForGraphs(histo2DPi0InvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                            0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

    TH2F * histo2DEtaInvMassDummy;
    histo2DEtaInvMassDummy             = new TH2F("histo2DEtaInvMassDummy","histo2DEtaInvMassDummy",11000,0.35,0.695,21000,-1000,600000);
    SetStyleHistoTH2ForGraphs(histo2DEtaInvMassDummy, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})","Counts",0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,
                            0.85*textsizeLabelsInvMass, textsizeLabelsInvMass,0.88, 0.115/(textsizeFacInvMass*marginInvMass));

    if (plotInvMassBins){
        for (Int_t i =0 ; i < 3; i++){
	  /******************************************************************************************/
	  /*********************************** Pi0 PCM **********************************************/
	  /******************************************************************************************/
            if (haveAllPi0InvMassPCM[i]){
                canvasInvMassSamplePlot->cd();
                histo2DPi0InvMassDummy->GetXaxis()->SetRangeUser(InvMassMinPi0,InvMassMaxPi0);
                histo2DPi0InvMassDummy->GetYaxis()->SetRangeUser(histoPi0InvMassSigRemBGSubPCM[i]->GetMinimum(), yMaxInvMassPi0PCM);
                histo2DPi0InvMassDummy->DrawCopy();

		TLatex *labelInvMassPtRangePi0PCM = new TLatex(0.945,0.9,Form("#pi^{0}: %0.1f GeV/#it{c} < #it{p}_{T} < %0.1f GeV/#it{c}", fBinsPi0HI5020GeVPt[exampleBinPi0PCM], fBinsPi0HI5020GeVPt[exampleBinPi0PCM+1]));

                DrawGammaSetMarker(histoPi0InvMassSigPlusBGPCM[i], markerStyleInvMassSGBG, markerSizeInvMassSGBG, markerColorInvMassSGBG, markerColorInvMassSGBG);
                histoPi0InvMassSigPlusBGPCM[i]->SetLineWidth(1);
                histoPi0InvMassSigPlusBGPCM[i]->Draw("hist,e,same");
                DrawGammaSetMarker(histoPi0InvMassBGTotPCM[i], markerStyleInvMassMBG, markerSizeInvMassMBG, markerColorInvMassMBG, markerColorInvMassMBG);
                histoPi0InvMassBGTotPCM[i]->Draw("same");

                DrawGammaSetMarker(histoPi0InvMassSigRemBGSubPCM[i], markerStyleInvMassSG, markerSizeInvMassSG, markerColorInvMassSG, markerColorInvMassSG);
                histoPi0InvMassSigRemBGSubPCM[i]->Draw("same");
                fitPi0InvMassSigPCM[i]->SetNpx(1000);
                fitPi0InvMassSigPCM[i]->SetRange(0,0.255);
                fitPi0InvMassSigPCM[i]->SetLineWidth(2);
                fitPi0InvMassSigPCM[i]->SetLineColor(fitColorInvMassSG);
                fitPi0InvMassSigPCM[i]->Draw("same");

                //
                TLatex *labelInvMassPerf      = new TLatex(0.135,0.9,"ALICE performance");
                SetStyleTLatex( labelInvMassPerf, 0.85*textSizeLabelsPixel,4);
                labelInvMassPerf->SetTextFont(43);
                labelInvMassPerf->Draw();

                TLatex *labelDate     = new TLatex(0.135,0.9-1*0.8*textsizeLabelsPP,date.Data());
                SetStyleTLatex( labelDate, 0.85*textSizeLabelsPixel,4);
                labelDate->SetTextFont(43);
                labelDate->Draw();

                TLatex *labelInvMassEnergy      = new TLatex(0.135,0.9-2*0.8*textsizeLabelsPP,collisionSystemPbPb5TeV.Data());
                SetStyleTLatex( labelInvMassEnergy, 0.85*textSizeLabelsPixel,4);
                labelInvMassEnergy->SetTextFont(43);
                labelInvMassEnergy->Draw();

                TLatex *labelInvMassTrigger      = new TLatex(0.135,0.9-3*0.8*textsizeLabelsPP,Form("%s triggered",nameTrigger[i].Data()));
                SetStyleTLatex( labelInvMassTrigger, 0.85*textSizeLabelsPixel,4);
                labelInvMassTrigger->SetTextFont(43);
                labelInvMassTrigger->Draw();

                TLatex *labelInvMassRecoPCM  = new TLatex(0.135,0.9-4*0.8*textsizeLabelsPP,"PCM");
                SetStyleTLatex( labelInvMassRecoPCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassRecoPCM->SetTextFont(43);
                labelInvMassRecoPCM->Draw();

                SetStyleTLatex( labelInvMassPtRangePi0PCM, 0.85*textSizeLabelsPixel,4);
                labelInvMassPtRangePi0PCM->SetTextAlign(31);
                labelInvMassPtRangePi0PCM->SetTextFont(43);
                labelInvMassPtRangePi0PCM->Draw();

                TLegend* legendInvMassPi0PCM  = GetAndSetLegend2(0.67, 0.88-5*0.75*textsizeLabelsPP, 0.9, 0.88, 0.85*textSizeLabelsPixel);
                legendInvMassPi0PCM->SetMargin(0.25);
                legendInvMassPi0PCM->AddEntry(histoPi0InvMassSigPlusBGPCM[i],"Raw real events","l");
                legendInvMassPi0PCM->AddEntry(histoPi0InvMassBGTotPCM[i],"Mixed event +","p");
                legendInvMassPi0PCM->AddEntry((TObject*)0,"corr. BG","");
                legendInvMassPi0PCM->AddEntry(histoPi0InvMassSigRemBGSubPCM[i],"BG subtracted","p");
                legendInvMassPi0PCM->AddEntry(fitPi0InvMassSigPCM[i], "Fit","l");
                legendInvMassPi0PCM->Draw();
                canvasInvMassSamplePlot->SaveAs(Form("%s/Pi0_InvMassBinPCM_%s_%s.%s",outputDir.Data(), nameTrigger[i].Data(), centralityString.Data(), suffix.Data()));
            } else {
                cout << "missing partial input for Pi0 invariant mass bin for PCM for trigger: " << nameTrigger[i].Data() << endl;
            }


	}
   }


     
 // **********************************************************************************************************************
 // ************************* Saving of final results ********************************************************************
 // **********************************************************************************************************************
 
    if(saveResults){
    TString nameOutputCommonFile    = Form("CombinedResultsPbPb5TeV_%s.root", dateForOutput.Data());    
    TFile fCombResults(nameOutputCommonFile.Data(), "RECREATE");

    fCombResults.mkdir("Pi0PbPb5TeV");
    TDirectoryFile* directoryPi0 = (TDirectoryFile*)fCombResults.Get("Pi0PbPb5TeV");
    fCombResults.cd("Pi0PbPb5TeV");
 
    histoPCMPi0Mass->Write("Pi0MassDataPCM");
    histoPCMPi0TrueMass->Write("Pi0MassMCPCM");
    histoPHOSPi0Mass->Write("Pi0MassDataPHOS");
    histoPHOSPi0TrueMass->Write("Pi0MassMCPHOS");
    graphEMCALPi0Mass->Write("Pi0MassDataEMCAL");
    graphEMCALPi0MassMC->Write("Pi0MassMCEMCAL");
    graphPCMEMCALPi0Mass->Write("Pi0MassDataPCMEMCAL");
    graphPCMEMCALPi0MassMC->Write("Pi0MassMCPCMEMCAL");

    histoPCMPi0FWHMMeV->Write("Pi0WidthDataPCM");
    histoPCMPi0TrueFWHMMeV->Write("Pi0WidthMCPCM");
    histoPHOSPi0FWHMMeV->Write("Pi0WidthDataPHOS");
    histoPHOSPi0TrueFWHMMeV->Write("Pi0WidthMCPHOS");
    graphEMCALPi0FWHM->Write("Pi0WidthDataEMCAL");
    graphEMCALPi0FWHMMC->Write("Pi0WidthMCEMCAL");
    graphPCMEMCALPi0FWHM->Write("Pi0WidthDataPCMEMCAL");
    graphPCMEMCALPi0FWHMMC->Write("Pi0WidthMCPCMEMCAL");
            
    fCombResults.mkdir("EtaPbPb5TeV");
    TDirectoryFile* directoryEta = (TDirectoryFile*)fCombResults.Get("EtaPbPb5TeV");
    fCombResults.cd("EtaPbPb5TeV");
              
    histoPCMEtaMass->Write("EtaMassDataPCM");
    histoPCMEtaTrueMass->Write("EtaMassMCPCM");
    graphEMCALEtaMass->Write("EtaMassDataEMCAL");
    graphEMCALEtaMassMC->Write("EtaMassMCEMCAL");
    graphPCMEMCALEtaMass->Write("EtaMassDataPCMEMCAL");
    graphPCMEMCALEtaMassMC->Write("EtaMassMCPCMEMCAL");
    
    histoPCMEtaFWHMMeV->Write("EtaWidthDataPCM");
    histoPCMEtaTrueFWHMMeV->Write("EtaWidthMCPCM");
    graphEMCALEtaFWHM->Write("EtaWidthDataEMCAL");
    graphEMCALEtaFWHMMC->Write("EtaWidthMCEMCAL");
    graphPCMEMCALEtaFWHM->Write("EtaWidthDataPCMEMCAL");
    graphPCMEMCALEtaFWHMMC->Write("EtaWidthMCPCMEMCAL");
        
    fCombResults.Close();

    }

}
