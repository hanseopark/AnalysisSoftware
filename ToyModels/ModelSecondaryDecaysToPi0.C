/****************************************************************************************************************************
******        provided by Gamma Conversion Group, PWGGA,                                                                *****
******        Friederike Bock, friederike.bock@cern.ch                                                                  *****
******        basics derived from Jason Kamin jason.kamin@cern.ch                                                       *****
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
#include "TRandom3.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGenPhaseSpace.h"
#include "TSpline.h"
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/FittingGammaConversion.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/ConversionFunctions.h"

//****************************************************************************
//*************** Draw 2D histogram ******************************************
//****************************************************************************
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
//******************* return minimum for 1 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestBin1DHist(TH1* hist, Double_t maxStart = 1e6 ){
    Double_t smallesContent     = maxStart;
    for (Int_t i= 0; i < hist->GetNbinsX(); i++){
        if (hist->GetBinContent(i) != 0 && smallesContent > hist->GetBinContent(i)){
            smallesContent = hist->GetBinContent(i);
        }    
    }    
    return smallesContent;
}

//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestBin2DHist(TH2* histo, Double_t maxStart = 1e6 ){
    Double_t minimum = maxStart;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        for (Int_t j = 1; j<histo->GetNbinsY(); j++){
            if (histo->GetBinContent(i,j) < minimum && histo->GetBinContent(i,j) > 0){
                minimum = histo->GetBinContent(i,j);
            }
        }
    }
    return minimum;
}


//****************************************************************************
//**************** Secondary pion decay simulation ***************************
//****************************************************************************
void ModelSecondaryDecaysToPi0(     Int_t nEvts             = 1000000, 
                                    Int_t particle          = 0,                // 0 = K0s, 1 = K0l, 2 = Lambda
                                    TString energy          = "2.76TeV",
                                    Double_t minPt          = 2,
                                    Double_t maxPt          = 50,
                                    TString filename        = "",
                                    TString suffix          = "eps",
                                    TString cutSelection    = "",
                                    Int_t mode              = 0,
                                    Bool_t kIsMC            = 0,
                                    TString filenameMC      = ""
                               ){

    //*************************************************************************************************
    //******************************** Set Style settings globally ************************************
    //*************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis(suffix);  
    SetPlotStyle();
    
    TString fCollisionSystenWrite       = ReturnCollisionEnergyOutputString(energy);
    
    TString fEventCutSelection          = "";
    TString fGammaCutSelection          = "";
    TString fClusterCutSelection        = "";
    TString fElectronCutSelection       = "";
    TString fMesonCutSelection          = "";
    ReturnSeparatedCutNumberAdvanced(cutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    
    Double_t fYMaxMother                = 2;  
    TString centralityString            = GetCentralityString(fEventCutSelection);
    cout << "centrality : "<< centralityString.Data() << endl;
    
    //*************************************************************************************************
    //*************************** Initialize pt spectra ***********************************************
    //*************************************************************************************************
    Double_t massParticle       = 0;
    
    TString outputlabel         = "";
    TString plotLabel           = "";
    TString daughterLabels[10]  = {"","","","","","","","","",""};
    Color_t colors[10]          = { kRed+1, kGreen+2, kOrange+7, kBlue+2, kAzure, 
                                    kYellow-8, kViolet, kCyan+2, kGray+1, kPink+2 };
    
    Int_t nDaughters            = 0;
    Double_t* masses            = NULL;
    Int_t* pdgCodesDaughters    = NULL;
    Double_t branchRatio        = 1;
    
    Double_t yRanges[10]        = { 0.1,    0.3,    0.5,    0.55,   0.6,
                                    0.65,   0.7,    0.75,   0.8,    0.85 };
    Double_t minDalitz[2]       = {0,0};
    Double_t maxDalitz[2]       = {2,2};
                                    
                                    
    // K0s
    if (particle == 0){
        massParticle                = TDatabasePDG::Instance()->GetParticle(310)->Mass();
        outputlabel                 = "K0S";
        plotLabel                   = "K^{0}_{s}";
        daughterLabels[0]           = "#pi^{0}_{1}";
        daughterLabels[1]           = "#pi^{0}_{2}";
        nDaughters                  = 2;
        pdgCodesDaughters           = new Int_t[nDaughters];
        pdgCodesDaughters[0]        = 111;
        pdgCodesDaughters[1]        = 111;
        masses                      = new Double_t[nDaughters];
        masses[0]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();;
        masses[1]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();;
        branchRatio                 = 0.3069;
    // K0l    
    } else if (particle == 1 ){   
        massParticle                = TDatabasePDG::Instance()->GetParticle(130)->Mass();
        outputlabel                 = "K0L";
        plotLabel                   = "K^{0}_{l}";
        daughterLabels[0]           = "#pi^{0}_{1}";
        daughterLabels[1]           = "#pi^{0}_{2}";
        daughterLabels[2]           = "#pi^{0}_{3}";
        nDaughters                  = 3;
        pdgCodesDaughters           = new Int_t[nDaughters];
        pdgCodesDaughters[0]        = 111;
        pdgCodesDaughters[1]        = 111;
        pdgCodesDaughters[2]        = 111;
        masses                      = new Double_t[nDaughters];
        masses[0]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();
        masses[1]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();
        masses[2]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[2])->Mass();
        branchRatio                 = (0.1952+1/3*0.1254); // not all of these are seen in the experiment, correction for that done in the task
        minDalitz[0]                = (masses[0]+masses[1])*0.8;
        minDalitz[1]                = (masses[0]+masses[2])*0.8;
        maxDalitz[0]                = (massParticle-masses[2])*1.2;
        maxDalitz[1]                = (massParticle-masses[1])*1.2;
        // Lambda    
    } else if (particle == 2 ){   
        massParticle                = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
        outputlabel                 = "Lambda";
        plotLabel                   = "#Lambda";
        daughterLabels[0]           = "n";
        daughterLabels[1]           = "#pi^{0}";
        nDaughters                  = 2;
        pdgCodesDaughters           = new Int_t[nDaughters];
        pdgCodesDaughters[0]        = 2112;
        pdgCodesDaughters[1]        = 111;
        masses                      = new Double_t[nDaughters];
        masses[0]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();;
        masses[1]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();;
        branchRatio                 = 0.3580;
    }
    // define output file name
    TString nameOutputRootFile          = Form("ToyMCOutput_%s_%s_%dMioEvt_%dMeV_%dMeV.root",outputlabel.Data(), fCollisionSystenWrite.Data(),(Int_t)(nEvts/1e6), (Int_t)(minPt*1000), (Int_t)(maxPt*1000));
    if (kIsMC) nameOutputRootFile       = Form("ToyMCOutput_MC_%s_%s_%dMioEvt_%dMeV_%dMeV.root",outputlabel.Data(), fCollisionSystenWrite.Data(),(Int_t)(nEvts/1e6), (Int_t)(minPt*1000), (Int_t)(maxPt*1000));
    cout << "searching for: " << nameOutputRootFile.Data() << endl;
    TFile* exisitingFile            = new TFile(nameOutputRootFile);
    if (!exisitingFile->IsZombie()){
        cout << "INFO: the file for this energy and rapidity range: " << nameOutputRootFile.Data() << "\t has been generated already" << endl;
        cout << "If you want to rerun the generation, please delete the file" << endl;
        cout << "exiting here!" << endl;
        
        cout << "adding output file name to ToyMCOutputs.txt" << endl;
        gSystem->Exec("echo '"+nameOutputRootFile+"' >> ToyMCOutputs.txt");
        return;
    }    
    delete exisitingFile; 

    // create output directory
    TString fOutputDir                  = Form("%s/%s_%dMioEvt_%dMeV_%dMeV/",fCollisionSystenWrite.Data(), outputlabel.Data(), (Int_t)(nEvts/1e6), (Int_t)(minPt*1000), (Int_t)(maxPt*1000));
    if (kIsMC) fOutputDir               = Form("%sMC/%s_%dMioEvt_%dMeV_%dMeV/",fCollisionSystenWrite.Data(), outputlabel.Data(), (Int_t)(nEvts/1e6), (Int_t)(minPt*1000), (Int_t)(maxPt*1000));
    gSystem->Exec("mkdir -p "+fOutputDir);


    //*************************************************************************************************
    //************************* Global variables for particle spectra from input **********************
    //*************************************************************************************************
    TFile* inputFile            = new TFile(filename.Data());
    TFile* inputFileMC          = NULL; 
    if (kIsMC){
        inputFileMC             = new TFile(filenameMC.Data());
    }    
    TH1D* partSpectrum          = NULL; 
    TH1D* partSpectrumAlter     = NULL; 
    TH1D* chHadDNdEta           = NULL;
    TF1* ptDistribution         = NULL; 
    TSpline3* partParam         = NULL; 
    TH1D* ratioPartSpecToFit    = NULL; 
    TH1D* ratioPartSpecToSpline = NULL; 
    TH1D* ratioPartSpecAlterToSpline        = NULL; 
    TH1D* ptXCheckDecayProds    = NULL;
    Bool_t haveMCxCheckInRap[10]= { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, 
                                    kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    
    TGraphErrors* partSpectrumExt           = NULL; 
    TSpline3* partParamExt                  = NULL; 
    TGraphErrors* ratioPartSpecExtToSpline  = NULL; 

    TH1F* partSpectrumFromParam     = new TH1F("partSpectrumFromParam","partSpectrumFromParam",7000,0,70); 
    TH1F* partSpectrumFromParamRest = new TH1F("partSpectrumFromParamRest","partSpectrumFromParamRest",7000,0,70); 
    // Data inputs
    TCanvas *canvasFitQA = new TCanvas("canvasFitQA","canvasFitQA",1000,800);
    DrawGammaCanvasSettings( canvasFitQA, 0.08, 0.015, 0.02, 0.08);
    canvasFitQA->cd();
    canvasFitQA->SetLogy(1);

    // careful what the input contains: 
    // if input is dN/dydpt 1/2pi 1/pt needs to be multiplied with pt
    // if input is dN/dydpt 1/2pi no need to be multiplied with pt
    
    // K0s & K0L
    if (particle == 0 || particle == 1 ){
        if (energy.CompareTo("2.76TeV") == 0){
            if (!kIsMC){
                partSpectrum            = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat2760GeV");
                partSpectrumAlter       = (TH1D*)inputFile->Get("histoNeutralKaonSpecStat2760GeV");
            } else {
                if (particle == 0){
                    partSpectrum            = (TH1D*)inputFileMC->Get("fHistoMCK0sPtRebinned2");
                    for (Int_t y=0; y< 10; y++){
                        if (!ptXCheckDecayProds){
                            ptXCheckDecayProds      = (TH1D*)inputFileMC->Get(Form("MCSecPi0FromK0s_%1.2f",yRanges[y]));
                            if (ptXCheckDecayProds){
                                haveMCxCheckInRap[y] = kTRUE;
                            }    
                        }    
                    }   
                } else {
                    partSpectrum            = (TH1D*)inputFileMC->Get("fHistoMCK0lPtRebinned2");
                    for (Int_t y=0; y< 10; y++){
                        if (!ptXCheckDecayProds){
                            ptXCheckDecayProds      = (TH1D*)inputFileMC->Get(Form("MCSecPi0FromK0l_%1.2f",yRanges[y]));
                            if (ptXCheckDecayProds){
                                haveMCxCheckInRap[y] = kTRUE;
                            }    
                        }    
                    }
                }
            }
            chHadDNdEta             = (TH1D*)inputFile->Get("histoChargedHadrondNdEtaALICEPP2760GeV");
            if (chHadDNdEta){
                cout << "found ch. had dist" << endl;
                chHadDNdEta->Sumw2();
                chHadDNdEta->Scale(1./chHadDNdEta->Integral());
            }    
            Double_t paramGraph[3]  = {partSpectrum->GetBinContent(1), 6., 0.5};
            ptDistribution          = FitObject("l","ptDistribution","K",partSpectrum,5.,20.,NULL,"QNRME");
            
            Double_t xValues[200];
            Double_t xErr[200];
            Double_t yValues[200];
            Double_t yErr[200];
            Int_t nBinsX            = 0;
            Int_t nBinsXOrg         = 0;
            Double_t maxOrSpec      = 12;
            if (kIsMC)maxOrSpec     = 8.;
            for (Int_t i = 1; (i< partSpectrum->GetNbinsX()+1 && partSpectrum->GetBinCenter(i) < maxOrSpec) ; i++){
//                 cout << partSpectrum->GetBinCenter(i) << "\t" << partSpectrum->GetBinContent(i) << endl ;
                if (partSpectrum->GetBinCenter(i) > 0.2){
//                     cout << "bin set at position: " << nBinsX << endl;
                    xValues[nBinsX] = partSpectrum->GetBinCenter(i);
                    xErr[nBinsX]    = partSpectrum->GetBinWidth(i)/2;
                    yValues[nBinsX] = partSpectrum->GetBinContent(i);
                    yErr[nBinsX]    = partSpectrum->GetBinError(i);
                    nBinsX++;
                }    
            }
            nBinsXOrg               = nBinsX;
            for (Int_t i = nBinsXOrg; (i < 200 && xValues[i-1] < 70); i++){
                xValues[i]          = maxOrSpec+(i-nBinsXOrg);    
                xErr[i]             = 0.5;    
                yValues[i]          = ptDistribution->Eval(xValues[i]);
                yErr[i]             = ptDistribution->Eval(xValues[i])*0.15;
                nBinsX++;
            }
            cout << nBinsX << endl;
            partSpectrumExt    = new TGraphErrors(nBinsX, xValues, yValues, xErr, yErr);
//             partSpectrumExt->Print();
            
            partParam               = new TSpline3(partSpectrum);
            partParamExt            = new TSpline3("extendPartSpline",xValues,yValues,nBinsX);
            for (Int_t i = 1; i<partSpectrumFromParam->GetNbinsX()+1; i++ ){
                partSpectrumFromParam->SetBinContent(i,partParamExt->Eval(partSpectrumFromParam->GetBinCenter(i)));
                partSpectrumFromParam->SetBinError(i,0);
                if (partSpectrumFromParamRest->GetBinCenter(i) > minPt && partSpectrumFromParamRest->GetBinCenter(i) < maxPt){
                    partSpectrumFromParamRest->SetBinContent(i,partParamExt->Eval(partSpectrumFromParamRest->GetBinCenter(i)));
                    partSpectrumFromParamRest->SetBinError(i,0);
                }    
//                 if (i%100 == 0){
//                     cout << partSpectrumFromParam->GetBinCenter(i) << "\t" << partSpectrumFromParam->GetBinContent(i) << endl;
//                 }    
            }
            
            TH2F* dummyDrawingHist  = new TH2F("dummyDrawingHist","dummyDrawingHist",5000,0,70,10000, 1e-14, 1); 
            SetStyleHistoTH2ForGraphs(  dummyDrawingHist, "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", 0.028, 0.04, 
                                        0.028, 0.04, 0.9, 0.95, 510, 505);
            dummyDrawingHist->Draw();
            
            DrawGammaSetMarker(partSpectrum, 20, 1.5, kAzure-6, kAzure-6);
            partSpectrum->Draw("samepe");

            if (partSpectrumAlter){
                DrawGammaSetMarker(partSpectrumAlter, 24, 1.5, kCyan+2, kCyan+2);
                partSpectrumAlter->Draw("samepe");
            }
            
            partParam->SetLineColor(kAzure+2);
            DrawGammaSetMarkerTGraphErr(partSpectrumExt, 24, 1.5, kRed+2, kRed+2);
            partSpectrumExt->Draw("same,pz");
            
            partParamExt->SetLineColor(kBlue+2);
            partParamExt->Draw("same");
            if (partSpectrumFromParam){
                partSpectrumFromParam->SetLineColor(kRed-6);
                partSpectrumFromParam->Draw("same,pe");
            }
            ptDistribution->SetRange(0,50);
            ptDistribution->SetLineColor(kRed+2);
            ptDistribution->Draw("same");
            
            ratioPartSpecToFit          = CalculateHistoRatioToFit (partSpectrum, ptDistribution); 
            ratioPartSpecToSpline       = CalculateHistoRatioToSpline (partSpectrum, partParamExt); 
            ratioPartSpecExtToSpline    = CalculateGraphErrRatioToSpline (partSpectrumExt, partParamExt);
            if (partSpectrumAlter){
                ratioPartSpecAlterToSpline  = CalculateHistoRatioToSpline (partSpectrumAlter, partParamExt); 
            }
//         } else if (energy.CompareTo("7TeV") == 0){
//             partSpectrum        = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat7TeV");
//             ptDistribution      = FitObject("tcm","ptDistribution","K",partSpectrum,0.4,20.);
//             ptDistribution->SetRange(minPt,maxPt);
        } else {
            cout << "ERROR: undefined energy for kaons" << endl;
            return;
        }    
    } else if (particle == 2){
        if (energy.CompareTo("2.76TeV") == 0){
            partSpectrum            = (TH1D*)inputFile->Get("histoLambda1115SpecStat2760GeV");
            chHadDNdEta             = (TH1D*)inputFile->Get("histoChargedHadrondNdEtaALICEPP2760GeV");
            if (chHadDNdEta){
                chHadDNdEta->Sumw2();
                chHadDNdEta->Scale(1./chHadDNdEta->Integral());
            }    

            Double_t paramGraph[3]  = {partSpectrum->GetBinContent(1), 6., 0.5};
            ptDistribution          = FitObject("l","ptDistribution","Lambda",partSpectrum,3.,15.,NULL,"QNRME");
            TF1* ptDistributionlow  = FitObject("l","ptDistribution","Lambda",partSpectrum,0.5,5,NULL,"QNRME");
            
            Double_t xValues[150];
            Double_t xErr[150];
            Double_t yValues[150];
            Double_t yErr[150];
            Int_t nBinsX            = 0;
            Int_t nBinsXOrg         = 0;
            for (Int_t i = 1; (i< partSpectrum->GetNbinsX()+1 && partSpectrum->GetBinCenter(i) < 10.) ; i++){
                cout << partSpectrum->GetBinCenter(i) << "\t" << partSpectrum->GetBinContent(i) << endl ;
                if (partSpectrum->GetBinCenter(i) > 0.2){
                    cout << "bin set at position: " << nBinsX << endl;
                    xValues[nBinsX] = partSpectrum->GetBinCenter(i);
                    xErr[nBinsX]    = partSpectrum->GetBinWidth(i)/2;
                    yValues[nBinsX] = partSpectrum->GetBinContent(i);
                    yErr[nBinsX]    = partSpectrum->GetBinError(i);
                    nBinsX++;
                }    
            }
            nBinsXOrg               = nBinsX;
            for (Int_t i = nBinsXOrg; (i < 150 && xValues[i-1] < 70); i++){
                xValues[i]          = 10+(i-nBinsXOrg);    
                xErr[i]             = 0.5;    
                yValues[i]          = ptDistribution->Eval(xValues[i]);
                yErr[i]             = ptDistribution->Eval(xValues[i])*0.15;
                nBinsX++;
            }
            cout << nBinsX << endl;
            partSpectrumExt    = new TGraphErrors(nBinsX, xValues, yValues, xErr, yErr);
//             partSpectrumExt->Print();
            
            partParam               = new TSpline3(partSpectrum);
            partParamExt            = new TSpline3("extendPartSpline",xValues,yValues,nBinsX);
            for (Int_t i = 1; i<partSpectrumFromParam->GetNbinsX()+1; i++ ){
                if (partSpectrumFromParam->GetBinCenter(i) < 0.5){
                    partSpectrumFromParam->SetBinContent(i,ptDistributionlow->Eval(partSpectrumFromParam->GetBinCenter(i)));
                    partSpectrumFromParam->SetBinError(i,0); 
                    if (partSpectrumFromParamRest->GetBinCenter(i) > minPt && partSpectrumFromParamRest->GetBinCenter(i) < maxPt){
                        partSpectrumFromParamRest->SetBinContent(i,ptDistributionlow->Eval(partSpectrumFromParamRest->GetBinCenter(i)));
                        partSpectrumFromParamRest->SetBinError(i,0);
                    }    
                } else {    
                    partSpectrumFromParam->SetBinContent(i,partParamExt->Eval(partSpectrumFromParam->GetBinCenter(i)));
                    partSpectrumFromParam->SetBinError(i,0);
                    if (partSpectrumFromParamRest->GetBinCenter(i) > minPt && partSpectrumFromParamRest->GetBinCenter(i) < maxPt){
                        partSpectrumFromParamRest->SetBinContent(i,partParamExt->Eval(partSpectrumFromParamRest->GetBinCenter(i)));
                        partSpectrumFromParamRest->SetBinError(i,0);
                    }
                }    
            }
            
            TH2F* dummyDrawingHist  = new TH2F("dummyDrawingHist","dummyDrawingHist",7000,0,70,10000, 1e-14, 1); 
            SetStyleHistoTH2ForGraphs(  dummyDrawingHist, "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", 0.028, 0.04, 
                                        0.028, 0.04, 0.9, 0.95, 510, 505);
            dummyDrawingHist->Draw();
            
            DrawGammaSetMarker(partSpectrum, 20, 1.5, kAzure-6, kAzure-6);
            partSpectrum->Draw("samepe");
            partParam->SetLineColor(kAzure+2);
            DrawGammaSetMarkerTGraphErr(partSpectrumExt, 24, 1.5, kRed+2, kRed+2);
            partSpectrumExt->Draw("same,pz");
            
            
            partParam->Draw("same");
            if (partSpectrumFromParam){
                partSpectrumFromParam->SetLineColor(kRed-6);
                partSpectrumFromParam->Draw("same,pe");
            }
            ptDistribution->SetRange(0,50);
            ptDistribution->SetLineColor(kRed+2);
            ptDistribution->Draw("same");
            
            ptDistributionlow->SetRange(0,50);
            ptDistributionlow->SetLineColor(kGreen+2);
            ptDistributionlow->Draw("same");
            
            ratioPartSpecToFit          = CalculateHistoRatioToFit (partSpectrum, ptDistribution); 
            ratioPartSpecToSpline       = CalculateHistoRatioToSpline (partSpectrum, partParamExt); 
            ratioPartSpecExtToSpline    = CalculateGraphErrRatioToSpline (partSpectrumExt, partParamExt);
            
//         } else if (energy.CompareTo("7TeV") == 0){
//             partSpectrum        = (TH1D*)inputFile->Get("histoProtonSpecPubStat7TeV");
//             ptDistribution      = FitObject("tcm","ptDistribution","P",partSpectrum,0.4,20.);
//             ptDistribution->SetRange(minPt,maxPt);
        } else {
            cout << "ERROR: undefined energy for lambda" << endl;
            return;
        }    
    }
        
    
    canvasFitQA->SaveAs(Form("%s%s_PtDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));    

    
    // plot ratio of input part to fit
    if (ratioPartSpecToFit){
        canvasFitQA->cd();
        canvasFitQA->SetLogy(0);
            TH2F* dummyDrawingRatio  = new TH2F("dummyDrawingRatio","dummyDrawingRatio",7000,0,70,1000, 0.5, 1.5); 
            SetStyleHistoTH2ForGraphs(  dummyDrawingRatio, "#it{p}_{T} (GeV/#it{c})", "#Data/Fit", 0.028, 0.04, 
                                        0.028, 0.04, 0.9, 0.95, 510, 505);
            dummyDrawingRatio->Draw();

            DrawGammaSetMarker(ratioPartSpecToFit, 20, 1.5, kAzure-6, kAzure-6);
            ratioPartSpecToFit->Draw("samepe");
            DrawGammaSetMarker(ratioPartSpecToSpline, 20, 1.5, kRed-6, kRed-6);
            ratioPartSpecToSpline->Draw("samepe");
            DrawGammaSetMarkerTGraphErr(ratioPartSpecExtToSpline, 24, 1.5, 807, 807);
            ratioPartSpecExtToSpline->Draw("same,pz");

            if (ratioPartSpecAlterToSpline){
                DrawGammaSetMarker(ratioPartSpecAlterToSpline, 24, 1.5, kCyan+2, kCyan+2);
                ratioPartSpecAlterToSpline->Draw("samepe");
            }    
            
            DrawGammaLines(0., 70,1., 1.,0.1);

        canvasFitQA->SaveAs(Form("%s%s_RatioToFitDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    }
    
    // find the correct factor due to generation in larger phi and y window
    Double_t minEtaRange            = 0;
    Double_t maxEtaRange            = 0;
    Double_t deltaEta               = 0;
    if (chHadDNdEta){
        minEtaRange                 = chHadDNdEta->GetBinCenter(1)- chHadDNdEta->GetBinWidth(1)/2;
        maxEtaRange                 = chHadDNdEta->GetBinCenter(chHadDNdEta->GetNbinsX())+ chHadDNdEta->GetBinWidth(chHadDNdEta->GetNbinsX())/2;
        deltaEta                    = maxEtaRange-minEtaRange;
    } else {
        minEtaRange                 = 0;
        maxEtaRange                 = 2*TMath::Pi();
        deltaEta                    = maxEtaRange-minEtaRange;
    }    
    Double_t scaleFactorDecayPart   = 2*TMath::Pi()* deltaEta;
    cout << "scale factor due to particle decay: " << scaleFactorDecayPart << endl;
    
    // find the scale factor due to the restriction in the pt range
    Double_t scaleFactor            = partSpectrumFromParam->Integral(partSpectrumFromParam->FindBin(minPt),partSpectrumFromParam->FindBin(maxPt))/partSpectrumFromParam->Integral();
    cout << scaleFactor << "\t" << partSpectrumFromParam->Integral(partSpectrumFromParam->FindBin(minPt),partSpectrumFromParam->FindBin(maxPt)) << "\t" << partSpectrumFromParam->Integral() << endl;    
    if (scaleFactor < 0){
        Double_t integralFullHist   = 0;
        for (Int_t iPt = 1; iPt< partSpectrumFromParam->GetNbinsX()+1; iPt++ ){
            if (partSpectrumFromParam->GetBinContent(iPt) > 0){
                integralFullHist    = integralFullHist+partSpectrumFromParam->GetBinContent(iPt);
            } else {
                cout << "bin content at " << partSpectrumFromParam->GetBinCenter(iPt) << " smaller 0." << endl;
            }    
        }    
        scaleFactor    = partSpectrumFromParam->Integral(partSpectrumFromParam->FindBin(minPt),partSpectrumFromParam->FindBin(maxPt))/integralFullHist;
        cout << "recalc scale fac, based on int by hand: "<<  scaleFactor << endl;
    }
    // put the input histo to the correct range from which it will be drawn
    ptDistribution->SetRange(minPt,maxPt);
    
    //*************************************************************************************************
    //********************* Initialize random number generators and TGenPhaseSpace ********************
    //*************************************************************************************************
    TRandom3 randy1;
    TRandom3 randy2;
    TRandom3 randy3;
    TGenPhaseSpace event;
    
    //*************************************************************************************************
    //********************** Create histograms for output *********************************************
    //*************************************************************************************************
    TH1D* h1_NEvt                       = new TH1D("h1_NEvt","", 1, 0, 1);
    h1_NEvt->SetBinContent(1,nEvts);
    TH2D *h2_ptMothervsDaughter         = new TH2D("h2_ptMothervsDaughter","", 500,0,50,700,0,70);  
    h2_ptMothervsDaughter->Sumw2();
    TH2D *h2_asym_geom                  = new TH2D("h2_asym_geom","", 200,-1,1, 700,0,70);
    h2_asym_geom->Sumw2();
    TH1D *h1_ptdistribution             = new TH1D("h1_ptdistribution","", 700,0,70);  
    h1_ptdistribution->Sumw2();
    TH1D *h1_ptdistributionInRap[10];
    for (Int_t y = 0; y < 10; y++){
        h1_ptdistributionInRap[y]    = new TH1D(Form("h1_ptdistributionInRap_%1.2f",yRanges[y]),"", 700,0,70);  
        h1_ptdistributionInRap[y]->Sumw2();
    }


    TH1D* h1_ptdistributionPerEv        = (TH1D*)h1_ptdistribution->Clone("h1_ptdistributionPerEv");
    h1_ptdistributionPerEv->Sumw2();
    TH1D *h1_phidistribution            = new TH1D("h1_phidistribution","", 100,0,2*TMath::Pi());  
    h1_phidistribution->Sumw2();
    TH1D *h1_etadistribution            = new TH1D("h1_etadistribution","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);  
    h1_etadistribution->Sumw2();
    TH1D *h1_ydistribution              = new TH1D("h1_ydistribution","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);  
    h1_ydistribution->Sumw2();
    TH1D *h1_ptDaughter[10];
    for (Int_t i = 0; i < nDaughters; i++){
        h1_ptDaughter[i]                = new TH1D(Form("h1_ptDaughter%d",i),"", 700,0,70);  
        h1_ptDaughter[i]->Sumw2();
    }
    
    TH1D *h1_ptDaughterInRap[10][10];
    for (Int_t i = 0; i < nDaughters; i++){
        for (Int_t y = 0; y < 10; y++){
            h1_ptDaughterInRap[y][i]    = new TH1D(Form("h1_ptDaughterInRap%d_%1.2f",i,yRanges[y]),"", 700,0,70);  
            h1_ptDaughterInRap[y][i]->Sumw2();
        }    
    }
    TH1D *h1_ptPiZeroDaughters          = new TH1D("h1_ptPiZeroDaughters","", 700,0,70);  
    h1_ptPiZeroDaughters->Sumw2();
    TH1D *h1_yPiZeroDaughters           = new TH1D("h1_yPiZeroDaughters","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);  
    h1_yPiZeroDaughters->Sumw2();
    TH1D *h1_ptPiZeroInRapDaughters[10];
    TH1D *h1_ptPiZeroToMotherInRap[10];
    TH1D *h1_ptPiZeroInRapDaughtersAlterNorm[10];
    for (Int_t y = 0; y < 10; y++){
        h1_ptPiZeroInRapDaughters[y]    = new TH1D(Form("h1_ptPiZeroInRapDaughters_%1.2f",yRanges[y]),"", 700,0,70);  
        h1_ptPiZeroInRapDaughters[y]->Sumw2();
        h1_ptPiZeroInRapDaughtersAlterNorm[y]    = new TH1D(Form("h1_ptPiZeroInRapDaughtersAlterNorm_%1.2f",yRanges[y]),"", 700,0,70);  
        h1_ptPiZeroInRapDaughtersAlterNorm[y]->Sumw2();
        h1_ptPiZeroToMotherInRap[y]     = new TH1D(Form("h1_ptPiZeroToMotherInRap_%1.2f",yRanges[y]),"", 700,0,70);  
        h1_ptPiZeroToMotherInRap[y]->Sumw2();
    }
    TH2D *h2_ptyPiZeroDaughters         = new TH2D("h2_ptyPiZeroDaughters","", 700,0,70,2*fYMaxMother*100,-fYMaxMother,fYMaxMother);  
    h2_ptyPiZeroDaughters->Sumw2();
    TH2D* h2_DalitzPlot                 = NULL;
    if (nDaughters == 3){
        h2_DalitzPlot                   = new TH2D("h2_DalitzPlot", "", (maxDalitz[0]-minDalitz[0])*1000, minDalitz[0],maxDalitz[0], (maxDalitz[1]-minDalitz[1])*1000, minDalitz[1], maxDalitz[1] ); 
    }    
    
    //*************************************************************************************************
    //**************************** Event loop *********************************************************
    //*************************************************************************************************
    for(Int_t n=0; n<nEvts; n++){ // this is the important loop (nEvents)
        // give a bit of stat in the printouts
        if (n%10000000 == 0) 
            cout << "generated " << (Double_t)n/1e6 << " Mio events" << endl;
        
        // draw pt to be generated for mother from constructed spectrum restricted to pt range
        Double_t ptcurrent      = partSpectrumFromParamRest->GetRandom();
        
        // asume phi is flat
        Double_t phiCurrent     = randy1.Uniform(2*TMath::Pi());
        
        // asume theta is flat
        Double_t thetaCurrent   = randy2.Uniform(2*TMath::Pi());
        Double_t etaCurrent     = -1000;
        // generate eta according to charged hadron distribution 
        if (chHadDNdEta){
            etaCurrent          = chHadDNdEta->GetRandom();
        // otherwise with the assumption theta is flat
        } else {    
            if (TMath::Cos(thetaCurrent)*TMath::Cos(thetaCurrent) < 1) 
                etaCurrent = -0.5* TMath::Log( (1.0-TMath::Cos(thetaCurrent))/(1.0+TMath::Cos(thetaCurrent)) );
            else 
                etaCurrent = -1000; 
        }    
        // weights
        Double_t weightFull     = 1./nEvts*scaleFactor*ptcurrent*scaleFactorDecayPart; // if input is dN/dydpt 1/2pi 1/pt needs to be multiplied with pt
        Double_t weightAlter    = 1./nEvts*scaleFactor*ptcurrent; // if input is dN/dydpt 1/2pi 1/pt needs to be multiplied with pt
        Double_t weightPartial  = 1./nEvts*scaleFactor;
        
        // create current particle
        TLorentzVector particle(0.0, 0.0, 0, massParticle); 
        particle.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
        // set decay properties
        event.SetDecay(particle, nDaughters, masses); 
     
        // filling of input distributions
        h1_ptdistribution->Fill(ptcurrent, weightFull);
        h1_ptdistributionPerEv->Fill(ptcurrent, weightPartial);
        h1_phidistribution->Fill(phiCurrent, weightFull);
        h1_etadistribution->Fill(etaCurrent, weightFull);
        h1_ydistribution->Fill(particle.Rapidity(), weightFull);
        for (Int_t y = 0; y < 10; y++){
            if (TMath::Abs(particle.Rapidity()) < yRanges[y]){
                h1_ptdistributionInRap[y]->Fill(ptcurrent, weightFull);
            }
        }    
        
        // let it decay
        Double_t weight         = event.Generate();
        
        
        // look at daughter quantites
        Double_t ptDaughter[nDaughters];  
        TLorentzVector* pDaughter[nDaughters];
        for (Int_t i = 0; i < nDaughters; i++){
            pDaughter[i]                = event.GetDecay(i); // these are my daughters !! 
            ptDaughter[i]               = pDaughter[i]->Pt(); 
            // look at the different daughters separately
            h1_ptDaughter[i]->Fill(ptDaughter[i],weightFull);
            // plot everything only for the pi0s in the decay chain
            if (pdgCodesDaughters[i] == 111){
                h1_ptPiZeroDaughters->Fill(ptDaughter[i], weightFull);            
                h2_ptyPiZeroDaughters->Fill(ptDaughter[i],pDaughter[i]->Rapidity(), weightFull);
                h1_yPiZeroDaughters->Fill(pDaughter[i]->Rapidity(), weightFull);
            }
            // look at relative quantites to mother
            if (i == 0)
                h2_ptMothervsDaughter->Fill(ptcurrent,ptDaughter[i],weightFull);
            // make cuts for different y ranges for later analysis
            for (Int_t y = 0; y < 10; y++){
                if (TMath::Abs(pDaughter[i]->Rapidity()) < yRanges[y]){
                    h1_ptDaughterInRap[y][i]->Fill(ptDaughter[i], weightFull);
                    if (pdgCodesDaughters[i] == 111){
                        h1_ptPiZeroInRapDaughters[y]->Fill(ptDaughter[i], weightFull);
                        h1_ptPiZeroInRapDaughtersAlterNorm[y]->Fill(ptDaughter[i], weightAlter);
                    }    
                }
            }    
        }
        // look at asymmetry for decays with more than 1 daughter
        if (nDaughters > 1){
            h2_asym_geom->Fill((ptDaughter[0]-ptDaughter[1])/(ptDaughter[0]+ptDaughter[1]), ptcurrent, weightFull);
        }
        // look at Dalitz plot for decays with 3 daug
        if (nDaughters == 3){
            TLorentzVector comb1(0.0, 0.0, 0, 0); 
            comb1.SetPxPyPzE(pDaughter[0]->Px()+pDaughter[1]->Px(),pDaughter[0]->Py()+pDaughter[1]->Py(),pDaughter[0]->Pz()+pDaughter[1]->Pz(),pDaughter[0]->E()+pDaughter[1]->E());
            TLorentzVector comb2(0.0, 0.0, 0, 0); 
            comb2.SetPxPyPzE(pDaughter[0]->Px()+pDaughter[2]->Px(),pDaughter[0]->Py()+pDaughter[2]->Py(),pDaughter[0]->Pz()+pDaughter[2]->Pz(),pDaughter[0]->E()+pDaughter[2]->E());
            h2_DalitzPlot->Fill(comb1.M(),comb2.M(), weightFull);
        }    
    }

    //*************************************************************************************************
    //******************** Multiply with correct branching Ratio for decay ****************************
    //*************************************************************************************************
    if (branchRatio != 1.){
        
        for (Int_t i = 0; i < nDaughters; i++){
            h1_ptDaughter[i]->Scale(branchRatio);
            for (Int_t y = 0; y < 10; y++){
                h1_ptDaughterInRap[y][i]->Scale(branchRatio);
            }    
        }
        h1_ptPiZeroDaughters->Scale(branchRatio);
        for (Int_t y = 0; y < 10; y++){
            h1_ptPiZeroInRapDaughters[y]->Scale(branchRatio);
        }    
        h2_ptMothervsDaughter->Scale(branchRatio);
        h2_asym_geom->Scale(branchRatio);
        h2_ptMothervsDaughter->Scale(branchRatio);
        h2_asym_geom->Scale(branchRatio);
    
    }

    //*************************************************************************************************
    //******************** Determine whether additional offset was accidentally created ***************
    //*************************************************************************************************    
    TF1* ptDistributionRefit            = NULL;
    if (particle == 0 || particle == 1 ){
        if (energy.CompareTo("2.76TeV") == 0){
            ptDistributionRefit         = FitObject("l","ptDistribution","K",h1_ptdistributionPerEv,5.,20.,NULL,"QNRME");
        } else { 
          cout << "ranges for second fit not defined" << endl;  
        }     
    } else if (particle == 2 ){
        if (energy.CompareTo("2.76TeV") == 0){
            ptDistributionRefit         = FitObject("l","ptDistribution","Lambda",h1_ptdistributionPerEv,5.,20.,NULL,"QNRME");
        } else { 
          cout << "ranges for second fit not defined" << endl;  
        }     
    }
    
    Double_t minPtForScaling            = 5;
    Double_t maxPtForScaling            = 10;
    Double_t additionalScalingFac       = 1.;
    if (ptDistributionRefit){
        additionalScalingFac            = ptDistribution->Integral(minPtForScaling,maxPtForScaling)/ptDistributionRefit->Integral(minPtForScaling,maxPtForScaling);
        if (additionalScalingFac != 1){
            cout << "somehow an offset has been generated of: " << additionalScalingFac << endl;
        }    
    }
    // scale if necessary all quantites
    if (additionalScalingFac != 1){
        cout << "rescaling the distribution according to offset" << endl;
        h1_ptdistribution->Scale(additionalScalingFac);
        h1_ptdistributionPerEv->Scale(additionalScalingFac);
        h1_phidistribution->Scale(additionalScalingFac);
        h1_etadistribution->Scale(additionalScalingFac);
        h1_ydistribution->Scale();
        for (Int_t i = 0; i < nDaughters; i++){
            h1_ptDaughter[i]->Scale(additionalScalingFac);
            for (Int_t y = 0; y < 10; y++){
                h1_ptDaughterInRap[y][i]->Scale(additionalScalingFac);
            }    
        }     
        h1_ptPiZeroDaughters->Scale(additionalScalingFac);
        for (Int_t y = 0; y < 10; y++){
            h1_ptPiZeroInRapDaughters[y]->Scale(additionalScalingFac);
        }
        h2_ptMothervsDaughter->Scale(additionalScalingFac);
        h2_ptyPiZeroDaughters->Scale(additionalScalingFac);
        h1_yPiZeroDaughters->Scale(additionalScalingFac);
        if (h2_DalitzPlot) h2_DalitzPlot->Scale(additionalScalingFac);
    
    }
    
    //*************************************************************************************************
    //******************************** Plot pt distribution   *****************************************
    //*************************************************************************************************
    
    TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1000,800);
    DrawGammaCanvasSettings( canvasQA, 0.07, 0.02, 0.02, 0.08);
    canvasQA->cd();
    canvasQA->SetLogy(1);
    TLegend* legendSpectra = GetAndSetLegend2(0.86, 0.70, 0.95, 0.93, 32,1); 
    
    DrawAutoGammaMesonHistos(   h1_ptdistribution, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_ptdistribution->GetYaxis()->SetTitleOffset(0.85);
    DrawGammaSetMarker(h1_ptdistribution, 20, 1.5, kAzure-6, kAzure-6);
    h1_ptdistribution->DrawClone("pe");
    legendSpectra->AddEntry(h1_ptdistribution,plotLabel.Data(),"pe");
    
    for (Int_t i = 0; i < nDaughters; i++){
        DrawGammaSetMarker(h1_ptDaughter[i], 21+i, 1., colors[i], colors[i]);
        h1_ptDaughter[i]->Draw("same,pe");
        legendSpectra->AddEntry(h1_ptDaughter[i],daughterLabels[i].Data(),"pe");
    }
    if (h1_ptPiZeroDaughters->GetEntries() > 0){
        DrawGammaSetMarker(h1_ptPiZeroDaughters, 24, 1.5, kGray+2, kGray+2);
        h1_ptPiZeroDaughters->Draw("same,pe");
        legendSpectra->AddEntry(h1_ptPiZeroDaughters,"#pi^{0}","pe");
    }
    if (ptXCheckDecayProds){
        DrawGammaSetMarker(ptXCheckDecayProds, 20, 1.5, kGray+2, kGray+2);
        ptXCheckDecayProds->Draw("same,pe");
//                 legendSpectra->AddEntry(ptXCheckDecayProds,"#pi^{0}, Full MC","pe");
    }    

    legendSpectra->Draw();
    canvasQA->SaveAs(Form("%s%s_PtDistribution_InputGenerated.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    
    for (Int_t y = 0;  y< 10; y++){
        canvasQA->cd();
        canvasQA->SetLogy(1);
        
        if (haveMCxCheckInRap[y]){
            DrawAutoGammaMesonHistos(   h1_ptPiZeroInRapDaughters[y], 
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{#pi^{0}}", // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
            h1_ptPiZeroInRapDaughters[y]->GetYaxis()->SetTitleOffset(0.85);
            DrawGammaSetMarker(h1_ptPiZeroInRapDaughters[y], 20, 1.5, kAzure-6, kAzure-6);
            h1_ptPiZeroInRapDaughters[y]->DrawClone("pe");

            DrawGammaSetMarker(h1_ptPiZeroInRapDaughtersAlterNorm[y], 24, 1.5, kRed-6, kRed-6);
            h1_ptPiZeroInRapDaughtersAlterNorm[y]->DrawClone("same,pe");
            
            if (ptXCheckDecayProds){
                DrawGammaSetMarker(ptXCheckDecayProds, 20, 1.5, kGray+2, kGray+2);
                ptXCheckDecayProds->Draw("same,pe");
//                 legendSpectra->AddEntry(ptXCheckDecayProds,"#pi^{0}, Full MC","pe");
            }    
            
            canvasQA->SaveAs(Form("%s%s_Pi0sInRap.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
            
        
            canvasQA->SetLogy(0);
            TH1D* ratioXCheckProd = (TH1D*)ptXCheckDecayProds->Clone("ratio");
            for (Int_t i = 1; i<ratioXCheckProd->GetNbinsX()+1; i++){
                if (ratioXCheckProd->GetBinContent(i) != 0){
                    Double_t value  = h1_ptPiZeroInRapDaughters[y]->GetBinContent(h1_ptPiZeroInRapDaughters[y]->FindBin(ratioXCheckProd->GetBinCenter(i)))/ratioXCheckProd->GetBinContent(i);
                    Double_t value1 = h1_ptPiZeroInRapDaughters[y]->GetBinContent(h1_ptPiZeroInRapDaughters[y]->FindBin(ratioXCheckProd->GetBinCenter(i)));
                    Double_t value2 = ratioXCheckProd->GetBinContent(i);
                    Double_t error1 = h1_ptPiZeroInRapDaughters[y]->GetBinError(h1_ptPiZeroInRapDaughters[y]->FindBin(ratioXCheckProd->GetBinCenter(i)));
                    Double_t error2 = ratioXCheckProd->GetBinError(i);
                    Double_t error  = TMath::Sqrt(error1*error1/(value1*value1)+error2*error2/(value2*value2))*value;
                    ratioXCheckProd->SetBinContent(i, value);
                    ratioXCheckProd->SetBinError(i, error);
                }
            }
            Double_t maxRatioXCheck = 20;
            if (particle == 1)
                maxRatioXCheck = 2000;
            TF1* constRatio     = new TF1("constRatio","[0]", 3,12);
            ratioXCheckProd->Fit(constRatio,"QRME","",1.5,12);
            cout << "offset seen off: " << constRatio->GetParameter(0) << endl;
            
            DrawAutoGammaMesonHistos(   ratioXCheckProd, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{#pi^{0}}/ reference", // (%)", 
                                kFALSE, 2, 0., kFALSE,
                                kTRUE, 0., maxRatioXCheck, 
                                kFALSE, 0., 10.);
            ratioXCheckProd->GetYaxis()->SetTitleOffset(0.85);
            DrawGammaSetMarker(ratioXCheckProd, 20, 1.5, kAzure-6, kAzure-6);
            ratioXCheckProd->DrawClone("pe");
            canvasQA->SaveAs(Form("%s%s_RatioPi0sInRap.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
            
        } else {
            cout << "no x-check histo available for |y|< " << yRanges[y] << endl;
        }    
    }
    canvasQA->cd();
    canvasQA->SetLogy(1);
    TLegend* legendSpectra2 = GetAndSetLegend2(0.7, 0.80, 0.95, 0.93, 32,1); 
    
    DrawAutoGammaMesonHistos(   h1_ptdistributionPerEv, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_ptdistributionPerEv->GetYaxis()->SetTitleOffset(0.85);
    DrawGammaSetMarker(h1_ptdistributionPerEv, 20, 1.5, kAzure-6, kAzure-6);
    h1_ptdistributionPerEv->DrawClone("pe");
    partSpectrumExt->Draw("same,p");
    
    legendSpectra2->AddEntry(h1_ptdistributionPerEv,Form("%s generated",plotLabel.Data()),"pe");
    legendSpectra2->AddEntry(partSpectrumExt,Form("%s input",plotLabel.Data()),"pe");
    legendSpectra2->Draw();
    
    canvasQA->SaveAs(Form("%s%s_PtDistribution_InputGeneratedCompInput.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    
    //*************************************************************************************************
    //******************************** Plot phi distribution   *****************************************
    //*************************************************************************************************
    
    DrawAutoGammaMesonHistos(   h1_phidistribution, 
                                "", "#varphi (rad)", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_phidistribution->GetYaxis()->SetTitleOffset(0.85);    
    h1_phidistribution->DrawClone();
    canvasQA->SaveAs(Form("%s%s_PhiDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    
    //*************************************************************************************************
    //******************************** Plot eta distribution   *****************************************
    //*************************************************************************************************

    TH1D* h1_etadistributionPlot    = (TH1D*)h1_etadistribution->Clone("Plot");
    if (chHadDNdEta){
        cout << "renormalizing" << endl;
        h1_etadistributionPlot->Scale(1/h1_etadistributionPlot->Integral());
    }
    DrawAutoGammaMesonHistos(   h1_etadistributionPlot, 
                                "", "#eta", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 1000, 1e-3, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_etadistributionPlot->GetYaxis()->SetTitleOffset(0.85);        
    h1_etadistributionPlot->DrawClone();
    if (chHadDNdEta){
        DrawGammaSetMarker(chHadDNdEta, 24, 1.5, kRed+2, kRed+2);
        chHadDNdEta->Draw("same,pe");
    }    
    canvasQA->SaveAs(Form("%s%s_EtaDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));

    
    //*************************************************************************************************
    //******************************** Plot y distribution   ******************************************
    //*************************************************************************************************

    DrawAutoGammaMesonHistos(   h1_ydistribution, 
                                "", "y", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_ydistribution->GetYaxis()->SetTitleOffset(0.85);        
    h1_ydistribution->DrawClone();
    if (h1_yPiZeroDaughters->GetEntries() > 0){
        DrawGammaSetMarker(h1_yPiZeroDaughters, 24, 1.5, kGray+2, kGray+2);
        h1_yPiZeroDaughters->Draw("same,pe");
    }

    canvasQA->SaveAs(Form("%s%s_YDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    
    
    
    //*************************************************************************************************
    //******************************** Plot asymmetry distribution   **********************************
    //*************************************************************************************************    
    TCanvas *canvasQA2D = new TCanvas("canvasQA2D","canvasQA2D",1000,800);
    DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.1, 0.02, 0.12);
    canvasQA2D->cd();    
    canvasQA2D->SetLogy(0);
    canvasQA2D->SetLogz(1);
    DrawAutoGammaHistoPaper2D(h2_ptMothervsDaughter,
                                "",
                                Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                                Form("#it{p}_{T,%s} (GeV/#it{c})",daughterLabels[0].Data()),
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);

    h2_ptMothervsDaughter->GetZaxis()->SetRangeUser( FindSmallestBin2DHist(h2_ptMothervsDaughter)*0.1, h2_ptMothervsDaughter->GetMaximum()*2);
    
    h2_ptMothervsDaughter->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_PtMothervsDaughter1.%s", fOutputDir.Data() ,outputlabel.Data(),suffix.Data()));
    DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.1, 0.02, 0.12);
    DrawAutoGammaHistoPaper2D(h2_asym_geom,
                                "",
                                "#it{A} = (#it{p}_{T,1}-#it{p}_{T,2})/(#it{p}_{T,1}+#it{p}_{T,2})",
                                Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);
    h2_asym_geom->GetZaxis()->SetRangeUser( FindSmallestBin2DHist(h2_asym_geom)*0.1, h2_asym_geom->GetMaximum()*2);
    h2_asym_geom->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_AsymmetryDaughtersVsMotherPt.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    
    if (h2_DalitzPlot){
        DrawGammaCanvasSettings( canvasQA2D, 0.12, 0.1, 0.02, 0.12);
        DrawAutoGammaHistoPaper2D(h2_DalitzPlot,
                                    "",
                                    Form("#it{M}_{%s%s} (GeV/#it{c})",daughterLabels[0].Data(),daughterLabels[1].Data()),
                                    Form("#it{M}_{%s%s} (GeV/#it{c})",daughterLabels[0].Data(),daughterLabels[2].Data()),
                                    0,0,0,
                                    0,0,5,
                                    0,0, 50,1,1.1);
        h2_DalitzPlot->GetZaxis()->SetRangeUser( FindSmallestBin2DHist(h2_DalitzPlot)*0.1, h2_DalitzPlot->GetMaximum()*2);
        h2_DalitzPlot->Draw("colz");
        canvasQA2D->SaveAs(Form("%s%s_DalitzPlot.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    }
    
    for (Int_t y = 0; y < 10; y++){
        TH1D* dummyForRatioDaughter = (TH1D*)h1_ptPiZeroInRapDaughters[y]->Clone(Form("h1_ptPiZeroForRatio_%1.2f",yRanges[y]));
        dummyForRatioDaughter->Rebin(4);
        TH1D* dummyForRatioMother   = (TH1D*)h1_ptdistributionInRap[y]->Clone(Form("h1_ptMotherForRatio_%1.2f",yRanges[y]));
        dummyForRatioMother->Rebin(4);
        
        h1_ptPiZeroToMotherInRap[y] = (TH1D*)dummyForRatioDaughter->Clone(Form("h1_ptPiZeroToMotherInRap_%1.2f",yRanges[y]));
        h1_ptPiZeroToMotherInRap[y]->Divide(h1_ptPiZeroToMotherInRap[y],dummyForRatioMother);
    
    }
    //*************************************************************************************************
    //********************** Write histograms to file *************************************************
    //*************************************************************************************************
    TFile* fileOutput = new TFile(nameOutputRootFile.Data(),"RECREATE");  
        h1_NEvt->Write();
        h1_ptdistribution->Write();
        h1_ptdistributionPerEv->Write();
        h1_phidistribution->Write();
        h1_etadistribution->Write();
        h1_ydistribution->Write();
        for (Int_t i = 0; i < nDaughters; i++){
            h1_ptDaughter[i]->Write();
            for (Int_t y = 0; y < 10; y++){
                h1_ptDaughterInRap[y][i]->Write();
            }    
        }     
        h1_ptPiZeroDaughters->Write();
        for (Int_t y = 0; y < 10; y++){
            h1_ptPiZeroInRapDaughters[y]->Write();
            h1_ptdistributionInRap[y]->Write();
            h1_ptPiZeroToMotherInRap[y]->Write();
        }
        h2_ptMothervsDaughter->Write();
        h2_ptyPiZeroDaughters->Write();
        h1_yPiZeroDaughters->Write();
        if (h2_DalitzPlot) h2_DalitzPlot->Write();
    fileOutput->Write();
    fileOutput->Close();

    cout << "adding output file name to ToyMCOutputs.txt" << endl;
    gSystem->Exec("echo '"+nameOutputRootFile+"' >> ToyMCOutputs.txt");

    
}