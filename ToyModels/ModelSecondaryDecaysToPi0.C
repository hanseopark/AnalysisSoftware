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
    
    Double_t fYMaxMother                = 10;
    TString centralityString            = GetCentralityString(fEventCutSelection);
    cout << "centrality : "<< centralityString.Data() << endl;
    
    //*************************************************************************************************
    //*************************** Initialize pt spectra ***********************************************
    //*************************************************************************************************
    Double_t massParticle       = 0;
    
    TString outputlabel         = "";
    TString plotLabel           = "";
    TString daughterLabels[10]  = { "","","","","","","","","",""};
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
    Double_t scaleFactorDecay[3]= {1, 30.971, 1};
                                    
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
    
    // Distribution in eta
    TH1D* chHadDNdEta           = NULL;
    if (inputFile){
        TString namedNdEta      = "";
        if (energy.CompareTo("900GeV") == 0){
            namedNdEta          = "histoChargedHadrondNdEtaALICEPP900GeV";
        } else if (energy.CompareTo("2.76TeV") == 0){
            namedNdEta          = "histoChargedHadrondNdEtaALICEPP2760GeV";
        } else if (energy.CompareTo("8TeV") == 0){
            namedNdEta          = "histoChargedHadrondNdEtaALICEPP8TeV";
        } else if (energy.CompareTo("7TeV") == 0){
            namedNdEta          = "histoChargedHadrondNdEtaALICEPP7TeV";
        }   
        chHadDNdEta             = (TH1D*)inputFile->Get(namedNdEta.Data());
        if (chHadDNdEta){
            cout << "found ch. had dist" << endl;
            chHadDNdEta->Sumw2();
        }
    }

    // Definition of histograms/ graphs and fits for input parametrisation
    TH1D* histoPartInputPt                  = NULL; 
    TH1D* histoPartInputPtAlter             = NULL; 
    TH1D* histoPartInputPtFromParam         = new TH1D("histoPartInputPtFromParam","histoPartInputPtFromParam",maxPt*100,0,maxPt); 
    TSpline3* splinePart                    = NULL; 
    TH1D* ratioPartSpecToFit                = NULL; 
    TH1D* ratioPartSpecToSpline             = NULL; 
    TH1D* ratioPartSpecAlterToSpline        = NULL; 
    TH1D* histoXCheckDecayProdsPt_InRap                = NULL;
    TH1D* histoXCheckMotherPt_InRap         = NULL;
    Bool_t haveMCxCheckInRap[10]            = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, 
                                                kFALSE, kFALSE, kFALSE, kFALSE, kFALSE};
    
    TGraphErrors* histoPartInputPtExt       = NULL; 
    TSpline3* splinePartExt                 = NULL; 
    TGraphErrors* ratioPartSpecExtToSpline  = NULL; 
    

    // Data inputs
    // careful what the input contains: 
    // if input is dN/dydpt 1/2pi 1/pt needs to be multiplied with pt
    // if input is dN/dydpt 1/2pi no need to be multiplied with pt
    TF1* fitPtPartInput                 = NULL; 
    Int_t nParam                        = 4;
    Double_t paramGraph[10]             = {0,0,0,0,0,0,0,0,0,0};
    Double_t fitRange[2]                = {0, 40};   
    Double_t maxOrSpec                  = 12;
    Double_t etaToBeGenerated           = 0.805;
    // K0s & K0L
    if (particle == 0 || particle == 1 ){
        etaToBeGenerated                = 0.805;
        if (!kIsMC){
            if (energy.CompareTo("2.76TeV") == 0){
                histoPartInputPt            = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat2760GeV");
                histoPartInputPtAlter       = (TH1D*)inputFile->Get("histoNeutralKaonSpecStat2760GeV");
                fitRange[0]             = 2;
                fitRange[1]             = 20;
                paramGraph[0]           = histoPartInputPt->GetBinContent(1);
                paramGraph[1]           = 6.; 
                paramGraph[2]           = 0.5;
                nParam                  = 3;
                maxOrSpec               = 12;
                fitPtPartInput          = FitObject("l","fitPtPartInput","K",NULL,fitRange[0],fitRange[1]);
            } else if (energy.CompareTo("7TeV") == 0){
              histoPartInputPt            = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat7TeV");
              histoPartInputPtAlter       = (TH1D*)inputFile->Get("histoNeutralKaonSpecStat7TeV");
              fitRange[0]             = 5;
              fitRange[1]             = 20;
              paramGraph[0]           = histoPartInputPt->GetBinContent(1);
              paramGraph[1]           = 6.;
              paramGraph[2]           = 0.5;
              nParam                  = 3;
              maxOrSpec               = 12;
              fitPtPartInput          = FitObject("l","fitPtPartInput","K",NULL,fitRange[0],fitRange[1]);
            }
        } else {
            if (particle == 0){
                histoPartInputPt            = (TH1D*)inputFileMC->Get("MC_K0s_Pt_Rebinned");
                for (Int_t y=0; y< 10; y++){
                    if (!histoXCheckDecayProdsPt_InRap){
                        histoXCheckMotherPt_InRap          = (TH1D*)inputFileMC->Get(Form("MCYield_K0s_Pt_%1.2f",yRanges[y]));
                        histoXCheckDecayProdsPt_InRap      = (TH1D*)inputFileMC->Get(Form("MCSecPi0FromK0s_%1.2f",yRanges[y]));
                        if (histoXCheckDecayProdsPt_InRap){
                            haveMCxCheckInRap[y] = kTRUE;
                        }
                    }
                }
            } else {
                histoPartInputPt            = (TH1D*)inputFileMC->Get("MC_K0l_Pt_Rebinned");
                for (Int_t y=0; y< 10; y++){
                    if (!histoXCheckDecayProdsPt_InRap){
                        histoXCheckMotherPt_InRap          = (TH1D*)inputFileMC->Get(Form("MCYield_K0l_Pt_%1.2f",yRanges[y]));
                        histoXCheckDecayProdsPt_InRap      = (TH1D*)inputFileMC->Get(Form("MCSecPi0FromK0l_%1.2f",yRanges[y]));
                        if (histoXCheckDecayProdsPt_InRap){
                            haveMCxCheckInRap[y] = kTRUE;
                        }
                    }
                }
            }
            if (energy.CompareTo("2.76TeV") == 0){
                fitRange[0]             = 5;
                fitRange[1]             = 20;
            } else if (energy.CompareTo("7TeV") == 0){
                fitRange[0]             = 5;
                fitRange[1]             = 20;
            } else if (energy.CompareTo("8TeV") == 0){
                fitRange[0]             = 5;
                fitRange[1]             = 50;
            }   
            maxOrSpec               = 8;
            nParam                  = 3;
            paramGraph[0]           = histoPartInputPt->GetBinContent(1);
            paramGraph[1]           = 6.; 
            paramGraph[2]           = 0.5;
            fitPtPartInput          = FitObject("l","fitPtPartInput","K",NULL,fitRange[0],fitRange[1]);
        }
        
        if (!histoPartInputPt){
            cout << "INFO: ToyMC not configured for kaons at this energy & data/MC, aborting here" << endl;
            return;
        }    

        for (Int_t i = 0; i<nParam; i++){
            fitPtPartInput->SetParameter(i, paramGraph[i]);
            if (i == 0)
                fitPtPartInput->SetParLimits(i, 0.01* paramGraph[i], 100*paramGraph[i]);
            else 
                fitPtPartInput->SetParLimits(i, 0.1* paramGraph[i], 10*paramGraph[i]);
        }
        histoPartInputPt->Fit(fitPtPartInput,"QNMRE","", fitRange[0],fitRange[1]);
            
        Double_t xValues[400];
        Double_t xErr[400];
        Double_t yValues[400];
        Double_t yErr[400];
        Int_t nBinsX            = 0;
        Int_t nBinsXOrg         = 0;
        Double_t widthBin       = 0.5;
        for (Int_t i = 1; (i< histoPartInputPt->GetNbinsX()+1 && histoPartInputPt->GetBinCenter(i) < maxOrSpec) ; i++){
            if (histoPartInputPt->GetBinCenter(i) > 0.2){
                xValues[nBinsX] = histoPartInputPt->GetBinCenter(i);
                xErr[nBinsX]    = histoPartInputPt->GetBinWidth(i)/2;
                yValues[nBinsX] = histoPartInputPt->GetBinContent(i);
                yErr[nBinsX]    = histoPartInputPt->GetBinError(i);
                nBinsX++;
            }
        }
        nBinsXOrg               = nBinsX;
        for (Int_t i = nBinsXOrg; (i < 400 && xValues[i-1] < 70); i++){
            xValues[i]          = xValues[i-1] + xErr[i-1]+ widthBin/2;
            xErr[i]             = widthBin/2;
            yValues[i]          = fitPtPartInput->Eval(xValues[i]);
            yErr[i]             = fitPtPartInput->Eval(xValues[i])*0.15;
            nBinsX++;
        }
        histoPartInputPtExt      = new TGraphErrors(nBinsX, xValues, yValues, xErr, yErr);
        splinePart               = new TSpline3(histoPartInputPt);
        splinePartExt            = new TSpline3("extendPartSpline",xValues,yValues,nBinsX);
        for (Int_t i = 1; i<histoPartInputPtFromParam->GetNbinsX()+1; i++ ){
            histoPartInputPtFromParam->SetBinContent(i,splinePartExt->Eval(histoPartInputPtFromParam->GetBinCenter(i)));
            histoPartInputPtFromParam->SetBinError(i,0);
        }
        ratioPartSpecToFit          = CalculateHistoRatioToFit (histoPartInputPt, fitPtPartInput); 
        ratioPartSpecToSpline       = CalculateHistoRatioToSpline (histoPartInputPt, splinePartExt); 
        ratioPartSpecExtToSpline    = CalculateGraphErrRatioToSpline (histoPartInputPtExt, splinePartExt);
        if (histoPartInputPtAlter){
            ratioPartSpecAlterToSpline  = CalculateHistoRatioToSpline (histoPartInputPtAlter, splinePartExt); 
        }
    } else if (particle == 2){
        etaToBeGenerated            = 1.6;
        TF1* fitPtPartInputlow      = NULL;
        Double_t fitRangeLow[2]     = {0.5, 5};
        if (!kIsMC){
            if (energy.CompareTo("2.76TeV") == 0){
                histoPartInputPt            = (TH1D*)inputFile->Get("histoLambda1115SpecStat2760GeV");
                fitRange[0]             = 3;
                fitRange[1]             = 15;
                paramGraph[0]           = histoPartInputPt->GetBinContent(1);
                paramGraph[1]           = 6.; 
                paramGraph[2]           = 0.5;
                nParam                  = 3;
                maxOrSpec               = 10;
                fitPtPartInput          = FitObject("l","fitPtPartInput","Lambda",NULL,fitRange[0],fitRange[1]);
                fitPtPartInputlow       = FitObject("l","fitPtPartInputLow","Lambda",NULL,fitRangeLow[0],fitRangeLow[1]);
            } else if (energy.CompareTo("7TeV") == 0){
// **************************************
// NEED TO ADD 7 TEV LAMBDA SPECTRA!!!!!*
// **************************************
              histoPartInputPt            = (TH1D*)inputFile->Get("histoLambda1115SpecStat2760GeV");
              fitRange[0]             = 3;
              fitRange[1]             = 20;
              paramGraph[0]           = histoPartInputPt->GetBinContent(1);
              paramGraph[1]           = 6.;
              paramGraph[2]           = 0.5;
              nParam                  = 3;
              maxOrSpec               = 10;
              fitPtPartInput          = FitObject("l","fitPtPartInput","Lambda",NULL,fitRange[0],fitRange[1]);
              fitPtPartInputlow       = FitObject("l","fitPtPartInputLow","Lambda",NULL,fitRangeLow[0],fitRangeLow[1]);
            }
        } else {
            histoPartInputPt            = (TH1D*)inputFileMC->Get("MC_Lambda_Pt_Rebinned");
            for (Int_t y=0; y< 10; y++){
                if (!histoXCheckDecayProdsPt_InRap){
                    histoXCheckMotherPt_InRap          = (TH1D*)inputFileMC->Get(Form("MCYield_Lambda_Pt_%1.2f",yRanges[y]));
                    histoXCheckDecayProdsPt_InRap      = (TH1D*)inputFileMC->Get(Form("MCSecPi0FromLambda_%1.2f",yRanges[y]));
                    if (histoXCheckDecayProdsPt_InRap){
                        haveMCxCheckInRap[y] = kTRUE;
                    }
                }
            }
            fitRange[0]             = 3;
            fitRange[1]             = 15;
            if (energy.CompareTo("7TeV") == 0){
               fitRange[0]             = 3;
               fitRange[1]             = 20;
            }
            paramGraph[0]           = histoPartInputPt->GetBinContent(1);
            paramGraph[1]           = 6.; 
            paramGraph[2]           = 0.5;
            nParam                  = 3;
            maxOrSpec               = 10;
            fitPtPartInput          = FitObject("l","fitPtPartInput","Lambda",NULL,fitRange[0],fitRange[1]);
            fitPtPartInputlow       = FitObject("l","fitPtPartInputLow","Lambda",NULL,fitRangeLow[0],fitRangeLow[1]);
        }    
        if (!histoPartInputPt){
            cout << "INFO: ToyMC not configured for Lambda at this energy & data/MC, aborting here"<< endl;
            return;
        }    

        for (Int_t i = 0; i<nParam; i++){
            fitPtPartInput->SetParameter(i, paramGraph[i]);
            if (i == 0)
                fitPtPartInput->SetParLimits(i, 0.01* paramGraph[i], 100*paramGraph[i]);
            else 
                fitPtPartInput->SetParLimits(i, 0.1* paramGraph[i], 10*paramGraph[i]);
        }
        histoPartInputPt->Fit(fitPtPartInput,"QNMRE","", fitRange[0],fitRange[1]);
        histoPartInputPt->Fit(fitPtPartInputlow,"QNMRE","", fitRangeLow[0],fitRangeLow[1]);    
            
        Double_t xValues[250];
        Double_t xErr[250];
        Double_t yValues[250];
        Double_t yErr[250];
        Int_t nBinsX            = 0;
        Int_t nBinsXOrg         = 0;
        for (Int_t i = 1; (i< histoPartInputPt->GetNbinsX()+1 && histoPartInputPt->GetBinCenter(i) < maxOrSpec) ; i++){
            cout << histoPartInputPt->GetBinCenter(i) << "\t" << histoPartInputPt->GetBinContent(i) << endl ;
            if (histoPartInputPt->GetBinCenter(i) > 0.2){
                cout << "bin set at position: " << nBinsX << endl;
                xValues[nBinsX] = histoPartInputPt->GetBinCenter(i);
                xErr[nBinsX]    = histoPartInputPt->GetBinWidth(i)/2;
                yValues[nBinsX] = histoPartInputPt->GetBinContent(i);
                yErr[nBinsX]    = histoPartInputPt->GetBinError(i);
                nBinsX++;
            }
        }
        nBinsXOrg               = nBinsX;
        for (Int_t i = nBinsXOrg; (i < 250 && xValues[i-1] < maxPt); i++){
            xValues[i]          = maxOrSpec+(i-nBinsXOrg);
            xErr[i]             = 0.5;
            yValues[i]          = fitPtPartInput->Eval(xValues[i]);
            yErr[i]             = fitPtPartInput->Eval(xValues[i])*0.15;
            nBinsX++;
        }
        cout << nBinsX << endl;
        histoPartInputPtExt         = new TGraphErrors(nBinsX, xValues, yValues, xErr, yErr);
        splinePart               = new TSpline3(histoPartInputPt);
        splinePartExt            = new TSpline3("extendPartSpline",xValues,yValues,nBinsX);
        for (Int_t i = 1; i<histoPartInputPtFromParam->GetNbinsX()+1; i++ ){
            if (histoPartInputPtFromParam->GetBinCenter(i) < 0.5){
                histoPartInputPtFromParam->SetBinContent(i,fitPtPartInputlow->Eval(histoPartInputPtFromParam->GetBinCenter(i)));
                histoPartInputPtFromParam->SetBinError(i,0); 
            } else {
                histoPartInputPtFromParam->SetBinContent(i,splinePartExt->Eval(histoPartInputPtFromParam->GetBinCenter(i)));
                histoPartInputPtFromParam->SetBinError(i,0);
            }
        }
        ratioPartSpecToFit          = CalculateHistoRatioToFit (histoPartInputPt, fitPtPartInput); 
        ratioPartSpecToSpline       = CalculateHistoRatioToSpline (histoPartInputPt, splinePartExt); 
        ratioPartSpecExtToSpline    = CalculateGraphErrRatioToSpline (histoPartInputPtExt, splinePartExt);
    }

    // plot input distribution and extended version
    TCanvas *canvasFitQA = new TCanvas("canvasFitQA","canvasFitQA",1000,800);
    DrawGammaCanvasSettings( canvasFitQA, 0.12, 0.015, 0.02, 0.08);
    canvasFitQA->cd();
    canvasFitQA->SetLogy(1);
    
    TH2F* dummyDrawingHist  = new TH2F("dummyDrawingHist","dummyDrawingHist",5000,0,70,10000, 1e-14, 1); 
    SetStyleHistoTH2ForGraphs(  dummyDrawingHist, "#it{p}_{T}(GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", 0.028, 0.04, 
                                0.028, 0.04, 0.9, 1.25, 510, 505);
    dummyDrawingHist->Draw();
  
        TString labelPlot[2]    = {"#Lambda","p"};
        if (particle < 2){
            labelPlot[0]        = "#frac{K^{+}+K^{-}}{2}";
            labelPlot[1]        = "K^{0}_{S}";
        }    
        TLegend* legendInput = GetAndSetLegend2(0.6, 0.50, 0.95, 0.93, 32,1);     
        legendInput->SetMargin(0.15);
        DrawGammaSetMarker(histoPartInputPt, 20, 1.5, kAzure-6, kAzure-6);
        histoPartInputPt->Draw("samepe");
        legendInput->AddEntry(histoPartInputPt,labelPlot[0].Data(),"p");
        if (histoPartInputPtAlter){
            DrawGammaSetMarker(histoPartInputPtAlter, 24, 1.5, kCyan+2, kCyan+2);
            histoPartInputPtAlter->Draw("samepe");
            legendInput->AddEntry(histoPartInputPtAlter,labelPlot[1].Data(),"p");
        }
        
        DrawGammaSetMarkerTGraphErr(histoPartInputPtExt, 24, 1.5, kRed+2, kRed+2);
        histoPartInputPtExt->Draw("same,pz");
        legendInput->AddEntry(histoPartInputPtExt,Form("%s, ext.",labelPlot[0].Data()),"p");
        
        if (histoPartInputPtFromParam){
            histoPartInputPtFromParam->SetLineColor(kRed-6);
            histoPartInputPtFromParam->Draw("same,pe");
            
        }
        splinePartExt->SetLineColor(kBlue+2);
        splinePartExt->SetLineStyle(5);
        splinePartExt->Draw("same");
        
        fitPtPartInput->SetRange(fitRange[0],fitRange[1]);
        fitPtPartInput->SetLineWidth(2);
        fitPtPartInput->SetLineStyle(7);
        fitPtPartInput->SetLineColor(kGreen+2);
        fitPtPartInput->Draw("same");

        legendInput->AddEntry(fitPtPartInput,Form("%s, fit.",labelPlot[0].Data()),"l");
        legendInput->AddEntry(splinePartExt,Form("%s, param.",labelPlot[0].Data()),"l");
        if (histoPartInputPtFromParam)
            legendInput->AddEntry(histoPartInputPtFromParam,Form("%s, from param.",labelPlot[0].Data()),"p");
        
        legendInput->Draw();
    canvasFitQA->SaveAs(Form("%s%s_PtDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    
    // plot ratio of input part to fit
    if (ratioPartSpecToFit){
        canvasFitQA->cd();
        canvasFitQA->SetLogy(0);
            TH2F* dummyDrawingRatio  = new TH2F("dummyDrawingRatio","dummyDrawingRatio",7000,0,70,1000, 0.5, 1.5); 
            SetStyleHistoTH2ForGraphs(  dummyDrawingRatio, "#it{p}_{T}(GeV/#it{c})", "#Data/Fit", 0.028, 0.04, 
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
    
    //*************************************************************************************************
    //********************* Initialize random number generators and TGenPhaseSpace ********************
    //*************************************************************************************************
    TRandom3* randy1    = new TRandom3(0);
    TRandom3* randy2    = new TRandom3(0);
    TRandom3* randy3    = new TRandom3(0);
    TGenPhaseSpace event;
    
    Int_t nBinsPt           = maxPt*10; //700;
    Double_t startBinPt     = 0.;
    Double_t endBinPt       = maxPt;
    
    //*************************************************************************************************
    //********************** Create histograms for output *********************************************
    //*************************************************************************************************
    TH1D* histoNEvents                  = new TH1D("histoNEvents","", 1, 0, 1);
    histoNEvents->SetBinContent(1,nEvts);
    TH2D *histoMothervsFirstDaughterPt  = new TH2D("histoMothervsFirstDaughterPt","", nBinsPt, startBinPt, endBinPt,nBinsPt, startBinPt, endBinPt);
    histoMothervsFirstDaughterPt->Sumw2();
    TH2D *histoAsymDaughters            = new TH2D("histoAsymDaughters","", 200,-1,1, nBinsPt, startBinPt, endBinPt);
    histoAsymDaughters->Sumw2();
    TH1D *histoMotherYieldPt            = new TH1D("histoMotherYieldPt","", nBinsPt, startBinPt, endBinPt);
    histoMotherYieldPt->Sumw2();
    TH1D *histoMotherPt_InRap[10];
    for (Int_t y = 0; y < 10; y++){
        histoMotherPt_InRap[y]          = new TH1D(Form("histoMotherPt_InRap_%1.2f",yRanges[y]),"", nBinsPt, startBinPt, endBinPt);
        histoMotherPt_InRap[y]->Sumw2();
    }


    TH1D* histoMotherInvYieldPt         = new TH1D("histoMotherInvYieldPt","", nBinsPt, startBinPt, endBinPt);
    histoMotherInvYieldPt->Sumw2();
    TH1D *histoPartPhi                  = new TH1D("histoPartPhi","", 100,0,2*TMath::Pi());
    histoPartPhi->Sumw2();
    TH1D *histMotherEta                 = new TH1D("histMotherEta","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);
    histMotherEta->Sumw2();
    TH1D *histoMotherY                  = new TH1D("histoMotherY","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);
    histoMotherY->Sumw2();
    TH1D *histoDaughterPt[10];
    for (Int_t i = 0; i < nDaughters; i++){
        histoDaughterPt[i]              = new TH1D(Form("histoDaughterPt%d",i),"", nBinsPt, startBinPt, endBinPt);
        histoDaughterPt[i]->Sumw2();
    }
    
    TH1D *histoDaughterPt_InRap[10][10];
    for (Int_t i = 0; i < nDaughters; i++){
        for (Int_t y = 0; y < 10; y++){
            histoDaughterPt_InRap[y][i] = new TH1D(Form("histoDaughterPt_InRap%d_%1.2f",i,yRanges[y]),"", nBinsPt, startBinPt, endBinPt);
            histoDaughterPt_InRap[y][i]->Sumw2();
        }
    }
    TH1D *histoPiZeroDaughtersPt        = new TH1D("histoPiZeroDaughtersPt","", nBinsPt, startBinPt, endBinPt);
    histoPiZeroDaughtersPt->Sumw2();
    TH1D *histoPi0ZeroDaughtersY        = new TH1D("histoPi0ZeroDaughtersY","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);
    histoPi0ZeroDaughtersY->Sumw2();
    TH1D *histoPiZeroDaughtersPt_InRap[10];
    TH1D *ratioPiZeroToMother_InRap[10];
    for (Int_t y = 0; y < 10; y++){
        histoPiZeroDaughtersPt_InRap[y] = new TH1D(Form("histoPiZeroDaughtersPt_InRap_%1.2f",yRanges[y]),"", nBinsPt, startBinPt, endBinPt);
        histoPiZeroDaughtersPt_InRap[y]->Sumw2();
        ratioPiZeroToMother_InRap[y]    = new TH1D(Form("ratioPiZeroToMother_InRap_%1.2f",yRanges[y]),"", nBinsPt, startBinPt, endBinPt);
        ratioPiZeroToMother_InRap[y]->Sumw2();
    }
    TH2D *histoPiZeroDaughtersPtY       = new TH2D("histoPiZeroDaughtersPtY","", nBinsPt, startBinPt, endBinPt,2*fYMaxMother*100,-fYMaxMother,fYMaxMother);
    histoPiZeroDaughtersPtY->Sumw2();
    TH2D* histoDalitzPlot               = NULL;
    if (nDaughters == 3){
        histoDalitzPlot                 = new TH2D("histoDalitzPlot", "", (maxDalitz[0]-minDalitz[0])*1000, minDalitz[0],maxDalitz[0], (maxDalitz[1]-minDalitz[1])*1000, minDalitz[1], maxDalitz[1] ); 
    }
    
    //*************************************************************************************************
    //**************************** Event loop *********************************************************
    //*************************************************************************************************
    Double_t deltaPt        = maxPt-minPt;
    
    for(Long_t n=0; n<nEvts; n++){// this is the important loop (nEvents)
        // give a bit of stat in the printouts
        if (n%1000000 == 0)
            cout << "generated " << (Double_t)n/1e6 << " Mio events" << endl;
        
        // draw pt to be generated for mother from constructed spectrum restricted to pt range
        Double_t ptcurrent      = randy3->Uniform(minPt, maxPt);
        Double_t weightPt       = splinePartExt->Eval(ptcurrent);
                        
        // asume phi is flat
        Double_t phiCurrent     = randy1->Uniform(2*TMath::Pi());
        
        // asume theta is flat
        Double_t thetaCurrent   = randy2->Uniform(2*TMath::Pi());
//         Double_t etaCurrent     = -1000;
        Double_t etaCurrent     = randy2->Uniform(-etaToBeGenerated,etaToBeGenerated);
        // generate eta according to charged hadron distribution 
//         if (chHadDNdEta){
//             etaCurrent          = chHadDNdEta->GetRandom();
//         // otherwise with the assumption theta is flat
//         } else {
//             if (TMath::Cos(thetaCurrent)*TMath::Cos(thetaCurrent) < 1) 
//                 etaCurrent = -0.5* TMath::Log( (1.0-TMath::Cos(thetaCurrent))/(1.0+TMath::Cos(thetaCurrent)) );
//             else 
//                 etaCurrent = -1000;
//         }
  
        Double_t weightPartial  = 1./nEvts*weightPt*deltaPt*1/histoMotherInvYieldPt->GetBinWidth(1);
        Double_t weightFull     = weightPartial*ptcurrent; // if input is dN/dydpt 1/2pi 1/pt needs to be multiplied with pt 
        
        // create current particle
        TLorentzVector particle(0.0, 0.0, 0, massParticle);
        particle.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
        // set decay properties
        event.SetDecay(particle, nDaughters, masses);

        // filling of input distributions
        // invariant yield
        histoMotherInvYieldPt->Fill(ptcurrent, weightPartial);
        // yield
        histoMotherYieldPt->Fill(ptcurrent, weightFull);
        
        histoPartPhi->Fill(phiCurrent, weightFull);
        histMotherEta->Fill(etaCurrent, weightFull);
        histoMotherY->Fill(particle.Rapidity(), weightFull);
        for (Int_t y = 0; y < 10; y++){
            if (TMath::Abs(particle.Rapidity()) < yRanges[y]){
                histoMotherPt_InRap[y]->Fill(ptcurrent, weightFull);
            }
        }
        
        // let it decay
        Double_t weight         = event.Generate();
        
        // look at daughter quantites
        Double_t ptDaughter[nDaughters];
        Double_t EDaughter[nDaughters];
        TLorentzVector* pDaughter[nDaughters];
        for (Int_t i = 0; i < nDaughters; i++){
            pDaughter[i]                = event.GetDecay(i);// these are my daughters !! 
            ptDaughter[i]               = pDaughter[i]->Pt();
            EDaughter[i]                = pDaughter[i]->E();
            // look at the different daughters separately
            histoDaughterPt[i]->Fill(ptDaughter[i],weightFull);
            // plot everything only for the pi0s in the decay chain
            if (pdgCodesDaughters[i] == 111){
                histoPiZeroDaughtersPt->Fill(ptDaughter[i], weightFull);
                histoPiZeroDaughtersPtY->Fill(ptDaughter[i],pDaughter[i]->Rapidity(), weightFull);
                histoPi0ZeroDaughtersY->Fill(pDaughter[i]->Rapidity(), weightFull);
            }
            // look at relative quantites to mother
            if (i == 0)
                histoMothervsFirstDaughterPt->Fill(ptcurrent,ptDaughter[i],weightFull);
            // make cuts for different y ranges for later analysis
            for (Int_t y = 0; y < 10; y++){
                if (TMath::Abs(pDaughter[i]->Rapidity()) < yRanges[y]){
                    histoDaughterPt_InRap[y][i]->Fill(ptDaughter[i], weightFull);
                    if (pdgCodesDaughters[i] == 111){
                        histoPiZeroDaughtersPt_InRap[y]->Fill(ptDaughter[i], weightFull);
                    }
                }
            }
        }
        // look at asymmetry for decays with more than 1 daughter
        if (nDaughters > 1){
            histoAsymDaughters->Fill((EDaughter[0]-EDaughter[1])/(EDaughter[0]+EDaughter[1]), ptcurrent, weightFull);
        }
        // look at Dalitz plot for decays with 3 daug
        if (nDaughters == 3){
            TLorentzVector comb1(0.0, 0.0, 0, 0); 
            comb1.SetPxPyPzE(pDaughter[0]->Px()+pDaughter[1]->Px(),pDaughter[0]->Py()+pDaughter[1]->Py(),pDaughter[0]->Pz()+pDaughter[1]->Pz(),pDaughter[0]->E()+pDaughter[1]->E());
            TLorentzVector comb2(0.0, 0.0, 0, 0); 
            comb2.SetPxPyPzE(pDaughter[0]->Px()+pDaughter[2]->Px(),pDaughter[0]->Py()+pDaughter[2]->Py(),pDaughter[0]->Pz()+pDaughter[2]->Pz(),pDaughter[0]->E()+pDaughter[2]->E());
            histoDalitzPlot->Fill(comb1.M(),comb2.M(), weightFull);
        }
    }

    //*************************************************************************************************
    //******************** Multiply with correct branching Ratio for decay ****************************
    //*************************************************************************************************
    if (branchRatio != 1.){        
        for (Int_t i = 0; i < nDaughters; i++){
            histoDaughterPt[i]->Scale(branchRatio);
            for (Int_t y = 0; y < 10; y++){
                histoDaughterPt_InRap[y][i]->Scale(branchRatio);
            }
        }
        histoPiZeroDaughtersPt->Scale(branchRatio/scaleFactorDecay[particle]);        
        for (Int_t y = 0; y < 10; y++){
            histoPiZeroDaughtersPt_InRap[y]->Scale(branchRatio/scaleFactorDecay[particle]);
        }    
    }

    //*************************************************************************************************
    //******************************** Plot pt distribution   *****************************************
    //*************************************************************************************************
    
    TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1000,800);
    DrawGammaCanvasSettings( canvasQA, 0.12, 0.02, 0.02, 0.08);
    canvasQA->cd();
    canvasQA->SetLogy(1);

    // cross check that input and generated pt distribution match
    DrawAutoGammaMesonHistos(   histoMotherInvYieldPt, 
                                "", "#it{p}_{T}(GeV/#it{c})", "#frac{1}{2#pi #it{N}_{ev}} #frac{d^{2}#it{N}}{#it{p}_{T}d#it{p}_{T}d#it{y}} (#it{c}/GeV)^{2}", // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    histoMotherInvYieldPt->GetYaxis()->SetTitleOffset(1.3);
    DrawGammaSetMarker(histoMotherInvYieldPt, 20, 1.5, kAzure-6, kAzure-6);
    histoMotherInvYieldPt->DrawClone("pe");
//     histoPartInputPtExt->Draw("same,p");
    histoPartInputPtFromParam->SetLineColor(kGreen+2);
    histoPartInputPtFromParam->SetMarkerColor(kGreen+2);
    histoPartInputPtFromParam->Draw("same,p");
    splinePartExt->SetLineColor(kRed+2);
    splinePartExt->SetLineStyle(7);
    splinePartExt->SetLineWidth(2);
    splinePartExt->Draw("same");
    
    DrawGammaSetMarker(histoPartInputPt, 24, 1.5, kGreen+2, kGreen+2);
    histoPartInputPt->Draw("same");
    
    TLegend* legendSpectra2 = GetAndSetLegend2(0.65, 0.71, 0.95, 0.93, 32,1);
    legendSpectra2->AddEntry(histoMotherInvYieldPt,Form("%s generated",plotLabel.Data()),"p");
    legendSpectra2->AddEntry(histoPartInputPtFromParam,Form("%s mod. input",plotLabel.Data()),"l");
    legendSpectra2->AddEntry(splinePartExt,Form("%s input spline",plotLabel.Data()),"l");
    legendSpectra2->AddEntry(histoPartInputPt,Form("%s orig. input hist",plotLabel.Data()),"p");
    legendSpectra2->Draw();
    
    canvasQA->SaveAs(Form("%s%s_PtDistribution_InputGeneratedCompInput.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));


    // check for MC if you have the respective MC input histograms as yields
    for (Int_t y = 0;y< 10; y++){
        canvasQA->cd();
        canvasQA->SetLogy(1);
        
        if (haveMCxCheckInRap[y]){
            // Do we have the same number of pi0s?
            DrawAutoGammaMesonHistos(   histoPiZeroDaughtersPt_InRap[y], 
                                "", "#it{p}_{T}(GeV/#it{c})",Form( "#frac{#it{N}_{#pi^0}}{#it{N}_{ev}} in |y| < %1.2f",yRanges[y]), // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
            histoPiZeroDaughtersPt_InRap[y]->GetYaxis()->SetTitleOffset(1.25);
            DrawGammaSetMarker(histoPiZeroDaughtersPt_InRap[y], 20, 1.5, kAzure-6, kAzure-6);
            histoPiZeroDaughtersPt_InRap[y]->DrawClone("pe");

            if (histoXCheckDecayProdsPt_InRap){
                DrawGammaSetMarker(histoXCheckDecayProdsPt_InRap, 20, 1.5, kGray+2, kGray+2);
                histoXCheckDecayProdsPt_InRap->Draw("same,pe");
            }
            
            canvasQA->SaveAs(Form("%s%s_Pi0sInRap.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

            // Do we have the same number of K0s?
            DrawAutoGammaMesonHistos(   histoMotherPt_InRap[y], 
                                "", "#it{p}_{T}(GeV/#it{c})", Form( "#frac{#it{N}_{K}}{#it{N}_{ev}} in |y| < %1.2f",yRanges[y] ),// (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
            histoMotherPt_InRap[y]->GetYaxis()->SetTitleOffset(1.25);
            DrawGammaSetMarker(histoMotherPt_InRap[y], 20, 1.5, kAzure-6, kAzure-6);
            histoMotherPt_InRap[y]->DrawClone("pe");

            if (histoXCheckMotherPt_InRap){
                DrawGammaSetMarker(histoXCheckMotherPt_InRap, 20, 1.5, kGray+2, kGray+2);
                histoXCheckMotherPt_InRap->Draw("same,pe");
            }
            
            canvasQA->SaveAs(Form("%s%s_InRap.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
            
            // Do we have the same number of K0s, quantify with ratio
            if (histoMotherPt_InRap[y]->GetBinWidth(1) == histoXCheckMotherPt_InRap->GetBinWidth(1)){
                canvasQA->SetLogy(0);
                TH1D* ratioXCheckMother = (TH1D*)histoXCheckMotherPt_InRap->Clone("ratio");
                for (Int_t i = 1; i<ratioXCheckMother->GetNbinsX()+1; i++){
                    if (ratioXCheckMother->GetBinContent(i) != 0){
                        Double_t value  = histoMotherPt_InRap[y]->GetBinContent(histoMotherPt_InRap[y]->FindBin(ratioXCheckMother->GetBinCenter(i)))/ratioXCheckMother->GetBinContent(i);
                        Double_t value1 = histoMotherPt_InRap[y]->GetBinContent(histoMotherPt_InRap[y]->FindBin(ratioXCheckMother->GetBinCenter(i)));
                        Double_t value2 = ratioXCheckMother->GetBinContent(i);
                        Double_t error1 = histoMotherPt_InRap[y]->GetBinError(histoMotherPt_InRap[y]->FindBin(ratioXCheckMother->GetBinCenter(i)));
                        Double_t error2 = ratioXCheckMother->GetBinError(i);
                        Double_t error  = TMath::Sqrt(error1*error1/(value1*value1)+error2*error2/(value2*value2))*value;
                        ratioXCheckMother->SetBinContent(i, value);
                        ratioXCheckMother->SetBinError(i, error);
                    }
                }
                Double_t maxRatioXCheck = 20;
                if (particle == 1)
                    maxRatioXCheck = 2000;
                TF1* constRatio     = new TF1("constRatio","[0]", 3,12);
                ratioXCheckMother->Fit(constRatio,"QRME","",1.5,12);
                cout << "Mother offset seen off: " << constRatio->GetParameter(0) << endl;
                
                DrawAutoGammaMesonHistos(   ratioXCheckMother, 
                                    "", "#it{p}_{T}(GeV/#it{c})", "#it{N}_{part}/ reference in same rap window", // (%)", 
                                    kTRUE, 2, 0.01, kFALSE,
                                    kFALSE, 0., maxRatioXCheck, 
                                    kFALSE, 0., 10.);
                ratioXCheckMother->GetYaxis()->SetTitleOffset(0.85);
                DrawGammaSetMarker(ratioXCheckMother, 20, 1.5, kAzure-6, kAzure-6);
                ratioXCheckMother->DrawClone("pe");
                canvasQA->SaveAs(Form("%s%s_RatioMotherInRap.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
            }
            if (histoPiZeroDaughtersPt_InRap[y]->GetBinWidth(1) == histoXCheckDecayProdsPt_InRap->GetBinWidth(1)){
                canvasQA->SetLogy(0);
                TH1D* ratioXCheckProd = (TH1D*)histoXCheckDecayProdsPt_InRap->Clone("ratio");
                for (Int_t i = 1; i<ratioXCheckProd->GetNbinsX()+1; i++){
                    if (ratioXCheckProd->GetBinContent(i) != 0){
                        Double_t value  = histoPiZeroDaughtersPt_InRap[y]->GetBinContent(histoPiZeroDaughtersPt_InRap[y]->FindBin(ratioXCheckProd->GetBinCenter(i)))/ratioXCheckProd->GetBinContent(i);
                        Double_t value1 = histoPiZeroDaughtersPt_InRap[y]->GetBinContent(histoPiZeroDaughtersPt_InRap[y]->FindBin(ratioXCheckProd->GetBinCenter(i)));
                        Double_t value2 = ratioXCheckProd->GetBinContent(i);
                        Double_t error1 = histoPiZeroDaughtersPt_InRap[y]->GetBinError(histoPiZeroDaughtersPt_InRap[y]->FindBin(ratioXCheckProd->GetBinCenter(i)));
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
                cout << "Pi0 offset seen off: " << constRatio->GetParameter(0) << endl;
                
                DrawAutoGammaMesonHistos(   ratioXCheckProd, 
                                    "", "#it{p}_{T}(GeV/#it{c})", "#it{N}_{#pi^{0}}/ reference in same rap window", // (%)", 
                                    kTRUE, 2, 0.01, kFALSE,
                                    kFALSE, 0., maxRatioXCheck, 
                                    kFALSE, 0., 10.);
                ratioXCheckProd->GetYaxis()->SetTitleOffset(0.85);
                DrawGammaSetMarker(ratioXCheckProd, 20, 1.5, kAzure-6, kAzure-6);
                ratioXCheckProd->DrawClone("pe");
                canvasQA->SaveAs(Form("%s%s_RatioPi0InRap.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
            }
        } else {
            cout << "no x-check histo available for |y|< " << yRanges[y] << endl;
        }
    }
    
    canvasQA->cd();
    canvasQA->SetLogy(1);

    TLegend* legendSpectra = GetAndSetLegend2(0.86, 0.70, 0.95, 0.93, 32,1);     
    DrawAutoGammaMesonHistos(   histoMotherYieldPt, 
                                "", "#it{p}_{T}(GeV/#it{c})", "#it{N}_{X}", // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    histoMotherYieldPt->GetYaxis()->SetTitleOffset(0.85);
    DrawGammaSetMarker(histoMotherYieldPt, 20, 1.5, kAzure-6, kAzure-6);
    histoMotherYieldPt->DrawClone("pe");
    legendSpectra->AddEntry(histoMotherYieldPt,plotLabel.Data(),"pe");
    
    for (Int_t i = 0; i < nDaughters; i++){
        DrawGammaSetMarker(histoDaughterPt[i], 21+i, 1., colors[i], colors[i]);
        histoDaughterPt[i]->Draw("same,pe");
        legendSpectra->AddEntry(histoDaughterPt[i],daughterLabels[i].Data(),"pe");
    }
    if (histoPiZeroDaughtersPt->GetEntries() > 0){
        DrawGammaSetMarker(histoPiZeroDaughtersPt, 24, 1.5, kGray+2, kGray+2);
        histoPiZeroDaughtersPt->Draw("same,pe");
        legendSpectra->AddEntry(histoPiZeroDaughtersPt,"#pi^{0}","pe");
    }
    legendSpectra->Draw();
    canvasQA->SaveAs(Form("%s%s_PtDistribution_InputGenerated.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    
    
    //*************************************************************************************************
    //******************************** Plot phi distribution   *****************************************
    //*************************************************************************************************
    
    DrawAutoGammaMesonHistos(   histoPartPhi, 
                                "", "#varphi (rad)", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    histoPartPhi->GetYaxis()->SetTitleOffset(0.85);
    histoPartPhi->DrawClone();
    canvasQA->SaveAs(Form("%s%s_PhiDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));

    //*************************************************************************************************
    //******************************** Plot eta distribution   *****************************************
    //*************************************************************************************************

    TH1D* histMotherEtaPlot    = (TH1D*)histMotherEta->Clone("Plot");
    if (chHadDNdEta){
        cout << "renormalizing" << endl;
        histMotherEtaPlot->Scale(1/histMotherEtaPlot->Integral());
    }
    DrawAutoGammaMesonHistos(   histMotherEtaPlot, 
                                "", "#eta", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 1000, 1e-3, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    histMotherEtaPlot->GetYaxis()->SetTitleOffset(0.85);
    histMotherEtaPlot->DrawClone();
    if (chHadDNdEta){
        DrawGammaSetMarker(chHadDNdEta, 24, 1.5, kRed+2, kRed+2);
        chHadDNdEta->Draw("same,pe");
    }
    canvasQA->SaveAs(Form("%s%s_EtaDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));

    
    //*************************************************************************************************
    //******************************** Plot y distribution   ******************************************
    //*************************************************************************************************

    DrawAutoGammaMesonHistos(   histoMotherY, 
                                "", "y", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-13, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    histoMotherY->GetYaxis()->SetTitleOffset(0.85);
    histoMotherY->DrawClone();
    if (histoPi0ZeroDaughtersY->GetEntries() > 0){
        DrawGammaSetMarker(histoPi0ZeroDaughtersY, 24, 1.5, kGray+2, kGray+2);
        histoPi0ZeroDaughtersY->Draw("same,pe");
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
    DrawAutoGammaHistoPaper2D(histoMothervsFirstDaughterPt,
                                "",
                                Form("#it{p}_{T,%s}(GeV/#it{c})",plotLabel.Data()),
                                Form("#it{p}_{T,%s}(GeV/#it{c})",daughterLabels[0].Data()),
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);

    histoMothervsFirstDaughterPt->GetZaxis()->SetRangeUser( FindSmallestBin2DHist(histoMothervsFirstDaughterPt)*0.1, histoMothervsFirstDaughterPt->GetMaximum()*2);
    
    histoMothervsFirstDaughterPt->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_PtMothervsDaughter1.%s", fOutputDir.Data() ,outputlabel.Data(),suffix.Data()));
    DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.1, 0.02, 0.12);
    DrawAutoGammaHistoPaper2D(histoAsymDaughters,
                                "",
//                                 "#it{A}= (#it{p}_{T,1}-#it{p}_{T,2})/(#it{p}_{T,1}+#it{p}_{T,2})",
                                "#it{A}= (#it{E}_{1}-#it{E}_{2})/(#it{E}_{1}+#it{E}_{2})",
                                Form("#it{p}_{T,%s}(GeV/#it{c})",plotLabel.Data()),
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);
    histoAsymDaughters->GetZaxis()->SetRangeUser( FindSmallestBin2DHist(histoAsymDaughters)*0.1, histoAsymDaughters->GetMaximum()*2);
    histoAsymDaughters->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_AsymmetryDaughtersVsMotherPt.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    
    if (histoDalitzPlot){
        DrawGammaCanvasSettings( canvasQA2D, 0.12, 0.1, 0.02, 0.12);
        DrawAutoGammaHistoPaper2D(histoDalitzPlot,
                                    "",
                                    Form("#it{M}_{%s%s}(GeV/#it{c})",daughterLabels[0].Data(),daughterLabels[1].Data()),
                                    Form("#it{M}_{%s%s}(GeV/#it{c})",daughterLabels[0].Data(),daughterLabels[2].Data()),
                                    0,0,0,
                                    0,0,5,
                                    0,0, 50,1,1.1);
        histoDalitzPlot->GetZaxis()->SetRangeUser( FindSmallestBin2DHist(histoDalitzPlot)*0.1, histoDalitzPlot->GetMaximum()*2);
        histoDalitzPlot->Draw("colz");
        canvasQA2D->SaveAs(Form("%s%s_DalitzPlot.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    }
    
    for (Int_t y = 0; y < 10; y++){
        TH1D* dummyForRatioDaughter = (TH1D*)histoPiZeroDaughtersPt_InRap[y]->Clone(Form("h1_ptPiZeroForRatio_%1.2f",yRanges[y]));
        dummyForRatioDaughter->Rebin(4);
        TH1D* dummyForRatioMother   = (TH1D*)histoMotherPt_InRap[y]->Clone(Form("h1_ptMotherForRatio_%1.2f",yRanges[y]));
        dummyForRatioMother->Rebin(4);
        
        ratioPiZeroToMother_InRap[y] = (TH1D*)dummyForRatioDaughter->Clone(Form("ratioPiZeroToMother_InRap_%1.2f",yRanges[y]));
        ratioPiZeroToMother_InRap[y]->Divide(ratioPiZeroToMother_InRap[y],dummyForRatioMother);
    
    }
    //*************************************************************************************************
    //********************** Write histograms to file *************************************************
    //*************************************************************************************************
    TFile* fileOutput = new TFile(nameOutputRootFile.Data(),"RECREATE");
        histoNEvents->Write();
        histoMotherYieldPt->Write();
        histoMotherInvYieldPt->Write();
        histoPartPhi->Write();
        histMotherEta->Write();
        histoMotherY->Write();
        for (Int_t i = 0; i < nDaughters; i++){
            histoDaughterPt[i]->Write();
            for (Int_t y = 0; y < 10; y++){
                histoDaughterPt_InRap[y][i]->Write();
            }
        }
        histoPiZeroDaughtersPt->Write();
        for (Int_t y = 0; y < 10; y++){
            histoPiZeroDaughtersPt_InRap[y]->Write();
            histoMotherPt_InRap[y]->Write();
            ratioPiZeroToMother_InRap[y]->Write();
        }
        histoMothervsFirstDaughterPt->Write();
        histoPiZeroDaughtersPtY->Write();
        histoPi0ZeroDaughtersY->Write();
        if (histoDalitzPlot) histoDalitzPlot->Write();
    fileOutput->Write();
    fileOutput->Close();

    cout << "adding output file name to ToyMCOutputs.txt" << endl;
    gSystem->Exec("echo '"+nameOutputRootFile+"' >> ToyMCOutputs.txt");

    
}
