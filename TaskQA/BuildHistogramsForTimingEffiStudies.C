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
#include "TH3F.h"
#include "TF1.h"
#include "TASImage.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TVectorT.h"
#include "TArc.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"


Bool_t debug = kFALSE;
typedef TVectorT<double> TVectorD;
typedef TVectorT<float> TVectorF;

using namespace std;

void SetLogBinningTH3(TH3* histoRebin){
    TAxis *axisafter    = histoRebin->GetZaxis();
    Int_t bins          = axisafter->GetNbins();
    Double_t from       = axisafter->GetXmin();
    Double_t to         = axisafter->GetXmax();
    Double_t *newbins   = new Double_t[bins+1];
    newbins[0]          = from;
    Double_t factor     = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
    axisafter->Set(bins, newbins);
    delete [] newbins;
}

void SetLogBinningTH2(TH2* histoRebin){
    TAxis *axisafter    = histoRebin->GetYaxis();
    Int_t bins          = axisafter->GetNbins();
    Double_t from       = axisafter->GetXmin();
    Double_t to         = axisafter->GetXmax();
    Double_t *newbins   = new Double_t[bins+1];
    newbins[0]          = from;
    Double_t factor     = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
    axisafter->Set(bins, newbins);
    delete [] newbins;
}

void SetLogBinningXTH2(TH2* histoRebin){
    TAxis *axisafter    = histoRebin->GetXaxis();
    Int_t bins          = axisafter->GetNbins();
    Double_t from       = axisafter->GetXmin();
    Double_t to         = axisafter->GetXmax();
    Double_t *newbins   = new Double_t[bins+1];
    newbins[0]          = from;
    Double_t factor     = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
    axisafter->Set(bins, newbins);
    delete [] newbins;
}

void SetLogBinningTH1(TH1* histoRebin){
    TAxis *axisafter    = histoRebin->GetXaxis();
    Int_t bins          = axisafter->GetNbins();
    Double_t from       = axisafter->GetXmin();
    Double_t to         = axisafter->GetXmax();
    Double_t *newbins   = new Double_t[bins+1];
    newbins[0]          = from;
    Double_t factor     = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
    axisafter->Set(bins, newbins);
    delete [] newbins;
}

Float_t FunctionM02(Float_t E, Float_t a, Float_t b, Float_t c, Float_t d, Float_t e){
  return ( exp( a+ b*E ) + c + d*E + e/E);
}


void BuildHistogramsForTimingEffiStudies(
    TString fileName                = "LHC18b8_fastandwoSDD_AnalysisResults_Calo.root",
    TString fCutSelection           = "000a2113_4117900010022700000",
    TString fOptionEnergy           = "XeXe_5.44TeV",
    TString nameOutputBase          = "TimingEffiTree",
    Int_t mode                      = 4,
    Long64_t maxNumMesonProcess     = -1,
    Double_t minETag                = 0.7,
    Double_t minInvMass             = 0.1,
    Double_t maxInvMass             = 0.15,
    TString optionCorrFrameworkDir  = "",
    TString suffix                  = "eps"
){


    //********************************************************************************
    //*            Definition of Cuts                                                *
    //********************************************************************************

    gROOT->Reset();

    // Set common default plot style
    StyleSettingsThesis();
    SetPlotStyle();

    //********************************************************************************
    //*            File definition/ loading                                          *
    //********************************************************************************

    TFile *f                        = (TFile*)gROOT->GetListOfFiles()->FindObject(fileName.Data());
    if (!f) {
        f                           = new TFile(fileName.Data());
    }
    if (!f) cout << "main List not found" << endl;
    if (f->IsZombie()) {
        cout <<fileName.Data() <<" file does not exist" << endl;
        f->Close();
        delete f;
        return;
    }
    TString outputDir           = Form("%s/%s", fCutSelection.Data(), suffix.Data());
    gSystem->Exec("mkdir -p "+outputDir);


    // Set collisions system
    TString collisionSystem     = ReturnFullCollisionsSystem(fOptionEnergy);
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    TString detectionProcess    = ReturnFullTextReconstructionProcess(mode);

    TString centralityString    = GetCentralityString(fCutSelection);
    if (centralityString.CompareTo("pp")==0){
        centralityString    = "";
    } else {
        if ( !centralityString.Contains("0-100%") )
            collisionSystem = Form("%s %s", centralityString.Data(), collisionSystem.Data());
    }


    //Declaration of leaves types
    TTree *TimingEffiTree = (TTree*)f->Get(Form("%s ClusterTimingEff",fCutSelection.Data()));
    if(!TimingEffiTree){ cout << "tree not found... returning!"<< endl; return;}

    Float_t timingCuts[9][2]    = { {-10000, 10000}, {-1000, 1000}, {-250, 250}, {-150, 150}, {-100, 100}, {-50, 50}, {-30, 50}, {-25, 25}, {-12.5,12.5}};
    Color_t colorTimeCut[9]     = { kBlack, kBlue+2, kViolet+2, 807, kGreen+2, kCyan+2, kBlue-6, kRed-6, 800 };
    Style_t styleMarker[9]      = { 24, 25, 42, 44, 29, 33, 20, 21, 27};

    // some variables have a seperate input variable to allow conversion of datatype
    Float_t         fInvMassMeson;                      //!<! array buffer
    Float_t         fPtMeson;                           //!<! array buffer
    Float_t         fClusterTimeTag;                    //!<! array buffer
    Float_t         fClusterTimeProbe;                  //!<! array buffer
    Float_t         fClusterETag;                       //!<! array buffer
    Float_t         fClusterEProbe;                     //!<! array buffer

    // Set branch adresses
    TimingEffiTree->SetBranchAddress("InvMass",                &fInvMassMeson);
    TimingEffiTree->SetBranchAddress("Pt",                     &fPtMeson);
    TimingEffiTree->SetBranchAddress("ClusterTimeTag",         &fClusterTimeTag);
    TimingEffiTree->SetBranchAddress("ClusterTimeProbe",       &fClusterTimeProbe);
    TimingEffiTree->SetBranchAddress("ClusterETag",            &fClusterETag);
    TimingEffiTree->SetBranchAddress("ClusterEProbe",          &fClusterEProbe);


    //********************************************************************************
    //*            Definition of Boundaries for Histograms                           *
    //********************************************************************************

    //Pt-plots
    Int_t nBinsE                    = 500;
    Int_t nBinsELog                 = 200;
    Double_t firstBinE              = 0.0;
    Double_t lastBinE               = 50.;
    Double_t firstBinELog           = 0.2;

    Int_t nBinsTime                 = 23000;
    Double_t firstBinTime           = -1500;
    Double_t lastBinTime            = 10000;
    if (mode == 4){
        nBinsTime                   = 14000;
        firstBinTime                = -1000;
        lastBinTime                 = 6000;
        firstBinELog                = 0.6;
    }


    //********************************************************************************
    //*      Definition of histograms for reconstructed Conversion Points            *
    //********************************************************************************
    TH2F* histoMesonInvMassvsPt[9]          = {NULL};
    TH2F* histoClusterTimevsETag[9]         = {NULL};
    TH2F* histoClusterTimevsEProbe[9]       = {NULL};
    TH1F* histoClusterETag[9]               = {NULL};
    TH1F* histoClusterEProbe[9]             = {NULL};
    TH1F* histoClusterEProbeWCut[9]         = {NULL};
    TH1F* histoClusterTimingEff[9]          = {NULL};

    for (Int_t k = 0; k<9; k++){
        cout << k << "\t" << timingCuts[k][0] << "\t"<< timingCuts[k][1] << endl;
        histoMesonInvMassvsPt[k]        = new TH2F(Form("histoMesonInvMassvsPt_%d", k),"",400, 0, 0.2, nBinsE, firstBinE, lastBinE);
        histoClusterTimevsETag[k]       = new TH2F(Form("histoClusterTimevsETag_%d",k),"",nBinsTime, firstBinTime, lastBinTime, nBinsE, firstBinE, lastBinE);
        histoClusterTimevsEProbe[k]     = new TH2F(Form("histoClusterTimevsEProbe_%d",k),"",nBinsTime, firstBinTime, lastBinTime, nBinsE, firstBinE, lastBinE);
        histoClusterETag[k]             = new TH1F(Form("histoClusterETag_%d",k),"",nBinsELog, firstBinELog, lastBinE);
        histoClusterEProbe[k]           = new TH1F(Form("histoClusterEProbe%d",k),"",nBinsELog, firstBinELog, lastBinE);
        histoClusterEProbeWCut[k]       = new TH1F(Form("histoClusterEProbeWCut%d",k),"",nBinsELog, firstBinELog, lastBinE);
        SetLogBinningTH1(histoClusterETag[k]);
        SetLogBinningTH1(histoClusterEProbe[k]);
        SetLogBinningTH1(histoClusterEProbeWCut[k]);
    }

    //********************************************************************************
    //*      Reading of Tree with reconstructed gammas/ filling of histograms        *
    //********************************************************************************

    Long64_t nEntriesRecMeson                 = TimingEffiTree->GetEntries();
    if(maxNumMesonProcess>0 && maxNumMesonProcess<nEntriesRecMeson)
        nEntriesRecMeson = maxNumMesonProcess;
    cout << "Number of mesons to be processed: " << nEntriesRecMeson << endl;

    for (Long64_t i=0; i<nEntriesRecMeson;i++) {
        TimingEffiTree->GetEntry(i);
        if (debug) cout << i << "\t" << fClusterETag << "\t" << fClusterTimeTag << "\t" << fClusterTimeProbe<< endl;
        if ( fClusterETag < minETag){ if (debug) cout << "failed E" << endl; continue;}
        for (Int_t k = 0; k < 9; k++){
            if (! (fClusterTimeTag > timingCuts[k][1]*1e-9 || fClusterTimeTag < timingCuts[k][0]*1e-9 )  && (fInvMassMeson < maxInvMass && fInvMassMeson > minInvMass)) {
                histoMesonInvMassvsPt[k]->Fill(fInvMassMeson,fPtMeson);
                histoClusterTimevsETag[k]->Fill(fClusterTimeTag*1e9,fClusterETag);
                histoClusterTimevsEProbe[k]->Fill(fClusterTimeProbe*1e9,fClusterEProbe);
                histoClusterETag[k]->Fill(fClusterETag);
                histoClusterEProbe[k]->Fill(fClusterEProbe);
                if (! (fClusterTimeProbe > timingCuts[k][1]*1e-9 || fClusterTimeProbe < timingCuts[k][0]*1e-9 ) ) {
                    histoClusterEProbeWCut[k]->Fill(fClusterEProbe);
                }
            }
        }
    }

    for (Int_t k = 0; k<9; k++){
        histoClusterETag[k]->Sumw2();
        histoClusterEProbe[k]->Sumw2();
        histoClusterEProbeWCut[k]->Sumw2();
        histoClusterTimingEff[k]    = (TH1F*)histoClusterEProbeWCut[k]->Clone(Form("histoClusterTimingEff_%d", k));
        histoClusterTimingEff[k]->Divide(histoClusterTimingEff[k],histoClusterEProbe[k],1.,1.,"B");
    }

    TCanvas* canvasCompEffSimple = new TCanvas("canvasCompEffSimple","",200,10,1350,900);// gives the page size
    DrawGammaCanvasSettings( canvasCompEffSimple, 0.07, 0.01, 0.01, 0.09);
    canvasCompEffSimple->SetLogx();
        TH1F* histo1DDummy       = new TH1F("histo1DDummy","histo1DDummy",1000, firstBinELog, 50);
        SetStyleHistoTH1ForGraphs(histo1DDummy, "#it{E} (GeV)", "#it{#varepsilon}_{t}", 0.035 ,0.04, 0.035,0.04, 0.9, 0.85, 510, 505);
        if (mode == 5)
            histo1DDummy->GetYaxis()->SetRangeUser(0,1.2);
        else if (mode == 4)
            histo1DDummy->GetYaxis()->SetRangeUser(0.95,1.015);
        histo1DDummy->GetXaxis()->SetLabelOffset(-0.005);
        histo1DDummy->GetXaxis()->SetMoreLogLabels();
        histo1DDummy->DrawCopy();
        Int_t nRows = 5;
        if (mode ==4)
            nRows   = 4;

        TLegend* legendEffi = GetAndSetLegend2(0.35,0.11,0.95,0.11+1.15*0.032*nRows, 1000*0.032, 2, "", 43, 0.1);

        for (Int_t k = 8; k > -1; k-- ){
            if (mode == 4 && k < 1 ) continue;
            DrawGammaSetMarker(histoClusterTimingEff[k], styleMarker[k], 1.5, colorTimeCut[k], colorTimeCut[k]);
            histoClusterTimingEff[k]->Draw("e1p,same");
            legendEffi->AddEntry(histoClusterTimingEff[k], Form("%4.1f < #it{t} (ns) < %4.1f", timingCuts[k][0],timingCuts[k][1]), "p" );
        }
        legendEffi->Draw();
        histo1DDummy->Draw("same,axis");

        TLatex *labelEnergy = new TLatex(0.11,0.93,Form("%s",collisionSystem.Data()));
        SetStyleTLatex( labelEnergy, 0.04,4);
        labelEnergy->Draw();
        TLatex *labelDet = new TLatex(0.11,0.88,Form("%s",detectionProcess.Data()));
        SetStyleTLatex( labelDet, 0.04,4);
        labelDet->Draw();

    canvasCompEffSimple->Update();
    canvasCompEffSimple->SaveAs(Form("%s/TimingEffi.%s",outputDir.Data(),suffix.Data()));


    //********************************************************************************
    //*                  Writing histograms to outputfile                            *
    //********************************************************************************
    TString fileNameOutput = Form("%s_Data.root",nameOutputBase.Data()) ;
    TFile* fileTimeEffiWrite = new TFile(fileNameOutput.Data(),"RECREATE");
        for (Int_t k = 0; k < 9; k++){
            fileTimeEffiWrite->mkdir(Form("%4.1fns_%4.1fns",timingCuts[k][0], timingCuts[k][1] ));
            fileTimeEffiWrite->cd(Form("%4.1fns_%4.1fns",timingCuts[k][0], timingCuts[k][1] ));
            histoMesonInvMassvsPt[k]->Write("histoMesonInvMassvsPt",TObject::kWriteDelete);
            histoClusterTimevsETag[k]->Write("histoClusterTimevsETag",TObject::kWriteDelete);
            histoClusterTimevsEProbe[k]->Write("histoClusterTimevsEProbe",TObject::kWriteDelete);
            histoClusterETag[k]->Write("histoClusterETag",TObject::kWriteDelete);
            histoClusterEProbe[k]->Write("histoClusterEProbe",TObject::kWriteDelete);
            histoClusterEProbeWCut[k]->Write("histoClusterEProbeWCut",TObject::kWriteDelete);
            histoClusterTimingEff[k]->Write("histoClusterTimingEff",TObject::kWriteDelete);
        }
    fileTimeEffiWrite->Write();
    fileTimeEffiWrite->Close();
    delete fileTimeEffiWrite;

    f->Close();
    delete f;

}
