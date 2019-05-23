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
#include "../CommonHeaders/ConversionFunctions.h"

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

Float_t FunctionM02(Float_t E, Float_t a, Float_t b, Float_t c, Float_t d, Float_t e){
  return ( exp( a+ b*E ) + c + d*E + e/E);
}


void BuildHistogramsForClusterQA(
    TString fileName                = "LHC18b8_fastandwoSDD_AnalysisResults_Calo.root",
    TString fCutSelection           = "000a2113_4117900010022700000",
    Bool_t kMC                      = 0,
    TString nameOutputBase          = "ClusterQA",
    Long64_t maxNumClusProcess      = -1,
    Int_t currPtHBin                = -1
){


    //********************************************************************************
    //*            Definition of Cuts                                                *
    //********************************************************************************

    gROOT->Reset();

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

    //Declaration of leaves types

    TTree *ClusterQA = (TTree*)f->Get(Form("ClusterQA_%s",fCutSelection.Data()));
    if(!ClusterQA){ cout << "tree not found... returning!"<< endl; return;}

    const Int_t kMaxActiveCells = 18000;
    const Int_t kMaxNTracks = 4010;

    // some variables have a seperate input variable to allow conversion of datatype
    Float_t         fBuffer_ClusterE;                    //!<! array buffer
    Float_t         fBuffer_ClusterPhi;                  //!<! array buffer
    Float_t         fBuffer_ClusterEta;                  //!<! array buffer
    Bool_t          fBuffer_ClusterIsEMCAL;              //!<! array buffer
    Int_t           fBuffer_ClusterSupMod;              //!<! array buffer
    Int_t           fBuffer_ClusterNumCells;             //!<! array buffer
    Int_t           fBuffer_LeadingCell_ID;        //!<! array buffer
    Float_t         fBuffer_LeadingCell_E;              //!<! array buffer
    Float_t         fBuffer_LeadingCell_Eta;              //!<! array buffer
    Float_t         fBuffer_LeadingCell_Phi;              //!<! array buffer
    Float_t         fBuffer_ClusterM02;              //!<! array buffer
    Float_t         fBuffer_ClusterM20;              //!<! array buffer

    Float_t         fBuffer_JJWeight;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_X;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Y;               //!<! array buffer
    Float_t         fBuffer_Event_Vertex_Z;               //!<! array buffer
    Float_t         fBuffer_Event_Multiplicity;             //!<! array buffer
    Int_t           fBuffer_Event_NumActiveCells;          //!<! array buffer

    Int_t*          fBuffer_Input_Cells_ID          = new Int_t[kMaxActiveCells];                      //!<! array buffer
    Float_t*        fBuffer_Input_Cells_E           = new Float_t[kMaxActiveCells];                      //!<! array buffer
    Float_t*        fBuffer_Input_Cells_RelativeEta = new Float_t[kMaxActiveCells];                      //!<! array buffer
    Float_t*        fBuffer_Input_Cells_RelativePhi = new Float_t[kMaxActiveCells];                      //!<! array buffer

    Int_t           fBuffer_Surrounding_NTracks;                //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_R           = new Float_t[kMaxNTracks];               //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_Pt          = new Float_t[kMaxNTracks];              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_P          = new Float_t[kMaxNTracks];              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_RelativeEta = new Float_t[kMaxNTracks];              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Tracks_RelativePhi = new Float_t[kMaxNTracks];              //!<! array buffer

    Int_t           fBuffer_Cluster_MC_Label;              //!<! array buffer
    Int_t           fBuffer_Mother_MC_Label;              //!<! array buffer
    Float_t        fBuffer_Cluster_MC_EFracFirstLabel;              //!<! array buffer
    Float_t        fBuffer_Cluster_MC_EFracLeadingPi0;              //!<! array buffer

    // Set branch adresses
    ClusterQA->SetBranchAddress("Cluster_E",                         &fBuffer_ClusterE);
    ClusterQA->SetBranchAddress("Cluster_Eta",                       &fBuffer_ClusterEta);
    ClusterQA->SetBranchAddress("Cluster_Phi",                       &fBuffer_ClusterPhi);
    ClusterQA->SetBranchAddress("Cluster_IsEMCAL",                   &fBuffer_ClusterIsEMCAL);
    ClusterQA->SetBranchAddress("Cluster_SM",                        &fBuffer_ClusterSupMod);
    ClusterQA->SetBranchAddress("Cluster_NumCells",                  &fBuffer_ClusterNumCells);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_ID",            &fBuffer_LeadingCell_ID);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_E",             &fBuffer_LeadingCell_E);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_Eta",           &fBuffer_LeadingCell_Eta);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_Phi",           &fBuffer_LeadingCell_Phi);
    ClusterQA->SetBranchAddress("Cluster_M02",                       &fBuffer_ClusterM02);
    ClusterQA->SetBranchAddress("Cluster_M20",                       &fBuffer_ClusterM20);
    // ClusterQA->SetBranchAddress("Event_Weight",                     &fBuffer_JJWeight);
    ClusterQA->SetBranchAddress("Event_Vertex_X",                  &fBuffer_Event_Vertex_X);
    ClusterQA->SetBranchAddress("Event_Vertex_Y",                  &fBuffer_Event_Vertex_Y);
    ClusterQA->SetBranchAddress("Event_Vertex_Z",                  &fBuffer_Event_Vertex_Z);
    ClusterQA->SetBranchAddress("Event_Multiplicity",              &fBuffer_Event_Multiplicity);
    ClusterQA->SetBranchAddress("Event_NumActiveCells",            &fBuffer_Event_NumActiveCells);
    ClusterQA->SetBranchAddress("Cluster_Cells_ID",                fBuffer_Input_Cells_ID);
    ClusterQA->SetBranchAddress("Cluster_Cells_E",                 fBuffer_Input_Cells_E);
    ClusterQA->SetBranchAddress("Cluster_Cells_RelativeEta",       fBuffer_Input_Cells_RelativeEta);
    ClusterQA->SetBranchAddress("Cluster_Cells_RelativePhi",       fBuffer_Input_Cells_RelativePhi);
    ClusterQA->SetBranchAddress("Surrounding_NTracks",             &fBuffer_Surrounding_NTracks);
    ClusterQA->SetBranchAddress("Surrounding_Tracks_R",            fBuffer_Surrounding_Tracks_R);
    ClusterQA->SetBranchAddress("Surrounding_Tracks_Pt",           fBuffer_Surrounding_Tracks_Pt);
    ClusterQA->SetBranchAddress("Surrounding_Tracks_P",            fBuffer_Surrounding_Tracks_P);
    ClusterQA->SetBranchAddress("Surrounding_Tracks_RelativeEta",  fBuffer_Surrounding_Tracks_RelativeEta);
    ClusterQA->SetBranchAddress("Surrounding_Tracks_RelativePhi",  fBuffer_Surrounding_Tracks_RelativePhi);

    if(kMC)
    {
      ClusterQA->SetBranchAddress("Cluster_MC_Label",                &fBuffer_Cluster_MC_Label);
      ClusterQA->SetBranchAddress("Mother_MC_Label",                 &fBuffer_Mother_MC_Label);
      ClusterQA->SetBranchAddress("Cluster_MC_EFracFirstLabel",                 &fBuffer_Cluster_MC_EFracFirstLabel);
      ClusterQA->SetBranchAddress("Cluster_MC_EFracLeadingPi0",                 &fBuffer_Cluster_MC_EFracLeadingPi0);
    }

    Double_t weightsBinsLHC18b8[20]      = { 16.1083,      4.60917,     2.15196,     0.782021,    0.26541,
                                           0.0978374,   0.0294286,   0.00989457,  0.0040615,   0.00135787,
                                           0.000531766, 0.000188772, 9.23331e-05, 4.30245e-05, 2.10196e-05,
                                           1.06695e-05, 5.78742e-06, 3.02897e-06, 1.62702e-06, 2.12118e-06 };

    if(currPtHBin>0 && kMC){
        fBuffer_JJWeight = weightsBinsLHC18b8[currPtHBin-1];
        cout << "analysing pT hard bin " << currPtHBin << " with weight " << fBuffer_JJWeight << endl;
    }else
        fBuffer_JJWeight = 1;

    //********************************************************************************
    //*            Definition of Boundaries for Histograms                           *
    //********************************************************************************

    //Pt-plots
    Int_t nBinsE                    = 200;
    Double_t firstBinE              = 0.01;
    Double_t lastBinE               = 200.;

    Int_t nBinsM02                   = 220;
    Double_t firstBinM02             = 0.0;
    Double_t lastBinM02              = 2.2;


    //********************************************************************************
    //*      Definition of histograms for reconstructed Conversion Points            *
    //********************************************************************************
    TH1F* histoM02perSM[20]             = {NULL};
    TH2F* histoM02perSMvsE[20]          = {NULL};
    TH1F* histoM20perSM[20]             = {NULL};
    TH2F* histoM20perSMvsE[20]          = {NULL};
    TH2F* histoM02vsM20perSM[20]        = {NULL};
    TH3F* histoM02vsM20perSMvsE[20]     = {NULL};

    TH2F* histoEFracFirstLabelvsE       = NULL;
    TH2F* histoEFracLeadPi0vsE          = NULL;

    TH1F* histoM02TrueMergedPi0perSM[20]             = {NULL};
    // TH1F* histoM02TrueMergedEtaperSM[20]             = {NULL};
    TH1F* histoM02TrueMergedPartConvPi0perSM[20]             = {NULL};
    // TH1F* histoM02TrueMergedPartConvEtaperSM[20]             = {NULL};
    TH1F* histoM02TrueGammaFromPi0perSM[20]             = {NULL};
    // TH1F* histoM02TrueGammaFromEtaperSM[20]             = {NULL};
    TH1F* histoM02TrueElectronFromPi0perSM[20]             = {NULL};
    // TH1F* histoM02TrueElectronFromEtaperSM[20]             = {NULL};
    TH1F* histoM02TrueBackgroundperSM[20]             = {NULL};

    TH2F* histoM02TrueMergedPi0perSMvsE[20]             = {NULL};
    // TH2F* histoM02TrueMergedEtaperSMvsE[20]             = {NULL};
    TH2F* histoM02TrueMergedPartConvPi0perSMvsE[20]             = {NULL};
    // TH2F* histoM02TrueMergedPartConvEtaperSMvsE[20]             = {NULL};
    TH2F* histoM02TrueGammaFromPi0perSMvsE[20]             = {NULL};
    // TH2F* histoM02TrueGammaFromEtaperSMvsE[20]             = {NULL};
    TH2F* histoM02TrueElectronFromPi0perSMvsE[20]             = {NULL};
    // TH2F* histoM02TrueElectronFromEtaperSMvsE[20]             = {NULL};
    TH2F* histoM02TrueBackgroundperSMvsE[20]             = {NULL};


    TH3F* histoTrueMergedPi0M02vsM20perSMvsE[20]     = {NULL};
    TH3F* histoTrueMergedPartConvPi0M02vsM20perSMvsE[20]     = {NULL};
    TH3F* histoTrueGammaFromPi0M02vsM20perSMvsE[20]     = {NULL};
    TH3F* histoTrueElectronFromPi0M02vsM20perSMvsE[20]     = {NULL};

    for(Int_t i=0;i<20;i++){
        histoM02perSM[i]                = new TH1F(Form("histoM02perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
        histoM02perSMvsE[i]             = new TH2F(Form("histoM02perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        SetLogBinningTH2(histoM02perSMvsE[i]);
        histoM20perSM[i]                = new TH1F(Form("histoM20perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
        histoM20perSMvsE[i]             = new TH2F(Form("histoM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        SetLogBinningTH2(histoM20perSMvsE[i]);
        histoM02vsM20perSM[i]           = new TH2F(Form("histoM02vsM20perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02);
        histoM02vsM20perSMvsE[i]        = new TH3F(Form("histoM02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        SetLogBinningTH3(histoM02vsM20perSMvsE[i]);
    }
    if (kMC){
        histoEFracFirstLabelvsE             = new TH2F("histoEFracFirstLabelvsE","",20, 0., 1., nBinsE, firstBinE, lastBinE);
        histoEFracLeadPi0vsE                = new TH2F("histoEFracLeadPi0vsE","",100, 0., 1., nBinsE, firstBinE, lastBinE);
        for(Int_t i=0;i<20;i++){
            histoM02TrueMergedPi0perSM[i]                   = new TH1F(Form("histoM02TrueMergedPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueMergedEtaperSM[i]                   = new TH1F(Form("histoM02TrueMergedEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            histoM02TrueMergedPartConvPi0perSM[i]           = new TH1F(Form("histoM02TrueMergedPartConvPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueMergedPartConvEtaperSM[i]           = new TH1F(Form("histoM02TrueMergedPartConvEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            histoM02TrueGammaFromPi0perSM[i]                = new TH1F(Form("histoM02TrueGammaFromPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueGammaFromEtaperSM[i]                = new TH1F(Form("histoM02TrueGammaFromEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            histoM02TrueElectronFromPi0perSM[i]             = new TH1F(Form("histoM02TrueElectronFromPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueElectronFromEtaperSM[i]             = new TH1F(Form("histoM02TrueElectronFromEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            histoM02TrueBackgroundperSM[i]                  = new TH1F(Form("histoM02TrueBackgroundperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);

            histoM02TrueMergedPi0perSMvsE[i]                = new TH2F(Form("histoM02TrueMergedPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueMergedEtaperSMvsE[i]                = new TH2F(Form("histoM02TrueMergedEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            histoM02TrueMergedPartConvPi0perSMvsE[i]        = new TH2F(Form("histoM02TrueMergedPartConvPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueMergedPartConvEtaperSMvsE[i]        = new TH2F(Form("histoM02TrueMergedPartConvEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            histoM02TrueGammaFromPi0perSMvsE[i]             = new TH2F(Form("histoM02TrueGammaFromPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueGammaFromEtaperSMvsE[i]             = new TH2F(Form("histoM02TrueGammaFromEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            histoM02TrueElectronFromPi0perSMvsE[i]          = new TH2F(Form("histoM02TrueElectronFromPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueElectronFromEtaperSMvsE[i]          = new TH2F(Form("histoM02TrueElectronFromEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            histoM02TrueBackgroundperSMvsE[i]               = new TH2F(Form("histoM02TrueBackgroundperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);

            histoTrueMergedPi0M02vsM20perSMvsE[i]           = new TH3F(Form("histoTrueMergedPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            histoTrueMergedPartConvPi0M02vsM20perSMvsE[i]   = new TH3F(Form("histoTrueMergedPartConvPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            histoTrueGammaFromPi0M02vsM20perSMvsE[i]        = new TH3F(Form("histoTrueGammaFromPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            histoTrueElectronFromPi0M02vsM20perSMvsE[i]     = new TH3F(Form("histoTrueElectronFromPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        }
    }

    for(Int_t i=0;i<20;i++){
        histoM02perSM[i]->Sumw2();
        histoM02perSMvsE[i]->Sumw2();
        histoM20perSM[i]->Sumw2();
        histoM20perSMvsE[i]->Sumw2();
        histoM02vsM20perSM[i]->Sumw2();
        histoM02vsM20perSMvsE[i]->Sumw2();
    }
    if (kMC){
        histoEFracFirstLabelvsE->Sumw2();
        histoEFracLeadPi0vsE->Sumw2();
        for(Int_t i=0;i<20;i++){
            histoM02TrueMergedPi0perSM[i]->Sumw2();
            // histoM02TrueMergedEtaperSM[i]->Sumw2();
            histoM02TrueMergedPartConvPi0perSM[i]->Sumw2();
            // histoM02TrueMergedPartConvEtaperSM[i]->Sumw2();
            histoM02TrueGammaFromPi0perSM[i]->Sumw2();
            // histoM02TrueGammaFromEtaperSM[i]->Sumw2();
            histoM02TrueElectronFromPi0perSM[i]->Sumw2();
            // histoM02TrueElectronFromEtaperSM[i]->Sumw2();
            histoM02TrueBackgroundperSM[i]->Sumw2();

            histoM02TrueMergedPi0perSMvsE[i]->Sumw2();
            // histoM02TrueMergedEtaperSMvsE[i]->Sumw2();
            histoM02TrueMergedPartConvPi0perSMvsE[i]->Sumw2();
            // histoM02TrueMergedPartConvEtaperSMvsE[i]->Sumw2();
            histoM02TrueGammaFromPi0perSMvsE[i]->Sumw2();
            // histoM02TrueGammaFromEtaperSMvsE[i]->Sumw2();
            histoM02TrueElectronFromPi0perSMvsE[i]->Sumw2();
            // histoM02TrueElectronFromEtaperSMvsE[i]->Sumw2();
            histoM02TrueBackgroundperSMvsE[i]->Sumw2();

            histoTrueMergedPi0M02vsM20perSMvsE[i]->Sumw2();
            histoTrueMergedPartConvPi0M02vsM20perSMvsE[i]->Sumw2();
            histoTrueGammaFromPi0M02vsM20perSMvsE[i]->Sumw2();
            histoTrueElectronFromPi0M02vsM20perSMvsE[i]->Sumw2();
        }
    }

    //********************************************************************************
    //*      Reading of Tree with reconstructed gammas/ filling of histograms        *
    //********************************************************************************

    Long64_t nEntriesRecClus                 = ClusterQA->GetEntries();
    if(maxNumClusProcess>0 && maxNumClusProcess<nEntriesRecClus)
        nEntriesRecClus = maxNumClusProcess;
    cout << "Number of Clusters to be processed: " << nEntriesRecClus << endl;

    for (Long64_t i=0; i<nEntriesRecClus;i++) {
        ClusterQA->GetEntry(i);
        // if(fBuffer_ClusterE<10) continue;
        // if(fBuffer_ClusterE>20) continue;
        if(fBuffer_ClusterM02<0.1) continue;
        Float_t isoTrackPt = 0;
        for(Int_t itr=0; itr<fBuffer_Surrounding_NTracks; itr++){
            if(fBuffer_Surrounding_Tracks_R[itr] < TMath::Power(0.4,2)){
                isoTrackPt+=fBuffer_Surrounding_Tracks_Pt[itr];
            }
        }
        if(isoTrackPt>2.0){
            // cout << "not iso" << endl;
             continue;
        }
        // cout << "\t ISOLATED" << endl;
        histoM02perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
        histoM02perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
        histoM20perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM20,fBuffer_JJWeight);
        histoM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
        histoM02vsM20perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_JJWeight);
        histoM02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);

        if (kMC){
            histoEFracFirstLabelvsE->Fill(fBuffer_Cluster_MC_EFracFirstLabel,fBuffer_ClusterE,fBuffer_JJWeight);
            histoEFracLeadPi0vsE->Fill(fBuffer_Cluster_MC_EFracLeadingPi0,fBuffer_ClusterE,fBuffer_JJWeight);
            switch(fBuffer_Cluster_MC_Label){
                case 10:  //
                case 11:  //
                    histoM02TrueMergedPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    histoM02TrueMergedPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    histoTrueMergedPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 12:  //
                    // histoM02TrueMergedEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    // histoM02TrueMergedEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 13:  //
                case 14:  //
                    histoM02TrueMergedPartConvPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    histoM02TrueMergedPartConvPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    histoTrueMergedPartConvPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 15:  //
                    // histoM02TrueMergedPartConvEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    // histoM02TrueMergedPartConvEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 20:  //
                case 21:  //
                    histoM02TrueGammaFromPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    histoM02TrueGammaFromPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    histoTrueGammaFromPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 22:  //
                    // histoM02TrueGammaFromEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    // histoM02TrueGammaFromEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 30:  //
                case 31:  //
                    histoM02TrueElectronFromPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    histoM02TrueElectronFromPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    histoTrueElectronFromPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 32:  //
                    // histoM02TrueElectronFromEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    // histoM02TrueElectronFromEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                default:
                    histoM02TrueBackgroundperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    histoM02TrueBackgroundperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                }
        }
    }

    //********************************************************************************
    //*                  Writing histograms to outputfile                            *
    //********************************************************************************
    TString fileNameOutput;
    if (kMC){
        fileNameOutput = Form("%s_MC%s.root",nameOutputBase.Data(),currPtHBin>0 ? Form("%d",currPtHBin) : "") ;
    } else {
        fileNameOutput = Form("%s%s_Data.root",nameOutputBase.Data(),currPtHBin>0 ? Form("%d",currPtHBin) : "") ;
    }
    TFile* filePhotoQAWrite = new TFile(fileNameOutput.Data(),"RECREATE");

    for(Int_t i=0;i<20;i++){
        histoM02perSM[i]->Write(Form("histoM02perSM%d",i),TObject::kWriteDelete);
        histoM02perSMvsE[i]->Write(Form("histoM02perSMvsE%d",i),TObject::kWriteDelete);
        histoM20perSM[i]->Write(Form("histoM20perSM%d",i),TObject::kWriteDelete);
        histoM20perSMvsE[i]->Write(Form("histoM20perSMvsE%d",i),TObject::kWriteDelete);
        histoM02vsM20perSM[i]->Write(Form("histoM02vsM20perSM%d",i),TObject::kWriteDelete);
        histoM02vsM20perSMvsE[i]->Write(Form("histoM02vsM20perSMvsE%d",i),TObject::kWriteDelete);
    }

    if (kMC){
        histoEFracFirstLabelvsE->Write("histoEFracFirstLabelvsE",TObject::kWriteDelete);
        histoEFracLeadPi0vsE->Write("histoEFracLeadPi0vsE",TObject::kWriteDelete);
        for(Int_t i=0;i<20;i++){
            // histoTrueMCKindQtPt[i]->Write(Form("histoTrueMCKindQtPt_kind%d",i),TObject::kWriteDelete);
            histoM02TrueMergedPi0perSM[i]->Write(Form("histoM02TrueMergedPi0perSM%d",i),TObject::kWriteDelete);
            // histoM02TrueMergedEtaperSM[i]->Write(Form("histoM02TrueMergedEtaperSM%d",i),TObject::kWriteDelete);
            histoM02TrueMergedPartConvPi0perSM[i]->Write(Form("histoM02TrueMergedPartConvPi0perSM%d",i),TObject::kWriteDelete);
            // histoM02TrueMergedPartConvEtaperSM[i]->Write(Form("histoM02TrueMergedPartConvEtaperSM%d",i),TObject::kWriteDelete);
            histoM02TrueGammaFromPi0perSM[i]->Write(Form("histoM02TrueGammaFromPi0perSM%d",i),TObject::kWriteDelete);
            // histoM02TrueGammaFromEtaperSM[i]->Write(Form("histoM02TrueGammaFromEtaperSM%d",i),TObject::kWriteDelete);
            histoM02TrueElectronFromPi0perSM[i]->Write(Form("histoM02TrueElectronFromPi0perSM%d",i),TObject::kWriteDelete);
            // histoM02TrueElectronFromEtaperSM[i]->Write(Form("histoM02TrueElectronFromEtaperSM%d",i),TObject::kWriteDelete);
            histoM02TrueBackgroundperSM[i]->Write(Form("histoM02TrueBackgroundperSM%d",i),TObject::kWriteDelete);

            histoM02TrueMergedPi0perSMvsE[i]->Write(Form("histoM02TrueMergedPi0perSMvsE%d",i),TObject::kWriteDelete);
            // histoM02TrueMergedEtaperSMvsE[i]->Write(Form("histoM02TrueMergedEtaperSMvsE%d",i),TObject::kWriteDelete);
            histoM02TrueMergedPartConvPi0perSMvsE[i]->Write(Form("histoM02TrueMergedPartConvPi0perSMvsE%d",i),TObject::kWriteDelete);
            // histoM02TrueMergedPartConvEtaperSMvsE[i]->Write(Form("histoM02TrueMergedPartConvEtaperSMvsE%d",i),TObject::kWriteDelete);
            histoM02TrueGammaFromPi0perSMvsE[i]->Write(Form("histoM02TrueGammaFromPi0perSMvsE%d",i),TObject::kWriteDelete);
            // histoM02TrueGammaFromEtaperSMvsE[i]->Write(Form("histoM02TrueGammaFromEtaperSMvsE%d",i),TObject::kWriteDelete);
            histoM02TrueElectronFromPi0perSMvsE[i]->Write(Form("histoM02TrueElectronFromPi0perSMvsE%d",i),TObject::kWriteDelete);
            // histoM02TrueElectronFromEtaperSMvsE[i]->Write(Form("histoM02TrueElectronFromEtaperSMvsE%d",i),TObject::kWriteDelete);
            histoM02TrueBackgroundperSMvsE[i]->Write(Form("histoM02TrueBackgroundperSMvsE%d",i),TObject::kWriteDelete);

            histoTrueMergedPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueMergedPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
            histoTrueMergedPartConvPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueMergedPartConvPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
            histoTrueGammaFromPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueGammaFromPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
            histoTrueElectronFromPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueElectronFromPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
        }
    }
    filePhotoQAWrite->Write();
    filePhotoQAWrite->Close();
    delete filePhotoQAWrite;

    if (kMC){
        delete   histoEFracFirstLabelvsE;
        delete   histoEFracLeadPi0vsE;
    }
    f->Close();
    delete f;

}
