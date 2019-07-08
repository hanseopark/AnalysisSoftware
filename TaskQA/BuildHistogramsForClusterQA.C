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
    Int_t currPtHBin                = -1,
    Bool_t doIsolation              = kTRUE
){
    if(doIsolation)
        cout << "isolating clusters" << endl;
    else
        cout << "using non-isolated clusters" << endl;

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
    Int_t           fBuffer_ClusterNLM;              //!<! array buffer
    Int_t*          fBuffer_ClusterNLM_ID  = new Int_t[kMaxActiveCells]; ;              //!<! array buffer
    Float_t*        fBuffer_ClusterNLM_E  = new Float_t[kMaxActiveCells]; ;              //!<! array buffer
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

    Int_t           fBuffer_Surrounding_NCells;                //!<! array buffer
    Int_t*          fBuffer_Surrounding_Cells_ID           = new Int_t[kMaxActiveCells];               //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_R           = new Float_t[kMaxActiveCells];               //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_E          = new Float_t[kMaxActiveCells];              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_RelativeEta = new Float_t[kMaxActiveCells];              //!<! array buffer
    Float_t*        fBuffer_Surrounding_Cells_RelativePhi = new Float_t[kMaxActiveCells];              //!<! array buffer

    Int_t           fBuffer_Cluster_MC_Label;              //!<! array buffer
    Int_t           fBuffer_Mother_MC_Label;              //!<! array buffer
    Float_t        fBuffer_Cluster_MC_EFracFirstLabel;              //!<! array buffer
    Float_t        fBuffer_Cluster_MC_EFracLeadingPi0;              //!<! array buffer
    Float_t        fBuffer_Cluster_MC_LeadingPi0_Pt;              //!<! array buffer
    Float_t        fBuffer_Cluster_MC_LeadingPi0_E;              //!<! array buffer


    // Set branch adresses
    ClusterQA->SetBranchAddress("Cluster_E",                         &fBuffer_ClusterE);
    ClusterQA->SetBranchAddress("Cluster_Eta",                       &fBuffer_ClusterEta);
    ClusterQA->SetBranchAddress("Cluster_Phi",                       &fBuffer_ClusterPhi);
    ClusterQA->SetBranchAddress("Cluster_IsEMCAL",                   &fBuffer_ClusterIsEMCAL);
    ClusterQA->SetBranchAddress("Cluster_SM",                        &fBuffer_ClusterSupMod);
    ClusterQA->SetBranchAddress("Cluster_NLM",                       &fBuffer_ClusterNLM);
    ClusterQA->SetBranchAddress("Cluster_NLM_ID",                    fBuffer_ClusterNLM_ID);
    ClusterQA->SetBranchAddress("Cluster_NLM_E",                     fBuffer_ClusterNLM_E);
    ClusterQA->SetBranchAddress("Cluster_NumCells",                  &fBuffer_ClusterNumCells);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_ID",            &fBuffer_LeadingCell_ID);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_E",             &fBuffer_LeadingCell_E);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_Eta",           &fBuffer_LeadingCell_Eta);
    ClusterQA->SetBranchAddress("Cluster_LeadingCell_Phi",           &fBuffer_LeadingCell_Phi);
    ClusterQA->SetBranchAddress("Cluster_M02",                       &fBuffer_ClusterM02);
    ClusterQA->SetBranchAddress("Cluster_M20",                       &fBuffer_ClusterM20);
    if(currPtHBin<1 && kMC){
        ClusterQA->SetBranchAddress("Event_Weight",                     &fBuffer_JJWeight);
    }
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

    ClusterQA->SetBranchAddress("Surrounding_NCells",              &fBuffer_Surrounding_NCells);
    ClusterQA->SetBranchAddress("Surrounding_Cells_ID",            fBuffer_Surrounding_Cells_ID);
    ClusterQA->SetBranchAddress("Surrounding_Cells_R",             fBuffer_Surrounding_Cells_R);
    ClusterQA->SetBranchAddress("Surrounding_Cells_E",             fBuffer_Surrounding_Cells_E);
    ClusterQA->SetBranchAddress("Surrounding_Cells_RelativeEta",   fBuffer_Surrounding_Cells_RelativeEta);
    ClusterQA->SetBranchAddress("Surrounding_Cells_RelativePhi",   fBuffer_Surrounding_Cells_RelativePhi);

    if(kMC)
    {
      ClusterQA->SetBranchAddress("Cluster_MC_Label",                &fBuffer_Cluster_MC_Label);
      ClusterQA->SetBranchAddress("Mother_MC_Label",                 &fBuffer_Mother_MC_Label);
      ClusterQA->SetBranchAddress("Cluster_MC_EFracFirstLabel",      &fBuffer_Cluster_MC_EFracFirstLabel);
      ClusterQA->SetBranchAddress("Cluster_MC_EFracLeadingPi0",      &fBuffer_Cluster_MC_EFracLeadingPi0);
      ClusterQA->SetBranchAddress("Cluster_MC_LeadingPi0_Pt",        &fBuffer_Cluster_MC_LeadingPi0_Pt);
      ClusterQA->SetBranchAddress("Cluster_MC_LeadingPi0_E",        &fBuffer_Cluster_MC_LeadingPi0_E);
    }


    Double_t weightsBinsLHC18b8[20]      = { 16.1083,      4.60917,     2.15196,     0.782021,    0.26541,
                                           0.0978374,   0.0294286,   0.00989457,  0.0040615,   0.00135787,
                                           0.000531766, 0.000188772, 9.23331e-05, 4.30245e-05, 2.10196e-05,
                                           1.06695e-05, 5.78742e-06, 3.02897e-06, 1.62702e-06, 2.12118e-06 };
    Double_t weightsBinsLHC18b10[6]      = { 0.000103227,      1.31851e-05,     1.94129e-06,     3.26392e-07,    6.51801e-08,
                                            2.21123e-08 };
    Double_t weightsBinsLHC16c2[20]      = { 28.3084, 8.43277, 4.07753, 1.54359, 0.543318,
                                            0.208394, 0.0652349, 0.0186904, 0.00834528, 0.00301414,
                                            0.00125939, 0.000474403, 0.000244052, 0.00011924, 6.09838e-05,
                                            3.24148e-05, 1.84314e-05, 1.00926e-05, 5.68632e-06, 8.38092e-06};
    Float_t eventNormalization = 1;
    if(kMC){
        TDirectory* fMainEventDir           = (TDirectory*)f->Get(Form("GammaCaloQA_%s",fCutSelection.Data()));
        if(fMainEventDir){
            TList* fMainEventList               = (TList*)fMainEventDir->Get(Form("GammaCaloQA_%s",fCutSelection.Data()));
            if(fMainEventList){
                TList* fEvenCutstList               = (TList*)fMainEventList->FindObject("ConvEventCuts_00010113");
                if(fEvenCutstList){
                    TH1F*  histoVertexZdummy            = (TH1F*) fEvenCutstList->FindObject("VertexZ 00010113");
                    if(histoVertexZdummy){
                        eventNormalization = histoVertexZdummy->GetEntries();
                        cout << "using " << eventNormalization << " events for current pTh bin" << endl;
                    }
                }
            }
        }
        if(currPtHBin>0){
            if(fileName.Contains("LHC18b8")){
                if(currPtHBin>0)
                    fBuffer_JJWeight = weightsBinsLHC18b8[currPtHBin-1];
                fBuffer_JJWeight /= eventNormalization;
            } else if(fileName.Contains("LHC18b10")){
                if(currPtHBin>0)
                    fBuffer_JJWeight = weightsBinsLHC18b10[currPtHBin-1];
                fBuffer_JJWeight /= eventNormalization;
            } else if(fileName.Contains("LHC16c2")){
                if(currPtHBin>0)
                    fBuffer_JJWeight = weightsBinsLHC16c2[currPtHBin-1];
                fBuffer_JJWeight /= eventNormalization;
            }
            cout << "analysing pT hard bin " << currPtHBin << " with weight " << fBuffer_JJWeight << endl;
        }
    }else
        fBuffer_JJWeight = 1;

    //********************************************************************************
    //*            Definition of Boundaries for Histograms                           *
    //********************************************************************************

    //Pt-plots
    Int_t nBinsE                    = 200;
    Double_t firstBinE              = 0.01;
    Double_t lastBinE               = 200.;

    Int_t nBinsM02                   = 300;
    Double_t firstBinM02             = 0.0;
    Double_t lastBinM02              = 1.5;


    //********************************************************************************
    //*      Definition of histograms for reconstructed Conversion Points            *
    //********************************************************************************
    TH1F* histoM02perSM[20]                             = {NULL};
    TH2F* histoM02perSMvsE[20]                          = {NULL};
    TH1F* histoM20perSM[20]                             = {NULL};
    TH2F* histoM20perSMvsE[20]                          = {NULL};
    TH2F* histoM02vsM20perSM[20]                        = {NULL};
    TH3F* histoM02vsM20perSMvsE[20]                     = {NULL};

    TH2F* histoEFracFirstLabelvsE                       = NULL;
    TH2F* histoEFracLeadPi0vsE                          = NULL;

    TH1F* histoM02TrueMergedPi0perSM[20]                = {NULL};
    TH1F* histoM02TrueMergedEtaperSM[20]                = {NULL};
    TH1F* histoM02TrueMergedPartConvPi0perSM[20]        = {NULL};
    TH1F* histoM02TrueMergedPartConvEtaperSM[20]        = {NULL};
    TH1F* histoM02TrueGammaFromPi0perSM[20]             = {NULL};
    TH1F* histoM02TrueGammaFromEtaperSM[20]             = {NULL};
    TH1F* histoM02TrueElectronFromPi0perSM[20]          = {NULL};
    TH1F* histoM02TrueElectronFromEtaperSM[20]          = {NULL};
    TH1F* histoM02TrueBackgroundperSM[20]               = {NULL};

    TH2F* histoM02TrueMergedPi0perSMvsE[20]             = {NULL};
    TH2F* histoM02TrueMergedEtaperSMvsE[20]             = {NULL};
    TH2F* histoM02TrueMergedPartConvPi0perSMvsE[20]     = {NULL};
    TH2F* histoM02TrueMergedPartConvEtaperSMvsE[20]     = {NULL};
    TH2F* histoM02TrueGammaFromPi0perSMvsE[20]          = {NULL};
    TH2F* histoM02TrueGammaFromEtaperSMvsE[20]          = {NULL};
    TH2F* histoM02TrueElectronFromPi0perSMvsE[20]       = {NULL};
    TH2F* histoM02TrueElectronFromEtaperSMvsE[20]       = {NULL};
    TH2F* histoM02TrueBackgroundperSMvsE[20]            = {NULL};


    TH1F* histoFracClusELeadPi0All                      = NULL;
    TH1F* histoFracClusELeadTrueMergedPi0               = NULL;
    TH1F* histoFracClusELeadTrueMergedPartConvPi0       = NULL;
    TH1F* histoFracClusELeadTrueGammaFromPi0            = NULL;

    TH1F* histoIsolationCurveAllClusR005[3]                     = {NULL};
    TH1F* histoIsolationCurveAllClusR01[3]                     = {NULL};
    TH1F* histoIsolationCurveAllClusR04[3]                     = {NULL};
    TH1F* histoIsolationCurveAllPi0ClusR04[3]                     = {NULL};
    TH1F* histoIsolationCurveTrueMergedPi0ClusR04[3]           = {NULL};
    TH1F* histoIsolationCurveTrueMergedPartConvPi0ClusR04[3]   = {NULL};
    TH1F* histoIsolationCurveGammaFromPi0ClusR04[3]            = {NULL};
    TH1F* histoIsolationCurveAllClusR02[3]                     = {NULL};
    TH1F* histoIsolationCurveAllPi0ClusR02[3]                     = {NULL};
    TH1F* histoIsolationCurveTrueMergedPi0ClusR02[3]           = {NULL};
    TH1F* histoIsolationCurveTrueMergedPartConvPi0ClusR02[3]   = {NULL};
    TH1F* histoIsolationCurveGammaFromPi0ClusR02[3]            = {NULL};

    TH2F* histoResolutionAllPi0Clus[4]                     = {NULL};
    TH2F* histoResolutionTrueMergedPi0Clus[4]              = {NULL};
    TH2F* histoResolutionTrueMergedPartConvPi0Clus[4]      = {NULL};
    TH2F* histoResolutionTrueGammaFromPi0Clus[4]           = {NULL};
    TH2F* histoResoMatrixAllPi0Clus[4]                     = {NULL};
    TH2F* histoResoMatrixTrueMergedPi0Clus[4]              = {NULL};
    TH2F* histoResoMatrixTrueMergedPartConvPi0Clus[4]      = {NULL};
    TH2F* histoResoMatrixTrueGammaFromPi0Clus[4]           = {NULL};

    TH3F* histoTrueMergedPi0M02vsM20perSMvsE[20]        = {NULL};
    TH3F* histoTrueMergedPartConvPi0M02vsM20perSMvsE[20]= {NULL};
    TH3F* histoTrueGammaFromPi0M02vsM20perSMvsE[20]     = {NULL};
    TH3F* histoTrueElectronFromPi0M02vsM20perSMvsE[20]  = {NULL};

    for(Int_t i=0;i<20;i++){
        // histoM02perSM[i]                = new TH1F(Form("histoM02perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
        // histoM02perSMvsE[i]             = new TH2F(Form("histoM02perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        // SetLogBinningTH2(histoM02perSMvsE[i]);
        // histoM20perSM[i]                = new TH1F(Form("histoM20perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
        // histoM20perSMvsE[i]             = new TH2F(Form("histoM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        // SetLogBinningTH2(histoM20perSMvsE[i]);
        // histoM02vsM20perSM[i]           = new TH2F(Form("histoM02vsM20perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02);
        // histoM02vsM20perSMvsE[i]        = new TH3F(Form("histoM02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        // SetLogBinningTH3(histoM02vsM20perSMvsE[i]);
    }
    for(Int_t i=0;i<3;i++){
        histoIsolationCurveAllClusR005[i]                   = new TH1F(Form("histoIsolationCurveAllClusR005_iso%d",i),"", 3000, 0, 300);
        histoIsolationCurveAllClusR01[i]                    = new TH1F(Form("histoIsolationCurveAllClusR01_iso%d",i),"", 3000, 0, 300);
        histoIsolationCurveAllClusR02[i]                    = new TH1F(Form("histoIsolationCurveAllClusR02_iso%d",i),"", 3000, 0, 300);
        histoIsolationCurveAllClusR04[i]                    = new TH1F(Form("histoIsolationCurveAllClusR04_iso%d",i),"", 3000, 0, 300);

    }
    if (kMC){
        histoEFracFirstLabelvsE                     = new TH2F("histoEFracFirstLabelvsE","",20, 0., 1., nBinsE, firstBinE, lastBinE);
        histoEFracLeadPi0vsE                        = new TH2F("histoEFracLeadPi0vsE","",100, 0., 1., nBinsE, firstBinE, lastBinE);

        histoFracClusELeadPi0All                    = new TH1F(Form("histoFracClusELeadPi0All"),"",100, -0.005, 1.005);
        histoFracClusELeadTrueMergedPi0             = new TH1F(Form("histoFracClusELeadTrueMergedPi0"),"",100, -0.005, 1.005);
        histoFracClusELeadTrueMergedPartConvPi0     = new TH1F(Form("histoFracClusELeadTrueMergedPartConvPi0"),"",100, -0.005, 1.005);
        histoFracClusELeadTrueGammaFromPi0          = new TH1F(Form("histoFracClusELeadTrueGammaFromPi0"),"",100, -0.005, 1.005);
        for(Int_t i=0;i<3;i++){
            histoIsolationCurveAllPi0ClusR04[i]                 = new TH1F(Form("histoIsolationCurveAllPi0ClusR04_iso%d",i),"", 1500, 0, 300);
            histoIsolationCurveTrueMergedPi0ClusR04[i]          = new TH1F(Form("histoIsolationCurveTrueMergedPi0ClusR04_iso%d",i),"", 1500, 0, 300);
            histoIsolationCurveTrueMergedPartConvPi0ClusR04[i]  = new TH1F(Form("histoIsolationCurveTrueMergedPartConvPi0ClusR04_iso%d",i),"", 1500, 0, 300);
            histoIsolationCurveGammaFromPi0ClusR04[i]           = new TH1F(Form("histoIsolationCurveGammaFromPi0ClusR04_iso%d",i),"", 1500, 0, 300);
            histoIsolationCurveAllPi0ClusR02[i]                 = new TH1F(Form("histoIsolationCurveAllPi0ClusR02_iso%d",i),"", 1500, 0, 300);
            histoIsolationCurveTrueMergedPi0ClusR02[i]          = new TH1F(Form("histoIsolationCurveTrueMergedPi0ClusR02_iso%d",i),"", 1500, 0, 300);
            histoIsolationCurveTrueMergedPartConvPi0ClusR02[i]  = new TH1F(Form("histoIsolationCurveTrueMergedPartConvPi0ClusR02_iso%d",i),"", 1500, 0, 300);
            histoIsolationCurveGammaFromPi0ClusR02[i]           = new TH1F(Form("histoIsolationCurveGammaFromPi0ClusR02_iso%d",i),"", 1500, 0, 300);
        }
        for(Int_t i=0;i<4;i++){
            histoResolutionAllPi0Clus[i]                    = new TH2F(Form("histoResolutionAllPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE,200,-1,1);
            histoResolutionTrueMergedPi0Clus[i]             = new TH2F(Form("histoResolutionTrueMergedPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE,200,-1,1);
            histoResolutionTrueMergedPartConvPi0Clus[i]     = new TH2F(Form("histoResolutionTrueMergedPartConvPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE,200,-1,1);
            histoResolutionTrueGammaFromPi0Clus[i]          = new TH2F(Form("histoResolutionTrueGammaFromPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE,200,-1,1);

            histoResoMatrixAllPi0Clus[i]                    = new TH2F(Form("histoResoMatrixAllPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE, nBinsE, firstBinE, lastBinE);
            histoResoMatrixTrueMergedPi0Clus[i]             = new TH2F(Form("histoResoMatrixTrueMergedPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE, nBinsE, firstBinE, lastBinE);
            histoResoMatrixTrueMergedPartConvPi0Clus[i]     = new TH2F(Form("histoResoMatrixTrueMergedPartConvPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE, nBinsE, firstBinE, lastBinE);
            histoResoMatrixTrueGammaFromPi0Clus[i]          = new TH2F(Form("histoResoMatrixTrueGammaFromPi0Clus_iso%d",i),"", nBinsE, firstBinE, lastBinE, nBinsE, firstBinE, lastBinE);
        }
        for(Int_t i=0;i<20;i++){
            // histoM02TrueMergedPi0perSM[i]                   = new TH1F(Form("histoM02TrueMergedPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueMergedEtaperSM[i]                   = new TH1F(Form("histoM02TrueMergedEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueMergedPartConvPi0perSM[i]           = new TH1F(Form("histoM02TrueMergedPartConvPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueMergedPartConvEtaperSM[i]           = new TH1F(Form("histoM02TrueMergedPartConvEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueGammaFromPi0perSM[i]                = new TH1F(Form("histoM02TrueGammaFromPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueGammaFromEtaperSM[i]                = new TH1F(Form("histoM02TrueGammaFromEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueElectronFromPi0perSM[i]             = new TH1F(Form("histoM02TrueElectronFromPi0perSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueElectronFromEtaperSM[i]             = new TH1F(Form("histoM02TrueElectronFromEtaperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);
            // histoM02TrueBackgroundperSM[i]                  = new TH1F(Form("histoM02TrueBackgroundperSM%d",i),"",nBinsM02, firstBinM02, lastBinM02);

            // histoM02TrueMergedPi0perSMvsE[i]                = new TH2F(Form("histoM02TrueMergedPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueMergedEtaperSMvsE[i]                = new TH2F(Form("histoM02TrueMergedEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueMergedPartConvPi0perSMvsE[i]        = new TH2F(Form("histoM02TrueMergedPartConvPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueMergedPartConvEtaperSMvsE[i]        = new TH2F(Form("histoM02TrueMergedPartConvEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueGammaFromPi0perSMvsE[i]             = new TH2F(Form("histoM02TrueGammaFromPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueGammaFromEtaperSMvsE[i]             = new TH2F(Form("histoM02TrueGammaFromEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueElectronFromPi0perSMvsE[i]          = new TH2F(Form("histoM02TrueElectronFromPi0perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueElectronFromEtaperSMvsE[i]          = new TH2F(Form("histoM02TrueElectronFromEtaperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoM02TrueBackgroundperSMvsE[i]               = new TH2F(Form("histoM02TrueBackgroundperSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);

            // histoTrueMergedPi0M02vsM20perSMvsE[i]           = new TH3F(Form("histoTrueMergedPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoTrueMergedPartConvPi0M02vsM20perSMvsE[i]   = new TH3F(Form("histoTrueMergedPartConvPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoTrueGammaFromPi0M02vsM20perSMvsE[i]        = new TH3F(Form("histoTrueGammaFromPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
            // histoTrueElectronFromPi0M02vsM20perSMvsE[i]     = new TH3F(Form("histoTrueElectronFromPi0M02vsM20perSMvsE%d",i),"",nBinsM02, firstBinM02, lastBinM02,nBinsM02, firstBinM02, lastBinM02, nBinsE, firstBinE, lastBinE);
        }
    }

    for(Int_t i=0;i<20;i++){
        if(histoM02perSM[i])histoM02perSM[i]->Sumw2();
        if(histoM02perSMvsE[i])histoM02perSMvsE[i]->Sumw2();
        if(histoM20perSM[i])histoM20perSM[i]->Sumw2();
        if(histoM20perSMvsE[i])histoM20perSMvsE[i]->Sumw2();
        if(histoM02vsM20perSM[i])histoM02vsM20perSM[i]->Sumw2();
        if(histoM02vsM20perSMvsE[i])histoM02vsM20perSMvsE[i]->Sumw2();
    }
    if (kMC){
        if(histoEFracFirstLabelvsE)histoEFracFirstLabelvsE->Sumw2();
        if(histoEFracLeadPi0vsE)histoEFracLeadPi0vsE->Sumw2();
        for(Int_t i=0;i<4;i++){
            if(histoResolutionAllPi0Clus[i])histoResolutionAllPi0Clus[i]->Sumw2();
            if(histoResolutionTrueMergedPi0Clus[i])histoResolutionTrueMergedPi0Clus[i]->Sumw2();
            if(histoResolutionTrueMergedPartConvPi0Clus[i])histoResolutionTrueMergedPartConvPi0Clus[i]->Sumw2();
            if(histoResolutionTrueGammaFromPi0Clus[i])histoResolutionTrueGammaFromPi0Clus[i]->Sumw2();

            if(histoResoMatrixAllPi0Clus[i])histoResoMatrixAllPi0Clus[i]->Sumw2();
            if(histoResoMatrixTrueMergedPi0Clus[i])histoResoMatrixTrueMergedPi0Clus[i]->Sumw2();
            if(histoResoMatrixTrueMergedPartConvPi0Clus[i])histoResoMatrixTrueMergedPartConvPi0Clus[i]->Sumw2();
            if(histoResoMatrixTrueGammaFromPi0Clus[i])histoResoMatrixTrueGammaFromPi0Clus[i]->Sumw2();
        }
        for(Int_t i=0;i<20;i++){
            if(histoM02TrueMergedPi0perSM[i])histoM02TrueMergedPi0perSM[i]->Sumw2();
            if(histoM02TrueMergedEtaperSM[i])histoM02TrueMergedEtaperSM[i]->Sumw2();
            if(histoM02TrueMergedPartConvPi0perSM[i])histoM02TrueMergedPartConvPi0perSM[i]->Sumw2();
            if(histoM02TrueMergedPartConvEtaperSM[i])histoM02TrueMergedPartConvEtaperSM[i]->Sumw2();
            if(histoM02TrueGammaFromPi0perSM[i])histoM02TrueGammaFromPi0perSM[i]->Sumw2();
            if(histoM02TrueGammaFromEtaperSM[i])histoM02TrueGammaFromEtaperSM[i]->Sumw2();
            if(histoM02TrueElectronFromPi0perSM[i])histoM02TrueElectronFromPi0perSM[i]->Sumw2();
            if(histoM02TrueElectronFromEtaperSM[i])histoM02TrueElectronFromEtaperSM[i]->Sumw2();
            if(histoM02TrueBackgroundperSM[i])histoM02TrueBackgroundperSM[i]->Sumw2();

            if(histoM02TrueMergedPi0perSMvsE[i])histoM02TrueMergedPi0perSMvsE[i]->Sumw2();
            if(histoM02TrueMergedEtaperSMvsE[i])histoM02TrueMergedEtaperSMvsE[i]->Sumw2();
            if(histoM02TrueMergedPartConvPi0perSMvsE[i])histoM02TrueMergedPartConvPi0perSMvsE[i]->Sumw2();
            if(histoM02TrueMergedPartConvEtaperSMvsE[i])histoM02TrueMergedPartConvEtaperSMvsE[i]->Sumw2();
            if(histoM02TrueGammaFromPi0perSMvsE[i])histoM02TrueGammaFromPi0perSMvsE[i]->Sumw2();
            if(histoM02TrueGammaFromEtaperSMvsE[i])histoM02TrueGammaFromEtaperSMvsE[i]->Sumw2();
            if(histoM02TrueElectronFromPi0perSMvsE[i])histoM02TrueElectronFromPi0perSMvsE[i]->Sumw2();
            if(histoM02TrueElectronFromEtaperSMvsE[i])histoM02TrueElectronFromEtaperSMvsE[i]->Sumw2();
            if(histoM02TrueBackgroundperSMvsE[i])histoM02TrueBackgroundperSMvsE[i]->Sumw2();

            if(histoTrueMergedPi0M02vsM20perSMvsE[i])histoTrueMergedPi0M02vsM20perSMvsE[i]->Sumw2();
            if(histoTrueMergedPartConvPi0M02vsM20perSMvsE[i])histoTrueMergedPartConvPi0M02vsM20perSMvsE[i]->Sumw2();
            if(histoTrueGammaFromPi0M02vsM20perSMvsE[i])histoTrueGammaFromPi0M02vsM20perSMvsE[i]->Sumw2();
            if(histoTrueElectronFromPi0M02vsM20perSMvsE[i])histoTrueElectronFromPi0M02vsM20perSMvsE[i]->Sumw2();
        }
    }

    //********************************************************************************
    //*      Reading of Tree with reconstructed gammas/ filling of histograms        *
    //********************************************************************************

    Long64_t nEntriesRecClus                 = ClusterQA->GetEntries();
    if(maxNumClusProcess>0 && maxNumClusProcess<nEntriesRecClus)
        nEntriesRecClus = maxNumClusProcess;
    Int_t fractionToAnalyse = 1;
    cout << "Number of Clusters to be processed: " << nEntriesRecClus/fractionToAnalyse << endl;
    Bool_t isIsolated[12]={kFALSE}; // 0: charged, 1: neutral, 2: charged+neutral
    for (Long64_t i=0; i<nEntriesRecClus;i++) {

        if(i%fractionToAnalyse > 0) continue;
        // cout << "cluster" << i << endl;
        if(i>0 && i%(nEntriesRecClus/(10*fractionToAnalyse)) ==0) cout << "//processed " << 100*(fractionToAnalyse*i)/nEntriesRecClus << "%"  << endl;
        ClusterQA->GetEntry(i);

        // minimun energy requirement for clusters
        // if(fBuffer_ClusterE<10) continue;

        // remove exotics with M02>0.1
        if(fBuffer_ClusterM02<0.1) continue;
        // set analysis M02 requirement
        // if(fBuffer_ClusterM02<0.27) continue;

        // if(0 && kMC && (fBuffer_Cluster_MC_Label==10||fBuffer_Cluster_MC_Label==11||fBuffer_Cluster_MC_Label==13||fBuffer_Cluster_MC_Label==14||fBuffer_Cluster_MC_Label==20||fBuffer_Cluster_MC_Label==21||fBuffer_Cluster_MC_Label==30||fBuffer_Cluster_MC_Label==31) && fBuffer_ClusterE > 60){
        if(kMC){
            cout<< "//_________________________________________________"<<endl;
            cout << "//MCLabel: " << fBuffer_Cluster_MC_Label << endl;
            cout << "//EvtMult: " << fBuffer_Event_Multiplicity << endl;
            cout << "//E: " << fBuffer_ClusterE << "\tM02: " << fBuffer_ClusterM02 << endl;
            cout << "//TruePi0E: " << fBuffer_Cluster_MC_LeadingPi0_E << "\tTruePi0Pt: " << fBuffer_Cluster_MC_LeadingPi0_Pt << endl;

            cout << "// surrounding cells for this cluster" << endl;
            cout << "\tconst Int_t ncellCLSSurr" << i << "=" <<fBuffer_Surrounding_NCells<<";"<<endl;
            cout << "\tInt_t IDsClsSurr" << i << "[ncellCLSSurr" << i << "] = {";
            for(Int_t cellI=0; cellI<fBuffer_Surrounding_NCells; cellI++){
                if(cellI==0)cout << fBuffer_Surrounding_Cells_ID[cellI];
                else cout << ", " << fBuffer_Surrounding_Cells_ID[cellI];
            }
            cout << "};"<<endl;
            cout << "\tDouble_t EClsSurr" << i << "[ncellCLSSurr" << i << "] = {";
            for(Int_t cellI=0; cellI<fBuffer_Surrounding_NCells; cellI++){
                if(cellI==0)cout << fBuffer_Surrounding_Cells_E[cellI];
                else cout << ", " << fBuffer_Surrounding_Cells_E[cellI];
            }
            cout << "};"<<endl;
            cout << "\tFill_histo(ncellCLSSurr" << i << ",IDsClsSurr" << i << ",EClsSurr" << i << ",hBadChannels,fEMCALGeo,fEMCALRecoUtils);"<<endl;

            cout << "// cells of this cluster" << endl;
            cout << "\tconst Int_t ncellCLS" << i << "=" <<fBuffer_ClusterNumCells<<";"<<endl;
            cout << "\tInt_t IDsCls" << i << "[ncellCLS" << i << "] = {";
            for(Int_t cellI=0; cellI<fBuffer_ClusterNumCells; cellI++){
                if(cellI==0)cout << fBuffer_Input_Cells_ID[cellI];
                else cout << ", " << fBuffer_Input_Cells_ID[cellI];
            }
            cout << "};"<<endl;
            cout << "\tDouble_t ECls" << i << "[ncellCLS" << i << "] = {";
            for(Int_t cellI=0; cellI<fBuffer_ClusterNumCells; cellI++){
                if(cellI==0)cout << fBuffer_Input_Cells_E[cellI];
                else cout << ", " << fBuffer_Input_Cells_E[cellI];
            }
            cout << "};"<<endl;
            cout << "\tFill_histo(ncellCLS" << i << ",IDsCls" << i << ",ECls" << i << ",hBadChannels,fEMCALGeo,fEMCALRecoUtils);"<<endl;

            cout << "// tracks around cluster" << endl;
            cout << "\tconst Int_t nTracksSur" << i << "=" <<fBuffer_Surrounding_NTracks<<";"<<endl;
            cout << "\tDouble_t EtaTracksSur" << i << "[nTracksSur" << i << "] = {";
            for(Int_t trackI=0; trackI<fBuffer_Surrounding_NTracks; trackI++){
                if(trackI==0)cout << fBuffer_LeadingCell_Eta+fBuffer_Surrounding_Tracks_RelativeEta[trackI];
                else cout << ", " << fBuffer_LeadingCell_Eta+fBuffer_Surrounding_Tracks_RelativeEta[trackI];
            }
            cout << "};"<<endl;
            cout << "\tDouble_t PhiTracksSur" << i << "[nTracksSur" << i << "] = {";
            for(Int_t trackI=0; trackI<fBuffer_Surrounding_NTracks; trackI++){
                if(trackI==0)cout << fBuffer_LeadingCell_Phi+fBuffer_Surrounding_Tracks_RelativePhi[trackI];
                else cout << ", " << fBuffer_LeadingCell_Phi+fBuffer_Surrounding_Tracks_RelativePhi[trackI];
            }
            cout << "};"<<endl;
            cout << "\tDouble_t ETracksSur" << i << "[nTracksSur" << i << "] = {";
            for(Int_t trackI=0; trackI<fBuffer_Surrounding_NTracks; trackI++){
                if(trackI==0)cout << fBuffer_Surrounding_Tracks_P[trackI];
                else cout << ", " << fBuffer_Surrounding_Tracks_P[trackI];
            }
            cout << "};"<<endl;
            cout << "\tFill_histoTracks(nTracksSur" << i << ",EtaTracksSur" << i << ",PhiTracksSur" << i << ",ETracksSur" << i << ",hBadChannels,fEMCALGeo,fEMCALRecoUtils);"<<endl;
            cout<< "//_________________________________________________"<<endl;
        }

            // cout << "analysing pT hard bin with weight " << fBuffer_JJWeight << endl;

        // apply isolation charged and/or neutral
        Float_t isoValuesR005[3]      = {0};
        Float_t isoTrackPtR005      = 0;
        Float_t isoClusterEtR005    = 0;
        Float_t isoValuesR01[3]      = {0};
        Float_t isoTrackPtR01      = 0;
        Float_t isoClusterEtR01    = 0;
        Float_t isoValuesR02[3]      = {0};
        Float_t isoTrackPtR02      = 0;
        Float_t isoClusterEtR02    = 0;
        Float_t isoValuesR04[3]      = {0};
        Float_t isoTrackPtR04      = 0;
        Float_t isoClusterEtR04    = 0;
        if(doIsolation){
            for(Int_t itr=0; itr<fBuffer_Surrounding_NTracks; itr++){
                if(fBuffer_Surrounding_Tracks_R[itr] < TMath::Power(0.05,2))
                    isoTrackPtR005+=fBuffer_Surrounding_Tracks_Pt[itr];
                if(fBuffer_Surrounding_Tracks_R[itr] < TMath::Power(0.1,2))
                    isoTrackPtR01+=fBuffer_Surrounding_Tracks_Pt[itr];
                if(fBuffer_Surrounding_Tracks_R[itr] < TMath::Power(0.2,2))
                    isoTrackPtR02+=fBuffer_Surrounding_Tracks_Pt[itr];
                if(fBuffer_Surrounding_Tracks_R[itr] < TMath::Power(0.4,2))
                    isoTrackPtR04+=fBuffer_Surrounding_Tracks_Pt[itr];
            }
            for(Int_t isurcell=0; isurcell<fBuffer_Surrounding_NCells; isurcell++){
                Bool_t cellInCls = kFALSE;
                for(Int_t icellcls=0; icellcls<fBuffer_ClusterNumCells; icellcls++){
                    if(fBuffer_Input_Cells_ID[icellcls]==fBuffer_Surrounding_Cells_ID[isurcell])
                        cellInCls = kTRUE;
                }
                if(!cellInCls){
                    if(fBuffer_Surrounding_Cells_R[isurcell] < TMath::Power(0.05,2))
                        isoClusterEtR005+=fBuffer_Surrounding_Cells_E[isurcell];
                    if(fBuffer_Surrounding_Cells_R[isurcell] < TMath::Power(0.1,2))
                        isoClusterEtR01+=fBuffer_Surrounding_Cells_E[isurcell];
                    if(fBuffer_Surrounding_Cells_R[isurcell] < TMath::Power(0.2,2))
                        isoClusterEtR02+=fBuffer_Surrounding_Cells_E[isurcell];
                    if(fBuffer_Surrounding_Cells_R[isurcell] < TMath::Power(0.4,2))
                        isoClusterEtR04+=fBuffer_Surrounding_Cells_E[isurcell];
                }
            }
            isoTrackPtR04<2.0                  ? isIsolated[0] = kTRUE :   isIsolated[0] = kFALSE;
            isoClusterEtR04<2.0                ? isIsolated[1] = kTRUE :   isIsolated[1] = kFALSE;
            (isoClusterEtR04+isoTrackPtR04)<2.0   ? isIsolated[2] = kTRUE :   isIsolated[2] = kFALSE;

            isoValuesR005[0] = isoTrackPtR005;
            isoValuesR005[1] = isoClusterEtR005;
            isoValuesR005[2] = isoClusterEtR005+isoTrackPtR005;
            isoValuesR01[0] = isoTrackPtR01;
            isoValuesR01[1] = isoClusterEtR01;
            isoValuesR01[2] = isoClusterEtR01+isoTrackPtR01;
            isoValuesR02[0] = isoTrackPtR02;
            isoValuesR02[1] = isoClusterEtR02;
            isoValuesR02[2] = isoClusterEtR02+isoTrackPtR02;
            isoValuesR04[0] = isoTrackPtR04;
            isoValuesR04[1] = isoClusterEtR04;
            isoValuesR04[2] = isoClusterEtR04+isoTrackPtR04;
        }
        if(histoM02perSM[fBuffer_ClusterSupMod])histoM02perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
        if(histoM02perSMvsE[fBuffer_ClusterSupMod])histoM02perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
        if(histoM20perSM[fBuffer_ClusterSupMod])histoM20perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM20,fBuffer_JJWeight);
        if(histoM20perSMvsE[fBuffer_ClusterSupMod])histoM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
        if(histoM02vsM20perSM[fBuffer_ClusterSupMod])histoM02vsM20perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_JJWeight);
        if(histoM02vsM20perSMvsE[fBuffer_ClusterSupMod])histoM02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
        for(Int_t iIso=0;iIso<3;iIso++){
            if(histoIsolationCurveAllClusR005[iIso])histoIsolationCurveAllClusR005[iIso]->Fill(isoValuesR005[iIso],fBuffer_JJWeight);
            if(histoIsolationCurveAllClusR01[iIso])histoIsolationCurveAllClusR01[iIso]->Fill(isoValuesR01[iIso],fBuffer_JJWeight);
            if(histoIsolationCurveAllClusR02[iIso])histoIsolationCurveAllClusR02[iIso]->Fill(isoValuesR02[iIso],fBuffer_JJWeight);
            if(histoIsolationCurveAllClusR04[iIso])histoIsolationCurveAllClusR04[iIso]->Fill(isoValuesR04[iIso],fBuffer_JJWeight);
        }
        if (kMC){
            if(histoEFracFirstLabelvsE)histoEFracFirstLabelvsE->Fill(fBuffer_Cluster_MC_EFracFirstLabel,fBuffer_ClusterE,fBuffer_JJWeight);
            if(histoEFracLeadPi0vsE)histoEFracLeadPi0vsE->Fill(fBuffer_Cluster_MC_EFracLeadingPi0,fBuffer_ClusterE,fBuffer_JJWeight);

            if(fBuffer_Cluster_MC_EFracLeadingPi0 > fBuffer_Cluster_MC_EFracFirstLabel){
                if(histoFracClusELeadPi0All)histoFracClusELeadPi0All->Fill(fBuffer_Cluster_MC_EFracLeadingPi0,fBuffer_JJWeight);

                if(histoResolutionAllPi0Clus[0])histoResolutionAllPi0Clus[0]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                if(histoResoMatrixAllPi0Clus[0])histoResoMatrixAllPi0Clus[0]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);

                for(Int_t iIso=1;iIso<4;iIso++){
                    if(isIsolated[iIso-1]){
                        if(histoResolutionAllPi0Clus[iIso])histoResolutionAllPi0Clus[iIso]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                        if(histoResoMatrixAllPi0Clus[iIso])histoResoMatrixAllPi0Clus[iIso]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);
                    }
                }
                for(Int_t iIso=0;iIso<3;iIso++){
                    if(histoIsolationCurveAllPi0ClusR02[iIso])histoIsolationCurveAllPi0ClusR02[iIso]->Fill(isoValuesR02[iIso],fBuffer_JJWeight);
                    if(histoIsolationCurveAllPi0ClusR04[iIso])histoIsolationCurveAllPi0ClusR04[iIso]->Fill(isoValuesR04[iIso],fBuffer_JJWeight);
                }
            }
            switch(fBuffer_Cluster_MC_Label){
                case 10:  //
                case 11:  //
                    if(histoM02TrueMergedPi0perSM[fBuffer_ClusterSupMod])histoM02TrueMergedPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueMergedPi0perSMvsE[fBuffer_ClusterSupMod])histoM02TrueMergedPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    if(histoTrueMergedPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod])histoTrueMergedPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    if(fBuffer_Cluster_MC_EFracLeadingPi0 > fBuffer_Cluster_MC_EFracFirstLabel){
                        if(histoFracClusELeadTrueMergedPi0)histoFracClusELeadTrueMergedPi0->Fill(fBuffer_Cluster_MC_EFracLeadingPi0,fBuffer_JJWeight);
                        if(histoResolutionTrueMergedPi0Clus[0])histoResolutionTrueMergedPi0Clus[0]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                        if(histoResoMatrixTrueMergedPi0Clus[0])histoResoMatrixTrueMergedPi0Clus[0]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);
                        for(Int_t iIso=1;iIso<4;iIso++){
                            if(isIsolated[iIso-1]){
                                if(histoResolutionTrueMergedPi0Clus[iIso])histoResolutionTrueMergedPi0Clus[iIso]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                                if(histoResoMatrixTrueMergedPi0Clus[iIso])histoResoMatrixTrueMergedPi0Clus[iIso]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);
                            }
                        }
                        for(Int_t iIso=0;iIso<3;iIso++){
                            if(histoIsolationCurveTrueMergedPi0ClusR02[iIso])histoIsolationCurveTrueMergedPi0ClusR02[iIso]->Fill(isoValuesR02[iIso],fBuffer_JJWeight);
                            if(histoIsolationCurveTrueMergedPi0ClusR04[iIso])histoIsolationCurveTrueMergedPi0ClusR04[iIso]->Fill(isoValuesR04[iIso],fBuffer_JJWeight);
                        }
                    }
                    break;
                case 12:  //
                    if(histoM02TrueMergedEtaperSM[fBuffer_ClusterSupMod])histoM02TrueMergedEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueMergedEtaperSMvsE[fBuffer_ClusterSupMod])histoM02TrueMergedEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 13:  //
                case 14:  //
                    if(histoM02TrueMergedPartConvPi0perSM[fBuffer_ClusterSupMod])histoM02TrueMergedPartConvPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueMergedPartConvPi0perSMvsE[fBuffer_ClusterSupMod])histoM02TrueMergedPartConvPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    if(histoTrueMergedPartConvPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod])histoTrueMergedPartConvPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    if(fBuffer_Cluster_MC_EFracLeadingPi0 > fBuffer_Cluster_MC_EFracFirstLabel){
                        if(histoFracClusELeadTrueMergedPartConvPi0)histoFracClusELeadTrueMergedPartConvPi0->Fill(fBuffer_Cluster_MC_EFracLeadingPi0,fBuffer_JJWeight);
                        if(histoResolutionTrueMergedPartConvPi0Clus[0])histoResolutionTrueMergedPartConvPi0Clus[0]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                        if(histoResoMatrixTrueMergedPartConvPi0Clus[0])histoResoMatrixTrueMergedPartConvPi0Clus[0]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);
                        for(Int_t iIso=1;iIso<4;iIso++){
                            if(isIsolated[iIso-1]){
                                if(histoResolutionTrueMergedPartConvPi0Clus[iIso])histoResolutionTrueMergedPartConvPi0Clus[iIso]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                                if(histoResoMatrixTrueMergedPartConvPi0Clus[iIso])histoResoMatrixTrueMergedPartConvPi0Clus[iIso]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);
                            }
                        }
                        for(Int_t iIso=0;iIso<3;iIso++){
                            if(histoIsolationCurveTrueMergedPartConvPi0ClusR02[iIso])histoIsolationCurveTrueMergedPartConvPi0ClusR02[iIso]->Fill(isoValuesR02[iIso],fBuffer_JJWeight);
                            if(histoIsolationCurveTrueMergedPartConvPi0ClusR04[iIso])histoIsolationCurveTrueMergedPartConvPi0ClusR04[iIso]->Fill(isoValuesR04[iIso],fBuffer_JJWeight);
                        }
                    }
                    break;
                case 15:  //
                    if(histoM02TrueMergedPartConvEtaperSM[fBuffer_ClusterSupMod])histoM02TrueMergedPartConvEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueMergedPartConvEtaperSMvsE[fBuffer_ClusterSupMod])histoM02TrueMergedPartConvEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 20:  //
                case 21:  //
                    if(histoM02TrueGammaFromPi0perSM[fBuffer_ClusterSupMod])histoM02TrueGammaFromPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueGammaFromPi0perSMvsE[fBuffer_ClusterSupMod])histoM02TrueGammaFromPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    if(histoTrueGammaFromPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod])histoTrueGammaFromPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    if(fBuffer_Cluster_MC_EFracLeadingPi0 > fBuffer_Cluster_MC_EFracFirstLabel){
                        if(histoFracClusELeadTrueGammaFromPi0)histoFracClusELeadTrueGammaFromPi0->Fill(fBuffer_Cluster_MC_EFracLeadingPi0,fBuffer_JJWeight);
                        if(histoResolutionTrueGammaFromPi0Clus[0])histoResolutionTrueGammaFromPi0Clus[0]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                        if(histoResoMatrixTrueGammaFromPi0Clus[0])histoResoMatrixTrueGammaFromPi0Clus[0]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);
                        for(Int_t iIso=1;iIso<4;iIso++){
                            if(isIsolated[iIso-1]){
                                if(histoResolutionTrueGammaFromPi0Clus[iIso])histoResolutionTrueGammaFromPi0Clus[iIso]->Fill(fBuffer_Cluster_MC_LeadingPi0_Pt,(fBuffer_ClusterE-fBuffer_Cluster_MC_LeadingPi0_Pt)/(fBuffer_Cluster_MC_LeadingPi0_Pt),fBuffer_JJWeight);
                                if(histoResoMatrixTrueGammaFromPi0Clus[iIso])histoResoMatrixTrueGammaFromPi0Clus[iIso]->Fill(fBuffer_ClusterE,fBuffer_Cluster_MC_LeadingPi0_Pt,fBuffer_JJWeight);
                            }
                        }
                        for(Int_t iIso=0;iIso<3;iIso++){
                            if(histoIsolationCurveGammaFromPi0ClusR02[iIso])histoIsolationCurveGammaFromPi0ClusR02[iIso]->Fill(isoValuesR02[iIso],fBuffer_JJWeight);
                            if(histoIsolationCurveGammaFromPi0ClusR04[iIso])histoIsolationCurveGammaFromPi0ClusR04[iIso]->Fill(isoValuesR04[iIso],fBuffer_JJWeight);
                        }
                    }
                    break;
                case 22:  //
                    if(histoM02TrueGammaFromEtaperSM[fBuffer_ClusterSupMod])histoM02TrueGammaFromEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueGammaFromEtaperSMvsE[fBuffer_ClusterSupMod])histoM02TrueGammaFromEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 30:  //
                case 31:  //
                    if(histoM02TrueElectronFromPi0perSM[fBuffer_ClusterSupMod])histoM02TrueElectronFromPi0perSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueElectronFromPi0perSMvsE[fBuffer_ClusterSupMod])histoM02TrueElectronFromPi0perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    if(histoTrueElectronFromPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod])histoTrueElectronFromPi0M02vsM20perSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterM20,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                case 32:  //
                    if(histoM02TrueElectronFromEtaperSM[fBuffer_ClusterSupMod])histoM02TrueElectronFromEtaperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueElectronFromEtaperSMvsE[fBuffer_ClusterSupMod])histoM02TrueElectronFromEtaperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
                    break;
                default:
                    if(histoM02TrueBackgroundperSM[fBuffer_ClusterSupMod])histoM02TrueBackgroundperSM[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_JJWeight);
                    if(histoM02TrueBackgroundperSMvsE[fBuffer_ClusterSupMod])histoM02TrueBackgroundperSMvsE[fBuffer_ClusterSupMod]->Fill(fBuffer_ClusterM02,fBuffer_ClusterE,fBuffer_JJWeight);
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
        if(histoM02perSM[i])histoM02perSM[i]->Write(Form("histoM02perSM%d",i),TObject::kWriteDelete);
        if(histoM02perSMvsE[i])histoM02perSMvsE[i]->Write(Form("histoM02perSMvsE%d",i),TObject::kWriteDelete);
        if(histoM20perSM[i])histoM20perSM[i]->Write(Form("histoM20perSM%d",i),TObject::kWriteDelete);
        if(histoM20perSMvsE[i])histoM20perSMvsE[i]->Write(Form("histoM20perSMvsE%d",i),TObject::kWriteDelete);
        if(histoM02vsM20perSM[i])histoM02vsM20perSM[i]->Write(Form("histoM02vsM20perSM%d",i),TObject::kWriteDelete);
        if(histoM02vsM20perSMvsE[i])histoM02vsM20perSMvsE[i]->Write(Form("histoM02vsM20perSMvsE%d",i),TObject::kWriteDelete);
    }
    for(Int_t i=0;i<3;i++){
        if(histoIsolationCurveAllClusR005[i])histoIsolationCurveAllClusR005[i]->Write(Form("histoIsolationCurveAllClusR005_iso%d",i),TObject::kWriteDelete);
        if(histoIsolationCurveAllClusR01[i])histoIsolationCurveAllClusR01[i]->Write(Form("histoIsolationCurveAllClusR01_iso%d",i),TObject::kWriteDelete);
        if(histoIsolationCurveAllClusR02[i])histoIsolationCurveAllClusR02[i]->Write(Form("histoIsolationCurveAllClusR02_iso%d",i),TObject::kWriteDelete);
        if(histoIsolationCurveAllClusR04[i])histoIsolationCurveAllClusR04[i]->Write(Form("histoIsolationCurveAllClusR04_iso%d",i),TObject::kWriteDelete);
    }

    if (kMC){
        if(histoEFracFirstLabelvsE)histoEFracFirstLabelvsE->Write("histoEFracFirstLabelvsE",TObject::kWriteDelete);
        if(histoEFracLeadPi0vsE)histoEFracLeadPi0vsE->Write("histoEFracLeadPi0vsE",TObject::kWriteDelete);
        if(histoFracClusELeadPi0All)histoFracClusELeadPi0All->Write("histoFracClusELeadPi0All",TObject::kWriteDelete);
        if(histoFracClusELeadTrueMergedPi0)histoFracClusELeadTrueMergedPi0->Write("histoFracClusELeadTrueMergedPi0",TObject::kWriteDelete);
        if(histoFracClusELeadTrueMergedPartConvPi0)histoFracClusELeadTrueMergedPartConvPi0->Write("histoFracClusELeadTrueMergedPartConvPi0",TObject::kWriteDelete);
        if(histoFracClusELeadTrueGammaFromPi0)histoFracClusELeadTrueGammaFromPi0->Write("histoFracClusELeadTrueGammaFromPi0",TObject::kWriteDelete);

        for(Int_t i=0;i<3;i++){
            if(histoIsolationCurveAllPi0ClusR02[i])histoIsolationCurveAllPi0ClusR02[i]->Write(Form("histoIsolationCurveAllPi0ClusR02_iso%d",i),TObject::kWriteDelete);
            if(histoIsolationCurveTrueMergedPi0ClusR02[i])histoIsolationCurveTrueMergedPi0ClusR02[i]->Write(Form("histoIsolationCurveTrueMergedPi0ClusR02_iso%d",i),TObject::kWriteDelete);
            if(histoIsolationCurveTrueMergedPartConvPi0ClusR02[i])histoIsolationCurveTrueMergedPartConvPi0ClusR02[i]->Write(Form("histoIsolationCurveTrueMergedPartConvPi0ClusR02_iso%d",i),TObject::kWriteDelete);
            if(histoIsolationCurveGammaFromPi0ClusR02[i])histoIsolationCurveGammaFromPi0ClusR02[i]->Write(Form("histoIsolationCurveGammaFromPi0ClusR02_iso%d",i),TObject::kWriteDelete);
            if(histoIsolationCurveAllPi0ClusR04[i])histoIsolationCurveAllPi0ClusR04[i]->Write(Form("histoIsolationCurveAllPi0ClusR04_iso%d",i),TObject::kWriteDelete);
            if(histoIsolationCurveTrueMergedPi0ClusR04[i])histoIsolationCurveTrueMergedPi0ClusR04[i]->Write(Form("histoIsolationCurveTrueMergedPi0ClusR04_iso%d",i),TObject::kWriteDelete);
            if(histoIsolationCurveTrueMergedPartConvPi0ClusR04[i])histoIsolationCurveTrueMergedPartConvPi0ClusR04[i]->Write(Form("histoIsolationCurveTrueMergedPartConvPi0ClusR04_iso%d",i),TObject::kWriteDelete);
            if(histoIsolationCurveGammaFromPi0ClusR04[i])histoIsolationCurveGammaFromPi0ClusR04[i]->Write(Form("histoIsolationCurveGammaFromPi0ClusR04_iso%d",i),TObject::kWriteDelete);
        }
        for(Int_t i=0;i<4;i++){
            if(histoResolutionAllPi0Clus[i])histoResolutionAllPi0Clus[i]->Write(Form("histoResolutionAllPi0Clus_iso%d",i),TObject::kWriteDelete);
            if(histoResolutionTrueMergedPi0Clus[i])histoResolutionTrueMergedPi0Clus[i]->Write(Form("histoResolutionTrueMergedPi0Clus_iso%d",i),TObject::kWriteDelete);
            if(histoResolutionTrueMergedPartConvPi0Clus[i])histoResolutionTrueMergedPartConvPi0Clus[i]->Write(Form("histoResolutionTrueMergedPartConvPi0Clus_iso%d",i),TObject::kWriteDelete);
            if(histoResolutionTrueGammaFromPi0Clus[i])histoResolutionTrueGammaFromPi0Clus[i]->Write(Form("histoResolutionTrueGammaFromPi0Clus_iso%d",i),TObject::kWriteDelete);

            if(histoResoMatrixAllPi0Clus[i])histoResoMatrixAllPi0Clus[i]->Write(Form("histoResoMatrixAllPi0Clus_iso%d",i),TObject::kWriteDelete);
            if(histoResoMatrixTrueMergedPi0Clus[i])histoResoMatrixTrueMergedPi0Clus[i]->Write(Form("histoResoMatrixTrueMergedPi0Clus_iso%d",i),TObject::kWriteDelete);
            if(histoResoMatrixTrueMergedPartConvPi0Clus[i])histoResoMatrixTrueMergedPartConvPi0Clus[i]->Write(Form("histoResoMatrixTrueMergedPartConvPi0Clus_iso%d",i),TObject::kWriteDelete);
            if(histoResoMatrixTrueGammaFromPi0Clus[i])histoResoMatrixTrueGammaFromPi0Clus[i]->Write(Form("histoResoMatrixTrueGammaFromPi0Clus_iso%d",i),TObject::kWriteDelete);
        }
        for(Int_t i=0;i<20;i++){
            if(histoM02TrueMergedPi0perSM[i])histoM02TrueMergedPi0perSM[i]->Write(Form("histoM02TrueMergedPi0perSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueMergedEtaperSM[i])histoM02TrueMergedEtaperSM[i]->Write(Form("histoM02TrueMergedEtaperSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueMergedPartConvPi0perSM[i])histoM02TrueMergedPartConvPi0perSM[i]->Write(Form("histoM02TrueMergedPartConvPi0perSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueMergedPartConvEtaperSM[i])histoM02TrueMergedPartConvEtaperSM[i]->Write(Form("histoM02TrueMergedPartConvEtaperSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueGammaFromPi0perSM[i])histoM02TrueGammaFromPi0perSM[i]->Write(Form("histoM02TrueGammaFromPi0perSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueGammaFromEtaperSM[i])histoM02TrueGammaFromEtaperSM[i]->Write(Form("histoM02TrueGammaFromEtaperSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueElectronFromPi0perSM[i])histoM02TrueElectronFromPi0perSM[i]->Write(Form("histoM02TrueElectronFromPi0perSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueElectronFromEtaperSM[i])histoM02TrueElectronFromEtaperSM[i]->Write(Form("histoM02TrueElectronFromEtaperSM%d",i),TObject::kWriteDelete);
            if(histoM02TrueBackgroundperSM[i])histoM02TrueBackgroundperSM[i]->Write(Form("histoM02TrueBackgroundperSM%d",i),TObject::kWriteDelete);

            if(histoM02TrueMergedPi0perSMvsE[i])histoM02TrueMergedPi0perSMvsE[i]->Write(Form("histoM02TrueMergedPi0perSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueMergedEtaperSMvsE[i])histoM02TrueMergedEtaperSMvsE[i]->Write(Form("histoM02TrueMergedEtaperSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueMergedPartConvPi0perSMvsE[i])histoM02TrueMergedPartConvPi0perSMvsE[i]->Write(Form("histoM02TrueMergedPartConvPi0perSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueMergedPartConvEtaperSMvsE[i])histoM02TrueMergedPartConvEtaperSMvsE[i]->Write(Form("histoM02TrueMergedPartConvEtaperSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueGammaFromPi0perSMvsE[i])histoM02TrueGammaFromPi0perSMvsE[i]->Write(Form("histoM02TrueGammaFromPi0perSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueGammaFromEtaperSMvsE[i])histoM02TrueGammaFromEtaperSMvsE[i]->Write(Form("histoM02TrueGammaFromEtaperSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueElectronFromPi0perSMvsE[i])histoM02TrueElectronFromPi0perSMvsE[i]->Write(Form("histoM02TrueElectronFromPi0perSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueElectronFromEtaperSMvsE[i])histoM02TrueElectronFromEtaperSMvsE[i]->Write(Form("histoM02TrueElectronFromEtaperSMvsE%d",i),TObject::kWriteDelete);
            if(histoM02TrueBackgroundperSMvsE[i])histoM02TrueBackgroundperSMvsE[i]->Write(Form("histoM02TrueBackgroundperSMvsE%d",i),TObject::kWriteDelete);

            if(histoTrueMergedPi0M02vsM20perSMvsE[i])histoTrueMergedPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueMergedPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
            if(histoTrueMergedPartConvPi0M02vsM20perSMvsE[i])histoTrueMergedPartConvPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueMergedPartConvPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
            if(histoTrueGammaFromPi0M02vsM20perSMvsE[i])histoTrueGammaFromPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueGammaFromPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
            if(histoTrueElectronFromPi0M02vsM20perSMvsE[i])histoTrueElectronFromPi0M02vsM20perSMvsE[i]->Write(Form("histoTrueElectronFromPi0M02vsM20perSMvsE%d",i),TObject::kWriteDelete);
        }
    }
    filePhotoQAWrite->Write();
    filePhotoQAWrite->Close();
    if(filePhotoQAWrite) delete filePhotoQAWrite;

    if (kMC){
        if(histoEFracFirstLabelvsE) delete   histoEFracFirstLabelvsE;
        if(histoEFracLeadPi0vsE) delete   histoEFracLeadPi0vsE;
        if(histoFracClusELeadPi0All) delete   histoFracClusELeadPi0All;
        if(histoFracClusELeadTrueMergedPi0) delete   histoFracClusELeadTrueMergedPi0;
        if(histoFracClusELeadTrueMergedPartConvPi0) delete   histoFracClusELeadTrueMergedPartConvPi0;
        if(histoFracClusELeadTrueGammaFromPi0) delete   histoFracClusELeadTrueGammaFromPi0;
    }
    f->Close();
    delete f;

}
