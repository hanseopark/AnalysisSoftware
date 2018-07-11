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

// ****************************************************************************************************************
// ****************************************************************************************************************
// ****************************************************************************************************************
TString AutoDetectTreeListTrigger(TList* fList, Int_t mode = 0){
    TString listName = "";
    TString nominalListName = "TriggerQA tree";
    TList *readList = (TList*)fList->At(1);
    listName = readList->GetName();
    if (listName.Contains(nominalListName) ) {
        cout << Form("-> found : %s",listName.Data()) << endl;
        return listName;
    } else {
        cout << "Could not find list named *TreeList* as last object of main list" << endl;
        return "";
    }
}


void BuildHistogramsTriggerQA(  TString fileName                = "GammaConvV1_QA_5460001022092970003190000000.root",
                                TString mainDirName             = "GammaTriggerQA",
                                TString fCutSelection           = "5460001022092970003190000000",
                                TString specificCutSelection    = "000000",
                                Bool_t merge                    = 0,
                                TString nameOutputBase          = "TriggerQA"
        ){

    //Reset ROOT and connect tree file
    gROOT->Reset();
    //********************************************************************************
    //*            Definition of Cuts                                                *
    //********************************************************************************
    TString CentralitySel               = specificCutSelection(0,1);
    TString centMinCutNumber            = specificCutSelection(1,2);
    TString centMaxCutNumber            = specificCutSelection(3,3);

    Int_t centSel                       = CentralitySel.Atoi();
    switch (centSel){
        case 0:
            cout << "CentralitySel: " <<  "VZERO cent sel" << endl;
        case 1:
            cout << "CentralitySel: " <<  "VZERO thresholdSel" << endl;
        case 2:
            cout << "CentralitySel: " <<  "TZERO semiCent" << endl;
        case 3:
            cout << "CentralitySel: " <<  "TZERO cent" << endl;
        default:
            cout << "CentralitySel not defined" << endl;
    }
    UInt_t centMinCut                 = centMinCutNumber.Atoi();
    UInt_t centMinBin                 = ((UInt_t)centMinCutNumber.Atoi())/10;
    UInt_t centMaxCut                 = centMaxCutNumber.Atoi();
    UInt_t centMaxBin                 = ((UInt_t)centMaxCutNumber.Atoi())/10;
                                        // 0,      10,     20,     30,   40,   50,   60,  70,  80,  90, 100
    UInt_t V0TriggTuned[11]           = {100000, 20800, 14770, 9865, 5000, 3545, 2000, 800, 100, 10, 0};

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
    TList* fTriggerMainDir          = (TList*)f->Get(mainDirName.Data());
    cout << mainDirName.Data() << endl;

    TString outputFolderName        = Form("Trigger_%s_%s",fCutSelection.Data(), specificCutSelection.Data());
    cout << "Will be written to " << outputFolderName.Data() << endl;

    TList* mainList                 = NULL;
    if (fTriggerMainDir){
        mainList                    = (TList*)fTriggerMainDir->FindObject(Form("Cut Number %s",fCutSelection.Data()));
    } else {
        cout << "main List not found" << endl;
    }
    if (!mainList){
        cout << "Main List not found, aboriting" << endl;
        return;
    }
    TString autoDetectedTreeList    = AutoDetectTreeListTrigger(mainList, 1);
    TList* treeList                 = (TList*)mainList->FindObject(autoDetectedTreeList.Data());
    if (!treeList) cout << "Tree List not found" << endl;


    //Declaration of leaves types
    Float_t fCent           = 0.;
    UShort_t fT0Trigg       = 0;
    UInt_t fV0Mult          = 0;
    UInt_t fV0Trigg         = 0;
    UInt_t fTPCMult         = 0;
    UInt_t fSPDHit          = 0;
    UInt_t fSPDTracklet     = 0;
    Float_t fZVertex        = 0;

    TTree *TriggerQA                 = (TTree*)treeList->FindObject("TriggerInfoTree");
    TriggerQA->SetBranchAddress("Cent",&fCent);
    TriggerQA->SetBranchAddress("T0Trigg",&fT0Trigg);
    TriggerQA->SetBranchAddress("V0Mult",&fV0Mult);
    TriggerQA->SetBranchAddress("V0Trigg",&fV0Trigg);
    TriggerQA->SetBranchAddress("TPCMult",&fTPCMult);
    TriggerQA->SetBranchAddress("SPDTrack",&fSPDHit);
    TriggerQA->SetBranchAddress("SPDTracklet",&fSPDTracklet);
    TriggerQA->SetBranchAddress("ZVertex",&fZVertex);

    //********************************************************************************
    //*            Definition of Boundaries for Histograms                           *
    //********************************************************************************
    Int_t nBinsCent              = 100;
    Double_t firstBinCent        = 0;
    Double_t lastBinCent         = 100;
    Int_t nBinsV0Mult            = 50000;
    Double_t firstBinV0Mult      = 0;
    Double_t lastBinV0Mult       = 50000;
    Int_t nBinsTPCMult           = 8000;
    Double_t firstBinTPCMult     = 0;
    Double_t lastBinTPCMult      = 8000;
    Int_t nBinsT0Trigg           = 100;
    Double_t firstBinT0Trigg     = 0;
    Double_t lastBinT0Trigg      = 100;

    //********************************************************************************
    //*      Definition of histograms for reconstructed Conversion Points            *
    //********************************************************************************
    TH1F* histoCent              = NULL;
    TH1F* histoV0Mult            = NULL;
    TH1F* histoV0Trigg           = NULL;
    TH1F* histoT0Trigg           = NULL;
    TH1F* histoTPCMult           = NULL;
    TString fileNameOutput       = Form("%s_Data.root",nameOutputBase.Data()) ;
    TFile* fileTriggerDetailed      = new TFile(fileNameOutput.Data());
    TDirectory*  directoryTrigger   = (TDirectory*)fileTriggerDetailed->Get(outputFolderName.Data());

    if (!merge || directoryTrigger==0) {
        histoCent                   = new TH1F("histoCent","", nBinsCent, firstBinCent, lastBinCent);
        histoV0Mult                 = new TH1F("histoV0Mult","", nBinsV0Mult, firstBinV0Mult, lastBinV0Mult);
        histoV0Trigg                = new TH1F("histoV0Trigg","", nBinsV0Mult, firstBinV0Mult, lastBinV0Mult);
        histoT0Trigg                = new TH1F("histoT0Trigg","", nBinsT0Trigg, firstBinT0Trigg, lastBinT0Trigg);
        histoTPCMult                = new TH1F("histoTPCMult","", nBinsTPCMult, firstBinTPCMult, lastBinTPCMult);
    } else {
        histoCent                   = (TH1F*)directoryTrigger->Get("histoCent");
        histoV0Mult                 = (TH1F*)directoryTrigger->Get("histoV0Mult");
        histoV0Trigg                = (TH1F*)directoryTrigger->Get("histoV0Trigg");
        histoTPCMult                = (TH1F*)directoryTrigger->Get("histoTPCMult");
        histoT0Trigg                = (TH1F*)directoryTrigger->Get("histoT0Trigg");
        // store sum of square of weights for the histograms
        histoCent->Sumw2();
        histoV0Mult->Sumw2();
        histoV0Trigg->Sumw2();
        histoTPCMult->Sumw2();
        histoT0Trigg->Sumw2();
    }


    //********************************************************************************
    //*      Reading of Tree with reconstructed gammas/ filling of histograms        *
    //********************************************************************************

    Long64_t nEntriesEvent                 = TriggerQA->GetEntries();
    Int_t nGammas                           = 0;
    cout << "Number of Events to be processed: " << nEntriesEvent << endl;
    //    (*ptGammaAssosiatedPi0)[k]
    Long64_t nbytesRecGam                   = 0;
    if (centSel==1) cout << centMinBin << "\t"<<  V0TriggTuned[centMinBin] << "\t" <<centMaxBin << "\t"<< V0TriggTuned[centMaxBin] << endl;
    for (Long64_t i=0; i<nEntriesEvent;i++) {
        nbytesRecGam                        += TriggerQA->GetEntry(i);
        if (centSel==0 && fCent > centMinCut && fCent < centMaxCut){
            histoCent->Fill(fCent);
            histoV0Mult->Fill(fV0Mult);
            histoTPCMult->Fill(fTPCMult);
            histoV0Trigg->Fill(fV0Trigg);
            histoT0Trigg->Fill(fT0Trigg);
        } else if (centSel==1 && fV0Trigg < V0TriggTuned[centMinBin] && fV0Trigg > V0TriggTuned[centMaxBin]){
            histoCent->Fill(fCent);
            histoV0Mult->Fill(fV0Mult);
            histoV0Trigg->Fill(fV0Trigg);
            histoTPCMult->Fill(fTPCMult);
            histoT0Trigg->Fill(fT0Trigg);
        } else if (centSel==2 && fT0Trigg > pow(2,3)){
            histoCent->Fill(fCent);
            histoV0Mult->Fill(fV0Mult);
            histoV0Trigg->Fill(fV0Trigg);
            histoTPCMult->Fill(fTPCMult);
            histoT0Trigg->Fill(fT0Trigg);
        } else if (centSel==3 && fT0Trigg > pow(2,4)){
            histoCent->Fill(fCent);
            histoV0Mult->Fill(fV0Mult);
            histoV0Trigg->Fill(fV0Trigg);
            histoTPCMult->Fill(fTPCMult);
            histoT0Trigg->Fill(fT0Trigg);
        }
    }

    //********************************************************************************
    //*                  Writing histograms to outputfile                            *
    //********************************************************************************
    TFile* fileTriggerQAWrite = new TFile(fileNameOutput.Data(),"UPDATE");

        fileTriggerQAWrite->mkdir(outputFolderName.Data());
        fileTriggerQAWrite->cd(outputFolderName.Data());

        histoCent->Write("histoCent",TObject::kWriteDelete);
        histoV0Mult->Write("histoV0Mult",TObject::kWriteDelete);
        histoV0Trigg->Write("histoV0Trigg",TObject::kWriteDelete);
        histoTPCMult->Write("histoTPCMult",TObject::kWriteDelete);
        histoT0Trigg->Write("histoT0Trigg",TObject::kWriteDelete);

    fileTriggerQAWrite->Write();
    fileTriggerQAWrite->Close();
    delete fileTriggerQAWrite;
    delete   histoCent;
    delete   histoV0Mult;

    if(directoryTrigger){
        directoryTrigger->Clear();
        delete directoryTrigger;
    }

    if(fileTriggerDetailed){
        fileTriggerDetailed->Close();
        delete fileTriggerDetailed;
    }

    if(treeList){
        treeList->Clear();
        delete treeList;
    }

    f->Close();
    delete f;

}
