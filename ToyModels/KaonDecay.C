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



//****************************************************************************
//*********************** Kaon decay simulation ******************************
//****************************************************************************
void KaonDecay(     Int_t nEvts         = 1000000, 
                    Int_t isMC          = 0,
                    TString energy      = "2.76TeV",
                    Double_t minPt      = 2,
                    Double_t maxPt      = 50,
                    TString filename    = "",
                    TString suffix      = "eps"
              ){

    //*************************************************************************************************
    //******************************** Set Style settings globally ************************************
    //*************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis(suffix);  
    SetPlotStyle();
    
    TString fCollisionSystenWrite       = ReturnCollisionEnergyOutputString(energy);
    TString fOutputDir                  = Form("%s/%s/",fCollisionSystenWrite.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+fOutputDir);
    
    //*************************************************************************************************
    //*************************** Initialize pt spectra ***********************************************
    //*************************************************************************************************
    Double_t massParticle       = 0;
    TF1* ptDistribution         = NULL; 
    TString outputlabel         = "";
    TString plotLabel           = "";
    TString daughterLabels[10]  = {"","","","","","","","","",""};
    Color_t colors[10]          = { kRed+1, kGreen+2, kOrange+7, kBlue+2, kAzure, 
                                    kYellow-8, kViolet, kCyan+2, kGray+1, kPink+2 };
    massParticle                = TDatabasePDG::Instance()->GetParticle(310)->Mass();
    Int_t nDaughters            = 0;
    Double_t* masses            = NULL;
    Int_t* pdgCodesDaughters    = NULL;
    outputlabel                 = "K";
    plotLabel                   = "K";
    daughterLabels[0]           = "#pi^{0}_{1}";
    daughterLabels[1]           = "#pi^{0}_{2}";
    nDaughters                  = 2;
    pdgCodesDaughters           = new Int_t[nDaughters];
    pdgCodesDaughters[0]        = 111;
    pdgCodesDaughters[1]        = 111;
    masses                      = new Double_t[nDaughters];
    masses[0]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();;
    masses[1]                   = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();;

    TFile* inputFile            = new TFile(filename.Data());
    TH1D* kaonSpectrum          = NULL; 
    // Data inputs
    TCanvas *canvasFitQA = new TCanvas("canvasFitQA","canvasFitQA",1000,800);
    DrawGammaCanvasSettings( canvasFitQA, 0.07, 0.02, 0.02, 0.08);
    canvasFitQA->cd();
    canvasFitQA->SetLogy(1);

    if (isMC == 0){
        if (energy.CompareTo("2.76TeV") == 0){
            kaonSpectrum            = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat2760GeV");
            Double_t paramGraph[3]  = {kaonSpectrum->GetBinContent(1), 6., 0.5};
            ptDistribution          = FitObject("l","ptDistribution","K",kaonSpectrum,0.1,20.,NULL,"QNRMEI");
            
            
            DrawAutoGammaMesonHistos(   kaonSpectrum, 
                                        "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", // (%)", 
                                        kTRUE, 10, 1e-10, kFALSE,
                                        kFALSE, 0., 0.7, 
                                        kFALSE, 0., 10.);
            kaonSpectrum->GetYaxis()->SetTitleOffset(0.85);
            DrawGammaSetMarker(kaonSpectrum, 20, 1.5, kAzure-6, kAzure-6);
            kaonSpectrum->DrawClone("pe");
            
            ptDistribution->SetLineColor(kRed+2);
            ptDistribution->Draw("same");
            
        } else if (energy.CompareTo("7TeV") == 0){
            kaonSpectrum        = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat7TeV");
            ptDistribution      = FitObject("tcm","ptDistribution","K",kaonSpectrum,0.4,20.);
            ptDistribution->SetRange(minPt,maxPt);
        } else {
            cout << "ERROR: undefined energy for data" << endl;
            return;
        }    
    // MC inputs    
    } else {
        cout << "decayChannel not defined, aborting" << endl;
    }    
    canvasFitQA->SaveAs(Form("%s%s_PtDistribution_Input_%d.%s", fOutputDir.Data(), outputlabel.Data(), isMC, suffix.Data()));    
        
    ptDistribution->SetRange(minPt,maxPt);
    
    //*************************************************************************************************
    //********************* Initialize random number generators and TGenPhaseSpace ********************
    //*************************************************************************************************
    TRandom3 randy1;
    TRandom3 randy2;
    TRandom3 randy3;
    TGenPhaseSpace event;
    

    TH2D *h2_ptMothervsDaughter = new TH2D("h2_ptMothervsDaughter","", 500,0,50,500,0,50);  
    h2_ptMothervsDaughter->Sumw2();
    TH2D *h2_asym_geom          = new TH2D("h2_asym_geom","", 500,-1,1, 500,0,50);
    h2_ptMothervsDaughter->Sumw2();
    
    TH1D *h1_ptdistribution     = new TH1D("h1_ptdistribution","", 500,0,50);  
    h1_ptdistribution->Sumw2();
    TH1D *h1_phidistribution    = new TH1D("h1_phidistribution","", 100,0,2*TMath::Pi());  
    h1_phidistribution->Sumw2();
    TH1D *h1_etadistribution    = new TH1D("h1_etadistribution","", 200,-1,1);  
    h1_etadistribution->Sumw2();
    TH1D *h1_ptDaughter[10];
    for (Int_t i = 0; i < nDaughters; i++){
        h1_ptDaughter[i]        = new TH1D(Form("h1_ptDaughter%d",i),"", 500,0,50);  
    }
    TH1D *h1_ptPiZeroDaughters  = new TH1D("h1_ptPiZeroDaughters","", 500,0,50);  

    //*************************************************************************************************
    //**************************** Event loop *********************************************************
    //*************************************************************************************************
    for(Int_t n=0; n<nEvts; n++){ // this is the important loop (nEvents)
//         cout << "event: " << n << endl;
        Double_t ptcurrent      = ptDistribution->GetRandom();
        
//         cout << "mother pt: " << ptcurrent << endl;

        Double_t phiCurrent     = randy1.Uniform(2*TMath::Pi());
        Double_t etaCurrent     = randy2.Uniform(-1,1);

        // assuming eta as a gaussian
//         Double_t etaCurrent     = randy2.Gaus(0,2);        
//         while (abs(etaCurrent) > 1){
//             etaCurrent          = randy2.Gaus(0,2);
//         }    
        h1_ptdistribution->Fill(ptcurrent);
        h1_phidistribution->Fill(phiCurrent);
        h1_etadistribution->Fill(etaCurrent);
        
        TLorentzVector particle(0.0, 0.0, 0, massParticle); 
        particle.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
        event.SetDecay(particle, nDaughters, masses); // ie. "set the pion to decay to 2 daughters with masses[0], masses[1]
        
//         cout << "mother pt vector: " << particle.Pt() << endl;
        
        
        Double_t weight         = event.Generate();
        Double_t ptDaughter[nDaughters];  
        for (Int_t i = 0; i < nDaughters; i++){
            TLorentzVector *pDaughter   = event.GetDecay(i); // these are my daughters !! 
            ptDaughter[i]               = pDaughter->Pt(); 
            h1_ptDaughter[i]->Fill(ptDaughter[i]);
            if (pdgCodesDaughters[i] == 111)
                h1_ptPiZeroDaughters->Fill(ptDaughter[i]);            
            if (i == 0)
                h2_ptMothervsDaughter->Fill(ptcurrent,ptDaughter[i]);
        }    
        if (nDaughters > 1){
            h2_asym_geom->Fill((ptDaughter[0]-ptDaughter[1])/(ptDaughter[0]+ptDaughter[1]), ptcurrent);
        }
    }

    TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1000,800);
    DrawGammaCanvasSettings( canvasQA, 0.07, 0.02, 0.02, 0.08);
    canvasQA->cd();
    canvasQA->SetLogy(1);
    TLegend* legendSpectra = GetAndSetLegend2(0.86, 0.70, 0.95, 0.93, 32,1); 
    
    DrawAutoGammaMesonHistos(   h1_ptdistribution, 
                                "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
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
    
    legendSpectra->Draw();
    canvasQA->SaveAs(Form("%s%s_PtDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    DrawAutoGammaMesonHistos(   h1_phidistribution, 
                                "", "#varphi (rad)", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_phidistribution->GetYaxis()->SetTitleOffset(0.85);    
    h1_phidistribution->DrawClone();
    canvasQA->SaveAs(Form("%s%s_PhiDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    
    DrawAutoGammaMesonHistos(   h1_etadistribution, 
                                "", "#eta", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_etadistribution->GetYaxis()->SetTitleOffset(0.85);        
    h1_etadistribution->DrawClone();
    canvasQA->SaveAs(Form("%s%s_EtaDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));

    TCanvas *canvasQA2D = new TCanvas("canvasQA2D","canvasQA2D",1000,800);
    DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
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

    h2_ptMothervsDaughter->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_PtMothervsDaughter1.%s", fOutputDir.Data() ,outputlabel.Data(),suffix.Data()));
    DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
    DrawAutoGammaHistoPaper2D(h2_asym_geom,
                                "",
                                "#it{A} = (#it{p}_{T,1}-#it{p}_{T,2})/(#it{p}_{T,1}+#it{p}_{T,2})",
                                Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);
    h2_asym_geom->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_AsymmetryDaughtersVsMotherPt.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    
    
    //*************************************************************************************************
    //********************** Write histograms to file *************************************************
    //*************************************************************************************************
    TFile* fileOutput = new TFile(Form("ToyMCOutput_%s_%d_%s.root",outputlabel.Data(), isMC, fCollisionSystenWrite.Data()),"RECREATE");  

        h1_ptdistribution->Write();
        h1_phidistribution->Write();
        h1_etadistribution->Write();
        for (Int_t i = 0; i < nDaughters; i++){
            h1_ptDaughter[i]->Write();
        }     
        h1_ptPiZeroDaughters->Write();
        h2_ptMothervsDaughter->Write();
    
    fileOutput->Write();
    fileOutput->Close();

    
}