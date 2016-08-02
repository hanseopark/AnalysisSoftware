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




void NeutralMesonDecayNonLin(   Int_t nEvts                 = 1000000, 
                                Int_t particle              = 1,
                                TString energy              = "2.76TeV",
                                Double_t minPt              = 2,
                                Double_t maxPt              = 50,
                                TString suffix              = "eps",
                                Double_t rEMC               = 440.,
                                Int_t mode                  = 10,
                                Int_t nonlin                = -1
                            ){

    //*************************************************************************************************
    //******************************** Set Style settings globally ************************************
    //*************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis(suffix);  
    SetPlotStyle();

    TString fCollisionSystenWrite       = ReturnCollisionEnergyOutputString(energy);
    TString fOutputDir                  = Form("%s/%s/NonLin/",fCollisionSystenWrite.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+fOutputDir);
    
    //*************************************************************************************************
    //*************************** Initialize pt spectra for particles**********************************
    //*************************************************************************************************
    Double_t massParticle       = 0;
    TF1* ptDistribution         = NULL; 
    TString outputlabel         = "";
    TString plotLabel           = "";
    TString daughterLabels[10]  = {"","","","","","","","","",""};
    Color_t colors[10]          = { kRed+1, kGreen+2, kOrange+7, kBlue+2, kAzure, 
                                    kYellow-8, kViolet, kCyan+2, kGray+1, kPink+2 };
    Int_t nDaughters            = 0;
    Double_t* masses            = NULL;
    Int_t* pdgCodesDaughters    = NULL;
    if (particle == 1){
        massParticle        = TDatabasePDG::Instance()->GetParticle(111)->Mass();
        if (energy.CompareTo("2.76TeV") == 0){
            ptDistribution      = FitObject("tcm","ptDistribution","Pi0",NULL,0.4,50.);
            ptDistribution->SetParameter(0,3.04120116616420654e+11);
            ptDistribution->SetParameter(1,1.33176621422733510e-01);
            ptDistribution->SetParameter(2,2.73953300141987877e+10);
            ptDistribution->SetParameter(3,5.64798456522145997e-01);
            ptDistribution->SetParameter(4,3.19909580984167752e+00);   
            ptDistribution->SetRange(minPt,maxPt);
        } else {
            cout << "ERROR: undefined energy for pi0" << endl;
            return;
        }    
        outputlabel         = "Pi0";
        plotLabel           = "#pi^{0}";
        daughterLabels[0]   = "#gamma_{1}";
        daughterLabels[1]   = "#gamma_{2}";
        nDaughters          = 2;
        pdgCodesDaughters   = new Int_t[nDaughters];
        pdgCodesDaughters[0]= 22;
        pdgCodesDaughters[1]= 22;
        masses              = new Double_t[nDaughters];
        masses[0]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();;
        masses[1]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();;
    } else if (particle == 2){
        massParticle        = TDatabasePDG::Instance()->GetParticle(221)->Mass();
        if (energy.CompareTo("2.76TeV") == 0){
            ptDistribution      = FitObject("tcm","ptDistribution","Eta",NULL,0.4,50.);
            ptDistribution->SetParameter(0,1.71597713217546654e+10);
            ptDistribution->SetParameter(1,1.51262659810865535e-01);
            ptDistribution->SetParameter(2,1.41711438002461696e+09);
            ptDistribution->SetParameter(3,8.49829150135323452e-01);
            ptDistribution->SetParameter(4,3.32752667783400913e+00);
            ptDistribution->SetRange(minPt,maxPt);
        } else {
            cout << "ERROR: undefined energy for eta" << endl;
            return;
        }    
        outputlabel         = "Eta";
        plotLabel           = "#eta";
        daughterLabels[0]   = "#gamma_{1}";
        daughterLabels[1]   = "#gamma_{2}";
        nDaughters          = 2;
        pdgCodesDaughters   = new Int_t[nDaughters];
        pdgCodesDaughters[0]= 22;
        pdgCodesDaughters[1]= 22;
        masses              = new Double_t[nDaughters];
        masses[0]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[0])->Mass();;
        masses[1]           = TDatabasePDG::Instance()->GetParticle(pdgCodesDaughters[1])->Mass();;
    } else {
        cout << "particle not defined, aborting" << endl;
    }    
        
    Color_t colorSmear[5]               = {kAzure, kRed+2, kGreen+2, kViolet+2, kCyan+2};
    Style_t markerSmear[5]              = {21, 20, 25, 24, 33};
    Double_t ptBinning[69]              = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                           2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                                           3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
                                           5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.5, 8.0, 8.5, 9.0,
                                           10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0,
                                           30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0};
    Int_t nBinsX                        = 68;    
    
    TF1* nonlinFunc         = NULL;
    Bool_t haveNonLin       = kFALSE;
    if (nonlin != -1){
        haveNonLin          = kTRUE;
        nonlinFunc          = new TF1("nonLin","1/((([0] +  [1] * TMath::Power(x,[2]))/([3] +  [4] * TMath::Power(x,[5])))+[6])");
        nonlinFunc->SetParameter(0,1.1100193881);
        nonlinFunc->SetParameter(1,-0.1389194936);
        nonlinFunc->SetParameter(2,-0.0800000242);
        nonlinFunc->SetParameter(3,1.1673716264);
        nonlinFunc->SetParameter(4,-0.1853095466);
        nonlinFunc->SetParameter(5,-0.0848801702);
        nonlinFunc->SetParameter(6,-0.017);
        nonlinFunc->SetRange(0,maxPt);
    }    
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
    TH2D *h2_OpenAngle          = new TH2D("h2_OpenAngle","", 500,0,50, 400,0,4);
    h2_OpenAngle->SetXTitle(Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()));
    h2_OpenAngle->SetYTitle("#theta_{open}");
    h2_OpenAngle->Sumw2();
    TH2D *h2_ClusterDistance    = new TH2D("h2_ClusterDistance","", 500,0,50, 5000,0,100);
    h2_ClusterDistance->SetXTitle(Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()));
    h2_ClusterDistance->SetYTitle("#it{R}_{#gamma's at r=440cm} (cm)");
    h2_ClusterDistance->Sumw2();
    TH2D *h2_E_EAfterNonlin     = new TH2D("h2_E_EAfterNonlin","", 500,0,50, 500,0,50);
    h2_E_EAfterNonlin->Sumw2();
    TH2D *h2_Pt_PtAfterNonlin   = new TH2D("h2_Pt_PtAfterNonlin","", 500,0,50, 500,0,50);
    h2_Pt_PtAfterNonlin->Sumw2();
    
    TH1D *h1_ptdistribution     = new TH1D("h1_ptdistribution","", 500,0,50);  
    h1_ptdistribution->Sumw2();
    TH1D *h1_ptdistributionReb  = new TH1D("h1_ptdistributionReb","", nBinsX, ptBinning);  
    h1_ptdistributionReb->Sumw2();
    
    TH1D *h1_phidistribution    = new TH1D("h1_phidistribution","", 100,0,2*TMath::Pi());  
    h1_phidistribution->Sumw2();
    TH1D *h1_etadistribution    = new TH1D("h1_etadistribution","", 200,-1,1);  
    h1_etadistribution->Sumw2();
    TH1D *h1_ptDaughter[10];
    for (Int_t i = 0; i < nDaughters; i++){
        h1_ptDaughter[i]        = new TH1D(Form("h1_ptDaughter%d",i),"", 500,0,50);  
    }
    TH1D *h1_ptGammaDaughters   = new TH1D("h1_ptGammaDaughters","", 500,0,50);  
    TH1D *h1_ptPiPlDaughters    = new TH1D("h1_ptPiPlDaughters","", 500,0,50);  
    TH1D *h1_ptPiMiDaughters    = new TH1D("h1_ptPiMiDaughters","", 500,0,50);  
    TH1D *h1_ptPiZeroDaughters  = new TH1D("h1_ptPiZeroDaughters","", 500,0,50);  
    
    //*************************************************************************************************
    //**************************** Event loop *********************************************************
    //*************************************************************************************************
    for(Int_t n=0; n<nEvts; n++){ // this is the important loop (nEvents)
        if (n%10000000 == 0) 
            cout << "generated " << (Double_t)n/1e6 << " Mio events" << endl;
//         cout << "event: " << n << endl;
        Double_t ptcurrent      = ptDistribution->GetRandom();        
//         cout << "mother pt: " << ptcurrent << endl;

        Double_t phiCurrent     = randy1.Uniform(2*TMath::Pi());
        Double_t etaCurrent     = randy2.Uniform(-1,1);

        h1_ptdistribution->Fill(ptcurrent);
        h1_ptdistributionReb->Fill(ptcurrent);
        h1_phidistribution->Fill(phiCurrent);
        h1_etadistribution->Fill(etaCurrent);
    
        // assuming eta as a gaussian
//         Double_t etaCurrent     = randy2.Gaus(0,2);        
//         while (abs(etaCurrent) > 1){
//             etaCurrent          = randy2.Gaus(0,2);
//         }    
        
        TLorentzVector particle(0.0, 0.0, 0, massParticle); 
        particle.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
        event.SetDecay(particle, nDaughters, masses); // ie. "set the pion to decay to 2 daughters with masses[0], masses[1]
        
//         cout << "mother pt vector: " << particle.Pt() << endl;
        
        if (haveNonLin){
            TLorentzVector particle2(0.0, 0.0, 0, massParticle); 
            particle2.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
//             cout << "before: " << particle2.Energy() << "\t after: " << particle2.Energy()*nonlinFunc->Eval(particle2.Energy())<< endl;
            particle2.SetE(particle2.Energy()*nonlinFunc->Eval(particle2.Energy()));
            
            
            
            h2_E_EAfterNonlin->Fill(particle.Energy(),particle2.Energy());
            h2_Pt_PtAfterNonlin->Fill(particle.Pt(),particle2.Pt());
//         event.SetDecay(particle, nDaughters, masses); // ie. "set the pion to decay to 2 daughters with masses[0], masses[1]
        }

        
        Double_t weight         = event.Generate();
        Double_t ptDaughter[nDaughters];  
        TLorentzVector* pDaughter[nDaughters];  
        for (Int_t i = 0; i < nDaughters; i++){
            pDaughter[i]                = event.GetDecay(i); // these are my daughters !! 
            ptDaughter[i]               = pDaughter[i]->Pt(); 
            
            h1_ptDaughter[i]->Fill(ptDaughter[i]);
            if (pdgCodesDaughters[i] == 22)
                h1_ptGammaDaughters->Fill(ptDaughter[i]);
            if (pdgCodesDaughters[i] == 111)
                h1_ptPiZeroDaughters->Fill(ptDaughter[i]);
            if (pdgCodesDaughters[i] == 211)
                h1_ptPiPlDaughters->Fill(ptDaughter[i]);
            if (pdgCodesDaughters[i] == -211)
                h1_ptPiMiDaughters->Fill(ptDaughter[i]);
            
            if (i == 0)
                h2_ptMothervsDaughter->Fill(ptcurrent,ptDaughter[i]);
        }    
        if (nDaughters > 1){
            h2_asym_geom->Fill((ptDaughter[0]-ptDaughter[1])/(ptDaughter[0]+ptDaughter[1]), ptcurrent);
            Double_t openangle      = pDaughter[0]->Angle(pDaughter[1]->Vect());
            Double_t distanceCl     = TMath::Tan(openangle/2)*rEMC*2;
            
            h2_OpenAngle->Fill(ptcurrent,openangle);
            h2_ClusterDistance->Fill(ptcurrent,distanceCl);
        }
    }

    
    //****************************** Plot input pT and decay distributions ***********************************************************************
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
    DrawGammaSetMarker(h1_ptdistribution, 20, 1.5, kBlack, kBlack);
    h1_ptdistribution->DrawClone("pe");
    legendSpectra->AddEntry(h1_ptdistribution,plotLabel.Data(),"pe");
    
    for (Int_t i = 0; i < nDaughters; i++){
        DrawGammaSetMarker(h1_ptDaughter[i], 21+i, 1., colors[i], colors[i]);
        h1_ptDaughter[i]->Draw("same,pe");
        legendSpectra->AddEntry(h1_ptDaughter[i],daughterLabels[i].Data(),"pe");
    }
    if (h1_ptGammaDaughters->GetEntries() > 0){
        DrawGammaSetMarker(h1_ptGammaDaughters, 24, 1.5, kGray+2, kGray+2);
        h1_ptGammaDaughters->Draw("same,pe");
        legendSpectra->AddEntry(h1_ptGammaDaughters,"#gamma's","pe");
    }
    
    legendSpectra->Draw();
    canvasQA->SaveAs(Form("%s%s_PtDistribution_Input.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    
    //******************************** Plot input phi phi distribution ***************************************************************
    canvasQA->SetLogy(1);
    DrawAutoGammaMesonHistos(   h1_phidistribution, 
                                "", "#varphi (rad)", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_phidistribution->GetYaxis()->SetTitleOffset(0.85);    
    h1_phidistribution->DrawClone();
    canvasQA->SaveAs(Form("%s_PhiDistribution_Input.%s",outputlabel.Data(),suffix.Data()));
    
    DrawAutoGammaMesonHistos(   h1_etadistribution, 
                                "", "#eta", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_etadistribution->GetYaxis()->SetTitleOffset(0.85);        
    h1_etadistribution->DrawClone();
    canvasQA->SaveAs(Form("%s%s_EtaDistribution_Input.%s",fOutputDir.Data(), outputlabel.Data(),suffix.Data()));

    if (haveNonLin){
        canvasQA->cd();
        canvasQA->SetLogy(0);
        TH2F* dummyDrawingHist  = new TH2F("dummyDrawingHist","dummyDrawingHist",5000,0,maxPt,10000, 0.95, 1.1); 
        SetStyleHistoTH2ForGraphs(  dummyDrawingHist, "#it{p}_{T} (GeV/#it{c})", "corr fac.", 0.028, 0.04, 
                                    0.028, 0.04, 0.86, 0.82, 510, 505);
        dummyDrawingHist->Draw();
        
        DrawGammaSetMarkerTF1( nonlinFunc, 7, 2, kBlue); 
        nonlinFunc->Draw("same");
 
        canvasQA->SaveAs(Form("%sNonlinearityCorr.%s",fOutputDir.Data(), suffix.Data()));
        
    }
    
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
    canvasQA2D->SaveAs(Form("%s%s_PtMothervsDaughter1.%s",fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
    DrawAutoGammaHistoPaper2D(h2_asym_geom,
                                "",
                                "#it{A} = (#it{p}_{T,1}-#it{p}_{T,2})/(#it{p}_{T,1}+#it{p}_{T,2})",
                                Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);
    h2_asym_geom->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_AsymmetryDaughtersVsMotherPt.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));
    
    DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
    DrawAutoGammaHistoPaper2D(h2_OpenAngle,
                                "",
                                Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                                "#theta_{open}",
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);
    h2_OpenAngle->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_MotherPtVsOpeningAngleDaughters.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));

    DrawAutoGammaHistoPaper2D(h2_ClusterDistance,
                                "",
                                Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                                "#it{R}_{#gamma's at r=440cm} (cm)",
                                0,0,0,
                                0,0,5,
                                0,0, 50,1,0.85);
    h2_ClusterDistance->Draw("colz");
    canvasQA2D->SaveAs(Form("%s%s_MotherPtVsClusterDistance.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));

    if (haveNonLin){
        DrawAutoGammaHistoPaper2D(h2_E_EAfterNonlin,
                                    "",
                                    Form("#it{E}_{%s} (GeV)",plotLabel.Data()),
                                    Form("#it{E}_{%s after NL} (GeV)",plotLabel.Data()),
                                    0,0,0,
                                    0,0,5,
                                    0,0, 50,1,0.85);
        h2_E_EAfterNonlin->Draw("colz");
        canvasQA2D->SaveAs(Form("%s%s_MotherEnergyVsEnergyAfterNonlin.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));

        DrawAutoGammaHistoPaper2D(h2_Pt_PtAfterNonlin,
                                    "",
                                    Form("#it{p}_{T,%s} (GeV/#it{c})",plotLabel.Data()),
                                    Form("#it{p}_{T,%s after NL} (GeV/#it{c})",plotLabel.Data()),
                                    0,0,0,
                                    0,0,5,
                                    0,0, 50,1,0.85);
        h2_Pt_PtAfterNonlin->Draw("colz");
        canvasQA2D->SaveAs(Form("%s%s_MotherPtVsPtAfterNonlin.%s",fOutputDir.Data(),outputlabel.Data(),suffix.Data()));
    }
    //*************************************************************************************************
    //********************** Write histograms to file *************************************************
    //*************************************************************************************************
    TFile* fileOutput = new TFile(Form("%s/ToyMCOutputNonLin_%s.root",fCollisionSystenWrite.Data(), outputlabel.Data()),"RECREATE");  

        h1_ptdistribution->Write();
        h1_phidistribution->Write();
        h1_etadistribution->Write();
        for (Int_t i = 0; i < nDaughters; i++){
            h1_ptDaughter[i]->Write();
        }     
        h2_ptMothervsDaughter->Write();
        h2_ClusterDistance->Write();
        h2_OpenAngle->Write();
        h2_asym_geom->Write();
        
    fileOutput->Write();
    fileOutput->Close();

    
}