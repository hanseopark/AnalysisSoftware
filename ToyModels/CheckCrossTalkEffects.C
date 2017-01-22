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




void CheckCrossTalkEffects( Int_t nEvts                 = 1000000, 
                            Int_t particle              = 1,
                            TString energy              = "2.76TeV",
                            Double_t minPt              = 2,
                            Double_t maxPt              = 50,
                            TString suffix              = "eps",
                            TString fileName            = "",
                            TString cutSting            = "",
                            Int_t mode                  = 10,
                            TString subMode             = "",
                            TString fileNameXTalk       = ""
                          ){

    //*************************************************************************************************
    //******************************** Set Style settings globally ************************************
    //*************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    
    StyleSettingsThesis(suffix);  
    SetPlotStyle();
    TString fCollisionSystem            = ReturnFullCollisionsSystem(energy);
    TString fCollisionSystenWrite       = ReturnCollisionEnergyOutputString(energy);
    TString fOutputDir                  = Form("%s/%s/CrossTalkCheck/",fCollisionSystenWrite.Data(),suffix.Data());
    gSystem->Exec("mkdir -p "+fOutputDir);
    
    //*************************************************************************************************
    //*************************** Initialize pt spectra for particles**********************************
    //*************************************************************************************************
    Double_t massParticle       = 0;
    TF1* ptDistribution         = NULL; 
    TString outputlabel         = "";
    TString plotLabel           = "";
    Color_t colors[10]          = { kRed+1, kGreen+2, kOrange+7, kBlue+2, kAzure, 
                                    kYellow-8, kViolet, kCyan+2, kGray+1, kPink+2 };
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
    } else {
        cout << "particle not defined, aborting" << endl;
    }    
        
    Color_t colorSmear[5]               = {kAzure, kRed+2, kGreen+2, kViolet+2, kCyan+2};
    Style_t markerSmear[5]              = {21, 20, 25, 24, 33};

    Double_t ptBinning[32]              = {3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,
                                           5.0, 5.4, 5.8, 6.2, 6.6, 7.0, 7.5, 8.0, 8.5, 9.0,
                                           10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 25.0,
                                           30.0, 50.0};
    Int_t nBinsX                        = 31;    
    Double_t energyFracLeadCell[3]      = {0.8, 0.6, 1};
    //*************************************************************************************************    
    //************************** Loading the resolutions **********************************************
    //*************************************************************************************************
    Bool_t haveResol                = kFALSE;
    Bool_t haveFractions            = kFALSE;
    TH2F* histoResolutionInputAll   = NULL;
    TString labelResolHistAll;
    TString labelResolHist[4];
    Int_t nResolHist                = 4;
    if (fileName.CompareTo("")!= 0 && cutSting.CompareTo("") != 0){
        TFile* resolutionFile           = new TFile(fileName.Data());
        TString autoDetectedMainDir     = AutoDetectMainTList(mode , resolutionFile);
        if (autoDetectedMainDir.CompareTo("") == 0){
            cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
            return;
        }
        
        //************************** Container Loading ********************************************************************
        TList* TopDir                   = (TList*)resolutionFile->Get(autoDetectedMainDir.Data());
        if(TopDir == NULL){
            cout<<"ERROR: TopDir not Found"<<endl;
            return;
        }
        TList *HistosGammaConversion    = (TList*)TopDir->FindObject(Form("Cut Number %s",cutSting.Data()));
        TList *TrueConversionContainer  = (TList*)HistosGammaConversion->FindObject(Form("%s True histograms",cutSting.Data()));

        //****************************** Resolution dPt vs Pt **********************************************
        TString nameResolHistAll           = Form("ESD_TruePrimary%s_MCPt_ResolPt",outputlabel.Data());
        labelResolHistAll                  = Form("%s", plotLabel.Data());
        
        histoResolutionInputAll           = (TH2F*)TrueConversionContainer->FindObject(nameResolHistAll.Data()); //
        if (histoResolutionInputAll){
            histoResolutionInputAll->Sumw2();
            haveResol               = kTRUE;            
            cout << "have found resolution hist" << endl;
        } else {
            haveResol               = kFALSE;            
        }    
    
        TString nameResolHist[5];       
        nameResolHist[0]                = Form("ESD_TruePrimary%s_MCPt_ResolPt",outputlabel.Data());
//         labelResolHist[0]               = Form("%s", plotLabel.Data());
        if (mode == 10){
            nameResolHist[0]            = Form("ESD_TruePrimary%sPureMerged_MCPt_ResolPt",outputlabel.Data());
            nameResolHist[1]            = Form("ESD_TruePrimary%sMergedPartConv_MCPt_ResolPt",outputlabel.Data());
            nameResolHist[2]            = Form("ESD_TruePrimary%s1Gamma_MCPt_ResolPt",outputlabel.Data());
            nameResolHist[3]            = Form("ESD_TruePrimary%s1Electron_MCPt_ResolPt",outputlabel.Data());
        }    
        TH2F* histoResolutionInput[5];
        for (Int_t i = 0; i < 4; i++){
            histoResolutionInput[i]     = (TH2F*)TrueConversionContainer->FindObject(nameResolHist[i].Data()); //
            if (histoResolutionInput[i]){
                histoResolutionInput[i]->Sumw2();
                if (i == 0){
                    histoResolutionInputAll = (TH2F*)histoResolutionInput[i]->Clone("resolAll");
                    histoResolutionInputAll->Sumw2();
                } else {    
                    histoResolutionInputAll->Add(histoResolutionInput[i]);
                }    
                haveResol               = kTRUE;            
            } else {
                haveResol               = kFALSE;            
            }    
        }
    }
    TH2F* histoXTalkSmear               = NULL;
    TH2F* histoXTalkSmearRescaled       = new TH2F("histoXTalkSmearRescaled","", nBinsX, ptBinning, 100, -1, 1);  
    Bool_t haveXTalk                    = kFALSE;
    TH2F* histoLeftRight[2]             = {NULL, NULL};
    TH2F* histoLeftRightRescaled[2]     = {NULL, NULL};
    histoLeftRightRescaled[0]           = new TH2F("histoLeftRightRescaled_0","", nBinsX, ptBinning, 100, -1, 1);  
    histoLeftRightRescaled[1]           = new TH2F("histoLeftRightRescaled_1","", nBinsX, ptBinning, 100, -1, 1);  
    if (fileNameXTalk.CompareTo("") != 0){    
        TFile* xTalkFile                = new TFile(fileNameXTalk.Data());
        histoLeftRight[0]               = (TH2F*)xTalkFile->Get("hLeftMinusRightVsEcell_0");
        histoLeftRight[1]               = (TH2F*)xTalkFile->Get("hLeftMinusRightVsEcell_1");
        histoXTalkSmear                 = (TH2F*)histoLeftRight[0]->Clone("histoXTalkSmear");
        histoXTalkSmear->Sumw2();
        histoXTalkSmear->Add(histoLeftRight[1]);
        haveXTalk                       = kTRUE;
        for (Int_t i = 1; i< histoXTalkSmear->GetNbinsX()+1; i++){ // x-axis
            for (Int_t j = 0; j < histoXTalkSmear->GetNbinsY()+1; j++) { // y-axis
                Double_t currDiff       = histoXTalkSmear->GetYaxis()->GetBinCenter(j)/histoXTalkSmear->GetXaxis()->GetBinCenter(i);
                Int_t amount            = histoXTalkSmear->GetBinContent(i,j);
                Int_t binXNew           = histoXTalkSmearRescaled->GetXaxis()->FindBin(histoXTalkSmear->GetXaxis()->GetBinCenter(i));
                Int_t binYNew           = histoXTalkSmearRescaled->GetYaxis()->FindBin(currDiff);
                if (amount > 0){
//                     cout << i << "\t" <<  j << "\t"<< histoXTalkSmear->GetXaxis()->Getroot BinCenter(i) << "\t" << histoXTalkSmear->GetYaxis()->GetBinCenter(j) << "\t" << currDiff << "\t" << amount << endl;
                    Double_t amountNew  = amount + histoXTalkSmearRescaled->GetBinContent(binXNew,binYNew);
                    histoXTalkSmearRescaled->SetBinContent(binXNew,binYNew,amountNew);
                }    
                amount                  = histoLeftRight[0]->GetBinContent(i,j);
                binXNew                 = histoLeftRightRescaled[0]->GetXaxis()->FindBin(histoLeftRight[0]->GetXaxis()->GetBinCenter(i));
                binYNew                 = histoLeftRightRescaled[0]->GetYaxis()->FindBin(currDiff);
                if (amount > 0){
                    Double_t amountNew  = amount + histoLeftRightRescaled[0]->GetBinContent(binXNew,binYNew);
                    histoLeftRightRescaled[0]->SetBinContent(binXNew,binYNew,amountNew);
                }    
                amount                  = histoLeftRight[1]->GetBinContent(i,j);
                binXNew                 = histoLeftRightRescaled[1]->GetXaxis()->FindBin(histoLeftRight[1]->GetXaxis()->GetBinCenter(i));
                binYNew                 = histoLeftRightRescaled[1]->GetYaxis()->FindBin(currDiff);
                if (amount > 0){
                    Double_t amountNew  = amount + histoLeftRightRescaled[1]->GetBinContent(binXNew,binYNew);
                    histoLeftRightRescaled[1]->SetBinContent(binXNew,binYNew,amountNew);
                }    
            }
        }
    }
//     return;
    //*************************************************************************************************
    //********************* Initialize random number generators and TGenPhaseSpace ********************
    //*************************************************************************************************
    TRandom3 randy1;
    TRandom3 randy2;
    TRandom3 randy3;
    TGenPhaseSpace event;
    
    TH1D *h1_ptdistribution     = new TH1D("h1_ptdistribution","", 500,0,50);  
    h1_ptdistribution->Sumw2();
    TH1D *h1_ptdistributionReb  = new TH1D("h1_ptdistributionReb","", nBinsX, ptBinning);  
    h1_ptdistributionReb->Sumw2();
    TH2D *h2_deltaRel           = new TH2D("h2_deltaRel","", nBinsX, ptBinning, 500, 0, 1);  
    h2_deltaRel->Sumw2();

    TH1D* h1_ptdistSmeared[nResolHist];
    TH1D* h1_ptdistSmearedReb[nResolHist];
    for (Int_t i = 0; i < nResolHist; i++){
        h1_ptdistSmeared[i]         = new TH1D(Form("h1_ptdistSmeared_%d",i),"", 500,0,50);  
        h1_ptdistSmeared[i]->Sumw2();
        h1_ptdistSmearedReb[i]      = new TH1D(Form("h1_ptdistSmearedReb_%d",i),"", nBinsX, ptBinning);  
        h1_ptdistSmearedReb[i]->Sumw2();
    }
    
    TH1D *h1_phidistribution    = new TH1D("h1_phidistribution","", 100,0,2*TMath::Pi());  
    h1_phidistribution->Sumw2();
    TH1D *h1_etadistribution    = new TH1D("h1_etadistribution","", 200,-1,1);  
    h1_etadistribution->Sumw2();

    TH1D* h1ResolProjection[1000];
    if (haveResol){
        for (Int_t i=0; (i < histoResolutionInputAll->GetNbinsX()+1 && i < 1000); i++){
            h1ResolProjection[i] = (TH1D*)histoResolutionInputAll->ProjectionY(Form("dummy_%d",i), i+1,i+1,"e");
        }
    }   
    TH1D* h1XTalkProjection[nBinsX];
    TH1D* h1XTalkLeftProjection[nBinsX];
    TH1D* h1XTalkRightProjection[nBinsX];
    if (haveXTalk){
        for (Int_t i=0; i < nBinsX; i++){
            h1XTalkProjection[i]        = (TH1D*)histoXTalkSmearRescaled->ProjectionY(Form("dummyXTalk_%d",i), i+1, i+1,"e");
            h1XTalkProjection[i]->Sumw2();
            h1XTalkProjection[i]->Scale(1./h1XTalkProjection[i]->GetEntries());
            h1XTalkLeftProjection[i]    = (TH1D*)histoLeftRightRescaled[0]->ProjectionY(Form("dummyXTalkLeft_%d",i), i+1, i+1,"e");
            h1XTalkRightProjection[i]   = (TH1D*)histoLeftRightRescaled[1]->ProjectionY(Form("dummyXTalkRight_%d",i), i+1, i+1,"e");
        }    
    }
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
    
        TLorentzVector particle(0.0, 0.0, 0, massParticle); 
        particle.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
        
        Double_t shift          = 0;
        
        Double_t ptcurrentDist      = ptcurrent;
        Double_t ptcurrentDistMod1  = ptcurrent;
        if (haveResol){
            shift                   = h1ResolProjection[histoResolutionInputAll->GetXaxis()->FindBin(ptcurrent)-1]->GetRandom();
            //             mesoncand->Pt()-mother->Pt())/mother->Pt()
            ptcurrentDist           = (shift*ptcurrent)+ptcurrent; 
//                 cout << ptcurrent << "\t" << shift <<"\t"<< ptcurrentDist<< endl; 
        }
        if (haveXTalk){
            for (Int_t i = 0; i < 3; i++){
                Double_t ReldeltaE      = 0;
                Int_t currBin           = 0;
                if (ptcurrent > 3. && ptcurrent < 50.){
    //                 cout << histoXTalkSmear->GetXaxis()->FindBin(particle.E()) << "\t" << particle.E() << endl;
                    currBin             = h2_deltaRel->GetXaxis()->FindBin(ptcurrent*energyFracLeadCell[i])-1;
                    ReldeltaE           = TMath::Abs(h1XTalkProjection[h2_deltaRel->GetXaxis()->FindBin(ptcurrent*energyFracLeadCell[i])-1]->GetRandom());
                }    
                h2_deltaRel->Fill(0.8*ptcurrent,ReldeltaE);
                Double_t amountLeft     =  h1XTalkLeftProjection[currBin]->GetBinContent(h1XTalkLeftProjection[currBin]->GetXaxis()->FindBin(ReldeltaE));
                Double_t amountRight    =  h1XTalkRightProjection[currBin]->GetBinContent(h1XTalkRightProjection[currBin]->GetXaxis()->FindBin(ReldeltaE));
    //             if (ReldeltaE > 0.7){
    //                 cout << ReldeltaE << endl;
    //             }    
                Double_t deltaEbef      = ReldeltaE;
                if (!((amountLeft+amountRight) > 0))
                    ReldeltaE           = 0.8*ReldeltaE;
                else 
                    ReldeltaE               = TMath::Abs(amountLeft-amountRight)/((amountLeft+amountRight)*2) * ReldeltaE;
                
    //             if (deltaEbef > 0.7){
    //                 cout << amountLeft << "\t" << amountRight << "\t" <<  "\t" << ReldeltaE<< "\t"<< deltaEbef<< endl;
    //             }    
    //             deltaE                  = 0.5;
                ptcurrentDistMod1       = (shift*(ptcurrent-ReldeltaE*ptcurrent)+(ptcurrent-ReldeltaE*ptcurrent));
                h1_ptdistSmeared[i+1]->Fill(ptcurrentDistMod1);
                h1_ptdistSmearedReb[i+1]->Fill(ptcurrentDistMod1);            
            }    
        }    

        h1_ptdistSmeared[0]->Fill(ptcurrentDist);
        h1_ptdistSmearedReb[0]->Fill(ptcurrentDist);            
        
        // assuming eta as a gaussian
//         Double_t etaCurrent     = randy2.Gaus(0,2);        
//         while (abs(etaCurrent) > 1){
//             etaCurrent          = randy2.Gaus(0,2);
//         }    
        
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
        
    legendSpectra->Draw();
    canvasQA->SaveAs(Form("%s%s_PtDistribution_Input.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    //****************************** Plot smeared and input pT ***********************************************************************
    canvasQA->cd();
    TLegend* legendSpectra2 = GetAndSetLegend2(0.45, 0.93-(nResolHist+1)*0.035, 0.65, 0.93, 32,1); 
    DrawGammaSetMarker(h1_ptdistribution, 20, 1.5, kBlack, kBlack);
    h1_ptdistribution->DrawClone("pe");
    legendSpectra2->AddEntry(h1_ptdistribution,plotLabel.Data(),"pe");
    for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
        DrawGammaSetMarker(h1_ptdistSmeared[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
        h1_ptdistSmeared[j]->Draw("same,pe");
        legendSpectra2->AddEntry(h1_ptdistSmeared[j],Form("%s smeared",labelResolHistAll.Data()),"pe");
    }
    legendSpectra2->Draw();
    
    canvasQA->SaveAs(Form("%s%s_PtDistribution_InputVsSmeared%s.%s",fOutputDir.Data(), outputlabel.Data(), subMode.Data(), suffix.Data()));

    
    //**************************** Calculate and draw ratio of smeared and input dist ************************************************
    
    TH1D* histoRatioSmearedDivInput[nResolHist];
    
    canvasQA->cd();
    canvasQA->SetLogy(0);

    TLegend* legendRatio1 = GetAndSetLegend2(0.15, 0.93-nResolHist*0.035, 0.35, 0.93, 32,1); 
    for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
        histoRatioSmearedDivInput[j]    = (TH1D*)h1_ptdistSmeared[j]->Clone(Form("histoRatioSmearedDivInput%d",j));
        histoRatioSmearedDivInput[j]->Divide(histoRatioSmearedDivInput[j],h1_ptdistribution,1,1,"B");
        if (j == 0){
            DrawAutoGammaMesonHistos(   histoRatioSmearedDivInput[j], 
                                        "", "#it{p}_{T} (GeV/#it{c})", "smeared/input", // (%)", 
                                        kFALSE, 10, 1e-1, kFALSE,
                                        kTRUE, 0., 2, 
                                        kFALSE, 0., 10.);
            histoRatioSmearedDivInput[j]->GetYaxis()->SetTitleOffset(0.85);
            DrawGammaSetMarker(histoRatioSmearedDivInput[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
            histoRatioSmearedDivInput[j]->DrawClone("pe");
            
        }  else {
            DrawGammaSetMarker(histoRatioSmearedDivInput[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
            histoRatioSmearedDivInput[j]->DrawClone("same,pe");            
        }   
        legendRatio1->AddEntry(histoRatioSmearedDivInput[j],labelResolHistAll.Data(),"pe");
    }
    legendRatio1->Draw();
    
    DrawGammaLines(0, maxPt , 1, 1 ,1, kGray, 7);   
    DrawGammaLines(0, maxPt , 1.2, 1.2 ,1, kGray, 8);   
    DrawGammaLines(0, maxPt , 0.8, 0.8 ,1, kGray, 8);   
    
    canvasQA->SaveAs(Form("%s%s_Ratio_SmearedDivInputVsPt.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    //**************************** Calculate and draw ratio of smeared and input dist ************************************************
    TH1D* histoRatioSmearedDivInputReb[nResolHist];
    canvasQA->cd();
    canvasQA->SetLogy(0);

    for (Int_t j = 0; (j < 5 && j < nResolHist); j++ ){
        histoRatioSmearedDivInputReb[j]    = (TH1D*)h1_ptdistSmearedReb[j]->Clone(Form("histoRatioSmearedDivInputReb%d",j));
        histoRatioSmearedDivInputReb[j]->Divide(histoRatioSmearedDivInputReb[j],h1_ptdistributionReb,1,1,"B");
        if (j == 0){
            DrawAutoGammaMesonHistos(   histoRatioSmearedDivInputReb[j], 
                                        "", "#it{p}_{T} (GeV/#it{c})", "smeared/input", // (%)", 
                                        kFALSE, 10, 1e-1, kFALSE,
                                        kTRUE, 0., 2, 
                                        kTRUE, 0., maxPt);
            histoRatioSmearedDivInputReb[j]->GetYaxis()->SetTitleOffset(0.85);
            DrawGammaSetMarker(histoRatioSmearedDivInputReb[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
            histoRatioSmearedDivInputReb[j]->DrawClone("pe");
        }  else {
            DrawGammaSetMarker(histoRatioSmearedDivInputReb[j], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
            histoRatioSmearedDivInputReb[j]->DrawClone("same,pe");            
        }    
    }
    legendRatio1->Draw();
    
    DrawGammaLines(0, maxPt , 1, 1 ,1, kGray, 7);   
    DrawGammaLines(0, maxPt , 1.2, 1.2 ,1, kGray, 8);   
    DrawGammaLines(0, maxPt , 0.8, 0.8 ,1, kGray, 8);   
    
    canvasQA->SaveAs(Form("%s%s_Ratio_SmearedDivInputVsPtRebined.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    
        TH1D* histoRatioSmearedDivSmeared[4];
        for (Int_t j = 1; j < nResolHist; j++){
            histoRatioSmearedDivSmeared[j-1]    = (TH1D*)h1_ptdistSmearedReb[j]->Clone(Form("histoRatioSmearedDivSmearedReb%d",j));
            histoRatioSmearedDivSmeared[j-1]->Divide(histoRatioSmearedDivSmeared[j-1],h1_ptdistSmearedReb[0],1,1,"B");
            if (j == 1){
                DrawAutoGammaMesonHistos(   histoRatioSmearedDivSmeared[j-1], 
                                            "", "#it{p}_{T} (GeV/#it{c})", "smeared mod/smeared", // (%)", 
                                            kFALSE, 10, 1e-1, kFALSE,
                                            kTRUE, 0.7, 1.3, 
                                            kTRUE, 10., maxPt);
                histoRatioSmearedDivSmeared[j-1]->GetYaxis()->SetTitleOffset(0.85);
                DrawGammaSetMarker(histoRatioSmearedDivSmeared[j-1], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
                histoRatioSmearedDivSmeared[j-1]->DrawClone("pe");
            }  else {
                DrawGammaSetMarker(histoRatioSmearedDivSmeared[j-1], markerSmear[j], 1.5, colorSmear[j], colorSmear[j]);
                histoRatioSmearedDivSmeared[j-1]->DrawClone("same,pe");
            }
        }   
        DrawGammaLines(10, maxPt , 1, 1 ,1, kGray, 7);   
        DrawGammaLines(10, maxPt , 1.2, 1.2 ,1, kGray, 8);   
        DrawGammaLines(10, maxPt , 0.8, 0.8 ,1, kGray, 8);   
    
    canvasQA->SaveAs(Form("%s%s_Ratio_SmearedDivSmearedVsPtRebined.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    
    
    
    //******************************** Plot input phi phi distribution ***************************************************************
    canvasQA->SetLogy(1);
    DrawAutoGammaMesonHistos(   h1_phidistribution, 
                                "", "#varphi (rad)", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_phidistribution->GetYaxis()->SetTitleOffset(0.85);    
    h1_phidistribution->DrawClone();
    canvasQA->SaveAs(Form("%s%s_PhiDistribution_Input.%s",fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
    
    DrawAutoGammaMesonHistos(   h1_etadistribution, 
                                "", "#eta", Form("#it{N}_{%s}",plotLabel.Data()), // (%)", 
                                kTRUE, 10, 1e-1, kFALSE,
                                kFALSE, 0., 0.7, 
                                kFALSE, 0., 10.);
    h1_etadistribution->GetYaxis()->SetTitleOffset(0.85);        
    h1_etadistribution->DrawClone();
    canvasQA->SaveAs(Form("%s%s_EtaDistribution_Input.%s",fOutputDir.Data(), outputlabel.Data(),suffix.Data()));

    
    if (haveResol){
        TCanvas *canvasQA2D = new TCanvas("canvasQA2D","canvasQA2D",1000,800);
        DrawGammaCanvasSettings( canvasQA2D, 0.1, 0.09, 0.02, 0.11);
        canvasQA2D->SetLogz(1);
        canvasQA2D->SetLeftMargin(0.115);
        canvasQA2D->SetRightMargin(0.11);

        DrawAutoGammaHistoPaper2D(  histoResolutionInputAll,
                                    "",
                                    "#it{p}_{T,MC} (GeV/#it{c})", 
                                    Form ("(#it{p}^{%s}_{T,rec} -#it{p}^{%s}_{T,MC})/#it{p}^{%s}_{T,MC}", plotLabel.Data(), plotLabel.Data(), plotLabel.Data()),
                                    0,0,0,
                                    0,0,5,
                                    0,0, 50,1,0.85);
        histoResolutionInputAll->GetZaxis()->SetRangeUser(1e-7, histoResolutionInputAll->GetMaximum());
        histoResolutionInputAll->GetYaxis()->SetTitleOffset(1.05);
        histoResolutionInputAll->Draw("colz");
        PutProcessLabelAndEnergyOnPlot(0.15, 0.25, 28, fCollisionSystem.Data(), "mEMC", labelResolHistAll.Data(), 63, 0.03);

        canvasQA2D->SaveAs(Form("%s%s_Resolutions.jpg",fOutputDir.Data(),outputlabel.Data()));
    }
    
    if (haveXTalk ){
        canvasQA->cd();
//         canvasQA->SetLogy(0);
        for (Int_t i=0; i < nBinsX; i++){
            if (i == 0){
                DrawAutoGammaMesonHistos(   h1XTalkProjection[i], 
                                            "", "#it{#delta E}/ #it{E}", "", // (%)", 
                                            kTRUE, 10, 1e-4, kFALSE,
                                            kFALSE, 0., 2, 
                                            kTRUE, -1., 1);
                h1XTalkProjection[i]->GetYaxis()->SetTitleOffset(0.85);
                DrawGammaSetMarker(h1XTalkProjection[i], markerSmear[i%10], 1.5, colorSmear[i%10], colorSmear[i%10]);
                h1XTalkProjection[i]->DrawClone("pe");
            }  else {
                DrawGammaSetMarker(h1XTalkProjection[i], markerSmear[i%10], 1.5, colorSmear[i%10], colorSmear[i%10]);
                h1XTalkProjection[i]->DrawClone("same,pe");            
            }    
        }
        
        canvasQA->SaveAs(Form("%s%s_Projections.%s",fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    }    
    //*************************************************************************************************
    //********************** Write histograms to file *************************************************
    //*************************************************************************************************
    TFile* fileOutput = new TFile(Form("%s/ToyMCOutput_CrossTalkCheck_%s.root",fCollisionSystenWrite.Data(), outputlabel.Data()),"RECREATE");  

        h1_ptdistribution->Write();
//         if(h1_ptdistSmeared) h1_ptdistSmeared->Write();
//         if(histoRatioSmearedDivInputReb) histoRatioSmearedDivInputReb->Write();
//         if(histoRatioSmearedDivInput) histoRatioSmearedDivInput->Write();
        h1_phidistribution->Write();
        h1_etadistribution->Write();
        histoXTalkSmearRescaled->Write();
        h2_deltaRel->Write();
        for (Int_t i=0; i < nBinsX; i++){
            h1XTalkProjection[i]->Write();
        }
    fileOutput->Write();
    fileOutput->Close();

    
}