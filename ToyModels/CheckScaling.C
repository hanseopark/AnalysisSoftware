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
//**************** Test toy for pion simulation ******************************
//****************************************************************************
void CheckScaling(      Int_t nEvts             = 1000000, 
                        TString energy          = "2.76TeV",
                        Double_t minPt          = 2,
                        Double_t maxPt          = 50,
                        TString suffix          = "eps",
                        TString cutSelection    = "",
                        Int_t mode              = 0,
                        TString externalFile    = ""
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
    TString rapidityRange               = "";
    Double_t deltaRapid                 =  ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);

    TString centralityString            = GetCentralityString(fEventCutSelection);
    cout << "centrality : "<< centralityString.Data() << endl;
    
    //*************************************************************************************************
    //*************************** Initialize pt spectra ***********************************************
    //*************************************************************************************************
    Double_t massParticle       = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    
    TString outputlabel         = "Pi0";
    TString plotLabel           = "#pi^{0}";
    
    Double_t yRanges[10]        = { 0.1,    0.3,    0.5,    0.55,   0.6,
                                    0.65,   0.7,    0.75,   0.8,    0.85 };
                                    
    Double_t branchRatio        = 0.98798; 
                                    
    // define output file name
    TString nameOutputRootFile          = Form("ToyMCOutput_%s_%s_%dMioEvt_%dMeV_%dMeV.root",outputlabel.Data(), fCollisionSystenWrite.Data(),(Int_t)(nEvts/1e6), (Int_t)(minPt*1000), (Int_t)(maxPt*1000));
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
    gSystem->Exec("mkdir -p "+fOutputDir);


    //*************************************************************************************************
    //************************* Global variables for particle spectra from input **********************
    //*************************************************************************************************
    TFile* inputFile            = new TFile(Form("%s/%s/Pi0_data_GammaConvV1Correction_%s.root",cutSelection.Data(), energy.Data(), cutSelection.Data()));
    TFile* inputFile2           = new TFile(Form("%s/%s/Pi0_MC_GammaConvV1CorrectionHistos_%s.root",cutSelection.Data(), energy.Data(), cutSelection.Data()));
    TFile* inputFile3           = new TFile(externalFile.Data());
    TH1D* partSpectrum          = NULL; 
    TH1D* partSpectrumInRap     = NULL;
    TH1D* histNEvtMC            = NULL;
    TH1D* chHadDNdEta           = NULL;
    TF1* ptDistribution         = NULL; 
    TSpline3* partParam         = NULL; 
    TH1D* ratioPartSpecToFit    = NULL; 
    TH1D* ratioPartSpecToSpline = NULL; 
    
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

    if (energy.CompareTo("2.76TeV") == 0){
        partSpectrum            = (TH1D*)inputFile->Get("MCYield_Meson_oldBin");
        partSpectrumInRap       = (TH1D*)inputFile2->Get("MC_Meson_genPt_oldBin");
        histNEvtMC              = (TH1D*)inputFile2->Get("NEvents");
        
        Float_t nEvtMC = 0;
        if ( fCollisionSystenWrite.Contains("pp") ){
            nEvtMC = GetNEvents(histNEvtMC);
        } else {
            nEvtMC = histNEvtMC->GetBinContent(1);
        }
        partSpectrumInRap->Scale(1./nEvtMC);
        
        chHadDNdEta             = (TH1D*)inputFile3->Get("histoChargedHadrondNdEtaALICEPP2760GeV");
        if (chHadDNdEta){
            cout << "found ch. had dist" << endl;
            chHadDNdEta->Sumw2();
            chHadDNdEta->Scale(1./chHadDNdEta->Integral());
        }    
        Double_t paramGraph[3]  = {partSpectrum->GetBinContent(1), 6., 0.5};
        ptDistribution          = FitObject("l","ptDistribution","Pi0",partSpectrum,5.,20.,NULL,"QNRME");
        
        Double_t xValues[150];
        Double_t xErr[150];
        Double_t yValues[150];
        Double_t yErr[150];
        Int_t nBinsX            = 0;
        Int_t nBinsXOrg         = 0;
        for (Int_t i = 1; (i< partSpectrum->GetNbinsX()+1 && partSpectrum->GetBinCenter(i) < 12.) ; i++){
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
        for (Int_t i = nBinsXOrg; (i < 150 && xValues[i-1] < 70); i++){
            xValues[i]          = 12+(i-nBinsXOrg);    
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
        }
        
        TH2F* dummyDrawingHist  = new TH2F("dummyDrawingHist","dummyDrawingHist",5000,0,70,10000, 1e-14, 1); 
        SetStyleHistoTH2ForGraphs(  dummyDrawingHist, "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", 0.028, 0.04, 
                                    0.028, 0.04, 0.9, 0.95, 510, 505);
        dummyDrawingHist->Draw();
        
        DrawGammaSetMarker(partSpectrum, 20, 1.5, kAzure-6, kAzure-6);
        partSpectrum->Draw("samepe");

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
        
    } else if (energy.CompareTo("7TeV") == 0){
        partSpectrum        = (TH1D*)inputFile->Get("histoChargedKaonSpecPubStat7TeV");
        ptDistribution      = FitObject("tcm","ptDistribution","K",partSpectrum,0.4,20.);
        ptDistribution->SetRange(minPt,maxPt);
    } else {
        cout << "ERROR: undefined energy for kaons" << endl;
        return;
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

            DrawGammaLines(0., 70,1., 1.,0.1);

        canvasFitQA->SaveAs(Form("%s%s_RatioToFitDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));
    }
    
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
    Double_t scaleFactorDecayPart   = 2*TMath::Pi()* (deltaEta);
    cout << "scale factor due to particle decay: " << scaleFactorDecayPart << endl;
    
    // find the scale factor due to the restriction in the pt range
    Double_t scaleFactor    = partSpectrumFromParam->Integral(partSpectrumFromParam->FindBin(minPt),partSpectrumFromParam->FindBin(maxPt))/partSpectrumFromParam->Integral();
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
    
    //*************************************************************************************************
    //********************** Create histograms for output *********************************************
    //*************************************************************************************************
    TH1D* h1_NEvt                       = new TH1D("h1_NEvt","", 1, 0, 1);
    h1_NEvt->SetBinContent(1,nEvts);
    TH1D *h1_ptdistribution             = new TH1D("h1_ptdistribution","", 700,0,70);  
    h1_ptdistribution->Sumw2();
    TH1D* h1_ptdistributionPerEv        = (TH1D*)h1_ptdistribution->Clone("h1_ptdistributionPerEv");
    h1_ptdistributionPerEv->Sumw2();
    TH1D *h1_phidistribution            = new TH1D("h1_phidistribution","", 100,0,2*TMath::Pi());  
    h1_phidistribution->Sumw2();
    TH1D *h1_etadistribution            = new TH1D("h1_etadistribution","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);  
    h1_etadistribution->Sumw2();
    TH1D *h1_ydistribution              = new TH1D("h1_ydistribution","", 2*fYMaxMother*100,-fYMaxMother,fYMaxMother);  
    h1_ydistribution->Sumw2();
    
    TH1D *h1_ptPiZeroInRap[10];
    for (Int_t y = 0; y < 10; y++){
        h1_ptPiZeroInRap[y]    = new TH1D(Form("h1_ptPiZeroInRap_%1.2f",yRanges[y]),"", 700,0,70);  
        h1_ptPiZeroInRap[y]->Sumw2();
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
        Double_t weightFull     = 1./nEvts*scaleFactor*ptcurrent; // if input is dN/dydpt 1/2pi 1/pt needs to be multiplied with pt
        Double_t weightPartial  = 1./nEvts*scaleFactor;
        
        // create current particle
        TLorentzVector particle(0.0, 0.0, 0, massParticle); 
        particle.SetPtEtaPhiM(ptcurrent, etaCurrent, phiCurrent, massParticle);
     
        // filling of input distributions
        h1_ptdistribution->Fill(ptcurrent, weightFull);
        h1_ptdistributionPerEv->Fill(ptcurrent, weightPartial);
        h1_phidistribution->Fill(phiCurrent, weightFull);
        h1_etadistribution->Fill(etaCurrent, weightFull);
        h1_ydistribution->Fill(particle.Rapidity(), weightFull);
        
        for (Int_t y = 0; y < 10; y++){
            if (TMath::Abs(particle.Rapidity()) < yRanges[y]){
                h1_ptPiZeroInRap[y]->Fill(ptcurrent, weightFull);
            }
        }    
    }

    //*************************************************************************************************
    //******************** Multiply with correct branching Ratio for decay ****************************
    //*************************************************************************************************
        
    for (Int_t y = 0; y < 10; y++){
        h1_ptPiZeroInRap[y]->Scale(branchRatio);
    }    
    
    //*************************************************************************************************
    //******************** Determine whether additional offset was accidentally created ***************
    //*************************************************************************************************    
    TF1* ptDistributionRefit            = NULL;    
    if (energy.CompareTo("2.76TeV") == 0){
        ptDistributionRefit         = FitObject("l","ptDistribution","Pi0",h1_ptdistributionPerEv,5.,20.,NULL,"QNRME");
    } else { 
        cout << "ranges for second fit not defined" << endl;  
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
        for (Int_t y = 0; y < 10; y++){
            h1_ptPiZeroInRap[y]->Scale(additionalScalingFac);
        }    
    }
    
    //*************************************************************************************************
    //******************************** Plot pt distribution   *****************************************
    //*************************************************************************************************
    
    TCanvas *canvasQA = new TCanvas("canvasQA","canvasQA",1000,800);
    DrawGammaCanvasSettings( canvasQA, 0.07, 0.02, 0.02, 0.08);
    canvasQA->cd();
    canvasQA->SetLogy(1);
    TLegend* legendSpectra = GetAndSetLegend2(0.76, 0.70, 0.95, 0.93, 32,1); 
    
    for (Int_t y = 0; y< 10; y++){
        if ( TMath::Abs(yRanges[y]-deltaRapid/2) < 0.00001){
            cout << "rapidity:" <<  y << "\t" << yRanges[y] << "\t" << deltaRapid/2 << endl;
            TH1D* histoInRapScaled      = (TH1D*)h1_ptPiZeroInRap[y]->Clone("histoInRapScaled");
            histoInRapScaled->Scale(scaleFactorDecayPart);
            DrawAutoGammaMesonHistos(   histoInRapScaled, 
                                        "", "#it{p}_{T} (GeV/#it{c})", "#it{N}_{X}", // (%)", 
                                        kTRUE, 10, 1e-13, kFALSE,
                                        kFALSE, 0., 0.7, 
                                        kFALSE, 0., 10.);
            histoInRapScaled->GetYaxis()->SetTitleOffset(0.85);
            DrawGammaSetMarker(histoInRapScaled, 20, 1.5, kAzure-6, kAzure-6);
            histoInRapScaled->DrawClone("pe");
            legendSpectra->AddEntry(histoInRapScaled,plotLabel.Data(),"pe");
        }
    }    
    DrawGammaSetMarker(partSpectrumInRap, 24, 1.5, kRed+2, kRed+2);
    partSpectrumInRap->DrawClone("same,pe");


    
    legendSpectra->AddEntry(partSpectrumInRap,Form("%s input", plotLabel.Data()),"pe");
    legendSpectra->Draw();
    canvasQA->SaveAs(Form("%s%s_PtDistribution_InputGeneratedInRap.%s", fOutputDir.Data(), outputlabel.Data(), suffix.Data()));

    
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

    canvasQA->SaveAs(Form("%s%s_YDistribution_Input.%s", fOutputDir.Data(), outputlabel.Data(),suffix.Data()));
        
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
        for (Int_t y = 0; y < 10; y++){
            h1_ptPiZeroInRap[y]->Write();
        }
    fileOutput->Write();
    fileOutput->Close();

    cout << "adding output file name to ToyMCOutputs.txt" << endl;
    gSystem->Exec("echo '"+nameOutputRootFile+"' >> ToyMCOutputs.txt");
   
}