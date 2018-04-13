/****************************************************************************************************************************
 ****** 	provided by Gamma Conversion Group, PWGGA, 																	*****
 ******		Friederike Bock, friederike.bock@cern.ch																	*****
 *****************************************************************************************************************************/

#include <Riostream.h>
#include <fstream>
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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TProfile.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

TString collisionSystem;
TString date;
TString fDetectionProcess ="";
TString ptHardRange = "";

extern TRandom*	gRandom;
extern TBenchmark*	gBenchmark;
extern TSystem*	gSystem;
extern TMinuit*  	gMinuit;

struct SysErrorConversion {
    Double_t value;
    Double_t error;
    //	TString name;
};

//**********************************************************************************
//******************* return minimum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn2D(TH2F* histo){
    if(!histo){cout << "WARNING: FindSmallestEntryIn2D, NULL pointer - returning 1... " << endl; return 1;}
    Double_t minimum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        for (Int_t j = 1; j<histo->GetNbinsY(); j++){
            if (histo->GetBinContent(i,j) < minimum && histo->GetBinContent(i,j) > 0){
                minimum = histo->GetBinContent(i,j);
            }
        }
    }
    return minimum;
}

//**********************************************************************************
//******************* return maximum for 2 D histo  ********************************
//**********************************************************************************
Double_t FindLargestEntryIn2D(TH2F* histo){
    if(!histo){cout << "WARNING: FindLargestEntry1D, NULL pointer - returning 1... " << endl; return 1;}
    Double_t maximum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        for (Int_t j = 1; j<histo->GetNbinsY(); j++){
            if (histo->GetBinContent(i,j) > maximum && histo->GetBinContent(i,j) > 0){
                maximum = histo->GetBinContent(i,j);
            }
        }
    }
    return maximum;
}


//**********************************************************************************
//******************* return minimum for 1 D histo  ********************************
//**********************************************************************************
Double_t FindSmallestEntryIn1D(TH1F* histo){
    if(!histo){cout << "WARNING: FindSmallestEntryIn1D, NULL pointer - returning 1... " << endl; return 1;}
    Double_t minimum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        if (histo->GetBinContent(i) < minimum && histo->GetBinContent(i) > 0){
            minimum = histo->GetBinContent(i);
        }
    }
    return minimum;
}

//**********************************************************************************
//******************* return maximum for 1 D histo  ********************************
//**********************************************************************************
Double_t FindLargestEntryIn1D(TH1F* histo){
    if(!histo){cout << "WARNING: FindLargestEntry1D, NULL pointer - returning 1... " << endl; return 1;}
    Double_t maximum = 1;
    for (Int_t i = 1; i<histo->GetNbinsX(); i++){
        if (histo->GetBinContent(i) > maximum ){
            maximum = histo->GetBinContent(i);
        }
    }
    return maximum;
}

//**********************************************************************************
//******************* Standardized plotting of 2D plots ****************************
//**********************************************************************************
void PlotStandard2D( TH2* histo2D, TString nameOutput, TString title, TString xTitle, TString yTitle, Bool_t kRangeY, Double_t startY, Double_t endY, Bool_t kRangeX, Double_t startX, Double_t endX, Int_t logX, Int_t logY, Int_t logZ, Float_t* floatLogo, Int_t canvasSizeX = 500, Int_t canvasSizeY = 500, TString generator ="" , TString period =""){
    TCanvas * canvasStandard = new TCanvas("canvasStandard","",10,10,canvasSizeX,canvasSizeY);  // gives the page size
    canvasStandard->SetLogx(logX);
    canvasStandard->SetLogy(logY);
    canvasStandard->SetLogz(logZ);
    canvasStandard->SetRightMargin(0.12);
    canvasStandard->SetLeftMargin(0.12);
    canvasStandard->SetBottomMargin(0.1);
    canvasStandard->SetTopMargin(0.04);
    canvasStandard->cd();
    histo2D->SetTitle("");
    DrawAutoGammaHisto2D(   histo2D,
                            title.Data(), xTitle.Data(), yTitle.Data(),"",kRangeY, startY, endY, kRangeX, startX, endX);
    histo2D->GetXaxis()->SetTitleOffset(1.05);
    //    cout << histo2D->GetYaxis()->GetTitleOffset() << endl;
    histo2D->GetYaxis()->SetTitleOffset(1.35);
    if (logX==1){
    //       cout << histo2D->GetXaxis()->GetLabelOffset() << endl;
        histo2D->GetXaxis()->SetLabelOffset(0.);
    }

    histo2D->Draw("colz");
    DrawLabelsEvents(floatLogo[0],floatLogo[1],floatLogo[2], 0.00, collisionSystem, generator, period);
    TLatex *detprocess = 	new TLatex(floatLogo[0], floatLogo[1] - 3.2*floatLogo[2], fDetectionProcess);
    detprocess->SetNDC();
    detprocess->SetTextColor(1);
    detprocess->SetTextSize(floatLogo[2]);
    detprocess->Draw();
    TLatex *ptHardBin = 	new TLatex(floatLogo[0], floatLogo[1] - 4.25*floatLogo[2], ptHardRange);
    ptHardBin->SetNDC();
    ptHardBin->SetTextColor(1);
    ptHardBin->SetTextSize(floatLogo[2]);
    ptHardBin->Draw();


    canvasStandard->Update();
    canvasStandard->SaveAs(nameOutput.Data());
    delete canvasStandard;
}

//**********************************************************************************
//************************ Scale the MC yield **************************************
//**********************************************************************************
void ScaleMCYield( TH1D* histoCorrectedToBeScaled, Double_t deltaRapid, Double_t scaling, Double_t nEvtMC, TString nameMeson){
    histoCorrectedToBeScaled->Sumw2();
    histoCorrectedToBeScaled->Scale(1./deltaRapid);
    histoCorrectedToBeScaled->Scale(scaling);
    histoCorrectedToBeScaled->Scale(1./nEvtMC);
    for (Int_t i = 1; i < histoCorrectedToBeScaled->GetNbinsX()+1 ; i++){
        Double_t newBinContent  = histoCorrectedToBeScaled->GetBinContent(i)/histoCorrectedToBeScaled->GetBinCenter(i)/histoCorrectedToBeScaled->GetBinWidth(i);
        Double_t newBinError    = histoCorrectedToBeScaled->GetBinError(i)/histoCorrectedToBeScaled->GetBinCenter(i)/histoCorrectedToBeScaled->GetBinWidth(i);
        histoCorrectedToBeScaled->SetBinContent(i,newBinContent);
        histoCorrectedToBeScaled->SetBinError(i,newBinError);
    }
    // no BR-ration correction needed
}


//**********************************************************************************
//********************** Main function for JJ properties ***************************
//**********************************************************************************
void  PlotJetJetMCProperties(   TString fileListInput   = "InputFile.txt",
                                TString cutSelection    = "",
                                Int_t mode              = 2,
                                Int_t numberOfBins      = 4,
                                TString suffix          = "eps",
                                TString optionEnergy    = "",
                                TString period          = "",
                                Bool_t additionalQA     = kFALSE,
                                Bool_t doubleLogplots   = kFALSE
                            ){

    //***************************************************************************************************************
    //************************************ Layouting preparations & general setup ***********************************
    //***************************************************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    TH1::AddDirectory(kFALSE);
    StyleSettingsThesis();
    SetPlotStyle();

    const Int_t MaxNumberOfFiles = 21;

    Color_t colorBins[21]       = { kBlack, kRed+2, kBlue+2, kGreen+2, kCyan+2,
                                    kViolet, kMagenta+2, kGray+1, kRed-2, kBlue-2,
                                    kViolet+2, kCyan-2, kTeal+2, kOrange+2, kSpring+2,
                                    kMagenta-2, kPink+2, kRed, kBlue, kTeal, kGreen};
    Color_t colorBinsShade[21]  = { kGray+1, kRed-6, kBlue-6, kGreen-8, kCyan-6,
                                    kViolet-8, kMagenta-8, kGray, kRed-8, kBlue-8,
                                    kViolet-6, kCyan-8, kTeal-6, 807, kSpring-6,
                                    kMagenta-9, kPink-8, kRed-9, kBlue-9, kTeal-7, kGreen-4 };
    Marker_t markerBins[21]     = { 20, 21, 33, 34, 29,
                                    24, 25, 27, 28, 30,
                                    20, 21, 33, 34, 29,
                                    24, 25, 27, 28, 30,
                                    20};

    TString outputDir =	Form("%s/%s/%s/JetJetMCProperties_%s", cutSelection.Data(), optionEnergy.Data(), suffix.Data(), period.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    TString fileNameInput                   [MaxNumberOfFiles];
    Float_t minPtHard                       [MaxNumberOfFiles];
    Float_t maxPtHard                       [MaxNumberOfFiles];

    //Variable defintion
    date                        = ReturnDateString();
    collisionSystem             = ReturnFullCollisionsSystem(optionEnergy);
    TString centralityString    = GetCentralityString(cutSelection.Data());
    fDetectionProcess           = ReturnFullTextReconstructionProcess(mode);

    Float_t floatLocationRightUp2D[4] = {0.45,0.95,0.035, 0.02};
    Float_t floatLocationLeftDown2D[4] = {0.15,0.25,0.035, 0.02};
    Float_t floatLocationRightDown2D[4] = {0.45,0.25,0.035, 0.02};

    Double_t maxPt      = 35;
    if (mode == 4 && period.Contains("LHC17g8")){
        maxPt           = 60;
    } else if (mode == 2){
        maxPt           = 30;
    } else if (mode == 10){
        maxPt           = 50;
    }
    cout << "going to maxPt: " << maxPt << endl;

    Double_t scaling                            = 1./(2.*TMath::Pi());
    Double_t recalcBarn                         = 1e-9;

    //***************************************************************************************************************
    //*************************** read setting from configuration file **********************************************
    //***************************************************************************************************************
    ifstream in(fileListInput.Data());
    cout<<"Available files:"<<endl;
    Int_t nrOfPtHardBins = 0;
    while(!in.eof() && nrOfPtHardBins<numberOfBins ){
        in >> fileNameInput[nrOfPtHardBins] >> minPtHard[nrOfPtHardBins] >> maxPtHard[nrOfPtHardBins];
        cout<< fileNameInput[nrOfPtHardBins]<< endl;
        nrOfPtHardBins++;
    }
// 	nrOfPtHardBins--;

    //***************************************************************************************************************
    //******************************** Load Pi0 histograms **********************************************************
    //***************************************************************************************************************
    TFile* fileInput                [MaxNumberOfFiles];
    TH1F* histoNEvents              [MaxNumberOfFiles];
    TH1F* histoNEventsWWeight       [MaxNumberOfFiles];
    TH1F* histoNTrials              [MaxNumberOfFiles];
    TProfile* profXSection          [MaxNumberOfFiles];
    TH1F* histoMCPi0Input           [MaxNumberOfFiles];
    TH1F* histoMCPi0InputW0EvtWeigth[MaxNumberOfFiles];
    TH1F* histoMCPi0InputAcc        [MaxNumberOfFiles];
    TH1F* histoMCEtaInput           [MaxNumberOfFiles];
    TH1F* histoMCEtaInputW0EvtWeigth[MaxNumberOfFiles];
    TH1F* histoMCEtaInputAcc        [MaxNumberOfFiles];
    TH2F* histoMCPi0vsJetPt         [MaxNumberOfFiles];
    TH2F* histoMCEtavsJetPt         [MaxNumberOfFiles];
    Double_t nTrials                [MaxNumberOfFiles];
    Double_t nGeneratedEvents       [MaxNumberOfFiles];
    Double_t nWeightedEvents        [MaxNumberOfFiles];
    Double_t xSection               [MaxNumberOfFiles];
    Double_t weight                 [MaxNumberOfFiles];
    Double_t weightApplied          [MaxNumberOfFiles];
    Double_t deltaRapid             [MaxNumberOfFiles];

    TString anchoredTo                  = "LHC11a";
    if (period.Contains("LHC13b4_fix"))
        anchoredTo                      = "LHC13[b-c]";
    if (period.Contains("LHC13b4_plus"))
        anchoredTo                      = "LHC13[d-e]";
    if (period.Contains("LHC15a3"))
        anchoredTo                      = "LHC13g";
    if (period.Contains("LHC16c2"))
        anchoredTo                      = "LHC12[c-i]";
    if (period.Contains("LHC16h3"))
        anchoredTo                      = "LHC15n";
    if (period.Contains("LHC17g8a"))
        anchoredTo                      = "LHC16qt";
    if (period.Contains("LHC17g8b"))
        anchoredTo                      = "LHC16r";
    if (period.Contains("LHC17g8c"))
        anchoredTo                      = "LHC16s";
    if (period.Contains("LHC18b8"))
        anchoredTo                      = "LHC17pq";

    TString acceptanceOf = "";
    if (mode == 0) acceptanceOf     = "|#eta_{#gamma}| < 0.9 (PCM acc.)";
    if (mode == 2) acceptanceOf     = "|#eta_{#gamma_{1}}| < 0.9,  #gamma_{2} in EMCal acc.";
    if (mode == 4) acceptanceOf     = "#gamma's in EMCal acc.";
    if (mode == 3) acceptanceOf     = "|#eta_{#gamma_{1}}| < 0.9,  #gamma_{2} in PHOS acc.";
    if (mode == 5) acceptanceOf     = "#gamma's in PHOS acc.";
    if (mode == 10) acceptanceOf    = "#gamma in EMCal acc.";
    if (mode == 11) acceptanceOf    = "#gamma in PHOS acc.";

    cout << __LINE__ << endl;
    for (Int_t i=0; i< nrOfPtHardBins; i++){
        // Define CutSelections
        TString fEventCutSelection                      = "";
        TString fGammaCutSelection                      = "";
        TString fClusterCutSelection                    = "";
        TString fElectronCutSelection                   = "";
        TString fMesonCutSelection                      = "";
        // disentangle cut selection
        ReturnSeparatedCutNumberAdvanced(cutSelection.Data(),fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);

        fileInput[i]                                = new TFile(fileNameInput[i]);
        if (fileInput[i]->IsZombie()) return;

        TString autoDetectedMainDir                 = AutoDetectMainTList(mode , fileInput[i]);
        if (autoDetectedMainDir.CompareTo("") == 0){
            cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
            return;
        }

        //************************** Container Loading ********************************************************************
        TList *TopDir                   = (TList*)fileInput[i]->Get(autoDetectedMainDir.Data());
        if(TopDir == NULL){
            cout<<"ERROR: TopDir not Found"<<endl;
            return;
        }
        TopDir->SetOwner();

        TString rapidityRange                       = "";
        deltaRapid[i]                               = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);


        TList *HistosGammaConversion                = (TList*)TopDir->FindObject(Form("Cut Number %s",cutSelection.Data()));
        if(HistosGammaConversion == NULL){
            cout<<"ERROR: " << Form("Cut Number %s",cutSelection.Data()) << " not Found in File"<<endl;
            return;
        }
        TList *ESDContainer                         = (TList*)HistosGammaConversion->FindObject(Form("%s ESD histograms",cutSelection.Data()));
        if(ESDContainer == NULL){
            cout<<"ERROR: " << Form("ESD histograms %s",cutSelection.Data()) << " not Found in File"<<endl;
            return;
        }
        histoNTrials[i]                             = (TH1F*)ESDContainer->FindObject("NTrials");
        histoNTrials[i]->SetName(Form("NTrials%d",i));
        nTrials[i]                                  = histoNTrials[i]->GetBinContent(1);
        histoNEvents[i]                             = (TH1F*)ESDContainer->FindObject("NEventsWOWeight");
        histoNEvents[i]->SetName(Form("NEventsWOWeight%d",i));
        histoNEventsWWeight[i]                      = (TH1F*)ESDContainer->FindObject("NEvents");
        histoNEventsWWeight[i]->SetName(Form("NEventsWWeight%d",i));

        nGeneratedEvents[i]                         = histoNEvents[i]->GetEntries();
        if (optionEnergy.Contains("PbPb") || optionEnergy.Contains("pPb"))
            nWeightedEvents[i]                      = histoNEventsWWeight[i]->GetBinContent(1);
        else
            nWeightedEvents[i]                      = GetNEvents(histoNEventsWWeight[i]);

        profXSection[i]                             = (TProfile*)ESDContainer->FindObject("XSection");
        profXSection[i]->SetName(Form("XSection%d",i));
        xSection[i]                                 = profXSection[i]->GetBinContent(1);

        weight[i]                                   = xSection[i]/(nTrials[i]/nGeneratedEvents[i]);
        weightApplied[i]                            = histoNEventsWWeight[i]->GetBinContent(1)/histoNEvents[i]->GetBinContent(1);

        TList *MCContainer                          = (TList*)HistosGammaConversion->FindObject(Form("%s MC histograms",cutSelection.Data()));
        if(MCContainer == NULL){
            cout<<"ERROR: " << Form("MC histograms %s",cutSelection.Data()) << " not Found in File"<<endl;
            return;
        }
        TH1F* temp = (TH1F*)MCContainer->FindObject("MC_Pi0_Pt");
        histoMCPi0Input[i]                          = new TH1F(*temp);
        histoMCPi0Input[i]->SetName(Form("MC_Pi0_Pt%d",i));
        temp = (TH1F*)MCContainer->FindObject("MC_Pi0InAcc_Pt");
        histoMCPi0InputAcc[i]                       = new TH1F(*temp);
        histoMCPi0InputAcc[i]->SetName(Form("MC_Pi0InAcc_Pt%d",i));
        temp = (TH1F*)MCContainer->FindObject("MC_Pi0_WOEventWeights_Pt");
        histoMCPi0InputW0EvtWeigth[i]               = new TH1F(*temp);
        histoMCPi0InputW0EvtWeigth[i]->SetName(Form("MC_Pi0_WOEventWeights_Pt%d",i));
        if (optionEnergy.CompareTo("pPb_5.023TeVRun2") == 0 || optionEnergy.CompareTo("13TeV") == 0){
            histoMCPi0Input[i]->Rebin(8);
            histoMCPi0InputAcc[i]->Rebin(8);
            histoMCPi0InputW0EvtWeigth[i]->Rebin(8);
        }
        if ( !(mode == 10 || mode == 11) ){
            temp = (TH1F*)MCContainer->FindObject("MC_Eta_Pt");
            histoMCEtaInput[i]                      = new TH1F(*temp);
            histoMCEtaInput[i]->SetName(Form("MC_Eta_Pt%d",i));
            temp = (TH1F*)MCContainer->FindObject("MC_EtaInAcc_Pt");
            histoMCEtaInputAcc[i]                   = new TH1F(*temp);
            histoMCEtaInputAcc[i]->SetName(Form("MC_EtaInAcc_Pt%d",i));
            temp = (TH1F*)MCContainer->FindObject("MC_Eta_WOEventWeights_Pt");
            histoMCEtaInputW0EvtWeigth[i]           = new TH1F(*temp);
            histoMCEtaInputW0EvtWeigth[i]->SetName(Form("MC_Eta_WOEventWeights_Pt%d",i));
            if (optionEnergy.CompareTo("pPb_5.023TeVRun2") == 0 || optionEnergy.CompareTo("13TeV") == 0){
                histoMCEtaInput[i]->Rebin(4);
                histoMCEtaInputAcc[i]->Rebin(4);
                histoMCEtaInputW0EvtWeigth[i]->Rebin(4);
            }


        }
        if (additionalQA){
            histoMCPi0vsJetPt[i] 							= (TH2F*)MCContainer->FindObject("MC_Pi0_Pt_JetPt");
            histoMCPi0vsJetPt[i]->SetName(Form("MC_Pi0_Pt_JetPt%d",i));
            if (!(mode == 10 || mode == 11) ){
                histoMCEtavsJetPt[i] 				= (TH2F*)MCContainer->FindObject("MC_Eta_Pt_JetPt");
                histoMCEtavsJetPt[i]->SetName(Form("MC_Eta_Pt_JetPt%d",i));
            }

            ptHardRange = Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[i], maxPtHard[i]);
            if (i==0) ptHardRange = "summed";

            Float_t maximumPi0 = FindLargestEntryIn2D(histoMCPi0vsJetPt[i]);
            Float_t minimumPi0 = FindSmallestEntryIn2D(histoMCPi0vsJetPt[i]);
            histoMCPi0vsJetPt[i]->GetZaxis()->SetRangeUser(minimumPi0,maximumPi0);
            PlotStandard2D( histoMCPi0vsJetPt[i] ,
                        Form("%s/MC_Pi0Pt_JetPt_%d.%s",outputDir.Data(),i,suffix.Data()),
                        "", "#it{p}_{#pi^{0},T} (GeV/#it{c})", "#it{p}_{jet,T} (GeV/#it{c})",
                        kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1,
                        floatLocationRightUp2D,500,500,"MC", period);
            if (!(mode == 10 || mode == 11)) {
                Float_t maximumEta = FindLargestEntryIn2D(histoMCEtavsJetPt[i]);
                Float_t minimumEta = FindSmallestEntryIn2D(histoMCEtavsJetPt[i]);
                histoMCEtavsJetPt[i]->GetZaxis()->SetRangeUser(minimumEta,maximumEta);
                PlotStandard2D( histoMCEtavsJetPt[i] ,
                            Form("%s/MC_EtaPt_JetPt_%d.%s",outputDir.Data(),i,suffix.Data()),
                            "", "#it{p}_{#eta},T} (GeV/#it{c})", "#it{p}_{jet,T} (GeV/#it{c})",
                            kFALSE, 30., 180., kFALSE, 0.01, 20., 0, 0, 1,
                            floatLocationRightUp2D,500,500,"MC", period);
            }
        }
        cout << ">>>>>>>> pT hard bin #: " << i << endl;
        cout << "ntrials: " <<  nTrials[i] << "\t xSection: " << xSection[i] << "\t number of generated events: " << nGeneratedEvents[i] << "\t weight: "
            << weight[i] << "\t weight applied: "<< weightApplied[i]<< endl;
        delete TopDir;

        fileInput[i]->Close();
        delete fileInput[i];

    }

    //***************************************************************************************************************
    //************************************Plotting unscaled inputs **************************************************
    //***************************************************************************************************************
    TCanvas* canvasInputUnscaled = new TCanvas("canvasInputUnscaled","",0,0,1350,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasInputUnscaled, 0.1, 0.015, 0.015, 0.07);
    canvasInputUnscaled->SetLogy();

    Float_t maximumPi0Unscaled = FindLargestEntryIn1D(histoMCPi0InputW0EvtWeigth[0])*10;
    if(period.Contains("LHC16c2")) maximumPi0Unscaled*=10;
    Float_t minimumPi0Unscaled = FindSmallestEntryIn1D(histoMCPi0InputW0EvtWeigth[nrOfPtHardBins-1]);

    TH2F * histo2DInputUnscaledPi0;
    histo2DInputUnscaledPi0 = new TH2F("histo2DInputUnscaledPi0","histo2DInputUnscaledPi0",1000,0., maxPt,10000,minimumPi0Unscaled,maximumPi0Unscaled);
    SetStyleHistoTH2ForGraphs(histo2DInputUnscaledPi0, "#it{p}_{T} (GeV/#it{c})","N_{#pi^{0}}",
                            0.032,0.04, 0.032,0.04, 0.8,1.1);
    histo2DInputUnscaledPi0->GetXaxis()->SetRangeUser(0,maxPt);
    histo2DInputUnscaledPi0->DrawCopy();

    TLegend* legendUnscaled = GetAndSetLegend2(0.2, 0.96-(1.*nrOfPtHardBins/2*0.028), 0.96, 0.97,22);
    legendUnscaled->SetMargin(0.12);
    legendUnscaled->SetNColumns(2);
    for (Int_t i = 0; i< nrOfPtHardBins; i++){
        DrawGammaSetMarker(histoMCPi0InputW0EvtWeigth[i], markerBins[i], 1., colorBins[i], colorBins[i]);
        histoMCPi0InputW0EvtWeigth[i]->DrawCopy("e1,same");
        if (i == 0) legendUnscaled->AddEntry(histoMCPi0InputW0EvtWeigth[i],"summed","p");
        else legendUnscaled->AddEntry(histoMCPi0InputW0EvtWeigth[i],Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[i], maxPtHard[i]),"p");
    }
    legendUnscaled->Draw();
    for (Int_t i = 0; i< nrOfPtHardBins; i++){
        histoMCPi0InputW0EvtWeigth[i]->DrawCopy("e1,same");
    }
    TLatex *labelMCName = new TLatex(0.45,0.99-(1.15*(nrOfPtHardBins/2+1)*0.032),Form("%s anchored to %s",period.Data(), anchoredTo.Data()));
    SetStyleTLatex( labelMCName, 32,4);
    labelMCName->SetTextFont(43);
    labelMCName->Draw();


    canvasInputUnscaled->Update();
    canvasInputUnscaled->SaveAs(Form("%s/Pi0_MC_InputUnscaled.%s",outputDir.Data(),suffix.Data()));


    if ( !(mode == 10 || mode == 11) ){
        Float_t maximumEtaUnscaled = FindLargestEntryIn1D(histoMCEtaInputW0EvtWeigth[0])*10;
        if(optionEnergy.CompareTo("8TeV")==0) maximumEtaUnscaled*=5;
        Float_t minimumEtaUnscaled = FindSmallestEntryIn1D(histoMCEtaInputW0EvtWeigth[nrOfPtHardBins-1]);

        TH2F * histo2DInputUnscaledEta;
        histo2DInputUnscaledEta = new TH2F("histo2DInputUnscaledEta","histo2DInputUnscaledEta",1000,0., maxPt,10000,minimumEtaUnscaled,maximumEtaUnscaled);
        SetStyleHistoTH2ForGraphs(histo2DInputUnscaledEta, "#it{p}_{T} (GeV/#it{c})","N_{#eta}",
                                0.032,0.04, 0.032,0.04, 0.8,1.1);
        histo2DInputUnscaledEta->GetXaxis()->SetRangeUser(0,maxPt);
        histo2DInputUnscaledEta->DrawCopy();

        legendUnscaled->Draw();
        for (Int_t i = 0; i< nrOfPtHardBins; i++){
            DrawGammaSetMarker(histoMCEtaInputW0EvtWeigth[i], markerBins[i], 1., colorBins[i], colorBins[i]);
            histoMCEtaInputW0EvtWeigth[i]->DrawCopy("e1,same");
        }

        labelMCName->Draw();

        canvasInputUnscaled->Update();
        canvasInputUnscaled->SaveAs(Form("%s/Eta_MC_InputUnscaled.%s",outputDir.Data(),suffix.Data()));
    }

    if(doubleLogplots){
        canvasInputUnscaled->cd();
        canvasInputUnscaled->SetLogy();
        canvasInputUnscaled->SetLogx();
        histo2DInputUnscaledPi0->GetYaxis()->SetRangeUser(1,2e6);
        histo2DInputUnscaledPi0->DrawCopy();

        TLegend* legendUnscaledDLog = GetAndSetLegend2(0.12, 0.1, 0.96, 0.1+(1.*nrOfPtHardBins*0.028),22);
        legendUnscaledDLog->SetMargin(0.12);
    //     legendUnscaledDLog->SetNColumns(2);
        for (Int_t i = 1; i< nrOfPtHardBins; i++){
            histoMCPi0InputW0EvtWeigth[i]->DrawCopy("e1,same");
            legendUnscaledDLog->AddEntry(histoMCPi0InputW0EvtWeigth[i],Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[i], maxPtHard[i]),"p");
        }
        legendUnscaledDLog->Draw();

        TLatex *labelMCNameDLog = new TLatex(0.5,0.93,Form("%s anchored to %s",period.Data(), anchoredTo.Data()));
        SetStyleTLatex( labelMCNameDLog, 32,4);
        labelMCNameDLog->SetTextFont(43);
        labelMCNameDLog->Draw();

        canvasInputUnscaled->Update();
        canvasInputUnscaled->SaveAs(Form("%s/Pi0_MC_InputUnscaledDLog.%s",outputDir.Data(),suffix.Data()));


        if ( !(mode == 10 || mode == 11) ){

            canvasInputUnscaled->cd();
            canvasInputUnscaled->SetLogy();
            canvasInputUnscaled->SetLogx();

            Float_t maximumEtaUnscaled = FindLargestEntryIn1D(histoMCEtaInputW0EvtWeigth[0])*10;
            if(optionEnergy.CompareTo("8TeV")==0) maximumEtaUnscaled*=5;
            Float_t minimumEtaUnscaled = FindSmallestEntryIn1D(histoMCEtaInputW0EvtWeigth[nrOfPtHardBins-1]);
            TH2F * histo2DInputUnscaledEta;
            histo2DInputUnscaledEta = new TH2F("histo2DInputUnscaledEta","histo2DInputUnscaledEta",1000,0., maxPt,10000,minimumEtaUnscaled,maximumEtaUnscaled);
            SetStyleHistoTH2ForGraphs(histo2DInputUnscaledEta, "#it{p}_{T} (GeV/#it{c})","N_{#eta}",
                                    0.032,0.04, 0.032,0.04, 0.8,1.1);
            histo2DInputUnscaledEta->GetXaxis()->SetRangeUser(0,maxPt);
            histo2DInputUnscaledEta->GetYaxis()->SetRangeUser(1,3e4);
            histo2DInputUnscaledEta->DrawCopy();

            for (Int_t i = 1; i< nrOfPtHardBins; i++){
                histoMCEtaInputW0EvtWeigth[i]->DrawCopy("e1,same");
            }
            legendUnscaledDLog->Draw();
            labelMCNameDLog->Draw();

            canvasInputUnscaled->Update();
            canvasInputUnscaled->SaveAs(Form("%s/Eta_MC_InputUnscaledDLog.%s",outputDir.Data(),suffix.Data()));
        }
    }
    delete canvasInputUnscaled;

    //***************************************************************************************************************
    //************************************Plotting scaled inputs **************************************************
    //***************************************************************************************************************
    TCanvas* canvasInputScaled = new TCanvas("canvasInputScaled","",0,0,1350,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasInputScaled, 0.1, 0.015, 0.015, 0.07);
    canvasInputScaled->SetLogy();
    Float_t maximumPi0Scaled = FindLargestEntryIn1D(histoMCPi0Input[0])*10;
    Float_t minimumPi0Scaled = FindSmallestEntryIn1D(histoMCPi0Input[nrOfPtHardBins-1]);

    TH2F * histo2DInputScaledPi0;
    histo2DInputScaledPi0 = new TH2F("histo2DInputScaledPi0","histo2DInputScaledPi0",1000,0., maxPt,10000,minimumPi0Scaled,maximumPi0Scaled);
    SetStyleHistoTH2ForGraphs(histo2DInputScaledPi0, "#it{p}_{T} (GeV/#it{c})","N_{#pi^{0}} reweighted",
                            0.032,0.04, 0.032,0.04, 0.8,1.2);
    histo2DInputScaledPi0->GetXaxis()->SetRangeUser(0,maxPt);
    histo2DInputScaledPi0->DrawCopy();

    TLegend* legendScaled = GetAndSetLegend2(0.2, 0.96-(1.*nrOfPtHardBins/2*0.028), 0.95, 0.96,22);
    legendScaled->SetMargin(0.12);
    legendScaled->SetNColumns(2);
    for (Int_t i = 0; i< nrOfPtHardBins; i++){
        DrawGammaSetMarker(histoMCPi0Input[i], markerBins[i], 1., colorBins[i], colorBins[i]);
        histoMCPi0Input[i]->DrawCopy("e1,same");
        if (i == 0) legendScaled->AddEntry(histoMCPi0Input[i],"summed","p");
        else legendScaled->AddEntry(histoMCPi0Input[i],Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[i], maxPtHard[i]),"p");
    }
    legendScaled->Draw();
    for (Int_t i = 0; i< nrOfPtHardBins; i++){
    histoMCPi0Input[i]->DrawCopy("e1,same");
    }
    labelMCName->Draw();


    canvasInputScaled->Update();
    canvasInputScaled->SaveAs(Form("%s/Pi0_MC_InputScaled.%s",outputDir.Data(),suffix.Data()));

    histo2DInputScaledPi0->DrawCopy();
    legendScaled->Draw();
    for (Int_t i = 0; i< nrOfPtHardBins; i++){
        DrawGammaSetMarker(histoMCPi0InputAcc[i], markerBins[i], 1., colorBins[i], colorBins[i]);
        histoMCPi0InputAcc[i]->DrawCopy("e1,same");
    }
    labelMCName->Draw();
    TLatex *labelAcceptance = new TLatex(0.45,0.99-(1.15*(nrOfPtHardBins/2+2)*0.032),Form("%s",acceptanceOf.Data()));
    SetStyleTLatex( labelAcceptance, 32,4);
    labelAcceptance->SetTextFont(43);
    labelAcceptance->Draw();


    canvasInputScaled->Update();
    canvasInputScaled->SaveAs(Form("%s/Pi0_MC_InputScaledInAcceptance.%s",outputDir.Data(),suffix.Data()));


    if ( !(mode == 10 || mode == 11) ){
        Float_t maximumEtaScaled = FindLargestEntryIn1D(histoMCEtaInput[0])*10;
        if(optionEnergy.CompareTo("8TeV")==0) maximumEtaScaled*=2;
        Float_t minimumEtaScaled = FindSmallestEntryIn1D(histoMCEtaInput[nrOfPtHardBins-1]);

        TH2F * histo2DInputScaledEta;
        histo2DInputScaledEta = new TH2F("histo2DInputScaledEta","histo2DInputScaledEta",1000,0., maxPt,10000,minimumEtaScaled,maximumEtaScaled);
        SetStyleHistoTH2ForGraphs(histo2DInputScaledEta, "#it{p}_{T} (GeV/#it{c})","N_{#eta} reweighted",
                                0.032,0.04, 0.032,0.04, 0.8,1.1);
        histo2DInputScaledEta->GetXaxis()->SetRangeUser(0,maxPt);
        histo2DInputScaledEta->DrawCopy();

        legendScaled->Draw();
        for (Int_t i = 0; i< nrOfPtHardBins; i++){
            DrawGammaSetMarker(histoMCEtaInput[i], markerBins[i], 1., colorBins[i], colorBins[i]);
            histoMCEtaInput[i]->DrawCopy("e1,same");
        }

        labelMCName->Draw();

        canvasInputScaled->Update();
        canvasInputScaled->SaveAs(Form("%s/Eta_MC_InputScaled.%s",outputDir.Data(),suffix.Data()));

        histo2DInputScaledEta->DrawCopy();
        legendScaled->Draw();
        for (Int_t i = 0; i< nrOfPtHardBins; i++){
            DrawGammaSetMarker(histoMCEtaInputAcc[i], markerBins[i], 1., colorBins[i], colorBins[i]);
            histoMCEtaInputAcc[i]->DrawCopy("e1,same");
        }

        labelMCName->Draw();
        labelAcceptance->Draw();

        canvasInputScaled->Update();
        canvasInputScaled->SaveAs(Form("%s/Eta_MC_InputScaledInAcceptance.%s",outputDir.Data(),suffix.Data()));
    }

    if(doubleLogplots){

        canvasInputScaled->cd();
        canvasInputScaled->SetLogy();
        canvasInputScaled->SetLogx();
//         histo2DInputScaledPi0->GetYaxis()->SetRangeUser(1e-3,1e8);
        histo2DInputScaledPi0->DrawCopy();

            TLegend* legendScaledDLog = GetAndSetLegend2(0.2, 0.96-(1.*nrOfPtHardBins/2*0.028), 0.95, 0.96,22);
            legendScaledDLog->SetNColumns(2);
            legendScaledDLog->SetMargin(0.12);
            for (Int_t i = 1; i< nrOfPtHardBins; i++){
                histoMCPi0Input[i]->DrawCopy("e1,same");
                legendScaledDLog->AddEntry(histoMCPi0Input[i],Form("%1.0f GeV/#it{c} < #it{p}^{hard}_{T} < %1.0f GeV/#it{c}", minPtHard[i], maxPtHard[i]),"p");
            }
//             legendScaledDLog->Draw();
            TLatex *labelMCNameDLog = new TLatex(0.5,0.93,Form("%s anchored to %s",period.Data(), anchoredTo.Data()));
            SetStyleTLatex( labelMCNameDLog, 32,4);
            labelMCNameDLog->SetTextFont(43);
            labelMCNameDLog->Draw();

        canvasInputScaled->Update();
        canvasInputScaled->SaveAs(Form("%s/Pi0_MC_InputScaledDLog.%s",outputDir.Data(),suffix.Data()));

        canvasInputScaled->cd();
        histo2DInputScaledPi0->DrawCopy();
//             legendScaledDLog->Draw();
            for (Int_t i = 1; i< nrOfPtHardBins; i++){
                histoMCPi0InputAcc[i]->DrawCopy("e1,same");
            }
            labelMCNameDLog->Draw();
    //         TLatex *labelAcceptance = new TLatex(0.45,0.99-(1.15*(nrOfPtHardBins/2+2)*0.032),Form("%s",acceptanceOf.Data()));
    //         SetStyleTLatex( labelAcceptance, 32,4);
    //         labelAcceptance->SetTextFont(43);
            labelAcceptance->Draw();

        canvasInputScaled->Update();
        canvasInputScaled->SaveAs(Form("%s/Pi0_MC_InputScaledInAcceptanceDLog.%s",outputDir.Data(),suffix.Data()));

        if ( !(mode == 10 || mode == 11) ){

            canvasInputScaled->cd();
            Float_t maximumEtaScaled = FindLargestEntryIn1D(histoMCEtaInput[0])*10;
            if(optionEnergy.CompareTo("8TeV")==0) maximumEtaScaled*=2;
            Float_t minimumEtaScaled = FindSmallestEntryIn1D(histoMCEtaInput[nrOfPtHardBins-1]);

            TH2F * histo2DInputScaledEta;
            histo2DInputScaledEta = new TH2F("histo2DInputScaledEta","histo2DInputScaledEta",1000,0., maxPt,10000,minimumEtaScaled,maximumEtaScaled);
            SetStyleHistoTH2ForGraphs(histo2DInputScaledEta, "#it{p}_{T} (GeV/#it{c})","N_{#eta} reweighted",
                                    0.032,0.04, 0.032,0.04, 0.8,1.1);
            histo2DInputScaledEta->GetXaxis()->SetRangeUser(0,maxPt);
//             histo2DInputScaledEta->GetYaxis()->SetRangeUser(1e-2,1e7);
            histo2DInputScaledEta->DrawCopy();

//                 legendScaledDLog->Draw();
                for (Int_t i = 1; i< nrOfPtHardBins; i++){
                    histoMCEtaInput[i]->DrawCopy("e1,same");
                }
                labelMCNameDLog->Draw();

                canvasInputScaled->Update();
                canvasInputScaled->SaveAs(Form("%s/Eta_MC_InputScaledDlog.%s",outputDir.Data(),suffix.Data()));

                histo2DInputScaledEta->DrawCopy();
//                 legendScaledDLog->Draw();
                for (Int_t i = 1; i< nrOfPtHardBins; i++){
                    histoMCEtaInputAcc[i]->DrawCopy("e1,same");
                }

                labelMCNameDLog->Draw();
                labelAcceptance->Draw();

            canvasInputScaled->Update();
            canvasInputScaled->SaveAs(Form("%s/Eta_MC_InputScaledInAcceptanceDlog.%s",outputDir.Data(),suffix.Data()));
        }
    }
    delete canvasInputScaled;

    //**************************************************************************************************************
    //*************************** Cross sections for pi0 and eta **************************************************
    //**************************************************************************************************************
    TH1D* histoPi0InvYield  [MaxNumberOfFiles];
    TH1D* histoEtaInvYield  [MaxNumberOfFiles];
    for (Int_t i=0; i< nrOfPtHardBins; i++){
        histoPi0InvYield[i]     = (TH1D*)histoMCPi0Input[i]->Clone(Form("MC_Pi0_InvYield_%i", i));
        ScaleMCYield( histoPi0InvYield[i], deltaRapid[i], scaling, nWeightedEvents[i], "Pi0");
        if ( !(mode == 10 || mode == 11) ){
            histoEtaInvYield[i]     = (TH1D*)histoMCEtaInput[i]->Clone(Form("MC_Eta_InvYield_%i", i));
            ScaleMCYield( histoEtaInvYield[i], deltaRapid[i], scaling, nWeightedEvents[i], "Eta");
        }
    }

    //**************************************************************************************************************
    //*************************** Write all pt hard bins in 1 file *************************************************
    //**************************************************************************************************************
    TString nameOutput  = Form("%s/CompilationJJMCBins.root",outputDir.Data());
    TFile* output       = new TFile(nameOutput.Data(),"RECREATE");
        for (Int_t i=0; i< nrOfPtHardBins; i++){
            histoMCPi0InputW0EvtWeigth[i]->Write();
            histoMCPi0Input[i]->Write();
            histoPi0InvYield[i]->Write();
            if ( !(mode == 10 || mode == 11) ){
                histoMCEtaInputW0EvtWeigth[i]->Write();
                histoMCEtaInput[i]->Write();
                histoEtaInvYield[i]->Write();
            }
        }
    output->Write();
    output->Close();


}
