// provided by Gamma Conversion Group, PWGGA/GammaConv
//A. Marin. July 2018. First maps shown by Hikari
// Takes the Maps created with the MaterialHistos Task and fits mean and width of the DeDx sigmam distribution for electrons and positrons.
// The mean is stored in eta vs p 2D histograms for 4 Radial bins

#include "TH3F.h"
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
#include <TLatex.h>
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "TMinuit.h"
#include "../CommonHeaders/PlottingMeson.h"
#include "TASImage.h"
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/ConversionFunctions.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"

TF1* FitSignal(TH1D* hSig,Color_t col);
void PlotdEdxSlices(TH1D** , TF1** , Double_t* , TString , TString , TString , Double_t* , TString , TString , Int_t , Int_t , Int_t, TString, TString );
void DrawSliceHisto(TH1* ,  TString , TString , TString , Float_t , Float_t , Size_t , Color_t );

TF1*     fGaus       = NULL;
Double_t Mean=0.0;
Double_t Width=0.0;
Double_t Chi2=0.0;
Double_t meanElectron[4][12][20];//eta pt
Double_t meanPositron[4][12][20];//eta pt
Double_t meanErrElectron[4][12][20];//eta pt
Double_t meanErrPositron[4][12][20];//eta pt
Double_t widthElectron[4][12][20];
Double_t widthPositron[4][12][20];
Double_t widthErrElectron[4][12][20];
Double_t widthErrPositron[4][12][20];
Double_t chi[20][20];

void CreateDEdxMaps(    TString fileNameWithMaps    ="" ,
                        TString cutSelection        = "",
                        TString optionMC            = "Data",
                        TString fEnergy             = "",
                        TString suffix              = "pdf",
                        TString period              = ""
)
{
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetEndErrorSize(0);
    StyleSettingsThesis();
    SetPlotStyle();

    TString outputDir = "";
    if(optionMC.CompareTo("Data")==0)    outputDir = "DeDxMapsData/";
    else if(optionMC.CompareTo("MC")==0) outputDir = "DeDxMapsMC/";
    if (period.CompareTo("") != 0) outputDir = outputDir+period+"/";
    gSystem->Exec(Form("mkdir -p %s/%s",outputDir.Data(),"DetailedFits"));

    TH3F *histoPositronDeDxPEta[4] = {NULL};
    TH3F *histoElectronDeDxPEta[4] = {NULL};

    TFile* fileMaterialHistos = new TFile(fileNameWithMaps.Data());
    TString fCutSelectionRead = cutSelection;
    TString nameMainDir = "GammaConvMaterial";

    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergy);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    TList *TopDir =(TList*)fileMaterialHistos->Get(nameMainDir.Data());
    if(TopDir != NULL){
        TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
        if(HistosGammaConversion == NULL){
            cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
            return;
        }

        TList *MapsContainer           = (TList*)HistosGammaConversion->FindObject(Form("%s  dEdx Maps",fCutSelectionRead.Data()));
        for (Int_t i = 0; i<4; i++){
            histoPositronDeDxPEta[i] = (TH3F*)MapsContainer->FindObject(Form("R%d positron sigma dEdx P Eta",i));
            histoElectronDeDxPEta[i] = (TH3F*)MapsContainer->FindObject(Form("R%d electron sigma dEdx P Eta",i));
        }
    } else if (TopDir == NULL){
        cout<<"WARNING: TopDir " << nameMainDir.Data() << " not Found... checking for PhotonQA output"<<endl;

        TString nameDirectory           = Form("GammaConvV1_QA_%s",  fCutSelectionRead.Data());
        TDirectory* directoryConv           = (TDirectory*)fileMaterialHistos->Get(nameDirectory.Data());
        if (directoryConv == NULL){
            cout<<"ERROR: PhotonQA directory " << nameDirectory.Data() << " not Found! Returning..."<<endl;
            return;
        } else {
            cout<<"INFO: PhotonQA directory " << nameDirectory.Data() << " found!"<<endl;
        }
        for (Int_t i = 0; i<4; i++){
            histoPositronDeDxPEta[i] = (TH3F*)directoryConv->Get(Form("R%d positron sigma dEdx P Eta",i));
            histoElectronDeDxPEta[i] = (TH3F*)directoryConv->Get(Form("R%d electron sigma dEdx P Eta",i));
        }
    }

    Color_t colorR[4]           = { kBlack, kBlue+2, kRed+2, kGreen+2};
    Color_t colorEta[18]        = { kBlack, kGray+2, kBlue-7, kBlue-4, kBlue+2, kCyan+2, kCyan-5, kGreen-7, kGreen+2, kYellow+1,
                                    kOrange-5, kOrange-3, kRed+2, kRed-4, kRed-6, kMagenta-5, kMagenta-2, kViolet-7};
    Color_t colorPt[12]         = { kBlack, kGray+2, kBlue-7, kBlue+2, kCyan+2, kCyan-5, kGreen-7, kGreen+2, kOrange-3, kRed+2,
                                    kRed-6, kMagenta-2};
    Marker_t markerR[8]         = { 20, 21, 33, 34, 24, 25, 27, 28};
    Marker_t markerEta[18]      = { 24, 25, 27, 28, 30, 29, 34, 33, 21, 20,
                                    42, 46, 40, 41, 27, 43, 24, 25};
    Marker_t markerPt[12]       = { 24, 25, 27, 28, 30, 29, 34, 33, 21, 20,
                                    42, 46};

    Int_t nPBins   = 12;
    Double_t *arrPBinning = new Double_t[13];
    for( Int_t i=0;i<nPBins+1;i++){
        if(i==0){
            arrPBinning[i]= 0.05;
        } else if(i>0 && i<11){
            arrPBinning[i]= 0.1*i;
        } else if(i==11){
            arrPBinning[i]= 2.0;
        } else if(i==12){
            arrPBinning[i]= 10.0;
        }
        cout << "p bins:: "<< i << " " <<  arrPBinning[i]<< endl;
    }

    Int_t nEtaBins      = 18;
    Int_t nEtaBinsOut   = 18;
    Double_t *arrEtaBinningOut = new Double_t[19];
    for( Int_t i=0;i<nEtaBinsOut+1;i++){
        arrEtaBinningOut[i]= -0.9+0.1*i;
        cout << "eta bins:: " << i << " " << arrEtaBinningOut[i] << endl;
    }

    Int_t nBinsR        = 4;
    Double_t arrRBinningOut[5] = {0, 33.5, 72, 145, 180};
    for (Int_t i=0;i<nBinsR+1;i++){
        cout << "R bins:: " << i << " " << arrRBinningOut[i] << endl;
    }

    TH2F* fhistoMeanPositron[4]     = {NULL};
    TH2F* fhistoMeanElectron[4]     = {NULL};
    TH2F* fhistoWidthElectron[4]    = {NULL};
    TH2F* fhistoWidthPositron[4]    = {NULL};

    for (Int_t i = 0; i<4; i++){
        fhistoMeanPositron[i]       = new TH2F(Form("MeanPosiR%d",i),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoMeanPositron[i]->Sumw2();
        fhistoMeanElectron[i]       = new TH2F(Form("MeanElecR%d",i),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoMeanElectron[i]->Sumw2();
        fhistoWidthPositron[i]      = new TH2F(Form("WidthPosiR%d",i),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoWidthPositron[i]->Sumw2();
        fhistoWidthElectron[i]      = new TH2F(Form("WidthElecR%d",i),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoWidthElectron[i]->Sumw2();
    }

    TH1D* histoPositronDeDx[4][12][20]  = {{{NULL}}};
    TH1D* histoElectronDeDx[4][12][20]  = {{{NULL}}};
    TF1* fitPositronDeDx[4][12][20]     = {{{NULL}}};
    TF1* fitElectronDeDx[4][12][20]     = {{{NULL}}};

    Int_t etaOff = 0;
//     if(fEnergy.CompareTo("5TeV2017")==0) etaOff = 0;

    //    cout<< histoPositronDeDxPEta0->GetNbinsY()<< " " << histoPositronDeDxPEta0->GetNbinsZ() << endl;
    for(Int_t i=0;i<nPBins;i++){

        Int_t minBinPt    = histoElectronDeDxPEta[0]->GetZaxis()->FindBin(arrPBinning[i]);
        Int_t maxBinPt    = histoElectronDeDxPEta[0]->GetZaxis()->FindBin(arrPBinning[i+1]-0.001);
        cout << "*****************\np = " << arrPBinning[i] << " GeV/c" << "< pT < " << arrPBinning[i+1] << " GeV/c"<< endl;
        cout << minBinPt << "\t" << maxBinPt << endl;
        for(Int_t j=0;j<nEtaBins;j++){
//             cout << "-> eta " << arrEtaBinningOut[j] << endl;
//             cout << "offset " << arrEtaBinningOut[j-etaOff] << endl;
//             cout << " j-etaOff " << j-etaOff << endl;
            Int_t minBinEta    = histoElectronDeDxPEta[0]->GetYaxis()->FindBin(arrEtaBinningOut[j]);
            Int_t maxBinEta    = histoElectronDeDxPEta[0]->GetYaxis()->FindBin(arrEtaBinningOut[j+1]-0.001);

            for (Int_t r = 0; r < 4; r++){
                // Electron, R bin 0
                histoElectronDeDx[r][i][j]  = (TH1D*)histoElectronDeDxPEta[r]->ProjectionX(Form("R%dElecSigdEdxP%dEta%d",r,i,j),minBinEta,maxBinEta,minBinPt,maxBinPt);
                if (histoElectronDeDx[r][i][j]->GetMaximum() < 200) histoElectronDeDx[r][i][j]->Rebin(2);
                if (histoElectronDeDx[r][i][j]->GetMaximum() < 100) histoElectronDeDx[r][i][j]->Rebin(2);
//                 fitElectronDeDx[r][i][j]    = FitSignal(histoElectronDeDx[r][i][j],kBlue);

                fitElectronDeDx[r][i][j]   = FitTH1DRecursivelyGaussian (histoElectronDeDx[r][i][j], 0.02, -3, 3, 2, 1.25);
                if (fitElectronDeDx[r][i][j]){
                    meanElectron[r][i][j]       = fitElectronDeDx[r][i][j]->GetParameter(1);
                    widthElectron[r][i][j]      = fitElectronDeDx[r][i][j]->GetParameter(2);
                    meanErrElectron[r][i][j]    = fitElectronDeDx[r][i][j]->GetParError(1);
                    widthErrElectron[r][i][j]   = fitElectronDeDx[r][i][j]->GetParError(2);
                } else {
                    meanElectron[r][i][j]       = 0;
                    widthElectron[r][i][j]      = 1;
                    meanErrElectron[r][i][j]    = 1;
                    widthErrElectron[r][i][j]   = 0.5;
                }
                if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                    fhistoMeanElectron[r]->SetBinContent(i+1,j+1-etaOff, meanElectron[r][i][j]);
                    fhistoMeanElectron[r]->SetBinError(i+1,j+1-etaOff, meanErrElectron[r][i][j]);
                    fhistoWidthElectron[r]->SetBinContent(i+1,j+1-etaOff, widthElectron[r][i][j]);
                    fhistoWidthElectron[r]->SetBinError(i+1,j+1-etaOff, widthErrElectron[r][i][j]);
                }

                // Positron, R bin 0
                histoPositronDeDx[r][i][j] = (TH1D*)histoPositronDeDxPEta[r]->ProjectionX(Form("R%dPosiSigdEdxP%dEta%d",r,i,j),minBinEta,maxBinEta,minBinPt,maxBinPt);
                if (histoPositronDeDx[r][i][j]->GetMaximum() < 200) histoPositronDeDx[r][i][j]->Rebin(2);
                if (histoPositronDeDx[r][i][j]->GetMaximum() < 100) histoPositronDeDx[r][i][j]->Rebin(2);
//                 fitPositronDeDx[r][i][j]    = FitSignal(histoPositronDeDx[r][i][j],kBlue);
                fitPositronDeDx[r][i][j]   = FitTH1DRecursivelyGaussian (histoPositronDeDx[r][i][j], 0.02, -3, 3, 2, 1.25 );
                if (fitElectronDeDx[r][i][j]){
                    meanPositron[r][i][j]       = fitPositronDeDx[r][i][j]->GetParameter(1);
                    widthPositron[r][i][j]      = fitPositronDeDx[r][i][j]->GetParameter(2);
                    meanErrPositron[r][i][j]    = fitPositronDeDx[r][i][j]->GetParError(1);
                    widthErrPositron[r][i][j]   = fitPositronDeDx[r][i][j]->GetParError(2);
                } else {
                    meanPositron[r][i][j]   = 0;
                    widthPositron[r][i][j]  = 1;
                    meanErrPositron[r][i][j]    = 1;
                    widthErrPositron[r][i][j]   = 0.5;
                }
                if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                    fhistoMeanPositron[r]->SetBinContent(i+1,j+1-etaOff, meanPositron[r][i][j]);
                    fhistoMeanPositron[r]->SetBinError(i+1,j+1-etaOff, meanErrPositron[r][i][j]);
                    fhistoWidthPositron[r]->SetBinContent(i+1,j+1-etaOff, widthPositron[r][i][j]);
                    fhistoWidthPositron[r]->SetBinError(i+1,j+1-etaOff, widthErrPositron[r][i][j]);
                }
            }
        }
    }



    Int_t textSizeLabelsPixel                   = 1000*0.04;
    TCanvas* canvasdEdxMeanStudies = new TCanvas("canvasdEdxMeanStudies","",200,10,1350,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasdEdxMeanStudies, 0.08, 0.02, 0.02, 0.09);
    canvasdEdxMeanStudies->SetLogx();
    TH2F * histo2DdEdxMeanStudies = new TH2F("histo2DdEdxMeanStudies","histo2DdEdxMeanStudies", 100, 0.04, 20, 200, -1.5, 1.5);
    SetStyleHistoTH2ForGraphs(histo2DdEdxMeanStudies, "#it{p} (GeV/#it{c})", "<n#sigma_{e}>", 0.035,0.04, 0.035,0.04, 0.98,1.0);

    TLegend* legendMeanEta  = GetAndSetLegend2(0.3, 0.11, 0.95, 0.11+(0.035*nEtaBinsOut/4.*1.05), 0.75*textSizeLabelsPixel,4,"",43,0.1);
    TLegend* legendMeanR    = GetAndSetLegend2(0.27, 0.11, 0.95, 0.11+(0.035*4.*1.05), 0.75*textSizeLabelsPixel,2,"",43,0.1);
    TCanvas* canvasdEdxWidthStudies = new TCanvas("canvasdEdxWidthStudies","",200,10,1350,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasdEdxWidthStudies, 0.08, 0.02, 0.02, 0.09);
    canvasdEdxWidthStudies->SetLogx();
    TH2F * histo2DdEdxWidthStudies = new TH2F("histo2DdEdxWidthStudies","histo2DdEdxWidthStudies", 100, 0.04, 20, 200, 0, 2);
    SetStyleHistoTH2ForGraphs(histo2DdEdxWidthStudies, "#it{p} (GeV/#it{c})", "#sigma(n#sigma_{e})", 0.035,0.04, 0.035,0.04, 0.98,1.0);

    TLegend* legendWidthEta  = GetAndSetLegend2(0.3, 0.11, 0.95, 0.11+(0.035*nEtaBinsOut/4.*1.05), 0.75*textSizeLabelsPixel,4,"",43,0.1);
    TLegend* legendWidthR    = GetAndSetLegend2(0.27, 0.11, 0.95, 0.11+(0.035*4.*1.05), 0.75*textSizeLabelsPixel,2,"",43,0.1);


    //   _____ _______                 ____
    //  |  ___|__   __|  /\           |  _ \
    //  | |_     | |    /  \          | |_| |
    //  |  _|    | |   / /\ \         |    /
    //  | |___   | |  / ___- \        | |\ \
    //  |_____|  |_| /_/    \_\       |_| \_\

    for (Int_t r = 0; r < 4; r++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanEta->Clear();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoMeanElectron[r]->ProjectionX(Form("%sEta%d",fhistoMeanElectron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
            legendMeanEta->AddEntry(histoProjection,Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]),"p");
        }
        legendMeanEta->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sElectronMeanNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoMeanPositron[r]->ProjectionX(Form("%sEta%d",fhistoMeanPositron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
        }
        legendMeanEta->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sPositronMeanNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    for (Int_t r = 0; r < 4; r++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthEta->Clear();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoWidthElectron[r]->ProjectionX(Form("%sEta%d",fhistoWidthElectron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
            legendWidthEta->AddEntry(histoProjection,Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]),"p");
        }
        legendWidthEta->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sElectronWidthNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoWidthPositron[r]->ProjectionX(Form("%sEta%d",fhistoWidthPositron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
        }
        legendWidthEta->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sPositronWidthNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    //   ____         _____ _______
    //  |  _ \       |  ___|__   __|  /\
    //  | |_| |      | |_     | |    /  \
    //  |    /       |  _|    | |   / /\ \
    //  | |\ \       | |___   | |  / ___- \
    //  |_| \_\      |_____|  |_| /_/    \_\


    for (Int_t j = 0; j< nEtaBinsOut; j++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]), fCollisionSystem, "", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sMeanNSigma_SummaryPtDepInRBinsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }
    for (Int_t j = 0; j< nEtaBinsOut; j++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]), fCollisionSystem, "", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sWidthNSigma_SummaryPtDepInRBinsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }

    for (Int_t j = 0; j< nEtaBinsOut/2; j++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanElectron[r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sElectronMeanNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));

        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanPositron[r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sPositronMeanNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }

    for (Int_t j = 0; j< nEtaBinsOut/2; j++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthElectron[r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sElectronWidthNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));

        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthPositron[r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sPositronWidthNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }

    //   ____       ____
    //  |  _ \     |  _ \
    //  | |_| |    | |_| |
    //  |  __/     |    /
    //  | |        | |\ \
    //  |_|        |_| \_\

    TH2F * histo2DdEdxMeanStudiesEta = new TH2F("histo2DdEdxMeanStudiesEta","histo2DdEdxMeanStudiesEta", 100, -1, 1, 200, -1.7, 1.7);
    SetStyleHistoTH2ForGraphs(histo2DdEdxMeanStudiesEta, "#eta", "<n#sigma_{e}>", 0.035,0.04, 0.035,0.04, 0.98,1.0);
    TH2F * histo2DdEdxWidthStudiesEta = new TH2F("histo2DdEdxWidthStudiesEta","histo2DdEdxWidthStudiesEta", 100, -1, 1, 200, 0, 2);
    SetStyleHistoTH2ForGraphs(histo2DdEdxWidthStudiesEta, "#eta", "#sigma(n#sigma_{e})", 0.035,0.04, 0.035,0.04, 0.98,1.0);

    TLegend* legendMeanPt  = GetAndSetLegend2(0.12, 0.11, 0.95, 0.11+(0.035*nPBins/3.*1.05), 0.75*textSizeLabelsPixel,3,"",43,0.1);
    TLegend* legendWidthPt  = GetAndSetLegend2(0.12, 0.11, 0.95, 0.11+(0.035*nPBins/3.*1.05), 0.75*textSizeLabelsPixel,3,"",43,0.1);

    canvasdEdxMeanStudies->SetLogx(kFALSE);
    for (Int_t r = 0; r < 4; r++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        legendMeanPt->Clear();
        for (Int_t i = 0; i < nPBins; i++){
            TH1D* histoProjection       = fhistoMeanElectron[r]->ProjectionY(Form("%sPt%d",fhistoMeanElectron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
            legendMeanPt->AddEntry(histoProjection,Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),"p");
        }
        legendMeanPt->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sElectronMeanNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        for (Int_t i = 0; i < nPBins; i++){
            TH1D* histoProjection       = fhistoMeanPositron[r]->ProjectionY(Form("%sPt%d",fhistoMeanPositron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection,  markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
        }
        legendMeanPt->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sPositronMeanNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    canvasdEdxWidthStudies->SetLogx(kFALSE);
    for (Int_t r = 0; r < 4; r++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        legendWidthPt->Clear();
        for (Int_t i = 0; i < nPBins; i++){
            TH1D* histoProjection       = fhistoWidthElectron[r]->ProjectionY(Form("%sPt%d",fhistoWidthElectron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
            legendWidthPt->AddEntry(histoProjection,Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),"p");
        }
        legendWidthPt->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sElectronWidthNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        for (Int_t i = 0; i < nPBins; i++){
            TH1D* histoProjection       = fhistoWidthPositron[r]->ProjectionY(Form("%sPt%d",fhistoWidthPositron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection,  markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
        }
        legendWidthPt->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sPositronWidthNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    //   ____         _____ _______
    //  |  _ \       |  ___|__   __|  /\
    //  | |_| |      | |_     | |    /  \
    //  |  __/       |  _|    | |   / /\ \
    //  | |          | |___   | |  / ___- \
    //  |_|          |_____|  |_| /_/    \_\

    for (Int_t i = 0; i< nPBins; i++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[r]->ProjectionY(Form("%sR%d",fhistoMeanElectron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[r]->ProjectionY(Form("%sR%d",fhistoMeanPositron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]), fCollisionSystem, "", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sMeanNSigma_SummaryEtaDepInRBinsPt%d.%s", outputDir.Data(), i, suffix.Data()));
    }

    for (Int_t i = 0; i< nPBins; i++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[r]->ProjectionY(Form("%sR%d",fhistoWidthElectron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[r]->ProjectionY(Form("%sR%d",fhistoWidthPositron[r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]), fCollisionSystem, "", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sWidthNSigma_SummaryEtaDepInRBinsPt%d.%s", outputDir.Data(), i, suffix.Data()));
    }

    Int_t fColumnPlot = 5;
    Int_t fRowPlot    = 4;
    TString nameCanvasElec = "";
    TString namePadElec    = "";
    TString nameCanvasPosit  = "";
    TString namePadPosit     = "";
    for(Int_t i=0;i<nPBins;i++){
        for (Int_t r = 0; r< 4; r++){
            // Plotting Electron
            nameCanvasElec = Form("Rbin%d_ElecNSigmaDeDx_Pbin%d",r,i);
            namePadElec    = Form("Rbin%d_ElecNSigmaDeDx_Pbin%d",r,i);
            PlotdEdxSlices( histoElectronDeDx[r][i],fitElectronDeDx[r][i],meanElectron[r][i],Form("%sDetailedFits/MapsElec_R%d_P%d.%s",outputDir.Data(),r,i,suffix.Data()),
                            nameCanvasElec,namePadElec,arrEtaBinningOut,optionMC.Data(),fCollisionSystem,nEtaBins,fColumnPlot,fRowPlot,
                            Form("%2.2f < #it{p}_{e^{-}} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),
                            Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1])
                          );

            // Ploting Positron
            nameCanvasPosit  = Form("Rbin%d_PositNSigmaDeDx_Pbin%d",r,i);
            namePadPosit     = Form("Rbin%d_PositNSigmaDeDx_Pbin%d",r,i);
            PlotdEdxSlices( histoPositronDeDx[r][i],fitElectronDeDx[r][i],meanPositron[r][i],Form("%sDetailedFits/MapsPos_R%d_P%d.%s",outputDir.Data(),r,i,suffix.Data()),
                            nameCanvasPosit,namePadPosit,arrEtaBinningOut,optionMC.Data(),fCollisionSystem,nEtaBins,fColumnPlot,fRowPlot,
                            Form("%2.2f < #it{p}_{e^{+}} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),
                            Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1])
                          );
        }
    }

    TFile outFileMonitoring(Form("%sDetailedFits/MonitoringDeDxMaps_%s.root",outputDir.Data(),fCutSelectionRead.Data()) ,"RECREATE");
        for(Int_t i=0;i<nPBins;i++){
            for(Int_t j=0;j<nEtaBins;j++){
                for (Int_t r = 0; r< 4; r++){
                    histoElectronDeDx[r][i][j]->Write();
                    histoPositronDeDx[r][i][j]->Write();
                }
            }
        }
    outFileMonitoring.Close();

    TFile outFileMaps(Form("%sDetailedFits/DeDxMaps_%s.root",outputDir.Data(),fCutSelectionRead.Data()) ,"RECREATE");
    for (Int_t r = 0; r< 4; r++){
        fhistoMeanPositron[r]->Write(Form("Pos_R%d_mean",r));
        fhistoMeanElectron[r]->Write(Form("Ele_R%d_mean",r));
        fhistoWidthPositron[r]->Write(Form("Pos_R%d_width",r));
        fhistoWidthElectron[r]->Write(Form("Ele_R%d_width",r));
    }
    outFileMaps.Close();
 }

//___________________________________________________
TF1* FitSignal(TH1D* hSig, Color_t col){
    TF1* fGaus = new TF1(Form("%sf1",hSig->GetName()),"gaus",-4,5);
    fGaus->SetParameters(hSig->GetMaximum(),0.,1);
    //  fGaus = new TF1("Gaussian","[0]*(exp(-0.5*((x-[1])/[2])^2)",-4,5);

    fGaus->SetLineColor(col);
    fGaus->SetLineWidth(1.);
    if(hSig->GetEntries()<40){
        Mean  = 0.;
        Width = 1.;
        fGaus = NULL;
    } else {

        if (hSig->GetEntries() >=40 && hSig->GetEntries()<150) hSig->Rebin(4);
        hSig->Fit(fGaus,"0RMEQF","",-3.5,+3.5);

        Double_t xMean = fGaus->GetParameter(1);
        fGaus->SetLineColor(kRed);
        hSig->Fit(fGaus,"0RMEQF","",xMean-2.,xMean+2.);
        Mean  = fGaus->GetParameter(1);
        Width = fGaus->GetParameter(2);
        Chi2 = fGaus->GetChisquare()/fGaus->GetNDF();
        if (Width > 2){
            cout << Mean << " " << Width << "  " << hSig->GetEntries() << " " << hSig->GetMaximum()<< endl;
        }
    }
    return fGaus;
}


//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotdEdxSlices(TH1D** fHistoMapsPlot,
                    TF1** fFitMapsPlot,
                    Double_t* fHistoMeanPlot,
                    TString namePlot,
                    TString nameCanvas,
                    TString namePad,
                    Double_t* arrEtaBinnin,
                    TString fMonteCarloInfo,
                    TString fEnergy,
                    Int_t nEtaBins,
                    Int_t fColumnPlot,
                    Int_t fRowPlot,
                    TString ptLegend    = "",
                    TString rLegend     = ""
){

    TGaxis::SetMaxDigits(3);
    TCanvas *canvasDataSpectra          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
    DrawGammaCanvasSettings( canvasDataSpectra, 0., 0.04, 0.04, 0.);
    canvasDataSpectra->cd();

    TPad * padDataSpectra               = new TPad(namePad.Data(),"",0,0,1.,1.,0);   // gives the size of the histo areas
    DrawGammaPadSettings( padDataSpectra, 0, 0, 0, 0);
    padDataSpectra->Divide(fColumnPlot,fRowPlot,0.0,0.0);
    padDataSpectra->Draw();


    Int_t place          = 0;
    Int_t legendPlace[2] = {fColumnPlot, fColumnPlot};
    if (fColumnPlot > 7) legendPlace[0]              = fColumnPlot-1;
    for(Int_t j=0;j<nEtaBins;j++){
        //         cout<<"Pt: "<<j<<" of "<<fNumberPtBins<<endl;
        Double_t startEta            = arrEtaBinnin[j];
        Double_t endEta              = arrEtaBinnin[j+1];

        place                       = place + 1; //give the right place in the page
        if ( place> legendPlace[0]-1 && place < legendPlace[1]+1 ){
            j--;
        } else {
            padDataSpectra->cd(place);
            padDataSpectra->cd(place)->SetTopMargin(0.12);
            padDataSpectra->cd(place)->SetBottomMargin(0.12);
            padDataSpectra->cd(place)->SetRightMargin(0.005);
            int remaining           = (int)((place-1)%fColumnPlot);
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.11);

            TString title = Form("%3.2f < | #eta | < %3.2f",startEta,endEta);
            DrawGammaHisto( fHistoMapsPlot[j],title,"n#sigma","Counts",-6.,6.,0.5,kBlack);
            DrawGammaLines(0.,0.,0.,fHistoMapsPlot[j]->GetMaximum(), 1, kGray+2, 7);
            DrawGammaLines(fHistoMeanPlot[j], fHistoMeanPlot[j],0.,fHistoMapsPlot[j]->GetMaximum(), 1, kBlue+2, 2);
            TLatex *meantext           = new TLatex(0.2, 0.75, Form("#chi^{2}/ndf = %2.2f",fFitMapsPlot[j]->GetChisquare()/fFitMapsPlot[j]->GetNDF()));
            SetStyleTLatex( meantext, 0.075, 1, 1, 42, kTRUE, 11);
            meantext->Draw();

            if (fFitMapsPlot[j]){
                fFitMapsPlot[j]->SetLineColor(kCyan+3);
                fFitMapsPlot[j]->SetLineWidth(2);
                fFitMapsPlot[j]->DrawCopy("same");
                TLatex *chi2text           = new TLatex(0.2, 0.75, Form("#chi^{2}/ndf = %2.2f",fFitMapsPlot[j]->GetChisquare()/fFitMapsPlot[j]->GetNDF()));
                SetStyleTLatex( chi2text, 0.075, 1, 1, 42, kTRUE, 11);
                chi2text->Draw();
            }
        }
    }

    canvasDataSpectra->cd();
    Double_t nPixels        = 13;
    Double_t textHeight     = 0.095;
    Double_t startTextX     = 0.10;
    Int_t columnsLegend     = 1;
    Double_t widthLegend    = 1./fColumnPlot;
    Double_t heightLegend   = 1./fRowPlot;
    Double_t marginWidthLeg = 0.15;
    Int_t exampleBin        = fColumnPlot-1;
    if (fColumnPlot > 7){
        startTextX          = 0.05;
        nPixels             = 12;
        widthLegend         = 2./fColumnPlot;
        marginWidthLeg      = 0.25;
    }

    // plotting Legend
    TPad * padLegend                = new TPad("dummyPad","",1-widthLegend,1-heightLegend,1.,1.,0);   // gives the size of the histo areas
    DrawGammaPadSettings( padLegend, 0, 0, 0, 0);
    padLegend->Draw();
    padLegend->cd();

    TString textAlice       = "ALICE performance";
    TString textEvents;
    if(fMonteCarloInfo.CompareTo("MC")==0){
        textEvents          = "MC";
    } else {
        textEvents          = "Data";
    }

    if (padLegend->XtoPixel(padLegend->GetX2()) < padLegend->YtoPixel(padLegend->GetY1())){
        textHeight          = (Double_t)nPixels/padLegend->XtoPixel(padLegend->GetX2()) ;
    } else {
        textHeight          = (Double_t)nPixels/padLegend->YtoPixel(padLegend->GetY1());
    }
    Double_t startTextY     = 0.9;
    Double_t differenceText = textHeight*1.05;
    // plot labels
    PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, ptLegend, rLegend, fEnergy, "", "");

    TLine *IdealMean = new TLine (0.,0.,0.,0.);
    IdealMean->SetLineColor(kGray+2);
    IdealMean->SetLineStyle(7);

    TLine *l1 = new TLine (0.,0.,0.,0.);
    l1->SetLineColor(kBlue+2);
    l1->SetLineStyle(2);

    TLegend* legendData     = GetAndSetLegend2(  startTextX, startTextY-5.75*differenceText, 0.85,  startTextY-(5.75+3/columnsLegend)*differenceText, nPixels, columnsLegend, "", 43, marginWidthLeg);
    Size_t markersize       = fHistoMapsPlot[exampleBin]->GetMarkerSize();
//     fHistoMapsPlot[exampleBin]->SetMarkerSize(3*markersize);
    legendData->AddEntry(fHistoMapsPlot[exampleBin],"n#sigma projection","ep");
    legendData->AddEntry(IdealMean,"ideal mean position","l");
    legendData->AddEntry(l1,"actual mean position","l");
    legendData->Draw();

    canvasDataSpectra->Print(namePlot.Data());
    delete padLegend;
    delete padDataSpectra;
    delete canvasDataSpectra;
}

void DrawSliceHisto(TH1* histo1,
                    TString Title,
                    TString XTitle,
                    TString YTitle,
                    Float_t xMin,
                    Float_t xMax,
                    Size_t  markerSize,
                    Color_t markerColor
){

    Double_t yMin = 0;
    Double_t yMax = 0;
    for (Int_t i = histo1->GetXaxis()->FindBin(xMin); i < histo1->GetXaxis()->FindBin(xMax); i++){
        if (histo1->GetBinContent(i) < yMin){
            yMin = histo1->GetBinContent(i);
        }
        if (histo1->GetBinContent(i) > yMax){
            yMax = histo1->GetBinContent(i);
        }
    }

    if (xMin > 0.2)  histo1->GetYaxis()->SetRangeUser(yMin, 1.5*yMax);
    else             histo1->GetYaxis()->SetRangeUser(yMin, 1.2*yMax);

    if(XTitle.Length() > 0){
        histo1->SetXTitle(XTitle.Data());
    }
    if(YTitle.Length() > 0){
        histo1->SetYTitle(YTitle.Data());
    }
    histo1->GetXaxis()->SetLabelSize(0.06);
    histo1->GetXaxis()->SetTitleSize(0.1);
    histo1->GetXaxis()->SetTitleOffset(-0.5);
    histo1->GetYaxis()->SetLabelSize(0.06);
    histo1->GetYaxis()->SetTitleSize(0.1);
    histo1->GetYaxis()->SetTitleOffset(-0.5);
    histo1->GetYaxis()->SetDecimals();
    histo1->GetXaxis()->SetNdivisions(507,kTRUE);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerColor(1);
    histo1->SetLineColor(1);
    histo1->SetLineWidth(1);
    histo1->SetMarkerSize(markerSize);
    histo1->SetTitleOffset(2.2,"xy");
    histo1->SetTitleSize(0.05,"xy");
    if(Title.Length() > 0){
        histo1->SetTitle("");
    }
    if(Title.Length() > 0){
        TLatex *alice = new TLatex(0.4,0.92,Form("%s",Title.Data()));
        alice->SetNDC();
        alice->SetTextColor(1);
        alice->SetTextSize(0.08);
        alice->Draw();
    }

    histo1->SetTitle(Title);
    histo1->Draw("same,e1,p");

}