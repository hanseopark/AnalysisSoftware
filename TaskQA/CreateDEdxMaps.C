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
Double_t meanElectron[2][4][14][20];//type r pt eta
Double_t meanPositron[2][4][14][20];//type r pt eta
Double_t meanErrElectron[2][4][14][20];//type r pt eta
Double_t meanErrPositron[2][4][14][20];//type r pt eta
Double_t widthElectron[2][4][14][20]; //type r pt eta
Double_t widthPositron[2][4][14][20]; //type r pt eta
Double_t widthErrElectron[2][4][14][20]; //type r pt eta
Double_t widthErrPositron[2][4][14][20]; //type r pt eta
Double_t chi[20][20];
Bool_t doFitting[2][4][14][20]; //type r pt eta
Double_t scalWidthFit[2][4][14][20]; //type r pt eta

void CreateDEdxMaps(    TString fileNameWithMaps    ="" ,
                        TString cutSelection        = "",
                        TString optionMC            = "Data",
                        TString fEnergy             = "",
                        TString suffix              = "pdf",
                        TString period              = "",
                        Bool_t isLowStat            = kFALSE
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
    gSystem->Exec(Form("mkdir -p %s/%s/%s",outputDir.Data(),"RConv","DetailedFits"));
    gSystem->Exec(Form("mkdir -p %s/%s/%s",outputDir.Data(),"TPCCl","DetailedFits"));

    TH3F *histoPositronDeDxPEta[2][4] = {NULL};
    TH3F *histoElectronDeDxPEta[2][4] = {NULL};

    TFile* fileMaterialHistos = new TFile(fileNameWithMaps.Data());
    TString fCutSelectionRead = cutSelection;
    TString nameMainDir = "GammaConvMaterial";

    TString fCollisionSystem = ReturnFullCollisionsSystem(fEnergy);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    TString centralityString    = GetCentralityString(cutSelection);
    if (centralityString.CompareTo("pp")!=0 && !centralityString.Contains("0-100%") ){
        fCollisionSystem    = Form("%s %s", centralityString.Data(), fCollisionSystem.Data());
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
            histoPositronDeDxPEta[0][i] = (TH3F*)MapsContainer->FindObject(Form("R%d positron sigma dEdx P Eta",i));
            histoElectronDeDxPEta[0][i] = (TH3F*)MapsContainer->FindObject(Form("R%d electron sigma dEdx P Eta",i));
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
        for (Int_t r = 0; r<4; r++){
            histoPositronDeDxPEta[0][r] = (TH3F*)directoryConv->Get(Form("R%d positron sigma dEdx P Eta",r));
            histoElectronDeDxPEta[0][r] = (TH3F*)directoryConv->Get(Form("R%d electron sigma dEdx P Eta",r));
            histoPositronDeDxPEta[1][r] = (TH3F*)directoryConv->Get(Form("Cl%d positron sigma dEdx P Eta",r));
            histoElectronDeDxPEta[1][r] = (TH3F*)directoryConv->Get(Form("Cl%d electron sigma dEdx P Eta",r));
        }
    }

    Color_t colorR[4]           = { kBlack, kBlue+2, kRed+2, kGreen+2};
    Color_t colorEta[18]        = { kBlack, kGray+2, kBlue-7, kBlue-4, kBlue+2, kCyan+2, kCyan-5, kGreen-7, kGreen+2, kYellow+1,
                                    kOrange-5, kOrange-3, kRed+2, kRed-4, kRed-6, kMagenta-5, kMagenta-2, kViolet-7};
    Color_t colorPt[14]         = { kBlack, kGray+2, kBlue-7, kBlue+2, kCyan+2, kCyan-5, kGreen-7, kGreen+2, kOrange-3, kRed+2,
                                    kRed-6, kMagenta-5, kMagenta-2, kViolet-7};
    Marker_t markerR[8]         = { 20, 21, 33, 34, 24, 25, 27, 28};
    Marker_t markerEta[18]      = { 24, 25, 27, 28, 30, 29, 34, 33, 21, 20,
                                    42, 46, 40, 41, 27, 43, 24, 25};
    Marker_t markerPt[14]       = { 24, 25, 27, 28, 30, 29, 34, 33, 21, 20,
                                    42, 46, 40, 41};

    Int_t nPBins   = 14;
    Double_t arrPBinning[15]    = { 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                                    2.0, 5.0, 10.0};
    if (isLowStat){
        arrPBinning[6]  = 0.6;
        arrPBinning[7]  = 0.8;
        arrPBinning[8]  = 1.0;
        arrPBinning[9]  = 2.0;
        arrPBinning[10] = 10.0;
        nPBins          = 10;
    }

    Int_t startBinP = 1;
    if (fEnergy.Contains("LowB") ) startBinP = 0;
    for( Int_t i=startBinP;i<nPBins+1;i++){
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

    Int_t nBinsTPCCl        = 4;
    Int_t arrTPCClBinningOut[5] = {0, 60, 100, 150, 180};
    for (Int_t i=0;i<nBinsTPCCl+1;i++){
        cout << "TPC Cl bins:: " << i << " " << arrTPCClBinningOut[i] << endl;
    }

    // fitting settings
    Double_t scalWidthFitStd    = 3;
    Double_t maxChi2            = 200;     // maximum allowed Chi2 of fit
    Bool_t fixSearchRange       = kFALSE;  // in case there is another, higher peak in the nSigma histos, should fix the range to search for bin with maximum content before fitting
    Double_t fitRange           = 3;       // -fitRange, fitRange
    Int_t firstRebinBoundary    = 200;     // If maximum counts is smaller than secRebinBoundary, do first rebin with factor 2
    Int_t secRebinBoundary      = 100;     // If maximum counts is smaller than secRebinBoundary, do second rebin with factor 2

    if(fEnergy.Contains("PbPb")){
        scalWidthFitStd    = 1.5;
        maxChi2            = 20;
        fixSearchRange     = kTRUE;
        firstRebinBoundary = 100;
        secRebinBoundary   = 0;

    } else if ((fEnergy.Contains("XeXe")) && period.Contains("Peri")){
        scalWidthFitStd    = 2;
        maxChi2            = 20;
    } else if (fEnergy.Contains("XeXe") && period.Contains("Cent")){
        scalWidthFitStd    = 1.5;
        maxChi2            = 20;
    }

    // enable fitting
    for (Int_t r = 0; r < 4; r++){
        for (Int_t j = 0; j < nPBins; j++){
            for (Int_t i = 0; i < nEtaBinsOut; i++){
                doFitting[0][r][j][i] = kTRUE;
                doFitting[1][r][j][i] = kTRUE;
                scalWidthFit[0][r][j][i] = scalWidthFitStd;
                scalWidthFit[1][r][j][i] = scalWidthFitStd;
            }
        }
    }


    if (fEnergy.Contains("XeXe")) {
        for (Int_t i = 0; i < nEtaBinsOut; i++){
            doFitting[1][0][3][i] = kFALSE;
            doFitting[1][0][4][i] = kFALSE;
        }
    } else if(fEnergy.Contains("PbPb")){

        for (Int_t i = 0; i < nEtaBinsOut; i++){

            //  [0][R][Pt][Eta]

            // R00
            doFitting[0][0][1][i]  = kFALSE;
            doFitting[0][0][9][i] = kFALSE;  // list2
            doFitting[0][0][10][i] = kFALSE;
            doFitting[0][0][11][i] = kFALSE;
            doFitting[0][0][12][i] = kFALSE;
            doFitting[0][0][13][i] = kFALSE;
            scalWidthFit[0][0][5][i]= 2.0;
            scalWidthFit[0][0][6][i]= 2.0;
            scalWidthFit[0][0][7][i]= 2.0;
            scalWidthFit[0][0][8][i]= 2.0;  // list2
            scalWidthFit[0][0][9][i]= 2.0;

             // R01
            //if( i < 2 || i > 15 ) doFitting[0][1][1][i] = kFALSE;   // list1
            if( i < 9 || i > 15 ) doFitting[0][1][1][i] = kFALSE;   // list2
            for(Int_t p = startBinP; p < nPBins; p++){
                if(p==1 || p==9) scalWidthFit[0][1][p][i] = 2.5;
                else if(p==6) scalWidthFit[0][1][p][i] = 3.0;
                else scalWidthFit[0][1][p][i] = 2.0;
            }
            doFitting[0][1][9][i] = kFALSE;  // list2
            doFitting[0][1][10][i] = kFALSE;  // list2
            doFitting[0][1][11][i] = kFALSE;  // list2
            doFitting[0][1][12][i] = kFALSE;
            doFitting[0][1][13][i] = kFALSE;

            // R02
            doFitting[0][2][9][i] = kFALSE;   // list2
            doFitting[0][2][10][i] = kFALSE;  // list2
            doFitting[0][2][11][i] = kFALSE;  // list2
            doFitting[0][2][12][i] = kFALSE;
            doFitting[0][2][13][i] = kFALSE;
            for(Int_t p = startBinP; p < nPBins; p++){
                if(p==10 || p==12) scalWidthFit[0][2][p][i] = 2.0;
                else scalWidthFit[0][2][p][i] = 3.0;
            }

            // R03
            if( i < 2 || i > 16 ) doFitting[0][3][1][i] = kFALSE;
            if( i < 2 || i > 16 ) doFitting[0][3][3][i] = kFALSE;  // list2
            for(Int_t p = 5; p < 8; p++){
                if( i < 1 || i > 16 ) doFitting[0][3][p][i] = kFALSE;
            }
            //for(Int_t p = 8; p < nPBins; p++){      // list1
            for(Int_t p = 3; p < nPBins; p++){        // list2
                doFitting[0][3][p][i] = kFALSE;
            }
            for(Int_t p = 0; p < nPBins; p++){
                scalWidthFit[0][3][p][i] = 2.0;
            }

            //  [1][TPCCl][Pt][Eta]

            // Cl00
            scalWidthFit[1][0][1][i] = 2.0;
            scalWidthFit[1][0][2][i] = 2.0;
            //for(Int_t p = 5; p < nPBins; p++){   // list1
            for(Int_t p = 3; p < nPBins; p++){     // list2
                doFitting[1][0][p][i]  = kFALSE;
            }

            // Cl01
            for(Int_t p = 8; p < nPBins; p++){
                doFitting[1][1][p][i]  = kFALSE;
            }
            if( i ==0 || i == 17 ) doFitting[1][1][1][i]  = kFALSE;
            scalWidthFit[1][1][1][i] = 2.5;
            scalWidthFit[1][1][2][i] = 2.5;
            scalWidthFit[1][1][5][i] = 2.2; // list1 2.0;
            scalWidthFit[1][1][6][i] = 2.2; // list1 2.0;
            scalWidthFit[1][1][7][i] = 2.2; // list1 2.0;

            // Cl02
            doFitting[1][2][1][i]  = kFALSE;
            doFitting[1][2][12][i] = kFALSE;
            doFitting[1][2][13][i] = kFALSE;
            scalWidthFit[1][2][2][i] = 2.0;
            for(Int_t p = 6; p < 11; p++){
                scalWidthFit[1][2][p][i] = 2.0;
            }

            // Cl03
            doFitting[1][3][1][i]  = kFALSE;
            if( i < 6 || i > 11 ) doFitting[1][3][2][i]  = kFALSE;
            if( i < 2 || i > 16 ) doFitting[1][3][3][i]  = kFALSE;
            for(Int_t p = startBinP; p < 10; p++){
                if((p>3) &&( i < 1 || i > 16) ) doFitting[1][3][p][i]  = kFALSE;
                else scalWidthFit[1][3][p][i] = 2.5;
            }
            doFitting[1][3][10][i] = kFALSE;
            doFitting[1][3][11][i] = kFALSE;
            doFitting[1][3][12][i] = kFALSE;
            doFitting[1][3][13][i] = kFALSE;
        }
    }

    TH2F* fhistoMeanPositron[2][4]     = {{NULL}};
    TH2F* fhistoMeanElectron[2][4]     = {{NULL}};
    TH2F* fhistoWidthElectron[2][4]    = {{NULL}};
    TH2F* fhistoWidthPositron[2][4]    = {{NULL}};

    TH2S* fhistoRecalibPositron[2][4]     = {{NULL}};
    TH2S* fhistoRecalibElectron[2][4]     = {{NULL}};

    for (Int_t r = 0; r<4; r++){
        fhistoMeanPositron[0][r]       = new TH2F(Form("MeanPosiR%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoMeanPositron[0][r]->Sumw2();
        fhistoMeanElectron[0][r]       = new TH2F(Form("MeanElecR%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoMeanElectron[0][r]->Sumw2();
        fhistoWidthPositron[0][r]      = new TH2F(Form("WidthPosiR%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoWidthPositron[0][r]->Sumw2();
        fhistoWidthElectron[0][r]      = new TH2F(Form("WidthElecR%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoWidthElectron[0][r]->Sumw2();
        fhistoMeanPositron[1][r]       = new TH2F(Form("MeanPosiCl%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoMeanPositron[1][r]->Sumw2();
        fhistoMeanElectron[1][r]       = new TH2F(Form("MeanElecCl%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoMeanElectron[1][r]->Sumw2();
        fhistoWidthPositron[1][r]      = new TH2F(Form("WidthPosiCl%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoWidthPositron[1][r]->Sumw2();
        fhistoWidthElectron[1][r]      = new TH2F(Form("WidthElecCl%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoWidthElectron[1][r]->Sumw2();
        fhistoRecalibPositron[0][r]       = new TH2S(Form("RecalibPosiR%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoRecalibPositron[0][r]->Sumw2();
        fhistoRecalibPositron[1][r]       = new TH2S(Form("RecalibPosiCl%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoRecalibPositron[1][r]->Sumw2();
        fhistoRecalibElectron[0][r]       = new TH2S(Form("RecalibElecR%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoRecalibElectron[0][r]->Sumw2();
        fhistoRecalibElectron[1][r]       = new TH2S(Form("RecalibElecCl%d",r),"",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
        fhistoRecalibElectron[1][r]->Sumw2();


    }


    TH1D* histoPositronDeDx[2][4][14][20]  = {{{{NULL}}}};
    TH1D* histoElectronDeDx[2][4][14][20]  = {{{{NULL}}}};
    TF1* fitPositronDeDx[2][4][14][20]     = {{{{NULL}}}};
    TF1* fitElectronDeDx[2][4][14][20]     = {{{{NULL}}}};

    for (Int_t pt = 1; pt < histoElectronDeDxPEta[0][0]->GetNbinsZ()+1; pt++){
        cout << histoElectronDeDxPEta[0][0]->GetZaxis()->GetBinCenter(pt) << "; "  ;
    }
    cout << endl;

    for(Int_t i=startBinP;i<nPBins;i++){
        Int_t minBinPt    = histoElectronDeDxPEta[0][0]->GetZaxis()->FindBin(arrPBinning[i]);
        Int_t maxBinPt    = histoElectronDeDxPEta[0][0]->GetZaxis()->FindBin(arrPBinning[i+1]-0.001);
        cout << "*****************\np = " << arrPBinning[i] << " GeV/c" << "< pT < " << arrPBinning[i+1] << " GeV/c"<< endl;
        cout << minBinPt << "\t" << maxBinPt << endl;
    }

    Int_t etaOff = 0;
//     if(fEnergy.CompareTo("5TeV2017")==0) etaOff = 0;

    //    cout<< histoPositronDeDxPEta0->GetNbinsY()<< " " << histoPositronDeDxPEta0->GetNbinsZ() << endl;
    for(Int_t i=startBinP;i<nPBins;i++){

        Int_t minBinPt    = histoElectronDeDxPEta[0][0]->GetZaxis()->FindBin(arrPBinning[i]);
        Int_t maxBinPt    = histoElectronDeDxPEta[0][0]->GetZaxis()->FindBin(arrPBinning[i+1]-0.001);
        cout << "*****************\np = " << arrPBinning[i] << " GeV/c" << "< pT < " << arrPBinning[i+1] << " GeV/c"<< endl;
        cout << minBinPt << "\t" << maxBinPt << endl;
        for(Int_t j=0;j<nEtaBins;j++){
//             cout << "-> eta " << arrEtaBinningOut[j] << endl;
//             cout << "offset " << arrEtaBinningOut[j-etaOff] << endl;
//             cout << " j-etaOff " << j-etaOff << endl;
            Int_t minBinEta    = histoElectronDeDxPEta[0][0]->GetYaxis()->FindBin(arrEtaBinningOut[j]);
            Int_t maxBinEta    = histoElectronDeDxPEta[0][0]->GetYaxis()->FindBin(arrEtaBinningOut[j+1]-0.001);

            for (Int_t r = 0; r < 4; r++){
                // initialize default values
                meanElectron[0][r][i][j]        = 0;
                widthElectron[0][r][i][j]       = 1;
                meanPositron[0][r][i][j]        = 0;
                widthPositron[0][r][i][j]       = 1;
                meanErrElectron[0][r][i][j]     = 0.2;
                widthErrElectron[0][r][i][j]    = 0.1;
                meanErrPositron[0][r][i][j]     = 0.2;
                widthErrPositron[0][r][i][j]    = 0.1;
                meanElectron[1][r][i][j]        = 0;
                widthElectron[1][r][i][j]       = 1;
                meanPositron[1][r][i][j]        = 0;
                widthPositron[1][r][i][j]       = 1;
                meanErrElectron[1][r][i][j]     = 0.2;
                widthErrElectron[1][r][i][j]    = 0.1;
                meanErrPositron[1][r][i][j]     = 0.2;
                widthErrPositron[1][r][i][j]    = 0.1;

                cout << r << "\t" << i << "\t" << j << "\t" << doFitting[0][r][i][j] << "\t"<< doFitting[1][r][i][j] << endl;
                // Electron, different R bins
                histoElectronDeDx[0][r][i][j]  = (TH1D*)histoElectronDeDxPEta[0][r]->ProjectionX(Form("R%dElecSigdEdxP%dEta%d",r,i,j),minBinEta,maxBinEta,minBinPt,maxBinPt);
                if (histoElectronDeDx[0][r][i][j]->GetMaximum() < firstRebinBoundary) histoElectronDeDx[0][r][i][j]->Rebin(2);
                if (histoElectronDeDx[0][r][i][j]->GetMaximum() < secRebinBoundary) histoElectronDeDx[0][r][i][j]->Rebin(2);
//                 fitElectronDeDx[0][r][i][j]    = FitSignal(histoElectronDeDx[0][r][i][j],kBlue);
                if (doFitting[0][r][i][j] && histoElectronDeDx[0][r][i][j]->GetMaximum() > 10)
                    fitElectronDeDx[0][r][i][j]   = FitTH1DRecursivelyGaussianWExp(histoElectronDeDx[0][r][i][j], 0.02, -fitRange, fitRange, scalWidthFit[0][r][i][j], maxChi2, fixSearchRange);
//                     fitElectronDeDx[0][r][i][j]   = FitTH1DRecursivelyGaussian (histoElectronDeDx[0][r][i][j], 0.02, -3, 3, 2, 1.25);

                if ( fitElectronDeDx[0][r][i][j]){
                    meanElectron[0][r][i][j]       = fitElectronDeDx[0][r][i][j]->GetParameter(1);
                    widthElectron[0][r][i][j]      = fitElectronDeDx[0][r][i][j]->GetParameter(2);
                    meanErrElectron[0][r][i][j]    = fitElectronDeDx[0][r][i][j]->GetParError(1);
                    widthErrElectron[0][r][i][j]   = fitElectronDeDx[0][r][i][j]->GetParError(2);
                } else if (doFitting[0][r][i][j] && histoElectronDeDx[0][r][i][j]->GetMaximum() > 5){
                    meanElectron[0][r][i][j]       = histoElectronDeDx[0][r][i][j]->GetMean();
                    widthElectron[0][r][i][j]      = histoElectronDeDx[0][r][i][j]->GetRMS();
                }
                if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                    fhistoMeanElectron[0][r]->SetBinContent(i+1,j+1-etaOff, meanElectron[0][r][i][j]);
                    fhistoMeanElectron[0][r]->SetBinError(i+1,j+1-etaOff, meanErrElectron[0][r][i][j]);
                    fhistoWidthElectron[0][r]->SetBinContent(i+1,j+1-etaOff, widthElectron[0][r][i][j]);
                    fhistoWidthElectron[0][r]->SetBinError(i+1,j+1-etaOff, widthErrElectron[0][r][i][j]);
                    fhistoRecalibElectron[0][r]->SetBinContent(i+1,j+1-etaOff, (Short_t)(meanElectron[0][r][i][j]*1000));
                    fhistoRecalibElectron[0][r]->SetBinError(i+1,j+1-etaOff, (Short_t)(widthElectron[0][r][i][j]*1000));
                }

                // Positron, different R bins
                histoPositronDeDx[0][r][i][j] = (TH1D*)histoPositronDeDxPEta[0][r]->ProjectionX(Form("R%dPosiSigdEdxP%dEta%d",r,i,j),minBinEta,maxBinEta,minBinPt,maxBinPt);
                if (histoPositronDeDx[0][r][i][j]->GetMaximum() < firstRebinBoundary) histoPositronDeDx[0][r][i][j]->Rebin(2);
                if (histoPositronDeDx[0][r][i][j]->GetMaximum() < secRebinBoundary) histoPositronDeDx[0][r][i][j]->Rebin(2);
//                 fitPositronDeDx[0][r][i][j]    = FitSignal(histoPositronDeDx[0][r][i][j],kBlue);
                if (doFitting[0][r][i][j] &&histoPositronDeDx[0][r][i][j]->GetMaximum() > 10)
                    fitPositronDeDx[0][r][i][j]   = FitTH1DRecursivelyGaussianWExp(histoPositronDeDx[0][r][i][j], 0.02, -fitRange, fitRange, scalWidthFit[0][r][i][j], maxChi2, fixSearchRange);
//                     fitPositronDeDx[0][r][i][j]   = FitTH1DRecursivelyGaussian (histoPositronDeDx[0][r][i][j], 0.02, -3, 3, 2, 1.25 );

                if (fitPositronDeDx[0][r][i][j]){
                    meanPositron[0][r][i][j]       = fitPositronDeDx[0][r][i][j]->GetParameter(1);
                    widthPositron[0][r][i][j]      = fitPositronDeDx[0][r][i][j]->GetParameter(2);
                    meanErrPositron[0][r][i][j]    = fitPositronDeDx[0][r][i][j]->GetParError(1);
                    widthErrPositron[0][r][i][j]   = fitPositronDeDx[0][r][i][j]->GetParError(2);
                } else if (doFitting[0][r][i][j] && histoPositronDeDx[0][r][i][j]->GetMaximum() > 5){
                    meanPositron[0][r][i][j]        = histoPositronDeDx[0][r][i][j]->GetMean();
                    widthPositron[0][r][i][j]       = histoPositronDeDx[0][r][i][j]->GetRMS();
                }
                if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                    fhistoMeanPositron[0][r]->SetBinContent(i+1,j+1-etaOff, meanPositron[0][r][i][j]);
                    fhistoMeanPositron[0][r]->SetBinError(i+1,j+1-etaOff, meanErrPositron[0][r][i][j]);
                    fhistoWidthPositron[0][r]->SetBinContent(i+1,j+1-etaOff, widthPositron[0][r][i][j]);
                    fhistoWidthPositron[0][r]->SetBinError(i+1,j+1-etaOff, widthErrPositron[0][r][i][j]);
                    fhistoRecalibPositron[0][r]->SetBinContent(i+1,j+1-etaOff, (Short_t)(meanPositron[0][r][i][j]*1000));
                    fhistoRecalibPositron[0][r]->SetBinError(i+1,j+1-etaOff, (Short_t)(widthPositron[0][r][i][j]*1000));
                }


                // Electron, different TPC Cl bins
                histoElectronDeDx[1][r][i][j]  = (TH1D*)histoElectronDeDxPEta[1][r]->ProjectionX(Form("Cl%dElecSigdEdxP%dEta%d",r,i,j),minBinEta,maxBinEta,minBinPt,maxBinPt);
                if (histoElectronDeDx[1][r][i][j]->GetMaximum() < firstRebinBoundary) histoElectronDeDx[1][r][i][j]->Rebin(2);
                if (histoElectronDeDx[1][r][i][j]->GetMaximum() < secRebinBoundary) histoElectronDeDx[1][r][i][j]->Rebin(2);
                //                 fitElectronDeDx[1][r][i][j]    = FitSignal(histoElectronDeDx[1][r][i][j],kBlue);
                if (doFitting[1][r][i][j] && histoElectronDeDx[1][r][i][j]->GetMaximum() > 10)
                    fitElectronDeDx[1][r][i][j]   = FitTH1DRecursivelyGaussianWExp(histoElectronDeDx[1][r][i][j], 0.02, -fitRange, fitRange, scalWidthFit[1][r][i][j], maxChi2, fixSearchRange);
//                     fitElectronDe[1][r][i][j]   = FitTH1DRecursivelyGaussian (histoElectronDeDx[1][r][i][j], 0.02, -3, 3, 2, 1.25);

                if (fitElectronDeDx[1][r][i][j]){
                    meanElectron[1][r][i][j]       = fitElectronDeDx[1][r][i][j]->GetParameter(1);
                    widthElectron[1][r][i][j]      = fitElectronDeDx[1][r][i][j]->GetParameter(2);
                    meanErrElectron[1][r][i][j]    = fitElectronDeDx[1][r][i][j]->GetParError(1);
                    widthErrElectron[1][r][i][j]   = fitElectronDeDx[1][r][i][j]->GetParError(2);
                } else if (doFitting[1][r][i][j] && histoElectronDeDx[1][r][i][j]->GetMaximum() > 5){
                    meanElectron[1][r][i][j]       = histoElectronDeDx[1][r][i][j]->GetMean();
                    widthElectron[1][r][i][j]      = histoElectronDeDx[1][r][i][j]->GetRMS();
                }
                if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                    fhistoMeanElectron[1][r]->SetBinContent(i+1,j+1-etaOff, meanElectron[1][r][i][j]);
                    fhistoMeanElectron[1][r]->SetBinError(i+1,j+1-etaOff, meanErrElectron[1][r][i][j]);
                    fhistoWidthElectron[1][r]->SetBinContent(i+1,j+1-etaOff, widthElectron[1][r][i][j]);
                    fhistoWidthElectron[1][r]->SetBinError(i+1,j+1-etaOff, widthErrElectron[1][r][i][j]);
                    fhistoRecalibElectron[1][r]->SetBinContent(i+1,j+1-etaOff, (Short_t)(meanElectron[1][r][i][j]*1000));
                    fhistoRecalibElectron[1][r]->SetBinError(i+1,j+1-etaOff,(Short_t)( widthElectron[1][r][i][j]*1000));
                }

                // Positron, different TPC Cl bins
                histoPositronDeDx[1][r][i][j] = (TH1D*)histoPositronDeDxPEta[1][r]->ProjectionX(Form("Cl%dPosiSigdEdxP%dEta%d",r,i,j),minBinEta,maxBinEta,minBinPt,maxBinPt);
                if (histoPositronDeDx[1][r][i][j]->GetMaximum() < firstRebinBoundary) histoPositronDeDx[1][r][i][j]->Rebin(2);
                if (histoPositronDeDx[1][r][i][j]->GetMaximum() < secRebinBoundary) histoPositronDeDx[1][r][i][j]->Rebin(2);
                //                 fitPositronDeDx[1][r][i][j]    = FitSignal(histoPositronDeDx[1][r][i][j],kBlue);
                if (doFitting[1][r][i][j] && histoPositronDeDx[1][r][i][j]->GetMaximum() > 10)
                    fitPositronDeDx[1][r][i][j]   = FitTH1DRecursivelyGaussianWExp(histoPositronDeDx[1][r][i][j], 0.02, -fitRange, fitRange, scalWidthFit[1][r][i][j], maxChi2, fixSearchRange);
            //  fitPositronDeDx[1][r][i][j]   = FitTH1DRecursivelyGaussian (histoPositronDeDx[1][r][i][j], 0.02, -3, 3, 2, 1.25 );

                if (fitPositronDeDx[1][r][i][j]){
                    meanPositron[1][r][i][j]        = fitPositronDeDx[1][r][i][j]->GetParameter(1);
                    widthPositron[1][r][i][j]       = fitPositronDeDx[1][r][i][j]->GetParameter(2);
                    meanErrPositron[1][r][i][j]     = fitPositronDeDx[1][r][i][j]->GetParError(1);
                    widthErrPositron[1][r][i][j]    = fitPositronDeDx[1][r][i][j]->GetParError(2);
                } else if (doFitting[1][r][i][j] && histoPositronDeDx[1][r][i][j]->GetMaximum() > 5){
                    meanPositron[1][r][i][j]        = histoPositronDeDx[1][r][i][j]->GetMean();
                    widthPositron[1][r][i][j]       = histoPositronDeDx[1][r][i][j]->GetRMS();
                }
                if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                    fhistoMeanPositron[1][r]->SetBinContent(i+1,j+1-etaOff, meanPositron[1][r][i][j]);
                    fhistoMeanPositron[1][r]->SetBinError(i+1,j+1-etaOff, meanErrPositron[1][r][i][j]);
                    fhistoWidthPositron[1][r]->SetBinContent(i+1,j+1-etaOff, widthPositron[1][r][i][j]);
                    fhistoWidthPositron[1][r]->SetBinError(i+1,j+1-etaOff, widthErrPositron[1][r][i][j]);
                    fhistoRecalibPositron[1][r]->SetBinContent(i+1,j+1-etaOff, (Short_t)(meanPositron[1][r][i][j]*1000));
                    fhistoRecalibPositron[1][r]->SetBinError(i+1,j+1-etaOff, (Short_t)(widthPositron[1][r][i][j]*1000));
                }

            }
        }
    }
    cout << "------------------------------------------------------------------------------------" << endl;
    cout << "------------------------------------------------------------------------------------" << endl;
    cout << "STARTED PLOTTING" << endl;
    cout << "------------------------------------------------------------------------------------" << endl;
    cout << "------------------------------------------------------------------------------------" << endl;

    Int_t textSizeLabelsPixel                   = 1000*0.04;
    TCanvas* canvasdEdxMeanStudies = new TCanvas("canvasdEdxMeanStudies","",200,10,1350,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasdEdxMeanStudies, 0.08, 0.02, 0.02, 0.09);
    canvasdEdxMeanStudies->SetLogx();
    TH2F * histo2DdEdxMeanStudies = new TH2F("histo2DdEdxMeanStudies","histo2DdEdxMeanStudies", 100, 0.01, 20, 200, -1.5, 1.5);
    SetStyleHistoTH2ForGraphs(histo2DdEdxMeanStudies, "#it{p} (GeV/#it{c})", "<n#sigma_{e}>", 0.035,0.04, 0.035,0.04, 0.98,1.0);

    TLegend* legendMeanEta  = GetAndSetLegend2(0.3, 0.11, 0.95, 0.11+(0.035*nEtaBinsOut/4.*1.05), 0.75*textSizeLabelsPixel,4,"",43,0.1);
    TLegend* legendMeanR    = GetAndSetLegend2(0.27, 0.11, 0.95, 0.11+(0.035*4.*1.05), 0.75*textSizeLabelsPixel,2,"",43,0.1);
    TCanvas* canvasdEdxWidthStudies = new TCanvas("canvasdEdxWidthStudies","",200,10,1350,1000);  // gives the page size
    DrawGammaCanvasSettings( canvasdEdxWidthStudies, 0.08, 0.02, 0.02, 0.09);
    canvasdEdxWidthStudies->SetLogx();
    TH2F * histo2DdEdxWidthStudies = new TH2F("histo2DdEdxWidthStudies","histo2DdEdxWidthStudies", 100, 0.01, 20, 200, 0, 2);
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
            TH1D* histoProjection       = fhistoMeanElectron[0][r]->ProjectionX(Form("%sEta%d",fhistoMeanElectron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
            legendMeanEta->AddEntry(histoProjection,Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]),"p");
        }
        legendMeanEta->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);

        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/ElectronMeanNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoMeanPositron[0][r]->ProjectionX(Form("%sEta%d",fhistoMeanPositron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
        }
        legendMeanEta->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);

        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/PositronMeanNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }


    for (Int_t r = 0; r < 4; r++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthEta->Clear();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoWidthElectron[0][r]->ProjectionX(Form("%sEta%d",fhistoWidthElectron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
            legendWidthEta->AddEntry(histoProjection,Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]),"p");
        }
        legendWidthEta->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);

        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/ElectronWidthNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoWidthPositron[0][r]->ProjectionX(Form("%sEta%d",fhistoWidthPositron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
        }
        legendWidthEta->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/PositronWidthNSigma_SummaryPtDepInEtaBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    //   _____ _______                 ____ _
    //  |  ___|__   __|  /\           / ___| |
    //  | |_     | |    /  \         | /   | |
    //  |  _|    | |   / /\ \        | |   | |
    //  | |___   | |  / ___- \       | \___| |___
    //  |_____|  |_| /_/    \_\       \____|_____|


    for (Int_t r = 0; r < 4; r++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanEta->Clear();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoMeanElectron[1][r]->ProjectionX(Form("%sEta%d",fhistoMeanElectron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
            legendMeanEta->AddEntry(histoProjection,Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]),"p");
        }
        legendMeanEta->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/ElectronMeanNSigma_SummaryPtDepInEtaBinsTPCCl%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoMeanPositron[1][r]->ProjectionX(Form("%sEta%d",fhistoMeanPositron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
        }
        legendMeanEta->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);

        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/PositronMeanNSigma_SummaryPtDepInEtaBinsTPCCl%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    for (Int_t r = 0; r < 4; r++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthEta->Clear();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoWidthElectron[1][r]->ProjectionX(Form("%sEta%d",fhistoWidthElectron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
            legendWidthEta->AddEntry(histoProjection,Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]),"p");
        }
        legendWidthEta->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/ElectronWidthNSigma_SummaryPtDepInEtaBinsTPCCl%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        for (Int_t j = 0; j < nEtaBinsOut; j++){
            TH1D* histoProjection       = fhistoWidthPositron[1][r]->ProjectionX(Form("%sEta%d",fhistoWidthPositron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerEta[j], 1.5, colorEta[j], colorEta[j]);
            histoProjection->Draw("same,p,e");
        }
        legendWidthEta->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/PositronWidthNSigma_SummaryPtDepInEtaBinsTPCCl%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    //    ____ _        _____ _______
    //   / ___| |      |  ___|__   __|  /\
    //  | /   |        | |_     | |    /  \
    //  | |   | |      |  _|    | |   / /\ \
    //  | \___| |___   | |___   | |  / ___- \
    //   \____|_____|  |_____|  |_| /_/    \_\

    for (Int_t j = 0; j< nEtaBinsOut; j++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[1][r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("e^{-}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[1][r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("e^{+}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]), fCollisionSystem, "", "", "",0, 31);

        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/MeanNSigma_SummaryPtDepInTPCClBinsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }
    for (Int_t j = 0; j< nEtaBinsOut; j++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[1][r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("e^{-}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[1][r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("e^{+}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]), fCollisionSystem, "", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/WidthNSigma_SummaryPtDepInTPCClBinsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }

    for (Int_t j = 0; j< nEtaBinsOut/2; j++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[1][r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("%d < N TPC Cl < %d, #eta < 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanElectron[1][r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[1][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("%d < N TPC Cl < %d, #eta > 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);

        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/ElectronMeanNSigma_SummaryPtDepInTPCClBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));

        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanPositron[1][r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("%d < N TPC Cl < %d, #eta < 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[1][r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[1][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("%d < N TPC Cl < %d, #eta > 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/PositronMeanNSigma_SummaryPtDepInTPCClBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }

    for (Int_t j = 0; j< nEtaBinsOut/2; j++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[1][r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("%d < N TPC Cl < %d, #eta < 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthElectron[1][r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[1][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("%d < N TPC Cl < %d, #eta > 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);

        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/ElectronWidthNSigma_SummaryPtDepInTPCClBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));

        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthPositron[1][r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[1][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("%d < N TPC Cl < %d, #eta < 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[1][r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[1][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("%d < N TPC Cl < %d, #eta > 0", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/PositronWidthNSigma_SummaryPtDepInTPCClBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));
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
            TH1D* histoProjection       = fhistoMeanElectron[0][r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[0][r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]), fCollisionSystem, "", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/MeanNSigma_SummaryPtDepInRBinsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }
    for (Int_t j = 0; j< nEtaBinsOut; j++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[0][r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[0][r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #eta < %2.1f", arrEtaBinningOut[j], arrEtaBinningOut[j+1]), fCollisionSystem, "", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/WidthNSigma_SummaryPtDepInRBinsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }

    for (Int_t j = 0; j< nEtaBinsOut/2; j++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[0][r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanElectron[0][r]->ProjectionX(Form("%sR%d",fhistoMeanElectron[0][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/ElectronMeanNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));


        histo2DdEdxMeanStudies->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanPositron[0][r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[0][r]->ProjectionX(Form("%sR%d",fhistoMeanPositron[0][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/PositronMeanNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));
    }

    for (Int_t j = 0; j< nEtaBinsOut/2; j++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[0][r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthElectron[0][r]->ProjectionX(Form("%sR%d",fhistoWidthElectron[0][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/ElectronWidthNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));

        histo2DdEdxWidthStudies->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthPositron[0][r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[0][r]->GetName(),j),j+1,j+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta < 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[0][r]->ProjectionX(Form("%sR%d",fhistoWidthPositron[0][r]->GetName(),nEtaBinsOut-j),nEtaBinsOut-j,nEtaBinsOut-j);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("%2.1f < #it{R}_{conv} (cm) < %2.1f, #eta > 0", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < |#eta| < %2.1f", arrEtaBinningOut[nEtaBinsOut-j-1], arrEtaBinningOut[nEtaBinsOut-j]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines(0.01, 20, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/PositronWidthNSigma_SummaryPtDepInRBinsAbsEta%d.%s", outputDir.Data(), j, suffix.Data()));
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

    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << " started with eta projections in " << nPBins << " P bins"<< endl;
    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;
    canvasdEdxMeanStudies->SetLogx(kFALSE);
    for (Int_t r = 0; r < 4; r++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        legendMeanPt->Clear();
        for(Int_t k=startBinP;k<nPBins;k++){
            TH1D* histoProjection       = fhistoMeanElectron[0][r]->ProjectionY(Form("%sPt%d",fhistoMeanElectron[0][r]->GetName(),k),k+1,k+1);
            DrawGammaSetMarker(  histoProjection, markerPt[k], 1.5, colorPt[k], colorPt[k]);
            histoProjection->Draw("same,p,e");
            legendMeanPt->AddEntry(histoProjection,Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[k], arrPBinning[k+1]),"p");
        }
        legendMeanPt->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/ElectronMeanNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        for(Int_t i=startBinP;i<nPBins;i++){
            TH1D* histoProjection       = fhistoMeanPositron[0][r]->ProjectionY(Form("%sPt%d",fhistoMeanPositron[0][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection,  markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
        }
        legendMeanPt->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines( -1, 1, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/PositronMeanNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    canvasdEdxWidthStudies->SetLogx(kFALSE);
    for (Int_t r = 0; r < 4; r++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        legendWidthPt->Clear();
        for (Int_t i = startBinP; i < nPBins; i++){
            TH1D* histoProjection       = fhistoWidthElectron[0][r]->ProjectionY(Form("%sPt%d",fhistoWidthElectron[0][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
            legendWidthPt->AddEntry(histoProjection,Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),"p");
        }
        legendWidthPt->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines( -1, 1, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/ElectronWidthNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        for (Int_t i = startBinP; i < nPBins; i++){
            TH1D* histoProjection       = fhistoWidthPositron[0][r]->ProjectionY(Form("%sPt%d",fhistoWidthPositron[0][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection,  markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
        }
        legendWidthPt->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines( -1, 1, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/PositronWidthNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    //   ____        ____ _
    //  |  _ \      / ___| |
    //  | |_| |    | /   | |
    //  |  __/     | |   | |
    //  | |        | \___| |___
    //  |_|         \____|_____|

    for (Int_t r = 0; r < 4; r++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        legendMeanPt->Clear();
        for(Int_t k=startBinP;k<nPBins;k++){
            TH1D* histoProjection       = fhistoMeanElectron[1][r]->ProjectionY(Form("%sPt%d",fhistoMeanElectron[1][r]->GetName(),k),k+1,k+1);
            DrawGammaSetMarker(  histoProjection, markerPt[k], 1.5, colorPt[k], colorPt[k]);
            histoProjection->Draw("same,p,e");
            legendMeanPt->AddEntry(histoProjection,Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[k], arrPBinning[k+1]),"p");
        }
        legendMeanPt->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines( -1, 1, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/ElectronMeanNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        for(Int_t i=startBinP;i<nPBins;i++){
            TH1D* histoProjection       = fhistoMeanPositron[1][r]->ProjectionY(Form("%sPt%d",fhistoMeanPositron[1][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection,  markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
        }
        legendMeanPt->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines( -1, 1, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/PositronMeanNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    for (Int_t r = 0; r < 4; r++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        legendWidthPt->Clear();
        for (Int_t i = startBinP; i < nPBins; i++){
            TH1D* histoProjection       = fhistoWidthElectron[1][r]->ProjectionY(Form("%sPt%d",fhistoWidthElectron[1][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
            legendWidthPt->AddEntry(histoProjection,Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),"p");
        }
        legendWidthPt->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);
        DrawGammaLines( -1, 1, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/ElectronWidthNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));

        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        for (Int_t i = startBinP; i < nPBins; i++){
            TH1D* histoProjection       = fhistoWidthPositron[1][r]->ProjectionY(Form("%sPt%d",fhistoWidthPositron[1][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection,  markerPt[i], 1.5, colorPt[i], colorPt[i]);
            histoProjection->Draw("same,p,e");
        }
        legendWidthPt->Draw();
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);
        DrawGammaLines( -1, 1, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/PositronWidthNSigma_SummaryEtaDepInPtBinsR%d.%s", outputDir.Data(), r, suffix.Data()));
    }

    //   ____         _____ _______
    //  |  _ \       |  ___|__   __|  /\
    //  | |_| |      | |_     | |    /  \
    //  |  __/       |  _|    | |   / /\ \
    //  | |          | |___   | |  / ___- \
    //  |_|          |_____|  |_| /_/    \_\

    for(Int_t i=startBinP;i<nPBins;i++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[0][r]->ProjectionY(Form("%sR%d",fhistoMeanElectron[0][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[0][r]->ProjectionY(Form("%sR%d",fhistoMeanPositron[0][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]), fCollisionSystem, "", "", "",0, 31);
        DrawGammaLines( -1, 1, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sRConv/MeanNSigma_SummaryEtaDepInRBinsPt%d.%s", outputDir.Data(), i, suffix.Data()));
    }

    for (Int_t i = startBinP; i< nPBins; i++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[0][r]->ProjectionY(Form("%sR%d",fhistoWidthElectron[0][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("e^{-}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[0][r]->ProjectionY(Form("%sR%d",fhistoWidthPositron[0][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("e^{+}, %2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]), fCollisionSystem, "", "", "",0, 31);
        DrawGammaLines( -1, 1, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sRConv/WidthNSigma_SummaryEtaDepInRBinsPt%d.%s", outputDir.Data(), i, suffix.Data()));
    }


    for(Int_t i=startBinP;i<nPBins;i++){
        canvasdEdxMeanStudies->cd();
        histo2DdEdxMeanStudiesEta->DrawCopy();
        legendMeanR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoMeanElectron[1][r]->ProjectionY(Form("%sR%d",fhistoMeanElectron[1][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendMeanR->AddEntry(histoProjection,Form("e^{-}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoMeanPositron[1][r]->ProjectionY(Form("%sR%d",fhistoMeanPositron[1][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendMeanR->AddEntry(histoProjection2,Form("e^{+}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendMeanR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]), fCollisionSystem, "", "", "",0, 31);
        DrawGammaLines( -1, 1, 0.0,0.0, 1, kGray+2 ,7);
        canvasdEdxMeanStudies->SaveAs(Form("%sTPCCl/MeanNSigma_SummaryEtaDepInTPClBinsPt%d.%s", outputDir.Data(), i, suffix.Data()));
    }

    for (Int_t i = startBinP; i< nPBins; i++){
        canvasdEdxWidthStudies->cd();
        histo2DdEdxWidthStudiesEta->DrawCopy();
        legendWidthR->Clear();
        for (Int_t r = 0; r < 4; r++){
            TH1D* histoProjection       = fhistoWidthElectron[1][r]->ProjectionY(Form("%sR%d",fhistoWidthElectron[1][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection, markerR[r], 1.5, colorR[r], colorR[r]);
            histoProjection->Draw("same,p,e");
            legendWidthR->AddEntry(histoProjection,Form("e^{-}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            TH1D* histoProjection2       = fhistoWidthPositron[1][r]->ProjectionY(Form("%sR%d",fhistoWidthPositron[1][r]->GetName(),i),i+1,i+1);
            DrawGammaSetMarker(  histoProjection2, markerR[r+4], 1.5, colorR[r], colorR[r]);
            legendWidthR->AddEntry(histoProjection2,Form("e^{+}, %d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]),"p");
            histoProjection2->Draw("same,p,e");
        }
        legendWidthR->Draw();
        // plot labels
        PlotLabelsInvMassInPtPlots ( 0.95, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.2f < #it{p}_{e} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]), fCollisionSystem, "", "", "",0, 31);
        DrawGammaLines( -1, 1, 1.0,1.0, 1, kGray+2 ,7);
        canvasdEdxWidthStudies->SaveAs(Form("%sTPCCl/WidthNSigma_SummaryEtaDepInTPCClBinsPt%d.%s", outputDir.Data(), i, suffix.Data()));
    }

    TCanvas* canvas2DStudies = new TCanvas("canvas2DStudies","",200,10,1350,1000);  // gives the page size
    DrawGammaCanvasSettings( canvas2DStudies, 0.08, 0.1, 0.02, 0.09);
    canvas2DStudies->SetLogx();

    TCanvas* canvas2DStudiesdEdx = new TCanvas("canvas2DStudiesdEdx","",200,10,1000,1000);  // gives the page size
    DrawGammaCanvasSettings( canvas2DStudiesdEdx, 0.08, 0.1, 0.02, 0.09);
    canvas2DStudiesdEdx->SetLogx();
    canvas2DStudiesdEdx->SetLogz();

    for (Int_t r = 0; r < 4; r++){
        canvas2DStudies->cd();
            SetStyleHistoTH2ForGraphs(fhistoMeanElectron[0][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoMeanElectron[0][r]->GetZaxis()->SetRangeUser(-1.5,1.5);
            fhistoMeanElectron[0][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                         Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "<n#sigma_{e}>, e^{-}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sRConv/Electron_MeanNSigma_2DMap_RBin%d.%s", outputDir.Data(), r, suffix.Data()));
            SetStyleHistoTH2ForGraphs(fhistoMeanPositron[0][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoMeanPositron[0][r]->GetZaxis()->SetRangeUser(-1.5,1.5);
            fhistoMeanPositron[0][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                         Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "<n#sigma_{e}>. e^{+}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sRConv/Positron_MeanNSigma_2DMap_RBin%d.%s", outputDir.Data(), r, suffix.Data()));
            SetStyleHistoTH2ForGraphs(fhistoWidthElectron[0][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoWidthElectron[0][r]->GetZaxis()->SetRangeUser(0,1.8);
            fhistoWidthElectron[0][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                        Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "#sigma(n#sigma_{e}), e^{-}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sRConv/Electron_WidthNSigma_2DMap_RBin%d.%s", outputDir.Data(), r, suffix.Data()));
            SetStyleHistoTH2ForGraphs(fhistoWidthPositron[0][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoWidthPositron[0][r]->GetZaxis()->SetRangeUser(0,1.8);
            fhistoWidthPositron[0][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                         Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "#sigma(n#sigma_{e}), e^{+}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sRConv/Positron_WidthNSigma_2DMap_RBin%d.%s", outputDir.Data(), r, suffix.Data()));

            SetStyleHistoTH2ForGraphs(fhistoMeanElectron[1][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoMeanElectron[1][r]->GetZaxis()->SetRangeUser(-1.5,1.5);
            fhistoMeanElectron[1][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                         Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "<n#sigma_{e}>, e^{-}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sTPCCl/Electron_MeanNSigma_2DMap_ClBin%d.%s", outputDir.Data(), r, suffix.Data()));
            SetStyleHistoTH2ForGraphs(fhistoMeanPositron[1][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoMeanPositron[1][r]->GetZaxis()->SetRangeUser(-1.5,1.5);
            fhistoMeanPositron[1][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                         Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "<n#sigma_{e}>. e^{+}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sTPCCl/Positron_MeanNSigma_2DMap_ClBin%d.%s", outputDir.Data(), r, suffix.Data()));
            SetStyleHistoTH2ForGraphs(fhistoWidthElectron[1][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoWidthElectron[1][r]->GetZaxis()->SetRangeUser(0,1.8);
            fhistoWidthElectron[1][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                         Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "#sigma(n#sigma_{e}), e^{-}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sTPCCl/Electron_WidthNSigma_2DMap_ClBin%d.%s", outputDir.Data(), r, suffix.Data()));
            SetStyleHistoTH2ForGraphs(fhistoWidthPositron[1][r], "#it{p} (GeV/#it{c})", "#eta", 0.035,0.04, 0.035,0.04, 0.98,1.0);
            fhistoWidthPositron[1][r]->GetZaxis()->SetRangeUser(0,1.8);
            fhistoWidthPositron[1][r]->Draw("colz");
            PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                         Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "#sigma(n#sigma_{e}), e^{+}", "", "",0, 31);
        canvas2DStudies->SaveAs(Form("%sTPCCl/Positron_WidthNSigma_2DMap_ClBin%d.%s", outputDir.Data(), r, suffix.Data()));

        canvas2DStudiesdEdx->cd();

        TH2D* histoElectronSigmadEdxP_R   = (TH2D*)histoElectronDeDxPEta[1][r]->Project3D("xz");
        SetStyleHistoTH2ForGraphs(histoElectronSigmadEdxP_R, "#it{p} (GeV/#it{c})", "n#sigma_{e} d#it{E}/d#it{x}", 0.035,0.04, 0.035,0.04, 0.98,1.0);
        histoElectronSigmadEdxP_R->Draw("colz");
        PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvas2DStudiesdEdx->SaveAs(Form("%sRConv/Electron_TPCNSigmadEdx_RBin%d.%s", outputDir.Data(), r, suffix.Data()));
        delete histoElectronSigmadEdxP_R;
        TH2D* histoPositronSigmadEdxP_R   = (TH2D*)histoPositronDeDxPEta[1][r]->Project3D("xz");
        SetStyleHistoTH2ForGraphs(histoPositronSigmadEdxP_R, "#it{p} (GeV/#it{c})", "n#sigma_{e} d#it{E}/d#it{x}", 0.035,0.04, 0.035,0.04, 0.98,1.0);
        histoPositronSigmadEdxP_R->Draw("colz");
        PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvas2DStudiesdEdx->SaveAs(Form("%sRConv/Positron_TPCNSigmadEdx_RBin%d.%s", outputDir.Data(), r, suffix.Data()));
        delete histoPositronSigmadEdxP_R;

        TH2D* histoElectronSigmadEdxP_Cl   = (TH2D*)histoElectronDeDxPEta[1][r]->Project3D("xz");
        SetStyleHistoTH2ForGraphs(histoElectronSigmadEdxP_Cl, "#it{p} (GeV/#it{c})", "n#sigma_{e} d#it{E}/d#it{x}", 0.035,0.04, 0.035,0.04, 0.98,1.0);
        histoElectronSigmadEdxP_Cl->Draw("colz");
        PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{-}", "", "",0, 31);

        canvas2DStudiesdEdx->SaveAs(Form("%sTPCCl/Electron_TPCNSigmadEdx_ClBin%d.%s", outputDir.Data(), r, suffix.Data()));
        delete histoElectronSigmadEdxP_Cl;
        TH2D* histoPositronSigmadEdxP_Cl   = (TH2D*)histoPositronDeDxPEta[1][r]->Project3D("xz");
        SetStyleHistoTH2ForGraphs(histoPositronSigmadEdxP_Cl, "#it{p} (GeV/#it{c})", "n#sigma_{e} d#it{E}/d#it{x}", 0.035,0.04, 0.035,0.04, 0.98,1.0);
        histoPositronSigmadEdxP_Cl->Draw("colz");
        PlotLabelsInvMassInPtPlots ( 0.85, 0.9, 0.035, 0.04*1.05, "ALICE performance",
                                     Form("%d < N TPC Cl < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1]), fCollisionSystem, "e^{+}", "", "",0, 31);

        canvas2DStudiesdEdx->SaveAs(Form("%sTPCCl/Positron_TPCNSigmadEdx_ClBin%d.%s", outputDir.Data(), r, suffix.Data()));
        delete histoPositronSigmadEdxP_Cl;

    }


    cout << "------------------------------------------------------------------------------------" << endl;
    cout << "------------------------------------------------------------------------------------" << endl;
    cout << "Finished with SUMMARY PLOTS, STARTED PLOTTING detailed fits" << endl;
    cout << "------------------------------------------------------------------------------------" << endl;
    cout << "------------------------------------------------------------------------------------" << endl;

    Int_t fColumnPlot = 5;
    Int_t fRowPlot    = 4;
    TString nameCanvasElec = "";
    TString namePadElec    = "";
    TString nameCanvasPosit  = "";
    TString namePadPosit     = "";
    for(Int_t i=startBinP;i<nPBins;i++){
        for (Int_t r = 0; r< 4; r++){
            cout << i << "\t" << r << endl;
            // Plotting Electron
            nameCanvasElec = Form("Rbin%d_ElecNSigmaDeDx_Pbin%d",r,i);
            namePadElec    = Form("Rbin%d_ElecNSigmaDeDx_Pbin%d",r,i);
            PlotdEdxSlices( histoElectronDeDx[0][r][i],fitElectronDeDx[0][r][i],meanElectron[0][r][i],Form("%sRConv/DetailedFits/MapsElec_R%02d_P%02d.%s",outputDir.Data(),r,i,suffix.Data()),
                            nameCanvasElec,namePadElec,arrEtaBinningOut,optionMC.Data(),fCollisionSystem,nEtaBins,fColumnPlot,fRowPlot,
                            Form("%2.2f < #it{p}_{e^{-}} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),
                            Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1])
                          );
            nameCanvasElec = Form("Clbin%d_ElecNSigmaDeDx_Pbin%d",r,i);
            namePadElec    = Form("Clbin%d_ElecNSigmaDeDx_Pbin%d",r,i);
            PlotdEdxSlices( histoElectronDeDx[1][r][i],fitElectronDeDx[1][r][i],meanElectron[1][r][i],Form("%sTPCCl/DetailedFits/MapsElec_Cl%02d_P%02d.%s",outputDir.Data(),r,i,suffix.Data()),
                            nameCanvasElec,namePadElec,arrEtaBinningOut,optionMC.Data(),fCollisionSystem,nEtaBins,fColumnPlot,fRowPlot,
                            Form("%2.2f < #it{p}_{e^{-}} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),
                            Form("%d < N TPC cl. < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1])
            );

            // Ploting Positron
            nameCanvasPosit  = Form("Rbin%d_PositNSigmaDeDx_Pbin%d",r,i);
            namePadPosit     = Form("Rbin%d_PositNSigmaDeDx_Pbin%d",r,i);
            PlotdEdxSlices( histoPositronDeDx[0][r][i],fitPositronDeDx[0][r][i],meanPositron[0][r][i],Form("%sRConv/DetailedFits/MapsPos_R%02d_P%02d.%s",outputDir.Data(),r,i,suffix.Data()),
                            nameCanvasPosit,namePadPosit,arrEtaBinningOut,optionMC.Data(),fCollisionSystem,nEtaBins,fColumnPlot,fRowPlot,
                            Form("%2.2f < #it{p}_{e^{+}} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),
                            Form("%2.1f < #it{R}_{conv} (cm) < %2.1f", arrRBinningOut[r], arrRBinningOut[r+1])
                          );
            nameCanvasPosit  = Form("Clbin%d_PositNSigmaDeDx_Pbin%d",r,i);
            namePadPosit     = Form("Clbin%d_PositNSigmaDeDx_Pbin%d",r,i);
            PlotdEdxSlices( histoPositronDeDx[1][r][i],fitPositronDeDx[1][r][i],meanPositron[1][r][i],Form("%sTPCCl/DetailedFits/MapsPos_Cl%02d_P%02d.%s",outputDir.Data(),r,i,suffix.Data()),
                            nameCanvasPosit,namePadPosit,arrEtaBinningOut,optionMC.Data(),fCollisionSystem,nEtaBins,fColumnPlot,fRowPlot,
                            Form("%2.2f < #it{p}_{e^{+}} (GeV/#it{c}) < %2.2f", arrPBinning[i], arrPBinning[i+1]),
                            Form("%d < N TPC cl. < %d", arrTPCClBinningOut[r], arrTPCClBinningOut[r+1])
            );
        }
    }

    TFile outFileMonitoring(Form("%sMonitoringDeDxMaps_%s.root",outputDir.Data(),fCutSelectionRead.Data()) ,"RECREATE");
        for(Int_t i=startBinP;i<nPBins;i++){
            for(Int_t j=0;j<nEtaBins;j++){
                for (Int_t r = 0; r< 4; r++){
                    histoElectronDeDx[0][r][i][j]->Write();
                    histoPositronDeDx[0][r][i][j]->Write();
                }
            }
        }
    outFileMonitoring.Close();

    TFile outFileMaps(Form("%sDeDxMaps_%s.root",outputDir.Data(),fCutSelectionRead.Data()) ,"RECREATE");
    for (Int_t r = 0; r< 4; r++){
        fhistoMeanPositron[0][r]->Write(Form("Pos_R%d_mean",r));
        fhistoMeanElectron[0][r]->Write(Form("Ele_R%d_mean",r));
        fhistoWidthPositron[0][r]->Write(Form("Pos_R%d_width",r));
        fhistoWidthElectron[0][r]->Write(Form("Ele_R%d_width",r));
        fhistoMeanPositron[1][r]->Write(Form("Pos_Cl%d_mean",r));
        fhistoMeanElectron[1][r]->Write(Form("Ele_Cl%d_mean",r));
        fhistoWidthPositron[1][r]->Write(Form("Pos_Cl%d_width",r));
        fhistoWidthElectron[1][r]->Write(Form("Ele_Cl%d_width",r));

        fhistoRecalibPositron[0][r]->Write(Form("Pos_R%d_recalib",r));
        fhistoRecalibElectron[0][r]->Write(Form("Ele_R%d_recalib",r));
        fhistoRecalibPositron[1][r]->Write(Form("Pos_Cl%d_recalib",r));
        fhistoRecalibElectron[1][r]->Write(Form("Ele_Cl%d_recalib",r));
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
            DrawGammaHisto( fHistoMapsPlot[j],title,"n#sigma","Counts",-6.,6.,0,0.7);
            DrawGammaLines(0.,0.,0.,fHistoMapsPlot[j]->GetMaximum(), 1, kGray+2, 7);
            DrawGammaLines(fHistoMeanPlot[j], fHistoMeanPlot[j],0.,fHistoMapsPlot[j]->GetMaximum(), 1, kBlue+2, 2);


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
