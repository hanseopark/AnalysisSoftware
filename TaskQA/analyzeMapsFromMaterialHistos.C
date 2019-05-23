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
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"

void FitSignal(TH1D* hSig,Color_t col);
void PlotdEdxSlices(TH1D** , Double_t* , TString , TString , TString , Double_t* , TString , TString , Int_t , Int_t , Int_t );
void DrawSliceHisto(TH1* ,  TString , TString , TString , Float_t , Float_t , Size_t , Color_t );

TF1*     fGaus       = NULL;
Double_t Mean=0.0;
Double_t Width=0.0;
Double_t Chi2=0.0;
Double_t meanElectronR0[12][20];//eta pt
Double_t meanElectronR1[12][20];//eta pt
Double_t meanElectronR2[12][20];//eta pt
Double_t meanElectronR3[12][20];//eta pt

Double_t meanPositronR0[12][20];//eta pt
Double_t meanPositronR1[12][20];//eta pt
Double_t meanPositronR2[12][20];//eta pt
Double_t meanPositronR3[12][20];//eta pt


Double_t widthElectronR0[12][20];
Double_t widthElectronR1[12][20];
Double_t widthElectronR2[12][20];
Double_t widthElectronR3[12][20];

Double_t widthPositronR0[12][20];
Double_t widthPositronR1[12][20];
Double_t widthPositronR2[12][20];
Double_t widthPositronR3[12][20];
Double_t chi[20][20];

void analyzeMapsFromMaterialHistos(TString fileNameWithMaps="" ,
                                   TString cutSelection = "",
                                   TString optionMC = "Data",
                                   TString fEnergy = ""
)
{
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetEndErrorSize(0);
    StyleSettingsThesis();
    SetPlotStyle();

    TString outputDir = "";
    if(optionMC.CompareTo("Data")==0)    outputDir = "DeDxMapsData";
    else if(optionMC.CompareTo("MC")==0) outputDir = "DeDxMapsMC";
    gSystem->Exec(Form("mkdir -p %s",outputDir.Data()));

    TH3F *histoPositronDeDxPEtaR0 = NULL;
    TH3F *histoPositronDeDxPEtaR1 = NULL;
    TH3F *histoPositronDeDxPEtaR2 = NULL;
    TH3F *histoPositronDeDxPEtaR3 = NULL;

    TH3F *histoElectronDeDxPEtaR0 = NULL;
    TH3F *histoElectronDeDxPEtaR1 = NULL;
    TH3F *histoElectronDeDxPEtaR2 = NULL;
    TH3F *histoElectronDeDxPEtaR3 = NULL;

    TFile* fileMaterialHistos = new TFile(fileNameWithMaps.Data());
    TString fCutSelectionRead = cutSelection;
    TString nameMainDir = "GammaConvMaterial";

    TList *TopDir =(TList*)fileMaterialHistos->Get(nameMainDir.Data());
    if(TopDir != NULL){
        TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
        if(HistosGammaConversion == NULL){
            cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
            return;
        }

        TList *MapsContainer           = (TList*)HistosGammaConversion->FindObject(Form("%s  dEdx Maps",fCutSelectionRead.Data()));
        histoPositronDeDxPEtaR0 = (TH3F*)MapsContainer->FindObject("R0 positron sigma dEdx P Eta");
        histoPositronDeDxPEtaR1 = (TH3F*)MapsContainer->FindObject("R1 positron sigma dEdx P Eta");
        histoPositronDeDxPEtaR2 = (TH3F*)MapsContainer->FindObject("R2 positron sigma dEdx P Eta");
        histoPositronDeDxPEtaR3 = (TH3F*)MapsContainer->FindObject("R3 positron sigma dEdx P Eta");

        histoElectronDeDxPEtaR0 = (TH3F*)MapsContainer->FindObject("R0 electron sigma dEdx P Eta");
        histoElectronDeDxPEtaR1 = (TH3F*)MapsContainer->FindObject("R1 electron sigma dEdx P Eta");
        histoElectronDeDxPEtaR2 = (TH3F*)MapsContainer->FindObject("R2 electron sigma dEdx P Eta");
        histoElectronDeDxPEtaR3 = (TH3F*)MapsContainer->FindObject("R3 electron sigma dEdx P Eta");
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
        histoPositronDeDxPEtaR0 = (TH3F*)directoryConv->Get("R0 positron sigma dEdx P Eta");
        histoPositronDeDxPEtaR1 = (TH3F*)directoryConv->Get("R1 positron sigma dEdx P Eta");
        histoPositronDeDxPEtaR2 = (TH3F*)directoryConv->Get("R2 positron sigma dEdx P Eta");
        histoPositronDeDxPEtaR3 = (TH3F*)directoryConv->Get("R3 positron sigma dEdx P Eta");

        histoElectronDeDxPEtaR0 = (TH3F*)directoryConv->Get("R0 electron sigma dEdx P Eta");
        histoElectronDeDxPEtaR1 = (TH3F*)directoryConv->Get("R1 electron sigma dEdx P Eta");
        histoElectronDeDxPEtaR2 = (TH3F*)directoryConv->Get("R2 electron sigma dEdx P Eta");
        histoElectronDeDxPEtaR3 = (TH3F*)directoryConv->Get("R3 electron sigma dEdx P Eta");
    }

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
//         cout << "p bins:: "<< i << " " <<  arrPBinning[i]<< endl;
    }

    Int_t nEtaBins = 19;
    Int_t nEtaBinsOut=18;
    Double_t *arrEtaBinningOut = new Double_t[19];
    for( Int_t i=0;i<nEtaBinsOut+1;i++){
        arrEtaBinningOut[i]= -0.9+0.1*i;
//         cout << "eta bins:: " << i << " " << arrEtaBinningOut[i] << endl;
    }


    TH2F* fhistoMeanPositronR0 = new TH2F("MeanPosiR0","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoMeanPositronR1 = new TH2F("MeanPosiR1","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoMeanPositronR2 = new TH2F("MeanPosiR2","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoMeanPositronR3 = new TH2F("MeanPosiR3","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);

    TH2F* fhistoMeanElectronR0 = new TH2F("MeanElecR0","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoMeanElectronR1 = new TH2F("MeanElecR1","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoMeanElectronR2 = new TH2F("MeanElecR2","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoMeanElectronR3 = new TH2F("MeanElecR3","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);

    TH2F* fhistoWidthPositronR0 = new TH2F("WidthPosiR0","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoWidthPositronR1 = new TH2F("WidthPosiR1","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoWidthPositronR2 = new TH2F("WidthPosiR2","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoWidthPositronR3 = new TH2F("WidthPosiR3","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);

    TH2F* fhistoWidthElectronR0 = new TH2F("WidthElecR0","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoWidthElectronR1 = new TH2F("WidthElecR1","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoWidthElectronR2 = new TH2F("WidthElecR2","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);
    TH2F* fhistoWidthElectronR3 = new TH2F("WidthElecR3","",nPBins,arrPBinning,nEtaBinsOut,arrEtaBinningOut);

    //  fhistomean->SetTitle("Mean;#it{p} (GeV/c);#eta ");

    TH1D*histoPositronR0DeDx[12][20];
    TH1D*histoPositronR1DeDx[12][20];
    TH1D*histoPositronR2DeDx[12][20];
    TH1D*histoPositronR3DeDx[12][20];

    TH1D*histoElectronR0DeDx[12][20];
    TH1D*histoElectronR1DeDx[12][20];
    TH1D*histoElectronR2DeDx[12][20];
    TH1D*histoElectronR3DeDx[12][20];

    Int_t etaOff = 1;
    if(fEnergy.CompareTo("5TeV2017")==0) etaOff = 0;

    //    cout<< histoPositronDeDxPEtaR0->GetNbinsY()<< " " << histoPositronDeDxPEtaR0->GetNbinsZ() << endl;
    for(Int_t i=0;i<nPBins;i++){
//         cout << "*****************\np = " << arrPBinning[i] << " GeV/c" << endl;
        for(Int_t j=0;j<nEtaBins;j++){
//             cout << "-> eta " << arrEtaBinningOut[j] << endl;
//             cout << "offset " << arrEtaBinningOut[j-etaOff] << endl;
//             cout << " j-etaOff " << j-etaOff << endl;

            // Electron, R bin 0
            histoElectronR0DeDx[i][j] = (TH1D*)histoElectronDeDxPEtaR0->ProjectionX(Form("R0ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoElectronR0DeDx[i][j],kBlue);
            meanElectronR0[i][j] = Mean;
            widthElectronR0[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanElectronR0->SetBinContent(i+1,j+1-etaOff, meanElectronR0[i][j]);
                fhistoWidthElectronR0->SetBinContent(i+1,j+1-etaOff, widthElectronR0[i][j]);
            }

            // Electron, R bin 1
            histoElectronR1DeDx[i][j] = (TH1D*)histoElectronDeDxPEtaR1->ProjectionX(Form("R1ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoElectronR1DeDx[i][j],kBlue);
            meanElectronR1[i][j] = Mean;
            widthElectronR1[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanElectronR1->SetBinContent(i+1,j+1-etaOff, meanElectronR1[i][j]);
                fhistoWidthElectronR1->SetBinContent(i+1,j+1-etaOff, widthElectronR1[i][j]);
            }

            // Electron, R bin 2
            histoElectronR2DeDx[i][j] = (TH1D*)histoElectronDeDxPEtaR2->ProjectionX(Form("R2ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoElectronR2DeDx[i][j],kBlue);
            meanElectronR2[i][j] = Mean;
            widthElectronR2[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanElectronR2->SetBinContent(i+1,j+1-etaOff, meanElectronR2[i][j]);
                fhistoWidthElectronR2->SetBinContent(i+1,j+1-etaOff, widthElectronR2[i][j]);
            }

            // Electron, R bin 3
            histoElectronR3DeDx[i][j] = (TH1D*)histoElectronDeDxPEtaR3->ProjectionX(Form("R3ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoElectronR3DeDx[i][j],kBlue);
            meanElectronR3[i][j] = Mean;
            widthElectronR3[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanElectronR3->SetBinContent(i+1,j+1-etaOff, meanElectronR3[i][j]);
                fhistoWidthElectronR3->SetBinContent(i+1,j+1-etaOff, widthElectronR3[i][j]);
            }


            // Positron, R bin 0
            histoPositronR0DeDx[i][j] = (TH1D*)histoPositronDeDxPEtaR0->ProjectionX(Form("R0PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoPositronR0DeDx[i][j],kBlue);
            meanPositronR0[i][j] = Mean;
            widthPositronR0[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanPositronR0->SetBinContent(i+1,j+1-etaOff, meanPositronR0[i][j]);
                fhistoWidthPositronR0->SetBinContent(i+1,j+1-etaOff, widthPositronR0[i][j]);
            }

            // Positron, R bin 1
            histoPositronR1DeDx[i][j] = (TH1D*)histoPositronDeDxPEtaR1->ProjectionX(Form("R1PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoPositronR1DeDx[i][j],kBlue);
            meanPositronR1[i][j] = Mean;
            widthPositronR1[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanPositronR1->SetBinContent(i+1,j+1-etaOff, meanPositronR1[i][j]);
                fhistoWidthPositronR1->SetBinContent(i+1,j+1-etaOff, widthPositronR1[i][j]);
            }

            // Positron, R bin 2
            histoPositronR2DeDx[i][j] = (TH1D*)histoPositronDeDxPEtaR2->ProjectionX(Form("R2PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoPositronR2DeDx[i][j],kBlue);
            meanPositronR2[i][j] = Mean;
            widthPositronR2[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanPositronR2->SetBinContent(i+1,j+1-etaOff, meanPositronR2[i][j]);
                fhistoWidthPositronR2->SetBinContent(i+1,j+1-etaOff, widthPositronR2[i][j]);
            }

            // Positron, R bin 3
            histoPositronR3DeDx[i][j] = (TH1D*)histoPositronDeDxPEtaR3->ProjectionX(Form("R3PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
            FitSignal(histoPositronR3DeDx[i][j],kBlue);
            meanPositronR3[i][j] = Mean;
            widthPositronR3[i][j] = Width;
            if( (j-etaOff) >= 0 && (j-etaOff) < nEtaBinsOut ){
                fhistoMeanPositronR3->SetBinContent(i+1,j+1-etaOff, meanPositronR3[i][j]);
                fhistoWidthPositronR3->SetBinContent(i+1,j+1-etaOff, widthPositronR3[i][j]);
            }
        }
    }

    Int_t fColumnPlot = 4;
    Int_t fRowPlot    = 5;
    TString nameCanvasElec = "";
    TString namePadElec    = "";
    TString nameCanvasPosit  = "";
    TString namePadPosit     = "";
    for(Int_t i=0;i<nPBins;i++){

        // Plotting Electron
        nameCanvasElec = Form("Rbin0_ElecNSigmaDeDx_Pbin%d",i);
        namePadElec    = Form("Rbin0_ElecNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoElectronR0DeDx[i],meanElectronR0[i],Form("%s/MapsElec_R0_P%d.pdf",outputDir.Data(),i),nameCanvasElec,namePadElec,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);

        nameCanvasElec  = Form("Rbin1_ElecNSigmaDeDx_Pbin%d",i);
        namePadElec     = Form("Rbin1_ElecNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoElectronR1DeDx[i],meanElectronR1[i],Form("%s/MapsElec_R1_P%d.pdf",outputDir.Data(),i),nameCanvasElec,namePadElec,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);

        nameCanvasElec  = Form("Rbin2_ElecNSigmaDeDx_Pbin%d",i);
        namePadElec     = Form("Rbin2_ElecNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoElectronR2DeDx[i],meanElectronR2[i],Form("%s/MapsElec_R2_P%d.pdf",outputDir.Data(),i),nameCanvasElec,namePadElec,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);

        nameCanvasElec  = Form("Rbin3_ElecNSigmaDeDx_Pbin%d",i);
        namePadElec     = Form("Rbin3_ElecNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoElectronR3DeDx[i],meanElectronR3[i],Form("%s/MapsElec_R3_P%d.pdf",outputDir.Data(),i),nameCanvasElec,namePadElec,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);


        // Ploting Positron
        nameCanvasPosit  = Form("Rbin0_PositNSigmaDeDx_Pbin%d",i);
        namePadPosit     = Form("Rbin0_PositNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoPositronR0DeDx[i],meanPositronR0[i],Form("%s/MapsPos_R0_P%d.pdf",outputDir.Data(),i),nameCanvasPosit,namePadPosit,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);

        nameCanvasPosit  = Form("Rbin1_PositNSigmaDeDx_Pbin%d",i);
        namePadPosit     = Form("Rbin1_PositNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoPositronR1DeDx[i],meanPositronR1[i],Form("%s/MapsPos_R1_P%d.pdf",outputDir.Data(),i),nameCanvasPosit,namePadPosit,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);

        nameCanvasPosit  = Form("Rbin2_PositNSigmaDeDx_Pbin%d",i);
        namePadPosit     = Form("Rbin2_PositNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoPositronR2DeDx[i],meanPositronR2[i],Form("%s/MapsPos_R2_P%d.pdf",outputDir.Data(),i),nameCanvasPosit,namePadPosit,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);

        nameCanvasPosit  = Form("Rbin3_PositNSigmaDeDx_Pbin%d",i);
        namePadPosit     = Form("Rbin3_PositNSigmaDeDx_Pbin%d",i);
        PlotdEdxSlices(histoPositronR3DeDx[i],meanPositronR3[i],Form("%s/MapsPos_R3_P%d.pdf",outputDir.Data(),i),nameCanvasPosit,namePadPosit,arrEtaBinningOut,optionMC.Data(),fEnergy,nEtaBins,fColumnPlot,fRowPlot);
    }

    TFile outFileMonitoring(Form("%s/MonitoringDeDxMaps_%s.root",outputDir.Data(),fCutSelectionRead.Data()) ,"RECREATE");
        for(Int_t i=0;i<nPBins;i++){
            for(Int_t j=0;j<nEtaBins;j++){
                histoElectronR0DeDx[i][j]->Write();
                histoElectronR1DeDx[i][j]->Write();
                histoElectronR2DeDx[i][j]->Write();
                histoElectronR3DeDx[i][j]->Write();
                histoPositronR0DeDx[i][j]->Write();
                histoPositronR1DeDx[i][j]->Write();
                histoPositronR2DeDx[i][j]->Write();
                histoPositronR3DeDx[i][j]->Write();
            }
        }
    outFileMonitoring.Close();

    TFile outFileMaps(Form("%s/DeDxMaps_%s.root",outputDir.Data(),fCutSelectionRead.Data()) ,"RECREATE");
        fhistoMeanPositronR0->Write("Pos_R0_mean");
        fhistoMeanPositronR1->Write("Pos_R1_mean");
        fhistoMeanPositronR2->Write("Pos_R2_mean");
        fhistoMeanPositronR3->Write("Pos_R3_mean");

        fhistoMeanElectronR0->Write("Ele_R0_mean");
        fhistoMeanElectronR1->Write("Ele_R1_mean");
        fhistoMeanElectronR2->Write("Ele_R2_mean");
        fhistoMeanElectronR3->Write("Ele_R3_mean");

        fhistoWidthPositronR0->Write("Pos_R0_width");
        fhistoWidthPositronR1->Write("Pos_R1_width");
        fhistoWidthPositronR2->Write("Pos_R2_width");
        fhistoWidthPositronR3->Write("Pos_R3_width");

        fhistoWidthElectronR0->Write("Ele_R0_width");
        fhistoWidthElectronR1->Write("Ele_R1_width");
        fhistoWidthElectronR2->Write("Ele_R2_width");
        fhistoWidthElectronR3->Write("Ele_R3_width");
    outFileMaps.Close();
 }

//___________________________________________________
void FitSignal(TH1D* hSig,Color_t col){

    fGaus= NULL;
    fGaus = new TF1("f1","gaus",-4,5);
    fGaus->SetParameters(hSig->GetMaximum(),0.,1);
    //  fGaus = new TF1("Gaussian","[0]*(exp(-0.5*((x-[1])/[2])^2)",-4,5);

    fGaus->SetLineColor(col);
    fGaus->SetLineWidth(1.);
    if(hSig->GetEntries()<40){
        Mean  = 0.;
        Width = 1.;
    } else {

        if (hSig->GetEntries() >=40 && hSig->GetEntries()<150) hSig->Rebin(4);
        hSig->Fit(fGaus,"QR+","",-3.5,+3.5);

        Double_t xMean = fGaus->GetParameter(1);
        fGaus->SetLineColor(kRed);
        hSig->Fit(fGaus,"QR+","",xMean-2.,xMean+2.);
        Mean  = fGaus->GetParameter(1);
        Width = fGaus->GetParameter(2);
        Chi2 = fGaus->GetChisquare()/fGaus->GetNDF();
        if (Width > 2){
            cout << Mean << " " << Width << "  " << hSig->GetEntries() << " " << hSig->GetMaximum()<< endl;
        }
    }
}


//__________________________________________ Plotting all Invariant Mass bins _______________________________________________
void PlotdEdxSlices(TH1D** fHistoMapsPlot,
                    Double_t* fHistoMeanPlot,
                    TString namePlot,
                    TString nameCanvas,
                    TString namePad,
                    Double_t* arrEtaBinnin,
                    TString fMonteCarloInfo,
                    TString fEnergy,
                    Int_t nEtaBins,
                    Int_t fColumnPlot,
                    Int_t fRowPlot
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
            padDataSpectra->cd(place)->SetBottomMargin(0.15);
            padDataSpectra->cd(place)->SetRightMargin(0.02);
            int remaining           = (int)((place-1)%fColumnPlot);
            if (remaining > 0) padDataSpectra->cd(place)->SetLeftMargin(0.15);
            else padDataSpectra->cd(place)->SetLeftMargin(0.25);

            TString title = Form("%3.2f < | #eta | < %3.2f GeV/#it{c}",startEta,endEta);

//             TH2F * histo2Ddummy = new TH2F("histo2Ddummy","histo2Ddummy",1000,-7, 7,1000,0.,1.3*fHistoMapsPlot[j]->GetMaximum());
//             SetStyleHistoTH2ForGraphs(histo2Ddummy, "n#sigma","",0.035,0.04, 0.035,0.04, 1.,1.);
//             histo2Ddummy->Draw("copy");
//             DrawSliceHisto( fHistoMapsPlot[j],title,"n#sigma","Counts",-5.,5.,0.5,kBlack);

            DrawGammaHisto( fHistoMapsPlot[j],title,"n#sigma","Counts",-6.,6.,0.5,kBlack);
            DrawGammaLines(0.,0.,0.,fHistoMapsPlot[j]->GetMaximum(), 1, kGray+2, 7);
            DrawGammaLines(fHistoMeanPlot[j], fHistoMeanPlot[j],0.,fHistoMapsPlot[j]->GetMaximum(), 1, kBlue+2, 2);

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
    TString dateDummy       = "";
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
    PlotLabelsInvMassInPtPlots ( startTextX, startTextY, textHeight, differenceText, textAlice, dateDummy, fEnergy, "", "", "");

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