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
#include "CommonHeaders/PlottingMeson.h"
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
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"

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
void FitSignal(TH1D* hSig,Color_t col);

Double_t chi[20][20];

void analyzeMapsFromMaterialHistos(TString fileNameWithMaps="" , TString cutSelection = "", TString optionMC = "Data")
{
  gSystem->Exec("mkdir -p DeDxMapsData");
  gSystem->Exec("mkdir -p DeDxMapsMC");
  //TString fileNameMaps = "",
  //  TString fileNameWithMaps = "/Users/marin/analysis/2018/Grid/OutputLegoTrains/pp/Legotrain-vAN-20180710-13TeV-MaterialHistos-Maps/GammaConv_Material_LHC16d_21.root";
  TFile* fileMaterialHistos = new TFile(fileNameWithMaps.Data());

  TString fCutSelectionRead = cutSelection;

  TString nameMainDir = "GammaConvMaterial";
  TList *TopDir =(TList*)fileMaterialHistos->Get(nameMainDir.Data());
    if(TopDir == NULL){
        cout<<"ERROR: TopDir not Found"<<endl;
        return;
    }
    TList *HistosGammaConversion = (TList*)TopDir->FindObject(Form("Cut Number %s",fCutSelectionRead.Data()));
    if(HistosGammaConversion == NULL){
        cout<<"ERROR: " << Form("Cut Number %s",fCutSelectionRead.Data()) << " not Found in File"<<endl;
        return;
    }


    TList *MapsContainer           = (TList*)HistosGammaConversion->FindObject(Form("%s  dEdx Maps",fCutSelectionRead.Data()));

 
    TH3F *histoPositronDeDxPEtaR0 = (TH3F*)MapsContainer->FindObject("R0 positron sigma dEdx P Eta");
    TH3F *histoPositronDeDxPEtaR1 = (TH3F*)MapsContainer->FindObject("R1 positron sigma dEdx P Eta");
    TH3F *histoPositronDeDxPEtaR2 = (TH3F*)MapsContainer->FindObject("R2 positron sigma dEdx P Eta");
    TH3F *histoPositronDeDxPEtaR3 = (TH3F*)MapsContainer->FindObject("R3 positron sigma dEdx P Eta");

    TH3F *histoElectronDeDxPEtaR0 = (TH3F*)MapsContainer->FindObject("R0 electron sigma dEdx P Eta");
    TH3F *histoElectronDeDxPEtaR1 = (TH3F*)MapsContainer->FindObject("R1 electron sigma dEdx P Eta");
    TH3F *histoElectronDeDxPEtaR2 = (TH3F*)MapsContainer->FindObject("R2 electron sigma dEdx P Eta");
    TH3F *histoElectronDeDxPEtaR3 = (TH3F*)MapsContainer->FindObject("R3 electron sigma dEdx P Eta");



    Int_t nPBins =12;
    Int_t nEtaBins =20;
    Double_t *arrPBinning = new Double_t[13]; 
    for( Int_t i=0;i<nPBins+1;i++){
      if(i==0){
	arrPBinning[i]= 0.05;
      }else if(i>0 && i<11){
	arrPBinning[i]= 0.1*i;
      }else if(i==11){
	arrPBinning[i]= 2.0;
      }else if(i==12){
	arrPBinning[i]= 10.0;
      }
      //cout<< "pbins::"<< i << " " <<  arrPBinning[i]<< endl;
    }
    Double_t *arrEtaBinning      = new Double_t[21]; 
    for( Int_t i=0;i<nEtaBins+1;i++){
      arrEtaBinning[i]= -1.+0.1*i;
      //cout<< "Etabins::"<< i << " " <<  arrEtaBinning[i]<< endl;
    }
    TH2F* fhistoMeanPositronR0 = new TH2F("MeanPosiR0","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoMeanPositronR1 = new TH2F("MeanPosiR1","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoMeanPositronR2 = new TH2F("MeanPosiR2","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoMeanPositronR3 = new TH2F("MeanPosiR3","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);

    TH2F* fhistoMeanElectronR0 = new TH2F("MeanElecR0","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoMeanElectronR1 = new TH2F("MeanElecR1","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoMeanElectronR2 = new TH2F("MeanElecR2","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoMeanElectronR3 = new TH2F("MeanElecR3","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);

    TH2F* fhistoWidthPositronR0 = new TH2F("WidthPosiR0","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoWidthPositronR1 = new TH2F("WidthPosiR1","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoWidthPositronR2 = new TH2F("WidthPosiR2","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoWidthPositronR3 = new TH2F("WidthPosiR3","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);

    TH2F* fhistoWidthElectronR0 = new TH2F("WidthElecR0","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoWidthElectronR1 = new TH2F("WidthElecR1","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoWidthElectronR2 = new TH2F("WidthElecR2","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);
    TH2F* fhistoWidthElectronR3 = new TH2F("WidthElecR3","",nPBins,arrPBinning,nEtaBins,arrEtaBinning);


    //  fhistomean->SetTitle("Mean;#it{p} (GeV/c);#eta ");



    TH1D*histoPositronR0DeDx[12][20];
    TH1D*histoPositronR1DeDx[12][20];
    TH1D*histoPositronR2DeDx[12][20];
    TH1D*histoPositronR3DeDx[12][20];

    TH1D*histoElectronR0DeDx[12][20];
    TH1D*histoElectronR1DeDx[12][20];
    TH1D*histoElectronR2DeDx[12][20];
    TH1D*histoElectronR3DeDx[12][20];



    //    cout<< histoPositronDeDxPEtaR0->GetNbinsY()<< " " << histoPositronDeDxPEtaR0->GetNbinsZ() << endl;
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	//	cout<< " i,j,R0::" << i << " "<< j<<endl;
	histoPositronR0DeDx[i][j] = new TH1D(Form("R0PosiSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoPositronR0DeDx[i][j] =
	  (TH1D*)histoPositronDeDxPEtaR0->ProjectionX(Form("R0PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoPositronR0DeDx[i][j],kBlack);
	meanPositronR0[i][j]=Mean;
	widthPositronR0[i][j]=Width;
	fhistoMeanPositronR0->SetBinContent(i+1,j+1, meanPositronR0[i][j]);	
	fhistoWidthPositronR0->SetBinContent(i+1,j+1, widthPositronR0[i][j]);	


	//	cout<< " i,j,R1::" << i << " "<< j<<endl;
	histoPositronR1DeDx[i][j] = new TH1D(Form("R1PosiSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoPositronR1DeDx[i][j] =
	  (TH1D*)histoPositronDeDxPEtaR1->ProjectionX(Form("R1PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoPositronR1DeDx[i][j],kBlack);
	meanPositronR1[i][j]=Mean;
	widthPositronR1[i][j]=Width;
	fhistoMeanPositronR1->SetBinContent(i+1,j+1, meanPositronR1[i][j]);	
	fhistoWidthPositronR1->SetBinContent(i+1,j+1, widthPositronR1[i][j]);	

	//	cout<< " i,j,R2::" << i << " "<< j<<endl;
	histoPositronR2DeDx[i][j] = new TH1D(Form("R2PosiSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoPositronR2DeDx[i][j] =
	  (TH1D*)histoPositronDeDxPEtaR2->ProjectionX(Form("R2PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoPositronR2DeDx[i][j],kBlack);
	meanPositronR2[i][j]=Mean;
	widthPositronR2[i][j]=Width;
	fhistoMeanPositronR2->SetBinContent(i+1,j+1, meanPositronR2[i][j]);	
	fhistoWidthPositronR2->SetBinContent(i+1,j+1, widthPositronR2[i][j]);	

	//	cout<< " i,j,R3::" << i << " "<< j<<endl;
	histoPositronR3DeDx[i][j] = new TH1D(Form("R3PosiSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoPositronR3DeDx[i][j] =
	  (TH1D*)histoPositronDeDxPEtaR3->ProjectionX(Form("R3PosiSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoPositronR3DeDx[i][j],kBlack);
	meanPositronR3[i][j]=Mean;
	widthPositronR3[i][j]=Width;
	fhistoMeanPositronR3->SetBinContent(i+1,j+1, meanPositronR3[i][j]);	
	fhistoWidthPositronR3->SetBinContent(i+1,j+1, widthPositronR3[i][j]);	


	//	cout<< " i,j,R0::" << i << " "<< j<<endl;
	histoElectronR0DeDx[i][j] = new TH1D(Form("R0ElecSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoElectronR0DeDx[i][j] =
	  (TH1D*)histoElectronDeDxPEtaR0->ProjectionX(Form("R0ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoElectronR0DeDx[i][j],kBlack);
	meanElectronR0[i][j]=Mean;
	widthElectronR0[i][j]=Width;
	fhistoMeanElectronR0->SetBinContent(i+1,j+1, meanElectronR0[i][j]);	
	fhistoWidthElectronR0->SetBinContent(i+1,j+1, widthElectronR0[i][j]);	

	//	cout<< " i,j,R1::" << i << " "<< j<<endl;
	histoElectronR1DeDx[i][j] = new TH1D(Form("R1ElecSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoElectronR1DeDx[i][j] =
	  (TH1D*)histoElectronDeDxPEtaR1->ProjectionX(Form("R1ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoElectronR1DeDx[i][j],kBlack);
	meanElectronR1[i][j] = Mean;
	widthElectronR1[i][j] = Width;
	fhistoMeanElectronR1->SetBinContent(i+1,j+1, meanElectronR1[i][j]);	
	fhistoWidthElectronR1->SetBinContent(i+1,j+1, widthElectronR1[i][j]);	

	//	cout<< " i,j,R2::" << i << " "<< j<<endl;
	histoElectronR2DeDx[i][j] = new TH1D(Form("R2ElecSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoElectronR2DeDx[i][j] =
	  (TH1D*)histoElectronDeDxPEtaR2->ProjectionX(Form("R2ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoElectronR2DeDx[i][j],kBlack);
	meanElectronR2[i][j]=Mean;
	widthElectronR2[i][j]=Width;
	fhistoMeanElectronR2->SetBinContent(i+1,j+1, meanElectronR2[i][j]);	
	fhistoWidthElectronR2->SetBinContent(i+1,j+1, widthElectronR2[i][j]);	

	//	cout<< " i,j,R3::" << i << " "<< j<<endl;
	histoElectronR3DeDx[i][j] = new TH1D(Form("R3ElecSigdEdxP%dEta%d",i,j),"",100,-5.,5.);
	histoElectronR3DeDx[i][j] =
	  (TH1D*)histoElectronDeDxPEtaR3->ProjectionX(Form("R3ElecSigdEdxP%dEta%d",i,j),j+1,j+1,i+1,i+1);
	FitSignal(histoElectronR3DeDx[i][j],kBlack);
	meanElectronR3[i][j]=Mean;
	widthElectronR3[i][j]=Width;
	fhistoMeanElectronR3->SetBinContent(i+1,j+1, meanElectronR3[i][j]);	
	fhistoWidthElectronR3->SetBinContent(i+1,j+1, widthElectronR3[i][j]);	
      }
    }

    // Ploting Electron

    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R0ElectronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R0ElectronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR0DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR0DeDx[i][j]->GetMaximum() , 1, kGray+2, 7);
	DrawGammaLines(meanElectronR0[i][j], meanElectronR0[i][j],.0,1.3*histoElectronR0DeDx[i][j]->GetMaximum() , 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }

    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R1ElectronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R1ElectronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR1DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR1DeDx[i][j]->GetMaximum() , 1, kGray+2, 7);
	DrawGammaLines(meanElectronR1[i][j], meanElectronR1[i][j],.0,1.3*histoElectronR1DeDx[i][j]->GetMaximum() , 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }
    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R2ElectronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R2ElectronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR2DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR2DeDx[i][j]->GetMaximum() , 1, kGray+2, 7);
	DrawGammaLines(meanElectronR2[i][j], meanElectronR2[i][j],.0,1.3*histoElectronR2DeDx[i][j]->GetMaximum() , 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }
    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R3ElectronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R3ElectronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR3DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR3DeDx[i][j]->GetMaximum(), 1, kGray+2, 7);
	DrawGammaLines(meanElectronR3[i][j], meanElectronR3[i][j],.0,1.3*histoElectronR3DeDx[i][j]->GetMaximum(), 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }



    // Ploting Positron

    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R0PositronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R0PositronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR0DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR0DeDx[i][j]->GetMaximum() , 1, kGray+2, 7);
	DrawGammaLines(meanElectronR0[i][j], meanElectronR0[i][j],.0,1.3*histoElectronR0DeDx[i][j]->GetMaximum() , 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }

    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R1PositronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R1PositronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR1DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR1DeDx[i][j]->GetMaximum() , 1, kGray+2, 7);
	DrawGammaLines(meanElectronR1[i][j], meanElectronR1[i][j],.0,1.3*histoElectronR1DeDx[i][j]->GetMaximum() , 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }
    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R2PositronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R2PositronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR2DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR2DeDx[i][j]->GetMaximum() , 1, kGray+2, 7);
	DrawGammaLines(meanElectronR2[i][j], meanElectronR2[i][j],.0,1.3*histoElectronR2DeDx[i][j]->GetMaximum() , 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }
    for(Int_t i=0;i<nPBins;i++){
      TString nameCanvas=Form("R3PositronSigmaDeDxForPbin%d",i);
      TCanvas *canvasDataSigmaDeDx          = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size
      DrawGammaCanvasSettings( canvasDataSigmaDeDx, 0, 0, 0, 0);
      canvasDataSigmaDeDx->cd();
      TString namePad=Form("R3PositronSigmaDeDxForPbin%d",i);
      TPad * padDataSigmaDeDx               = new TPad(namePad.Data(),"",-0.0,0.0,1.,1.,0);   // gives the size of the histo areas
      DrawGammaPadSettings( padDataSigmaDeDx, 0, 0, 0, 0);
      padDataSigmaDeDx->Divide(5,5,0.0,0.0);
      padDataSigmaDeDx->Draw();
      for(Int_t j=0;j<nEtaBins;j++){
	padDataSigmaDeDx->cd(j+1);
	histoElectronR3DeDx[i][j]->DrawCopy();
	DrawGammaLines(0.0, 0,.0,1.3*histoElectronR3DeDx[i][j]->GetMaximum(), 1, kGray+2, 7);
	DrawGammaLines(meanElectronR3[i][j], meanElectronR3[i][j],.0,1.3*histoElectronR3DeDx[i][j]->GetMaximum(), 1, kBlue+2, 2);
      }
      canvasDataSigmaDeDx->Print(Form("./DeDxMaps%s/Maps%s.pdf",optionMC.Data(),nameCanvas.Data()));
    }

    TFile outFileMonitoring(Form("./DeDxMaps%s/MonitoringDeDxMaps_%s.root",optionMC.Data(),fCutSelectionRead.Data()) ,"RECREATE");
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoPositronR0DeDx[i][j]->Write();
      }
    }
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoPositronR1DeDx[i][j]->Write();
      }
    }
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoPositronR2DeDx[i][j]->Write();
     }
    }
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoPositronR3DeDx[i][j]->Write();
      }
    }

    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoElectronR0DeDx[i][j]->Write();
    }
    }
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoElectronR1DeDx[i][j]->Write();
    }
    }
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoElectronR2DeDx[i][j]->Write();
    }
    }
    for(Int_t i=0;i<nPBins;i++){
      for(Int_t j=0;j<nEtaBins;j++){
	histoElectronR3DeDx[i][j]->Write();
      }
    }
    outFileMonitoring.Close();
    TFile outFileMaps(Form("./DeDxMaps%s/DeDxMaps_%s.root",optionMC.Data(),fCutSelectionRead.Data()) ,"RECREATE");
 
    fhistoMeanPositronR0->Write();
    fhistoMeanPositronR1->Write();
    fhistoMeanPositronR2->Write();
    fhistoMeanPositronR3->Write();


    fhistoMeanElectronR0->Write();
    fhistoMeanElectronR1->Write();
    fhistoMeanElectronR2->Write();
    fhistoMeanElectronR3->Write();

    fhistoWidthPositronR0->Write();
    fhistoWidthPositronR1->Write();
    fhistoWidthPositronR2->Write();
    fhistoWidthPositronR3->Write();


    fhistoWidthElectronR0->Write();
    fhistoWidthElectronR1->Write();
    fhistoWidthElectronR2->Write();
    fhistoWidthElectronR3->Write();

    outFileMaps.Close();
 }

//___________________________________________________
void FitSignal(TH1D* hSig,Color_t col){
  fGaus= NULL;
  fGaus = new TF1("f1","gaus",-4,5);
  fGaus->SetParameters(hSig->GetMaximum(),0.0,1);
  //  fGaus = new TF1("Gaussian","[0]*(exp(-0.5*((x-[1])/[2])^2)",-4,5);

  Int_t binmax = hSig->GetMaximumBin();
  Double_t x = hSig->GetXaxis()->GetBinCenter(binmax);
  fGaus->SetLineColor(col);
  fGaus->SetParameter(0,hSig->GetMaximum());
  fGaus->SetParameter(1,0.);
  if(  hSig->GetEntries()<40){  
    Mean=0.;
    Width=1.;
  }else if (hSig->GetEntries() >=40 && hSig->GetEntries()<150){
    hSig->Rebin(4);
    hSig->Fit(fGaus,"QR+","",-3.5,+3.5);
    Double_t xMean=fGaus->GetParameter(1);
    fGaus->SetLineColor(kRed); 
    hSig->Fit(fGaus,"QR+","",xMean-2.,xMean+2.);
    Mean  = fGaus->GetParameter(1);
    Width = fGaus->GetParameter(2);
    Chi2 = fGaus->GetChisquare()/fGaus->GetNDF();
    if (Width > 2){
      cout << Mean << " " << Width << "  " << hSig->GetEntries() << " " << hSig->GetMaximum()<< endl;
    }
  }else{
    hSig->Fit(fGaus,"QR+","",-3.5,+3.5);
    Double_t xMean=fGaus->GetParameter(1);
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
