/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/


#include <iostream>
#include <vector>
#include <fstream>

#include <TApplication.h>
#include <TArrow.h>
#include <TASImage.h>
#include <TAxis.h>
#include <TAttAxis.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TFrame.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TList.h>
#include <TMath.h>
#include <TMarker.h>
#include <TMinuit.h>
#include <TMultiGraph.h>
#include <TObject.h>
#include <TPaletteAxis.h>
#include <TPaveLabel.h>
#include <TPostScript.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TVirtualFitter.h>
#include <TPoint.h>

#include "../CommonHeaders/PlottingGammaConversionHistos.h"

struct pi0cand{
    pi0cand(){
      cluster_NClusters.clear();
      cluster_phi.clear();
      cluster_eta.clear();
      cluster_cellNCells.clear();
      cluster_cellAbsId.clear();
      cluster_cellNSupMod.clear();
      cluster_cellNPhi.clear();
      cluster_cellNEta.clear();
      cluster_cellEnergy.clear();

      pi0cand_openAngle.clear();
      pi0cand_pT.clear();
      pi0cand_clusterInvMass.clear();
    }

    std::vector<Int_t> cluster_NClusters;
    std::vector<Float_t> cluster_phi;
    std::vector<Float_t> cluster_eta;
    std::vector<Int_t> cluster_cellAbsId;
    std::vector<Int_t> cluster_cellNSupMod;
    std::vector<Int_t> cluster_cellNPhi;
    std::vector<Int_t> cluster_cellNEta;
    std::vector<Int_t> cluster_cellNCells;
    std::vector<Float_t> cluster_cellEnergy;

    std::vector<Float_t> pi0cand_openAngle;
    std::vector<Float_t> pi0cand_pT;
    std::vector<Float_t> pi0cand_clusterInvMass;
};

struct sepClusters{
    Float_t lowX;
    Float_t upX;
    Float_t lowY;
    Float_t upY;
};

void DetermineGlobalColRow(Int_t nSM, Int_t nPhi, Int_t nEta, Int_t &phi, Int_t &eta);

void GetMinMaxBin(TH1* hist, Int_t &iMin, Int_t &iMax);
void GetMinMaxBin(TH2* hist2D, Int_t &iMin, Int_t &iMax);
void GetMinMaxBinY(TH2* hist2D, Int_t &iMin, Int_t &iMax);
void SetXRange(TH2* hist, Int_t min, Int_t max);
void SetYRange(TH2* hist, Int_t min, Int_t max);
void SetZMinMaxTH2(TH2* hist, Int_t minX, Int_t maxX, Int_t minY, Int_t maxY, Bool_t setZero = kFALSE);
void DrawPeriodQAHistoTH2(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
                    TH2* fHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset);
void SaveCanvasAndWriteHistogram(TCanvas* canvas, TObject* hist, TString output);

void Grid_PlotGammaCaloDebug(TString filePath = "/home/daniel/data/work/debugOutput.txt", TString outputDir = "DebugPlots", Bool_t plotFullEMCal = kFALSE, Int_t debugOption = 0)
{
  gROOT->Reset();
  StyleSettingsThesis();
  SetPlotStyle();

  pi0cand cand;
  Int_t nEvent = 0;

  cout << Form("Reading from %s...", filePath.Data()) << endl;
  fstream file;
  TString fVar;
  file.open(filePath.Data(), ios::in);
  if(file.good())
  {
      file.seekg(0L, ios::beg);
      Int_t iClus = 0;
      while(!file.eof())
      {
          file >> fVar;
          if(fVar.Sizeof()>1)
          {
            if(fVar.BeginsWith("--pi0cand")){
                if(debugOption>1) cout << "--------------------------" << endl;
                file >> fVar;
                file >> fVar;
                cand.pi0cand_openAngle.push_back(fVar.Atof());
                file >> fVar;
                file >> fVar;
                cand.pi0cand_pT.push_back(fVar.Atof());
                file >> fVar;
                file >> fVar;
                cand.pi0cand_clusterInvMass.push_back(fVar.Atof());
                if(debugOption>1) cout << "pi0cand - opening angle: '" << cand.pi0cand_openAngle.at(cand.pi0cand_openAngle.size()-1) << "' - pT: '" << cand.pi0cand_pT.at(cand.pi0cand_openAngle.size()-1) << "'" << endl;
            }else if(fVar.BeginsWith("--cluster")){
                if(debugOption>1) cout << "cluster" << iClus;
                file >> fVar;
                Int_t nCells = 0;
                do{
                    cand.cluster_cellAbsId.push_back(fVar.Atoi());
                    file >> fVar;
                    cand.cluster_cellNSupMod.push_back(fVar.Atoi());
                    file >> fVar;
                    cand.cluster_cellNPhi.push_back(fVar.Atoi());
                    file >> fVar;
                    cand.cluster_cellNEta.push_back(fVar.Atoi());
                    file >> fVar;
                    cand.cluster_cellEnergy.push_back(fVar.Atof());
                    file >> fVar;
                    if(debugOption>1) cout << ".";
                    nCells++;
                }while(!fVar.Contains("phi"));
                if(debugOption>1) cout << " - " << nCells << " cells";
                cand.cluster_cellNCells.push_back(nCells);
                file >> fVar;
                cand.cluster_phi.push_back(fVar.Atof());
                file >> fVar;
                file >> fVar;
                cand.cluster_eta.push_back(fVar.Atof());
                iClus++;
                if(debugOption>1) cout << endl;
            } else if(fVar.BeginsWith("--event")){
              file >> fVar;
              file >> fVar;
              if(debugOption>1) cout << "--------------------------" << endl;
              if(debugOption>1) cout << "#event: " << nEvent++ << endl;
              if(debugOption>1) cout << "#clusters in event: " << fVar.Atoi() << endl;
              cand.pi0cand_openAngle.push_back(0.);
              cand.pi0cand_pT.push_back(0.);
              cand.pi0cand_clusterInvMass.push_back(0.);
            } else if(fVar.BeginsWith("---")){
              cand.cluster_NClusters.push_back(iClus);
              iClus = 0;
            }

          }
      }
  }
  cout << "reading done!" << endl;
  file.close();

  TCanvas* canvas         = new TCanvas("canvas","",10,10,750,750);  // gives the page size
  Double_t leftMargin     = 0.15;
  Double_t rightMargin    = 0.15;
  Double_t topMargin      = 0.15;
  Double_t bottomMargin   = 0.15;

  gSystem->Exec("mkdir -p "+outputDir);

//  const Int_t nEmcalEtaBins             = 96;
//  const Int_t nEmcalPhiBins             = 124;
//  Float_t EmcalEtaBins[nEmcalEtaBins+1] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
//  Float_t EmcalPhiBins[nEmcalPhiBins+1] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};
  const Int_t nEmcalColBins             = 96;
  const Int_t nEmcalRowBins             = 120;

  Int_t iCell = 0;
  Int_t iCellMax = 0;

  Int_t iClus = 0;
  Int_t iClusMax = 0;

  TLine *line             = new TLine();
  line->SetLineStyle(2);
  line->SetLineWidth(0.5);
  line->SetLineColor(1);

  TLine *lineClusters             = new TLine();
  lineClusters->SetLineStyle(9);
  lineClusters->SetLineWidth(6);
  lineClusters->SetLineColor(1);

  TLine *borderSM             = new TLine();
  borderSM->SetLineStyle(9);
  borderSM->SetLineWidth(3);
  borderSM->SetLineColor(1);

  TLine *borderEMCal             = new TLine();
  borderEMCal->SetLineStyle(3);
  borderEMCal->SetLineWidth(4);
  borderEMCal->SetLineColor(2);

  TMarker *point             = new TMarker();
  point->SetMarkerStyle(34);
  point->SetMarkerSize(2);
  point->SetMarkerColor(1);

  std::vector<sepClusters> vecSepClus;
  //std::vector<sepClusters> vecCenterClus;

  for(Int_t iEvent = 0; iEvent<(Int_t)cand.cluster_NClusters.size(); iEvent++){
    //vecCenterClus.clear();
    vecSepClus.clear();
    if(debugOption>1) cout << "--------" << endl;
    if(debugOption>1) cout << "#event: " << iEvent << endl;
    if(debugOption>1) cout << cand.cluster_NClusters.at(iEvent) << endl;
    if((Int_t)cand.cluster_NClusters.at(iEvent)==0) continue;
    iClusMax += (Int_t)cand.cluster_NClusters.at(iEvent);
    TH2F* fHistClusterEtavsPhiAfterQA     = new TH2F(Form("clusters_%i",iEvent),"EtaPhi_afterClusterQA",nEmcalRowBins,0,nEmcalRowBins,nEmcalColBins,0,nEmcalColBins);
    if((Int_t)cand.cluster_NClusters.at(iEvent)==1){
      for(; iClus<iClusMax; iClus++){
//        sepClusters centerTemp;
//        centerTemp.lowX = cand.cluster_phi.at(iClus);
//        centerTemp.lowY = cand.cluster_eta.at(iClus);
//        if(debugOption>1) cout << "clus" << iClus << " - " << centerTemp.lowX << "/" << centerTemp.lowY << endl;
//        vecCenterClus.push_back(centerTemp);

        iCellMax+=(Int_t)cand.cluster_cellNCells.at(iClus);
        if(debugOption>1) cout << iCell << "/" << iCellMax << endl;

        for(; iCell<iCellMax; iCell++){
          if(debugOption>1) cout << iCell << " - " << cand.cluster_cellAbsId.at(iCell) << endl;
          Int_t tempX, tempY;
          DetermineGlobalColRow(cand.cluster_cellNSupMod.at(iCell),cand.cluster_cellNPhi.at(iCell),cand.cluster_cellNEta.at(iCell),tempX,tempY);
          fHistClusterEtavsPhiAfterQA->Fill(tempX,tempY,cand.cluster_cellEnergy.at(iCell));
        }
      }
    }else{
      for(; iClus<iClusMax; iClus++){
//        sepClusters centerTemp;
//        centerTemp.lowX = cand.cluster_phi.at(iClus);
//        centerTemp.lowY = cand.cluster_eta.at(iClus);
//        if(debugOption>1) cout << "clus" << iClus << " - " << centerTemp.lowX << "/" << centerTemp.lowY << endl;
//        vecCenterClus.push_back(centerTemp);

        iCellMax+=(Int_t)cand.cluster_cellNCells.at(iClus);
        if(debugOption>1) cout << iCell << "/" << iCellMax  << endl;

        std::vector<Int_t> tempBins;
        for(; iCell<iCellMax; iCell++){
          if(debugOption>1) cout << iCell<< " - " << cand.cluster_cellAbsId.at(iCell) << endl;
          Int_t tempX, tempY;
          DetermineGlobalColRow(cand.cluster_cellNSupMod.at(iCell),cand.cluster_cellNPhi.at(iCell),cand.cluster_cellNEta.at(iCell),tempX,tempY);
          tempBins.push_back(fHistClusterEtavsPhiAfterQA->Fill(tempX,tempY,cand.cluster_cellEnergy.at(iCell)));
        }

        Int_t iCell_sec=iCell;
        Int_t iCellMax_sec=iCellMax;
        for(Int_t iClus2=iClus+1; iClus2<iClusMax; iClus2++){
          iCellMax_sec+=(Int_t)cand.cluster_cellNCells.at(iClus2);

          for(Int_t iCell2=iCell_sec; iCell2<iCellMax_sec; iCell2++){
            Int_t tempX, tempY;
            DetermineGlobalColRow(cand.cluster_cellNSupMod.at(iCell2),cand.cluster_cellNPhi.at(iCell2),cand.cluster_cellNEta.at(iCell2),tempX,tempY);
            Int_t tempBin2 = fHistClusterEtavsPhiAfterQA->FindBin(tempX,tempY);
            for(Int_t iB = 0; iB < (Int_t)tempBins.size(); iB++){
              Int_t diff = tempBins.at(iB)-tempBin2;
              sepClusters sepTemp;
              if(diff == 1){
                sepTemp.lowX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.lowY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinLowEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                vecSepClus.push_back(sepTemp);
              }else if(diff == -1){
                sepTemp.lowX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinLowEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.lowY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinLowEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinLowEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                vecSepClus.push_back(sepTemp);
              }else if (diff == 122){
                sepTemp.lowX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinLowEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.lowY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                vecSepClus.push_back(sepTemp);
              }else if (diff == -122){
                sepTemp.lowX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinLowEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.lowY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinLowEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upX = fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(tempBin2%(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                sepTemp.upY = fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinLowEdge(tempBin2/(fHistClusterEtavsPhiAfterQA->GetNbinsX()+2));
                vecSepClus.push_back(sepTemp);
              }
            }
          }
          iCell_sec+=(Int_t)cand.cluster_cellNCells.at(iClus2);
        }
        tempBins.clear();
      }
    }
    Int_t minB = 1;
    Int_t maxB = fHistClusterEtavsPhiAfterQA->GetNbinsX();
    Int_t minYB = 1;
    Int_t maxYB = fHistClusterEtavsPhiAfterQA->GetNbinsY();
    if(!plotFullEMCal){
      GetMinMaxBin(fHistClusterEtavsPhiAfterQA,minB,maxB);
      if(minB-2>=1) minB-=2;
      if(maxB+2<fHistClusterEtavsPhiAfterQA->GetNbinsX()) maxB+=2;
      GetMinMaxBinY(fHistClusterEtavsPhiAfterQA,minYB,maxYB);
      if(minYB-2>=1) minYB-=2;
      if(maxYB+2<fHistClusterEtavsPhiAfterQA->GetNbinsY()) maxYB+=2;
    }
    SetXRange(fHistClusterEtavsPhiAfterQA,minB,maxB);
    SetYRange(fHistClusterEtavsPhiAfterQA,minYB,maxYB);
    SetZMinMaxTH2(fHistClusterEtavsPhiAfterQA,minB,maxB,minYB,maxYB);

    fHistClusterEtavsPhiAfterQA->GetZaxis()->SetMoreLogLabels();
    fHistClusterEtavsPhiAfterQA->GetZaxis()->SetNoExponent();
    fHistClusterEtavsPhiAfterQA->GetZaxis()->SetTitle("GeV");
    fHistClusterEtavsPhiAfterQA->GetZaxis()->SetTitleOffset(1.1);
    DrawPeriodQAHistoTH2(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                         fHistClusterEtavsPhiAfterQA,"","col (0-119 in #phi)","row (0-95 in #eta)",1,1.4);

    if(!plotFullEMCal){
      for(Int_t iBx=minB; iBx<maxB; iBx++){
        if(minYB-1 < 1) line->DrawLine(fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(iBx), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(minYB),fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(iBx), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(maxYB));
        else line->DrawLine(fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(iBx), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(minYB-1),fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(iBx), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(maxYB));
      }
      for(Int_t iBy=minYB; iBy<maxYB; iBy++){
        if(minB-1 < 1) line->DrawLine(fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(minB), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(iBy),fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(maxB), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(iBy));
        else line->DrawLine(fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(minB-1), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(iBy),fHistClusterEtavsPhiAfterQA->GetXaxis()->GetBinUpEdge(maxB), fHistClusterEtavsPhiAfterQA->GetYaxis()->GetBinUpEdge(iBy));
      }
    }

//    for(Int_t iVec = 0; iVec<(Int_t) vecCenterClus.size(); iVec++){
//      sepClusters temp = vecCenterClus.at(iVec);
//      point->DrawMarker(temp.lowX,temp.lowY);
//    }
    for(Int_t iVec = 0; iVec<(Int_t) vecSepClus.size(); iVec++){
      sepClusters temp = vecSepClus.at(iVec);
      lineClusters->DrawLine(temp.lowX,
                             temp.lowY,
                             temp.upX,
                             temp.upY);
    }

    if(minB < 24 && maxB > 24) borderSM->DrawLine(24,minYB-1,24,maxYB);
    if(minB < 48 && maxB > 48) borderSM->DrawLine(48,minYB-1,48,maxYB);
    if(minB < 72 && maxB > 72) borderSM->DrawLine(72,minYB-1,72,maxYB);
    if(minB < 96 && maxB > 96) borderSM->DrawLine(96,minYB-1,96,maxYB);
    if(minYB < 48 && maxYB > 48) borderSM->DrawLine(minB-1,48,maxB,48);

    //cout << minB << ", " << maxB << ", " << minYB << ", " << maxYB << endl;
    if(minB<=1) borderEMCal->DrawLine(0,minYB-1,0,maxYB);
    if(maxB>=119) borderEMCal->DrawLine(120,minYB-1,120,maxYB);
    if(minYB<=1) borderEMCal->DrawLine(minB-1,0,maxB,0);
    if(maxYB>=95) borderEMCal->DrawLine(minB-1,96,maxB,96);

    SaveCanvasAndWriteHistogram(canvas, fHistClusterEtavsPhiAfterQA, Form("%s/output_%i.%s", outputDir.Data(), iEvent, "eps"));

    if(debugOption>1) cout << "end:" << iCell << "/" << iCellMax << endl;
    delete fHistClusterEtavsPhiAfterQA;
  }

  return;
}

void DetermineGlobalColRow(Int_t nSM, Int_t nPhi, Int_t nEta, Int_t &phi, Int_t &eta){

  //cout << "nSM: " << nSM << " - nPhi: " << nPhi << " - nEta: " << nEta << endl;
  phi = -1;
  eta = -1;

  if(nSM%2 == 1){
    eta = nEta + 48;
  }else{
    eta = nEta;
  }

  Int_t tempRest = nSM / 2;
  phi = nPhi + tempRest * 24;

  //cout << "eta: " << eta << " - phi: " << phi << endl;
  return;
}

void SetZMinMaxTH2(TH2* hist, Int_t minX, Int_t maxX, Int_t minY, Int_t maxY, Bool_t setZero){
    if(hist->GetEntries()<2) return;

    if(setZero){
      for(Int_t iX=minX; iX<=maxX; iX++){
        for(Int_t iY=minY; iY<=maxY; iY++){
          Double_t temp = hist->GetBinContent(iX,iY);
          if(temp<=0.) hist->SetBinContent(iX,iY,0.);
        }
      }
    }

    Double_t min = 0;
    Double_t max = 0;
    Bool_t bSet = kTRUE;

    for(Int_t iX=minX; iX<=maxX; iX++){
        for(Int_t iY=minY; iY<=maxY; iY++){
            Double_t temp = hist->GetBinContent(iX,iY);
            if(temp!=0.){
              if(bSet){
                min = temp;
                max = temp;
                bSet = kFALSE;
              }
              if(temp > max) max = temp;
              if(temp < min) min = temp;
            }
        }
    }

    hist->GetZaxis()->SetRangeUser(min,max);

    return;
}

void DrawPeriodQAHistoTH2(TCanvas* canvas,Double_t leftMargin,Double_t rightMargin,Double_t topMargin,Double_t bottomMargin,Bool_t logX, Bool_t logY, Bool_t logZ,
                    TH2* fHist, TString title,TString xTitle,TString yTitle, Double_t xOffset, Double_t yOffset)
{
    canvas->cd();
    canvas->SetLeftMargin(leftMargin);
    canvas->SetRightMargin(rightMargin);
    canvas->SetTopMargin(topMargin);
    canvas->SetBottomMargin(bottomMargin);
    //x-axis
    fHist->GetXaxis()->SetTitleOffset(xOffset);
    fHist->SetXTitle(xTitle.Data());
    fHist->GetXaxis()->SetLabelFont(42);
    fHist->GetXaxis()->SetTitleFont(62);
    fHist->GetXaxis()->SetTitleSize(0.04);
    fHist->GetXaxis()->SetLabelSize(0.035);
    //y-axis
    fHist->GetYaxis()->SetTitleOffset(yOffset);
    fHist->SetYTitle(yTitle.Data());
    fHist->GetYaxis()->SetLabelFont(42);
    fHist->GetYaxis()->SetTitleFont(62);
    fHist->GetYaxis()->SetTitleSize(0.04);
    fHist->GetYaxis()->SetLabelSize(0.035);
    fHist->GetYaxis()->SetDecimals();
    //z-axis
    fHist->GetZaxis()->SetLabelFont(42);
    fHist->GetZaxis()->SetTitleFont(62);
    fHist->GetZaxis()->SetTitleSize(0.04);
    fHist->GetZaxis()->SetLabelSize(0.035);
    //draw
    if(title.IsNull()) fHist->SetTitle("");
    else{
      if(topMargin<0.06) canvas->SetTopMargin(0.06);
      fHist->SetTitle(title.Data());
    }
    fHist->DrawCopy("colz");

    canvas->SetLogx(logX); canvas->SetLogy(logY); canvas->SetLogz(logZ);
    return;
}

void SaveCanvasAndWriteHistogram(TCanvas* canvas, TObject* hist, TString output){
    canvas->Update();
    canvas->SaveAs(output.Data());
    canvas->Clear();
    return;
}

void SetXRange(TH2* hist, Int_t min, Int_t max){

    if( max < min ) cout << "ERROR: SetXRangeTH2, max (" << max << ") < min (" << min << ")" << endl;
    if( min < 1 ) min = 1;
    if( max > hist->GetNbinsX() ) max = hist->GetNbinsX();

    hist->GetXaxis()->SetRange(min,max);

    return;
}

void SetYRange(TH2* hist, Int_t min, Int_t max){

    if( max < min ) cout << "ERROR: SetYRangeTH2, max (" << max << ") < min (" << min << ")" << endl;
    if( min < 1 ) min = 1;
    if( max > hist->GetNbinsY() ) max = hist->GetNbinsY();

    hist->GetYaxis()->SetRange(min,max);

    return;
}

void GetMinMaxBin(TH1* hist, Int_t &iMin, Int_t &iMax){
    if(!hist) {cout << "INFO: NULL pointer given to GetMinMaxBin!'" << endl; return;}
    iMin = 1;
    iMax = hist->GetNbinsX();
    Double_t min = hist->GetBinContent(iMin);
    Double_t max = hist->GetBinContent(iMax);

    while(min==0.){
        min = hist->GetBinContent(++iMin);
        if(iMin==iMax){
            iMin=1;
            break;
        }
    }

    while(max==0.){
        max = hist->GetBinContent(--iMax);
        if(iMax==1 && max==0.){
            iMax = hist->GetNbinsX();
            cout << "INFO: Empty histogram given to GetMinMaxBin: '" << hist->GetName() << "'" << endl;
            break;
        }
    }
    return;
}

void GetMinMaxBin(TH2* hist2D, Int_t &iMin, Int_t &iMax){
    if(!hist2D){cout << "INFO: NULL pointer given to GetMinMaxBin!'" << endl; return;}
    //if(hist2D->IsA()!=TH2::Class()){cout << Form("INFO: No TH2 given to GetMinMaxBin(TH2*): %s!'",hist2D->GetName()) << endl; return;}
    TH1* hist = (TH1*) hist2D->ProjectionX(Form("%s_project",hist2D->GetName()),1,hist2D->GetNbinsY());
    GetMinMaxBin(hist, iMin, iMax);
    delete hist;
    return;
}

void GetMinMaxBinY(TH2* hist2D, Int_t &iMin, Int_t &iMax){
    if(!hist2D){cout << "INFO: NULL pointer given to GetMinMaxBin!'" << endl; return;}
    //if(hist2D->IsA()!=TH2::Class()){cout << Form("INFO: No TH2 given to GetMinMaxBin(TH2*): %s!'",hist2D->GetName()) << endl; return;}
    TH1* hist = (TH1*) hist2D->ProjectionY(Form("%s_project",hist2D->GetName()),1,hist2D->GetNbinsX());
    GetMinMaxBin(hist, iMin, iMax);
    delete hist;
    return;
}
