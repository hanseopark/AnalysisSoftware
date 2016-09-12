/****************************************************************************************************************************
 ******     provided by Gamma Conversion Group, PWGGA,                                                                  *****
 ******        Daniel Mühlheim, d.muehlheim@cern.ch                                                                     *****
 *****************************************************************************************************************************/
#include <TString.h>
#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TFrame.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TF1.h>
#include <TASImage.h>
#include <TMarker.h>
#include <TDatime.h>
#include <TArrow.h>
#include <TROOT.h>
#include <TObjString.h>
#include <TFile.h>
#include "TMath.h"


#include <Riostream.h>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <fstream>
#include <math.h>

#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

using namespace std;

Int_t GetPosInVec( std::vector<TString> vec, TString lookup){
  for(Int_t i=0; i<(Int_t)vec.size(); i++){
    if(vec.at(i).CompareTo(lookup) == 0) return i;
  }
  cout << "ERROR in GetPosInVec, did not find '" << lookup.Data() << "'" << endl;
  return -1;
}

void ComputeCorrelationFactors(
                                    TString fileInput = "input.txt",
                                    TString combMode = "triggers",
                                    TString meson = "",
                                    TString energy = "",
                                    TString mode = "",
                                    TString suffix = "eps"
                                   ){

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  TH1::AddDirectory(kFALSE);

  StyleSettingsThesis();
  SetPlotStyle();

  TString mesonPlot = "";
  if(meson.CompareTo("Pi0") == 0) mesonPlot = "#pi^{0}";
  else if(meson.CompareTo("Eta") == 0) mesonPlot = "#eta";
  else if(meson.CompareTo("Pi0EtaBinning") == 0 || meson.CompareTo("EtaToPi0") == 0) mesonPlot = "#eta/#pi^{0}";

  if(combMode.CompareTo("systems")==0) mode = "Systems";

  TString dateForOutput   = ReturnDateStringForOutput();
  TString collisionSystem = ReturnFullCollisionsSystem(energy);
  TString energyForOutput = energy;
  TString outputDir       = Form("%s/%s/ComputeCorrelationFactors_%s",suffix.Data(),dateForOutput.Data(), energyForOutput.Data());

  fstream fLog;
  fLog.open(Form("%s/%s_%s_%s_mode%s.log",outputDir.Data(),meson.Data(),combMode.Data(),energyForOutput.Data(),mode.Data()), ios::out);
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "Energy: " << energyForOutput.Data() << endl;
  fLog << collisionSystem.Data() << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

  cout << "creating outpur dir: " << outputDir.Data() << endl;
  gSystem->Exec("mkdir -p "+outputDir);
  gSystem->Exec("mkdir -p "+outputDir+"/sys");
  gSystem->Exec("mkdir -p "+outputDir+"/corrFactors");
  gSystem->Exec("mkdir -p "+outputDir+"/corrPlotting");
  gSystem->Exec(Form("cp %s %s/%s_%s_%s_mode%s.config", fileInput.Data(), outputDir.Data(), meson.Data(), combMode.Data(), energyForOutput.Data(), mode.Data()));

  cout << "running combination mode: " << combMode.Data() << endl;
  fLog << "running combination mode: " << combMode.Data() << endl;
  if(combMode.CompareTo("triggers") && combMode.CompareTo("systems")){ cout << "ERROR in line '" << __LINE__ << "' combMode can only be 'triggers' or 'systems', returning..."<< endl; return;}

  Int_t nRead = 0;
  vector<TString>  vecComb;
  vector<TString>  vecSysFiles;
  vector<Int_t>    vecNBins;
  vector<Double_t> vecUseMinPt;
  vector<Double_t> vecUseMaxPt;

  //---------------------------------------------------------------------------------------------------------------
  //--------------------------------- read config file
  //---------------------------------------------------------------------------------------------------------------

  cout << "reading input file" << endl;

  ifstream in(fileInput.Data());
  in >> nRead;
  cout << Form("number of %s available: ",combMode.Data()) << nRead << "\n" << endl;
  fLog << Form("number of %s available: ",combMode.Data()) << nRead << "\n" << endl;
  if(nRead<1){ cout << "ERROR in line '" << __LINE__ << "', returning..."<< endl; return;}

  for(Int_t iR=0; iR<nRead; iR++) {
    TString temp = "";
    TString tempSysFile = "";
    Int_t nBinsTemp = 0;
    Double_t minTemp = 0;
    Double_t maxTemp = 0;
    in >> temp >> nBinsTemp >> minTemp >> maxTemp >> tempSysFile;
    if(!temp.IsNull() && !tempSysFile.IsNull() && nBinsTemp > 0 && minTemp > 0 && maxTemp > 0){
      vecComb.push_back(temp);
      vecSysFiles.push_back(tempSysFile);
      vecNBins.push_back(nBinsTemp);
      vecUseMinPt.push_back(minTemp);
      vecUseMaxPt.push_back(maxTemp);
      TObjArray *rToken = tempSysFile.Tokenize("/");
      TObjString *rString = (TObjString*)rToken->At(rToken->GetLast());
      gSystem->Exec(Form("cp %s %s/sys/%s", tempSysFile.Data(), outputDir.Data(), ((TString)rString->GetString()).Data()));
      fLog << "read: " << temp.Data() << ", " << minTemp << ", " << maxTemp << ", " << tempSysFile.Data() << endl;
    }else{ cout << "ERROR in line '" << __LINE__ << "', returning..."<< endl; return;}
  }
  fLog << endl;

  vector<TString>** corrName = new vector<TString>*[nRead];
  vector<TString>** corr = new vector<TString>*[nRead];
  vector<Double_t> **corrFactorsBins = new vector<Double_t>*[nRead];
  vector<Double_t> **corrFactors = new vector<Double_t>*[nRead];
  for(Int_t iR=0; iR<nRead; iR++){
    corrName[iR] = new vector<TString>[nRead];
    corr[iR] = new vector<TString>[nRead];
    corrFactorsBins[iR] = new vector<Double_t>[nRead];
    corrFactors[iR] = new vector<Double_t>[nRead];
  }

  while(!in.eof()){
    TString tempA = ""; TString tempB = ""; TString tempC = ""; TString tempD = "";
    in >> tempA >> tempB >> tempC >> tempD;
    if(!tempA.IsNull() && !tempB.IsNull() && !tempC.IsNull() && !tempD.IsNull()){
      corrName[GetPosInVec(vecComb,tempA)][GetPosInVec(vecComb,tempB)].push_back(tempC);
      corr    [GetPosInVec(vecComb,tempA)][GetPosInVec(vecComb,tempB)].push_back(tempD);
    }
  }
  in.close();

  //---------------------------------------------------------------------------------------------------------------
  //--------------------------------- read detailed systematics for each input
  //---------------------------------------------------------------------------------------------------------------

  vector<TString>** ptSys = new vector<TString>*[nRead];
  for(Int_t iR=0; iR<nRead; iR++) ptSys[iR] = new vector<TString>[vecNBins.at(iR)+1];

  for(Int_t iR=0; iR<nRead; iR++){
    ifstream fileSysErr;
    fileSysErr.open(vecSysFiles.at(iR).Data(),ios_base::in);
    cout << "opening: " << vecSysFiles.at(iR).Data() << endl;
    fLog << "opening: " << vecSysFiles.at(iR).Data() << endl;
    Int_t counter = 0;
    string line;
    while (getline(fileSysErr, line) && counter < 100) {
      istringstream ss(line);
      TString temp="";
      TString tempBin="";
      Int_t counterColumn = 0;
      fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      fLog << "\t*** reading line " << counter << " ***" << endl;
      while(ss && counterColumn < 100){
        ss >> temp;
        if(counter > 0 && counterColumn == 0){
          tempBin = temp;
          if( !tempBin.IsNull() && (tempBin.Atof() < vecUseMinPt.at(iR) || tempBin.Atof() > vecUseMaxPt.at(iR)) ){
            cout << "INFO: skipping pT " << tempBin.Atof() << ", outside of specified interval of " << vecUseMinPt.at(iR) << "-" << vecUseMaxPt.at(iR) << endl;
            break;
          }
        }
        if( !(counter==0 && temp.CompareTo("bin")==0) && !temp.IsNull() && counter < (vecNBins.at(iR)+1)){
          ptSys[iR][counter].push_back(temp);
          fLog << temp.Data() << ", ";
          counterColumn++;
        }
      }

      if(counter == 0){
        ptSys[iR][counter++].push_back("TotalError");
        fLog << "TotalError";
        counterColumn++;
      }else if(!tempBin.IsNull() && (tempBin.Atof() >= vecUseMinPt.at(iR) && tempBin.Atof() <= vecUseMaxPt.at(iR))) counter++;
      fLog << "\n\t***" << counterColumn << " columns read ***" << endl;
      fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << endl;
    }
    vecNBins.at(iR) = --counter;
    cout << "\t---" << counter << " pT bins read ---" << endl;
    fLog << "\t---" << counter << " pT bins read ---" << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    if(counter < 1){
      cout << "ERROR: read " << counter << " pT-bins, returning!" << endl;
      fLog << "ERROR: read " << counter << " pT-bins, returning!" << endl;
      return;
    }
    fileSysErr.close();
  }

  //---------------------------------------------------------------------------------------------------------------
  //--------------------------------- calculate correlation factors
  //---------------------------------------------------------------------------------------------------------------

  Int_t textSizeLabelsPixel           = 900*0.04;

  TCanvas* canvasWeights = new TCanvas("canvasWeights","",200,10,1350,900);  // gives the page size
  DrawGammaCanvasSettings( canvasWeights, 0.08, 0.02, 0.035, 0.09);
  canvasWeights->SetLogx();

  Int_t plotErr = 0;
  TH2F * histo2DPi0Weights;
  Double_t minCorrYaxis = 0.45;
  if(combMode.CompareTo("systems") == 0) minCorrYaxis = 0.35;
  histo2DPi0Weights = new TH2F("histo2DPi0Weights","histo2DPi0Weights",11000,vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1),1000,minCorrYaxis,1.05);
  SetStyleHistoTH2ForGraphs(histo2DPi0Weights, "#it{p}_{T} (GeV/#it{c})","#rho_{A_B}",0.035,0.04, 0.035,0.04, 1.,1.);
  histo2DPi0Weights->GetXaxis()->SetMoreLogLabels();
  histo2DPi0Weights->GetXaxis()->SetLabelOffset(-0.01);
  histo2DPi0Weights->GetYaxis()->SetTitleOffset(0.8);
  histo2DPi0Weights->Draw("copy");

  TLatex *labelWeightsEnergy      = new TLatex(0.75,0.20,collisionSystem.Data());
  SetStyleTLatex( labelWeightsEnergy, 0.85*textSizeLabelsPixel,4);
  labelWeightsEnergy->SetTextFont(43);
  labelWeightsEnergy->Draw();
  TLatex *labelWeightsPi0 = 0x0;
  if(meson.CompareTo("Pi0EtaBinning") == 0 || meson.CompareTo("EtaToPi0") == 0) labelWeightsPi0 = new TLatex(0.75,0.16,Form("%s",mesonPlot.Data()));
  else labelWeightsPi0 = new TLatex(0.75,0.16,Form("%s #rightarrow #gamma#gamma",mesonPlot.Data()));
  SetStyleTLatex( labelWeightsPi0, 0.85*textSizeLabelsPixel,4);
  labelWeightsPi0->SetTextFont(43);
  labelWeightsPi0->Draw();

  Color_t colorTrigg      [10] = {kBlack, kGray+1, kRed+2, kBlue+2, kGreen+3, kCyan+2, kViolet, kMagenta+2,  kRed-2, kBlue-2};
  Marker_t markerTrigg    [10] = {20, 20, 21, 34, 29, 33, 21, 27, 28, 30 };
  Int_t iTrigg = 0;

  TLegend* legendWeights   = GetAndSetLegend2(0.11, 0.11, 0.42, 0.1+(0.035*8), 32);
  legendWeights->SetMargin(0.15);
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" << endl;

  fstream fCorr;
  fCorr.open(Form("%s/corrFactors/%s_%s_%s_mode%s.log",outputDir.Data(),meson.Data(),combMode.Data(),energyForOutput.Data(),mode.Data()), ios::out);
  TFile *fOutput = new TFile(Form("%s/%s.root",outputDir.Data(),energyForOutput.Data()),"UPDATE");

  for(Int_t iC=0; iC<nRead; iC++){
    for(Int_t iC2=0; iC2<nRead; iC2++){
      if(iC==iC2) continue;
      Int_t binC  = 1;
      Int_t binC2 = 1;
      Double_t pT  = -1;
      Double_t pT2 = -1;

      Double_t factor = 0;
      TString tempCorr = "";
      TString tempPlot = "";
      if(iC<iC2){
        tempCorr = Form("%s_%s-%s",vecComb.at(iC).Data(),vecComb.at(iC).Data(),vecComb.at(iC2).Data());
        tempPlot = Form("%s_%s",vecComb.at(iC).Data(),vecComb.at(iC2).Data());
        cout << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "correlation factors " << vecComb.at(iC) << "_" << vecComb.at(iC) << "-" << vecComb.at(iC2) << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "correlation factors " << vecComb.at(iC) << "_" << vecComb.at(iC) << "-" << vecComb.at(iC2) << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      }else{
        tempCorr = Form("%s_%s-%s",vecComb.at(iC).Data(),vecComb.at(iC2).Data(),vecComb.at(iC).Data());
        tempPlot = Form("%s_%s",vecComb.at(iC).Data(),vecComb.at(iC2).Data());
        cout << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "correlation factors " << vecComb.at(iC) << "_" << vecComb.at(iC2) << "-" << vecComb.at(iC) << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "correlation factors " << vecComb.at(iC) << "_" << vecComb.at(iC2) << "-" << vecComb.at(iC) << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
      }
      fCorr << endl;
      fCorr << tempCorr.Data() << endl;
      TH1D* histoCorr = new TH1D(Form("%s_%s_%s",mode.Data(),meson.Data(),tempCorr.Data()),Form("%s_%s_%s",mode.Data(),meson.Data(),tempCorr.Data()),1000,0,100);
      for(Int_t iBin=1; iBin<histoCorr->GetNbinsX()+1; iBin++) histoCorr->SetBinContent(iBin,0.);

      for(;; binC++, binC2++){
        if( binC == vecNBins.at(iC)+1 || binC2 == vecNBins.at(iC2)+1 ) break;

        pT  = ((TString)ptSys[iC][binC].at(0)).Atof();
        pT2 = ((TString)ptSys[iC2][binC2].at(0)).Atof();
        if( binC < vecNBins.at(iC) && binC2 < vecNBins.at(iC2) ){
          do{
            if(pT < pT2){
              pT = ((TString)ptSys[iC][++binC].at(0)).Atof();
            }else if(pT2 < pT){
              pT2 = ((TString)ptSys[iC2][++binC2].at(0)).Atof();
            }
          }while(pT!=pT2 && binC < vecNBins.at(iC) && binC2 < vecNBins.at(iC2));

        }
        if( pT!=pT2 ) break;

        cout << "\n\t---------- pT: " << pT << " ----------" << endl;
        fLog << "\n\t---------- pT: " << pT << " ----------" << endl;

        if( (Int_t) corrName[iC][iC2].size() > 0){
          Double_t totalErr = ((TString)ptSys[iC][binC].at(GetPosInVec(ptSys[iC][0],"TotalError"))).Atof();
          cout << vecComb.at(iC) << " totErr: " << totalErr << endl;
          fLog << vecComb.at(iC) << " totErr: " << totalErr << endl;
          Double_t unCorrErr = 0;
          for(Int_t iCorr=0; iCorr<(Int_t) corrName[iC][iC2].size(); iCorr++){
            TString tempName = corrName[iC][iC2].at(iCorr);
            TString tempD = corr[iC][iC2].at(iCorr);
            Double_t tempUnCorrErr = ((TString)ptSys[iC][binC].at(GetPosInVec(ptSys[iC][0],tempName))).Atof();
            Double_t tempErr = 0;
            if(tempD.Contains("%")){
              tempD.Resize(tempD.Length()-1);
              tempErr = tempUnCorrErr*tempD.Atof()/100.;
            }else tempErr = tempD.Atof();
            unCorrErr += tempErr*tempErr;
            cout << corrName[iC][iC2].at(iCorr) << "(" << tempUnCorrErr << "), unCorrErr: " << tempErr << endl;
            fLog << corrName[iC][iC2].at(iCorr) << "(" << tempUnCorrErr << "), unCorrErr: " << tempErr << endl;
          }
          factor = sqrt(totalErr*totalErr-unCorrErr)/totalErr;
        }
        corrFactorsBins[iC][iC2].push_back(pT);
        corrFactors[iC][iC2].push_back(factor);

        cout << "\n\t|-----------" << endl;
        cout << "\t|corrFactor:\t" << factor << "|" << endl;
        cout << "\t|------------" << endl;
        fLog << "\n\t|------------" << endl;
        fLog << "\n\tcorrFactor:\t" << factor << "|" << endl;
        fLog << "\t|------------" << endl;
        fCorr << pT << " \t-\t " << factor << endl;
        histoCorr->SetBinContent(histoCorr->FindBin(pT),factor);
        if(factor>0 && factor<=plotErr){
          cout << "\n\n\tWARNING: point out of range (" << factor <<  ") for generated plots, please adjust y-range of histograms!!!\n\n" << endl;
          fLog << "\n\n\tWARNING: point out of range (" << factor <<  ") for generated plots, please adjust y-range of histograms!!!\n\n" << endl;
          plotErr++;
        }
      }
      histoCorr->SetMarkerSize(2.);
      histoCorr->SetMarkerColor(colorTrigg[iTrigg]);
      histoCorr->SetMarkerStyle(markerTrigg[iTrigg++]);
      histoCorr->Draw("p,same");
      legendWeights->AddEntry(histoCorr,tempPlot.Data(),"p");
      histoCorr->Write(Form("%s_%s_%s",mode.Data(),meson.Data(),tempCorr.Data()),TObject::kOverwrite);
    }
  }
  legendWeights->Draw();
  canvasWeights->RedrawAxis();

  labelWeightsEnergy->Draw();
  labelWeightsPi0->Draw();

  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 1., 1.,0.1, kGray, 1);
  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 0.98, 0.98,0.1, kGray, 3);
  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 0.96, 0.96,0.1, kGray, 7);
  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 0.94, 0.94,0.1, kGray, 3);
  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 0.92, 0.92,0.1, kGray, 7);
  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 0.9, 0.9,0.1, kGray, 1);
  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 0.8, 0.8,0.1, kGray, 3);
  DrawGammaLines(vecUseMinPt.at(0),vecUseMaxPt.at(nRead-1), 0.7, 0.7,0.1, kGray, 7);

  canvasWeights->SaveAs(Form("%s/corrPlotting/%s_%s_%s_corrFactors.%s",outputDir.Data(),energyForOutput.Data(),mode.Data(),meson.Data(),suffix.Data()));

  if(plotErr>0){
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "\n\t|-----------" << endl;
    cout << "\t|total plotting errors:\t" << plotErr << "| -> adjust y-range in ComputeCorrelationFactors!" << endl;
    cout << "\t|------------" << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "\n\t|------------" << endl;
    fLog << "\n\ttotal plotting errors:\t" << plotErr << "| -> adjust y-range in ComputeCorrelationFactors!" << endl;
    fLog << "\t|------------" << endl;
  }

  fOutput->Write();
  fOutput->Close();
  fCorr.close();
  fLog.close();

  return;
}
