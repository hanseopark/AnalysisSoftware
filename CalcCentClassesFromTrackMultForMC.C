#include <Riostream.h>
#include <fstream>
#include "TMath.h"
#include <stdio.h>
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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include <vector>
#include "CommonHeaders/PlottingGammaConversionHistos.h"
#include "CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "CommonHeaders/PlottingGammaConversionAdditional.h"

using std::cout;
using std::endl;

/* works only with c++11 i.e. root6 because of l.50
purpose: calculate multiplicity of central barrel tracks in centrality classes from data in order to adjust centrality classes accordingly for MC
input: data train output GammaConvV1 root file ( path to and names of files)
       photon and meson cut number
have to run ChangeStructureToStandard first    
 */


void CalcCentClassesFromTrackMultForMC(){

  StyleSettingsThesis();
  SetPlotStyle();

  std::vector<TString> fileNames{"/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_248_list1_train352.root","/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_250_list1_train352.root","/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_252_list1_train352.root","/home/meike/analysis/data/GridOutput/GammaConv/PbPb/ESD/LHC15o/GammaConvV1_254_list1_train352.root"};
  TString photonMesonCuts = "00200009247602008250404000_0652501500000000";  
  const Int_t maxX=3500;         // max x value for plotting (maximum possible:4000)
  const Int_t maxY=1000000;      // max y value for plotting
  const Bool_t drawLegend = kTRUE;
  const Bool_t debug = kFALSE;
  const Bool_t logX = kFALSE;
  const Int_t minBinContent = 100;        
  const Int_t maxBin = 2000;
  const TString outputDir = "eps";

  // read files
  std::vector<TList*> lists;
  for (unsigned int nFile=0; nFile<fileNames.size(); nFile++){
    TFile file(fileNames.at(nFile));
    TList *list = (TList*)file.Get("GammaConvV1");
    lists.push_back(list);
  }

  // vector containing all cutnumbers
  std::vector<TString> cutNumbers;
  for(Int_t i=0; i<9; i++){
    cutNumbers.push_back(Form("3%d%d10013_%s",i,i+1,photonMesonCuts.Data()));   // 0-45% in steps of 5%
  }
  for(Int_t i=0; i<9; i++){
    cutNumbers.push_back(Form("4%d%d10013_%s",i,i+1,photonMesonCuts.Data()));   // 45-90% in steps of 5%
  }

  // plotting
  TCanvas canvas("canvas","",1200,600);
  canvas.cd();
  canvas.SetLogy();  
  if(logX) canvas.SetLogx();

  gStyle->SetOptStat(0);
  if(!logX) canvas.SetRightMargin(0.04);
  canvas.SetLeftMargin(0.08);

  TLegend legend = TLegend(0.85,0.13,0.95,0.93);
  legend.SetLineColor(0);
  std::vector<TH1F*> histos;
  Int_t fileNo=0;
  cout << "found " << cutNumbers.size() << " cut selections: " << endl;
  for(unsigned int nCut=0; nCut<cutNumbers.size(); nCut++){
    cout << nCut << " " << cutNumbers.at(nCut).Data() << endl;
    TList *cutList=(TList*)(lists.at(fileNo)->FindObject(Form("Cut Number %s",cutNumbers.at(nCut).Data())));
    if(cutList){
      TList *subList=(TList*)cutList->FindObject(Form("%s ESD histograms",cutNumbers.at(nCut).Data()));
      TH1F *histoTracks=(TH1F*)subList->FindObject("GoodESDTracks");
      histoTracks->SetLineColor(46-nCut);
      histoTracks->SetTitle("");
      histoTracks->GetXaxis()->SetTitle("N_{Good Tracks}");
      histoTracks->GetYaxis()->SetTitle("#frac{dN}{dN_{Good Tracks}}");
      histoTracks->GetXaxis()->SetRangeUser(1,maxX);
      histoTracks->GetYaxis()->SetRangeUser(1,maxY);
      legend.AddEntry(histoTracks, GetCentralityString(cutNumbers.at(nCut)).Data());
      histoTracks->Draw("same");
      histos.push_back(histoTracks);
    } else{
      fileNo++;
      TList *cutList=(TList*)(lists.at(fileNo)->FindObject(Form("Cut Number %s",cutNumbers.at(nCut).Data())));
      TList *subList=(TList*)cutList->FindObject(Form("%s ESD histograms",cutNumbers.at(nCut).Data()));
      TH1F *histoTracks=(TH1F*)subList->FindObject("GoodESDTracks");
      histoTracks->SetLineColor(46-nCut);
      histoTracks->SetTitle("");
      histoTracks->GetXaxis()->SetTitle("N_{Good Tracks}");
      histoTracks->GetYaxis()->SetTitle("#frac{dN}{dN_{Good Tracks}}");
      histoTracks->GetXaxis()->SetRangeUser(1,maxX);
      histoTracks->GetYaxis()->SetRangeUser(1,maxY);
      legend.AddEntry(histoTracks, GetCentralityString(cutNumbers.at(nCut)).Data());
      histoTracks->Draw("same");
      histos.push_back(histoTracks);
    }
  }

  if(drawLegend) legend.Draw("same");

  // calculate intersections between histograms
  Int_t nBinsX = 0;  
  std::vector<Int_t> intersections;
  for (unsigned int nHisto=0; nHisto<histos.size()-1; nHisto++){  // compare every histogram to the next histogram, up to the second-to-last
    Int_t nBinsXtmp = histos.at(nHisto)->GetNbinsX();
    if(nHisto==0) nBinsX=nBinsXtmp;
    else {
      if(nBinsXtmp!=nBinsX) { // nBinsXtmp should be the same number for all histograms
	cout << "ERROR: different number of bins!" << endl;
	return;
      }
    }
    if (debug) cout << "histo " << nHisto << " " << nBinsX << " bins" << endl;
    for(Int_t nBin=maxBin; nBin>0; nBin--){  // go through all bins and compare bin content
      Double_t content = histos.at(nHisto)->GetBinContent(nBin);
      Double_t contentNext = histos.at(nHisto+1)->GetBinContent(nBin);
      if( (content > minBinContent) && (contentNext > minBinContent) && (content < contentNext) ){
	Double_t x = histos.at(nHisto)->GetXaxis()->GetBinUpEdge(nBin);
	if(debug) cout << "intersection at bin " << nBin << " corresponding to x = " << x << endl; 
	intersections.push_back(x);
	break;
      }
    }
  }

  cout << "found " << intersections.size() << " intersetions between neighbouring ones of the " << histos.size() << " histograms: " << endl;
  // plot intersections and write to stdout
  std::vector<TLine*> lines;
  for(unsigned int nLine=0; nLine<intersections.size(); nLine++){
    Int_t x = intersections.at(nLine);
    cout << nLine << "\t" << x << endl;
    TLine *line = new TLine(x,1,x,maxY); // x1, y1, x2, y2
    lines.push_back(line);
    line->SetLineStyle(2);
    line->Draw("same");
  }

  gSystem->Exec("mkdir -p "+outputDir);
  if(logX) canvas.SaveAs(Form("%s/TrackMultInCentClassesLog.eps",outputDir.Data()));
  else canvas.SaveAs(Form("%s/TrackMultInCentClasses.eps",outputDir.Data()));

  for(unsigned int i=0; i<lines.size(); i++){
   delete lines.at(i);
  }

}


