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

TH2D *h2_asym_dphi;
TH2D *h2_asym_dphi_pt[11];
TH2D *h2_asym_geom;
TH2D *h2_asym_arth;

Color_t COLORS[11]      = { kRed+1, kOrange+7, kYellow+1, kSpring+3, kGreen, 
                            kTeal-6, kCyan-7, kAzure-2, kBlue, kViolet-5, 
                            kMagenta};
Double_t pt[11]         = { 0.1,  0.20, 0.40, 0.60, 0.80, 
                            0.90, 2.0,  1.25, 1.5,  2.0,  
                            3.5};
Double_t gammaFac[11]   = { 1.02, 1.1,  1.3,  1.6,  1.90, 
                            2.01, 2.2,  2.5,  2.8,  4.0,  
                            8.0};

char saythis[500];
char saythis2[500];

void DecayExample(){
  
  
  Int_t nEvts = 1000000;
  Double_t mpi = 0.139;
  Double_t mom = 1.0;
  Double_t E   = sqrt(mom*mom + mpi*mpi);
  TRandom3 randy;

  TGenPhaseSpace event;
  TLorentzVector pion(0.0, 0.0, mom, mpi);  
  Double_t masses[2] = {0.0, 0.0}; //here are the masses of the decay particle (ie. 0 for photons!)
  
  
  sprintf(saythis2,"#pi^{0} #rightarrow 2(m=%2.3f GeV/c^{2});asymmetry (0=symmetric);opening angle [deg]",masses[0]);
  h2_asym_dphi = new TH2D("h2_asym_dphi",saythis2, 500,0,1, 360,0,180);  
  h2_asym_geom = new TH2D("h2_asym_geom","h2_asym_geom;asymmetry (0=symmetric);geometric mean [GeV]", 500,0,1, 500,0,20);
  h2_asym_arth = new TH2D("h2_asym_arth","h2_asym_arth;asymmetry (0=symmetric);arithmetic mean [GeV]", 500,0,1, 500,0,20);
  
  TCanvas *c3 = new TCanvas("c3","c3",1200,750);
  c3->SetLeftMargin(0.08);
  c3->SetRightMargin(0.15);
  TLegend *leg1 = new TLegend(0.85,0.10,0.998,0.90);
  leg1->SetFillColor(0);
  
  for(Int_t i=0; i<11; i++){ // you can remove this i loop; I was doing some specific thing. 

    //mom = pt[i];
    //E   = sqrt(mom*mom + mpi*mpi);
    E   = mpi*gammaFac[i];
    mom = sqrt(E*E-mpi*mpi);
    sprintf(saythis,"h2_asym_dphi_pt_%d",i);
    sprintf(saythis2,"#pi^{0} #rightarrow 2(m=%2.3f GeV/c^{2});asymmetry (0=symmetric);opening angle [deg]",masses[0]);
    h2_asym_dphi_pt[i] = new TH2D(saythis,saythis2, 500,0,1, 360,0,180);  
    h2_asym_dphi_pt[i]->SetMarkerStyle(20);
    h2_asym_dphi_pt[i]->SetMarkerSize(1.5);
    h2_asym_dphi_pt[i]->SetMarkerColor(COLORS[i]);
    sprintf(saythis,"#gamma_{#pi^{0}}=%1.2f, #gamma_{e}=%1.2f",E/mpi,0.5*E/masses[0]);    
    leg1->AddEntry(h2_asym_dphi_pt[i]->Clone(),saythis,"P");
    h2_asym_dphi_pt[i]->SetMarkerSize(0.75);


    for(Int_t n=0; n<nEvts; n++){ // this is the important loop (nEvents)
      
      E   = mpi*gammaFac[i];
      E   = mpi*(1.01+randy.Rndm()*8.099);

      mom = sqrt(E*E-mpi*mpi);
      TLorentzVector pion(0.0, 0.0, mom, mpi); // in this simple case, I don't care about eta,phi
      event.SetDecay(pion, 2, masses); // ie. "set the pion to decay to 2 daughters with masses[0], masses[1]

      Double_t weight = event.Generate();
      TLorentzVector *p1    = event.GetDecay(0); // these are my daughters !! 
      TLorentzVector *p2    = event.GetDecay(1); // these are my daughters !! 

      Double_t pt1  = p1->P(); // grab whatever you want here; yes, I'm substituting P for pT 
      Double_t pt2  = p2->P(); // (although I think it's the same since eta=0 )

      h2_asym_dphi_pt[i]->Fill(fabs(pt1-pt2)/(pt1+pt2), 57.3*p1->Angle(p2->Vect()));
      h2_asym_dphi      ->Fill(fabs(pt1-pt2)/(pt1+pt2), 57.3*p1->Angle(p2->Vect()));
      h2_asym_geom      ->Fill(fabs(pt1-pt2)/(pt1+pt2), sqrt(pt1*pt2));
      h2_asym_arth      ->Fill(fabs(pt1-pt2)/(pt1+pt2), 0.5*(pt1+pt2));

    }

    c3->cd();
    if(i==0) h2_asym_dphi_pt[i]->DrawClone();
    else     h2_asym_dphi_pt[i]->DrawClone("same");
    h2_asym_dphi_pt[i]->SetMarkerStyle(1);
    h2_asym_dphi_pt[i]->SetMarkerColor(1);


  }//i loop (pt bins)

  c3->cd();
  leg1->Draw();

  TCanvas *c1 = new TCanvas("c1","c1",1200,750);
  c1->cd()->SetLogz();
  h2_asym_dphi->Draw("colz");
  for(Int_t i=0; i<11; i++)
    h2_asym_dphi_pt[i]->DrawClone("same");


  /*
  TCanvas *c2 = new TCanvas("c2","c2",1050,750);
  //h2_asym_geom->SetMarkerStyle(20);
  //h2_asym_arth->SetMarkerStyle(20);;
  h2_asym_geom->SetMarkerColor(2);
  h2_asym_arth->SetMarkerColor(4);;
  h2_asym_geom->Draw();
  h2_asym_arth->Draw("same");
  */  
  
}