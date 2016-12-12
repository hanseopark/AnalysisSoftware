/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/

using namespace std;

#include <Riostream.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <math.h>
#include <queue>
#include <algorithm>

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

#include "AliOADBContainer.h"
#include "AliEMCALGeometry.h"

void PrepareEMCalAcceptanceMap(TString str = "mod_acc"){

  TH1::AddDirectory(kFALSE);

  //********************
  // general definitions
  //********************
  TString pathOutput = "~/Desktop/out";
  gSystem->mkdir(pathOutput.Data());

  AliEMCALGeometry* rec = new AliEMCALGeometry("EMCAL_COMPLETE12SMV1"); // EMCAL_COMPLETE12SMV1

  Int_t nMaxCellsEMCAL  = 10*48*24;
  Int_t imod = -1;Int_t iTower = -1, iIphi = -1, iIeta = -1;
  Int_t icol = -1;Int_t irow = -1;

  //********************
  // histos
  //********************
  TH1S* output = new TH1S(str.Data(),str.Data(),nMaxCellsEMCAL,0,nMaxCellsEMCAL);

  TH2S* histCells[10];
  TH2S* histModules[10];
  for(Int_t i=0; i<10; i++){
    histCells[i] = new TH2S(Form("cellID - SM%i",i),Form("cellID - SM%i",i),48,0,48,24,0,24);
    histModules[i] = new TH2S(Form("moduleID - SM%i",i),Form("moduleID - SM%i",i),24,0,24,12,0,12);
  }

  //********************
  // processing
  //********************
  for(Int_t i=0; i<nMaxCellsEMCAL; i++){
    rec->GetCellIndex(i,imod,iTower,iIphi,iIeta);
    rec->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi,iIeta,irow,icol);

    //examples
    //cout << i << " - " << imod << ":" << irow << "/" << icol << endl;
    histCells[imod]->SetBinContent(histCells[imod]->FindBin(icol,irow),i);
    //cout << i << " - " << imod << " - " << iTower << ":" << iIphi << "/" << iIeta << endl;
    histModules[imod]->SetBinContent(histModules[imod]->FindBin(icol/2,irow/2),iTower);

    //output
    if(str.CompareTo("mod_acc") == 0){
      if(iTower%2 == 1) output->SetBinContent(i+1,1);
      else output->SetBinContent(i+1,0);
    }else if(str.BeginsWith("SM-")){
      TString temp = str;
      temp.Replace(0,3,"");
      Int_t SMnumber = temp.Atoi();
      if(imod == SMnumber){
        output->SetBinContent(i+1,1);
        histModules[imod]->SetBinContent(histModules[imod]->FindBin(icol/2,irow/2),iTower);
      }else{
        output->SetBinContent(i+1,0);
        histModules[imod]->SetBinContent(histModules[imod]->FindBin(icol/2,irow/2),0);
      }
    }else if(str.CompareTo("mod_acc-EMCAL") == 0){
      const Int_t nAcc = 116;
      Int_t accept[nAcc] = {
        13,14,15,16,17,18,19,20,21,22,
        25,26,27,28,29,30,31,32,33,34,
        37,38,41,42,45,46,49,50,53,54,
        57,58,61,62,65,66,69,70,
        109,110,111,112,113,114,115,116,117,118,
        121,122,123,124,125,126,127,128,129,130,
        140,141,151,152,164,165,
        169,170,171,172,173,174,175,176,177,178,
        181,182,183,184,185,186,187,188,189,190,
        229,230,231,232,233,234,235,236,237,238,
        241,242,243,244,245,246,247,248,249,250,
        261,262,273,274,285,286,253,254,265,266,
        277,278
      };
      Bool_t isFound = kFALSE;
      for(Int_t i=0; i<nAcc; i++){if(accept[i]==iTower) isFound = kTRUE;}
      if(isFound){
        output->SetBinContent(i+1,1);
        histModules[imod]->SetBinContent(histModules[imod]->FindBin(icol/2,irow/2),iTower);
      }else{
        output->SetBinContent(i+1,0);
        histModules[imod]->SetBinContent(histModules[imod]->FindBin(icol/2,irow/2),0);
      }
    }
  }

  //********************
  // write output file
  //********************
  TFile* fOutput = new TFile(Form("%s/cellID.root",pathOutput.Data()),"UPDATE");
  TDirectory *tmpDir;
  if(!fOutput->cd("example")) tmpDir = fOutput->mkdir("example");
  else tmpDir = fOutput->GetDirectory("example");

  fOutput->cd();
  output->Write(str.Data(),TObject::kOverwrite);

  tmpDir->cd();
  for(Int_t i=0; i<10; i++){
    histCells[i]->Write(Form("CellID_SM%i",i),TObject::kOverwrite);
    histModules[i]->Write(Form("ModuleID_SM%i",i),TObject::kOverwrite);
  }
  fOutput->Close();
  delete fOutput;

  TH1::AddDirectory(kTRUE);

  return;
}
