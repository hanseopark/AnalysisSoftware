
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

TObjArray* fEMCALBadChannelMap;

void ClusterQA_ReadEMCalBadCells(){

  TString pathFile = "/home/daniel/Desktop/EMCALBadChannels.root";
  TString pathOutput = "/home/daniel/Desktop/out/";

  gSystem->mkdir(pathOutput.Data());

  fstream fLog;
  AliEMCALGeometry* rec = new AliEMCALGeometry("EMCAL_FIRSTYEARV1"); // EMCAL_COMPLETE12SMV1

  TFile* fFile = new TFile(pathFile.Data(),"READ");
  TString nameMainDir="";
  if(fFile->IsZombie()){cout << "ERROR: File " << pathFile.Data() << " could not be openend! Returning..." << endl; return;}
  else{
      cout << "Processing file: " << pathFile.Data();
      TKey *key;
      TIter next(fFile->GetListOfKeys());
      while ((key=(TKey*)next())){
          cout << Form(" - found TopDir: %s",key->GetName());
          nameMainDir = key->GetName();
      }
      cout << endl;

      AliOADBContainer* cont = 0;
      fFile->GetObject("AliEMCALBadChannels",cont);
      for(Int_t i=0; i<cont->GetNumberOfEntries(); i++){
          TObjArray* obj = (TObjArray*)cont->GetObjectByIndex(i);
          TString str = obj->GetName();
          cout << obj->GetName() << endl;
          if(!str.Contains("10") && !str.Contains("BadChannelsbc_1")) continue;
          fLog.open(Form("%s/%s.txt",pathOutput.Data(),obj->GetName()), ios::out);
          for(Int_t j=0; j<obj->GetEntries(); j++){
            TH2I* h=(TH2I*)obj->FindObject(Form("EMCALBadChannelMap_Mod%i",j));
            //cout << j << ", " << h->GetNbinsX() << ", " << h->GetNbinsY() << endl;
            for(Int_t hX=0; hX<h->GetNbinsX(); hX++){
              for(Int_t hY=0; hY<h->GetNbinsY(); hY++){
                //cout << hX << ", " << hY << ", Content: " << h->GetBinContent(hX,hY) << endl;
                if(h->GetBinContent(hX,hY)>0){
                  fLog << rec->GetAbsCellIdFromCellIndexes(j,hY,hX) << endl;
                  cout << rec->GetAbsCellIdFromCellIndexes(j,hY,hX) << ", ";
                }
              }
            }
            cout << endl;
           }
         gSystem->Exec(Form("sort -k 1 %s/%s.txt > %s/%s_sorted.txt",pathOutput.Data(),obj->GetName(),pathOutput.Data(),obj->GetName()));
         fLog.close();
      }


  }

  return;
}
