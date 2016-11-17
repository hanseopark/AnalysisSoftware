#include <Riostream.h>
#include <fstream>
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
#include "TLatex.h"
#include "TASImage.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TMarker.h"
#include "TGraphAsymmErrors.h" 
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"

void ChangeStructureToStandardCalo(TString nameInputFile, TString namefileOutput, TString nameInputList){

   TFile *fileInput = NULL;
   fileInput = new TFile(nameInputFile.Data());
   cout << fileInput << endl;
  
   TString autodetectedString   = AutoDetectMainTList(2 , fileInput);
   cout << "auto: " << autodetectedString << endl;
   
   TList *listInput =(TList*)fileInput->Get(nameInputList.Data());
   if (listInput == NULL){ 
      return;
   }   
   TFile *fileOutput = new TFile(namefileOutput,"RECREATE");
   TList *listOutput =(TList*)fileOutput->Get("GammaCalo");
   Bool_t kNewList = kFALSE;
   if (!listOutput){
      kNewList = kTRUE;
      listOutput = new TList();
      listOutput->SetName("GammaCalo");
   }

   for(Int_t i = 0; i<listInput->GetSize(); i++){
      TList *listToSave = (TList*)listInput->At(i);
      TString dirname = listToSave->GetName();
      cout<<dirname<<endl;
      if(listToSave){
         cout<<"found"<<endl;
         listOutput->Add(listToSave);
      }
   }
   
   listOutput->Write("",TObject::kSingleKey);
   fileOutput->Close();
}
