// provided by Gamma Conversion Group, PWGGA
// Friederike Bock: fbock@cern.ch

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
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ExtractSignalBinning.h"
// #include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
using namespace std;


void CombineDifferentFolders(   TString nameInputFile1, 
                                TString nameInputFile2, 
                                TString namefileOutput, 
                                Int_t mode = 0
                            ){

    TFile *fileInput1 = NULL;
    fileInput1 = new TFile(nameInputFile1.Data());
    cout << fileInput1 << endl;
    TFile *fileInput2 = NULL;
    fileInput2 = new TFile(nameInputFile2.Data());
    cout << fileInput2 << endl;

    TString autoDetectedMainDir1     = AutoDetectMainTList(mode , fileInput1);
    TString autoDetectedMainDir2     = AutoDetectMainTList(mode , fileInput2);
    TString defaultMainDir           = GetDefaultMainTListName(mode);
    
    if (autoDetectedMainDir1.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }
    if (autoDetectedMainDir2.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }
  
    TList *listInput1 =(TList*)fileInput1->Get(autoDetectedMainDir1.Data());
    if (listInput1 == NULL){ 
        return;
    }   
    TList *listInput2 =(TList*)fileInput2->Get(autoDetectedMainDir2.Data());
    if (listInput2 == NULL){ 
        return;
    }   
    
    TFile *fileOutput = new TFile(namefileOutput,"RECREATE");
    TList *listOutput =(TList*)fileOutput->Get(defaultMainDir.Data());
    Bool_t kNewList = kFALSE;
    if (!listOutput){
        kNewList = kTRUE;
        listOutput = new TList();
        listOutput->SetName(defaultMainDir.Data());
    }

    for(Int_t i = 0; i<listInput1->GetSize(); i++){
        TList *listToSave = (TList*)listInput1->At(i);
        TString dirname = listToSave->GetName();
        cout<<dirname<<endl;
        if(listToSave){
            cout<<"found"<<endl;
            listOutput->Add(listToSave);
        }
    }
    for(Int_t i = 0; i<listInput2->GetSize(); i++){
        TList *listToSave = (TList*)listInput2->At(i);
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
