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


void SeparateDifferentCutnumbers(   TString nameInputFile1          = "",
                                    TString nameOutputFileBase      = "",
                                    Int_t mode                      = 0,
                                    Int_t splitAll                  = kFALSE,
                                    Bool_t removeDCAtree            = kFALSE
                                ){

    TFile *fileInput1 = new TFile(nameInputFile1.Data());
    if(fileInput1->IsZombie()){
      cout << "ERROR: file is zombie, returning..." << endl;
      return;
    }

    TString autoDetectedMainDir1     = AutoDetectMainTList(mode , fileInput1);
    TString defaultMainDir           = GetDefaultMainTListName(mode);

    if (autoDetectedMainDir1.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }

    TList *listInput1 =(TList*)fileInput1->Get(autoDetectedMainDir1.Data());
    if (listInput1 == NULL){
        return;
    }

    Int_t nDiffCutNumbers = 0;
    for(Int_t i = 0; i<listInput1->GetSize(); i++){
        TList *listToSave       = (TList*)listInput1->At(i);
        TString dirname         = listToSave->GetName();
        if (dirname.Contains("Cut Number ") && listToSave){
            nDiffCutNumbers++;
        }
    }
    if (nDiffCutNumbers < 2 && !splitAll){
        cout << "only one cut present + additional folders for basis cuts available, no need to split" << endl;
        return;
    }

    TString additionOutputName[20] = {"A", "B", "C", "D", "E",
                                      "F", "G", "H", "I", "J",
                                      "K", "L", "M", "N", "O",
                                      "P", "Q", "R", "S", "T" };

    Int_t n = 0;
    for(Int_t i = 0; i<listInput1->GetSize(); i++){
        TList *listToSave       = (TList*)listInput1->At(i);
        TString dirname         = listToSave->GetName();
        if (dirname.Contains("Cut Number ") && listToSave){
            cout<< "found:" << dirname<<endl;
            TString nameOutputFile  = Form("%s_%s.root", nameOutputFileBase.Data(), additionOutputName[n].Data());
            cout << "writing to file: " << nameOutputFile.Data() << endl;
            TFile *fileOutput       = new TFile(nameOutputFile.Data(),"RECREATE");
                TList *listOutput       = (TList*)fileOutput->Get(defaultMainDir.Data());
                Bool_t kNewList         = kFALSE;
                if (!listOutput){
                    kNewList            = kTRUE;
                    listOutput          = new TList();
                    listOutput->SetName(defaultMainDir.Data());
                }
                listOutput->Add(listToSave);
                listOutput->Write("",TObject::kSingleKey);
            fileOutput->Close();
            delete fileOutput;
            n++;
        } else {
            cout<< "found:" << dirname<<endl;
            TString nameOutputFile  = Form("%s_Basic.root", nameOutputFileBase.Data());
            TFile *fileOutput       = new TFile(nameOutputFile.Data(),"UPDATE");
                TList *listOutput       = (TList*)fileOutput->Get(defaultMainDir.Data());
                Bool_t kNewList         = kFALSE;
                if (!listOutput){
                    kNewList            = kTRUE;
                    listOutput          = new TList();
                    listOutput->SetName(defaultMainDir.Data());
                }
                listOutput->Add(listToSave);
                listOutput->Write("",TObject::kSingleKey);
            fileOutput->Close();
            delete fileOutput;

        }
    }
    delete fileInput1;
}
