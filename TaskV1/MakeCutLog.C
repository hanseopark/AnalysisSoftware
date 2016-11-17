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
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"
#include "../CommonHeaders/ConversionFunctions.h"
using namespace std;

void MakeCutLog(const char *inputRootFile = "GammaConvV1.root", const char *InputName, Int_t mode = 0){

    TString filename = inputRootFile; 
    TFile *file = new TFile(filename.Data());  
    if (file->IsZombie()) return;
    
    fstream outputFile(InputName,ios::out);
    if(!outputFile.is_open()){
        cout<<"Problem opening file"<<endl;
        return;
    }

    //  Char_t filename_input1[200] = (Form("%s%s",path,input1));	
    TKey *key;
    file->ls();

    cout<<file<<endl;
    //TIter nextkey(l.GetListOfKeys());

    TString autoDetectedMainDir     = AutoDetectMainTList(mode , file);
    if (autoDetectedMainDir.CompareTo("") == 0){
        cout << "ERROR: trying to read file, which is incompatible with mode selected" << endl;;
        return;
    }
    TList *list = (TList*) file->Get(autoDetectedMainDir.Data());

    for(Int_t i = 0; i<list->GetEntries(); i++){
        TList *l2 = (TList*) list->At(i);
        TString dirname = l2->GetName();
        if(dirname.BeginsWith("Cut") == 1){
            TString CutNumber(dirname(11,dirname.Length()-11));
            outputFile << CutNumber.Data() <<endl;
            cout<<CutNumber<<endl;
        }
    }
    outputFile.close();

}
