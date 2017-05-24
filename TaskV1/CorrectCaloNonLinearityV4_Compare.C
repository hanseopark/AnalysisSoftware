#include <stdlib.h>
#include <iostream>
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
#include "TProfile2D.h"
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
#include "TMath.h"
#include "TPostScript.h"
#include "TGraphErrors.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "../CommonHeaders/PlottingGammaConversionHistos.h"
#include "../CommonHeaders/PlottingGammaConversionAdditional.h"
#include "../CommonHeaders/ConversionFunctionsBasicsAndLabeling.h"

Float_t FunctionNL_kTestBeamv3(Float_t e){
    return ( 0.9615 / ( 0.976941 *( 1. / ( 1. + 0.162310 * exp( -e / 1.08689 ) ) * 1. / ( 1. + 0.0819592 * exp( ( e - 152.338 ) / 30.9594 ) ) ) ) );
}
Float_t FunctionNL_kPi0MCv3(Float_t e){
    return ( 1.0 / ( 0.981039 * ( 1. / ( 1. + 0.113508 * exp( -e / 1.00173 ) ) * 1. / ( 1. + 0.0967998 * exp( ( e - 219.381 ) / 63.1604 ) ) ) ) );
}
Float_t FunctionNL_kPi0MC(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5, Float_t p6){
    return ( p6 / ( p0 * ( 1. / ( 1. + p1 * exp( -e / p2 ) ) * 1. / ( 1. + p3 * exp( ( e - p4 ) / p5 ) ) ) ) );
}
Float_t FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2){
    return ( p0 + exp( p1 + ( p2 * e ) ) );
}
Float_t FunctionNL_DPOW(Float_t e, Float_t p0, Float_t p1, Float_t p2, Float_t p3, Float_t p4, Float_t p5){
  Float_t ret = 1;
  if ((p3 +  p4 * TMath::Power(e,p5 ) ) != 0)
    ret = ( (p0 +  p1 * TMath::Power(e,p2 ) )/(p3 +  p4 * TMath::Power(e,p5 ) ) );
  if (ret != 0.)
    return ret;
  else
    return 1.;
}

template<class ForwardIt>
void SetLogBinningXTH(ForwardIt* histoRebin){
    TAxis *axisafter    = histoRebin->GetXaxis();
    Int_t bins          = axisafter->GetNbins();
    Double_t from       = axisafter->GetXmin();
    Double_t to         = axisafter->GetXmax();
    Double_t *newbins   = new Double_t[bins+1];
    newbins[0]          = from;
    Double_t factor     = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) 
        newbins[i]      = factor * newbins[i-1];
    axisafter->Set(bins, newbins);
    delete [] newbins;
    return;
}


//================================================================================================================
//Multiply two TF1s
//================================================================================================================
TF1* MultiplyTF1(TF1* f1, TF1* f2, TString name) {

    if (!f1 || !f2) return NULL;

    cout << "Multiply functions" << endl;
    Double_t xmin, xmax;
    f1->GetRange(xmin, xmax);
    Int_t nPar1                         = f1->GetNpar();
    Int_t nPar2                         = f2->GetNpar();
    TString formula1                    = f1->GetExpFormula();
    TString formula2                    = f2->GetExpFormula();

    for (Int_t i = 0; i< nPar2; i++){
        formula2.ReplaceAll(Form("[%d]",i), Form("[%d]",i+nPar1));
    }

    TF1* result = new TF1(name.Data(),Form("(%s)*(%s)",formula1.Data(), formula2.Data()), xmin, xmax);
    for (Int_t i = 0; i < nPar1; i++ ){
        result->SetParameter(i, f1->GetParameter(i));
    }
    for (Int_t j = 0; j < nPar2; j++ ){
        result->SetParameter(nPar1+j, f2->GetParameter(j));
    }

    return result;
}



//****************************************************************************
//************** Function to compare different CaloNonLinearities ************
//****************************************************************************
void CorrectCaloNonLinearityV4_Compare( TString configFileName  = "config.txt",
                                        TString suffix          = "eps",
                                        Bool_t enableAddCouts   = kFALSE
                                      ) {
    
    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();

    TString inputDir                = "CorrectCaloNonLinearity";
    TString outputDir               = Form("%s/Compare",inputDir.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    // Setup general arrays and variables
    TString select                  = "";
    Int_t nNL                       = 0;
    Int_t nNLData                   = 0;
    Int_t nNLMC                     = 0;
    Int_t nNLRatio                  = 0;
    TString addPathFiles[10]        = {"", "", "", "", "",  "", "", "", "", "" };
    TString addPathFilesAdd[3][10]  = { {"", "", "", "", "",  "", "", "", "", "" }, 
                                        {"", "", "", "", "",  "", "", "", "", "" },
                                        {"", "", "", "", "",  "", "", "", "", "" } };
    TString inputFileNames[10]      = {"", "", "", "", "",  "", "", "", "", "" };
    TString inputFileNamesAdd[3][10]= { {"", "", "", "", "",  "", "", "", "", "" }, 
                                        {"", "", "", "", "",  "", "", "", "", "" },
                                        {"", "", "", "", "",  "", "", "", "", "" } };
    TString plotLabelsData[10]      = {"", "", "", "", "",  "", "", "", "", "" };
    TString plotLabelsMC[10]        = {"", "", "", "", "",  "", "", "", "", "" };
    TString plotLabelsRatio[10]     = {"", "", "", "", "",  "", "", "", "", "" };
    TString inputFilePaths[10]      = {"", "", "", "", "",  "", "", "", "", "" };
    TString inputFilePathsAdd[3][10]= { {"", "", "", "", "",  "", "", "", "", "" },
                                        {"", "", "", "", "",  "", "", "", "", "" }, 
                                        {"", "", "", "", "",  "", "", "", "", "" } };
    
    TString nameUsedCorr[10]        = {"", "", "", "", "",  "", "", "", "", "" };
    TString nameUsedCorrAdd[3][10]  = { {"", "", "", "", "",  "", "", "", "", "" },
                                        {"", "", "", "", "",  "", "", "", "", "" }, 
                                        {"", "", "", "", "",  "", "", "", "", "" } };
    Bool_t plotMassData[10]         = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    Bool_t plotMassMC[10]           = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    Bool_t plotMassRatio[10]        = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    TString xTitle                  = "#it{E}_{Cluster} (GeV)";
    TFile *inputFiles[10]           = { 0x0, 0x0, 0x0, 0x0, 0x0,
                                        0x0, 0x0, 0x0, 0x0, 0x0 };
    TFile *inputFilesAdd[3][10]     = { { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 }, 
                                        { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 },
                                        { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 } };
    Bool_t hasAddInput[3][10]       = { { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE }, 
                                        { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE }, 
                                        { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE } };
    Bool_t isRecursAddInput[3][10]  = { { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE }, 
                                        { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE }, 
                                        { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE, kFALSE } };
    
    // Setup colors
    const Int_t nColor              = 11;
    const Int_t nMarkerStyle        = 10;
    Color_t color[nColor]           = { kBlack, 633, 807, /*800,*/ 418, 
                                        kYellow+1, 601, 879, 806, 852,
                                        kCyan+3, 426};
    Int_t markerStyle[nMarkerStyle] = {24, 25, 27, 28, 29, 
                                       30, 31, 20, 33, 34};
    Int_t lineStyle[nMarkerStyle]   = {1, 2, 7, 8, 9, 
                                       3, 4, 5, 10, 1};

    Bool_t isDataVariation          = kFALSE;
    TString optionEnergy            = "";
    TString optionPeriod            = "";
    Int_t mode                      = -1;
    Int_t nLinesLabel               = 0;
    //**************************************************************************************************************
    //******************************* Read config file for detailed settings ***************************************
    //**************************************************************************************************************
    // ATTENTION: The data set has to be separated with either tabs or spaces a mixture of 
    //            both will most likely lead to misconfigurations
    //**************************************************************************************************************
    
    cout << "INFO: You have chosen the following config file: " << configFileName.Data() << endl;
    ifstream fileConfigNonLin;
    fileConfigNonLin.open(configFileName,ios_base::in);
    if (!fileConfigNonLin) {
        cout << "ERROR: settings " << configFileName.Data() << " not found!" << endl;
        return;
    }
    
    // read settings from file
    for( TString tempLine; tempLine.ReadLine(fileConfigNonLin, kTRUE); ) {        
        // check if line should be considered
        if (tempLine.BeginsWith("%") || tempLine.BeginsWith("#")){
            continue;
        } 
        if (enableAddCouts) cout << tempLine.Data() << endl;
     
        // Separate the string according to tabulators
        TObjArray *tempArr  = tempLine.Tokenize("\t");
        if(tempArr->GetEntries()<1){
            cout << "nothing to be done" << endl;
            delete tempArr;
            continue;
        } else if (tempArr->GetEntries() == 1){
            // Separate the string according to space
            tempArr       = tempLine.Tokenize(" ");
            if(tempArr->GetEntries()<1){
                cout << "nothing to be done" << endl;
                delete tempArr;
                continue;
            } else if (tempArr->GetEntries() == 1 ) {
                cout << ((TString)((TObjString*)tempArr->At(0))->GetString()).Data() << " has not be reset, no value given!" << endl;    
                delete tempArr;
                continue;
            }
        }

        // Put them to the correct variables
        TString tempValue   = (TString)((TObjString*)tempArr->At(0))->GetString();
        // reading selection string
        if (tempValue.BeginsWith("select",TString::kIgnoreCase)){
            select          = (TString)((TObjString*)tempArr->At(1))->GetString();
        // setting options for general plot labeling     
        } else if (tempValue.BeginsWith("energy",TString::kIgnoreCase)){
            optionEnergy    = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("period",TString::kIgnoreCase)){
            optionPeriod    = (TString)((TObjString*)tempArr->At(1))->GetString();
        } else if (tempValue.BeginsWith("mode",TString::kIgnoreCase)){
            mode            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        } else if (tempValue.BeginsWith("nNL",TString::kIgnoreCase)){
            nNL             = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // read input names
        } else if (tempValue.BeginsWith("inputFileNamesAdd",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting input file names" << endl;
            Int_t currArray             = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
            for(Int_t i = 2; i<tempArr->GetEntries() && i < 12 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    inputFileNamesAdd[currArray][i-2]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                                   = tempArr->GetEntries();
            }    
        } else if (tempValue.BeginsWith("inputFileNames",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting input file names" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    inputFileNames[i-1] = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        // additional path files
        } else if (tempValue.BeginsWith("addPathFilesAdd",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting input file names" << endl;
            Int_t currArray             = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
            for(Int_t i = 2; i<tempArr->GetEntries() && i < 12 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    addPathFilesAdd[currArray][i-2]     = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                                   = tempArr->GetEntries();
            }    
        // additional path files
        } else if (tempValue.BeginsWith("addPathFiles",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting input file names" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    addPathFiles[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }   
        // read input labels data
        } else if (tempValue.BeginsWith("plotLabelsData",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting plot Labels Data" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    plotLabelsData[i-1]    = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        // read input labels MC
        } else if (tempValue.BeginsWith("plotLabelsMC",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting plot Labels MC" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    plotLabelsMC[i-1]    = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        // read input labels ratio
        } else if (tempValue.BeginsWith("plotLabelsRatio",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting plot Labels Ratio" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    plotLabelsRatio[i-1]= ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        } else if (tempValue.BeginsWith("nameUsedCorrAdd",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting name of used correction" << endl;
            Int_t currArray             = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
            for(Int_t i = 2; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    nameUsedCorrAdd[currArray][i-2]     = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                                   = tempArr->GetEntries();
            }    
        } else if (tempValue.BeginsWith("isRecursAddInput",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Check whether its recursion or not" << endl;
            Int_t currArray             = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
            for(Int_t i = 2; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase)){
                    if ( ((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("1") == 0){
                        isRecursAddInput[currArray][i-2]= kTRUE;
                        cout << "set recursion for " << i-2 << endl;
                    }
                } else {
                    i                                   = tempArr->GetEntries();
                }
            }
        } else if (tempValue.BeginsWith("nameUsedCorr",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting name of used correction" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    nameUsedCorr[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }
        }  else if (tempValue.BeginsWith("plotMassData",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "enable mass plotting" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase)){
                    if ( ((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("1") == 0)
                        plotMassData[i-1]   = kTRUE;
                } else {
                    i                   = tempArr->GetEntries();
                }
            }
        }  else if (tempValue.BeginsWith("plotMassMC",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "enable mass plotting" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase)){
                    if ( ((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("1") == 0)
                        plotMassMC[i-1] = kTRUE;
                } else {
                    i                   = tempArr->GetEntries();
                }
            }
        }  else if (tempValue.BeginsWith("plotMassRatio",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "enable mass plotting" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase)){
                    if ( ((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("1") == 0)
                        plotMassRatio[i-1] = kTRUE;
                } else {
                    i                   = tempArr->GetEntries();
                }
            }
        }  else if (tempValue.BeginsWith("color",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting color" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 12 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase)){
                    color[i-1]          = ((TString)((TObjString*)tempArr->At(i))->GetString()).Atoi();
                } else {
                    i                   = tempArr->GetEntries();
                }
            }
        }  else if (tempValue.BeginsWith("lineStyle",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting lineStyle" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 12 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase)){
                    lineStyle[i-1]      = ((TString)((TObjString*)tempArr->At(i))->GetString()).Atoi();
                } else {
                    i                   = tempArr->GetEntries();
                }
            }
        }    
        delete tempArr;
    }
    
    
    //**************************************************************************************************************
    //******************************* Check wether settings were valid *********************************************
    //**************************************************************************************************************
    
    cout << "**************************************************************************" << endl;
    cout << "**************** Settings found in config file ***************************" << endl;
    cout << "**************************************************************************" << endl;
    cout << "select:\t"<< select.Data() << endl;
    cout << "nNL:\t"<< nNL << endl;
    if (nNL > 10){
        nNL = 10;
        cout << "attention reset to maximum 10 non lins" << endl;
    }
    if (nNL == 0){
        cout << "ABORTING: You are missing the nNL setting, can't continue like that" << endl;
        return;
    }
    TString fCollisionSystem    = "";
    TString recGamma            = "";
    if (optionEnergy.CompareTo("") != 0){
        cout << "energyFlag:\t"<< optionEnergy.Data() << endl;
        fCollisionSystem    = ReturnFullCollisionsSystem(optionEnergy);
        nLinesLabel++;
    }    
    if (optionPeriod.CompareTo("") != 0){
        cout << "periodFlag:\t"<< optionPeriod.Data() << endl;
        nLinesLabel++;
    }    
    if (mode != -1 ){
        cout << "mode:\t"<< mode << endl;
        recGamma            = ReturnFullTextReconstructionProcess(mode);
        nLinesLabel++;
    }
    
    cout << "**************************************************************************" << endl;
    cout << "Data set setup: " << endl;
    for (Int_t i = 0; i < nNL; i++){
        if ( inputFileNames[i].CompareTo("") != 0 ){
            inputFilePaths[i] = Form("%s%s/CorrectCaloNonLinearity_%s.root", inputDir.Data(), addPathFiles[i].Data(), inputFileNames[i].Data());
            cout << i << "\t" << inputFileNames[i].Data() << "\t" << inputFilePaths[i].Data() ;
            if ( plotMassData[i] || plotMassMC[i] || plotMassRatio[i] ) cout << "\t plotting mass";
            if ( plotMassData[i] ){
              cout << "\t Data";
              nNLData++;
            }  
            if ( plotMassMC[i] ){
              cout << "\t MC";
              nNLMC++;
            }  
            if ( plotMassRatio[i] ){
              cout << "\t Ratio";
              nNLRatio++;
            }  
            cout << endl;
        } else {
            cout << "ERROR: no correct data set was set, aborting" << endl;
            cout << inputFileNames[i].Data() << endl;
            return;
        }
    }
    cout << "**************************************************************************" << endl;
    cout << "**************************************************************************" << endl;

    cout << "**************************************************************************" << endl;
    cout << "Data set setup: " << endl;
    for (Int_t i = 0; i < nNL; i++){
        for (Int_t k = 0; k < 3; k++ ){
//             cout << inputFileNamesAdd[k][i].Data() << endl;
            if ( inputFileNamesAdd[k][i].CompareTo("") != 0 && inputFileNamesAdd[k][i].CompareTo("-") != 0 ){
                inputFilePathsAdd[k][i] = Form("%s%s/CorrectCaloNonLinearity_%s.root", inputDir.Data(), addPathFilesAdd[k][i].Data(), inputFileNamesAdd[k][i].Data());
                cout << k << "\t" << i << "\t" << inputFileNamesAdd[k][i].Data() << "\t" << inputFilePathsAdd[k][i].Data() ;
                cout << endl;
            }
        }
    }    
    cout << "**************************************************************************" << endl;
    cout << "**************************************************************************" << endl;
    

    //*******************************************************************************
    // Input

    for(Int_t i=0; i<nNL; i++){
        inputFiles[i] = new TFile(inputFilePaths[i].Data(),"READ");
        if(inputFiles[i]->IsZombie()) {cout << "Info: ROOT file '" << inputFilePaths[i].Data() << "' could not be openend, return!" << endl; return;}
        for (Int_t k = 0; k < 3; k++){
            if (inputFilePathsAdd[k][i].CompareTo("") != 0){
                inputFilesAdd[k][i]     = new TFile(inputFilePathsAdd[k][i].Data(),"READ");
                if(inputFilesAdd[k][i]->IsZombie()) {cout << "Info: ROOT file '" << inputFilePathsAdd[k][i].Data() << "' could not be openend, return!" << endl; return;}
                hasAddInput[k][i]       = kTRUE;
            }  
        }    
    }

    //*******************************************************************************
    // Output

    TCanvas *canvas     = new TCanvas("canvas","",200,0,1350,900);  // gives the page size
    Double_t leftMargin = 0.1; Double_t rightMargin = 0.02; Double_t topMargin = 0.06; Double_t bottomMargin = 0.1;
    DrawGammaCanvasSettings(canvas, 0.1, 0.02, 0.02, 0.1);
    canvas->SetLogx(1); canvas->SetLogy(0);
    
    Int_t nColumns      = 1;
    if (nNL > 3) 
        nColumns        = 2;
    TLegend *legendData = GetAndSetLegend2(0.13, 0.95-(nNLData/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42, 0.1);
    TLegend *legendMC   = GetAndSetLegend2(0.13, 0.95-(nNLMC/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42, 0.1);
    TLegend *legendRatio= GetAndSetLegend2(0.13, 0.95-(nNLRatio/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42, 0.1);
    //*******************************************************************************
    // plotting masses MC+Data
    TH1D* histoMassMC[nNL];
    for(Int_t i=0; i<nNL; i++){
        if(!plotMassMC[i]) continue;
        TString drawOption = (i==0)?"p":"p, same";
        histoMassMC[i] = (TH1D*)inputFiles[i]->Get("Mean mass MC");
        if(!histoMassMC[i]){
            cout << "ERROR: Could not find histogram 'Mean mass MC' in " << inputFilePaths[i].Data() << endl;
            continue;
        }
        DrawGammaSetMarker(histoMassMC[i], markerStyle[i], 1, color[i], color[i]);
        histoMassMC[i]->GetYaxis()->SetTitleOffset(1.1);
        histoMassMC[i]->GetXaxis()->SetLabelOffset(-0.01);
        histoMassMC[i]->Draw(drawOption.Data());
        legendMC->AddEntry(histoMassMC[i],plotLabelsMC[i].Data());
    }

    PutProcessLabelAndEnergyOnPlot(0.94, 0.16+nLinesLabel*0.035, 0.035, fCollisionSystem.Data(), optionPeriod.Data(), recGamma.Data(), 42, 0.035, "", 1, 1.25, 31);
    
    legendMC->Draw("same");
    canvas->Update();
    canvas->SaveAs(Form("%s/Mass_MC_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvas->Clear();
    legendMC->Clear();
    
    TH1D* histoMassData[nNL];
    for(Int_t i=0; i<nNL; i++){
        if(!plotMassData[i]) continue;
        TString drawOption = (i==0)?"p":"p, same";
        histoMassData[i] = (TH1D*)inputFiles[i]->Get("Mean mass Data");

        if(!histoMassData[i]){
            cout << "ERROR: Could not find histogram 'Mean mass Data' in " << inputFilePaths[i].Data() << endl;
            continue;
        }
        DrawGammaSetMarker(histoMassData[i], markerStyle[i], 1, color[i], color[i]);
        histoMassData[i]->GetXaxis()->SetLabelOffset(-0.01);
        histoMassData[i]->GetYaxis()->SetTitleOffset(1.1);
        histoMassData[i]->Draw(drawOption.Data());
        legendData->AddEntry(histoMassData[i],plotLabelsData[i].Data());
    }

    PutProcessLabelAndEnergyOnPlot(0.94, 0.16+nLinesLabel*0.035, 0.035, fCollisionSystem.Data(), optionPeriod.Data(), recGamma.Data(), 42, 0.035, "", 1, 1.25, 31);
    
    legendData->Draw("same");
    canvas->Update();
    canvas->SaveAs(Form("%s/Mass_Data_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvas->Clear();
    legendData->Clear();

    //*******************************************************************************
    // plotting mass ratios
    TCanvas *canvasRatio     = new TCanvas("canvasRatio","",0,0,750,500);  // gives the page size
    DrawGammaCanvasSettings(canvasRatio, 0.08, 0.02, 0.02, 0.08);
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0);

    TH1D* histoMeanMassRatio[nNL];
    for(Int_t i=0; i<nNL; i++){
        if (!plotMassRatio[i]) continue;
        TString drawOption = (i==0)?"p":"p, same";
        histoMeanMassRatio[i] = (TH1D*)inputFiles[i]->Get("MeanMassRatioMCData");
        if(!histoMeanMassRatio[i]){
        cout << "ERROR: Could not find histogram 'MeanMassRatioMCData' in " << inputFilePaths[i].Data() << endl;
        continue;
        }
        histoMeanMassRatio[i]->SetTitle(" ");
        DrawGammaSetMarker(histoMeanMassRatio[i], markerStyle[i], 1, color[i], color[i]);
        if(select.Contains("LHC10-")) histoMeanMassRatio[i]->GetYaxis()->SetRangeUser(0.9,1.05);
//         histoMeanMassRatio[i]->GetYaxis()->SetTitleOffset(1.1);
        histoMeanMassRatio[i]->GetXaxis()->SetLabelOffset(-0.01);
        histoMeanMassRatio[i]->GetXaxis()->SetTitleOffset(1.1);
        histoMeanMassRatio[i]->Draw(drawOption.Data());
        legendRatio->AddEntry(histoMeanMassRatio[i],plotLabelsRatio[i].Data());
    }

    PutProcessLabelAndEnergyOnPlot(0.13, 0.15+nLinesLabel*0.035, 0.035, fCollisionSystem.Data(), optionPeriod.Data(), recGamma.Data(), 42, 0.035, "", 1, 1.25, 11);
    
    legendRatio->Draw("same");
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
    canvasRatio->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasRatio->Clear();
    legendRatio->Clear();

    //*******************************************************************************
    TCanvas *canvasRatioC     = new TCanvas("canvasRatioC","",0,0,750,500);  // gives the page size
    DrawGammaCanvasSettings(canvasRatioC, 0.065, 0.012, 0.02, 0.08);
    canvasRatioC->SetLogx(1); canvasRatioC->SetLogy(0);

    // plotting total correction
    Double_t minPlotY = 0.93;
    if(select.Contains("LHC10-Calo")) minPlotY = 0.9;
    if(select.Contains("LHC12-JetJet")) minPlotY = 0.93;

    TH1D* totalCorrection = new TH1D("Total Correction","",1000,0.3,50);
    SetStyleHistoTH1ForGraphs(totalCorrection, "#it{E}_{Cluster} (GeV)","correction factor",0.035,0.043, 0.035,0.043, 0.82,0.8);
    totalCorrection->GetYaxis()->SetRangeUser(minPlotY+0.031,1.1);
    totalCorrection->GetXaxis()->SetMoreLogLabels();
    totalCorrection->GetXaxis()->SetLabelOffset(-0.01);
    SetLogBinningXTH(totalCorrection);
    totalCorrection->DrawCopy("p");

    TLegend *legend     = GetAndSetLegend2(0.1, 0.945-(nNL/nColumns+1)*0.035, 0.62, 0.945, 0.035, nColumns, "", 42, 0.15);

    TF1* fitTotalCorrection[nNL];
    TF1* fitTotalCorrectionAdd[3][nNL];
    TH1D* totalCorrectionDiffNL[nNL]; 
    for(Int_t i=0; i<nNL; i++){
        totalCorrectionDiffNL[i]        = new TH1D(Form("%s_TotalCorrHist",nameUsedCorr[i].Data()),"",7000,0.0,70);
        SetStyleHisto(totalCorrectionDiffNL[i], 2, lineStyle[i], color[i]);
        
        fitTotalCorrection[i]           = (TF1*)inputFiles[i]->Get(Form("%s_TotalCorr",nameUsedCorr[i].Data()));
        if(!fitTotalCorrection[i]){
            cout << "ERROR: Could not find histogram 'Total Correction' in " << inputFilePaths[i].Data() << endl;
            continue;
        } else {
            for (Int_t l = 1; l < totalCorrectionDiffNL[i]->GetNbinsX()+1; l++){
                totalCorrectionDiffNL[i]->SetBinContent(l, fitTotalCorrection[i]->Eval(totalCorrectionDiffNL[i]->GetBinCenter(l)));
            }
            for (Int_t k = 0; k < 3; k++){
                if (hasAddInput[k][i]){
                    fitTotalCorrectionAdd[k][i]     = (TF1*)inputFilesAdd[k][i]->Get(Form("%s_TotalCorr",nameUsedCorrAdd[k][i].Data()));
                    if (fitTotalCorrectionAdd[k][i] && !isRecursAddInput[k][i]){    // needs to be multiplied
//                         cout << k<< "\t" << i << "\t" << fitTotalCorrectionAdd[k][i]->GetParameter(0) << endl;
                        fitTotalCorrection[i]       = MultiplyTF1(fitTotalCorrection[i], fitTotalCorrectionAdd[k][i], Form("%s_TotalCorr",nameUsedCorr[i].Data()));
                        for (Int_t l = 1; l < totalCorrectionDiffNL[i]->GetNbinsX()+1; l++){
                            totalCorrectionDiffNL[i]->SetBinContent(l, fitTotalCorrection[i]->Eval(totalCorrectionDiffNL[i]->GetBinCenter(l)));
                        }    
                    } else if (fitTotalCorrectionAdd[k][i]){                        // is a recursion and needs to be calculated
                        TH1D* dummyCorr         = (TH1D*)totalCorrectionDiffNL[i]->Clone("dummy"); 
                        for (Int_t l = 1; l < totalCorrectionDiffNL[i]->GetNbinsX()+1; l++){
                            Double_t ptOld              = totalCorrectionDiffNL[i]->GetBinCenter(l);
                            Double_t ptInt              = ptOld * dummyCorr->GetBinContent(l);
                            Double_t ptNew              = ptInt * fitTotalCorrectionAdd[k][i]->Eval(ptInt);
                            Double_t corrFactorTotal    = ptNew/ptOld; 
//                             if ( l%100 == 0 ) cout << ptOld   << "\t" << ptInt << "\t" << ptNew << "\t" << corrFactorTotal << endl;
                            totalCorrectionDiffNL[i]->SetBinContent(l, corrFactorTotal);
                        }    
                    }
                }
            }
        }
        totalCorrectionDiffNL[i]->Draw("hist,same,c");
        legend->AddEntry(totalCorrectionDiffNL[i],plotLabelsRatio[i].Data(),"l");
    }

    TH1D* testBeam = (TH1D*) totalCorrection->Clone("testbeam");
    testBeam->Reset("ICE");
    SetStyleHisto(testBeam, 2, 1, color[nNL]);
    for(Int_t iBin = 1; iBin <= testBeam->GetNbinsX()+1; iBin++) {
      Float_t e = testBeam->GetXaxis()->GetBinCenter(iBin);
      Float_t factor = 1;
      factor *= FunctionNL_kPi0MCv3(e);
      factor /= FunctionNL_kTestBeamv3(e);
      testBeam->SetBinContent(iBin,factor);
    }
    legend->AddEntry(testBeam,"kPi0MCv3 / kTestBeamv3","l");
    testBeam->Draw("hist,c,same");

    DrawGammaLines(0.3, 50.,1.0, 1.0, 1, kGray+2, 2);
    
    PutProcessLabelAndEnergyOnPlot(0.958, 0.15+nLinesLabel*0.035, 0.035, fCollisionSystem.Data(), optionPeriod.Data(), recGamma.Data(), 42, 0.035, "", 1, 1.25, 31);
    
    legend->Draw("same");
    canvasRatioC->Update();
    canvasRatioC->SaveAs(Form("%s/TotalCorrection_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    delete canvasRatioC;
    
    return;
}
