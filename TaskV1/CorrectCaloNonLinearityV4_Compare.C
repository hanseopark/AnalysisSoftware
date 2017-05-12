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
    TString addPathFilesAdd[10]     = {"", "", "", "", "",  "", "", "", "", "" };
    TString inputFileNames[10]      = {"", "", "", "", "",  "", "", "", "", "" };
    TString inputFileNamesAdd[10]   = {"", "", "", "", "",  "", "", "", "", "" };
    TString plotLabelsData[10]      = {"", "", "", "", "",  "", "", "", "", "" };
    TString plotLabelsMC[10]        = {"", "", "", "", "",  "", "", "", "", "" };
    TString plotLabelsRatio[10]     = {"", "", "", "", "",  "", "", "", "", "" };
    TString inputFilePaths[10]      = {"", "", "", "", "",  "", "", "", "", "" };
    TString inputFilePathsAdd[10]   = {"", "", "", "", "",  "", "", "", "", "" };
    
    TString nameUsedCorr[10]        = {"", "", "", "", "",  "", "", "", "", "" };
    TString nameUsedCorrAdd[10]     = {"", "", "", "", "",  "", "", "", "", "" };
    Bool_t plotMassData[10]         = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    Bool_t plotMassMC[10]           = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    Bool_t plotMassRatio[10]        = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    TString xTitle                  = "#it{E}_{Cluster} (GeV)";
    TFile *inputFiles[10]           = { 0x0, 0x0, 0x0, 0x0, 0x0,
                                        0x0, 0x0, 0x0, 0x0, 0x0 };
    TFile *inputFilesAdd[10]        = { 0x0, 0x0, 0x0, 0x0, 0x0,
                                        0x0, 0x0, 0x0, 0x0, 0x0 };
    Bool_t hasAddInput[10]          = { kFALSE, kFALSE, kFALSE, kFALSE, kFALSE,
                                        kFALSE, kFALSE, kFALSE, kFALSE, kFALSE };
    
                                        
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
        } else if (tempValue.BeginsWith("nNL",TString::kIgnoreCase)){
            nNL             = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // read input names
        } else if (tempValue.BeginsWith("inputFileNamesAdd",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "Setting input file names" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    inputFileNamesAdd[i-1] = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                       = tempArr->GetEntries();
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
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    addPathFilesAdd[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
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
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 11 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    nameUsedCorrAdd[i-1]  = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                     = tempArr->GetEntries();
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
//     cout << "energyFlag:\t"<< optionEnergy.Data() << endl;
//     cout << "mode:\t"<< mode << endl;
    
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
        cout << inputFileNamesAdd[i].Data() << endl;
        if ( inputFileNamesAdd[i].CompareTo("") != 0 ){
            inputFilePathsAdd[i] = Form("%s%s/CorrectCaloNonLinearity_%s.root", inputDir.Data(), addPathFilesAdd[i].Data(), inputFileNamesAdd[i].Data());
            cout << i << "\t" << inputFileNamesAdd[i].Data() << "\t" << inputFilePathsAdd[i].Data() ;
            cout << endl;
        }
    }
    cout << "**************************************************************************" << endl;
    cout << "**************************************************************************" << endl;
    

    //*******************************************************************************
    // Input

    for(Int_t i=0; i<nNL; i++){
        inputFiles[i] = new TFile(inputFilePaths[i].Data(),"READ");
        if(inputFiles[i]->IsZombie()) {cout << "Info: ROOT file '" << inputFilePaths[i].Data() << "' could not be openend, return!" << endl; return;}
        if (inputFilePathsAdd[i].CompareTo("") != 0){
            inputFilesAdd[i]  = new TFile(inputFilePathsAdd[i].Data(),"READ");
            if(inputFilesAdd[i]->IsZombie()) {cout << "Info: ROOT file '" << inputFilePathsAdd[i].Data() << "' could not be openend, return!" << endl; return;}
            hasAddInput[i]    = kTRUE;
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
    TLegend *legendData = GetAndSetLegend2(0.13, 0.95-(nNLData/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42);
    legendData->SetMargin(0.1);
    TLegend *legendMC   = GetAndSetLegend2(0.13, 0.95-(nNLMC/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42);
    legendMC->SetMargin(0.1);
    TLegend *legendRatio= GetAndSetLegend2(0.13, 0.95-(nNLRatio/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42);
    legendRatio->SetMargin(0.1);
    TLegend *legend     = GetAndSetLegend2(0.13, 0.95-(nNL/nColumns+1)*0.035, 0.75, 0.95, 0.035, nColumns, "", 42);
    legend->SetMargin(0.1);

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
        histoMeanMassRatio[i] = (TH1D*)inputFiles[i]->Get("MeanMassRatioMCData-noFit");
        if(!histoMeanMassRatio[i]){
        cout << "ERROR: Could not find histogram 'MeanMassRatioMCData-noFit' in " << inputFilePaths[i].Data() << endl;
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

    legendRatio->Draw("same");
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
    canvasRatio->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasRatio->Clear();
    legendRatio->Clear();

    //*******************************************************************************
    // plotting total correction
    Double_t minPlotY = 0.93;
    if(select.Contains("LHC10-Calo")) minPlotY = 0.9;
    if(select.Contains("LHC12-JetJet")) minPlotY = 0.93;

    TH1D* totalCorrection = new TH1D("Total Correction","; #it{E}_{Cluster} (GeV); correction factor",1000,0.3,50);
    SetStyleHistoTH1ForGraphs(totalCorrection, "#it{E}_{Cluster} (GeV)","correction factor",0.035,0.043, 0.035,0.043, 0.82,0.9);
    totalCorrection->GetYaxis()->SetRangeUser(minPlotY+0.031,1.1);
    totalCorrection->GetXaxis()->SetMoreLogLabels();
    totalCorrection->GetXaxis()->SetLabelOffset(-0.01);
    SetLogBinningXTH(totalCorrection);
    totalCorrection->DrawCopy("p");

    TF1* fitTotalCorrection[nNL];
    TF1* fitTotalCorrectionAdd[nNL];
    for(Int_t i=0; i<nNL; i++){
        TString drawOption = (i==0)?"p":"p, same";
        fitTotalCorrection[i] = (TF1*)inputFiles[i]->Get(Form("%s_TotalCorr",nameUsedCorr[i].Data()));
        if(!fitTotalCorrection[i]){
            cout << "ERROR: Could not find histogram 'Total Correction' in " << inputFilePaths[i].Data() << endl;
            continue;
        } else {
            if (hasAddInput[i]){
              fitTotalCorrectionAdd[i]  = (TF1*)inputFilesAdd[i]->Get(Form("%s_TotalCorr",nameUsedCorrAdd[i].Data()));
              if (fitTotalCorrectionAdd[i]){
                fitTotalCorrection[i]   = MultiplyTF1(fitTotalCorrection[i], fitTotalCorrectionAdd[i], Form("TotalCorr%d",i));
              } 
            }
        }  
        DrawGammaSetMarkerTF1( fitTotalCorrection[i], lineStyle[i], 2, color[i]);
        fitTotalCorrection[i]->Draw("same");
        legend->AddEntry(fitTotalCorrection[i],plotLabelsRatio[i].Data(),"l");
    }

    TH1D* testBeam = (TH1D*) totalCorrection->Clone("testbeam");
    testBeam->Reset("ICE");
    DrawGammaSetMarker(testBeam, markerStyle[nNL], 0.2, color[nNL], color[nNL]);
    testBeam->SetLineWidth(2);
    for(Int_t iBin = 1; iBin <= testBeam->GetNbinsX()+1; iBin++) {
      Float_t e = testBeam->GetXaxis()->GetBinCenter(iBin);
      Float_t factor = 1;
      factor *= FunctionNL_kPi0MCv3(e);
      factor /= FunctionNL_kTestBeamv3(e);
      testBeam->SetBinContent(iBin,factor);
    }
    legend->AddEntry(testBeam,"kPi0MCv3 / kTestBeamv3","l");
    testBeam->Draw("p,hist,c,same");

    legend->Draw("same");
    canvasRatio->SetLogx(1); canvasRatio->SetLogy(0); canvasRatio->SetLogz(0); canvasRatio->Update();
    canvasRatio->SaveAs(Form("%s/TotalCorrection_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasRatio->Clear();
    legend->Clear();

    return;
}
