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
#include "../CommonHeaders/FittingGammaConversion.h"
#include "../CommonHeaders/ConversionFunctions.h"

TF1* FitRecursiveGaussian (TH1* histo, Double_t precision, Double_t correctRange, Double_t fitRangeMin, Double_t fitRangeMax);
TF1* FitExpPlusGaussian(TH1D* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t mode, Double_t ptcenter);
TF1* FitBckg(TH1* fHisto, Double_t minFit, Double_t maxFit);
TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, TString selection, Double_t constPar = -1);
Float_t    FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2);
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

//****************************************************************************
//************** Function to Correct CaloNonLinearity4 ***********************
//****************************************************************************
void CorrectCaloNonLinearityV4( 
                                TString configFileName  = "config.txt",
                                TString suffix          = "eps",
                                Bool_t enableAddCouts   = kFALSE
                              ){
    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();


    // General options
    TString select              = "";
    TString optionEnergy        = "";
    Int_t mode                  = -1;
    
    // variables for data set indentifiers and maximum number of sets
    Int_t nSets                 = 0;
    TString fPlotLabelsRatio[6] = {"", "", "", "", "", ""};
    TString strDataFile[6]      = {"", "", "", "", "", ""};
    TString strMCFile[6]        = {"", "", "", "", "", ""};
    TString dataCut[6]          = {"", "", "", "", "", ""};
    TString mcCut[6]            = {"", "", "", "", "", ""};
    
    // color settings
    const Int_t nColor          = 13;
    const Int_t nStyle          = 7;
    Color_t color[nColor]       = {kBlack,633,807,/*800,*/418,/*kGreen+4,*/435,601,879,806,852,kCyan+3,426};
    Int_t markerStyle[nStyle]   = {24,25,27,28,29,30,31};
    Double_t massPi0            = 0.1349766;

    // pT range for mass fitting
    Int_t ptBinRange[2]         = {0, 17};
    Int_t firstPtBinSet[6]      = { -1, -1, -1, -1, -1,     -1};
    Double_t ptBinsForRebin[10] = { -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1 };
    Int_t rebin[10]             = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1};
    Int_t exampleBin[30]        = { -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1 };
    Int_t nExampleBins          = 0;
    Int_t fixedOffSet           = -1;
    Double_t startFit           = 3;
    
    // pT binning general initialization
    Int_t fNBinsPt              = 50;
    Double_t fBinsPt[51];
    for (Int_t i = 0; i< 51; i++){
        fBinsPt[i]              = -1.;
    }    
    
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
        // reading energy setup
        } else if (tempValue.BeginsWith("energy",TString::kIgnoreCase)){
            optionEnergy    = (TString)((TObjString*)tempArr->At(1))->GetString();
        // reading mode    
        } else if (tempValue.BeginsWith("mode",TString::kIgnoreCase)){    
            mode            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // reading number of data sets to be compared
        } else if (tempValue.BeginsWith("nSets",TString::kIgnoreCase)){    
            nSets            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // reding maximum number of pt bins
        } else if (tempValue.BeginsWith("maxNPtBins",TString::kIgnoreCase)){    
            fNBinsPt        = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // reading ptBins
        } else if (tempValue.BeginsWith("fBinsPt",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting ptBins" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 52 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    fBinsPt[i-1]        = ((TString)((TObjString*)tempArr->At(i))->GetString()).Atof();
                else 
                    i                   = tempArr->GetEntries();
            }    
        // reading min and max bin for pt
        } else if (tempValue.BeginsWith("ptRangeBins",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting ptBin bin Range" << endl;
            ptBinRange[0]               = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
            ptBinRange[1]               = ((TString)((TObjString*)tempArr->At(2))->GetString()).Atoi();        
        // read data file names
        } else if (tempValue.BeginsWith("dataFileNames",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting data file paths" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 6+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    strDataFile[i-1]    = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        // read cut string for data files    
        } else if (tempValue.BeginsWith("cutData",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting cutstring in data files" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 6+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    dataCut[i-1]        = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        // read MC file names
        } else if (tempValue.BeginsWith("MCFileNames",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting MC file paths" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 6+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    strMCFile[i-1]    = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        // read cut string for MC files    
        } else if (tempValue.BeginsWith("cutMC",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting cutstring in MC files" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 6+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    mcCut[i-1]        = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                   = tempArr->GetEntries();
            }    
        // read plot labels    
        } else if (tempValue.BeginsWith("plotLabels",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting cutstring in plot labels for different data combinations" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 6+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    fPlotLabelsRatio[i-1]   = ((TString)((TObjString*)tempArr->At(i))->GetString());
                else 
                    i                       = tempArr->GetEntries();
            }  
        // read example bins    
        } else if (tempValue.BeginsWith("exampleBins",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting exampleBins" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 30+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    exampleBin[i-1]         = (((TString)((TObjString*)tempArr->At(i))->GetString())).Atoi();
                else 
                    i                       = tempArr->GetEntries();
            }  
        // read rebin pt boundaries    
        } else if (tempValue.BeginsWith("rebinPt",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting rebinPt boundaries for rebins" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 10+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    ptBinsForRebin[i-1]     = (((TString)((TObjString*)tempArr->At(i))->GetString())).Atof();
                else 
                    i                       = tempArr->GetEntries();
            }  
        // read rebin values    
        } else if (tempValue.BeginsWith("rebinValue",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting value for rebin" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 10+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    rebin[i-1]              = (((TString)((TObjString*)tempArr->At(i))->GetString())).Atoi();
                else 
                    i                       = tempArr->GetEntries();
            }  
        // read first bin of every trigger set    
        } else if (tempValue.BeginsWith("firstPtBinSet",TString::kIgnoreCase)){    
            if (enableAddCouts) cout << "setting value for switches for trigger bins" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 6+1 ; i++){
                if (enableAddCouts) cout << i << "\t" <<((TString)((TObjString*)tempArr->At(i))->GetString()).Data() << endl;
                if (((TString)((TObjString*)tempArr->At(i))->GetString()).CompareTo("stop",TString::kIgnoreCase))
                    firstPtBinSet[i-1]      = (((TString)((TObjString*)tempArr->At(i))->GetString())).Atoi();
                else 
                    i                       = tempArr->GetEntries();
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
    cout << "energyFlag:\t"<< optionEnergy.Data() << endl;
    cout << "mode:\t"<< mode << endl;
    cout << "fNBinsPt:\t" << fNBinsPt << endl;
    cout << "ptBinning:" << endl;
    for (Int_t i = 0; i < fNBinsPt+1; i++){
        cout <<  fBinsPt[i] << "\t" ;
    }
    cout << endl;
    if (ptBinsForRebin[0] != -1 ){
        cout << "setting rebin values" << endl;
        for (Int_t i = 0; (ptBinsForRebin[i] != -1 && i < 10); i++ ){
            cout << "pt range:\t"<< ptBinsForRebin[i] ;
            if (ptBinsForRebin[i+1] != -1)
                cout << " - " << ptBinsForRebin[i+1] << "\t" ;
            else 
                cout << " - " << fBinsPt[fNBinsPt] << "\t" ;
            cout << "rebin invMass: \t" << rebin[i] << endl;
        }
    }
    cout << "pt bin range set: " << ptBinRange[0] << " - " << ptBinRange[1] << endl;
    
    cout << "exampleBins:" << endl;
    for (Int_t i = 0; exampleBin[i] != -1; i++){
        cout <<  exampleBin[i] << "\t" ;
        nExampleBins++;
    }
    cout << endl;
    cout << "nSets:\t" << nSets << endl;
    if ( !optionEnergy.CompareTo("") || mode == -1 || nSets == 0 || !select.CompareTo("") ){
        cout << "ABORTING: You are missing the select, nSets, energy or mode setting, can't continue like that" << endl;
        return;
    }
    cout << "**************************************************************************" << endl;
    cout << "Data set setup: " << endl;
    for (Int_t i = 0; i < nSets; i++){
        if ( strDataFile[i].CompareTo("") != 0 && strMCFile[i].CompareTo("") != 0 ){
            if (dataCut[i].CompareTo("") == 0 ){
                cout << "ERROR: you forgot to tell us which cut, aborting" << endl;
                return;
            }    
            if (mcCut[i].CompareTo("") == 0){
                cout << "----------" << endl;
                cout << "ATTENTION: setting same cut as data for MC" << endl;
                cout << "----------" << endl;
                mcCut[i]    = dataCut[i];
            }    
            if (fPlotLabelsRatio[i].CompareTo("") == 0){
                cout << "----------" << endl;
                cout << "ATTENTION: no label for plotting has been set for this combination" << endl;
                cout << "----------" << endl;
            }   
            if (firstPtBinSet[i] == -1 || firstPtBinSet[i] > fNBinsPt){
                cout << "----------" << endl;
                cout << "ATTENTION: will never switch to this trigger" << endl;
                cout << "----------" << endl;
            }    
            cout << i << "\t" << fPlotLabelsRatio[i].Data() << "\t" << strDataFile[i].Data() << "\t" << dataCut[i].Data() << "\t" << strMCFile[i].Data() << "\t" << mcCut[i].Data() << "\t switch at:\t" << fBinsPt[firstPtBinSet[i]] << endl;
        } else {
            cout << "ERROR: no correct data set combi was set, aborting" << endl;
            cout << strDataFile[i].Data() << "\t" <<  strMCFile[i].Data() << endl;
            return;
        }
    }
    cout << "**************************************************************************" << endl;
    cout << "**************************************************************************" << endl;


    //*******************************************************************************
    //***************** create outputfolder and root file ***************************
    //*******************************************************************************
    TString outputDir           = "CorrectCaloNonLinearity";
    gSystem->Exec("mkdir -p "+outputDir);
    TString outputDirSample     = Form("CorrectCaloNonLinearity/%s",select.Data());
    gSystem->Exec("mkdir -p "+outputDirSample);
    // Output - root file
    TString nameOutput          = Form("%s/CorrectCaloNonLinearity_%s.root",outputDir.Data(),select.Data());
    TFile* fOutput              = new TFile(nameOutput,"RECREATE");
    
    
    //*******************************************************************************
    //***************** set up proper labeling **************************************
    //*******************************************************************************
    for (Int_t i = 0; i < nSets; i++){
        if (firstPtBinSet[i+1] != -1)
            fPlotLabelsRatio[i] = Form("%0.1f < #it{E}_{Cluster} < %0.1f GeV : %s",fBinsPt[firstPtBinSet[i]],fBinsPt[firstPtBinSet[i+1]],fPlotLabelsRatio[i].Data());
        else 
            fPlotLabelsRatio[i] = Form("#it{E}_{Cluster} #geq %0.1f GeV : %s",fBinsPt[firstPtBinSet[i]],fPlotLabelsRatio[i].Data());
    }    
    
    TString recGamma            = ReturnFullTextReconstructionProcess(mode);
    TString fTextMeasurement    = Form("#pi^{0} #rightarrow #gamma#gamma");
    TString fCollisionSystem    = ReturnFullCollisionsSystem(optionEnergy);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    
    //*******************************************************************************
    //********************* setting hist names to be loaded *************************
    //*******************************************************************************
    TString nameSignalHisto     = "";
    TString nameBGHisto         = "";
    if(mode==2||mode==3){
        nameSignalHisto         = "ESD_Mother_InvMass_E_Calib";
        nameBGHisto             = "ESD_Background_InvMass_E_Calib";
    } else {
        nameSignalHisto         = "ESD_Mother_InvMass_vs_Pt_Alpha";
        nameBGHisto             = "ESD_Background_InvMass_vs_Pt_Alpha";
    }
    
    //*******************************************************************************
    //**********************  create output histos **********************************
    //*******************************************************************************
    TH1D* histMCResults         = new TH1D("Mean mass MC","; #it{E}_{Cluster} (GeV); mean mass MC",fNBinsPt,fBinsPt);
    TH1D* histDataResults       = new TH1D("Mean mass Data","; #it{E}_{Cluster} (GeV); mean mass Data",fNBinsPt,fBinsPt);
    TH1D* histDataMCResults     = new TH1D("Mean mass ratio MC/Data","; #it{E}_{Cluster} (GeV); mean mass ratio (MC/Data)",fNBinsPt,fBinsPt);
    histMCResults->SetDirectory(0);
    histDataResults->SetDirectory(0);
    histMCResults->GetXaxis()->SetRangeUser(fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    histDataResults->GetXaxis()->SetRangeUser(fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    histMCResults->GetYaxis()->SetRangeUser(0.12,0.142);
    histDataResults->GetYaxis()->SetRangeUser(0.12,0.142);
    TF1* fFitReco;
    TF1* fFitMassPos;
    TCanvas *canvas             = new TCanvas("canvas","",200,0,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvas, 0.1, 0.02, 0.06, 0.1);

    
    //*******************************************************************************
    //**************************** read Input ***************************************
    //*******************************************************************************

    TFile* dataFile             = new TFile(strDataFile[0].Data(),"READ"); 
    if(dataFile->IsZombie()) {cout << "Info: ROOT file '" << strDataFile[0].Data() << "' could not be openend, return!" << endl; return;}
    TString mainDirNameData     =  AutoDetectMainTList(mode , dataFile);
    TList* dataTopDir           = (TList*) dataFile->Get(mainDirNameData.Data()); 
    if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
    TList* dataTopContainer     = (TList*) dataTopDir->FindObject(Form("Cut Number %s",dataCut[0].Data())); 
    if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",dataCut[0].Data()) << " not found in Data-File" << endl; return;}
    TList* dataESDContainer     = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",dataCut[0].Data())); 
    if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",dataCut[0].Data()) << " not found in Data-File" << endl; return;}

    TFile* mcFile               = new TFile(strMCFile[0].Data(),"READ"); 
    if(mcFile->IsZombie()) {cout << "Info: ROOT file '" << strMCFile[0].Data() << "' could not be openend, return!" << endl; return;}
    TString mainDirNameMC       =  AutoDetectMainTList(mode , mcFile);
    TList* mcTopDir             = (TList*) mcFile->Get(mainDirNameMC.Data()); 
    if(mcTopDir == NULL) {cout << "ERROR: mcTopDir not Found"<<endl; return;}
    TList* mcTopContainer       = (TList*) mcTopDir->FindObject(Form("Cut Number %s",mcCut[0].Data())); 
    if(mcTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",mcCut[0].Data()) << " not found in MC-File" << endl; return;}
    TList* mcESDContainer       = (TList*) mcTopContainer->FindObject(Form("%s ESD histograms",mcCut[0].Data())); 
    if(mcESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",mcCut[0].Data()) << " not found in MC-File" << endl; return;}

    // Getting correct hist input
    TH2F* dataInvMassPt    = (TH2F*) dataESDContainer->FindObject(nameSignalHisto.Data());
    TH2F* dataBGInvMassPt  = (TH2F*) dataESDContainer->FindObject(nameBGHisto.Data());
    TH2F* mcInvMassPt      = (TH2F*) mcESDContainer->FindObject(nameSignalHisto.Data());
    TH2F* mcBGInvMassPt    = (TH2F*) mcESDContainer->FindObject(nameBGHisto.Data());
    if(!dataInvMassPt){cout << "did not find ESD_Mother_InvMass_E_Calib in data" << endl; return;}
    if(!dataBGInvMassPt){cout << "did not find ESD_Background_InvMass_E_Calib in data" << endl; return;}
    if(!mcInvMassPt){cout << "did not find ESD_Mother_InvMass_E_Calib in mc" << endl; return;}
    if(!mcBGInvMassPt){cout << "did not find ESD_Background_InvMass_E_Calib in mc" << endl; return;}

    
    // storing 2D hists in output root file
    fOutput->cd();
    dataInvMassPt->Write("Data - ESD_Mother_InvMass");
    dataBGInvMassPt->Write("Data - ESD_Background_InvMass");
    mcInvMassPt->Write("MC - ESD_Mother_InvMass");
    mcBGInvMassPt->Write("MC - ESD_Background_InvMass");
        
    TString dataMC[2]           = {"Data","MC"};
    Int_t triggerSel            = 0;
    Int_t rebinWindow           = 0;
    for(Int_t iClusterPt=ptBinRange[0]; iClusterPt<ptBinRange[1]; iClusterPt++){
        if(iClusterPt==firstPtBinSet[triggerSel+1]){
            // switching to trigger
            cout << endl;
            cout << "-----------------------------------------------------" << endl;
            cout << "\t Closing open files, switching to Trigger!" << endl;
            cout << "bin: " << firstPtBinSet[triggerSel] << endl;
            cout << "-----------------------------------------------------" << endl;

            if(!strDataFile[triggerSel+1].IsNull()){
                dataESDContainer->Clear(); dataTopContainer->Clear(); dataTopDir->Clear(); dataFile->Delete();
                dataFile            = new TFile(strDataFile[triggerSel+1].Data(),"READ"); 
                if(dataFile->IsZombie()) {cout << "Info: ROOT file '" << strDataFile[triggerSel+1].Data() << "' could not be openend, return!" << endl; return;}
                mainDirNameData     =  AutoDetectMainTList(mode , dataFile);
                dataTopDir          = (TList*) dataFile->Get(mainDirNameData.Data()); 
                if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
                dataTopContainer    = (TList*) dataTopDir->FindObject(Form("Cut Number %s",dataCut[triggerSel+1].Data())); 
                if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",dataCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
                dataESDContainer    = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",dataCut[triggerSel+1].Data())); 
                if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",dataCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
            }
            if(!strMCFile[triggerSel+1].IsNull()){
                mcESDContainer->Clear(); mcTopContainer->Clear(); mcTopDir->Clear(); mcFile->Delete();
                mcFile              = new TFile(strMCFile[triggerSel+1].Data(),"READ"); 
                if(mcFile->IsZombie()) {cout << "Info: ROOT file '" << strMCFile[triggerSel+1].Data() << "' could not be openend, return!" << endl; return;}
                mainDirNameMC       =  AutoDetectMainTList(mode , mcFile);
                mcTopDir            = (TList*) mcFile->Get(mainDirNameMC.Data()); 
                if(mcTopDir == NULL) {cout << "ERROR: mcTopDir not Found"<<endl; return;}
                mcTopContainer      = (TList*) mcTopDir->FindObject(Form("Cut Number %s",mcCut[triggerSel+1].Data())); 
                if(mcTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",mcCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
                mcESDContainer      = (TList*) mcTopContainer->FindObject(Form("%s ESD histograms",mcCut[triggerSel+1].Data())); 
                if(mcESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",mcCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
            }
            
            dataInvMassPt      = (TH2F*) dataESDContainer->FindObject(nameSignalHisto.Data());
            dataBGInvMassPt    = (TH2F*) dataESDContainer->FindObject(nameBGHisto.Data());
            mcInvMassPt        = (TH2F*) mcESDContainer->FindObject(nameSignalHisto.Data());
            mcBGInvMassPt      = (TH2F*) mcESDContainer->FindObject(nameBGHisto.Data());
            if(!dataInvMassPt){cout << "did not find ESD_Mother_InvMass_E_Calib in triggered data" << endl; return;}
            if(!dataBGInvMassPt){cout << "did not find ESD_Background_InvMass_E_Calib in triggered data" << endl; return;}
            if(!mcInvMassPt){cout << "did not find ESD_Mother_InvMass_E_Calib in trigger mc" << endl; return;}
            if(!mcBGInvMassPt){cout << "did not find ESD_Background_InvMass_E_Calib in trigger mc" << endl; return;}

            triggerSel++;
            fOutput->cd();
        }

        cout << endl;
        cout << "-----------------------------------------------------" << endl;
        cout << "\t MC/Data Fitting mass positions" << endl;
        cout << "loop: " << iClusterPt << ", " << fBinsPt[iClusterPt] << " - " << fBinsPt[iClusterPt+1] << " GeV" << endl;
        cout << "-----------------------------------------------------" << endl;

        TH2* Hist2D;
        TH2* HistBG2D;
        
        for(Int_t iDataMC = 0; iDataMC < 2; iDataMC++){
            if(iDataMC==0) {
                Hist2D              = dataInvMassPt;
                HistBG2D            = dataBGInvMassPt;
            } else if(iDataMC==1) {
                Hist2D              = mcInvMassPt;
                HistBG2D            = mcBGInvMassPt;
            } else {
                cout << "ERROR: data/mc loop, returning..." << endl; return;
            }
            Hist2D->Sumw2();
            HistBG2D->Sumw2();
            
            Double_t projectMin;
            Double_t projectMax;
            if(mode==2||mode==3){
              projectMin            = Hist2D->GetYaxis()->FindBin(fBinsPt[iClusterPt]+0.001);
              projectMax            = Hist2D->GetYaxis()->FindBin(fBinsPt[iClusterPt+1]-0.001);
            }else{
              projectMin            = Hist2D->GetYaxis()->FindBin((fBinsPt[iClusterPt]*2)+0.001);
              projectMax            = Hist2D->GetYaxis()->FindBin((fBinsPt[iClusterPt+1]*2)-0.001);
            }
            TH1D* sliceHist         = (TH1D*) Hist2D->ProjectionX(Form("slice%sAlpha_%f-%f",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]),projectMin,projectMax);
            sliceHist->SetDirectory(0);
            sliceHist->SetTitle(Form("%s - %.01f < #it{E}_{Cluster} < %.01f (GeV)",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
            sliceHist->GetYaxis()->SetTitle("#frac{d#it{M}_{inv}}{dN}");
            sliceHist->Sumw2();
            TH1D* sliceBGHist       = (TH1D*) HistBG2D->ProjectionX(Form("sliceBG%sAlpha_%f-%f",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]),projectMin,projectMax);
            sliceBGHist->SetDirectory(0);
            sliceBGHist->SetTitle(Form("%s - %.01f < #it{E}_{Cluster} < %.01f (GeV)",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
            sliceBGHist->GetYaxis()->SetTitle("#frac{d#it{M}_{inv}}{dN}");
            sliceBGHist->Sumw2();
            
            //*******************************************************************************
            // Rebin in mass hist
            if (ptBinsForRebin[rebinWindow+1] != -1 && rebinWindow < 10){
                cout << ptBinsForRebin[rebinWindow+1] << "\t" << fBinsPt[iClusterPt] << endl;
                if (ptBinsForRebin[rebinWindow+1] < fBinsPt[iClusterPt]){
                    rebinWindow++;
                    cout << "changed rebin factor: \t" << rebin[rebinWindow-1] << "\t -> \t" << rebin[rebinWindow] << endl;
                }
            }
            sliceHist->Rebin(rebin[rebinWindow]);
            sliceBGHist->Rebin(rebin[rebinWindow]);
        
            //*******************************************************************************
            // Background subtraction ranges
            Double_t integralSigAndBG   = sliceHist->Integral(sliceHist->FindBin(0.17), sliceHist->FindBin(0.3));
            Double_t integralBG         = sliceBGHist->Integral(sliceBGHist->FindBin(0.17), sliceBGHist->FindBin(0.3));
            cout << integralSigAndBG << "\t" << integralBG << "\t" << integralSigAndBG/ integralBG << endl;
            
            if (integralBG != 0){
                sliceBGHist->Scale( integralSigAndBG/ integralBG );
            }    
            for (Int_t i = 1; i< sliceHist->GetNbinsX()+1; i++){
                if (sliceHist->GetBinContent(i) == 0)
                    sliceHist->SetBinError(i,1.);
            }
            
            TH1D* sliceHistCopy         = (TH1D*)sliceHist->Clone("SliceCopy");
            if (integralBG != 0 && fBinsPt[iClusterPt]<16.0){
                sliceHist->Add( sliceBGHist, -1);
            }                

            sliceHistCopy->GetXaxis()->SetRangeUser(0.0,0.3);
            sliceHistCopy->Write(Form("slice%sAlpha_%f-%f-withWithBG",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
            sliceBGHist->SetLineColor(kGreen+2);
            sliceBGHist->GetXaxis()->SetRangeUser(0.0,0.3);
            sliceBGHist->Write(Form("slice%sAlpha_%f-%f-BG",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));                
            sliceHist->SetLineColor(kRed+2);
            sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);
            sliceHist->Write(Form("slice%sAlpha_%f-%f-withRemainingBckg",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
            
            sliceHistCopy->GetYaxis()->SetRangeUser(sliceHist->GetMinimum(),sliceHistCopy->GetMaximum());
            sliceHistCopy->GetXaxis()->SetRangeUser(0.0,0.3);
            sliceHistCopy->DrawCopy();
            sliceBGHist->DrawCopy("same");
            sliceHist->DrawCopy("same");
            DrawGammaLines(0., 0.3, 0, 0, 1, kGray+1, 7);
            canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
            
            for (Int_t i = 0; i < nExampleBins; i++ ){
                if(iClusterPt==exampleBin[i] ){
                    canvas->SaveAs(Form("%s/ExampleBin_%sAlpha_%.01f-%.01f-withBckgAndFit.%s",outputDirSample.Data(),dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1],suffix.Data()));
                }
            }    
            canvas->Write(Form("canvas_slice%sAlpha_%.01f-%.01f-withBckgAndFit",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
            
            canvas->Clear();
            sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);

            
            //*******************************************************************************
            // adjust min and max fitting range for inv mass fits
            //*******************************************************************************
            Double_t minMax[2]={0.04,0.3};
            // special setting for PCM-EMC
            if( mode == 2 ){
                if (fBinsPt[iClusterPt] < 1)
                    minMax[1]       = 0.2;
                if (fBinsPt[iClusterPt] < 2)
                    minMax[1]       = 0.22;
                if (fBinsPt[iClusterPt] < 3)
                    minMax[1]       = 0.25;
                minMax[0]       = 0.02;
                Double_t min    = 0.005*fBinsPt[iClusterPt] - 0.001;
                if (min > minMax[0])
                    minMax[0]   = min;
            // special setting for PCM-PHOS
            } else if( mode == 3 ){
                if (fBinsPt[iClusterPt] < 1)
                    minMax[1]       = 0.20;
                else 
                    minMax[1]       = 0.25;
                minMax[0]       = 0.03;
            // special setting for EMC
            } else if( mode == 4 ){
                minMax[1]       = 0.25;
                Double_t min    = 0.02*fBinsPt[iClusterPt] - 0.001;
                if (min > minMax[0])
                    minMax[0]   = min;
            // special setting for PHOS    
            } else if( mode == 5){
                minMax[1]       = 0.25;
                Double_t min    = 0.005*fBinsPt[iClusterPt] - 0.001;
                if (min > minMax[0])
                    minMax[0]   = min;
            }
            cout << "invMass fit range: \t" << minMax[0] << "\t" << minMax[1] << endl;
            //*******************************************************************************
            // Fit
            //*******************************************************************************
            fFitReco = FitExpPlusGaussian (sliceHist, minMax[0], minMax[1], mode, (fBinsPt[iClusterPt]+fBinsPt[iClusterPt+1])/2);

            if(iDataMC==0) {
                histDataResults->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
                histDataResults->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
            }
            else if(iDataMC==1) {
                histMCResults->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
                histMCResults->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
            }
            fFitReco->SetNpx(100000); 
            
            sliceHist->GetListOfFunctions()->Add(fFitReco);
            sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);
            
            sliceHist->DrawCopy();
            DrawGammaLines(0., 0.3, 0, 0, 1, kGray+1, 7);
            canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
            for (Int_t i = 0; i < nExampleBins; i++ ){
                if(iClusterPt==exampleBin[i] ){
                    canvas->SaveAs(Form("%s/ExampleBin_%sAlpha_%.01f-%.01f.%s",outputDirSample.Data(),dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1],suffix.Data()));
                }
            }
            canvas->Write(Form("canvas_slice%sAlpha_%.01f-%.01f",dataMC[iDataMC].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
            canvas->Clear();
        }
    }
    delete canvas;

    cout << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;

    histDataMCResults->Divide(histMCResults,histDataResults,1,1);
    DrawGammaSetMarker(histDataMCResults, 24, 2, kBlack, kBlack);

    Double_t minPlotY = 0.95;
    if(select.Contains("LHC10-Calo")) minPlotY = 0.9;
    if(select.Contains("LHC12-JetJet")) minPlotY = 0.93;


    //*********************************************************************************************************************************
    //************************************ Write mean mass for MC and data into output file *******************************************
    //*********************************************************************************************************************************    
    SetStyleHistoTH1ForGraphs(histMCResults, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC)} #GT",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histMCResults, markerStyle[1], 1, color[1], color[1]);
    histMCResults->Write("Mean mass MC");
    SetStyleHistoTH1ForGraphs(histDataResults, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (data)} #GT",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histDataResults, markerStyle[0], 1, color[0], color[0]);
    histDataResults->Write("Mean mass Data");
    
    //*********************************************************************************************************************************
    //*********************************** Plotting Mean mass for data and MC vs PDG value *********************************************
    //*********************************************************************************************************************************
    TCanvas *canvasMassPDG = new TCanvas("canvasMassPDG","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvasMassPDG, 0.1, 0.02, 0.06, 0.1);
    canvasMassPDG->SetLogx(1); 
    canvasMassPDG->SetLogy(0); 
  
    TH2F * histoDummyMeanMassVsPDG;
    histoDummyMeanMassVsPDG = new TH2F("histoDummyMeanMassVsPDG","histoDummyMeanMassVsPDG",11000,fBinsPt[ptBinRange[0]]/1.5,fBinsPt[ptBinRange[1]]*1.5,1000,0.89,1.1);
    SetStyleHistoTH2ForGraphs(histoDummyMeanMassVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC/data)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
    histoDummyMeanMassVsPDG->GetXaxis()->SetMoreLogLabels();
    histoDummyMeanMassVsPDG->GetXaxis()->SetLabelOffset(-0.01);
    histoDummyMeanMassVsPDG->DrawCopy("");
    
    Double_t rangeExponent[2]   = {-0.5, -0.08};
    Double_t rangeMult[2]       = {-0.2, -0.001};
    
    TLegend *legend = GetAndSetLegend2(0.15, 0.95, 0.95, 0.99, 0.043, 2, "", 42);
   
    TH1D* histDataResultsVsPDG =  (TH1D*)histDataResults->Clone("Mean mass data / mass PDG Pi0");
    histDataResultsVsPDG->Scale(1/massPi0);
    SetStyleHistoTH1ForGraphs(histDataResultsVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (data)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histDataResultsVsPDG, markerStyle[0], 1, color[0], color[0]);
    TF1* fitMassDataVsPDG = new TF1("fitMassDataVsPDG", "[0] + [1]*pow(x,[2])" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    fitMassDataVsPDG->SetParLimits(1, rangeMult[0], rangeMult[1]);
    fitMassDataVsPDG->SetParLimits(2, rangeExponent[0], rangeExponent[1]);
    
    histDataResultsVsPDG->Fit(fitMassDataVsPDG,"QRME0");
    fitMassDataVsPDG->SetLineColor(kBlack);
    fitMassDataVsPDG->Draw("same");
    cout << WriteParameterToFile(fitMassDataVsPDG) << endl;
    
    histDataResultsVsPDG->Write();
    histDataResultsVsPDG->DrawCopy("same");
    legend->AddEntry(histDataResultsVsPDG,"Data");

    TH1D* histMCResultsVsPDG =  (TH1D*)histMCResults->Clone("Mean mass MC / mass PDG Pi0");
    histMCResultsVsPDG->Scale(1/massPi0);
    SetStyleHistoTH1ForGraphs(histMCResultsVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
    DrawGammaSetMarker(histMCResultsVsPDG, markerStyle[1], 1, color[1], color[1]);
    TF1* fitMassMCVsPDG = new TF1("fitMassMCVsPDG", "[0] + [1]*pow(x,[2])" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    fitMassMCVsPDG->SetParLimits(1, rangeMult[0], rangeMult[1]);
    fitMassMCVsPDG->SetParLimits(2, rangeExponent[0], rangeExponent[1]);

    histMCResultsVsPDG->Fit(fitMassMCVsPDG,"QRME0");
    fitMassMCVsPDG->SetLineColor(kRed+2);
    fitMassMCVsPDG->Draw("same");
    cout << WriteParameterToFile(fitMassMCVsPDG) << endl;
    
    histMCResultsVsPDG->Write();
    histMCResultsVsPDG->DrawCopy("same");
    legend->AddEntry(histMCResultsVsPDG,"MC");
    
    PutProcessLabelAndEnergyOnPlot(0.94, 0.915, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 62, 0.03, "", 1, 1.25, 31);
    for (Int_t i = 0; i < nSets; i++){
       PutProcessLabelAndEnergyOnPlot(0.15, 0.915-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "");
    }
        
    legend->Draw("same");
    canvasMassPDG->Update();
    canvasMassPDG->SaveAs(Form("%s/MeanMass_Pi0_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasMassPDG->Clear();
    delete canvasMassPDG;

    //*********************************************************************************************************************************
    //****************************** Fitting ratio of mean mass position in MC/data ***************************************************
    //*********************************************************************************************************************************
    TF1* fFitConst = new TF1("DataMCConst", "[0]" ,startFit,fBinsPt[ptBinRange[1]]);
    histDataMCResults->Fit(fFitConst,"QRME0");
    Double_t highPtConst            = fixedOffSet;
    if (highPtConst == -1){
        highPtConst = fFitConst->GetParameter(0);
    }
    
    fFitMassPos = FitDataMC(histDataMCResults, fBinsPt[ptBinRange[0]], fBinsPt[ptBinRange[1]], select,  highPtConst);
    TF1* fFitComposit = new TF1("fFitComposit", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    fFitComposit->SetParameter(0, fitMassMCVsPDG->GetParameter(0) );
    fFitComposit->SetParameter(1, fitMassMCVsPDG->GetParameter(1) );
    fFitComposit->SetParameter(2, fitMassMCVsPDG->GetParameter(2) );
    fFitComposit->SetParameter(3, fitMassDataVsPDG->GetParameter(0) );
    fFitComposit->SetParameter(4, fitMassDataVsPDG->GetParameter(1) );
    fFitComposit->SetParameter(5, fitMassDataVsPDG->GetParameter(2) );

    
    TF1* fFitCompositFitted = new TF1("fFitCompositFitted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]-2]);
    fFitCompositFitted->SetParameter(0, fitMassMCVsPDG->GetParameter(0) );
    fFitCompositFitted->SetParameter(1, fitMassMCVsPDG->GetParameter(1) );
    fFitCompositFitted->FixParameter(2, fitMassMCVsPDG->GetParameter(2) );
    fFitCompositFitted->SetParameter(3, fitMassDataVsPDG->GetParameter(0) );
    fFitCompositFitted->SetParameter(4, fitMassDataVsPDG->GetParameter(1) );
    fFitCompositFitted->FixParameter(5, fitMassDataVsPDG->GetParameter(2) );
    histDataMCResults->Fit(fFitCompositFitted,"QRME0");
    
    TF1* fFitCompositInverted = new TF1("fFitCompositInverted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    fFitCompositInverted->SetParameter(0, fitMassDataVsPDG->GetParameter(0) );
    fFitCompositInverted->SetParameter(1, fitMassDataVsPDG->GetParameter(1) );
    fFitCompositInverted->SetParameter(2, fitMassDataVsPDG->GetParameter(2) );
    fFitCompositInverted->SetParameter(3, fitMassMCVsPDG->GetParameter(0) );
    fFitCompositInverted->SetParameter(4, fitMassMCVsPDG->GetParameter(1) );
    fFitCompositInverted->SetParameter(5, fitMassMCVsPDG->GetParameter(2) );
    
    TF1* fFitCompositInvertedFitted = new TF1("fFitCompositInvertedFitted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    fFitCompositInvertedFitted->SetParameter(0, fFitCompositFitted->GetParameter(3) );
    fFitCompositInvertedFitted->SetParameter(1, fFitCompositFitted->GetParameter(4) );
    fFitCompositInvertedFitted->SetParameter(2, fFitCompositFitted->GetParameter(5) );
    fFitCompositInvertedFitted->SetParameter(3, fFitCompositFitted->GetParameter(0) );
    fFitCompositInvertedFitted->SetParameter(4, fFitCompositFitted->GetParameter(1) );
    fFitCompositInvertedFitted->SetParameter(5, fFitCompositFitted->GetParameter(2) );
    
    histDataMCResults->GetYaxis()->SetRangeUser(minPlotY,1.05);
    histDataMCResults->GetXaxis()->SetRangeUser(fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);

    histDataMCResults->Write("MeanMassRatioMCData-noFit");
    fFitMassPos->Write("MeanMassRatioMCData-Fit");
    histDataMCResults->GetListOfFunctions()->Add(fFitMassPos);

    fstream fLog;
    fLog.open(Form("%s/CorrectCaloNonLinearity3_%s.log",outputDir.Data(),select.Data()), ios::out);
    fLog << "FitDataMC results:" << endl;
    for(Int_t i=0;i<=2;i++) fLog << "Par " << i << ": " << fFitMassPos->GetParameter(i) << " +- " << fFitMassPos->GetParError(i) << endl;
    
    fLog << WriteParameterToFile(fFitComposit) << endl;
    fLog << WriteParameterToFile(fFitCompositFitted) << endl;
    fLog << WriteParameterToFile(fFitCompositInverted) << endl;
    fLog.close();
    
    //*******************************************************************************
    // plotting mass ratios
    //*******************************************************************************
    TCanvas *canvasMassRatioMCData = new TCanvas("canvasMassPDG","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvasMassRatioMCData, 0.1, 0.02, 0.02, 0.1);
    canvasMassRatioMCData->SetLogx(1); 
    canvasMassRatioMCData->SetLogy(0); 
    
    SetStyleHistoTH1ForGraphs(histDataMCResults, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC)} #GT / #LT M_{#pi^{0} (data)} #GT",0.035,0.043, 0.035,0.043, 1.,0.9);
    DrawGammaSetMarker(histDataMCResults, markerStyle[0], 1, color[0], color[0]);
    histDataMCResults->Draw();
    fFitComposit->SetLineColor(kGreen+2);
    fFitComposit->Draw("same");

    fFitMassPos->Draw("same");
    PutProcessLabelAndEnergyOnPlot(0.94, 0.96, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 62, 0.03, "", 1, 1.25, 31);
    for (Int_t i = 0; i < nSets; i++){
       PutProcessLabelAndEnergyOnPlot(0.15, 0.945-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "");
    }
    
    canvasMassRatioMCData->Update();
    canvasMassRatioMCData->SaveAs(Form("%s/MeanMassRatio_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasMassRatioMCData->Clear();
  
    //*******************************************************************************
    // plotting total correction
    //*******************************************************************************
    canvasMassRatioMCData->cd();
    TH1D* totalCorrection = new TH1D("Total Correction","; #it{E}_{Cluster} (GeV); correction factor",1000,0.3,50);

    
    SetStyleHistoTH1ForGraphs(totalCorrection, "#it{E}_{Cluster} (GeV)","correction factor",0.035,0.043, 0.035,0.043, 1.,0.9);
    DrawGammaSetMarker(totalCorrection, 8, 1, kBlack, kBlack);
    totalCorrection->GetYaxis()->SetRangeUser(minPlotY,1.1);
    SetLogBinningXTH(totalCorrection);
 
    for(Int_t iBin = 1; iBin < totalCorrection->GetNbinsX()+1; iBin++){
        Float_t e = totalCorrection->GetXaxis()->GetBinCenter(iBin);
        Float_t factor = 1;
        Float_t p0 = fFitMassPos->GetParameter(0);
        factor /= FunctionNL_kSDM(e,p0,fFitMassPos->GetParameter(1),fFitMassPos->GetParameter(2));
        totalCorrection->SetBinContent(iBin,factor);
    }
    totalCorrection->DrawCopy("p");
    fFitCompositInverted->SetRange(0.3,50);
    fFitCompositInverted->SetLineColor(kGreen+2);
    fFitCompositInverted->Draw("same");
    
    TLegend *legend2 = GetAndSetLegend2(0.2, 0.2, 0.4, 0.29, 0.03, 1, "", 42);
    legend2->AddEntry(totalCorrection,"Correction factor for MC, direct ratio fit","l");
    legend2->AddEntry(fFitCompositInverted,"Correction factor for MC from mass fits","l");
    legend2->Draw("same");

    PutProcessLabelAndEnergyOnPlot(0.94, 0.96, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 62, 0.03, "", 1, 1.25, 31);
    for (Int_t i = 0; i < nSets; i++){
       PutProcessLabelAndEnergyOnPlot(0.15, 0.945-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "");
    }   

    canvasMassRatioMCData->Update();
    canvasMassRatioMCData->SaveAs(Form("%s/TotalCorrection_%s.%s", outputDir.Data(), select.Data(), suffix.Data()));
    canvasMassRatioMCData->Clear();
    delete canvasMassRatioMCData;

    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;

    fOutput->Write();
    fOutput->Close();
    return;
}

//*******************************************************************************
//****************** functional form for SDM correction *************************
//*******************************************************************************

Float_t FunctionNL_kSDM(Float_t e, Float_t p0, Float_t p1, Float_t p2){
    return ( p0 + exp( p1 + ( p2 * e ) ) );
}

//*******************************************************************************
//************** define exclusion function for BG fitting ***********************
//*******************************************************************************
Double_t fitExcludeSignal(Double_t *x, Double_t *par)
{
    if (x[0] > 0.08 && x[0] < 0.17) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0]; //+ par[2]*x[0]*x[0]
}

//*******************************************************************************
//************** fit BG only ****************************************************
//*******************************************************************************
TF1* FitBckg(TH1* fHisto, Double_t minFit, Double_t maxFit){
    TF1* fFitBckg = new TF1("fFitBckg",fitExcludeSignal,minFit,maxFit,3);
    fFitBckg->SetLineColor(kBlue);
    fFitBckg->SetLineWidth(2);
    fFitBckg->SetLineStyle(1);
    fHisto->Fit(fFitBckg,"QRME0");
    return fFitBckg;
}

//*******************************************************************************
//************** fit ratio of data and MC mass positions ************************
//*******************************************************************************
TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, TString selection, Double_t constPar){

    cout << "running standard fit from " <<  minFit << "\t"<<  maxFit << endl;
    TF1* fFitReco = new TF1("DataMC", "[0]+exp([1]+([2]*x))" ,minFit,maxFit);

    fFitReco->SetParameter(0,1.);
    fFitReco->SetParameter(1,-1.);
    fFitReco->SetParameter(2,-0.5);
    if(constPar!=-1) fFitReco->FixParameter(0,constPar);
    fHisto->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(kRed);
    fFitReco->SetLineWidth(2);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        cout << "Parameters for DataMC: " << endl;
        for(Int_t i=0;i<=2;i++) cout << "Par " << i << ": " << fFitReco->GetParameter(i) << " +- " << fFitReco->GetParError(i) << endl;
    } else {
        cout << "DataMC fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }

    return fFitReco;
}

//*******************************************************************************
//************** Definition of fitting with pure Gaussian ***********************
//*******************************************************************************
TF1* FitExpPlusGaussian(TH1D* histo, Double_t fitRangeMin, Double_t fitRangeMax, Int_t mode, Double_t ptcenter ){

    Double_t mesonAmplitude =histo->GetMaximum();
    Double_t mesonAmplitudeMin;
    Double_t mesonAmplitudeMax;

    // special setting for PCM-EMC
    if (mode == 2){
        if (ptcenter > 1.5)
            mesonAmplitudeMin = mesonAmplitude*95./100.;
        else 
            mesonAmplitudeMin = mesonAmplitude*90./100.;
    // special setting for PCM-PHOS
    } else if (mode == 3){
        if (ptcenter > 1.)
            mesonAmplitudeMin = mesonAmplitude*95./100.;
        else 
            mesonAmplitudeMin = mesonAmplitude*90./100.;
    // special setting for EMC
    } else if (mode == 4){
        mesonAmplitudeMin = mesonAmplitude*80./100.;
    // special setting for PHOS
    } else if (mode == 5){
        if (ptcenter > 1.)
            mesonAmplitudeMin = mesonAmplitude*80./100.;
        else 
            mesonAmplitudeMin = mesonAmplitude*10./100.;
    } else {
        mesonAmplitudeMin = mesonAmplitude*50./100.;
    }    
    mesonAmplitudeMax = mesonAmplitude*400./100.;

    
    TF1* fFitReco    = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1.-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)", 
                               fitRangeMin, fitRangeMax);
    Double_t fMesonMassExpect = TDatabasePDG::Instance()->GetParticle(111)->Mass();
    
    fFitReco->SetParameter(0, mesonAmplitude);
    fFitReco->SetParameter(1, fMesonMassExpect);
    fFitReco->SetParameter(2, 0.01);
    fFitReco->SetParameter(3, 0.012);

    fFitReco->SetParLimits(0, mesonAmplitudeMin, mesonAmplitudeMax);
    if (mode == 4 || mode == 5){
        fFitReco->SetParLimits(1, fMesonMassExpect*0.7, fMesonMassExpect*1.3);
    } else {
        fFitReco->SetParLimits(1, fMesonMassExpect*0.9, fMesonMassExpect*1.1);
    }    
    fFitReco->SetParLimits(2, 0.001, 0.1);
    fFitReco->SetParLimits(3, 0.001, 0.09);
    
    histo->Fit(fFitReco,"QRME0");
    
    fFitReco->SetLineColor(kRed+1);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        cout << "Parameter for exponential+Gaussian "<< endl;
        cout << gMinuit->fCstatu.Data() << endl;
        cout << "Gausexp: \t" << fFitReco->GetParameter(0) <<"+-" << fFitReco->GetParError(0) << "\t " << fFitReco->GetParameter(1)<<"+-" << fFitReco->GetParError(1) << "\t "<< fFitReco->GetParameter(2) <<"+-" << fFitReco->GetParError(2)<<endl;
        cout << "chi2/ndf:" << fFitReco->GetChisquare()/fFitReco->GetNDF() << endl;
    } else {
        cout << "Exp+Gaussian fitting failed in with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }
    return fFitReco;
}

//*******************************************************************************
//************** Definition of fitting with Gaussian + Exp tail *****************
//*******************************************************************************
TF1* FitRecursiveGaussian (TH1* histo, Double_t precision, Double_t correctRange, Double_t fitRangeMin, Double_t fitRangeMax ){
    TF1 *f0             = new TF1("f0", "gaus", fitRangeMin,fitRangeMax);
    histo->Fit(f0,"0RMEQ");
    Double_t rp         = f0->GetParameter(2);
    Double_t mp         = f0->GetParameter(1);
    Double_t ymin       = mp -(rp * correctRange);
    Double_t ymax       = mp + (rp * correctRange);
    Double_t deviation  = 100;
    Int_t counter       = 0;
    TF1* f1             = new TF1 ("f1", "gaus", ymin, ymax);
    while(deviation > precision && counter < 100){
        f1->SetRange(ymin,ymax);
        histo->Fit(f1,"0RMEQ");
        Double_t rp2    = f1->GetParameter(2);
        if (rp2>rp){ 
            deviation   = rp2-rp;
        } else { 
            deviation   = rp-rp2 ;
        }
        rp              = rp2 ;
        mp              = f1->GetParameter(1);
        ymin            = mp -(rp * correctRange);
        ymax            = mp +(rp * correctRange);
        counter++;
    }
    delete f0;

    if(TString(gMinuit->fCstatu.Data()).CompareTo("CONVERGED") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("SUCCESSFUL") == 0 || TString(gMinuit->fCstatu.Data()).CompareTo("PROBLEMS") == 0){
        cout << "Parameters for FitRecursiveGaussian: " << endl;
        for(Int_t i=0;i<=2;i++) cout << "Par " << i << ": " << f1->GetParameter(i) << " +- " << f1->GetParError(i) << endl;
    } else {
        cout << "FitRecursiveGaussian fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }

    return f1;
}
