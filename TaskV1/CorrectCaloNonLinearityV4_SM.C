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
TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, TString selection, Double_t constPar = -1, Int_t mode = 2);
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
void CorrectCaloNonLinearityV4_SM(
                                TString configFileName  = "/home/gustav/AnalysisSoftware/AnalysisSoftware/TaskV1/ExampleConfigs/NonLin_SM_PCMEDC_1617.txt",
                                TString suffix          = "png",
                                Bool_t enableAddCouts   = kFALSE
                              ){
    gROOT->Reset();

    StyleSettingsThesis();
    SetPlotStyle();


    // General options
    TString select              = "";
    TString optionEnergy        = "";
    Int_t mode                  = -1;
    TString FittingFunction    ="";

    // variables for data set indentifiers and maximum number of sets
    Int_t nSets                 = 0;
    TString fPlotLabelsRatio[6] = {"", "", "", "", "", ""};
    TString strDataFile[6]      = {"", "", "", "", "", ""};
    TString strMCFile[6]        = {"", "", "", "", "", ""};
    TString dataCut[6]          = {"", "", "", "", "", ""};
    TString mcCut[6]            = {"", "", "", "", "", ""};

    // color settings
    const Int_t nColor          = 21;
    const Int_t nStyle          = 7;
    Color_t color[nColor]       = {kBlack, kRed, kBlue, kGreen + 2, kYellow + 2, kMagenta + 1, kCyan + 1, kTeal - 7, kRed - 3, kAzure - 7, kPink + 7, kGreen - 2, kRed - 2, kBlue - 2, kOrange - 3, kTeal - 5, kMagenta + 3, kCyan - 6, kViolet - 5, kOrange + 7, kGreen + 2};//{kBlack,633,807,/*800,*/418,/*kGreen+4,*/435,601,879,806,852,kCyan+3,426};
    Int_t markerStyle[nStyle]   = {24,25,27,28,29,30,31};
    Double_t massPi0            = 0.1349766;

    // pT range for mass fitting
    Int_t ptBinRange[2]         = {0, 17};
    Int_t firstPtBinSet[6]      = { -1, -1, -1, -1, -1,     -1};
    Double_t ptBinsForRebin[10] = { -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1 };
    Int_t rebin[10]             = { 1, 1, 1, 1, 1,  1, 1, 1, 1, 1};
    Int_t exampleBin[40]        = { -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1,     -1, -1, -1, -1, -1 };
    Int_t nExampleBins          = 0;
    Int_t fixedOffSet           = -1;
    Int_t numOfSM               = 20;


    // pT binning general initialization
    Int_t fNBinsPt              = 300;
    Double_t fBinsPt[301];
    for (Int_t i = 0; i< 301; i++){
        fBinsPt[i]              = -1.;
    }
    Double_t rangeHighPtFitMass[4]  = {5, 10, 5, 10};
    Double_t rangeHighPtFitRatio[2] = {3, 10};
    Bool_t isNotFirstIte        = kFALSE;
    Bool_t onlyConstFit         = kFALSE;
    Bool_t doLightOutput       = kTRUE;

    TString dataMC[21]           = {"ALL SMs", "SM0", "SM1","SM2","SM3","SM4","SM5","SM6","SM7","SM8","SM9","SM10","SM11","SM12","SM13","SM14","SM15","SM16","SM17","SM18","SM19"};

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
        // Reading the number of Supermodules
        } else if (tempValue.BeginsWith("numOfSM",TString::kIgnoreCase)){
            numOfSM            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // reading if only const fit is applied
        } else if (tempValue.BeginsWith("onlyConstFit",TString::kIgnoreCase)){
            onlyConstFit            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // reading if only light output (only inv. Mass for All & SM 0 is plotted)
        } else if (tempValue.BeginsWith("doLightOutput",TString::kIgnoreCase)){
            doLightOutput            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // reading number of data sets to be compared
        } else if (tempValue.BeginsWith("nSets",TString::kIgnoreCase)){
            nSets            = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // enable add plotting if not first iteration
        } else if (tempValue.BeginsWith("notFirstIteration",TString::kIgnoreCase)){
            if (((TString)((TObjString*)tempArr->At(1))->GetString()).CompareTo("1") == 0)
                isNotFirstIte    = kTRUE;
        // reading maximum number of pt bins
        } else if (tempValue.BeginsWith("maxNPtBins",TString::kIgnoreCase)){
            fNBinsPt        = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atoi();
        // reading ptBins
        } else if (tempValue.BeginsWith("fBinsPt",TString::kIgnoreCase)){
            if (enableAddCouts) cout << "setting ptBins" << endl;
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 302 ; i++){
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
        // reading range for high pt const fit to ratio
        } else if (tempValue.BeginsWith("rangeHighPtFitRatio",TString::kIgnoreCase)){
            rangeHighPtFitRatio[0]      = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
            rangeHighPtFitRatio[1]      = ((TString)((TObjString*)tempArr->At(2))->GetString()).Atof();
        // reading min and max pt range for mass fits
        } else if (tempValue.BeginsWith("rangeHighPtFitMass",TString::kIgnoreCase)){
            if (enableAddCouts) cout << "setting min and max pt range for mass fits" << endl;
            rangeHighPtFitMass[0]       = ((TString)((TObjString*)tempArr->At(1))->GetString()).Atof();
            rangeHighPtFitMass[1]       = ((TString)((TObjString*)tempArr->At(2))->GetString()).Atof();
            rangeHighPtFitMass[2]       = ((TString)((TObjString*)tempArr->At(3))->GetString()).Atof();
            rangeHighPtFitMass[3]       = ((TString)((TObjString*)tempArr->At(4))->GetString()).Atof();
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
            for(Int_t i = 1; i<tempArr->GetEntries() && i < 40+1 ; i++){
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
    if (isNotFirstIte) cout << "This is not the first iteration" << endl;
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
        if ( strDataFile[i].CompareTo("") != 0 ){
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
    TString outputDir               = "CorrectCaloNonLinearity_SM";
    gSystem->Exec("mkdir -p "+outputDir);
    TString outputDirSampleSummary  = Form("%s/%s",outputDir.Data(),((TString)ReturnTextReconstructionProcess(mode)).Data());
    TString outputDirSample         = Form("%s/%s",outputDirSampleSummary.Data(), select.Data());
    gSystem->Exec("mkdir -p "+outputDirSampleSummary);
    gSystem->Exec("mkdir -p "+outputDirSample);
    // Output - root file
    TString nameOutput              = Form("%s/CorrectCaloNonLinearity_%s.root",outputDirSampleSummary.Data(),select.Data());
    TFile* fOutput                = new TFile(nameOutput,"RECREATE");

    fstream fLogOffsets;
        fLogOffsets.open(Form("%s/SuperModuleOffsets_%s.log",outputDirSampleSummary.Data(),select.Data()), ios::out);
        fLogOffsets << "------------------------------------------------------------------------------------" << endl;
        fLogOffsets << "------------------------------------------------------------------------------------" << endl;
        fLogOffsets << "List of offsets for each SM compared to the average of all" << endl;
        fLogOffsets << "Energy of clusters in each supermodule needs to be divided by the respective offset:" << endl;
        fLogOffsets << "------------------------------------------------------------------------------------" << endl;
        fLogOffsets << "------------------------------------------------------------------------------------" << endl;
    //*******************************************************************************
    //***************** set up proper labeling **************************************
    //*******************************************************************************
    for (Int_t i = 0; i < nSets; i++){
        if (firstPtBinSet[i+1] != -1)
            fPlotLabelsRatio[i] = Form("%0.2f < #it{E}_{Cluster} < %0.2f GeV : %s",fBinsPt[firstPtBinSet[i]],fBinsPt[firstPtBinSet[i+1]],fPlotLabelsRatio[i].Data());
        else
            fPlotLabelsRatio[i] = Form("#it{E}_{Cluster} #geq %0.2f GeV : %s",fBinsPt[firstPtBinSet[i]],fPlotLabelsRatio[i].Data());
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
    TString nameSignalHisto      = "ESD_Mother_InvMass_E_Calib";
    TString nameBGHisto          = "ESD_Back_InvMass_E_Calib";
    TString nameSignalHisto2     = "ESD_Mother_InvMass_vs_Pt_Alpha"; // for old files
    TString nameBGHisto2         = "ESD_Back_InvMass_vs_Pt_Alpha"; // for old files

    //*******************************************************************************
    //**********************  create output histos **********************************
    //*******************************************************************************
    TH1D* histAllResults         = new TH1D("Mean mass MC","; #it{E}_{Cluster} (GeV); mean mass MC",fNBinsPt,fBinsPt);
    TH1D** histSMResults       = new TH1D*[numOfSM];
    for(Int_t iSM = 0; iSM < numOfSM; iSM++){
      histSMResults[iSM]       = new TH1D(Form("Mean mass Data SM %i",iSM), Form("Mean mass Data SM %i; #it{E}_{Cluster} (GeV); mean mass Data",iSM),fNBinsPt,fBinsPt);

      histSMResults[iSM]->SetDirectory(0);
      histSMResults[iSM]->GetXaxis()->SetRangeUser(fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
      histSMResults[iSM]->GetYaxis()->SetRangeUser(0.12,0.142);
    }
    TH1D** histSMWiseRatioResults     = new TH1D*[numOfSM];//

    histAllResults->SetDirectory(0);
    histAllResults->GetXaxis()->SetRangeUser(fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
    histAllResults->GetYaxis()->SetRangeUser(0.12,0.142);
    TF1* fFitReco;
    TF1* fFitMassPos = nullptr;
    TCanvas *canvas             = new TCanvas("canvas","",200,0,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvas, 0.1, 0.02, 0.06, 0.1);


    //*******************************************************************************
    //**************************** read Input ***************************************
    //*******************************************************************************

    TFile* dataFile             = new TFile(strDataFile[0].Data(),"READ");
    if(dataFile->IsZombie()) {cout << "Info: ROOT file '" << strDataFile[0].Data() << "' could not be openend, return!" << endl; return;}
    TString mainDirNameData     = AutoDetectMainTList(91 , dataFile); //ConvCaloCalibration_1_0_1   mode!!!
    TList* dataTopDir           = (TList*) dataFile->Get(mainDirNameData.Data());
    if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
    TList* dataTopContainer     = (TList*) dataTopDir->FindObject(Form("Cut Number %s",dataCut[0].Data()));
    if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",dataCut[0].Data()) << " not found in Data-File" << endl; return;}
    TList* dataESDContainer     = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",dataCut[0].Data()));
    if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",dataCut[0].Data()) << " not found in Data-File" << endl; return;}

  // Getting correct hist input
    TH2F** InvMassPt_SM = new TH2F*[numOfSM];
    TH2F** BGInvMassPt_SM = new TH2F*[numOfSM];
    TH2F* InvMassPt_All      = (TH2F*) dataESDContainer->FindObject(nameSignalHisto.Data());
    TH2F* BGInvMassPt_All    = (TH2F*) dataESDContainer->FindObject(nameBGHisto.Data());

    for(Int_t iSM = 0; iSM < numOfSM; iSM++){
      InvMassPt_SM[iSM]    = (TH2F*) dataESDContainer->FindObject(Form("%s_SM%i",nameSignalHisto.Data(),iSM));
      BGInvMassPt_SM[iSM]  = (TH2F*) dataESDContainer->FindObject(Form("%s_SM%i",nameBGHisto.Data(),iSM));
    }

    // storing 2D hists in output root file
    fOutput->cd();
    // dataInvMassPt->Write("Data - ESD_Mother_InvMass");
    // dataBGInvMassPt->Write("Data - ESD_Background_InvMass");
    // mcInvMassPt->Write("MC - ESD_Mother_InvMass");
    // mcBGInvMassPt->Write("MC - ESD_Background_InvMass");

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
                mainDirNameData     =  AutoDetectMainTList(91 , dataFile); //mode!!
                dataTopDir          = (TList*) dataFile->Get(mainDirNameData.Data());
                if(dataTopDir == NULL) {cout << "ERROR: dataTopDir not Found"<<endl; return;}
                dataTopContainer    = (TList*) dataTopDir->FindObject(Form("Cut Number %s",dataCut[triggerSel+1].Data()));
                if(dataTopContainer == NULL) {cout << "ERROR: " << Form("Cut Number '%s'",dataCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
                dataESDContainer    = (TList*) dataTopContainer->FindObject(Form("%s ESD histograms",dataCut[triggerSel+1].Data()));
                if(dataESDContainer == NULL) {cout << "ERROR: " << Form("'%s' ESD histograms",dataCut[triggerSel+1].Data()) << " not found in File" << endl; return;}
            }


            InvMassPt_All      = (TH2F*) dataESDContainer->FindObject(nameSignalHisto.Data());
            BGInvMassPt_All    = (TH2F*) dataESDContainer->FindObject(nameBGHisto.Data());

            for(Int_t iSM = 0; iSM < numOfSM; iSM++){
              InvMassPt_SM[iSM]    = (TH2F*) dataESDContainer->FindObject(Form("%s_SM%i",nameSignalHisto.Data(),iSM));
              BGInvMassPt_SM[iSM]  = (TH2F*) dataESDContainer->FindObject(Form("%s_SM%i",nameBGHisto.Data(),iSM));
            }

            triggerSel++;
            fOutput->cd();
        }

        cout << endl;
        cout << "-----------------------------------------------------" << endl;
        cout << "\t SM/Data Fitting mass positions" << endl;
        cout << "loop: " << iClusterPt << ", " << fBinsPt[iClusterPt] << " - " << fBinsPt[iClusterPt+1] << " GeV" << endl;
        cout << "-----------------------------------------------------" << endl;

        TH2* Hist2D;
        TH2* HistBG2D;

        for(Int_t iSM = 0; iSM < numOfSM + 1; iSM++){
            if(iSM==0) {
                Hist2D              = InvMassPt_All;
                HistBG2D            = BGInvMassPt_All;
            } else if(iSM > 0 && iSM <= numOfSM) {
                Hist2D              = InvMassPt_SM[iSM - 1];
                HistBG2D            = BGInvMassPt_SM[iSM - 1];
            } else {
                cout << "ERROR: data/SM loop, returning..." << endl; return;
            }
            cout<<"iSM: "<<iSM<<endl;
            Hist2D->Sumw2();
            HistBG2D->Sumw2();

            Double_t projectMin;
            Double_t projectMax;
            if(mode==14){
              projectMin            = Hist2D->GetYaxis()->FindBin(fBinsPt[iClusterPt]+0.001);
              projectMax            = Hist2D->GetYaxis()->FindBin(fBinsPt[iClusterPt+1]-0.001);
            }else{
              projectMin            = Hist2D->GetYaxis()->FindBin((fBinsPt[iClusterPt]*2)+0.001);
              projectMax            = Hist2D->GetYaxis()->FindBin((fBinsPt[iClusterPt+1]*2)-0.001);
            }
            TH1D* sliceHist         = (TH1D*) Hist2D->ProjectionX(Form("slice%sAlpha_%f-%f",dataMC[iSM].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]),projectMin,projectMax);
            sliceHist->SetDirectory(0);
            sliceHist->SetTitle(Form("%s - %.02f < #it{E}_{Cluster} < %.02f (GeV)",dataMC[iSM].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
            sliceHist->GetYaxis()->SetTitle("#frac{d#it{M}_{inv}}{dN}");
            sliceHist->Sumw2();
            TH1D* sliceBGHist       = (TH1D*) HistBG2D->ProjectionX(Form("sliceBG%sAlpha_%f-%f",dataMC[iSM].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]), projectMin,projectMax);
            sliceBGHist->SetDirectory(0);
            sliceBGHist->SetTitle(Form("%s - %.02f < #it{E}_{Cluster} < %.02f (GeV)",dataMC[iSM].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1]));
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
            Double_t minMGGBG           = 0.22;
            Double_t integralSigAndBG   = sliceHist->Integral(sliceHist->FindBin(minMGGBG), sliceHist->FindBin(0.29));
            Double_t integralBG         = sliceBGHist->Integral(sliceBGHist->FindBin(minMGGBG), sliceBGHist->FindBin(0.29));
            cout << integralSigAndBG << "\t" << integralBG << "\t" << integralSigAndBG/ integralBG << endl;

            if (integralBG != 0){
                sliceBGHist->Scale( integralSigAndBG/ integralBG );
            }
            for (Int_t i = 1; i< sliceHist->GetNbinsX()+1; i++){
                if (sliceHist->GetBinContent(i) == 0)
                    sliceHist->SetBinError(i,1.);
            }

            TH1D* sliceHistCopy         = (TH1D*)sliceHist->Clone("SliceCopy");
            if (integralBG != 0 /*&& fBinsPt[iClusterPt]<16.0*/){
                sliceHist->Add( sliceBGHist, -1);
            }

            sliceHistCopy->GetYaxis()->SetRangeUser(sliceHist->GetMinimum(),sliceHistCopy->GetMaximum());
            sliceHistCopy->GetXaxis()->SetRangeUser(0.0,0.3);
            sliceHistCopy->DrawCopy();
            sliceBGHist->DrawCopy("same");
            sliceHist->DrawCopy("same");
            DrawGammaLines(0., 0.3, 0, 0, 1, kGray+1, 7);
            canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
            if(iSM > 2 && doLightOutput){
              cout<<"light output, skipping draw..."<<endl;
            } else{
              for (Int_t i = 0; i < nExampleBins; i++ ){
                  if(iClusterPt==exampleBin[i] ){
                      canvas->SaveAs(Form("%s/ExampleBin_%sAlpha_%.02f-%.02f-withBckgAndFit.%s",outputDirSample.Data(),dataMC[iSM].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1],suffix.Data()));
                  }
              }
            }
            canvas->Clear();
            sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);


            //*******************************************************************************
            // adjust min and max fitting range for inv mass fits
            //*******************************************************************************
            Double_t minMax[2]={0.04,0.3};
            // special setting for PCM-EDC
            if(mode == 14){
                if (optionEnergy.Contains("PbPb") || optionEnergy.Contains("XeXe")){
                    minMax[1]   = 0.3;
                } else if (optionEnergy.Contains("pPb_5.023TeVRun2") ){
                    if (fBinsPt[iClusterPt] < 1)
                        minMax[1]       = 0.2;
                    else if (fBinsPt[iClusterPt] < 2)
                        minMax[1]       = 0.23;
                    else if (fBinsPt[iClusterPt] < 3)
                        minMax[1]       = 0.26;
                }
                minMax[0]       = 0.09;
                Double_t min    = 0.008*fBinsPt[iClusterPt] - 0.001;
                if (min > minMax[0] && min < 0.1){
                    minMax[0]   = min;
                } else {
                  minMax[0] = 0.1;
                }

            } else if( mode == 15){
                minMax[1]       = 0.25;
                Double_t min    = 0.02*fBinsPt[iClusterPt] - 0.001;
                if (min > minMax[0])
                    minMax[0]   = min;

                if (optionEnergy.Contains("PbPb") || optionEnergy.Contains("XeXe")){
                    minMax[1]   = 0.3;
                }
              }
            cout << "invMass fit range: \t" << minMax[0] << "\t" << minMax[1] << endl;
            //*******************************************************************************
            // Fit
            //*******************************************************************************
            fFitReco = FitExpPlusGaussian (sliceHist, minMax[0], minMax[1], mode, (fBinsPt[iClusterPt]+fBinsPt[iClusterPt+1])/2);

            if(iSM==0) {
                histAllResults->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
                histAllResults->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
                cout<<"\n\nGetBinContent ALL: "<< iClusterPt+1<<"  "<<fFitReco->GetParameter(1)<<"\n\n";
            }
            else if(iSM>0) {
                cout<<"\n\nGetBinContent: "<< iClusterPt+1<<"  "<<fFitReco->GetParameter(1)<<"\n\n";
                histSMResults[iSM - 1]->SetBinContent(iClusterPt+1,fFitReco->GetParameter(1));
                histSMResults[iSM - 1]->SetBinError(iClusterPt+1,fFitReco->GetParError(1));
            }
            fFitReco->SetNpx(100000);

            sliceHist->GetListOfFunctions()->Add(fFitReco);
            sliceHist->GetXaxis()->SetRangeUser(0.0,0.3);

            sliceHist->DrawCopy();
            DrawGammaLines(0., 0.3, 0, 0, 1, kGray+1, 7);
            canvas->SetLogx(0); canvas->SetLogy(0); canvas->SetLogz(0); canvas->Update();
            if(iSM > 2 && doLightOutput){
              cout<<"light output, skipping draw..."<<endl;
            } else{
              for (Int_t i = 0; i < nExampleBins; i++ ){
                  if(iClusterPt==exampleBin[i] ){
                      canvas->SaveAs(Form("%s/ExampleBin_%sAlpha_%.02f-%.02f.%s",outputDirSample.Data(),dataMC[iSM].Data(),fBinsPt[iClusterPt],fBinsPt[iClusterPt+1],suffix.Data()));
                  }
              }
            }
            canvas->Clear();

            cout<<"iSM: "<<iSM<<endl;
        }
    }
    delete canvas;

    cout << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << "-----------------------------------------------------" << endl;

    //*********************************************************************************************************************************
    //************************* Initialize histos, legends, function ... needed for every Supermodule *********************************
    //*********************************************************************************************************************************

    TCanvas *canvas_All = new TCanvas("canvasMassRatioMCData_All","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvas_All, 0.082, 0.15, 0.02, 0.1);
    canvas_All->SetLogx(1);
    canvas_All->SetLogy(0);
    TCanvas *canvas_All_Fitted = new TCanvas("canvasMassRatioMCData_All_Fitted","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvas_All_Fitted, 0.082, 0.15, 0.02, 0.1);
    canvas_All_Fitted->SetLogx(1);
    canvas_All_Fitted->SetLogy(0);
    TCanvas *canvas_All_FitOnly = new TCanvas("canvasMassRatioMCData_All_FitOnly","",200,10,1350,900);  // gives the page size
    DrawGammaCanvasSettings(canvas_All_FitOnly, 0.082, 0.15, 0.02, 0.1);
    canvas_All_FitOnly->SetLogx(1);
    canvas_All_FitOnly->SetLogy(0);
    Bool_t fitconst = kTRUE;

    TH1D* histSMResultsRatio = nullptr;
    TH1D* histAllResultsVsPDG = nullptr;
    TH2F * histoDummyMeanMassVsPDG = nullptr;
    TH1D* histSMResultsVsPDG = nullptr;

    TF1* fitMassDataVsPDG = nullptr;
    TF1* fitMassDataVsPDGConst = nullptr;
    TF1* fitMassMCVsPDG = nullptr;
    TF1* fitMassDataVsPDG2 = nullptr;
    TF1* fitMassMCVsPDGConst = nullptr;
    TF1* fitMassMCVsPDG2 = nullptr;
    TF1* fFitConst = nullptr;
    TF1* fFitConstFull = nullptr;
    TF1* fFitConstFullInv = nullptr;
    TF1* fFitMassPosInverted = nullptr;
    TF1* fFitCompositFitted = nullptr;
    TF1* fFitComposit = nullptr;
    TF1* fFitCompositInverted = nullptr;
    TF1* fFitCompositInvertedFitted = nullptr;
    TF1* fFitExpComb = nullptr;
    TF1* fFitExpCombInverted = nullptr;
    TH2F * histoDummyDataMCRatio = nullptr;
    TH1D* totalCorrection = nullptr;


    // TCanvas *canvasMassPDG = nullptr;
    TCanvas *canvasMassRatioMCData = nullptr;

    TLegend* legend_All = GetAndSetLegend2(0.87,0.1, 0.98,0.98, 0.04, 1, "", 42, 0.);
    TLegend *legendPDG = nullptr;
    TLegend *legendFits = nullptr;
    TLegend* legendCorrectionFunctions = nullptr;
    TLegend *legend2 = nullptr;



    Double_t minPlotY = 0.93;
    if (select.Contains("LHC16"))           minPlotY = 0.94;
    if (select.Contains("LHC17"))           minPlotY = 0.94;

    Double_t maxPlotY = 1.06;

    Double_t minMass  = 0.89;
    Double_t maxMass  = 1.1;
    if (mode == 15){
        minMass  = 0.87;
        maxMass  = 1.25;
    }
    Bool_t isFirstLoop = kTRUE;
    //*********************************************************************************************************************************
    //********************************************* Start loop over all Supermodules **************************************************
    //*********************************************************************************************************************************

    for(Int_t iSM = numOfSM - 1; iSM >= 0; iSM-=1){


      SetStyleHistoTH1ForGraphs(histSMResults[iSM], "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (SM)} #GT",0.035,0.043, 0.035,0.043, 1.,1.);
      DrawGammaSetMarker(histSMResults[iSM], markerStyle[1], 1, color[iSM + 1], color[iSM + 1]);
      SetStyleHistoTH1ForGraphs(histAllResults, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (EDC)} #GT",0.035,0.043, 0.035,0.043, 1.,1.);
      DrawGammaSetMarker(histAllResults, markerStyle[0], 1, color[0], color[0]);

      //*********************************************************************************************************************************
      //************************************ Crate Ratio of SM wise vs. EDC integrated **************************************************
      //*********************************************************************************************************************************

      histSMResultsRatio = (TH1D*) histSMResults[iSM]->Clone("SM / EDC");
      cout<<"\n\n "<<"Before: " <<histSMResults[iSM]->GetBinContent(3) << "  "<< histAllResults->GetBinContent(3) <<endl<<endl;
      histSMResultsRatio->Divide(histSMResults[iSM],histAllResults,1,1);
      if(mode == 14) histSMResultsRatio->Multiply(histSMResultsRatio,histSMResultsRatio,1.,1.,"B");
      DrawGammaSetMarker(histSMResultsRatio, 24, 2, kBlack, kBlack);



      //*********************************************************************************************************************************
      //*********************************** Plotting Mean mass for data and SM vs PDG value *********************************************
      //*********************************************************************************************************************************
      cout<<"Plotting Mean mass for EDC and SM vs PDG value"<<endl;
      TCanvas *canvasMassPDG = new TCanvas("canvasMassPDG","",200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings(canvasMassPDG, 0.08, 0.02, 0.055, 0.08);
      canvasMassPDG->SetLogx(1);
      canvasMassPDG->SetLogy(0);

      histoDummyMeanMassVsPDG = new TH2F("histoDummyMeanMassVsPDG", "histoDummyMeanMassVsPDG", 11000, fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5, 1000, minMass, maxMass);
      SetStyleHistoTH2ForGraphs(histoDummyMeanMassVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (SM/EDC)} #GT / M_{#pi^{0} (PDG)}", 0.035, 0.043, 0.035, 0.043, 0.82, 0.9);
      histoDummyMeanMassVsPDG->GetXaxis()->SetMoreLogLabels();
      histoDummyMeanMassVsPDG->GetXaxis()->SetLabelOffset(-0.01);
      histoDummyMeanMassVsPDG->DrawCopy("");

      Double_t standardrangeExponent0 =   -0.5;
      Double_t standardrangeExponent1 =   -0.08;
      Double_t rangeMult0             =   -0.2;
      Double_t rangeMult1             =   -0.001;
      Double_t rangeExponent[2]       =   {standardrangeExponent0, standardrangeExponent1};
      Double_t rangeMult[2]           =   {rangeMult0, rangeMult1};

      legendPDG = GetAndSetLegend2(0.15, 0.95, 0.95, 0.99, 0.043, 2, "", 42);

      // create scaled mass vs pt histos
      histAllResultsVsPDG =  (TH1D*)histAllResults->Clone("Mean mass EDC / mass PDG Pi0");
      histAllResultsVsPDG->Scale(1/massPi0);
      SetStyleHistoTH1ForGraphs(histAllResultsVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (EDC)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
      if(mode == 14) {
        histAllResultsVsPDG->Multiply(histAllResultsVsPDG,histAllResultsVsPDG,1.,1.,"B");
        histAllResultsVsPDG->SetXTitle("#it{E}_{Cluster} (GeV)");
        histAllResultsVsPDG->SetYTitle("#LT M^{2}_{#pi^{0} (EDC)} #GT / M^{2}_{#pi^{0} (PDG)}");
      }
      DrawGammaSetMarker(histAllResultsVsPDG, markerStyle[0], 1, color[0], color[0]);

      histSMResultsVsPDG =  (TH1D*)histSMResults[iSM]->Clone("Mean mass MC / mass PDG Pi0");
      histSMResultsVsPDG->Scale(1/massPi0);
      SetStyleHistoTH1ForGraphs(histSMResultsVsPDG, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (MC)} #GT / M_{#pi^{0} (PDG)}",0.035,0.043, 0.035,0.043, 1.,1.);
      if(mode == 14) {
        histSMResultsVsPDG->Multiply(histSMResultsVsPDG,histSMResultsVsPDG,1.,1.,"B");
        histSMResultsVsPDG->SetXTitle("#it{E}_{Cluster} (GeV)");
        histSMResultsVsPDG->SetYTitle("#LT M^{2}_{#pi^{0} (MC)} #GT / M^{2}_{#pi^{0} (PDG)}");
      }
      DrawGammaSetMarker(histSMResultsVsPDG, markerStyle[1], 1, color[1], color[1]);

      // fitting data mass positions
      FittingFunction="[0]";
      fitMassDataVsPDGConst  = new TF1("fitMassDataVsPDGConst", FittingFunction ,rangeHighPtFitMass[0],rangeHighPtFitMass[1]);
      histAllResultsVsPDG->Fit(fitMassDataVsPDGConst,"QRME0");
      cout << WriteParameterToFile(fitMassDataVsPDGConst) << endl;


      FittingFunction="[0] + [1]*pow(x,[2])";
      fitMassDataVsPDG       = new TF1("fitMassDataVsPDG", FittingFunction ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
      fitMassDataVsPDG->SetParLimits(1, rangeMult[0], rangeMult[1]);
      fitMassDataVsPDG->SetParLimits(2, rangeExponent[0], rangeExponent[1]);

      histAllResultsVsPDG->Fit(fitMassDataVsPDG,"QRME0");
      cout << WriteParameterToFile(fitMassDataVsPDG) << endl;

      FittingFunction="[0]-[3]*TMath::Exp(-[1]*x+[2])";

      fitMassDataVsPDG2      = new TF1("fitMassDataVsPDG2", FittingFunction ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
      fitMassDataVsPDG2->SetParameter(0, fitMassDataVsPDGConst->GetParameter(0));
      if ( (mode ==  14) && (optionEnergy.Contains("PbPb") || optionEnergy.Contains("XeXe")  ) ){
          fitMassDataVsPDG2->SetParLimits(0, fitMassDataVsPDGConst->GetParameter(0)-0.1*fitMassDataVsPDGConst->GetParError(0), fitMassDataVsPDGConst->GetParameter(0)+0.1*fitMassDataVsPDGConst->GetParError(0));
      } else if (mode == 14){
          fitMassDataVsPDG2->SetParLimits(0, fitMassDataVsPDGConst->GetParameter(0)-3*fitMassDataVsPDGConst->GetParError(0), fitMassDataVsPDGConst->GetParameter(0)+3*fitMassDataVsPDGConst->GetParError(0));
      } else {
        fitMassDataVsPDG2->SetParLimits(0, fitMassDataVsPDGConst->GetParameter(0)-0.5*fitMassDataVsPDGConst->GetParError(0), fitMassDataVsPDGConst->GetParameter(0)+0.5*fitMassDataVsPDGConst->GetParError(0));
      }
      fitMassDataVsPDG2->FixParameter(3,1.);

      histAllResultsVsPDG->Fit(fitMassDataVsPDG2,"QRME0");

      // fitting MC mass positions
      FittingFunction="[0] + [1]*pow(x,[2])";
      fitMassMCVsPDG = new TF1("fitMassMCVsPDG", FittingFunction ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
      histSMResultsVsPDG->Fit(fitMassMCVsPDG,"QRME0");
      cout << WriteParameterToFile(fitMassMCVsPDG) << endl;

      FittingFunction="[0]";
      fitMassMCVsPDGConst  = new TF1("fitMassMCVsPDGConst", FittingFunction ,rangeHighPtFitMass[2],rangeHighPtFitMass[3]);
      histSMResultsVsPDG->Fit(fitMassMCVsPDGConst,"QRME0");
      cout << WriteParameterToFile(fitMassMCVsPDGConst) << endl;

      FittingFunction="[0]-[3]*TMath::Exp(-[1]*x+[2])";
      fitMassMCVsPDG2 = new TF1("fitMassMCVsPDG2", FittingFunction ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
      fitMassMCVsPDG2->SetParameter(0, fitMassMCVsPDGConst->GetParameter(0));
      if ( ( mode == 14) && (optionEnergy.Contains("PbPb") || optionEnergy.Contains("XeXe")  ) ){
          fitMassMCVsPDG2->SetParLimits(0, fitMassMCVsPDGConst->GetParameter(0)-0.1*fitMassMCVsPDGConst->GetParError(0), fitMassMCVsPDGConst->GetParameter(0)+0.1*fitMassMCVsPDGConst->GetParError(0));
      } else if ( mode == 14){
          fitMassMCVsPDG2->SetParLimits(0, fitMassMCVsPDGConst->GetParameter(0)-2*fitMassMCVsPDGConst->GetParError(0), fitMassMCVsPDGConst->GetParameter(0)+2*fitMassMCVsPDGConst->GetParError(0));
      } else {
          fitMassMCVsPDG2->SetParLimits(0, fitMassMCVsPDGConst->GetParameter(0)-0.5*fitMassMCVsPDGConst->GetParError(0), fitMassMCVsPDGConst->GetParameter(0)+0.5*fitMassMCVsPDGConst->GetParError(0));
          cout << fitMassMCVsPDGConst->GetParameter(0)-0.5*fitMassMCVsPDGConst->GetParError(0) << " " << fitMassMCVsPDGConst->GetParameter(0)+0.5*fitMassMCVsPDGConst->GetParError(0) << endl;
      }
      fitMassMCVsPDG2->FixParameter(3,1.);


      histSMResultsVsPDG->Fit(fitMassMCVsPDG2,"QRME0");
      cout << WriteParameterToFile(fitMassMCVsPDG2) << endl;


      // draw data graphs and fits
      DrawGammaSetMarkerTF1( fitMassDataVsPDGConst, 7, 2, kGray);
      fitMassDataVsPDGConst->Draw("same");
      DrawGammaSetMarkerTF1( fitMassDataVsPDG, 1, 2, kBlack);
      fitMassDataVsPDG->Draw("same");
      DrawGammaSetMarkerTF1( fitMassDataVsPDG2, 7, 2, kGray+1);
      fitMassDataVsPDG2->Draw("same");
      histAllResultsVsPDG->DrawCopy("same");
      legendPDG->AddEntry(histAllResultsVsPDG,"EDC");

      DrawGammaSetMarkerTF1( fitMassMCVsPDG, 1, 2, kRed+2);
      fitMassMCVsPDG->Draw("same");
      DrawGammaSetMarkerTF1( fitMassMCVsPDGConst, 7, 2, kRed-8);
      fitMassMCVsPDGConst->Draw("same");
      DrawGammaSetMarkerTF1( fitMassMCVsPDG2, 7, 2, kRed-6);
      fitMassMCVsPDG2->Draw("same");
      histSMResultsVsPDG->DrawCopy("same");
      legendPDG->AddEntry(histSMResultsVsPDG,Form("SM_%i",iSM));

      PutProcessLabelAndEnergyOnPlot(0.94, 0.915, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 42, 0.03, "", 1, 1.25, 31);
      for (Int_t i = 0; i < nSets; i++){
         PutProcessLabelAndEnergyOnPlot(0.12, 0.915-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "", 42, 0.03, "", 1, 1.25, 11);
      }

      legendFits   = GetAndSetLegend2(0.12, 0.12 , 0.37, 0.12 + 3*0.03, 0.03, 2, "", 42, 0.35);
      legendFits->AddEntry(fitMassDataVsPDG, " ", "l");
      legendFits->AddEntry(fitMassMCVsPDG, "powerlaw fit","l" );
      legendFits->AddEntry(fitMassDataVsPDG2, " ", "l");
      legendFits->AddEntry(fitMassMCVsPDG2, "exponential fit","l" );
      legendFits->AddEntry(fitMassDataVsPDGConst, " ", "l");
      legendFits->AddEntry(fitMassMCVsPDGConst, "high #it{p}_{T} const.","l" );
      legendFits->Draw("same");

      legendPDG->Draw("same");
      canvasMassPDG->Update();
      canvasMassPDG->SaveAs(Form("%s/MeanMass_Pi0_SM%i_%s.%s", outputDirSampleSummary.Data(),iSM, select.Data(), suffix.Data()));
      canvasMassPDG->Clear();
      delete canvasMassPDG;

      //*********************************************************************************************************************************
      //****************************** Fitting ratio of mean mass position in MC/data ***************************************************
      //*********************************************************************************************************************************
      fFitConst          = new TF1("DataMCConst", "[0]" ,rangeHighPtFitRatio[0],rangeHighPtFitRatio[1]);
      fFitConstFull      = new TF1("ConstFullPtRange", "[0]" ,fBinsPt[ptBinRange[0]], fBinsPt[ptBinRange[1]]);
      histSMResultsRatio->Fit(fFitConst,"QRME0");
      histSMResultsRatio->Fit(fFitConstFull,"QRME0");
      cout << "Offset SM" << iSM << ": " << fFitConstFull->GetParameter(0) << endl;
      Double_t highPtConst            = fixedOffSet;
      if (highPtConst == -1){
          highPtConst = fFitConst->GetParameter(0);
      }

      fFitConstFullInv   = new TF1("ConstFullPtRangeInverted", "[0]" ,fBinsPt[ptBinRange[0]], fBinsPt[ptBinRange[1]]);
      fFitConstFullInv->SetParameter(0, 1./fFitConstFull->GetParameter(0));

      // creating real fit functions
      if(!onlyConstFit){
        fFitMassPos = FitDataMC(histSMResultsRatio, fBinsPt[ptBinRange[0]], fBinsPt[ptBinRange[1]], select,  highPtConst, mode);
        fFitMassPosInverted        =   new TF1("fFitMassPosInverted", "1./([0]+exp([1]+([2]*x)))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]-2]);
        fFitMassPosInverted->SetParameter(0, fFitMassPos->GetParameter(0));
        fFitMassPosInverted->SetParameter(1, fFitMassPos->GetParameter(1));
        fFitMassPosInverted->SetParameter(2, fFitMassPos->GetParameter(2));

        fFitCompositFitted         = new TF1("fFitCompositFitted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]-2]);
        fFitCompositFitted->SetParameter(0, fitMassMCVsPDG->GetParameter(0) );
        fFitCompositFitted->SetParameter(1, fitMassMCVsPDG->GetParameter(1) );
        fFitCompositFitted->FixParameter(2, fitMassMCVsPDG->GetParameter(2) );
        fFitCompositFitted->SetParameter(3, fitMassDataVsPDG->GetParameter(0) );
        fFitCompositFitted->SetParameter(4, fitMassDataVsPDG->GetParameter(1) );
        fFitCompositFitted->FixParameter(5, fitMassDataVsPDG->GetParameter(2) );
        histSMResults[iSM]->Fit(fFitCompositFitted,"QRME0");


        // calculating fit functions based on mass fits with powerlaw like mass functions
        fFitComposit               = new TF1("fFitComposit", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
        fFitComposit->SetParameter(0, fitMassMCVsPDG->GetParameter(0) );
        fFitComposit->SetParameter(1, fitMassMCVsPDG->GetParameter(1) );
        fFitComposit->SetParameter(2, fitMassMCVsPDG->GetParameter(2) );
        fFitComposit->SetParameter(3, fitMassDataVsPDG->GetParameter(0) );
        fFitComposit->SetParameter(4, fitMassDataVsPDG->GetParameter(1) );
        fFitComposit->SetParameter(5, fitMassDataVsPDG->GetParameter(2) );

        fFitCompositInverted       = new TF1("fFitCompositInverted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
        fFitCompositInverted->SetParameter(0, fitMassDataVsPDG->GetParameter(0) );
        fFitCompositInverted->SetParameter(1, fitMassDataVsPDG->GetParameter(1) );
        fFitCompositInverted->SetParameter(2, fitMassDataVsPDG->GetParameter(2) );
        fFitCompositInverted->SetParameter(3, fitMassMCVsPDG->GetParameter(0) );
        fFitCompositInverted->SetParameter(4, fitMassMCVsPDG->GetParameter(1) );
        fFitCompositInverted->SetParameter(5, fitMassMCVsPDG->GetParameter(2) );

        fFitCompositInvertedFitted = new TF1("fFitCompositInvertedFitted", "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
        fFitCompositInvertedFitted->SetParameter(0, fFitCompositFitted->GetParameter(3) );
        fFitCompositInvertedFitted->SetParameter(1, fFitCompositFitted->GetParameter(4) );
        fFitCompositInvertedFitted->SetParameter(2, fFitCompositFitted->GetParameter(5) );
        fFitCompositInvertedFitted->SetParameter(3, fFitCompositFitted->GetParameter(0) );
        fFitCompositInvertedFitted->SetParameter(4, fFitCompositFitted->GetParameter(1) );
        fFitCompositInvertedFitted->SetParameter(5, fFitCompositFitted->GetParameter(2) );


        // calculating fit functions with exponential like mass functions
        fFitExpComb = new TF1("fFitExpComb", "([0]-[6]*TMath::Exp(-[1]*x+[2]))/([3]-[7]*TMath::Exp(-[4]*x+[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
        fFitExpComb->SetParameter(0, fitMassMCVsPDG2->GetParameter(0) );
        fFitExpComb->SetParameter(1, fitMassMCVsPDG2->GetParameter(1) );
        fFitExpComb->SetParameter(2, fitMassMCVsPDG2->GetParameter(2) );
        fFitExpComb->SetParameter(6, fitMassMCVsPDG2->GetParameter(3) );
        fFitExpComb->SetParameter(3, fitMassDataVsPDG2->GetParameter(0) );
        fFitExpComb->SetParameter(4, fitMassDataVsPDG2->GetParameter(1) );
        fFitExpComb->SetParameter(5, fitMassDataVsPDG2->GetParameter(2) );
        fFitExpComb->SetParameter(7, fitMassDataVsPDG2->GetParameter(3) );

        // calculating fit functions with exponential like mass functions
        fFitExpCombInverted = new TF1("fFitExpCombInverted", "([0]-[6]*TMath::Exp(-[1]*x+[2]))/([3]-[7]*TMath::Exp(-[4]*x+[5]))" ,fBinsPt[ptBinRange[0]],fBinsPt[ptBinRange[1]]);
        fFitExpCombInverted->SetParameter(3, fitMassMCVsPDG2->GetParameter(0) );
        fFitExpCombInverted->SetParameter(4, fitMassMCVsPDG2->GetParameter(1) );
        fFitExpCombInverted->SetParameter(5, fitMassMCVsPDG2->GetParameter(2) );
        fFitExpCombInverted->SetParameter(7, fitMassMCVsPDG2->GetParameter(3) );
        fFitExpCombInverted->SetParameter(0, fitMassDataVsPDG2->GetParameter(0) );
        fFitExpCombInverted->SetParameter(1, fitMassDataVsPDG2->GetParameter(1) );
        fFitExpCombInverted->SetParameter(2, fitMassDataVsPDG2->GetParameter(2) );
        fFitExpCombInverted->SetParameter(6, fitMassDataVsPDG2->GetParameter(3) );
      }
        fLogOffsets << "SM" << iSM << ":" << fFitConstFull->GetParameter(0) << endl;
        // fLogOffsets << WriteParameterToFile(fFitConstFull) << endl;

        fstream fLog;
        fLog.open(Form("%s/CorrectCaloNonLinearity_%s_SM%i.log",outputDirSampleSummary.Data(),select.Data(), iSM), ios::out);
        fLog << "use non-inverted with =/ in AliCaloPhotonCuts, and the inverted with *=" << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "Fit EDC SM results ([0]), AliCaloPhotonCuts::FunctionNL_const, use with /=):" << endl;
        fLog << WriteParameterToFile(fFitConst) << endl;
        fLog << WriteParameterToFile(fFitConstFull) << endl;
        fLog << WriteParameterToFile(fFitConstFullInv) << endl;
        fLog << "-----------------------------------" << endl;
      if(!onlyConstFit){
        fLog << "Exponential function fitted" << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "FitDataMC results ([0]+[3]*exp([1]+([2]*x)), AliCaloPhotonCuts::FunctionNL_kSDM, use with /=):" << endl;
        for(Int_t i=0;i<=3;i++) fLog << "Par " << i << ": " << fFitMassPos->GetParameter(i) << " +- " << fFitMassPos->GetParError(i) << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "Ind. Mass fitted with powerlaws" << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "([0] + [1]*pow(x,[2]))/([3] + [4]*pow(x,[5]))" << endl;
        fLog << "AliCaloPhotonCuts::FunctionNL_DPOW" << endl;
        fLog << WriteParameterToFile(fFitComposit) << endl;
        fLog << WriteParameterToFile(fFitCompositFitted) << endl;
        fLog << WriteParameterToFile(fFitCompositInverted) << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "Ind. Mass fitted with exponentials" << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "([0]-[3]*TMath::Exp(-[1]*x+[2]))/([4]-[7]*TMath::Exp(-[5]*x+[6]))" << endl;
        fLog << "AliCaloPhotonCuts::FunctionNL_DExp" << endl;
        fLog << WriteParameterToFile(fFitExpComb) << endl;
        fLog << WriteParameterToFile(fFitExpCombInverted) << endl;
        fLog << "-----------------------------------" << endl;
        fLog << "-----------------------------------" << endl;
        fLog.close();
      }

      //*******************************************************************************
      // plotting mass ratios
      //*******************************************************************************
      canvasMassRatioMCData = new TCanvas("canvasMassRatioMCData","",200,10,1350,900);  // gives the page size
      DrawGammaCanvasSettings(canvasMassRatioMCData, 0.082, 0.012, 0.02, 0.1);
      canvasMassRatioMCData->SetLogx(1);
      canvasMassRatioMCData->SetLogy(0);


      histoDummyDataMCRatio = new TH2F("histoDummyDataMCRatio","histoDummyDataMCRatio", 11000, fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5, 1000, minPlotY, maxPlotY);
      if(mode == 14) SetStyleHistoTH2ForGraphs(histoDummyDataMCRatio, "#it{E}_{Cluster} (GeV)","#LT M^{2}_{#pi^{0} (SM)} #GT / #LT M^{2}_{#pi^{0} (mean)} #GT", 0.035, 0.043, 0.035, 0.043, 0.82, 0.9);
      else SetStyleHistoTH2ForGraphs(histoDummyDataMCRatio, "#it{E}_{Cluster} (GeV)","#LT M_{#pi^{0} (SM)} #GT / #LT M_{#pi^{0} (mean)} #GT", 0.035, 0.043, 0.035, 0.043, 0.82, 0.9);
      histoDummyDataMCRatio->GetXaxis()->SetMoreLogLabels();
      histoDummyDataMCRatio->GetXaxis()->SetLabelOffset(-0.01);
      histoDummyDataMCRatio->DrawCopy("");


      Int_t nCorrections      = 3;
      DrawGammaSetMarkerTF1( fFitConstFull, 9, 2, 807);
      fFitConstFull->Draw("same");
      nCorrections++;

      if(!onlyConstFit){
        fFitComposit->SetRange(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5);
        fFitExpComb->SetRange(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5);
        fFitMassPos->SetRange(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5);

        DrawGammaSetMarkerTF1( fFitComposit, 7, 2, kGreen+2);
        DrawGammaSetMarkerTF1( fFitExpComb, 8, 2, kBlue+2);
        DrawGammaSetMarkerTF1( fFitMassPos, 1, 2, kRed+1);

        fFitMassPos->Draw("same");
        fFitComposit->Draw("same");
        fFitExpComb->Draw("same");
      }
      DrawGammaSetMarker(histSMResultsRatio, 20, 1, color[0], color[0]);
      histSMResultsRatio->DrawCopy("same,pe");

      PutProcessLabelAndEnergyOnPlot(0.94, 0.95, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 42, 0.03, "", 1, 1.25, 31);
      PutProcessLabelAndEnergyOnPlot(0.62, 0.95, 0.03, Form("Supermodule: %i",iSM), "", "", 42, 0.038, "", 1, 1.25, 31);
      for (Int_t i = 0; i < nSets; i++){
         PutProcessLabelAndEnergyOnPlot(0.12, 0.945-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "", 42, 0.03, "", 1, 1.25, 11);
      }

      legendCorrectionFunctions = GetAndSetLegend2(0.125,0.15, 0.4,0.15+nCorrections*1.1*0.03, 0.03, 1, "", 42, 0.15);
      if(!onlyConstFit){
        legendCorrectionFunctions->AddEntry(fFitMassPos,"Exponential function fitted","l");
        legendCorrectionFunctions->AddEntry(fFitComposit,"Ind. Mass fitted with powerlaws","l");
        legendCorrectionFunctions->AddEntry(fFitExpComb,"Ind. Mass fitted with exponentials","l");
      }
      legendCorrectionFunctions->AddEntry(fFitConstFull,"Constant fitted","l");
      legendCorrectionFunctions->Draw();

      DrawGammaLines(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5, 1.0, 1.0, 1, kGray+2, 2);

      canvasMassRatioMCData->Update();
      canvasMassRatioMCData->SaveAs(Form("%s/MeanMassRatio_%i_%s.%s", outputDirSampleSummary.Data(), iSM, select.Data(), suffix.Data()));
      cout<<"saving plot for SM: "<<iSM<<endl;
      canvasMassRatioMCData->Clear();

      //*******************************************************************************
      // plotting mass ratios for all SM
      //*******************************************************************************
      canvas_All->cd();

      if(isFirstLoop){
        histoDummyDataMCRatio->DrawCopy("");
        PutProcessLabelAndEnergyOnPlot(0.75, 0.96, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 42, 0.03, "", 1, 1.25, 31);
        DrawGammaLines(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5, 1.0, 1.0, 1, kGray+2, 2);
        for (Int_t i = 0; i < nSets; i++){
           PutProcessLabelAndEnergyOnPlot(0.12, 0.945-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "", 42, 0.03, "", 1, 1.25, 11);
        }
      }
      legend_All->AddEntry(histSMResults[iSM],Form("%s",dataMC[iSM + 1].Data()),"p");
      DrawGammaSetMarker(histSMResultsRatio, 20, 1, color[iSM], color[iSM]);
      histSMResultsRatio->DrawCopy("same");
      if(iSM == 0) legend_All->Draw("same");
      canvas_All->Update();
      canvas_All->SaveAs(Form("%s/MeanMassRatio_All.%s",outputDir.Data(),suffix.Data()));


      //*******************************************************************************
      // plotting total correction
      //*******************************************************************************
      canvasMassRatioMCData->cd();
      totalCorrection = new TH1D("Total Correction","; #it{E}_{Cluster} (GeV); correction factor",1000,0.3,50);
      SetStyleHistoTH1ForGraphs(totalCorrection, "#it{E}_{Cluster} (GeV)","correction factor",0.035,0.043, 0.035,0.043, 1.,0.9);
      totalCorrection->GetYaxis()->SetRangeUser(minPlotY+0.031,1.1);
      totalCorrection->GetXaxis()->SetMoreLogLabels();
      totalCorrection->GetXaxis()->SetLabelOffset(-0.01);
      SetLogBinningXTH(totalCorrection);
      totalCorrection->DrawCopy("p");

      if(!onlyConstFit){
        fFitMassPosInverted->SetRange(0.3,50);
        fFitCompositInverted->SetRange(0.3,50);
        fFitExpCombInverted->SetRange(0.3,50);
      }
      fFitConstFullInv->SetRange(0.3,50);

      if(!onlyConstFit){
        DrawGammaSetMarkerTF1( fFitMassPosInverted, 1, 2, kRed+2);
        DrawGammaSetMarkerTF1( fFitCompositInverted, 7, 2, kGreen+2);
        DrawGammaSetMarkerTF1( fFitExpCombInverted, 8, 2, kBlue+2);
      }
      DrawGammaSetMarkerTF1( fFitConstFullInv, 9, 2, 807);
      fFitConstFullInv->Draw("same");
      if(!onlyConstFit){
        fFitMassPosInverted->Draw("same");
        fFitCompositInverted->Draw("same");
        fFitExpCombInverted->Draw("same");
      }


      legend2 = GetAndSetLegend2(0.45,0.15, 0.725,0.15+nCorrections*1.1*0.03, 0.03, 1, "", 42, 0.15);//GetAndSetLegend2(0.2, 0.2, 0.4, 0.29, 0.03, 1, "", 42);
      if(!onlyConstFit){
        legend2->AddEntry(fFitMassPosInverted,"Correction factor for MC, exponential function (ratio fit)","l");
        legend2->AddEntry(fFitCompositInverted,"Correction factor for MC from mass fits (powerlaws)","l");
        legend2->AddEntry(fFitExpCombInverted,"Correction factor for MC from mass fits (exponentials)","l");
      }
      legend2->AddEntry(fFitConstFullInv,"Constant fitted","l");
      legend2->Draw("same");

      DrawGammaLines(0.3, 50.,1.0, 1.0, 1, kGray+2, 2);

      PutProcessLabelAndEnergyOnPlot(0.94, 0.96, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 42, 0.03, "", 1, 1.25, 31);
      for (Int_t i = 0; i < nSets; i++){
         PutProcessLabelAndEnergyOnPlot(0.12, 0.945-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "", 42, 0.03, "", 1, 1.25, 11);
      }

      canvasMassRatioMCData->Update();
      canvasMassRatioMCData->SaveAs(Form("%s/TotalCorrection_%s_SM%i.%s", outputDir.Data(), select.Data(), iSM, suffix.Data()));
      canvasMassRatioMCData->Clear();
      delete canvasMassRatioMCData;

      //*******************************************************************************
      // plotting mass ratios for all SM
      //*******************************************************************************
      canvas_All_Fitted->cd();

      if(isFirstLoop){
        histoDummyDataMCRatio->GetYaxis()->SetRangeUser(0.971,1.049);
        histoDummyDataMCRatio->DrawCopy("");
        PutProcessLabelAndEnergyOnPlot(0.75, 0.96, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 42, 0.03, "", 1, 1.25, 31);
        DrawGammaLines(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5, 1.0, 1.0, 1, kGray+2, 2);
        for (Int_t i = 0; i < nSets; i++){
           PutProcessLabelAndEnergyOnPlot(0.12, 0.945-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "", 42, 0.03, "", 1, 1.25, 11);
        }
      }
      fFitConstFull->SetRange(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5);
      DrawGammaSetMarkerTF1( fFitConstFull, 1, 1, color[iSM+1]);
      fFitConstFull->Draw("same");
    //   legend_All->AddEntry(histSMResults[iSM],Form("%s",dataMC[iSM + 1].Data()),"p");
      DrawGammaSetMarker(histSMResultsRatio, 24, 1, color[iSM+1], color[iSM+1]);
      histSMResultsRatio->DrawCopy("same");
      if(iSM == 0) legend_All->Draw("same");
      canvas_All_Fitted->Update();
      canvas_All_Fitted->SaveAs(Form("%s/MeanMassRatio_All_Fitted.%s",outputDir.Data(),suffix.Data()));
      //*******************************************************************************
      // plotting mass ratio fits for all SM
      //*******************************************************************************
      canvas_All_FitOnly->cd();

      if(isFirstLoop){
        histoDummyDataMCRatio->GetYaxis()->SetRangeUser(0.991,1.019);
        histoDummyDataMCRatio->DrawCopy("");
        PutProcessLabelAndEnergyOnPlot(0.75, 0.96, 0.03, fCollisionSystem.Data(), fTextMeasurement.Data(), recGamma.Data(), 42, 0.03, "", 1, 1.25, 31);
        DrawGammaLines(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5, 1.0, 1.0, 1, kGray+2, 2);
        for (Int_t i = 0; i < nSets; i++){
           PutProcessLabelAndEnergyOnPlot(0.12, 0.945-2*0.03*(i), 0.03, fPlotLabelsRatio[i].Data(),"", "", 42, 0.03, "", 1, 1.25, 11);
        }
      }
      fFitConstFull->SetRange(fBinsPt[ptBinRange[0]]/1.5, fBinsPt[ptBinRange[1]]*1.5);
      DrawGammaSetMarkerTF1( fFitConstFull, 2, 2, color[iSM+1]);
      fFitConstFull->Draw("same");
      if(iSM == 0) legend_All->Draw("same");
      canvas_All_FitOnly->Update();
      canvas_All_FitOnly->SaveAs(Form("%s/MeanMassRatio_All_FitOnly.%s",outputDir.Data(),suffix.Data()));

      cout << "-----------------------------------------------------" << endl;
      cout << "-----------------------------------------------------" << endl;
      fOutput->cd();
      histAllResults->Write("Mean mass Data");
      histSMResults[iSM]->Write(Form("Mean mass SM %i",iSM));
      histAllResultsVsPDG->Write();
      histSMResultsVsPDG->Write();
      histSMResults[iSM]->Write("MeanMassRatioMCData");
      if(!onlyConstFit){
        fFitMassPosInverted->Write("REXP_TotalCorr");
        fFitCompositInverted->Write("DPOW_TotalCorr");
        fFitExpCombInverted->Write("DEXP_TotalCorr");
      }
      fFitConstFullInv->Write("CONST_TotalCorr");

      fOutput->Write();
      fOutput->Close();
      isFirstLoop = kFALSE;
    }
    fLogOffsets.close();
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
TF1* FitDataMC(TH1* fHisto, Double_t minFit, Double_t maxFit, TString selection, Double_t constPar, Int_t mode){

    cout << "running standard fit from " <<  minFit << "\t"<<  maxFit << endl;
    TF1* fFitReco = new TF1("DataMC", "[0]+[3]*exp([1]+([2]*x))" ,minFit,maxFit);

    fFitReco->SetParameter(0,1.);
    fFitReco->SetParameter(2,-0.5);
    fFitReco->SetParameter(1,-1.);
    fFitReco->FixParameter(3, 1.);

    // if(constPar!=-1 && (mode != 3 || mode != 5 || mode != 13)) fFitReco->FixParameter(0,constPar);
    fHisto->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(kRed);
    fFitReco->SetLineWidth(2);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
        cout << "Parameters for DataMC: " << endl;
        for(Int_t i=0;i<=3;i++) cout << "Par " << i << ": " << fFitReco->GetParameter(i) << " +- " << fFitReco->GetParError(i) << endl;
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
    if (mode == 14){
        if (ptcenter > 1.5)
            mesonAmplitudeMin = mesonAmplitude*90./100.;
        else
            mesonAmplitudeMin = mesonAmplitude*85./100.;
    // special setting for PCM-PHOS
    } else if (mode == 15){
        mesonAmplitudeMin = mesonAmplitude*20./100.;
        fitRangeMin = 0.08;
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
    if (mode == 15){
        fFitReco->SetParLimits(1, fMesonMassExpect*0.65, fMesonMassExpect*1.45);
    } else {
        fFitReco->SetParLimits(1, fMesonMassExpect*0.9, fMesonMassExpect*1.1);
    }
    fFitReco->SetParLimits(2, 0.001, 0.1);
    fFitReco->SetParLimits(3, 0.001, 0.09);

    histo->Fit(fFitReco,"QRME0");

    fFitReco->SetLineColor(kRed+1);
    fFitReco->SetLineWidth(1);
    fFitReco->SetLineStyle(1);

    if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
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

    if(TString(gMinuit->fCstatu.Data()).Contains("CONVERGED") == 1 || TString(gMinuit->fCstatu.Data()).Contains("SUCCESSFUL") == 1 || TString(gMinuit->fCstatu.Data()).Contains("PROBLEMS") == 1){
        cout << "Parameters for FitRecursiveGaussian: " << endl;
        for(Int_t i=0;i<=2;i++) cout << "Par " << i << ": " << f1->GetParameter(i) << " +- " << f1->GetParError(i) << endl;
    } else {
        cout << "FitRecursiveGaussian fitting failed with status " << gMinuit->fCstatu.Data() <<endl << endl;
    }

    return f1;
}
