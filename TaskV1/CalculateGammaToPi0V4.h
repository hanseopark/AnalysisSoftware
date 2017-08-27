#ifndef GAMMACONV_CalculateGammaToPi0
#define GAMMACONV_CalculateGammaToPi0

    #include "../CommonHeaders/ExtractSignalBinning.h"

    //***********************************************************************************************
    //**************************  variables for global configuration ********************************
    //***********************************************************************************************
    TString nameRec;
    TString nameOutputLabel;

    Int_t centCutNumberI                                        = 0;
    TString centrality                                          = "";
    Double_t fNcoll                                             = 0;
    Double_t fNcollErr                                          = 0;
    TString triggerCutNumber                                    = "";
    TString subTriggerCutNumber                                 = "";

    Bool_t kDoPileup                                            = kFALSE;
    Double_t textSizeSpectra                                    = 0.035;

    Double_t fitMinPt                                           = 0.4;
    Double_t fitMaxPt                                           = 0;
    TString forOutput                                           = "";
    TString nameIntRanges[6]                                    = {"","Wide", "Narrow", "Left", "LeftWide", "LeftNarrow"};
    TString fEventCutSelection                                  = "";
    TString fGammaCutSelection                                  = "";
    TString fClusterCutSelection                                = "";
    TString fElectronCutSelection                               = "";
    TString fMesonCutSelection                                  = "";

    //***********************************************************************************************
    //********************************  Inputs for pi0 spectrum  ************************************
    //***********************************************************************************************
    TFile *filePi0                                              = NULL;
    TH1D *histoCorrectedPi0Yield[6]                             = {NULL, NULL, NULL, NULL, NULL, NULL};
    TGraphAsymmErrors* graphSysPi0PileUpOptions                 = NULL;
    TGraphAsymmErrors* graphSysPi0PileUpIterations              = NULL;
    TH1D *histoMCYieldMeson                                     = NULL;
    TH1D *histoMCYieldMesonOldBin                               = NULL;
    TFile *fileEta                                              = NULL;
    TH1D *histoCorrectedEtaYield                                = NULL;

    //***********************************************************************************************
    //*******************************  Different fits for pi0 yield *********************************
    //***********************************************************************************************
    TF1 *fitPi0YieldA                                           = NULL;
    TF1 *fitPi0YieldB                                           = NULL;
    TF1 *fitPi0YieldC                                           = NULL;
    TH1D *histoRatioFitPi0A                                     = NULL;
    TH1D *histoRatioFitPi0B                                     = NULL;
    TH1D *histoRatioFitPi0C                                     = NULL;

    //***********************************************************************************************
    //***************************  input for gamma spectrum  ****************************************
    //***********************************************************************************************
    TFile *fileGamma                                            = NULL;
    TDirectoryFile* directoryGamma                              = NULL;
    TH1D *histoGammaSpecCorrPurity                              = NULL;
    TH1D *histoGammaSpecMCAll                                   = NULL;
    TH1D *histoIncRatioFitPurity[6]                             = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *histoIncRatioPurityTrueEff[6]                         = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *histoMCIncRatio                                       = NULL;
    TH1D *histoMCDecaySumGammaPt                                = NULL;
    //***********************************************************************************************
    //***************************  Different fits to inclusive gamma  *******************************
    //***********************************************************************************************
    TF1 *fitGammaA                                              = NULL;
    TF1 *fitGammaB                                              = NULL;
    TH1D *histoRatioFitGammaA                                   = NULL;
    TH1D *histoRatioFitGammaB                                   = NULL;

    //***********************************************************************************************
    //***************************  cocktail definition **********************************************
    //***********************************************************************************************
    TFile *cocktailFile                                         = NULL;
    TDirectoryFile *cocktailDir                                 = NULL;
    TH1D *cocktailAllGamma                                      = NULL;
    TH1D *cocktailPi0                                           = NULL;
    TH1D *cocktailEta                                           = NULL;
    TH1D *cocktailAllGammaPi0                                   = NULL;
    TH1D *cocktailPi0Rebinned                                   = NULL;
    TH1D *cocktailAllGammaRebinned                              = NULL;

    //***********************************************************************************************
    //***************************  calculated double Ratios  ****************************************
    //***********************************************************************************************
    TH1D *histoDoubleRatioTrueEffPurity[6]                      = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *histoDoubleRatioFitPi0YieldPurity[6]                  = {NULL, NULL, NULL, NULL, NULL, NULL};
    TH1D *histoDoubleRatioUpperLimits                           = NULL;

    //***********************************************************************************************
    //**************************  output file  ******************************************************
    //***********************************************************************************************
    TFile *fileCorrectedOutput                                  = NULL;

    //***********************************************************************************************
    //************* Wrapper function for separated cut string and centrality settings****************
    //**    - sets global variables for the different cut strings                                  **
    //**    - sets trigger cutnumbers                                                              **
    //**    - sets centrality name and number                                                      **
    //**    - sets Ncoll values + errors                                                           **
    //***********************************************************************************************
    void SeparateCutnumberString(TString cutSelStr, Int_t mode, TString energy){
        TString fCutSelection                       = cutSelStr;
        if (mode == 9){
            ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
            fEventCutSelection                      = fGammaCutSelection(0,7);
            fGammaCutSelection                      = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
            cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
        } else {
            ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
        }
        centCutNumberI                     = ((TString)fEventCutSelection(GetEventCentralityMinCutPosition(),1)).Atoi();
        centrality                         = GetCentralityString(fEventCutSelection);
        triggerCutNumber                   = fEventCutSelection(3,1);
        subTriggerCutNumber                = fEventCutSelection(4,1);
        fNcoll                             = GetNCollFromCutNumber(fEventCutSelection,energy);
        fNcollErr                          = GetNCollErrFromCutNumber(fEventCutSelection,energy);
    }

    //***********************************************************************************************
    //************* Wrapper function to create output names and proper mesonType ********************
    //***********************************************************************************************
    void CreateNamingStrings(TString nameMeson, TString isMC){
    nameOutputLabel                      = Form("Gamma_%s",nameMeson.Data());
    if (isMC.CompareTo("kTRUE") == 0){
        nameRec                           = "recMC";
    } else {
        nameRec                           = "data";
    }

    TString mesonType;
    if ((nameMeson.CompareTo("Pi0") == 0) || (nameMeson.CompareTo("Pi0EtaBinning") == 0)){
        mesonType                               = "Pi0";
    } else {
        mesonType                               = "Eta";
    }
    }

    //***********************************************************************************************
    //******** Creates proper labeling with collision system, detection system and particle  ********
    //***********************************************************************************************
    void PlotLatexLegend(   Double_t xPosition        = 0.5,
                            Double_t yPosition        = 0.5,
                            Double_t textSizeSpectra  = 0.04,
                            TString collisionSystem   = "",
                            TString detectionProcess  = "",
                            Int_t labelNums           = 3,
                            Int_t align               = 11
                        )
    {
        TLatex* labelPi0Plot    = NULL;
        TLatex* labelEnergyPlot = NULL;
        TLatex* labelDetProcPlot= NULL;

        if(labelNums<3){
            labelEnergyPlot = new TLatex(xPosition, yPosition+0.85*textSizeSpectra,collisionSystem.Data());
            labelDetProcPlot= new TLatex(xPosition, yPosition,detectionProcess.Data());
        }else{
            labelEnergyPlot = new TLatex(xPosition, yPosition+2*0.85*textSizeSpectra,collisionSystem.Data());
            labelPi0Plot    = new TLatex(xPosition, yPosition+0.85*textSizeSpectra,"#pi^{0} #rightarrow #gamma#gamma");
            labelDetProcPlot= new TLatex(xPosition, yPosition,detectionProcess.Data());
        }
        SetStyleTLatex( labelEnergyPlot, 0.85*textSizeSpectra,4, 1, 42, kTRUE, align);
        if (labelPi0Plot) SetStyleTLatex( labelPi0Plot, 0.85*textSizeSpectra,4, 1, 42, kTRUE, align);
        SetStyleTLatex( labelDetProcPlot, 0.85*textSizeSpectra,4, 1, 42, kTRUE, align);

        labelEnergyPlot->Draw();
        if (labelPi0Plot) labelPi0Plot->Draw();
        labelDetProcPlot->Draw();
    }

#endif