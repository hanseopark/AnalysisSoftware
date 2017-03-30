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


//***********************************************************************************************
//**********************  NLO file reading  *****************************************************
//***********************************************************************************************
TFile* fileTheoryCompilation                                = NULL;
TH1D *cocktailAllGammaNLO                                   = NULL;
TGraphAsymmErrors *graphDirectPhotonNLO                     = NULL;
TGraphAsymmErrors *graphPromptPhotonNLO                     = NULL;
TGraphAsymmErrors *graphFragmentationPhotonNLO              = NULL;
TGraphAsymmErrors* graphDirectPhotonNLOCopy                 = NULL;
TGraphAsymmErrors* graphNLODoubleRatio                      = NULL;
TGraphAsymmErrors* graphNLODirGammaSpectra                  = NULL;
TF1 *fitNLODoubleRatio                                      = NULL;
TF1* fitNLODirectPhoton                                     = NULL;
TF1* fitNLOPromptPhoton                                     = NULL;
TF1* fitNLOFragmentationPhoton                              = NULL;
TH1D *histoRatioNLODirectPhoton                             = NULL;
TF1 *fitCocktailAllGammaForNLO                              = NULL;

//***********************************************************************************************
//********************************  Inputs for pi0 spectrum  ************************************
//***********************************************************************************************
TFile *filePi0                                              = NULL;
TH1D *histoCorrectedPi0Yield[6]                             = {NULL, NULL, NULL, NULL, NULL, NULL};
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
//**************************  calculated direct gamma spectra  **********************************
//***********************************************************************************************
TH1D *histoDirectPhotonSpectrum                             = NULL;
TH1D *histoThermalPhotonSpectrum                            = NULL;
TH1D *histoPromptPhotonSpectrum                             = NULL;

//***********************************************************************************************
//**************************  output file  ******************************************************
//***********************************************************************************************
TFile *fileCorrectedOutput                                  = NULL;


//***********************************************************************************************
//****************** variable setting for the different error graphs ****************************
//***********************************************************************************************
ifstream fileSysErrGamma;
Int_t nPointsGamma                                          = 0;
Double_t relSystErrorGammaUp[50];
Double_t relSystErrorGammaDown[50];
Double_t relSystErrorWOMaterialGammaUp[50];
Double_t relSystErrorWOMaterialGammaDown[50];
Double_t systErrorGammaUp[50];
Double_t systErrorGammaDown[50];

ifstream fileSysErrInclRatio;
Int_t nPointsInclRatio                                      = 0;
Double_t relSystErrorInclRatioUp[50];
Double_t relSystErrorInclRatioDown[50];
Double_t relSystErrorWOMaterialInclRatioUp[50];
Double_t relSystErrorWOMaterialInclRatioDown[50];
Double_t systErrorInclRatioUp[50];
Double_t systErrorInclRatioDown[50];

ifstream fileSysErrDoubleRatio;
Int_t nPointsDoubleRatio                                    = 0;
Double_t relSystErrorDoubleRatioUp[50];
Double_t relSystErrorDoubleRatioDown[50];
Double_t relSystErrorWOMaterialDoubleRatioUp[50];
Double_t relSystErrorWOMaterialDoubleRatioDown[50];
Double_t systErrorDoubleRatioUp[50];
Double_t systErrorDoubleRatioDown[50];

TGraphAsymmErrors*  graphGammaYieldSysErr                   = NULL;
TGraphAsymmErrors*  graphInclRatioSysErr                    = NULL;
TGraphAsymmErrors*  graphDoubleRatioSysErr                  = NULL;
TGraphAsymmErrors*  graphDoubleRatioFitSysErr               = NULL;

//***********************************************************************************************
//************* Wrapper function for separated cut string and centrality settings****************
//**    - sets global variables for the different cut strings                                  **
//**    - sets trigger cutnumbers                                                              **
//**    - sets centrality name and number                                                      **
//**    - sets Ncoll values + errors                                                           **
//***********************************************************************************************
void SeparateCutnumberString(TString cutSelStr, Int_t mode){
    TString fCutSelection                       = cutSelStr;
    TString fEventCutSelection                  = "";
    TString fGammaCutSelection                  = "";
    TString fClusterCutSelection                = "";
    TString fElectronCutSelection               = "";
    TString fMesonCutSelection                  = "";
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
    fNcoll                             = GetNCollFromCutNumber(fEventCutSelection);
    fNcollErr                          = GetNCollErrFromCutNumber(fEventCutSelection);
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
void PlotLatexLegend( Double_t xPosition        = 0.5,
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

//***********************************************************************************************
//**  function to return upper limit on photon excess, using a Bayesian approach               **
//**  with the heaviside function used as prior (excluding R_gamma < 1)                        **
//***********************************************************************************************
Double_t GetUpperLimit( Double_t mean, Double_t statErr, Double_t sysErr,
                        Double_t confidenceLevel, Double_t& confidenceLevelReached, 
                        Double_t accuracy = 1e-9, Int_t maxNIterations = 1e8
                      ) {
    
    // R_gamma limits
    Double_t    minRGamma           = 1.;
    Double_t    maxRGamma           = 10.;
    
    // total uncertainty
    Double_t    sigmaTot            = TMath::Sqrt(statErr*statErr + sysErr*sysErr);
    
    // cond. prob norm
    Double_t    condProbNorm        = TMath::Erf( (mean - 1)/(TMath::Sqrt(2)*sigmaTot) ) + 1;   // - 1 term from limit erf(-inf)
    
    // cond. prob
    TF1         condProb("condProb", Form("[0] * ( TMath::Erf( ([1] - 1)/(TMath::Sqrt(2)*[2]) ) - TMath::Erf( ([1] - x)/(TMath::Sqrt(2)*[2]) ) )"), minRGamma, maxRGamma);
    condProb.SetParameter(0, 1./condProbNorm);
    condProb.SetParameter(1, mean);
    condProb.SetParameter(2, sigmaTot);

    // iteratively find upper limit (interval bisection)
    Double_t    upperLimit          = (maxRGamma-1)/2;
    Double_t    upperLimitPrev      = upperLimit;
    Double_t    step                = 0.;
    Int_t       nIterations         = 0;
    while (((condProb.Eval(upperLimit) < (confidenceLevel-accuracy)) || condProb.Eval(upperLimit) > (confidenceLevel+accuracy)) && nIterations < maxNIterations) {
        
        if (condProb.Eval(upperLimit) > confidenceLevel)
            step                    = - TMath::Abs(upperLimit-1)/2;
        else
            step                    = TMath::Abs(upperLimitPrev-upperLimit)/2;
        upperLimitPrev              = upperLimit;
        upperLimit                  = upperLimit + step;
        
        //if ( !(nIterations%10) ) cout << "   condProb.Eval( " << upperLimit << ") = " << condProb.Eval(upperLimit) << endl;
        
        nIterations++;
    }
    
    confidenceLevelReached          = condProb.Eval(upperLimit);
    return upperLimit;
}

//***********************************************************************************************
//************ Wrapper function to extract the direct photon spectrum upper limit ***************
//**    - extracts the direct photon signal upper limits based on statistical error            **
//**      given in a histogram and systematic errors in a graphs                               **
//**    - confidence level can be varied as well as the accuracy and the number of iterations  **
//***********************************************************************************************
TH1D* GetUpperLimitsHisto(TH1D* histo, 
                          TGraphAsymmErrors* sysErrGraph, 
                          Double_t confidenceLevel          = 0.95, 
                          Double_t accuracy                 = 0.004, 
                          Int_t maxNIterations              = 1e3
                         ) {

    cout << endl;
    cout << "*************************************************************" << endl;
    cout << "**                                                         **" << endl;
    cout << "**      STARTING UPPER LIMIT CALCULATION                   **" << endl;
    cout << "**                                                         **" << endl;
    cout << "*************************************************************" << endl;
    cout << endl;

    // upper limits histo
    TH1D*       upperLimits             = (TH1D*)histo->Clone("upperLimits");

    // get graph quantities
    Int_t       nBinsGraph              = sysErrGraph->GetN();
    Double_t*   xValueGraph             = sysErrGraph->GetX();
    Double_t*   yErrorLowGraph          = sysErrGraph->GetEYlow();
    Double_t*   yErrorHighGraph         = sysErrGraph->GetEYhigh();

    // fill upper limits histo
    for (Int_t i=1; i<histo->GetNbinsX()+1; i++) {
        if (!histo->GetBinContent(i)) {
            upperLimits->SetBinContent( i, 0);
            upperLimits->SetBinError(   i, 0);
        } else {
            for (Int_t j=0; j<sysErrGraph->GetN(); j++) {
                if (xValueGraph[j] == histo->GetBinCenter(i)) {

                    Double_t reached    = 0.;

                    upperLimits->SetBinContent( i, GetUpperLimit(histo->GetBinContent(i),histo->GetBinError(i),(yErrorLowGraph[j]+yErrorHighGraph[j])/2,confidenceLevel,reached,accuracy,maxNIterations));
                    upperLimits->SetBinError(   i, 0);

                    cout << "p_T = " << histo->GetBinCenter(i) << ":\t" << histo->GetBinContent(i) << " ( +/- " << TMath::Sqrt(histo->GetBinError(i)*histo->GetBinError(i) + (yErrorLowGraph[j]+yErrorHighGraph[j])/2*(yErrorLowGraph[j]+yErrorHighGraph[j])/2) << " )\t->\t" << upperLimits->GetBinContent(i) << "\tat CL = " << reached*100 << "%" << endl;
                }
            }
        }
    }

    cout << endl;
    cout << "*************************************************************" << endl;
    cout << "**                                                         **" << endl;
    cout << "**      DONE WITH UPPER LIMIT CALCULATION                  **" << endl;
    cout << "**                                                         **" << endl;
    cout << "*************************************************************" << endl;
    cout << endl;
    
    return upperLimits;
}

//***********************************************************************************************
//*******  functions to calculate uncertainty from fit (variation with parameter errors)  *******
//***********************************************************************************************
TF1* VaryFunctionWithParameterError(TF1* func, Bool_t up = kTRUE) {
    
    Double_t    sign            = 1.;
    if (!up)    sign            = -1.;
    
    Int_t nPar                  = func->GetNpar();
    std::vector<Double_t> parameter(nPar);
    std::vector<Double_t> parameterError(nPar);
    std::vector<Double_t> parameterModified(nPar);
    for (Int_t i=0; i<nPar; i++) {
        parameter[i]            = func->GetParameter(i);
        parameterError[i]       = func->GetParError(i);
        parameterModified[i]    = parameter[i] + sign*parameterError[i];
    }
    
    TF1* returnFunc             = (TF1*)func->Clone(Form("%s_up", func->GetName()));
    for (Int_t i=0; i<nPar; i++) returnFunc->SetParameter(i, parameterModified[i]);
    
    return returnFunc;
}

//***********************************************************************************************
//****************** Calculate systematic errors with the fitted pi0 ****************************
//***********************************************************************************************
TGraphAsymmErrors* ProduceTotalSystematicUncertaintyWithFit( TH1D* doubleRatio, 
                                                             TGraphAsymmErrors* systUncert, 
                                                             TF1* pi0Fit
                                                           ) {
    
    // syst. uncert. graph
    Int_t       nPoints             = systUncert->GetN();
    Double_t*   xVal                = systUncert->GetX();
    Double_t*   xErrUp              = systUncert->GetEXhigh();
    Double_t*   xErrDown            = systUncert->GetEXlow();
    Double_t*   yVal                = systUncert->GetY();
    Double_t*   yErrUp              = systUncert->GetEYhigh();
    Double_t*   yErrDown            = systUncert->GetEYlow();

    // declare total syst. uncert.
    Double_t*   yErrUpTot           = new Double_t[nPoints];
    Double_t*   yErrDownTot         = new Double_t[nPoints];
    
    // vary fit using parameter errors
    TF1*        pi0FitUp            = VaryFunctionWithParameterError(pi0Fit);
    TF1*        pi0FitDown          = VaryFunctionWithParameterError(pi0Fit, kFALSE);
    
    // calculate total syst. uncert.
    Double_t    totSystUncertUp     = 0.;
    Double_t    totSystUncertDown   = 0.;
    Double_t    systUncertFit       = 0.;
    for (Int_t i=0; i<nPoints; i++) {
        for (Int_t j=1; j<doubleRatio->GetNbinsX()+1; j++) {
            
            if (doubleRatio->GetBinCenter(j)!=xVal[i]) continue;
            
            systUncertFit           = TMath::Abs(pi0FitUp->Eval(xVal[i])-pi0FitDown->Eval(xVal[i])) / 2;
            systUncertFit           = systUncertFit / pi0Fit->Eval(xVal[i]);
            systUncertFit           = systUncertFit * doubleRatio->GetBinContent(j);
            
            totSystUncertUp         = TMath::Sqrt( yErrUp[j]*yErrUp[j] + systUncertFit*systUncertFit );
            totSystUncertDown       = TMath::Sqrt( yErrDown[j]*yErrDown[j] + systUncertFit*systUncertFit );
            
            yErrUpTot[i]            = totSystUncertUp;
            yErrDownTot[i]          = totSystUncertDown;
        }
    }
    
    TGraphAsymmErrors* systUncertTotal = new TGraphAsymmErrors(nPoints,xVal,yVal,xErrDown,xErrUp,yErrDownTot,yErrUpTot);
    systUncertTotal->SetName(Form("%s_total", systUncert->GetName()));
    return systUncertTotal;
}
