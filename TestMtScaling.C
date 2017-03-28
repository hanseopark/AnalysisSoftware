//  **********************************************************************************
//  ******     provided by Gamma Conversion Group, PWGGA,                        *****
//  ******     Lucas Altenkämper, lucas.altenkamper@cern.ch                      *****
//  **********************************************************************************
#include "TestMtScaling.h"

Bool_t kIsMC = kFALSE;

//  **********************************************************************************
//  ******     main function                                                     *****
//  **********************************************************************************
void TestMtScaling(     TString     fileNamePi0                 = "",
                        TString     fileNamePi0EtaBinning       = "",
                        TString     fileNameEta                 = "",
                        TString     fileNameCocktailInput       = "",
                        TString     fileNameCocktailProduction  = "",
                        TString     fCutSelection               = "",
                        TString     suffix                      = "",
                        TString     isMC                        = "",
                        TString     optionEnergy                = "",
                        TString     optionPeriod                = "",
                        Bool_t      optDalitz                   = kFALSE,
                        Int_t       mode                        = 0,
                        Bool_t      recalcScalingFactors        = kTRUE,
                        Double_t    rapidityCocktail            = 0.85
                     ) {
    
    //  ******************************************************************************
    //  ******     general style settings                                        *****
    //  ******************************************************************************
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    StyleSettingsThesis(suffix);
    SetPlotStyle();

    TString fDalitz                                                     = "";
    Bool_t kDalitz                                                      = optDalitz;
    if (kDalitz)    fDalitz                                             = "Dalitz";
    else            fDalitz                                             = "";
    
    //  ******************************************************************************
    //  ******     set global variables                                          *****
    //  ******************************************************************************
    outputDir                                                           = Form("%s/%s/%s/%s/TestMtScaling%s",fCutSelection.Data(),optionEnergy.Data(),optionPeriod.Data(),suffix.Data(),fDalitz.Data());
    gSystem->Exec("mkdir -p "+outputDir);

    TString date                                                        = ReturnDateString();
    TString prefix2                                                     = "";

    // plot labeling
    fDetectionProcess                                                   = ReturnFullTextReconstructionProcess(mode);
    collisionSystem                                                     = ReturnFullCollisionsSystem(optionEnergy);

    // cut strings
    TString fEventCutSelection                                          = "";
    TString fGammaCutSelection                                          = "";
    TString fClusterCutSelection                                        = "";
    TString fElectronCutSelection                                       = "";
    TString fMesonCutSelection                                          = "";

    if (mode == 9){
        if (kDalitz)    ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection,kTRUE);
        else            ReturnSeparatedCutNumber(fCutSelection, fGammaCutSelection, fElectronCutSelection,fMesonCutSelection);
        
        fEventCutSelection                                              = fGammaCutSelection(0,7);
        fGammaCutSelection                                              = fGammaCutSelection(7,fGammaCutSelection.Length()-7);
        cout << fEventCutSelection.Data() << "\t" << fGammaCutSelection.Data() << endl;
        
    } else {
        ReturnSeparatedCutNumberAdvanced(fCutSelection,fEventCutSelection, fGammaCutSelection, fClusterCutSelection, fElectronCutSelection, fMesonCutSelection, mode);
    }
    
    // Set flags for MC case
    if (isMC.CompareTo("kTRUE") == 0){
        prefix2                                                         = "MC";
        kIsMC                                                           = kTRUE;
    } else {
        prefix2                                                         = "data";
    }
    
    // flags for collisions sytem
    energyString                                                        = optionEnergy;
    
    if (collisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    Int_t kCollisionSystem                                              = 0; // 0 : pp, 1: PbPb, 2: pPb
    if (optionEnergy.CompareTo("PbPb_2.76TeV") == 0) kCollisionSystem   = 1;
    if (optionEnergy.CompareTo("pPb_5.023TeV") == 0) kCollisionSystem   = 2;
    
    TString centralityString                                            = GetCentralityString(fEventCutSelection);
    TString fTextCent                                                   = "";
    if (centralityString.CompareTo("pp")==0){
        fTextCent                                                       = "MinBias";
        centralityString                                                = "";
    } else {
        fTextCent                                                       = Form("%s central", centralityString.Data());
        if ( !centralityString.Contains("0-100%") )
            collisionSystem                                             = Form("%s", collisionSystem.Data());
    }
    if (optionPeriod.CompareTo("") != 0 ){
        collisionSystem                                                 = Form("%s, %s",collisionSystem.Data(),optionPeriod.Data());
    }
    
    TString rapidityRange                                               = "";
    Double_t deltaRapid                                                 = ReturnRapidityStringAndDouble(fMesonCutSelection, rapidityRange);
    Double_t rapidityRangeDouble                                        = rapidityRange.Atof();
    if (rapidityCocktail) rapidityRangeDouble                           = rapidityCocktail;     // only used to read in from cocktail production file
    
    TString trigger                                                     = fEventCutSelection(GetEventSelectSpecialTriggerCutPosition(),2);
    
    // event normalization scaling factor for cocktail quantities (spectra + params)
    eventNormScalingFactor                                              = ReturnCocktailNormalization(optionEnergy, fEventCutSelection);

    //  ******************************************************************************
    //  ******     read in from files                                            *****
    //  ******************************************************************************
    Bool_t hasFilePi0                                                   = kTRUE;
    Bool_t hasFilePi0EtaBinning                                         = kTRUE;
    Bool_t hasFileEta                                                   = kTRUE;
    Bool_t hasFileCocktailProduction                                    = kTRUE;
    Bool_t hasFileCocktailInput                                         = kTRUE;
    
    // pi0 ***************************************************************************
    TFile filePi0(fileNamePi0.Data());
    if (filePi0.IsZombie()) {
        cout << "No pi0 file found, skipping pi0 related actions!" << endl;
        hasFilePi0                                                      = kFALSE;
    }
    
    if (hasFilePi0) {
        pi0InvYield                                                     = (TH1D*)filePi0.Get("CorrectedYieldTrueEff");
        pi0InvYieldWOResFeedDown                                        = (TH1D*)filePi0.Get("FeedDownCorrectedYieldTrueEff");
        
        if (pi0InvYield)                pi0InvYield->SetName(               Form("%d_invYield",                 particlePDGCode[0]));
        if (pi0InvYieldWOResFeedDown)   pi0InvYieldWOResFeedDown->SetName(  Form("%d_invYield_resFeedDownCorr", particlePDGCode[0]));
    }
    
    // pi0 in eta binning ************************************************************
    TFile filePi0EtaBinning(fileNamePi0EtaBinning.Data());
    if (filePi0EtaBinning.IsZombie()) {
        cout << "No pi0 eta binning file found, skipping pi0 eta binning related actions!" << endl;
        hasFilePi0EtaBinning                                            = kFALSE;
    }
    
    if (hasFilePi0EtaBinning) {
        pi0EtaBinningInvYield                                           = (TH1D*)filePi0EtaBinning.Get("CorrectedYieldTrueEff");
        pi0EtaBinningInvYieldWOResFeedDown                              = (TH1D*)filePi0EtaBinning.Get("FeedDownCorrectedYieldTrueEff");
        
        if (pi0EtaBinningInvYield)                pi0EtaBinningInvYield->SetName(               Form("%d_invYield_etaBinning",                 particlePDGCode[0]));
        if (pi0EtaBinningInvYieldWOResFeedDown)   pi0EtaBinningInvYieldWOResFeedDown->SetName(  Form("%d_invYield_etaBinning_resFeedDownCorr", particlePDGCode[0]));
    }
    
    // eta ***************************************************************************
    TFile fileEta(fileNameEta.Data());
    if (fileEta.IsZombie()) {
        cout << "No eta file found, skipping eta related actions!" << endl;
        hasFileEta                                                      = kFALSE;
    }
    
    if (hasFileEta) {
        etaInvYield                                                     = (TH1D*)fileEta.Get("CorrectedYieldTrueEff");
        
        if (etaInvYield) etaInvYield->SetName(Form("%d_invYield", particlePDGCode[1]));
    }
    
    // cocktail input ****************************************************************
    TFile fileCocktailInput(fileNameCocktailInput.Data());
    if (fileCocktailInput.IsZombie()) {
        cout << "No cocktail input file found, skipping comparisons to measured particle ratios!" << endl;
        hasFileCocktailInput                                            = kFALSE;
    }

    // list name in cocktail input file
    TString fileCocktailInputListName                                   = "";
    if (kCollisionSystem == 0) fileCocktailInputListName                = Form("pp_%s", optionEnergy.Data());
    else if (kCollisionSystem == 1) {
        cout << "Macro not yet final for kCollisionSystem = " << kCollisionSystem << ", loading MB list." << endl;
        fileCocktailInputListName                                       = "2.76TeV_MB";
    } else if (kCollisionSystem == 2) {
        cout << "Macro not yet final for kCollisionSystem = " << kCollisionSystem << ", loading MB list." << endl;
        fileCocktailInputListName                                       = "5TeV_MB";
    } else {
        cout << "Collision system not recognized, skipping read in from cocktail input file!" << endl;
        hasFileCocktailInput                                            = kFALSE;
    }

    TList* fileCocktailInputList                                        = NULL;
    if (hasFileCocktailInput) {
        fileCocktailInputList                                           = (TList*)fileCocktailInput.Get(fileCocktailInputListName.Data());
        if (!fileCocktailInputList) {
            cout << "List " << fileCocktailInputListName.Data() << " not found in " << fileNameCocktailInput.Data() << ", skipping" << endl;
            hasFileCocktailInput                                        = kFALSE;
        }
    }
    
    if (hasFileCocktailInput) {
        
        TObject*    tempSpectrum                                        = NULL;
        TObject*    tempRatio                                           = NULL;
        Int_t       tempParticleRatioIter                               = 0;
        for (Int_t i = 0; i < nParticles; i++) {
            // spectra
            TString nameSpectrum                                        = "";
            if (i == 0 || i == 1) {     // pi0 or eta
                nameSpectrum                                            = Form("%sPCMStat",     particleNameCocktailInput[i].Data());
            } else if (i == 3) {        // omega
                nameSpectrum                                            = Form("%sPHOSStat",    particleNameCocktailInput[i].Data());
            } else if (i == 6) {        // phi
                nameSpectrum                                            = Form("%sCombStat",    particleNameCocktailInput[i].Data());
            } else {
                nameSpectrum                                            = Form("%sStat",        particleNameCocktailInput[i].Data());
            }
            tempSpectrum                                                = (TObject*)fileCocktailInputList->FindObject(nameSpectrum.Data());
            
            if (tempSpectrum && i != 0 && i != 1) {
                particleYield[i]                                        = (TH1D*)TransformToTH1D(tempSpectrum, Form("%d_yield", particlePDGCode[i]));
                particleYield[i]->Sumw2();
                particleYield[i]->Scale(eventNormScalingFactor);
            } else
                particleYield[i]                                        = NULL;
            
            for (Int_t j = 0; j < nParticles; j++) {
                // ratios
                TString nameRatio                                       = "";
                if (i == 0 && j == 1) {             // pi0 or eta
                    nameRatio                                           = Form("%sTo%sPCMStat",     particleNameCocktailInput[j].Data(), particleNameCocktailInput[i].Data());
                } else if (i == 0 && j == 3) {      // omega
                    nameRatio                                           = Form("%sTo%sPHOSStat",    particleNameCocktailInput[j].Data(), particleNameCocktailInput[i].Data());
                } else if (i == 0 && j == 6) {      // phi
                    nameRatio                                           = Form("%sTo%sCombStat",    particleNameCocktailInput[j].Data(), particleNameCocktailInput[i].Data());
                } else {
                    nameRatio                                           = Form("%sTo%sStat",        particleNameCocktailInput[j].Data(), particleNameCocktailInput[i].Data());
                }
                tempRatio                                               = (TObject*)fileCocktailInputList->FindObject(nameRatio.Data());
                tempParticleRatioIter                                   = GetParticleRatioIter(particlePDGCode[j], particlePDGCode[i]);
                
                if (tempRatio && !(i == 0 && j == 1))
                    particleRatio[tempParticleRatioIter]                = (TH1D*)TransformToTH1D(tempRatio, Form("%d-%d_ratio", particlePDGCode[j], particlePDGCode[i]));
                else
                    particleRatio[tempParticleRatioIter]                = NULL;
            }
        }
    }
    
    // cocktail train output *********************************************************
    TFile fileCocktailProduction(fileNameCocktailProduction.Data());
    if (fileCocktailProduction.IsZombie()) {
        cout << "No cocktail production file found, skipping cocktail related actions!" << endl;
        hasFileCocktailProduction                                       = kFALSE;
    }
    
    TDirectoryFile* fileCocktailProductionDir                           = NULL;
    if (hasFileCocktailProduction) {
        fileCocktailProductionDir                                       = (TDirectoryFile*)fileCocktailProduction.Get("HadronicCocktailMC");
        if (!fileCocktailProductionDir) {
            cout << "ERROR: HadronicCocktailMC not found!" << endl;
            hasFileCocktailProduction                                   = kFALSE;
        }
    }
    
    TList*  fileCocktailProductionDirList                               = NULL;
    TString fileCocktailProductionDirListName                           = Form("HadronicCocktailMC_pi0_%.2f", rapidityRangeDouble);
    if (hasFileCocktailProduction) {
        fileCocktailProductionDirList                                   = (TList*)fileCocktailProductionDir->Get(fileCocktailProductionDirListName);
        if (!fileCocktailProductionDirList) {
            cout << "ERROR: " << fileCocktailProductionDirListName.Data() << " not found!" << endl;
            hasFileCocktailProduction                                   = kFALSE;
        }
    }
    
    if (hasFileCocktailProduction) {
        TTree* cocktailSettingsTree                                     = (TTree*)fileCocktailProductionDirList->FindObject("cocktailSettings");
        TList* cocktailSettingsList                                     = NULL;
        if (cocktailSettingsTree) {
            cocktailSettingsList                                        = (TList*)cocktailSettingsTree->GetUserInfo();
            if (cocktailSettingsList) {
                
                // get mt scaling factors from file
                TH1F* histMtScalingFactors                              = (TH1F*)cocktailSettingsList->FindObject("histoMtScaleFactor");
                TString tempBinLabel                                    = "";
                if (histMtScalingFactors) {
                    for (Int_t i=1; i<histMtScalingFactors->GetNbinsX()+1; i++) {
                        tempBinLabel                                    = (TString)histMtScalingFactors->GetXaxis()->GetBinLabel(i);
                        for (Int_t j=0; j<nParticles; j++) {
                            if (tempBinLabel.CompareTo(Form("%d",particlePDGCode[j])) == 0)
                                particleMtScalingFactor[j]             = histMtScalingFactors->GetBinContent(i);
                        }
                    }
                }
                
                // get parametrizations from file
                TF1* paramTemp                                          = NULL;
                TF1* paramMtTemp                                        = NULL;
                for (Int_t i=0; i<nParticles; i++) {
                    paramTemp                                           = (TF1*)cocktailSettingsList->FindObject(Form("%d_pt",          particlePDGCode[i]));
                    //paramMtTemp                                         = (TF1*)cocktailSettingsList->FindObject(Form("%d_pt_mtScaled", particlePDGCode[i]));
                    
                    if (paramTemp) {
                        if (eventNormScalingFactor == 1.)
                            yieldParametrizations[i]                    = (TF1*)paramTemp->Clone(Form("%d_yieldParam", particlePDGCode[i]));
                        else
                            yieldParametrizations[i]                    = (TF1*)ScaleTF1(paramTemp, eventNormScalingFactor, Form("%d_yieldParam", particlePDGCode[i]));
                    } else
                        yieldParametrizations[i]                        = NULL;
                    
                    // no need to load... will be (re-)calculated anyway
                    //if (paramMtTemp)
                    //    yieldParametrizationsMtScaled[i]                = (TF1*)paramMtTemp->Clone(Form("%d_yieldParam_mtScaled", particlePDGCode[i]));
                    //else
                    //    yieldParametrizationsMtScaled[i]                = NULL;
                }
            }
        }
    }
    
    // fill mt scaling factors from cocktail train input in histo
    mtScalingFactorCocktail                                             = FillMtScalingFactorsHisto("mtScalingFactors_Cocktail", particleMtScalingFactor, particleMtScalingFactorErr);
    
    // skip rest of macro if no files were provided
    if (!hasFilePi0 && !hasFilePi0EtaBinning && !hasFileEta && !hasFileCocktailInput && !hasFileCocktailProduction) {
        cout << "No files provided, returning!" << endl;
        return;
    }
    
    //  ******************************************************************************
    //  ******     transform invariant yields to yields                          *****
    //  ******************************************************************************
    if (pi0InvYield)                        pi0Yield                        = (TH1D*)TransformInvYield(pi0InvYield, Form("%d_yield", particlePDGCode[0]));
    if (pi0InvYieldWOResFeedDown)           pi0YieldWOResFeedDown           = (TH1D*)TransformInvYield(pi0InvYieldWOResFeedDown, Form("%d_yield_resFeedDownCorr", particlePDGCode[0]));
    if (pi0EtaBinningInvYield)              pi0EtaBinningYield              = (TH1D*)TransformInvYield(pi0EtaBinningInvYield, Form("%d_yield_etaBinning", particlePDGCode[0]));
    if (pi0EtaBinningInvYieldWOResFeedDown) pi0EtaBinningYieldWOResFeedDown = (TH1D*)TransformInvYield(pi0EtaBinningInvYieldWOResFeedDown, Form("%d_yield_etaBinning_resFeedDownCorr", particlePDGCode[0]));
    if (etaInvYield)                        etaYield                        = (TH1D*)TransformInvYield(etaInvYield, Form("%d_yield", particlePDGCode[1]));
    
    //  ******************************************************************************
    //  ******     calculate eta to pi0 ratio                                    *****
    //  ******************************************************************************
    if (etaYield && pi0EtaBinningYield) {
        etaToPi0Ratio                                                       = (TH1D*)etaYield->Clone(Form("%d-%d_ratio", particlePDGCode[1], particlePDGCode[0]));
        etaToPi0Ratio->Sumw2();
        etaToPi0Ratio->Divide(pi0EtaBinningYield);
    }
    
    if (etaYield && pi0EtaBinningYieldWOResFeedDown) {
        etaToPi0RatioWOResFeedDown                                          = (TH1D*)etaYield->Clone(Form("%d-%d_ratio_resFeedDownCorr", particlePDGCode[1], particlePDGCode[0]));
        etaToPi0RatioWOResFeedDown->Sumw2();
        etaToPi0RatioWOResFeedDown->Divide(pi0EtaBinningYieldWOResFeedDown);
    }
    
    //  ******************************************************************************
    //  ******     fit resonance feed down corrected pi0 yield                   *****
    //  ******************************************************************************
    if (pi0YieldWOResFeedDown && yieldParametrizations[0]) {
        yieldParametrizationPi0WOResFeedDown                                = (TF1*)yieldParametrizations[0]->Clone(Form("%d_yieldParam_resFeedDownCorr", particlePDGCode[0]));
        TFitResultPtr   fitResultPtr                                        = pi0YieldWOResFeedDown->Fit(yieldParametrizationPi0WOResFeedDown,"QMNRE+","",pi0YieldWOResFeedDown->GetXaxis()->GetXmin(),pi0YieldWOResFeedDown->GetXaxis()->GetXmax());
        Int_t           fitStatus                                           = fitResultPtr;
        
        if (fitStatus != 0 && fitStatus < 1000)
            yieldParametrizationPi0WOResFeedDown                            = NULL;
        
        if (yieldParametrizationPi0WOResFeedDown)
            yieldParametrizationPi0WOResFeedDownRatioToData                 = CalculateRatioToFit(pi0YieldWOResFeedDown, yieldParametrizationPi0WOResFeedDown, Form("%s_RatioToFit", pi0YieldWOResFeedDown->GetName()));
    }
    
    if (pi0YieldWOResFeedDown && yieldParametrizationPi0WOResFeedDown && yieldParametrizationPi0WOResFeedDownRatioToData) {
        PlotYield(pi0YieldWOResFeedDown, yieldParametrizationPi0WOResFeedDown, yieldParametrizationPi0WOResFeedDownRatioToData, NULL, NULL, NULL, NULL, NULL, NULL, 0, "Pi0ResFeedDownCorrect_Param", suffix);
    }

    //  ******************************************************************************
    //  ******     calculate mt scaling factors from particle ratios             *****
    //  ******************************************************************************
    for (Int_t i = 0; i < nParticles*nParticles; i++) {
        constantFitParticleRatio[i]                                         = NULL;
    }

    if (recalcScalingFactors) {
        Double_t    constantXToPi0                                          = 1.;
        Double_t    constantXToPi0Err                                       = 0.;
        Int_t       tempIter                                                = 0;
        
        // take mt scaling factor for eta from measured eta/pi0 ratio if available
        if (etaToPi0Ratio) {
            cout << "using " << particleName[1] << " to " << particleName[0] << " ratio for recalculation of " << particleName[1].Data() << " mt scaling factor" << endl;
            tempIter                                                        = GetParticleRatioIter(particlePDGCode[1], particlePDGCode[0]);
            constantFitParticleRatio[tempIter]                              = FitPlateau(etaToPi0Ratio);
            constantFitParticleRatio[tempIter]->SetName("221-111_ratio_constFit");
            constantXToPi0                                                  = constantFitParticleRatio[tempIter]->GetParameter(0);
            constantXToPi0Err                                               = constantFitParticleRatio[tempIter]->GetParError(0);
            particleMtScalingFactor[1]                                      = constantXToPi0;
            particleMtScalingFactorErr[1]                                   = constantXToPi0Err;
        }
        
        // calculate scaling factors from particle ratios
        for (Int_t i = 1; i < nParticles; i++) {
                
            // X to pi0 ratio
            tempIter                                                        = GetParticleRatioIter(particlePDGCode[i], particlePDGCode[0]);
            if (particleRatio[tempIter]) {
                cout << "using " << particleName[i] << " to " << particleName[0] << " ratio for recalculation of " << particleName[i] << " mt scaling factor" << endl;
                constantFitParticleRatio[tempIter]                          = FitPlateau(particleRatio[tempIter]);
                constantFitParticleRatio[tempIter]->SetName(Form("%d-%d_ratio_constFit",particlePDGCode[i],particlePDGCode[0]));
                constantXToPi0                                              = constantFitParticleRatio[tempIter]->GetParameter(0);
                constantXToPi0Err                                           = constantFitParticleRatio[tempIter]->GetParError(0);
                particleMtScalingFactor[i]                                  = constantXToPi0;
                particleMtScalingFactorErr[i]                               = constantXToPi0Err;
            }
            
            
            // X to pi+- ratio
            tempIter                                                        = GetParticleRatioIter(particlePDGCode[i], particlePDGCode[4]);
            if (particleRatio[tempIter]) {
                cout << "using " << particleName[i] << " to " << particleName[4] << " ratio for recalculation of " << particleName[i] << " mt scaling factor, be careful with the overall normalization, check input!" << endl;
                constantFitParticleRatio[tempIter]                          = FitPlateau(particleRatio[tempIter]);
                constantFitParticleRatio[tempIter]->SetName(Form("%d-%d_ratio_constFit",particlePDGCode[i],particlePDGCode[4]));
                constantXToPi0                                              = constantFitParticleRatio[tempIter]->GetParameter(0);
                constantXToPi0Err                                           = constantFitParticleRatio[tempIter]->GetParError(0);
                particleMtScalingFactor[i]                                  = constantXToPi0;
                particleMtScalingFactorErr[i]                               = constantXToPi0Err;
            }
        }
    }
    
    // fill calculated mt scaling factors in histo
    mtScalingFactorRecalc                                                   = FillMtScalingFactorsHisto("mtScalingFactors_Recalc", particleMtScalingFactor, particleMtScalingFactorErr);

    //  ******************************************************************************
    //  ******     calculate mt scaling factors (res. feed down corr. pi0)       *****
    //  ******************************************************************************
    if (etaToPi0RatioWOResFeedDown) {

        Double_t constantEtaToPi0WOResFeedDown                              = 0.;
        Double_t constantEtaToPi0WOResFeedDownErr                           = 0.;

        constantFitEtaToPi0WOResFeedDownRatio                               = FitPlateau(etaToPi0RatioWOResFeedDown);
        constantFitEtaToPi0WOResFeedDownRatio->SetName("221-111_ratio_resFeedDownCorr_constFit");
        constantEtaToPi0WOResFeedDown                                       = constantFitEtaToPi0WOResFeedDownRatio->GetParameter(0);
        constantEtaToPi0WOResFeedDownErr                                    = constantFitEtaToPi0WOResFeedDownRatio->GetParError(0);

        particleMtScalingFactorResFeedDownCorr[1]                           = constantEtaToPi0WOResFeedDown;
        particleMtScalingFactorErrResFeedDownCorr[1]                        = constantEtaToPi0WOResFeedDownErr;

        for (Int_t i = 2; i < nParticles; i++) {
            particleMtScalingFactorResFeedDownCorr[i]                       = particleMtScalingFactor[i] * constantEtaToPi0WOResFeedDown / particleMtScalingFactor[1];
            particleMtScalingFactorErrResFeedDownCorr[i]                    = TMath::Sqrt( TMath::Power(particleMtScalingFactorErr[i]*constantEtaToPi0WOResFeedDown/particleMtScalingFactor[1],2)
                                                                                          + TMath::Power(particleMtScalingFactor[i]*constantEtaToPi0WOResFeedDownErr/particleMtScalingFactor[1],2)
                                                                                          + TMath::Power(particleMtScalingFactor[i]*constantEtaToPi0WOResFeedDown*particleMtScalingFactorErr[1]/particleMtScalingFactor[1]/particleMtScalingFactor[1],2));
        }

        // fill calculated mt scaling factors in histo
        mtScalingFactorResFeedDownCorr                                      = FillMtScalingFactorsHisto("mtScalingFactors_resFeedDownCorr", particleMtScalingFactorResFeedDownCorr, particleMtScalingFactorErrResFeedDownCorr);
    }
    
    //  ******************************************************************************
    //  ******     calculate scaling factors for scaling from phi                *****
    //  ******************************************************************************
    Double_t    constantPhiToPi0                                            = 0.11;
    Double_t    constantPhiToPi0Err                                         = 0.;
    TString     particleRatioName                                           = "";
    for (Int_t i = 0; i < nParticles*nParticles; i++) {
        if (!particleRatio[i]) continue;
        
        particleRatioName                                                   = particleRatio[i]->GetName();
        if (particleRatioName.CompareTo("333-111_ratio") == 0) {
            constantFitPhiToPi0Ratio                                        = FitPlateau(particleRatio[i]);
            constantFitPhiToPi0Ratio->SetName("333-111_ratio_constFit");
            constantPhiToPi0                                                = constantFitPhiToPi0Ratio->GetParameter(0);
            constantPhiToPi0Err                                             = constantFitPhiToPi0Ratio->GetParError(0);
            break;
        }
    }
    
    for (Int_t i = 0; i < nParticles; i++) {
        particleMtScalingFactorPhi[i]                                       = particleMtScalingFactor[i] / constantPhiToPi0;
        particleMtScalingFactorErrPhi[i]                                    = TMath::Sqrt( TMath::Power( particleMtScalingFactorErr[i] / constantPhiToPi0, 2 ) + TMath::Power( (constantPhiToPi0Err * particleMtScalingFactor[i])/(constantPhiToPi0 * constantPhiToPi0), 2 ) );
    }
    
    // fill calculated mt scaling factors in histo
    mtScalingFactorPhi                                                      = FillMtScalingFactorsHisto("mtScalingFactors_phi", particleMtScalingFactorPhi, particleMtScalingFactorErrPhi);
    
    //  ******************************************************************************
    //  ******     calculate mt scaled parametrizations                          *****
    //  ******************************************************************************
    
    // scaling from resonance feed down corrected pi0
    if (yieldParametrizationPi0WOResFeedDown) {
        for (Int_t i = 0; i < nParticles; i++) {
            if (i > 0)  {
                yieldParametrizationsMtScaledWOResFeedDown[i]               = MtScaledParam(yieldParametrizationPi0WOResFeedDown, particlePDGCode[i], particlePDGCode[0], particleMtScalingFactor[i], kFALSE, kTRUE);
                yieldParametrizationsMtScaledWOResFeedDown[i]->SetName(Form("%d_yieldParam_mtScaled_resFeedDownCorr", particlePDGCode[i]));
            } else
                yieldParametrizationsMtScaledWOResFeedDown[i]               = NULL;
        }
    } else {
        for (Int_t i = 0; i < nParticles; i++)
            yieldParametrizationsMtScaledWOResFeedDown[i]                   = NULL;
    }

    // scaling from pi0
    if (yieldParametrizations[0]) {
        for (Int_t i = 0; i < nParticles; i++) {
            if (i > 0 /*&& !yieldParametrizationsMtScaled[i]*/) {   // recalc all mt scaled params (scaling factors might have been updated!)
                yieldParametrizationsMtScaled[i]                            = MtScaledParam(yieldParametrizations[0], particlePDGCode[i], particlePDGCode[0], particleMtScalingFactor[i], kFALSE, kTRUE);
                yieldParametrizationsMtScaled[i]->SetName(Form("%d_yieldParam_mtScaled", particlePDGCode[i]));
            } else
                yieldParametrizationsMtScaled[i]                            = NULL;
        }
    } else {
        for (Int_t i = 0; i < nParticles; i++)
            yieldParametrizationsMtScaled[i]                                = NULL;
    }
    
    // scaling from phi
    if (yieldParametrizations[6]) {
        for (Int_t i = 0; i < nParticles; i++) {
            if (i != 6) {
                yieldParametrizationsMtScaledPhi[i]                         = MtScaledParam(yieldParametrizations[6], particlePDGCode[i], particlePDGCode[6], particleMtScalingFactorPhi[i], kFALSE, kTRUE);
                yieldParametrizationsMtScaledPhi[i]->SetName(Form("%d_yieldParam_mtScaled_fromPhi", particlePDGCode[i]));
            } else
                yieldParametrizationsMtScaledPhi[i]                         = NULL;
        }
    } else {
        for (Int_t i = 0; i < nParticles; i++)
            yieldParametrizationsMtScaledPhi[i]                             = NULL;
    }
    
    //  ******************************************************************************
    //  ******     calculate ratio to param                                      *****
    //  ******************************************************************************
    if (pi0Yield) {
        if (yieldParametrizations[0])
            yieldParametrizationsRatioToData[0]                             = CalculateRatioToFit(pi0Yield, yieldParametrizations[0], Form("%s_RatioToFit", yieldParametrizations[0]->GetName()));
        else
            yieldParametrizationsRatioToData[0]                             = NULL;
        
        if (yieldParametrizationsMtScaled[0])
            yieldParametrizationsMtScaledRatioToData[0]                     = CalculateRatioToFit(pi0Yield, yieldParametrizationsMtScaled[0], Form("%s_RatioToFit", yieldParametrizationsMtScaled[0]->GetName()));
        else
            yieldParametrizationsMtScaledRatioToData[0]                     = NULL;
        
        if (yieldParametrizationsMtScaledPhi[0])
            yieldParametrizationsMtScaledPhiRatioToData[0]                  = CalculateRatioToFit(pi0Yield, yieldParametrizationsMtScaledPhi[0], Form("%s_RatioToFit", yieldParametrizationsMtScaledPhi[0]->GetName()));
        else
            yieldParametrizationsMtScaledPhiRatioToData[0]                  = NULL;
    }
    
    if (etaYield) {
        if (yieldParametrizations[1])
            yieldParametrizationsRatioToData[1]                             = CalculateRatioToFit(etaYield, yieldParametrizations[1], Form("%s_RatioToFit", yieldParametrizations[1]->GetName()));
        else
            yieldParametrizationsRatioToData[1]                             = NULL;
        
        if (yieldParametrizationsMtScaled[1])
            yieldParametrizationsMtScaledRatioToData[1]                     = CalculateRatioToFit(etaYield, yieldParametrizationsMtScaled[1], Form("%s_RatioToFit", yieldParametrizationsMtScaled[1]->GetName()));
        else
            yieldParametrizationsMtScaledRatioToData[1]                     = NULL;
        
        if (yieldParametrizationsMtScaledWOResFeedDown[1])
            yieldParametrizationsMtScaledWOResFeedDownRatioToData[1]        = CalculateRatioToFit(etaYield, yieldParametrizationsMtScaledWOResFeedDown[1], Form("%s_RatioToFit", yieldParametrizationsMtScaledWOResFeedDown[1]->GetName()));
        else
            yieldParametrizationsMtScaledWOResFeedDownRatioToData[1]        = NULL;
        
        if (yieldParametrizationsMtScaledPhi[1])
            yieldParametrizationsMtScaledPhiRatioToData[1]                  = CalculateRatioToFit(etaYield, yieldParametrizationsMtScaledPhi[1], Form("%s_RatioToFit", yieldParametrizationsMtScaledPhi[1]->GetName()));
        else
            yieldParametrizationsMtScaledPhiRatioToData[1]                  = NULL;
    }
    
    for (Int_t i = 2; i < nParticles; i++) {
        if (particleYield[i] && yieldParametrizations[i])
            yieldParametrizationsRatioToData[i]                             = CalculateRatioToFit(particleYield[i], yieldParametrizations[i], Form("%s_RatioToFit", yieldParametrizations[i]->GetName()));
        else
            yieldParametrizationsRatioToData[i]                             = NULL;

        if (particleYield[i] && yieldParametrizationsMtScaled[i])
            yieldParametrizationsMtScaledRatioToData[i]                     = CalculateRatioToFit(particleYield[i], yieldParametrizationsMtScaled[i], Form("%s_RatioToFit", yieldParametrizationsMtScaled[i]->GetName()));
        else
            yieldParametrizationsMtScaledRatioToData[i]                     = NULL;

        if (particleYield[i] && yieldParametrizationsMtScaledWOResFeedDown[i])
            yieldParametrizationsMtScaledWOResFeedDownRatioToData[i]        = CalculateRatioToFit(particleYield[i], yieldParametrizationsMtScaledWOResFeedDown[i], Form("%s_RatioToFit", yieldParametrizationsMtScaledWOResFeedDown[i]->GetName()));
        else
            yieldParametrizationsMtScaledWOResFeedDownRatioToData[i]        = NULL;
        
        if (particleYield[i] && yieldParametrizationsMtScaledPhi[i])
            yieldParametrizationsMtScaledPhiRatioToData[i]                  = CalculateRatioToFit(particleYield[i], yieldParametrizationsMtScaledPhi[i], Form("%s_RatioToFit", yieldParametrizationsMtScaledPhi[i]->GetName()));
        else
            yieldParametrizationsMtScaledPhiRatioToData[i]                  = NULL;
    }
    
    //  ******************************************************************************
    //  ******     comparisons...                                                *****
    //  ******************************************************************************
    if (pi0Yield) {
        PlotYield(pi0Yield, yieldParametrizations[0], yieldParametrizationsRatioToData[0], NULL, NULL, NULL, NULL, yieldParametrizationsMtScaledPhi[0], yieldParametrizationsMtScaledPhiRatioToData[0], 0, "Pi0Yield_MtScalingFromPhi", suffix);
    }

    if (etaYield) {
        PlotYield(etaYield, yieldParametrizations[1], yieldParametrizationsRatioToData[1], yieldParametrizationsMtScaled[1], yieldParametrizationsMtScaledRatioToData[1], yieldParametrizationsMtScaledWOResFeedDown[1], yieldParametrizationsMtScaledWOResFeedDownRatioToData[1], yieldParametrizationsMtScaledPhi[1], yieldParametrizationsMtScaledPhiRatioToData[1], 1, "EtaYield_Comparisons", suffix);
    }

    for (Int_t i = 0; i < nParticles; i++) {
        if (particleYield[i]) {
            PlotYield(particleYield[i], yieldParametrizations[i], yieldParametrizationsRatioToData[i], yieldParametrizationsMtScaled[i], yieldParametrizationsMtScaledRatioToData[i], yieldParametrizationsMtScaledWOResFeedDown[i], yieldParametrizationsMtScaledWOResFeedDownRatioToData[i], yieldParametrizationsMtScaledPhi[i], yieldParametrizationsMtScaledPhiRatioToData[i], i, Form("%sYield_Comparisons",particleName[i].Data()), suffix);
        }
    }

    //  ******************************************************************************
    //  ******     write to file                                                 *****
    //  ******************************************************************************
    TString nameOutputFile                                                  = Form("%s/%s/TestMtScaling%s_%s.root", fCutSelection.Data(), optionEnergy.Data(), optionPeriod.Data(), fCutSelection.Data());
    TFile*  outputFile                                                      = new TFile(nameOutputFile.Data(),"UPDATE");
    cout << "INFO: writing into: " << nameOutputFile << endl;

    if (pi0InvYield)                        pi0InvYield->Write(                             pi0InvYield->GetName(),                             TObject::kOverwrite);
    if (pi0Yield)                           pi0Yield->Write(                                pi0Yield->GetName(),                                TObject::kOverwrite);
    if (pi0InvYieldWOResFeedDown)           pi0InvYieldWOResFeedDown->Write(                pi0InvYieldWOResFeedDown->GetName(),                TObject::kOverwrite);
    if (pi0YieldWOResFeedDown)              pi0YieldWOResFeedDown->Write(                   pi0YieldWOResFeedDown->GetName(),                   TObject::kOverwrite);
    
    if (pi0EtaBinningInvYield)              pi0EtaBinningInvYield->Write(                   pi0EtaBinningInvYield->GetName(),                   TObject::kOverwrite);
    if (pi0EtaBinningYield)                 pi0EtaBinningYield->Write(                      pi0EtaBinningYield->GetName(),                      TObject::kOverwrite);
    if (pi0EtaBinningInvYieldWOResFeedDown) pi0EtaBinningInvYieldWOResFeedDown->Write(      pi0EtaBinningInvYieldWOResFeedDown->GetName(),      TObject::kOverwrite);
    if (pi0EtaBinningYieldWOResFeedDown)    pi0EtaBinningYieldWOResFeedDown->Write(         pi0EtaBinningYieldWOResFeedDown->GetName(),         TObject::kOverwrite);
    
    if (etaInvYield)                        etaInvYield->Write(                             etaInvYield->GetName(),                             TObject::kOverwrite);
    if (etaYield)                           etaYield->Write(                                etaYield->GetName(),                                TObject::kOverwrite);
    
    if (etaToPi0Ratio)                      etaToPi0Ratio->Write(                           etaToPi0Ratio->GetName(),                           TObject::kOverwrite);
    if (etaToPi0RatioWOResFeedDown)         etaToPi0RatioWOResFeedDown->Write(              etaToPi0RatioWOResFeedDown->GetName(),              TObject::kOverwrite);

    if (constantFitPhiToPi0Ratio)           constantFitPhiToPi0Ratio->Write(                constantFitPhiToPi0Ratio->GetName(),                TObject::kOverwrite);
    if (constantFitEtaToPi0WOResFeedDownRatio)
                                            constantFitEtaToPi0WOResFeedDownRatio->Write(   constantFitEtaToPi0WOResFeedDownRatio->GetName(),   TObject::kOverwrite);
    
    for (Int_t i = 0; i < nParticles; i++) {
        if (particleYield[i])               particleYield[i]->Write(                        particleYield[i]->GetName(),                        TObject::kOverwrite);
    }
    
    for (Int_t i = 0; i < nParticles*nParticles; i++) {
        if (particleRatio[i])               particleRatio[i]->Write(                        particleRatio[i]->GetName(),                        TObject::kOverwrite);
        if (constantFitParticleRatio[i])    constantFitParticleRatio[i]->Write(             constantFitParticleRatio[i]->GetName(),             TObject::kOverwrite);
    }
    
    if (mtScalingFactorCocktail)            mtScalingFactorCocktail->Write(                 mtScalingFactorCocktail->GetName(),                 TObject::kOverwrite);
    if (mtScalingFactorRecalc)              mtScalingFactorRecalc->Write(                   mtScalingFactorRecalc->GetName(),                   TObject::kOverwrite);
    if (mtScalingFactorPhi)                 mtScalingFactorPhi->Write(                      mtScalingFactorPhi->GetName(),                      TObject::kOverwrite);
    if (mtScalingFactorResFeedDownCorr)     mtScalingFactorResFeedDownCorr->Write(          mtScalingFactorResFeedDownCorr->GetName(),          TObject::kOverwrite);

    if (yieldParametrizationPi0WOResFeedDown)
        yieldParametrizationPi0WOResFeedDown->Write(            yieldParametrizationPi0WOResFeedDown->GetName(),                        TObject::kOverwrite);
    if (yieldParametrizationPi0WOResFeedDownRatioToData)
        yieldParametrizationPi0WOResFeedDownRatioToData->Write( yieldParametrizationPi0WOResFeedDownRatioToData->GetName(),             TObject::kOverwrite);
    
    for (Int_t i = 0; i < nParticles; i++) {
        if (yieldParametrizations[i])
            yieldParametrizations[i]->Write(                                yieldParametrizations[i]->GetName(),                                    TObject::kOverwrite);
        if (yieldParametrizationsRatioToData[i])
            yieldParametrizationsRatioToData[i]->Write(                     yieldParametrizationsRatioToData[i]->GetName(),                         TObject::kOverwrite);
        
        if (yieldParametrizationsMtScaled[i])
            yieldParametrizationsMtScaled[i]->Write(                        yieldParametrizationsMtScaled[i]->GetName(),                            TObject::kOverwrite);
        if (yieldParametrizationsMtScaledRatioToData[i])
            yieldParametrizationsMtScaledRatioToData[i]->Write(             yieldParametrizationsMtScaledRatioToData[i]->GetName(),                 TObject::kOverwrite);

        if (yieldParametrizationsMtScaledWOResFeedDown[i])
            yieldParametrizationsMtScaledWOResFeedDown[i]->Write(           yieldParametrizationsMtScaledWOResFeedDown[i]->GetName(),               TObject::kOverwrite);
        if (yieldParametrizationsMtScaledWOResFeedDownRatioToData[i])
            yieldParametrizationsMtScaledWOResFeedDownRatioToData[i]->Write(yieldParametrizationsMtScaledWOResFeedDownRatioToData[i]->GetName(),    TObject::kOverwrite);

        if (yieldParametrizationsMtScaledPhi[i])
            yieldParametrizationsMtScaledPhi[i]->Write(                     yieldParametrizationsMtScaledPhi[i]->GetName(),                         TObject::kOverwrite);
        if (yieldParametrizationsMtScaledPhiRatioToData[i])
            yieldParametrizationsMtScaledPhiRatioToData[i]->Write(          yieldParametrizationsMtScaledPhiRatioToData[i]->GetName(),              TObject::kOverwrite);
    }
    
    outputFile->Write();
    outputFile->Close();
}
