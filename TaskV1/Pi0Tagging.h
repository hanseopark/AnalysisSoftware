#ifndef GAMMACONV_Pi0Tagging
#define GAMMACONV_Pi0Tagging

    #include "../CommonHeaders/ExtractSignalBinning.h"

    //***********************************************************************************************
    //**************************  variables for global configuration ********************************
    //***********************************************************************************************
    TString nameRec;
    TString nameOutputLabel;

    TString triggerCutNumber                                    = "";
    TString subTriggerCutNumber                                 = "";
    TString fEventCutSelection                                  = "";
    TString fGammaCutSelection                                  = "";
    TString fClusterCutSelection                                = "";
    TString fElectronCutSelection                               = "";
    TString fMesonCutSelection                                  = "";

    //***********************************************************************************************
    //********************************  Inputs for pi0 spectrum  ************************************
    //***********************************************************************************************
    TFile *filePi0                                              = NULL;
    TFile *filePi0Effi                                          = NULL;
    TH1D* histoInputTaggedPi0                                   = NULL;
    TH1D* histoTaggedPi0_SecAndEffiCorr                         = NULL;
    TH1D* histoPi0TaggingEffi                                   = NULL;

    TH1D* histoRGammaMC                                         = NULL;

    //***********************************************************************************************
    //***************************  input for gamma spectrum  ****************************************
    //***********************************************************************************************
    TFile *fileGamma                                            = NULL;
    TH1D *histoInputGammaSpec                                   = NULL;

    //***********************************************************************************************
    //***************************  cocktail definition **********************************************
    //***********************************************************************************************
    TFile *fileCocktail                                         = NULL;
    TH1D *cocktailAllGamma                                      = NULL;
    TH1D *cocktailAllGammaPi0                                   = NULL;

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
        triggerCutNumber                   = fEventCutSelection(3,1);
        subTriggerCutNumber                = fEventCutSelection(4,1);
    }

    //***********************************************************************************************
    //************* Wrapper function to create output names and proper mesonType ********************
    //***********************************************************************************************
    void CreateNamingStrings(TString isMC){
        nameOutputLabel                      = "Pi0Tagging";
        if (isMC.CompareTo("kTRUE") == 0){
            nameRec                           = "recMC";
        } else {
            nameRec                           = "data";
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
