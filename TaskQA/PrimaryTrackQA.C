/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 ******     Friederike Bock, fbock@cern.ch                                 *****
 *******************************************************************************/
#include "QA.h"
//**************************************************************************************************************
//***************************** Main routine *******************************************************************
//**************************************************************************************************************
void PrimaryTrackQA(
        Int_t nSetsIn,                          // number of data sets to be analysed
        TString fEnergyFlag,                    // energy flag
        TString* DataSets,                      // technical names of data sets for output
        TString* plotDataSets,                  // labels of data sets in plots
        TString* pathDataSets,                  // path for data sets
        Int_t mode              = 2,            // standard mode for analysis
        Int_t cutNr             = -1,           // if -1: you have to choose number at runtime
        TString suffix          = "eps",        // output format of plots
        TString labelData       = "Data",       // Label for data
        Bool_t addSubfolder     = kFALSE        // flag to enable subdirectory creation for primary cut
        )
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "PrimaryTrackQA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    //**************************************************************************************************************
    //**************************** Setting general plotting style **************************************************
    //**************************************************************************************************************
    gROOT->Reset();
    TH1::AddDirectory(kFALSE);
    StyleSettingsThesis();
    SetPlotStyle();
    //**************************************************************************************************************
    //****************************** Setting common variables ******************************************************
    //**************************************************************************************************************
    const Int_t nSets                   = nSetsIn;
    // set correct mode and enable respective flags
    Int_t fMode                         = mode;
    Bool_t isCalo                       = kFALSE;
    Bool_t isConv                       = kFALSE;
    Bool_t isMC                         = kFALSE;
    Int_t iParticleType                 =0;
    // mode:    40 // new output PCM-PCM
    //          41 // new output PCM EMCAL
    //          42 // new output PCM-PHOS
    //          43 // new output PCM-DCAL
    //          44 // new output EMCAL-EMCAL
    //          45 // new output PHOS-PHOS
    //          46 // old output DCAL-DCAL
    //          47 // new output PCM-DALITZ
    //          48 // new output EMCAL-DALITZ
    //          49 // new output PHOS-DALITZ
    //          50 // new output DCAL-DALITZ
    if (fMode == 40 || fMode == 41 || fMode == 42 || fMode == 43 || fMode == 47)
        isConv                          = kTRUE;
    if (fMode == 41 || fMode == 42 || fMode == 43 || fMode == 44 || fMode == 45 || fMode == 46)
        isCalo                          = kTRUE;
    TString fCollisionSystem            = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    TString fDetectionProcess           = ReturnFullTextReconstructionProcess(fMode);
    TString fDate                       = ReturnDateString();
    TString fTextMeasurement            = Form("#omega #rightarrow #pi^{0} #pi^{+} #pi^{-}, #omega #rightarrow #pi^{0} #gamma");
    TString fTextMeasurementEta         = Form("#eta #rightarrow #pi^{0} #pi^{+} #pi^{-}");
    TString fTextMeasurementMeson[2]    = {fTextMeasurement, fTextMeasurementEta};
    const Int_t maxSets                 = 12;
    //nSets == 0 is always data!
    if(nSets>maxSets){
        cout << "Maximum hardcoded number of Data Sets in PrimaryTrackQA: " << maxSets << endl;
        cout << "You have chosen: " << nSets << ", returning!" << endl;
        return;
    }
    Color_t colorCompare[maxSets]       = {kBlack, kRed+1, kMagenta+2, 807, 800, kGreen+2, kCyan+2, kBlue+1, kOrange+2, kAzure, kViolet, kGray+1};
    TString nameMainDir[maxSets];
    Double_t processLabelOffsetX1       = 0.9;
    Double_t processLabelOffsetX2       = 0.83;
    TString* fCutSelection              = new TString[nSets];
    TString* fTypeCutSelection          = new TString[nSets];
    TString* fEventCutSelection         = new TString[nSets];
    TString* fGammaCutSelection         = new TString[nSets];
    TString* fClusterCutSelection       = new TString[nSets];
    TString* fPionCutSelection          = new TString[nSets];
    TString* fNeutralPionCutSelection   = new TString[nSets];
    TString* fMesonCutSelection         = new TString[nSets];
    Int_t iNdivisions=510;//Stores TAttAxis; histo->GetXaxis()->SetNdivisions(value,kTRUE);
    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Determine which cut to process ***************************************
    //*****************************************************************************************************
    //*****************************************************************************************************
    vector <TString> cutsTemp;
    map<TString,Int_t> mapCuts;
    TString fPionCuts2ContainerCutString;
    TString* DataSetsMC=new TString[nSets];
    for(Int_t i=0; i<nSets; i++){
        DataSetsMC[i]="True";
        TFile *fFile = new TFile(pathDataSets[i].Data(),"READ");
        if(fFile->IsZombie()){cout << "ERROR: File " << pathDataSets[i].Data() << " could not be openend! Returning..." << endl; return;}
        else{
            cout << "Processing file: " << pathDataSets[i].Data();
            TKey *key;
            TIter next(fFile->GetListOfKeys());
            while ((key=(TKey*)next())){
                cout << Form(" - found TopDir: %s",key->GetName());
                nameMainDir[i] = key->GetName();
            }
            cout << endl;
            if(nameMainDir[i].IsNull() || !nameMainDir[i].BeginsWith("Gamma")){cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir[i].Data() << "|, running in mode: " << fMode << endl; return;}
            TList *listInput =(TList*)fFile->Get(nameMainDir[i].Data());
            listInput->SetOwner(kTRUE);
            if(!listInput) {cout << "ERROR: Could not find main dir: " << nameMainDir[i].Data() << " in file! Returning..." << endl; return;}
            for(Int_t j = 0; j<listInput->GetSize(); j++){
                TList *listCuts = (TList*)listInput->At(j);
                TString nameCuts = listCuts->GetName();
                if(nameCuts.BeginsWith("Cut Number")){
                    nameCuts.Replace(0,11,"");
                    if(i==0) {
                        cutsTemp.push_back(nameCuts);
                        mapCuts[nameCuts]=1;
                    } else if( mapCuts.find(nameCuts) != mapCuts.end() ) mapCuts[nameCuts] += 1;
                }
                if(nameCuts.BeginsWith("PionCuts_")){
                    nameCuts.Replace(0,9,"");
                    fPionCuts2ContainerCutString=nameCuts.Data();
                }
            }
            delete listInput;
        }
        fFile->Close();
        delete fFile;
    }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "The following cuts are available:" << endl;
    Int_t nCuts = 0;
    vector <TString> cuts;
    for(Int_t i = 0; i < (Int_t) cutsTemp.size(); i++) {
        if(mapCuts[cutsTemp.at(i)]==nSets){ cout << Form("(%i) -- %s", i, cutsTemp.at(i).Data()) << endl; cuts.push_back(cutsTemp.at(i)); nCuts++;}
    }
    if(nCuts==0){ cout << "ERROR: No cut is available in all given data/MC sets! Returning..." << endl; return;}
    if(cutNr == -1){
        do{ cin >> cutNr;}
        while( (cutNr < 0) || (cutNr > (Int_t) cuts.size()) );
    }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Processing Cut Number: " << cutNr << endl;
    cout << cuts.at(cutNr) << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Pictures are saved as " << suffix.Data() << "!" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(Int_t i=0; i<nSets; i++){
        fCutSelection[i]            = cuts.at(cutNr);
        fTypeCutSelection[i]        = "";
        fEventCutSelection[i]       = "";
        fGammaCutSelection[i]       = "";
        fClusterCutSelection[i]     = "";
        fPionCutSelection[i]        = "";
        fNeutralPionCutSelection[i] = "";
        fMesonCutSelection[i]       = "";
        // Dirty fix, because cut string has wrong ordering!!
        fMode=ReturnSeparatedCutNumberPiPlPiMiPiZero(fCutSelection[i], fTypeCutSelection[i], fEventCutSelection[i], fGammaCutSelection[i], fClusterCutSelection[i], fNeutralPionCutSelection[i],fPionCutSelection[i], fMesonCutSelection[i]);
        if(fMode!=mode){
            cout << "ERROR: Chosen mode ("<<mode<<") is not identical to mode extracted from cutstring ("<<fMode<<")! " << endl;
            return;
        } // check if extracted mode from cutnumber is the same mode given by user
    }
    fMode                  = mode;
    cout << "\t MODE = " << mode << endl;
    Int_t fIsHeavyIonInt = -1;
    Bool_t fIsPbPb = kFALSE;
    Bool_t fIspPb = kFALSE;
    TString fIsHeavyIon = fEventCutSelection[0](0,1);
    fIsHeavyIonInt = fIsHeavyIon.Atoi();
    if(fIsHeavyIonInt > 0 & fIsHeavyIonInt < 8){
        cout << "Detected from event cuts that dataset is PbPb" << endl;
        //cout << "Will produce centrality and event plane angle histograms" << endl;
        fIsPbPb = kTRUE;
        processLabelOffsetX1 = 0.76;
        processLabelOffsetX2 = 0.58;
    } else if(fIsHeavyIonInt == 8 || fIsHeavyIonInt == 9){
        fIspPb = kTRUE;
    } else if (fIsHeavyIonInt == -1){
        cout << "ERROR detecting collision system" << endl;
    }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Obtaining trigger - ";
    TString* fTrigger       = new TString[nSets];
    for(Int_t iT=0;iT<nSets;iT++){
        fTrigger[iT] = "";
    }
    TString fTriggerCut     = fEventCutSelection[0](3,2);
    fTrigger[0]             = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), DataSets[0].Data());
    cout  << "'" << fTrigger[0].Data() << "' - was found!" << endl;
    if(fTrigger[0].Contains("not defined")){
        fTrigger[0]         = "";
        cout << "INFO: Trigger cut not defined!" << endl;
    }
    //*****************************************************************************************************
    //************************** Define output directories*************************************************
    //*****************************************************************************************************
    TString outputDir                   = Form("%s/%s/PrimaryTrackQA/%s",cuts.at(cutNr).Data(),fEnergyFlag.Data(),suffix.Data());
    TString outputDirRootFile           = Form("%s/%s/PrimaryTrackQA",cuts.at(cutNr).Data(),fEnergyFlag.Data());
    if(addSubfolder) outputDir          +=Form("/%s",DataSets[0].Data());
    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir+"/Comparison");
    gSystem->Exec("mkdir -p "+outputDir+"/Comparison/Ratios");
    gSystem->Exec("mkdir -p "+outputDir+"/Debug");
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    //*****************************************************************************************************
    //************************** Set proper cluster nomenclature ******************************************
    //*****************************************************************************************************
    TString calo            = "";
    TString fClusters       = "";
    if(isCalo){
        if(fClusterCutSelection[0].BeginsWith('1')){
            calo            = "EMCal";
            fClusters       = Form("%s clusters", calo.Data());
        } else if(fClusterCutSelection[0].BeginsWith('2')){
            calo            = "PHOS";
            fClusters       = Form("%s clusters", calo.Data());
        } else {cout << "No correct calorimeter type found: " << calo.Data() << ", returning..." << endl; return;}
    }
    //*****************************************************************************************************
    //******************************* create log file *****************************************************
    //*****************************************************************************************************
    fstream fLog;
    fLog.open(Form("%s/A-%s.log",outputDirRootFile.Data(),DataSets[0].Data()), ios::out);
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(Int_t i=0; i<nSets; i++) fLog << "Using file: " << pathDataSets[i].Data() << endl;
    fLog << fCollisionSystem.Data() << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "processed cut: " << cuts.at(cutNr).Data() << endl;
    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Vectors for Histograms ***********************************************
    //*****************************************************************************************************
    //*****************************************************************************************************
    std::vector<TH1D*> vecESD_PrimaryNegPions_Pt;
    std::vector<TH1D*> vecESD_PrimaryPosPions_Pt;
    std::vector<TH1D*> vecESD_PrimaryNegPions_Phi;
    std::vector<TH1D*> vecESD_PrimaryPosPions_Phi;
    std::vector<TH1D*> vecESD_PrimaryNegPions_Eta;
    std::vector<TH1D*> vecESD_PrimaryPosPions_Eta;
    std::vector<TH2D*> vecESD_PrimaryNegPions_ClsTPC;
    std::vector<TH2D*> vecESD_PrimaryPosPions_ClsTPC;
    std::vector<TH2D*> vecESD_PrimaryPions_DCAxy;
    std::vector<TH2D*> vecESD_PrimaryPions_DCAz;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdx;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdxSignal;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdx_LowPt;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdxSignal_LowPt;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdx_MidPt;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdxSignal_MidPt;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdx_HighPt;
    std::vector<TH2D*> vecESD_PrimaryPions_TPCdEdxSignal_HighPt;
    //AfterQA=Cut Number Pion Cut (Pion Cut folder in Cut number directory); Pre Selection=Main Dir Pion Cut
    //AfterQA
    std::vector<TH1D*> vecIsPionSelected_AfterQA;
    std::vector<TH1D*> vecdEdxCuts_AfterQA;
    std::vector<TH2D*> vecPion_ITS_after_AfterQA;
    std::vector<TH2D*> vecPion_ITS_after_AfterQA_LowPt;
    std::vector<TH2D*> vecPion_ITS_after_AfterQA_MidPt;
    std::vector<TH2D*> vecPion_ITS_after_AfterQA_HighPt;
    std::vector<TH2D*> vecPion_dEdx_after_AfterQA;
    std::vector<TH2D*> vecPion_dEdx_after_AfterQA_LowPt;
    std::vector<TH2D*> vecPion_dEdx_after_AfterQA_MidPt;
    std::vector<TH2D*> vecPion_dEdx_after_AfterQA_HighPt;
    std::vector<TH2D*> vecPion_dEdxSignal_after_AfterQA;
    std::vector<TH2D*> vecPion_dEdxSignal_after_AfterQA_LowPt;
    std::vector<TH2D*> vecPion_dEdxSignal_after_AfterQA_MidPt;
    std::vector<TH2D*> vecPion_dEdxSignal_after_AfterQA_HighPt;
    std::vector<TH2D*> vecPion_TOF_after_AfterQA;
    std::vector<TH2D*> vecPion_TOF_after_AfterQA_LowPt;
    std::vector<TH2D*> vecPion_TOF_after_AfterQA_MidPt;
    std::vector<TH2D*> vecPion_TOF_after_AfterQA_HighPt;
    std::vector<TH2D*> vechTrack_DCAxy_Pt_after_AfterQA;
    std::vector<TH2D*> vechTrack_DCAz_Pt_after_AfterQA;
    std::vector<TH2D*> vechTrack_NFindCls_Pt_TPC_after_AfterQA;
    //Pre Selection
    std::vector<TH1D*> vecIsPionSelected_PreSel;
    std::vector<TH1D*> vecdEdxCuts_PreSel;
    std::vector<TH2D*> vecPion_ITS_before_PreSel;
    std::vector<TH2D*> vecPion_ITS_before_PreSel_LowPt;
    std::vector<TH2D*> vecPion_ITS_before_PreSel_MidPt;
    std::vector<TH2D*> vecPion_ITS_before_PreSel_HighPt;
    std::vector<TH2D*> vecPion_dEdx_before_PreSel;
    std::vector<TH2D*> vecPion_dEdx_before_PreSel_LowPt;
    std::vector<TH2D*> vecPion_dEdx_before_PreSel_MidPt;
    std::vector<TH2D*> vecPion_dEdx_before_PreSel_HighPt;
    std::vector<TH2D*> vecPion_dEdxSignal_before_PreSel;
    std::vector<TH2D*> vecPion_dEdxSignal_before_PreSel_LowPt;
    std::vector<TH2D*> vecPion_dEdxSignal_before_PreSel_MidPt;
    std::vector<TH2D*> vecPion_dEdxSignal_before_PreSel_HighPt;
    std::vector<TH2D*> vecPion_TOF_before_PreSel;
    std::vector<TH2D*> vecPion_TOF_before_PreSel_LowPt;
    std::vector<TH2D*> vecPion_TOF_before_PreSel_MidPt;
    std::vector<TH2D*> vecPion_TOF_before_PreSel_HighPt;
    std::vector<TH2D*> vechTrack_DCAxy_Pt_before_PreSel;
    std::vector<TH2D*> vechTrack_DCAz_Pt_before_PreSel;
    std::vector<TH2D*> vechTrack_NFindCls_Pt_TPC_before_PreSel;
    std::vector<TH2D*> vecPion_ITS_after_PreSel;
    std::vector<TH2D*> vecPion_ITS_after_PreSel_LowPt;
    std::vector<TH2D*> vecPion_ITS_after_PreSel_MidPt;
    std::vector<TH2D*> vecPion_ITS_after_PreSel_HighPt;
    std::vector<TH2D*> vecPion_dEdx_after_PreSel;
    std::vector<TH2D*> vecPion_dEdx_after_PreSel_LowPt;
    std::vector<TH2D*> vecPion_dEdx_after_PreSel_MidPt;
    std::vector<TH2D*> vecPion_dEdx_after_PreSel_HighPt;
    std::vector<TH2D*> vecPion_dEdxSignal_after_PreSel;
    std::vector<TH2D*> vecPion_dEdxSignal_after_PreSel_LowPt;
    std::vector<TH2D*> vecPion_dEdxSignal_after_PreSel_MidPt;
    std::vector<TH2D*> vecPion_dEdxSignal_after_PreSel_HighPt;
    std::vector<TH2D*> vecPion_TOF_after_PreSel;
    std::vector<TH2D*> vecPion_TOF_after_PreSel_LowPt;
    std::vector<TH2D*> vecPion_TOF_after_PreSel_MidPt;
    std::vector<TH2D*> vecPion_TOF_after_PreSel_HighPt;
    std::vector<TH2D*> vechTrack_DCAxy_Pt_after_PreSel;
    std::vector<TH2D*> vechTrack_DCAz_Pt_after_PreSel;
    std::vector<TH2D*> vechTrack_NFindCls_Pt_TPC_after_PreSel;
    //MC histograms
    std::vector<TH1D*> vecMC_AllPosPions_Pt;
    std::vector<TH1D*> vecMC_AllNegPions_Pt;
    std::vector<TH1D*> vecMC_PosPionsFromNeutralMeson_Pt;
    std::vector<TH1D*> vecMC_NegPionsFromNeutralMeson_Pt;
    //True histograms
    std::vector<TH1D*> vecESD_TruePosPion_Pt;
    std::vector<TH1D*> vecESD_TrueNegPion_Pt;
    std::vector<TH1D*> vecESD_TruePosPionFromNeutralMeson_Pt;
    std::vector<TH1D*> vecESD_TrueNegPionFromNeutralMeson_Pt;
    //MC vs True
    std::vector<std::vector <TH1D*>> vecMCvsTrueAllPosPions_Pt(nSetsIn,std::vector<TH1D*>(2));
    std::vector<std::vector <TH1D*>> vecMCvsTrueAllNegPions_Pt(nSetsIn,std::vector<TH1D*>(2));
    std::vector<std::vector <TH1D*>> vecMCvsTruePosPionsFromNeutralMeson_Pt(nSetsIn,std::vector<TH1D*>(2));
    std::vector<std::vector <TH1D*>> vecMCvsTrueNegPionsFromNeutralMeson_Pt(nSetsIn,std::vector<TH1D*>(2));
    //---------------------------------------------Projections-----------------------------------------------------------------------
    std::vector<TH1D*> vecESD_PrimaryNegPions_ClsTPC_ProjPt;                //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPosPions_ClsTPC_ProjPt;                //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_DCAxy_ProjPt;                    //Pt was y
    std::vector<TH1D*> vecESD_PrimaryPions_DCAz_ProjPt;                     //Pt was y
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdx_ProjPt;                  //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt;            //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt;            //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt;           //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdxSignal_ProjPt;            //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt;      //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt;      //Pt was x
    std::vector<TH1D*> vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt;     //Pt was x
    //AfterQA
    std::vector<TH1D*> vecPion_ITS_after_AfterQA_ProjPt;                    //Pt was x
    std::vector<TH1D*> vecPion_ITS_after_AfterQA_LowPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_ITS_after_AfterQA_MidPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_ITS_after_AfterQA_HighPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_AfterQA_ProjPt;                   //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_AfterQA_LowPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_AfterQA_MidPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_AfterQA_HighPt_ProjPt;            //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_AfterQA_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt;       //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt;       //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt;      //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_AfterQA_ProjPt;                    //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_AfterQA_LowPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_AfterQA_MidPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_AfterQA_HighPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vechTrack_DCAxy_Pt_after_AfterQA_ProjPt;             //Pt was y
    std::vector<TH1D*> vechTrack_DCAz_Pt_after_AfterQA_ProjPt;              //Pt was y
    std::vector<TH1D*> vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt;      //Pt was x
    //Pre Selection
    std::vector<TH1D*> vecPion_ITS_before_PreSel_ProjPt;                    //Pt was x
    std::vector<TH1D*> vecPion_ITS_before_PreSel_LowPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_ITS_before_PreSel_MidPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_ITS_before_PreSel_HighPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdx_before_PreSel_ProjPt;                   //Pt was x
    std::vector<TH1D*> vecPion_dEdx_before_PreSel_LowPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdx_before_PreSel_MidPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdx_before_PreSel_HighPt_ProjPt;            //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_before_PreSel_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt;       //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt;       //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt;      //Pt was x
    std::vector<TH1D*> vecPion_TOF_before_PreSel_ProjPt;                    //Pt was x
    std::vector<TH1D*> vecPion_TOF_before_PreSel_LowPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_TOF_before_PreSel_MidPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_TOF_before_PreSel_HighPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vechTrack_DCAxy_Pt_before_PreSel_ProjPt;             //Pt was y
    std::vector<TH1D*> vechTrack_DCAz_Pt_before_PreSel_ProjPt;              //Pt was y
    std::vector<TH1D*> vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt;      //Pt was x
    std::vector<TH1D*> vecPion_ITS_after_PreSel_ProjPt;                     //Pt was x
    std::vector<TH1D*> vecPion_ITS_after_PreSel_LowPt_ProjPt;               //Pt was x
    std::vector<TH1D*> vecPion_ITS_after_PreSel_MidPt_ProjPt;               //Pt was x
    std::vector<TH1D*> vecPion_ITS_after_PreSel_HighPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_PreSel_ProjPt;                    //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_PreSel_LowPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_PreSel_MidPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_dEdx_after_PreSel_HighPt_ProjPt;             //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_PreSel_ProjPt;              //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt;        //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt;        //Pt was x
    std::vector<TH1D*> vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt;       //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_PreSel_ProjPt;                     //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_PreSel_LowPt_ProjPt;               //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_PreSel_MidPt_ProjPt;               //Pt was x
    std::vector<TH1D*> vecPion_TOF_after_PreSel_HighPt_ProjPt;              //Pt was x
    std::vector<TH1D*> vechTrack_DCAxy_Pt_after_PreSel_ProjPt;              //Pt was y
    std::vector<TH1D*> vechTrack_DCAz_Pt_after_PreSel_ProjPt;               //Pt was y
    std::vector<TH1D*> vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt;       //Pt was x
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    Double_t* nEventsAll    = new Double_t[nSets];
    Double_t* nEvents       = new Double_t[nSets];
    Int_t minB          = 0;    Int_t maxB          = 0;
    Int_t minYB         = 0;    Int_t maxYB         = 0;
    Int_t minB_SPD      = 0;    Int_t maxB_SPD      = 0;
    Int_t minYB_SPD     = 0;    Int_t maxYB_SPD     = 0;
    Bool_t isMinMaxSPD  = kTRUE;
    MesonFit fitter;
    TString StrNameOfHistogram;

    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    Double_t NmbSelectedMCPosPions     = 0.;
    Double_t NmbSelectedMCNegPions     = 0.;
    Double_t NmbSelectedPosPions       = 0.;
    Double_t NmbSelectedNegPions       = 0.;
    Double_t NmbOutPions               = 0.;
    Double_t NmbSelectedTruePosPions   = 0.;
    Double_t NmbSelectedTrueNegPions   = 0.;
    Double_t efficiency                = 0.;
    Double_t validatedefficiency       = 0.;
    Double_t purity                    = 0.;
    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Looping over DataSets ************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************
    // canvas definition
    TCanvas* canvas         = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    TCanvas* cvsQuadratic   = new TCanvas("cvsQuadratic","",10,10,500,500);  // gives the page size
    Double_t leftMargin     = 0.09;
    Double_t rightMargin    = 0.02;
    Double_t topMargin      = 0.04;
    Double_t bottomMargin   = 0.09;
    DrawGammaCanvasSettings(canvas,leftMargin,rightMargin,topMargin,bottomMargin);
    for(Int_t i=0; i<nSets; i++) {
        TFile *fFile = new TFile(pathDataSets[i].Data(),"READ");
        if(fFile->IsZombie()) {
            cout << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
            fLog << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
            return;
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // reading respective containers
        //TopDir; There should only be one; It is the main directory
        TList* TopDir                   = (TList*) fFile->Get(nameMainDir[i].Data());
        if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
        else TopDir->SetOwner(kTRUE);
        //TopContainer; The directory with the chosen cut number
        TList* TopContainer             = (TList*) TopDir->FindObject(Form("Cut Number %s",fCutSelection[i].Data()));
        if(TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",fCutSelection[i].Data()) << " not found in File" << endl; return;}
        else TopContainer->SetOwner(kTRUE);
        //ESDContainer; ESD histogram directory in the "Cut Number" directory with the chosen cut number
        TList* ESDContainer             = (TList*) TopContainer->FindObject(Form("%s ESD histograms",fCutSelection[i].Data()));
        if(ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",fCutSelection[i].Data()) << " not found in File" << endl; return;}
        else ESDContainer->SetOwner(kTRUE);
        //MCContainer not needed; MC histograms directory in the "Cut Number" directory with the chosen cut number
        TList* MCContainer            = (TList*) TopContainer->FindObject(Form("%s MC histograms",fCutSelection[i].Data()));
        if(MCContainer == NULL) {cout << "INFO: " << Form("%s MC histograms",fCutSelection[i].Data()) << " not found in File, processing data?" << endl;}
        else {
            MCContainer->SetOwner(kTRUE);
            isMC=kTRUE;
        }
        //TrueContainer; True histograms  directory in the "Cut Number" directory with the chosen cut number, only if file is MC
        TList* TrueContainer            = (TList*) TopContainer->FindObject(Form("%s True histograms",fCutSelection[i].Data()));
        if(TrueContainer == NULL) {cout << "INFO: " << Form("%s True histograms",fCutSelection[i].Data()) << " not found in File, processing data?" << endl;}
        else TrueContainer->SetOwner(kTRUE);
        //ConvEventCutsContainer; ConvEventCuts directory in the "Cut Number" directory with the chosen cut number
        TList* ConvEventCutsContainer   = (TList*) TopContainer->FindObject(Form("ConvEventCuts_%s",fEventCutSelection[i].Data()));
        if(ConvEventCutsContainer == NULL) {cout << "ERROR: " << Form("ConvEventCuts_%s",fEventCutSelection[i].Data()) << " not found in File" << endl; return;}
        else if(ConvEventCutsContainer) ConvEventCutsContainer->SetOwner(kTRUE);
        //PionCutsContainer; PionCuts directory in the "Cut Number" directory with the chosen cut number
        TString fPionCutsContainerCutString=Form("%s_%s_%s_%s_%s",fEventCutSelection[i].Data(),fGammaCutSelection[i].Data(),fPionCutSelection[i].Data(),fNeutralPionCutSelection[i].Data(),fMesonCutSelection[i].Data()); // default for mode 40
        if((fMode==41)||(fMode==42)){
            fPionCutsContainerCutString=Form("%s_%s_%s_%s_%s_%s",fEventCutSelection[i].Data(),fGammaCutSelection[i].Data(),fClusterCutSelection[i].Data(),fPionCutSelection[i].Data(),fNeutralPionCutSelection[i].Data(),fMesonCutSelection[i].Data());
        } else if ((fMode==44)||(fMode==45)){
            fPionCutsContainerCutString=Form("%s_%s_%s_%s_%s",fEventCutSelection[i].Data(),fClusterCutSelection[i].Data(),fPionCutSelection[i].Data(),fNeutralPionCutSelection[i].Data(),fMesonCutSelection[i].Data());
        }
        TList* PionCutsContainer       = (TList*) TopContainer->FindObject(Form("PionCuts_%s",fPionCutsContainerCutString.Data()));
        if(PionCutsContainer == NULL) {cout << "ERROR: " << Form("PionCuts_%s",fPionCutsContainerCutString.Data()) << " not found in File" << endl; return;}
        else if(PionCutsContainer) PionCutsContainer->SetOwner(kTRUE);
        //ConvCutsContainer; ConvCuts directory in the "Cut Number" directory with the chosen cut number
        TList* ConvCutsContainer        = (TList*) TopContainer->FindObject(Form("ConvCuts_%s",fGammaCutSelection[i].Data()));
        if(isConv && ConvCutsContainer == NULL) {cout << "ERROR: " << Form("ConvCuts_%s",fGammaCutSelection[i].Data()) << " not found in File" << endl; return;}
        else if(ConvCutsContainer) ConvCutsContainer->SetOwner(kTRUE);
        //MesonCutsContainer; ConvMesonCuts directory in the main directory
        TList* MesonCutsContainer       = (TList*) TopContainer->FindObject(Form("ConvMesonCuts_%s",fMesonCutSelection[i].Data()));
        if(MesonCutsContainer == NULL) {cout << "ERROR: " << Form("ConvMesonCuts_%s",fMesonCutSelection[i].Data()) << " not found in File" << endl; return;}
        else if(MesonCutsContainer) MesonCutsContainer->SetOwner(kTRUE);
        //MesonCutsContainer2; ConvMesonCuts directory in the main directory
        //CaloCutsContainer; CaloCuts directory in the main directory
        TList* CaloCutsContainer        = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fClusterCutSelection[i].Data()));
        if(isCalo && CaloCutsContainer == NULL) {cout << "ERROR: " << Form("CaloCuts_%s",fClusterCutSelection[i].Data()) << " not found in File" << endl; return;}
        else if(CaloCutsContainer) CaloCutsContainer->SetOwner(kTRUE);
        TList* TopContainerGamma        = NULL;
        TString ContainerGammaCut       = "";
        //PionCuts2: Main directory pion cuts
        //TString fPionCuts2ContainerCutString=Form("%s_%s_%s_%s_%s",fEventCutSelection[i].Data(),fGammaCutSelection[i].Data(),fNeutralPionCutSelection[i].Data(),fPionCutSelection[i].Data(),fMesonCutSelection[i].Data());
        //fPionCuts2ContainerCutString=Form("000000200");
        TList* PionCuts2Container       = (TList*) TopDir->FindObject(Form("PionCuts_%s",fPionCuts2ContainerCutString.Data()));
        if(PionCuts2Container == NULL) {cout << "ERROR: " << Form("PionCuts_%s",fPionCuts2ContainerCutString.Data()) << " not found in File" << endl; return;}
        else if(PionCuts2Container) PionCuts2Container->SetOwner(kTRUE);
        if(isConv){
            TList *listCuts=NULL; TString name=""; Int_t j=0;
            do{ listCuts = (TList*)TopDir->At(j); name = listCuts->GetName(); } while(!name.BeginsWith("ConvCuts_") && ++j<TopDir->GetSize());
            if(j>=TopDir->GetSize()) TopContainerGamma = NULL;
            else TopContainerGamma  = listCuts;
            if(TopContainerGamma == NULL) {cout << "ERROR: " << "ConvCuts_*" << " not found in File" << endl; return;}
            else TopContainerGamma->SetOwner(kTRUE);
            ContainerGammaCut = TopContainerGamma->GetName();
            ContainerGammaCut.Remove(0,9);
        }
        TList* TopContainerEvent = NULL;
        TString ContainerEventCut = "";
        TList* listCuts=NULL; TString name=""; Int_t j=0;
        do{ listCuts = (TList*)TopDir->At(j); name = listCuts->GetName(); } while(!name.BeginsWith("ConvEventCuts_") && ++j<TopDir->GetSize());
        if(j>=TopDir->GetSize()) TopContainerEvent = NULL;
        else TopContainerEvent = listCuts;

        if(TopContainerEvent == NULL) {cout << "ERROR: " << "ConvEventCuts_*" << " not found in File" << endl; return;}
        else TopContainerEvent->SetOwner(kTRUE);
        ContainerEventCut = TopContainerEvent->GetName();
        ContainerEventCut.Remove(0,14);
        //-------------------------------------------------------------------------------------------------------------------------------
        nEvents[i]          = 0;
        nEventsAll[i]       = 0;
        TH1D* fHistNEvents  = NULL;
        fHistNEvents        = (TH1D*)ESDContainer->FindObject("NEventsWOWeight");
        if (fHistNEvents ){
            cout << "INFO: Output contains event weights" << endl;
        } else {
            fHistNEvents     = (TH1D*)ESDContainer->FindObject("NEvents");
        }
        //if(plotDataSets[i].Contains("JetJet") || plotDataSets[i].Contains("jetjet")) fHistNEvents = (TH1D*)ESDContainer->FindObject("NEventsWOWeight");
        if(fHistNEvents){
            if(fIsPbPb || fIspPb )
                nEvents[i]     = (Double_t) fHistNEvents->GetBinContent(1);
            else
                nEvents[i]     = (Double_t) GetNEvents(fHistNEvents,kFALSE);
            nEventsAll[i]      = fHistNEvents->GetEntries() - fHistNEvents->GetBinContent(4);
        }
        else{
            cout << "ERROR: Object |fHistNEvents| could not be found in File '" << pathDataSets[i].Data() << "'! Returning..." << endl;
            return;
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        cout << endl;
        cout << "----------------------------------------------------------------------------" << endl;
        cout << "Processing file: " << pathDataSets[i].Data() << endl;
        cout << "Set: " << plotDataSets[i].Data() << endl;
        cout << "NEventsAll: '" << nEventsAll[i] << "', NEvents (for normalization): '" << nEvents[i] << "'" << endl;
        cout << "----------------------------------------------------------------------------" << endl;
        cout << endl;
        fLog << "----------------------------------------------------------------------------" << endl;
        fLog << "Processing file: " << pathDataSets[i].Data() << endl;
        fLog << "Set: " << plotDataSets[i].Data() << endl;
        fLog << "NEventsAll: '" << nEventsAll[i] << "', NEvents (for normalization): '" << nEvents[i] << "'" << endl;
        fLog << "----------------------------------------------------------------------------" << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        if( nEvents[i] < 1. ){cout << "ERROR: number of accepted events in data set is <1: " << nEvents[i] << "! Returning..." << endl; return;}
        //-------------------------------------------------------------------------------------------------------------------------------
        const char* nameOutput = Form("%s/PrimaryTrackQA_%s.root",outputDirRootFile.Data(),DataSets[i].Data());
        TFile* fOutput = new TFile(nameOutput,"RECREATE");
        cout << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Output file: " << nameOutput << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "Output file: " << nameOutput << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------|Get Histograms: ESD-Histograms|----------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_PrimaryNegPions_Pt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryNegPions_Pt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistESD_PrimaryNegPions_Pt = (TH1D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryNegPions_Pt){
            GetMinMaxBin(fHistESD_PrimaryNegPions_Pt,minB,maxB);
            SetXRange(fHistESD_PrimaryNegPions_Pt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                 fHistESD_PrimaryNegPions_Pt,"","#it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistESD_PrimaryNegPions_Pt);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryNegPions_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryNegPions_Pt.push_back(new TH1D(*fHistESD_PrimaryNegPions_Pt));
            NmbSelectedNegPions = fHistESD_PrimaryNegPions_Pt->GetEntries();
        } else cout << "INFO: Object |fHistESD_PrimaryNegPions_Pt| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_PrimaryPosPions_Pt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPosPions_Pt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistESD_PrimaryPosPions_Pt = (TH1D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPosPions_Pt){
            GetMinMaxBin(fHistESD_PrimaryPosPions_Pt,minB,maxB);
            SetXRange(fHistESD_PrimaryPosPions_Pt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPosPions_Pt,"","#it{p}_{T, #it{#pi^{+}}} (GeV/#it{c})","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistESD_PrimaryPosPions_Pt);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPosPions_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPosPions_Pt.push_back(new TH1D(*fHistESD_PrimaryPosPions_Pt));
            NmbSelectedPosPions = fHistESD_PrimaryPosPions_Pt->GetEntries();
        } else cout << "INFO: Object |fHistESD_PrimaryPosPions_Pt| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_PrimaryNegPions_Phi
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryNegPions_Phi";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistESD_PrimaryNegPions_Phi = (TH1D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryNegPions_Phi){
            GetMinMaxBin(fHistESD_PrimaryNegPions_Phi,minB,maxB);
            SetXRange(fHistESD_PrimaryNegPions_Phi,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistESD_PrimaryNegPions_Phi,"","#it{#varphi}_{#it{#pi^{-}}}","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistESD_PrimaryNegPions_Phi);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryNegPions_Phi, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryNegPions_Phi.push_back(new TH1D(*fHistESD_PrimaryNegPions_Phi));
        } else cout << "INFO: Object |fHistESD_PrimaryNegPions_Phi| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_PrimaryPosPions_Phi
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPosPions_Phi";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistESD_PrimaryPosPions_Phi = (TH1D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPosPions_Phi){
            GetMinMaxBin(fHistESD_PrimaryPosPions_Phi,minB,maxB);
            SetXRange(fHistESD_PrimaryPosPions_Phi,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistESD_PrimaryPosPions_Phi,"","#it{#varphi}_{#it{#pi^{+}}}","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistESD_PrimaryPosPions_Phi);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPosPions_Phi, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPosPions_Phi.push_back(new TH1D(*fHistESD_PrimaryPosPions_Phi));
        } else cout << "INFO: Object |fHistESD_PrimaryPosPions_Phi| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_PrimaryNegPions_Eta
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryNegPions_Eta";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistESD_PrimaryNegPions_Eta = (TH1D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryNegPions_Eta){
            GetMinMaxBin(fHistESD_PrimaryNegPions_Eta,minB,maxB);
            SetXRange(fHistESD_PrimaryNegPions_Eta,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistESD_PrimaryNegPions_Eta,"","#it{#eta}_{#it{#pi^{-}}}","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistESD_PrimaryNegPions_Eta);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryNegPions_Eta, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryNegPions_Eta.push_back(new TH1D(*fHistESD_PrimaryNegPions_Eta));
        } else cout << "INFO: Object |fHistESD_PrimaryNegPions_Eta| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_PrimaryPosPions_Eta
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPosPions_Eta";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistESD_PrimaryPosPions_Eta = (TH1D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPosPions_Eta){
            GetMinMaxBin(fHistESD_PrimaryPosPions_Eta,minB,maxB);
            SetXRange(fHistESD_PrimaryPosPions_Eta,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistESD_PrimaryPosPions_Eta,"","#it{#eta}_{#it{#pi^{+}}}","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistESD_PrimaryPosPions_Eta);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPosPions_Eta, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPosPions_Eta.push_back(new TH1D(*fHistESD_PrimaryPosPions_Eta));
        } else cout << "INFO: Object |fHistESD_PrimaryPosPions_Eta| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryNegPions_ClsTPC
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryNegPions_ClsTPC";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryNegPions_ClsTPC = (TH2D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryNegPions_ClsTPC){
            TH2D* fHistESD_PrimaryNegPions_ClsTPC_switchedAxis=SwitchTF2DAxis(fHistESD_PrimaryNegPions_ClsTPC);
            GetMinMaxBin(fHistESD_PrimaryNegPions_ClsTPC_switchedAxis,minB,maxB);
            SetXRange(fHistESD_PrimaryNegPions_ClsTPC_switchedAxis,minB-1,maxB+1);
            GetMinMaxBinY(fHistESD_PrimaryNegPions_ClsTPC_switchedAxis,minYB,maxYB);
            SetYRange(fHistESD_PrimaryNegPions_ClsTPC_switchedAxis,minYB,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryNegPions_ClsTPC_switchedAxis,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistESD_PrimaryNegPions_ClsTPC_switchedAxis,"", "#it{p}_{T, #pi^{-}} (GeV/#it{c}) TPC",
                                 "Findable Clusters",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryNegPions_ClsTPC_switchedAxis, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryNegPions_ClsTPC.push_back(new TH2D(*fHistESD_PrimaryNegPions_ClsTPC_switchedAxis));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt= (TH1D*)fHistESD_PrimaryNegPions_ClsTPC_switchedAxis->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryNegPions_ClsTPC_switchedAxis->GetYaxis()->GetFirst(),fHistESD_PrimaryNegPions_ClsTPC_switchedAxis->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt,"",fHistESD_PrimaryNegPions_ClsTPC_switchedAxis->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryNegPions_ClsTPC_ProjPt.push_back(new TH1D(*fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt));
            delete fHistESD_PrimaryNegPions_ClsTPC_switchedAxis;
            fHistESD_PrimaryNegPions_ClsTPC_switchedAxis=NULL;
            delete fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt;
            fHistESD_PrimaryNegPions_ClsTPC_switchedAxis_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryNegPions_ClsTPC %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPosPions_ClsTPC
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPosPions_ClsTPC";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPosPions_ClsTPC = (TH2D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPosPions_ClsTPC){
            TH2D* fHistESD_PrimaryPosPions_ClsTPC_switchedAxis=SwitchTF2DAxis(fHistESD_PrimaryPosPions_ClsTPC);
            GetMinMaxBin(fHistESD_PrimaryPosPions_ClsTPC_switchedAxis,minB,maxB);
            SetXRange(fHistESD_PrimaryPosPions_ClsTPC_switchedAxis,minB-1,maxB+1);
            GetMinMaxBinY(fHistESD_PrimaryPosPions_ClsTPC_switchedAxis,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPosPions_ClsTPC_switchedAxis,minYB,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPosPions_ClsTPC_switchedAxis,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPosPions_ClsTPC_switchedAxis,"","#it{p}_{T, #pi^{+}} (GeV/#it{c}) TPC",
                                 "Findable Clusters",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPosPions_ClsTPC_switchedAxis, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPosPions_ClsTPC.push_back(new TH2D(*fHistESD_PrimaryPosPions_ClsTPC_switchedAxis));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt= (TH1D*)fHistESD_PrimaryPosPions_ClsTPC_switchedAxis->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPosPions_ClsTPC_switchedAxis->GetYaxis()->GetFirst(),fHistESD_PrimaryPosPions_ClsTPC_switchedAxis->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt,"",fHistESD_PrimaryPosPions_ClsTPC_switchedAxis->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPosPions_ClsTPC_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt));
            delete fHistESD_PrimaryPosPions_ClsTPC_switchedAxis;
            fHistESD_PrimaryPosPions_ClsTPC_switchedAxis=NULL;
            delete fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt;
            fHistESD_PrimaryPosPions_ClsTPC_switchedAxis_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPosPions_ClsTPC %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_DCAxy
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_DCAxy";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_DCAxy = (TH2D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPions_DCAxy){
            GetMinMaxBin(fHistESD_PrimaryPions_DCAxy,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_DCAxy,minB-1,maxB+1);
            GetMinMaxBinY(fHistESD_PrimaryPions_DCAxy,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_DCAxy,minYB,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_DCAxy,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_DCAxy,"","DCA_{#it{xy}} (cm)",
                                 "#it{p}_{T, #pi} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_DCAxy, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_DCAxy.push_back(new TH2D(*fHistESD_PrimaryPions_DCAxy));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_DCAxy_ProjPt= (TH1D*)fHistESD_PrimaryPions_DCAxy->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_DCAxy->GetYaxis()->GetFirst(),fHistESD_PrimaryPions_DCAxy->GetYaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_DCAxy_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_DCAxy_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_DCAxy_ProjPt,"",fHistESD_PrimaryPions_DCAxy->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_DCAxy_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_DCAxy_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_DCAxy_ProjPt));
            delete fHistESD_PrimaryPions_DCAxy_ProjPt;
            fHistESD_PrimaryPions_DCAxy_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPions_DCAxy %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_DCAz
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_DCAz";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_DCAz = (TH2D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPions_DCAz){
            GetMinMaxBin(fHistESD_PrimaryPions_DCAz,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_DCAz,minB-1,maxB+1);
            GetMinMaxBinY(fHistESD_PrimaryPions_DCAz,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_DCAz,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_DCAz,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_DCAz,"","DCA_{#it{z}} (cm)",
                                 "#it{p}_{T, #pi} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_DCAz, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_DCAz.push_back(new TH2D(*fHistESD_PrimaryPions_DCAz));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_DCAz_ProjPt= (TH1D*)fHistESD_PrimaryPions_DCAz->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_DCAz->GetYaxis()->GetFirst(),fHistESD_PrimaryPions_DCAz->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_DCAz_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_DCAz_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_DCAz_ProjPt,"",fHistESD_PrimaryPions_DCAz->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_DCAz_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_DCAz_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_DCAz_ProjPt));
            delete fHistESD_PrimaryPions_DCAz_ProjPt;
            fHistESD_PrimaryPions_DCAz_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPions_DCAz %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdx = (TH2D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPions_TPCdEdx){
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx,minB-1,maxB+1);
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdx,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdx,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdx,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdx,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdx, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdx));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdx_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdx->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdx_ProjPt,"",fHistESD_PrimaryPions_TPCdEdx->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdx_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdx_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdx_ProjPt;
            fHistESD_PrimaryPions_TPCdEdx_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPions_TPCdEdx %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx_LowPt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdx_LowPt =new TH2D(*fHistESD_PrimaryPions_TPCdEdx);
        if(fHistESD_PrimaryPions_TPCdEdx_LowPt){
            fHistESD_PrimaryPions_TPCdEdx_LowPt->SetName(StrNameOfHistogram.Data());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_LowPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx_LowPt,0,fHistESD_PrimaryPions_TPCdEdx_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdx_LowPt,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdx_LowPt,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdx_LowPt,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdx_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdx_LowPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx_LowPt.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdx_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdx_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdx_LowPt->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt,"",fHistESD_PrimaryPions_TPCdEdx_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt;
            fHistESD_PrimaryPions_TPCdEdx_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object | ESD_PrimaryPions_TPCdEdx %s| could not be found! Skipping Draw for LowPt...", fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx_MidPt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdx_MidPt =new TH2D(*fHistESD_PrimaryPions_TPCdEdx);
        if(fHistESD_PrimaryPions_TPCdEdx_MidPt){
            fHistESD_PrimaryPions_TPCdEdx_MidPt->SetName(StrNameOfHistogram.Data());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_MidPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx_MidPt,fHistESD_PrimaryPions_TPCdEdx_MidPt->GetXaxis()->FindBin(2),fHistESD_PrimaryPions_TPCdEdx_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdx_MidPt,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdx_MidPt,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdx_MidPt,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdx_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdx_MidPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx_MidPt.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdx_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdx_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdx_MidPt->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt,"",fHistESD_PrimaryPions_TPCdEdx_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt;
            fHistESD_PrimaryPions_TPCdEdx_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object | ESD_PrimaryPions_TPCdEdx %s| could not be found! Skipping Draw for MidPt...", fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx_HighPt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdx_HighPt =new TH2D(*fHistESD_PrimaryPions_TPCdEdx);
        if(fHistESD_PrimaryPions_TPCdEdx_HighPt){
            fHistESD_PrimaryPions_TPCdEdx_HighPt->SetName(StrNameOfHistogram.Data());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_HighPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx_HighPt,fHistESD_PrimaryPions_TPCdEdx_HighPt->GetXaxis()->FindBin(3),fHistESD_PrimaryPions_TPCdEdx_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdx_HighPt,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdx_HighPt,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdx_HighPt,1,maxB+1,minB-1,maxB+1);
            iNdivisions=fHistESD_PrimaryPions_TPCdEdx_HighPt->GetNdivisions();
            fHistESD_PrimaryPions_TPCdEdx_HighPt->GetXaxis()->SetNdivisions(300000,kFALSE);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdx_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdx_HighPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx_HighPt.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdx_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdx_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdx_HighPt->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt,"",fHistESD_PrimaryPions_TPCdEdx_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt;
            fHistESD_PrimaryPions_TPCdEdx_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object | ESD_PrimaryPions_TPCdEdx %s| could not be found! Skipping Draw for HighPt...", fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdxSignal = (TH2D*)ESDContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
        if(fHistESD_PrimaryPions_TPCdEdxSignal){
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal,minB-1,maxB+1);
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdxSignal,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdxSignal,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdxSignal,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdxSignal, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdxSignal));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdxSignal->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt,"",fHistESD_PrimaryPions_TPCdEdxSignal->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt;
            fHistESD_PrimaryPions_TPCdEdxSignal_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPions_TPCdEdxSignal %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal_LowPt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdxSignal_LowPt = new TH2D(*fHistESD_PrimaryPions_TPCdEdxSignal);
        if(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt){
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt,0,fHistESD_PrimaryPions_TPCdEdxSignal_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdxSignal_LowPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal_LowPt.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdxSignal_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdxSignal_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdxSignal_LowPt->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt,"",fHistESD_PrimaryPions_TPCdEdxSignal_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt;
            fHistESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPions_TPCdEdxSignal_LowPt %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal_MidPt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdxSignal_MidPt = new TH2D(*fHistESD_PrimaryPions_TPCdEdxSignal);
        if(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt){
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt,fHistESD_PrimaryPions_TPCdEdxSignal_MidPt->GetXaxis()->FindBin(2),fHistESD_PrimaryPions_TPCdEdxSignal_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdxSignal_MidPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal_MidPt.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdxSignal_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdxSignal_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdxSignal_MidPt->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt,"",fHistESD_PrimaryPions_TPCdEdxSignal_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt;
            fHistESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPions_TPCdEdxSignal_MidPt %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal_HighPt
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistESD_PrimaryPions_TPCdEdxSignal_HighPt = new TH2D(*fHistESD_PrimaryPions_TPCdEdxSignal);
        if(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt){
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt,fHistESD_PrimaryPions_TPCdEdxSignal_HighPt->GetXaxis()->FindBin(3),fHistESD_PrimaryPions_TPCdEdxSignal_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt,minYB,maxYB);
            SetYRange(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt,minYB-1,maxYB+1);
            SetZMinMaxTH2(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt,1,maxB+1,minB-1,maxB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistESD_PrimaryPions_TPCdEdxSignal_HighPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal_HighPt.push_back(new TH2D(*fHistESD_PrimaryPions_TPCdEdxSignal_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt= (TH1D*)fHistESD_PrimaryPions_TPCdEdxSignal_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistESD_PrimaryPions_TPCdEdxSignal_HighPt->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt,"",fHistESD_PrimaryPions_TPCdEdxSignal_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.push_back(new TH1D(*fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt));
            delete fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt;
            fHistESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |ESD_PrimaryPions_TPCdEdxSignal_HighPt %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------|Get Histograms: MC-Histograms|-----------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        if (MCContainer){
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Drawing MC Histograms" << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // MC_AllPosPions_Pt
            if (iParticleType==0){StrNameOfHistogram="MC_AllPosPions_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistMC_AllPosPions_Pt = (TH1D*)MCContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistMC_AllPosPions_Pt){
                GetMinMaxBin(fHistMC_AllPosPions_Pt,minB,maxB);
                SetXRange(fHistMC_AllPosPions_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kFALSE,kFALSE,
                                     fHistMC_AllPosPions_Pt,"","#it{p}_{T, #pi^{+}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistMC_AllPosPions_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistMC_AllPosPions_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecMC_AllPosPions_Pt.push_back(new TH1D(*fHistMC_AllPosPions_Pt));
                NmbSelectedMCPosPions = fHistMC_AllPosPions_Pt->GetEntries();
            } else cout << "INFO: Object |: MC_AllPosPions_Pt| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // MC_AllNegPions_Pt
            if (iParticleType==0){StrNameOfHistogram="MC_AllNegPions_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistMC_AllNegPions_Pt = (TH1D*)MCContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistMC_AllNegPions_Pt){
                GetMinMaxBin(fHistMC_AllNegPions_Pt,minB,maxB);
                SetXRange(fHistMC_AllNegPions_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kFALSE,kFALSE,
                                     fHistMC_AllNegPions_Pt,"","#it{p}_{T, #pi^{-}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistMC_AllNegPions_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistMC_AllNegPions_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecMC_AllNegPions_Pt.push_back(new TH1D(*fHistMC_AllNegPions_Pt));
                NmbSelectedMCNegPions = fHistMC_AllNegPions_Pt->GetEntries();
            } else cout << "INFO: Object |: MC_AllNegPions_Pt| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // MC_PosPionsFromNeutralMeson_Pt
            if (iParticleType==0){StrNameOfHistogram="MC_PosPionsFromNeutralMeson_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistMC_PosPionsFromNeutralMeson_Pt = (TH1D*)MCContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistMC_PosPionsFromNeutralMeson_Pt){
                GetMinMaxBin(fHistMC_PosPionsFromNeutralMeson_Pt,minB,maxB);
                SetXRange(fHistMC_PosPionsFromNeutralMeson_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kFALSE,kFALSE,
                                     fHistMC_PosPionsFromNeutralMeson_Pt,"","#it{p}_{T, #pi^{+}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistMC_PosPionsFromNeutralMeson_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistMC_PosPionsFromNeutralMeson_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecMC_PosPionsFromNeutralMeson_Pt.push_back(new TH1D(*fHistMC_PosPionsFromNeutralMeson_Pt));
            } else cout << "INFO: Object |: MC_PosPionsFromNeutralMeson_Pt| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // MC_NegPionsFromNeutralMeson_Pt
            if (iParticleType==0){StrNameOfHistogram="MC_NegPionsFromNeutralMeson_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistMC_NegPionsFromNeutralMeson_Pt = (TH1D*)MCContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistMC_NegPionsFromNeutralMeson_Pt){
                GetMinMaxBin(fHistMC_NegPionsFromNeutralMeson_Pt,minB,maxB);
                SetXRange(fHistMC_NegPionsFromNeutralMeson_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kFALSE,kFALSE,
                                     fHistMC_NegPionsFromNeutralMeson_Pt,"","#it{p}_{T, #pi^{-}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistMC_NegPionsFromNeutralMeson_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistMC_NegPionsFromNeutralMeson_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecMC_NegPionsFromNeutralMeson_Pt.push_back(new TH1D(*fHistMC_NegPionsFromNeutralMeson_Pt));
            } else cout << "INFO: Object |: MC_NegPionsFromNeutralMeson_Pt| could not be found! Skipping Draw..." << endl;
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------|Get Histograms: True Histograms|---------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        if (TrueContainer){
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "Drawing True Histograms" << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // ESD_TruePosPion_Pt
            if (iParticleType==0){StrNameOfHistogram="ESD_TruePosPion_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistESD_TruePosPion_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistESD_TruePosPion_Pt){
                GetMinMaxBin(fHistESD_TruePosPion_Pt,minB,maxB);
                SetXRange(fHistESD_TruePosPion_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                     fHistESD_TruePosPion_Pt,"","#it{p}_{T, #pi^{+}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistESD_TruePosPion_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistESD_TruePosPion_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecESD_TruePosPion_Pt.push_back(new TH1D(*fHistESD_TruePosPion_Pt));
                NmbSelectedTruePosPions = fHistESD_TruePosPion_Pt->GetEntries();
            } else cout << "INFO: Object |: ESD_TruePosPion_Pt| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // ESD_TrueNegPion_Pt
            if (iParticleType==0){StrNameOfHistogram="ESD_TrueNegPion_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistESD_TrueNegPion_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistESD_TrueNegPion_Pt){
                GetMinMaxBin(fHistESD_TrueNegPion_Pt,minB,maxB);
                SetXRange(fHistESD_TrueNegPion_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                     fHistESD_TrueNegPion_Pt,"","#it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistESD_TrueNegPion_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistESD_TrueNegPion_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecESD_TrueNegPion_Pt.push_back(new TH1D(*fHistESD_TrueNegPion_Pt));
                NmbSelectedTrueNegPions = fHistESD_TrueNegPion_Pt->GetEntries();
            } else cout << "INFO: Object |: ESD_TrueNegPion_Pt| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // ESD_TruePosPionFromNeutralMeson_Pt
            if (iParticleType==0){StrNameOfHistogram="ESD_TruePosPionFromNeutralMeson_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistESD_TruePosPionFromNeutralMeson_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistESD_TruePosPionFromNeutralMeson_Pt){
                GetMinMaxBin(fHistESD_TruePosPionFromNeutralMeson_Pt,minB,maxB);
                SetXRange(fHistESD_TruePosPionFromNeutralMeson_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                     fHistESD_TruePosPionFromNeutralMeson_Pt,"","#it{p}_{T, #it{#pi^{+}}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistESD_TruePosPionFromNeutralMeson_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistESD_TruePosPionFromNeutralMeson_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecESD_TruePosPionFromNeutralMeson_Pt.push_back(new TH1D(*fHistESD_TruePosPionFromNeutralMeson_Pt));
            } else cout << "INFO: Object |: ESD_TruePosPionFromNeutralMeson_Pt| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
            // ESD_TrueNegPionFromNeutralMeson_Pt
            if (iParticleType==0){StrNameOfHistogram="ESD_TrueNegPionFromNeutralMeson_Pt";}
            if (iParticleType==1){StrNameOfHistogram="";}
            TH1D* fHistESD_TrueNegPionFromNeutralMeson_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",StrNameOfHistogram.Data()));
            if(fHistESD_TrueNegPionFromNeutralMeson_Pt){
                GetMinMaxBin(fHistESD_TrueNegPionFromNeutralMeson_Pt,minB,maxB);
                SetXRange(fHistESD_TrueNegPionFromNeutralMeson_Pt,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                     fHistESD_TrueNegPionFromNeutralMeson_Pt,"","#it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})","# Entries",1,1,
                                     processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                //WriteHistogram(fHistESD_TrueNegPionFromNeutralMeson_Pt);
                SaveCanvasAndWriteHistogram(canvas, fHistESD_TrueNegPionFromNeutralMeson_Pt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
                vecESD_TrueNegPionFromNeutralMeson_Pt.push_back(new TH1D(*fHistESD_TrueNegPionFromNeutralMeson_Pt));
            } else cout << "INFO: Object |: ESD_TrueNegPionFromNeutralMeson_Pt| could not be found! Skipping Draw..." << endl;
            //-------------------------------------------------------------------------------------------------------------------------------
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        //-------------------------|Get Histograms: AfterQA (Pion Cut folder in Cut number directory)|--------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Drawing Pion Cut Histograms AfterQA" << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // PionCutsContainer: IsPionSelected
        if (iParticleType==0){StrNameOfHistogram="IsPionSelected";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistIsPionSelected_AfterQA = (TH1D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHistIsPionSelected_AfterQA){
            GetMinMaxBin(fHistIsPionSelected_AfterQA,minB,maxB);
            SetXRange(fHistIsPionSelected_AfterQA,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistIsPionSelected_AfterQA,"","","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistIsPionSelected_AfterQA);
            SaveCanvasAndWriteHistogram(canvas, fHistIsPionSelected_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecIsPionSelected_AfterQA.push_back(new TH1D(*fHistIsPionSelected_AfterQA));
            NmbOutPions = fHistIsPionSelected_AfterQA->GetBinContent(4);
        } else cout << "INFO: Object |AfterQA: fHistIsPionSelected_AfterQA| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // PionCutsContainer: dEdxCuts
        if (iParticleType==0){StrNameOfHistogram="dEdxCuts";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistdEdxCuts_AfterQA = (TH1D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHistdEdxCuts_AfterQA){
            GetMinMaxBin(fHistdEdxCuts_AfterQA,minB,maxB);
            SetXRange(fHistdEdxCuts_AfterQA,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistdEdxCuts_AfterQA,"","","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistdEdxCuts_AfterQA);
            SaveCanvasAndWriteHistogram(canvas, fHistdEdxCuts_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecdEdxCuts_AfterQA.push_back(new TH1D(*fHistdEdxCuts_AfterQA));
        } else cout << "INFO: Object |AfterQA: fHistdEdxCuts_AfterQA| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_ITS_after
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_AfterQA = (TH2D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHistPion_ITS_after_AfterQA){
            GetMinMaxBin(fHistPion_ITS_after_AfterQA,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_ITS_after_AfterQA,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_AfterQA,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_AfterQA,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_AfterQA,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA.push_back(new TH2D(*fHistPion_ITS_after_AfterQA));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_AfterQA_ProjPt= (TH1D*)fHistPion_ITS_after_AfterQA->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_AfterQA_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_AfterQA_ProjPt,"",fHistPion_ITS_after_AfterQA->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_AfterQA_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_AfterQA_ProjPt));
            delete fHistPion_ITS_after_AfterQA_ProjPt;
            fHistPion_ITS_after_AfterQA_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_ITS_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_ITS_after_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_AfterQA_LowPt = new TH2D (*fHistPion_ITS_after_AfterQA);
        if(fHistPion_ITS_after_AfterQA_LowPt){
            fHistPion_ITS_after_AfterQA_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_ITS_after_AfterQA_LowPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA_LowPt,0,fHistPion_ITS_after_AfterQA_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_ITS_after_AfterQA_LowPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_AfterQA_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_AfterQA_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_AfterQA_LowPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_AfterQA_LowPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA_LowPt.push_back(new TH2D(*fHistPion_ITS_after_AfterQA_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_AfterQA_LowPt_ProjPt= (TH1D*)fHistPion_ITS_after_AfterQA_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_AfterQA_LowPt->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_AfterQA_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_AfterQA_LowPt_ProjPt,"",fHistPion_ITS_after_AfterQA_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_AfterQA_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA_LowPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_AfterQA_LowPt_ProjPt));
            delete fHistPion_ITS_after_AfterQA_LowPt_ProjPt;
            fHistPion_ITS_after_AfterQA_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_ITS_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_ITS_after_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_AfterQA_MidPt = new TH2D (*fHistPion_ITS_after_AfterQA);
        if(fHistPion_ITS_after_AfterQA_MidPt){
            fHistPion_ITS_after_AfterQA_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_ITS_after_AfterQA_MidPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA_MidPt,fHistPion_ITS_after_AfterQA_MidPt->GetXaxis()->FindBin(2),fHistPion_ITS_after_AfterQA_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_ITS_after_AfterQA_MidPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_AfterQA_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_AfterQA_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_AfterQA_MidPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_AfterQA_MidPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA_MidPt.push_back(new TH2D(*fHistPion_ITS_after_AfterQA_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_AfterQA_MidPt_ProjPt= (TH1D*)fHistPion_ITS_after_AfterQA_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_AfterQA_MidPt->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_AfterQA_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_AfterQA_MidPt_ProjPt,"",fHistPion_ITS_after_AfterQA_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_AfterQA_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA_MidPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_AfterQA_MidPt_ProjPt));
            delete fHistPion_ITS_after_AfterQA_MidPt_ProjPt;
            fHistPion_ITS_after_AfterQA_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_ITS_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_ITS_after_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_AfterQA_HighPt = new TH2D (*fHistPion_ITS_after_AfterQA);
        if(fHistPion_ITS_after_AfterQA_HighPt){
            fHistPion_ITS_after_AfterQA_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_ITS_after_AfterQA_HighPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA_HighPt,fHistPion_ITS_after_AfterQA_MidPt->GetXaxis()->FindBin(3),fHistPion_ITS_after_AfterQA_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_ITS_after_AfterQA_HighPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_AfterQA_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_AfterQA_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_AfterQA_HighPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_AfterQA_HighPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA_HighPt.push_back(new TH2D(*fHistPion_ITS_after_AfterQA_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_AfterQA_HighPt_ProjPt= (TH1D*)fHistPion_ITS_after_AfterQA_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_AfterQA_HighPt->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_AfterQA_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_AfterQA_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_AfterQA_HighPt_ProjPt,"",fHistPion_ITS_after_AfterQA_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_AfterQA_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_AfterQA_HighPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_AfterQA_HighPt_ProjPt));
            delete fHistPion_ITS_after_AfterQA_HighPt_ProjPt;
            fHistPion_ITS_after_AfterQA_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_ITS_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_AfterQA = (TH2D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHistPion_dEdx_after_AfterQA){
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_dEdx_after_AfterQA,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_AfterQA,minYB-1,maxYB+1);
            //cout<<"minYB: "<<minYB<<"; maxYB: "<<maxYB<<endl;
            //SetYRange(fHistPion_dEdx_after_AfterQA,-2,2);
            SetZMinMaxTH2(fHistPion_dEdx_after_AfterQA,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_AfterQA,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA.push_back(new TH2D(*fHistPion_dEdx_after_AfterQA));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_AfterQA_ProjPt= (TH1D*)fHistPion_dEdx_after_AfterQA->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_AfterQA_ProjPt,"",fHistPion_dEdx_after_AfterQA->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_AfterQA_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_AfterQA_ProjPt));
            delete fHistPion_dEdx_after_AfterQA_ProjPt;
            fHistPion_dEdx_after_AfterQA_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdx_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_AfterQA_LowPt = new TH2D (*fHistPion_dEdx_after_AfterQA);
        if(fHistPion_dEdx_after_AfterQA_LowPt){
            fHistPion_dEdx_after_AfterQA_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA_LowPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA_LowPt,0,fHistPion_dEdx_after_AfterQA_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_dEdx_after_AfterQA_LowPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_AfterQA_LowPt,minYB-1,maxYB+1);
            //cout<<"minYB: "<<minYB<<"; maxYB: "<<maxYB<<endl;
            //SetYRange(fHistPion_dEdx_after_AfterQA_LowPt,-2,2);
            SetZMinMaxTH2(fHistPion_dEdx_after_AfterQA_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_AfterQA_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_AfterQA_LowPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA_LowPt.push_back(new TH2D(*fHistPion_dEdx_after_AfterQA_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_AfterQA_LowPt_ProjPt= (TH1D*)fHistPion_dEdx_after_AfterQA_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_AfterQA_LowPt->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_AfterQA_LowPt_ProjPt,"",fHistPion_dEdx_after_AfterQA_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_AfterQA_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA_LowPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_AfterQA_LowPt_ProjPt));
            delete fHistPion_dEdx_after_AfterQA_LowPt_ProjPt;
            fHistPion_dEdx_after_AfterQA_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdx_after_LowPt %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_AfterQA_MidPt = new TH2D (*fHistPion_dEdx_after_AfterQA);
        if(fHistPion_dEdx_after_AfterQA_MidPt){
            fHistPion_dEdx_after_AfterQA_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA_MidPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA_MidPt,fHistPion_dEdx_after_AfterQA_MidPt->GetXaxis()->FindBin(2),fHistPion_dEdx_after_AfterQA_MidPt->GetXaxis()->FindBin(2));
            GetMinMaxBinY(fHistPion_dEdx_after_AfterQA_MidPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_AfterQA_MidPt,minYB-1,maxYB+1);
            //cout<<"minYB: "<<minYB<<"; maxYB: "<<maxYB<<endl;
            //SetYRange(fHistPion_dEdx_after_AfterQA_MidPt,-2,2);
            SetZMinMaxTH2(fHistPion_dEdx_after_AfterQA_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_AfterQA_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_AfterQA_MidPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA_MidPt.push_back(new TH2D(*fHistPion_dEdx_after_AfterQA_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_AfterQA_MidPt_ProjPt= (TH1D*)fHistPion_dEdx_after_AfterQA_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_AfterQA_MidPt->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_AfterQA_MidPt_ProjPt,"",fHistPion_dEdx_after_AfterQA_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_AfterQA_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA_MidPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_AfterQA_MidPt_ProjPt));
            delete fHistPion_dEdx_after_AfterQA_MidPt_ProjPt;
            fHistPion_dEdx_after_AfterQA_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdx_after_MidPt %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_AfterQA_HighPt = new TH2D (*fHistPion_dEdx_after_AfterQA);
        if(fHistPion_dEdx_after_AfterQA_HighPt){
            fHistPion_dEdx_after_AfterQA_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA_HighPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA_HighPt,fHistPion_dEdx_after_AfterQA_HighPt->GetXaxis()->FindBin(3),fHistPion_dEdx_after_AfterQA_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_dEdx_after_AfterQA_HighPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_AfterQA_HighPt,minYB-1,maxYB+1);
            //cout<<"minYB: "<<minYB<<"; maxYB: "<<maxYB<<endl;
            //SetYRange(fHistPion_dEdx_after_AfterQA_HighPt,-2,2);
            SetZMinMaxTH2(fHistPion_dEdx_after_AfterQA_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_AfterQA_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_AfterQA_HighPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA_HighPt.push_back(new TH2D(*fHistPion_dEdx_after_AfterQA_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_AfterQA_HighPt_ProjPt= (TH1D*)fHistPion_dEdx_after_AfterQA_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_AfterQA_HighPt->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_AfterQA_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_AfterQA_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_AfterQA_HighPt_ProjPt,"",fHistPion_dEdx_after_AfterQA_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_AfterQA_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_AfterQA_HighPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_AfterQA_HighPt_ProjPt));
            delete fHistPion_dEdx_after_AfterQA_HighPt_ProjPt;
            fHistPion_dEdx_after_AfterQA_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdx_after_HighPt %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_AfterQA = (TH2D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHistPion_dEdxSignal_after_AfterQA){
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_dEdxSignal_after_AfterQA,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_AfterQA,fHistPion_dEdxSignal_after_AfterQA->FindBin(0.),fHistPion_dEdxSignal_after_AfterQA->FindBin(200.));
            fHistPion_dEdxSignal_after_AfterQA->SetAxisRange(0.,200,"Y");
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_AfterQA,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_AfterQA,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA.push_back(new TH2D(*fHistPion_dEdxSignal_after_AfterQA));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_AfterQA_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_AfterQA->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_AfterQA_ProjPt,"",fHistPion_dEdxSignal_after_AfterQA->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_AfterQA_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_AfterQA_ProjPt));
            delete fHistPion_dEdxSignal_after_AfterQA_ProjPt;
            fHistPion_dEdxSignal_after_AfterQA_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdxSignal_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_AfterQA_LowPt = new TH2D (*fHistPion_dEdxSignal_after_AfterQA);
        if(fHistPion_dEdxSignal_after_AfterQA_LowPt){
            fHistPion_dEdxSignal_after_AfterQA_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_LowPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA_LowPt,0,fHistPion_dEdxSignal_after_AfterQA_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_dEdxSignal_after_AfterQA_LowPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_AfterQA_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_AfterQA_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_AfterQA_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_AfterQA_LowPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA_LowPt.push_back(new TH2D(*fHistPion_dEdxSignal_after_AfterQA_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_AfterQA_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_AfterQA_LowPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt,"",fHistPion_dEdxSignal_after_AfterQA_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt));
            delete fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt;
            fHistPion_dEdxSignal_after_AfterQA_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdxSignal_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_AfterQA_MidPt = new TH2D (*fHistPion_dEdxSignal_after_AfterQA);
        if(fHistPion_dEdxSignal_after_AfterQA_MidPt){
            fHistPion_dEdxSignal_after_AfterQA_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_MidPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA_MidPt,fHistPion_dEdxSignal_after_AfterQA_MidPt->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_after_AfterQA_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_dEdxSignal_after_AfterQA_MidPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_AfterQA_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_AfterQA_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_AfterQA_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_AfterQA_MidPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA_MidPt.push_back(new TH2D(*fHistPion_dEdxSignal_after_AfterQA_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_AfterQA_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_AfterQA_MidPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt,"",fHistPion_dEdxSignal_after_AfterQA_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt));
            delete fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt;
            fHistPion_dEdxSignal_after_AfterQA_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdxSignal_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_AfterQA_HighPt = new TH2D (*fHistPion_dEdxSignal_after_AfterQA);
        if(fHistPion_dEdxSignal_after_AfterQA_HighPt){
            fHistPion_dEdxSignal_after_AfterQA_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_HighPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA_HighPt,fHistPion_dEdxSignal_after_AfterQA_HighPt->GetXaxis()->FindBin(3),fHistPion_dEdxSignal_after_AfterQA_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_dEdxSignal_after_AfterQA_HighPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_AfterQA_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_AfterQA_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_AfterQA_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_AfterQA_HighPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA_HighPt.push_back(new TH2D(*fHistPion_dEdxSignal_after_AfterQA_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_AfterQA_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_AfterQA_HighPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt,"",fHistPion_dEdxSignal_after_AfterQA_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt));
            delete fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt;
            fHistPion_dEdxSignal_after_AfterQA_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdxSignal_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_TOF_after
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_AfterQA = (TH2D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHistPion_TOF_after_AfterQA){
            GetMinMaxBin(fHistPion_TOF_after_AfterQA,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_TOF_after_AfterQA,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_AfterQA,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_AfterQA,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_AfterQA,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA.push_back(new TH2D(*fHistPion_TOF_after_AfterQA));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_AfterQA_ProjPt= (TH1D*)fHistPion_TOF_after_AfterQA->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_AfterQA_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_AfterQA_ProjPt,"",fHistPion_TOF_after_AfterQA->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_AfterQA_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_AfterQA_ProjPt));
            delete fHistPion_TOF_after_AfterQA_ProjPt;
            fHistPion_TOF_after_AfterQA_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_TOF_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_TOF_after_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_AfterQA_LowPt = new TH2D (*fHistPion_TOF_after_AfterQA);
        if(fHistPion_TOF_after_AfterQA_LowPt){
            fHistPion_TOF_after_AfterQA_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_TOF_after_AfterQA_LowPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA_LowPt,0,fHistPion_TOF_after_AfterQA_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_TOF_after_AfterQA_LowPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_AfterQA_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_AfterQA_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_AfterQA_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_AfterQA_LowPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA_LowPt.push_back(new TH2D(*fHistPion_TOF_after_AfterQA_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_AfterQA_LowPt_ProjPt= (TH1D*)fHistPion_TOF_after_AfterQA_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_AfterQA_LowPt->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_AfterQA_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_AfterQA_LowPt_ProjPt,"",fHistPion_TOF_after_AfterQA_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_AfterQA_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA_LowPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_AfterQA_LowPt_ProjPt));
            delete fHistPion_TOF_after_AfterQA_LowPt_ProjPt;
            fHistPion_TOF_after_AfterQA_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_TOF_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_AfterQA_MidPt = new TH2D (*fHistPion_TOF_after_AfterQA);
        if(fHistPion_TOF_after_AfterQA_MidPt){
            fHistPion_TOF_after_AfterQA_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_TOF_after_AfterQA_MidPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA_MidPt,fHistPion_TOF_after_AfterQA_MidPt->GetXaxis()->FindBin(2),fHistPion_TOF_after_AfterQA_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_TOF_after_AfterQA_MidPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_AfterQA_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_AfterQA_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_AfterQA_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_AfterQA_MidPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA_MidPt.push_back(new TH2D(*fHistPion_TOF_after_AfterQA_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_AfterQA_MidPt_ProjPt= (TH1D*)fHistPion_TOF_after_AfterQA_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_AfterQA_MidPt->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_AfterQA_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_AfterQA_MidPt_ProjPt,"",fHistPion_TOF_after_AfterQA_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_AfterQA_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA_MidPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_AfterQA_MidPt_ProjPt));
            delete fHistPion_TOF_after_AfterQA_MidPt_ProjPt;
            fHistPion_TOF_after_AfterQA_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdxSignal_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_AfterQA_HighPt = new TH2D (*fHistPion_TOF_after_AfterQA);
        if(fHistPion_TOF_after_AfterQA_HighPt){
            fHistPion_TOF_after_AfterQA_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_TOF_after_AfterQA_HighPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA_HighPt,fHistPion_TOF_after_AfterQA_HighPt->GetXaxis()->FindBin(3),fHistPion_TOF_after_AfterQA_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_TOF_after_AfterQA_HighPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_AfterQA_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_AfterQA_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_AfterQA_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_AfterQA_HighPt, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA_HighPt.push_back(new TH2D(*fHistPion_TOF_after_AfterQA_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_AfterQA_HighPt_ProjPt= (TH1D*)fHistPion_TOF_after_AfterQA_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_AfterQA_HighPt->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_AfterQA_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_AfterQA_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_AfterQA_HighPt_ProjPt,"",fHistPion_TOF_after_AfterQA_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_AfterQA_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_AfterQA_HighPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_AfterQA_HighPt_ProjPt));
            delete fHistPion_TOF_after_AfterQA_HighPt_ProjPt;
            fHistPion_TOF_after_AfterQA_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: Pion_dEdxSignal_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: hTrack_DCAxy_Pt_after
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAxy_Pt_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_DCAxy_Pt_after_AfterQA = (TH2D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHisthTrack_DCAxy_Pt_after_AfterQA){
            GetMinMaxBin(fHisthTrack_DCAxy_Pt_after_AfterQA,minB,maxB);
            SetXRange(fHisthTrack_DCAxy_Pt_after_AfterQA,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_DCAxy_Pt_after_AfterQA,minYB,maxYB);
            SetYRange(fHisthTrack_DCAxy_Pt_after_AfterQA,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_DCAxy_Pt_after_AfterQA,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_DCAxy_Pt_after_AfterQA,"",
                                 "DCA_{#it{xy}} (cm)","#it{p}_{T} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_DCAxy_Pt_after_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAxy_Pt_after_AfterQA.push_back(new TH2D(*fHisthTrack_DCAxy_Pt_after_AfterQA));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt= (TH1D*)fHisthTrack_DCAxy_Pt_after_AfterQA->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_DCAxy_Pt_after_AfterQA->GetXaxis()->GetFirst(),fHisthTrack_DCAxy_Pt_after_AfterQA->GetYaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt,"",fHisthTrack_DCAxy_Pt_after_AfterQA->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAxy_Pt_after_AfterQA_ProjPt.push_back(new TH1D(*fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt));
            delete fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt;
            fHisthTrack_DCAxy_Pt_after_AfterQA_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: hTrack_DCAxy_Pt_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: hTrack_DCAz_Pt_after
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAz_Pt_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_DCAz_Pt_after_AfterQA = (TH2D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHisthTrack_DCAz_Pt_after_AfterQA){
            GetMinMaxBin(fHisthTrack_DCAz_Pt_after_AfterQA,minB,maxB);
            SetXRange(fHisthTrack_DCAz_Pt_after_AfterQA,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_DCAz_Pt_after_AfterQA,minYB,maxYB);
            SetYRange(fHisthTrack_DCAz_Pt_after_AfterQA,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_DCAz_Pt_after_AfterQA,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_DCAz_Pt_after_AfterQA,"",
                                 "DCA_{#it{z}} (cm)","#it{p}_{T} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_DCAz_Pt_after_AfterQA, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAz_Pt_after_AfterQA.push_back(new TH2D(*fHisthTrack_DCAz_Pt_after_AfterQA));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt= (TH1D*)fHisthTrack_DCAz_Pt_after_AfterQA->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_DCAz_Pt_after_AfterQA->GetXaxis()->GetFirst(),fHisthTrack_DCAz_Pt_after_AfterQA->GetYaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt,"",fHisthTrack_DCAz_Pt_after_AfterQA->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAz_Pt_after_AfterQA_ProjPt.push_back(new TH1D(*fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt));
            delete fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt;
            fHisthTrack_DCAz_Pt_after_AfterQA_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: hTrack_DCAz_Pt_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: hTrack_NFindCls_Pt_TPC_after
        if (iParticleType==0){StrNameOfHistogram="hTrack_NFindCls_Pt_TPC_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_NFindCls_Pt_TPC_after_AfterQA = (TH2D*)PionCutsContainer->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCutsContainerCutString.Data()));
        if(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA){
            TH2D* fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis = SwitchTF2DAxis(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA);
            GetMinMaxBin(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis,minB,maxB);
            SetXRange(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis,minYB,maxYB);
            SetYRange(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis,"","#it{p}_{T} (GeV/#it{c}) TPC",
                                 "Findable Clusters after Cut",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis, Form("%s/%s_AfterQA_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_NFindCls_Pt_TPC_after_AfterQA.push_back(new TH2D(*fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt;
            fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt= (TH1D*)fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis->GetYaxis()->GetFirst(),fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis->GetXaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt,"",fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.push_back(new TH1D(*fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt));
            delete fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis;
            fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis=NULL;
            delete fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt;
            fHisthTrack_NFindCls_Pt_TPC_after_AfterQA_switchedAxis_ProjPt=NULL;
        } else cout << Form("INFO: Object |AfterQA: hTrack_NFindCls_Pt_TPC_after %s| could not be found! Skipping Draw...",fPionCutsContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //-----------------------------|Get Histograms: Pre Selection (Pion Cut folder in Main Directory)|----------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Drawing Pion Cut Histograms Pre Selection" << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // PionCuts2Container: IsPionSelected
        if (iParticleType==0){StrNameOfHistogram="IsPionSelected";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistIsPionSelected_PreSel = (TH1D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistIsPionSelected_PreSel){
            GetMinMaxBin(fHistIsPionSelected_PreSel,minB,maxB);
            SetXRange(fHistIsPionSelected_PreSel,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistIsPionSelected_PreSel,"","","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistIsPionSelected_PreSel);
            SaveCanvasAndWriteHistogram(canvas, fHistIsPionSelected_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecIsPionSelected_PreSel.push_back(new TH1D(*fHistIsPionSelected_PreSel));
        } else cout << "INFO: Object |Pre Selection: fHistIsPionSelected_PreSel| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        // PionCuts2Container: dEdxCuts
        if (iParticleType==0){StrNameOfHistogram="dEdxCuts";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH1D* fHistdEdxCuts_PreSel = (TH1D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistdEdxCuts_PreSel){
            GetMinMaxBin(fHistdEdxCuts_PreSel,minB,maxB);
            SetXRange(fHistdEdxCuts_PreSel,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                 fHistdEdxCuts_PreSel,"","","# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            //WriteHistogram(fHistdEdxCuts_PreSel);
            SaveCanvasAndWriteHistogram(canvas, fHistdEdxCuts_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecdEdxCuts_PreSel.push_back(new TH1D(*fHistdEdxCuts_PreSel));
        } else cout << "INFO: Object |Pre Selection: fHistdEdxCuts_PreSel| could not be found! Skipping Draw..." << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //--------------------------|Get Histograms: Pre Selection (Pion Cut folder in Main Directory) before|-------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_before
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_before_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_ITS_before_PreSel){
            GetMinMaxBin(fHistPion_ITS_before_PreSel,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_ITS_before_PreSel,minYB,maxYB);
            SetYRange(fHistPion_ITS_before_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_before_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_before_PreSel,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_before_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel.push_back(new TH2D(*fHistPion_ITS_before_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_before_PreSel_ProjPt= (TH1D*)fHistPion_ITS_before_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_before_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_before_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_before_PreSel_ProjPt,"",fHistPion_ITS_before_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_before_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel_ProjPt.push_back(new TH1D(*fHistPion_ITS_before_PreSel_ProjPt));
            delete fHistPion_ITS_before_PreSel_ProjPt;
            fHistPion_ITS_before_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_before_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_before_PreSel_LowPt = new TH2D (*fHistPion_ITS_before_PreSel);
        if(fHistPion_ITS_before_PreSel_LowPt){
            fHistPion_ITS_before_PreSel_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_ITS_before_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel_LowPt,0,fHistPion_ITS_before_PreSel_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_ITS_before_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_before_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_before_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_before_PreSel_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_before_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel.push_back(new TH2D(*fHistPion_ITS_before_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_before_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_ITS_before_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_before_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_before_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_before_PreSel_LowPt_ProjPt,"",fHistPion_ITS_before_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_before_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_before_PreSel_LowPt_ProjPt));
            delete fHistPion_ITS_before_PreSel_LowPt_ProjPt;
            fHistPion_ITS_before_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_before_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_before_PreSel_MidPt = new TH2D (*fHistPion_ITS_before_PreSel);
        if(fHistPion_ITS_before_PreSel_MidPt){
            fHistPion_ITS_before_PreSel_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_ITS_before_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel_MidPt,fHistPion_ITS_before_PreSel_MidPt->GetXaxis()->FindBin(2),fHistPion_ITS_before_PreSel_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_ITS_before_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_before_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_before_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_before_PreSel_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_before_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel.push_back(new TH2D(*fHistPion_ITS_before_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_before_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_ITS_before_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_before_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_before_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_before_PreSel_MidPt_ProjPt,"",fHistPion_ITS_before_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_before_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_before_PreSel_MidPt_ProjPt));
            delete fHistPion_ITS_before_PreSel_MidPt_ProjPt;
            fHistPion_ITS_before_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_before_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_before_PreSel_HighPt = new TH2D (*fHistPion_ITS_before_PreSel);
        if(fHistPion_ITS_before_PreSel_HighPt){
            fHistPion_ITS_before_PreSel_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_ITS_before_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel_HighPt,fHistPion_ITS_before_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_ITS_before_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_ITS_before_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_before_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_before_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_before_PreSel_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_before_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel.push_back(new TH2D(*fHistPion_ITS_before_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_before_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_ITS_before_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_before_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_before_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_before_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_before_PreSel_HighPt_ProjPt,"",fHistPion_ITS_before_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_before_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_before_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_before_PreSel_HighPt_ProjPt));
            delete fHistPion_ITS_before_PreSel_HighPt_ProjPt;
            fHistPion_ITS_before_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_before_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_dEdx_before_PreSel){
            GetMinMaxBin(fHistPion_dEdx_before_PreSel,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_dEdx_before_PreSel,minYB,maxYB);
            SetYRange(fHistPion_dEdx_before_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_before_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_before_PreSel,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_before_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel.push_back(new TH2D(*fHistPion_dEdx_before_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_before_PreSel_ProjPt= (TH1D*)fHistPion_dEdx_before_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_before_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_before_PreSel_ProjPt,"",fHistPion_dEdx_before_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_before_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel_ProjPt.push_back(new TH1D(*fHistPion_dEdx_before_PreSel_ProjPt));
            delete fHistPion_dEdx_before_PreSel_ProjPt;
            fHistPion_dEdx_before_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_before_PreSel_LowPt = new TH2D (*fHistPion_dEdx_before_PreSel);
        if(fHistPion_dEdx_before_PreSel_LowPt){
            fHistPion_dEdx_before_PreSel_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdx_before_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel_LowPt,0,fHistPion_dEdx_before_PreSel_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_dEdx_before_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_before_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_before_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_before_PreSel_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_before_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel_LowPt.push_back(new TH2D(*fHistPion_dEdx_before_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_before_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_dEdx_before_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_before_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_before_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_before_PreSel_LowPt_ProjPt,"",fHistPion_dEdx_before_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_before_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_before_PreSel_LowPt_ProjPt));
            delete fHistPion_dEdx_before_PreSel_LowPt_ProjPt;
            fHistPion_dEdx_before_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_before_PreSel_MidPt = new TH2D (*fHistPion_dEdx_before_PreSel);
        if(fHistPion_dEdx_before_PreSel_MidPt){
            fHistPion_dEdx_before_PreSel_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdx_before_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel_MidPt,fHistPion_dEdx_before_PreSel_MidPt->GetXaxis()->FindBin(2),fHistPion_dEdx_before_PreSel_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_dEdx_before_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_before_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_before_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_before_PreSel_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_before_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel_MidPt.push_back(new TH2D(*fHistPion_dEdx_before_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_before_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_dEdx_before_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_before_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_before_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_before_PreSel_MidPt_ProjPt,"",fHistPion_dEdx_before_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_before_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_before_PreSel_MidPt_ProjPt));
            delete fHistPion_dEdx_before_PreSel_MidPt_ProjPt;
            fHistPion_dEdx_before_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_before_PreSel_HighPt = new TH2D (*fHistPion_dEdx_before_PreSel);
        if(fHistPion_dEdx_before_PreSel_HighPt){
            fHistPion_dEdx_before_PreSel_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdx_before_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel_HighPt,fHistPion_dEdx_before_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_dEdx_before_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_dEdx_before_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_before_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_before_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_before_PreSel_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_before_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel_HighPt.push_back(new TH2D(*fHistPion_dEdx_before_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_before_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_dEdx_before_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_before_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_before_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_before_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_before_PreSel_HighPt_ProjPt,"",fHistPion_dEdx_before_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_before_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_before_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_before_PreSel_HighPt_ProjPt));
            delete fHistPion_dEdx_before_PreSel_HighPt_ProjPt;
            fHistPion_dEdx_before_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_before_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_dEdxSignal_before_PreSel){
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_before_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_before_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_before_PreSel,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_before_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel.push_back(new TH2D(*fHistPion_dEdxSignal_before_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_before_PreSel_ProjPt= (TH1D*)fHistPion_dEdxSignal_before_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_before_PreSel_ProjPt,"",fHistPion_dEdxSignal_before_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_before_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_before_PreSel_ProjPt));
            delete fHistPion_dEdxSignal_before_PreSel_ProjPt;
            fHistPion_dEdxSignal_before_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_before_PreSel_LowPt = new TH2D (*fHistPion_dEdxSignal_before_PreSel);
        if(fHistPion_dEdxSignal_before_PreSel_LowPt){
            fHistPion_dEdxSignal_before_PreSel_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel_LowPt,0,fHistPion_dEdxSignal_before_PreSel_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_before_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_before_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_before_PreSel_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_before_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel_LowPt.push_back(new TH2D(*fHistPion_dEdxSignal_before_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_before_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_before_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt,"",fHistPion_dEdxSignal_before_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt));
            delete fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt;
            fHistPion_dEdxSignal_before_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_before_PreSel_MidPt = new TH2D (*fHistPion_dEdxSignal_before_PreSel);
        if(fHistPion_dEdxSignal_before_PreSel_MidPt){
            fHistPion_dEdxSignal_before_PreSel_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel_MidPt,fHistPion_dEdxSignal_before_PreSel_MidPt->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_before_PreSel_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_before_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_before_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_before_PreSel_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_before_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel_MidPt.push_back(new TH2D(*fHistPion_dEdxSignal_before_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_before_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_before_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt,"",fHistPion_dEdxSignal_before_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt));
            delete fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt;
            fHistPion_dEdxSignal_before_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_before_PreSel_HighPt = new TH2D (*fHistPion_dEdxSignal_before_PreSel);
        if(fHistPion_dEdxSignal_before_PreSel_HighPt){
            fHistPion_dEdxSignal_before_PreSel_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel_HighPt,fHistPion_dEdxSignal_before_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_dEdxSignal_before_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_before_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_before_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_before_PreSel_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_before_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel_HighPt.push_back(new TH2D(*fHistPion_dEdxSignal_before_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_before_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_before_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt,"",fHistPion_dEdxSignal_before_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt));
            delete fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt;
            fHistPion_dEdxSignal_before_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_before_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_TOF_before_PreSel){
            GetMinMaxBin(fHistPion_TOF_before_PreSel,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_TOF_before_PreSel,minYB,maxYB);
            SetYRange(fHistPion_TOF_before_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_before_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_before_PreSel,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_before_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel.push_back(new TH2D(*fHistPion_TOF_before_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_before_PreSel_ProjPt= (TH1D*)fHistPion_TOF_before_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_before_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_before_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_before_PreSel_ProjPt,"",fHistPion_TOF_before_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_before_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel_ProjPt.push_back(new TH1D(*fHistPion_TOF_before_PreSel_ProjPt));
            delete fHistPion_TOF_before_PreSel_ProjPt;
            fHistPion_TOF_before_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_before_PreSel_LowPt = new TH2D (*fHistPion_TOF_before_PreSel);
        if(fHistPion_TOF_before_PreSel_LowPt){
            fHistPion_TOF_before_PreSel_LowPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_TOF_before_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel_LowPt,0,fHistPion_TOF_before_PreSel_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_TOF_before_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_before_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_before_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_before_PreSel_LowPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_before_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel_LowPt.push_back(new TH2D(*fHistPion_TOF_before_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_before_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_TOF_before_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_before_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_before_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_before_PreSel_LowPt_ProjPt,"",fHistPion_TOF_before_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_before_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_before_PreSel_LowPt_ProjPt));
            delete fHistPion_TOF_before_PreSel_LowPt_ProjPt;
            fHistPion_TOF_before_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_before_PreSel_MidPt = new TH2D (*fHistPion_TOF_before_PreSel);
        if(fHistPion_TOF_before_PreSel_MidPt){
            fHistPion_TOF_before_PreSel_MidPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_TOF_before_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel_MidPt,fHistPion_TOF_before_PreSel_MidPt->GetXaxis()->FindBin(2),fHistPion_TOF_before_PreSel_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_TOF_before_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_before_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_before_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_before_PreSel_MidPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_before_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel_MidPt.push_back(new TH2D(*fHistPion_TOF_before_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_before_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_TOF_before_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_before_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_before_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_before_PreSel_MidPt_ProjPt,"",fHistPion_TOF_before_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_before_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_before_PreSel_MidPt_ProjPt));
            delete fHistPion_TOF_before_PreSel_MidPt_ProjPt;
            fHistPion_TOF_before_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_before_PreSel_HighPt = new TH2D (*fHistPion_TOF_before_PreSel);
        if(fHistPion_TOF_before_PreSel_HighPt){
            fHistPion_TOF_before_PreSel_HighPt->SetName(StrNameOfHistogram);
            GetMinMaxBin(fHistPion_TOF_before_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel_HighPt,fHistPion_TOF_before_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_TOF_before_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_TOF_before_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_before_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_before_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_before_PreSel_HighPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_before_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel_HighPt.push_back(new TH2D(*fHistPion_TOF_before_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_before_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_TOF_before_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_before_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_before_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_before_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_before_PreSel_HighPt_ProjPt,"",fHistPion_TOF_before_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_before_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_before_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_before_PreSel_HighPt_ProjPt));
            delete fHistPion_TOF_before_PreSel_HighPt_ProjPt;
            fHistPion_TOF_before_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAxy_Pt_before
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAxy_Pt_before";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_DCAxy_Pt_before_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHisthTrack_DCAxy_Pt_before_PreSel){
            GetMinMaxBin(fHisthTrack_DCAxy_Pt_before_PreSel,minB,maxB);
            SetXRange(fHisthTrack_DCAxy_Pt_before_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_DCAxy_Pt_before_PreSel,minYB,maxYB);
            SetYRange(fHisthTrack_DCAxy_Pt_before_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_DCAxy_Pt_before_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_DCAxy_Pt_before_PreSel,"",
                                 "DCA_{#it{xy}} (cm)","#it{p}_{T} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_DCAxy_Pt_before_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAxy_Pt_before_PreSel.push_back(new TH2D(*fHisthTrack_DCAxy_Pt_before_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt= (TH1D*)fHisthTrack_DCAxy_Pt_before_PreSel->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_DCAxy_Pt_before_PreSel->GetXaxis()->GetFirst(),fHisthTrack_DCAxy_Pt_before_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt,"",fHisthTrack_DCAxy_Pt_before_PreSel->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAxy_Pt_before_PreSel_ProjPt.push_back(new TH1D(*fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt));
            delete fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt;
            fHisthTrack_DCAxy_Pt_before_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: hTrack_DCAxy_Pt_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAz_Pt_before
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAz_Pt_before";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_DCAz_Pt_before_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHisthTrack_DCAz_Pt_before_PreSel){
            GetMinMaxBin(fHisthTrack_DCAz_Pt_before_PreSel,minB,maxB);
            SetXRange(fHisthTrack_DCAz_Pt_before_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_DCAz_Pt_before_PreSel,minYB,maxYB);
            SetYRange(fHisthTrack_DCAz_Pt_before_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_DCAz_Pt_before_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_DCAz_Pt_before_PreSel,"",
                                 "DCA_{#it{z}} (cm)","#it{p}_{T} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_DCAz_Pt_before_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAz_Pt_before_PreSel.push_back(new TH2D(*fHisthTrack_DCAz_Pt_before_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_DCAz_Pt_before_PreSel_ProjPt= (TH1D*)fHisthTrack_DCAz_Pt_before_PreSel->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_DCAz_Pt_before_PreSel->GetXaxis()->GetFirst(),fHisthTrack_DCAz_Pt_before_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_DCAz_Pt_before_PreSel_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_DCAz_Pt_before_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_DCAz_Pt_before_PreSel_ProjPt,"",fHisthTrack_DCAz_Pt_before_PreSel->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_DCAz_Pt_before_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAz_Pt_before_PreSel_ProjPt.push_back(new TH1D(*fHisthTrack_DCAz_Pt_before_PreSel_ProjPt));
            delete fHisthTrack_DCAz_Pt_before_PreSel_ProjPt;
            fHisthTrack_DCAz_Pt_before_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: hTrack_DCAz_Pt_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_NFindCls_Pt_TPC_before
        if (iParticleType==0){StrNameOfHistogram="hTrack_NFindCls_Pt_TPC_before";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_NFindCls_Pt_TPC_before_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHisthTrack_NFindCls_Pt_TPC_before_PreSel){
            TH2D* fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis=SwitchTF2DAxis(fHisthTrack_NFindCls_Pt_TPC_before_PreSel);
            GetMinMaxBin(fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis,minB,maxB);
            SetXRange(fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis,minYB,maxYB);
            SetYRange(fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis,"","#it{p}_{T} (GeV/#it{c}) TPC",
                                 "Findable Clusters before Cut",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_NFindCls_Pt_TPC_before_PreSel.push_back(new TH2D(*fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt= (TH1D*)fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis->GetXaxis()->GetFirst(),fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis->GetXaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt,"",fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.push_back(new TH1D(*fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt));
            delete fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis;
            fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis=NULL;
            delete fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt;
            fHisthTrack_NFindCls_Pt_TPC_before_PreSel_switchedAxis_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: hTrack_NFindCls_Pt_TPC_before %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //--------------------------|Get Histograms: Pre Selection (Pion Cut folder in Main Directory) after|-------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_after
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_ITS_after_PreSel){
            GetMinMaxBin(fHistPion_ITS_after_PreSel,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_ITS_after_PreSel,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_PreSel,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel.push_back(new TH2D(*fHistPion_ITS_after_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_PreSel_ProjPt= (TH1D*)fHistPion_ITS_after_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_after_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_PreSel_ProjPt,"",fHistPion_ITS_after_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_PreSel_ProjPt));
            delete fHistPion_ITS_after_PreSel_ProjPt;
            fHistPion_ITS_after_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_after %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_after_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_PreSel_LowPt = new TH2D (*fHistPion_ITS_after_PreSel);
        if(fHistPion_ITS_after_PreSel_LowPt){
            GetMinMaxBin(fHistPion_ITS_after_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel_LowPt,0,fHistPion_ITS_after_PreSel_LowPt->GetXaxis()->FindBin(0.2));
            GetMinMaxBinY(fHistPion_ITS_after_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_PreSel_LowPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel_LowPt.push_back(new TH2D(*fHistPion_ITS_after_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_ITS_after_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_ITS_after_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_PreSel_LowPt_ProjPt,"",fHistPion_ITS_after_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_PreSel_LowPt_ProjPt));
            delete fHistPion_ITS_after_PreSel_LowPt_ProjPt;
            fHistPion_ITS_after_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_after_LowPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_after_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_PreSel_MidPt = new TH2D (*fHistPion_ITS_after_PreSel);
        if(fHistPion_ITS_after_PreSel_MidPt){
            GetMinMaxBin(fHistPion_ITS_after_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel_MidPt,fHistPion_ITS_after_PreSel_MidPt->GetXaxis()->FindBin(2),fHistPion_ITS_after_PreSel_MidPt->GetXaxis()->FindBin(3));
            GetMinMaxBinY(fHistPion_ITS_after_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_PreSel_MidPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel_MidPt.push_back(new TH2D(*fHistPion_ITS_after_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_ITS_after_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_ITS_after_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_PreSel_MidPt_ProjPt,"",fHistPion_ITS_after_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_PreSel_MidPt_ProjPt));
            delete fHistPion_ITS_after_PreSel_MidPt_ProjPt;
            fHistPion_ITS_after_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_after_MidPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_after_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_ITS_after_PreSel_HighPt = new TH2D (*fHistPion_ITS_after_PreSel);
        if(fHistPion_ITS_after_PreSel_HighPt){
            GetMinMaxBin(fHistPion_ITS_after_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel_HighPt,fHistPion_ITS_after_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_ITS_after_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_ITS_after_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_ITS_after_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_ITS_after_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_ITS_after_PreSel_HighPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_ITS_after_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel_HighPt.push_back(new TH2D(*fHistPion_ITS_after_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_ITS_after_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_ITS_after_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_ITS_after_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_ITS_after_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_ITS_after_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_ITS_after_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_ITS_after_PreSel_HighPt_ProjPt,"",fHistPion_ITS_after_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_ITS_after_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_ITS_after_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_ITS_after_PreSel_HighPt_ProjPt));
            delete fHistPion_ITS_after_PreSel_HighPt_ProjPt;
            fHistPion_ITS_after_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_ITS_after_HighPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_dEdx_after_PreSel){
            GetMinMaxBin(fHistPion_dEdx_after_PreSel,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_dEdx_after_PreSel,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_after_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_PreSel,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel.push_back(new TH2D(*fHistPion_dEdx_after_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_PreSel_ProjPt= (TH1D*)fHistPion_dEdx_after_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_PreSel_ProjPt,"",fHistPion_dEdx_after_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_PreSel_ProjPt));
            delete fHistPion_dEdx_after_PreSel_ProjPt;
            fHistPion_dEdx_after_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_after %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_PreSel_LowPt = new TH2D (*fHistPion_dEdx_after_PreSel);
        if(fHistPion_dEdx_after_PreSel_LowPt){
            GetMinMaxBin(fHistPion_dEdx_after_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel_LowPt,0,fHistPion_dEdx_after_PreSel_LowPt->FindBin(3));
            GetMinMaxBinY(fHistPion_dEdx_after_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_after_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_PreSel_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel_LowPt.push_back(new TH2D(*fHistPion_dEdx_after_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_dEdx_after_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_PreSel_LowPt_ProjPt,"",fHistPion_dEdx_after_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_PreSel_LowPt_ProjPt));
            delete fHistPion_dEdx_after_PreSel_LowPt_ProjPt;
            fHistPion_dEdx_after_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_after_LowPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_PreSel_MidPt = new TH2D (*fHistPion_dEdx_after_PreSel);
        if(fHistPion_dEdx_after_PreSel_MidPt){
            GetMinMaxBin(fHistPion_dEdx_after_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel_MidPt,fHistPion_dEdx_after_PreSel_MidPt->GetXaxis()->FindBin(3),fHistPion_dEdx_after_PreSel_MidPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_dEdx_after_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_after_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_PreSel_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel_MidPt.push_back(new TH2D(*fHistPion_dEdx_after_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_dEdx_after_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_PreSel_MidPt_ProjPt,"",fHistPion_dEdx_after_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_PreSel_MidPt_ProjPt));
            delete fHistPion_dEdx_after_PreSel_MidPt_ProjPt;
            fHistPion_dEdx_after_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_after_MidPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdx_after_PreSel_HighPt = new TH2D (*fHistPion_dEdx_after_PreSel);
        if(fHistPion_dEdx_after_PreSel_HighPt){
            GetMinMaxBin(fHistPion_dEdx_after_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel_HighPt,fHistPion_dEdx_after_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_dEdx_after_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_dEdx_after_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_dEdx_after_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdx_after_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdx_after_PreSel_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdx_after_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel_HighPt.push_back(new TH2D(*fHistPion_dEdx_after_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdx_after_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_dEdx_after_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdx_after_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdx_after_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdx_after_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdx_after_PreSel_HighPt_ProjPt,"",fHistPion_dEdx_after_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdx_after_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdx_after_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_dEdx_after_PreSel_HighPt_ProjPt));
            delete fHistPion_dEdx_after_PreSel_HighPt_ProjPt;
            fHistPion_dEdx_after_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdx_after_HighPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_dEdxSignal_after_PreSel){
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_dEdxSignal_after_PreSel,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_PreSel,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel.push_back(new TH2D(*fHistPion_dEdxSignal_after_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_PreSel_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_PreSel_ProjPt,"",fHistPion_dEdxSignal_after_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_PreSel_ProjPt));
            delete fHistPion_dEdxSignal_after_PreSel_ProjPt;
            fHistPion_dEdxSignal_after_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_after %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_PreSel_LowPt = new TH2D (*fHistPion_dEdxSignal_after_PreSel);
        if(fHistPion_dEdxSignal_after_PreSel_LowPt){
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel_LowPt,0,fHistPion_dEdxSignal_after_PreSel_LowPt->FindBin(0.2));
            GetMinMaxBinY(fHistPion_dEdxSignal_after_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_PreSel_LowPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel_LowPt.push_back(new TH2D(*fHistPion_dEdxSignal_after_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt,"",fHistPion_dEdxSignal_after_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt));
            delete fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt;
            fHistPion_dEdxSignal_after_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_after_LowPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_PreSel_MidPt = new TH2D (*fHistPion_dEdxSignal_after_PreSel);
        if(fHistPion_dEdxSignal_after_PreSel_MidPt){
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel_MidPt,fHistPion_dEdxSignal_after_PreSel_MidPt->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_after_PreSel_MidPt->FindBin(3));
            GetMinMaxBinY(fHistPion_dEdxSignal_after_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_PreSel_MidPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel_MidPt.push_back(new TH2D(*fHistPion_dEdxSignal_after_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt,"",fHistPion_dEdxSignal_after_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt));
            delete fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt;
            fHistPion_dEdxSignal_after_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_after_MidPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_dEdxSignal_after_PreSel_HighPt = new TH2D (*fHistPion_dEdxSignal_after_PreSel);
        if(fHistPion_dEdxSignal_after_PreSel_HighPt){
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel_HighPt,fHistPion_dEdxSignal_after_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_dEdxSignal_after_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_dEdxSignal_after_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_dEdxSignal_after_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_dEdxSignal_after_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_dEdxSignal_after_PreSel_HighPt,"",
                                 "#it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_dEdxSignal_after_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel_HighPt.push_back(new TH2D(*fHistPion_dEdxSignal_after_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_dEdxSignal_after_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_dEdxSignal_after_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt,"",fHistPion_dEdxSignal_after_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt));
            delete fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt;
            fHistPion_dEdxSignal_after_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_dEdxSignal_after_HighPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHistPion_TOF_after_PreSel){
            GetMinMaxBin(fHistPion_TOF_after_PreSel,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHistPion_TOF_after_PreSel,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_PreSel,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel.push_back(new TH2D(*fHistPion_TOF_after_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_PreSel_ProjPt= (TH1D*)fHistPion_TOF_after_PreSel->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_after_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_PreSel_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_PreSel_ProjPt,"",fHistPion_TOF_after_PreSel->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_PreSel_ProjPt));
            delete fHistPion_TOF_after_PreSel_ProjPt;
            fHistPion_TOF_after_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_after %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel_LowPt
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_LowPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_PreSel_LowPt = new TH2D (*fHistPion_TOF_after_PreSel);
        if(fHistPion_TOF_after_PreSel_LowPt){
            GetMinMaxBin(fHistPion_TOF_after_PreSel_LowPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel_LowPt,0,fHistPion_TOF_after_PreSel_LowPt->FindBin(0.2));
            GetMinMaxBinY(fHistPion_TOF_after_PreSel_LowPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_PreSel_LowPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_PreSel_LowPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_PreSel_LowPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_PreSel_LowPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel_LowPt.push_back(new TH2D(*fHistPion_TOF_after_PreSel_LowPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_PreSel_LowPt_ProjPt= (TH1D*)fHistPion_TOF_after_PreSel_LowPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_PreSel_LowPt->GetXaxis()->GetFirst(),fHistPion_TOF_after_PreSel_LowPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_PreSel_LowPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel_LowPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_PreSel_LowPt_ProjPt,"",fHistPion_TOF_after_PreSel_LowPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_PreSel_LowPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel_LowPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_PreSel_LowPt_ProjPt));
            delete fHistPion_TOF_after_PreSel_LowPt_ProjPt;
            fHistPion_TOF_after_PreSel_LowPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_after_LowPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel_MidPt
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_MidPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_PreSel_MidPt = new TH2D (*fHistPion_TOF_after_PreSel);
        if(fHistPion_TOF_after_PreSel_MidPt){
            GetMinMaxBin(fHistPion_TOF_after_PreSel_MidPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel_MidPt,fHistPion_TOF_after_PreSel_MidPt->GetXaxis()->FindBin(2),fHistPion_TOF_after_PreSel_MidPt->FindBin(3));
            GetMinMaxBinY(fHistPion_TOF_after_PreSel_MidPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_PreSel_MidPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_PreSel_MidPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_PreSel_MidPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_PreSel_MidPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel_MidPt.push_back(new TH2D(*fHistPion_TOF_after_PreSel_MidPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_PreSel_MidPt_ProjPt= (TH1D*)fHistPion_TOF_after_PreSel_MidPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_PreSel_MidPt->GetXaxis()->GetFirst(),fHistPion_TOF_after_PreSel_MidPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_PreSel_MidPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel_MidPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_PreSel_MidPt_ProjPt,"",fHistPion_TOF_after_PreSel_MidPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_PreSel_MidPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel_MidPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_PreSel_MidPt_ProjPt));
            delete fHistPion_TOF_after_PreSel_MidPt_ProjPt;
            fHistPion_TOF_after_PreSel_MidPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_after_MidPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel_HighPt
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_HighPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHistPion_TOF_after_PreSel_HighPt = new TH2D (*fHistPion_TOF_after_PreSel);
        if(fHistPion_TOF_after_PreSel_HighPt){
            GetMinMaxBin(fHistPion_TOF_after_PreSel_HighPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel_HighPt,fHistPion_TOF_after_PreSel_HighPt->GetXaxis()->FindBin(3),fHistPion_TOF_after_PreSel_HighPt->GetNbinsX());
            GetMinMaxBinY(fHistPion_TOF_after_PreSel_HighPt,minYB,maxYB);
            SetYRange(fHistPion_TOF_after_PreSel_HighPt,minYB,maxYB+1);
            SetZMinMaxTH2(fHistPion_TOF_after_PreSel_HighPt,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                 fHistPion_TOF_after_PreSel_HighPt,"",
                                 "#it{p}_{T, #pi} (GeV/#it{c})","#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHistPion_TOF_after_PreSel_HighPt, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel_HighPt.push_back(new TH2D(*fHistPion_TOF_after_PreSel_HighPt));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHistPion_TOF_after_PreSel_HighPt_ProjPt= (TH1D*)fHistPion_TOF_after_PreSel_HighPt->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHistPion_TOF_after_PreSel_HighPt->GetXaxis()->GetFirst(),fHistPion_TOF_after_PreSel_HighPt->GetXaxis()->GetLast());
            GetMinMaxBin(fHistPion_TOF_after_PreSel_HighPt_ProjPt,minB,maxB);
            SetXRange(fHistPion_TOF_after_PreSel_HighPt_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistPion_TOF_after_PreSel_HighPt_ProjPt,"",fHistPion_TOF_after_PreSel_HighPt->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHistPion_TOF_after_PreSel_HighPt_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vecPion_TOF_after_PreSel_HighPt_ProjPt.push_back(new TH1D(*fHistPion_TOF_after_PreSel_HighPt_ProjPt));
            delete fHistPion_TOF_after_PreSel_HighPt_ProjPt;
            fHistPion_TOF_after_PreSel_HighPt_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: Pion_TOF_after_HighPt %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAxy_Pt_after
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAxy_Pt_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_DCAxy_Pt_after_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHisthTrack_DCAxy_Pt_after_PreSel){
            GetMinMaxBin(fHisthTrack_DCAxy_Pt_after_PreSel,minB,maxB);
            SetXRange(fHisthTrack_DCAxy_Pt_after_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_DCAxy_Pt_after_PreSel,minYB,maxYB);
            SetYRange(fHisthTrack_DCAxy_Pt_after_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_DCAxy_Pt_after_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_DCAxy_Pt_after_PreSel,"",
                                 "DCA_{#it{xy}} (cm)","#it{p}_{T} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_DCAxy_Pt_after_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAxy_Pt_after_PreSel.push_back(new TH2D(*fHisthTrack_DCAxy_Pt_after_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt= (TH1D*)fHisthTrack_DCAxy_Pt_after_PreSel->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_DCAxy_Pt_after_PreSel->GetXaxis()->GetFirst(),fHisthTrack_DCAxy_Pt_after_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt,"",fHisthTrack_DCAxy_Pt_after_PreSel->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAxy_Pt_after_PreSel_ProjPt.push_back(new TH1D(*fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt));
            delete fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt;
            fHisthTrack_DCAxy_Pt_after_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: hTrack_DCAxy_Pt_after %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAz_Pt_after
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAz_Pt_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_DCAz_Pt_after_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHisthTrack_DCAz_Pt_after_PreSel){
            GetMinMaxBin(fHisthTrack_DCAz_Pt_after_PreSel,minB,maxB);
            SetXRange(fHisthTrack_DCAz_Pt_after_PreSel,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_DCAz_Pt_after_PreSel,minYB,maxYB);
            SetYRange(fHisthTrack_DCAz_Pt_after_PreSel,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_DCAz_Pt_after_PreSel,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_DCAz_Pt_after_PreSel,"",
                                 "DCA_{#it{z}} (cm)","#it{p}_{T} (GeV/#it{c})",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_DCAz_Pt_after_PreSel, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAz_Pt_after_PreSel.push_back(new TH2D(*fHisthTrack_DCAz_Pt_after_PreSel));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_DCAz_Pt_after_PreSel_ProjPt= (TH1D*)fHisthTrack_DCAz_Pt_after_PreSel->ProjectionX(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_DCAz_Pt_after_PreSel->GetXaxis()->GetFirst(),fHisthTrack_DCAz_Pt_after_PreSel->GetXaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_DCAz_Pt_after_PreSel_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_DCAz_Pt_after_PreSel_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_DCAz_Pt_after_PreSel_ProjPt,"",fHisthTrack_DCAz_Pt_after_PreSel->GetXaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_DCAz_Pt_after_PreSel_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_DCAz_Pt_after_PreSel_ProjPt.push_back(new TH1D(*fHisthTrack_DCAz_Pt_after_PreSel_ProjPt));
            delete fHisthTrack_DCAz_Pt_after_PreSel_ProjPt;
            fHisthTrack_DCAz_Pt_after_PreSel_ProjPt=NULL;
        } else cout << Form("INFO: Object |Pre Selection: hTrack_DCAz_Pt_after %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_NFindCls_Pt_TPC_after
        if (iParticleType==0){StrNameOfHistogram="hTrack_NFindCls_Pt_TPC_after";}
        if (iParticleType==1){StrNameOfHistogram="";}
        TH2D* fHisthTrack_NFindCls_Pt_TPC_after_PreSel = (TH2D*)PionCuts2Container->FindObject(Form("%s %s",StrNameOfHistogram.Data(),fPionCuts2ContainerCutString.Data()));
        if(fHisthTrack_NFindCls_Pt_TPC_after_PreSel){
            TH2D* fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis=SwitchTF2DAxis(fHisthTrack_NFindCls_Pt_TPC_after_PreSel);
            GetMinMaxBin(fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis,minB,maxB);
            SetXRange(fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis,minB-1,maxB+1);
            GetMinMaxBinY(fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis,minYB,maxYB);
            SetYRange(fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis,minYB,maxYB+1);
            SetZMinMaxTH2(fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis,1,maxB+1,minYB,maxYB+1);
            DrawPeriodQAHistoTH2(cvsQuadratic,0.12,0.12,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                 fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis,"","#it{p}_{T} (GeV/#it{c}) TPC",
                                 "Findable Clusters after Cut",1,1.4,
                                 processLabelOffsetX2,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(cvsQuadratic, fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis, Form("%s/%s_PreSel_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_NFindCls_Pt_TPC_after_PreSel.push_back(new TH2D(*fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis));
            StrNameOfHistogram+="_ProjPt";
            TH1D* fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt= (TH1D*)fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis->ProjectionY(Form("%s",StrNameOfHistogram.Data()),fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis->GetXaxis()->GetFirst(),fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis->GetXaxis()->GetLast());
            GetMinMaxBin(fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt,minB,maxB);
            SetXRange(fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt,"",fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis->GetYaxis()->GetTitle(),"# Entries",1,1,
                                 processLabelOffsetX1,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
            SaveCanvasAndWriteHistogram(canvas, fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt, Form("%s/%s_%s.%s", outputDir.Data(),StrNameOfHistogram.Data(), DataSets[i].Data(), suffix.Data()));
            vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.push_back(new TH1D(*fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt));
            delete fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt;
            fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis_ProjPt=NULL;
            delete fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis;
            fHisthTrack_NFindCls_Pt_TPC_after_PreSel_switchedAxis=NULL;
        } else cout << Form("INFO: Object |Pre Selection: hTrack_NFindCls_Pt_TPC_after %s| could not be found! Skipping Draw...",fPionCuts2ContainerCutString.Data()) << endl;
        //-------------------------------------------------------------------------------------------------------------------------------
        //-------------------------------------------- clean up -------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------------------------------
        fFile->Close();
        delete fFile;
        fFile=NULL;
        delete TopDir;
        TopDir=NULL;
        fOutput->Close();
        delete fOutput;
        fOutput=NULL;

        //-------------------------------------------------------------------
        // display some useful values, to get a rough idea what is going on

        if(MCContainer&&TrueContainer){
            validatedefficiency = (NmbSelectedTruePosPions+NmbSelectedTrueNegPions)/(NmbSelectedMCPosPions+NmbSelectedMCNegPions);
            purity = (NmbSelectedTrueNegPions+NmbSelectedTruePosPions)/(NmbSelectedPosPions+NmbSelectedNegPions);
            cout << "--------------------------------" << endl;
            cout << "Selected Pos Pions = " << NmbSelectedPosPions << endl;
            cout << "Selected MC Pos Pions = " << NmbSelectedMCPosPions << endl;
            cout << "Selected True Pos Pions = " << NmbSelectedTruePosPions << endl;
            cout << "Selected Neg Pions = " << NmbSelectedNegPions << endl;
            cout << "Selected MC Neg Pions = " << NmbSelectedMCNegPions << endl;
            cout << "Selected True Neg Pions = " << NmbSelectedTrueNegPions << endl;
            cout << "Selected MC Pos + MC Neg Pions = " << NmbSelectedMCNegPions + NmbSelectedMCPosPions<< endl;
            cout << "Selected Pos + Pions = " << NmbSelectedNegPions + NmbSelectedPosPions<< endl;
            cout << "Pions Out = " << NmbOutPions << endl;
            cout << "--------------------------------" << endl;
            cout << "Integrated validated efficiency = " << validatedefficiency << endl;
            cout << "Integrated purity = " << purity << endl;
        }


    }
    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Comparison Histograms ************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************
    std::vector<TH1D*> temp;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms" << endl;
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------|Compare Histograms: ESD-Histograms|------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryNegPions_Pt
    if ((Int_t)vecESD_PrimaryNegPions_Pt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryNegPions_Pt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryNegPions_Pt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryNegPions_Pt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryNegPions_Pt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryNegPions_Pt,"",vecESD_PrimaryNegPions_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryNegPions_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryNegPions_Pt,"",vecESD_PrimaryNegPions_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryNegPions_Pt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPosPions_Pt
    if ((Int_t)vecESD_PrimaryPosPions_Pt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPosPions_Pt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPosPions_Pt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPosPions_Pt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPosPions_Pt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPosPions_Pt,"",vecESD_PrimaryPosPions_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPosPions_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPosPions_Pt,"",vecESD_PrimaryPosPions_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPosPions_Pt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryNegPions_Phi
    if ((Int_t)vecESD_PrimaryNegPions_Phi.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryNegPions_Phi";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryNegPions_Phi,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryNegPions_Phi.size(); iVec++){
            TH1D* temp = vecESD_PrimaryNegPions_Phi.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryNegPions_Phi,"",vecESD_PrimaryNegPions_Phi.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryNegPions_Phi.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryNegPions_Phi,"",vecESD_PrimaryNegPions_Phi.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryNegPions_Phi.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPosPions_Phi
    if ((Int_t)vecESD_PrimaryPosPions_Phi.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPosPions_Phi";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPosPions_Phi,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPosPions_Phi.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPosPions_Phi.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPosPions_Phi,"",vecESD_PrimaryPosPions_Phi.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPosPions_Phi.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPosPions_Phi,"",vecESD_PrimaryPosPions_Phi.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPosPions_Phi.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryNegPions_Eta
    if ((Int_t)vecESD_PrimaryNegPions_Eta.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryNegPions_Eta";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryNegPions_Eta,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryNegPions_Eta.size(); iVec++){
            TH1D* temp = vecESD_PrimaryNegPions_Eta.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryNegPions_Eta,"",vecESD_PrimaryNegPions_Eta.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryNegPions_Eta.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryNegPions_Eta,"",vecESD_PrimaryNegPions_Eta.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryNegPions_Eta.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPosPions_Eta
    if ((Int_t)vecESD_PrimaryPosPions_Eta.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPosPions_Eta";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPosPions_Eta,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPosPions_Eta.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPosPions_Eta.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPosPions_Eta,"",vecESD_PrimaryPosPions_Eta.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPosPions_Eta.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPosPions_Eta,"",vecESD_PrimaryPosPions_Eta.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPosPions_Eta.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdx_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdx_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdx_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdx_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdx_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdx_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdx_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdx_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdx_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdx_LowPt_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdx_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdx_MidPt_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdx_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdx_HighPt_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdx_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdx_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdxSignal_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdxSignal_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdxSignal_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdxSignal_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdxSignal_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdxSignal_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdxSignal_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdxSignal_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdxSignal_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // ESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt
    if ((Int_t)vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="ESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt,"",vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    //---------------------------------------|Compare Histograms: MC vs True Histograms|---------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    if (isMC){
        for (Int_t iVec=0; iVec<(Int_t)vecMCvsTrueAllPosPions_Pt.size();iVec++){
            vecMCvsTrueAllPosPions_Pt[iVec].at(0)=vecMC_AllPosPions_Pt.at(iVec);
            vecMCvsTrueAllPosPions_Pt[iVec].at(1)=(TH1D*)vecESD_TruePosPion_Pt.at(iVec)->Clone("ESD_TruePosPion_PtClone");
            vecMCvsTrueAllNegPions_Pt[iVec].at(0)=vecMC_AllNegPions_Pt.at(iVec);
            vecMCvsTrueAllNegPions_Pt[iVec].at(1)=(TH1D*)vecESD_TrueNegPion_Pt.at(iVec)->Clone("ESD_TrueNegPion_PtClone");
            vecMCvsTruePosPionsFromNeutralMeson_Pt[iVec].at(0)=vecMC_PosPionsFromNeutralMeson_Pt.at(iVec);
            vecMCvsTruePosPionsFromNeutralMeson_Pt[iVec].at(1)=(TH1D*)vecESD_TruePosPionFromNeutralMeson_Pt.at(iVec)->Clone("ESD_TruePosPionFromNeutralMeson_PtClone");
            vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVec].at(0)=vecMC_NegPionsFromNeutralMeson_Pt.at(iVec);
            vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVec].at(1)=(TH1D*)vecESD_TrueNegPionFromNeutralMeson_Pt.at(iVec)->Clone("ESD_TrueNegPionFromNeutralMeson_PtClone");
        }
        for (Int_t iVecLoop=0; iVecLoop<(Int_t)vecMCvsTrueAllPosPions_Pt.size();iVecLoop++){
            //-------------------------------------------------------------------------------------------------------------------------------
            // MCvsTrueAllPosPions_Pt
            if ((Int_t)vecMCvsTrueAllPosPions_Pt[iVecLoop].size()>0){
                GetMinMaxBin(vecMCvsTrueAllPosPions_Pt[iVecLoop],minB,maxB);
                for(Int_t iVec=0; iVec<(Int_t)vecMCvsTrueAllPosPions_Pt[iVecLoop].size(); iVec++){
                    TH1D* temp = vecMCvsTrueAllPosPions_Pt[iVecLoop].at(iVec);
                    temp->Sumw2();
                    //temp->Scale(1./temp->Integral());
                    //temp->Scale(1./nEvents[iVec]);
                    SetXRange(temp,minB,maxB);
                }
                TGaxis::SetExponentOffset(-0.03, -0.04, "x");
                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMCvsTrueAllPosPions_Pt[iVecLoop],"",vecMCvsTrueAllPosPions_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Counts",1,1.1,
                                            "MC", colorCompare, kTRUE, 5, 5, kFALSE,
                                            0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/MCvsTrueAllPosPions_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()), kFALSE, kTRUE);
                TGaxis::SetExponentOffset(0, 0, "x");
                DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                                 vecMCvsTrueAllPosPions_Pt[iVecLoop],"",vecMCvsTrueAllPosPions_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                                 "MC", colorCompare, kTRUE, 5, 5, kTRUE,
                                                 0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MCvsTrueAllPosPions_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()));
            }
            //-------------------------------------------------------------------------------------------------------------------------------
            // MCvsTrueAllNegPions_Pt
            if ((Int_t)vecMCvsTrueAllNegPions_Pt[iVecLoop].size()>0){
                GetMinMaxBin(vecMCvsTrueAllNegPions_Pt[iVecLoop],minB,maxB);
                for(Int_t iVec=0; iVec<(Int_t)vecMCvsTrueAllNegPions_Pt[iVecLoop].size(); iVec++){
                    TH1D* temp = vecMCvsTrueAllNegPions_Pt[iVecLoop].at(iVec);
                    temp->Sumw2();
                    //temp->Scale(1./temp->Integral());
                    //temp->Scale(1./nEvents[iVec]);
                    SetXRange(temp,minB,maxB);
                }
                TGaxis::SetExponentOffset(-0.03, -0.04, "x");
                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMCvsTrueAllNegPions_Pt[iVecLoop],"",vecMCvsTrueAllNegPions_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Counts",1,1.1,
                                            "MC", colorCompare, kTRUE, 5, 5, kFALSE,
                                            0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/MCvsTrueAllNegPions_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()), kFALSE, kTRUE);
                TGaxis::SetExponentOffset(0, 0, "x");
                DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                                 vecMCvsTrueAllNegPions_Pt[iVecLoop],"",vecMCvsTrueAllNegPions_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                                 "MC", colorCompare, kTRUE, 5, 5, kTRUE,
                                                 0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MCvsTrueAllNegPions_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()));
            }
            //-------------------------------------------------------------------------------------------------------------------------------
            // MCvsTruePosPionsFromNeutralMeson_Pt
            if ((Int_t)vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop].size()>0){
                GetMinMaxBin(vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop],minB,maxB);
                for(Int_t iVec=0; iVec<(Int_t)vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop].size(); iVec++){
                    TH1D* temp = vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop].at(iVec);
                    temp->Sumw2();
                    //->Scale(1./temp->Integral());
                    //temp->Scale(1./nEvents[iVec]);
                    SetXRange(temp,minB,maxB);
                }
                TGaxis::SetExponentOffset(-0.03, -0.04, "x");
                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop],"",vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Counts",1,1.1,
                                            "MC", colorCompare, kTRUE, 5, 5, kFALSE,
                                            0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/MCvsTruePosPionsFromNeutralMeson_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()), kFALSE, kTRUE);
                TGaxis::SetExponentOffset(0, 0, "x");
                DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                                 vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop],"",vecMCvsTruePosPionsFromNeutralMeson_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                                 "MC", colorCompare, kTRUE, 5, 5, kTRUE,
                                                 0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MCvsTruePosPionsFromNeutralMeson_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()));
            }
            //-------------------------------------------------------------------------------------------------------------------------------
            // MCvsTrueNegPionsFromNeutralMeson_Pt
            if ((Int_t)vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop].size()>0){
                GetMinMaxBin(vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop],minB,maxB);
                for(Int_t iVec=0; iVec<(Int_t)vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop].size(); iVec++){
                    TH1D* temp = vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop].at(iVec);
                    temp->Sumw2();
                    //temp->Scale(1./temp->Integral());
                    //temp->Scale(1./nEvents[iVec]);
                    SetXRange(temp,minB,maxB);
                }
                TGaxis::SetExponentOffset(-0.03, -0.04, "x");
                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop],"",vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Counts",1,1.1,
                                            "MC", colorCompare, kTRUE, 5, 5, kFALSE,
                                            0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/MCvsTrueNegPionsFromNeutralMeson_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()), kFALSE, kTRUE);
                TGaxis::SetExponentOffset(0, 0, "x");
                DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                                 vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop],"",vecMCvsTrueNegPionsFromNeutralMeson_Pt[iVecLoop].at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                                 "MC", colorCompare, kTRUE, 5, 5, kTRUE,
                                                 0.95,0.92,0.03,fCollisionSystem,DataSetsMC,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MCvsTrueNegPionsFromNeutralMeson_Pt_%s.%s", outputDir.Data(), DataSets[iVecLoop].Data(), suffix.Data()));
            }
        }
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------|Compare Histograms: MC-Histograms|-------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    if(isMC){
        // MC_AllPosPions_Pt
        if ((Int_t)vecMC_AllPosPions_Pt.size()>0){
            GetMinMaxBin(vecMC_AllPosPions_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecMC_AllPosPions_Pt.size(); iVec++){
                TH1D* temp = vecMC_AllPosPions_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMC_AllPosPions_Pt,"",vecMC_AllPosPions_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/MC_AllPosPions_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecMC_AllPosPions_Pt,"",vecMC_AllPosPions_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MC_AllPosPions_Pt.%s", outputDir.Data(), suffix.Data()));
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // MC_AllNegPions_Pt
        if ((Int_t)vecMC_AllNegPions_Pt.size()>0){
            GetMinMaxBin(vecMC_AllNegPions_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecMC_AllNegPions_Pt.size(); iVec++){
                TH1D* temp = vecMC_AllNegPions_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMC_AllNegPions_Pt,"",vecMC_AllNegPions_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/MC_AllNegPions_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecMC_AllNegPions_Pt,"",vecMC_AllNegPions_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MC_AllNegPions_Pt.%s", outputDir.Data(), suffix.Data()));
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // MC_PosPionsFromNeutralMeson_Pt
        if ((Int_t)vecMC_PosPionsFromNeutralMeson_Pt.size()>0){
            GetMinMaxBin(vecMC_PosPionsFromNeutralMeson_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecMC_PosPionsFromNeutralMeson_Pt.size(); iVec++){
                TH1D* temp = vecMC_PosPionsFromNeutralMeson_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMC_PosPionsFromNeutralMeson_Pt,"",vecMC_PosPionsFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/MC_PosPionsFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecMC_PosPionsFromNeutralMeson_Pt,"",vecMC_PosPionsFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MC_PosPionsFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()));
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // MC_NegPionsFromNeutralMeson_Pt
        if ((Int_t)vecMC_NegPionsFromNeutralMeson_Pt.size()>0){
            GetMinMaxBin(vecMC_NegPionsFromNeutralMeson_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecMC_NegPionsFromNeutralMeson_Pt.size(); iVec++){
                TH1D* temp = vecMC_NegPionsFromNeutralMeson_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMC_NegPionsFromNeutralMeson_Pt,"",vecMC_NegPionsFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/MC_NegPionsFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecMC_NegPionsFromNeutralMeson_Pt,"",vecMC_NegPionsFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_MC_NegPionsFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()));
        }
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------|Compare Histograms: True Histograms|-------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    if (isMC){
        // ESD_TruePosPion_Pt
        if ((Int_t)vecESD_TruePosPion_Pt.size()>0){
            GetMinMaxBin(vecESD_TruePosPion_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecESD_TruePosPion_Pt.size(); iVec++){
                TH1D* temp = vecESD_TruePosPion_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecESD_TruePosPion_Pt,"",vecESD_TruePosPion_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/ESD_TruePosPion_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecESD_TruePosPion_Pt,"",vecESD_TruePosPion_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_TruePosPion_Pt.%s", outputDir.Data(), suffix.Data()));
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_TrueNegPion_Pt
        if ((Int_t)vecESD_TrueNegPion_Pt.size()>0){
            GetMinMaxBin(vecESD_TrueNegPion_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecESD_TrueNegPion_Pt.size(); iVec++){
                TH1D* temp = vecESD_TrueNegPion_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecESD_TrueNegPion_Pt,"",vecESD_TrueNegPion_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/ESD_TrueNegPion_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecESD_TrueNegPion_Pt,"",vecESD_TrueNegPion_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_TrueNegPion_Pt.%s", outputDir.Data(), suffix.Data()));
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_TruePosPionFromNeutralMeson_Pt
        if ((Int_t)vecESD_TruePosPionFromNeutralMeson_Pt.size()>0){
            GetMinMaxBin(vecESD_TruePosPionFromNeutralMeson_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecESD_TruePosPionFromNeutralMeson_Pt.size(); iVec++){
                TH1D* temp = vecESD_TruePosPionFromNeutralMeson_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecESD_TruePosPionFromNeutralMeson_Pt,"",vecESD_TruePosPionFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/ESD_TruePosPionFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecESD_TruePosPionFromNeutralMeson_Pt,"",vecESD_TruePosPionFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_TruePosPionFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()));
        }
        //-------------------------------------------------------------------------------------------------------------------------------
        // ESD_TrueNegPionFromNeutralMeson_Pt
        if ((Int_t)vecESD_TrueNegPionFromNeutralMeson_Pt.size()>0){
            GetMinMaxBin(vecESD_TrueNegPionFromNeutralMeson_Pt,minB,maxB);
            for(Int_t iVec=0; iVec<(Int_t)vecESD_TrueNegPionFromNeutralMeson_Pt.size(); iVec++){
                TH1D* temp = vecESD_TrueNegPionFromNeutralMeson_Pt.at(iVec);
                temp->Sumw2();
                temp->Scale(1./temp->Integral());
                //temp->Scale(1./nEvents[iVec]);
                SetXRange(temp,minB,maxB);
            }
            TGaxis::SetExponentOffset(-0.03, -0.04, "x");
            DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecESD_TrueNegPionFromNeutralMeson_Pt,"",vecESD_TrueNegPionFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                        labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/ESD_TrueNegPionFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
            TGaxis::SetExponentOffset(0, 0, "x");
            DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                             vecESD_TrueNegPionFromNeutralMeson_Pt,"",vecESD_TrueNegPionFromNeutralMeson_Pt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                             labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                             0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
            SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ESD_TrueNegPionFromNeutralMeson_Pt.%s", outputDir.Data(), suffix.Data()));
        }
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------|Compare Histograms: AfterQA-Histograms|-----------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    //AfterQA
    // AfterQA: IsPionSelected
    if ((Int_t)vecIsPionSelected_AfterQA.size()>0){
        if (iParticleType==0){StrNameOfHistogram="IsPionSelected_AfterQA";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecIsPionSelected_AfterQA,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecIsPionSelected_AfterQA.size(); iVec++){
            TH1D* temp = vecIsPionSelected_AfterQA.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecIsPionSelected_AfterQA,"",vecIsPionSelected_AfterQA.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    plotDataSets[0], colorCompare, kTRUE, 5, 5, kFALSE,
                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/IsPionSelected_AfterQA.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecIsPionSelected_AfterQA,"",vecIsPionSelected_AfterQA.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         plotDataSets[0], colorCompare, kTRUE, 5, 5, kTRUE,
                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_IsPionSelected_AfterQA.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: dEdxCuts
    if ((Int_t)vecdEdxCuts_AfterQA.size()>0){
        if (iParticleType==0){StrNameOfHistogram="dEdxCuts_AfterQA";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecdEdxCuts_AfterQA,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecdEdxCuts_AfterQA.size(); iVec++){
            TH1D* temp = vecdEdxCuts_AfterQA.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecdEdxCuts_AfterQA,"",vecdEdxCuts_AfterQA.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/dEdxCuts_AfterQA.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecdEdxCuts_AfterQA,"",vecdEdxCuts_AfterQA.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_dEdxCuts_AfterQA.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_ITS_after
    if ((Int_t)vecPion_ITS_after_AfterQA_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_AfterQA_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_AfterQA_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_AfterQA_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_AfterQA_ProjPt,"",vecPion_ITS_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_AfterQA_ProjPt,"",vecPion_ITS_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_ITS_after_LowPt
    if ((Int_t)vecPion_ITS_after_AfterQA_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_AfterQA_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_AfterQA_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_AfterQA_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_AfterQA_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_AfterQA_LowPt_ProjPt,"",vecPion_ITS_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_AfterQA_LowPt_ProjPt,"",vecPion_ITS_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_ITS_after_MidPt
    if ((Int_t)vecPion_ITS_after_AfterQA_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_AfterQA_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_AfterQA_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_AfterQA_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_AfterQA_MidPt_ProjPt,"",vecPion_ITS_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_AfterQA_MidPt_ProjPt,"",vecPion_ITS_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_ITS_after_HighPt
    if ((Int_t)vecPion_ITS_after_AfterQA_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_AfterQA_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_AfterQA_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_AfterQA_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_AfterQA_HighPt_ProjPt,"",vecPion_ITS_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_AfterQA_HighPt_ProjPt,"",vecPion_ITS_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdx_after
    if ((Int_t)vecPion_dEdx_after_AfterQA_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_AfterQA_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_AfterQA_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_AfterQA_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_AfterQA_ProjPt,"",vecPion_dEdx_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_AfterQA_ProjPt,"",vecPion_dEdx_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdx_after_LowPt
    if ((Int_t)vecPion_dEdx_after_AfterQA_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_AfterQA_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_AfterQA_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_AfterQA_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_AfterQA_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_AfterQA_LowPt_ProjPt,"",vecPion_dEdx_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_AfterQA_LowPt_ProjPt,"",vecPion_dEdx_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdx_after_MidPt
    if ((Int_t)vecPion_dEdx_after_AfterQA_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_AfterQA_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_AfterQA_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_AfterQA_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_AfterQA_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_AfterQA_MidPt_ProjPt,"",vecPion_dEdx_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_AfterQA_MidPt_ProjPt,"",vecPion_dEdx_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdx_after_HighPt
    if ((Int_t)vecPion_dEdx_after_AfterQA_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_AfterQA_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_AfterQA_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_AfterQA_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_AfterQA_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_AfterQA_HighPt_ProjPt,"",vecPion_dEdx_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_AfterQA_HighPt_ProjPt,"",vecPion_dEdx_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdxSignal_after
    if ((Int_t)vecPion_dEdxSignal_after_AfterQA_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_AfterQA_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_AfterQA_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_AfterQA_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_AfterQA_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_AfterQA_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdxSignal_after_LowPt
    if ((Int_t)vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_AfterQA_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdxSignal_after_MidPt
    if ((Int_t)vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_AfterQA_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_dEdxSignal_after_HighPt
    if ((Int_t)vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_AfterQA_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt,"",vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_TOF_after
    if ((Int_t)vecPion_TOF_after_AfterQA_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_AfterQA_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_AfterQA_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_AfterQA_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_AfterQA_ProjPt,"",vecPion_TOF_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_AfterQA_ProjPt,"",vecPion_TOF_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_TOF_after_LowPt
    if ((Int_t)vecPion_TOF_after_AfterQA_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_AfterQA_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_AfterQA_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_AfterQA_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_AfterQA_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_AfterQA_LowPt_ProjPt,"",vecPion_TOF_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_AfterQA_LowPt_ProjPt,"",vecPion_TOF_after_AfterQA_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_AfterQA_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_TOF_after_MidPt
    if ((Int_t)vecPion_TOF_after_AfterQA_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_AfterQA_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_AfterQA_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_AfterQA_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_AfterQA_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_AfterQA_MidPt_ProjPt,"",vecPion_TOF_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_AfterQA_MidPt_ProjPt,"",vecPion_TOF_after_AfterQA_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_AfterQA_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: Pion_TOF_after_HighPt
    if ((Int_t)vecPion_TOF_after_AfterQA_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_AfterQA_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_AfterQA_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_AfterQA_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_AfterQA_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_AfterQA_HighPt_ProjPt,"",vecPion_TOF_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_AfterQA_HighPt_ProjPt,"",vecPion_TOF_after_AfterQA_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_AfterQA_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: hTrack_DCAxy_Pt_after
    if ((Int_t)vechTrack_DCAxy_Pt_after_AfterQA_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAxy_Pt_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_DCAxy_Pt_after_AfterQA_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_DCAxy_Pt_after_AfterQA_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_DCAxy_Pt_after_AfterQA_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_DCAxy_Pt_after_AfterQA_ProjPt,"",vechTrack_DCAxy_Pt_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_DCAxy_Pt_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_DCAxy_Pt_after_AfterQA_ProjPt,"",vechTrack_DCAxy_Pt_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_DCAxy_Pt_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: hTrack_DCAz_Pt_after
    if ((Int_t)vechTrack_DCAz_Pt_after_AfterQA_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAz_Pt_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_DCAz_Pt_after_AfterQA_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_DCAz_Pt_after_AfterQA_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_DCAz_Pt_after_AfterQA_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_DCAz_Pt_after_AfterQA_ProjPt,"",vechTrack_DCAz_Pt_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_DCAz_Pt_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_DCAz_Pt_after_AfterQA_ProjPt,"",vechTrack_DCAz_Pt_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_DCAz_Pt_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // AfterQA: hTrack_NFindCls_Pt_TPC_after
    if ((Int_t)vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt,"",vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt,"",vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------|Compare Histograms: Pre Selection-Histograms|-------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: IsPionSelected
    if ((Int_t)vecIsPionSelected_PreSel.size()>0){
        if (iParticleType==0){StrNameOfHistogram="IsPionSelected_PreSel";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecIsPionSelected_PreSel,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecIsPionSelected_PreSel.size(); iVec++){
            TH1D* temp = vecIsPionSelected_PreSel.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecIsPionSelected_PreSel,"",vecIsPionSelected_PreSel.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/IsPionSelected_PreSel.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecIsPionSelected_PreSel,"",vecIsPionSelected_PreSel.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_IsPionSelected_PreSel.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: dEdxCuts
    if ((Int_t)vecdEdxCuts_PreSel.size()>0){
        if (iParticleType==0){StrNameOfHistogram="dEdxCuts_PreSel";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecdEdxCuts_PreSel,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecdEdxCuts_PreSel.size(); iVec++){
            TH1D* temp = vecdEdxCuts_PreSel.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecdEdxCuts_PreSel,"",vecdEdxCuts_PreSel.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/dEdxCuts_PreSel.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecdEdxCuts_PreSel,"",vecdEdxCuts_PreSel.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_dEdxCuts_PreSel.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_before
    if ((Int_t)vecPion_ITS_before_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_before_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_before_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_before_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_before_PreSel_ProjPt,"",vecPion_ITS_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_before_PreSel_ProjPt,"",vecPion_ITS_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_before_LowPt
    if ((Int_t)vecPion_ITS_before_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_before_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_before_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_before_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_before_PreSel_LowPt_ProjPt,"",vecPion_ITS_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_before_PreSel_LowPt_ProjPt,"",vecPion_ITS_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_before_MidPt
    if ((Int_t)vecPion_ITS_before_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_before_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_before_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_before_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_before_PreSel_MidPt_ProjPt,"",vecPion_ITS_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_before_PreSel_MidPt_ProjPt,"",vecPion_ITS_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_before_HighPt
    if ((Int_t)vecPion_ITS_before_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_before_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_before_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_before_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_before_PreSel_HighPt_ProjPt,"",vecPion_ITS_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_before_PreSel_HighPt_ProjPt,"",vecPion_ITS_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_before
    if ((Int_t)vecPion_dEdx_before_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_before_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_before_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_before_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_before_PreSel_ProjPt,"",vecPion_dEdx_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_before_PreSel_ProjPt,"",vecPion_dEdx_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_before_LowPt
    if ((Int_t)vecPion_dEdx_before_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_before_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_before_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_before_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_before_PreSel_LowPt_ProjPt,"",vecPion_dEdx_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_before_PreSel_LowPt_ProjPt,"",vecPion_dEdx_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_before_MidPt
    if ((Int_t)vecPion_dEdx_before_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_before_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_before_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_before_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_before_PreSel_MidPt_ProjPt,"",vecPion_dEdx_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_before_PreSel_MidPt_ProjPt,"",vecPion_dEdx_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_before_HighPt
    if ((Int_t)vecPion_dEdx_before_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_before_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_before_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_before_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_before_PreSel_HighPt_ProjPt,"",vecPion_dEdx_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_before_PreSel_HighPt_ProjPt,"",vecPion_dEdx_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_before
    if ((Int_t)vecPion_dEdxSignal_before_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_before_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_before_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_before_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_before_PreSel_ProjPt,"",vecPion_dEdxSignal_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_before_PreSel_ProjPt,"",vecPion_dEdxSignal_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_before_LowPt
    if ((Int_t)vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt,"",vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt,"",vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_before
    if ((Int_t)vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt,"",vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt,"",vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_before_HighPt
    if ((Int_t)vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt,"",vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt,"",vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_before
    if ((Int_t)vecPion_TOF_before_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_before_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_before_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_before_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_before_PreSel_ProjPt,"",vecPion_TOF_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_before_PreSel_ProjPt,"",vecPion_TOF_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_before_LowPt
    if ((Int_t)vecPion_TOF_before_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_before_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_before_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_before_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_before_PreSel_LowPt_ProjPt,"",vecPion_TOF_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_before_PreSel_LowPt_ProjPt,"",vecPion_TOF_before_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_before_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_before_MidPt
    if ((Int_t)vecPion_TOF_before_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_before_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_before_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_before_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_before_PreSel_MidPt_ProjPt,"",vecPion_TOF_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_before_PreSel_MidPt_ProjPt,"",vecPion_TOF_before_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_before_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_before_HighPt
    if ((Int_t)vecPion_TOF_before_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_before_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_before_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_before_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_before_PreSel_HighPt_ProjPt,"",vecPion_TOF_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_before_PreSel_HighPt_ProjPt,"",vecPion_TOF_before_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_before_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: hTrack_DCAxy_Pt_before
    if ((Int_t)vechTrack_DCAxy_Pt_before_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAxy_Pt_before_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_DCAxy_Pt_before_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_DCAxy_Pt_before_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_DCAxy_Pt_before_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_DCAxy_Pt_before_PreSel_ProjPt,"",vechTrack_DCAxy_Pt_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_DCAxy_Pt_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_DCAxy_Pt_before_PreSel_ProjPt,"",vechTrack_DCAxy_Pt_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_DCAxy_Pt_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: hTrack_DCAz_Pt_before
    if ((Int_t)vechTrack_DCAz_Pt_before_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAz_Pt_before_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_DCAz_Pt_before_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_DCAz_Pt_before_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_DCAz_Pt_before_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_DCAz_Pt_before_PreSel_ProjPt,"",vechTrack_DCAz_Pt_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_DCAz_Pt_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_DCAz_Pt_before_PreSel_ProjPt,"",vechTrack_DCAz_Pt_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_DCAz_Pt_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: hTrack_NFindCls_Pt_TPC_before
    if ((Int_t)vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt,"",vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt,"",vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_after
    if ((Int_t)vecPion_ITS_after_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_PreSel_ProjPt,"",vecPion_ITS_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_PreSel_ProjPt,"",vecPion_ITS_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_after_LowPt
    if ((Int_t)vecPion_ITS_after_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_PreSel_LowPt_ProjPt,"",vecPion_ITS_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_PreSel_LowPt_ProjPt,"",vecPion_ITS_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_after_MidPt
    if ((Int_t)vecPion_ITS_after_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_PreSel_MidPt_ProjPt,"",vecPion_ITS_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_PreSel_MidPt_ProjPt,"",vecPion_ITS_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_ITS_after_HighPt
    if ((Int_t)vecPion_ITS_after_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_ITS_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_ITS_after_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_ITS_after_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_ITS_after_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_ITS_after_PreSel_HighPt_ProjPt,"",vecPion_ITS_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_ITS_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_ITS_after_PreSel_HighPt_ProjPt,"",vecPion_ITS_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_ITS_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_after
    if ((Int_t)vecPion_dEdx_after_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_PreSel_ProjPt,"",vecPion_dEdx_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_PreSel_ProjPt,"",vecPion_dEdx_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_after_LowPt
    if ((Int_t)vecPion_dEdx_after_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_PreSel_LowPt_ProjPt,"",vecPion_dEdx_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_PreSel_LowPt_ProjPt,"",vecPion_dEdx_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_after_MidPt
    if ((Int_t)vecPion_dEdx_after_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_PreSel_MidPt_ProjPt,"",vecPion_dEdx_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_PreSel_MidPt_ProjPt,"",vecPion_dEdx_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdx_after_HighPt
    if ((Int_t)vecPion_dEdx_after_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdx_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdx_after_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdx_after_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdx_after_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdx_after_PreSel_HighPt_ProjPt,"",vecPion_dEdx_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdx_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdx_after_PreSel_HighPt_ProjPt,"",vecPion_dEdx_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdx_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_after
    if ((Int_t)vecPion_dEdxSignal_after_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_PreSel_ProjPt,"",vecPion_dEdxSignal_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_PreSel_ProjPt,"",vecPion_dEdxSignal_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_after_LowPt
    if ((Int_t)vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt,"",vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt,"",vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_after_MidPt
    if ((Int_t)vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt,"",vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt,"",vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_dEdxSignal_after_HighPt
    if ((Int_t)vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_dEdxSignal_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt,"",vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_dEdxSignal_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt,"",vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_dEdxSignal_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_after
    if ((Int_t)vecPion_TOF_after_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_PreSel_ProjPt,"",vecPion_TOF_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_PreSel_ProjPt,"",vecPion_TOF_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_after_LowPt
    if ((Int_t)vecPion_TOF_after_PreSel_LowPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_PreSel_LowPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_PreSel_LowPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_PreSel_LowPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_PreSel_LowPt_ProjPt,"",vecPion_TOF_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_PreSel_LowPt_ProjPt,"",vecPion_TOF_after_PreSel_LowPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_PreSel_LowPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_after_MidPt
    if ((Int_t)vecPion_TOF_after_PreSel_MidPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_PreSel_MidPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_PreSel_MidPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_PreSel_MidPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_PreSel_MidPt_ProjPt,"",vecPion_TOF_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_PreSel_MidPt_ProjPt,"",vecPion_TOF_after_PreSel_MidPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_PreSel_MidPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: Pion_TOF_after_HighPt
    if ((Int_t)vecPion_TOF_after_PreSel_HighPt_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="Pion_TOF_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vecPion_TOF_after_PreSel_HighPt_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vecPion_TOF_after_PreSel_HighPt_ProjPt.size(); iVec++){
            TH1D* temp = vecPion_TOF_after_PreSel_HighPt_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecPion_TOF_after_PreSel_HighPt_ProjPt,"",vecPion_TOF_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Pion_TOF_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vecPion_TOF_after_PreSel_HighPt_ProjPt,"",vecPion_TOF_after_PreSel_HighPt_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Pion_TOF_after_PreSel_HighPt_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: hTrack_DCAxy_Pt_after
    if ((Int_t)vechTrack_DCAxy_Pt_after_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAxy_Pt_after_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_DCAxy_Pt_after_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_DCAxy_Pt_after_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_DCAxy_Pt_after_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_DCAxy_Pt_after_PreSel_ProjPt,"",vechTrack_DCAxy_Pt_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_DCAxy_Pt_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_DCAxy_Pt_after_PreSel_ProjPt,"",vechTrack_DCAxy_Pt_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_DCAxy_Pt_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: hTrack_DCAz_Pt_after
    if ((Int_t)vechTrack_DCAz_Pt_after_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_DCAz_Pt_after_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_DCAz_Pt_after_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_DCAz_Pt_after_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_DCAz_Pt_after_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_DCAz_Pt_after_PreSel_ProjPt,"",vechTrack_DCAz_Pt_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_DCAz_Pt_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_DCAz_Pt_after_PreSel_ProjPt,"",vechTrack_DCAz_Pt_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_DCAz_Pt_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //-------------------------------------------------------------------------------------------------------------------------------
    // Pre Selection: hTrack_NFindCls_Pt_TPC_after
    if ((Int_t)vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.size()>0){
        if (iParticleType==0){StrNameOfHistogram="hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt";}
        if (iParticleType==1){StrNameOfHistogram="";}
        GetMinMaxBin(vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt,minB,maxB);
        for(Int_t iVec=0; iVec<(Int_t)vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.size(); iVec++){
            TH1D* temp = vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.at(iVec);
            temp->Sumw2();
            temp->Scale(1./temp->Integral());
            //temp->Scale(1./nEvents[iVec]);
            SetXRange(temp,minB,maxB);
        }
        TGaxis::SetExponentOffset(-0.03, -0.04, "x");
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt,"",vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Normalized Counts",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        TGaxis::SetExponentOffset(0, 0, "x");
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                         vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt,"",vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.at(0)->GetXaxis()->GetTitle(),"Ratio",1,1.1,
                                         labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                         0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt.%s", outputDir.Data(), suffix.Data()));
    }
    //*****************************************************************************************************
    //***************************** Cleanup vectors *******************************************************
    //*****************************************************************************************************
    DeleteVecTH1D(vecESD_PrimaryNegPions_Pt);
    DeleteVecTH1D(vecESD_PrimaryPosPions_Pt);
    DeleteVecTH1D(vecESD_PrimaryNegPions_Phi);
    DeleteVecTH1D(vecESD_PrimaryPosPions_Phi);
    DeleteVecTH1D(vecESD_PrimaryNegPions_Eta);
    DeleteVecTH1D(vecESD_PrimaryPosPions_Eta);
    DeleteVecTH2D(vecESD_PrimaryNegPions_ClsTPC);
    DeleteVecTH2D(vecESD_PrimaryPosPions_ClsTPC);
    DeleteVecTH2D(vecESD_PrimaryPions_DCAxy);
    DeleteVecTH2D(vecESD_PrimaryPions_DCAz);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdx);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdx_LowPt);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdx_MidPt);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdx_HighPt);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdxSignal);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdxSignal_LowPt);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdxSignal_MidPt);
    DeleteVecTH2D(vecESD_PrimaryPions_TPCdEdxSignal_HighPt);
    DeleteVecTH1D(vecESD_PrimaryNegPions_ClsTPC_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPosPions_ClsTPC_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_DCAxy_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_DCAz_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt);
    DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt);
    //-------------------------------------------------------------------------------------------------------------------------------
    DeleteVecTH1D(vecIsPionSelected_AfterQA);
    DeleteVecTH1D(vecdEdxCuts_AfterQA);
    DeleteVecTH2D(vecPion_ITS_after_AfterQA);
    DeleteVecTH2D(vecPion_ITS_after_AfterQA_LowPt);
    DeleteVecTH2D(vecPion_ITS_after_AfterQA_MidPt);
    DeleteVecTH2D(vecPion_ITS_after_AfterQA_HighPt);
    DeleteVecTH2D(vecPion_dEdx_after_AfterQA);
    DeleteVecTH2D(vecPion_dEdx_after_AfterQA_LowPt);
    DeleteVecTH2D(vecPion_dEdx_after_AfterQA_MidPt);
    DeleteVecTH2D(vecPion_dEdx_after_AfterQA_HighPt);
    DeleteVecTH2D(vecPion_dEdxSignal_after_AfterQA);
    DeleteVecTH2D(vecPion_dEdxSignal_after_AfterQA_LowPt);
    DeleteVecTH2D(vecPion_dEdxSignal_after_AfterQA_MidPt);
    DeleteVecTH2D(vecPion_dEdxSignal_after_AfterQA_HighPt);
    DeleteVecTH2D(vecPion_TOF_after_AfterQA);
    DeleteVecTH2D(vecPion_TOF_after_AfterQA_LowPt);
    DeleteVecTH2D(vecPion_TOF_after_AfterQA_MidPt);
    DeleteVecTH2D(vecPion_TOF_after_AfterQA_HighPt);
    DeleteVecTH2D(vechTrack_DCAxy_Pt_after_AfterQA);
    DeleteVecTH2D(vechTrack_DCAz_Pt_after_AfterQA);
    DeleteVecTH2D(vechTrack_NFindCls_Pt_TPC_after_AfterQA);
    DeleteVecTH1D(vecPion_ITS_after_AfterQA_ProjPt);
    DeleteVecTH1D(vecPion_ITS_after_AfterQA_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_ITS_after_AfterQA_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_ITS_after_AfterQA_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_AfterQA_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_AfterQA_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_AfterQA_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_AfterQA_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_AfterQA_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_AfterQA_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_AfterQA_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_AfterQA_HighPt_ProjPt);
    DeleteVecTH1D(vechTrack_DCAxy_Pt_after_AfterQA_ProjPt);
    DeleteVecTH1D(vechTrack_DCAz_Pt_after_AfterQA_ProjPt);
    DeleteVecTH1D(vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt);
    //-------------------------------------------------------------------------------------------------------------------------------
    DeleteVecTH1D(vecIsPionSelected_PreSel);
    DeleteVecTH1D(vecdEdxCuts_PreSel);
    DeleteVecTH2D(vecPion_ITS_before_PreSel);
    DeleteVecTH2D(vecPion_ITS_before_PreSel_LowPt);
    DeleteVecTH2D(vecPion_ITS_before_PreSel_MidPt);
    DeleteVecTH2D(vecPion_ITS_before_PreSel_HighPt);
    DeleteVecTH2D(vecPion_dEdx_before_PreSel);
    DeleteVecTH2D(vecPion_dEdx_before_PreSel_LowPt);
    DeleteVecTH2D(vecPion_dEdx_before_PreSel_MidPt);
    DeleteVecTH2D(vecPion_dEdx_before_PreSel_HighPt);
    DeleteVecTH2D(vecPion_dEdxSignal_before_PreSel);
    DeleteVecTH2D(vecPion_dEdxSignal_before_PreSel_LowPt);
    DeleteVecTH2D(vecPion_dEdxSignal_before_PreSel_MidPt);
    DeleteVecTH2D(vecPion_dEdxSignal_before_PreSel_HighPt);
    DeleteVecTH2D(vecPion_TOF_before_PreSel);
    DeleteVecTH2D(vecPion_TOF_before_PreSel_LowPt);
    DeleteVecTH2D(vecPion_TOF_before_PreSel_MidPt);
    DeleteVecTH2D(vecPion_TOF_before_PreSel_HighPt);
    DeleteVecTH2D(vechTrack_DCAxy_Pt_before_PreSel);
    DeleteVecTH2D(vechTrack_DCAz_Pt_before_PreSel);
    DeleteVecTH2D(vechTrack_NFindCls_Pt_TPC_before_PreSel);
    DeleteVecTH2D(vecPion_ITS_after_PreSel);
    DeleteVecTH2D(vecPion_ITS_after_PreSel_LowPt);
    DeleteVecTH2D(vecPion_ITS_after_PreSel_MidPt);
    DeleteVecTH2D(vecPion_ITS_after_PreSel_HighPt);
    DeleteVecTH2D(vecPion_dEdx_after_PreSel);
    DeleteVecTH2D(vecPion_dEdx_after_PreSel_LowPt);
    DeleteVecTH2D(vecPion_dEdx_after_PreSel_MidPt);
    DeleteVecTH2D(vecPion_dEdx_after_PreSel_HighPt);
    DeleteVecTH2D(vecPion_dEdxSignal_after_PreSel);
    DeleteVecTH2D(vecPion_dEdxSignal_after_PreSel_LowPt);
    DeleteVecTH2D(vecPion_dEdxSignal_after_PreSel_MidPt);
    DeleteVecTH2D(vecPion_dEdxSignal_after_PreSel_HighPt);
    DeleteVecTH2D(vecPion_TOF_after_PreSel);
    DeleteVecTH2D(vecPion_TOF_after_PreSel_LowPt);
    DeleteVecTH2D(vecPion_TOF_after_PreSel_MidPt);
    DeleteVecTH2D(vecPion_TOF_after_PreSel_HighPt);
    DeleteVecTH2D(vechTrack_DCAxy_Pt_after_PreSel);
    DeleteVecTH2D(vechTrack_DCAz_Pt_after_PreSel);
    DeleteVecTH2D(vechTrack_NFindCls_Pt_TPC_after_PreSel);
    DeleteVecTH1D(vecPion_ITS_before_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_ITS_before_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_ITS_before_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_ITS_before_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_before_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_before_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_before_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_before_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_before_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_TOF_before_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_before_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_before_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vechTrack_DCAxy_Pt_before_PreSel_ProjPt);
    DeleteVecTH1D(vechTrack_DCAz_Pt_before_PreSel_ProjPt);
    DeleteVecTH1D(vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_ITS_after_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_ITS_after_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_ITS_after_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_ITS_after_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdx_after_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_PreSel_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_PreSel_LowPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_PreSel_MidPt_ProjPt);
    DeleteVecTH1D(vecPion_TOF_after_PreSel_HighPt_ProjPt);
    DeleteVecTH1D(vechTrack_DCAxy_Pt_after_PreSel_ProjPt);
    DeleteVecTH1D(vechTrack_DCAz_Pt_after_PreSel_ProjPt);
    DeleteVecTH1D(vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt);
    //-------------------------------------------------------------------------------------------------------------------------------
    fLog.close();
    //-------------------------------------------------------------------------------------------------------------------------------
    delete cvsQuadratic;
    cvsQuadratic=NULL;
    delete canvas;
    canvas=NULL;
    //-------------------------------------------------------------------------------------------------------------------------------
    delete[] nEvents;
    delete[] nEventsAll;
    //-------------------------------------------------------------------------------------------------------------------------------
    delete[] fTrigger;
    delete[] fCutSelection;
    delete[] fTypeCutSelection;
    delete[] fEventCutSelection;
    delete[] fGammaCutSelection;
    delete[] fClusterCutSelection;
    delete[] fPionCutSelection;
    delete[] fNeutralPionCutSelection;
    delete[] fMesonCutSelection;
    //-------------------------------------------------------------------------------------------------------------------------------
    TH1::AddDirectory(kTRUE);
    //-------------------------------------------------------------------------------------------------------------------------------
    cout << "Done with PrimaryTrackQA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    //-------------------------------------------------------------------------------------------------------------------------------
    return;
}//end

