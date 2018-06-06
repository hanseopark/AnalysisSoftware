/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 ******     Friederike Bock, fbock@cern.ch                                 *****
 *******************************************************************************/
#include "QA.h"
//**************************************************************************************************************
//***************************** Main routine *******************************************************************
//**************************************************************************************************************
void PrimaryTrackQA_Runwise(
        Int_t nSetsIn,                                      // number of sets to be analysed
        Int_t nDataIn,                                      // number of real data sets to be analysed
        TString fEnergyFlag,                                // energy flag
        TString filePath,                                   // path to the data
        TString fileName,                                   // file name of the data
        TString* DataSets,                                  // technical names of data sets for output
        TString* plotDataSets,                              // labels of data sets in plots
        Int_t mode                      = 2,                // standard mode for analysis
        Int_t cutNr                     = -1,               // if -1: you have to choose number at runtime
        Int_t doExtQA                   = 2,                // 0: switched off, 1: normal extQA, 2: with Cell level plots, 3: with mean value calculations
        Bool_t doEquidistantXaxis       = kFALSE,           // kTRUE: each run in runlist corresponds to 1 bin in X in histogram,
        // kFALSE: histograms contain the complete specified run number range, where each run represents a bin - even if it is not specified
        Bool_t doTrigger                = kTRUE,            // enables trigger analysis
        Bool_t doHistsForEverySet       = kTRUE,            // kTRUE: output done for each set separately as well
        // kFALSE: only full run range output is produced
        Bool_t addSubFolder             = kFALSE,           // kTRUE: adds another subfolder for QA output fo reach DataSet[i]
        // kFALSE: stores the runwise output all together
        Bool_t useDataRunListForMC      = kFALSE,           // kTRUE: use the same run list for data and MC
        // kFALSE: use specified
        Size_t markerSize               = 1,                // how large should the markers be?
        TString suffix                  = "eps",            // output format of plots
        TString folderRunlists          = "",               // path to the runlists
        TString addLabelRunList         = ""                // additional name for runlist
        )
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "PrimaryTrackQA_Runwise" << endl;
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
    
    const Int_t nSets   = nSetsIn;
    const Int_t nData   = nDataIn;
    
    const Int_t maxSets = 20;
    if (nSets>maxSets){
        cout << "Maximum hardcoded number of Data Sets: " << maxSets << endl;
        cout << "You have chosen: " << nSets << ", returning!" << endl;
        return;
    }
    
    Int_t fMode         = mode;
    Bool_t isCalo       = kFALSE;
    Bool_t isConv       = kFALSE;
    Bool_t isMC     = kFALSE;
    Int_t iParticleType       =0;
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
    
    
    TString fDate               = ReturnDateString();
    TString fTextMeasurement    = Form("#omega #rightarrow #pi^{0} #pi^{+} #pi^{-}, #omega #rightarrow #pi^{0} #gamma");
    TString fCentrality[30];
    for(Int_t i=0; i<nSets; i++) {
        if (fEnergyFlag.Contains("PbPb")){
            if (plotDataSets[i].Contains("0-10%")) fCentrality[i] = "0-10%";
            else if (plotDataSets[i].Contains("10-20%")) fCentrality[i] = "10-20%";
            else if (plotDataSets[i].Contains("20-50%")) fCentrality[i] = "20-50%";
            else if (plotDataSets[i].Contains("50-90%")) fCentrality[i] = "50-90%";
            else fCentrality[i] = "";
        } else {
            fCentrality[i] = "";
        }
    }
    TString fCollisionSystem    = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }
    
    TString fDetectionProcess   = ReturnFullTextReconstructionProcess(fMode);
    TString nameMainDir;
    
    //**************************************************************************************************************
    //****************************** Define plotting settings ******************************************************
    //**************************************************************************************************************
    std::vector<TString> vecDataSet;
    Style_t hMarkerStyle[maxSets];
    Size_t hMarkerSize[maxSets];
    Color_t hMarkerColor[maxSets];
    Color_t hLineColor[maxSets];
    
    for(Int_t i=0; i<nSets; i++){
        vecDataSet.push_back(DataSets[i].Data());
        hMarkerStyle[i]         = GetDefaultMarkerStyle(fEnergyFlag.Data(),DataSets[i].Data(),fCentrality[i].Data());
        hMarkerColor[i]         = GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),fCentrality[i].Data());
        hLineColor[i]           = GetColorDefaultColor(fEnergyFlag.Data(),DataSets[i].Data(),fCentrality[i].Data());
        hMarkerSize[i]          = markerSize;
    }
    
    Float_t xPosLabel = 0.8;
    Bool_t drawVerticalLines = kFALSE;
    Int_t nLines;        // number of vertical lines
    Int_t runRanges[10]; // array of bin numbers where to draw vertical lines
    TLine* verticalLines[10] = {NULL};
    if (fEnergyFlag.CompareTo("pPb_5.023TeV") == 0)
        xPosLabel = 0.75;
    if (fEnergyFlag.Contains("PbPb")){
        xPosLabel = 0.75;
        if (fEnergyFlag.CompareTo("PbPb_5.02TeV") == 0){
            drawVerticalLines = kTRUE;
            nLines = 4;
            runRanges[0] = 5; runRanges[1] = 37; runRanges[2] = 70; runRanges[3] = 72;
        }
    }
    if (nLines > 10) cout << "ERROR: nLines cannot be larger than 10. Increase size of runRanges[10] and verticalLines[10]" << endl;
    
    //*************************************************************************************************************
    // runNumbers
    std::vector<TString> vecRuns;
    TString fileRuns[maxSets];
    // bad QA runs
    std::vector<TString> vecRunsBad;
    TString fileRunsBad[maxSets];
    
    for(Int_t i=0; i<nSets; i++){
        fileRuns[i]             = Form("%s/runNumbers%s%s.txt", folderRunlists.Data(), (vecDataSet.at(i)).Data(),addLabelRunList.Data());
        fileRunsBad[i]          = Form("%s/runNumbers%s%sBadQA.txt", folderRunlists.Data(), (vecDataSet.at(i)).Data(),addLabelRunList.Data());
        
        if (useDataRunListForMC && i>=nData) {
            fileRuns[i]         = Form("%s/runNumbers%s%s-%s.txt", folderRunlists.Data(), vecDataSet.at(i).Data(), addLabelRunList.Data(),vecDataSet.at(0).Data());
            cout << "Switch useDataRunListForMC is true, reading runs from: " << fileRuns[i].Data() << endl;
        }
        cout << "trying to read: " << fileRuns[i].Data() << endl;
        if (!readin(fileRuns[i], vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}
    }
    
    //*************************************************************************************************************
    //****************************** Determine which cut to process ***********************************************
    //*************************************************************************************************************
    TFile* fCutFile             = new TFile(Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data()));
    if (fCutFile->IsZombie()) {cout << "ERROR: ROOT file '" << Form("%s/%s/%s/%s", filePath.Data(), ((TString)vecDataSet.at(0)).Data(), ((TString)vecRuns.at(0)).Data(), fileName.Data()) << "' could not be openend, return!" << endl; return;}
    
    TKey *key;
    TIter next(fCutFile->GetListOfKeys());
    while ((key=(TKey*)next())){
        cout << Form("Found TopDir: '%s' ",key->GetName());
        nameMainDir             = key->GetName();
    }
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    if (nameMainDir.IsNull() || !nameMainDir.BeginsWith("Gamma")){cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir.Data() << "|, running in mode: " << fMode << endl; return;}
    
    TList *listInput            = (TList*)fCutFile->Get(nameMainDir.Data());
    if (!listInput) {cout << "ERROR: Could not find main dir: " << nameMainDir.Data() << " in file! Returning..." << endl; return;}
    listInput->SetOwner(kTRUE);
    vector <TString> cuts;
    TString fPionCuts2ContainerCutString;
    for(Int_t i = 0; i<listInput->GetSize(); i++){
        TList *listCuts         = (TList*)listInput->At(i);
        TString nameCuts        = listCuts->GetName();
        if (nameCuts.BeginsWith("Cut Number")){
            nameCuts.Replace(0,11,"");
            cuts.push_back(nameCuts);
        }
        if (nameCuts.BeginsWith("PionCuts_")){
            nameCuts.Replace(0,9,"");
            fPionCuts2ContainerCutString=nameCuts.Data();
            cout<<"fPionCuts2ContainerCutString: "<<fPionCuts2ContainerCutString<<endl;
        }
    }
    delete listInput;
    
    cout << "The following cuts are available:" << endl;
    for(Int_t i = 0; i < (Int_t) cuts.size(); i++) {cout << Form("(%i) -- %s", i, cuts[i].Data()) << endl;}
    
    if (cutNr == -1){
        do{ cin >> cutNr;}
        while( (cutNr < 0) || (cutNr > (Int_t) cuts.size()) );
    }
    
    fCutFile->Close();
    delete fCutFile;
    
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Processing Cut Number: " << cutNr << endl;
    cout << cuts.at(cutNr) << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Pictures are saved as " << suffix.Data() << "!" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    
    TString fCutSelection            = cuts.at(cutNr);
    TString fTypeCutSelection        = "";
    TString fEventCutSelection       = "";
    TString fGammaCutSelection       = "";
    TString fClusterCutSelection     = "";
    TString fMClusterCutSelection    = "";
    TString fPionCutSelection        = "";
    TString fNeutralPionCutSelection = "";
    TString fMesonCutSelection       = "";
    fMode=ReturnSeparatedCutNumberPiPlPiMiPiZero(fCutSelection, fTypeCutSelection, fEventCutSelection, fGammaCutSelection, fClusterCutSelection,fPionCutSelection, fNeutralPionCutSelection, fMesonCutSelection);
    cout<<"fCutSelection: "<<fCutSelection<<endl<<"fTypeCutSelection: "<<fTypeCutSelection<<endl<<"fEventCutSelection: "<<fEventCutSelection<<endl<<"fGammaCutSelection: "<<fGammaCutSelection<<endl<<"fClusterCutSelection: "<<fClusterCutSelection<<endl<<"fMClusterCutSelection: "<<fMClusterCutSelection<<endl<<"fPionCutSelection:"<<fPionCutSelection<<endl<<"fNeutralPionCutSelection: "<<fNeutralPionCutSelection<<endl<<"fMesonCutSelection: "<<fMesonCutSelection<<endl;
    if (fMode!=mode){
        cout << "ERROR: Chosen mode ("<<mode<<") is not identical to mode extracted from cutstring ("<<fMode<<")! " << endl;
        return;
    } // check if extracted mode from cutnumber is the same mode given by user
    fMode                  = mode;
    cout << "\t MODE = " << mode << endl;
    
    
    //*****************************************************************************************************
    //************************** Set proper cluster nomenclature ******************************************
    //*****************************************************************************************************
    
    TString calo                = "";
    TString fClusters           = "";
    
    if (isCalo){
        if (fClusterCutSelection.BeginsWith('1')){
            calo                = "EMCal";
            fClusters           = Form("%s clusters", calo.Data());
        }else if (fClusterCutSelection.BeginsWith('2')){
            calo                = "PHOS";
            fClusters           = Form("%s clusters", calo.Data());
        } else if (fClusterCutSelection.BeginsWith('3')){
            calo                = "DCal";
            fClusters           = Form("%s clusters", calo.Data());
        }else {cout << "No correct calorimeter type found: " << calo.Data() << ", returning..." << endl; return;}
    }
    
    //*****************************************************************************************************
    //************************** Define output directories*************************************************
    //*****************************************************************************************************
    TString outputDir           = Form("%s/%s/PrimaryTrackQA/%s/Runwise",fCutSelection.Data(),fEnergyFlag.Data(),suffix.Data());
    if (addSubFolder){
        outputDir               += "/";
        outputDir               += DataSets[0];
    }
    gSystem->Exec("mkdir -p "+outputDir);
    
    
    //*****************************************************************************************************
    //**************************** Determine global run list **********************************************
    //*****************************************************************************************************
    std::vector<TString> globalRuns;
    Float_t rangesRuns[nSets][2];
    
    for(Int_t i=0; i<nSets; i++) {
        vecRuns.clear();
        if (!readin(fileRuns[i], vecRuns, kFALSE)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}
        
        for(Int_t j=0; j<(Int_t) vecRuns.size();j++){
            if ( i==0 ) globalRuns.push_back(vecRuns.at(j));
            else {
                Bool_t bFound = kFALSE;
                for(Int_t k=0; k<(Int_t) globalRuns.size();k++){ if (globalRuns.at(k)==vecRuns.at(j)) bFound=kTRUE;}
                if (!bFound) globalRuns.push_back(vecRuns.at(j));
            }
        }
        
        if ( !doEquidistantXaxis && ((Int_t) vecRuns.size())>0 ){
            if (nSets>2){
                rangesRuns[i][0]=((TString)vecRuns.front()).Atof() - 500.;
                rangesRuns[i][1]=((TString)vecRuns.back()).Atof() + 500.;
            }else{
                rangesRuns[i][0]=((TString)vecRuns.front()).Atof() - 25.;
                rangesRuns[i][1]=((TString)vecRuns.back()).Atof() + 25.;
            }
        }else{ rangesRuns[i][0]=0; rangesRuns[i][1]=0; }
    }
    
    selection_sort(globalRuns.begin(), globalRuns.end());
    
    map<TString,Int_t> mapBin;
    
    
    //*****************************************************************************************************
    //********************* Create histograms for plotting ************************************************
    //*****************************************************************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Processing following list of " << globalRuns.size() << " Runs:";
    for(Int_t i=0; i<(Int_t) globalRuns.size(); i++) {
        mapBin[globalRuns.at(i)]=i+1;
        if (i%10==0) cout << endl;
        cout << globalRuns.at(i) << ", ";
    }
    cout << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    
    TH1D* hESD_PrimaryNegPions_Pt[nSets];
    TH1D* hESD_PrimaryPosPions_Pt[nSets];
    TH1D* hESD_PrimaryNegPions_Phi[nSets];
    TH1D* hESD_PrimaryPosPions_Phi[nSets];
    TH1D* hESD_PrimaryNegPions_Eta[nSets];
    TH1D* hESD_PrimaryPosPions_Eta[nSets];
    //AfterQA=Cut Number Pion Cut (Pion Cut folder in Cut number directory); Pre Selection=Main Dir Pion Cut
    //AfterQA
    TH1D* hIsPionSelected_AfterQA[nSets];
    TH1D* hdEdxCuts_AfterQA[nSets];
    //Pre Selection
    TH1D* hIsPionSelected_PreSel[nSets];
    TH1D* hdEdxCuts_PreSel[nSets];
    //MC histograms
    TH1D* hMC_AllPosPions_Pt[nSets];
    TH1D* hMC_AllNegPions_Pt[nSets];
    TH1D* hMC_PosPionsFromNeutralMeson_Pt[nSets];
    TH1D* hMC_NegPionsFromNeutralMeson_Pt[nSets];
    //True histograms
    TH1D* hESD_TruePosPion_Pt[nSets];
    TH1D* hESD_TrueNegPion_Pt[nSets];
    TH1D* hESD_TruePosPionFromNeutralMeson_Pt[nSets];
    TH1D* hESD_TrueNegPionFromNeutralMeson_Pt[nSets];
    //---------------------------------------------Projections-----------------------------------------------------------------------
    TH1D* hESD_PrimaryNegPions_ClsTPC_ProjPt[nSets];                //Pt was x
    TH1D* hESD_PrimaryPosPions_ClsTPC_ProjPt[nSets];                //Pt was x
    TH1D* hESD_PrimaryPions_DCAxy_ProjPt[nSets];                    //Pt was y
    TH1D* hESD_PrimaryPions_DCAz_ProjPt[nSets];                     //Pt was y
    TH1D* hESD_PrimaryPions_TPCdEdx_ProjPt[nSets];                  //Pt was x
    TH1D* hESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[nSets];            //Pt was x
    TH1D* hESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[nSets];            //Pt was x
    TH1D* hESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[nSets];           //Pt was x
    TH1D* hESD_PrimaryPions_TPCdEdxSignal_ProjPt[nSets];            //Pt was x
    TH1D* hESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[nSets];      //Pt was x
    TH1D* hESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[nSets];      //Pt was x
    TH1D* hESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[nSets];     //Pt was x
    //AfterQA
    TH1D* hPion_ITS_after_AfterQA_ProjPt[nSets];                    //Pt was x
    TH1D* hPion_ITS_after_AfterQA_LowPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_ITS_after_AfterQA_MidPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_ITS_after_AfterQA_HighPt_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdx_after_AfterQA_ProjPt[nSets];                   //Pt was x
    TH1D* hPion_dEdx_after_AfterQA_LowPt_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdx_after_AfterQA_MidPt_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdx_after_AfterQA_HighPt_ProjPt[nSets];            //Pt was x
    TH1D* hPion_dEdxSignal_after_AfterQA_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[nSets];       //Pt was x
    TH1D* hPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[nSets];       //Pt was x
    TH1D* hPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[nSets];      //Pt was x
    TH1D* hPion_TOF_after_AfterQA_ProjPt[nSets];                    //Pt was x
    TH1D* hPion_TOF_after_AfterQA_LowPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_TOF_after_AfterQA_MidPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_TOF_after_AfterQA_HighPt_ProjPt[nSets];             //Pt was x
    TH1D* hTrack_DCAxy_Pt_after_AfterQA_ProjPt[nSets];             //Pt was y
    TH1D* hTrack_DCAz_Pt_after_AfterQA_ProjPt[nSets];              //Pt was y
    TH1D* hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[nSets];      //Pt was x
    //Pre Selection
    TH1D* hPion_ITS_before_PreSel_ProjPt[nSets];                    //Pt was x
    TH1D* hPion_ITS_before_PreSel_LowPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_ITS_before_PreSel_MidPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_ITS_before_PreSel_HighPt_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdx_before_PreSel_ProjPt[nSets];                   //Pt was x
    TH1D* hPion_dEdx_before_PreSel_LowPt_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdx_before_PreSel_MidPt_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdx_before_PreSel_HighPt_ProjPt[nSets];            //Pt was x
    TH1D* hPion_dEdxSignal_before_PreSel_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdxSignal_before_PreSel_LowPt_ProjPt[nSets];       //Pt was x
    TH1D* hPion_dEdxSignal_before_PreSel_MidPt_ProjPt[nSets];       //Pt was x
    TH1D* hPion_dEdxSignal_before_PreSel_HighPt_ProjPt[nSets];      //Pt was x
    TH1D* hPion_TOF_before_PreSel_ProjPt[nSets];                    //Pt was x
    TH1D* hPion_TOF_before_PreSel_LowPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_TOF_before_PreSel_MidPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_TOF_before_PreSel_HighPt_ProjPt[nSets];             //Pt was x
    TH1D* hTrack_DCAxy_Pt_before_PreSel_ProjPt[nSets];             //Pt was y
    TH1D* hTrack_DCAz_Pt_before_PreSel_ProjPt[nSets];              //Pt was y
    TH1D* hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[nSets];      //Pt was x
    TH1D* hPion_ITS_after_PreSel_ProjPt[nSets];                     //Pt was x
    TH1D* hPion_ITS_after_PreSel_LowPt_ProjPt[nSets];               //Pt was x
    TH1D* hPion_ITS_after_PreSel_MidPt_ProjPt[nSets];               //Pt was x
    TH1D* hPion_ITS_after_PreSel_HighPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_dEdx_after_PreSel_ProjPt[nSets];                    //Pt was x
    TH1D* hPion_dEdx_after_PreSel_LowPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_dEdx_after_PreSel_MidPt_ProjPt[nSets];              //Pt was x
    TH1D* hPion_dEdx_after_PreSel_HighPt_ProjPt[nSets];             //Pt was x
    TH1D* hPion_dEdxSignal_after_PreSel_ProjPt[nSets];              //Pt was x
    TH1D* hPion_dEdxSignal_after_PreSel_LowPt_ProjPt[nSets];        //Pt was x
    TH1D* hPion_dEdxSignal_after_PreSel_MidPt_ProjPt[nSets];        //Pt was x
    TH1D* hPion_dEdxSignal_after_PreSel_HighPt_ProjPt[nSets];       //Pt was x
    TH1D* hPion_TOF_after_PreSel_ProjPt[nSets];                     //Pt was x
    TH1D* hPion_TOF_after_PreSel_LowPt_ProjPt[nSets];               //Pt was x
    TH1D* hPion_TOF_after_PreSel_MidPt_ProjPt[nSets];               //Pt was x
    TH1D* hPion_TOF_after_PreSel_HighPt_ProjPt[nSets];              //Pt was x
    TH1D* hTrack_DCAxy_Pt_after_PreSel_ProjPt[nSets];              //Pt was y
    TH1D* hTrack_DCAz_Pt_after_PreSel_ProjPt[nSets];               //Pt was y
    TH1D* hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[nSets];       //Pt was x
    
    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Vectors for Histograms ***********************************************
    //*****************************************************************************************************
    //*****************************************************************************************************
    
    std::vector<TH1D*>* vecESD_PrimaryNegPions_Pt               = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_PrimaryPosPions_Pt               = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_PrimaryNegPions_Phi              = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_PrimaryPosPions_Phi              = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_PrimaryNegPions_Eta              = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_PrimaryPosPions_Eta              = new std::vector<TH1D*>[nSets];
    //AfterQA=Cut Number Pion Cut (Pion Cut folder in Cut number directory); Pre Selection=Main Dir Pion Cut
    //AfterQA
    std::vector<TH1D*>* vecIsPionSelected_AfterQA               = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecdEdxCuts_AfterQA                     = new std::vector<TH1D*>[nSets];
    //Pre Selection
    std::vector<TH1D*>* vecIsPionSelected_PreSel                = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecdEdxCuts_PreSel                      = new std::vector<TH1D*>[nSets];
    //MC histograms
    std::vector<TH1D*>* vecMC_AllPosPions_Pt                    = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecMC_AllNegPions_Pt                    = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecMC_PosPionsFromNeutralMeson_Pt       = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecMC_NegPionsFromNeutralMeson_Pt       = new std::vector<TH1D*>[nSets];
    //True histograms
    std::vector<TH1D*>* vecESD_TruePosPion_Pt                   = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_TrueNegPion_Pt                   = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_TruePosPionFromNeutralMeson_Pt   = new std::vector<TH1D*>[nSets];
    std::vector<TH1D*>* vecESD_TrueNegPionFromNeutralMeson_Pt   = new std::vector<TH1D*>[nSets];
    //---------------------------------------------Projections-----------------------------------------------------------------------
    std::vector<TH1D*>* vecESD_PrimaryNegPions_ClsTPC_ProjPt                = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPosPions_ClsTPC_ProjPt                = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_DCAxy_ProjPt                    = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vecESD_PrimaryPions_DCAz_ProjPt                     = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdx_ProjPt                  = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt            = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt            = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt           = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdxSignal_ProjPt            = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt      = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt      = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt     = new std::vector<TH1D*>[nSets];//Pt was x
    
    //AfterQA
    std::vector<TH1D*>* vecPion_ITS_after_AfterQA_ProjPt                    = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_after_AfterQA_LowPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_after_AfterQA_MidPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_after_AfterQA_HighPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_AfterQA_ProjPt                   = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_AfterQA_LowPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_AfterQA_MidPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_AfterQA_HighPt_ProjPt            = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_AfterQA_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt       = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt       = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt      = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_AfterQA_ProjPt                    = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_AfterQA_LowPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_AfterQA_MidPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_AfterQA_HighPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vechTrack_DCAxy_Pt_after_AfterQA_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vechTrack_DCAz_Pt_after_AfterQA_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt      = new std::vector<TH1D*>[nSets];//Pt was x
    //Pre Selection
    std::vector<TH1D*>* vecPion_ITS_before_PreSel_ProjPt                    = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_before_PreSel_LowPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_before_PreSel_MidPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_before_PreSel_HighPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_before_PreSel_ProjPt                   = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_before_PreSel_LowPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_before_PreSel_MidPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_before_PreSel_HighPt_ProjPt            = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_before_PreSel_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt       = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt       = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt      = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_before_PreSel_ProjPt                    = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_before_PreSel_LowPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_before_PreSel_MidPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_before_PreSel_HighPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vechTrack_DCAxy_Pt_before_PreSel_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vechTrack_DCAz_Pt_before_PreSel_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt      = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_after_PreSel_ProjPt                     = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_after_PreSel_LowPt_ProjPt               = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_after_PreSel_MidPt_ProjPt               = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_ITS_after_PreSel_HighPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_PreSel_ProjPt                    = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_PreSel_LowPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_PreSel_MidPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdx_after_PreSel_HighPt_ProjPt             = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_PreSel_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt        = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt        = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt       = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_PreSel_ProjPt                     = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_PreSel_LowPt_ProjPt               = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_PreSel_MidPt_ProjPt               = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vecPion_TOF_after_PreSel_HighPt_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was x
    std::vector<TH1D*>* vechTrack_DCAxy_Pt_after_PreSel_ProjPt              = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vechTrack_DCAz_Pt_after_PreSel_ProjPt               = new std::vector<TH1D*>[nSets];//Pt was y
    std::vector<TH1D*>* vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt       = new std::vector<TH1D*>[nSets];//Pt was x
    //-------------------------------------------------------------------------------------------------------------------------------
    std::vector<TH1D*>* vecHistos                               = new std::vector<TH1D*>[nSets];
    std::vector<TString> vecHistosName;
    TString histoName;
    TCanvas* canvas1                 = new TCanvas("canvas1","",10,10,750,500);//Just for Debugging!
    
    
    Int_t hFBin;
    Int_t hLBin;
    Int_t hNBin;
    Int_t minB          = 0;    Int_t maxB          = 0;
    Int_t minYB         = 0;    Int_t maxYB         = 0;
    
    if (doEquidistantXaxis)    {
        hFBin       = 0;
        hLBin       = globalRuns.size();
        hNBin       = globalRuns.size();
    } else {
        if (nSets>2){
            hFBin   = ((TString) globalRuns.at(0)).Atoi() - 500;
            hLBin   = ((TString) globalRuns.at(globalRuns.size()-1)).Atoi()  + 500;
            hNBin   = hLBin - hFBin;
        }else{
            hFBin   = ((TString) globalRuns.at(0)).Atoi() - 25;
            hLBin   = ((TString) globalRuns.at(globalRuns.size()-1)).Atoi()  + 25;
            hNBin   = hLBin - hFBin;
        }
    }
    
    for(Int_t i=0; i<nSets; i++){
        //*******************************************************************************************************************************
        //-------------------------------------------------------------------------------------------------------------------------------
        //*******************************************************************************************************************************
        if (iParticleType==0){histoName="ESD_PrimaryNegPions_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryNegPions_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryNegPions_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryNegPions_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_PrimaryPosPions_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPosPions_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #it{#pi^{+}}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPosPions_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPosPions_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_PrimaryNegPions_Phi";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryNegPions_Phi[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{#varphi}_{#it{#pi^{-}}}",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryNegPions_Phi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryNegPions_Phi[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_PrimaryPosPions_Phi";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPosPions_Phi[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Mean #it{#varphi}_{#it{#pi^{+}}}",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPosPions_Phi[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPosPions_Phi[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_PrimaryNegPions_Eta";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryNegPions_Eta[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{#eta}_{#it{#pi^{-}}}",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryNegPions_Eta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryNegPions_Eta[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_PrimaryPosPions_Eta";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPosPions_Eta[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{#eta}_{#it{#pi^{+}}}",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPosPions_Eta[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPosPions_Eta[i]);
        //*******************************************************************************************************************************
        //AfterQA=Cut Number Pion Cut (Pion Cut folder in Cut number directory); Pre Selection=Main Dir Pion Cut
        //AfterQA
        if (iParticleType==0){histoName="IsPionSelected_AfterQA_OutDivIn";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hIsPionSelected_AfterQA[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Number",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hIsPionSelected_AfterQA[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hIsPionSelected_AfterQA[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="dEdxCuts_AfterQA_OutDivIn";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hdEdxCuts_AfterQA[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Number",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hdEdxCuts_AfterQA[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hdEdxCuts_AfterQA[i]);
        //*******************************************************************************************************************************
        //Pre Selection=Main Dir Pion Cut
        if (iParticleType==0){histoName="IsPionSelected_PreSel_OutDivIn";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hIsPionSelected_PreSel[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Number",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hIsPionSelected_PreSel[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hIsPionSelected_PreSel[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="dEdxCuts_PreSel_OutDivIn";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hdEdxCuts_PreSel[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Number",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hdEdxCuts_PreSel[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hdEdxCuts_PreSel[i]);
        //*******************************************************************************************************************************
        //MC histograms
        if (iParticleType==0){histoName="MC_AllPosPions_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hMC_AllPosPions_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #pi^{+}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hMC_AllPosPions_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hMC_AllPosPions_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="MC_AllNegPions_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hMC_AllNegPions_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #pi^{-}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hMC_AllNegPions_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hMC_AllNegPions_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="MC_PosPionsFromNeutralMeson_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hMC_PosPionsFromNeutralMeson_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #pi^{+}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hMC_PosPionsFromNeutralMeson_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hMC_PosPionsFromNeutralMeson_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="MC_NegPionsFromNeutralMeson_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hMC_NegPionsFromNeutralMeson_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #pi^{-}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hMC_NegPionsFromNeutralMeson_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hMC_NegPionsFromNeutralMeson_Pt[i]);
        //*******************************************************************************************************************************
        //True histograms
        if (iParticleType==0){histoName="ESD_TruePosPion_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_TruePosPion_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #pi^{+}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_TruePosPion_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_TruePosPion_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_TrueNegPion_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_TrueNegPion_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_TrueNegPion_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_TrueNegPion_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_TruePosPionFromNeutralMeson_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_TruePosPionFromNeutralMeson_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #it{#pi^{+}}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_TruePosPionFromNeutralMeson_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_TruePosPionFromNeutralMeson_Pt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        if (iParticleType==0){histoName="ESD_TrueNegPionFromNeutralMeson_Pt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_TrueNegPionFromNeutralMeson_Pt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_TrueNegPionFromNeutralMeson_Pt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_TrueNegPionFromNeutralMeson_Pt[i]);
        //*******************************************************************************************************************************
        //--------------------------------------------------------Projections------------------------------------------------------------
        //*******************************************************************************************************************************
        //ESD_PrimaryNegPions_ClsTPC
        if (iParticleType==0){histoName="ESD_PrimaryNegPions_ClsTPC_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryNegPions_ClsTPC_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Findable Clusters",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryNegPions_ClsTPC_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryNegPions_ClsTPC_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPosPions_ClsTPC
        if (iParticleType==0){histoName="ESD_PrimaryPosPions_ClsTPC_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPosPions_ClsTPC_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Findable Clusters",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPosPions_ClsTPC_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPosPions_ClsTPC_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_DCAxy
        if (iParticleType==0){histoName="ESD_PrimaryPions_DCAxy_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_DCAxy_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{xy}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_DCAxy_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_DCAxy_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_DCAz
        if (iParticleType==0){histoName="ESD_PrimaryPions_DCAz_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_DCAz_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{z}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_DCAz_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_DCAz_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdx_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdx_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdx_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx_LowPt
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx_MidPt
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdx_HighPt
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdxSignal_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdxSignal_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdxSignal_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal_LowPt
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal_MidPt
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //ESD_PrimaryPions_TPCdEdxSignal_HighPt
        if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i]);
        //*******************************************************************************************************************************
        //AfterQA: Pion_ITS_after
        if (iParticleType==0){histoName="Pion_ITS_after_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_AfterQA_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_AfterQA_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_AfterQA_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_ITS_after_LowPt
        if (iParticleType==0){histoName="Pion_ITS_after_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_AfterQA_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_AfterQA_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_AfterQA_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_ITS_after_MidPt
        if (iParticleType==0){histoName="Pion_ITS_after_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_AfterQA_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_AfterQA_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_AfterQA_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_ITS_after_HighPt
        if (iParticleType==0){histoName="Pion_ITS_after_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_AfterQA_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_AfterQA_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_AfterQA_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after
        if (iParticleType==0){histoName="Pion_dEdx_after_AfterQA_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_AfterQA_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_AfterQA_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_AfterQA_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after_LowPt
        if (iParticleType==0){histoName="Pion_dEdx_after_AfterQA_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_AfterQA_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_AfterQA_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_AfterQA_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after_MidPt
        if (iParticleType==0){histoName="Pion_dEdx_after_AfterQA_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_AfterQA_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_AfterQA_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_AfterQA_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdx_after_HighPt
        if (iParticleType==0){histoName="Pion_dEdx_after_AfterQA_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_AfterQA_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_AfterQA_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_AfterQA_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_AfterQA_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_AfterQA_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_AfterQA_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_AfterQA_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_LowPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_AfterQA_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_MidPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_AfterQA_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_dEdxSignal_after_HighPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_AfterQA_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_TOF_after
        if (iParticleType==0){histoName="Pion_TOF_after_AfterQA_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_AfterQA_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_AfterQA_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_AfterQA_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_TOF_after_LowPt
        if (iParticleType==0){histoName="Pion_TOF_after_AfterQA_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_AfterQA_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_AfterQA_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_AfterQA_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_TOF_after_MidPt
        if (iParticleType==0){histoName="Pion_TOF_after_AfterQA_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_AfterQA_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_AfterQA_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_AfterQA_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: Pion_TOF_after_HighPt
        if (iParticleType==0){histoName="Pion_TOF_after_AfterQA_HighPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_AfterQA_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_AfterQA_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_AfterQA_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: hTrack_DCAxy_Pt_after
        if (iParticleType==0){histoName="hTrack_DCAxy_Pt_after_AfterQA";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_DCAxy_Pt_after_AfterQA_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{xy}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_DCAxy_Pt_after_AfterQA_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_DCAxy_Pt_after_AfterQA_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: hTrack_DCAz_Pt_after
        if (iParticleType==0){histoName="hTrack_DCAz_Pt_after_AfterQA";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_DCAz_Pt_after_AfterQA_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{z}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_DCAz_Pt_after_AfterQA_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_DCAz_Pt_after_AfterQA_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //AfterQA: hTrack_NFindCls_Pt_TPC_after
        if (iParticleType==0){histoName="hTrack_NFindCls_Pt_TPC_after_AfterQA";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Findable Clusters after Cut",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i]);
        
        //*******************************************************************************************************************************
        //BEFORE (In Pre Selection)
        //Pre Selection: Pion_ITS_before
        if (iParticleType==0){histoName="Pion_ITS_before_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_before_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_before_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_before_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_before_LowPt
        if (iParticleType==0){histoName="Pion_ITS_before_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_before_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_before_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_before_PreSel_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_before_MidPt
        if (iParticleType==0){histoName="Pion_ITS_before_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_before_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_before_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_before_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_before_HighPt
        if (iParticleType==0){histoName="Pion_ITS_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_before_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_before_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_before_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before
        if (iParticleType==0){histoName="Pion_dEdx_before_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_before_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_before_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_before_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before_LowPt
        if (iParticleType==0){histoName="Pion_dEdx_before_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        hPion_dEdx_before_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_before_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_before_PreSel_LowPt_ProjPt[i]);
        if (i==0) vecHistosName.push_back(histoName);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before_MidPt
        if (iParticleType==0){histoName="Pion_dEdx_before_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_before_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_before_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_before_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_before_HighPt
        if (iParticleType==0){histoName="Pion_dEdx_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_before_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_before_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_before_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before
        if (iParticleType==0){histoName="Pion_dEdxSignal_before_PreSel";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_before_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_before_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_before_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before_LowPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_before_PreSel_LowPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before_MidPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_before_PreSel_MidPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_before_HighPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{p}_{#pi} (GeV/#it{c})","d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel
        if (iParticleType==0){histoName="Pion_TOF_before_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_before_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_before_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_before_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel_LowPt
        if (iParticleType==0){histoName="Pion_TOF_before_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_before_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_before_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_before_PreSel_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel_MidPt
        if (iParticleType==0){histoName="Pion_TOF_before_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_before_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_before_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_before_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_before_PreSel_HighPt
        if (iParticleType==0){histoName="Pion_TOF_before_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_before_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_before_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_before_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAxy_Pt_before
        if (iParticleType==0){histoName="hTrack_DCAxy_Pt_before_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_DCAxy_Pt_before_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{xy}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_DCAxy_Pt_before_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_DCAxy_Pt_before_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAz_Pt_before
        if (iParticleType==0){histoName="hTrack_DCAz_Pt_before_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_DCAz_Pt_before_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{z}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_DCAz_Pt_before_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_DCAz_Pt_before_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_NFindCls_Pt_TPC_before
        if (iParticleType==0){histoName="hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Findable Clusters before Cut",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i]);
        //*******************************************************************************************************************************
        //AFTER (In Pre Selection)
        //Pre Selection: Pion_ITS_after
        if (iParticleType==0){histoName="Pion_ITS_after_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_after_LowPt
        if (iParticleType==0){histoName="Pion_ITS_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_PreSel_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_after_MidPt
        if (iParticleType==0){histoName="Pion_ITS_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_ITS_after_HighPt
        if (iParticleType==0){histoName="Pion_ITS_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_ITS_after_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_ITS_after_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_ITS_after_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after
        if (iParticleType==0){histoName="Pion_dEdx_after_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after_LowPt
        if (iParticleType==0){histoName="Pion_dEdx_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_PreSel_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after_MidPt
        if (iParticleType==0){histoName="Pion_dEdx_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdx_after_HighPt
        if (iParticleType==0){histoName="Pion_dEdx_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdx_after_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdx_after_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdx_after_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA_LowPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA_MidPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_dEdxSignal_after_AfterQA_HighPt
        if (iParticleType==0){histoName="Pion_dEdxSignal_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean d#it{E}/d#it{x} TPC",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel
        if (iParticleType==0){histoName="Pion_TOF_after_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel_LowPt
        if (iParticleType==0){histoName="Pion_TOF_after_PreSel_LowPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_PreSel_LowPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_PreSel_LowPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_PreSel_LowPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel_MidPt
        if (iParticleType==0){histoName="Pion_TOF_after_PreSel_MidPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_PreSel_MidPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_PreSel_MidPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_PreSel_MidPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: Pion_TOF_after_PreSel_HighPt
        if (iParticleType==0){histoName="Pion_TOF_after_PreSel_HighPt_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hPion_TOF_after_PreSel_HighPt_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean #it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hPion_TOF_after_PreSel_HighPt_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hPion_TOF_after_PreSel_HighPt_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAxy_Pt_after
        if (iParticleType==0){histoName="hTrack_DCAxy_Pt_after_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_DCAxy_Pt_after_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{xy}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_DCAxy_Pt_after_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_DCAxy_Pt_after_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_DCAz_Pt_after
        if (iParticleType==0){histoName="hTrack_DCAz_Pt_after_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_DCAz_Pt_after_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean DCA_{#it{z}} (cm)",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_DCAz_Pt_after_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_DCAz_Pt_after_PreSel_ProjPt[i]);
        //-------------------------------------------------------------------------------------------------------------------------------
        //Pre Selection: hTrack_NFindCls_Pt_TPC_after
        if (iParticleType==0){histoName="hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt";}
        if (iParticleType==1){histoName="";}
        if (i==0) vecHistosName.push_back(histoName);
        hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i]                 = new TH1D(Form("%s_%s", histoName.Data(), DataSets[i].Data()),Form("h%s; Run Number; Mean Findable Clusters after Cut",histoName.Data()),hNBin,hFBin,hLBin);
        EditTH1(globalRuns, doEquidistantXaxis, hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i], hMarkerStyle[i], hMarkerSize[i], hMarkerColor[i], hLineColor[i]);
        vecHistos[i].push_back(hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i]);
    } // end of loop over datasets
    
    //*****************************************************************************************************
    //******************************* create log files *****************************************************
    //*****************************************************************************************************
    fstream* fLog           = new fstream[nSets];
    fstream* fEventLog      = new fstream[nSets];
    for(Int_t iStr=0; iStr<nSets; iStr++){
        if (useDataRunListForMC && iStr>=nData) fLog[iStr].open(Form("%s/A-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
        else fLog[iStr].open(Form("%s/A-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << DataSets[iStr].Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << fCollisionSystem.Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[iStr] << "processed cut: " << fCutSelection.Data() << endl;
        fLog[iStr] << calo.Data() << endl;
        fLog[iStr] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        
        if (useDataRunListForMC && iStr>=nData) fEventLog[iStr].open(Form("%s/A-NEvents-%s-%s.log",outputDir.Data(),DataSets[iStr].Data(),DataSets[0].Data()), ios::out);
        else fEventLog[iStr].open(Form("%s/A-NEvents-%s.log",outputDir.Data(),DataSets[iStr].Data()), ios::out);
        fEventLog[iStr] << Form("Run\t%s",DataSets[iStr].Data()) << endl;
    }
    
    //*****************************************************************************************************
    //********************************** Supporting variables *********************************************
    //*****************************************************************************************************
    TString fRootFile;
    TString fDataSet;
    TString fRunNumber;
    Int_t bin                               = -1;
    TCanvas* cvsQuadratic   = new TCanvas("cvsQuadratic","",10,10,500,500);  // gives the page size
    
    std::vector<TString>* vecMissingRuns    = new std::vector<TString>[nSets];
    
    
    //*****************************************************************************************************
    //****************************** Looping over DataSets ************************************************
    //*****************************************************************************************************
    for(Int_t i=0; i<nSets; i++) {
        vecRuns.clear();
        fDataSet                      = vecDataSet.at(i);
        if (!readin(fileRuns[i], vecRuns)) {cout << "ERROR, no Run Numbers could be found! Returning..." << endl; return;}
        
        //*****************************************************************************************************
        //****************************** Looping over Runs ****************************************************
        //*****************************************************************************************************
        cout << endl;
        cout << "\t----------------------------------------------------------------------------" << endl;
        cout << "\tLooping over Runs of DataSet |" << (vecDataSet.at(i)).Data() << "|" << endl;
        cout << "\t----------------------------------------------------------------------------" << endl;
        cout << endl;
        fLog[i] << "Looping over Runs:" << endl;
        
        vecMissingRuns[i].clear();
        for(Int_t j=0; j<(Int_t) vecRuns.size(); j++){
            //--------------------------------------------------------------------------------------------------------
            //------------------------- Read in individual files -----------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            fRunNumber                  = vecRuns.at(j);
            fRootFile                   = Form("%s/%s/%s/%s", filePath.Data(), fDataSet.Data(), fRunNumber.Data(), fileName.Data());
            TFile* RootFile             = new TFile(fRootFile.Data(),"READ");
            if (RootFile->IsZombie()) {vecMissingRuns[i].push_back(fRunNumber); cout << "INFO: ROOT file '" << fRootFile.Data() << "' could not be openend, continue!" << endl; continue;}
            
            cout << endl;
            cout << "\t\t----------------------------------------------------------------------------" << endl;
            cout << Form("\t\tRun %s", fRunNumber.Data()) << endl;
            cout << "\t\tProcessing file: " << fRootFile.Data() << endl;
            cout << "\t\t----------------------------------------------------------------------------" << endl;
            cout << endl;
            
            // reading respective containers
            //TopDir; There should only be one; It is the main directory
            TList* TopDir                   = (TList*) RootFile->Get(nameMainDir.Data());
            if (TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
            else TopDir->SetOwner(kTRUE);
            //TopContainer; The directory with the chosen cut number
            TList* TopContainer             = (TList*) TopDir->FindObject(Form("Cut Number %s",fCutSelection.Data()));
            if (TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",fCutSelection.Data()) << " not found in File" << endl; return;}
            else TopContainer->SetOwner(kTRUE);
            //ESDContainer; ESD histogram directory in the "Cut Number" directory with the chosen cut number
            TList* ESDContainer             = (TList*) TopContainer->FindObject(Form("%s ESD histograms",fCutSelection.Data()));
            if (ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",fCutSelection.Data()) << " not found in File" << endl; return;}
            else ESDContainer->SetOwner(kTRUE);
            //MCContainer not needed; MC histograms directory in the "Cut Number" directory with the chosen cut number
            TList* MCContainer            = (TList*) TopContainer->FindObject(Form("%s MC histograms",fCutSelection.Data()));
            if (MCContainer == NULL) {cout << "INFO: " << Form("%s MC histograms",fCutSelection.Data()) << " not found in File, processing data?" << endl;}
            else {
                MCContainer->SetOwner(kTRUE);
                isMC=kTRUE;
            }
            //TrueContainer; True histograms  directory in the "Cut Number" directory with the chosen cut number, only if file is MC
            TList* TrueContainer            = (TList*) TopContainer->FindObject(Form("%s True histograms",fCutSelection.Data()));
            if (TrueContainer == NULL) {cout << "INFO: " << Form("%s True histograms",fCutSelection.Data()) << " not found in File, processing data?" << endl;}
            else TrueContainer->SetOwner(kTRUE);
            //ConvEventCutsContainer; ConvEventCuts directory in the "Cut Number" directory with the chosen cut number
            TList* ConvEventCutsContainer   = (TList*) TopContainer->FindObject(Form("ConvEventCuts_%s",fEventCutSelection.Data()));
            if (ConvEventCutsContainer == NULL) {cout << "ERROR: " << Form("ConvEventCuts_%s",fEventCutSelection.Data()) << " not found in File" << endl; return;}
            else if (ConvEventCutsContainer) ConvEventCutsContainer->SetOwner(kTRUE);
            //PionCutsContainer; PionCuts directory in the "Cut Number" directory with the chosen cut number
            TString fPionCutsContainerCutString=Form("%s_%s_%s_%s_%s",fEventCutSelection.Data(),fGammaCutSelection.Data(),fNeutralPionCutSelection.Data(),fPionCutSelection.Data(),fMesonCutSelection.Data());
            TList* PionCutsContainer       = (TList*) TopContainer->FindObject(Form("PionCuts_%s",fPionCutsContainerCutString.Data()));
            if (PionCutsContainer == NULL) {
                cout << "INFO: " << Form("PionCuts_%s",fPionCutsContainerCutString.Data()) << " not found in File; Trying other Cut String" << endl;
                fPionCutsContainerCutString=Form("%s_%s_%s_%s_%s_%s",fEventCutSelection.Data(),fGammaCutSelection.Data(),fClusterCutSelection.Data(),fNeutralPionCutSelection.Data(),fPionCutSelection.Data(),fMesonCutSelection.Data());
                PionCutsContainer       = (TList*) TopContainer->FindObject(Form("PionCuts_%s",fPionCutsContainerCutString.Data()));
                if (PionCutsContainer == NULL) {
                    cout << "ERROR: " << Form("PionCuts_%s",fPionCutsContainerCutString.Data()) << " not found in File" << endl; return;
                }
            }
            if (PionCutsContainer != NULL) PionCutsContainer->SetOwner(kTRUE);
            //ConvCutsContainer; ConvCuts directory in the "Cut Number" directory with the chosen cut number
            TList* ConvCutsContainer        = (TList*) TopContainer->FindObject(Form("ConvCuts_%s",fGammaCutSelection.Data()));
            if (isConv && ConvCutsContainer == NULL) {cout << "ERROR: " << Form("ConvCuts_%s",fGammaCutSelection.Data()) << " not found in File" << endl; return;}
            else if (ConvCutsContainer) ConvCutsContainer->SetOwner(kTRUE);
            //MesonCutsContainer; ConvMesonCuts directory in the main directory
            TList* MesonCutsContainer       = (TList*) TopContainer->FindObject(Form("ConvMesonCuts_%s",fMesonCutSelection.Data()));
            if (MesonCutsContainer == NULL) {cout << "ERROR: " << Form("ConvMesonCuts_%s",fMesonCutSelection.Data()) << " not found in File" << endl; return;}
            else if (MesonCutsContainer) MesonCutsContainer->SetOwner(kTRUE);
            //MesonCutsContainer2; ConvMesonCuts directory in the main directory
            //CaloCutsContainer; CaloCuts directory in the main directory
            TList* CaloCutsContainer        = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fClusterCutSelection.Data()));
            if (isCalo && CaloCutsContainer == NULL) {cout << "ERROR: " << Form("CaloCuts_%s",fClusterCutSelection.Data()) << " not found in File" << endl; return;}
            else if (CaloCutsContainer) CaloCutsContainer->SetOwner(kTRUE);
            TList* TopContainerGamma        = NULL;
            TString ContainerGammaCut       = "";
            //PionCuts2: Main directory pion cuts
            //TString fPionCuts2ContainerCutString=Form("%s_%s_%s_%s_%s",fEventCutSelection.Data(),fGammaCutSelection.Data(),fNeutralPionCutSelection.Data(),fPionCutSelection.Data(),fMesonCutSelection.Data());
            //fPionCuts2ContainerCutString=Form("000000200");
            TList* PionCuts2Container       = (TList*) TopDir->FindObject(Form("PionCuts_%s",fPionCuts2ContainerCutString.Data()));
            if (PionCuts2Container == NULL) {cout << "ERROR: " << Form("PionCuts_%s",fPionCuts2ContainerCutString.Data()) << " not found in File" << endl; return;}
            else if (PionCuts2Container) PionCuts2Container->SetOwner(kTRUE);
            if (isConv){
                TList *listCuts=NULL; TString name=""; Int_t k=0;
                do{ listCuts = (TList*)TopDir->At(k); name = listCuts->GetName(); } while(!name.BeginsWith("ConvCuts_") && ++k<TopDir->GetSize());
                if (k>=TopDir->GetSize()) TopContainerGamma = NULL;
                else TopContainerGamma  = listCuts;
                if (TopContainerGamma == NULL) {cout << "ERROR: " << "ConvCuts_*" << " not found in File" << endl; return;}
                else TopContainerGamma->SetOwner(kTRUE);
                ContainerGammaCut = TopContainerGamma->GetName();
                ContainerGammaCut.Remove(0,9);
            }
            TList* TopContainerEvent = NULL;
            TString ContainerEventCut = "";
            TList* listCuts=NULL; TString name=""; Int_t k=0;
            do{ listCuts = (TList*)TopDir->At(k); name = listCuts->GetName(); } while(!name.BeginsWith("ConvEventCuts_") && ++k<TopDir->GetSize());
            if (k>=TopDir->GetSize()) TopContainerEvent = NULL;
            else TopContainerEvent = listCuts;
            
            if (TopContainerEvent == NULL) {cout << "ERROR: " << "ConvEventCuts_*" << " not found in File" << endl; return;}
            else TopContainerEvent->SetOwner(kTRUE);
            ContainerEventCut = TopContainerEvent->GetName();
            ContainerEventCut.Remove(0,14);
            //--------------------------------------------------------------------------------------------------------
            if (doEquidistantXaxis) bin  = mapBin[fRunNumber];
            else bin = fRunNumber.Atoi() - hFBin;
            
            //Get Histograms From Container
            //-------------------------------------------------------------------------------------------------------------------------------
            //-------------------------------------------|Get Histograms: ESD-Histograms|----------------------------------------------------
            //-------------------------------------------------------------------------------------------------------------------------------
            if (ESDContainer!=NULL){
                cout<<"ESDContainer"<<endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_PrimaryNegPions_Pt
                if (iParticleType==0){histoName="ESD_PrimaryNegPions_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_PrimaryNegPions_Pt = (TH1D*)ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryNegPions_Pt)&&(fHistESD_PrimaryNegPions_Pt->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryNegPions_Pt__temp = new TH1D(*fHistESD_PrimaryNegPions_Pt);
                    hESD_PrimaryNegPions_Pt[i]->SetBinContent(bin, fHistESD_PrimaryNegPions_Pt__temp->GetMean());
                    hESD_PrimaryNegPions_Pt[i]->SetBinError(bin, fHistESD_PrimaryNegPions_Pt__temp->GetMeanError());
                    fHistESD_PrimaryNegPions_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})");
                    fHistESD_PrimaryNegPions_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryNegPions_Pt__temp->Sumw2();
                    vecESD_PrimaryNegPions_Pt[i].push_back(fHistESD_PrimaryNegPions_Pt__temp);
                }//else cout << "INFO: Object |fHistESD_PrimaryNegPions_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPosPions_Pt
                if (iParticleType==0){histoName="ESD_PrimaryPosPions_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_PrimaryPosPions_Pt = (TH1D*)ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPosPions_Pt)&&(fHistESD_PrimaryPosPions_Pt->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPosPions_Pt__temp = new TH1D(*fHistESD_PrimaryPosPions_Pt);
                    hESD_PrimaryPosPions_Pt[i]->SetBinContent(bin, fHistESD_PrimaryPosPions_Pt__temp->GetMean());
                    hESD_PrimaryPosPions_Pt[i]->SetBinError(bin, fHistESD_PrimaryPosPions_Pt__temp->GetMeanError());
                    fHistESD_PrimaryPosPions_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #it{#pi^{+}}} (GeV/#it{c})");
                    fHistESD_PrimaryPosPions_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPosPions_Pt__temp->Sumw2();
                    vecESD_PrimaryPosPions_Pt[i].push_back(fHistESD_PrimaryPosPions_Pt__temp);
                }//else cout << "INFO: Object |fHistESD_PrimaryPosPions_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_PrimaryNegPions_Phi
                if (iParticleType==0){histoName="ESD_PrimaryNegPions_Phi";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_PrimaryNegPions_Phi = (TH1D*)ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryNegPions_Phi)&&(fHistESD_PrimaryNegPions_Phi->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryNegPions_Phi__temp = new TH1D(*fHistESD_PrimaryNegPions_Phi);
                    hESD_PrimaryNegPions_Phi[i]->SetBinContent(bin, fHistESD_PrimaryNegPions_Phi__temp->GetMean());
                    hESD_PrimaryNegPions_Phi[i]->SetBinError(bin, fHistESD_PrimaryNegPions_Phi__temp->GetMeanError());
                    fHistESD_PrimaryNegPions_Phi__temp->GetXaxis()->SetTitle("#it{#varphi}_{#it{#pi^{-}}}");
                    fHistESD_PrimaryNegPions_Phi__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryNegPions_Phi__temp->Sumw2();
                    vecESD_PrimaryNegPions_Phi[i].push_back(fHistESD_PrimaryNegPions_Phi__temp);
                }//else cout << "INFO: Object |fHistESD_PrimaryNegPions_Phi| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPosPions_Phi
                if (iParticleType==0){histoName="ESD_PrimaryPosPions_Phi";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_PrimaryPosPions_Phi = (TH1D*)ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPosPions_Phi)&&(fHistESD_PrimaryPosPions_Phi->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPosPions_Phi__temp = new TH1D(*fHistESD_PrimaryPosPions_Phi);
                    hESD_PrimaryPosPions_Phi[i]->SetBinContent(bin, fHistESD_PrimaryPosPions_Phi__temp->GetMean());
                    hESD_PrimaryPosPions_Phi[i]->SetBinError(bin, fHistESD_PrimaryPosPions_Phi__temp->GetMeanError());
                    fHistESD_PrimaryPosPions_Phi__temp->GetXaxis()->SetTitle("#it{#varphi}_{#it{#pi^{+}}}");
                    fHistESD_PrimaryPosPions_Phi__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPosPions_Phi__temp->Sumw2();
                    vecESD_PrimaryPosPions_Phi[i].push_back(fHistESD_PrimaryPosPions_Phi__temp);
                }//else cout << "INFO: Object |fHistESD_PrimaryPosPions_Phi| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_PrimaryNegPions_Eta
                if (iParticleType==0){histoName="ESD_PrimaryNegPions_Eta";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_PrimaryNegPions_Eta = (TH1D*)ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryNegPions_Eta)&&(fHistESD_PrimaryNegPions_Eta->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryNegPions_Eta__temp = new TH1D(*fHistESD_PrimaryNegPions_Eta);
                    hESD_PrimaryNegPions_Eta[i]->SetBinContent(bin, fHistESD_PrimaryNegPions_Eta__temp->GetMean());
                    hESD_PrimaryNegPions_Eta[i]->SetBinError(bin, fHistESD_PrimaryNegPions_Eta__temp->GetMeanError());
                    fHistESD_PrimaryNegPions_Eta__temp->GetXaxis()->SetTitle("#it{#eta}_{#it{#pi^{-}}}");
                    fHistESD_PrimaryNegPions_Eta__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryNegPions_Eta__temp->Sumw2();
                    vecESD_PrimaryNegPions_Eta[i].push_back(fHistESD_PrimaryNegPions_Eta__temp);
                }//else cout << "INFO: Object |fHistESD_PrimaryNegPions_Eta| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPosPions_Eta
                if (iParticleType==0){histoName="ESD_PrimaryPosPions_Eta";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_PrimaryPosPions_Eta = (TH1D*)ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPosPions_Eta)&&(fHistESD_PrimaryPosPions_Eta->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPosPions_Eta__temp = new TH1D(*fHistESD_PrimaryPosPions_Eta);
                    hESD_PrimaryPosPions_Eta[i]->SetBinContent(bin, fHistESD_PrimaryPosPions_Eta__temp->GetMean());
                    hESD_PrimaryPosPions_Eta[i]->SetBinError(bin, fHistESD_PrimaryPosPions_Eta__temp->GetMeanError());
                    fHistESD_PrimaryPosPions_Eta__temp->GetXaxis()->SetTitle("#it{#eta}_{#it{#pi^{+}}}");
                    fHistESD_PrimaryPosPions_Eta__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPosPions_Eta__temp->Sumw2();
                    vecESD_PrimaryPosPions_Eta[i].push_back(fHistESD_PrimaryPosPions_Eta__temp);
                }//else cout << "INFO: Object |fHistESD_PrimaryPosPions_Eta| could not be found!" << endl;
            } else {cout<<"no ESDContainer"<<endl;}
            //*******************************************************************************************************************************
            //-------------------------------------------------------------------------------------------------------------------------------
            //*******************************************************************************************************************************
            if (PionCutsContainer!=NULL){
                //AfterQA=Cut Number Pion Cut (Pion Cut folder in Cut number directory); Pre Selection=Main Dir Pion Cut
                //AfterQA
                cout<<"PionCutsContainer"<<endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                // PionCutsContainer: IsPionSelected
                if (iParticleType==0){histoName="IsPionSelected";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistIsPionSelected_AfterQA = (TH1D*)PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHistIsPionSelected_AfterQA)&&(fHistIsPionSelected_AfterQA->GetEntries()>0)){
                    TH1D* fHistIsPionSelected_AfterQA__temp = new TH1D(*fHistIsPionSelected_AfterQA);
                    hIsPionSelected_AfterQA[i]->SetBinContent(bin, (fHistIsPionSelected_AfterQA__temp->GetBinContent(5)/fHistIsPionSelected_AfterQA__temp->GetBinContent(1)));
                    hIsPionSelected_AfterQA[i]->SetBinError(bin, 0);
                    fHistIsPionSelected_AfterQA__temp->GetXaxis()->SetTitle("");
                    fHistIsPionSelected_AfterQA__temp->GetYaxis()->SetTitle("# Entries");
                    fHistIsPionSelected_AfterQA__temp->Sumw2();
                    vecIsPionSelected_AfterQA[i].push_back(fHistIsPionSelected_AfterQA__temp);
                }//else cout << "INFO: Object |fHistIsPionSelected_AfterQA| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                // PionCutsContainer: dEdxCuts
                if (iParticleType==0){histoName="dEdxCuts";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistdEdxCuts_AfterQA = (TH1D*)PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHistdEdxCuts_AfterQA)&&(fHistdEdxCuts_AfterQA->GetEntries()>0)){
                    TH1D* fHistdEdxCuts_AfterQA__temp = new TH1D(*fHistdEdxCuts_AfterQA);
                    hdEdxCuts_AfterQA[i]->SetBinContent(bin, (fHistdEdxCuts_AfterQA__temp->GetBinContent(5)/fHistdEdxCuts_AfterQA__temp->GetBinContent(1)));
                    hdEdxCuts_AfterQA[i]->SetBinError(bin,0);
                    fHistdEdxCuts_AfterQA__temp->GetXaxis()->SetTitle("");
                    fHistdEdxCuts_AfterQA__temp->GetYaxis()->SetTitle("# Entries");
                    fHistdEdxCuts_AfterQA__temp->Sumw2();
                    vecdEdxCuts_AfterQA[i].push_back(fHistdEdxCuts_AfterQA__temp);
                }//else cout << "INFO: Object |fHistdEdxCuts_AfterQA| could not be found!" << endl;
            }else {cout<<"no PionCutsContainer"<<endl;}
            //*******************************************************************************************************************************
            //-------------------------------------------------------------------------------------------------------------------------------
            //*******************************************************************************************************************************
            if (PionCuts2Container!=NULL){
                //Pre Selection (Pion Cut folder in Main Directory)
                cout<<"PionCuts2Container"<<endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //IsPionSelected_PreSel
                if (iParticleType==0){histoName="IsPionSelected";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistIsPionSelected_PreSel = (TH1D*)PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistIsPionSelected_PreSel)&&(fHistIsPionSelected_PreSel->GetEntries()>0)){
                    TH1D* fHistIsPionSelected_PreSel__temp = new TH1D(*fHistIsPionSelected_PreSel);
                    hIsPionSelected_PreSel[i]->SetBinContent(bin, (fHistIsPionSelected_PreSel__temp->GetBinContent(5)/fHistIsPionSelected_PreSel__temp->GetBinContent(1)));
                    //cout<<fHistIsPionSelected_PreSel__temp->GetBinContent(1)<<endl;
                    //cout<<fHistIsPionSelected_PreSel__temp->GetBinContent(5)<<endl;
                    hIsPionSelected_PreSel[i]->SetBinError(bin, 0);
                    fHistIsPionSelected_PreSel__temp->GetXaxis()->SetTitle("");
                    fHistIsPionSelected_PreSel__temp->GetYaxis()->SetTitle("# Entries");
                    fHistIsPionSelected_PreSel__temp->Sumw2();
                    vecIsPionSelected_PreSel[i].push_back(fHistIsPionSelected_PreSel__temp);
                }//else cout << "INFO: Object |fHistIsPionSelected_PreSel| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //dEdxCuts_PreSel
                if (iParticleType==0){histoName="dEdxCuts";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistdEdxCuts_PreSel = (TH1D*)PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistdEdxCuts_PreSel)&&(fHistdEdxCuts_PreSel->GetEntries()>0)){
                    TH1D* fHistdEdxCuts_PreSel__temp = new TH1D(*fHistdEdxCuts_PreSel);
                    hdEdxCuts_PreSel[i]->SetBinContent(bin, (fHistdEdxCuts_PreSel__temp->GetBinContent(5)/fHistdEdxCuts_PreSel__temp->GetBinContent(1)));
                    hdEdxCuts_PreSel[i]->SetBinError(bin, 0);
                    fHistdEdxCuts_PreSel__temp->GetXaxis()->SetTitle("");
                    fHistdEdxCuts_PreSel__temp->GetYaxis()->SetTitle("# Entries");
                    //fHistdEdxCuts_PreSel__temp->Sumw2();
                    vecdEdxCuts_PreSel[i].push_back(fHistdEdxCuts_PreSel__temp);
                }//else cout << "INFO: Object |fHistdEdxCuts_PreSel| could not be found!" << endl;
            }else {cout<<"no PionCuts2Container"<<endl;}
            //*******************************************************************************************************************************
            //-------------------------------------------------------------------------------------------------------------------------------
            //*******************************************************************************************************************************
            if (MCContainer!=NULL){
                //MC histograms
                cout<<"MCContainer"<<endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //MC_AllPosPions_Pt
                if (iParticleType==0){histoName="MC_AllPosPions_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistMC_AllPosPions_Pt = (TH1D*)MCContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistMC_AllPosPions_Pt)&&(fHistMC_AllPosPions_Pt->GetEntries()>0)){
                    TH1D* fHistMC_AllPosPions_Pt__temp = new TH1D(*fHistMC_AllPosPions_Pt);
                    hMC_AllPosPions_Pt[i]->SetBinContent(bin, fHistMC_AllPosPions_Pt__temp->GetMean());
                    hMC_AllPosPions_Pt[i]->SetBinError(bin, fHistMC_AllPosPions_Pt__temp->GetMeanError());
                    fHistMC_AllPosPions_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #pi^{+}} (GeV/#it{c})");
                    fHistMC_AllPosPions_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistMC_AllPosPions_Pt__temp->Sumw2();
                    vecMC_AllPosPions_Pt[i].push_back(fHistMC_AllPosPions_Pt__temp);
                }//else cout << "INFO: Object |fHistMC_AllPosPions_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //MC_AllNegPions_Pt
                if (iParticleType==0){histoName="MC_AllNegPions_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistMC_AllNegPions_Pt = (TH1D*)MCContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistMC_AllNegPions_Pt)&&(fHistMC_AllNegPions_Pt->GetEntries()>0)){
                    TH1D* fHistMC_AllNegPions_Pt__temp = new TH1D(*fHistMC_AllNegPions_Pt);
                    hMC_AllNegPions_Pt[i]->SetBinContent(bin, fHistMC_AllNegPions_Pt__temp->GetMean());
                    hMC_AllNegPions_Pt[i]->SetBinError(bin, fHistMC_AllNegPions_Pt__temp->GetMeanError());
                    fHistMC_AllNegPions_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #pi^{-}} (GeV/#it{c})");
                    fHistMC_AllNegPions_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistMC_AllNegPions_Pt__temp->Sumw2();
                    vecMC_AllNegPions_Pt[i].push_back(fHistMC_AllNegPions_Pt__temp);
                }//else cout << "INFO: Object |fHistMC_AllNegPions_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //MC_PosPionsFromNeutralMeson_Pt
                if (iParticleType==0){histoName="MC_PosPionsFromNeutralMeson_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistMC_PosPionsFromNeutralMeson_Pt = (TH1D*)MCContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistMC_PosPionsFromNeutralMeson_Pt)&&(fHistMC_PosPionsFromNeutralMeson_Pt->GetEntries()>0)){
                    TH1D* fHistMC_PosPionsFromNeutralMeson_Pt__temp = new TH1D(*fHistMC_PosPionsFromNeutralMeson_Pt);
                    hMC_PosPionsFromNeutralMeson_Pt[i]->SetBinContent(bin, fHistMC_PosPionsFromNeutralMeson_Pt__temp->GetMean());
                    hMC_PosPionsFromNeutralMeson_Pt[i]->SetBinError(bin, fHistMC_PosPionsFromNeutralMeson_Pt__temp->GetMeanError());
                    fHistMC_PosPionsFromNeutralMeson_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #pi^{+}} (GeV/#it{c})");
                    fHistMC_PosPionsFromNeutralMeson_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistMC_PosPionsFromNeutralMeson_Pt__temp->Sumw2();
                    vecMC_PosPionsFromNeutralMeson_Pt[i].push_back(fHistMC_PosPionsFromNeutralMeson_Pt__temp);
                }//else cout << "INFO: Object |fHistMC_PosPionsFromNeutralMeson_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //MC_NegPionsFromNeutralMeson_Pt
                if (iParticleType==0){histoName="MC_NegPionsFromNeutralMeson_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistMC_NegPionsFromNeutralMeson_Pt = (TH1D*)MCContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistMC_NegPionsFromNeutralMeson_Pt)&&(fHistMC_NegPionsFromNeutralMeson_Pt->GetEntries()>0)){
                    TH1D* fHistMC_NegPionsFromNeutralMeson_Pt__temp = new TH1D(*fHistMC_NegPionsFromNeutralMeson_Pt);
                    hMC_NegPionsFromNeutralMeson_Pt[i]->SetBinContent(bin, fHistMC_NegPionsFromNeutralMeson_Pt__temp->GetMean());
                    hMC_NegPionsFromNeutralMeson_Pt[i]->SetBinError(bin, fHistMC_NegPionsFromNeutralMeson_Pt__temp->GetMeanError());
                    fHistMC_NegPionsFromNeutralMeson_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #pi^{-}} (GeV/#it{c})");
                    fHistMC_NegPionsFromNeutralMeson_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistMC_NegPionsFromNeutralMeson_Pt__temp->Sumw2();
                    vecMC_NegPionsFromNeutralMeson_Pt[i].push_back(fHistMC_NegPionsFromNeutralMeson_Pt__temp);
                }//else cout << "INFO: Object |fHistMC_NegPionsFromNeutralMeson_Pt| could not be found!" << endl;
            }else {cout<<"no MCContainer"<<endl;}
            //*******************************************************************************************************************************
            //-------------------------------------------------------------------------------------------------------------------------------
            //*******************************************************************************************************************************
            if (TrueContainer != NULL){
                //True histograms
                cout<<"TrueContainer"<<endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_TruePosPion_Pt
                if (iParticleType==0){histoName="ESD_TruePosPion_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_TruePosPion_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_TruePosPion_Pt)&&(fHistESD_TruePosPion_Pt->GetEntries()>0)){
                    TH1D* fHistESD_TruePosPion_Pt__temp = new TH1D(*fHistESD_TruePosPion_Pt);
                    hESD_TruePosPion_Pt[i]->SetBinContent(bin, fHistESD_TruePosPion_Pt__temp->GetMean());
                    hESD_TruePosPion_Pt[i]->SetBinError(bin, fHistESD_TruePosPion_Pt__temp->GetMeanError());
                    fHistESD_TruePosPion_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #pi^{+}} (GeV/#it{c})");
                    fHistESD_TruePosPion_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_TruePosPion_Pt__temp->Sumw2();
                    vecESD_TruePosPion_Pt[i].push_back(fHistESD_TruePosPion_Pt__temp);
                }//else cout << "INFO: Object |fHistESD_TruePosPion_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_TrueNegPion_Pt
                if (iParticleType==0){histoName="ESD_TrueNegPion_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_TrueNegPion_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_TrueNegPion_Pt)&&(fHistESD_TrueNegPion_Pt->GetEntries()>0)){
                    TH1D* fHistESD_TrueNegPion_Pt__temp = new TH1D(*fHistESD_TrueNegPion_Pt);
                    hESD_TrueNegPion_Pt[i]->SetBinContent(bin, fHistESD_TrueNegPion_Pt__temp->GetMean());
                    hESD_TrueNegPion_Pt[i]->SetBinError(bin, fHistESD_TrueNegPion_Pt__temp->GetMeanError());
                    fHistESD_TrueNegPion_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})");
                    fHistESD_TrueNegPion_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_TrueNegPion_Pt__temp->Sumw2();
                    vecESD_TrueNegPion_Pt[i].push_back(fHistESD_TrueNegPion_Pt__temp);
                }//else cout << "INFO: Object |fHistESD_TrueNegPion_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_TruePosPionFromNeutralMeson_Pt
                if (iParticleType==0){histoName="ESD_TruePosPionFromNeutralMeson_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_TruePosPionFromNeutralMeson_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_TruePosPionFromNeutralMeson_Pt)&&(fHistESD_TruePosPionFromNeutralMeson_Pt->GetEntries()>0)){
                    TH1D* fHistESD_TruePosPionFromNeutralMeson_Pt__temp = new TH1D(*fHistESD_TruePosPionFromNeutralMeson_Pt);
                    hESD_TruePosPionFromNeutralMeson_Pt[i]->SetBinContent(bin, fHistESD_TruePosPionFromNeutralMeson_Pt__temp->GetMean());
                    hESD_TruePosPionFromNeutralMeson_Pt[i]->SetBinError(bin, fHistESD_TruePosPionFromNeutralMeson_Pt__temp->GetMeanError());
                    fHistESD_TruePosPionFromNeutralMeson_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #it{#pi^{+}}} (GeV/#it{c})");
                    fHistESD_TruePosPionFromNeutralMeson_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_TruePosPionFromNeutralMeson_Pt__temp->Sumw2();
                    vecESD_TruePosPionFromNeutralMeson_Pt[i].push_back(fHistESD_TruePosPionFromNeutralMeson_Pt__temp);
                }//else cout << "INFO: Object |fHistESD_TruePosPionFromNeutralMeson_Pt| could not be found!" << endl;
                //--------------------------------------------------------------------------------------------------------
                //ESD_TrueNegPionFromNeutralMeson_Pt
                if (iParticleType==0){histoName="ESD_TrueNegPionFromNeutralMeson_Pt";}
                if (iParticleType==1){histoName="";}
                TH1D* fHistESD_TrueNegPionFromNeutralMeson_Pt = (TH1D*)TrueContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_TrueNegPionFromNeutralMeson_Pt)&&(fHistESD_TrueNegPionFromNeutralMeson_Pt->GetEntries()>0)){
                    TH1D* fHistESD_TrueNegPionFromNeutralMeson_Pt__temp = new TH1D(*fHistESD_TrueNegPionFromNeutralMeson_Pt);
                    hESD_TrueNegPionFromNeutralMeson_Pt[i]->SetBinContent(bin, fHistESD_TrueNegPionFromNeutralMeson_Pt__temp->GetMean());
                    hESD_TrueNegPionFromNeutralMeson_Pt[i]->SetBinError(bin, fHistESD_TrueNegPionFromNeutralMeson_Pt__temp->GetMeanError());
                    fHistESD_TrueNegPionFromNeutralMeson_Pt__temp->GetXaxis()->SetTitle("#it{p}_{T, #it{#pi^{-}}} (GeV/#it{c})");
                    fHistESD_TrueNegPionFromNeutralMeson_Pt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_TrueNegPionFromNeutralMeson_Pt__temp->Sumw2();
                    vecESD_TrueNegPionFromNeutralMeson_Pt[i].push_back(fHistESD_TrueNegPionFromNeutralMeson_Pt__temp);
                }//else cout << "INFO: Object |fHistESD_TrueNegPionFromNeutralMeson_Pt| could not be found!" << endl;
            }else {cout<<"no TrueContainer"<<endl;}
            //*******************************************************************************************************************************
            //--------------------------------------------------------Projections------------------------------------------------------------
            //*******************************************************************************************************************************
            if (ESDContainer!=NULL){
                cout<<"ESDContainer Projections"<<endl;
                //ESD_PrimaryNegPions_ClsTPC
                if (iParticleType==0){histoName="ESD_PrimaryNegPions_ClsTPC";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistESD_PrimaryNegPions_ClsTPC = (TH2D*) ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryNegPions_ClsTPC)&&(fHistESD_PrimaryNegPions_ClsTPC->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryNegPions_ClsTPC__temp = (TH1D*) fHistESD_PrimaryNegPions_ClsTPC->ProjectionX(Form("%s",histoName.Data()),fHistESD_PrimaryNegPions_ClsTPC->GetYaxis()->GetFirst(),fHistESD_PrimaryNegPions_ClsTPC->GetYaxis()->GetLast());
                    hESD_PrimaryNegPions_ClsTPC_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryNegPions_ClsTPC__temp->GetMean());
                    hESD_PrimaryNegPions_ClsTPC_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryNegPions_ClsTPC__temp->GetMeanError());
                    fHistESD_PrimaryNegPions_ClsTPC__temp->GetXaxis()->SetTitle("Findable Clusters");
                    fHistESD_PrimaryNegPions_ClsTPC__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryNegPions_ClsTPC__temp->Sumw2();
                    vecESD_PrimaryNegPions_ClsTPC_ProjPt[i].push_back(fHistESD_PrimaryNegPions_ClsTPC__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPosPions_ClsTPC
                if (iParticleType==0){histoName="ESD_PrimaryPosPions_ClsTPC";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistESD_PrimaryPosPions_ClsTPC = (TH2D*) ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPosPions_ClsTPC)&&(fHistESD_PrimaryPosPions_ClsTPC->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPosPions_ClsTPC__temp = (TH1D*) fHistESD_PrimaryPosPions_ClsTPC->ProjectionX(Form("%s",histoName.Data()),fHistESD_PrimaryPosPions_ClsTPC->GetYaxis()->GetFirst(),fHistESD_PrimaryPosPions_ClsTPC->GetYaxis()->GetLast());
                    hESD_PrimaryPosPions_ClsTPC_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPosPions_ClsTPC__temp->GetMean());
                    hESD_PrimaryPosPions_ClsTPC_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPosPions_ClsTPC__temp->GetMeanError());
                    fHistESD_PrimaryPosPions_ClsTPC__temp->GetXaxis()->SetTitle("Findable Clusters");
                    fHistESD_PrimaryPosPions_ClsTPC__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPosPions_ClsTPC__temp->Sumw2();
                    vecESD_PrimaryPosPions_ClsTPC_ProjPt[i].push_back(fHistESD_PrimaryPosPions_ClsTPC__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_DCAxy
                if (iParticleType==0){histoName="ESD_PrimaryPions_DCAxy";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistESD_PrimaryPions_DCAxy = (TH2D*) ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPions_DCAxy)&&(fHistESD_PrimaryPions_DCAxy->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPions_DCAxy__temp = (TH1D*) fHistESD_PrimaryPions_DCAxy->ProjectionX(Form("%s",histoName.Data()),fHistESD_PrimaryPions_DCAxy->GetYaxis()->GetFirst(),fHistESD_PrimaryPions_DCAxy->GetYaxis()->GetLast());
                    hESD_PrimaryPions_DCAxy_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_DCAxy__temp->GetMean());
                    hESD_PrimaryPions_DCAxy_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_DCAxy__temp->GetMeanError());
                    fHistESD_PrimaryPions_DCAxy__temp->GetXaxis()->SetTitle("DCA_{#it{xy}} (cm)");
                    fHistESD_PrimaryPions_DCAxy__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_DCAxy__temp->Sumw2();
                    vecESD_PrimaryPions_DCAxy_ProjPt[i].push_back(fHistESD_PrimaryPions_DCAxy__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_DCAz
                if (iParticleType==0){histoName="ESD_PrimaryPions_DCAz";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistESD_PrimaryPions_DCAz = (TH2D*) ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPions_DCAz)&&(fHistESD_PrimaryPions_DCAz->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPions_DCAz__temp = (TH1D*) fHistESD_PrimaryPions_DCAz->ProjectionX(Form("%s",histoName.Data()),fHistESD_PrimaryPions_DCAz->GetYaxis()->GetFirst(),fHistESD_PrimaryPions_DCAz->GetYaxis()->GetLast());
                    hESD_PrimaryPions_DCAz_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_DCAz__temp->GetMean());
                    hESD_PrimaryPions_DCAz_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_DCAz__temp->GetMeanError());
                    fHistESD_PrimaryPions_DCAz__temp->GetXaxis()->SetTitle("DCA_{#it{z}} (cm)");
                    fHistESD_PrimaryPions_DCAz__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_DCAz__temp->Sumw2();
                    vecESD_PrimaryPions_DCAz_ProjPt[i].push_back(fHistESD_PrimaryPions_DCAz__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdx
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistESD_PrimaryPions_TPCdEdx = (TH2D*) ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPions_TPCdEdx)&&(fHistESD_PrimaryPions_TPCdEdx->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPions_TPCdEdx__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdx->ProjectionY(Form("%s",histoName.Data()),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetLast());
                    hESD_PrimaryPions_TPCdEdx_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdx__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdx_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdx__temp->GetMeanError());
                    fHistESD_PrimaryPions_TPCdEdx__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdx__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdx__temp->Sumw2();
                    vecESD_PrimaryPions_TPCdEdx_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdx__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdx_LowPt
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistESD_PrimaryPions_TPCdEdx,fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetLast());
                if ((fHistESD_PrimaryPions_TPCdEdx)&&(fHistESD_PrimaryPions_TPCdEdx->GetEntries()>0)){
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdx,0,fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->FindBin(0.2));
                    //SetXRange(fHistESD_PrimaryPions_TPCdEdx,minB-1,maxB+1);
                    GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdx,minYB,maxYB);
                    SetYRange(fHistESD_PrimaryPions_TPCdEdx,minYB-1,maxYB+1);
                    TH1D* fHistESD_PrimaryPions_TPCdEdx_LowPt__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdx->ProjectionY(Form("%s",histoName.Data()));
                    hESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdx_LowPt__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdx_LowPt__temp->GetMeanError());
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_LowPt__temp,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdx_LowPt__temp,minB-1,maxB+1);
                    fHistESD_PrimaryPions_TPCdEdx_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdx_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdx_LowPt__temp->Sumw2();
                    vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdx_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdx_MidPt
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistESD_PrimaryPions_TPCdEdx,fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetLast());
                if ((fHistESD_PrimaryPions_TPCdEdx)&&(fHistESD_PrimaryPions_TPCdEdx->GetEntries()>0)){
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdx,fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->FindBin(2),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdx,minYB,maxYB);
                    SetYRange(fHistESD_PrimaryPions_TPCdEdx,minYB-1,maxYB+1);
                    TH1D* fHistESD_PrimaryPions_TPCdEdx_MidPt__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdx->ProjectionY(Form("%s",histoName.Data()));
                    hESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdx_MidPt__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdx_MidPt__temp->GetMeanError());
                    fHistESD_PrimaryPions_TPCdEdx_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdx_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdx_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_MidPt__temp,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdx_MidPt__temp,minB-1,maxB+1);
                    vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdx_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdx_HighPt
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdx_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistESD_PrimaryPions_TPCdEdx,fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetLast());
                if ((fHistESD_PrimaryPions_TPCdEdx)&&(fHistESD_PrimaryPions_TPCdEdx->GetEntries()>0)){
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdx,fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->FindBin(3),fHistESD_PrimaryPions_TPCdEdx->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdx,minYB,maxYB);
                    SetYRange(fHistESD_PrimaryPions_TPCdEdx,minYB-1,maxYB+1);
                    TH1D* fHistESD_PrimaryPions_TPCdEdx_HighPt__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdx->ProjectionY(Form("%s",histoName.Data()));
                    hESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdx_HighPt__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdx_HighPt__temp->GetMeanError());
                    fHistESD_PrimaryPions_TPCdEdx_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdx_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdx_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdx_HighPt__temp,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdx_HighPt__temp,minB-1,maxB+1);
                    vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdx_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdxSignal
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistESD_PrimaryPions_TPCdEdxSignal = (TH2D*) ESDContainer->FindObject(Form("%s",histoName.Data()));
                if ((fHistESD_PrimaryPions_TPCdEdxSignal)&&(fHistESD_PrimaryPions_TPCdEdxSignal->GetEntries()>0)){
                    TH1D* fHistESD_PrimaryPions_TPCdEdxSignal__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdxSignal->ProjectionY(Form("%s",histoName.Data(),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetLast()));
                    hESD_PrimaryPions_TPCdEdxSignal_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdxSignal__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdxSignal_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdxSignal__temp->GetMeanError());
                    fHistESD_PrimaryPions_TPCdEdxSignal__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdxSignal__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdxSignal__temp->Sumw2();
                    vecESD_PrimaryPions_TPCdEdxSignal_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdxSignal__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdxSignal_LowPt
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal,fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetLast());
                if ((fHistESD_PrimaryPions_TPCdEdxSignal)&&(fHistESD_PrimaryPions_TPCdEdxSignal->GetEntries()>0)){
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal,0,fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdxSignal,minYB,maxYB);
                    SetYRange(fHistESD_PrimaryPions_TPCdEdxSignal,minYB-1,maxYB+1);
                    TH1D* fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdxSignal->ProjectionY(Form("%s",histoName.Data()));
                    hESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp->GetMeanError());
                    fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp,minB-1,maxB+1);
                    vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdxSignal_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdxSignal_MidPt
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal,fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetLast());
                if ((fHistESD_PrimaryPions_TPCdEdxSignal)&&(fHistESD_PrimaryPions_TPCdEdxSignal->GetEntries()>0)){
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal,fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->FindBin(2),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdxSignal,minYB,maxYB);
                    SetYRange(fHistESD_PrimaryPions_TPCdEdxSignal,minYB-1,maxYB+1);
                    TH1D* fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdxSignal->ProjectionY(Form("%s",histoName.Data()));
                    hESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp->GetMeanError());
                    fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp,minB-1,maxB+1);
                    vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdxSignal_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //ESD_PrimaryPions_TPCdEdxSignal_HighPt
                if (iParticleType==0){histoName="ESD_PrimaryPions_TPCdEdxSignal_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal,fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetFirst(),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetLast());
                if ((fHistESD_PrimaryPions_TPCdEdxSignal)&&(fHistESD_PrimaryPions_TPCdEdxSignal->GetEntries()>0)){
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal,fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->FindBin(3),fHistESD_PrimaryPions_TPCdEdxSignal->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistESD_PrimaryPions_TPCdEdxSignal,minYB,maxYB);
                    SetYRange(fHistESD_PrimaryPions_TPCdEdxSignal,minYB-1,maxYB+1);
                    TH1D* fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp = (TH1D*) fHistESD_PrimaryPions_TPCdEdxSignal->ProjectionY(Form("%s",histoName.Data()));
                    hESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i]->SetBinContent(bin, fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp->GetMean());
                    hESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i]->SetBinError(bin, fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp->GetMeanError());
                    fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp,minB,maxB);
                    SetXRange(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp,minB-1,maxB+1);
                    vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i].push_back(fHistESD_PrimaryPions_TPCdEdxSignal_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
            }else {cout<<"no ESDContainer Projections"<<endl;}
            //*******************************************************************************************************************************
            //After QA
            if (PionCutsContainer!=NULL){
                cout<<"PionCutsContainer Projections, After QA: "<<endl<<fPionCutsContainerCutString.Data()<<endl;
                //AfterQA: Pion_ITS_after
                if (iParticleType==0){histoName="Pion_ITS_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_ITS_after_AfterQA = (TH2D*) PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHistPion_ITS_after_AfterQA)&&(fHistPion_ITS_after_AfterQA->GetEntries()>0)){
                    TH1D* fHistPion_ITS_after_AfterQA__temp = (TH1D*) fHistPion_ITS_after_AfterQA->ProjectionY(Form("%s",histoName.Data(),fHistPion_ITS_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA->GetXaxis()->GetLast()));
                    hPion_ITS_after_AfterQA_ProjPt[i]->SetBinContent(bin, fHistPion_ITS_after_AfterQA__temp->GetMean());
                    hPion_ITS_after_AfterQA_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_AfterQA__temp->GetMeanError());
                    fHistPion_ITS_after_AfterQA__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_AfterQA__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_AfterQA__temp->Sumw2();
                    vecPion_ITS_after_AfterQA_ProjPt[i].push_back(fHistPion_ITS_after_AfterQA__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_ITS_after_LowPt
                if (iParticleType==0){histoName="Pion_ITS_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_after_AfterQA,fHistPion_ITS_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_ITS_after_AfterQA)&&(fHistPion_ITS_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_ITS_after_AfterQA,0,fHistPion_ITS_after_AfterQA->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_ITS_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_ITS_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_after_AfterQA_LowPt__temp = (TH1D*) fHistPion_ITS_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_after_AfterQA_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_ITS_after_AfterQA_LowPt__temp->GetMean());
                    hPion_ITS_after_AfterQA_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_AfterQA_LowPt__temp->GetMeanError());
                    fHistPion_ITS_after_AfterQA_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_AfterQA_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_AfterQA_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_after_AfterQA_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_after_AfterQA_LowPt__temp,minB-1,maxB+1);
                    vecPion_ITS_after_AfterQA_LowPt_ProjPt[i].push_back(fHistPion_ITS_after_AfterQA_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_ITS_after_MidPt
                if (iParticleType==0){histoName="Pion_ITS_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_after_AfterQA,fHistPion_ITS_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_ITS_after_AfterQA)&&(fHistPion_ITS_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_ITS_after_AfterQA,fHistPion_ITS_after_AfterQA->GetXaxis()->FindBin(2),fHistPion_ITS_after_AfterQA->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_ITS_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_ITS_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_after_AfterQA_MidPt__temp = (TH1D*) fHistPion_ITS_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_after_AfterQA_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_ITS_after_AfterQA_MidPt__temp->GetMean());
                    hPion_ITS_after_AfterQA_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_AfterQA_MidPt__temp->GetMeanError());
                    fHistPion_ITS_after_AfterQA_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_AfterQA_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_AfterQA_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_after_AfterQA_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_after_AfterQA_MidPt__temp,minB-1,maxB+1);
                    vecPion_ITS_after_AfterQA_MidPt_ProjPt[i].push_back(fHistPion_ITS_after_AfterQA_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_ITS_after_HighPt
                if (iParticleType==0){histoName="Pion_ITS_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_after_AfterQA,fHistPion_ITS_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_ITS_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_ITS_after_AfterQA)&&(fHistPion_ITS_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_ITS_after_AfterQA,fHistPion_ITS_after_AfterQA->GetXaxis()->FindBin(2),fHistPion_ITS_after_AfterQA->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_ITS_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_ITS_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_after_AfterQA_HighPt__temp = (TH1D*) fHistPion_ITS_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_after_AfterQA_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_ITS_after_AfterQA_HighPt__temp->GetMean());
                    hPion_ITS_after_AfterQA_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_AfterQA_HighPt__temp->GetMeanError());
                    fHistPion_ITS_after_AfterQA_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_AfterQA_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_AfterQA_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_after_AfterQA_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_after_AfterQA_HighPt__temp,minB-1,maxB+1);
                    vecPion_ITS_after_AfterQA_HighPt_ProjPt[i].push_back(fHistPion_ITS_after_AfterQA_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdx_after
                if (iParticleType==0){histoName="Pion_dEdx_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_dEdx_after_AfterQA = (TH2D*) PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHistPion_dEdx_after_AfterQA)&&(fHistPion_dEdx_after_AfterQA->GetEntries()>0)){
                    TH1D* fHistPion_dEdx_after_AfterQA__temp = (TH1D*) fHistPion_dEdx_after_AfterQA->ProjectionY(Form("%s",histoName.Data(),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetLast()));
                    hPion_dEdx_after_AfterQA_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_AfterQA__temp->GetMean());
                    hPion_dEdx_after_AfterQA_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_AfterQA__temp->GetMeanError());
                    fHistPion_dEdx_after_AfterQA__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_AfterQA__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_AfterQA__temp->Sumw2();
                    vecPion_dEdx_after_AfterQA_ProjPt[i].push_back(fHistPion_dEdx_after_AfterQA__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdx_after_LowPt
                if (iParticleType==0){histoName="Pion_dEdx_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_after_AfterQA,fHistPion_dEdx_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_after_AfterQA)&&(fHistPion_dEdx_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_AfterQA,0,fHistPion_dEdx_after_AfterQA->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_dEdx_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_after_AfterQA_LowPt__temp = (TH1D*) fHistPion_dEdx_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_after_AfterQA_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_AfterQA_LowPt__temp->GetMean());
                    hPion_dEdx_after_AfterQA_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_AfterQA_LowPt__temp->GetMeanError());
                    fHistPion_dEdx_after_AfterQA_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_AfterQA_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_AfterQA_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_after_AfterQA_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_AfterQA_LowPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_after_AfterQA_LowPt_ProjPt[i].push_back(fHistPion_dEdx_after_AfterQA_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdx_after_MidPt
                if (iParticleType==0){histoName="Pion_dEdx_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_after_AfterQA,fHistPion_dEdx_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_after_AfterQA)&&(fHistPion_dEdx_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_AfterQA,fHistPion_dEdx_after_AfterQA->GetXaxis()->FindBin(2),fHistPion_dEdx_after_AfterQA->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_dEdx_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_after_AfterQA_MidPt__temp = (TH1D*) fHistPion_dEdx_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_after_AfterQA_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_AfterQA_MidPt__temp->GetMean());
                    hPion_dEdx_after_AfterQA_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_AfterQA_MidPt__temp->GetMeanError());
                    fHistPion_dEdx_after_AfterQA_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_AfterQA_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_AfterQA_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_after_AfterQA_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_AfterQA_MidPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_after_AfterQA_MidPt_ProjPt[i].push_back(fHistPion_dEdx_after_AfterQA_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdx_after_HighPt
                if (iParticleType==0){histoName="Pion_dEdx_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_after_AfterQA,fHistPion_dEdx_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_after_AfterQA)&&(fHistPion_dEdx_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_AfterQA,fHistPion_dEdx_after_AfterQA->GetXaxis()->FindBin(2),fHistPion_dEdx_after_AfterQA->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_dEdx_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_after_AfterQA_HighPt__temp = (TH1D*) fHistPion_dEdx_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_after_AfterQA_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_AfterQA_HighPt__temp->GetMean());
                    hPion_dEdx_after_AfterQA_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_AfterQA_HighPt__temp->GetMeanError());
                    fHistPion_dEdx_after_AfterQA_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_AfterQA_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_AfterQA_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_after_AfterQA_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_AfterQA_HighPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_after_AfterQA_HighPt_ProjPt[i].push_back(fHistPion_dEdx_after_AfterQA_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdxSignal_after
                if (iParticleType==0){histoName="Pion_dEdxSignal_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_dEdxSignal_after_AfterQA = (TH2D*) PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHistPion_dEdxSignal_after_AfterQA)&&(fHistPion_dEdxSignal_after_AfterQA->GetEntries()>0)){
                    TH1D* fHistPion_dEdxSignal_after_AfterQA__temp = (TH1D*) fHistPion_dEdxSignal_after_AfterQA->ProjectionY(Form("%s",histoName.Data(),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetLast()));
                    hPion_dEdxSignal_after_AfterQA_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_AfterQA__temp->GetMean());
                    hPion_dEdxSignal_after_AfterQA_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_AfterQA__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_AfterQA__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_AfterQA__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_AfterQA__temp->Sumw2();
                    vecPion_dEdxSignal_after_AfterQA_ProjPt[i].push_back(fHistPion_dEdxSignal_after_AfterQA__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdxSignal_after_LowPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_after_AfterQA,fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_after_AfterQA)&&(fHistPion_dEdxSignal_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_AfterQA,0,fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_dEdxSignal_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_after_AfterQA_LowPt__temp = (TH1D*) fHistPion_dEdxSignal_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_AfterQA_LowPt__temp->GetMean());
                    hPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_AfterQA_LowPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_AfterQA_LowPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_AfterQA_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_AfterQA_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_AfterQA_LowPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i].push_back(fHistPion_dEdxSignal_after_AfterQA_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdxSignal_after_MidPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_after_AfterQA,fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_after_AfterQA)&&(fHistPion_dEdxSignal_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_AfterQA,fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_dEdxSignal_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_after_AfterQA_MidPt__temp = (TH1D*) fHistPion_dEdxSignal_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_AfterQA_MidPt__temp->GetMean());
                    hPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_AfterQA_MidPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_AfterQA_MidPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_AfterQA_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_AfterQA_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_AfterQA_MidPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i].push_back(fHistPion_dEdxSignal_after_AfterQA_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_dEdxSignal_after_HighPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_after_AfterQA,fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_after_AfterQA)&&(fHistPion_dEdxSignal_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_AfterQA,fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->FindBin(3),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_dEdxSignal_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_after_AfterQA_HighPt__temp = (TH1D*) fHistPion_dEdxSignal_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_AfterQA_HighPt__temp->GetMean());
                    hPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_AfterQA_HighPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_AfterQA_HighPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_AfterQA_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_AfterQA_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_after_AfterQA_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_AfterQA_HighPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i].push_back(fHistPion_dEdxSignal_after_AfterQA_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_TOF_after
                if (iParticleType==0){histoName="Pion_TOF_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_TOF_after_AfterQA = (TH2D*) PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHistPion_TOF_after_AfterQA)&&(fHistPion_TOF_after_AfterQA->GetEntries()>0)){
                    TH1D* fHistPion_TOF_after_AfterQA__temp = (TH1D*) fHistPion_TOF_after_AfterQA->ProjectionY(Form("%s",histoName.Data(),fHistPion_TOF_after_AfterQA__temp->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA__temp->GetXaxis()->GetLast()));
                    hPion_TOF_after_AfterQA_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_AfterQA__temp->GetMean());
                    hPion_TOF_after_AfterQA_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_AfterQA__temp->GetMeanError());
                    fHistPion_TOF_after_AfterQA__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_AfterQA__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_AfterQA__temp->Sumw2();
                    vecPion_TOF_after_AfterQA_ProjPt[i].push_back(fHistPion_TOF_after_AfterQA__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_TOF_after_LowPt
                if (iParticleType==0){histoName="Pion_TOF_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_after_AfterQA,fHistPion_TOF_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_TOF_after_AfterQA)&&(fHistPion_TOF_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_TOF_after_AfterQA,0,fHistPion_TOF_after_AfterQA->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_TOF_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_TOF_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_after_AfterQA_LowPt__temp = (TH1D*) fHistPion_TOF_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_after_AfterQA_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_AfterQA_LowPt__temp->GetMean());
                    hPion_TOF_after_AfterQA_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_AfterQA_LowPt__temp->GetMeanError());
                    fHistPion_TOF_after_AfterQA_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_AfterQA_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_AfterQA_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_after_AfterQA_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_after_AfterQA_LowPt__temp,minB-1,maxB+1);
                    vecPion_TOF_after_AfterQA_LowPt_ProjPt[i].push_back(fHistPion_TOF_after_AfterQA_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_TOF_after_MidPt
                if (iParticleType==0){histoName="Pion_TOF_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_after_AfterQA,fHistPion_TOF_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_TOF_after_AfterQA)&&(fHistPion_TOF_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_TOF_after_AfterQA,fHistPion_TOF_after_AfterQA->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_after_AfterQA->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_TOF_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_TOF_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_after_AfterQA_MidPt__temp = (TH1D*) fHistPion_TOF_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_after_AfterQA_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_AfterQA_MidPt__temp->GetMean());
                    hPion_TOF_after_AfterQA_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_AfterQA_MidPt__temp->GetMeanError());
                    fHistPion_TOF_after_AfterQA_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_AfterQA_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_AfterQA_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_after_AfterQA_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_after_AfterQA_MidPt__temp,minB-1,maxB+1);
                    vecPion_TOF_after_AfterQA_MidPt_ProjPt[i].push_back(fHistPion_TOF_after_AfterQA_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: Pion_TOF_after_HighPt
                if (iParticleType==0){histoName="Pion_TOF_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_after_AfterQA,fHistPion_TOF_after_AfterQA->GetXaxis()->GetFirst(),fHistPion_TOF_after_AfterQA->GetXaxis()->GetLast());
                if ((fHistPion_TOF_after_AfterQA)&&(fHistPion_TOF_after_AfterQA->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_after_AfterQA,minB,maxB);
                    SetXRange(fHistPion_TOF_after_AfterQA,fHistPion_TOF_after_AfterQA->GetXaxis()->FindBin(3),fHistPion_TOF_after_AfterQA->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_TOF_after_AfterQA,minYB,maxYB);
                    SetYRange(fHistPion_TOF_after_AfterQA,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_after_AfterQA_HighPt__temp = (TH1D*) fHistPion_TOF_after_AfterQA->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_after_AfterQA_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_AfterQA_HighPt__temp->GetMean());
                    hPion_TOF_after_AfterQA_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_AfterQA_HighPt__temp->GetMeanError());
                    fHistPion_TOF_after_AfterQA_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_AfterQA_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_AfterQA_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_after_AfterQA_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_after_AfterQA_HighPt__temp,minB-1,maxB+1);
                    vecPion_TOF_after_AfterQA_HighPt_ProjPt[i].push_back(fHistPion_TOF_after_AfterQA_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: hTrack_DCAxy_Pt_after
                if (iParticleType==0){histoName="hTrack_DCAxy_Pt_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_DCAxy_Pt_after_AfterQA = (TH2D*) PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHisthTrack_DCAxy_Pt_after_AfterQA)&&(fHisthTrack_DCAxy_Pt_after_AfterQA->GetEntries()>0)){
                    TH1D* fHisthTrack_DCAxy_Pt_after_AfterQA__temp = (TH1D*) fHisthTrack_DCAxy_Pt_after_AfterQA->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_DCAxy_Pt_after_AfterQA->GetYaxis()->GetFirst(),fHisthTrack_DCAxy_Pt_after_AfterQA->GetYaxis()->GetLast()));
                    hTrack_DCAxy_Pt_after_AfterQA_ProjPt[i]->SetBinContent(bin, fHisthTrack_DCAxy_Pt_after_AfterQA__temp->GetMean());
                    hTrack_DCAxy_Pt_after_AfterQA_ProjPt[i]->SetBinError(bin, fHisthTrack_DCAxy_Pt_after_AfterQA__temp->GetMeanError());
                    fHisthTrack_DCAxy_Pt_after_AfterQA__temp->GetYaxis()->SetTitle("DCA_{#it{xy}} (cm)");
                    fHisthTrack_DCAxy_Pt_after_AfterQA__temp->Sumw2();
                    vechTrack_DCAxy_Pt_after_AfterQA_ProjPt[i].push_back(fHisthTrack_DCAxy_Pt_after_AfterQA__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: hTrack_DCAz_Pt_after
                if (iParticleType==0){histoName="hTrack_DCAz_Pt_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_DCAz_Pt_after_AfterQA = (TH2D*) PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHisthTrack_DCAz_Pt_after_AfterQA)&&(fHisthTrack_DCAz_Pt_after_AfterQA->GetEntries()>0)){
                    TH1D* fHisthTrack_DCAz_Pt_after_AfterQA__temp = (TH1D*) fHisthTrack_DCAz_Pt_after_AfterQA->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_DCAz_Pt_after_AfterQA->GetYaxis()->GetFirst(),fHisthTrack_DCAz_Pt_after_AfterQA->GetYaxis()->GetLast()));
                    hTrack_DCAz_Pt_after_AfterQA_ProjPt[i]->SetBinContent(bin, fHisthTrack_DCAz_Pt_after_AfterQA__temp->GetMean());
                    hTrack_DCAz_Pt_after_AfterQA_ProjPt[i]->SetBinError(bin, fHisthTrack_DCAz_Pt_after_AfterQA__temp->GetMeanError());
                    fHisthTrack_DCAz_Pt_after_AfterQA__temp->GetYaxis()->SetTitle("DCA_{#it{z}} (cm)");
                    fHisthTrack_DCAz_Pt_after_AfterQA__temp->Sumw2();
                    vechTrack_DCAz_Pt_after_AfterQA_ProjPt[i].push_back(fHisthTrack_DCAz_Pt_after_AfterQA__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //AfterQA: hTrack_NFindCls_Pt_TPC_after
                if (iParticleType==0){histoName="hTrack_NFindCls_Pt_TPC_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_NFindCls_Pt_TPC_after_AfterQA = (TH2D*) PionCutsContainer->FindObject(Form("%s %s",histoName.Data(),fPionCutsContainerCutString.Data()));
                if ((fHisthTrack_NFindCls_Pt_TPC_after_AfterQA)&&(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA->GetEntries()>0)){
                    TH1D* fHisthTrack_NFindCls_Pt_TPC_after_AfterQA__temp = (TH1D*) fHisthTrack_NFindCls_Pt_TPC_after_AfterQA->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_NFindCls_Pt_TPC_after_AfterQA->GetYaxis()->GetFirst(),fHisthTrack_NFindCls_Pt_TPC_after_AfterQA->GetYaxis()->GetLast()));
                    hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i]->SetBinContent(bin, fHisthTrack_NFindCls_Pt_TPC_after_AfterQA__temp->GetMean());
                    hTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i]->SetBinError(bin, fHisthTrack_NFindCls_Pt_TPC_after_AfterQA__temp->GetMeanError());
                    fHisthTrack_NFindCls_Pt_TPC_after_AfterQA__temp->GetXaxis()->SetTitle("Findable Clusters after Cut");
                    fHisthTrack_NFindCls_Pt_TPC_after_AfterQA__temp->Sumw2();
                    vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i].push_back(fHisthTrack_NFindCls_Pt_TPC_after_AfterQA__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
            }else {cout<<"no PionCutsContainer Projections, After QA"<<endl<<fPionCutsContainerCutString.Data()<<endl;}
            //*******************************************************************************************************************************
            //BEFORE (In Pre Selection)
            if (PionCuts2Container!=NULL){
                cout<<"PionCuts2Container Projections, BEFORE (In Pre Selection): "<<endl<<fPionCuts2ContainerCutString.Data()<<endl;;
                //Pre Selection: Pion_ITS_before
                if (iParticleType==0){histoName="Pion_ITS_before";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_ITS_before_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_ITS_before_PreSel)&&(fHistPion_ITS_before_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_ITS_before_PreSel__temp = (TH1D*) fHistPion_ITS_before_PreSel->ProjectionY(Form("%s",histoName.Data(),fHistPion_ITS_before_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel->GetXaxis()->GetLast()));
                    hPion_ITS_before_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_ITS_before_PreSel__temp->GetMean());
                    hPion_ITS_before_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_ITS_before_PreSel__temp->GetMeanError());
                    fHistPion_ITS_before_PreSel__temp->GetXaxis()->SetTitle("#it{p}_{#pi} (GeV/#it{c})");
                    fHistPion_ITS_before_PreSel__temp->GetYaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_before_PreSel__temp->Sumw2();
                    vecPion_ITS_before_PreSel_ProjPt[i].push_back(fHistPion_ITS_before_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_ITS_before_LowPt
                if (iParticleType==0){histoName="Pion_ITS_before_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_before_PreSel,fHistPion_ITS_before_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_ITS_before_PreSel)&&(fHistPion_ITS_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_ITS_before_PreSel,0,fHistPion_ITS_before_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_ITS_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_ITS_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_before_PreSel_LowPt__temp = (TH1D*) fHistPion_ITS_before_PreSel->ProjectionY(Form("%s",histoName.Data(),fHistPion_ITS_before_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel->GetXaxis()->GetLast()));
                    hPion_ITS_before_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_before_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_ITS_before_PreSel_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_before_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_before_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_before_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_before_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_ITS_before_PreSel_LowPt_ProjPt[i].push_back(fHistPion_ITS_before_PreSel_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_ITS_before_MidPt
                if (iParticleType==0){histoName="Pion_ITS_before_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_before_PreSel,fHistPion_ITS_before_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_ITS_before_PreSel)&&(fHistPion_ITS_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_ITS_before_PreSel,fHistPion_ITS_before_PreSel->GetXaxis()->FindBin(2),fHistPion_ITS_before_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_ITS_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_ITS_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_before_PreSel_MidPt__temp = (TH1D*) fHistPion_ITS_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_before_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_before_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_ITS_before_PreSel_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_before_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_before_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_before_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_before_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_ITS_before_PreSel_MidPt_ProjPt[i].push_back(fHistPion_ITS_before_PreSel_MidPt__temp);
                }//else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_ITS_before_HighPt
                if (iParticleType==0){histoName="Pion_ITS_before_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_before_PreSel,fHistPion_ITS_before_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_ITS_before_PreSel)&&(fHistPion_ITS_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_ITS_before_PreSel,fHistPion_ITS_before_PreSel->GetXaxis()->FindBin(2),fHistPion_ITS_before_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_ITS_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_ITS_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_before_PreSel_HighPt__temp = (TH1D*) fHistPion_ITS_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_before_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_before_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_ITS_before_PreSel_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_before_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_before_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_before_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_before_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_ITS_before_PreSel_HighPt_ProjPt[i].push_back(fHistPion_ITS_before_PreSel_HighPt__temp);
                }//else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_before
                if (iParticleType==0){histoName="Pion_dEdx_before";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_dEdx_before_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_dEdx_before_PreSel)&&(fHistPion_dEdx_before_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_dEdx_before_PreSel__temp = (TH1D*) fHistPion_dEdx_before_PreSel->ProjectionY(Form("%s",histoName.Data()),fHistPion_dEdx_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel->GetXaxis()->GetLast());
                    hPion_dEdx_before_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_before_PreSel__temp->GetMean());
                    hPion_dEdx_before_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_before_PreSel__temp->GetMeanError());
                    fHistPion_dEdx_before_PreSel__temp->GetXaxis()->SetTitle("#it{p}_{#pi} (GeV/#it{c})");
                    fHistPion_dEdx_before_PreSel__temp->GetYaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_before_PreSel__temp->Sumw2();
                    vecPion_dEdx_before_PreSel_ProjPt[i].push_back(fHistPion_dEdx_before_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_before_LowPt
                if (iParticleType==0){histoName="Pion_dEdx_before_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_before_PreSel,fHistPion_dEdx_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_before_PreSel)&&(fHistPion_dEdx_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdx_before_PreSel,0,fHistPion_dEdx_before_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_dEdx_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_before_PreSel_LowPt__temp = (TH1D*) fHistPion_dEdx_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_before_PreSel_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_before_PreSel_LowPt__temp->GetMean());
                    hPion_dEdx_before_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_before_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_dEdx_before_PreSel_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_before_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_before_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_before_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_before_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_before_PreSel_LowPt_ProjPt[i].push_back(fHistPion_dEdx_before_PreSel_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_before_MidPt
                if (iParticleType==0){histoName="Pion_dEdx_before_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_before_PreSel,fHistPion_dEdx_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_before_PreSel)&&(fHistPion_dEdx_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdx_before_PreSel,fHistPion_dEdx_before_PreSel->GetXaxis()->FindBin(2),fHistPion_dEdx_before_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_dEdx_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_before_PreSel_MidPt__temp = (TH1D*) fHistPion_dEdx_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_before_PreSel_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_before_PreSel_MidPt__temp->GetMean());
                    hPion_dEdx_before_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_before_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_dEdx_before_PreSel_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_before_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_before_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_before_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_before_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_before_PreSel_MidPt_ProjPt[i].push_back(fHistPion_dEdx_before_PreSel_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_before_HighPt
                if (iParticleType==0){histoName="Pion_dEdx_before_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_before_PreSel,fHistPion_dEdx_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_before_PreSel)&&(fHistPion_dEdx_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdx_before_PreSel,fHistPion_dEdx_before_PreSel->GetXaxis()->FindBin(2),fHistPion_dEdx_before_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_dEdx_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_before_PreSel_HighPt__temp = (TH1D*) fHistPion_dEdx_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_before_PreSel_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_before_PreSel_HighPt__temp->GetMean());
                    hPion_dEdx_before_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_before_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_dEdx_before_PreSel_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_before_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_before_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_before_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_before_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_before_PreSel_HighPt_ProjPt[i].push_back(fHistPion_dEdx_before_PreSel_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_before
                if (iParticleType==0){histoName="Pion_dEdxSignal_before";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_dEdxSignal_before_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_dEdxSignal_before_PreSel)&&(fHistPion_dEdxSignal_before_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_dEdxSignal_before_PreSel__temp = (TH1D*) fHistPion_dEdxSignal_before_PreSel->ProjectionY(Form("%s",histoName.Data()),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetLast());
                    hPion_dEdxSignal_before_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_before_PreSel__temp->GetMean());
                    hPion_dEdxSignal_before_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_before_PreSel__temp->GetMeanError());
                    fHistPion_dEdxSignal_before_PreSel__temp->GetXaxis()->SetTitle("#it{p}_{#pi} (GeV/#it{c})");
                    fHistPion_dEdxSignal_before_PreSel__temp->GetYaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_before_PreSel__temp->Sumw2();
                    vecPion_dEdxSignal_before_PreSel_ProjPt[i].push_back(fHistPion_dEdxSignal_before_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_before_LowPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_before_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_before_PreSel,fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_before_PreSel)&&(fHistPion_dEdxSignal_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_before_PreSel,0,fHistPion_dEdxSignal_before_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_before_PreSel_LowPt__temp = (TH1D*) fHistPion_dEdxSignal_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_before_PreSel_LowPt__temp->GetMean());
                    hPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_before_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_before_PreSel_LowPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_before_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_before_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_before_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i].push_back(fHistPion_dEdxSignal_before_PreSel_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_before_MidPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_before_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_before_PreSel,fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_before_PreSel)&&(fHistPion_dEdxSignal_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_before_PreSel,fHistPion_dEdxSignal_before_PreSel->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_before_PreSel_MidPt__temp = (TH1D*) fHistPion_dEdxSignal_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_before_PreSel_MidPt__temp->GetMean());
                    hPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_before_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_before_PreSel_MidPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_before_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_before_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_before_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i].push_back(fHistPion_dEdxSignal_before_PreSel_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_before_HighPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_before_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_before_PreSel,fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_before_PreSel)&&(fHistPion_dEdxSignal_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_before_PreSel,fHistPion_dEdxSignal_before_PreSel->GetXaxis()->FindBin(3),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_before_PreSel_HighPt__temp = (TH1D*) fHistPion_dEdxSignal_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_before_PreSel_HighPt__temp->GetMean());
                    hPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_before_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_before_PreSel_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_dEdxSignal_before_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_before_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_before_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_before_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i].push_back(fHistPion_dEdxSignal_before_PreSel_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_before_PreSel
                if (iParticleType==0){histoName="Pion_TOF_before";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_TOF_before_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_TOF_before_PreSel)&&(fHistPion_TOF_before_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_TOF_before_PreSel__temp = (TH1D*) fHistPion_TOF_before_PreSel->ProjectionY(Form("%s",histoName.Data()),fHistPion_TOF_before_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel->GetXaxis()->GetLast());
                    hPion_TOF_before_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_before_PreSel__temp->GetMean());
                    hPion_TOF_before_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_TOF_before_PreSel__temp->GetMeanError());
                    fHistPion_TOF_before_PreSel__temp->GetXaxis()->SetTitle("#it{p}_{T, #pi} (GeV/#it{c})");
                    fHistPion_TOF_before_PreSel__temp->GetYaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_before_PreSel__temp->Sumw2();
                    vecPion_TOF_before_PreSel_ProjPt[i].push_back(fHistPion_TOF_before_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_before_PreSel_LowPt
                if (iParticleType==0){histoName="Pion_TOF_before_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_before_PreSel,fHistPion_TOF_before_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_TOF_before_PreSel)&&(fHistPion_TOF_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_TOF_before_PreSel,0,fHistPion_TOF_before_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_TOF_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_TOF_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_before_PreSel_LowPt__temp = (TH1D*) fHistPion_TOF_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_before_PreSel_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_before_PreSel_LowPt__temp->GetMean());
                    hPion_TOF_before_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_before_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_TOF_before_PreSel_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_before_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_before_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_before_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_before_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_TOF_before_PreSel_LowPt_ProjPt[i].push_back(fHistPion_TOF_before_PreSel_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_before_PreSel_MidPt
                if (iParticleType==0){histoName="Pion_TOF_before_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_before_PreSel,fHistPion_TOF_before_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_TOF_before_PreSel)&&(fHistPion_TOF_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_TOF_before_PreSel,fHistPion_TOF_before_PreSel->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_before_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_dEdxSignal_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_before_PreSel_MidPt__temp = (TH1D*) fHistPion_TOF_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_before_PreSel_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_before_PreSel_MidPt__temp->GetMean());
                    hPion_TOF_before_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_before_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_TOF_before_PreSel_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_before_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_before_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_before_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_before_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_TOF_before_PreSel_MidPt_ProjPt[i].push_back(fHistPion_TOF_before_PreSel_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_before_PreSel_HighPt
                if (iParticleType==0){histoName="Pion_TOF_before_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_before_PreSel,fHistPion_TOF_before_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_before_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_TOF_before_PreSel)&&(fHistPion_TOF_before_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_before_PreSel,minB,maxB);
                    SetXRange(fHistPion_TOF_before_PreSel,fHistPion_TOF_before_PreSel->GetXaxis()->FindBin(3),fHistPion_TOF_before_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_TOF_before_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_TOF_before_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_before_PreSel_HighPt__temp = (TH1D*) fHistPion_TOF_before_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_before_PreSel_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_before_PreSel_HighPt__temp->GetMean());
                    hPion_TOF_before_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_before_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_TOF_before_PreSel_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_before_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_before_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_before_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_before_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_TOF_before_PreSel_HighPt_ProjPt[i].push_back(fHistPion_TOF_before_PreSel_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: hTrack_DCAxy_Pt_before
                if (iParticleType==0){histoName="hTrack_DCAxy_Pt_before";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_DCAxy_Pt_before_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHisthTrack_DCAxy_Pt_before_PreSel)&&(fHisthTrack_DCAxy_Pt_before_PreSel->GetEntries()>0)){
                    TH1D* fHisthTrack_DCAxy_Pt_before_PreSel__temp = (TH1D*) fHisthTrack_DCAxy_Pt_before_PreSel->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_DCAxy_Pt_before_PreSel->GetYaxis()->GetFirst(),fHisthTrack_DCAxy_Pt_before_PreSel->GetYaxis()->GetLast()));
                    hTrack_DCAxy_Pt_before_PreSel_ProjPt[i]->SetBinContent(bin, fHisthTrack_DCAxy_Pt_before_PreSel__temp->GetMean());
                    hTrack_DCAxy_Pt_before_PreSel_ProjPt[i]->SetBinError(bin, fHisthTrack_DCAxy_Pt_before_PreSel__temp->GetMeanError());
                    fHisthTrack_DCAxy_Pt_before_PreSel__temp->GetXaxis()->SetTitle("DCA_{#it{xy}} (cm)");
                    fHisthTrack_DCAxy_Pt_before_PreSel__temp->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    fHisthTrack_DCAxy_Pt_before_PreSel__temp->Sumw2();
                    //DrawPeriodQAHistoTH1(canvas1,0.09,0.02,0.04,0.09,kTRUE,kTRUE,kFALSE,
                    //                     fHisthTrack_DCAxy_Pt_before_PreSel__temp,"","X","Y",1,1,
                    //                     0.82,0.94,0.03,fCollisionSystem,plotDataSets[i],"");
                    //WriteHistogram(fHistESD_PrimaryNegPions_Pt);
                    //SaveCanvasAndWriteHistogram(canvas1, fHisthTrack_DCAxy_Pt_before_PreSel__temp, Form("%s/Test01.%s", outputDir.Data(), suffix.Data()));
                    vechTrack_DCAxy_Pt_before_PreSel_ProjPt[i].push_back(fHisthTrack_DCAxy_Pt_before_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: hTrack_DCAz_Pt_before
                if (iParticleType==0){histoName="hTrack_DCAz_Pt_before";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_DCAz_Pt_before_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHisthTrack_DCAz_Pt_before_PreSel)&&(fHisthTrack_DCAz_Pt_before_PreSel->GetEntries()>0)){
                    TH1D* fHisthTrack_DCAz_Pt_before_PreSel__temp = (TH1D*) fHisthTrack_DCAz_Pt_before_PreSel->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_DCAz_Pt_before_PreSel->GetYaxis()->GetFirst(),fHisthTrack_DCAz_Pt_before_PreSel->GetYaxis()->GetLast()));
                    hTrack_DCAz_Pt_before_PreSel_ProjPt[i]->SetBinContent(bin, fHisthTrack_DCAz_Pt_before_PreSel__temp->GetMean());
                    hTrack_DCAz_Pt_before_PreSel_ProjPt[i]->SetBinError(bin, fHisthTrack_DCAz_Pt_before_PreSel__temp->GetMeanError());
                    fHisthTrack_DCAz_Pt_before_PreSel__temp->GetXaxis()->SetTitle("DCA_{#it{z}} (cm)");
                    fHisthTrack_DCAz_Pt_before_PreSel__temp->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    fHisthTrack_DCAz_Pt_before_PreSel__temp->Sumw2();
                    vechTrack_DCAz_Pt_before_PreSel_ProjPt[i].push_back(fHisthTrack_DCAz_Pt_before_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: hTrack_NFindCls_Pt_TPC_before
                if (iParticleType==0){histoName="hTrack_NFindCls_Pt_TPC_before";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_NFindCls_Pt_TPC_before_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHisthTrack_NFindCls_Pt_TPC_before_PreSel)&&(fHisthTrack_NFindCls_Pt_TPC_before_PreSel->GetEntries()>0)){
                    TH1D* fHisthTrack_NFindCls_Pt_TPC_before_PreSel__temp = (TH1D*) fHisthTrack_NFindCls_Pt_TPC_before_PreSel->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_NFindCls_Pt_TPC_before_PreSel->GetYaxis()->GetFirst(),fHisthTrack_NFindCls_Pt_TPC_before_PreSel->GetYaxis()->GetLast()));
                    hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i]->SetBinContent(bin, fHisthTrack_NFindCls_Pt_TPC_before_PreSel__temp->GetMean());
                    hTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i]->SetBinError(bin, fHisthTrack_NFindCls_Pt_TPC_before_PreSel__temp->GetMeanError());
                    fHisthTrack_NFindCls_Pt_TPC_before_PreSel__temp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c}) TPC");
                    fHisthTrack_NFindCls_Pt_TPC_before_PreSel__temp->GetYaxis()->SetTitle("Findable Clusters before Cut");
                    fHisthTrack_NFindCls_Pt_TPC_before_PreSel__temp->Sumw2();
                    vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i].push_back(fHisthTrack_NFindCls_Pt_TPC_before_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
            }else {cout<<"no PionCuts2Container Projections, BEFORE (In Pre Selection): "<<endl<<fPionCuts2ContainerCutString.Data()<<endl;}
            //*******************************************************************************************************************************
            //AFTER (In Pre Selection)
            if (PionCuts2Container!=NULL){
                cout<<"PionCuts2Container Projections, AFTER (In Pre Selection)"<<endl<<fPionCuts2ContainerCutString.Data()<<endl;
                //Pre Selection: Pion_ITS_after
                if (iParticleType==0){histoName="Pion_ITS_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_ITS_after_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_ITS_after_PreSel)&&(fHistPion_ITS_after_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_ITS_after_PreSel__temp = (TH1D*) fHistPion_ITS_after_PreSel->ProjectionY(Form("%s",histoName.Data()),fHistPion_ITS_after_PreSel->GetYaxis()->GetFirst(),fHistPion_ITS_after_PreSel->GetXaxis()->GetLast());
                    hPion_ITS_after_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_ITS_after_PreSel__temp->GetMean());
                    hPion_ITS_after_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_PreSel__temp->GetMeanError());
                    fHistPion_ITS_after_PreSel__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_PreSel__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_PreSel__temp->Sumw2();
                    vecPion_ITS_after_PreSel_ProjPt[i].push_back(fHistPion_ITS_after_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_ITS_after_LowPt
                if (iParticleType==0){histoName="Pion_ITS_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_after_PreSel,fHistPion_ITS_after_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_ITS_after_PreSel)&&(fHistPion_ITS_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_ITS_after_PreSel,0,fHistPion_ITS_after_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_ITS_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_ITS_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_after_PreSel_LowPt__temp = (TH1D*) fHistPion_ITS_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_after_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_ITS_after_PreSel_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_after_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_after_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_ITS_after_PreSel_LowPt_ProjPt[i].push_back(fHistPion_ITS_after_PreSel_LowPt__temp);
                }//else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_ITS_after_MidPt
                if (iParticleType==0){histoName="Pion_ITS_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_after_PreSel,fHistPion_ITS_after_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_ITS_after_PreSel)&&(fHistPion_ITS_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_ITS_after_PreSel,fHistPion_ITS_after_PreSel->GetXaxis()->FindBin(2),fHistPion_ITS_after_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_ITS_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_ITS_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_after_PreSel_MidPt__temp = (TH1D*) fHistPion_ITS_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_after_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_ITS_after_PreSel_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_after_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_after_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_ITS_after_PreSel_MidPt_ProjPt[i].push_back(fHistPion_ITS_after_PreSel_MidPt__temp);
                }//else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_ITS_after_HighPt
                if (iParticleType==0){histoName="Pion_ITS_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_ITS_after_PreSel,fHistPion_ITS_after_PreSel->GetXaxis()->GetFirst(),fHistPion_ITS_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_ITS_after_PreSel)&&(fHistPion_ITS_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_ITS_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_ITS_after_PreSel,fHistPion_ITS_after_PreSel->GetXaxis()->FindBin(3),fHistPion_ITS_after_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_ITS_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_ITS_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_ITS_after_PreSel_HighPt__temp = (TH1D*) fHistPion_ITS_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_ITS_after_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_ITS_after_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_ITS_after_PreSel_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} ITS");
                    fHistPion_ITS_after_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_ITS_after_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_ITS_after_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_ITS_after_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_ITS_after_PreSel_HighPt_ProjPt[i].push_back(fHistPion_ITS_after_PreSel_HighPt__temp);
                }//else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_after
                if (iParticleType==0){histoName="Pion_dEdx_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_dEdx_after_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_dEdx_after_PreSel)&&(fHistPion_dEdx_after_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_dEdx_after_PreSel__temp = (TH1D*) fHistPion_dEdx_after_PreSel->ProjectionY(Form("%s",histoName.Data()),fHistPion_dEdx_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel->GetXaxis()->GetLast());
                    hPion_dEdx_after_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_PreSel__temp->GetMean());
                    hPion_dEdx_after_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_PreSel__temp->GetMeanError());
                    fHistPion_dEdx_after_PreSel__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_PreSel__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_PreSel__temp->Sumw2();
                    vecPion_dEdx_after_PreSel_ProjPt[i].push_back(fHistPion_dEdx_after_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_after_LowPt
                if (iParticleType==0){histoName="Pion_dEdx_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_after_PreSel,fHistPion_dEdx_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_after_PreSel)&&(fHistPion_dEdx_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_PreSel,0,fHistPion_dEdx_after_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_dEdx_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_after_PreSel_LowPt__temp = (TH1D*) fHistPion_dEdx_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_after_PreSel_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_PreSel_LowPt__temp->GetMean());
                    hPion_dEdx_after_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_dEdx_after_PreSel_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_after_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_after_PreSel_LowPt_ProjPt[i].push_back(fHistPion_dEdx_after_PreSel_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_after_MidPt
                if (iParticleType==0){histoName="Pion_dEdx_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_after_PreSel,fHistPion_dEdx_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_after_PreSel)&&(fHistPion_dEdx_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_PreSel,fHistPion_dEdx_after_PreSel->GetXaxis()->FindBin(2),fHistPion_dEdx_after_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_dEdx_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_after_PreSel_MidPt__temp = (TH1D*) fHistPion_dEdx_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_after_PreSel_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_PreSel_MidPt__temp->GetMean());
                    hPion_dEdx_after_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_dEdx_after_PreSel_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_after_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_after_PreSel_MidPt_ProjPt[i].push_back(fHistPion_dEdx_after_PreSel_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdx_after_HighPt
                if (iParticleType==0){histoName="Pion_dEdx_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdx_after_PreSel,fHistPion_dEdx_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdx_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdx_after_PreSel)&&(fHistPion_dEdx_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdx_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_PreSel,fHistPion_dEdx_after_PreSel->GetXaxis()->FindBin(3),fHistPion_dEdx_after_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_dEdx_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdx_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdx_after_PreSel_HighPt__temp = (TH1D*) fHistPion_dEdx_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdx_after_PreSel_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdx_after_PreSel_HighPt__temp->GetMean());
                    hPion_dEdx_after_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdx_after_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_dEdx_after_PreSel_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TPC");
                    fHistPion_dEdx_after_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdx_after_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdx_after_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdx_after_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_dEdx_after_PreSel_HighPt_ProjPt[i].push_back(fHistPion_dEdx_after_PreSel_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_after
                if (iParticleType==0){histoName="Pion_dEdxSignal_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_dEdxSignal_after_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_dEdxSignal_after_PreSel)&&(fHistPion_dEdxSignal_after_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_dEdxSignal_after_PreSel__temp = (TH1D*) fHistPion_dEdxSignal_after_PreSel->ProjectionY(Form("%s",histoName.Data()),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetLast());
                    hPion_dEdxSignal_after_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_PreSel__temp->GetMean());
                    hPion_dEdxSignal_after_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_PreSel__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_PreSel__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_PreSel__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_PreSel__temp->Sumw2();
                    vecPion_dEdxSignal_after_PreSel_ProjPt[i].push_back(fHistPion_dEdxSignal_after_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_after_LowPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_after_PreSel,fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_after_PreSel)&&(fHistPion_dEdxSignal_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_PreSel,0,fHistPion_dEdxSignal_after_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_dEdxSignal_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_after_PreSel_LowPt__temp = (TH1D*) fHistPion_dEdxSignal_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_PreSel_LowPt__temp->GetMean());
                    hPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_PreSel_LowPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i].push_back(fHistPion_dEdxSignal_after_PreSel_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_after_MidPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_after_PreSel,fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_after_PreSel)&&(fHistPion_dEdxSignal_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_PreSel,fHistPion_dEdxSignal_after_PreSel->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_dEdxSignal_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_after_PreSel_MidPt__temp = (TH1D*) fHistPion_dEdxSignal_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_PreSel_MidPt__temp->GetMean());
                    hPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_PreSel_MidPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i].push_back(fHistPion_dEdxSignal_after_PreSel_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_dEdxSignal_after_HighPt
                if (iParticleType==0){histoName="Pion_dEdxSignal_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_dEdxSignal_after_PreSel,fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetFirst(),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_dEdxSignal_after_PreSel)&&(fHistPion_dEdxSignal_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_PreSel,fHistPion_dEdxSignal_after_PreSel->GetXaxis()->FindBin(3),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_dEdxSignal_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_dEdxSignal_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_dEdxSignal_after_PreSel_HighPt__temp = (TH1D*) fHistPion_dEdxSignal_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_dEdxSignal_after_PreSel_HighPt__temp->GetMean());
                    hPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_dEdxSignal_after_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_dEdxSignal_after_PreSel_HighPt__temp->GetXaxis()->SetTitle("d#it{E}/d#it{x} TPC");
                    fHistPion_dEdxSignal_after_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_dEdxSignal_after_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_dEdxSignal_after_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_dEdxSignal_after_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_TOF_after_PreSel_HighPt_ProjPt[i].push_back(fHistPion_dEdxSignal_after_PreSel_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_after_PreSel
                if (iParticleType==0){histoName="Pion_TOF_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHistPion_TOF_after_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHistPion_TOF_after_PreSel)&&(fHistPion_TOF_after_PreSel->GetEntries()>0)){
                    TH1D* fHistPion_TOF_after_PreSel__temp = (TH1D*) fHistPion_TOF_after_PreSel->ProjectionY(Form("%s",histoName.Data()),fHistPion_TOF_after_PreSel->GetYaxis()->GetFirst(),fHistPion_TOF_after_PreSel->GetXaxis()->GetLast());
                    hPion_TOF_after_PreSel_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_PreSel__temp->GetMean());
                    hPion_TOF_after_PreSel_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_PreSel__temp->GetMeanError());
                    fHistPion_TOF_after_PreSel__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_PreSel__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_PreSel__temp->Sumw2();
                    vecPion_TOF_after_PreSel_ProjPt[i].push_back(fHistPion_TOF_after_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_after_PreSel_LowPt
                if (iParticleType==0){histoName="Pion_TOF_after_LowPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_after_PreSel,fHistPion_TOF_after_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_TOF_after_PreSel)&&(fHistPion_TOF_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_TOF_after_PreSel,0,fHistPion_TOF_after_PreSel->GetXaxis()->FindBin(0.2));
                    GetMinMaxBinY(fHistPion_TOF_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_TOF_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_after_PreSel_LowPt__temp = (TH1D*) fHistPion_TOF_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_after_PreSel_LowPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_PreSel_LowPt__temp->GetMean());
                    hPion_TOF_after_PreSel_LowPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_PreSel_LowPt__temp->GetMeanError());
                    fHistPion_TOF_after_PreSel_LowPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_PreSel_LowPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_PreSel_LowPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_after_PreSel_LowPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_after_PreSel_LowPt__temp,minB-1,maxB+1);
                    vecPion_TOF_after_PreSel_LowPt_ProjPt[i].push_back(fHistPion_TOF_after_PreSel_LowPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_after_PreSel_MidPt
                if (iParticleType==0){histoName="Pion_TOF_after_MidPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_after_PreSel,fHistPion_TOF_after_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_TOF_after_PreSel)&&(fHistPion_TOF_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_TOF_after_PreSel,fHistPion_TOF_after_PreSel->GetXaxis()->FindBin(2),fHistPion_dEdxSignal_after_PreSel->GetXaxis()->FindBin(3));
                    GetMinMaxBinY(fHistPion_TOF_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_TOF_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_after_PreSel_MidPt__temp = (TH1D*) fHistPion_TOF_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_after_PreSel_MidPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_PreSel_MidPt__temp->GetMean());
                    hPion_TOF_after_PreSel_MidPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_PreSel_MidPt__temp->GetMeanError());
                    fHistPion_TOF_after_PreSel_MidPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_PreSel_MidPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_PreSel_MidPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_after_PreSel_MidPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_after_PreSel_MidPt__temp,minB-1,maxB+1);
                    vecPion_TOF_after_PreSel_MidPt_ProjPt[i].push_back(fHistPion_TOF_after_PreSel_MidPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: Pion_TOF_after_PreSel_HighPt
                if (iParticleType==0){histoName="Pion_TOF_after_HighPt";}
                if (iParticleType==1){histoName="";}
                SetXRange(fHistPion_TOF_after_PreSel,fHistPion_TOF_after_PreSel->GetXaxis()->GetFirst(),fHistPion_TOF_after_PreSel->GetXaxis()->GetLast());
                if ((fHistPion_TOF_after_PreSel)&&(fHistPion_TOF_after_PreSel->GetEntries()>0)){
                    GetMinMaxBin(fHistPion_TOF_after_PreSel,minB,maxB);
                    SetXRange(fHistPion_TOF_after_PreSel,fHistPion_TOF_after_PreSel->GetXaxis()->FindBin(3),fHistPion_TOF_after_PreSel->GetXaxis()->GetLast());
                    GetMinMaxBinY(fHistPion_TOF_after_PreSel,minYB,maxYB);
                    SetYRange(fHistPion_TOF_after_PreSel,minYB-1,maxYB+1);
                    TH1D* fHistPion_TOF_after_PreSel_HighPt__temp = (TH1D*) fHistPion_TOF_after_PreSel->ProjectionY(Form("%s",histoName.Data()));
                    hPion_TOF_after_PreSel_HighPt_ProjPt[i]->SetBinContent(bin, fHistPion_TOF_after_PreSel_HighPt__temp->GetMean());
                    hPion_TOF_after_PreSel_HighPt_ProjPt[i]->SetBinError(bin, fHistPion_TOF_after_PreSel_HighPt__temp->GetMeanError());
                    fHistPion_TOF_after_PreSel_HighPt__temp->GetXaxis()->SetTitle("#it{n} #sigma_{#pi} d#it{E}/d#it{x} TOF");
                    fHistPion_TOF_after_PreSel_HighPt__temp->GetYaxis()->SetTitle("# Entries");
                    fHistPion_TOF_after_PreSel_HighPt__temp->Sumw2();
                    GetMinMaxBin(fHistPion_TOF_after_PreSel_HighPt__temp,minB,maxB);
                    SetXRange(fHistPion_TOF_after_PreSel_HighPt__temp,minB-1,maxB+1);
                    vecPion_TOF_after_PreSel_HighPt_ProjPt[i].push_back(fHistPion_TOF_after_PreSel_HighPt__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: hTrack_DCAxy_Pt_after
                if (iParticleType==0){histoName="hTrack_DCAxy_Pt_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_DCAxy_Pt_after_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHisthTrack_DCAxy_Pt_after_PreSel)&&(fHisthTrack_DCAxy_Pt_after_PreSel->GetEntries()>0)){
                    TH1D* fHisthTrack_DCAxy_Pt_after_PreSel__temp = (TH1D*) fHisthTrack_DCAxy_Pt_after_PreSel->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_DCAxy_Pt_after_PreSel->GetYaxis()->GetFirst(),fHisthTrack_DCAxy_Pt_after_PreSel->GetYaxis()->GetLast()));
                    hTrack_DCAxy_Pt_after_PreSel_ProjPt[i]->SetBinContent(bin, fHisthTrack_DCAxy_Pt_after_PreSel__temp->GetMean());
                    hTrack_DCAxy_Pt_after_PreSel_ProjPt[i]->SetBinError(bin, fHisthTrack_DCAxy_Pt_after_PreSel__temp->GetMeanError());
                    fHisthTrack_DCAxy_Pt_after_PreSel__temp->GetXaxis()->SetTitle("DCA_{#it{xy}} (cm)");
                    fHisthTrack_DCAxy_Pt_after_PreSel__temp->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    fHisthTrack_DCAxy_Pt_after_PreSel__temp->Sumw2();
                    vechTrack_DCAxy_Pt_after_PreSel_ProjPt[i].push_back(fHisthTrack_DCAxy_Pt_after_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: hTrack_DCAz_Pt_after
                if (iParticleType==0){histoName="hTrack_DCAz_Pt_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_DCAz_Pt_after_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHisthTrack_DCAz_Pt_after_PreSel)&&(fHisthTrack_DCAz_Pt_after_PreSel->GetEntries()>0)){
                    TH1D* fHisthTrack_DCAz_Pt_after_PreSel__temp = (TH1D*) fHisthTrack_DCAz_Pt_after_PreSel->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_DCAz_Pt_after_PreSel->GetYaxis()->GetFirst(),fHisthTrack_DCAz_Pt_after_PreSel->GetYaxis()->GetLast()));
                    hTrack_DCAz_Pt_after_PreSel_ProjPt[i]->SetBinContent(bin, fHisthTrack_DCAz_Pt_after_PreSel__temp->GetMean());
                    hTrack_DCAz_Pt_after_PreSel_ProjPt[i]->SetBinError(bin, fHisthTrack_DCAz_Pt_after_PreSel__temp->GetMeanError());
                    fHisthTrack_DCAz_Pt_after_PreSel__temp->GetXaxis()->SetTitle("DCA_{#it{z}} (cm)");
                    fHisthTrack_DCAz_Pt_after_PreSel__temp->GetYaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
                    fHisthTrack_DCAz_Pt_after_PreSel__temp->Sumw2();
                    vechTrack_DCAz_Pt_after_PreSel_ProjPt[i].push_back(fHisthTrack_DCAz_Pt_after_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
                //-------------------------------------------------------------------------------------------------------------------------------
                //Pre Selection: hTrack_NFindCls_Pt_TPC_after
                if (iParticleType==0){histoName="hTrack_NFindCls_Pt_TPC_after";}
                if (iParticleType==1){histoName="";}
                TH2D* fHisthTrack_NFindCls_Pt_TPC_after_PreSel = (TH2D*) PionCuts2Container->FindObject(Form("%s %s",histoName.Data(),fPionCuts2ContainerCutString.Data()));
                if ((fHisthTrack_NFindCls_Pt_TPC_after_PreSel)&&(fHisthTrack_NFindCls_Pt_TPC_after_PreSel->GetEntries()>0)){
                    TH1D* fHisthTrack_NFindCls_Pt_TPC_after_PreSel__temp = (TH1D*) fHisthTrack_NFindCls_Pt_TPC_after_PreSel->ProjectionX(Form("%s",histoName.Data(),fHisthTrack_NFindCls_Pt_TPC_after_PreSel->GetYaxis()->GetFirst(),fHisthTrack_NFindCls_Pt_TPC_after_PreSel->GetYaxis()->GetLast()));
                    hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i]->SetBinContent(bin, fHisthTrack_NFindCls_Pt_TPC_after_PreSel__temp->GetMean());
                    hTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i]->SetBinError(bin, fHisthTrack_NFindCls_Pt_TPC_after_PreSel__temp->GetMeanError());
                    fHisthTrack_NFindCls_Pt_TPC_after_PreSel__temp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c}) TPC");
                    fHisthTrack_NFindCls_Pt_TPC_after_PreSel__temp->GetYaxis()->SetTitle("Findable Clusters before Cut");
                    fHisthTrack_NFindCls_Pt_TPC_after_PreSel__temp->Sumw2();
                    vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i].push_back(fHisthTrack_NFindCls_Pt_TPC_after_PreSel__temp);
                } //else cout << "INFO: Object |fHist"<<histoName<<"| could not be found!" << endl;
            }else {cout<<"no PionCuts2Container Projections, AFTER (In Pre Selection): "<<endl<<fPionCuts2ContainerCutString.Data()<<endl;}
            
            std::vector<TH1D*>* vecHistos                               = new std::vector<TH1D*>[nSets];
            
            //--------------------------------------------------------------------------------------------------------
            //---------------------------- Cleanup -------------------------------------------------------------------
            //--------------------------------------------------------------------------------------------------------
            TopDir->Clear();
            delete TopDir;
            
            RootFile->Close();
            delete RootFile;
        } // end of loop over runs
    } // end of loop over datasets
    
    //****************************** Drawing Histograms ************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Histograms" << endl;
    
    TCanvas* canvas                 = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    Double_t leftMar                = 0.09;
    Double_t rightMar               = 0.025;
    Double_t topMargin              = 0.04;
    Double_t bottomMargin           = 0.09;
    DrawGammaCanvasSettings(canvas, leftMar, rightMar, topMargin, bottomMargin);
    
    if (doHistsForEverySet) {
        TBox *boxLabel              = new TBox(1.37,0.7,1.78,0.83);
        boxLabel->SetFillStyle(0);boxLabel->SetFillColor(0);boxLabel->SetLineColor(1);boxLabel->SetLineWidth(0.6);
        TBox *boxLabel2             = new TBox(-0.4,51,6.5,56.5);
        boxLabel2->SetFillStyle(0);boxLabel2->SetFillColor(0);boxLabel2->SetLineColor(1);boxLabel2->SetLineWidth(0.6);
        
        TString outputDirDataSet;
        
        for(Int_t i=0; i<nSets; i++) {
            cout << "DataSet: " << DataSets[i].Data() << endl;
            outputDirDataSet        = Form("%s/%s",outputDir.Data(), DataSets[i].Data());
            
            if (useDataRunListForMC && i>=nData) {
                outputDirDataSet    = Form("%s/%s-%s", outputDir.Data(), DataSets[i].Data(),DataSets[0].Data());
                cout << "Switch useDataRunListForMC is true, output to: " << outputDirDataSet.Data() << endl;
            }
            gSystem->Exec("mkdir -p "+outputDirDataSet);
            gSystem->Exec("mkdir -p "+outputDirDataSet+"/ExtQA");
            
            TString fTrigger        = "";
            if (i<nData){
                TString fTriggerCut = fEventCutSelection(3,2);
                fTrigger            = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                cout << "Trigger: '" << fTrigger.Data() << "'" << endl;
                if (fTrigger.Contains("not defined")){
                    fTrigger        = "";
                    cout << "INFO: Trigger cut not defined!" << endl;
                }
            }
            if (vecESD_PrimaryNegPions_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryNegPions_Pt[i], vecESD_PrimaryNegPions_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPosPions_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPosPions_Pt[i], vecESD_PrimaryPosPions_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryNegPions_Phi[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryNegPions_Phi[i], vecESD_PrimaryNegPions_Phi[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPosPions_Phi[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPosPions_Phi[i], vecESD_PrimaryPosPions_Phi[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryNegPions_Eta[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryNegPions_Eta[i], vecESD_PrimaryNegPions_Eta[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPosPions_Eta[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPosPions_Eta[i], vecESD_PrimaryPosPions_Eta[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecIsPionSelected_AfterQA[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecIsPionSelected_AfterQA[i], vecIsPionSelected_AfterQA[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecdEdxCuts_AfterQA[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecdEdxCuts_AfterQA[i], vecdEdxCuts_AfterQA[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecIsPionSelected_PreSel[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecIsPionSelected_PreSel[i], vecIsPionSelected_PreSel[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecdEdxCuts_PreSel[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecdEdxCuts_PreSel[i], vecdEdxCuts_PreSel[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecMC_AllPosPions_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecMC_AllPosPions_Pt[i], vecMC_AllPosPions_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecMC_AllNegPions_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecMC_AllNegPions_Pt[i], vecMC_AllNegPions_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecMC_PosPionsFromNeutralMeson_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecMC_PosPionsFromNeutralMeson_Pt[i], vecMC_PosPionsFromNeutralMeson_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecMC_NegPionsFromNeutralMeson_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecMC_NegPionsFromNeutralMeson_Pt[i], vecMC_NegPionsFromNeutralMeson_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_TruePosPion_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_TruePosPion_Pt[i], vecESD_TruePosPion_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_TrueNegPion_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_TrueNegPion_Pt[i], vecESD_TrueNegPion_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_TruePosPionFromNeutralMeson_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_TruePosPionFromNeutralMeson_Pt[i], vecESD_TruePosPionFromNeutralMeson_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_TrueNegPionFromNeutralMeson_Pt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_TrueNegPionFromNeutralMeson_Pt[i], vecESD_TrueNegPionFromNeutralMeson_Pt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            //Projections
            if (vecESD_PrimaryNegPions_ClsTPC_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryNegPions_ClsTPC_ProjPt[i], vecESD_PrimaryNegPions_ClsTPC_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPosPions_ClsTPC_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPosPions_ClsTPC_ProjPt[i], vecESD_PrimaryPosPions_ClsTPC_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_DCAxy_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_DCAxy_ProjPt[i], vecESD_PrimaryPions_DCAxy_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_DCAz_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_DCAz_ProjPt[i], vecESD_PrimaryPions_DCAz_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdx_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdx_ProjPt[i], vecESD_PrimaryPions_TPCdEdx_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i], vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i], vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i], vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdxSignal_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdxSignal_ProjPt[i], vecESD_PrimaryPions_TPCdEdxSignal_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i], vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i], vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i], vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            //After QA
            if (vecPion_ITS_after_AfterQA_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_AfterQA_ProjPt[i], vecPion_ITS_after_AfterQA_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_after_AfterQA_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_AfterQA_LowPt_ProjPt[i], vecPion_ITS_after_AfterQA_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_after_AfterQA_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_AfterQA_MidPt_ProjPt[i], vecPion_ITS_after_AfterQA_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_after_AfterQA_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_AfterQA_HighPt_ProjPt[i], vecPion_ITS_after_AfterQA_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_AfterQA_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_AfterQA_ProjPt[i], vecPion_dEdx_after_AfterQA_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_AfterQA_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_AfterQA_LowPt_ProjPt[i], vecPion_dEdx_after_AfterQA_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_AfterQA_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_AfterQA_MidPt_ProjPt[i], vecPion_dEdx_after_AfterQA_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_AfterQA_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_AfterQA_HighPt_ProjPt[i], vecPion_dEdx_after_AfterQA_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_AfterQA_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_AfterQA_ProjPt[i], vecPion_dEdxSignal_after_AfterQA_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i], vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i], vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i], vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_AfterQA_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_AfterQA_ProjPt[i], vecPion_TOF_after_AfterQA_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_AfterQA_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_AfterQA_LowPt_ProjPt[i], vecPion_TOF_after_AfterQA_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_AfterQA_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_AfterQA_MidPt_ProjPt[i], vecPion_TOF_after_AfterQA_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_AfterQA_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_AfterQA_HighPt_ProjPt[i], vecPion_TOF_after_AfterQA_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_DCAxy_Pt_after_AfterQA_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_DCAxy_Pt_after_AfterQA_ProjPt[i], vechTrack_DCAxy_Pt_after_AfterQA_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_DCAz_Pt_after_AfterQA_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_DCAz_Pt_after_AfterQA_ProjPt[i], vechTrack_DCAz_Pt_after_AfterQA_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i], vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            //BEFORE (In Pre Selection)
            if (vecPion_ITS_before_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_before_PreSel_ProjPt[i], vecPion_ITS_before_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_before_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_before_PreSel_LowPt_ProjPt[i], vecPion_ITS_before_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_before_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_before_PreSel_MidPt_ProjPt[i], vecPion_ITS_before_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_before_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_before_PreSel_HighPt_ProjPt[i], vecPion_ITS_before_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_before_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_before_PreSel_ProjPt[i], vecPion_dEdx_before_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_before_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_before_PreSel_LowPt_ProjPt[i], vecPion_dEdx_before_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_before_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_before_PreSel_MidPt_ProjPt[i], vecPion_dEdx_before_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_before_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_before_PreSel_HighPt_ProjPt[i], vecPion_dEdx_before_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_before_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_before_PreSel_ProjPt[i], vecPion_dEdxSignal_before_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i], vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i], vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i], vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_before_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_before_PreSel_ProjPt[i], vecPion_TOF_before_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_before_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_before_PreSel_LowPt_ProjPt[i], vecPion_TOF_before_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_before_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_before_PreSel_MidPt_ProjPt[i], vecPion_TOF_before_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_before_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_before_PreSel_HighPt_ProjPt[i], vecPion_TOF_before_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_DCAxy_Pt_before_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_DCAxy_Pt_before_PreSel_ProjPt[i], vechTrack_DCAxy_Pt_before_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_DCAz_Pt_before_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_DCAz_Pt_before_PreSel_ProjPt[i], vechTrack_DCAz_Pt_before_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i], vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            //AFTER (In Pre Selection)
            if (vecPion_ITS_after_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_PreSel_ProjPt[i], vecPion_ITS_after_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_after_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_PreSel_LowPt_ProjPt[i], vecPion_ITS_after_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_after_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_PreSel_MidPt_ProjPt[i], vecPion_ITS_after_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_ITS_after_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_ITS_after_PreSel_HighPt_ProjPt[i], vecPion_ITS_after_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_PreSel_ProjPt[i], vecPion_dEdx_after_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_PreSel_LowPt_ProjPt[i], vecPion_dEdx_after_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_PreSel_MidPt_ProjPt[i], vecPion_dEdx_after_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdx_after_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdx_after_PreSel_HighPt_ProjPt[i], vecPion_dEdx_after_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_PreSel_ProjPt[i], vecPion_dEdxSignal_after_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i], vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i], vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i], vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_PreSel_ProjPt[i], vecPion_TOF_after_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_PreSel_LowPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_PreSel_LowPt_ProjPt[i], vecPion_TOF_after_PreSel_LowPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_PreSel_MidPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_PreSel_MidPt_ProjPt[i], vecPion_TOF_after_PreSel_MidPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vecPion_TOF_after_PreSel_HighPt_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vecPion_TOF_after_PreSel_HighPt_ProjPt[i], vecPion_TOF_after_PreSel_HighPt_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_DCAxy_Pt_after_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_DCAxy_Pt_after_PreSel_ProjPt[i], vechTrack_DCAxy_Pt_after_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_DCAz_Pt_after_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_DCAz_Pt_after_PreSel_ProjPt[i], vechTrack_DCAz_Pt_after_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
            if (vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i].size()>0){
                DrawVectorOverviewTH1D( canvas, vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i], vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i].at(0)->GetTitle(), outputDirDataSet, suffix, 0.13, 0.15, 0.1, 0.1, 0.6, 0.8, 0.12, 0.93, 0x0, kFALSE, kFALSE);
            }
        }
        delete boxLabel;
        delete boxLabel2;
    }
    //****************************** Combined Trending Histograms ************************************************
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Trending Histograms" << endl;
    //---------
    if (useDataRunListForMC) cout << "WARNING: useDataRunListForMC is true, overwriting histograms for DataSet!" << endl;
    
    TLegend *legend = new TLegend(0.15,0.95,0.95,0.98);
    legend->SetNColumns(nSets);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->SetTextSize(0.03);
    legend->SetTextFont(42);
    //---------
    if (nSets > 1)
    {
        for(Int_t h=0; h<(Int_t) vecHistos[0].size(); h++)
        {
            TString fTrigger = "";
            TString fTriggerCut = fEventCutSelection(3,2);
            for(Int_t i=0; i<nSets; i++)
            {
                if (doTrigger && i<nData && i==0){
                    fTrigger = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), plotDataSets[i]);
                    if (fTrigger.Contains("not defined")) fTrigger = "";
                }
                legend->AddEntry(((TH1D*) vecHistos[i].at(h)),plotDataSets[i].Data(),"p");
            }
            AdjustHistRange(vecHistos,1.1,1.1,h,nSets,kTRUE);
            for(Int_t i=nSets-1; i>=0; i--)
            {
                TString draw = (i==nSets-1)?"px0e1":"px0e1, same";
                ((TH1D*) vecHistos[i].at(h))->SetTitle("");
                ((TH1D*) vecHistos[i].at(h))->Draw(draw.Data());
            }
            legend->Draw();
            
            if (doTrigger) PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), fTrigger.Data(), "");
            else PutProcessLabelAndEnergyOnPlot(xPosLabel, 0.92, 0.03, fCollisionSystem.Data(), "", "");
            
            if (canvas->GetTopMargin()!=0.06) {
                canvas->SetTopMargin(0.06);}
            if (useDataRunListForMC && !addSubFolder){
                SaveCanvas(canvas, Form("%s/%s/%s.%s", outputDir.Data(), DataSets[0].Data(),vecHistosName.at(h).Data(),suffix.Data()));
            }else{
                SaveCanvas(canvas, Form("%s/%s.%s", outputDir.Data(), vecHistosName.at(h).Data(), suffix.Data()));
            }
            legend->Clear();
        }
    }

    //****************************** Create Output ROOT-File ************************************************

    const char* nameOutput;

    for(Int_t i=0; i<nSets; i++)
    {
        fDataSet = vecDataSet.at(i);

        if (useDataRunListForMC && i>=nData) nameOutput = Form("%s/%s/PrimaryTrackQA/%s-%s_PrimaryTrackQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),fDataSet.Data(),vecDataSet.at(0).Data());
        else nameOutput = Form("%s/%s/PrimaryTrackQA/%s_PrimaryTrackQARunwise.root",fCutSelection.Data(),fEnergyFlag.Data(),fDataSet.Data());

        TFile* fOutput = new TFile(nameOutput,"RECREATE");
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Output file: " << nameOutput << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog[i] << "Output file: " << nameOutput << endl;
        fLog[i] << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        for(Int_t h=0; h<(Int_t) vecHistos[i].size(); h++) WriteHistogram(((TH1D*) vecHistos[i].at(h)));

        DeleteVecTH1D(vecESD_PrimaryNegPions_Pt[i]);
        DeleteVecTH1D(vecESD_PrimaryPosPions_Pt[i]);
        DeleteVecTH1D(vecESD_PrimaryNegPions_Phi[i]);
        DeleteVecTH1D(vecESD_PrimaryPosPions_Phi[i]);
        DeleteVecTH1D(vecESD_PrimaryNegPions_Eta[i]);
        DeleteVecTH1D(vecESD_PrimaryPosPions_Eta[i]);
        //AfterQA=Cut Number Pion Cut (Pion Cut folder in Cut number directory); Pre Selection=Main Dir Pion Cut
        //AfterQA
        DeleteVecTH1D(vecIsPionSelected_AfterQA[i]);
        DeleteVecTH1D(vecdEdxCuts_AfterQA[i]);
        //Pre Selection
        DeleteVecTH1D(vecIsPionSelected_PreSel[i]);
        DeleteVecTH1D(vecdEdxCuts_PreSel[i]);
        //MC histograms
        DeleteVecTH1D(vecMC_AllPosPions_Pt[i]);
        DeleteVecTH1D(vecMC_AllNegPions_Pt[i]);
        DeleteVecTH1D(vecMC_PosPionsFromNeutralMeson_Pt[i]);
        DeleteVecTH1D(vecMC_NegPionsFromNeutralMeson_Pt[i]);
        //True histograms
        DeleteVecTH1D(vecESD_TruePosPion_Pt[i]);
        DeleteVecTH1D(vecESD_TrueNegPion_Pt[i]);
        DeleteVecTH1D(vecESD_TruePosPionFromNeutralMeson_Pt[i]);
        DeleteVecTH1D(vecESD_TrueNegPionFromNeutralMeson_Pt[i]);
        //---------------------------------------------Projections-----------------------------------------------------------------------
        DeleteVecTH1D(vecESD_PrimaryNegPions_ClsTPC_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPosPions_ClsTPC_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_DCAxy_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_DCAz_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt[i]);
        //AfterQA
        DeleteVecTH1D(vecPion_ITS_after_AfterQA_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_after_AfterQA_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_after_AfterQA_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_after_AfterQA_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_AfterQA_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_AfterQA_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_AfterQA_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_AfterQA_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_AfterQA_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_AfterQA_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_AfterQA_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_AfterQA_HighPt_ProjPt[i]);
        DeleteVecTH1D(vechTrack_DCAxy_Pt_after_AfterQA_ProjPt[i]);
        DeleteVecTH1D(vechTrack_DCAz_Pt_after_AfterQA_ProjPt[i]);
        DeleteVecTH1D(vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt[i]);
        //Pre Selection
        DeleteVecTH1D(vecPion_ITS_before_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_before_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_before_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_before_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_before_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_before_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_before_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_before_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_before_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_before_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_before_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_before_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vechTrack_DCAxy_Pt_before_PreSel_ProjPt[i]);
        DeleteVecTH1D(vechTrack_DCAz_Pt_before_PreSel_ProjPt[i]);
        DeleteVecTH1D(vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_after_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_after_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_after_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_ITS_after_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdx_after_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_PreSel_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_PreSel_LowPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_PreSel_MidPt_ProjPt[i]);
        DeleteVecTH1D(vecPion_TOF_after_PreSel_HighPt_ProjPt[i]);
        DeleteVecTH1D(vechTrack_DCAxy_Pt_after_PreSel_ProjPt[i]);
        DeleteVecTH1D(vechTrack_DCAz_Pt_after_PreSel_ProjPt[i]);
        DeleteVecTH1D(vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt[i]);
        //--------------------------------------------------------------------------------------------------------
        /*WriteHistogramTH2DVec(fOutput,vecESD_PrimaryPions_TPCdEdx_ProjPt[i],"ESD_PrimaryPions_TPCdEdx_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecESD_PrimaryPions_TPCdEdxSignal_ProjPt[i],"ESD_PrimaryPions_TPCdEdxSignal_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_ITS_after_AfterQA_ProjPt[i],"Pion_ITS_after_AfterQA_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_dEdx_after_AfterQA_ProjPt[i],"Pion_dEdx_after_AfterQA_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_dEdxSignal_after_AfterQA_ProjPt[i],"Pion_dEdxSignal_after_AfterQA_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_TOF_after_AfterQA_ProjPt[i],"Pion_TOF_after_AfterQA_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_ITS_before_PreSel_ProjPt[i],"Pion_ITS_before_PreSel_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_dEdx_before_PreSel_ProjPt[i],"Pion_dEdx_before_PreSel_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_dEdxSignal_before_PreSel_ProjPt[i],"Pion_dEdxSignal_before_PreSel_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_TOF_before_PreSel_ProjPt[i],"Pion_TOF_before_PreSel_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_ITS_after_PreSel_ProjPt[i],"Pion_ITS_after_PreSel_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_dEdx_after_PreSel_ProjPt[i],"Pion_dEdx_after_PreSel_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_dEdxSignal_after_PreSel_ProjPt[i],"Pion_dEdxSignal_after_PreSel_ProjPt");
        WriteHistogramTH2DVec(fOutput,vecPion_TOF_after_PreSel_ProjPt[i],"Pion_TOF_after_PreSel_ProjPt");*/
        fOutput->Write();
        fOutput->Close();
        delete fOutput;
        fLog[i].close();
    }
    
    delete[] vecESD_PrimaryNegPions_Pt;
    delete[] vecESD_PrimaryPosPions_Pt;
    delete[] vecESD_PrimaryNegPions_Phi;
    delete[] vecESD_PrimaryPosPions_Phi;
    delete[] vecESD_PrimaryNegPions_Eta;
    delete[] vecESD_PrimaryPosPions_Eta;
    //AfterQA=Cut Number Pion Cut (Pion Cut folder in Cut number directory); Pre Selection=Main Dir Pion Cut
    //AfterQA
    delete[] vecIsPionSelected_AfterQA;
    delete[] vecdEdxCuts_AfterQA;
    //Pre Selection
    delete[] vecIsPionSelected_PreSel;
    delete[] vecdEdxCuts_PreSel;
    //MC histograms
    delete[] vecMC_AllPosPions_Pt;
    delete[] vecMC_AllNegPions_Pt;
    delete[] vecMC_PosPionsFromNeutralMeson_Pt;
    delete[] vecMC_NegPionsFromNeutralMeson_Pt;
    //True histograms
    delete[] vecESD_TruePosPion_Pt;
    delete[] vecESD_TrueNegPion_Pt;
    delete[] vecESD_TruePosPionFromNeutralMeson_Pt;
    delete[] vecESD_TrueNegPionFromNeutralMeson_Pt;
    //---------------------------------------------Projections-----------------------------------------------------------------------
    delete[] vecESD_PrimaryNegPions_ClsTPC_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPosPions_ClsTPC_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_DCAxy_ProjPt;//Pt was y
    delete[] vecESD_PrimaryPions_DCAz_ProjPt;//Pt was y
    delete[] vecESD_PrimaryPions_TPCdEdx_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_TPCdEdx_LowPt_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_TPCdEdx_MidPt_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_TPCdEdx_HighPt_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_TPCdEdxSignal_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_TPCdEdxSignal_LowPt_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_TPCdEdxSignal_MidPt_ProjPt;//Pt was x
    delete[] vecESD_PrimaryPions_TPCdEdxSignal_HighPt_ProjPt;//Pt was x
    //AfterQA
    delete[] vecPion_ITS_after_AfterQA_ProjPt;//Pt was x
    delete[] vecPion_ITS_after_AfterQA_LowPt_ProjPt;//Pt was x
    delete[] vecPion_ITS_after_AfterQA_MidPt_ProjPt;//Pt was x
    delete[] vecPion_ITS_after_AfterQA_HighPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_AfterQA_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_AfterQA_LowPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_AfterQA_MidPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_AfterQA_HighPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_AfterQA_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_AfterQA_LowPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_AfterQA_MidPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_AfterQA_HighPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_AfterQA_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_AfterQA_LowPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_AfterQA_MidPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_AfterQA_HighPt_ProjPt;//Pt was x
    delete[] vechTrack_DCAxy_Pt_after_AfterQA_ProjPt;//Pt was y
    delete[] vechTrack_DCAz_Pt_after_AfterQA_ProjPt;//Pt was y
    delete[] vechTrack_NFindCls_Pt_TPC_after_AfterQA_ProjPt;//Pt was x
    //Pre Selection
    delete[] vecPion_ITS_before_PreSel_ProjPt;//Pt was x
    delete[] vecPion_ITS_before_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_ITS_before_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_ITS_before_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_before_PreSel_ProjPt;//Pt was x
    delete[] vecPion_dEdx_before_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_before_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_before_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_before_PreSel_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_before_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_before_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_before_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_before_PreSel_ProjPt;//Pt was x
    delete[] vecPion_TOF_before_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_before_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_before_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vechTrack_DCAxy_Pt_before_PreSel_ProjPt;//Pt was y
    delete[] vechTrack_DCAz_Pt_before_PreSel_ProjPt;//Pt was y
    delete[] vechTrack_NFindCls_Pt_TPC_before_PreSel_ProjPt;//Pt was x
    delete[] vecPion_ITS_after_PreSel_ProjPt;//Pt was x
    delete[] vecPion_ITS_after_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_ITS_after_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_ITS_after_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_PreSel_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_dEdx_after_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_PreSel_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_dEdxSignal_after_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_PreSel_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_PreSel_LowPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_PreSel_MidPt_ProjPt;//Pt was x
    delete[] vecPion_TOF_after_PreSel_HighPt_ProjPt;//Pt was x
    delete[] vechTrack_DCAxy_Pt_after_PreSel_ProjPt;//Pt was y
    delete[] vechTrack_DCAz_Pt_after_PreSel_ProjPt;//Pt was y
    delete[] vechTrack_NFindCls_Pt_TPC_after_PreSel_ProjPt;//Pt was x
    
    cout << "Done with PrimaryTrackQA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    return;
}
