/*******************************************************************************
 ******  provided by Gamma Conversion Group, PWGGA,                        *****
 ******     Daniel Muehlheim, d.muehlheim@cern.ch                          *****
 *******************************************************************************/


#include "QA.h"

void ClusterQA(
                Int_t nSetsIn,                          // number of data sets to be analysed
                TString fEnergyFlag,                    // energy flag
                TString* DataSets,                      // technical names of data sets for output
                TString* plotDataSets,                  // labels of data sets in plots
                TString* pathDataSets,                  // path for data sets
                Int_t mode              = 2,            // standard mode for analysis
                Int_t cutNr             = -1,           // if -1: you have to choose number at runtime
                Int_t doExtQA           = 2,            // 0: switched off, 1: normal extQA / trackmatching, 2: with Cell level plots / trackmatching+extQA
                TString suffix          = "eps",        // output format of plots
                TString labelData       = "Data",       // Label for data
                Bool_t addSubfolder     = kFALSE,       // flag to enable subdirectory creation for primary cut
                Bool_t runMergedClust   = kFALSE,        // flag to switch to merged cluster cuts instead of standard cluster cut
                Int_t useCellQAcutsExt  = 0,        //
                Double_t* arrSetQAEnergy= NULL,
                Double_t* arrSetQATime  = NULL,
                Double_t* arrSetQAHotCells1D= NULL,
                Double_t* arrSetMin2D   = NULL,
                Double_t* arrSetMax2D   = NULL
              )
{
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "ClusterQA" << endl;
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
    const Int_t nSets               = nSetsIn;

    // set correct mode and enable respective flags
    Int_t fMode                     = mode;
    Bool_t isPCMCalo                = kTRUE;
    Bool_t isMerged                 = kFALSE;
    // mode:    0 // new output PCM-PCM
    //          1 // new output PCM dalitz
    //          2 // new output PCM-EMCal
    //          3 // new output PCM-PHOS
    //          4 // new output EMCal-EMCal
    //          5 // new output PHOS-PHOS
    //          9 // old output PCM-PCM
    //          10 // merged EMCal
    //          11 // merged PHOS
    if(fMode == 0 || fMode == 1 || fMode == 9){
        cout << "Returning, given mode contains no calo information: " << fMode << endl; return;}
    if(fMode != 2 && fMode != 3)
        isPCMCalo                   = kFALSE;
    if(fMode == 10 || fMode == 11)
        isMerged                    = kTRUE;

    TString fDetectionProcess       = ReturnFullTextReconstructionProcess(fMode);
    TString fDate                   = ReturnDateString();
    TString fTextMeasurement        = Form("#pi^{0} #rightarrow #gamma#gamma");
    TString fTextMeasurementEta     = Form("#eta #rightarrow #gamma#gamma");
    TString fTextMeasurementMeson[2]= {fTextMeasurement, fTextMeasurementEta};

    TString fCollisionSystem        = ReturnFullCollisionsSystem(fEnergyFlag);
    if (fCollisionSystem.CompareTo("") == 0){
        cout << "No correct collision system specification, has been given" << endl;
        return;
    }

    const Int_t maxSets             = 12;
    //nSets == 0 is always data!

    if(nSets>maxSets){
        cout << "Maximum hardcoded number of Data Sets in ClusterQA: " << maxSets << endl;
        cout << "You have chosen: " << nSets << ", returning!" << endl;
        return;
    }

    Color_t colorCompare[maxSets]   = {kBlack, kRed+1, kMagenta+2, 807, 800, kGreen+2, kCyan+2, kBlue+1, kOrange+2, kAzure, kViolet, kGray+1};
    TString nameMainDir[maxSets];


    TString* fCutSelection          = new TString[nSets];
    TString* fEventCutSelection     = new TString[nSets];
    TString* fGammaCutSelection     = new TString[nSets];
    TString* fClusterCutSelection   = new TString[nSets];
    TString* fElectronCutSelection  = new TString[nSets];
    TString* fMesonCutSelection     = new TString[nSets];
    TString* fMClusterCutSelection  = new TString[nSets];

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Determine which cut to process ***************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    vector <TString> cutsTemp;
    map<TString,Int_t> mapCuts;

    for(Int_t i=0; i<nSets; i++){
        TFile* fFile        = new TFile(pathDataSets[i].Data(),"READ");
        if(fFile->IsZombie()){
            cout << "ERROR: File " << pathDataSets[i].Data() << " could not be openend! Returning..." << endl;
            return;
        } else {
            cout << "Processing file: " << pathDataSets[i].Data();
            TKey *key;
            TIter next(fFile->GetListOfKeys());
            while ((key=(TKey*)next())){
                cout << Form(" - found TopDir: %s",key->GetName());
                nameMainDir[i]  = key->GetName();
                if (nameMainDir[i].Contains("Gamma")) break;
            }
            cout << endl;
            if(nameMainDir[i].IsNull() || !nameMainDir[i].BeginsWith("Gamma")){
                cout << "ERROR, Unable to obtain valid name of MainDir:|" << nameMainDir[i].Data() << "|, running in mode: " << fMode << endl; return;
            }

            TList *listInput    = (TList*)fFile->Get(nameMainDir[i].Data());
            listInput->SetOwner(kTRUE);
            if(!listInput) {cout << "ERROR: Could not find main dir: " << nameMainDir[i].Data() << " in file! Returning..." << endl; return;}
            for(Int_t j = 0; j<listInput->GetSize(); j++){
                TList *listCuts     = (TList*)listInput->At(j);
                TString nameCuts    = listCuts->GetName();
                if(nameCuts.BeginsWith("Cut Number")){
                nameCuts.Replace(0,11,"");
                if(i==0) {
                    cutsTemp.push_back(nameCuts);
                    mapCuts[nameCuts]   =1;
                } else if( mapCuts.find(nameCuts) != mapCuts.end() )
                    mapCuts[nameCuts]   += 1;
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

    TString TrackMatchingEtaMax[10] = {"0","0.008", "0.012", "0.016", "0.018", "0.02", "0.022", "p_{T} dependent", "p_{T} dependent", "p_{T} dependent"};
    TString TrackMatchingPhiMin[10] = {"0","-0.03", "-0.05", "-0.09", "-0.11", "-0.13", "-0.15", "p_{T} dependent", "p_{T} dependent", "p_{T} dependent"};
    TString TrackMatchingPhiMax[10] = {"0","0.03", "0.04", "0.06", "0.07", "0.08", "0.1", "p_{T} dependent", "p_{T} dependent", "p_{T} dependent"};

    Bool_t isTrackMatching          = kTRUE;
    Int_t trackMatch                = 0;
    for(Int_t i=0; i<nSets; i++){
        fCutSelection[i]            = cuts.at(cutNr);
        fEventCutSelection[i]       = "";
        fGammaCutSelection[i]       = "";
        fClusterCutSelection[i]     = "";
        fMClusterCutSelection[i]    = "";
        fElectronCutSelection[i]    = "";
        fMesonCutSelection[i]       = "";
        if (!isMerged){
            ReturnSeparatedCutNumberAdvanced(fCutSelection[i], fEventCutSelection[i], fGammaCutSelection[i], fClusterCutSelection[i], fElectronCutSelection[i], fMesonCutSelection[i], fMode);
        } else {
            ReturnSeparatedCutNumberAdvanced(fCutSelection[i], fEventCutSelection[i], fClusterCutSelection[i], fMClusterCutSelection[i], fElectronCutSelection[i], fMesonCutSelection[i], fMode);
        }
        // reset cluster cut to merged cluster cut if running QA for merged clus
        if (runMergedClust){
            fClusterCutSelection[i] = fMClusterCutSelection[i];
        }
        TString trackMatchingCut(fClusterCutSelection[i](GetClusterTrackMatchingCutPosition(fClusterCutSelection[i]),1));
        trackMatch                  = trackMatchingCut.Atoi();
        if(trackMatch == 0 || runMergedClust ){
            cout << "INFO: TrackMatching cut found to be '0' for '" << DataSets[i] << "', deactivating generation of track matching histograms!" << endl;
            isTrackMatching         = kFALSE;
        }
    }

    TString centralityString        = GetCentralityString(fEventCutSelection[0]);
    if (centralityString.CompareTo("pp")!=0 && !centralityString.Contains("0-100%") ){
      fCollisionSystem              = Form("%s %s", centralityString.Data(), fCollisionSystem.Data());
      cout << "added cent to Collision system string" << endl;
    }

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    cout << "Obtaining trigger - ";
    TString* fTrigger               = new TString[nSets];
    for(Int_t iT=0;iT<nSets;iT++){
        fTrigger[iT]                = "";
    }
    TString fTriggerCut             = fEventCutSelection[0](3,2);
    fTrigger[0]                     = AnalyseSpecialTriggerCut(fTriggerCut.Atoi(), DataSets[0].Data()); //fTrigger[0]+=" ";
    cout  << "'" << fTrigger[0].Data() << "' - was found!" << endl;
    if(fTrigger[0].Contains("not defined")){
        fTrigger[0]                 = "";
        cout << "INFO: Trigger cut not defined!" << endl;
    }
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    //*****************************************************************************************************
    //************************** Define output directories*************************************************
    //*****************************************************************************************************
    TString outputDir           = "";
    TString outputDirRootFiles  = "";
    if (!runMergedClust){
        outputDir               = Form("%s/%s/ClusterQA/%s",cuts.at(cutNr).Data(),fEnergyFlag.Data(),suffix.Data());
        outputDirRootFiles      = Form("%s/%s/ClusterQA",cuts.at(cutNr).Data(),fEnergyFlag.Data());
    } else {
        outputDir               = Form("%s/%s/MergedClusterQA/%s",cuts.at(cutNr).Data(),fEnergyFlag.Data(),suffix.Data());
        outputDirRootFiles      = Form("%s/%s/MergedClusterQA",cuts.at(cutNr).Data(),fEnergyFlag.Data());
    }
    if(addSubfolder)
        outputDir       +=Form("/%s",DataSets[0].Data());

    gSystem->Exec("mkdir -p "+outputDir);
    gSystem->Exec("mkdir -p "+outputDir+"/Comparison/Ratios");
    if (doExtQA > 0 && !runMergedClust){
        gSystem->Exec("mkdir -p "+outputDir+"/Cells/Detailed/Good");
    }
    //*****************************************************************************************************
    //************************** Set proper cluster nomenclature ******************************************
    //*****************************************************************************************************

    TString calo            = "";
    TString fClusters       = "";
    Int_t iCalo             = 0;
    Int_t nCaloModules      = 0;
    Int_t nCaloCells        = 0;

    if(fClusterCutSelection[0].BeginsWith('1')){
        calo                = "EMCal";
        iCalo               = 1;
        fClusters           = Form("%s clusters", calo.Data());
        nCaloModules        = 10;
        nCaloCells          = 11520;
        if(DataSets[0].Contains("LHC10")){
            nCaloModules    = 4;
            nCaloCells      = 4608;
        }
    } else if(fClusterCutSelection[0].BeginsWith('2')){
        calo                = "PHOS";
        iCalo               = 2;
        fClusters           = Form("%s clusters", calo.Data());
        nCaloModules        = 5;
        nCaloCells          = 10700;
        if(DataSets[0].Contains("LHC10")){
            nCaloModules    = 5;
            nCaloCells      = 10700;
        }
        if(DataSets[0].Contains("LHC15") || DataSets[0].Contains("LHC16") || DataSets[0].Contains("LHC17")|| DataSets[0].Contains("LHC18")){
            nCaloModules    = 5;
            nCaloCells      = 14336;
        }
    } else if(fClusterCutSelection[0].BeginsWith('3')){
        calo                = "DCal";
        iCalo               = 1;
        fClusters           = Form("%s clusters", calo.Data());
        nCaloModules        = 6;
        nCaloCells          = 19000;
    } else {
        cout << "No correct calorimeter type found: " << calo.Data() << ", returning..." << endl; return;
    }

    //*****************************************************************************************************
    //******************************* create log file *****************************************************
    //*****************************************************************************************************

    fstream fLog;
    fLog.open(Form("%s/A-%s.log",outputDirRootFiles.Data(),DataSets[0].Data()), ios::out);
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    for(Int_t i=0; i<nSets; i++) fLog << "Using file: " << pathDataSets[i].Data() << endl;
    fLog << "Energy: " << fEnergyFlag.Data() << endl;
    fLog << fCollisionSystem.Data() << endl;
    fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    fLog << "processed cut: " << cuts.at(cutNr).Data() << endl;
    fLog << calo.Data() << ", Modules: " << nCaloModules << ", Cells: " << nCaloCells << endl;

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Histograms for MC/Data ***********************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    // Binning for cluster energy histograms

    const Int_t fNBinsClusterPt                     = 60;
    Double_t fBinsClusterPt[fNBinsClusterPt+1]      = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                                        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                                        2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8,
                                                        4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8,
                                                        6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.4, 7.8, 8.2, 8.6,
                                                        9.0, 9.5, 10, 11, 12, 14, 16, 18, 20, 25, 30};

    TH1D* fDeltaPt = new TH1D("deltaPt","",fNBinsClusterPt,fBinsClusterPt);
    for(Int_t iPt=1;iPt<fNBinsClusterPt+1;iPt++){
        fDeltaPt->SetBinContent(iPt,fBinsClusterPt[iPt]-fBinsClusterPt[iPt-1]);
        fDeltaPt->SetBinError(iPt,0);
    }

    const Int_t fNBins2760GeVPt                     = 16;
    Double_t fBins2760GeVPt[fNBins2760GeVPt+1]      = { 0.6, 0.8, 1.0,
                                                        1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 3.0, 3.5,
                                                        4.0, 5.0, 6.0, 8.0};

    const Int_t fNBins8000GeVPt                     = 26;
    Double_t fBins8000GeVPt[fNBins8000GeVPt+1]      = { 0.6, 0.8, 1.0,
                                                        1.2, 1.4, 1.6, 1.8, 2.0,
                                                        2.2, 2.4, 2.6, 2.8, 3.0,
                                                        3.2, 3.4, 3.6, 3.8, 4.0,
                                                        4.5, 5.0, 5.5, 6.0, 7.0,
                                                        8.0, 10.0, 12.0, 16.0};
    const Int_t fNBins7000GeVPt                     = 27;
    Double_t fBins7000GeVPt[fNBins7000GeVPt+1]      = { 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,
                                                        1.8, 2.0, 2.2, 2.4, 2.6,
                                                        2.8, 3.0, 3.2, 3.4, 3.6,
                                                        3.8, 4.0, 4.5, 5.0, 5.5,
                                                        6.0, 7.0, 8.0, 10.0, 12.0,
                                                        16.0, 20.0};

    const Int_t fNBins8000GeVTriggerPt                          = 31;
    Double_t fBins8000GeVTriggerPt[fNBins8000GeVTriggerPt+1]    = { 0.6, 0.8, 1.0,
                                                                    1.2, 1.4, 1.6, 1.8, 2.0,
                                                                    2.2, 2.4, 2.6, 2.8, 3.0,
                                                                    3.2, 3.4, 3.6, 3.8, 4.0,
                                                                    4.5, 5.0, 5.5, 6.0, 7.0,
                                                                    8.0, 10.0, 12.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0};

    Int_t fNBinsAnalysisPt          = 0;
    Double_t* fBinsAnalysisPt       = 0x0;
    if(DataSets[0].Contains("LHC11a") || DataSets[0].Contains("LHC13g")){
        fNBinsAnalysisPt            = fNBins2760GeVPt;
        fBinsAnalysisPt             = new Double_t[fNBinsAnalysisPt+1];
        for(Int_t iB=0; iB<fNBinsAnalysisPt+1; iB++) {
            fBinsAnalysisPt[iB]     = fBins2760GeVPt[iB];
        }
    } else if(CheckForData8TeV(DataSets[0])){
        fNBinsAnalysisPt            = fNBins8000GeVPt;
        fBinsAnalysisPt             = new Double_t[fNBinsAnalysisPt+1];
        for(Int_t iB=0; iB<fNBinsAnalysisPt+1; iB++) {
            fBinsAnalysisPt[iB]     = fBins8000GeVPt[iB];
        }
    } else if(CheckForTriggerData8TeV(DataSets[0])){
        fNBinsAnalysisPt            = fNBins8000GeVTriggerPt;
        fBinsAnalysisPt             = new Double_t[fNBinsAnalysisPt+1];
        for(Int_t iB=0; iB<fNBinsAnalysisPt+1; iB++) {
            fBinsAnalysisPt[iB] = fBins8000GeVTriggerPt[iB];
        }
    } else if(CheckForData7TeV(DataSets[0],kTRUE)){
        fNBinsAnalysisPt            = fNBins7000GeVPt;
        fBinsAnalysisPt             = new Double_t[fNBinsAnalysisPt+1];
        for(Int_t iB=0; iB<fNBinsAnalysisPt+1; iB++) {
            fBinsAnalysisPt[iB]     = fBins7000GeVPt[iB];
        }
    } else{
        cout << "WARNING: No Binning loaded, nothing defined for " << DataSets[0].Data() << ". Loading standard 0.6-30 GeV..." << endl;
        fLog << "WARNING: No Binning loaded, nothing defined for " << DataSets[0].Data() << ". Loading standard 0.6-30 GeV..." << endl;
        fNBinsAnalysisPt            = fNBins8000GeVTriggerPt;
        fBinsAnalysisPt             = new Double_t[fNBinsAnalysisPt+1];
        for(Int_t iB=0; iB<fNBinsAnalysisPt+1; iB++) {
            fBinsAnalysisPt[iB]     = fBins8000GeVTriggerPt[iB];
        }
    }

    cout << "Number of Bins: " << fNBinsAnalysisPt << ", analyzing pT-Bins: " << endl;
    fLog << "Number of Bins: " << fNBinsAnalysisPt << ", analyzing pT-Bins: " << endl;
    for(Int_t iB=0; iB<fNBinsAnalysisPt; iB++) {
        cout << fBinsAnalysisPt[iB] << "-" << fBinsAnalysisPt[iB+1] << ", ";
        fLog << fBinsAnalysisPt[iB] << "-" << fBinsAnalysisPt[iB+1] << ", ";
    }
    cout << endl;
    fLog << endl;

    //*****************************************************************************************************
    //*****************************************************************************************************
    //************************************* Cell QA *******************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    // Setup of CellQA+defining cuts
    CellQAObj* cellQAData   = 0x0;
    Bool_t doCellQA         = kFALSE;
    Bool_t useGoodRuns      = kFALSE;
    TLine *line             = new TLine();
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(1);
    if(useCellQAcutsExt==1){
        cout << "setting cuts externally" << endl;
        doCellQA                = kTRUE;
        cellQAData              = new CellQAObj();
        cout << arrSetQAEnergy[0] << "\t"<< arrSetQAEnergy[1] << "\t"<< arrSetQAEnergy[2] << "\t"<< arrSetQAEnergy[3] << endl;
        setQAEnergy(cellQAData,arrSetQAEnergy[0],arrSetQAEnergy[1],arrSetQAEnergy[2],arrSetQAEnergy[3]);
        cout << arrSetQATime[0] << "\t"<< arrSetQATime[1] << "\t"<< arrSetQATime[2] << "\t"<< arrSetQATime[3] << endl;
        setQATime(cellQAData,arrSetQATime[0],arrSetQATime[1],arrSetQATime[2],arrSetQATime[3]);
        cout << arrSetQAHotCells1D[0] << "\t"<< arrSetQAHotCells1D[1] << "\t"<< arrSetQAHotCells1D[2] << "\t"<< arrSetQAHotCells1D[3] << endl;
        setQAHotCells1D(cellQAData,arrSetQAHotCells1D[0],arrSetQAHotCells1D[1],arrSetQAHotCells1D[2],arrSetQAHotCells1D[3]);
        const Int_t dim2D       = 9;
        Double_t min2D[dim2D];
        Double_t max2D[dim2D];
        for(Int_t i=0;i<9;i++){
            min2D[i]=arrSetMin2D[i];
            max2D[i]=arrSetMax2D[i];
        }
        setQAHotCells2D(cellQAData,dim2D,min2D,max2D);
    }else{
        if(CheckForData8TeV(DataSets[0])){
            if(!DataSets[0].Contains("-kEMC") && DataSets[0].CompareTo("LHC12")!=0 && DataSets[0].CompareTo("LHC12-periods")!=0){
                doCellQA                = kTRUE;
                cellQAData              = new CellQAObj();
                setQAEnergy(cellQAData,0,0.3,0,0.25);
                setQATime(cellQAData,-0.025E-6,0.05E-6,0.002E-6,0.15E-6);
                const Int_t dim2D       = 9;
                Double_t min2D[dim2D]   = {0,2.5,1.5,0.8,0.4,0.25,0.15,0.1,0.08};
                Double_t max2D[dim2D]   = {0,105,105,105,105,105,105,105,105};
                setQAHotCells2D(cellQAData,dim2D,min2D,max2D);

                setQAHotCells1D(cellQAData,0,4E4,0,1.8);
                if(useGoodRuns && DataSets[0].CompareTo("LHC12")==0){
                    //            const Int_t nGoodCells = ;
                    //            Int_t goodCells[nGoodCells]={129,760,869,1039,1367,1377,};
                    //            FillGoodCells(cellQAData,nGoodCells,goodCells);
                }
                if(DataSets[0].CompareTo("LHC12a")==0) setQAHotCells1D(cellQAData,0,9E3,0,1.8);
                if(DataSets[0].CompareTo("LHC12b")==0) setQAHotCells1D(cellQAData,0,5E3,0,1.8);
                if(DataSets[0].CompareTo("LHC12c")==0) setQAHotCells1D(cellQAData,0,1.5E3,0,1.8);
                if(DataSets[0].CompareTo("LHC12d")==0) setQAHotCells1D(cellQAData,0,1.5E4,0,1.8);
                if(DataSets[0].CompareTo("LHC12f")==0) setQAHotCells1D(cellQAData,0,4E3,0,1.8);
                if(DataSets[0].CompareTo("LHC12g")==0){
                    setQAEnergy(cellQAData,0,0.3,0,0.4);
                    setQATime(cellQAData,-0.1E-6,0.1E-6,-0.002E-6,0.3E-6);
                    setQAHotCells1D(cellQAData,0,50,0,10);
                    min2D[8]=0.15; min2D[7]=0.2; min2D[6]=0.25;min2D[5]=0.3;
                    setQAHotCells2D(cellQAData,dim2D,min2D,max2D);
                }
                if(DataSets[0].CompareTo("LHC12h")==0) {
                    setQATime(cellQAData,0,0.05E-6,0.05E-6,0.2E-6);
                    setQAHotCells1D(cellQAData,0,1E4,0,1.8);
                }
                if(DataSets[0].CompareTo("LHC12i")==0) {
                    setQATime(cellQAData,0,0.07E-6,0.05E-6,0.2E-6);
                    setQAHotCells1D(cellQAData,0,1E3,0,1.8);
                }
            }
            //    } else if(DataSets[0].BeginsWith("LHC12-kEMC")){//CheckForTriggerData8TeV(DataSet)){
            //        doCellQA = kTRUE;
            //        cellQAData = new CellQAObj();
            //        setQAEnergy(cellQAData,0,0.4,0,0.35);
            //        setQATime(cellQAData,-0.1E-6,0.1E-6,0,0.4E-6);
            //        setQAHotCells1D(cellQAData,1E4,2E6,0,10);
            //        const Int_t dim2D= 9;
            //        Double_t min2D[dim2D]={0,4,2.5,1.5,1,0.75,0.5,0.3,0.2};
            //        Double_t max2D[dim2D]={0,105,105,105,105,105,105,105,105};
            //        setQAHotCells2D(cellQAData,dim2D,min2D,max2D);
        } else if(DataSets[0].CompareTo("LHC11a_p4_wSDD")==0){
            doCellQA                = kTRUE;
            cellQAData              = new CellQAObj();
            setQAEnergy(cellQAData,0,0.25,0,0.25);
            setQATime(cellQAData,-0.005E-6,0.025E-6,0.005E-6,0.04E-6);
            setQAHotCells1D(cellQAData,0,6E3,0,1.5);
            const Int_t dim2D       = 9;
            Double_t min2D[dim2D]   = {0,2.5,1.5,0.8,0.4,0.25,0.15,0.1,0.08};
            Double_t max2D[dim2D]   = {0,105,105,105,105,105,105,105,105};
            setQAHotCells2D(cellQAData,dim2D,min2D,max2D);
        } else if(DataSets[0].CompareTo("LHC11a_p4_wSDD-kEMC1")==0){
            doCellQA                = kTRUE;
            cellQAData              = new CellQAObj();
            setQAEnergy(cellQAData,0,0.4,0,0.5);
            setQATime(cellQAData,-0.01E-6,0.025E-6,0,0.04E-6);
            setQAHotCells1D(cellQAData,0,1E3,0,1.5);
            const Int_t dim2D       = 9;
            Double_t min2D[dim2D]   = {0,2.5,1.5,0.8,0.7,0.6,0.5,0.4,0.3};
            Double_t max2D[dim2D]   = {0,105,105,105,105,105,105,105,105};
        } else if(DataSets[0].CompareTo("LHC13bc")==0){
            doCellQA                = kTRUE;
            cellQAData              = new CellQAObj();
            setQAEnergy(cellQAData,0.065,0.205,0.045,0.23);
            setQATime(cellQAData,0.05E-6,0.245E-6,0.4E-6,0.82E-6);
            setQAHotCells1D(cellQAData,-500,14E3,0,1.5);
            const Int_t dim2D       = 9;
            Double_t min2D[dim2D]   = {0,2.0,1,0.7,0.6,0.5,0.3,0.3,0.25};
            Double_t max2D[dim2D]   = {0,105,105,105,105,105,105,105,105};
            setQAHotCells2D(cellQAData,dim2D,min2D,max2D);
        } else if(CheckForData7TeV(DataSets[0])){
            doCellQA                = kTRUE;
            cellQAData              = new CellQAObj();
            setQAEnergy(cellQAData,0,0.2,0,0.2);
            const Int_t dim2D       = 9;
            Double_t min2D[dim2D]   = {0,2.5,1.5,0.8,0.4,0.25,0.15,0.1,0.08};
            Double_t max2D[dim2D]   = {0,105,105,105,105,105,105,105,105};
            setQAHotCells2D(cellQAData,dim2D,min2D,max2D);

            if(DataSets[0].CompareTo("LHC10b_pass4")==0){
                setQAHotCells1D(cellQAData,1E3,7E3,0,4);
                setQATime(cellQAData,6E-7,7E-7,4E-8,6E-8);
            } else if(DataSets[0].CompareTo("LHC10c_pass4")==0){
                setQAHotCells1D(cellQAData,2E3,2E4,0,2);
                setQATime(cellQAData,6E-7,6.8E-7,2E-8,5E-8);
            } else if(DataSets[0].CompareTo("LHC10d_pass4")==0){
                setQAHotCells1D(cellQAData,5E3,2.5E4,0,2);
                setQATime(cellQAData,6E-7,6.8E-7,0,3E-8);
            } else if(DataSets[0].CompareTo("LHC10e_pass4")==0){
                setQAHotCells1D(cellQAData,3E3,3E4,0,2);
                setQATime(cellQAData,5.5E-7,6.5E-7,0,3E-8);
            } else if(DataSets[0].CompareTo("LHC10f_pass4")==0){
                setQAHotCells1D(cellQAData,2E3,1.2E4,0,2);
                setQATime(cellQAData,6E-7,7E-7,8E-8,1.4E-7);
            } else{
                setQAHotCells1D(cellQAData,0,1E5,0,4);
                setQATime(cellQAData,6E-7,7E-7,0,1E-7);
            }
        } else if(DataSets[0].CompareTo("LHC15o")==0){
            doCellQA                = kTRUE;
            cellQAData              = new CellQAObj();
            setQAEnergy(cellQAData,0.06,0.22,0.06,0.22);
            setQATime(cellQAData,-0.005E-6,0.015E-6,0.02E-6,0.05E-6);
            setQAHotCells1D(cellQAData,0.1E6,1.7E6,0,1.7);
            const Int_t dim2D       = 9;
            Double_t min2D[dim2D]   = {0.7,0.6,0.5,0.4,0.3,0.25,0.2,0.2,0.2};
            Double_t max2D[dim2D]   = {1E3,1E3,105,105,105,105,105,105,105};
            setQAHotCells2D(cellQAData,dim2D,min2D,max2D);
        } else{
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "\tWARNING: For DataSet " << DataSets[0].Data() << " there are no cell cuts defined yet!" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "\tWARNING: For DataSet " << DataSets[0].Data() << " there are no cell cuts defined yet!" << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        }
    }

    TString charge[2]       = {"pos","neg"};
    TString meson[2]        = {"Pi0","Eta"};

    Bool_t isMC             = kFALSE;
    CellQAObj* cellQA       = cellQAData;

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Vectors for Histograms ***********************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    std::vector<TH1D*> vecNEvents;
    std::vector<TH1D*> vecClusterTimeBefore;
    std::vector<TH1D*> vecClusterTimeAfter;
    std::vector<TH2D*> vecEtaPhiAfterQA;
    std::vector<TH2D*> vecCellEnergyForComparison;
    std::vector<TH2D*> vecCellTimingForComparison;
    std::vector<TH1D*> vecClusterR;
    std::vector<TH1D*> vecClusterEnergy;
    std::vector<TH1D*> vecClusterM02;
    std::vector<TH1D*> vecClusterM20;
    std::vector<TH1D*> vecClusterDispersion;
    std::vector<TH1D*> vecClusterNCells;
    std::vector<TH1D*> vecMatchingDeltaEtaPhi_matched[2];
    std::vector<TH1D*> vecMatchingDeltaEtaPhi_allTracks[2];
    std::vector<TH1D*> vecMatchingDeltaEtaPhiTracksAll[2][2];
    std::vector<TH1D*> vecMatchingDeltaEtaPhiTracksMatched[2][2];
    std::vector<TH2D*> vecMatchingDeltaEtaPhiTracksAllPt[3][2];
    std::vector<TH1D*> vecConvPhotonEtaPhi_Pi0[2];
    std::vector<TH1D*> vecConvPhotonEtaPhi_Eta[2];
    std::vector<TH1D*> vecClusterEnergyFracCells;
    std::vector<TH1D*> vecClusterIncludedCells;
    std::vector<TH1D*> vecClusterNLM;
    std::vector<TH1D*> vecClusterTimingInCluster;
    std::vector<TH1D*> vecClusterDistanceRowWithin;
    std::vector<TH1D*> vecClusterDistanceColumnWithin;

    const Int_t nEM02_bins = 8;
    Double_t EM02_bins[nEM02_bins]={0.7,1.0,1.5,2.0,4.0,8.0,12.0,16.0};
    std::vector<TH2D*> vecClusterEM02;
    std::vector<TH1D*>* vecClusterEM02_TH1D = new std::vector<TH1D*>[nEM02_bins];

    std::vector<TH2D*> vecESDMother;
    std::vector<TH2D*> vecESDMotherMatched;
    std::vector<TH2D*> vecESDTrueMeson[2];
    std::vector<TH2D*> vecESDTrueMesonMatched[2];
    std::vector<TH2D*> vecESDTrueMesonCaloConvertedPhotonMatched[2];
    std::vector<TH2D*> vecESDTrueMesonCaloConvertedPhoton[2];
    std::vector<TH2D*> vecESDTrueMesonCaloPhoton[2];
    std::vector<TH2D*> vecESDTrueMesonCaloElectron[2];

    Double_t* nEvents       = new Double_t[nSets];
    Bool_t* isTrueContainer = new Bool_t[nSets];
    Int_t minB              = 0;
    Int_t maxB              = 0;

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Looping over DataSets ************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    TCanvas* canvas         = new TCanvas("canvas","",10,10,750,500);  // gives the page size
    TCanvas* canvas2NTEO    = new TCanvas("canvas2NTEO","",10,10,1500,500);  // gives the page size
    TCanvas* canvasJPG      = new TCanvas("canvasJPG","",10,10,1500,1000);  // gives the page size
    TCanvas* cvsQuadratic   = new TCanvas("cvsQuadratic","",10,10,500,500);  // gives the page size

    Double_t leftMargin         = 0.09;
    Double_t rightMargin        = 0.02;
    Double_t topMargin          = 0.04;
    Double_t bottomMargin       = 0.09;
    Double_t topMarginQuad      = 0.014;
    Double_t bottomMarginQuad   = 0.092;
    Double_t leftMarginQuad     = 0.117;
    Double_t rightMarginQuad    = 0.117;
    Double_t leftMargin2NTEO    = 0.075;
    Double_t rightMargin2NTEO   = 0.005;
    Double_t topMargin2NTEO     = 0.005;
    Double_t bottomMargin2NTEO  = 0.085;

    DrawGammaCanvasSettings(canvas,leftMargin,rightMargin,topMargin,bottomMargin);
    DrawGammaCanvasSettings(canvasJPG,leftMargin,rightMargin,topMargin,bottomMargin);
    DrawGammaCanvasSettings(canvas2NTEO,leftMargin2NTEO,rightMargin2NTEO,topMargin2NTEO,bottomMargin2NTEO);
    canvas2NTEO->cd();
    TPad* padL2NTEO          = new TPad("padL2NTEO", "", 0., 0., 0.5, 1.,-1, -1, -2);
    TPad* padR2NTEO          = new TPad("padR2NTEO", "", 0.5, 0., 1., 1.,-1, -1, -2);
    DrawGammaPadSettings( padL2NTEO, leftMargin2NTEO, rightMargin2NTEO, topMargin2NTEO, bottomMargin2NTEO);
    DrawGammaPadSettings( padR2NTEO, leftMargin2NTEO, rightMargin2NTEO, topMargin2NTEO, bottomMargin2NTEO);
    padL2NTEO->Draw();
    padR2NTEO->Draw();



    for(Int_t i=0; i<nSets; i++){

        TFile* fFile = new TFile(pathDataSets[i].Data(),"READ");
        if(fFile->IsZombie()) {
            cout << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
            fLog << "ERROR: ROOT file '" << pathDataSets[i].Data() << "' could not be openend, returning!" << endl;
            return;
        }
        //---------------------------------------------------------------------------------------------------------------
        if(i>0){isMC = kTRUE; cellQA = 0x0;}

        //-------------------------------------------------------------------------------------------------------------------------------
        // reading respective containers
        //---------------------------------------------------------------------------------------------------------------
        TList* TopDir                   = (TList*) fFile->Get(nameMainDir[i].Data());
        if(TopDir == NULL) {cout << "ERROR: TopDir not Found"<<endl; return;}
        else TopDir->SetOwner(kTRUE);
        TList* TopContainer             = (TList*) TopDir->FindObject(Form("Cut Number %s",fCutSelection[i].Data()));
        if(TopContainer == NULL) {cout << "ERROR: " << Form("Cut Number %s",fCutSelection[i].Data()) << " not found in File" << endl; return;}
            else TopContainer->SetOwner(kTRUE);
        TList* ESDContainer             = (TList*) TopContainer->FindObject(Form("%s ESD histograms",fCutSelection[i].Data()));
        if(ESDContainer == NULL) {cout << "ERROR: " << Form("%s ESD histograms",fCutSelection[i].Data()) << " not found in File" << endl; return;}
            else ESDContainer->SetOwner(kTRUE);
        TList* CaloCutsContainer        = (TList*) TopContainer->FindObject(Form("CaloCuts_%s",fClusterCutSelection[i].Data()));
        if(CaloCutsContainer == NULL) {cout << "ERROR: " << Form("CaloCuts_%s",fClusterCutSelection[i].Data()) << " not found in File" << endl; return;}
            else CaloCutsContainer->SetOwner(kTRUE);
        TList* ClusterContainer         = (TList*) TopContainer->FindObject(Form("%s Cluster Output",fCutSelection[i].Data()));
        if(isPCMCalo && ClusterContainer == NULL) {cout << "ERROR: " << Form("%s Cluster Output",fCutSelection[i].Data()) << " not found in File" << endl; return;}
            else if(ClusterContainer) ClusterContainer->SetOwner(kTRUE);
        TList* CaloExtQAContainer       = (TList*) TopContainer->FindObject(Form("CaloExtQA_%s",fClusterCutSelection[i].Data()));
        if(CaloExtQAContainer == NULL) {
            cout << "WARNING: " << Form("CaloExtQA_%s",fClusterCutSelection[i].Data()) << " not found in File, using CaloCuts-Container" << endl;
            cout << "INFO: Setting doExtQA to 1, reducing drastically output size!" << endl;
            doExtQA                     = 1;
            CaloExtQAContainer          = CaloCutsContainer;
        } else CaloExtQAContainer->SetOwner(kTRUE);
        isTrueContainer[i]              = kTRUE;
        TList* TrueContainer            = (TList*) TopContainer->FindObject(Form("%s True histograms",fCutSelection[i].Data()));
        if(TrueContainer == NULL) {isTrueContainer[i] = kFALSE; cout << "INFO: " << Form("%s True histograms",fCutSelection[i].Data()) << " not found in File, processing data?" << endl;}
            else TrueContainer->SetOwner(kTRUE);
        //---------------------------------------------------------------------------------------------------------------
        nEvents[i]                      = 0;
        TH1D* fHistNEvents              = (TH1D*)ESDContainer->FindObject("NEventsWOWeight");
        if (fHistNEvents){
            cout << "INFO: Output contains event weights" << endl;
        } else {
            fHistNEvents                = (TH1D*)ESDContainer->FindObject("NEvents");
        }

        if(fHistNEvents){
            if (fEnergyFlag.Contains("PbPb") || fEnergyFlag.Contains("pPb"))
                nEvents[i]     = (Double_t) fHistNEvents->GetBinContent(1);
            else
                nEvents[i]     = (Double_t) GetNEvents(fHistNEvents,kFALSE);
        } else cout << "INFO: Object |fHistNEvents| could not be found!" << endl;
        //---------------------------------------------------------------------------------------------------------------
        cout << endl;
        cout << "----------------------------------------------------------------------------" << endl;
        cout << "Processing file: " << pathDataSets[i].Data() << endl;
        cout << "Set: " << plotDataSets[i].Data() << " - NEvents: " << nEvents[i] << endl;
        cout << "----------------------------------------------------------------------------" << endl;
        cout << endl;
        fLog << "----------------------------------------------------------------------------" << endl;
        fLog << "Processing file: " << pathDataSets[i].Data() << endl;
        fLog << "Set: " << plotDataSets[i].Data() << " - NEvents: " << nEvents[i] << endl;
        fLog << "----------------------------------------------------------------------------" << endl;
        //---------------------------------------------------------------------------------------------------------------
        if( nEvents[i] < 1. ){cout << "ERROR: number of accepted events in data set is <1: " << nEvents[i] << "! Returning..." << endl; return;}

        //---------------------------------------------------------------------------------------------------------------
        const char* nameOutput = Form("%s/ClusterQA_%s.root",outputDirRootFiles.Data(),DataSets[i].Data());
        TFile* fOutput = new TFile(nameOutput,"RECREATE");
        cout << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Output file: " << nameOutput << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "Output file: " << nameOutput << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        //---------------------------------------------------------------------------------------------------------------
        // event histo
        if(fHistNEvents){
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                fHistNEvents,"","","# of Entries",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistNEvents, Form("%s/NEvents_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecNEvents.push_back(new TH1D(*fHistNEvents));
        } else cout << "INFO: Object |fHistNEvents| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // cluster distribution on detector before acceptance cuts
        TH2D* fHistClusterEtavsPhiBeforeAcc = (TH2D*)CaloCutsContainer->FindObject(Form("EtaPhi_beforeAcceptance %s", fClusterCutSelection[i].Data()));
        if(fHistClusterEtavsPhiBeforeAcc){
            fHistClusterEtavsPhiBeforeAcc->Sumw2();
            fHistClusterEtavsPhiBeforeAcc->Scale(1/nEvents[i]);
            fHistClusterEtavsPhiBeforeAcc->Scale(1/GetMeanTH2(fHistClusterEtavsPhiBeforeAcc));
            fHistClusterEtavsPhiBeforeAcc->GetZaxis()->SetRangeUser(0,2);
            DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                fHistClusterEtavsPhiBeforeAcc,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),"#phi_{Cluster}","#eta_{Cluster}",1,1);
            SaveCanvasAndWriteHistogram(canvas, fHistClusterEtavsPhiBeforeAcc, Form("%s/EtaVsPhi_%s_beforeAcceptance.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |fHistClusterEtavsPhiBeforeAcc| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // cluster distribution on detector after acceptance cuts
        TH2D* fHistClusterEtavsPhiAfterAcc = (TH2D*)CaloCutsContainer->FindObject(Form("EtaPhi_afterAcceptance %s", fClusterCutSelection[i].Data()));
        if(fHistClusterEtavsPhiAfterAcc){
            fHistClusterEtavsPhiAfterAcc->Sumw2();
            fHistClusterEtavsPhiAfterAcc->Scale(1/nEvents[i]);
            fHistClusterEtavsPhiAfterAcc->Scale(1/GetMeanTH2(fHistClusterEtavsPhiAfterAcc));
            fHistClusterEtavsPhiAfterAcc->GetZaxis()->SetRangeUser(0,2);
            DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                fHistClusterEtavsPhiAfterAcc,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),"#phi_{Cluster}","#eta_{Cluster}",1,1);
            SaveCanvasAndWriteHistogram(canvas, fHistClusterEtavsPhiAfterAcc, Form("%s/EtaVsPhi_%s_afterAcceptance.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |fHistClusterEtavsPhiAfterAcc| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // cluster distribution on detector after cluster QA cuts
        TH2D* fHistClusterEtavsPhiAfterQA = (TH2D*)CaloCutsContainer->FindObject(Form("EtaPhi_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistClusterEtavsPhiAfterQA){
            fHistClusterEtavsPhiAfterQA->Sumw2();
            fHistClusterEtavsPhiAfterQA->Scale(1/nEvents[i]);
            fHistClusterEtavsPhiAfterQA->Scale(1/GetMeanTH2(fHistClusterEtavsPhiAfterQA));
            fHistClusterEtavsPhiAfterQA->GetZaxis()->SetRangeUser(0,2);
            DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                fHistClusterEtavsPhiAfterQA,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),"#phi_{Cluster}","#eta_{Cluster}",1,1);
            SaveCanvasAndWriteHistogram(canvas, fHistClusterEtavsPhiAfterQA, Form("%s/EtaVsPhi_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecEtaPhiAfterQA.push_back(new TH2D(*fHistClusterEtavsPhiAfterQA));
        } else cout << "INFO: Object |fHistClusterEtavsPhiAfterQA| could not be found! Skipping Draw..." << endl;

        //---------------------------------------------------------------------------------------------------------------
        //----------------------------- Cluster properties - timing  ----------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        TGaxis::SetExponentOffset(-0.08, 0.01, "x");
        //---------------------------------------------------------------------------------------------------------------
        // timing distribution before cuts
        TH2D* fHistClusterTimevsEBeforeQA   = (TH2D*)CaloCutsContainer->FindObject(Form("ClusterTimeVsE_beforeClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistClusterTimevsEBeforeQA){
          TH1D* fHistClusterTimeBeforeQA      = (TH1D*)fHistClusterTimevsEBeforeQA->ProjectionX("fHistClusterTimeBeforeQA",1,fHistClusterTimevsEBeforeQA->GetNbinsY());
          GetMinMaxBin(fHistClusterTimeBeforeQA,minB,maxB); minB-=10; maxB+=10;

          if(fHistClusterTimevsEBeforeQA){
              SetXRange(fHistClusterTimevsEBeforeQA,minB,maxB);
              SetZMinMaxTH2(fHistClusterTimevsEBeforeQA,minB,maxB,1,fHistClusterTimevsEBeforeQA->GetNbinsY());
              DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                  fHistClusterTimevsEBeforeQA,"","#it{t}_{Cluster} (s)","#it{E}_{Cluster} (GeV)",1,1,
                                  0.85,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
              SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterTimevsEBeforeQA, Form("%s/ClusterTimeVsE_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
          } else cout << "INFO: Object |fHistClusterTimevsEBeforeQA| could not be found! Skipping Draw..." << endl;

          if(fHistClusterTimeBeforeQA){
              SetXRange(fHistClusterTimeBeforeQA,minB,maxB);
              DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                  fHistClusterTimeBeforeQA,"","#it{t}_{Cluster} (s)","#frac{d#it{N}}{d#it{t}_{Cluster}}",1,1,
                                  0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
              SaveCanvasAndWriteHistogram(canvas, fHistClusterTimeBeforeQA,Form("%s/ClusterTime_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()) );
              vecClusterTimeBefore.push_back(fHistClusterTimeBeforeQA);
          } else cout << "INFO: Object |fHistClusterTimeBeforeQA| could not be found! Skipping Draw..." << endl;
        }
        //---------------------------------------------------------------------------------------------------------------
        // timing distribution after cuts
        TH2D* fHistClusterTimevsEAfterQA    = (TH2D*)CaloCutsContainer->FindObject(Form("ClusterTimeVsE_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistClusterTimevsEAfterQA){
          TH1D* fHistClusterTimeAfterQA       = (TH1D*)fHistClusterTimevsEAfterQA->ProjectionX("fHistClusterTimeAfterQA",1,fHistClusterTimevsEAfterQA->GetNbinsY());
          GetMinMaxBin(fHistClusterTimeAfterQA,minB,maxB); minB-=10; maxB+=10;

          if(fHistClusterTimevsEAfterQA){
              SetXRange(fHistClusterTimevsEAfterQA,minB,maxB);
              SetZMinMaxTH2(fHistClusterTimevsEAfterQA,minB,maxB,1,fHistClusterTimevsEAfterQA->GetNbinsY());
              DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                  fHistClusterTimevsEAfterQA,"","#it{t}_{Cluster} (s)","#it{E}_{Cluster} (GeV)",1,1,
                                  0.85,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
              SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterTimevsEAfterQA, Form("%s/ClusterTimeVsE_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
          } else cout << "INFO: Object |fHistClusterTimevsEAfterQA| could not be found! Skipping Draw..." << endl;

          if(fHistClusterTimeAfterQA){
              SetXRange(fHistClusterTimeAfterQA,minB,maxB);
              DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                  fHistClusterTimeAfterQA,"","#it{t}_{Cluster} (s)","#frac{d#it{N}}{d#it{t}_{Cluster}}",1,1,
                                  0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
              SaveCanvasAndWriteHistogram(canvas, fHistClusterTimeAfterQA, Form("%s/ClusterTime_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
              vecClusterTimeAfter.push_back(fHistClusterTimeAfterQA);
          } else cout << "INFO: Object |fHistClusterTimeAfterQA| could not be found! Skipping Draw..." << endl;
        }
        //---------------------------------------------------------------------------------------------------------------
        //-------------------------- Detailed look at cell timing distributions in clusters------------------------------
        //---------------------------------------------------------------------------------------------------------------
        // only do if doExtQA>1
        if(doExtQA>1){
            //---------------------------------------------------------------------------------------------------------------
            // cell time for cells contained in a cluster after QA versus cluster energy
            TH2D* fHistClusterTimingInCluster = (TH2D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCellsTiming_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterTimingInCluster){
                GetMinMaxBin(fHistClusterTimingInCluster,minB,maxB);
                SetXRange(fHistClusterTimingInCluster,minB,fHistClusterTimingInCluster->GetXaxis()->FindBin(30.));
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                    fHistClusterTimingInCluster,"",
                                    "Cluster Energy (GeV)","Cell Time in Cluster (ns)",1,1,
                                    0.85,0.92,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterTimingInCluster, Form("%s/ClusterIncludedCellsTiming_vs_E_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                //---------------------------------------------------------------------------------------------------------------
                // show cell timing integrated in energy after QA cuts
                TH1D* fHistClusterTimingInClusterY = (TH1D*) fHistClusterTimingInCluster->ProjectionY("fHistClusterTimingInClusterY",1,fHistClusterTimingInCluster->GetNbinsX());
                if(fHistClusterTimingInClusterY){
            //          GetMinMaxBin(fHistClusterTimingInClusterY,minB,maxB);
                    SetXRange(fHistClusterTimingInClusterY,1,fHistClusterTimingInCluster->GetNbinsY());
                    DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                        fHistClusterTimingInClusterY,"","Cell Time in Cluster (ns)","#frac{d#it{N}}{d#it{t_{cell}}}",1,1,
                                        0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, fHistClusterTimingInClusterY, Form("%s/ClusterIncludedCellsTiming_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    vecClusterTimingInCluster.push_back(fHistClusterTimingInClusterY);
                }
                //---------------------------------------------------------------------------------------------------------------
                // show cell timing for different cluster energies after QA cuts
                canvas->cd(); DrawGammaCanvasSettings(canvas, 130.5/(1350.+(1*108.)), (40.5+(1*108.))/(1350.+(1*108.)), topMargin, bottomMargin);
                const Int_t nB      = 6;
                Double_t range[nB]  = {0.6,1.0,2.0,4.0,10.0,16.0};
                Double_t rangeM[nB] = {0.7,1.2,2.4,5.0,16.0,30.0};
                std::vector<TH1D*> vecCellTimingBins;
                for(Int_t iBins=0; iBins<nB; iBins++){
                    TH1D* proj      = (TH1D*) fHistClusterTimingInCluster->ProjectionY(Form("proj_%i",iBins),fHistClusterTimingInCluster->GetXaxis()->FindBin(range[iBins]),
                                                                                       fHistClusterTimingInCluster->GetXaxis()->FindBin(rangeM[iBins]));
                    if(isMC) proj->Scale(1./proj->GetBinContent(proj->FindBin(60.)));
                    else proj->Scale(1./proj->GetBinContent(proj->FindBin(-1.)));
                    //proj->Scale(1./proj->Integral(1,proj->GetNbinsX()));
                    vecCellTimingBins.push_back(proj);
                }
                TLegend *legendRuns = new TLegend(0.9,0.09,0.98,0.94);
                legendRuns->SetNColumns(1);legendRuns->SetFillColor(0);legendRuns->SetLineColor(0);legendRuns->SetTextSize(0.03);legendRuns->SetTextFont(42);
                //GetMinMaxBin(vecCellTimingBins,minB,maxB);
                AdjustHistRange(vecCellTimingBins,2,5,kFALSE);
                for(Int_t h=0; h<(Int_t) vecCellTimingBins.size(); h++){
                    TString draw = (h==0)?"p":"p,same";
                    EditRunwiseHists(((TH1D*) vecCellTimingBins.at(h)), h, "");
                    //SetXRange(((TH1D*) vecCellTimingBins.at(h)), minB, maxB);
                    if(isMC) ((TH1D*) vecCellTimingBins.at(h))->GetXaxis()->SetRangeUser(40,80);
                    else ((TH1D*) vecCellTimingBins.at(h))->GetXaxis()->SetRangeUser(-100,100);
                    ((TH1D*) vecCellTimingBins.at(h))->Draw(draw.Data());
                    legendRuns->AddEntry(((TH1D*) vecCellTimingBins.at(h)),Form("%.01f - %.01f",range[h],rangeM[h]),"p");
                }
                legendRuns->Draw();
                PutProcessLabelAndEnergyOnPlot(0.95, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), 62, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/ClusterIncludedCellsTiming_vs_E_projected_%s.%s", outputDir.Data(), DataSets[i].Data(),suffix.Data()),kFALSE,kTRUE);

                DeleteVecTH1D(vecCellTimingBins);
            } else cout << "INFO: Object |ClusterIncludedCellsTiming_afterClusterQA| could not be found! Skipping Draw..." << endl;

            //---------------------------------------------------------------------------------------------------------------
            // cell timing of cells included in a cluster vs. cell energy
            TH2D* fHistClusterTimingEnergyInCluster = (TH2D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCellsTimingEnergy_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterTimingEnergyInCluster){
                GetMinMaxBin(fHistClusterTimingEnergyInCluster,minB,maxB);
                SetXRange(fHistClusterTimingEnergyInCluster,minB,maxB);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kTRUE,kFALSE,kTRUE,
                                    fHistClusterTimingEnergyInCluster,"",
                                    "Cell Energy in Cluster (GeV)","Cell Time in Cluster (ns)",1,1,
                                    0.85,0.92,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterTimingEnergyInCluster, Form("%s/ClusterIncludedCellsTiming_vs_CellEnergy_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |ClusterIncludedCellsTimingEnergy_afterClusterQA| could not be found! Skipping Draw..." << endl;

            //---------------------------------------------------------------------------------------------------------------
            // distance of cluster to other clusters withing the same time window in the same event - reflects cluster overlap
            TH2D* fHistClusterDistanceTo_within = (TH2D*)CaloExtQAContainer->FindObject(Form("ClusterDistanceTo_withinTimingCut %s", fClusterCutSelection[i].Data()));
            if(fHistClusterDistanceTo_within){
                AdjustHistRange(fHistClusterDistanceTo_within,1,1.1,kTRUE,1,1);
                GetMinMaxBin(fHistClusterDistanceTo_within,minB,maxB);
                SetXRange(fHistClusterDistanceTo_within,1,maxB+5);
                SetZMinMaxTH2(fHistClusterDistanceTo_within,1,maxB+5,1,fHistClusterDistanceTo_within->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterDistanceTo_within,"",
                                    "#Delta row","#Delta column",1,1,
                                    0.85,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterDistanceTo_within, Form("%s/ClusterDistanceTo_withinTimingCut_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                //---------------------------------------------------------------------------------------------------------------
                // distance of cluster to other clusters withing the same time window in the same event projection to the rows
                TH1D* ClusterDistanceToRow = (TH1D*) fHistClusterDistanceTo_within->ProjectionX("ClusterDistanceToRow",1,fHistClusterDistanceTo_within->GetNbinsY());
                if(ClusterDistanceToRow){
                GetMinMaxBin(ClusterDistanceToRow,minB,maxB);
                SetXRange(ClusterDistanceToRow,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    ClusterDistanceToRow,"","#Delta row","# of entries",1,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, ClusterDistanceToRow, Form("%s/ClusterDistanceTo_withinTimingCut_Row_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecClusterDistanceRowWithin.push_back(ClusterDistanceToRow);
                }
                //---------------------------------------------------------------------------------------------------------------
                // distance of cluster to other clusters withing the same time window in the same event projections to the columns
                TH1D* ClusterDistanceToColumn = (TH1D*) fHistClusterDistanceTo_within->ProjectionY("ClusterDistanceToColumn",1,fHistClusterDistanceTo_within->GetNbinsX());
                if(ClusterDistanceToColumn){
                GetMinMaxBin(ClusterDistanceToColumn,minB,maxB);
                SetXRange(ClusterDistanceToColumn,minB,maxB);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    ClusterDistanceToColumn,"","#Delta column","# of entries",1,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, ClusterDistanceToColumn, Form("%s/ClusterDistanceTo_withinTimingCut_Column_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecClusterDistanceColumnWithin.push_back(ClusterDistanceToColumn);
                }
            } else cout << "INFO: Object |ClusterDistanceTo_withinTimingCut| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // distance of cluster to other clusters withing the time window in the same event - reflects probability for contribution of cells from other bunches to accepted cluster
            TH2D* fHistClusterDistanceTo_outside = (TH2D*)CaloExtQAContainer->FindObject(Form("ClusterDistanceTo_outsideTimingCut %s", fClusterCutSelection[i].Data()));
            if(fHistClusterDistanceTo_outside){
                AdjustHistRange(fHistClusterDistanceTo_outside,1,1.1,kTRUE,1,1);
                GetMinMaxBin(fHistClusterDistanceTo_outside,minB,maxB);
                SetXRange(fHistClusterDistanceTo_outside,1,maxB+5);
                SetZMinMaxTH2(fHistClusterDistanceTo_outside,1,maxB+5,1,fHistClusterDistanceTo_outside->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterDistanceTo_outside,"",
                                    "#Delta row","#Delta column",1,1,
                                    0.85,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterDistanceTo_outside, Form("%s/ClusterDistanceTo_outsideTimingCut_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |ClusterDistanceTo_outsideTimingCut| could not be found! Skipping Draw..." << endl;
        }

        //---------------------------------------------------------------------------------------------------------------
        //---------------------------- Cluster properties energy --------------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        TGaxis::SetExponentOffset(0, 0, "x");
        //---------------------------------------------------------------------------------------------------------------
        // cluster energy before cluster quality cuts
        TH1D* fHistEnergyOfClusterBeforeQA  = (TH1D*)CaloCutsContainer->FindObject(Form("EnergyOfCluster_beforeClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistEnergyOfClusterBeforeQA){
            fHistEnergyOfClusterBeforeQA    = (TH1D*)fHistEnergyOfClusterBeforeQA->Rebin(fNBinsClusterPt, "energyBefore", fBinsClusterPt);
            fHistEnergyOfClusterBeforeQA->Divide(fDeltaPt);
            GetMinMaxBin(fHistEnergyOfClusterBeforeQA,minB,maxB);
            SetXRange(fHistEnergyOfClusterBeforeQA,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                fHistEnergyOfClusterBeforeQA,"","#it{E}_{Cluster} (GeV)","#frac{d#it{N}}{d#it{E}_{Cluster}}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistEnergyOfClusterBeforeQA, Form("%s/EnergyOfCluster_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |fHistEnergyOfClusterBeforeQA| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // cluster energy after cluster quality cuts
        TH1D* fHistEnergyOfClusterAfterQA   = (TH1D*)CaloCutsContainer->FindObject(Form("EnergyOfCluster_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistEnergyOfClusterAfterQA){
            fHistEnergyOfClusterAfterQA     = (TH1D*)fHistEnergyOfClusterAfterQA->Rebin(fNBinsClusterPt, "energyAfter", fBinsClusterPt);
            fHistEnergyOfClusterAfterQA->Sumw2();
            fHistEnergyOfClusterAfterQA->Divide(fDeltaPt);
            GetMinMaxBin(fHistEnergyOfClusterAfterQA,minB,maxB);
            SetXRange(fHistEnergyOfClusterAfterQA,minB,maxB);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kTRUE,kTRUE,kFALSE,
                                fHistEnergyOfClusterAfterQA,"","#it{E}_{Cluster} (GeV)","#frac{d#it{N}}{d#it{E}_{Cluster}}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistEnergyOfClusterAfterQA, Form("%s/EnergyOfCluster_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecClusterEnergy.push_back(new TH1D(*fHistEnergyOfClusterAfterQA));
        } else cout << "INFO: Object |fHistEnergyOfClusterAfterQA| could not be found! Skipping Draw..." << endl;


        //---------------------------------------------------------------------------------------------------------------
        //---------------------------- Cluster properties NCells --------------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        // Cluster nCells before cluster QA cuts
        TH1D* fHistNCellsBeforeQA = (TH1D*)CaloCutsContainer->FindObject(Form("NCellPerCluster_beforeClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistNCellsBeforeQA){
            GetMinMaxBin(fHistNCellsBeforeQA,minB,maxB);
            SetXRange(fHistNCellsBeforeQA,minB,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                 fHistNCellsBeforeQA,"","#it{N}_{cells} per Cluster","",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistNCellsBeforeQA, Form("%s/NCellPerCluster_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |fHistNCellsBeforeQA| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // Cluster nCells after cluster QA cuts
        TH1D* fHistNCellsAfterQA = (TH1D*)CaloCutsContainer->FindObject(Form("NCellPerCluster_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistNCellsAfterQA){
            GetMinMaxBin(fHistNCellsAfterQA,minB,maxB);
            SetXRange(fHistNCellsAfterQA,minB,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistNCellsAfterQA,"","#it{N}_{cells} per Cluster","",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistNCellsAfterQA, Form("%s/NCellPerCluster_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecClusterNCells.push_back(new TH1D(*fHistNCellsAfterQA));
        } else cout << "INFO: Object |fHistNCellsAfterQA| could not be found! Skipping Draw..." << endl;

        if (doExtQA > 1){
            TGaxis::SetExponentOffset(0, 0, "x");
            //---------------------------------------------------------------------------------------------------------------
            // NCells versus cluster energy - reflects size of clusters vs energy
            TH2D* fHistClusterEnergyVsNCells = (TH2D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyVsNCells_afterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterEnergyVsNCells){
                AdjustHistRange(fHistClusterEnergyVsNCells,1,1.1,kTRUE,1,1);
                GetMinMaxBin(fHistClusterEnergyVsNCells,minB,maxB);
                SetXRange(fHistClusterEnergyVsNCells,1,maxB+5);
                SetZMinMaxTH2(fHistClusterEnergyVsNCells,1,maxB+5,1,fHistClusterEnergyVsNCells->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterEnergyVsNCells,"",
                                    "Cluster Energy (GeV)","#it{N}_{Cells} per Cluster",1,1,
                                    0.85,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterEnergyVsNCells, Form("%s/ClusterEnergyVsNCells_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |ClusterEnergyVsNCells_afterQA| could not be found! Skipping Draw..." << endl;
        }

        //---------------------------------------------------------------------------------------------------------------
        //---------------------------- Cluster properties shower shape parameters ---------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        // Cluster M02 before cluster QA (long axis)
        TH1D* fHistM02BeforeQA = (TH1D*)CaloCutsContainer->FindObject(Form("M02_beforeClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistM02BeforeQA){
            GetMinMaxBin(fHistM02BeforeQA,minB,maxB);
            SetXRange(fHistM02BeforeQA,minB-1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistM02BeforeQA,"","#lambda_{0}^{2}","#frac{d#it{N}}{d#lambda_{0}^{2}}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistM02BeforeQA, Form("%s/M02_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |fHistM02BeforeQA| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // Cluster M02 after cluster QA (long axis)
        TH1D* fHistM02AfterQA = (TH1D*)CaloCutsContainer->FindObject(Form("M02_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistM02AfterQA){
            GetMinMaxBin(fHistM02AfterQA,minB,maxB);
            SetXRange(fHistM02AfterQA,minB-1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistM02AfterQA,"","#lambda_{0}^{2}","#frac{d#it{N}}{d#lambda_{0}^{2}}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistM02AfterQA, Form("%s/M02_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecClusterM02.push_back(new TH1D(*fHistM02AfterQA));
        } else cout << "INFO: Object |fHistM02AfterQA| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // Cluster M20 before cluster QA (short axis)
        TH1D* fHistM20BeforeQA = (TH1D*)CaloCutsContainer->FindObject(Form("M20_beforeClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistM20BeforeQA){
            GetMinMaxBin(fHistM20BeforeQA,minB,maxB);
            SetXRange(fHistM20BeforeQA,1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistM20BeforeQA,"","#lambda_{0}^{2}","#frac{d#it{N}}{d#lambda_{1}^{2}}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistM20BeforeQA, Form("%s/M20_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |fHistM20BeforeQA| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // Cluster M20 after cluster QA (short axis)
        TH1D* fHistM20AfterQA = (TH1D*)CaloCutsContainer->FindObject(Form("M20_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistM20AfterQA){
            GetMinMaxBin(fHistM20AfterQA,minB,maxB);
            SetXRange(fHistM20AfterQA,1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistM20AfterQA,"","#lambda_{0}^{2}","#frac{d#it{N}}{d#lambda_{1}^{2}}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistM20AfterQA, Form("%s/M20_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecClusterM20.push_back(new TH1D(*fHistM20AfterQA));
        } else cout << "INFO: Object |fHistM20AfterQA| could not be found! Skipping Draw..." << endl;

        //---------------------------------------------------------------------------------------------------------------
        // Cluster dispersion before cluster QA cuts
        TH1D* fHistDispersionBeforeQA = (TH1D*)CaloCutsContainer->FindObject(Form("Dispersion_beforeClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistDispersionBeforeQA){
            GetMinMaxBin(fHistDispersionBeforeQA,minB,maxB);
            SetXRange(fHistDispersionBeforeQA,minB-1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistDispersionBeforeQA,"","Cluster Dispersion","#frac{d#it{N}}{dDisp}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistDispersionBeforeQA, Form("%s/Dispersion_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << "INFO: Object |fHistDispersionBeforeQA| could not be found! Skipping Draw..." << endl;
        //---------------------------------------------------------------------------------------------------------------
        // Cluster dispersion before cluster QA cuts
        TH1D* fHistDispersionAfterQA = (TH1D*)CaloCutsContainer->FindObject(Form("Dispersion_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistDispersionAfterQA){
            GetMinMaxBin(fHistDispersionAfterQA,minB,maxB);
            SetXRange(fHistDispersionAfterQA,minB-1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistDispersionAfterQA,"","Cluster Dispersion","#frac{d#it{N}}{dDisp}",1,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistDispersionAfterQA, Form("%s/Dispersion_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecClusterDispersion.push_back(new TH1D(*fHistDispersionAfterQA));
        } else cout << "INFO: Object |fHistDispersionAfterQA| could not be found! Skipping Draw..." << endl;

        //---------------------------------------------------------------------------------------------------------------
        //--------------------------------- Detailed shower shape QA  ---------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        //only plot when doExtQA>0
        if(doExtQA>0){
            //---------------------------------------------------------------------------------------------------------------
            // cluster shape: M02 (long axis) vs energy of cluster before QA cuts
            TH2D* fHistClusterEM02BeforeQA = (TH2D*)CaloCutsContainer->FindObject(Form("EVsM02_beforeClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterEM02BeforeQA){
                SetZMinMaxTH2(fHistClusterEM02BeforeQA,1,fHistClusterEM02BeforeQA->GetNbinsX(),1,fHistClusterEM02BeforeQA->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterEM02BeforeQA,"",
                                    "E (GeV)","#lambda_{0}^{2}",1,1,
                                    0.85,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterEM02BeforeQA, Form("%s/EVsM02_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |fHistClusterEM02BeforeQA| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // cluster shape: M02 (long axis) vs energy of cluster after QA cuts
            TH2D* fHistClusterEM02AfterQA = (TH2D*)CaloCutsContainer->FindObject(Form("EVsM02_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterEM02AfterQA){
                SetZMinMaxTH2(fHistClusterEM02AfterQA,1,fHistClusterEM02AfterQA->GetNbinsX(),1,fHistClusterEM02AfterQA->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterEM02AfterQA,"",
                                    "E (GeV)","#lambda_{0}^{2}",1,1,
                                    0.85,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterEM02AfterQA, Form("%s/EVsM02_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecClusterEM02.push_back(new TH2D(*fHistClusterEM02AfterQA));
            } else cout << "INFO: Object |fHistClusterEM02AfterQA| could not be found! Skipping Draw..." << endl;

            if ( !runMergedClust ){ // only do if its not a merged cluster cut
                //---------------------------------------------------------------------------------------------------------------
                // cluster shape: correlatio of M02 (long axis) and M20 (short axis) before cluster QA cuts
                TH2D* fHistClusterM20M02BeforeQA = (TH2D*)CaloCutsContainer->FindObject(Form("M20VsM02_beforeClusterQA %s", fClusterCutSelection[i].Data()));
                if(fHistClusterM20M02BeforeQA){
                    fHistClusterM20M02BeforeQA->GetXaxis()->SetDecimals();
                    fHistClusterM20M02BeforeQA->GetXaxis()->SetLabelOffset(0.01);
                    GetMinMaxBin(fHistClusterM20M02BeforeQA,minB,maxB);
                    SetXRange(fHistClusterM20M02BeforeQA,1,maxB+1);
                    Int_t maxBX = maxB;
                    fHistClusterM20M02BeforeQA->GetYaxis()->SetLabelOffset(0.01);
                    GetMinMaxBinY(fHistClusterM20M02BeforeQA,minB,maxB);
                    SetYRange(fHistClusterM20M02BeforeQA,minB-1,maxB+1);
                    SetZMinMaxTH2(fHistClusterM20M02BeforeQA,1,maxBX+1,minB-1,maxB+1);
                    DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                        fHistClusterM20M02BeforeQA,"",
                                        "#lambda_{1}^{2}","#lambda_{0}^{2}",1,1.2,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterM20M02BeforeQA, Form("%s/M20VsM02_all_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterM20M02BeforeQA| could not be found! Skipping Draw..." << endl;
                //---------------------------------------------------------------------------------------------------------------
                // cluster shape: correlatio of M02 (long axis) and M20 (short axis) after cluster QA cuts
                TH2D* fHistClusterM20M02AfterQA = (TH2D*)CaloCutsContainer->FindObject(Form("M20VsM02_afterClusterQA %s", fClusterCutSelection[i].Data()));
                if(fHistClusterM20M02AfterQA){
                    fHistClusterM20M02AfterQA->GetXaxis()->SetDecimals();
                    fHistClusterM20M02AfterQA->GetXaxis()->SetLabelOffset(0.01);
                    GetMinMaxBin(fHistClusterM20M02AfterQA,minB,maxB);
                    SetXRange(fHistClusterM20M02AfterQA,1,maxB+1);
                    Int_t maxBX = maxB;
                    fHistClusterM20M02AfterQA->GetYaxis()->SetLabelOffset(0.01);
                    GetMinMaxBinY(fHistClusterM20M02AfterQA,minB,maxB);
                    SetYRange(fHistClusterM20M02AfterQA,minB-1,maxB+1);
                    SetZMinMaxTH2(fHistClusterM20M02AfterQA,1,maxBX+1,minB-1,maxB+1);
                    DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                        fHistClusterM20M02AfterQA,"",
                                        "#lambda_{1}^{2}","#lambda_{0}^{2}",1,1.2,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterM20M02AfterQA, Form("%s/M20VsM02_accepted_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterM20M02AfterQA| could not be found! Skipping Draw..." << endl;
                //---------------------------------------------------------------------------------------------------------------
                if(fHistClusterM20M02BeforeQA && fHistClusterM20M02AfterQA){
                    TH2D* fHistClusterM20M02MatchedAfterQA = (TH2D*) fHistClusterM20M02BeforeQA->Clone();
                    fHistClusterM20M02MatchedAfterQA->Add(fHistClusterM20M02AfterQA,-1);
                    DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                        fHistClusterM20M02MatchedAfterQA,"",
                                        "#lambda_{1}^{2}","#lambda_{0}^{2}",1,1.2,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterM20M02MatchedAfterQA, Form("%s/M20VsM02_rejected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    delete fHistClusterM20M02MatchedAfterQA;
                }
            }
        }//end only plot when doExtQA>0

        //---------------------------------------------------------------------------------------------------------------
        //---------------------------- Drawing shower shape parameter M02 vs E for identified particles -----------------
        //---------------------------------------------------------------------------------------------------------------
        if (doExtQA > 1 || isMerged){
            cout << "entered plotting of merged cluster paper plot" << endl;
            TGaxis::SetExponentOffset(0, 0, "x");
            TGaxis::SetExponentOffset(0, 0, "y");
            //---------------------------------------------------------------------------------------------------------------
            if( ((isMerged && TrueContainer) || ClusterContainer) && isMC){
                //---------------------------------------------------------------------------------------------------------------
                // create canvas
                TCanvas* cvsEM02 = new TCanvas("cvs2pads","",2250,1500);
                DrawGammaCanvasSettings(cvsEM02, 0.09, 0.13, 0.01, 0.1);
                cvsEM02->cd();
                TF1 *min = new TF1("min","exp(2.135-0.245*x)",4.95,13.63);
                min->SetLineColor(kBlack);
                min->SetLineStyle(2);
                TF1 *cnst = new TF1("cnst","0.3",13.63,50.05);
                cnst->SetLineColor(kBlack);
                cnst->SetLineStyle(2);
                TF1 *max = new TF1("max","exp(0.353-0.0264*x)-0.524+0.00559*x+21.9/x",9.5,50.05);
                max->SetLineColor(kBlack);
                max->SetLineStyle(2);
                //---------------------------------------------------------------------------------------------------------------
                // validated Gamma plot
                TH2D* fHistClusterTrueGammaEM02 = NULL;
                if (isMerged)
                    fHistClusterTrueGammaEM02   = (TH2D*)TrueContainer->FindObject("TrueClusGammaEM02");
                else
                    fHistClusterTrueGammaEM02   = (TH2D*)ClusterContainer->FindObject("TrueClusGammaEM02");

                if(fHistClusterTrueGammaEM02){
                    fHistClusterTrueGammaEM02->Scale(1/nEvents[i]);
                    DrawAutoGammaHistoPaper2D(fHistClusterTrueGammaEM02,
                                                "",
                                                "E (GeV)",
                                                "#lambda_{0}^{2}",
                                                0,0,0,
                                                1,0.1,2.95,
                                                1,4.95,50.05,0.8,0.65);
                    fHistClusterTrueGammaEM02->GetXaxis()->SetMoreLogLabels();
                    fHistClusterTrueGammaEM02->GetXaxis()->SetLabelOffset(-0.02);
                    SetZMinMaxTH2(fHistClusterTrueGammaEM02,fHistClusterTrueGammaEM02->GetXaxis()->FindBin(4.95),fHistClusterTrueGammaEM02->GetXaxis()->FindBin(50.05),
                                  fHistClusterTrueGammaEM02->GetYaxis()->FindBin(0.1),fHistClusterTrueGammaEM02->GetYaxis()->FindBin(2.95));
                    //fHistClusterTrueGammaEM02->GetZaxis()->SetRangeUser(1E-8,5E-5);
                    fHistClusterTrueGammaEM02->GetZaxis()->SetLabelSize(0.051);
                    fHistClusterTrueGammaEM02->GetXaxis()->SetTickLength(0.05);
                    fHistClusterTrueGammaEM02->DrawCopy("COLZ");
                    min->Draw("same");
                    cnst->Draw("same");
                    max->Draw("same");
                    WriteHistogram(fHistClusterTrueGammaEM02);
                    PutProcessLabelAndEnergyOnPlot(0.55, 0.99, 0.06, "ALICE simulation", fCollisionSystem.Data(), "", 42, 0.03, "");
                    PutProcessLabelAndEnergyOnPlot(0.12, 0.2, 0.06, "#gamma", "", "", 42, 0.03, "");
                    cvsEM02->SetLogx(1); cvsEM02->SetLogy(0); cvsEM02->SetLogz(1); cvsEM02->Update();
                    cvsEM02->SaveAs(Form("%s/EVsM02_TrueGamma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    cvsEM02->Clear();
                } else cout << "INFO: Object |TrueClusGammaEM02| could not be found! Skipping Draw..." << endl;

                //---------------------------------------------------------------------------------------------------------------
                // validated Pi0 plot
                TH2D* fHistClusterTruePi0EM02   = NULL;
                if (isMerged)
                    fHistClusterTruePi0EM02     = (TH2D*)TrueContainer->FindObject("TrueClusPi0EM02");
                else
                    fHistClusterTruePi0EM02     = (TH2D*)ClusterContainer->FindObject("TrueClusPi0EM02");

                if(fHistClusterTruePi0EM02){
                    fHistClusterTruePi0EM02->Scale(1/nEvents[i]);
                    DrawAutoGammaHistoPaper2D(fHistClusterTruePi0EM02,
                                                "",
                                                "E (GeV)",
                                                "#lambda_{0}^{2}",
                                                0,0,0,
                                                1,0.1,2.95,
                                                1,4.95,50.05,0.8,0.65);
                    fHistClusterTruePi0EM02->GetXaxis()->SetMoreLogLabels();
                    fHistClusterTruePi0EM02->GetXaxis()->SetLabelOffset(-0.02);
                    SetZMinMaxTH2(fHistClusterTruePi0EM02,fHistClusterTruePi0EM02->GetXaxis()->FindBin(4.95),fHistClusterTruePi0EM02->GetXaxis()->FindBin(50.05),
                                  fHistClusterTruePi0EM02->GetYaxis()->FindBin(0.1),fHistClusterTruePi0EM02->GetYaxis()->FindBin(2.95));
                    //fHistClusterTruePi0EM02->GetZaxis()->SetRangeUser(1E-8,5E-5);
                    fHistClusterTruePi0EM02->GetZaxis()->SetLabelSize(0.051);
                    fHistClusterTruePi0EM02->GetXaxis()->SetTickLength(0.05);
                    fHistClusterTruePi0EM02->DrawCopy("COLZ");
                    min->Draw("same");
                    cnst->Draw("same");
                    max->Draw("same");

                    WriteHistogram(fHistClusterTruePi0EM02);
                    PutProcessLabelAndEnergyOnPlot(0.55, 0.99, 0.06,  "ALICE simulation", fCollisionSystem.Data(), "", 42, 0.03, "");
                    PutProcessLabelAndEnergyOnPlot(0.3, 0.6, 0.06, "#pi^{0}", "", "", 42, 0.03, "");
                    cvsEM02->SetLogx(1); cvsEM02->SetLogy(0); cvsEM02->SetLogz(1); cvsEM02->Update();
                    cvsEM02->SaveAs(Form("%s/EVsM02_TruePi0_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    cvsEM02->Clear();
                } else cout << "INFO: Object |TrueClusPi0EM02| could not be found! Skipping Draw..." << endl;

                //---------------------------------------------------------------------------------------------------------------
                // paper plot for merged analysis
                if(fHistClusterTrueGammaEM02&&fHistClusterTruePi0EM02){
                    fHistClusterTrueGammaEM02->Add(fHistClusterTruePi0EM02,1);
                    fHistClusterTrueGammaEM02->SetName("TrueClusGammaPi0EM02");
                    SetZMinMaxTH2(fHistClusterTrueGammaEM02,fHistClusterTrueGammaEM02->GetXaxis()->FindBin(4.95),fHistClusterTrueGammaEM02->GetXaxis()->FindBin(50.05),
                                 fHistClusterTrueGammaEM02->GetYaxis()->FindBin(0.1),fHistClusterTrueGammaEM02->GetYaxis()->FindBin(2.95));
//                     fHistClusterTrueGammaEM02->GetZaxis()->SetRangeUser(5E-8,9E-5);
                    fHistClusterTrueGammaEM02->GetZaxis()->SetLabelSize(0.051);
                    fHistClusterTrueGammaEM02->GetXaxis()->SetTickLength(0.05);
                    fHistClusterTrueGammaEM02->DrawCopy("COLZ");

                    min->Draw("same");
                    cnst->Draw("same");
                    max->Draw("same");
                    PutProcessLabelAndEnergyOnPlot(0.5, 0.99, 0.06,  "ALICE simulation", fCollisionSystem.Data(), "", 42, 0.03, "");
                    PutProcessLabelAndEnergyOnPlot(0.28, 0.6, 0.08, "#pi^{0}", "", "", 42, 0.03, "",0);
                    PutProcessLabelAndEnergyOnPlot(0.12, 0.23, 0.08, "#gamma", "", "", 42, 0.03, "",0);
                    cvsEM02->SetLogx(1); cvsEM02->SetLogy(0); cvsEM02->SetLogz(1); cvsEM02->Update();
                    cvsEM02->SaveAs(Form("%s/EVsM02_TrueGamma_Pi0_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    cvsEM02->Clear();
                }
                if (isMerged){
                    //---------------------------------------------------------------------------------------------------------------
                    // validated Eta plot
                    TH2D* fHistClusterTrueEtaEM02   = NULL;
                    if (isMerged)
                        fHistClusterTrueEtaEM02     = (TH2D*)TrueContainer->FindObject("TrueClusEtaEM02");
                    else
                        fHistClusterTrueEtaEM02     = (TH2D*)ClusterContainer->FindObject("TrueClusEtaEM02");

                    if(fHistClusterTrueEtaEM02){
                        fHistClusterTrueEtaEM02->Scale(1/nEvents[i]);
                        DrawAutoGammaHistoPaper2D(fHistClusterTrueEtaEM02,
                                                    "",
                                                    "E (GeV)",
                                                    "#lambda_{0}^{2}",
                                                    0,0,0,
                                                    1,0.1,2.95,
                                                    1,4.95,50.05,0.8,0.65);
                        fHistClusterTrueEtaEM02->GetXaxis()->SetMoreLogLabels();
                        fHistClusterTrueEtaEM02->GetXaxis()->SetLabelOffset(-0.02);
                        SetZMinMaxTH2(fHistClusterTrueEtaEM02,fHistClusterTrueEtaEM02->GetXaxis()->FindBin(4.95),fHistClusterTrueEtaEM02->GetXaxis()->FindBin(50.05),
                                    fHistClusterTrueEtaEM02->GetYaxis()->FindBin(0.1),fHistClusterTrueEtaEM02->GetYaxis()->FindBin(2.95));
                        //fHistClusterTrueEtaEM02->GetZaxis()->SetRangeUser(1E-8,5E-5);
                        fHistClusterTrueEtaEM02->GetZaxis()->SetLabelSize(0.051);
                        fHistClusterTrueEtaEM02->GetXaxis()->SetTickLength(0.05);
                        fHistClusterTrueEtaEM02->DrawCopy("COLZ");
                        min->Draw("same");
                        cnst->Draw("same");
                        max->Draw("same");
                        WriteHistogram(fHistClusterTrueEtaEM02);
                        PutProcessLabelAndEnergyOnPlot(0.55, 0.99, 0.06,  "ALICE simulation", fCollisionSystem.Data(), "", 42, 0.03, "");
                        PutProcessLabelAndEnergyOnPlot(0.3, 0.6, 0.06, "#eta", "", "", 42, 0.03, "");
                        cvsEM02->SetLogx(1); cvsEM02->SetLogy(0); cvsEM02->SetLogz(1); cvsEM02->Update();
                        cvsEM02->SaveAs(Form("%s/EVsM02_TrueEta_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                        cvsEM02->Clear();
                    } else cout << "INFO: Object |TrueClusEtaEM02| could not be found! Skipping Draw..." << endl;
                    //---------------------------------------------------------------------------------------------------------------
                    // validated Electron plot
                    TH2D* fHistClusterTrueElectronEM02   = NULL;
                    if (isMerged)
                        fHistClusterTrueElectronEM02     = (TH2D*)TrueContainer->FindObject("TrueClusElectronEM02");
                    else
                        fHistClusterTrueElectronEM02     = (TH2D*)ClusterContainer->FindObject("TrueClusElectronEM02");

                    if(fHistClusterTrueElectronEM02){
                        fHistClusterTrueElectronEM02->Scale(1/nEvents[i]);
                        DrawAutoGammaHistoPaper2D(fHistClusterTrueElectronEM02,
                                                    "",
                                                    "E (GeV)",
                                                    "#lambda_{0}^{2}",
                                                    0,0,0,
                                                    1,0.1,2.95,
                                                    1,4.95,50.05,0.8,0.65);
                        fHistClusterTrueElectronEM02->GetXaxis()->SetMoreLogLabels();
                        fHistClusterTrueElectronEM02->GetXaxis()->SetLabelOffset(-0.02);
                        SetZMinMaxTH2(fHistClusterTrueElectronEM02,fHistClusterTrueElectronEM02->GetXaxis()->FindBin(4.95),fHistClusterTrueElectronEM02->GetXaxis()->FindBin(50.05),
                                    fHistClusterTrueElectronEM02->GetYaxis()->FindBin(0.1),fHistClusterTrueElectronEM02->GetYaxis()->FindBin(2.95));
                        //fHistClusterTrueElectronEM02->GetZaxis()->SetRangeUser(1E-8,5E-5);
                        fHistClusterTrueElectronEM02->GetZaxis()->SetLabelSize(0.051);
                        fHistClusterTrueElectronEM02->GetXaxis()->SetTickLength(0.05);
                        fHistClusterTrueElectronEM02->DrawCopy("COLZ");
                        min->Draw("same");
                        cnst->Draw("same");
                        max->Draw("same");
                        WriteHistogram(fHistClusterTrueElectronEM02);
                        PutProcessLabelAndEnergyOnPlot(0.55, 0.99, 0.06,  "ALICE simulation", fCollisionSystem.Data(), "", 42, 0.03, "");
                        PutProcessLabelAndEnergyOnPlot(0.12, 0.23, 0.06, "e^{#pm}", "", "", 42, 0.03, "");
                        cvsEM02->SetLogx(1); cvsEM02->SetLogy(0); cvsEM02->SetLogz(1); cvsEM02->Update();
                        cvsEM02->SaveAs(Form("%s/EVsM02_TrueElectron_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                        cvsEM02->Clear();
                    } else cout << "INFO: Object |TrueClusElectronEM02| could not be found! Skipping Draw..." << endl;
                }
                delete cvsEM02;
                delete min;
                delete cnst;
                delete max;

            }
        }

        //---------------------------------------------------------------------------------------------------------------
        //--------------------------- Cluster property: number of local maxima in cluster -------------------------------
        //---------------------------------------------------------------------------------------------------------------
        // number of local maxima in cluster before QA cuts
        TH1D* fHistNLMBeforeQA = (TH1D*)CaloCutsContainer->FindObject(Form("NLM_beforeClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistNLMBeforeQA){
            GetMinMaxBin(fHistNLMBeforeQA,minB,maxB);
            SetXRange(fHistNLMBeforeQA,minB-1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistNLMBeforeQA,"","Number of Local Maxima","#frac{1}{N} #frac{dN}{dNLM}",0.9,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistNLMBeforeQA, Form("%s/NLM_beforeQA_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << Form("INFO: Object |NLM_beforeClusterQA %s| could not be found! Skipping Draw...",fClusterCutSelection[i].Data()) << endl;
        //---------------------------------------------------------------------------------------------------------------
        // number of local maxima in cluster after QA cuts
        TH1D* fHistNLMAfterQA = (TH1D*)CaloCutsContainer->FindObject(Form("NLM_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistNLMAfterQA){
            GetMinMaxBin(fHistNLMAfterQA,minB,maxB);
            SetXRange(fHistNLMAfterQA,minB-1,maxB+1);
            DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                fHistNLMAfterQA,"","Number of Local Maxima","#frac{1}{N} #frac{dN}{dNLM}",0.9,1,
                                0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistNLMAfterQA, Form("%s/NLM_afterQA_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            vecClusterNLM.push_back(new TH1D(*fHistNLMAfterQA));
        } else cout << Form("INFO: Object |NLM_afterClusterQA %s| could not be found! Skipping Draw...",fClusterCutSelection[i].Data()) << endl;
        //---------------------------------------------------------------------------------------------------------------
        // number of local maxima in cluster vs NCells in cluster
        TH2D* fHistNLMNCells = (TH2D*)CaloCutsContainer->FindObject(Form("NLM_NCells_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistNLMNCells){
            GetMinMaxBin(fHistNLMNCells,minB,maxB); minB-=1; maxB+=1;
            SetXRange(fHistNLMNCells,minB,maxB);
            Int_t minYB=0; Int_t maxYB=1;
            GetMinMaxBinY(fHistNLMNCells,minYB,maxYB); maxYB+=1;
            SetYRange(fHistNLMNCells,minYB,maxYB);
            SetZMinMaxTH2(fHistNLMNCells,minB,maxB,minYB,maxYB);
            DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                fHistNLMNCells,"",
                                "Number of Local Maxima","#it{N}_{Cells} per Cluster",1,1,
                                0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
            SaveCanvasAndWriteHistogram(canvas, fHistNLMNCells, Form("%s/NLM_NCells_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

        } else cout << Form("INFO: Object |NLM_NCells_afterClusterQA %s| could not be found! Skipping Draw...",fClusterCutSelection[i].Data()) << endl;
        //---------------------------------------------------------------------------------------------------------------
        // number of local maxima in cluster vs cluster energy
        TH2D* fHistNLME = (TH2D*)CaloCutsContainer->FindObject(Form("NLM_E_afterClusterQA %s", fClusterCutSelection[i].Data()));
        if(fHistNLME){
        GetMinMaxBin(fHistNLME,minB,maxB); minB-=1; maxB+=1;
        SetXRange(fHistNLME,minB,maxB);
        Int_t minYB=0; Int_t maxYB=1;
        GetMinMaxBinY(fHistNLME,minYB,maxYB);
        SetYRange(fHistNLME,minYB,maxYB);
        SetZMinMaxTH2(fHistNLME,minB,maxB,minYB,maxYB);
        DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                            fHistNLME,"",
                            "Number of Local Maxima","Cluster Energy (GeV)",1,1,
                            0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
        SaveCanvasAndWriteHistogram(canvas, fHistNLME, Form("%s/NLM_E_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
        } else cout << Form("INFO: Object |NLM_E_afterClusterQA %s| could not be found! Skipping Draw...",fClusterCutSelection[i].Data()) << endl;




        //---------------------------------------------------------------------------------------------------------------
        //---------------------------------- Cluster track matching detailed distributions ------------------------------
        //---------------------------------------------------------------------------------------------------------------
        if(isTrackMatching){
            TBox *cutPos    = new TBox(-0.016,-0.09,0.016,0.06);
            cutPos->SetFillStyle(0);
            cutPos->SetFillColor(0);
            cutPos->SetLineColor(kRed+3);
            cutPos->SetLineStyle(3);
            cutPos->SetLineWidth(4);
            TBox *cutNeg    = new TBox(-0.016,-0.06,0.016,0.09);
            cutNeg->SetFillStyle(0);
            cutNeg->SetFillColor(0);
            cutNeg->SetLineColor(kRed+3);
            cutNeg->SetLineStyle(5);
            cutNeg->SetLineWidth(4);

            //---------------------------------------------------------------------------------------------------------------
            //cluster track matching in dEta, dPhi before cluster QA cuts
            TH2D* fHistClusterdEtadPhiBeforeQA = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_beforeClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterdEtadPhiBeforeQA){
                SetZMinMaxTH2(fHistClusterdEtadPhiBeforeQA,1,fHistClusterdEtadPhiBeforeQA->GetNbinsX(),1,fHistClusterdEtadPhiBeforeQA->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiBeforeQA,"",
                                    "#Delta#eta_{cluster - all tracks}","#Delta#phi_{cluster - all tracks}",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiBeforeQA, Form("%s/dEtaVsdPhi_allTracks_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiBeforeQA,"",
                                    "#Delta#eta_{cluster - all tracks}","#Delta#phi_{cluster - all tracks}",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                if (trackMatch == 3){
                    cutPos->Draw();
                    cutNeg->Draw();
                }
                SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_allTracks_WithCuts_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                // projection of dEta, dPhi before cluster QA cuts to dEta axis
                TH1D* fHistClusterdEtaBeforeQA = (TH1D*) fHistClusterdEtadPhiBeforeQA->ProjectionX("fHistClusterdEtaBeforeQA",1,fHistClusterdEtadPhiBeforeQA->GetNbinsY());
                SetXRange(fHistClusterdEtaBeforeQA,1,300);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterdEtaBeforeQA,"","#Delta#eta_{cluster - all tracks}","#frac{d#it{N}}{d#Delta#eta}",0.9,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and charged tracks", "");
                SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtaBeforeQA, Form("%s/dEta_allTracks_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecMatchingDeltaEtaPhi_allTracks[0].push_back(fHistClusterdEtaBeforeQA);

                // projection of dEta, dPhi before cluster QA cuts to dPhi axis
                TH1D* fHistClusterdPhiBeforeQA = (TH1D*) fHistClusterdEtadPhiBeforeQA->ProjectionY("fHistClusterdPhiBeforeQA",1,fHistClusterdEtadPhiBeforeQA->GetNbinsX());
                SetXRange(fHistClusterdPhiBeforeQA,1,300);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterdPhiBeforeQA,"","#Delta#phi_{cluster - all tracks}","#frac{d#it{N}}{d#Delta#phi}",0.9,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and charged tracks", "");
                SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhiBeforeQA, Form("%s/dPhi_allTracks_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecMatchingDeltaEtaPhi_allTracks[1].push_back(fHistClusterdPhiBeforeQA);

                // |axis|<0.1 - Zoom in for notes
                SetXRange(fHistClusterdEtadPhiBeforeQA,fHistClusterdEtadPhiBeforeQA->GetXaxis()->FindBin(-0.1),fHistClusterdEtadPhiBeforeQA->GetXaxis()->FindBin(0.1));
                SetYRange(fHistClusterdEtadPhiBeforeQA,fHistClusterdEtadPhiBeforeQA->GetYaxis()->FindBin(-0.1),fHistClusterdEtadPhiBeforeQA->GetYaxis()->FindBin(0.1));
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiBeforeQA,"",
                                    "#Delta#eta","#Delta#phi",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and charged tracks", "");
                SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_allTracks_zoom_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiBeforeQA,"",
                                    "#Delta#eta","#Delta#phi",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                if (trackMatch == 3){
                    cutPos->Draw();
                    cutNeg->Draw();
                }
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and charged tracks", "");
                SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_allTracks_WithCuts_zoom_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |fHistClusterdEtadPhiBeforeQA| could not be found! Skipping Draw..." << endl;
            delete cutPos;
            delete cutNeg;

            //---------------------------------------------------------------------------------------------------------------
            //cluster track matching in dEta, dPhi after cluster QA cuts
            TH2D* fHistClusterdEtadPhiAfterQA = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterdEtadPhiAfterQA){
                SetZMinMaxTH2(fHistClusterdEtadPhiAfterQA,1,fHistClusterdEtadPhiAfterQA->GetNbinsX(),1,fHistClusterdEtadPhiAfterQA->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiAfterQA,"",
                                    "#Delta#eta","#Delta#phi",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and unmatched tracks", "");
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiAfterQA, Form("%s/dEtaVsdPhi_unmatched_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                // projection of dEta, dPhi after cluster QA cuts to dEta axis
                TH1D* fHistClusterdEtaAfterQA = (TH1D*) fHistClusterdEtadPhiAfterQA->ProjectionX("fHistClusterdEtaAfterQA",1,fHistClusterdEtadPhiAfterQA->GetNbinsY());
                SetXRange(fHistClusterdEtaAfterQA,1,300);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterdEtaAfterQA,"","#Delta#eta_{cluster - unmatched tracks}","#frac{d#it{N}}{d#Delta#eta}",0.9,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtaAfterQA, Form("%s/dEta_unmatched_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                delete fHistClusterdEtaAfterQA;

                // projection of dEta, dPhi after cluster QA cuts to dPhi axis
                TH1D* fHistClusterdPhiAfterQA = (TH1D*) fHistClusterdEtadPhiAfterQA->ProjectionY("fHistClusterdPhiAfterQA",1,fHistClusterdEtadPhiAfterQA->GetNbinsX());
                SetXRange(fHistClusterdPhiAfterQA,1,300);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterdPhiAfterQA,"","#Delta#phi_{cluster - unmatched tracks}","#frac{d#it{N}}{d#Delta#phi}",0.9,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhiAfterQA, Form("%s/dPhi_unmatched_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                delete fHistClusterdPhiAfterQA;

                // |axis|<0.1
                SetXRange(fHistClusterdEtadPhiAfterQA,fHistClusterdEtadPhiAfterQA->GetXaxis()->FindBin(-0.1),fHistClusterdEtadPhiAfterQA->GetXaxis()->FindBin(0.1));
                SetYRange(fHistClusterdEtadPhiAfterQA,fHistClusterdEtadPhiAfterQA->GetYaxis()->FindBin(-0.1),fHistClusterdEtadPhiAfterQA->GetYaxis()->FindBin(0.1));
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiAfterQA,"",
                                    "#Delta#eta","#Delta#phi",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and unmatched tracks", "");
                SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_unmatched_zoom_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |fHistClusterdEtadPhiAfterQA| could not be found! Skipping Draw..." << endl;

            //---------------------------------------------------------------------------------------------------------------
            // check which clusters have been thrown out by the track matching cuts in dPhi & dEta
            if(fHistClusterdEtadPhiBeforeQA && fHistClusterdEtadPhiAfterQA){
                TH2D* fHistClusterdEtadPhiMatchedAfterQA = (TH2D*) fHistClusterdEtadPhiBeforeQA->Clone();
                fHistClusterdEtadPhiMatchedAfterQA->Add(fHistClusterdEtadPhiAfterQA,-1);
                SetZMinMaxTH2(fHistClusterdEtadPhiMatchedAfterQA,1,fHistClusterdEtadPhiMatchedAfterQA->GetNbinsX(),1,fHistClusterdEtadPhiMatchedAfterQA->GetNbinsY());
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiMatchedAfterQA,"",
                                    "#Delta#eta","#Delta#phi",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and matched tracks", "");
                SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiMatchedAfterQA, Form("%s/dEtaVsdPhi_matched_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                // check which clusters have been thrown out by the track matching cuts in dEta
                TGaxis::SetExponentOffset(-0.04, -0.04, "x");
                TH1D* fHistClusterdEtaMatchedAfterQA = (TH1D*) fHistClusterdEtadPhiMatchedAfterQA->ProjectionX("fHistClusterdEtaMatchedAfterQA",1,fHistClusterdEtadPhiMatchedAfterQA->GetNbinsY());
                GetMinMaxBin(fHistClusterdEtaMatchedAfterQA,minB,maxB);
                SetXRange(fHistClusterdEtaMatchedAfterQA,minB-1,maxB+1);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterdEtaMatchedAfterQA,"","#Delta#eta_{cluster - matched tracks}","#frac{d#it{N}}{d#Delta#eta}",0.9,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtaMatchedAfterQA, Form("%s/dEta_matched_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecMatchingDeltaEtaPhi_matched[0].push_back(fHistClusterdEtaMatchedAfterQA);

                // check which clusters have been thrown out by the track matching cuts in dPhi
                TH1D* fHistClusterdPhiMatchedAfterQA = (TH1D*) fHistClusterdEtadPhiMatchedAfterQA->ProjectionY("fHistClusterdPhiMatchedAfterQA",1,fHistClusterdEtadPhiMatchedAfterQA->GetNbinsX());
                GetMinMaxBin(fHistClusterdPhiMatchedAfterQA,minB,maxB);
                SetXRange(fHistClusterdPhiMatchedAfterQA,minB-1,maxB+1);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterdPhiMatchedAfterQA,"","#Delta#phi_{cluster - matched tracks}","#frac{d#it{N}}{d#Delta#phi}",0.9,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhiMatchedAfterQA, Form("%s/dPhi_matched_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                TGaxis::SetExponentOffset(0, 0, "x");
                vecMatchingDeltaEtaPhi_matched[1].push_back(fHistClusterdPhiMatchedAfterQA);

                // |axis|<0.1 reduced range for note
                SetXRange(fHistClusterdEtadPhiMatchedAfterQA,fHistClusterdEtadPhiMatchedAfterQA->GetXaxis()->FindBin(-0.1),fHistClusterdEtadPhiMatchedAfterQA->GetXaxis()->FindBin(0.1));
                SetYRange(fHistClusterdEtadPhiMatchedAfterQA,fHistClusterdEtadPhiMatchedAfterQA->GetYaxis()->FindBin(-0.1),fHistClusterdEtadPhiMatchedAfterQA->GetYaxis()->FindBin(0.1));
                DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                    fHistClusterdEtadPhiMatchedAfterQA,"",
                                    "#Delta#eta","#Delta#phi",1,1.2,
                                    0.85,0.9,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                PutProcessLabelAndEnergyOnPlot(0.15,0.9,0.03, fClusters.Data(), "and matched tracks", "");
                SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_matched_zoom_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                delete fHistClusterdEtadPhiMatchedAfterQA;
            }

            //---------------------------------------------------------------------------------------------------------------
            // distance of track to cluster before QA cuts
            TH1D* fHistClusterRBeforeQA = (TH1D*)CaloCutsContainer->FindObject(Form("R_Cluster_beforeClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterRBeforeQA){
                GetMinMaxBin(fHistClusterRBeforeQA,minB,maxB);
                SetXRange(fHistClusterRBeforeQA,minB-5,maxB+5);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterRBeforeQA,"","#it{R}_{Cluster} (cm)","#frac{d#it{N}}{d#it{R}}",1,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterRBeforeQA, Form("%s/R_Cluster_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |fHistClusterRBeforeQA| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // distance of track to cluster after QA cuts
            TH1D* fHistClusterRAfterQA = (TH1D*)CaloCutsContainer->FindObject(Form("R_Cluster_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterRAfterQA){
                GetMinMaxBin(fHistClusterRAfterQA,minB,maxB);
                SetXRange(fHistClusterRAfterQA,minB-5,maxB+5);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterRAfterQA,"","#it{R}_{Cluster} (cm)","#frac{d#it{N}}{d#it{R}}",1,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterRAfterQA, Form("%s/R_Cluster_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecClusterR.push_back(new TH1D(*fHistClusterRAfterQA));
            } else cout << "INFO: Object |fHistClusterRAfterQA| could not be found! Skipping Draw..." << endl;


            //---------------------------------------------------------------------------------------------------------------
            // distance of track to cluster in units of dphi and dEta before QA cuts
            TH1D* fHistDistanceTrackToClusterBeforeQA = (TH1D*)CaloCutsContainer->FindObject(Form("DistanceToTrack_beforeClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistDistanceTrackToClusterBeforeQA){
                AdjustHistRange(fHistDistanceTrackToClusterBeforeQA,5,5,kTRUE,1,1.);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,0.15,kFALSE,kTRUE,kFALSE,
                                    fHistDistanceTrackToClusterBeforeQA,"","r = #sqrt{#Delta#eta_{cluster - all tracks}^{2}+#Delta#phi_{cluster - all tracks}^{2}}","#frac{d#it{N}}{d#it{r}}",1.3,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistDistanceTrackToClusterBeforeQA, Form("%s/DistanceToTrack_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |fHistDistanceTrackToClusterBeforeQA| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // distance of track to cluster in units of dphi and dEta after QA cuts
            TH1D* fHistDistanceTrackToClusterAfterQA = (TH1D*)CaloCutsContainer->FindObject(Form("DistanceToTrack_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistDistanceTrackToClusterAfterQA){
                AdjustHistRange(fHistDistanceTrackToClusterAfterQA,5,5,kTRUE,1,1.);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,0.15,kFALSE,kTRUE,kFALSE,
                                    fHistDistanceTrackToClusterAfterQA,"","r = #sqrt{#Delta#eta_{cluster - unmatched tracks}^{2}+#Delta#phi_{cluster - unmatched tracks}^{2}}","#frac{d#it{N}}{d#it{r}}",1.3,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistDistanceTrackToClusterAfterQA, Form("%s/DistanceToTrack_unmatched_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |fHistDistanceTrackToClusterAfterQA| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // distance of matched cluster to the tracks
            if(fHistDistanceTrackToClusterBeforeQA && fHistDistanceTrackToClusterAfterQA){
                TH1D* fHistDistanceTrackToClusterMatchedAfterQA = (TH1D*) fHistDistanceTrackToClusterBeforeQA->Clone();
                fHistDistanceTrackToClusterMatchedAfterQA->Add(fHistDistanceTrackToClusterAfterQA, -1);
                AdjustHistRange(fHistDistanceTrackToClusterMatchedAfterQA,5,5,kTRUE,1,1.);
                GetMinMaxBin(fHistDistanceTrackToClusterMatchedAfterQA,minB,maxB);
                SetXRange(fHistDistanceTrackToClusterMatchedAfterQA,minB,maxB+5);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,0.15,kFALSE,kTRUE,kFALSE,
                                    fHistDistanceTrackToClusterMatchedAfterQA,"","r = #sqrt{#Delta#eta_{cluster - matched tracks}^{2}+#Delta#phi_{cluster - matched tracks}^{2}}","#frac{d#it{N}}{d#it{r}}",1.3,1,
                                    0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                SaveCanvasAndWriteHistogram(canvas, fHistDistanceTrackToClusterMatchedAfterQA, Form("%s/DistanceToTrack_matched_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                delete fHistDistanceTrackToClusterMatchedAfterQA;
            }

            //---------------------------------------------------------------------------------------------------------------
            //--------------------------------- Pt dependent matching properties --------------------------------------------
            //---------------------------------------------------------------------------------------------------------------
            if(doExtQA>0){
                //---------------------------------------------------------------------------------------------------------------
                // Pt dependent matching variable dEta
                TH2D* fHistClusterdEtadPtBeforeQA   = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsPt_beforeClusterQA %s", fClusterCutSelection[i].Data()));
                if(fHistClusterdEtadPtBeforeQA){
                    TH1D* fHistClusterdEtadPtBeforeQAPt = (TH1D*) fHistClusterdEtadPtBeforeQA->ProjectionY("dEtaVsPtOnPt_beforeClusterQA",1,fHistClusterdEtadPtBeforeQA->GetNbinsX());
                    GetMinMaxBin(fHistClusterdEtadPtBeforeQAPt,minB,maxB);
                    delete fHistClusterdEtadPtBeforeQAPt;
                    SetYRange(fHistClusterdEtadPtBeforeQA,minB,maxB+5);
                    SetZMinMaxTH2(fHistClusterdEtadPtBeforeQA,1,fHistClusterdEtadPtBeforeQA->GetNbinsX(),minB,maxB+5);
                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                        fHistClusterdEtadPtBeforeQA,"",
                                        "#Delta#eta_{cluster - charged tracks}","#it{p}_{T} (GeV/#it{c})",1,1,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtadPtBeforeQA, Form("%s/dEtaVsPt_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterdEtadPtBeforeQA| could not be found! Skipping Draw..." << endl;
                //---------------------------------------------------------------------------------------------------------------
                // Pt dependent matching variable dPhi
                TH2D* fHistClusterdPhidPt_pos_BeforeQA   = (TH2D*)CaloCutsContainer->FindObject(Form("dPhiVsPt_posTracks_beforeClusterQA %s", fClusterCutSelection[i].Data()));
                TH2D* fHistClusterdPhidPt_neg_BeforeQA   = (TH2D*)CaloCutsContainer->FindObject(Form("dPhiVsPt_negTracks_beforeClusterQA %s", fClusterCutSelection[i].Data()));
                TH2D* fHistClusterdPhidPtBeforeQA = 0x0;
                if(fHistClusterdPhidPt_neg_BeforeQA && fHistClusterdPhidPt_pos_BeforeQA){
                  fHistClusterdPhidPtBeforeQA = fHistClusterdPhidPt_pos_BeforeQA;
                  fHistClusterdPhidPtBeforeQA->Add(fHistClusterdPhidPt_neg_BeforeQA);
                }
                if(fHistClusterdPhidPtBeforeQA){
                    TH1D* fHistClusterdPhidPtBeforeQAPt = (TH1D*) fHistClusterdPhidPtBeforeQA->ProjectionY("dPhiVsPtOnPt_beforeClusterQA",1,fHistClusterdPhidPtBeforeQA->GetNbinsX());
                    GetMinMaxBin(fHistClusterdPhidPtBeforeQAPt,minB,maxB);
                    delete fHistClusterdPhidPtBeforeQAPt;
                    SetYRange(fHistClusterdPhidPtBeforeQA,minB,maxB+5);
                    SetZMinMaxTH2(fHistClusterdPhidPtBeforeQA,1,fHistClusterdPhidPtBeforeQA->GetNbinsX(),minB,maxB+5);
                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                        fHistClusterdPhidPtBeforeQA,"",
                                        "#Delta#phi_{cluster - charged tracks}","#it{p}_{T} (GeV/#it{c})",1,1,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhidPtBeforeQA, Form("%s/dPhiVsPt_%s_beforeClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterdPhidPtBeforeQA| could not be found! Skipping Draw..." << endl;

                // Pt dependent matching variable dEta
                TH2D* fHistClusterdEtadPtAfterQA   = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsPt_afterClusterQA %s", fClusterCutSelection[i].Data()));
                if(fHistClusterdEtadPtAfterQA){
                    TH1D* fHistClusterdEtadPtAfterQAPt = (TH1D*) fHistClusterdEtadPtAfterQA->ProjectionY("dEtaVsPtOnPt_afterClusterQA",1,fHistClusterdEtadPtAfterQA->GetNbinsX());
                    GetMinMaxBin(fHistClusterdEtadPtAfterQAPt,minB,maxB);
                    delete fHistClusterdEtadPtAfterQAPt;
                    SetYRange(fHistClusterdEtadPtAfterQA,minB,maxB+5);
                    SetZMinMaxTH2(fHistClusterdEtadPtAfterQA,1,fHistClusterdEtadPtAfterQA->GetNbinsX(),minB,maxB+5);
                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                        fHistClusterdEtadPtAfterQA,"",
                                        "#Delta#eta_{cluster - charged tracks}","#it{p}_{T} (GeV/#it{c})",1,1,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtadPtAfterQA, Form("%s/dEtaVsPt_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterdEtadPtAfterQA| could not be found! Skipping Draw..." << endl;
                //---------------------------------------------------------------------------------------------------------------
                // Pt dependent matching variable dPhi
                TH2D* fHistClusterdPhidPtAfterQA   = (TH2D*)CaloCutsContainer->FindObject(Form("dPhiVsPt_afterClusterQA %s", fClusterCutSelection[i].Data()));
                if(fHistClusterdPhidPtAfterQA){
                    TH1D* fHistClusterdPhidPtAfterQAPt = (TH1D*) fHistClusterdPhidPtAfterQA->ProjectionY("dPhiVsPtOnPt_afterClusterQA",1,fHistClusterdPhidPtAfterQA->GetNbinsX());
                    GetMinMaxBin(fHistClusterdPhidPtAfterQAPt,minB,maxB);
                    delete fHistClusterdPhidPtAfterQAPt;
                    SetYRange(fHistClusterdPhidPtAfterQA,minB,maxB+5);
                    SetZMinMaxTH2(fHistClusterdPhidPtAfterQA,1,fHistClusterdPhidPtAfterQA->GetNbinsX(),minB,maxB+5);
                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                        fHistClusterdPhidPtAfterQA,"",
                                        "#Delta#phi_{cluster - charged tracks}","#it{p}_{T} (GeV/#it{c})",1,1,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhidPtAfterQA, Form("%s/dPhiVsPt_%s_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterdPhidPtAfterQA| could not be found! Skipping Draw..." << endl;


                // Pt dependent matching variable dEta
                if(fHistClusterdEtadPtBeforeQA && fHistClusterdEtadPtAfterQA){
                    fHistClusterdEtadPtBeforeQA->Add(fHistClusterdEtadPtAfterQA,-1.);
                    TH1D* fHistClusterdEtadPtAfterQAPt = (TH1D*) fHistClusterdEtadPtBeforeQA->ProjectionY("dEtaVsPtOnPt_unmatched_afterClusterQA",1,fHistClusterdEtadPtAfterQA->GetNbinsX());
                    GetMinMaxBin(fHistClusterdEtadPtAfterQAPt,minB,maxB);
                    delete fHistClusterdEtadPtAfterQAPt;
                    SetYRange(fHistClusterdEtadPtBeforeQA,minB,maxB+5);
                    SetZMinMaxTH2(fHistClusterdEtadPtBeforeQA,1,fHistClusterdEtadPtBeforeQA->GetNbinsX(),minB,maxB+5);
                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                        fHistClusterdEtadPtBeforeQA,"",
                                        "#Delta#eta_{cluster - charged tracks (not matched)}","#it{p}_{T} (GeV/#it{c})",1,1,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtadPtBeforeQA, Form("%s/dEtaVsPt_%s__unmatched_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterdEtadPtAfterQA| could not be found! Skipping Draw..." << endl;
                //---------------------------------------------------------------------------------------------------------------
                // Pt dependent matching variable dPhi
                if(fHistClusterdPhidPtBeforeQA && fHistClusterdPhidPtAfterQA){
                    fHistClusterdPhidPtBeforeQA->Add(fHistClusterdPhidPtAfterQA,-1.);
                    TH1D* fHistClusterdPhidPtAfterQAPt = (TH1D*) fHistClusterdPhidPtBeforeQA->ProjectionY("dPhiVsPtOnPt_unmatched_afterClusterQA",1,fHistClusterdPhidPtAfterQA->GetNbinsX());
                    GetMinMaxBin(fHistClusterdPhidPtAfterQAPt,minB,maxB);
                    delete fHistClusterdPhidPtAfterQAPt;
                    SetYRange(fHistClusterdPhidPtBeforeQA,minB,maxB+5);
                    SetZMinMaxTH2(fHistClusterdPhidPtBeforeQA,1,fHistClusterdPhidPtBeforeQA->GetNbinsX(),minB,maxB+5);
                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                        fHistClusterdPhidPtBeforeQA,"",
                                        "#Delta#phi_{cluster - charged tracks (not matched)}","#it{p}_{T} (GeV/#it{c})",1,1,
                                        0.85,0.25,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhidPtBeforeQA, Form("%s/dPhiVsPt_%s_unmatched_afterClusterQA.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                } else cout << "INFO: Object |fHistClusterdPhidPtAfterQA| could not be found! Skipping Draw..." << endl;

                // separating the matching parameters according to their charges
                for(Int_t iCharge=0;iCharge<2;iCharge++){
                    //---------------------------------------------------------------------------------------------------------------
                    // dPhi dEta before QA cuts for a single charge
                    TH2D* fHistClusterdEtadPhiTracksBeforeQA = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_%sTracks_beforeClusterQA %s", charge[iCharge].Data(),fClusterCutSelection[i].Data()));
                    if(fHistClusterdEtadPhiTracksBeforeQA){
                        SetZMinMaxTH2(fHistClusterdEtadPhiTracksBeforeQA,1,fHistClusterdEtadPhiTracksBeforeQA->GetNbinsX(),1,fHistClusterdEtadPhiTracksBeforeQA->GetNbinsY());
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksBeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "");
                        SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiTracksBeforeQA, Form("%s/dEtaVsdPhi_all_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(),
                                                    DataSets[i].Data(), suffix.Data()));

                        TH1D* fHistClusterdEtaBeforeQA = (TH1D*) fHistClusterdEtadPhiTracksBeforeQA->ProjectionX(Form("fHistClusterdEta%sBeforeQA", charge[iCharge].Data()),1,
                                                                                                                 fHistClusterdEtadPhiTracksBeforeQA->GetNbinsY());
                        SetXRange(fHistClusterdEtaBeforeQA,1,300);
                        DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                            fHistClusterdEtaBeforeQA,"",Form("#Delta#eta_{Cluster - %s. Tracks}", charge[iCharge].Data()),"#frac{d#it{N}}{d#Delta#eta}",0.9,1,
                                            0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtaBeforeQA, Form("%s/dEta_all_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(),
                                                    suffix.Data()));
                        vecMatchingDeltaEtaPhiTracksAll[0][iCharge].push_back(fHistClusterdEtaBeforeQA);

                        TH1D* fHistClusterdPhiBeforeQA = (TH1D*) fHistClusterdEtadPhiTracksBeforeQA->ProjectionY(Form("fHistClusterdPhi%sBeforeQA", charge[iCharge].Data()),1,
                                                                                                                 fHistClusterdEtadPhiTracksBeforeQA->GetNbinsX());
                        SetXRange(fHistClusterdPhiBeforeQA,1,300);
                        DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                            fHistClusterdPhiBeforeQA,"",Form("#Delta#phi_{Cluster - %s. Tracks}", charge[iCharge].Data()),"#frac{d#it{N}}{d#Delta#phi}",0.9,1,
                                            0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhiBeforeQA, Form("%s/dPhi_all_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                        vecMatchingDeltaEtaPhiTracksAll[1][iCharge].push_back(fHistClusterdPhiBeforeQA);

                        // |axis|<0.1
                        SetXRange(fHistClusterdEtadPhiTracksBeforeQA,fHistClusterdEtadPhiTracksBeforeQA->GetXaxis()->FindBin(-0.1),fHistClusterdEtadPhiTracksBeforeQA->GetXaxis()->FindBin(0.1));
                        SetYRange(fHistClusterdEtadPhiTracksBeforeQA,fHistClusterdEtadPhiTracksBeforeQA->GetYaxis()->FindBin(-0.1),fHistClusterdEtadPhiTracksBeforeQA->GetYaxis()->FindBin(0.1));
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksBeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "");
                        SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_all_%sTracks_zoom_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                    } else cout << "INFO: Object |fHistClusterdEtadPhiTracksBeforeQA| could not be found! Skipping Draw..." << endl;
                    //---------------------------------------------------------------------------------------------------------------
                    // dPhi dEta after QA cuts for a single charge
                    TH2D* fHistClusterdEtadPhiTracksAfterQA = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_%sTracks_afterClusterQA %s", charge[iCharge].Data(),fClusterCutSelection[i].Data()));
                    if(fHistClusterdEtadPhiTracksAfterQA){
                        SetZMinMaxTH2(fHistClusterdEtadPhiTracksAfterQA,1,fHistClusterdEtadPhiTracksAfterQA->GetNbinsX(),1,fHistClusterdEtadPhiTracksAfterQA->GetNbinsY());
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksAfterQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and unmatched %s. tracks",charge[iCharge].Data()), "");
                        SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiTracksAfterQA, Form("%s/dEtaVsdPhi_unmatched_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(),
                                                    DataSets[i].Data(), suffix.Data()));

                        TH1D* fHistClusterdEtaAfterQA = (TH1D*) fHistClusterdEtadPhiTracksAfterQA->ProjectionX(Form("fHistClusterdEta%sAfterQA", charge[iCharge].Data()),1,
                                                                                                               fHistClusterdEtadPhiTracksAfterQA->GetNbinsY());
                        SetXRange(fHistClusterdEtaAfterQA,1,300);
                        DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                            fHistClusterdEtaAfterQA,"",Form("#Delta#eta_{Cluster - unmatched %s. Tracks}", charge[iCharge].Data()),"#frac{d#it{N}}{d#Delta#eta}",0.9,1,
                                            0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtaAfterQA, Form("%s/dEta_unmatched_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                        delete fHistClusterdEtaAfterQA;

                        TH1D* fHistClusterdPhiAfterQA = (TH1D*) fHistClusterdEtadPhiTracksAfterQA->ProjectionY(Form("fHistClusterdPhi%sAfterQA", charge[iCharge].Data()),1,
                                                                                                               fHistClusterdEtadPhiTracksAfterQA->GetNbinsX());
                        SetXRange(fHistClusterdPhiAfterQA,1,300);
                        DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                            fHistClusterdPhiAfterQA,"",Form("#Delta#phi_{Cluster - unmatched %s. Tracks}", charge[iCharge].Data()),"#frac{d#it{N}}{d#Delta#phi}",0.9,1,
                                            0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhiAfterQA, Form("%s/dPhi_unmatched_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                        delete fHistClusterdPhiAfterQA;

                        // |axis|<0.1
                        SetXRange(fHistClusterdEtadPhiTracksAfterQA,fHistClusterdEtadPhiTracksAfterQA->GetXaxis()->FindBin(-0.1),fHistClusterdEtadPhiTracksAfterQA->GetXaxis()->FindBin(0.1));
                        SetYRange(fHistClusterdEtadPhiTracksAfterQA,fHistClusterdEtadPhiTracksAfterQA->GetYaxis()->FindBin(-0.1),fHistClusterdEtadPhiTracksAfterQA->GetYaxis()->FindBin(0.1));
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksAfterQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and unmatched %s. tracks",charge[iCharge].Data()), "");
                        SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_unmatched_%sTracks_zoom_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));

                    } else cout << "INFO: Object |fHistClusterdEtadPhiTracksAfterQA| could not be found! Skipping Draw..." << endl;
                    //---------------------------------------------------------------------------------------------------------------
                    // dPhi dEta for a single charges matched to a cluster
                    if(fHistClusterdEtadPhiTracksBeforeQA && fHistClusterdEtadPhiTracksAfterQA){
                        TH2D* fHistClusterdEtadPhiTracksMatchedAfterQA = (TH2D*) fHistClusterdEtadPhiTracksBeforeQA->Clone();
                        fHistClusterdEtadPhiTracksMatchedAfterQA->Add(fHistClusterdEtadPhiTracksAfterQA,-1);
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksMatchedAfterQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and matched %s. tracks",charge[iCharge].Data()), "");
                        SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiTracksMatchedAfterQA, Form("%s/dEtaVsdPhi_matched_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(),
                                                    DataSets[i].Data(), suffix.Data()));

                        TH1D* fHistClusterdEtaMatchedAfterQA = (TH1D*) fHistClusterdEtadPhiTracksMatchedAfterQA->ProjectionX(Form("fHistClusterdEta%sAfterQA", charge[iCharge].Data()),
                                                                                                                             1, fHistClusterdEtadPhiTracksMatchedAfterQA->GetNbinsY());
                        GetMinMaxBin(fHistClusterdEtaMatchedAfterQA,minB,maxB);
                        SetXRange(fHistClusterdEtaMatchedAfterQA,minB-1,maxB+1);
                        DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                            fHistClusterdEtaMatchedAfterQA,"",Form("#Delta#eta_{Cluster - matched %s. Tracks}", charge[iCharge].Data()),"#frac{d#it{N}}{d#Delta#eta}",0.9,1,
                                            0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        SaveCanvasAndWriteHistogram(canvas, fHistClusterdEtaMatchedAfterQA, Form("%s/dEta_matched_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                        vecMatchingDeltaEtaPhiTracksMatched[0][iCharge].push_back(fHistClusterdEtaMatchedAfterQA);

                        TH1D* fHistClusterdPhiMatchedAfterQA = (TH1D*) fHistClusterdEtadPhiTracksMatchedAfterQA->ProjectionY(Form("fHistClusterdPhi%sAfterQA", charge[iCharge].Data()), 1,
                                                                                                                             fHistClusterdEtadPhiTracksMatchedAfterQA->GetNbinsX());
                        GetMinMaxBin(fHistClusterdPhiMatchedAfterQA,minB,maxB);
                        SetXRange(fHistClusterdPhiMatchedAfterQA,minB-1,maxB+1);
                        DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                            fHistClusterdPhiMatchedAfterQA,"",Form("#Delta#phi_{Cluster - matched %s. Tracks}", charge[iCharge].Data()),"#frac{d#it{N}}{d#Delta#phi}",0.9,1,
                                            0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        SaveCanvasAndWriteHistogram(canvas, fHistClusterdPhiMatchedAfterQA, Form("%s/dPhi_matched_%sTracks_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                        vecMatchingDeltaEtaPhiTracksMatched[1][iCharge].push_back(fHistClusterdPhiMatchedAfterQA);

                        // |axis|<0.1
                        SetXRange(fHistClusterdEtadPhiTracksMatchedAfterQA,fHistClusterdEtadPhiTracksMatchedAfterQA->GetXaxis()->FindBin(-0.1),fHistClusterdEtadPhiTracksMatchedAfterQA->GetXaxis()->FindBin(0.1));
                        SetYRange(fHistClusterdEtadPhiTracksMatchedAfterQA,fHistClusterdEtadPhiTracksMatchedAfterQA->GetYaxis()->FindBin(-0.1),fHistClusterdEtadPhiTracksMatchedAfterQA->GetYaxis()->FindBin(0.1));
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksMatchedAfterQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and matched %s. tracks",charge[iCharge].Data()), "");
                        SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_matched_%sTracks_zoom_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));

                        delete fHistClusterdEtadPhiTracksMatchedAfterQA;
                    }
                    //---------------------------------------------------------------------------------------------------------------
                    // dPhi dEta before QA cuts for a single charge tracks at low track momentum (pT < 0.75)
                    TH2D* fHistClusterdEtadPhiTracksP_000_075BeforeQA = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_%sTracks_P<0.75_beforeClusterQA %s", charge[iCharge].Data(),
                                                                                                             fClusterCutSelection[i].Data()));
                    if(fHistClusterdEtadPhiTracksP_000_075BeforeQA){
                        SetZMinMaxTH2(fHistClusterdEtadPhiTracksP_000_075BeforeQA,1,fHistClusterdEtadPhiTracksP_000_075BeforeQA->GetNbinsX(),1,fHistClusterdEtadPhiTracksP_000_075BeforeQA->GetNbinsY());
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksP_000_075BeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "with #it{p}_{T}<0.75 GeV/#it{c}");
                        SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiTracksP_000_075BeforeQA, Form("%s/dEtaVsdPhi_%sTracks_000_075_%s.%s", outputDir.Data(), charge[iCharge].Data(),
                                                    DataSets[i].Data(), suffix.Data()));
                        vecMatchingDeltaEtaPhiTracksAllPt[0][iCharge].push_back(new TH2D(*fHistClusterdEtadPhiTracksP_000_075BeforeQA));

                        // |axis|<0.1
                        SetXRange(fHistClusterdEtadPhiTracksP_000_075BeforeQA,fHistClusterdEtadPhiTracksP_000_075BeforeQA->GetXaxis()->FindBin(-0.1),
                                  fHistClusterdEtadPhiTracksP_000_075BeforeQA->GetXaxis()->FindBin(0.1));
                        SetYRange(fHistClusterdEtadPhiTracksP_000_075BeforeQA,fHistClusterdEtadPhiTracksP_000_075BeforeQA->GetYaxis()->FindBin(-0.1),
                                  fHistClusterdEtadPhiTracksP_000_075BeforeQA->GetYaxis()->FindBin(0.1));
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksP_000_075BeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "with #it{p}_{T}<0.75 GeV/#it{c}");
                        SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_%sTracks_000_075_zoom_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                    } else cout << "INFO: Object |fHistClusterdEtadPhiTracksP_000_075BeforeQA| could not be found! Skipping Draw..." << endl;
                    //---------------------------------------------------------------------------------------------------------------
                    // dPhi dEta before QA cuts for a single charge tracks at medium track momentum (0.75< pT < 1.25)
                    TH2D* fHistClusterdEtadPhiTracksP_075_125BeforeQA = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_%sTracks_0.75<P<1.25_beforeClusterQA %s",
                                                                                                             charge[iCharge].Data(),fClusterCutSelection[i].Data()));
                    if(fHistClusterdEtadPhiTracksP_075_125BeforeQA){
                        SetZMinMaxTH2(fHistClusterdEtadPhiTracksP_075_125BeforeQA,1,fHistClusterdEtadPhiTracksP_075_125BeforeQA->GetNbinsX(),1,fHistClusterdEtadPhiTracksP_075_125BeforeQA->GetNbinsY());
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksP_075_125BeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "with 0.75<#it{p}_{T}<1.25 GeV/#it{c}");
                        SaveCanvasAndWriteHistogram( cvsQuadratic, fHistClusterdEtadPhiTracksP_075_125BeforeQA, Form("%s/dEtaVsdPhi_%sTracks_075_125_%s.%s", outputDir.Data(), charge[iCharge].Data(),
                                                     DataSets[i].Data(), suffix.Data()));
                        vecMatchingDeltaEtaPhiTracksAllPt[1][iCharge].push_back(new TH2D(*fHistClusterdEtadPhiTracksP_075_125BeforeQA));

                        // |axis|<0.1
                        SetXRange(fHistClusterdEtadPhiTracksP_075_125BeforeQA,fHistClusterdEtadPhiTracksP_075_125BeforeQA->GetXaxis()->FindBin(-0.1),
                                  fHistClusterdEtadPhiTracksP_075_125BeforeQA->GetXaxis()->FindBin(0.1));
                        SetYRange(fHistClusterdEtadPhiTracksP_075_125BeforeQA,fHistClusterdEtadPhiTracksP_075_125BeforeQA->GetYaxis()->FindBin(-0.1),
                                  fHistClusterdEtadPhiTracksP_075_125BeforeQA->GetYaxis()->FindBin(0.1));
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksP_075_125BeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "with 0.75<#it{p}_{T}<1.25 GeV/#it{c}");
                        SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_%sTracks_075_125_zoom_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                    } else cout << "INFO: Object |fHistClusterdEtadPhiTracksP_075_125BeforeQA| could not be found! Skipping Draw..." << endl;
                    //---------------------------------------------------------------------------------------------------------------
                    // dPhi dEta before QA cuts for a single charge tracks at medium high track momentum (1.25< pT < \infty)
                    TH2D* fHistClusterdEtadPhiTracksP_125_999BeforeQA = (TH2D*)CaloCutsContainer->FindObject(Form("dEtaVsdPhi_%sTracks_P>1.25_beforeClusterQA %s", charge[iCharge].Data(),
                                                                                                             fClusterCutSelection[i].Data()));
                    if(fHistClusterdEtadPhiTracksP_125_999BeforeQA){
                        SetZMinMaxTH2(fHistClusterdEtadPhiTracksP_125_999BeforeQA,1,fHistClusterdEtadPhiTracksP_125_999BeforeQA->GetNbinsX(),1,fHistClusterdEtadPhiTracksP_125_999BeforeQA->GetNbinsY());
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksP_125_999BeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "with #it{p}_{T}>1.25 GeV/#it{c}");
                        SaveCanvasAndWriteHistogram(cvsQuadratic, fHistClusterdEtadPhiTracksP_125_999BeforeQA, Form("%s/dEtaVsdPhi_%sTracks_125_999_%s.%s", outputDir.Data(), charge[iCharge].Data(),
                                                    DataSets[i].Data(), suffix.Data()));
                        vecMatchingDeltaEtaPhiTracksAllPt[2][iCharge].push_back(new TH2D(*fHistClusterdEtadPhiTracksP_125_999BeforeQA));

                        // |axis|<0.1
                        SetXRange(fHistClusterdEtadPhiTracksP_125_999BeforeQA,fHistClusterdEtadPhiTracksP_125_999BeforeQA->GetXaxis()->FindBin(-0.1),
                                  fHistClusterdEtadPhiTracksP_125_999BeforeQA->GetXaxis()->FindBin(0.1));
                        SetYRange(fHistClusterdEtadPhiTracksP_125_999BeforeQA,fHistClusterdEtadPhiTracksP_125_999BeforeQA->GetYaxis()->FindBin(-0.1),
                                  fHistClusterdEtadPhiTracksP_125_999BeforeQA->GetYaxis()->FindBin(0.1));
                        DrawPeriodQAHistoTH2(cvsQuadratic,leftMarginQuad,rightMarginQuad,topMarginQuad,bottomMarginQuad,kFALSE,kFALSE,kTRUE,
                                            fHistClusterdEtadPhiTracksP_125_999BeforeQA,"",
                                            "#Delta#eta","#Delta#phi",1,1.2,
                                            0.85,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                        PutProcessLabelAndEnergyOnPlot(0.15,0.94,0.03, fClusters.Data(), Form("and %s. tracks",charge[iCharge].Data()), "with #it{p}_{T}>1.25 GeV/#it{c}");
                        SaveCanvasOnly(cvsQuadratic, Form("%s/dEtaVsdPhi_%sTracks_125_999_zoom_%s.%s", outputDir.Data(), charge[iCharge].Data(), DataSets[i].Data(), suffix.Data()));
                    } else cout << "INFO: Object |fHistClusterdEtadPhiTracksP_125_999BeforeQA| could not be found! Skipping Draw..." << endl;

                }
             }
        }


        //---------------------------------------------------------------------------------------------------------------
        //--------------------------------------- Cell QA and Bad channel detection -------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        if(doExtQA>1){
            //---------------------------------------------------------------------------------------------------------------
            // frequency of cells to be included in a cluster before cluster QA
            TGaxis::SetExponentOffset(-0.06, -0.04, "x");
            TH1D* fHistClusterCellIDsBeforeQA = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_beforeClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterCellIDsBeforeQA){
                fHistClusterCellIDsBeforeQA->GetXaxis()->SetRangeUser(0,nCaloCells);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                    fHistClusterCellIDsBeforeQA,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                     "CellID in all Clusters","#frac{d#it{N}}{d#it{CellID}}",1,1);
                WriteHistogram(fHistClusterCellIDsBeforeQA);
            } else cout << "INFO: Object |ClusterIncludedCells_beforeClusterQA| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // frequency of cells to be included in a cluster after cluster QA
            TH1D* fHistClusterCellIDsAfterQA = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterIncludedCells_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterCellIDsAfterQA){
                fHistClusterCellIDsAfterQA->GetXaxis()->SetRangeUser(0,nCaloCells);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kTRUE,kFALSE,
                                     fHistClusterCellIDsAfterQA,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                     "Cell ID in accepted Clusters","#frac{d#it{N}}{d#it{CellID}}", 1, 1);
                WriteHistogram(fHistClusterCellIDsAfterQA);
                vecClusterIncludedCells.push_back(new TH1D(*fHistClusterCellIDsAfterQA));
            } else cout << "INFO: Object |ClusterIncludedCells_afterClusterQA| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // fraction energy of cells with respect to cluster energy before cluster QA
            TH1D* fHistClusterCellEFracBeforeQA = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyFracCells_beforeClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterCellEFracBeforeQA){
                fHistClusterCellEFracBeforeQA->Sumw2();
                fHistClusterCellIDsBeforeQA->Sumw2();
                fHistClusterCellEFracBeforeQA->Divide(fHistClusterCellIDsBeforeQA);
                fHistClusterCellEFracBeforeQA->GetXaxis()->SetRangeUser(0,nCaloCells);
                fHistClusterCellEFracBeforeQA->GetYaxis()->SetRangeUser(0,1);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterCellEFracBeforeQA,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                     "Energy Fraction of CellID in all Clusters","#frac{d#it{N}}{d#it{CellID}}",1,1);
                WriteHistogram(fHistClusterCellEFracBeforeQA);
            } else cout << "INFO: Object |ClusterEnergyFracCells_beforeClusterQA| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // fraction energy of cells with respect to cluster energy after cluster QA
            TH1D* fHistClusterCellEFracAfterQA = (TH1D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyFracCells_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterCellEFracAfterQA){
                fHistClusterCellEFracAfterQA->Sumw2();
                fHistClusterCellIDsAfterQA->Sumw2();
                fHistClusterCellEFracAfterQA->Divide(fHistClusterCellIDsAfterQA);
                fHistClusterCellEFracAfterQA->GetXaxis()->SetRangeUser(0,nCaloCells);
                fHistClusterCellEFracAfterQA->GetYaxis()->SetRangeUser(0,1);
                DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistClusterCellEFracAfterQA,Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()),
                                     "Energy Fraction of CellID in accepted Clusters","#frac{d#it{N}}{d#it{CellID}}",1,1);
                WriteHistogram(fHistClusterCellEFracAfterQA);
                vecClusterEnergyFracCells.push_back(new TH1D(*fHistClusterCellEFracAfterQA));
            } else cout << "INFO: Object |ClusterEnergyFracCells_afterClusterQA| could not be found! Skipping Draw..." << endl;


            TGaxis::SetExponentOffset(0, 0, "x");
            //---------------------------------------------------------------------------------------------------------------
            // cell ID vs energy and time of the cells
            TH2D* fHistCellEnergyVsCellID   =  (TH2D*)CaloExtQAContainer->FindObject(Form("CellEnergyVsCellID %s", fClusterCutSelection[i].Data()));
            TH2D* fHistCellTimeVsCellID     = (TH2D*)CaloExtQAContainer->FindObject(Form("CellTimeVsCellID %s", fClusterCutSelection[i].Data()));

            if(fHistCellEnergyVsCellID){
                // plot cell ID vs Energy
                fHistCellEnergyVsCellID->GetYaxis()->SetRangeUser(0,nCaloCells);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistCellEnergyVsCellID,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),
                                    "Cell Energy (GeV)","CellID",1,1);
                SaveCanvasAndWriteHistogram(canvas, fHistCellEnergyVsCellID, Form("%s/CellEnergyVsCellID_%s.%s", outputDir.Data(), DataSets[i].Data(),suffix.Data()));
                vecCellEnergyForComparison.push_back(new TH2D(*fHistCellEnergyVsCellID));

                // plot cell ID vs Energy - zoom in on low energies
                fHistCellEnergyVsCellID->GetXaxis()->SetRangeUser(0,2);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistCellEnergyVsCellID,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),
                                    "Cell Energy (GeV)","CellID",1,1);
                SaveCanvasAndWriteHistogram(canvas, fHistCellEnergyVsCellID, Form("%s/CellEnergyVsCellID_LowEnergy_%s.%s", outputDir.Data(), DataSets[i].Data(),suffix.Data()));

                // Calculate correlatio between mean and sigma of energy distribution per cell
                // if limits are set in cellQA object also the bad cell candidates for 2D energy vs sigma will be evaluated
                PlotCellMeanVsSigma(cellQA,nCaloCells,fHistCellEnergyVsCellID,
                                    "Mean Cell Energy (GeV)",
                                    "#sigma_{Cell Energy} (GeV)",
                                    0,0,0,
                                    0,0,0,kTRUE,isMC);
                if(!isMC && doCellQA){
                    line->DrawLine(cellQA->EnergyMean[1],cellQA->EnergySigma[0],cellQA->EnergyMean[1],cellQA->EnergySigma[1]);
                    line->DrawLine(cellQA->EnergyMean[0],cellQA->EnergySigma[1],cellQA->EnergyMean[1],cellQA->EnergySigma[1]);
                }
                PutProcessLabelAndEnergyOnPlot(0.87, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(),42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/CellEnergyVsSigma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kFALSE, kTRUE);

                // estimate hot cells from energy distribution
                if(mode == 3){
                  PlotHotCells(cellQA,1,nCaloCells,fHistCellEnergyVsCellID,
                            "N_{Cell fired with Energy L to 30 GeV} / N_{Cell fired with Energy 0 to L GeV}",
                            "Integration Limit 'L' (GeV)",
                            1,1E-8,1,
                            0,0,0,isMC,0.95);
                }else{
                  PlotHotCells(cellQA,1,nCaloCells,fHistCellEnergyVsCellID,
                            "N_{Cell fired with Energy L to 30 GeV} / N_{Cell fired with Energy 0 to L GeV}",
                            "Integration Limit 'L' (GeV)",
                            0,0,0,
                            0,0,0,isMC,0.95);
                }
                if(!isMC && doCellQA){
                    const Int_t dim2D= 9;
                    for(Int_t ii=0; ii<dim2D; ii++){
                        cout << ii << "\t" << cellQA->HotCells2D[ii][0] << "\t" << cellQA->HotCells2D[ii][1] << endl;
                        line->DrawLine(cellQA->HotCells2D[ii][0],0.2+ii*0.1,cellQA->HotCells2D[ii][0],0.3+ii*0.1);
                        line->DrawLine(cellQA->HotCells2D[ii][1],0.2+ii*0.1,cellQA->HotCells2D[ii][1],0.3+ii*0.1);
                    }
                }
                PutProcessLabelAndEnergyOnPlot(0.87, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), 42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/CellHotCells2D_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kFALSE, kTRUE);

                if (mode==5 || mode == 4 || mode == 3) canvas->SetRightMargin(rightMargin*3.);
                else canvas->SetRightMargin(rightMargin);
                PlotHotCells(cellQA,0,nCaloCells,fHistCellTimeVsCellID,
                            "N_{Cell fired} / N_{Events}",
                            "#frac{dN_{Cell fired with E>0.2 GeV}}{dN}",
                            0,0,0,
                            0,0,0,isMC,1,1,kFALSE,nEvents[i]);
                if(!isMC && doCellQA){
                    line->DrawLine(cellQA->HotCells1D[0],0,cellQA->HotCells1D[0],10);
                    line->DrawLine(cellQA->HotCells1D[1],0,cellQA->HotCells1D[1],10);
                }
                PutProcessLabelAndEnergyOnPlot(0.9, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(),42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/CellHotCells_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE);

                if (mode==5 || mode == 3) canvas->SetRightMargin(rightMargin*3.);
                else canvas->SetRightMargin(rightMargin);
                PlotHotCells(cellQA,0,nCaloCells,fHistCellTimeVsCellID,
                             "N_{Cell fired} / N_{Events} ",
                             "#frac{dN_{Cell fired with E>0.2 GeV}}{dN}",
                             0,0,0,
                             0,0,0,isMC,
                             1,1,kTRUE, nEvents[i]
                            );
                if(!isMC && doCellQA){
                  line->DrawLine(cellQA->HotCells1D[0],0,cellQA->HotCells1D[0],10);
                  line->DrawLine(cellQA->HotCells1D[1],0,cellQA->HotCells1D[1],10);
                }
                PutProcessLabelAndEnergyOnPlot(0.9, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(),42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/CellHotCellsRescaled_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE);

                TGaxis::SetMaxDigits(3);

            } else cout << "INFO: Object |CellEnergyVsCellID| could not be found! Skipping Draw..." << endl;

            // saved as .jpg per default due to size of histogram
            TGaxis::SetExponentOffset(0.5, 0, "x");
            //---------------------------------------------------------------------------------------------------------------
            // cell ID vs energy and time of the cells
            if(fHistCellTimeVsCellID){
                // cell ID vs energy and time of the cells
                fHistCellTimeVsCellID->GetYaxis()->SetRangeUser(0,nCaloCells);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistCellTimeVsCellID,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),
                                    "Cell Time (#mus)","CellID",1,1);
                SaveCanvasAndWriteHistogram(canvas, fHistCellTimeVsCellID, Form("%s/CellTimeVsCellID_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                vecCellTimingForComparison.push_back(new TH2D(*fHistCellTimeVsCellID));

                // Calculate correlatio between mean and sigma of time distribution per cell
                // if limits are set in cellQA object also the bad cell candidates for 2D energy vs sigma will be evaluated
                PlotCellMeanVsSigma(cellQA,nCaloCells,fHistCellTimeVsCellID,
                                    "Mean Cell Time (#mus)",
                                    "#sigma_{Cell Time} (#mus)",
                                    0,0,0,
                                    0,0,0,kFALSE,isMC);
                if(!isMC && doCellQA){
                    line->DrawLine(cellQA->TimeMean[0],cellQA->TimeSigma[0],cellQA->TimeMean[1],cellQA->TimeSigma[0]);
                    line->DrawLine(cellQA->TimeMean[0],cellQA->TimeSigma[0],cellQA->TimeMean[0],cellQA->TimeSigma[1]);
                    line->DrawLine(cellQA->TimeMean[1],cellQA->TimeSigma[0],cellQA->TimeMean[1],cellQA->TimeSigma[1]);
                    line->DrawLine(cellQA->TimeMean[0],cellQA->TimeSigma[1],cellQA->TimeMean[1],cellQA->TimeSigma[1]);
                }
                TGaxis::SetExponentOffset(0, 0.5, "y");
                PutProcessLabelAndEnergyOnPlot(0.87, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(),42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/CellTimeVsSigma_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kFALSE, kTRUE);

                // estimate Hot cells based on timing information
                canvas->SetRightMargin(rightMargin);
                PlotHotCells(cellQA,2,nCaloCells,fHistCellTimeVsCellID,
                            "N_{Cell fired} / N_{Cell fired at |t|< 0.02 (#mus)} - E_{Cell}>0.2 GeV",
                            "# of Entries",
                            0,0,0,
                            0,0,0,isMC,0.95);
                if(!isMC && doCellQA){
                    line->DrawLine(cellQA->HotCellsTime1D[0],0,cellQA->HotCellsTime1D[0],10);
                    line->DrawLine(cellQA->HotCellsTime1D[1],0,cellQA->HotCellsTime1D[1],10);
                }
                PutProcessLabelAndEnergyOnPlot(0.95, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(),42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/CellHotCellsTime1D_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE, kFALSE);
            } else cout << "INFO: Object |CellTimeVsCellID| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // plot bad channels masked in the data already
            TProfile* fEMCalBadChannels = (TProfile*)CaloExtQAContainer->FindObject(Form("%s - Bad Channels",calo.Data()));
            if(fEMCalBadChannels){
                fEMCalBadChannels->SetTitle(Form("%s - %s - %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fClusters.Data()));
                fEMCalBadChannels->GetXaxis()->SetTitle("Cell ID");
                fEMCalBadChannels->GetYaxis()->SetTitle("Cell ID Bad in Fraction of Events");
                fEMCalBadChannels->GetYaxis()->SetLabelFont(42);
                fEMCalBadChannels->GetXaxis()->SetLabelFont(42);
                fEMCalBadChannels->GetYaxis()->SetTitleFont(62);
                fEMCalBadChannels->GetXaxis()->SetTitleFont(62);
                fEMCalBadChannels->GetYaxis()->SetLabelSize(0.035);
                fEMCalBadChannels->GetYaxis()->SetTitleSize(0.043);
                fEMCalBadChannels->GetYaxis()->SetDecimals();
                fEMCalBadChannels->GetXaxis()->SetTitleSize(0.043);
                fEMCalBadChannels->GetXaxis()->SetLabelSize(0.035);
                fEMCalBadChannels->DrawCopy();
                WriteHistogram(fEMCalBadChannels);
                SaveCanvas(canvas, Form("%s/BadChannels_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
            } else cout << "INFO: Object |Bad Channels| could not be found! Skipping Draw..." << endl;
        } //end only plotting for doExtQA>1

        //---------------------------------------------------------------------------------------------------------------
        //--------------------------------- SuperModule QA --------------------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        if (doExtQA > 1){  // only do if doExtQA > 1
            TGaxis::SetExponentOffset(0, 0, "x");
            //---------------------------------------------------------------------------------------------------------------
            //--------------------------------- SuperModule QA --------------------------------------------------------------
            //---------------------------------------------------------------------------------------------------------------
            // cluster energy per super module number
            TH2D* fHistClusterEnergyVsModule = (TH2D*)CaloExtQAContainer->FindObject(Form("ClusterEnergyVsModule_afterClusterQA %s", fClusterCutSelection[i].Data()));
            if(fHistClusterEnergyVsModule){
                GetMinMaxBin(fHistClusterEnergyVsModule,minB,maxB);
                SetXRange(fHistClusterEnergyVsModule,minB,maxB+5);
                fHistClusterEnergyVsModule->GetYaxis()->SetRangeUser(0,nCaloModules);
                SetZMinMaxTH2(fHistClusterEnergyVsModule,minB,maxB+5,1,nCaloModules+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistClusterEnergyVsModule,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),
                                    "Cluster Energy (GeV)","SuperModule Number",1,1);
                SaveCanvasAndWriteHistogram(canvas, fHistClusterEnergyVsModule, Form("%s/ClusterEnergyVsModule_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                canvas->SetRightMargin(rightMargin);
                PlotCaloQAModule(fHistClusterEnergyVsModule,
                                nCaloModules,
                                "Cluster Energy (GeV)",
                                "#frac{dE}{dN}",
                                1,10,1,0.1,kTRUE,minB,maxB,kTRUE);
                PutProcessLabelAndEnergyOnPlot(0.95, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(),42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/ClusterEnergyVsModule_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
            } else cout << "INFO: Object |ClusterEnergyVsModule_afterClusterQA| could not be found! Skipping Draw..." << endl;


            //---------------------------------------------------------------------------------------------------------------
            // summed energy per super module
            TH2D* fHistModuleEnergyVsModule = (TH2D*)CaloExtQAContainer->FindObject(Form("ModuleEnergyVsModule %s", fClusterCutSelection[i].Data()));
            if(fHistModuleEnergyVsModule){
                GetMinMaxBin(fHistModuleEnergyVsModule,minB,maxB);
                SetXRange(fHistModuleEnergyVsModule,minB,maxB+5);
                fHistModuleEnergyVsModule->GetYaxis()->SetRangeUser(0,nCaloModules);
                SetZMinMaxTH2(fHistModuleEnergyVsModule,minB,maxB+5,1,nCaloModules+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistModuleEnergyVsModule,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),
                                    "Total SuperModule Energy per Event (GeV)","SuperModule Number",1,1);
                SaveCanvasAndWriteHistogram(canvas, fHistModuleEnergyVsModule, Form("%s/ModuleEnergyVsModule_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                canvas->SetRightMargin(rightMargin);
                PlotCaloQAModule(fHistModuleEnergyVsModule,
                                nCaloModules,
                                "Total SuperModule Energy per Event (GeV)",
                                "#frac{dE}{dN}",
                                1,10,1,0.1,kTRUE,minB,maxB+5,kFALSE);
                PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data());
                SaveCanvas(canvas, Form("%s/ModuleEnergyVsModule_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE);
                PlotCaloQAModule(fHistModuleEnergyVsModule,
                                nCaloModules,
                                "Total SuperModule Energy per Event (GeV)",
                                "#frac{dE}{dN}",
                                1,10,1,0.1,kTRUE,2,maxB+5,kFALSE);
                PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data());
                SaveCanvas(canvas, Form("%s/ModuleEnergyVsModule_Projected_%s_LOG.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kTRUE, kTRUE);
            } else cout << "INFO: Object |ModuleEnergyVsModule_| could not be found! Skipping Draw..." << endl;

            //---------------------------------------------------------------------------------------------------------------
            // number of cells above 100 MeV in a super module
            TH2D* fHistNCellsAbove100VsModule = (TH2D*)CaloExtQAContainer->FindObject(Form("NCellsAbove100VsModule %s", fClusterCutSelection[i].Data()));
            if(fHistNCellsAbove100VsModule){
                GetMinMaxBin(fHistNCellsAbove100VsModule,minB,maxB);
                SetXRange(fHistNCellsAbove100VsModule,minB,maxB+5);
                fHistNCellsAbove100VsModule->GetYaxis()->SetRangeUser(0,nCaloModules);
                SetZMinMaxTH2(fHistNCellsAbove100VsModule,minB,maxB+5,1,nCaloModules+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistNCellsAbove100VsModule,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),
                                    "#it{N}_{Cells}>100 MeV in SM per Event","SuperModule Number",1,1);
                SaveCanvasAndWriteHistogram(canvas, fHistNCellsAbove100VsModule, Form("%s/NCellsAbove100VsModule_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                canvas->SetRightMargin(rightMargin);
                PlotCaloQAModule(fHistNCellsAbove100VsModule,
                                nCaloModules,
                                "#it{N}_{Cells}>100 MeV in SM per Event",
                                "#frac{d#it{N}_{Cells}>100 MeV}{dN}",
                                1,10,1,0.1,kTRUE,minB,maxB,kFALSE);
                PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data());
                SaveCanvas(canvas, Form("%s/NCellsAbove100VsModule_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE);
            } else cout << "INFO: Object |NCellsAbove100VsModule| could not be found! Skipping Draw..." << endl;
            //---------------------------------------------------------------------------------------------------------------
            // number of cells above 1.5 GeV in a super module
            TH2D* fHistNCellsAbove1500VsModule = (TH2D*)CaloExtQAContainer->FindObject(Form("NCellsAbove1500VsModule %s", fClusterCutSelection[i].Data()));
            if(fHistNCellsAbove1500VsModule){
                GetMinMaxBin(fHistNCellsAbove1500VsModule,minB,maxB);
                SetXRange(fHistNCellsAbove1500VsModule,minB,maxB+5);
                fHistNCellsAbove1500VsModule->GetYaxis()->SetRangeUser(0,nCaloModules);
                SetZMinMaxTH2(fHistNCellsAbove1500VsModule,minB,maxB+5,1,nCaloModules+1);
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kTRUE,
                                    fHistNCellsAbove1500VsModule,Form("%s - %s %s- %s",fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data(), fClusters.Data()),
                                    "#it{N}_{Cells}>1500 MeV in SM per Event","SuperModule Number",1,1);
                SaveCanvasAndWriteHistogram(canvas, fHistNCellsAbove1500VsModule, Form("%s/NCellsAbove1500VsModule_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                canvas->SetRightMargin(rightMargin);
                PlotCaloQAModule(fHistNCellsAbove1500VsModule,
                                nCaloModules,
                                "#it{N}_{Cells}>1500 MeV in SM per Event",
                                "#frac{d#it{N}_{Cells}>1500 MeV}{dN}",
                                1,10,1,0.1,kTRUE,minB,maxB,kFALSE);
                PutProcessLabelAndEnergyOnPlot(0.8, 0.92, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTrigger[i].Data());
                SaveCanvas(canvas, Form("%s/NCellsAbove1500VsModule_Projected_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE);
            } else cout << "INFO: Object |NCellsAbove1500VsModule| could not be found! Skipping Draw..." << endl;
        } //end only plotting for doExtQA>1


        //---------------------------------------------------------------------------------------------------------------
        //------------------------- detailed QA for conv-calo properties ------------------------------------------------
        //---------------------------------------------------------------------------------------------------------------
        if(isPCMCalo){
            //---------------------------------------------------------------------------------------------------------------
            //------------------- distribution of conversion photons contribution to meson peaks ----------------------------
            //---------------------------------------------------------------------------------------------------------------
            TH2D* MotherPi0Conv_Eta_Phi_ESD = (TH2D*) ESDContainer->FindObject("ESD_MotherPi0ConvPhoton_Eta_Phi");
            TH2D* MotherEtaConv_Eta_Phi_ESD = (TH2D*) ESDContainer->FindObject("ESD_MotherEtaConvPhoton_Eta_Phi");

            //---------------------------------------------------------------------------------------------------------------
            //------------------- distribution of conversion photons contribution to pi0 peaks ------------------------------
            //---------------------------------------------------------------------------------------------------------------
            if(MotherPi0Conv_Eta_Phi_ESD){
                //---------------------------------------------------------------------------------------------------------------
                // 2 dimensional distribution of conversion photons contribution to pi0 peak in eta, phi
                TH2D* MotherPi0Conv_Eta_Phi = (TH2D*) MotherPi0Conv_Eta_Phi_ESD->Clone();
                GetMinMaxBin(MotherPi0Conv_Eta_Phi,minB,maxB);
                SetXRange(MotherPi0Conv_Eta_Phi,minB-5,maxB+5);
                //MotherPi0Conv_Eta_Phi->GetXaxis()->SetRangeUser(0.8,3.8);
                SetZMinMaxTH2(MotherPi0Conv_Eta_Phi,MotherPi0Conv_Eta_Phi->GetXaxis()->FindBin(0.8),MotherPi0Conv_Eta_Phi->GetXaxis()->FindBin(3.8),1,MotherPi0Conv_Eta_Phi->GetNbinsY());
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    MotherPi0Conv_Eta_Phi,"",
                                    "#phi_{#gamma_{conv} under #pi^{0}-peak}","#eta_{#gamma_{conv} under #pi^{0}-peak}",1,1,
                                    0.75,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(canvas, MotherPi0Conv_Eta_Phi, Form("%s/ConvPhotonPi0_Eta_Phi_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                //---------------------------------------------------------------------------------------------------------------
                // distribution of conversion photons contribution to pi0 peak in phi
                TH1D* MotherPi0Conv_Phi = (TH1D*) MotherPi0Conv_Eta_Phi->ProjectionX("MotherPi0Conv_Phi",1,MotherPi0Conv_Eta_Phi->GetNbinsY());
                if(MotherPi0Conv_Phi){
                    GetMinMaxBin(MotherPi0Conv_Phi,minB,maxB);
                    SetXRange(MotherPi0Conv_Phi,minB-5,maxB+5);
                    DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                        MotherPi0Conv_Phi,"","#phi_{#gamma_{conv} under #pi^{0}-peak}","#frac{d#it{N}}{d#eta}",1,1,
                                        0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, MotherPi0Conv_Phi, Form("%s/ConvPhotonPi0_Phi_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    vecConvPhotonEtaPhi_Pi0[0].push_back(MotherPi0Conv_Phi);
                }
                //---------------------------------------------------------------------------------------------------------------
                // distribution of conversion photons contribution to pi0 peak in phi
                TH1D* MotherPi0Conv_Eta = (TH1D*) MotherPi0Conv_Eta_Phi->ProjectionY("MotherPi0Conv_Eta",1,MotherPi0Conv_Eta_Phi->GetNbinsX());
                if(MotherPi0Conv_Eta){
                    MotherPi0Conv_Eta->GetXaxis()->SetRangeUser(-1.,1.);
                    DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                        MotherPi0Conv_Eta,"","#eta_{#gamma_{conv} under #pi^{0}-peak}","#frac{d#it{N}}{d#phi}",1,1,
                                        0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, MotherPi0Conv_Eta, Form("%s/ConvPhotonPi0_Eta_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    vecConvPhotonEtaPhi_Pi0[1].push_back(MotherPi0Conv_Eta);
                }
                delete MotherPi0Conv_Eta_Phi;
            } else cout << "INFO: Object |MotherPi0Conv_Eta| could not be found! Skipping Draw..." << endl;

            //---------------------------------------------------------------------------------------------------------------
            //------------------- distribution of conversion photons contribution to eta peaks ------------------------------
            //---------------------------------------------------------------------------------------------------------------
            if(MotherEtaConv_Eta_Phi_ESD){
                //---------------------------------------------------------------------------------------------------------------
                // 2 dimensional distribution of conversion photons contribution to eta peak in eta, phi
                TH2D* MotherEtaConv_Eta_Phi = (TH2D*) MotherEtaConv_Eta_Phi_ESD->Clone();
                SetZMinMaxTH2(MotherEtaConv_Eta_Phi,1,MotherEtaConv_Eta_Phi->GetNbinsX(),1,MotherEtaConv_Eta_Phi->GetNbinsY());
                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    MotherEtaConv_Eta_Phi,"",
                                    "#phi_{#gamma_{conv} under #eta-peak}","#eta_{#gamma_{conv} under #eta-peak}",1,1,
                                    0.75,0.95,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i]);
                SaveCanvasAndWriteHistogram(canvas, MotherEtaConv_Eta_Phi, Form("%s/ConvPhotonEta_Eta_Phi_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                //---------------------------------------------------------------------------------------------------------------
                // distribution of conversion photons contribution to eta peak in phi
                TH1D* MotherEtaConv_Phi = (TH1D*) MotherEtaConv_Eta_Phi->ProjectionX("MotherEtaConv_Phi",1,MotherEtaConv_Eta_Phi->GetNbinsY());
                if(MotherEtaConv_Phi){
                    DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                        MotherEtaConv_Phi,"","#phi_{#gamma_{conv} under #eta-peak}","#frac{d#it{N}}{d#phi}",1,1,
                                        0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, MotherEtaConv_Phi, Form("%s/ConvPhotonEta_Phi_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    vecConvPhotonEtaPhi_Eta[0].push_back(MotherEtaConv_Phi);
                }
                //---------------------------------------------------------------------------------------------------------------
                // distribution of conversion photons contribution to eta peak in phi
                TH1D* MotherEtaConv_Eta = (TH1D*) MotherEtaConv_Eta_Phi->ProjectionY("MotherEtaConv_Eta",1,MotherEtaConv_Eta_Phi->GetNbinsX());
                if(MotherEtaConv_Eta){
                    MotherEtaConv_Eta->GetXaxis()->SetRangeUser(-1.,1.);
                    DrawPeriodQAHistoTH1(canvas,leftMargin,rightMargin,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                        MotherEtaConv_Eta,"","#eta_{#gamma_{conv} under #eta-peak}","#frac{d#it{N}}{d#eta}",1,1,
                                        0.95,0.94,0.03,fCollisionSystem,plotDataSets[i],fTrigger[i],31);
                    SaveCanvasAndWriteHistogram(canvas, MotherEtaConv_Eta, Form("%s/ConvPhotonEta_Eta_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));
                    vecConvPhotonEtaPhi_Eta[1].push_back(MotherEtaConv_Eta);
                }
                delete MotherEtaConv_Eta_Phi;
            } else cout << "INFO: Object |MotherEtaConv_Eta| could not be found! Skipping Draw..." << endl;

            //---------------------------------------------------------------------------------------------------------------
            //--------------------------- Read histograms for matching efficiency & fake rate for V0 matching ---------------
            //---------------------------------------------------------------------------------------------------------------
            // reconstructed mass spectrum
            TH2D* ESDMother = (TH2D*) ESDContainer->FindObject("ESD_Mother_InvMass_Pt");
            if(ESDMother){
                WriteHistogram(ESDMother);
                vecESDMother.push_back(new TH2D(*ESDMother));
            } else {cout << "ERROR: Object |ESD_Mother_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; return;}
            //---------------------------------------------------------------------------------------------------------------
            // rejected meson candidates due to V0 matching
            TH2D* ESDMotherMatched = (TH2D*) ESDContainer->FindObject("ESD_MotherMatched_InvMass_Pt");
            if(ESDMotherMatched){
                WriteHistogram(ESDMotherMatched);
                vecESDMotherMatched.push_back(new TH2D(*ESDMotherMatched));
            } else {cout << "ERROR: Object |ESD_MotherMatched_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; return;}

            //---------------------------------------------------------------------------------------------------------------
            //reconstructed validated
            if(isTrueContainer[i]){
                for(Int_t iM=0; iM<2; iM++){
                    //---------------------------------------------------------------------------------------------------------------
                    // validated meson
                    TH2D* ESDTrueMeson = (TH2D*) TrueContainer->FindObject(Form("ESD_True%s_InvMass_Pt",meson[iM].Data()));
                    if(ESDTrueMeson){
                        WriteHistogram(ESDTrueMeson);
                        vecESDTrueMeson[iM].push_back(new TH2D(*ESDTrueMeson));
                    } else {cout << "INFO: Object |ESD_True%s_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; isTrueContainer[i] = kFALSE;}
                    //---------------------------------------------------------------------------------------------------------------
                    // falsely matched validated mesons
                    TH2D* ESDTrueMesonMatched = (TH2D*) TrueContainer->FindObject(Form("ESD_True%s_Matched_InvMass_Pt",meson[iM].Data()));
                    if(ESDTrueMesonMatched){
                        WriteHistogram(ESDTrueMesonMatched);
                        vecESDTrueMesonMatched[iM].push_back(new TH2D(*ESDTrueMesonMatched));
                    } else {cout << "INFO: Object |ESD_True%s_Matched_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; isTrueContainer[i] = kFALSE;}
                    //---------------------------------------------------------------------------------------------------------------
                    //
                    TH2D* ESDTrueMesonCaloConvertedPhotonMatched = (TH2D*) TrueContainer->FindObject(Form("ESD_True%sCaloConvertedPhotonMatched_InvMass_Pt",meson[iM].Data()));
                    if(ESDTrueMesonCaloConvertedPhotonMatched){
                        WriteHistogram(ESDTrueMesonCaloConvertedPhotonMatched);
                        vecESDTrueMesonCaloConvertedPhotonMatched[iM].push_back(new TH2D(*ESDTrueMesonCaloConvertedPhotonMatched));
                    } else {cout << "INFO: Object |ESD_True%sCaloConvertedPhotonMatched_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; isTrueContainer[i] = kFALSE;}
                    //---------------------------------------------------------------------------------------------------------------
                    // meson reconstructed with cluster from conversion
                    TH2D* ESDTrueMesonCaloConvertedPhoton = (TH2D*) TrueContainer->FindObject(Form("ESD_True%sCaloConvertedPhoton_InvMass_Pt",meson[iM].Data()));
                    if(ESDTrueMesonCaloConvertedPhoton){
                        WriteHistogram(ESDTrueMesonCaloConvertedPhoton);
                        vecESDTrueMesonCaloConvertedPhoton[iM].push_back(new TH2D(*ESDTrueMesonCaloConvertedPhoton));
                    } else {cout << "INFO: Object |ESD_True%sCaloConvertedPhoton_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; isTrueContainer[i] = kFALSE;}
                    //---------------------------------------------------------------------------------------------------------------
                    // meson reconstructed with calorimeter photon
                    TH2D* ESDTrueMesonCaloPhoton = (TH2D*) TrueContainer->FindObject(Form("ESD_True%sCaloPhoton_InvMass_Pt",meson[iM].Data()));
                    if(ESDTrueMesonCaloPhoton){
                        WriteHistogram(ESDTrueMesonCaloPhoton);
                        vecESDTrueMesonCaloPhoton[iM].push_back(new TH2D(*ESDTrueMesonCaloPhoton));
                    } else {cout << "INFO: Object |ESD_True%sCaloPhoton_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; isTrueContainer[i] = kFALSE;}
                    //---------------------------------------------------------------------------------------------------------------
                    // meson reconstructed with calorimeter electron
                    TH2D* ESDTrueMesonCaloElectron = (TH2D*) TrueContainer->FindObject(Form("ESD_True%sCaloElectron_InvMass_Pt",meson[iM].Data()));
                    if(ESDTrueMesonCaloElectron){
                        WriteHistogram(ESDTrueMesonCaloElectron);
                        vecESDTrueMesonCaloElectron[iM].push_back(new TH2D(*ESDTrueMesonCaloElectron));
                    } else {cout << "INFO: Object |ESD_True%sCaloElectron_InvMass_Pt| could not be found! Skipping Draw & return..." << endl; isTrueContainer[i] = kFALSE;}
                }
            }
            if(!isTrueContainer[i]){
                for(Int_t iM=0; iM<2; iM++){
                    vecESDTrueMeson[iM].push_back(0x0);
                    vecESDTrueMesonMatched[iM].push_back(0x0);
                    vecESDTrueMesonCaloConvertedPhotonMatched[iM].push_back(0x0);
                    vecESDTrueMesonCaloConvertedPhoton[iM].push_back(0x0);
                    vecESDTrueMesonCaloPhoton[iM].push_back(0x0);
                    vecESDTrueMesonCaloElectron[iM].push_back(0x0);
                }
            }
        }

        fFile->Close();
        delete fFile;

        delete TopDir;

        fOutput->Close();
        delete fOutput;
    }

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Comparison Histograms ************************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms" << endl;

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms for Event qualities"<< endl;
    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------- compare event quality histos -----------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    for(Int_t iVec=0; iVec<(Int_t)vecNEvents.size(); iVec++){
        TH1D* temp = vecNEvents.at(iVec);
        temp->Sumw2();
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                vecNEvents,"","","N_{Events}",1,1.1,
                                labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/NEvents.%s", outputDir.Data(), suffix.Data()));

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecNEvents,"","","N_{Events}",1,1.1,
                                    labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_NEvents.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecNEvents);

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms for Timing"  << endl;
    //---------------------------------------------------------------------------------------------------------------
    //-------------------------------- Cluster properties: timing ---------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // shift MC to data
    TGaxis::SetExponentOffset(-0.08, 0.01, "x");
    TH1D* mcClusterTime[nSets];
    std::vector<TH1D*> vecClusterTimeAfterCorrected;
    for(Int_t iVec=0; iVec<nSets; iVec++){
        TH1D* temp = vecClusterTimeAfter.at(iVec);
        if(iVec==0){
            for(Int_t ii=1; ii<nSets; ii++){
                mcClusterTime[ii] = (TH1D*) temp->Clone(Form("clone_mc%iClusterTimeAfter",ii));
                mcClusterTime[ii]->Reset("ICE");
            }
            temp->Sumw2();
            Double_t nMax = temp->GetBinContent(temp->GetMaximumBin());
            temp->Scale(1./nMax);
            vecClusterTimeAfterCorrected.push_back(new TH1D(*temp));
        } else {
            Double_t meanMC = GetMaximumBinValueTH1(temp);//temp->GetMean(1);
            Double_t meanData = GetMaximumBinValueTH1(vecClusterTimeAfter.at(0));
            for(Int_t iBin=1; iBin<temp->GetXaxis()->GetNbins(); iBin++){
                mcClusterTime[iVec]->SetBinContent(temp->GetXaxis()->FindBin(temp->GetBinCenter(iBin) - meanMC + meanData),
                                                temp->GetBinContent(iBin));
            }
            mcClusterTime[iVec]->Sumw2();
            Double_t nMax = mcClusterTime[iVec]->GetBinContent(mcClusterTime[iVec]->GetMaximumBin());
            mcClusterTime[iVec]->Scale(1/nMax);
            vecClusterTimeAfterCorrected.push_back(mcClusterTime[iVec]);
        }
    }
    //---------------------------------------------------------------------------------------------------------------
    // cluster time after cuts
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterTimeAfterCorrected,"","#it{t}_{Cluster} (s)","#frac{1}{N_{max}} #frac{dN}{dt}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/ClusterTimeAfter.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
    DeleteVecTH1D(vecClusterTimeAfterCorrected);
    DeleteVecTH1D(vecClusterTimeAfter);

    //---------------------------------------------------------------------------------------------------------------
    // shift MC to data
    std::vector<TH1D*> vecClusterTimeBeforeCorrected;
    for(Int_t iVec=0; iVec<(Int_t) vecClusterTimeBefore.size(); iVec++){
        TH1D* temp = vecClusterTimeBefore.at(iVec);
        if(iVec==0){
        for(Int_t ii=1; ii<nSets; ii++){
            mcClusterTime[ii] = (TH1D*) temp->Clone(Form("clone_mc%iClusterTimeBefore",ii));
            mcClusterTime[ii]->Reset("ICE");
        }
        temp->Sumw2();
        Double_t nMax = temp->GetBinContent(temp->GetMaximumBin());
        temp->Scale(1./nMax);
        vecClusterTimeBeforeCorrected.push_back(new TH1D(*temp));
        } else{
        Double_t meanMC = GetMaximumBinValueTH1(temp); //temp->GetMean(1);
        Double_t meanData = GetMaximumBinValueTH1(vecClusterTimeBefore.at(0));
        for(Int_t iBin=1; iBin<temp->GetXaxis()->GetNbins(); iBin++){
            mcClusterTime[iVec]->SetBinContent(temp->GetXaxis()->FindBin(temp->GetBinCenter(iBin) - meanMC + meanData),
                                            temp->GetBinContent(iBin));
        }
        mcClusterTime[iVec]->Sumw2();
        Double_t nMax = mcClusterTime[iVec]->GetBinContent(mcClusterTime[iVec]->GetMaximumBin());
        mcClusterTime[iVec]->Scale(1/nMax);
        vecClusterTimeBeforeCorrected.push_back(mcClusterTime[iVec]);
        }
    }
    //---------------------------------------------------------------------------------------------------------------
    // cluster time before cuts
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterTimeBeforeCorrected,"","#it{t}_{Cluster} (s)","#frac{1}{N_{max}} #frac{dN}{dt}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/ClusterTimeBefore.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
    DeleteVecTH1D(vecClusterTimeBeforeCorrected);
    DeleteVecTH1D(vecClusterTimeBefore);
    TGaxis::SetExponentOffset(0, 0, "x");

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------- detailed information on cell/cluster time ----------------------------------
    //---------------------------------------------------------------------------------------------------------------
    if(doExtQA>1){
        TGaxis::SetExponentOffset(0, 0, "x");

        // cell time in cluster
        Double_t yStart=0; Double_t yEnd=0;
        for(Int_t iVec=0; iVec<(Int_t)vecClusterTimingInCluster.size(); iVec++){
            TH1D* temp = vecClusterTimingInCluster.at(iVec);
            temp->Sumw2();
            Double_t nEntries = temp->Integral();
            temp->Scale(1./nEntries);
            if(iVec==0){
                yStart = temp->GetMinimum();
                yEnd = temp->GetMaximum();
            }
        }
        //---------------------------------------------------------------------------------------------------------------
        // cell time in cluster
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterTimingInCluster,"","Cell Time in Clusters (ns)","#frac{1}{N} #frac{dN}{d#it{t_{cell}}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ClusterIncludedCellTime.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterTimingInCluster,"","Cell Time in Clusters (ns)","#frac{1}{N} #frac{dN}{d#it{t_{cell}}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        line->DrawLine(-50, yStart, -50, yEnd);
        line->DrawLine(50, yStart, 50, yEnd);
        SaveCanvas(canvas, Form("%s/Comparison/ClusterIncludedCellTime_50ns.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecClusterTimingInCluster,"","Cell Time in Clusters (ns)","#frac{1}{N} #frac{dN}{d#it{t_{cell}}}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ClusterIncludedCellTime.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecClusterTimingInCluster);

        //---------------------------------------------------------------------------------------------------------------
        // distance of clusters within same time window
        for(Int_t iVec=0; iVec<(Int_t)vecClusterDistanceRowWithin.size(); iVec++){
            TH1D* temp = vecClusterDistanceRowWithin.at(iVec);
            temp->Sumw2();
            Double_t nEntries = temp->Integral();
            temp->Scale(1./nEntries);
        }
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecClusterDistanceRowWithin,"","#Delta row","#frac{1}{N} #frac{dN}{d#Delta row}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 1.1, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ClusterDistanceRowWithin.%s", outputDir.Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        vecClusterDistanceRowWithin,"","#Delta row","#frac{1}{N} #frac{dN}{d#Delta row}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 1.1, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ClusterDistanceRowWithin.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecClusterDistanceRowWithin);

        for(Int_t iVec=0; iVec<(Int_t)vecClusterDistanceColumnWithin.size(); iVec++){
            TH1D* temp = vecClusterDistanceColumnWithin.at(iVec);
            temp->Sumw2();
            Double_t nEntries = temp->Integral();
            temp->Scale(1./nEntries);
        }
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                    vecClusterDistanceColumnWithin,"","#Delta column","#frac{1}{N} #frac{dN}{d#Delta column}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 1.1, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ClusterDistanceColumnWithin.%s", outputDir.Data(), suffix.Data()));

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                        vecClusterDistanceColumnWithin,"","#Delta column","#frac{1}{N} #frac{dN}{d#Delta column}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 1.1, kFALSE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ClusterDistanceColumnWithin.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecClusterDistanceColumnWithin);

        for(Int_t iEM02 = 0; iEM02<nEM02_bins-1; iEM02++){
          for(Int_t iVec=0; iVec<(Int_t)vecClusterEM02.size(); iVec++){
              TH2D* temp = vecClusterEM02.at(iVec);
              temp->GetXaxis()->SetRangeUser(0.1,0.5);
              temp->GetYaxis()->SetRangeUser(0.,30.);
              TH1D* temp1D = (TH1D*) temp->ProjectionY(Form("projEM02_%i_%i",iVec,iEM02),temp->GetXaxis()->FindBin(EM02_bins[iEM02]+0.001),temp->GetXaxis()->FindBin(EM02_bins[iEM02+1]-0.001));
              Double_t nEntries = temp1D->Integral(1,temp1D->GetNbinsX());
              temp1D->Sumw2();
              temp1D->Scale(1./nEntries);
              temp1D->GetXaxis()->SetRangeUser(0.1,0.5);
              vecClusterEM02_TH1D[iEM02].push_back(temp1D);
          }
          Int_t Xmin, Xmax = 0;
          GetMinMaxBin(vecClusterEM02_TH1D[iEM02],Xmin,Xmax);
          for(Int_t iVec=0; iVec<(Int_t)vecClusterEM02_TH1D[iEM02].size(); iVec++){
            ((TH1D*)vecClusterEM02_TH1D[iEM02].at(iVec))->GetXaxis()->SetRange(Xmin,Xmax);
          }
          DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                      vecClusterEM02_TH1D[iEM02],"",Form("#lambda_{0}^{2} for %.01f < E_{#mbox{cluster}} < %.01f GeV",EM02_bins[iEM02],EM02_bins[iEM02+1]),"#frac{1}{N} #frac{dN}{d#lambda_{0}^{2}}",1,1.1,
                                      labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                      0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
          SaveCanvas(canvas, Form("%s/Comparison/ClusterEM02_bin%i.%s", outputDir.Data(), iEM02, suffix.Data()),kFALSE,kTRUE);

          DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                          vecClusterEM02_TH1D[iEM02],"","#lambda_{0}^{2}","#frac{1}{N} #frac{dN}{d#lambda_{0}^{2}}",1,1.1,
                                          labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                          0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
          SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ClusterEM02_bin%i.%s", outputDir.Data(), iEM02, suffix.Data()));
        }
        DeleteVecTH2D(vecClusterEM02);
    }

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms for cluster energy" << endl;
    //---------------------------------------------------------------------------------------------------------------
    //------------------------------Cluster properties: energy ------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    for(Int_t iVec=0; iVec<(Int_t)vecClusterEnergy.size(); iVec++){
        TH1D* temp = vecClusterEnergy.at(iVec);
        temp->Scale(1./temp->Integral());
        //temp->GetXaxis()->SetRangeUser(minClusE,fBinsClusterPt[fNBinsClusterPt+1]);
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kTRUE,kTRUE,kFALSE,
                                vecClusterEnergy,"","Cluster Energy (GeV)","#frac{1}{N_{cl.}} #frac{dN}{dE}",1,1.1,
                                labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Energy_Cluster_afterQA.%s", outputDir.Data(), suffix.Data()), kTRUE, kTRUE);

    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterEnergy,"","Cluster Energy (GeV)","#frac{1}{N_{cl.}} #frac{dN}{dE}",1,1.1,
                                labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Energy_Cluster_afterQA_lin.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterEnergy,"","Cluster Energy (GeV)","#frac{1}{N_{cl.}} #frac{dN}{dE}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Energy_Cluster_afterQA.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecClusterEnergy);


    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms for cluster shapes"   << endl;
    //---------------------------------------------------------------------------------------------------------------
    //--------------------------------------Cluster properties: Cluster shape ---------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    for(Int_t iVec=0; iVec<(Int_t)vecClusterM02.size(); iVec++){
        TH1D* temp = vecClusterM02.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
    }
    //---------------------------------------------------------------------------------------------------------------
    // M02 (long axis)
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterM02,"","#lambda_{0}^{2}","#frac{1}{N_{cl.}} #frac{dN}{d#lambda_{0}^{2}}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/M02_afterQA.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterM02,"","#lambda_{0}^{2}","#frac{1}{N_{cl.}} #frac{dN}{d#lambda_{0}^{2}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_M02_afterQA.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecClusterM02);
    //---------------------------------------------------------------------------------------------------------------
    // M20 (short axis)
    for(Int_t iVec=0; iVec<(Int_t)vecClusterM20.size(); iVec++){
        TH1D* temp = vecClusterM20.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterM20,"","#lambda_{1}^{2}","#frac{1}{N_{cl.}} #frac{dN}{d#lambda_{1}^{2}}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/M20_afterQA.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterM20,"","#lambda_{1}^{2}","#frac{1}{N_{cl.}} #frac{dN}{d#lambda_{1}^{2}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_M20_afterQA.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecClusterM20);
    //---------------------------------------------------------------------------------------------------------------
    // dispersion
    for(Int_t iVec=0; iVec<(Int_t)vecClusterDispersion.size(); iVec++){
        TH1D* temp = vecClusterDispersion.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterDispersion,"","Dispersion","#frac{1}{N_{cl}} #frac{dN}{dDisp}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Dispersion_afterQA.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterDispersion,"","Dispersion","#frac{1}{N_{cl}} #frac{dN}{dDisp}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_Dispersion_afterQA.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecClusterDispersion);

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms for cluster NCells"  << endl;
    //---------------------------------------------------------------------------------------------------------------
    //----------------------------- Cluster properties: NCells ------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    for(Int_t iVec=0; iVec<(Int_t)vecClusterNCells.size(); iVec++){
        TH1D* temp = vecClusterNCells.at(iVec);
        temp->Sumw2();
        temp->Scale(1./temp->Integral());
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterNCells,"","N_{Cells} per Cluster","#frac{1}{N_{cl}} #frac{dN}{dN_{Cells}}",1,1.1,
                                labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/NCells_afterQA.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterNCells,"","N_{Cells} per Cluster","#frac{1}{N_{cl}} #frac{dN}{dN_{Cells}}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_NCells_afterQA.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecClusterNCells);

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Comparison Histograms for cluster NLMs" << endl;
    //---------------------------------------------------------------------------------------------------------------
    //------------------------------ Cluster properties: number of local maxima -------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    for(Int_t iVec=0; iVec<(Int_t)vecClusterNLM.size(); iVec++){
        TH1D* temp = vecClusterNLM.at(iVec);
        temp->Sumw2();
        Double_t nEntries = temp->Integral();
        temp->Scale(1./temp->Integral());
    }
    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                vecClusterNLM,"","Number of Local Maxima","#frac{1}{N_{cl}} #frac{dN}{dNLM}",1,1.1,
                                labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/NLM_Cluster_afterQA.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

    DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterNLM,"","Number of Local Maxima","#frac{1}{N_{cl}} #frac{dN}{dNLM}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 5, kFALSE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
    SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_NLM_Cluster_afterQA.%s", outputDir.Data(), suffix.Data()));
    DeleteVecTH1D(vecClusterNLM);


    //---------------------------------------------------------------------------------------------------------------
    //--------------------------------- Cluster properties: track matching ------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    if(isTrackMatching){

        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Drawing Comparison Histograms for cluster track matching" << endl;

        // distance of tracks to clusters
        for(Int_t iVec=0; iVec<(Int_t)vecClusterR.size(); iVec++){
            TH1D* temp = vecClusterR.at(iVec);
            temp->Sumw2();
            Double_t nEntries = temp->Integral();
            temp->Scale(1./nEntries);
        }
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecClusterR,"","R_{Cluster} (cm)","#frac{1}{N} #frac{dN}{dR}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/R_Cluster_afterQA.%s", outputDir.Data(), suffix.Data()), kFALSE , kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecClusterR,"","R_{Cluster} (cm)","#frac{1}{N} #frac{dN}{dR}",1,1.1,
                                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_R_Cluster_afterQA.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecClusterR);


        TGaxis::SetExponentOffset(-0.06, 0.01, "x");

        for(Int_t iProj=0; iProj<2; iProj++){
            for(Int_t iVec=0; iVec<(Int_t)vecMatchingDeltaEtaPhi_matched[iProj].size(); iVec++){
                TH1D* temp = vecMatchingDeltaEtaPhi_matched[iProj].at(iVec);
                temp->Sumw2();
                Double_t nEntries = temp->Integral();
                temp->Scale(1./nEntries);
            }
        }
        //---------------------------------------------------------------------------------------------------------------
        // distance of tracks and matched clusters
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecMatchingDeltaEtaPhi_matched[0],"","#Delta#eta_{Cluster - matched Tracks}","#frac{1}{N} #frac{dN}{d#Delta#eta}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 2, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/DeltaEta.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecMatchingDeltaEtaPhi_matched[1],"","#Delta#phi_{Cluster - matched Tracks}","#frac{1}{N} #frac{dN}{d#Delta#phi}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 2, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/DeltaPhi.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMatchingDeltaEtaPhi_matched[0],"","#Delta#eta_{Cluster - matched Tracks}","#frac{1}{N} #frac{dN}{d#Delta#eta}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 2, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_DeltaEta.%s", outputDir.Data(), suffix.Data()));
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMatchingDeltaEtaPhi_matched[1],"","#Delta#phi_{Cluster - matched Tracks}","#frac{1}{N} #frac{dN}{d#Delta#phi}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 2, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_DeltaPhi.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecMatchingDeltaEtaPhi_matched[0]);
        DeleteVecTH1D(vecMatchingDeltaEtaPhi_matched[1]);
        //---------------------------------------------------------------------------------------------------------------
        // distance of all tracks and clusters
        for(Int_t iProj=0; iProj<2; iProj++){
            for(Int_t iVec=0; iVec<(Int_t)vecMatchingDeltaEtaPhi_allTracks[iProj].size(); iVec++){
                TH1D* temp = vecMatchingDeltaEtaPhi_allTracks[iProj].at(iVec);
                temp->Sumw2();
                Double_t nEntries = temp->Integral();
                temp->Scale(1./nEntries);
            }
        }
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecMatchingDeltaEtaPhi_allTracks[0],"","#Delta#eta_{Cluster - all Tracks}","#frac{1}{N} #frac{dN}{d#Delta#eta}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/DeltaEta_allTracks.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecMatchingDeltaEtaPhi_allTracks[1],"","#Delta#phi_{Cluster - all Tracks}","#frac{1}{N} #frac{dN}{d#Delta#phi}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/DeltaPhi_allTracks.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMatchingDeltaEtaPhi_allTracks[0],"","#Delta#eta_{Cluster - all Tracks}","#frac{1}{N} #frac{dN}{d#Delta#eta}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_DeltaEta_allTracks.%s", outputDir.Data(), suffix.Data()));
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecMatchingDeltaEtaPhi_allTracks[1],"","#Delta#phi_{Cluster - all Tracks}","#frac{1}{N} #frac{dN}{d#Delta#phi}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_DeltaPhi_allTracks.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecMatchingDeltaEtaPhi_allTracks[0]);
        DeleteVecTH1D(vecMatchingDeltaEtaPhi_allTracks[1]);

        //---------------------------------------------------------------------------------------------------------------
        // further detailed TM cluster properties separated in
        if(doExtQA>0){
            for(Int_t iCharge=0; iCharge<2; iCharge++){
                // distance of all tracks and clusters for different charges
                for(Int_t iProj=0; iProj<2; iProj++){
                    for(Int_t iVec=0; iVec<(Int_t)vecMatchingDeltaEtaPhiTracksAll[iProj][iCharge].size(); iVec++){
                        TH1D* temp = vecMatchingDeltaEtaPhiTracksAll[iProj][iCharge].at(iVec);
                        temp->Sumw2();
                        Double_t nEntries = temp->Integral();
                        temp->Scale(1./nEntries);
                    }
                }
                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMatchingDeltaEtaPhiTracksAll[0][iCharge],"",Form("#Delta#eta_{Cluster - %s. Tracks}", charge[iCharge].Data()),"#frac{1}{N} #frac{dN}{d#Delta#eta}",1,1.1,
                                            labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                            0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/DeltaEta_%s.%s", outputDir.Data(), charge[iCharge].Data(), suffix.Data()), kFALSE, kTRUE);
                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMatchingDeltaEtaPhiTracksAll[1][iCharge],"",Form("#Delta#phi_{Cluster - %s. Tracks}", charge[iCharge].Data()),"#frac{1}{N} #frac{dN}{d#Delta#phi}",1,1.1,
                                            labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                            0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/DeltaPhi_%s.%s", outputDir.Data(), charge[iCharge].Data(), suffix.Data()), kFALSE, kTRUE);
                DeleteVecTH1D(vecMatchingDeltaEtaPhiTracksAll[0][iCharge]);
                DeleteVecTH1D(vecMatchingDeltaEtaPhiTracksAll[1][iCharge]);
                //---------------------------------------------------------------------------------------------------------------
                // distance of matched tracks and clusters for different charges
                for(Int_t iProj=0; iProj<2; iProj++){
                    for(Int_t iVec=0; iVec<(Int_t)vecMatchingDeltaEtaPhiTracksMatched[iProj][iCharge].size(); iVec++){
                        TH1D* temp = vecMatchingDeltaEtaPhiTracksMatched[iProj][iCharge].at(iVec);
                        temp->Sumw2();
                        Double_t nEntries = temp->Integral();
                        temp->Scale(1./nEntries);
                    }
                }

                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMatchingDeltaEtaPhiTracksMatched[0][iCharge],"",Form("#Delta#eta_{Cluster - matched %s. Tracks}", charge[iCharge].Data()),"#frac{1}{N} #frac{dN}{d#Delta#eta}",1,1.1,
                                            labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                            0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/DeltaEta_%s_matched.%s", outputDir.Data(), charge[iCharge].Data(), suffix.Data()), kFALSE, kTRUE);

                DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                            vecMatchingDeltaEtaPhiTracksMatched[1][iCharge],"",Form("#Delta#phi_{Cluster - matched %s. Tracks}", charge[iCharge].Data()),"#frac{1}{N} #frac{dN}{d#Delta#phi}",1,1.1,
                                            labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                            0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
                SaveCanvas(canvas, Form("%s/Comparison/DeltaPhi_%s_matched.%s", outputDir.Data(), charge[iCharge].Data(), suffix.Data()), kFALSE, kTRUE);
                DeleteVecTH1D(vecMatchingDeltaEtaPhiTracksMatched[0][iCharge]);
                DeleteVecTH1D(vecMatchingDeltaEtaPhiTracksMatched[1][iCharge]);
                //---------------------------------------------------------------------------------------------------------------
                // distance of matched trakcs for different charges and pt ranges
                for(Int_t iPt=0; iPt<3; iPt++){
                    std::vector<TH1D*> tempProjections[2];
                    for(Int_t iProj=0; iProj<2; iProj++){
                        tempProjections[iProj].clear();
                        for(Int_t iVec=0; iVec<(Int_t)vecMatchingDeltaEtaPhiTracksAllPt[iPt][iCharge].size(); iVec++){
                            TH2D* temp = (TH2D*) vecMatchingDeltaEtaPhiTracksAllPt[iPt][iCharge].at(iVec);

                            if(iProj==0){
                                tempProjections[iProj].push_back((TH1D*) temp->ProjectionX(Form("ProjEtaCharge%i-%s-%i", iVec, charge[iCharge].Data(), iPt),1,temp->GetNbinsY()));
                                SetXRange(tempProjections[iProj].at(iVec),1,temp->GetNbinsY());
                            } else if(iProj==1){
                                tempProjections[iProj].push_back((TH1D*) temp->ProjectionY(Form("ProjPhiCharge%i-%s-%i", iVec, charge[iCharge].Data(), iPt),1,temp->GetNbinsX()));
                                SetXRange(tempProjections[iProj].at(iVec),1,temp->GetNbinsX());
                            }
                            tempProjections[iProj].at(iVec)->Sumw2();
                            Double_t nEntries = tempProjections[iProj].at(iVec)->Integral();
                            tempProjections[iProj].at(iVec)->Scale(1/nEntries);
                        }
                    }

                    TString xTitle;
                    if(iPt==0)      xTitle = Form("_{Cluster - %s. Tracks with #it{p}_{T}<0.75 GeV/#it{c}}", charge[iCharge].Data());
                    else if(iPt==1) xTitle = Form("_{Cluster - %s. Tracks with 0.75<#it{p}_{T}<1.25 GeV/#it{c}}", charge[iCharge].Data());
                    else if(iPt==2) xTitle = Form("_{Cluster - %s. Tracks with #it{p}_{T}>1.25 GeV/#it{c}}", charge[iCharge].Data());

                    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                                tempProjections[0],"",Form("#Delta#eta%s",xTitle.Data()),"#frac{1}{N} #frac{dN}{d#Delta#eta}",1,1.1,
                                                labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
                    SaveCanvas(canvas, Form("%s/Comparison/DeltaEta_%s_%i.%s", outputDir.Data(), charge[iCharge].Data(), iPt, suffix.Data()), kFALSE, kTRUE);
                    DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                                tempProjections[1],"",Form("#Delta#phi%s",xTitle.Data()),"#frac{1}{N} #frac{dN}{d#Delta#phi}",1,1.1,
                                                labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                                0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
                    SaveCanvas(canvas, Form("%s/Comparison/DeltaPhi_%s_%i.%s", outputDir.Data(), charge[iCharge].Data(), iPt, suffix.Data()), kFALSE, kTRUE);
                    DeleteVecTH1D(tempProjections[0]);
                    DeleteVecTH1D(tempProjections[1]);
                    DeleteVecTH2D(vecMatchingDeltaEtaPhiTracksAllPt[iPt][iCharge]);
                }
            }
            TGaxis::SetExponentOffset(0, 0, "x");
        }
    }


    //---------------------------------------------------------------------------------------------------------------
    //---------------------------------- Cell properties: energy ----------------------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    if(doExtQA>1 && (Int_t)vecCellEnergyForComparison.size()>1){
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Drawing Comparison Histograms for cell properties"  << endl;

        for(Int_t iMC=1; iMC<nSets; iMC++){
        canvas->SetRightMargin(0.1);
        CheckCellsDataMC(cellQAData,
                        vecCellEnergyForComparison.at(0),
                        vecCellEnergyForComparison.at(iMC),
                        Form("Compare Cells: Data and %s",plotDataSets[iMC].Data()),
                        "CellID",
                        nCaloCells,
                        "Data",
                        plotDataSets[iMC].Data());
        PutProcessLabelAndEnergyOnPlot(0.7, 0.9, 0.03, fCollisionSystem.Data(), plotDataSets[iMC].Data(), fClusters.Data());
        SaveCanvas(canvas, Form("%s/Comparison/CompareCells_Data_%s.%s", outputDir.Data(), plotDataSets[iMC].Data(), suffix.Data()));
        canvas->SetRightMargin(rightMargin);
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //------------------------------------ Cell properties: frequency of usage and energy ---------------------------
    //---------------------------------------------------------------------------------------------------------------
    if(doExtQA>1){
        TGaxis::SetExponentOffset(-0.06, -0.04, "x");
        // fractional cluster energy in cells
        for(Int_t iVec=1; iVec<(Int_t)vecClusterEnergyFracCells.size(); iVec++){
            TH1D* temp = vecClusterEnergyFracCells.at(iVec);
            temp->Sumw2();
            temp->Add(vecClusterEnergyFracCells.at(0),-1);
            temp->GetXaxis()->SetRangeUser(0,nCaloCells);
            AdjustHistRange(temp,1.2,1.2,kTRUE);
            DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                temp,"","Cell ID",Form("Difference between %s and Data #frac{dN}{dCellID}",plotDataSets[iVec].Data()),1,1,
                                0.75,0.92,0.03,fCollisionSystem,fTrigger[0],"");
            SaveCanvas(canvas, Form("%s/Comparison/ClusterEnergyFracCells_%s.%s", outputDir.Data(), DataSets[iVec].Data(), suffix.Data()));
        }
        // frequency of cells included in cluster
        for(Int_t iVec=1; iVec<(Int_t)vecClusterIncludedCells.size(); iVec++){
            TH1D* temp = vecClusterIncludedCells.at(iVec);
            temp->Sumw2();
            temp->Add(vecClusterIncludedCells.at(0),-1);
            temp->GetXaxis()->SetRangeUser(0,nCaloCells);
            AdjustHistRange(temp,1.2,1.2,kTRUE);
            DrawPeriodQAHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kFALSE,kFALSE,
                                temp,"","Cell ID in accepted Clusters",Form("Difference between %s and Data #frac{d#it{N}}{d#it{Cell ID}}",plotDataSets[iVec].Data()),1,1,
                                0.75,0.92,0.03,fCollisionSystem,fTrigger[0],"");
            SaveCanvas(canvas, Form("%s/Comparison/ClusterIncludedCells_%s.%s", outputDir.Data(), DataSets[iVec].Data(), suffix.Data()));
        }
        DeleteVecTH1D(vecClusterEnergyFracCells);
        DeleteVecTH1D(vecClusterIncludedCells);
    }


    //---------------------------------------------------------------------------------------------------------------
    //--------------------- Comparison of conversion photon spatial distribution under meson peaks ------------------
    //---------------------------------------------------------------------------------------------------------------
    if(isPCMCalo){
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Drawing Comparison Histograms for conversion photons under meson peaks"<< endl;
        //---------------------------------------------------------------------------------------------------------------
        // Comparison of conversion photon spatial distribution under pi0 peak
        for(Int_t iProj=0; iProj<2; iProj++){
            for(Int_t iVec=0; iVec<(Int_t)vecConvPhotonEtaPhi_Pi0[iProj].size(); iVec++){
                TH1D* temp = vecConvPhotonEtaPhi_Pi0[iProj].at(iVec);
                temp->Sumw2();
                Double_t nEntries = temp->Integral();
                temp->Scale(1./nEntries);
            }
        }
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecConvPhotonEtaPhi_Pi0[1],"","#eta_{#gamma_{conv} under #pi^{0}-peak}","#frac{1}{N} #frac{dN}{d#eta}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 8, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ConvPhotonPi0_Eta.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecConvPhotonEtaPhi_Pi0[0],"","#phi_{#gamma_{conv} under #pi^{0}-peak}","#frac{1}{N} #frac{dN}{d#phi}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 16, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ConvPhotonPi0_Phi.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecConvPhotonEtaPhi_Pi0[1],"","#eta_{#gamma_{conv} under #pi^{0}-peak}","#frac{1}{N} #frac{dN}{d#eta}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 8, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ConvPhotonPi0_Eta.%s", outputDir.Data(), suffix.Data()));
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecConvPhotonEtaPhi_Pi0[0],"","#phi_{#gamma_{conv} under #pi^{0}-peak}","#frac{1}{N} #frac{dN}{d#phi}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 16, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ConvPhotonPi0_Phi.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecConvPhotonEtaPhi_Pi0[0]);
        DeleteVecTH1D(vecConvPhotonEtaPhi_Pi0[1]);
                //---------------------------------------------------------------------------------------------------------------
        // Comparison of conversion photon spatial distribution under pi0 peak

        for(Int_t iProj=0; iProj<2; iProj++){
            for(Int_t iVec=0; iVec<(Int_t)vecConvPhotonEtaPhi_Eta[iProj].size(); iVec++){
                TH1D* temp = vecConvPhotonEtaPhi_Eta[iProj].at(iVec);
                temp->Sumw2();
                Double_t nEntries = temp->Integral();
                temp->Scale(1./nEntries);
            }
        }
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecConvPhotonEtaPhi_Eta[1],"","#eta_{#gamma_{conv} under #eta-peak}","#frac{1}{N} #frac{dN}{d#eta}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ConvPhotonEta_Eta.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecConvPhotonEtaPhi_Eta[0],"","#phi_{#gamma_{conv} under #eta-peak}","#frac{1}{N} #frac{dN}{d#phi}",1,1.1,
                                    labelData, colorCompare, kTRUE, 2, 8, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/ConvPhotonEta_Phi.%s", outputDir.Data(), suffix.Data()), kFALSE, kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecConvPhotonEtaPhi_Eta[1],"","#eta_{#gamma_{conv} under #eta-peak}","#frac{1}{N} #frac{dN}{d#eta}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 4, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ConvPhotonEta_Eta.%s", outputDir.Data(), suffix.Data()));
        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecConvPhotonEtaPhi_Eta[0],"","#phi_{#gamma_{conv} under #eta-peak}","#frac{1}{N} #frac{dN}{d#phi}",1,1.1,
                                        labelData, colorCompare, kTRUE, 2, 8, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_ConvPhotonEta_Phi.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecConvPhotonEtaPhi_Eta[0]);
        DeleteVecTH1D(vecConvPhotonEtaPhi_Eta[1]);
    }


    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Drawing Special Histograms *******************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Drawing Special Histograms" << endl;

    const char* nameOutput = Form("%s/ClusterQA_%s.root",outputDirRootFiles.Data(),DataSets[0].Data());
    TFile* fOutput = new TFile(nameOutput,"UPDATE");

    for(Int_t iEM02 = 0; iEM02<nEM02_bins-1; iEM02++){
      WriteHistogramTH1DVec(fOutput,vecClusterEM02_TH1D[iEM02],Form("EM02_%i",iEM02));
    }

    //---------------------------------------------------------------------------------------------------------------
    //---------------------------- Cluster distribution on detector MC vs Data --------------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // calculate "asymetry of data and MC" for eta phi distribution
    for(Int_t iMC=1; iMC<nSets; iMC++){
        TH2D* temp = (TH2D*) vecEtaPhiAfterQA.at(iMC)->Clone(Form("CloneEtaVsPhi-%i",iMC));
        for(Int_t iX=1; iX<=temp->GetNbinsX(); iX++){
            for(Int_t iY=1; iY<=temp->GetNbinsY(); iY++){
                Double_t tempBinMC = temp->GetBinContent(iX,iY);
                Double_t tempBinData = vecEtaPhiAfterQA.at(0)->GetBinContent(iX,iY);
                if( (tempBinMC+tempBinData) > 0) temp->SetBinContent(iX,iY,(tempBinData-tempBinMC)/(tempBinData+tempBinMC));
                else temp->SetBinContent(iX,iY,-2);
            }
        }
        temp->GetZaxis()->SetRangeUser(-1,1);
        DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                            temp,Form("%s - %s %s& %s - (Data-MC) / (Data+MC) - %s",fCollisionSystem.Data(), plotDataSets[0].Data(), fTrigger[0].Data(), plotDataSets[iMC].Data(),
                                      fClusters.Data()),"#phi_{Cluster}","#eta_{Cluster}",1,1);
        SaveCanvasAndWriteHistogram(canvas, temp, Form("%s/EtaVsPhi_Asym%s_%s_afterClusterQA.%s", outputDir.Data(), DataSets[iMC].Data(), DataSets[0].Data(), suffix.Data()));
        delete temp;
    }
    //---------------------------------------------------------------------------------------------------------------
    // calculate ratio of data and MC
    for(Int_t iMC=1; iMC<nSets; iMC++){
        TH2D* temp = (TH2D*) vecEtaPhiAfterQA.at(iMC)->Clone(Form("CloneEtaVsPhi-%i",iMC));
        for(Int_t iX=1; iX<=temp->GetNbinsX(); iX++){
            for(Int_t iY=1; iY<=temp->GetNbinsY(); iY++){
                Double_t tempBinMC = temp->GetBinContent(iX,iY);
                Double_t tempBinData = vecEtaPhiAfterQA.at(0)->GetBinContent(iX,iY);
                if( tempBinMC != 0) temp->SetBinContent(iX,iY,tempBinData/tempBinMC);
                else temp->SetBinContent(iX,iY,-2);
            }
        }
        temp->GetZaxis()->SetRangeUser(0,2);
        DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                            temp,Form("%s - %s %s& %s - Data / MC - %s",fCollisionSystem.Data(), plotDataSets[0].Data(), fTrigger[0].Data(), plotDataSets[iMC].Data(),
                                      fClusters.Data()),"#phi_{Cluster}","#eta_{Cluster}",1,1);
        SaveCanvasAndWriteHistogram(canvas, temp, Form("%s/EtaVsPhi_Ratio%s_%s_afterClusterQA.%s", outputDir.Data(), DataSets[iMC].Data(), DataSets[0].Data(), suffix.Data()));
        delete temp;
    }

    //---------------------------------------------------------------------------------------------------------------
    //------------- Compare matching efficiencies in MC and data or for different data sets for V0 candidates -------
    //---------------------------------------------------------------------------------------------------------------
    if(isPCMCalo){
        for(Int_t i=0; i<nSets; i++){
            //---------------------------------------------------------------------------------------------------------------
            // Compare pure mass distribution and rejected candidates due to V0 matching
            cout << "Set: " << DataSets[i].Data() << endl;
            TH2D* mother = (TH2D*) vecESDMother.at(i);
            TH2D* motherMatched = (TH2D*) vecESDMotherMatched.at(i);

            if(mother->GetNbinsX() != motherMatched->GetNbinsX()) {cout << "ERROR: NBinsX of ESD_Mother and ESD_Mother_Matched do not agree, return..." << endl; return;}
            if(mother->GetNbinsY() != motherMatched->GetNbinsY()) {cout << "ERROR: NBinsY of ESD_Mother and ESD_Mother_Matched do not agree, return..." << endl; return;}

            if(mother && motherMatched){
                // two dimenstional comparison
                TH2D* fHistInvMassPtMotherMotherMatched = (TH2D*) mother->Clone("Mother_InvMass/(MotherMatched_InvMass + Mother_InvMass)");
                //TH2D* fHistInvMassPtMotherMotherMatched = new TH2D("Mother_InvMass/(MotherMatched_InvMass + Mother_InvMass)","Mother_InvMass/(MotherMatched_InvMass + Mother_InvMass)", 800, 0, 0.8, 250, 0, 25);
                fHistInvMassPtMotherMotherMatched->Reset("ICE");
                for(Int_t x = 1; x <= fHistInvMassPtMotherMotherMatched->GetNbinsX(); x++) {
                    for(Int_t y = 1; y <= fHistInvMassPtMotherMotherMatched->GetNbinsY(); y++) {
                        Double_t a = mother->GetBinContent(x,y);
                        Double_t b = motherMatched->GetBinContent(x,y);
                        Double_t r = 0;
                        if ( (a+b)!= 0 ) r = a/(a+b);
                        fHistInvMassPtMotherMotherMatched->SetBinContent(x,y,r);
                    }
                }

                DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                    fHistInvMassPtMotherMotherMatched,"","#it{M}_{inv, Mother} (GeV/#it{c}^{2})","#it{p}_{T, Mother} (GeV/#it{c})",1,0.8);

                TString whichTrMatch = fClusterCutSelection[i].Copy();
                whichTrMatch.Remove(10,9); whichTrMatch.Remove(0,9);

                TLatex* label[4];
                label[0] = new TLatex(0.5, 0.86, Form("%s: #frac{Accepted Mother-Candidates}{All Mother-Candidates}",plotDataSets[i].Data()));
                label[1] = new TLatex(0.5, 0.81, Form( "PCM, %s TrackMatchingCuts:",calo.Data()));
                if(whichTrMatch.Atoi() < 7){
                  label[2] = new TLatex(0.5, 0.77, Form( "- |#Delta#eta|<%s", TrackMatchingEtaMax[whichTrMatch.Atoi()].Data()));
                  label[3] = new TLatex(0.5, 0.73, Form( "- e^{+}: %s<#Delta#phi<%s; e^{-}: %.2f<#Delta#phi<%.2f", TrackMatchingPhiMin[whichTrMatch.Atoi()].Data(), TrackMatchingPhiMax[whichTrMatch.Atoi()].Data(),
                                                         (-1* TrackMatchingPhiMax[whichTrMatch.Atoi()].Atof()) , (-1* TrackMatchingPhiMin[whichTrMatch.Atoi()].Atof()) ));
                }else{
                  label[2] = new TLatex(0.5, 0.77, Form( "- #Delta#eta: %s", TrackMatchingEtaMax[whichTrMatch.Atoi()].Data()));
                  label[3] = new TLatex(0.5, 0.73, Form( "- #Delta#phi: %s", TrackMatchingPhiMin[whichTrMatch.Atoi()].Data()));
                }
                for(Int_t iL=0; iL<4; iL++){
                    SetStyleTLatex( label[iL], 0.03,4);
                    label[iL]->Draw();
                }

                SaveCanvasAndWriteHistogram(canvas, fHistInvMassPtMotherMotherMatched, Form("%s/%s_%s_cut%i.%s", outputDir.Data(),"TrackMatchingCut", DataSets[i].Data(), whichTrMatch.Atoi(), suffix.Data()));
                delete fHistInvMassPtMotherMotherMatched;

                //---------------------------------------------------------------------------------------------------------------
                // check if MC is available as well
                TH2D* mcPi0                 = 0x0;
                TH2D* mcEta                 = 0x0;
                TH2D* mcPi0Matched_MC       = 0x0;
                TH2D* mcEtaMatched_MC       = 0x0;
                TH2D* mcPi0Matched_True     = 0x0;
                TH2D* mcEtaMatched_True     = 0x0;

                Double_t nDataMatches       = 0;
                Double_t nMCMatches         = 0;
                Double_t nTrueMCMatches     = 0;

                Double_t nMCMatches_Missing = 0;
                Double_t nMCMatches_TooMany = 0;

                if(isTrueContainer[i]){
                    mcPi0                   = vecESDTrueMeson[0].at(i);
                    mcEta                   = vecESDTrueMeson[1].at(i);
                    mcPi0Matched_MC         = vecESDTrueMesonMatched[0].at(i);
                    mcEtaMatched_MC         = vecESDTrueMesonMatched[1].at(i);
                    mcPi0Matched_True       = vecESDTrueMesonCaloConvertedPhotonMatched[0].at(i);
                    mcEtaMatched_True       = vecESDTrueMesonCaloConvertedPhotonMatched[1].at(i);
                }

                //---------------------------------------------------------------------------------------------------------------
                //
                if( mcPi0Matched_True && mcEtaMatched_True){
                    TH2D* matched_True = (TH2D*) mcPi0Matched_True->Clone("matched_True");
                    matched_True->Add(mcEtaMatched_True,1);

                    TH2D* fHistInvMassPtMotherMatchedCompare_Matched_True = (TH2D*) mother->Clone("(MotherMatched_InvMass-TrueMotherMatched_InvMass)/(MotherMatched_InvMass+TrueMotherMatched_InvMass)");
                    fHistInvMassPtMotherMatchedCompare_Matched_True->Reset("ICE");
                    for(Int_t x = 1; x <= fHistInvMassPtMotherMatchedCompare_Matched_True->GetNbinsX(); x++) {
                        for(Int_t y = 1; y <= fHistInvMassPtMotherMatchedCompare_Matched_True->GetNbinsY(); y++) {
                            Double_t a = motherMatched->GetBinContent(x,y);
                            Double_t b = matched_True->GetBinContent(x,y);
                            nDataMatches += a;
                            Double_t r = 0;
                            if( a == 0 && b == 0) r = -2;
                            else r = (a-b)/(a+b);
                            fHistInvMassPtMotherMatchedCompare_Matched_True->SetBinContent(x,y,r);
                        }
                    }
                    fHistInvMassPtMotherMatchedCompare_Matched_True->GetZaxis()->SetRangeUser(-1,1);

                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                        fHistInvMassPtMotherMatchedCompare_Matched_True,"","#it{M}_{inv, Mother} (GeV/#it{c}^{2})","#it{p}_{T, Mother} (GeV/#it{c})",1,0.8);

                    TLatex* label1 = new TLatex(0.5, 0.88, Form("%s: #frac{V0^{full MC}_{rejected} - V0^{TrueMC}_{rejected}}{V0^{full MC}_{rejected} + V0^{TrueMC}_{rejected}}",plotDataSets[i].Data()));
                    SetStyleTLatex( label1, 0.03,4);
                    label1->Draw(); label[1]->Draw(); label[2]->Draw(); label[3]->Draw();
                    SaveCanvasAndWriteHistogram(canvas, fHistInvMassPtMotherMatchedCompare_Matched_True, Form("%s/%s_%s_cut%i.%s", outputDir.Data(),"TrackMatchingRatio_fullMC_TrueMC", DataSets[i].Data(),
                                                whichTrMatch.Atoi(), suffix.Data()));
                    delete fHistInvMassPtMotherMatchedCompare_Matched_True;
                    delete matched_True;
                    delete label1;
                }

                //---------------------------------------------------------------------------------------------------------------
                //
                if( mcPi0Matched_MC && mcEtaMatched_MC){
                    TH2D* matched_MC = (TH2D*) mcPi0Matched_MC->Clone("matched_MC");
                    matched_MC->Add(mcEtaMatched_MC,1);

                    TH2D* fHistInvMassPtMotherMatchedCompare_Matched_MC = (TH2D*) mother->Clone("(MotherMatched_InvMass-MCMatched_InvMass)/(MotherMatched_InvMass+MCMatched_InvMass)");
                    fHistInvMassPtMotherMatchedCompare_Matched_MC->Reset("ICE");
                    for(Int_t x = 1; x <= fHistInvMassPtMotherMatchedCompare_Matched_MC->GetNbinsX(); x++) {
                        for(Int_t y = 1; y <= fHistInvMassPtMotherMatchedCompare_Matched_MC->GetNbinsY(); y++) {
                            Double_t a = motherMatched->GetBinContent(x,y);
                            Double_t b = matched_MC->GetBinContent(x,y);
                            Double_t r = 0;
                            if( a == 0 && b == 0) r = -2;
                            else r = (a-b)/(a+b);
                            fHistInvMassPtMotherMatchedCompare_Matched_MC->SetBinContent(x,y,r);
                        }
                    }
                    fHistInvMassPtMotherMatchedCompare_Matched_MC->GetZaxis()->SetRangeUser(-1,1);

                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                        fHistInvMassPtMotherMatchedCompare_Matched_MC,"","#it{M}_{inv, Mother} (GeV/#it{c}^{2})","#it{p}_{T, Mother} (GeV/#it{c})",1,0.8);

                    TLatex* label2 = new TLatex(0.5, 0.88, Form("%s: #frac{V0^{full MC}_{rejected} - V0^{MC}_{rejected}}{V0^{full MC}_{rejected} + V0^{MC}_{rejected}}",plotDataSets[i].Data()));
                    label2->SetNDC();
                    label2->SetTextColor(1);
                    label2->SetTextSize(0.03);
                    label2->Draw(); label[1]->Draw(); label[2]->Draw(); label[3]->Draw();
                    SaveCanvasAndWriteHistogram(canvas, fHistInvMassPtMotherMatchedCompare_Matched_MC, Form("%s/%s_%s_cut%i.%s", outputDir.Data(),"TrackMatchingRatio_fullMC_MC", DataSets[i].Data(),
                                                whichTrMatch.Atoi(), suffix.Data()));
                    delete fHistInvMassPtMotherMatchedCompare_Matched_MC;
                    delete matched_MC;
                    delete label2;
                }

                //---------------------------------------------------------------------------------------------------------------
                if( mcPi0Matched_True && mcEtaMatched_True && mcPi0Matched_MC && mcEtaMatched_MC){
                    TH2D* matched_True = (TH2D*) mcPi0Matched_True->Clone("matched_True2");
                    matched_True->Add(mcEtaMatched_True,1);
                    TH2D* matched_MC = (TH2D*) mcPi0Matched_MC->Clone("matched_MC2");
                    matched_MC->Add(mcEtaMatched_MC,1);

                    TH2D* fHistInvMassPtMotherMatchedCompare_True_MC = (TH2D*) mother->Clone("(TrueMotherMatched_InvMass-MCMatched_InvMass)/(TrueMotherMatched_InvMass+MCMatched_InvMass)");
                    fHistInvMassPtMotherMatchedCompare_True_MC->Reset("ICE");
                    for(Int_t x = 1; x <= fHistInvMassPtMotherMatchedCompare_True_MC->GetNbinsX(); x++) {
                        for(Int_t y = 1; y <= fHistInvMassPtMotherMatchedCompare_True_MC->GetNbinsY(); y++) {
                            Double_t a = matched_True->GetBinContent(x,y);
                            Double_t b = matched_MC->GetBinContent(x,y);
                            nTrueMCMatches += a;
                            nMCMatches += b;
                            Double_t r = 0;
                            if( a == 0 && b == 0) r = -2;
                            else r = (a-b)/(a+b);
                            fHistInvMassPtMotherMatchedCompare_True_MC->SetBinContent(x,y,r);
                            if( a > b ) nMCMatches_Missing += abs(a-b);
                            if( a < b ) nMCMatches_TooMany += abs(a-b);
                        }
                    }
                    fHistInvMassPtMotherMatchedCompare_True_MC->GetZaxis()->SetRangeUser(-1,1);

                    DrawPeriodQAHistoTH2(canvas,leftMargin,0.1,topMargin,bottomMargin,kFALSE,kFALSE,kFALSE,
                                        fHistInvMassPtMotherMatchedCompare_True_MC,"","#it{M}_{inv, Mother} (GeV/#it{c}^{2})","#it{p}_{T, Mother} (GeV/#it{c})",1,0.8);

                    TLatex* label3 = new TLatex(0.5, 0.88, Form("%s: #frac{V0^{TrueMC}_{rejected} - V0^{MC}_{rejected}}{V0^{TrueMC}_{rejected} + V0^{MC}_{rejected}}",plotDataSets[i].Data()));
                    SetStyleTLatex( label3, 0.03,4);
                    label3->Draw(); label[1]->Draw(); label[2]->Draw(); label[3]->Draw();
                    SaveCanvasAndWriteHistogram(canvas, fHistInvMassPtMotherMatchedCompare_True_MC, Form("%s/%s_%s_cut%i.%s", outputDir.Data(),"TrackMatchingRatio_TrueMC_MC", DataSets[i].Data(),
                                                whichTrMatch.Atoi(), suffix.Data()));
                    delete fHistInvMassPtMotherMatchedCompare_True_MC;
                    delete matched_True;
                    delete matched_MC;
                    delete label3;
                }
                delete label[0]; delete label[1]; delete label[2]; delete label[3];

                //---------------------------------------------------------------------------------------------------------------
                // Calculate how many we accidently rejected and how many we didn't get
                if( mcPi0Matched_True && mcEtaMatched_True && mcPi0Matched_MC && mcEtaMatched_MC ){
                    TH2D* mcCompare = (TH2D*) mcPi0->Clone("mcPi0Eta_TrueTH1D");
                    mcCompare->Add(mcEta,1);
                    TH2D* matched_True = (TH2D*) mcPi0Matched_True->Clone("matched_TrueTH1D");
                    matched_True->Add(mcEtaMatched_True,1);
                    TH2D* matched_MC = (TH2D*) mcPi0Matched_MC->Clone("matched_MCTH1D");
                    matched_MC->Add(mcEtaMatched_MC,1);

                    Int_t rangeMin[2]={103,488};
                    Int_t rangeMax[2]={157,598};

                    for(Int_t iMeson=0; iMeson<2; iMeson++){
                        TH1D* mcFractionMatchesMissingPt = new TH1D(Form("fractionMissing%i",iMeson),"",fNBinsAnalysisPt,fBinsAnalysisPt);
                        TH1D* mcFractionMatchesTooManyPt = new TH1D(Form("fractionTooMany%i",iMeson),"",fNBinsAnalysisPt,fBinsAnalysisPt);

                        for(Int_t iPt = 0; iPt<fNBinsAnalysisPt; iPt++){
                            Double_t bMatchesMissing = 0;
                            Double_t cMatchesTooMany = 0;
                            Double_t tMatches = 0;
                            for(Int_t x = rangeMin[iMeson]; x <= rangeMax[iMeson]; x++) {
                                for(Int_t y = (fBinsAnalysisPt[iPt]*10)+1; y <= (fBinsAnalysisPt[iPt+1]*10); y++) {
                                    tMatches += mcCompare->GetBinContent(x,y);
                                    Double_t a = matched_True->GetBinContent(x,y);
                                    Double_t b = matched_MC->GetBinContent(x,y);
                                    if( a > b ) bMatchesMissing += abs(a-b);
                                    if( a < b) cMatchesTooMany += abs(a-b);
                                }
                            }
                            Double_t rMiss, rTooMany;
                            if(tMatches > 0 ) {
                                rMiss = (bMatchesMissing / tMatches) * 100;
                                rTooMany = (cMatchesTooMany / tMatches) * 100;
                            } else {
                                rMiss = 0;
                                rTooMany = 0;
                            }
                            //cout << iMeson << "" << bMatchesMissing << "" << cMatchesTooMany << "" << tMatches << endl;
                            mcFractionMatchesMissingPt->SetBinContent(iPt+1,rMiss);
                            mcFractionMatchesTooManyPt->SetBinContent(iPt+1,rTooMany);
                            //cout << fBinsAnalysisPt[iPt] << " - " << fBinsAnalysisPt[iPt+1] << ", rMiss:" << rMiss << ", rTooMany: " << rTooMany << ", rTotal: " << rTotal
                            //     << ", rTooMany*rTotal*0.01: " << rTooMany*rTotal*0.01 <<endl;
                        }

                        EditTH1NoRunwise(mcFractionMatchesTooManyPt,24,1,kAzure,kAzure,1.2,0.9);
                        mcFractionMatchesTooManyPt->SetLineWidth(1);
                        mcFractionMatchesTooManyPt->GetXaxis()->SetRangeUser(fBinsAnalysisPt[0],fBinsAnalysisPt[fNBinsAnalysisPt]+2);
                        mcFractionMatchesTooManyPt->GetXaxis()->SetTitle("#it{p}_{T, #gamma#gamma} (GeV/#it{c})");
                        mcFractionMatchesTooManyPt->GetYaxis()->SetRangeUser(1E-3,120);
                        mcFractionMatchesTooManyPt->GetYaxis()->SetTitle("fraction of all true candidates (%)");
                        mcFractionMatchesTooManyPt->Sumw2();
                        mcFractionMatchesTooManyPt->Draw("p");

                        EditTH1NoRunwise(mcFractionMatchesMissingPt,20,1,kRed,kRed,1.2,0.9);
                        mcFractionMatchesMissingPt->SetLineWidth(1);
                        mcFractionMatchesMissingPt->Sumw2();
                        mcFractionMatchesMissingPt->Draw("p, same");

                        TLegend *legend = GetAndSetLegend2(0.35, 0.82, 0.75, 0.9, 0.03, 1,"",42);
                        legend->AddEntry(mcFractionMatchesMissingPt,"... misses true matching");
                        legend->AddEntry(mcFractionMatchesTooManyPt,"... removes true candidate by mistake");
                        legend->Draw("same");

                        TLatex* clusVZero = new TLatex(0.45, 0.91, "Cluster - V^{0}-track matching...");
                        SetStyleTLatex( clusVZero, 0.03,4);
                        clusVZero->Draw("SAME");

                        PutProcessLabelAndEnergyOnPlot(0.15, 0.95, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTextMeasurementMeson[iMeson].Data());
                        PutProcessLabelAndEnergyOnPlot(0.15, 0.84, 0.03, Form("%i < #it{M}_{#gamma#gamma} < %i MeV/#it{c}^{2}",rangeMin[iMeson],rangeMax[iMeson]), Form("#gamma's rec. with PCM, %s",
                                                       calo.Data()), "");

                        canvas->SetRightMargin(0.02); canvas->SetLeftMargin(0.08); canvas->SetBottomMargin(0.12); canvas->SetTopMargin(0.02);
                        SaveWriteCanvas(canvas, Form("%s/%s_%s_cut%i_%i.%s", outputDir.Data(),"TrackMatchingMissing_TrueMC_MC_Pt", DataSets[i].Data(), whichTrMatch.Atoi(), iMeson, suffix.Data()), kTRUE, kTRUE,
                                        kFALSE);
                        canvas->SetRightMargin(rightMargin); canvas->SetLeftMargin(leftMargin); canvas->SetBottomMargin(bottomMargin); canvas->SetTopMargin(topMargin);

                        delete mcFractionMatchesMissingPt;
                        delete mcFractionMatchesTooManyPt;
                        delete legend;
                        delete clusVZero;
                    }

                    delete mcCompare;
                    delete matched_MC;
                    delete matched_True;
                }

                //---------------------------------------------------------------------------------------------------------------
                //
                if( mcPi0 && mcPi0Matched_True && mcPi0Matched_MC){
                    TLegend *legend     = GetAndSetLegend2(0.3, 0.8, 0.7, 0.9, 0.03, 1,"",42);
                    TH2D* mcCompare     = (TH2D*) mcPi0->Clone("mcCompare");
                    for(Int_t x = 1; x <= mcCompare->GetNbinsX(); x++) {
                        for(Int_t y = 1; y <= mcCompare->GetNbinsY(); y++) {
                            Double_t a = mcCompare->GetBinContent(x,y);
                            Double_t b = mcPi0Matched_MC->GetBinContent(x,y);
                            mcCompare->SetBinContent(x,y,(a+b));
                        }
                    }

                    TH1D* projectTrue = (TH1D*) mcCompare->ProjectionX("projecttrue",1,mcCompare->GetNbinsY());
                    legend->AddEntry(projectTrue,"True MC");
                    projectTrue->Rebin(5);
                    projectTrue->SetLineColor(1);
                    DrawAutoGammaHisto(projectTrue,
                                        "",
                                        "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma#gamma}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,0.7);

                    TH1D* projectCompare = (TH1D*) mcPi0Matched_True->ProjectionX("projectcompare",1,mcPi0Matched_True->GetNbinsY());
                    legend->AddEntry(projectCompare,"Cluster - V^{0}-Track Matching");
                    projectCompare->Rebin(5);
                    projectCompare->SetLineColor(2);
                    projectCompare->SetFillColor(2);
                    projectCompare->SetFillStyle(3003);
                    projectCompare->DrawCopy("same,e,hist");

                    legend->Draw("same");

                    PutProcessLabelAndEnergyOnPlot(0.94, 0.9, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), Form("%s",fTextMeasurement.Data()), 42, 0.03, "", 1, 1.25, 31);
                    PutProcessLabelAndEnergyOnPlot(0.94, 0.79, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "", "", 42, 0.03, "", 1, 1.25, 31);
                    SaveWriteCanvas(canvas, Form("%s/TrackMatching_invMassPi0_withTrueMatches_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE , kFALSE);
                    delete projectTrue;
                    delete projectCompare;

                    mcCompare->Add(mcPi0Matched_True,-1);
                    TH1D* projectCompareW = (TH1D*) mcCompare->ProjectionX("projectcomparew",1,mcCompare->GetNbinsY());
                    projectCompareW->Rebin(5);
                    projectCompareW->SetLineColor(1);
                    DrawAutoGammaHisto(projectCompareW,
                                        "",
                                        "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma#gamma}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,0.7);
                    PutProcessLabelAndEnergyOnPlot(0.94, 0.9, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), Form("%s",fTextMeasurement.Data()), 42, 0.03, "", 1, 1.25, 31);
                    PutProcessLabelAndEnergyOnPlot(0.94, 0.79, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "", "", 42, 0.03, "", 1, 1.25, 31);
                    SaveWriteCanvas(canvas, Form("%s/TrackMatching_invMassPi0_withoutTrueMatches_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()), kFALSE, kTRUE , kFALSE);
                    delete projectCompareW;
                    delete legend;
                    delete mcCompare;
                }

                //---------------------------------------------------------------------------------------------------------------
                //
                if( mcPi0 && mcEta && mcPi0Matched_True && mcEtaMatched_True && mcPi0Matched_MC && mcEtaMatched_MC) {
                    TCanvas* cvs2pads   = new TCanvas("cvs2pads","",1350,1500);  // gives the page size
                    TPad* pad1          = new TPad("pad1", "", 0., 0.54, 1., 1.,-1, -1, -2);
                    DrawGammaPadSettings( pad1, 0.08, 0.02, 0.02, 0);
                    pad1->SetLogy(1);
                    pad1->Draw();
                    TPad* pad2          = new TPad("pad2", "", 0., 0., 1., 0.54,-1, -1, -2);
                    DrawGammaPadSettings( pad2, 0.08, 0.02, 0.00, 0.14);
                    pad2->SetLogy(1);
                    pad2->Draw();

                    const Int_t nBin    = 3;
                    Int_t lowBin[nBin]  = {0,16,40};
                    Int_t upBin[nBin]   = {300,18,50};
                    Int_t lastBinning[2];
                    lastBinning[0]      = 0;
                    lastBinning[1]      = 0;
                    for(Int_t iB=0; iB<nBin; iB++){
                        TLegend *legend = GetAndSetLegend2(0.3, 0.8, 0.7, 0.88, 0.03, 1,"",42);

                        TH2D* mcCompare = (TH2D*) mcPi0->Clone(Form("mcComparePi0Eta_Bins%i-%i",lowBin[iB],upBin[iB]));
                        mcCompare->Add(mcEta,1);
                        for(Int_t x = 1; x <= mcCompare->GetNbinsX(); x++) {
                            for(Int_t y = 1; y <= mcCompare->GetNbinsY(); y++) {
                                Double_t a = mcCompare->GetBinContent(x,y);
                                Double_t b = mcPi0Matched_MC->GetBinContent(x,y);
                                Double_t c = mcEtaMatched_MC->GetBinContent(x,y);
                                mcCompare->SetBinContent(x,y,(a+b+c));
                            }
                        }

                        TH1D* projectTrue = (TH1D*) mcCompare->ProjectionX(Form("projecttruePi0Eta_Bins%i-%i",lowBin[iB],upBin[iB]),lowBin[iB],upBin[iB]);
                        legend->AddEntry(projectTrue,"true validated Monte Carlo");
                        projectTrue->Rebin(10);
                        projectTrue->SetLineColor(1);
                        canvas->cd();
                        DrawAutoGammaHisto(projectTrue,
                                        "",
                                        "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma#gamma}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,0.7);

                        canvas->SetRightMargin(rightMargin);
                        TH1D* projectCompare = (TH1D*) mcPi0Matched_True->ProjectionX(Form("projectcomparePi0_Bins%i-%i",lowBin[iB],upBin[iB]),lowBin[iB],upBin[iB]);
                        TH1D* projectCompareEta = (TH1D*) mcEtaMatched_True->ProjectionX(Form("projectcompareEta_Bins%i-%i",lowBin[iB],upBin[iB]),lowBin[iB],upBin[iB]);
                        projectCompare->Add(projectCompareEta,1);
                        legend->AddEntry(projectCompare,"cluster - V^{0}-track matching");
                        projectCompare->Rebin(10);
                        projectCompare->SetLineColor(2);
                        projectCompare->SetFillColor(2);
                        projectCompare->SetFillStyle(3003);
                        projectCompare->DrawCopy("same,e,hist");
                        legend->Draw("same");

                        PutProcessLabelAndEnergyOnPlot(0.94, 0.9, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), Form("%s & %s",fTextMeasurement.Data(),fTextMeasurementEta.Data()),
                                                       42, 0.03, "", 1, 1.25, 31);
                        PutProcessLabelAndEnergyOnPlot(0.94, 0.785, 0.03, Form("%.1f < #it{p}_{T, #gamma#gamma} < %.1f GeV/#it{c}",mcCompare->GetYaxis()->GetBinUpEdge(lowBin[iB]),
                                                                               mcCompare->GetYaxis()->GetBinUpEdge(upBin[iB])),"","", 42, 0.03, "", 1, 1.25, 31);
                        PutProcessLabelAndEnergyOnPlot(0.94, 0.74, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "","", 42, 0.03, "", 1, 1.25, 31);
                        SaveWriteCanvas(canvas, Form("%s/TrackMatching_invMassPi0Eta_withTrueMatches_%s_Bins%i-%i.%s", outputDir.Data(), DataSets[i].Data(), lowBin[iB],upBin[iB], suffix.Data()), kFALSE,
                                        kTRUE , kFALSE);

                        if(iB == nBin-1){
                            pad1->cd();
                            AdjustHistRange(projectTrue,1,10,kTRUE,1,0.7);
                            if(mode == 2 && fCollisionSystem.Contains("8 TeV")) projectTrue->GetYaxis()->SetRangeUser(0.7,7E4);
                            //projectTrue->GetYaxis()->SetRangeUser(0.7,9.5E3);
                            //projectTrue->GetYaxis()->SetTitleSize(0.09);
                            DrawAutoGammaHistoPaper(projectTrue,
                                                    "",
                                                    "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                                    "d#it{N}/d#it{M}_{#gamma#gamma}",
                                                    0,0,0,
                                                    0,0,0,
                                                    1,0,0.8,0.5,0.07);
                            projectCompare->DrawCopy("same,e,hist");
                            PutProcessLabelAndEnergyOnPlot(0.94, 0.98, 0.063, "ALICE simulation", "","", 42, 0.03, "", 1, 1.25, 31);
                            PutProcessLabelAndEnergyOnPlot(0.94, 0.9, 0.063, fCollisionSystem.Data(), plotDataSets[i].Data(), Form("%s & %s",fTextMeasurement.Data(),fTextMeasurementEta.Data()), 42, 0.03, "", 1, 1.25, 31);
                            PutProcessLabelAndEnergyOnPlot(0.1, 0.98, 0.063, Form("%.1f < #it{p}_{T, #gamma#gamma} < %.1f GeV/#it{c}",mcCompare->GetYaxis()->GetBinUpEdge(lowBin[iB]),
                                                                                  mcCompare->GetYaxis()->GetBinUpEdge(upBin[iB])),"","", 42, 0.03, "", 1, 1.25, 11);
                            PutProcessLabelAndEnergyOnPlot(0.94, 0.66, 0.063, Form("#gamma's rec. with PCM, %s",calo.Data()), "" ,"", 42, 0.03, "", 1, 1.25, 31);
                        } else if(iB == nBin-2){
                            pad2->cd();
                            AdjustHistRange(projectTrue,1,10,kTRUE,1,0.7);
                            if(mode == 2 && fCollisionSystem.Contains("8 TeV")) projectTrue->GetYaxis()->SetRangeUser(0.7,7E4);
                            //projectTrue->GetYaxis()->SetRangeUser(0.7,9.5E3);
                            DrawAutoGammaHistoPaper(projectTrue,
                                                    "",
                                                    "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                                    "d#it{N}/d#it{M}_{#gamma#gamma}",
                                                    0,0,0,
                                                    0,0,0,
                                                    1,0,0.8,0.6);
                            projectCompare->DrawCopy("same,e,hist");
                            TLegend* leg2 = (TLegend*) legend->Clone();
                            leg2->SetX1(0.46);
                            leg2->SetX2(0.86);
                            leg2->SetY1(0.847);
                            leg2->SetY2(0.967);
                            leg2->SetTextSize(0.054);
                            leg2->Draw("same");
                            PutProcessLabelAndEnergyOnPlot(0.1, 1, 0.054, Form("%.1f < #it{p}_{T, #gamma#gamma} < %.1f GeV/#it{c}",mcCompare->GetYaxis()->GetBinUpEdge(lowBin[iB]),
                                                                               mcCompare->GetYaxis()->GetBinUpEdge(upBin[iB])),"","", 42, 0.03, "", 1, 1.25, 11);
                        }

                        canvas->cd();
                        mcCompare->Add(mcPi0Matched_True,-1);
                        mcCompare->Add(mcEtaMatched_True,-1);
                        TH1D* projectCompareW = (TH1D*) mcCompare->ProjectionX(Form("projectcomparew_Bins%i-%i",lowBin[iB],upBin[iB]),lowBin[iB],upBin[iB]);
                        projectCompareW->Rebin(5);
                        projectCompareW->SetLineColor(1);
                        DrawAutoGammaHisto(projectCompareW,
                                        "",
                                        "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma#gamma}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8);
                        PutProcessLabelAndEnergyOnPlot(0.94, 0.9, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), Form("%s & %s",fTextMeasurement.Data(),fTextMeasurementEta.Data()), 42, 0.03, "", 1, 1.25, 31);
                        PutProcessLabelAndEnergyOnPlot(0.94, 0.785, 0.03, Form("%.1f < #it{p}_{T, #gamma#gamma} < %.1f GeV/#it{c}",mcCompare->GetYaxis()->GetBinUpEdge(lowBin[iB]),
                                                                               mcCompare->GetYaxis()->GetBinUpEdge(upBin[iB])),"","", 42, 0.03, "", 1, 1.25, 31);
                        PutProcessLabelAndEnergyOnPlot(0.94, 0.74, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "","", 42, 0.03, "", 1, 1.25, 31);
                        SaveWriteCanvas(canvas, Form("%s/TrackMatching_invMassPi0Eta_withoutTrueMatches_%s_Bins%i-%i.%s", outputDir.Data(), DataSets[i].Data(),lowBin[iB],upBin[iB], suffix.Data()), kFALSE,
                                        kTRUE , kFALSE);

                        delete projectTrue;
                        delete projectCompare;
                        delete projectCompareEta;
                        delete projectCompareW;
                        delete mcCompare;
                        delete legend;
                    }
                    SaveCanvas(cvs2pads, Form("%s/TrackMatching_invMassPi0Eta_withTrueMatches_%s_Bins%i-%i_%i-%i.%s", outputDir.Data(), DataSets[i].Data(),lowBin[nBin-2],upBin[nBin-2],
                               lowBin[nBin-1],upBin[nBin-1], suffix.Data()), kFALSE, kFALSE , kFALSE);
                    delete cvs2pads;
                }

                if(isTrueContainer[i]){
                    if(mcPi0Matched_MC && mcEtaMatched_MC && mcPi0Matched_True && mcEtaMatched_True){
                        if(nTrueMCMatches>0 && nMCMatches>0){
                        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        cout << "Matching efficiencies: "<< endl;
                        cout << "nDataMatches: " << nDataMatches << endl;
                        cout << "nMCMatches: " << nMCMatches << endl;
                        cout << "nTrueMCMatches: " << nTrueMCMatches << endl;
                        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        cout << "nData/nMC: " << nDataMatches/nMCMatches*100 << "%" << endl;
                        cout << "nData/nTrueMC: " << nDataMatches/nTrueMCMatches*100 << "%" << endl;
                        cout << "nMC/nTrueMC: " << (nMCMatches-nMCMatches_TooMany)/nTrueMCMatches*100  << "%, missing matches in MC: " << nMCMatches_Missing
                            << "(" << nMCMatches_Missing/nTrueMCMatches*100 << "%) - too many matches in MC(already subtracted): " << nMCMatches_TooMany
                            << "(" << nMCMatches_TooMany/nTrueMCMatches*100 << "%)" << endl;
                        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        fLog << "Matching efficiencies: "<< endl;
                        fLog << "nDataMatches: " << nDataMatches << endl;
                        fLog << "nMCMatches: " << nMCMatches << endl;
                        fLog << "nTrueMCMatches: " << nTrueMCMatches << endl;
                        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        fLog << "nData/nMC: " << nDataMatches/nMCMatches*100 << "%" << endl;
                        fLog << "nData/nTrueMC: " << nDataMatches/nTrueMCMatches*100 << "%" << endl;
                        fLog << "nMC/nTrueMC: " << (nMCMatches-nMCMatches_TooMany)/nTrueMCMatches*100  << "%, missing matches in MC: " << nMCMatches_Missing
                            << "(" << nMCMatches_Missing/nTrueMCMatches*100 << "%) - too many matches in MC(already subtracted): " << nMCMatches_TooMany
                            << "(" << nMCMatches_TooMany/nTrueMCMatches*100 << "%)" << endl;
                        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
                        } else{cout << "ERROR: nTrueMCMatches == 0 (" << nTrueMCMatches << ") or nMCMatches == 0 (" << nMCMatches << "), returning" << endl; return;}
                    } else {cout << "INFO: Could not find 'ESD_TruePi0_Matched_InvMass_Pt' and/or 'ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt'" << endl;}
                }
            } else {cout << "WARNING: ESD_Mother_InvMass_Pt or ESD_MotherMatched_InvMass_Pt not found for " << DataSets[i].Data() << " ! Skipping drawing of related histograms..." << endl;}
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //----------------------------- Further detailed QA histos for ConvCalo concerning matching----------------------
    //---------------------------------------------------------------------------------------------------------------
    canvas->cd();
    DrawGammaCanvasSettings( canvas, leftMargin, rightMargin, topMargin, bottomMargin);

    if(isPCMCalo){
        Int_t nESDRebin = 5;
        Int_t iMinPt[5]={0,6,10,20,50};
        std::vector<TH1D*> vecRange;
        std::vector<TH1D*> vecESD_InvMass_Matched_Sum;
        for(Int_t i=0; i<nSets; i++){
            for(Int_t iMin=0; iMin<5; iMin++){
                //---------------------------------------------------------------------------------------------------------------
                //
                TH2D* ESD_InvMass_Pt = (TH2D*) vecESDMother.at(i)->Clone();
                TH2D* ESD_Matched_InvMass_Pt = (TH2D*) vecESDMotherMatched.at(i)->Clone();

                TH1D* ESD_InvMass = (TH1D*) ESD_InvMass_Pt->ProjectionX("ESD_InvMass",iMinPt[iMin],ESD_InvMass_Pt->GetNbinsY());
                TH1D* ESD_Matched_InvMass = (TH1D*) ESD_Matched_InvMass_Pt->ProjectionX("ESD_InvMass_Matched",iMinPt[iMin],ESD_Matched_InvMass_Pt->GetNbinsY());
                TH1D* ESD_InvMass_Matched_Sum = (TH1D*) ESD_InvMass->Clone();
                ESD_InvMass->Sumw2();
                vecRange.push_back(ESD_InvMass);
                ESD_Matched_InvMass->Sumw2();
                vecRange.push_back(ESD_Matched_InvMass);
                ESD_InvMass_Matched_Sum->Sumw2();
                vecRange.push_back(ESD_InvMass_Matched_Sum);
                ESD_InvMass_Matched_Sum->Add(ESD_Matched_InvMass,1);

                ESD_InvMass->Rebin(nESDRebin);
                ESD_Matched_InvMass->Rebin(nESDRebin);
                ESD_InvMass_Matched_Sum->Rebin(nESDRebin);
                AdjustHistRange(vecRange,1.,10.,kTRUE,1,1.);
                DrawAutoGammaHistMatch3H(ESD_InvMass_Matched_Sum, ESD_InvMass, ESD_Matched_InvMass,
                                        "",
                                        "#it{M}_{inv, Mother} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{inv}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,
                                        1,1.2,"Sum", "Accepted", "Rejected", kFALSE);
                PutProcessLabelAndEnergyOnPlot(0.13, 0.93, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(),
                                               Form("Excluded M_{inv}^{Mother}-candidates by V^{0}-track matching, #it{p}^{Mother}_{T} > %.01f GeV/#it{c}",((Float_t)iMinPt[iMin])/10));
                SaveCanvas(canvas, Form("%s/TrackMatching_ESD_IntegratedPt_fromPtBin_%i_%s.%s", outputDir.Data(), iMinPt[iMin], DataSets[i].Data(),suffix.Data()), kFALSE, kTRUE);

                delete ESD_InvMass_Pt;
                delete ESD_Matched_InvMass_Pt;
                delete ESD_InvMass;
                if(iMin == 0) vecESD_InvMass_Matched_Sum.push_back(new TH1D(*ESD_Matched_InvMass));
                delete ESD_Matched_InvMass;
                delete ESD_InvMass_Matched_Sum;

                vecRange.clear();
            }
        }

        for(Int_t iVec=0; iVec<(Int_t)vecESD_InvMass_Matched_Sum.size(); iVec++){
            TH1D* temp = vecESD_InvMass_Matched_Sum.at(iVec);
            temp->Sumw2();
            Double_t nEntries = temp->Integral(); //temp->FindBin(0.05),temp->FindBin(0.1)); // normalization on left side of peak
            temp->Scale(1./nEntries);
        }
        DrawPeriodQACompareHistoTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                    vecESD_InvMass_Matched_Sum,"","#it{M}_{inv, Mother} (GeV/#it{c}^{2})","#frac{1}{N} d#it{N}/d#it{M}_{inv}",1,1.1,
                                    labelData, colorCompare, kTRUE, 5, 5, kTRUE,
                                    0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/InvMass_matched_afterQA.%s", outputDir.Data(), suffix.Data()), kFALSE , kTRUE);

        DrawPeriodQACompareHistoRatioTH1(canvas,0.11, 0.02, 0.05, 0.11,kFALSE,kTRUE,kFALSE,
                                        vecESD_InvMass_Matched_Sum,"","#it{M}_{inv, Mother} (GeV/#it{c}^{2})","#frac{1}{N} d#it{N}/d#it{M}_{inv}",1,1.1,
                                        labelData, colorCompare, kTRUE, 1.1, 1.1, kTRUE,
                                        0.95,0.92,0.03,fCollisionSystem,plotDataSets,fTrigger[0],31);
        SaveCanvas(canvas, Form("%s/Comparison/Ratios/ratio_InvMass_matched_afterQA.%s", outputDir.Data(), suffix.Data()));
        DeleteVecTH1D(vecESD_InvMass_Matched_Sum);

        //---------------------------------------------------------------------------------------------------------------
        //
        for(Int_t i=1; i<nSets; i++){
            if(isTrueContainer[i]){
                Int_t nPi0Rebin = 5;
                Int_t nEtaRebin = 5;

                TH2D* TruePi0_InvMass_Pt = vecESDTrueMesonCaloConvertedPhoton[0].at(i);
                TH2D* TruePi0_InvMass_Pt_matched = vecESDTrueMesonCaloConvertedPhotonMatched[0].at(i);
                TH2D* TrueEta_InvMass_Pt = vecESDTrueMesonCaloConvertedPhoton[1].at(i);
                TH2D* TrueEta_InvMass_Pt_matched = vecESDTrueMesonCaloConvertedPhotonMatched[1].at(i);

                if(TruePi0_InvMass_Pt && TruePi0_InvMass_Pt_matched && TrueEta_InvMass_Pt && TrueEta_InvMass_Pt_matched){
                TH1D* TruePi0_InvMass = TruePi0_InvMass_Pt->ProjectionX("TruePi0_InvMass",1,TruePi0_InvMass_Pt->GetNbinsY());
                TH1D* TruePi0_InvMass_matched = TruePi0_InvMass_Pt_matched->ProjectionX("TruePi0_InvMass_matched",1,TruePi0_InvMass_Pt_matched->GetNbinsY());
                TH1D* TrueEta_InvMass = TrueEta_InvMass_Pt->ProjectionX("TrueEta_InvMass",1,TrueEta_InvMass_Pt->GetNbinsY());
                TH1D* TrueEta_InvMass_matched = TrueEta_InvMass_Pt_matched->ProjectionX("TrueEta_InvMass_matched",1,TrueEta_InvMass_Pt_matched->GetNbinsY());

                TruePi0_InvMass->Rebin(nPi0Rebin);
                TruePi0_InvMass_matched->Rebin(nPi0Rebin);
                TrueEta_InvMass->Rebin(nEtaRebin);
                TrueEta_InvMass_matched->Rebin(nEtaRebin);

                TruePi0_InvMass->Sumw2();
                TruePi0_InvMass_matched->Sumw2();
                TrueEta_InvMass->Sumw2();
                TrueEta_InvMass_matched->Sumw2();

                TH1D* TruePi0_InvMass_sum = (TH1D*) TruePi0_InvMass->Clone();
                TH1D* TrueEta_InvMass_sum = (TH1D*) TrueEta_InvMass->Clone();

                TruePi0_InvMass_sum->Add(TruePi0_InvMass_matched,1);
                TrueEta_InvMass_sum->Add(TrueEta_InvMass_matched,1);
                TruePi0_InvMass_sum->Sumw2();
                TrueEta_InvMass_sum->Sumw2();
                //---------
                vecRange.push_back(TruePi0_InvMass_sum);
                vecRange.push_back(TruePi0_InvMass_matched);
                vecRange.push_back(TruePi0_InvMass);
                AdjustHistRange(vecRange,1.,1.1,kTRUE,1,0.);
                DrawAutoGammaHistMatch3H(TruePi0_InvMass_sum, TruePi0_InvMass_matched, TruePi0_InvMass,
                                        "",
                                        "#it{M}_{#gamma e^{+/-}} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma e^{+/-}}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,
                                        1,1.2, "Sum", "Rejected", "Accepted");
                PutProcessLabelAndEnergyOnPlot(0.94, 0.78, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTextMeasurement.Data(), 42, 0.03, "", 1, 1.25, 31);
                PutProcessLabelAndEnergyOnPlot(0.94, 0.665, 0.03, "#gamma's rec. with PCM and", Form("%s, e^{+/-} leads cluster",calo.Data()), "", 42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/TrackMatching_TrueMCPi0_IntegratedPt_%s.%s", outputDir.Data(), DataSets[i].Data(), suffix.Data()));

                delete TruePi0_InvMass;
                delete TruePi0_InvMass_matched;
                delete TruePi0_InvMass_sum;
                vecRange.clear();
                //---------------------------------------------------------------------------------------------------------------
                vecRange.push_back(TrueEta_InvMass_sum);
                vecRange.push_back(TrueEta_InvMass_matched);
                vecRange.push_back(TrueEta_InvMass);
                AdjustHistRange(vecRange,1.,1.1,kTRUE,1,0.);
                DrawAutoGammaHistMatch3H(TrueEta_InvMass_sum, TrueEta_InvMass_matched, TrueEta_InvMass,
                                        "",
                                        "#it{M}_{#gamma e^{+/-}} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma e^{+/-}}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,
                                        1,1.2, "Sum", "Rejected", "Accepted");
                PutProcessLabelAndEnergyOnPlot(0.94, 0.78, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTextMeasurement.Data(), 42, 0.03, "", 1, 1.25, 31);
                PutProcessLabelAndEnergyOnPlot(0.94, 0.665, 0.03, "#gamma's rec. with PCM and", Form("%s, e^{+/-} leads cluster",calo.Data()), "", 42, 0.03, "", 1, 1.25, 31);
                SaveCanvas(canvas, Form("%s/TrackMatching_TrueMCEta_IntegratedPt_%s.%s", outputDir.Data(), DataSets[i].Data(),suffix.Data()));

                delete TrueEta_InvMass;
                delete TrueEta_InvMass_matched;
                delete TrueEta_InvMass_sum;
                vecRange.clear();
                //---------------------------------------------------------------------------------------------------------------
                const Int_t nPtExample          = 2;
                Int_t ptExampleMin[nPtExample]  = {10,40};
                Int_t ptExampleMax[nPtExample]  = {12,50};
                for(Int_t iPt=0; iPt<nPtExample; iPt++){
                    Double_t ptMin = TruePi0_InvMass_Pt->GetYaxis()->GetBinUpEdge(ptExampleMin[iPt]);
                    Double_t ptMax = TruePi0_InvMass_Pt->GetYaxis()->GetBinUpEdge(ptExampleMax[iPt]);

                    TH1D* ExampleTruePi0_InvMass = TruePi0_InvMass_Pt->ProjectionX(Form("ExampleTruePi0_InvMass%01i",ptExampleMin[iPt]),ptExampleMin[iPt],ptExampleMax[iPt]);
                    TH1D* ExampleTruePi0_InvMass_matched = TruePi0_InvMass_Pt_matched->ProjectionX(Form("ExampleTruePi0_InvMass_matched%01i",ptExampleMin[iPt]),ptExampleMin[iPt],ptExampleMax[iPt]);
                    TH1D* ExampleTrueEta_InvMass = TrueEta_InvMass_Pt->ProjectionX(Form("ExampleTrueEta_InvMass%01i",ptExampleMin[iPt]),ptExampleMin[iPt],ptExampleMax[iPt]);
                    TH1D* ExampleTrueEta_InvMass_matched = TrueEta_InvMass_Pt_matched->ProjectionX(Form("ExampleTrueEta_InvMass_matched%01i",ptExampleMin[iPt]),ptExampleMin[iPt],ptExampleMax[iPt]);

                    ExampleTruePi0_InvMass->Rebin(nPi0Rebin);
                    ExampleTruePi0_InvMass_matched->Rebin(nPi0Rebin);
                    ExampleTrueEta_InvMass->Rebin(nEtaRebin);
                    ExampleTrueEta_InvMass_matched->Rebin(nEtaRebin);

                    ExampleTruePi0_InvMass->Sumw2();
                    ExampleTruePi0_InvMass_matched->Sumw2();
                    ExampleTrueEta_InvMass->Sumw2();
                    ExampleTrueEta_InvMass_matched->Sumw2();

                    TH1D* ExampleTruePi0_InvMass_sum = (TH1D*) ExampleTruePi0_InvMass->Clone();
                    TH1D* ExampleTrueEta_InvMass_sum = (TH1D*) ExampleTrueEta_InvMass->Clone();

                    ExampleTruePi0_InvMass_sum->Add(ExampleTruePi0_InvMass_matched,1);
                    ExampleTrueEta_InvMass_sum->Add(ExampleTrueEta_InvMass_matched,1);
                    ExampleTruePi0_InvMass_sum->Sumw2();
                    ExampleTrueEta_InvMass_sum->Sumw2();
                    //---------
                    vecRange.push_back(ExampleTruePi0_InvMass_sum);
                    vecRange.push_back(ExampleTruePi0_InvMass_matched);
                    vecRange.push_back(ExampleTruePi0_InvMass);
                    AdjustHistRange(vecRange,1.,1.1,kTRUE,1,0.);
                    DrawAutoGammaHistMatch3H(ExampleTruePi0_InvMass_sum, ExampleTruePi0_InvMass_matched, ExampleTruePi0_InvMass,
                                            "",
                                            "#it{M}_{#gamma e^{+/-}} (GeV/#it{c}^{2})",
                                            "d#it{N}/d#it{M}_{#gamma e^{+/-}}",
                                            0,0,0,
                                            0,0,0,
                                            1,0,0.8,
                                            1,1.2, "Sum", "Rejected", "Accepted");
                    PutProcessLabelAndEnergyOnPlot(0.94, 0.78, 0.03, fCollisionSystem.Data(), Form("%s, %s",plotDataSets[i].Data(),fTextMeasurement.Data()), Form("%.01f < #it{p}_{T} < %.01f GeV/#it{c}",ptMin, ptMax), 42, 0.03, "", 1, 1.25, 31);
                    PutProcessLabelAndEnergyOnPlot(0.94, 0.665, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "e^{+/-} leads cluster", "", 42, 0.03, "", 1, 1.25, 31);
                    SaveCanvas(canvas, Form("%s/TrackMatching_TrueMCPi0_PtBins%01i-%01i_%s.%s", outputDir.Data(), ptExampleMin[iPt], ptExampleMax[iPt], DataSets[i].Data(),suffix.Data()));

                    delete ExampleTruePi0_InvMass;
                    delete ExampleTruePi0_InvMass_matched;
                    delete ExampleTruePi0_InvMass_sum;
                    vecRange.clear();
                    //---------
                    vecRange.push_back(ExampleTrueEta_InvMass_sum);
                    vecRange.push_back(ExampleTrueEta_InvMass_matched);
                    vecRange.push_back(ExampleTrueEta_InvMass);
                    AdjustHistRange(vecRange,1.,1.1,kTRUE,1,0.);
                    DrawAutoGammaHistMatch3H(ExampleTrueEta_InvMass_sum, ExampleTrueEta_InvMass_matched, ExampleTrueEta_InvMass,
                                            "",
                                            "#it{M}_{#gamma e^{+/-}} (GeV/#it{c}^{2})",
                                            "d#it{N}/d#it{M}_{#gamma e^{+/-}}",
                                            0,0,0,
                                            0,0,0,
                                            1,0,0.8,
                                            1,1.2, "Sum", "Rejected", "Accepted");
                    PutProcessLabelAndEnergyOnPlot(0.94, 0.78, 0.03, fCollisionSystem.Data(), Form("%s, %s",plotDataSets[i].Data(),fTextMeasurementEta.Data()),
                                                   Form("%.01f < #it{p}_{T} < %.01f GeV/#it{c}",ptMin, ptMax), 42, 0.03, "", 1, 1.25, 31);
                    PutProcessLabelAndEnergyOnPlot(0.94, 0.665, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "e^{+/-} leads cluster", "", 42, 0.03, "", 1, 1.25, 31);
                    SaveCanvas(canvas, Form("%s/TrackMatching_TrueMCEta_PtBins%01i-%01i_%s.%s", outputDir.Data(), ptExampleMin[iPt], ptExampleMax[iPt], DataSets[i].Data(),suffix.Data()));

                    delete ExampleTrueEta_InvMass;
                    delete ExampleTrueEta_InvMass_matched;
                    delete ExampleTrueEta_InvMass_sum;
                    vecRange.clear();
                }
                //---------
                TH2D* TruePi0Photon_InvMass_Pt = vecESDTrueMesonCaloPhoton[0].at(i);
                TH2D* TruePi0Electron_InvMass_Pt = vecESDTrueMesonCaloElectron[0].at(i);
                TH2D* TrueEtaPhoton_InvMass_Pt = vecESDTrueMesonCaloPhoton[1].at(i);
                TH2D* TrueEtaElectron_InvMass_Pt = vecESDTrueMesonCaloElectron[1].at(i);

                TH1D* TruePi0CaloPhoton_InvMass = (TH1D*) TruePi0Photon_InvMass_Pt->ProjectionX("TruePi0CaloPhoton_InvMass",1,TruePi0Photon_InvMass_Pt->GetNbinsY());
                TH1D* TruePi0CaloElectron_InvMass = (TH1D*) TruePi0Electron_InvMass_Pt->ProjectionX("TruePi0CaloElectron_InvMass",1,TruePi0Electron_InvMass_Pt->GetNbinsY());
                TH1D* TrueEtaCaloPhoton_InvMass = (TH1D*) TrueEtaPhoton_InvMass_Pt->ProjectionX("TrueEtaCaloPhoton_InvMass",1,TrueEtaPhoton_InvMass_Pt->GetNbinsY());
                TH1D* TrueEtaCaloElectron_InvMass = (TH1D*) TrueEtaElectron_InvMass_Pt->ProjectionX("TrueEtaCaloElectron_InvMass",1,TrueEtaElectron_InvMass_Pt->GetNbinsY());

                TruePi0CaloPhoton_InvMass->Rebin(nPi0Rebin);
                TruePi0CaloElectron_InvMass->Rebin(nPi0Rebin);
                TrueEtaCaloPhoton_InvMass->Rebin(nEtaRebin);
                TrueEtaCaloElectron_InvMass->Rebin(nEtaRebin);

                TruePi0CaloPhoton_InvMass->Sumw2();
                TruePi0CaloElectron_InvMass->Sumw2();
                TrueEtaCaloPhoton_InvMass->Sumw2();
                TrueEtaCaloElectron_InvMass->Sumw2();

                TH1D* TruePi0Calo_InvMass_sum = (TH1D*) TruePi0CaloPhoton_InvMass->Clone();
                TH1D* TrueEtaCalo_InvMass_sum = (TH1D*) TrueEtaCaloPhoton_InvMass->Clone();

                TruePi0Calo_InvMass_sum->Add(TruePi0CaloElectron_InvMass,1);
                TrueEtaCalo_InvMass_sum->Add(TrueEtaCaloElectron_InvMass,1);
                TruePi0Calo_InvMass_sum->Sumw2();
                TrueEtaCalo_InvMass_sum->Sumw2();
                //---------
                vecRange.push_back(TruePi0Calo_InvMass_sum);
                vecRange.push_back(TruePi0CaloPhoton_InvMass);
                vecRange.push_back(TruePi0CaloElectron_InvMass);
                AdjustHistRange(vecRange,1.,10.,kTRUE,1,0.5);
                DrawAutoGammaHistMatch3H(TruePi0Calo_InvMass_sum, TruePi0CaloPhoton_InvMass, TruePi0CaloElectron_InvMass,
                                        "",
                                        "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma#gamma}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,
                                        1,1.2,"Sum","TrueMC: #gamma_{EMCAL}","TrueMC: e^{+/-}_{EMCAL}");
                PutProcessLabelAndEnergyOnPlot(0.75, 0.72, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTextMeasurement.Data());
                PutProcessLabelAndEnergyOnPlot(0.75, 0.61, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "", "");
                SaveCanvas(canvas, Form("%s/TrueMCPi0_GammaGamma_GammaElectron_%s.%s", outputDir.Data(), DataSets[i].Data(),suffix.Data()), kFALSE, kTRUE);
                vecRange.clear();
                //---------
                vecRange.push_back(TrueEtaCalo_InvMass_sum);
                vecRange.push_back(TrueEtaCaloPhoton_InvMass);
                vecRange.push_back(TrueEtaCaloElectron_InvMass);
                AdjustHistRange(vecRange,1.,10.,kTRUE,1,0.5);
                DrawAutoGammaHistMatch3H(TrueEtaCalo_InvMass_sum, TrueEtaCaloPhoton_InvMass, TrueEtaCaloElectron_InvMass,
                                        "",
                                        "#it{M}_{#gamma#gamma}  (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma#gamma}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,
                                        1,1.2,"Sum","TrueMC: #gamma_{EMCAL}","TrueMC: e^{+/-}_{EMCAL}");
                PutProcessLabelAndEnergyOnPlot(0.2, 0.9, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), fTextMeasurementEta.Data());
                PutProcessLabelAndEnergyOnPlot(0.2, 0.79, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "", "");
                SaveCanvas(canvas, Form("%s/TrueMCEta_GammaGamma_GammaElectron_%s.%s", outputDir.Data(), DataSets[i].Data(),suffix.Data()), kFALSE, kTRUE);
                vecRange.clear();
                //---------
                vecRange.push_back(TruePi0Calo_InvMass_sum);
                vecRange.push_back(TruePi0CaloPhoton_InvMass);
                vecRange.push_back(TruePi0CaloElectron_InvMass);
                AdjustHistRange(vecRange,1.,10.,kTRUE,1,0.5);
                TruePi0Calo_InvMass_sum->Add(TrueEtaCalo_InvMass_sum,1);
                TruePi0CaloPhoton_InvMass->Add(TrueEtaCaloPhoton_InvMass,1);
                TruePi0CaloElectron_InvMass->Add(TrueEtaCaloElectron_InvMass,1);
                DrawAutoGammaHistMatch3H(TruePi0Calo_InvMass_sum, TruePi0CaloPhoton_InvMass, TruePi0CaloElectron_InvMass,
                                        "",
                                        "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})",
                                        "d#it{N}/d#it{M}_{#gamma#gamma}",
                                        0,0,0,
                                        0,0,0,
                                        1,0,0.8,
                                        1,1.2,"Sum","TrueMC: #gamma_{EMCAL}","TrueMC: e^{+/-}_{EMCAL}");
                PutProcessLabelAndEnergyOnPlot(0.4, 0.9, 0.03, fCollisionSystem.Data(), plotDataSets[i].Data(), Form("%s & %s", fTextMeasurement.Data(), fTextMeasurementEta.Data()));
                PutProcessLabelAndEnergyOnPlot(0.4, 0.79, 0.03, Form("#gamma's rec. with PCM, %s",calo.Data()), "", "");
                SaveCanvas(canvas, Form("%s/TrueMCPi0Eta_GammaGamma_GammaElectron_%s.%s", outputDir.Data(), DataSets[i].Data(),suffix.Data()), kFALSE, kTRUE);

                delete TruePi0CaloPhoton_InvMass;
                delete TruePi0CaloElectron_InvMass;
                delete TruePi0Calo_InvMass_sum;
                delete TrueEtaCaloPhoton_InvMass;
                delete TrueEtaCaloElectron_InvMass;
                delete TrueEtaCalo_InvMass_sum;
                vecRange.clear();
                //---------
                } else {cout << "INFO: Could not find 'ESD_True*CaloConvertedPhoton_InvMass_Pt' and/or 'ESD_TruePi0CaloConvertedPhotonMatched_InvMass_Pt'" << endl;}
            }
        }
    }

    DeleteVecTH2D(vecESDMother);
    DeleteVecTH2D(vecESDMotherMatched);
    for(Int_t i=0; i<2; i++){
        DeleteVecTH2D(vecESDTrueMeson[i]);
        DeleteVecTH2D(vecESDTrueMesonMatched[i]);
        DeleteVecTH2D(vecESDTrueMesonCaloConvertedPhotonMatched[i]);
        DeleteVecTH2D(vecESDTrueMesonCaloConvertedPhoton[i]);
        DeleteVecTH2D(vecESDTrueMesonCaloPhoton[i]);
        DeleteVecTH2D(vecESDTrueMesonCaloElectron[i]);
    }


    fOutput->Write();
    fOutput->Close();
    delete fOutput;

    //*****************************************************************************************************
    //*****************************************************************************************************
    //****************************** Create Output ROOT-File **********************************************
    //*****************************************************************************************************
    //*****************************************************************************************************

    if(doExtQA>1 && doCellQA){
        cellQA = cellQAData;
        std::vector<TH2D*> DataMCHists;
        std::vector<TH2D*> DataMCHistsTime;
        for(Int_t j=0; j<nSets; j++){
            TH2D* fHistDataMCCell = (TH2D*) vecCellEnergyForComparison.at(j)->Clone(Form("vecEnergyComp_%i",j));
            fHistDataMCCell->Sumw2();
            TH2D* fHistDataMCCellTime = (TH2D*) vecCellTimingForComparison.at(j)->Clone(Form("vecTimeComp_%i",j));
            fHistDataMCCellTime->Sumw2();
            DataMCHists.push_back(fHistDataMCCell);
            DataMCHistsTime.push_back(fHistDataMCCellTime);
        }

        vector<Int_t> allCellsBad;
        vector<Int_t> allCellsGood;

        CollectAndPlotBadCellCandidates(fLog, allCellsBad, cellQA->cellIDsEnergy,"Energy - Mean/Sigma");
        CollectAndPlotBadCellCandidates(fLog, allCellsBad, cellQA->cellIDsTime,"Time - Mean/Sigma");
        CollectAndPlotBadCellCandidates(fLog, allCellsBad, cellQA->cellIDsHotCells1D,"HotCells1D");
        CheckBadCellCandidatesVec(DataMCHistsTime, cellQA->cellIDsHotCellsTime1D,"HotCellsTime1D");
        CollectAndPlotBadCellCandidates(fLog, allCellsBad, cellQA->cellIDsHotCellsTime1D,"HotCellsTime1D");
        CollectAndPlotBadCellCandidates(fLog, allCellsBad, cellQA->cellIDsHotCells2D,"HotCells2D");
        CollectAndPlotBadCellCandidates(fLog, allCellsBad, cellQA->cellIDsMissing,"Missing MC-Data");

        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "GOOD cells:" << endl;
        Int_t emptyCells = 0;
        DetermineGoodCellCandidates(allCellsGood, cellQA, nCaloCells, DataMCHists, emptyCells);
        cout << "found " << allCellsGood.size() << " good cells and " << emptyCells << " empty cells" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "allCellsBad.size before sort and unique: " << allCellsBad.size() << ".";
        fLog << "allCellsBad.size before sort and unique: " << allCellsBad.size() << ".";
        Int_t iMaximumAmountOfBadCells=3000;
        if((Int_t)allCellsBad.size()>iMaximumAmountOfBadCells){
            cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            cout << "ERROR: allCellsBad.size() too big " << allCellsBad.size() << ", check cuts!, Maximum is "<<iMaximumAmountOfBadCells <<  endl;
            cout << "RETURNING..." << endl;
            fLog << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
            fLog << "ERROR: allCellsBad.size() too big " << allCellsBad.size() << ", check cuts!, Maximum is "<<iMaximumAmountOfBadCells << endl;
            fLog << "RETURNING..." << endl;
            return;
        }
        selection_sort(allCellsBad.begin(),allCellsBad.end());
        vector<Int_t>::iterator it;
        it = unique(allCellsBad.begin(),allCellsBad.end());
        allCellsBad.resize( distance(allCellsBad.begin(),it) );
        cout << "Finally " << allCellsBad.size() << " different cells found!" << endl;
        fLog << "Finally " << allCellsBad.size() << " different cells found!" << endl;

        fstream fBadCells;
        fBadCells.open(Form("%s/Cells/%s.log",outputDir.Data(),DataSets[0].Data()), ios::out);
        for(Int_t iC=0; iC<(Int_t)allCellsBad.size(); iC++){
            cout << allCellsBad.at(iC) << ", ";
            fLog << allCellsBad.at(iC) << ", ";
            fBadCells << allCellsBad.at(iC) << endl;
        }
        fBadCells.close();
        cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        canvas->SetLeftMargin(0.1);canvas->SetRightMargin(0.1);canvas->SetTopMargin(0.06);canvas->SetBottomMargin(0.1);
        if((Int_t)allCellsBad.size()>0){
            cout << "start plotting bad cells" << endl;
            PlotBadCellReasons(cellQA,allCellsBad,canvas,outputDir,suffix,fClusters,plotDataSets[0],DataSets[0],fCollisionSystem);
            PlotBadCellOverview(kTRUE,kFALSE,DataMCHists.at(0),allCellsBad,canvas,outputDir,suffix,fClusters,plotDataSets[0],DataSets[0],fCollisionSystem);
            PlotBadCellOverview(kFALSE,kFALSE,DataMCHistsTime.at(0),allCellsBad,canvas,outputDir,suffix,fClusters,plotDataSets[0],DataSets[0],fCollisionSystem);
            for(Int_t j=1; j<nSets; j++){
                PlotBadCellOverview(kTRUE,kTRUE,DataMCHists.at(j),allCellsBad,canvas,outputDir,suffix,fClusters,plotDataSets[j],DataSets[j],fCollisionSystem);
                PlotBadCellOverview(kFALSE,kTRUE,DataMCHistsTime.at(j),allCellsBad,canvas,outputDir,suffix,fClusters,plotDataSets[j],DataSets[j],fCollisionSystem);
            }
            // plot for each bad cell the energy distribution together with the MC one as well as the time distribution together with that of a good cell
            PlotBadCellComparisonVecBoth(DataMCHists,DataMCHistsTime,colorCompare,allCellsBad, allCellsGood, nEvents, canvas2NTEO, padL2NTEO, padR2NTEO, outputDir, suffix, fClusters, plotDataSets,
                                         fCollisionSystem, kFALSE, 50,100);
        }

        if ((Int_t)allCellsGood.size()>0){
            cout << "start plotting good cells" << endl;
            // plot for 50 good cells the energy distribution together with the MC one as well as the time distribution together with that of another good cell
            PlotBadCellComparisonVecBoth(DataMCHists,DataMCHistsTime,colorCompare,allCellsGood, allCellsGood, nEvents, canvas2NTEO, padL2NTEO, padR2NTEO, outputDir, suffix, fClusters, plotDataSets,
                                         fCollisionSystem, kTRUE, 50,100);
        }
        char* nameOutput = Form("%s/ClusterQA_%s.root",outputDirRootFiles.Data(),DataSets[0].Data());
        TFile* fOutput = new TFile(nameOutput,"UPDATE");
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "Reopening Output file: " << nameOutput << " to store BadCellCandidates" << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        fLog << "Reopening Output file: " << nameOutput << " to store BadCellCandidates" << endl;
        fLog << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        SaveBadCellCandidates(cellQA->cellIDsEnergy,"Energy");
        SaveBadCellCandidates(cellQA->cellIDsTime,"Time");
        SaveBadCellCandidates(cellQA->cellIDsHotCells1D,"HotCells1D");
        SaveBadCellCandidates(cellQA->cellIDsHotCellsTime1D,"HotCellsTime1D");
        SaveBadCellCandidates(cellQA->cellIDsHotCells2D,"HotCells2D");
        SaveBadCellCandidates(cellQA->cellIDsMissing,"Missing");
        SaveBadCellCandidates(allCellsBad,"allCellsBad");

        fOutput->Write();
        fOutput->Close();
        delete fOutput;

        for(Int_t j=0; j<nSets; j++){
        delete DataMCHists.at(j);
        delete DataMCHistsTime.at(j);
        }
        DataMCHists.clear();
        DataMCHistsTime.clear();
    }

    DeleteVecTH2D(vecCellEnergyForComparison);
    DeleteVecTH2D(vecCellTimingForComparison);

    fLog.close();

    delete line;
    delete fDeltaPt;
    delete fBinsAnalysisPt;
    if(cellQAData) delete cellQAData;
    delete cvsQuadratic;
    delete canvas;

    delete[] nEvents;
    delete[] isTrueContainer;

    delete[] fTrigger;
    delete[] fCutSelection;
    delete[] fEventCutSelection;
    delete[] fGammaCutSelection;
    delete[] fClusterCutSelection;
    delete[] fElectronCutSelection;
    delete[] fMesonCutSelection;

    TH1::AddDirectory(kTRUE);

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "Done with ClusterQA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    return;

}//end
